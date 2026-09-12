import json

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.library import (
    Library,
    export_bundle,
    export_profile_table,
    import_bundle,
    import_profile_table,
    parse_profile_table,
)
from wmlstudio.project import Project


def table_file(tmp_path, text):
    path = tmp_path / "profiles.tsv"
    path.write_text(text)
    return path


@pytest.mark.parametrize("token", ["LNF", "PLOT3", "NIPH", "NIPHEM", "PAMA", "ASM", "ALM", "?", "NA", "000", "-0", "+0.0", "INF", "١٢", "１２", "new-tool-token"])
def test_external_missing_tokens_never_become_shared_alleles(token):
    parsed = parse_profile_table(f"FILE\tl1\tl2\nA\t1\t{token}\nB\t1\t{token}\n", "scheme", "verified")
    assert parsed["results"][0]["alleles"]["l2"] is None
    pair = pairwise_distances(parsed["results"])[0]
    assert pair["shared_loci"] == 1 and pair["total_loci"] == 2
    assert pair["distance"] is None


def test_numeric_padding_is_normalized_and_novel_identity_preserved_separately():
    parsed = parse_profile_table("FILE\ta\tb\tc\td\nA\t007\tINF-12\t6~ABCDEF\t12~\n", "scheme")
    result = parsed["results"][0]
    assert result["alleles"] == {"a": "7", "b": None, "c": None, "d": None}
    assert result["external_alleles"] == {"a": "7", "b": "12~", "c": "6~abcdef", "d": "12~"}
    assert any("inferred/novel" in warning for warning in parsed["warnings"])


@pytest.mark.parametrize("text", [
    "FILE\ta\tb\nA\t1\t2\nB\t1\n",
    "FILE\ta\nA\t1\nA\t2\n",
    "FILE\ta\ta\nA\t1\t2\n",
    "FILE\t\nA\t1\n",
])
def test_malformed_or_ambiguous_table_is_rejected_before_any_write(tmp_path, text):
    with Project(tmp_path / "study.wmlstudio") as project:
        with pytest.raises(ValueError):
            import_profile_table(project, table_file(tmp_path, text), "scheme")
        assert project.samples() == []


def test_repeated_column_merges_matching_calls_and_keeps_full_scheme_denominator():
    parsed = parse_profile_table("FILE\ta\ta\nb\t1\tLNF\n", "scheme", "verified", ["a", "b", "c"])
    assert parsed["results"][0]["alleles"] == {"a": "1", "b": None, "c": None}
    assert len(parsed["loci"]) == 3
    assert parsed["warnings"]


def test_table_roundtrip_preserves_novel_markers_without_manufacturing_fasta(tmp_path):
    source = table_file(tmp_path, "FILE\tl1\tl2\tl3\nA\t7\t6~abcdef\t0\nB\t2\t12~\t3\n")
    first_output, second_output = tmp_path / "first.tsv", tmp_path / "second.tsv"
    with Project(tmp_path / "first.wmlstudio") as project:
        imported = import_profile_table(project, source, "scheme", "verified")
        assert imported["imported"] == 2
        assert all(sample["profile_only"] for sample in project.samples())
        export_profile_table(project, first_output)
    with Project(tmp_path / "second.wmlstudio") as project:
        import_profile_table(project, first_output, "scheme", "verified")
        export_profile_table(project, second_output)
    assert first_output.read_bytes() == second_output.read_bytes()
    assert "6~abcdef" in second_output.read_text(encoding="utf-8-sig")


def test_unverified_tables_do_not_claim_shared_scheme_revision():
    a = parse_profile_table("FILE\tl1\nA\t1\n", "same-name")["results"][0]
    b = parse_profile_table("FILE\tl1\nB\t1\n", "same-name")["results"][0]
    assert a["scheme_digest"] != b["scheme_digest"]
    assert pairwise_distances([a, b])[0]["distance"] is None


def test_portable_bundle_keeps_metadata_collections_secondary_profiles_and_ids(tmp_path):
    bundle = tmp_path / "evidence.json"
    with Project(tmp_path / "original.wmlstudio") as project:
        sid = project.add_profile("Isolate", {"scheme": "MLST", "scheme_digest": "mlst", "st": "131", "alleles": {"a": "1"}},
                                  {"organism": {"genus": "Escherichia", "species": "coli"}, "ward": "ICU"})
        project.set_analysis(sid, {"scheme": "cgMLST", "scheme_digest": "cgmlst", "st": None, "alleles": {"b": "7"}})
        collection = project.create_collection("Investigation")
        project.set_collection_members(collection, [sid])
        export_bundle(project, bundle)
    assert json.loads(bundle.read_text())["contains_sequences"] is False
    with Project(tmp_path / "imported.wmlstudio") as project:
        imported = import_bundle(project, bundle)
        assert imported["sample_ids"] == [sid]
        restored = project.get_sample(sid)
        assert restored["profile_only"] and not restored["missing_input"]
        assert restored["metadata"]["ward"] == "ICU"
        assert restored["result"]["st"] == "131"
        assert len(project.analysis_results(sid)) == 2
        assert project.collections()[0]["sample_ids"] == [sid]
        repeated = import_bundle(project, bundle)
        assert repeated["sample_ids"][0] != sid
        assert len(project.samples()) == 2
        assert project.get_sample(sid)["result"]["st"] == "131"


def test_bad_bundle_collection_rolls_back_all_samples(tmp_path):
    source = tmp_path / "bad.json"
    source.write_text(json.dumps({"format": "WMLSTudio-evidence-bundle", "format_version": 1,
                                  "samples": [{"id": "id", "name": "Sample", "result": {"alleles": {}}}],
                                  "collections": [{"name": "Broken", "sample_ids": ["absent"]}]}))
    with Project(tmp_path / "study.wmlstudio") as project:
        with pytest.raises(ValueError, match="outside"):
            import_bundle(project, source)
        assert project.samples() == [] and project.collections() == []


def test_cross_project_library_remains_searchable_after_reopen_and_source_removal(tmp_path):
    database = tmp_path / "library.sqlite"
    with Library(database) as library:
        for number in range(2):
            with Project(tmp_path / f"run{number}.wmlstudio") as project:
                project.add_profile(f"isolate-{number}", {"scheme": "MLST", "scheme_digest": "mlst", "st": str(7 + number), "alleles": {"a": "1"}},
                                    {"organism": {"genus": "Listeria", "species": "monocytogenes"},
                                     "hydra": {"hits": [{"gene": "fosX", "element_type": "AMR", "primary": True}]}})
                assert library.index_project(project) == 1
    (tmp_path / "run0.wmlstudio").unlink()
    with Library(database) as library:
        assert len(library.search(genus="listeria", amr_gene="FOSX")) == 2
        assert len(library.search(st=7)) == 1
        assert len(library.search(q="isolate-1", scheme_digest="mlst")) == 1
        assert library.search(species="wrong") == []


@pytest.mark.parametrize('prefix', ['NOVEL_', 'SHA256_'])
def test_full_sequence_id_table_roundtrip_retains_token_without_inventing_validation(tmp_path, prefix):
    token = prefix + 'A' * 64
    source = table_file(tmp_path, 'FILE\tlocus\nA\t' + token + '\nB\t' + token + '\n')
    exported = tmp_path / 'roundtrip.tsv'
    with Project(tmp_path / 'study.sqlite') as project:
        import_profile_table(project, source, 'scheme', 'verified')
        records = [sample['result'] for sample in project.samples()]
        assert records[0]['alleles']['locus'] is None
        assert records[0]['external_alleles']['locus'] == prefix + 'a' * 64
        assert records[0]['calls'][0]['status'] == 'external_sequence_unverified'
        assert pairwise_distances(records)[0]['comparable'] is False
        export_profile_table(project, exported)
    parsed = parse_profile_table(exported.read_text(encoding='utf-8-sig'), 'scheme', 'verified')
    assert parsed['results'][0]['external_alleles']['locus'] == prefix + 'a' * 64
    assert 'NOVEL_' not in str(parsed['results'][0]['alleles'])


@pytest.mark.parametrize('status', ['failed', 'running', 'interrupted', 'queued'])
def test_library_reuse_refuses_noncompleted_current_state(tmp_path, status):
    from wmlstudio.library_dialog import reuse_library_profiles
    with Project(tmp_path / 'source.sqlite') as source, Project(tmp_path / 'target.sqlite') as target:
        sid = source.add_profile('isolate', {'scheme': 's', 'scheme_digest': 'd', 'alleles': {'a': '1'}})
        source.set_status(sid, status)
        snapshot = {**source.get_sample(sid), 'project_path': str(source.path),
                    'analyses': source.analysis_results(sid)}
        with pytest.raises(ValueError, match='not a completed current sample'):
            reuse_library_profiles(target, [snapshot])
        assert target.samples() == []


def test_library_reuse_refuses_old_secondary_after_assembly_replacement(tmp_path):
    from wmlstudio.library_dialog import reuse_library_profiles
    with Project(tmp_path / 'target.sqlite') as target:
        snapshot = {'id': 'source-id', 'name': 'isolate', 'project_path': '/old/project.sqlite',
                    'status': 'completed', 'result': None,
                    'metadata': {'assembly': {'provenance': {'assembly_sha256': 'b' * 64}}},
                    'analyses': [{'scheme': 's', 'scheme_digest': 'd', 'input_sha256': 'a' * 64,
                                  'alleles': {'a': '1'}}]}
        with pytest.raises(ValueError, match='no saved allelic profile'):
            reuse_library_profiles(target, [snapshot])
        assert target.samples() == []
