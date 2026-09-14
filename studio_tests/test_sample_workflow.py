import json

import pytest

from wmlstudio.export import export_results
from wmlstudio.project import Project
from wmlstudio.sample_workflow import (
    combined_feature_rows,
    current_hydra_evidence,
    feature_fields,
    hydra_evidence_status,
    link_hydra,
    set_cluster,
    suggest_hydra_links,
)


def add(project, path, name):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(">a\nACGT\n")
    return project.add_sample(path, name)


def report(*names):
    return {"import_provenance": {"sha256": "abc", "source_path": "hydra.json"},
            "samples": [{"sample": name, "summary": {"amr_genes": 1}, "hits": [
                {"gene": "blaKPC-2", "element_type": "AMR", "class": "BETA-LACTAM", "primary": True},
                {"gene": "secondary", "element_type": "AMR", "primary": False},
                {"gene": "blaKPC-2", "element_type": "AMR", "class": "BETA-LACTAM", "primary": True},
            ]} for name in names]}


def test_hydra_same_stems_require_explicit_disambiguation(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        first = add(project, tmp_path / "a" / "same.fasta", "same")
        second = add(project, tmp_path / "b" / "same.fasta", "same")
        upstream = report("same", "missing")
        suggestions = suggest_hydra_links(project, upstream)
        assert suggestions["mapping"] == {}
        assert set(suggestions["ambiguous"]["same"]) == {first, second}
        assert suggestions["missing"] == ["missing"]
        link_hydra(project, upstream, {"same": first})
        assert project.get_sample(first)["metadata"]["hydra"]["source_sample"] == "same"
        assert "hydra" not in project.get_sample(second)["metadata"]


def test_hydra_invalid_mapping_cannot_partially_attach_evidence(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        with pytest.raises(KeyError):
            link_hydra(project, report("a", "b"), {"a": sid, "b": "missing-id"})
        assert project.get_sample(sid)["metadata"] == {}
        with pytest.raises(ValueError, match="Multiple"):
            link_hydra(project, report("a", "b"), {"a": sid, "b": sid})


def test_hydra_execution_provenance_survives_reopen_and_bundle(tmp_path):
    from wmlstudio.library import export_bundle, import_bundle

    project_path = tmp_path / "study.wmlstudio"
    with Project(project_path) as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        upstream = report("a")
        upstream["execution_provenance"] = {
            "tool_version": "1.0.0", "command": ["hydra", "run"],
            "database_sha256": "database123", "inputs": [{"sha256": "input123"}],
        }
        link_hydra(project, upstream, {"a": sid})
        upstream["execution_provenance"]["inputs"][0]["sha256"] = "mutated"
    bundle = tmp_path / "portable.json"
    with Project(project_path) as project:
        assert project.get_sample(sid)["metadata"]["hydra"]["execution_provenance"]["inputs"][0]["sha256"] == "input123"
        export_bundle(project, bundle)
    with Project(tmp_path / "imported.wmlstudio") as project:
        import_bundle(project, bundle)
        evidence = project.get_sample(sid)["metadata"]["hydra"]
        assert evidence["execution_provenance"]["database_sha256"] == "database123"


@pytest.mark.parametrize("format", ["html", "json", "csv", "tsv"])
def test_selected_report_contains_only_selected_evidence_and_highlights(tmp_path, format):
    with Project(tmp_path / "study.wmlstudio") as project:
        first = add(project, tmp_path / "first.fasta", "selected_isolate")
        second = add(project, tmp_path / "second.fasta", "excluded_isolate")
        project.set_metadata(first, {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}})
        link_hydra(project, report("selected_isolate"), {"selected_isolate": first})
        set_cluster(project, [first], "Ward <A>", "#123456")
        features = combined_feature_rows(project.samples(), selected_ids=[first])
        assert features[0]["amr_genes"] == ["blaKPC-2"]
        assert features[0]["organism"] == "Klebsiella pneumoniae"
        assert features[0]["cluster_highlight"] is True
        output = tmp_path / f"selected.{format}"
        export_results(project.samples(), output, selected_ids=[first])
        text = output.read_text(encoding="utf-8-sig")
        assert "selected_isolate" in text and "excluded_isolate" not in text
        assert "blaKPC-2" in text
        if format == "html":
            assert "Ward &lt;A&gt;" in text and "border-left:5px solid #123456" in text
        elif format == "json":
            assert json.loads(text)["report_scope"] == {"mode": "selected", "sample_ids": [first]}
        assert second not in [row["sample_id"] for row in features]


def test_unknown_selected_id_is_error_and_cannot_replace_existing_report(tmp_path):
    path = tmp_path / "report.json"
    path.write_text("keep")
    with pytest.raises(KeyError):
        export_results([{"sample_id": "one", "sample_name": "one"}], path, selected_ids=["two"])
    assert path.read_text() == "keep"


def test_managed_originals_are_protected_even_when_record_is_not_selected(tmp_path):
    original = tmp_path / "original.fasta"
    original.write_text(">original\nACGT\n")
    records = [{"id": "one", "name": "one", "result": None, "input_path": "copy.fasta",
                "metadata": {"workflow": {"source_path": str(original)}}},
               {"id": "two", "name": "two", "result": None, "input_path": "other.fasta"}]
    with pytest.raises(ValueError, match="different output"):
        export_results(records, original, "html", selected_ids=["two"])
    assert original.read_text() == ">original\nACGT\n"


def test_hashed_hydra_becomes_archived_not_current_when_typing_input_changes(tmp_path):
    with Project(tmp_path / 'study.sqlite') as project:
        sid = add(project, tmp_path / 'assembly.fa', 'isolate')
        project.set_result(sid, {'input_sha256': 'a' * 64, 'st': '7'})
        upstream = report('isolate')
        upstream['execution_provenance'] = {'inputs': [{'sha256': 'a' * 64}]}
        link_hydra(project, upstream, {'isolate': sid})
        sample = project.get_sample(sid)
        assert hydra_evidence_status(sample)['status'] == 'current'
        assert feature_fields(sample)['amr_genes'] == ['blaKPC-2']
        project.set_result(sid, {'input_sha256': 'b' * 64, 'st': '8'})
        sample = project.get_sample(sid)
        assert hydra_evidence_status(sample)['status'] == 'stale'
        assert current_hydra_evidence(sample) == {}
        assert feature_fields(sample)['amr_genes'] == []
        assert feature_fields(sample)['hydra_amr_genes'] is None
        assert sample['metadata']['hydra']['hits'][0]['gene'] == 'blaKPC-2'
        assert sample['result']['st'] == '8'


def test_hashless_import_stays_unverified_even_with_matching_link_baseline(tmp_path):
    with Project(tmp_path / 'study.sqlite') as project:
        sid = add(project, tmp_path / 'assembly.fa', 'isolate')
        project.set_result(sid, {'input_sha256': 'a' * 64})
        link_hydra(project, report('isolate'), {'isolate': sid})
        sample = project.get_sample(sid)
        assert hydra_evidence_status(sample)['status'] == 'unverified'
        assert hydra_evidence_status(sample)['evidence_input_sha256'] is None
        assert sample['metadata']['hydra']['linked_input_sha256'] == 'a' * 64
        assert feature_fields(sample)['amr_genes'] == ['blaKPC-2']
        project.set_result(sid, {'input_sha256': 'b' * 64})
        assert hydra_evidence_status(project.get_sample(sid))['status'] == 'stale'


def test_multi_sample_report_does_not_assign_input_hashes_by_unverified_list_order(tmp_path):
    with Project(tmp_path / 'study.sqlite') as project:
        ids = [add(project, tmp_path / f'{name}.fa', name) for name in ('one', 'two')]
        upstream = report('one', 'two')
        upstream['execution_provenance'] = {'inputs': [{'sha256': 'a' * 64}, {'sha256': 'b' * 64}]}
        for sid in ids:
            project.set_result(sid, {'input_sha256': 'a' * 64})
        link_hydra(project, upstream, dict(zip(('one', 'two'), ids)))
        assert all(hydra_evidence_status(project.get_sample(sid))['status'] == 'unverified' for sid in ids)


def test_current_assembly_hash_is_used_before_primary_typing_exists_and_old_amr_is_retained(tmp_path):
    with Project(tmp_path / 'study.sqlite') as project:
        sid = add(project, tmp_path / 'assembly.fa', 'isolate')
        project.update_metadata(sid, {'assembly': {'provenance': {'assembly_sha256': 'a' * 64}}})
        upstream = report('isolate')
        upstream['execution_provenance'] = {'inputs': [{'sha256': 'a' * 64}]}
        link_hydra(project, upstream, {'isolate': sid})
        # The assembly hash agrees, but its typing job is still queued.
        assert hydra_evidence_status(project.get_sample(sid))['status'] == 'unverified'
        replacement = report('isolate')
        replacement['import_provenance']['sha256'] = 'new-report'
        link_hydra(project, replacement, {'isolate': sid})
        archived = [event for event in project.history(sid) if event['action'] == 'hydra_superseded']
        assert len(archived) == 1
        assert archived[0]['details']['evidence']['execution_provenance']['inputs'][0]['sha256'] == 'a' * 64


def test_failed_rerun_does_not_upgrade_saved_hash_match_to_current_evidence(tmp_path):
    with Project(tmp_path / 'study.sqlite') as project:
        sid = add(project, tmp_path / 'assembly.fa', 'isolate')
        project.set_result(sid, {'input_sha256': 'a' * 64})
        upstream = report('isolate')
        upstream['execution_provenance'] = {'inputs': [{'sha256': 'a' * 64}]}
        link_hydra(project, upstream, {'isolate': sid})
        project.set_status(sid, 'failed', 'Current typing attempt failed')
        assert hydra_evidence_status(project.get_sample(sid))['status'] == 'unverified'


def test_organism_evidence_does_not_displace_the_hydra_provenance_fields(tmp_path):
    """Both records reach the same row, so neither may overwrite the other's fields.

    The organism decision and the HYDRA evidence are separate provenance: one says
    how a label was chosen, the other says which report an AMR gene came from and
    against which input hash. A row that quietly dropped the second would make an
    unverifiable AMR call look like a sourced one.
    """
    with Project(tmp_path / 'study.sqlite') as project:
        sid = add(project, tmp_path / 'assembly.fa', 'isolate')
        project.set_result(sid, {'input_sha256': 'a' * 64, 'st': '7'})
        upstream = report('isolate')
        upstream['execution_provenance'] = {'inputs': [{'sha256': 'a' * 64}]}
        link_hydra(project, upstream, {'isolate': sid})
        project.update_metadata(sid, {'organism_evidence': {
            'basis': 'reference_ani', 'confidence': 'strong', 'status': 'resolved'}})
        fields = feature_fields(project.get_sample(sid))
        assert fields['organism_basis'] == 'reference_ani'
        assert fields['organism_confidence'] == 'strong'
        assert fields['organism_status'] == 'resolved'
        assert fields['hydra_source_sample'] == 'isolate'
        assert fields['hydra_report_sha256'] == 'abc'
        assert fields['hydra_amr_genes'] == 1


def test_organism_typing_column_is_empty_when_the_assembly_changed(tmp_path):
    """Typing from an earlier assembly must not be exported as this one's.

    The flat column is a convenience for someone reading a spreadsheet rather
    than the drill-down, which is exactly the reader least able to notice that
    a call belonged to a different input.
    """
    from wmlstudio.characterization import persist_characterization
    from wmlstudio.sequence import file_sha256
    with Project(tmp_path / 'study.sqlite') as project:
        source = tmp_path / 'assembly.fa'
        sid = add(project, source, 'isolate')
        digest = file_sha256(source)
        project.set_result(sid, {'input_sha256': digest, 'st': '7'})
        persist_characterization(project, sid, {
            'format_version': 1, 'input_sha256': digest, 'input_path': str(source),
            'sccmec': {'status': 'typed', 'candidate_types': ['IV'],
                       'summary': 'SCCmec candidate type IV'}})
        current = feature_fields(project.get_sample(sid))['organism_typing']
        # The module reports its own conservative wording; the column carries it
        # verbatim rather than restating a type the assay declined to assert.
        assert current and 'withheld' in current
        project.set_result(sid, {'input_sha256': 'b' * 64, 'st': '8'})
        assert feature_fields(project.get_sample(sid))['organism_typing'] == ''


def test_a_core_genome_run_never_empties_the_classical_st_column(tmp_path):
    from wmlstudio.sample_workflow import typing_profiles
    mlst = {"scheme": "Klebsiella MLST", "scheme_digest": "mlst-v1", "st": "258",
            "status": "complete", "alleles": {f"gene{i}": "1" for i in range(7)}}
    cgmlst = {"scheme": "Klebsiella cgMLST", "scheme_digest": "cg-v1", "st": None,
              "status": "complete", "parameters": {"method": "full-cds-cgmlst-v2"},
              "alleles": {f"locus{i:04d}": ("3" if i % 2 else None) for i in range(400)}}
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        project.set_result(sid, mlst)
        project.set_result(sid, cgmlst)
        record = dict(project.get_sample(sid), analyses=project.analysis_results(sid))

        profiles = typing_profiles(record)
        assert profiles["mlst"]["st"] == "258" and profiles["cgmlst"]["scheme_digest"] == "cg-v1"
        fields = feature_fields(record)
        assert fields["mlst_st"] == "258" and fields["mlst_scheme"] == "Klebsiella MLST"
        assert fields["cgmlst_scheme"] == "Klebsiella cgMLST"
        assert fields["typing_kinds"] == ["mlst", "cgmlst"]
        # The scheme size and the number of loci actually called stay separate.
        assert fields["cgmlst_loci"] == 400 and fields["cgmlst_called"] == 200


def test_a_sample_with_only_an_st_reports_no_core_genome_profile(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        project.set_result(sid, {"scheme": "MLST", "scheme_digest": "mlst-v1", "st": "11",
                                 "status": "complete", "alleles": {"adk": "1"}})
        fields = feature_fields(dict(project.get_sample(sid),
                                     analyses=project.analysis_results(sid)))
        assert fields["mlst_st"] == "11" and fields["typing_kinds"] == ["mlst"]
        assert fields["cgmlst_scheme"] == "" and fields["cgmlst_loci"] == 0
        assert fields["cgmlst_called"] == 0
