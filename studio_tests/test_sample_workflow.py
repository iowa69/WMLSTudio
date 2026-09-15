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


MLST_RESULT = {"scheme": "Klebsiella MLST", "scheme_digest": "mlst-v1", "st": "258",
               "status": "complete", "input_sha256": "a" * 64,
               "alleles": {f"gene{index}": "1" for index in range(7)},
               "scheme_metadata": {"type": "MLST", "last_updated": "2026-01-05"},
               "parameters": {"method": "exact-nucleotide", "strands": "both"}}
CG_RESULT = {"scheme": "Klebsiella cgMLST", "scheme_digest": "cg-v1", "st": None,
             "status": "complete", "input_sha256": "a" * 64,
             "alleles": {f"locus{index:04d}": "3" for index in range(400)},
             "parameters": {"method": "full-cds-cgmlst-v2", "genetic_code": 11}}


def test_the_whole_sample_configuration_is_read_from_one_place(tmp_path):
    """Every menu reads these names, so no surface has to keep a second copy.

    A project written before the single home stored some of the same choices at the
    top of metadata; those are still read, and the workflow record wins when both
    are present, so an old project keeps working without being rewritten.
    """
    from wmlstudio.sample_workflow import sample_configuration
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        project.set_metadata(sid, {
            "organism": {"genus": "Klebsiella", "species": "pneumoniae"},
            "typing_mode": "auto", "scheme_path": "/legacy/scheme",
            "workflow": {"typing_mode": "manual", "scheme_path": "/schemes/kpneumoniae",
                         "calling_mode": "full_cds", "genetic_code": 11,
                         "cgmlst_scheme_key": "cgmlst.org:kpneumoniae-2358",
                         "cgmlst_scheme_path": "/cgmlst/Klebsiella", "run_cgmlst": True}})
        configuration = sample_configuration(project.get_sample(sid))
        assert configuration["genus"] == "Klebsiella" and configuration["species"] == "pneumoniae"
        assert configuration["typing_mode"] == "manual"
        assert configuration["scheme_path"] == "/schemes/kpneumoniae"
        assert configuration["calling_mode"] == "full_cds" and configuration["genetic_code"] == 11
        assert configuration["cgmlst_scheme_key"] == "cgmlst.org:kpneumoniae-2358"
        assert configuration["run_cgmlst"] is True and configuration["run_hydra"] is False

    # A sample nothing has configured reports the absence rather than a guess.
    blank = sample_configuration({"id": "x", "name": "x"})
    assert blank["typing_mode"] == "" and blank["scheme_path"] is None
    assert blank["cgmlst_scheme_key"] is None and blank["run_cgmlst"] is False
    assert blank["organism_source"] == "unknown" and blank["set_by"] == {}


def test_an_mlst_scheme_match_fills_the_organism_in_but_never_as_taxonomy(tmp_path):
    """The lineage a panel matched is corroborating evidence, and says so wherever it goes."""
    from wmlstudio.sample_workflow import MLST_ORGANISM_CAVEAT, sample_organism
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        project.set_result(sid, {**MLST_RESULT, "identification": {
            "identification_status": "assigned",
            "organism": {"genus": "Klebsiella", "species": "pneumoniae"}}})
        identity = sample_organism(project.get_sample(sid))
        assert identity["organism"] == "Klebsiella pneumoniae"
        assert identity["source"] == "scheme_match" and identity["basis"] == "mlst_panel"
        assert identity["confidence"] == "panel_compatibility"
        assert identity["note"] == MLST_ORGANISM_CAVEAT
        # The flat row calls it a detection, never an assignment somebody made.
        assert feature_fields(project.get_sample(sid))["organism_source"] == "local_scheme_detection"

        project.update_metadata(sid, {"organism": {"genus": "Klebsiella", "species": "variicola"}})
        assigned = sample_organism(project.get_sample(sid))
        assert assigned["species"] == "variicola" and assigned["source"] == "assigned"
        assert feature_fields(project.get_sample(sid))["organism_source"] == "assigned"


def test_typing_is_not_repeated_when_nothing_it_was_run_against_changed():
    from wmlstudio.sample_workflow import rerun_decision
    request = {"input_sha256": "a" * 64, "scheme_digest": "mlst-v1", "scheme": "Klebsiella MLST",
               "scheme_version": "2026-01-05", "caller": {"method": "exact-nucleotide"}}
    decision = rerun_decision(MLST_RESULT, request)
    assert decision["rerun"] is False and decision["status"] == "current"
    assert decision["changes"] == [] and decision["unverified"] == []
    assert set(decision["compared"]) == {"input_sha256", "scheme_digest", "scheme",
                                         "scheme_version", "caller.method"}
    assert "stands" in decision["summary"]


def test_a_rerun_names_which_of_the_four_things_changed():
    from wmlstudio.sample_workflow import rerun_decision
    request = {"input_sha256": "a" * 64, "scheme_digest": "mlst-v1", "scheme": "Klebsiella MLST",
               "scheme_version": "2026-01-05", "caller": {"method": "exact-nucleotide"}}
    for field, value, expected in [("input_sha256", "b" * 64, "the sequence input"),
                                   ("scheme_digest", "mlst-v2", "the scheme's contents"),
                                   ("scheme", "Other MLST", "the scheme"),
                                   ("scheme_version", "2026-06-01", "the scheme version")]:
        decision = rerun_decision(MLST_RESULT, {**request, field: value})
        assert decision["rerun"] is True and decision["status"] == "changed"
        assert [change["field"] for change in decision["changes"]] == [field]
        assert expected in decision["summary"]
    caller = rerun_decision(MLST_RESULT, {**request, "caller": {"method": "full-cds-cgmlst-v2"}})
    assert [change["field"] for change in caller["changes"]] == ["caller.method"]
    assert caller["changes"][0]["was"] == "exact-nucleotide"
    # Nothing stored at all is its own answer, not a change.
    assert rerun_decision(None, request)["status"] == "never_run"
    assert rerun_decision({"status": "qc_only", "scheme_digest": None}, request)["rerun"] is True


def test_a_result_that_cannot_be_shown_to_be_current_is_never_called_current():
    """Unverifiable is a third answer. It is not currency and it is not staleness."""
    from wmlstudio.sample_workflow import rerun_decision
    request = {"input_sha256": "a" * 64, "scheme_digest": "mlst-v1"}
    hashless = rerun_decision({**MLST_RESULT, "input_sha256": None}, request)
    assert hashless["rerun"] is True and hashless["status"] == "unverifiable"
    assert hashless["changes"] == []
    assert hashless["unverified"][0]["field"] == "input_sha256"
    assert "the stored result does not record it" in hashless["unverified"][0]["reason"]
    # A setting the request asks about that the stored result never recorded is the
    # same kind of gap: it is reported, not assumed to have matched.
    asked = rerun_decision(MLST_RESULT, {**request, "caller": {"min_identity": 0.9}})
    assert asked["status"] == "unverifiable"
    assert [entry["field"] for entry in asked["unverified"]] == ["caller.min_identity"]
    # A request that does not state the scheme version is not asking about it.
    assert rerun_decision(MLST_RESULT, request)["status"] == "current"


def test_each_typing_kind_answers_its_own_rerun_question(tmp_path):
    """A current cgMLST profile never excuses an MLST rerun, or the other way round."""
    from wmlstudio.sample_workflow import rerun_decision
    with Project(tmp_path / "study.wmlstudio") as project:
        sid = add(project, tmp_path / "a.fasta", "a")
        project.set_result(sid, MLST_RESULT)
        project.set_result(sid, CG_RESULT)
        classical = rerun_decision(project.latest_analysis(sid, "mlst"),
                                   {"input_sha256": "a" * 64, "scheme_digest": "mlst-v1"})
        core = rerun_decision(project.latest_analysis(sid, "cgmlst"),
                              {"input_sha256": "b" * 64, "scheme_digest": "cg-v1"})
        assert classical["status"] == "current" and core["status"] == "changed"
        assert core["changes"][0]["field"] == "input_sha256"


KPNEUMONIAE = {"key": "cgmlst.org:kpneumoniae-2358", "genus": "Klebsiella",
               "species": "pneumoniae", "organism": "Klebsiella pneumoniae",
               "scheme_name": "Klebsiella pneumoniae cgMLST", "locus_count": 2358,
               "path": "/cgmlst/Klebsiella_pneumoniae__cgmlst_org_2358"}


def named(name, genus="", species=""):
    organism = {"genus": genus, "species": species}
    return {"id": name, "name": name, "metadata": {"organism": organism}, "result": None}


def test_a_cgmlst_scheme_is_refused_for_a_cohort_of_another_organism():
    from wmlstudio.sample_workflow import cgmlst_scheme_applies
    cohort = [named("kp1", "Klebsiella", "pneumoniae"), named("ec1", "Escherichia", "coli")]
    verdict = cgmlst_scheme_applies(cohort, KPNEUMONIAE)
    assert verdict["applies"] is False and verdict["status"] == "organism_mismatch"
    assert [entry["name"] for entry in verdict["mismatched"]] == ["ec1"]
    assert "Klebsiella pneumoniae" in verdict["message"]
    assert "ec1 (Escherichia coli)" in verdict["message"]
    assert "never crosses organisms" in verdict["message"]


def test_a_cgmlst_scheme_is_not_confirmed_for_samples_whose_organism_is_unknown():
    from wmlstudio.sample_workflow import cgmlst_scheme_applies
    verdict = cgmlst_scheme_applies([named("unknown1"), named("unknown2")], KPNEUMONIAE)
    assert verdict["applies"] is False and verdict["status"] == "organism_unknown"
    assert [entry["name"] for entry in verdict["unknown"]] == ["unknown1", "unknown2"]
    assert "Identify or assign the organism first" in verdict["message"]
    # One identified sample is enough to run, and the unidentified one is named.
    mixed = cgmlst_scheme_applies([named("kp1", "Klebsiella", "pneumoniae"), named("u")],
                                  KPNEUMONIAE)
    assert mixed["applies"] is True
    assert any("were not checked" in warning for warning in mixed["warnings"])


def test_a_complex_member_may_run_the_scheme_and_is_told_what_it_was_defined_on():
    from wmlstudio.sample_workflow import cgmlst_scheme_applies
    verdict = cgmlst_scheme_applies([named("kv1", "Klebsiella", "variicola")], KPNEUMONIAE)
    assert verdict["applies"] is True and verdict["status"] == "applies"
    assert "Klebsiella variicola" in verdict["warnings"][0]
    assert "defined on Klebsiella pneumoniae" in verdict["warnings"][0]


def test_a_seven_locus_scheme_is_never_offered_as_a_core_genome_run():
    from wmlstudio.sample_workflow import cgmlst_scheme_applies
    classical = {"key": "kp-mlst", "genus": "Klebsiella", "species": "pneumoniae",
                 "scheme_name": "Klebsiella MLST", "locus_count": 7, "path": "/schemes/kp"}
    verdict = cgmlst_scheme_applies([named("kp1", "Klebsiella", "pneumoniae")], classical)
    assert verdict["applies"] is False and verdict["status"] == "not_a_cgmlst_scheme"
    assert "cannot be run as a cgMLST scheme" in verdict["message"]
    # A scheme that names no organism cannot be shown to apply to anything.
    unnamed = cgmlst_scheme_applies([named("kp1", "Klebsiella", "pneumoniae")],
                                    {"key": "x", "locus_count": 2000, "path": "/x"})
    assert unnamed["applies"] is False and unnamed["status"] == "scheme_organism_unknown"


# The launch dialog is the surface cgmlst_scheme_applies exists for, so the refusal
# and the plan it produces are checked here beside the rule they come from.
def test_the_run_plan_asks_for_cgmlst_beside_mlst_and_names_the_target_set(qtbot):
    from wmlstudio.analysis_plan import RunPlanDialog
    dialog = RunPlanDialog([named("kp1", "Klebsiella", "pneumoniae")],
                           cgmlst_schemes=[{**KPNEUMONIAE, "ready": True}])
    qtbot.addWidget(dialog)
    dialog.hydra.setChecked(False)
    dialog.cgmlst.setChecked(True)
    dialog.accept()
    assert dialog.plan["cgmlst"] is True and dialog.plan["hydra"] is False
    assert dialog.plan["cgmlst_scheme"]["key"] == "cgmlst.org:kpneumoniae-2358"
    assert dialog.plan["cgmlst_scheme"]["path"] == KPNEUMONIAE["path"]
    assert dialog.plan["cgmlst_scheme"]["locus_count"] == 2358


def test_the_run_plan_refuses_a_cgmlst_scheme_the_cohort_is_not(qtbot):
    from PySide6.QtWidgets import QDialog

    from wmlstudio.analysis_plan import RunPlanDialog
    dialog = RunPlanDialog([named("ec1", "Escherichia", "coli")],
                           cgmlst_schemes=[{**KPNEUMONIAE, "ready": True}])
    qtbot.addWidget(dialog)
    dialog.cgmlst.setChecked(True)
    dialog.accept()
    assert dialog.result() != QDialog.DialogCode.Accepted
    assert "never crosses organisms" in dialog.feedback.text()
    assert "Nothing was run" in dialog.feedback.text()
    dialog.cgmlst.setChecked(False)
    dialog.accept()
    assert dialog.plan["cgmlst"] is False and dialog.plan["cgmlst_scheme"] is None


def test_the_run_plan_offers_nothing_when_no_cgmlst_scheme_is_installed(qtbot):
    from wmlstudio.analysis_plan import RunPlanDialog
    dialog = RunPlanDialog([named("kp1", "Klebsiella", "pneumoniae")],
                           cgmlst_schemes=[{**KPNEUMONIAE, "path": "", "ready": False}])
    qtbot.addWidget(dialog)
    assert dialog.cgmlst_choices == []
    assert dialog.cgmlst.isEnabled() is False
    assert "cgMLST schemes tab" in dialog.cgmlst.toolTip()


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
