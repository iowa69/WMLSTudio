"""Launch-plan, assembly/typing/AMR chaining and frozen evidence reuse in real Qt."""

import hashlib
import json
import time

import pytest
from PySide6.QtWidgets import QDialog

from wmlstudio.analysis_plan import RunPlanDialog
from wmlstudio.app import MainWindow
from wmlstudio.hydra import load_hydra_report
from wmlstudio.library import Library
from wmlstudio.library_dialog import reuse_library_profiles
from wmlstudio.pairing_dialog import PairReadsDialog, mate_hint
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, file_sha256

ARC, GYR = "AACCGTACGTTAG", "TTGGCATACCTGA"


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "app")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    if widget.worker and widget.worker.isRunning():
        widget.cancel_analysis()
        qtbot.waitUntil(lambda: not widget.worker.isRunning(), timeout=15000)
    qtbot.waitUntil(lambda: not widget.worker_role, timeout=15000)
    widget.close()


def idle(qtbot, window):
    qtbot.waitUntil(lambda: not window.worker_role and (not window.worker or not window.worker.isRunning())
                   and window.run_button.isEnabled(), timeout=20000)


def scheme(root):
    folder = root / "schemes" / "two_loci"
    folder.mkdir(parents=True)
    (folder / "arcA.tfa").write_text(f">arcA_1\n{ARC}\n")
    (folder / "gyrB.tfa").write_text(f">gyrB_1\n{GYR}\n")
    (folder / "profiles.tsv").write_text("ST\tarcA\tgyrB\n17\t1\t1\n")
    return folder


def hydra_report(tmp_path, identifier):
    path = tmp_path / f"{identifier}.json"
    path.write_text(json.dumps({"hydra_version": "1.4.0", "command": "test engine adapter",
        "parameters": {}, "databases": ["ncbi"], "samples": [{"sample": identifier,
        "input_type": "assembly", "inputs": [], "hits": [{"gene": "blaTEST", "database": "ncbi",
        "element_type": "AMR", "method": "blastn", "primary": True}]}]}))
    report = load_hydra_report(path)
    report["execution_provenance"] = {"engine": "test adapter", "version": "test"}
    return report


def test_pairing_suggestions_are_explicit_and_duplicate_mates_refused(qtbot):
    samples = [{"id": f"r{i}", "name": f"isolate_R{i}", "input_path": f"/incoming/isolate_R{i}_001.fastq.gz"} for i in (1, 2)]
    dialog = PairReadsDialog(samples)
    qtbot.addWidget(dialog)
    assert mate_hint(samples[0]["input_path"]) == ("isolate", 1)
    assert dialog.combos[0][1].currentData() == "r2"
    assert dialog.combos[1][1].currentData() is None
    dialog.combos[1][1].setCurrentIndex(dialog.combos[1][1].findData("r1"))
    dialog.accept()
    assert not dialog.assignments
    assert "only one pair" in dialog.feedback.text()
    dialog.combos[1][1].setCurrentIndex(0)
    dialog.accept()
    assert dialog.assignments == [{"primary_id": "r1", "mate_id": "r2"}]


def test_launch_review_does_not_allow_unavailable_amr_runtime(qtbot, monkeypatch):
    import wmlstudio.hydra_runtime as runtime
    monkeypatch.setattr(runtime, "runtime_capabilities",
                        lambda db_root=None: {"available": False, "message": "Missing native engine",
                                              "tools": {}, "databases": {}})
    dialog = RunPlanDialog([{"id": "one", "name": "Unknown isolate"}])
    qtbot.addWidget(dialog)
    dialog.hydra.setChecked(True)
    dialog.accept()
    assert dialog.result() != QDialog.DialogCode.Accepted
    # The refusal now names the missing piece and the one action that fixes it,
    # rather than echoing a runtime string the user cannot act on.
    refusal = dialog.feedback.text()
    assert "HYDRA cannot start" in refusal and "nothing was run" in refusal
    assert "Install" in refusal
    dialog.hydra.setChecked(False)
    dialog.accept()
    assert dialog.plan["hydra"] is False
    assert dialog.plan["thresholds"]["protein_min_coverage"] == 90


def test_launch_review_keeps_resource_controls_in_collapsed_advanced(qtbot):
    dialog = RunPlanDialog([{"id": "one", "name": "Reviewed isolate"}])
    qtbot.addWidget(dialog)
    dialog.show()
    assert not dialog.advanced.isVisible()
    assert "sample" in dialog.resource_summary.text()
    assert "threads each" in dialog.resource_summary.text()
    dialog.advanced_button.setChecked(True)
    assert dialog.advanced.isVisible()
    assert dialog.threads.isVisible() and dialog.memory.isVisible()


def test_typing_storage_and_hydra_execute_in_order_with_stable_ids(window, qtbot, tmp_path, monkeypatch):
    import wmlstudio.hydra_runtime as runtime
    source = tmp_path / "isolate.fasta"
    source.write_text(f">a\n{ARC}\n>b\n{GYR}\n")
    digest = file_sha256(source)
    reference = scheme(window.root)
    window.populate_schemes()
    window.import_assignments([{"path": str(source), "typing_mode": "manual", "scheme_path": str(reference),
                                "genus": "Escherichia", "species": "coli"}], {"managed": True, "append_st": True})
    idle(qtbot, window)
    identifier = window.project.samples()[0]["id"]
    calls = []

    def execute(paths, db_root, databases, **options):
        calls.append((paths, options))
        assert "ST_17" in paths[0]
        assert options["sample_names"] == [identifier]
        assert window.project.get_sample(identifier)["result"]["st"] == "17"
        return hydra_report(tmp_path, identifier)

    monkeypatch.setattr(runtime, "run_assemblies", execute)
    monkeypatch.setattr(window, "review_run_plan", lambda *args, **kwargs: {"hydra": True, "db_root": "test snapshot"})
    window.start_analysis(confirm=True)
    idle(qtbot, window)
    sample = window.project.get_sample(identifier)
    assert len(calls) == 1
    assert sample["metadata"]["hydra"]["source_sample"] == identifier
    assert sample["metadata"]["hydra"]["execution_provenance"]["engine"] == "test adapter"
    assert sample["result"]["st"] == "17"
    assert file_sha256(source) == digest
    assert window.test_errors == []


def test_cancelled_typing_never_launches_optional_hydra(window, qtbot, tmp_path, monkeypatch):
    import wmlstudio.hydra_runtime as runtime
    import wmlstudio.jobs as jobs
    source = tmp_path / "isolate.fasta"
    source.write_text(f">a\n{ARC}\n>b\n{GYR}\n")
    reference = scheme(window.root)
    window.import_paths([source])
    window.populate_schemes()
    window.scheme_combo.setCurrentIndex(window.scheme_combo.findData(str(reference)))
    calls = []

    def slow(path, scheme, cancelled=None, **kwargs):
        while not cancelled():
            time.sleep(0.005)
        raise AnalysisCancelled("Test cancellation")

    monkeypatch.setattr(jobs, "call_assembly", slow)
    monkeypatch.setattr(runtime, "run_assemblies", lambda *args, **kwargs: calls.append(True))
    monkeypatch.setattr(window, "review_run_plan", lambda *args, **kwargs: {"hydra": True, "db_root": "test snapshot"})
    window.start_analysis(confirm=True)
    qtbot.waitUntil(lambda: window.project.samples()[0]["status"] == "running", timeout=10000)
    window.cancel_analysis()
    idle(qtbot, window)
    assert not calls
    assert window.project.samples()[0]["status"] == "interrupted"


def test_assembly_association_then_typing_retains_both_read_records(window, qtbot, tmp_path, monkeypatch):
    import wmlstudio.assembly as assembly
    reference = scheme(window.root)
    reads = []
    for number in (1, 2):
        path = tmp_path / f"isolate_R{number}.fastq"
        path.write_text(f"@read1/{number}\nACGT\n+\nIIII\n")
        reads.append(path)
    originals = {path: path.read_bytes() for path in reads}
    window.import_paths(reads)
    window.populate_schemes()
    window.scheme_combo.setCurrentIndex(window.scheme_combo.findData(str(reference)))
    primary, mate = window.project.samples()

    def assemble(read1, read2, output_dir, **options):
        output_dir.mkdir(parents=True)
        path = output_dir / "contigs.fasta"
        path.write_text(f">a\n{ARC}\n>b\n{GYR}\n")
        return {"assembly_path": str(path), "provenance": {"assembly_sha256": file_sha256(path),
            "inputs": [{"path": str(source), "sha256": file_sha256(source)} for source in (read1, read2)]},
            "read_qc": [{}, {}], "pairing": {}, "qc": {}}

    def accept_pairing(dialog):
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(assembly, "run_skesa", assemble)
    monkeypatch.setattr(PairReadsDialog, "exec", accept_pairing)
    monkeypatch.setattr(window, "review_run_plan", lambda *args, **kwargs: {"assemble": True})
    window.start_analysis(confirm=True)
    idle(qtbot, window)
    result = window.project.get_sample(primary["id"])
    linked_mate = window.project.get_sample(mate["id"])
    assert result["result"]["st"] == "17"
    assert linked_mate["metadata"]["workflow"]["source_kind"] == "read_mate"
    assert linked_mate["result"] is None
    assert len(window.project.samples()) == 2
    for path, before in originals.items():
        assert path.read_bytes() == before
    assert window.test_errors == []


def test_reuse_frozen_cross_project_profiles_preserves_origin_and_secondary_results(tmp_path):
    profile = {"scheme": "classical", "scheme_digest": "classical-digest", "status": "profile_imported",
               "alleles": {"arcA": "1"}, "st": "17", "input_sha256": hashlib.sha256(b"genome").hexdigest()}
    with Project(tmp_path / "source.wmlstudio") as source, Project(tmp_path / "target.wmlstudio") as target:
        identifier = source.add_profile("saved isolate", profile, {"annotations": {"ward": "ICU"}})
        source.set_analysis(identifier, {**profile, "scheme": "cgMLST", "scheme_digest": "cg-snapshot"})
        with Library(tmp_path / "index.sqlite") as library:
            library.index_project(source)
            snapshots = library.search(st="17")
            imported = reuse_library_profiles(target, snapshots)
            assert imported == reuse_library_profiles(target, snapshots)
            assert len(target.samples()) == 1
            assert len(target.analysis_results(imported[0])) == 2
            assert target.get_sample(imported[0])["metadata"]["library_origin"]["sample_id"] == identifier
            assert target.get_sample(imported[0])["input_path"] == ""
            assert source.get_sample(identifier)["metadata"]["annotations"]["ward"] == "ICU"


def test_stale_amr_is_visible_as_archived_not_current_matrix_or_report(window, tmp_path):
    from wmlstudio.export import export_results
    from wmlstudio.sample_workflow import link_hydra
    source = tmp_path / "isolate.fa"
    source.write_text(">contig\nACGTACGT\n")
    identifier = window.project.add_sample(source)
    original_hash = file_sha256(source)
    window.project.set_result(identifier, {"status": "qc_only", "kind": "fasta", "input_sha256": original_hash})
    report = hydra_report(tmp_path, identifier)
    report["execution_provenance"]["inputs"] = [{"path": str(source), "sha256": original_hash}]
    link_hydra(window.project, report, {identifier: identifier})
    window.feature_ids = {identifier}
    window.refresh()
    assert window.amr_model.rows[0]["blaTEST"] == "Present"
    source.write_text(">contig\nTTTTCCCC\n")
    window.project.set_result(identifier, {"status": "qc_only", "kind": "fasta", "input_sha256": file_sha256(source)})
    window.refresh()
    assert window.amr_model.rows[0]["Evidence state"] == "stale"
    assert window.feature_model.rows[0]["AMR genes"] == ""
    assert window.project.get_sample(identifier)["metadata"]["hydra"]["hits"][0]["gene"] == "blaTEST"
    destination = tmp_path / "stale.html"
    export_results(window.report_records(), destination, "html")
    html = destination.read_text()
    assert "Archived metadata (old AMR evidence, not current calls)" in html
    assert "earlier or different input" in html
