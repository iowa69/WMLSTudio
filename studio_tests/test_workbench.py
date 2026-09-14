"""End-to-end native microbiology workflows across the sample-centric workbench."""

import hashlib
import json
from pathlib import Path

import pytest
from PySide6.QtCore import QItemSelectionModel, Qt
from PySide6.QtWidgets import QComboBox, QDialog, QMessageBox, QTableWidget

from wmlstudio.app import MainWindow
from wmlstudio.library import export_bundle
from wmlstudio.project import Project
from wmlstudio.sample_workflow import set_cluster
from wmlstudio.typing import load_scheme
from wmlstudio.workflow_dialogs import BatchAssignmentDialog, ImportSamplesDialog

ARC = "AACCGTACGTTAG"
GYR = "TTGGCATACCTGA"


@pytest.fixture
def workbench(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    window = MainWindow(storage_root=tmp_path / "application")
    window.test_errors = errors
    qtbot.addWidget(window)
    window.show()
    yield window
    if window.worker and window.worker.isRunning():
        window.worker.cancel()
        qtbot.waitUntil(lambda: not window.worker.isRunning(), timeout=15000)
    qtbot.waitUntil(lambda: not window.worker_role, timeout=15000)
    window.close()


def idle(qtbot, window):
    qtbot.waitUntil(lambda: (window.worker is None or not window.worker.isRunning())
                   and not window.worker_role and window.run_button.isEnabled(), timeout=20000)


def make_scheme(path, locus="arcA", st="17"):
    path.mkdir(parents=True)
    (path / f"{locus}.tfa").write_text(f">{locus}_1\n{ARC}\n")
    (path / "gyrB.tfa").write_text(f">gyrB_1\n{GYR}\n")
    (path / "profiles.tsv").write_text(f"ST\t{locus}\tgyrB\n{st}\t1\t1\n")
    return path


def assembly(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">one\n{ARC}\n>two\n{GYR}\n")
    return path


def add_profile(window, name, genus="Klebsiella", species="pneumoniae", st="17", digest="MLST"):
    sid = window.project.add_profile(name, {
        "sample_name": name, "kind": "profile", "scheme": "MLST", "scheme_digest": digest,
        "status": "profile_imported", "st": st, "alleles": {"arcA": "1", "gyrB": "1"},
        "calls": [{"locus": "arcA", "allele": "1", "status": "external_exact"},
                  {"locus": "gyrB", "allele": "1", "status": "external_exact"}],
    }, {"organism": {"genus": genus, "species": species}})
    window.cohort_ids = set(window.cohort_ids or ()) | {sid}
    return sid


def select_ids(window, identifiers):
    window.navigate(1)
    window.sample_table.clearSelection()
    for row in range(window.sample_table.rowCount()):
        if window.sample_table.item(row, 0).data(Qt.ItemDataRole.UserRole) in identifiers:
            index = window.sample_table.model().index(row, 0)
            window.sample_table.selectionModel().select(
                index, QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows)
    window.sample_selection_changed()


def test_secondary_snapshot_from_changed_input_is_not_comparable(workbench, tmp_path):
    path = assembly(tmp_path / "changing.fasta")
    sid = workbench.project.add_sample(path)
    primary = {"scheme": "MLST", "scheme_digest": "primary", "status": "complete",
               "alleles": {"arcA": "1"}, "input_sha256": "old-input", "st": "17"}
    workbench.project.set_result(sid, primary)
    secondary = {**primary, "scheme": "cgMLST", "scheme_digest": "secondary", "st": None}
    workbench.project.set_analysis(sid, secondary)
    workbench.project.set_result(sid, {**primary, "input_sha256": "new-input", "st": "42"})
    available = workbench.available_profiles(workbench.project.get_sample(sid))
    assert [r["scheme_digest"] for r in available] == ["primary"]
    archived = workbench.project.analysis_results(sid)
    assert next(r for r in archived if r["scheme_digest"] == "secondary")["input_sha256"] == "old-input"


def test_pair_dialog_rejects_reusing_reads(qtbot):
    from wmlstudio.pairing_dialog import PairReadsDialog

    dialog = PairReadsDialog([
        {"id": "one", "name": "one", "input_path": "isolate_R1.fastq.gz"},
        {"id": "two", "name": "two", "input_path": "isolate_R2.fastq.gz"},
        {"id": "three", "name": "three", "input_path": "other.fastq"},
    ])
    qtbot.addWidget(dialog)
    assert dialog.combos[0][1].currentData() == "two"
    assert dialog.combos[1][1].currentData() is None
    third = dialog.combos[2][1]
    third.setCurrentIndex(third.findData("two"))
    dialog.accept()
    assert not dialog.assignments and "only one pair" in dialog.feedback.text()
    third.setCurrentIndex(0)
    dialog.accept()
    assert dialog.assignments == [{"primary_id": "one", "mate_id": "two"}]


def test_assembly_queue_associates_reads_then_runs_real_typing(workbench, qtbot, tmp_path, monkeypatch):
    """Native engine itself is CI-tested; use an explicit assembly stub to test Qt chaining."""
    import wmlstudio.assembly as assembler

    scheme = make_scheme(tmp_path / "scheme")
    identifiers = []
    for mate in (1, 2):
        path = tmp_path / f"isolate_R{mate}.fastq"
        path.write_text(f"@pair/{mate}\nACGTTGCAACGT\n+\nIIIIIIIIIIII\n")
        sid = workbench.project.add_sample(path)
        workbench.project.set_metadata(sid, {"workflow": {"typing_mode": "manual", "scheme_path": str(scheme)}})
        identifiers.append(sid)
    primary_id, mate_id = identifiers
    inputs_before = [workbench.project.get_sample(sid)["input_path"] for sid in identifiers]

    def test_only_assembly(read1, read2, destination, **kwargs):
        path = assembly(Path(destination) / "contigs.fasta")
        return {"assembly_path": str(path), "read_qc": [{"records": 1}, {"records": 1}],
                "provenance": {"engine": "TEST-ONLY assembly fixture", "inputs": [
                    {"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()}
                    for p in (read1, read2)],
                    "assembly_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}}

    monkeypatch.setattr(assembler, "run_skesa", test_only_assembly)
    workbench._run_ids = set(identifiers)
    workbench._run_plan = {}
    workbench._run_cancelled = False
    workbench._typing_override = str(scheme)
    workbench.assemble_pairs([{"primary_id": primary_id, "mate_id": mate_id}], {})
    idle(qtbot, workbench)
    primary, mate = [workbench.project.get_sample(sid) for sid in identifiers]
    assert primary["result"]["st"] == "17" and primary["status"] == "completed"
    assert primary["metadata"]["assembly"]["read_qc"][0]["records"] == 1
    assert mate["metadata"]["workflow"]["source_kind"] == "read_mate"
    assert mate["input_path"] == inputs_before[1]
    assert all(Path(path).is_file() for path in inputs_before)
    assert not workbench.test_errors


def test_import_dialog_actual_worker_copy_typing_and_st_filename(workbench, qtbot, tmp_path, monkeypatch):
    scheme = make_scheme(workbench.root / "schemes" / "local_mlst")
    workbench.populate_schemes()
    originals = [assembly(tmp_path / "incoming" / f"sample{i}.fasta") for i in range(2)]
    before = {path: path.read_bytes() for path in originals}
    managed = tmp_path / "managed"

    def configure(dialog):
        dialog.table.selectRow(0)
        dialog.mode.setCurrentIndex(dialog.mode.findData("manual"))
        dialog.genus.setText("Escherichia")
        dialog.species.setText("coli")
        dialog.scheme.setCurrentIndex(dialog.scheme.findData(str(scheme)))
        dialog.apply_selected()
        dialog.table.clearSelection()
        dialog.table.selectRow(1)
        dialog.mode.setCurrentIndex(dialog.mode.findData("unknown"))
        dialog.apply_selected()
        dialog.storage_root.setText(str(managed))
        dialog.append_st.setChecked(True)
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ImportSamplesDialog, "exec", configure)
    workbench.import_paths(originals, configure=True)
    idle(qtbot, workbench)
    samples = workbench.project.samples()
    assert len(samples) == 2
    assert {sample["id"] for sample in samples} == workbench.selection_ids
    assert all(sample["metadata"]["workflow"]["managed"] for sample in samples)
    assert samples[0]["metadata"]["workflow"]["scheme_path"] == str(scheme)
    assert samples[1]["metadata"]["workflow"]["typing_mode"] == "unknown"
    workbench.start_analysis()
    idle(qtbot, workbench)
    samples = {sample["name"]: sample for sample in workbench.project.samples()}
    typed = samples["sample0"]
    assert typed["result"]["st"] == "17"
    assert Path(typed["input_path"]).name == "sample0_ST_17.fasta"
    assert "Escherichia/coli/ST_17" in Path(typed["input_path"]).as_posix()
    assert samples["sample1"]["result"]["status"] == "qc_only"
    for original in originals:
        assert original.read_bytes() == before[original]
    assert typed["result"]["input_sha256"] == hashlib.sha256(before[originals[0]]).hexdigest()
    assert workbench.test_errors == []


def test_batch_assignment_ui_keeps_annotations_and_stable_selection(workbench, monkeypatch):
    first = add_profile(workbench, "first")
    second = add_profile(workbench, "second")
    workbench.project.update_metadata(first, {"annotations": {"ward": "ICU"}})
    workbench.refresh()
    select_ids(workbench, {first, second})

    def assign(dialog):
        dialog.mode.setCurrentIndex(dialog.mode.findData("manual"))
        dialog.genus.setText("Enterococcus")
        dialog.species.setText("faecium")
        dialog.apply_all()
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(BatchAssignmentDialog, "exec", assign)
    workbench.assign_selected()
    assert workbench.selection_ids == {first, second}
    assert all(sample["metadata"]["organism"]["genus"] == "Enterococcus"
               for sample in workbench.project.samples())
    assert workbench.project.get_sample(first)["metadata"]["annotations"]["ward"] == "ICU"
    assert len(workbench.project.analysis_results(first)) == 1
    assert workbench.test_errors == []


def test_library_filters_sorting_and_navigation_keep_hidden_selection(workbench):
    first = add_profile(workbench, "Zulu", "Escherichia", "coli", "2")
    second = add_profile(workbench, "Alpha", "Klebsiella", "pneumoniae", "10")
    third = add_profile(workbench, "Bravo", "Escherichia", "coli", "131")
    workbench.refresh()
    select_ids(workbench, {first, second})
    workbench.sample_table.sortItems(0, Qt.SortOrder.AscendingOrder)
    assert workbench.sample_table.item(0, 0).text() == "Alpha"
    assert workbench.selection_ids == {first, second}
    workbench.genus_filter.setCurrentIndex(workbench.genus_filter.findData("Escherichia"))
    assert workbench.sample_table.rowCount() == 2
    assert workbench.selection_ids == {first, second}
    workbench.select_visible_samples()
    assert workbench.selection_ids == {first, second, third}
    workbench.navigate(4)
    workbench.navigate(1)
    assert workbench.selection_ids == {first, second, third}
    workbench.clear_filters()
    assert workbench.sample_table.rowCount() == 3
    workbench.clear_sample_selection()
    assert workbench.selection_ids == set()


def test_selected_html_and_actual_pdf_use_same_highlighted_cohort(workbench, tmp_path, monkeypatch):
    chosen = add_profile(workbench, "CHOSEN_ISOLATE")
    add_profile(workbench, "EXCLUDED_ISOLATE")
    set_cluster(workbench.project, [chosen], "Investigation <A>", "#886644")
    workbench.refresh()
    select_ids(workbench, {chosen})
    workbench.report_selected()
    assert workbench.report_sample_ids() == {chosen}
    destination = tmp_path / "selected.html"
    monkeypatch.setattr("wmlstudio.ui_reports.QFileDialog.getSaveFileName",
                        lambda *args: (str(destination), "HTML"))
    workbench.export_report("html")
    report = destination.read_text()
    assert "CHOSEN_ISOLATE" in report and "EXCLUDED_ISOLATE" not in report
    assert "Investigation &lt;A&gt;" in report
    captured = []
    from wmlstudio import ui_reports
    real_report = ui_reports.review_report_html

    def audited_report(samples, **kwargs):
        captured.extend(sample["id"] for sample in samples)
        return real_report(samples, **kwargs)

    monkeypatch.setattr(ui_reports, "review_report_html", audited_report)
    destination = tmp_path / "selected.pdf"
    workbench.export_report("pdf")
    assert destination.read_bytes().startswith(b"%PDF")
    assert captured == [chosen]
    assert workbench.test_errors == []


def test_hydra_mapping_ui_keeps_ambiguous_names_unlinked_and_primary_matrix_correct(workbench, monkeypatch):
    first = add_profile(workbench, "duplicate")
    second = add_profile(workbench, "duplicate")
    report = {"import_provenance": {"sha256": "validated"}, "samples": [{
        "sample": "duplicate", "hits": [
            {"gene": "blaKPC-2", "element_type": "AMR", "primary": True},
            {"gene": "SECONDARY_ONLY", "element_type": "AMR", "primary": False},
        ], "summary": {"amr_genes": 1},
    }]}
    workbench.project.set_setting("hydra_report", report)
    defaults = []

    def map_dialog(dialog):
        table = dialog.findChild(QTableWidget)
        combo = table.cellWidget(0, 1)
        assert isinstance(combo, QComboBox)
        defaults.append(combo.currentData())
        combo.setCurrentIndex(combo.findData(first))
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(QDialog, "exec", map_dialog)
    workbench.feature_ids = {first, second}
    workbench.link_hydra_samples()
    assert defaults == [None]
    assert workbench.project.get_sample(first)["metadata"]["hydra"]["source_sample"] == "duplicate"
    assert "hydra" not in workbench.project.get_sample(second)["metadata"]
    assert "blaKPC-2" in workbench.amr_model.headers
    assert "SECONDARY_ONLY" not in workbench.amr_model.headers
    assert {row["blaKPC-2"] for row in workbench.amr_model.rows} == {"Present", "No report"}
    assert workbench.project.get_sample(first)["result"]["st"] == "17"
    assert workbench.test_errors == []


def test_bundle_import_without_fasta_can_compare_and_reexport(workbench, qtbot, tmp_path, monkeypatch):
    bundle = tmp_path / "shared.json"
    source_input = assembly(tmp_path / "original.fasta")
    with Project(tmp_path / "source.wmlstudio") as project:
        for name in ("Imported A", "Imported B"):
            sid = project.add_sample(source_input, name)
            project.set_result(sid, {"sample_name": name, "kind": "fasta", "scheme": "shared",
                                     "scheme_digest": "shared-fingerprint", "st": "7", "status": "complete",
                                     "alleles": {"l1": "1", "l2": "2"}, "calls": []})
        export_bundle(project, bundle)
    source_input.unlink()
    monkeypatch.setattr("wmlstudio.ui_reports.QFileDialog.getOpenFileName",
                        lambda *args: (str(bundle), "JSON"))
    workbench.import_sample_bundle()
    idle(qtbot, workbench)
    assert len(workbench.project.samples()) == 2
    assert all(sample["profile_only"] and not sample["missing_input"] for sample in workbench.project.samples())
    workbench.compare_selected()
    workbench.refresh_comparison()
    assert len(workbench.distance_rows) == 1
    assert workbench.distance_rows[0]["distance"] == 0
    workbench.start_analysis(all_samples=True)
    assert not workbench.worker_role
    output = tmp_path / "again.json"
    monkeypatch.setattr("wmlstudio.ui_reports.QFileDialog.getSaveFileName",
                        lambda *args: (str(output), "JSON"))
    workbench.report_selected()
    workbench.export_sample_bundle()
    assert len(json.loads(output.read_text())["samples"]) == 2
    assert workbench.test_errors == []


def test_secondary_scheme_worker_preserves_primary_mlst_and_builds_saved_matrix(workbench, qtbot, tmp_path, monkeypatch):
    cg_scheme = make_scheme(workbench.root / "schemes" / "cg_schema", locus="core001", st="99")
    source = assembly(tmp_path / "isolate.fasta")
    ids = []
    for name in ("A", "B"):
        sid = workbench.project.add_sample(source, name)
        workbench.project.set_result(sid, {"sample_name": name, "kind": "fasta", "scheme": "MLST",
                                          "scheme_digest": "original-mlst", "st": "17", "status": "complete",
                                          "alleles": {"arcA": "1", "gyrB": "1"}, "calls": []})
        ids.append(sid)
    workbench.populate_schemes()
    workbench.refresh()
    select_ids(workbench, set(ids))
    workbench.compare_selected()
    workbench.compare_scheme.setCurrentIndex(workbench.compare_scheme.findData(str(cg_scheme)))
    monkeypatch.setattr(QMessageBox, "question", lambda *args: QMessageBox.StandardButton.Yes)
    workbench.type_comparison_scheme()
    idle(qtbot, workbench)
    digest = load_scheme(cg_scheme).digest
    assert workbench.compare_scheme.currentData() == "digest:" + digest
    for sid in ids:
        assert workbench.project.get_sample(sid)["result"]["st"] == "17"
        assert {result["scheme_digest"] for result in workbench.project.analysis_results(sid)} == {"original-mlst", digest}
    assert workbench.profile_model.columnCount() == 5
    assert "core001" in workbench.profile_model.headers
    assert len(workbench.distance_rows) == 1 and workbench.distance_rows[0]["distance"] == 0
    assert workbench.test_errors == []
