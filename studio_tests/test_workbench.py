"""End-to-end native microbiology workflows across the sample-centric workbench."""

import hashlib
import json
from pathlib import Path

import pytest
from PySide6.QtCore import QItemSelectionModel, Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QInputDialog,
    QMenu,
    QMessageBox,
    QTableWidget,
)

from wmlstudio import archive
from wmlstudio.app import MainWindow
from wmlstudio.context_menus import SEPARATOR, Selection
from wmlstudio.library import export_bundle
from wmlstudio.project import Project
from wmlstudio.sample_workflow import set_cluster
from wmlstudio.storage import confirm_organism, import_samples, quarantine_samples
from wmlstudio.typing import load_scheme
from wmlstudio.ui_workbench import (
    REVIEW_STYLES,
    identification_support,
    organism_review,
)
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


def assembly(path, tail=""):
    """A typeable assembly. `tail` makes two isolates genuinely different bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">one\n{ARC}\n>two\n{GYR}\n" + (f">tail\n{tail}\n" if tail else ""))
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
    originals = [assembly(tmp_path / "incoming" / f"sample{i}.fasta", "ACGT" * (i + 1))
                 for i in range(2)]
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


def test_batch_assignment_ui_keeps_annotations_and_stable_selection(workbench, qtbot, monkeypatch):
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
    idle(qtbot, workbench)
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


def library_rows(window):
    return {window.sample_table.item(row, 0).text() for row in range(window.sample_table.rowCount())}


def tree_labels(window):
    return [window.library_tree.topLevelItem(index).text(0)
            for index in range(window.library_tree.topLevelItemCount())]


def test_right_click_archive_hides_an_isolate_and_restores_it_with_its_evidence(workbench, monkeypatch):
    add_profile(workbench, "kept")
    stored = add_profile(workbench, "stored")
    workbench.refresh()
    monkeypatch.setattr(QInputDialog, "getText", lambda *args, **kwargs: ("Closed investigation", True))

    workbench.context_archive(Selection("library", (stored,)))
    assert library_rows(workbench) == {"kept"}
    assert "Archived  (1)" in tree_labels(workbench)
    assert "1 archived" in workbench.selection_label.text()
    archived = workbench.project.get_sample(stored)
    assert archived["result"]["st"] == "17", "archiving keeps every stored result"
    assert archive.archive_record(archived)["reason"] == "Closed investigation"
    assert {entry["action"] for entry in workbench.project.history(stored)} >= {"sample_archived"}

    workbench.context_restore(Selection("library", (stored,)))
    assert library_rows(workbench) == {"kept", "stored"}
    assert "Archived  (1)" not in tree_labels(workbench)
    assert archive.archive_record(workbench.project.get_sample(stored)) is None
    assert workbench.test_errors == []


def test_right_click_assign_organism_moves_the_managed_copy_and_keeps_the_original(workbench, qtbot, tmp_path):
    source = assembly(tmp_path / "inbox" / "isolate.fasta")
    before = source.read_bytes()
    managed = tmp_path / "managed"
    sid = import_samples(workbench.project, [{"path": str(source), "typing_mode": "manual",
                                              "genus": "Klebsiella", "species": "pneumoniae"}],
                         str(managed))[0]
    workbench.refresh()

    workbench.context_assign_organism_quick(Selection("library", (sid,)), ("Enterobacter", "cloacae"))
    idle(qtbot, workbench)

    filed = Path(workbench.project.get_sample(sid)["input_path"])
    assert filed.relative_to(managed).parts[:3] == ("Enterobacter", "cloacae", "ST_unassigned")
    assert sid in filed.parts and filed.is_file()
    assert not (managed / "Klebsiella").exists(), "the vacated organism folder is pruned"
    assert source.read_bytes() == before and source.is_file()
    evidence = workbench.project.get_sample(sid)["metadata"]["organism_evidence"]
    assert evidence["basis"] == "user_assigned" and evidence["confidence"] == "unresolved"
    assert evidence["accepted"] == {"genus": "Enterobacter", "species": "cloacae"}
    assert workbench.test_errors == []


def test_removing_an_isolate_keeps_the_original_file_and_can_be_undone(workbench, tmp_path, monkeypatch):
    source = assembly(tmp_path / "inbox" / "gone.fasta")
    sid = workbench.project.add_sample(source, "gone")
    workbench.project.set_result(sid, {"sample_name": "gone", "kind": "fasta", "scheme": "MLST",
                                       "scheme_digest": "digest", "status": "complete", "st": "17",
                                       "alleles": {"arcA": "1"}, "calls": []})
    workbench.refresh()
    monkeypatch.setattr(QMessageBox, "exec", lambda self: QMessageBox.StandardButton.Yes)

    workbench.context_remove(Selection("library", (sid,)))
    assert workbench.project.samples() == []
    assert source.is_file(), "removal never deletes the user's own file"

    menu = QMenu(workbench)
    workbench.fill_recently_removed(menu)
    entries = [action for action in menu.actions() if action.isEnabled()]
    assert len(entries) == 1 and entries[0].text().startswith("gone · removed ")
    entries[0].trigger()
    restored = workbench.project.samples()[0]
    assert restored["name"] == "gone" and restored["result"]["st"] == "17"
    assert library_rows(workbench) == {"gone"}
    assert workbench.test_errors == []


def test_the_library_right_click_offers_only_commands_this_window_can_run(workbench):
    selection = Selection("library", ("first", "second"))
    plan = [entry for entry in workbench.context_menu_plan(selection) if entry is not SEPARATOR]
    titles = [entry.format_title(selection) for entry in plan]

    assert all(callable(getattr(workbench, entry.handler, None)) for entry in plan)
    assert "Archive 2 isolates (keeps all evidence)…" in titles
    assert titles.index("Archive 2 isolates (keeps all evidence)…") < \
        titles.index("Remove 2 isolates from this project…"), "removal is last, never the default"
    assert set(workbench.context_menu_report()) <= {"context_select_group"}
    for view in ("library", "library.tree", "overview.recent"):
        assert workbench._context_adapters[view].view_id == view


def test_the_cohort_picker_folder_nodes_resolve_to_isolates_and_leave_archived_out(workbench, qtbot):
    from wmlstudio.cohort_picker import PICKER_HANDLERS, CohortPickerDialog
    first = add_profile(workbench, "kp1", "Klebsiella", "pneumoniae")
    second = add_profile(workbench, "kp2", "Klebsiella", "oxytoca")
    third = add_profile(workbench, "ec", "Escherichia", "coli")
    archive.archive_samples(workbench.project, [third], "not this investigation")
    workbench.refresh()

    dialog = CohortPickerDialog(workbench.project.samples(), workbench.project,
                                "Which isolates?", None, workbench)
    qtbot.addWidget(dialog)
    assert {sample["id"] for sample in dialog.samples} == {first, second}
    assert sorted(dialog.folder_members("Klebsiella")) == sorted([first, second])
    assert dialog.folder_members("Klebsiella", "oxytoca") == [second]
    assert dialog.folder_members("Escherichia") == []
    assert all(callable(getattr(workbench, name, None)) for name in PICKER_HANDLERS)
    assert dialog.table.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu
    assert dialog.folders.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu


def test_refiling_shows_where_each_managed_copy_goes_before_moving_it(workbench, qtbot, tmp_path, monkeypatch):
    from wmlstudio.storage import confirm_organism
    source = assembly(tmp_path / "inbox" / "isolate.fasta")
    managed = tmp_path / "managed"
    sid = import_samples(workbench.project, [{"path": str(source), "typing_mode": "manual",
                                              "genus": "Klebsiella", "species": "pneumoniae"}],
                         str(managed))[0]
    confirm_organism(workbench.project, [sid], "Enterobacter", "cloacae")
    workbench.refresh()
    workbench.selection_ids = {sid}
    asked = {}

    def confirm(parent, title, text, *args, **kwargs):
        asked["text"] = text
        return QMessageBox.StandardButton.Yes

    monkeypatch.setattr(QMessageBox, "question", confirm)
    workbench.refile_selected()
    idle(qtbot, workbench)

    assert "Enterobacter/cloacae/ST_unassigned" in asked["text"]
    assert "filing decision rather than a laboratory identification" in asked["text"]
    filed = Path(workbench.project.get_sample(sid)["input_path"])
    assert filed.relative_to(managed).parts[0] == "Enterobacter" and filed.is_file()
    assert source.is_file()
    assert workbench.test_errors == []


# --- the Samples hub: one load, four angles, and a gate before every run -----


def column_of(table, title):
    """The column a header names, so a reordered view never silently changes a test."""
    return next(index for index in range(table.columnCount())
                if table.horizontalHeaderItem(index).text() == title)


def row_of(table, name):
    return next(row for row in range(table.rowCount()) if table.item(row, 0).text() == name)


def headers(table):
    return [table.horizontalHeaderItem(i).text() for i in range(table.columnCount())]


def local_scheme(window, folder, genus, species="", *, loci=2, name=None):
    """An installed reference that records the organism it was built from."""
    path = window.root / "schemes" / folder
    path.mkdir(parents=True, exist_ok=True)
    for index in range(loci):
        (path / f"locus{index:03d}.tfa").write_text(f">locus{index:03d}_1\n{ARC}\n")
    (path / "scheme.json").write_text(json.dumps(
        {"name": name or folder, "organism": {"genus": genus, "species": species}}))
    window.populate_schemes()
    return path


def evidence_for(genus, species, *, basis="genomic_ani", confidence="genomic_reference_supported",
                 status="confirmed"):
    accepted = {"genus": genus, "species": species} if status == "confirmed" else {"genus": "", "species": ""}
    return {"format_version": 1, "status": status, "basis": basis, "confidence": confidence,
            "proposed": {"genus": genus, "species": species}, "accepted": accepted,
            "quarantine_reason": None}


def identified(window, sid, genus, species, **kwargs):
    window.project.update_metadata(sid, {"organism": {"genus": genus, "species": species},
                                         "organism_evidence": evidence_for(genus, species, **kwargs)})


def select_only(window, sid):
    window.selection_ids = {sid}
    window.refresh()
    window.update_step_gates()


def test_st_and_cgmlst_are_separate_sub_tabs_with_their_own_schemes_and_denominators(workbench, tmp_path):
    source = assembly(tmp_path / "inbox" / "both.fasta", "AAAA")
    sid = workbench.project.add_sample(source, "both")
    identified(workbench, sid, "Klebsiella", "pneumoniae")
    workbench.project.set_result(sid, {"sample_name": "both", "kind": "fasta", "status": "complete",
                                       "scheme": "klebsiella_mlst", "scheme_digest": "mlst-digest",
                                       "st": "17", "alleles": {"arcA": "1", "gyrB": "1"}, "calls": []},
                                 kind="mlst")
    targets = {f"core{index:04d}": (str(index + 1) if index < 35 else None) for index in range(40)}
    workbench.project.set_result(sid, {"sample_name": "both", "kind": "fasta", "status": "incomplete",
                                       "scheme": "klebsiella_cgmlst", "scheme_digest": "cg-digest",
                                       "st": None, "alleles": targets, "calls": []}, kind="cgmlst")
    workbench.refresh()

    assert [workbench.sample_tabs.tabText(i) for i in range(workbench.sample_tabs.count())] == \
        ["Assembly", "ST", "cgMLST", "HYDRA"]

    roster = workbench.sample_table
    assert roster.item(row_of(roster, "both"), column_of(roster, "ST (7-locus)")).text() == "ST 17", \
        "a later cgMLST run must not take over the seven-locus ST column"

    st_table = workbench.step_tables["st"]
    st_row = row_of(st_table, "both")
    assert st_table.item(st_row, column_of(st_table, "7-locus scheme")).text() == "klebsiella_mlst"
    assert st_table.item(st_row, column_of(st_table, "ST")).text() == "ST 17"
    assert st_table.item(st_row, column_of(st_table, "Loci called")).text() == "2 / 2"

    cg_table = workbench.step_tables["cgmlst"]
    cg_row = row_of(cg_table, "both")
    assert cg_table.item(cg_row, column_of(cg_table, "cgMLST scheme")).text() == "klebsiella_cgmlst"
    assert cg_table.item(cg_row, column_of(cg_table, "Targets in scheme")).text() == "40"
    assert cg_table.item(cg_row, column_of(cg_table, "Loci called")).text() == "35"
    assert cg_table.item(cg_row, column_of(cg_table, "Missing")).text() == "5"

    # The two quantities never share a column name, so no reader can line them up.
    assert "ST" in headers(st_table) and "ST" not in headers(cg_table)
    assert "Targets in scheme" not in headers(st_table)
    assert workbench.test_errors == []


def test_a_genome_comparison_and_a_typing_panel_that_disagree_are_flagged_in_colour_and_in_words(workbench, tmp_path):
    local_scheme(workbench, "disagreeing_ecoli", "Escherichia", "coli")
    source = assembly(tmp_path / "inbox" / "mixed.fasta", "CCCC")
    sid = workbench.project.add_sample(source, "mixed")
    identified(workbench, sid, "Klebsiella", "pneumoniae")
    workbench.project.set_result(sid, {"sample_name": "mixed", "kind": "fasta", "status": "complete",
                                       "scheme": "disagreeing_ecoli", "scheme_digest": "ec-digest",
                                       "st": "11", "alleles": {"locus000": "1"}, "calls": []},
                                 kind="mlst")
    workbench.refresh()

    roster = workbench.sample_table
    cell = roster.item(row_of(roster, "mixed"), column_of(roster, "Organism evidence"))
    assert "Needs organism review" in cell.text()
    assert cell.background().color().name() == REVIEW_STYLES["conflict"][0].casefold()
    assert "disagree" in cell.toolTip() and "Klebsiella" in cell.toolTip()
    st_table = workbench.step_tables["st"]
    marker = st_table.item(row_of(st_table, "mixed"), column_of(st_table, "Organism review"))
    assert "Needs organism review" in marker.text(), "the words carry it, not only the colour"
    assert marker.background().color().name() == REVIEW_STYLES["conflict"][0].casefold()

    # Agreement is the quiet case: no amber, and the panel is still only corroboration.
    local_scheme(workbench, "agreeing_kp", "Klebsiella", "pneumoniae")
    workbench.project.set_result(sid, {"sample_name": "mixed", "kind": "fasta", "status": "complete",
                                       "scheme": "agreeing_kp", "scheme_digest": "kp-digest",
                                       "st": "258", "alleles": {"locus000": "1"}, "calls": []},
                                 kind="mlst")
    workbench.refresh()
    settled = roster.item(row_of(roster, "mixed"), column_of(roster, "Organism evidence"))
    assert "Needs organism review" not in settled.text()
    assert "corroborates" in settled.toolTip()
    assert "not an independent identification" in settled.toolTip()
    assert workbench.test_errors == []


def test_a_species_complex_or_an_undecided_organism_is_marked_for_review_in_a_second_colour(workbench, tmp_path):
    complex_source = assembly(tmp_path / "inbox" / "complex.fasta", "GGGG")
    complex_id = workbench.project.add_sample(complex_source, "complex")
    identified(workbench, complex_id, "Escherichia", "coli", confidence="complex_only")
    waiting_source = assembly(tmp_path / "inbox" / "waiting.fasta", "TTTT")
    waiting_id = workbench.project.add_sample(waiting_source, "waiting")
    quarantine_samples(workbench.project, [waiting_id], "not_in_reference_panel")
    workbench.refresh()

    roster = workbench.sample_table
    column = column_of(roster, "Organism evidence")
    complex_cell = roster.item(row_of(roster, "complex"), column)
    waiting_cell = roster.item(row_of(roster, "waiting"), column)
    assert "species complex" in complex_cell.text()
    assert complex_cell.background().color().name() == REVIEW_STYLES["conflict"][0].casefold()
    assert "not in reference panel" in waiting_cell.text()
    assert waiting_cell.background().color().name() == REVIEW_STYLES["undecided"][0].casefold()
    assert "Needs organism review" in waiting_cell.text()
    assert workbench.test_errors == []


def test_typing_is_disabled_until_an_organism_and_a_matching_scheme_are_both_set(workbench, tmp_path):
    source = assembly(tmp_path / "inbox" / "gated.fasta", "AACC")
    sid = workbench.project.add_sample(source, "gated")
    quarantine_samples(workbench.project, [sid], "awaiting_identification")
    select_only(workbench, sid)

    assert not workbench.step_buttons["st"].isEnabled()
    assert "an accepted organism" in workbench.step_notes["st"].text()
    assert "Assign organism" in workbench.step_notes["st"].text()

    confirm_organism(workbench.project, [sid], "Testarella", "ficta")
    select_only(workbench, sid)
    assert not workbench.step_buttons["st"].isEnabled()
    assert "seven-locus MLST scheme for Testarella ficta" in workbench.step_notes["st"].text()
    assert "Schemes" in workbench.step_notes["st"].text()

    local_scheme(workbench, "testarella_mlst", "Testarella", "ficta")
    select_only(workbench, sid)
    assert workbench.step_buttons["st"].isEnabled()
    assert "1 of 1 selected isolates are ready for ST" in workbench.step_notes["st"].text()
    # A seven-locus scheme is not a core-genome scheme, so cgMLST stays blocked.
    assert not workbench.step_buttons["cgmlst"].isEnabled()
    assert "cgMLST scheme for Testarella ficta" in workbench.step_notes["cgmlst"].text()
    assert workbench.test_errors == []


def test_raw_reads_and_profile_only_isolates_say_what_they_need_instead_of_failing_a_run(workbench, tmp_path):
    reads = tmp_path / "inbox" / "isolate_R1.fastq"
    reads.parent.mkdir(parents=True, exist_ok=True)
    reads.write_text("@read1\nACGT\n+\nIIII\n")
    read_id = workbench.project.add_sample(reads, "reads")
    identified(workbench, read_id, "Klebsiella", "pneumoniae")
    local_scheme(workbench, "kp_reads_scheme", "Klebsiella", "pneumoniae")
    select_only(workbench, read_id)
    assert not workbench.step_buttons["st"].isEnabled()
    assert "an assembly" in workbench.step_notes["st"].text()
    assert "Assemble them on the Assembly tab first" in workbench.step_notes["st"].text()
    assert workbench.step_buttons["assembly"].isEnabled(), "reads are what the Assembly tab is for"

    profile_id = add_profile(workbench, "profile only")
    select_only(workbench, profile_id)
    assert not workbench.step_buttons["assembly"].isEnabled()
    assert "no sequence file" in workbench.step_notes["assembly"].text()
    assert workbench.test_errors == []


def test_dropping_files_identifies_them_files_them_and_never_makes_a_second_isolate(workbench, qtbot, tmp_path, monkeypatch):
    inbox = tmp_path / "inbox"
    first = assembly(inbox / "KPNIH1.fasta", "AAAA")
    stranger = assembly(inbox / "stranger.fasta", "CCCC")
    # A panel is present for this test whatever this machine has installed; the
    # identification itself is stubbed, so the panel is named but never read.
    panel = tmp_path / "panel"
    panel.mkdir()
    monkeypatch.setattr(MainWindow, "installed_species_panel", lambda self: panel)

    def verdict(path, genus="", species="", **extra):
        record = {"format_version": 1, "input_path": str(path), "input_sha256": "", "kind": "fasta",
                  "basis": "genomic_ani" if genus else "none",
                  "confidence": "genomic_reference_supported" if genus else "unresolved",
                  "status": "proposed", "proposed": {"genus": genus, "species": species},
                  "accepted": {"genus": "", "species": ""}, "quarantine_reason": None,
                  "confirmed_by": None, "confirmed_utc": None, "panel": {}, "detail": {},
                  "margin_ani": None, "runner_up": None, "reason": "", "notes": [],
                  "limitations": [], "errors": []}
        record.update(extra)
        return record

    monkeypatch.setattr("wmlstudio.organism_id.identify_batch", lambda paths, **kwargs: [
        verdict(path, "Klebsiella", "pneumoniae") if Path(path).name == "KPNIH1.fasta"
        else verdict(path, quarantine_reason="not_in_reference_panel", status="quarantined")
        for path in paths])

    workbench.intake_paths([str(inbox)])
    idle(qtbot, workbench)

    root = workbench.project_path.with_suffix(".files")
    samples = {sample["name"]: sample for sample in workbench.project.samples()}
    assert set(samples) == {"KPNIH1", "stranger"}
    assert Path(samples["KPNIH1"]["input_path"]).relative_to(root).parts[:2] == ("Klebsiella", "pneumoniae")
    assert Path(samples["stranger"]["input_path"]).relative_to(root).parts[:2] == ("_Unresolved", "Not_in_reference_panel")
    assert first.is_file() and stranger.is_file(), "the user's own files are never moved"
    assert "1 are marked 'Needs organism review'" in workbench.progress_text.text()

    roster = workbench.sample_table
    flagged = roster.item(row_of(roster, "stranger"), column_of(roster, "Organism evidence"))
    assert "Needs organism review" in flagged.text()
    assert flagged.background().color().name() == REVIEW_STYLES["undecided"][0].casefold()

    workbench.intake_paths([str(inbox)])
    idle(qtbot, workbench)
    assert len(workbench.project.samples()) == 2, "the same bytes stay one isolate"
    assert "already in this project" in workbench.progress_text.text()
    assert workbench.test_errors == []


def test_the_hydra_organism_is_your_own_assignment_and_gates_the_organism_specific_tools(workbench, tmp_path, monkeypatch):
    from wmlstudio.provisioning import Requirement

    def store(*args, **kwargs):
        return Requirement(key="hydra_database", title="HYDRA AMR database store",
                           capability="HYDRA resistance screening", blocks="every HYDRA run stops",
                           required=True, state="ready", reason="Test store.", probe="Test probe.",
                           detail={"organisms": ["Klebsiella pneumoniae", "Staphylococcus aureus"]})

    monkeypatch.setattr("wmlstudio.provisioning.hydra_database_requirement", store)
    source = assembly(tmp_path / "inbox" / "amr.fasta", "AAAA")
    sid = workbench.project.add_sample(source, "amr")
    identified(workbench, sid, "Klebsiella", "pneumoniae")
    # Choosing a snapshot is what makes the application re-read the installed store.
    workbench.project.set_setting("hydra_database_root", str(tmp_path / "snapshot"))
    select_only(workbench, sid)

    assert workbench.hydra_organism.findData("Staphylococcus aureus") > 0
    assert "point-mutation catalogues" in workbench.hydra_organism_note.text()
    assert workbench.hydra_organism_for(workbench.project.get_sample(sid)) == "Klebsiella pneumoniae"

    workbench.hydra_organism.setCurrentIndex(workbench.hydra_organism.findData("Staphylococcus aureus"))
    workbench.apply_hydra_organism("selected")
    sample = workbench.project.get_sample(sid)
    assert sample["metadata"]["workflow"]["hydra_organism"] == "Staphylococcus aureus"
    assert sample["metadata"]["organism"] == {"genus": "Klebsiella", "species": "pneumoniae"}, \
        "choosing a HYDRA organism never rewrites the organism the evidence supported"
    assert "your assignment" in workbench.progress_text.text()
    hydra_table = workbench.step_tables["hydra"]
    assert hydra_table.item(row_of(hydra_table, "amr"),
                            column_of(hydra_table, "HYDRA organism")).text() == "Staphylococcus aureus"

    workbench.hydra_organism.setCurrentIndex(workbench.hydra_organism.findData(""))
    workbench.apply_hydra_organism("project")
    assert workbench.project.get_sample(sid)["metadata"]["workflow"]["hydra_organism"] == ""
    assert "stay unavailable" in workbench.progress_text.text()
    assert workbench.test_errors == []


def test_hydra_is_blocked_while_no_amr_database_is_installed_and_says_where_to_get_one(workbench, tmp_path, monkeypatch):
    from wmlstudio.provisioning import Requirement

    def empty(*args, **kwargs):
        return Requirement(key="hydra_database", title="HYDRA AMR database store",
                           capability="HYDRA resistance screening", blocks="every HYDRA run stops",
                           required=True, state="missing",
                           reason="No AMR database is installed in any store this application "
                                  "looks at, so HYDRA has nothing to search.",
                           probe="Test probe.", detail={"organisms": []})

    monkeypatch.setattr("wmlstudio.provisioning.hydra_database_requirement", empty)
    source = assembly(tmp_path / "inbox" / "amr.fasta", "AAAA")
    sid = workbench.project.add_sample(source, "amr")
    identified(workbench, sid, "Klebsiella", "pneumoniae")
    workbench.project.set_setting("hydra_database_root", str(tmp_path / "empty-store"))
    select_only(workbench, sid)

    assert not workbench.step_buttons["hydra"].isEnabled()
    assert "an installed AMR database" in workbench.step_notes["hydra"].text()
    assert "nothing to search" in workbench.step_notes["hydra"].text()
    assert "organism rules are unavailable" in workbench.hydra_organism_note.text()
    assert "no point mutation can be called" in workbench.hydra_organism_note.text()
    assert workbench.test_errors == []


def test_the_rarely_used_controls_stay_behind_the_advanced_disclosure(workbench):
    workbench.navigate(1)
    assert not workbench.samples_advanced.isVisible()
    assert not workbench.samples_advanced_button.isChecked()
    for column in ("ST (7-locus)", "AMR genes", "Scheme"):
        assert workbench.sample_table.isColumnHidden(column_of(workbench.sample_table, column))
    for column in ("Sample", "Organism evidence", "Storage"):
        assert not workbench.sample_table.isColumnHidden(column_of(workbench.sample_table, column))

    workbench.samples_advanced_button.setChecked(True)
    assert workbench.samples_advanced.isVisible()
    workbench.sample_columns_button.setChecked(True)
    assert not workbench.sample_table.isColumnHidden(column_of(workbench.sample_table, "AMR genes"))
    assert workbench.test_errors == []


def test_a_typing_panel_match_corroborates_an_organism_but_never_becomes_the_identification():
    panel_only = {"id": "a", "name": "a", "metadata": {"organism_evidence": {
        "status": "proposed", "basis": "mlst_panel", "confidence": "panel_compatibility",
        "proposed": {"genus": "Enterobacter", "species": ""},
        "accepted": {"genus": "", "species": ""},
        "detail": {"mlst_top": [{"scheme": "ecloacae", "organism_label": "Enterobacter cloacae"}]}}}}
    support = identification_support(panel_only, {})
    assert support["agreement"] == "panel_only"
    assert "not a species identification" in support["sentence"]
    assert organism_review(panel_only, {})["needs_review"] is True

    genomic = {"id": "b", "name": "b", "metadata": {
        "organism": {"genus": "Klebsiella", "species": "pneumoniae"},
        "organism_evidence": {"status": "confirmed", "basis": "genomic_ani",
                              "confidence": "genomic_reference_supported",
                              "proposed": {"genus": "Klebsiella", "species": "pneumoniae"},
                              "accepted": {"genus": "Klebsiella", "species": "pneumoniae"}}}}
    alone = identification_support(genomic, {})
    assert alone["agreement"] == "claim_only"
    assert "No classical typing panel has corroborated it yet" in alone["sentence"]
    assert organism_review(genomic, {})["needs_review"] is False

    # An organism a person typed in is a claim too, so a panel that contradicts it
    # is still a conflict worth a second look — and it is never silently corrected.
    assigned = {"id": "c", "name": "c", "metadata": {
        "organism": {"genus": "Enterobacter", "species": "cloacae"},
        "organism_evidence": {"status": "confirmed", "basis": "user_assigned",
                              "confidence": "unresolved",
                              "proposed": {"genus": "Enterobacter", "species": "cloacae"},
                              "accepted": {"genus": "Enterobacter", "species": "cloacae"}}}}
    schemes = {"kp_panel": {"name": "kp_panel", "id": "kp_panel", "genus": "Klebsiella",
                            "species": "pneumoniae", "kind": "mlst"}}
    disputed = identification_support(assigned, schemes, {"scheme": "kp_panel", "st": "258"})
    assert disputed["agreement"] == "conflict"
    assert "You assigned it" in disputed["sentence"] and "Klebsiella" in disputed["sentence"]
    review = organism_review(assigned, schemes, {"scheme": "kp_panel", "st": "258"})
    assert review["needs_review"] is True and review["level"] == "conflict"
    assert assigned["metadata"]["organism"] == {"genus": "Enterobacter", "species": "cloacae"}


def test_the_cohort_picker_says_which_isolates_still_need_an_organism(workbench, qtbot, tmp_path):
    from wmlstudio.cohort_picker import CohortPickerDialog
    ready = add_profile(workbench, "settled")
    # A later core-genome result must not empty the seven-locus column of the picker.
    targets = {f"core{index:04d}": str(index + 1) for index in range(40)}
    workbench.project.set_result(ready, {"sample_name": "settled", "kind": "profile",
                                         "scheme": "kp_cgmlst", "scheme_digest": "cg-digest",
                                         "status": "complete", "st": None, "alleles": targets,
                                         "calls": []}, kind="cgmlst")
    source = assembly(tmp_path / "inbox" / "unsettled.fasta", "GGGG")
    unsettled = workbench.project.add_sample(source, "unsettled")
    quarantine_samples(workbench.project, [unsettled], "low_confidence")
    workbench.refresh()

    dialog = CohortPickerDialog(workbench.project.samples(), workbench.project,
                                "Which isolates?", None, workbench)
    qtbot.addWidget(dialog)
    review = column_of(dialog.table, "Organism review")
    rows = {dialog.table.item(row, 1).text(): row for row in range(dialog.table.rowCount())}
    assert dialog.table.item(rows["settled"], review).text() == "Ready"
    flagged = dialog.table.item(rows["unsettled"], review)
    assert "Needs organism review" in flagged.text()
    assert flagged.background().color().name() == REVIEW_STYLES["undecided"][0].casefold()
    assert dialog.table.item(rows["settled"], column_of(dialog.table, "ST (7-locus)")).text() == "17"


def test_loading_without_a_reference_panel_says_so_instead_of_quietly_quarantining(workbench, qtbot, tmp_path, monkeypatch):
    # Pinned to the state the test is about, in both directions: no installed panel
    # and no installed schemes, whatever this machine happens to have staged.
    monkeypatch.setattr(MainWindow, "installed_species_panel", lambda self: None)
    monkeypatch.setattr("wmlstudio.characterization_refs.bundled_reference_root", lambda: None)
    workbench.scheme_paths = []
    source = assembly(tmp_path / "inbox" / "orphan.fasta", "AAAA")

    workbench.intake_paths([str(source)])
    idle(qtbot, workbench)

    sample = workbench.project.samples()[0]
    root = workbench.project_path.with_suffix(".files")
    assert Path(sample["input_path"]).relative_to(root).parts[0] == "_Unresolved"
    assert "No species reference panel is installed" in workbench.progress_text.text()
    assert "Install the species reference panel" in workbench.progress_text.text()
    assert source.is_file()
    assert workbench.test_errors == []


def test_running_the_st_step_types_against_a_seven_locus_scheme_and_records_it_as_mlst(workbench, qtbot, tmp_path, monkeypatch):
    scheme = make_scheme(workbench.root / "schemes" / "kp_seven")
    (scheme / "scheme.json").write_text(json.dumps(
        {"name": "kp_seven", "organism": {"genus": "Klebsiella", "species": "pneumoniae"}}))
    workbench.populate_schemes()
    source = assembly(tmp_path / "inbox" / "kp.fasta")
    sid = workbench.project.add_sample(source, "kp")
    identified(workbench, sid, "Klebsiella", "pneumoniae")
    monkeypatch.setattr(MainWindow, "review_run_plan", lambda self, *args, **kwargs: {})
    select_only(workbench, sid)
    assert workbench.step_buttons["st"].isEnabled()

    workbench.run_step("st")
    idle(qtbot, workbench)

    assert workbench.project.latest_analysis(sid, "mlst")["st"] == "17"
    assert workbench.project.latest_analysis(sid, "cgmlst") is None, \
        "typing seven loci must never be recorded as a core-genome profile"
    st_table, cg_table = workbench.step_tables["st"], workbench.step_tables["cgmlst"]
    assert st_table.item(row_of(st_table, "kp"), column_of(st_table, "ST")).text() == "ST 17"
    assert cg_table.item(row_of(cg_table, "kp"), column_of(cg_table, "cgMLST scheme")).text() == "—"
    assert cg_table.item(row_of(cg_table, "kp"), column_of(cg_table, "Targets in scheme")).text() == "—"
    stored = workbench.project.get_sample(sid)["metadata"].get("workflow") or {}
    assert stored.get("scheme_path") is None, \
        "a step run chooses a scheme for that run; it does not repin the isolate"
    assert source.is_file()
    assert workbench.test_errors == []
