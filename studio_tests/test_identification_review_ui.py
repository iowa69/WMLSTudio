"""Import that identifies the user's own files first, then files only what was accepted."""

import gzip

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog

from wmlstudio import practice_cohorts, storage
from wmlstudio.app import MainWindow
from wmlstudio.workflow_dialogs import (
    IdentificationReviewDialog,
    ImportSamplesDialog,
    PracticeCohortDialog,
)

ARC = "AACCGTACGTTAG"
GYR = "TTGGCATACCTGA"


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    if widget.worker and widget.worker.isRunning():
        widget.worker.cancel()
        qtbot.waitUntil(lambda: not widget.worker.isRunning(), timeout=15000)
    qtbot.waitUntil(lambda: not widget.worker_role, timeout=15000)
    widget.close()


def idle(qtbot, window):
    qtbot.waitUntil(lambda: (window.worker is None or not window.worker.isRunning())
                    and not window.worker_role and window.run_button.isEnabled(), timeout=20000)


# Four loci, because the identification index deliberately ignores a panel with
# fewer than four (identification._build_index), and this test types in auto mode.
LOCI = {"arcA": ARC, "gyrB": GYR, "mdh": "GGTTACGATCCAAT", "purA": "CCATTGACGGTACA"}


def make_scheme(path, st="17"):
    path.mkdir(parents=True)
    for locus, sequence in LOCI.items():
        (path / f"{locus}.tfa").write_text(f">{locus}_1\n{sequence}\n")
    (path / "profiles.tsv").write_text("ST\t" + "\t".join(LOCI) + f"\n{st}\t" +
                                       "\t".join("1" for _ in LOCI) + "\n")
    return path


def typeable(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{locus}\n{sequence}\n" for locus, sequence in LOCI.items()))
    return path


def assembly(path, seed=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">one{seed}\n{ARC}\n>two\n{GYR}\n")
    return path


def verdict(path, genus="", species="", *, confidence="unresolved", basis="none",
            quarantine_reason=None, **extra):
    """A verdict shaped exactly as organism_id.identify_input returns one."""
    record = {"format_version": 1, "input_path": str(path), "input_sha256": "a" * 64,
              "kind": "fasta", "basis": basis, "confidence": confidence, "status": "proposed",
              "proposed": {"genus": genus, "species": species},
              "accepted": {"genus": "", "species": ""}, "quarantine_reason": quarantine_reason,
              "confirmed_by": None, "confirmed_utc": None, "panel": {}, "detail": {},
              "margin_ani": None, "runner_up": None, "reason": "", "notes": [],
              "limitations": [], "errors": []}
    record.update(extra)
    return record


def accept_import(dialog):
    dialog.accept()
    return QDialog.DialogCode.Accepted


def test_identification_reads_the_original_files_before_any_copy_is_written(window, qtbot, tmp_path,
                                                                           monkeypatch):
    inbox = tmp_path / "inbox"
    first = assembly(inbox / "KPNIH1.fasta", "kp")
    second = assembly(inbox / "mystery.fasta", "mystery")
    duplicate = assembly(inbox / "also_KPNIH1.fasta", "kp")
    root = window.project_path.with_suffix(".files")
    seen = {}

    def identify(paths, **kwargs):
        seen["paths"] = list(paths)
        seen["copies"] = sorted(str(p) for p in root.rglob("*") if p.is_file())
        seen["scheme_paths"] = list(kwargs.get("scheme_paths") or ())
        return [verdict(first, "Klebsiella", "pneumoniae", basis="genomic_ani",
                        confidence="genomic_reference_supported"),
                verdict(duplicate, "Klebsiella", "pneumoniae", basis="genomic_ani",
                        confidence="genomic_reference_supported"),
                verdict(second, quarantine_reason="not_in_reference_panel")]

    reviews = []

    def review(dialog):
        reviews.append(dialog)
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr("wmlstudio.organism_id.identify_batch", identify)
    monkeypatch.setattr(ImportSamplesDialog, "exec", accept_import)
    monkeypatch.setattr(IdentificationReviewDialog, "exec", review)
    window.import_paths([str(inbox)], configure=True)
    idle(qtbot, window)

    assert seen["paths"] == [str(first), str(duplicate), str(second)]
    assert seen["copies"] == [], "identification must run before anything is copied"
    samples = {sample["name"]: sample for sample in window.project.samples()}
    assert set(samples) == {"KPNIH1", "mystery"}, "byte-identical files are imported once"
    assert "not copied again" in window.progress_text.text()
    assert duplicate.is_file(), "a skipped duplicate is the user's file and is never touched"
    filed = storage.Path(samples["KPNIH1"]["input_path"]).relative_to(root).parts
    assert filed[:3] == ("Klebsiella", "pneumoniae", "ST_unassigned")
    quarantined = storage.Path(samples["mystery"]["input_path"]).relative_to(root).parts
    assert quarantined[:2] == ("_Unresolved", "Not_in_reference_panel")
    assert first.read_text().startswith(">one"), "the user's original file is never moved"
    assert second.is_file()
    assert (root / "_Unresolved" / "README.txt").is_file()
    evidence = samples["KPNIH1"]["metadata"]["organism_evidence"]
    assert evidence["status"] == "confirmed" and evidence["confirmed_by"] == "user"
    assert evidence["basis"] == "genomic_ani"
    assert samples["KPNIH1"]["metadata"]["workflow"]["typing_mode"] == "auto", \
        "accepting an organism decides the folder, not how the typing scheme is chosen"
    assert window.project.get_sample(samples["mystery"]["id"])["metadata"]["organism"]["genus"] == ""
    assert window.test_errors == []


def test_needs_review_isolates_are_shown_as_a_branch_of_their_own(window, qtbot, tmp_path,
                                                                 monkeypatch):
    source = assembly(tmp_path / "inbox" / "unknown.fasta")
    monkeypatch.setattr("wmlstudio.organism_id.identify_batch",
                        lambda paths, **kwargs: [verdict(source,
                                                         quarantine_reason="not_in_reference_panel")])
    monkeypatch.setattr(ImportSamplesDialog, "exec", accept_import)
    monkeypatch.setattr(IdentificationReviewDialog, "exec", accept_import)
    window.import_paths([str(source)], configure=True)
    idle(qtbot, window)

    labels = [window.library_tree.topLevelItem(i).text(0)
              for i in range(window.library_tree.topLevelItemCount())]
    assert "Needs review  (1)" in labels
    branch = next(window.library_tree.topLevelItem(i)
                  for i in range(window.library_tree.topLevelItemCount())
                  if window.library_tree.topLevelItem(i).text(0) == "Needs review  (1)")
    assert branch.child(0).text(0) == "Not in reference panel  (1)"
    assert "declined to decide" in branch.toolTip(0)
    assert "Needs review" in window.selection_label.text()
    assert window.test_errors == []


def test_cancelling_the_review_imports_nothing_and_creates_no_folder(window, qtbot, tmp_path,
                                                                    monkeypatch):
    source = assembly(tmp_path / "inbox" / "isolate.fasta")
    monkeypatch.setattr("wmlstudio.organism_id.identify_batch",
                        lambda paths, **kwargs: [verdict(source, "Klebsiella", "pneumoniae",
                                                         basis="genomic_ani",
                                                         confidence="genomic_reference_supported")])
    monkeypatch.setattr(ImportSamplesDialog, "exec", accept_import)
    monkeypatch.setattr(IdentificationReviewDialog, "exec",
                        lambda dialog: QDialog.DialogCode.Rejected)
    window.import_paths([str(source)], configure=True)
    idle(qtbot, window)

    assert window.project.samples() == []
    assert not window.project_path.with_suffix(".files").exists()
    assert source.is_file()
    assert window.test_errors == []


def test_a_deferred_proposal_is_filed_as_your_choice_not_as_missing_evidence(window, qtbot,
                                                                            tmp_path, monkeypatch):
    source = assembly(tmp_path / "inbox" / "deferred.fasta")
    monkeypatch.setattr("wmlstudio.organism_id.identify_batch",
                        lambda paths, **kwargs: [verdict(source, "Klebsiella", "pneumoniae",
                                                         basis="genomic_ani",
                                                         confidence="genomic_reference_supported")])

    def review(dialog):
        dialog.quarantine_rest()
        dialog._set_action([0], "quarantine")
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ImportSamplesDialog, "exec", accept_import)
    monkeypatch.setattr(IdentificationReviewDialog, "exec", review)
    window.import_paths([str(source)], configure=True)
    idle(qtbot, window)

    sample = window.project.samples()[0]
    root = window.project_path.with_suffix(".files")
    assert storage.Path(sample["input_path"]).relative_to(root).parts[:2] == ("_Unresolved",
                                                                              "User_deferred")
    assert sample["metadata"]["organism_evidence"]["quarantine_reason"] == "user_deferred"
    assert window.test_errors == []


def test_an_accepted_organism_is_typed_and_then_refiled_under_its_st(window, qtbot, tmp_path,
                                                                     monkeypatch):
    scheme = make_scheme(window.root / "schemes" / "local_mlst")
    window.populate_schemes()
    # Only this panel is installed for the purposes of the run, so identification
    # is a real assignment rather than a search of every bundled scheme.
    window.scheme_paths = [scheme]
    source = typeable(tmp_path / "inbox" / "typed.fasta")
    monkeypatch.setattr("wmlstudio.organism_id.identify_batch",
                        lambda paths, **kwargs: [verdict(source, "Klebsiella", "pneumoniae",
                                                         basis="genomic_ani",
                                                         confidence="genomic_reference_supported")])
    monkeypatch.setattr(ImportSamplesDialog, "exec", accept_import)
    monkeypatch.setattr(IdentificationReviewDialog, "exec", accept_import)
    window.import_paths([str(source)], configure=True)
    idle(qtbot, window)

    window.start_analysis(all_samples=True)
    idle(qtbot, window)
    sample = window.project.samples()[0]
    root = window.project_path.with_suffix(".files")
    assert sample["result"]["st"] == "17", "an accepted organism still lets typing choose a scheme"
    assert storage.Path(sample["input_path"]).relative_to(root).parts[:3] == (
        "Klebsiella", "pneumoniae", "ST_17")
    assert source.is_file()
    assert window.test_errors == []


def review_dialog(qtbot, tmp_path, verdicts):
    dialog = IdentificationReviewDialog(verdicts, tmp_path / "managed")
    qtbot.addWidget(dialog)
    return dialog


def test_the_goes_to_column_is_computed_by_the_function_the_import_itself_calls(qtbot, tmp_path):
    root = tmp_path / "managed"
    verdicts = [
        verdict(tmp_path / "strong.fasta", "Klebsiella", "pneumoniae", basis="genomic_ani",
                confidence="genomic_reference_supported"),
        verdict(tmp_path / "complex.fasta", "Escherichia", "coli", basis="genomic_ani",
                confidence="complex_only"),
        verdict(tmp_path / "panel.fasta", "Enterobacter", basis="mlst_panel",
                confidence="panel_compatibility"),
        verdict(tmp_path / "none.fasta", quarantine_reason="not_in_reference_panel"),
    ]
    dialog = review_dialog(qtbot, tmp_path, verdicts)

    assert [dialog.row_action(row) for row in range(4)] == ["accept", "quarantine", "quarantine",
                                                            "quarantine"]
    expected = [
        storage.preview_target(root, dialog.PLACEHOLDER, tmp_path / "strong.fasta",
                               "Klebsiella", "pneumoniae"),
        storage.preview_target(root, dialog.PLACEHOLDER, tmp_path / "complex.fasta",
                               quarantine="user_deferred"),
        storage.preview_target(root, dialog.PLACEHOLDER, tmp_path / "panel.fasta",
                               quarantine="user_deferred"),
        storage.preview_target(root, dialog.PLACEHOLDER, tmp_path / "none.fasta",
                               quarantine="not_in_reference_panel"),
    ]
    for row, target in enumerate(expected):
        shown = dialog.table.item(row, 5).text()
        assert shown == target.parent.parent.relative_to(root.resolve()).as_posix() + "/"
        assert shown == dialog.destination_text(row)
    assert dialog.table.item(0, 4).text() == "Strong"
    assert dialog.table.item(2, 4).text() == "Panel match only"
    assert "not a species confirmation" in dialog.table.item(2, 4).toolTip()


def test_accepting_all_strong_proposals_leaves_every_weaker_one_for_a_person(qtbot, tmp_path):
    verdicts = [
        verdict(tmp_path / "a.fasta", "Klebsiella", "pneumoniae", basis="genomic_ani",
                confidence="genomic_reference_supported"),
        verdict(tmp_path / "b.fasta", "Escherichia", "coli", basis="genomic_ani",
                confidence="complex_only"),
        verdict(tmp_path / "c.fasta", "Enterobacter", basis="mlst_panel",
                confidence="panel_compatibility"),
    ]
    dialog = review_dialog(qtbot, tmp_path, verdicts)
    assert dialog.table.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu
    dialog._set_action([0, 1, 2], "quarantine")
    dialog.accept_strong()

    assert [dialog.row_action(row) for row in range(3)] == ["accept", "quarantine", "quarantine"]
    assert dialog.row_verdict(1)["status"] == "quarantined"
    dialog.accept_proposed()
    assert [dialog.row_action(row) for row in range(3)] == ["accept", "accept", "accept"]
    assert dialog.row_verdict(2)["confidence"] == "panel_compatibility", \
        "accepting a weak proposal must not strengthen its recorded confidence"


def test_editing_one_row_changes_only_that_row_and_drops_the_engine_basis(qtbot, tmp_path):
    verdicts = [
        verdict(tmp_path / "a.fasta", "Klebsiella", "pneumoniae", basis="genomic_ani",
                confidence="genomic_reference_supported"),
        verdict(tmp_path / "b.fasta", "Klebsiella", "pneumoniae", basis="genomic_ani",
                confidence="genomic_reference_supported"),
    ]
    dialog = review_dialog(qtbot, tmp_path, verdicts)
    unchanged = dialog.destination_text(0)
    dialog.rows[1][0].setCurrentText("Enterobacter")
    dialog.rows[1][1].setCurrentText("cloacae")

    assert dialog.destination_text(0) == unchanged
    assert dialog.destination_text(1).startswith("Enterobacter/cloacae/")
    assert dialog.table.item(1, 5).text().startswith("Enterobacter/cloacae/")
    corrected = dialog.row_verdict(1)
    assert corrected["basis"] == "user_assigned" and corrected["confidence"] == "unresolved"
    assert dialog.row_verdict(0)["basis"] == "genomic_ani"


def test_the_review_dialog_repeats_the_engines_own_sentences_and_the_folder_caveat(qtbot, tmp_path):
    limitation = ("ANI cannot reliably distinguish Escherichia coli from Shigella; report the "
                  "E. coli/Shigella complex.")
    note = "Nearest reference is within the configured species threshold."
    dialog = review_dialog(qtbot, tmp_path, [
        verdict(tmp_path / "a.fasta", "Escherichia", "coli", basis="genomic_ani",
                confidence="complex_only", reason="Nearest species meets the configured gates.",
                notes=[note], limitations=[limitation])])
    dialog.table.selectRow(0)
    shown = dialog.evidence.toPlainText()

    assert limitation in shown and note in shown
    assert "Nearest species meets the configured gates." in shown
    assert dialog.FOOTER == ("A folder name is where the file is stored. It is not a laboratory "
                             "identification. Confirm the organism before clinical interpretation.")
    texts = [child.text() for child in dialog.findChildren(type(dialog.feedback))]
    assert dialog.FOOTER in texts


def test_the_review_dialog_fits_a_small_screen_without_horizontal_scrolling(qtbot, tmp_path):
    verdicts = [verdict(tmp_path / f"{index}.fasta", "Klebsiella", "pneumoniae",
                        basis="genomic_ani", confidence="genomic_reference_supported")
                for index in range(12)]
    dialog = review_dialog(qtbot, tmp_path, verdicts)
    dialog.resize(1080, 720)
    dialog.show()
    qtbot.waitUntil(lambda: dialog.table.width() > 0, timeout=5000)

    assert dialog.minimumSizeHint().width() <= 1080
    assert dialog.minimumSizeHint().height() <= 720
    assert dialog.table.horizontalScrollBar().maximum() == 0
    dialog.hide()


def test_the_practice_download_dialog_repeats_every_caveat_word_for_word(qtbot, tmp_path):
    cohorts = practice_cohorts.describe_cohorts()
    destinations = {entry["name"]: tmp_path / entry["name"] for entry in cohorts}
    dialog = PracticeCohortDialog(cohorts, destinations)
    qtbot.addWidget(dialog)

    for index, entry in enumerate(cohorts):
        dialog.list.setCurrentRow(index)
        shown = dialog.detail.toPlainText()
        assert dialog.chosen == entry["name"]
        for caveat in entry["caveats"]:
            assert caveat in shown
        assert str(destinations[entry["name"]]) in shown
        assert "not a validation set" in shown


def test_downloading_practice_data_runs_in_the_cancellable_worker(window, qtbot, tmp_path,
                                                                  monkeypatch):
    destination = tmp_path / "practice"
    destination.mkdir()
    with gzip.open(destination / "Klebsiella_pneumoniae_HS11286_GCF_000240185.1.fna.gz", "wt") as handle:
        handle.write(f">one\n{ARC}\n")
    seen = {}

    def download(name, path, *, cancelled=None, progress=None, **kwargs):
        seen["name"] = name
        seen["cancellable"] = callable(cancelled) and callable(progress)
        progress(1, 1, "Verifying published checksums…")
        return {"path": str(destination), "cohort": name, "genomes": 1, "downloaded": 1,
                "reused": 0, "total_bytes": 1, "installed": True}

    monkeypatch.setattr(practice_cohorts, "download_cohort", download)
    monkeypatch.setattr(PracticeCohortDialog, "exec", accept_import)
    monkeypatch.setattr(ImportSamplesDialog, "exec", lambda dialog: QDialog.DialogCode.Rejected)
    window.download_practice_cohort()
    idle(qtbot, window)

    assert seen["name"] in practice_cohorts.cohort_names()
    assert seen["cancellable"], "download must run through the worker's cancel and progress hooks"
    assert str(destination) in window.progress_text.text()
    assert window.project.samples() == [], "the download itself imports nothing"
    assert window.test_errors == []
