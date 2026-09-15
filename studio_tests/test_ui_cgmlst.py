"""The dedicated cgMLST page: the table of calls, its denominator, and the read recheck."""

import json

import pytest
from PySide6.QtCore import Qt

from wmlstudio import ui_cgmlst
from wmlstudio.app import MainWindow
from wmlstudio.cgmlst_calls import (
    CALL_STATES,
    RECHECK_LIMITATIONS,
    TABLE_LIMITATIONS,
    interpret_target_recheck,
    plan_target_recheck,
)
from wmlstudio.ui_cgmlst import CgmlstCallsPanel

DIGEST = "c" * 64
OTHER_DIGEST = "e" * 64
# Above project.CGMLST_LOCUS_FLOOR, so a stored profile of this size is classified as
# core-genome typing on its own evidence rather than on anything this test asserts.
LOCI = [f"target{index:03d}" for index in range(1, 41)]


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.resize(1280, 860)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


@pytest.fixture
def panel(window):
    return window.cgmlst_calls


def call_row(locus, status, allele=None):
    return {"locus": locus, "status": status, "allele": allele, "reason": f"{status} evidence"}


def cgmlst_result(states=None, *, scheme="Demo cgMLST", digest=DIGEST):
    """A stored cgMLST analysis: every target called exactly unless `states` says otherwise."""
    states = states or {}
    calls = []
    for index, locus in enumerate(LOCI, 1):
        status, allele = states.get(locus, ("exact", str(index)))
        calls.append(call_row(locus, status, allele))
    return {"scheme": scheme, "scheme_digest": digest, "status": "incomplete",
            "input_sha256": "a" * 64, "analysis_kind": "cgmlst",
            "alleles": {row["locus"]: row["allele"] for row in calls}, "calls": calls}


def typed(window, name, states=None, **kwargs):
    sid = window.project.add_profile(name, dict(cgmlst_result(states, **kwargs), sample_name=name),
                                     {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}})
    window.cohort_ids = set(window.cohort_ids) | {sid}
    return sid


def attach_reads(window, sample_id):
    """A verified original FASTQ pair, recorded exactly as the read attachment records it."""
    window.project.update_metadata(sample_id, {"reads": {"reads": [
        {"mate": 1, "path": "/reads/R1.fastq.gz", "sha256": "1" * 64},
        {"mate": 2, "path": "/reads/R2.fastq.gz", "sha256": "2" * 64}]}})


def build(panel, qtbot, window):
    assert panel.refresh_references()
    assert panel.build_table() is not False
    qtbot.waitUntil(lambda: panel.calls is not None and not window.worker_role, timeout=20000)
    assert window.test_errors == []
    return panel.calls


def support_locus(locus, status, compatible=()):
    return {"locus": locus, "status": status, "candidates": [], "assigned_allele": None,
            "compatible_candidates": list(compatible)}


def support_payload(loci, *, complete=True):
    return {"format_version": 1, "status": "completed", "scheme": "Demo cgMLST",
            "scheme_digest": DIGEST, "loci": list(loci),
            "sampling": {"complete_files": complete, "sampled": not complete, "pairs_checked": 5000},
            "panel": {"loci": [row["locus"] for row in loci]}, "parameters": {"min_depth": 3},
            "reads": [], "input_sha256": "a" * 64, "provenance": {}}


def test_the_table_of_calls_states_the_schemes_full_target_count_beside_every_number(
        panel, window, qtbot):
    typed(window, "Isolate A")
    typed(window, "Isolate B", {"target003": ("missing", None), "target004": ("partial", None)})
    table = build(panel, qtbot, window)
    assert table["target_count"] == len(LOCI) and table["isolate_count"] == 2
    assert panel.model.rowCount() == len(LOCI)
    assert panel.model.columnCount() == 3  # the target column, then one per isolate
    assert f"{len(LOCI)} targets" in panel.denominator.text()
    assert f"Showing {len(LOCI)} of {len(LOCI)} targets" in panel.denominator.text()
    summaries = [panel.isolate_table.item(row, 1).text() for row in range(panel.isolate_table.rowCount())]
    assert sorted(summaries) == [f"38 of {len(LOCI)} targets called (95.0% of the scheme)",
                                 f"40 of {len(LOCI)} targets called (100.0% of the scheme)"]
    assert panel.limitations.text() == "\n".join(TABLE_LIMITATIONS)
    assert f"{len(LOCI)} targets × 2 isolate(s)" in panel.status.text()


def test_every_reason_there_is_no_allele_keeps_its_own_cell_and_never_reads_as_a_call(
        panel, window, qtbot):
    typed(window, "Isolate A", {
        "target002": ("novel_validated", "NOVEL_" + "b" * 64),
        "target003": ("missing", None), "target004": ("partial", None),
        "target005": ("low_similarity", None), "target006": ("ambiguous", None),
        "target007": ("mixed", None), "target008": ("a_state_from_the_future", None)})
    build(panel, qtbot, window)
    states, displays = {}, {}
    for row in range(panel.model.rowCount()):
        index = panel.model.index(row, 1)
        call = panel.model.cell_at(row, 1)
        states[panel.model.data(panel.model.index(row, 0))] = call["state"]
        displays[call["state"]] = panel.model.data(index)
    assert states["target001"] == "called" and states["target002"] == "called_novel"
    assert [states[locus] for locus in ("target003", "target004", "target005", "target006",
                                        "target007", "target008")] == [
        "missing", "partial", "low_similarity", "ambiguous", "mixed", "unknown"]
    # No state that assigns no allele may ever render as one, and the eight states
    # that are on screen are eight distinct pieces of text.
    for state, text in displays.items():
        if not CALL_STATES[state]["has_allele"]:
            assert text == CALL_STATES[state]["label"]
    assert len(set(displays.values())) == len(displays)
    tooltip = panel.model.data(panel.model.index(2, 1), Qt.ItemDataRole.ToolTipRole)
    assert CALL_STATES["missing"]["meaning"] in tooltip and "target003" in tooltip
    assert CALL_STATES["missing"]["meaning"] in panel.legend.toolTip()


def test_filtering_the_rows_never_moves_the_denominator_under_them(panel, window, qtbot):
    typed(window, "Isolate A", {"target003": ("missing", None)})
    typed(window, "Isolate B", {"target003": ("missing", None), "target009": ("partial", None)})
    build(panel, qtbot, window)
    panel.row_filter.setCurrentIndex(1)
    assert [row["locus"] for row in panel.visible_rows()] == ["target003", "target009"]
    assert panel.model.rowCount() == 2
    assert f"{len(LOCI)} targets this scheme defines" in panel.denominator.text()
    assert f"Showing 2 of {len(LOCI)} targets" in panel.denominator.text()
    panel.row_filter.setCurrentIndex(2)
    assert [row["locus"] for row in panel.visible_rows()] == ["target003"]
    panel.row_filter.setCurrentIndex(0)
    panel.locus_search.setText("target01")
    assert panel.model.rowCount() == 10
    assert f"Showing 10 of {len(LOCI)} targets" in panel.denominator.text()


def test_two_references_stay_two_tables_and_are_never_tabulated_together(panel, window, qtbot):
    first, second = typed(window, "Isolate A"), typed(window, "Isolate B")
    other = typed(window, "Isolate C", scheme="Other cgMLST", digest=OTHER_DIGEST)
    references = panel.refresh_references()
    assert [entry["scheme_digest"] for entry in references] == [DIGEST, OTHER_DIGEST]
    assert [entry["isolates"] for entry in references] == [2, 1]
    table = build(panel, qtbot, window)
    assert table["scheme_digest"] == DIGEST
    assert {row["sample_id"] for row in table["samples"]} == {first, second}
    assert other not in {row["sample_id"] for row in table["samples"]}
    panel.reference.setCurrentIndex(1)
    second = build(panel, qtbot, window)
    assert second["scheme_digest"] == OTHER_DIGEST and second["isolate_count"] == 1


def test_an_isolate_with_no_core_genome_profile_is_absent_rather_than_shown_uncalled(
        panel, window, qtbot):
    typed(window, "Isolate A")
    classical = window.project.add_profile("Classical only", {
        "sample_name": "Classical only", "scheme": "Study MLST", "scheme_digest": "m" * 64,
        "status": "complete", "alleles": dict(zip("abcdefg", "1111111")), "calls": [],
        "input_sha256": "a" * 64})
    window.cohort_ids = set(window.cohort_ids) | {classical}
    table = build(panel, qtbot, window)
    assert [row["sample_name"] for row in table["samples"]] == ["Isolate A"]
    assert classical not in {row["sample_id"] for row in table["samples"]}


def test_the_page_reads_the_comparison_cohort_and_can_be_widened_to_the_project(panel, window):
    inside = typed(window, "Isolate A")
    outside = typed(window, "Isolate B")
    window.cohort_ids = {inside}
    assert [entry["isolates"] for entry in panel.refresh_references()] == [1]
    panel.scope.setCurrentIndex(1)
    assert [entry["isolates"] for entry in panel.references] == [2]
    assert outside in {row["sample_id"] for row in panel.references[0]["records"]}


def test_clearing_the_page_empties_it_and_leaves_every_stored_profile_in_place(
        panel, window, qtbot):
    typed(window, "Isolate A", {"target003": ("missing", None)})
    build(panel, qtbot, window)
    assert panel.model.rowCount() == len(LOCI)
    panel.clear()
    assert panel.calls is None and panel.recheck is None and panel.plan is None
    assert panel.model.rowCount() == 0 and panel.model.columnCount() == 0
    assert panel.isolate_table.rowCount() == 0 and panel.recheck_table.rowCount() == 0
    assert panel.legend.text() == "" and panel.limitations.text() == ""
    assert "Nothing stored changed" in panel.status.text()
    assert len(window.project.samples()) == 1
    assert window.project.samples()[0]["result"]["alleles"]["target003"] is None
    assert build(panel, qtbot, window)["target_count"] == len(LOCI)


def test_a_recheck_without_a_verified_read_pair_is_refused_in_the_users_own_words(
        panel, window, qtbot):
    typed(window, "Isolate A", {"target003": ("missing", None)})
    build(panel, qtbot, window)
    assert panel.recheck_missing() is not True
    assert "Attach a verified original FASTQ pair" in panel.status.text()
    assert panel.recheck is None and not window.worker_role


def test_a_recheck_runs_the_bounded_assay_and_reports_it_as_evidence_not_as_a_call(
        panel, window, qtbot, monkeypatch, tmp_path):
    sid = typed(window, "Isolate A", {"target003": ("missing", None), "target004": ("partial", None)})
    attach_reads(window, sid)
    table = build(panel, qtbot, window)
    plan = plan_target_recheck(table, sid, sample=window.project.get_sample(sid), control_count=2)
    assert plan["targets"] == ["target003", "target004"] and len(plan["controls"]) == 2
    interpreted = interpret_target_recheck(support_payload(
        [support_locus("target003", "supported", compatible=["1", "7"]),
         support_locus("target004", "no_support"),
         *[support_locus(locus, "supported") for locus in plan["controls"]]]), plan, table=table)
    folder = tmp_path / "reference"
    folder.mkdir()
    monkeypatch.setattr(ui_cgmlst, "cached_scheme", lambda path, cancelled=None: object())
    monkeypatch.setattr(ui_cgmlst, "recheck_missing_targets",
                        lambda *args, **kwargs: interpreted)
    panel.scheme_folder.addItem("Demo cgMLST", str(folder))
    panel.scheme_folder.setCurrentIndex(panel.scheme_folder.count() - 1)
    assert panel.recheck_missing() is True
    qtbot.waitUntil(lambda: panel.recheck is not None and not window.worker_role, timeout=20000)
    assert window.test_errors == []
    rows = {panel.recheck_table.item(row, 0).text():
            [panel.recheck_table.item(row, column).text() for column in range(1, 6)]
            for row in range(panel.recheck_table.rowCount())}
    assert set(rows) == {"target003", "target004"}
    assert rows["target003"][0] == "Present in the reads" and rows["target003"][1] == "present"
    assert rows["target004"][0] == "Consistent with absence" and rows["target004"][1] == "absent"
    assert rows["target003"][2] == CALL_STATES["missing"]["label"]
    # Compatible reference alleles are never offered as a call, in any column.
    assert rows["target003"][3] == "2 compatible, none assigned"
    assert all(row["assigned_allele"] is None for row in panel.recheck["targets"])
    assert panel.recheck_headline.text() == interpreted["headline"]
    assert "control target" in panel.recheck_controls.text()
    for sentence in RECHECK_LIMITATIONS:
        assert sentence in panel.recheck_limits.text()
    # Nothing the reads showed is written back onto the stored profile.
    assert window.project.get_sample(sid)["result"]["alleles"]["target003"] is None


def test_a_recheck_names_the_targets_it_deferred_rather_than_dropping_them(
        panel, window, qtbot, monkeypatch, tmp_path):
    missing = {locus: ("missing", None) for locus in LOCI[:30]}
    sid = typed(window, "Isolate A", missing)
    attach_reads(window, sid)
    table = build(panel, qtbot, window)
    plan = plan_target_recheck(table, sid, sample=window.project.get_sample(sid), control_count=2)
    assert len(plan["targets"]) + len(plan["controls"]) == 20
    assert len(plan["deferred_targets"]) == 30 - len(plan["targets"])
    notices = []
    monkeypatch.setattr(type(window), "notify", lambda self, message: notices.append(str(message)))
    interpreted = interpret_target_recheck(support_payload(
        [support_locus(locus, "no_support") for locus in plan["targets"]]
        + [support_locus(locus, "supported") for locus in plan["controls"]]), plan, table=table)
    monkeypatch.setattr(ui_cgmlst, "cached_scheme", lambda path, cancelled=None: object())
    monkeypatch.setattr(ui_cgmlst, "recheck_missing_targets", lambda *args, **kwargs: interpreted)
    folder = tmp_path / "reference"
    folder.mkdir()
    panel.scheme_folder.addItem("Demo cgMLST", str(folder))
    panel.scheme_folder.setCurrentIndex(panel.scheme_folder.count() - 1)
    assert panel.recheck_missing() is True
    qtbot.waitUntil(lambda: panel.recheck is not None and not window.worker_role, timeout=20000)
    assert any("are not checked until you run them" in notice for notice in notices)
    assert "have not been checked" in panel.recheck_controls.text()
    assert len(panel.recheck["deferred_targets"]) == len(plan["deferred_targets"])


def test_the_page_says_so_when_nothing_in_the_selection_has_a_core_genome_profile(panel, window):
    window.project.add_profile("Classical only", {
        "sample_name": "Classical only", "scheme": "Study MLST", "scheme_digest": "m" * 64,
        "status": "complete", "alleles": dict(zip("abcdefg", "1111111")), "calls": []})
    panel.scope.setCurrentIndex(1)
    assert panel.references == []
    assert not panel.build_button.isEnabled()
    assert "unknown evidence here" in panel.status.text()
    assert panel.build_table() is not True


def test_a_panel_needs_a_workspace_window_and_says_which_part_is_missing(qtbot):
    class Incomplete:
        project = None

    with pytest.raises(TypeError, match="launch_task"):
        CgmlstCallsPanel(Incomplete())


def installed_reference(window, folder_name, scheme_name, loci=2):
    """An installed cgMLST folder, declared as such in its own metadata."""
    folder = window.root / "schemes" / folder_name
    folder.mkdir(parents=True)
    (folder / "scheme.json").write_text(
        json.dumps({"name": scheme_name, "type": "cgmlst", "organism": "Klebsiella pneumoniae"}))
    for locus in LOCI[:loci]:
        (folder / f"{locus}.tfa").write_text(f">{locus}_1\nACGTACGTACGT\n")
    window.scheme_paths = [*window.scheme_paths, folder]
    return folder


def test_the_recheck_offers_installed_references_and_preselects_the_matching_one(
        panel, window, qtbot):
    demo = installed_reference(window, "demo", "Demo cgMLST")
    installed_reference(window, "other", "Another cgMLST")
    typed(window, "Isolate A", {"target003": ("missing", None)})
    build(panel, qtbot, window)
    titles = [panel.scheme_folder.itemText(index) for index in range(panel.scheme_folder.count())]
    assert len(titles) == 2 and any("Demo cgMLST" in title for title in titles)
    assert panel.scheme_folder.currentData() == str(demo)
    assert panel.recheck_button.isEnabled()


def test_the_calls_page_hands_its_isolates_to_the_cgmlst_tree_on_its_own_scale(
        panel, window, qtbot):
    first, second = typed(window, "Isolate A"), typed(window, "Isolate B")
    window.cohort_ids = {first, second}
    table = build(panel, qtbot, window)
    assert panel.tree_button.isVisibleTo(panel) and panel.tree_button.isEnabled()
    panel.show_tree()
    assert window.typing_kind == "cgmlst"
    assert window.cohort_ids == {first, second}
    assert sorted(window.tree.nodes) == sorted([first, second])
    assert window.graph_tabs.currentWidget() is window.graph_split
    # The table's own quantity travels with the tree; nothing else does.
    assert window.tree.scale["targets"] == table["target_count"]
    assert window.tree.scale["kind"] == "cgmlst"
    assert "share no scale" in panel.status.text()
