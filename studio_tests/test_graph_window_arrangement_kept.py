"""An arrangement made in a detached tree window survives closing that window.

The window is a copy of the page's forest, so for a long time nothing it was
arranged into ever came back: the picture somebody had spent ten minutes laying
out was thrown away by the close button. These tests hold the arrangement to the
tree it was opened from, to that tree's own project setting, and to the sentence
the window prints about what is kept and what is not.
"""

import pytest

from wmlstudio.graph_window import EDIT_NOTE, GraphIdentity, GraphWindow
from wmlstudio.investigation import InvestigationStore
from wmlstudio.widgets import TreeView


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "application")
    errors = []
    monkeypatch.setattr(widget, "error", lambda message: errors.append(str(message)))
    widget.test_errors = errors
    qtbot.addWidget(widget)
    # Laid out with resize, never shown: repainting a second TreeView after many
    # torn-down offscreen widgets crashes inside the graph label painter.
    widget.resize(1280, 860)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


def imported(window, name, vector="1111", ward="ICU"):
    return window.project.add_profile(name, {
        "sample_name": name, "scheme": "Study MLST", "scheme_digest": "study-reference",
        "status": "profile_imported", "alleles": dict(zip("abcd", vector)), "calls": [],
        "st": "20" if vector.startswith("1") else "44", "input_sha256": "a" * 64,
    }, {"organism": {"genus": "Staphylococcus", "species": "aureus"}, "annotations": {"ward": ward}})


def build(window, ids, investigation_id=None, name="Ward A"):
    """Save (or re-save) the investigation and build it, which stores a snapshot."""
    window.refresh_cohort_table()
    plan = InvestigationStore(window.project).save(
        name, list(ids), investigation_id=investigation_id, scheme="Study MLST",
        scheme_digest="study-reference", threshold=1, min_overlap=0.95,
        protocol="Synthetic test protocol; not a clinical cutoff")
    window.cohort_ids = set(ids)
    window.select_investigation(plan["id"])
    assert window._current_snapshot, window.tree_status.text()
    return plan["id"]


def grown_investigation(window):
    """A baseline of two isolates and a current comparison of three."""
    a, b = imported(window, "A"), imported(window, "B", "2111")
    investigation_id = build(window, [a, b])
    c = imported(window, "C", "2211", ward="Ward 2")
    build(window, [a, b, c], investigation_id=investigation_id)
    window.dual_toggle.setChecked(True)
    return investigation_id, a, b, c


def test_arranging_a_detached_window_and_closing_it_leaves_the_page_and_the_project_arranged(window):
    """Prevents the reported loss: an arrangement made in the popup was discarded on close.

    Nothing listened to the window's stateChanged, so ten minutes of laying out a
    forest ended at the close button. The page's own tree and the project setting it
    already saves to must both carry what the window was left showing.
    """
    _, a, b, c = grown_investigation(window)
    detached = window.open_graph_window('current')
    detached.view.nodes[a].setPos(211, 132)
    detached.view.set_node_colors({b: "#ff0055"})
    detached.view.set_node_label(c, "Ward 2 index case")
    detached.close()
    assert window.tree.nodes[a].pos().x() == pytest.approx(211)
    assert window.tree.nodes[a].pos().y() == pytest.approx(132)
    assert window.tree.nodes[b].slices == [("#ff0055", 1)]
    assert window.tree.node_label(c) == "Ward 2 index case"
    stored = window.project.get_setting("graph_style", {})
    assert stored["positions"][a] == [211.0, 132.0]
    assert stored["colors"][b] == "#ff0055"
    assert stored["labels"][c] == "Ward 2 index case"
    # Presentation only: the stored sample name and the drawn evidence are untouched.
    assert [row["name"] for row in window.project.samples() if row["id"] == c] == ["C"]
    assert set(window.tree._results) == {a, b, c}


def test_opening_the_window_a_second_time_starts_from_the_arrangement_the_first_one_left(window):
    """Prevents a fix that saves the arrangement but still reopens on the automatic layout.

    A window that forgets what the last window was arranged into is the same
    complaint one step later.
    """
    _, a, _b, _c = grown_investigation(window)
    first = window.open_graph_window('current')
    first.view.nodes[a].setPos(305, 44)
    first.close()
    second = window.open_graph_window('current')
    try:
        assert second.view.nodes[a].pos().x() == pytest.approx(305)
        assert second.view.nodes[a].pos().y() == pytest.approx(44)
    finally:
        second.close()


def test_an_arrangement_made_on_the_baseline_window_never_lands_on_the_current_tree(window):
    """Prevents one stored arrangement being written over another tree's.

    The baseline is a frozen snapshot of a different cohort. Its layout going into
    the current tree's setting would move nodes that stand for other isolates.
    """
    _, a, _b, _c = grown_investigation(window)
    window.tree.nodes[a].setPos(10, 20)
    window.tree.update_edges()
    detached = window.open_graph_window('baseline')
    detached.view.nodes[a].setPos(400, 410)
    detached.close()
    assert window.baseline_tree.nodes[a].pos().x() == pytest.approx(400)
    assert window.tree.nodes[a].pos().x() == pytest.approx(10)
    assert window.project.get_setting("graph_style.baseline", {})["positions"][a] == [400.0, 410.0]
    assert window.project.get_setting("graph_style", {})["positions"][a] == [10.0, 20.0]


def test_the_detached_window_says_the_arrangement_is_kept_and_that_zoom_is_not(window):
    """Prevents a silent save-back, and prevents implying the zoom comes back too.

    Zoom and pan are not part of the exported presentation state, so a window that
    only said "changes are kept" would be promising something it does not carry.
    """
    grown_investigation(window)
    detached = window.open_graph_window('current')
    try:
        note = detached.keep_note.text()
        assert detached.keep_note.isVisibleTo(detached) is True
        assert "kept" in note and "current tree on the comparison page" in note
        assert "zoom and pan do not" in note
        assert "Allele calls, distances and group membership are untouched" in note
        assert note in detached.statusBar().currentMessage()
    finally:
        detached.close()


def test_a_window_nobody_is_keeping_promises_nothing(qtbot):
    """Prevents the window claiming a save when no host connected stateChanged.

    The same window class is opened from pages that do not write the arrangement
    anywhere, and an unperformed save must not be announced as a performed one.
    """
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(900, 600)
    window = GraphWindow(GraphIdentity(kind="mlst", title="Classical MLST"), view=view)
    qtbot.addWidget(window)
    assert window.keep_note.text() == ""
    assert window.keep_note.isVisibleTo(window) is False
    assert window.statusBar().currentMessage() == EDIT_NOTE
    window.keep_arrangement_in("the current tree on the comparison page")
    assert "kept" in window.keep_note.text()
    assert window.keep_note.isVisibleTo(window) is True
