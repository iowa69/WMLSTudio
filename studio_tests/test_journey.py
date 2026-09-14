"""Workflow maps stay scoped and preserve unknown/stale evidence states."""

import copy

from PySide6.QtCore import QUrl

from wmlstudio.journey import characterization_state, journey_summary
from wmlstudio.workflow_guide import WorkflowGuide, guide_chapters, guide_path


def sample(identifier, *, status="completed", result=None, metadata=None):
    return {"id": identifier, "name": identifier, "input_path": "/not/opened.fasta",
            "status": status, "result": result, "metadata": metadata or {}}


def test_map_empty_does_not_invent_analysis():
    result = journey_summary([])
    assert result["counts"]["total"] == 0
    assert len(result["steps"]) == 6
    assert result["steps"][0].state == "empty"
    assert result["counts"].get("typed", 0) == 0


def test_map_counts_isolates_not_mates_and_preserves_scope():
    records = [sample("a", result={"status": "complete", "st": "7", "alleles": {"a": "1"}}),
               sample("b", status="queued"), sample("c", status="failed"),
               sample("mate", metadata={"workflow": {"source_kind": "read_mate"}})]
    original = copy.deepcopy(records)
    result = journey_summary(records, selected_ids={"a", "mate"}, comparison_ids={"a", "b"})
    assert result["counts"]["total"] == 3
    assert result["counts"]["linked_mates"] == 1
    assert result["counts"]["selected"] == 1
    assert result["counts"]["comparison"] == 2
    assert result["counts"]["typed"] == 1
    assert result["pending_ids"] == ["b"]
    assert result["review_ids"] == ["c"]
    assert records == original


def test_characterization_stale_and_hashless_are_not_current():
    record = sample("a", result={"input_sha256": "a" * 64}, metadata={
        "characterization": {"input_sha256": "b" * 64, "status": "completed"}})
    assert characterization_state(record)[0] == "stale"
    assert journey_summary([record])["counts"]["characterized"] == 0
    record["metadata"]["characterization"].pop("input_sha256")
    assert characterization_state(record)[0] == "unverified"
    record["metadata"]["characterization"]["input_sha256"] = "a" * 64
    assert characterization_state(record)[0] == "completed"
    record["result"] = None
    assert characterization_state(record)[0] == "unverified"


def test_offline_guide_has_all_problem_topics():
    text = guide_path().read_text(encoding="utf-8")
    chapters = guide_chapters(text)
    assert len(chapters) >= 15
    assert any("New samples" in title for title, _ in chapters)
    assert any("Missing targets" in title for title, _ in chapters)
    assert "not measured antimicrobial susceptibility" in text


def test_guide_search_and_safe_action_routes(qtbot, monkeypatch):
    guide = WorkflowGuide()
    qtbot.addWidget(guide)
    actions = []
    opened = []
    guide.actionRequested.connect(actions.append)
    monkeypatch.setattr("wmlstudio.workflow_guide.QDesktopServices.openUrl", opened.append)
    guide.search.setText("plasmid")
    assert 0 < guide.contents.count() < len(guide.chapters)
    guide.search.setText("a-topic-that-is-not-present")
    assert guide.contents.count() == 0
    guide.follow_link(QUrl("file:///private/sequence.fasta"))
    guide.follow_link(QUrl("wmlstudio:delete_everything"))
    assert not actions and not opened
    guide.follow_link(QUrl("wmlstudio:compare"))
    assert actions == ["compare"]
    guide.close()


def test_map_connected_to_native_navigation(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    window = MainWindow(storage_root=tmp_path)
    qtbot.addWidget(window)
    assert window.overview_tabs.tabText(0) == "Investigation map"
    # The map routes by index; the workspace now also names each page, so a
    # renumbering would be caught here instead of silently landing elsewhere.
    for action, index, key in [("samples", 1, "isolates"), ("compare", 2, "compare"),
                               ("reports", 5, "reports")]:
        window.journey_action(action)
        assert window.pages.currentIndex() == index
        assert window.pages.current_key() == key
        assert window.page_index[key] == index
    window.open_workflow_guide("plasmid")
    assert window.workflow_guide.isVisible()
    window.workflow_guide.close()
    window.close()


def test_stage_button_labels_do_not_turn_an_ampersand_into_a_keyboard_shortcut():
    """Qt reads '&' in button text as a mnemonic, so 'samples & reads' drew as 'samples _reads'.

    The investigation map is the first thing a new user reads, so a mangled word
    there costs more trust than it looks like it should.
    """
    from PySide6.QtWidgets import QPushButton

    from wmlstudio.journey import journey_summary
    for step in journey_summary([])["steps"]:
        # Only the action label becomes a QPushButton (journey_widgets.py). The
        # titles are QLabels, which render '&' literally, so they stay unescaped.
        assert "&" not in step.action_label.replace("&&", ""), \
            f"unescaped ampersand in button label {step.action_label!r}"
        # Qt strips the escape when it renders, so the user reads a real ampersand.
        assert QPushButton(step.action_label).text() == step.action_label
