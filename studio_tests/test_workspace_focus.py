"""Shared focus connects the tabs; only an explicit click ever writes a cohort."""

import pytest
from PySide6.QtWidgets import QWidget

from wmlstudio.workspace_focus import (
    COHORT_ATTRIBUTES,
    EMPTY_COHORT_TEXT,
    FOCUS_ACCENT,
    FOCUS_LEGEND,
    SAVED_HIGHLIGHT_ACCENT,
    CohortBar,
    CohortLedger,
    FocusBar,
    FocusBus,
    adopt_focus,
    cohort_bar,
    cohort_ids,
    send_selection,
)


class FakeProject:
    """Only the two calls workspace_focus is allowed to make on a project."""

    def __init__(self):
        self.settings = {}

    def set_setting(self, key, value):
        self.settings[key] = value


class FakeWindow(QWidget):
    """The window duck-type these primitives are wired into, without importing the app."""

    def __init__(self):
        super().__init__()
        self.selection_ids = set()
        self.cohort_ids = set()
        self.feature_ids = set()
        self.report_ids = set()
        self.project = FakeProject()
        self.focus = FocusBus(self)
        self.cohort_origins = CohortLedger(self)
        self.page_index = {"overview": 0, "isolates": 1, "compare": 2, "schemes": 3,
                           "evidence": 4, "reports": 5, "settings": 6}
        self.messages = []
        self.visited = []
        self.refreshed = []

    def notify(self, message):
        self.messages.append(message)

    def navigate(self, index):
        self.visited.append(index)

    def refresh_tables(self):
        self.refreshed.append("isolates")

    def refresh_comparison(self):
        self.refreshed.append("compare")

    def refresh_features(self):
        self.refreshed.append("evidence")

    def refresh_report_table(self):
        self.refreshed.append("reports")

    def choose_comparison_cohort(self):
        self.refreshed.append("chooser")

    def cohorts(self):
        return {key: set(getattr(self, attribute))
                for key, attribute in COHORT_ATTRIBUTES.items()}


@pytest.fixture
def window(qtbot):
    widget = FakeWindow()
    qtbot.addWidget(widget)
    return widget


def test_focusing_isolates_in_one_tab_never_widens_another_tabs_cohort(window):
    """The honesty invariant: focus is a highlight, not a cohort."""
    changes = []
    window.focus.focusChanged.connect(lambda: changes.append(set(window.focus.ids)))
    window.focus.set_focus(["a", "b", "c"], "Isolate table selection")
    assert window.focus.ids == {"a", "b", "c"}
    assert window.focus.origin == "Isolate table selection"
    assert window.focus.at
    assert changes == [{"a", "b", "c"}]
    assert window.cohorts() == {"isolates": set(), "compare": set(),
                                "evidence": set(), "reports": set()}
    assert window.project.settings == {}, "focus is never persisted"
    assert window.refreshed == []


def test_repeating_the_same_focus_does_not_churn_the_interface(window):
    changes = []
    window.focus.focusChanged.connect(lambda: changes.append(1))
    window.focus.set_focus(["a"], "Graph selection")
    window.focus.set_focus(["a"], "Graph selection")
    assert len(changes) == 1
    window.focus.set_focus(["a"], "Group: Ward A")
    assert len(changes) == 2


def test_clearing_focus_forgets_the_isolates_and_where_they_came_from(window):
    window.focus.set_focus(["a"], "Graph selection")
    window.focus.clear()
    assert window.focus.ids == set()
    assert window.focus.origin == ""
    assert window.focus.at == ""
    assert window.focus.describe() == "No isolates focused."


def test_focus_drops_isolates_that_are_no_longer_in_the_project(window):
    window.focus.set_focus(["a", "b", "c"], "Isolate table selection")
    window.focus.prune({"a", "c"})
    assert window.focus.ids == {"a", "c"}
    window.focus.prune(set())
    assert window.focus.ids == set()
    assert window.focus.origin == ""


def test_adopting_focus_writes_one_cohort_and_records_where_it_came_from(window):
    window.focus.set_focus(["a", "b"], "Graph selection")
    adopted = adopt_focus(window, "reports")
    assert adopted == {"a", "b"}
    assert window.report_ids == {"a", "b"}
    assert window.cohort_ids == set() and window.feature_ids == set()
    assert window.refreshed == ["reports"]
    entry = window.cohort_origins.entry("reports")
    assert entry["origin"] == "Graph selection"
    assert entry["count"] == 2
    described = window.cohort_origins.describe("reports", 2)
    assert "Graph selection" in described
    assert "belongs to this tab only" in described
    assert "Graph selection" in window.messages[-1]


def test_adopting_nothing_reports_it_and_leaves_every_cohort_empty(window):
    assert adopt_focus(window, "compare") == set()
    assert window.cohorts() == {"isolates": set(), "compare": set(),
                                "evidence": set(), "reports": set()}
    assert "Nothing is included automatically." in window.messages[-1]
    assert window.refreshed == []


def test_a_tab_without_its_own_cohort_cannot_be_given_one(window):
    window.focus.set_focus(["a"], "Isolate table selection")
    with pytest.raises(ValueError):
        adopt_focus(window, "schemes")
    with pytest.raises(ValueError):
        send_selection(window, "settings", ["a"])


def test_sending_a_selection_sets_that_tabs_cohort_and_goes_there(window):
    sent = send_selection(window, "compare", ["b", "a"], origin="Isolate table selection")
    assert sent == {"a", "b"}
    assert window.cohort_ids == {"a", "b"}
    assert window.report_ids == set() and window.feature_ids == set()
    assert window.project.settings["comparison_cohort"] == ["a", "b"]
    assert window.visited == [window.page_index["compare"]]
    assert window.refreshed == ["compare"]
    assert "belongs to that tab only" in window.messages[-1]


def test_sending_an_empty_selection_changes_nothing(window):
    assert send_selection(window, "reports", []) == set()
    assert window.report_ids == set()
    assert window.visited == []
    assert window.project.settings == {}
    assert "Select isolates first." in window.messages[-1]


def test_only_the_comparison_cohort_is_persisted_as_a_project_setting(window):
    send_selection(window, "reports", ["a"], navigate=False)
    send_selection(window, "evidence", ["a"], navigate=False)
    assert window.project.settings == {}
    send_selection(window, "compare", ["a"], navigate=False)
    assert list(window.project.settings) == ["comparison_cohort"]


def test_the_cohort_bar_says_plainly_that_nothing_is_inherited(window, qtbot):
    bar = cohort_bar(window, "reports")
    qtbot.addWidget(bar)
    assert bar.text.text() == EMPTY_COHORT_TEXT
    assert "nothing here is inherited from another tab" in bar.text.text()
    assert bar.adopt_button.isEnabled() is False

    window.focus.set_focus(["a", "b", "c"], "Folder: Klebsiella / pneumoniae")
    assert bar.adopt_button.text() == "Use current focus (3)"
    assert bar.adopt_button.isEnabled() is True
    bar.adopt()
    assert window.report_ids == {"a", "b", "c"}
    assert "Folder: Klebsiella / pneumoniae" in bar.text.text()
    assert "3 isolates" in bar.text.text()


def test_the_cohort_bar_offers_the_tabs_own_chooser_when_the_window_has_one(window, qtbot):
    bar = CohortBar(window, "compare")
    qtbot.addWidget(bar)
    buttons = [child.text() for child in bar.findChildren(type(bar.adopt_button))]
    assert "Choose isolates…" in buttons
    evidence = CohortBar(window, "evidence")
    qtbot.addWidget(evidence)
    labels = [child.text() for child in evidence.findChildren(type(bar.adopt_button))]
    assert "Choose isolates…" not in labels, "the window has no evidence chooser here"


def test_the_focus_bar_shows_the_project_count_until_something_is_focused(window, qtbot):
    bar = FocusBar(window.focus)
    qtbot.addWidget(bar)
    bar.show()
    qtbot.waitUntil(bar.isVisible, timeout=2000)
    bar.set_project_count(12)
    assert bar.text.text() == "12 isolates in project"
    assert bar.clear_button.isVisible() is False

    window.focus.set_focus(["a", "b"], "Graph selection")
    assert "2 isolates" in bar.text.text()
    assert "Graph selection" in bar.text.text()
    assert bar.show_button.isVisible() is True
    assert bar.clear_button.isVisible() is True
    window.focus.clear()
    assert bar.text.text() == "12 isolates in project"
    assert bar.show_button.isVisible() is False


def test_focus_colour_and_legend_cannot_be_mistaken_for_a_saved_highlight():
    assert FOCUS_ACCENT != SAVED_HIGHLIGHT_ACCENT
    assert "temporary" in FOCUS_LEGEND
    assert "saved" in FOCUS_LEGEND


def test_cohort_ids_reads_each_tabs_own_attribute(window):
    window.feature_ids = {"x"}
    assert cohort_ids(window, "evidence") == {"x"}
    assert cohort_ids(window, "compare") == set()
    assert cohort_ids(window, "schemes") == set()


# ---------------------------------------------------------------------------
# The real window: focus is installed, shared and visible, and still writes
# nothing into the cohorts this revision keeps separate.
# ---------------------------------------------------------------------------


@pytest.fixture
def studio(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "workspace")
    qtbot.addWidget(widget)
    yield widget
    widget.close()


def two_isolates(studio, tmp_path):
    ids = []
    for name in ("alpha", "bravo"):
        path = tmp_path / f"{name}.fasta"
        path.write_text(">c\nACGTACGTACGT\n", encoding="utf-8")
        ids.append(path)
    studio.import_paths(ids)
    return [sample["id"] for sample in studio.project.samples()]


def test_the_window_carries_one_focus_bus_and_one_cohort_ledger(studio):
    assert isinstance(studio.focus, FocusBus)
    assert isinstance(studio.cohort_origins, CohortLedger)
    assert studio.focus.ids == set()


def test_the_top_strip_shows_the_project_count_and_then_what_is_focused(studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    assert isinstance(studio.focus_bar, FocusBar)
    # refresh_journey still writes the project count through the old attribute name.
    assert studio.scope_label is studio.focus_bar
    assert studio.focus_bar.text.text() == "2 isolates in project"
    studio.focus.set_focus(ids, "Isolate table selection")
    assert "2 isolates" in studio.focus_bar.text.text()
    assert "Isolate table selection" in studio.focus_bar.text.text()
    studio.focus.clear()
    assert studio.focus_bar.text.text() == "2 isolates in project"


def test_focusing_isolates_on_the_real_window_writes_no_cohort(studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    studio.focus.set_focus(ids, "Graph selection")
    assert studio.focus.ids == set(ids)
    assert studio.cohort_ids == set()
    assert studio.feature_ids == set()
    assert studio.report_ids == set()
    assert studio.project.get_setting("comparison_cohort", None) in (None, [])


def test_adopting_focus_into_the_report_cohort_records_where_it_came_from(studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    studio.focus.set_focus(ids, "Graph selection")
    adopted = adopt_focus(studio, "reports")
    assert adopted == set(ids)
    assert studio.report_ids == set(ids)
    assert studio.cohort_ids == set()
    assert studio.feature_ids == set()
    entry = studio.cohort_origins.entry("reports")
    assert entry["origin"] == "Graph selection"
    assert "Graph selection" in studio.cohort_origins.describe("reports", len(ids))
    assert "belongs to this tab only" in studio.cohort_origins.describe("reports", len(ids))


def test_sending_a_selection_to_compare_goes_there_and_persists_only_that_cohort(
        studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    send_selection(studio, "compare", ids, origin="Isolate table selection")
    assert studio.cohort_ids == set(ids)
    assert studio.pages.currentIndex() == studio.page_index["compare"]
    assert sorted(studio.project.get_setting("comparison_cohort", [])) == sorted(ids)
    assert studio.report_ids == set()
    assert studio.feature_ids == set()


def test_showing_the_focus_takes_you_to_the_isolates_tab_without_changing_anything(
        studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    studio.focus.set_focus(ids, "Graph selection")
    studio.navigate(studio.page_index["compare"])
    studio.show_focus_in_library()
    assert studio.pages.currentIndex() == studio.page_index["isolates"]
    assert studio.cohort_ids == set()
    assert studio.feature_ids == set()
    assert studio.report_ids == set()


def test_removing_an_isolate_drops_it_from_the_shared_focus(studio, tmp_path):
    ids = two_isolates(studio, tmp_path)
    studio.focus.set_focus(ids, "Isolate table selection")
    studio.project.remove_sample(ids[0])
    studio.refresh()
    studio.focus.prune({sample["id"] for sample in studio.current_samples})
    assert studio.focus.ids == {ids[1]}


def test_each_cohort_tab_says_what_it_is_reviewing_and_offers_only_a_click(studio, tmp_path):
    from wmlstudio.workspace_focus import COHORT_ATTRIBUTES
    for key in COHORT_ATTRIBUTES:
        bar = studio.page_headers[key].cohort_bar
        assert bar.text.text() == "Nothing chosen here yet"
        assert EMPTY_COHORT_TEXT in bar.text.toolTip()
        assert bar.adopt_button.isEnabled() is False


def test_the_adopt_button_is_the_only_thing_that_turns_focus_into_a_cohort(studio, tmp_path,
                                                                          qtbot):
    ids = two_isolates(studio, tmp_path)
    bar = studio.page_headers["reports"].cohort_bar
    studio.focus.set_focus(ids, "Isolate table selection")
    assert bar.adopt_button.text() == "Use current focus (2)"
    assert bar.adopt_button.isEnabled() is True
    # Focusing alone changed nothing.
    assert studio.report_ids == set()
    bar.adopt_button.click()
    assert studio.report_ids == set(ids)
    assert "Isolate table selection" in bar.text.text()
    assert "belongs to this tab only" in bar.text.toolTip()
    # And only that tab's cohort moved.
    assert studio.cohort_ids == set()
    assert studio.feature_ids == set()
