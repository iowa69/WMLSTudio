"""Key-addressed workspace tabs: stable keys, working integer indices, sub-tab registry."""

import pytest
from PySide6.QtWidgets import QTabWidget, QWidget

from wmlstudio.ui_tabs import (
    NAV_SYMBOLS,
    NEXT_STEP,
    PAGE_KEYS,
    PAGE_PURPOSE,
    TAB_LABELS,
    WorkspaceTabs,
    page_key,
    tab_title,
)


@pytest.fixture
def tabs(qtbot):
    widget = WorkspaceTabs()
    qtbot.addWidget(widget)
    widget.resize(1032, 640)
    for key in PAGE_KEYS:
        page = QWidget()
        page.setObjectName(f"page-{key}")
        widget.add_page(page)
    widget.show()
    return widget


def test_pages_are_added_in_build_order_and_carry_their_stable_key(tabs):
    assert tabs.count() == len(PAGE_KEYS) == 7
    assert tabs.keys() == PAGE_KEYS
    assert [tabs.tabBar().tabData(index) for index in range(tabs.count())] == list(PAGE_KEYS)
    assert tabs.page_index() == {key: index for index, key in enumerate(PAGE_KEYS)}
    for index, key in enumerate(PAGE_KEYS):
        assert tabs.widget(index).objectName() == f"page-{key}"
        assert tabs.key_at(index) == key


def test_every_tab_shows_its_glyph_label_and_a_purpose_tooltip(tabs):
    for index, key in enumerate(PAGE_KEYS):
        assert TAB_LABELS[key] in tabs.tabText(index)
        assert NAV_SYMBOLS[index] in tabs.tabText(index)
        assert tabs.tabToolTip(index) == PAGE_PURPOSE[key]
        assert tabs.tabToolTip(index).strip()


def test_the_old_integer_indices_still_address_the_same_pages(tabs):
    """The window navigates by number in seventeen places; the swap must not move a page."""
    assert tabs.show_page(2) is True
    assert tabs.currentIndex() == 2
    assert tabs.current_key() == "compare"
    assert tabs.show_page("reports") is True
    assert tabs.currentIndex() == PAGE_KEYS.index("reports")
    assert tabs.index_of(5) == 5
    assert tabs.index_of("settings") == 6
    assert tabs.page_for(1) is tabs.widget(1)
    assert tabs.page_for("compare") is tabs.widget(2)


def test_unknown_keys_and_out_of_range_indices_are_refused_instead_of_guessed(tabs):
    tabs.show_page("overview")
    assert tabs.show_page("does-not-exist") is False
    assert tabs.show_page(99) is False
    assert tabs.show_page(-1) is False
    assert tabs.index_of("does-not-exist") == -1
    assert tabs.page_for(42) is None
    assert tabs.currentIndex() == 0
    assert page_key(42) == ""


def test_stacked_widget_call_sites_keep_working_through_add_widget(qtbot):
    widget = WorkspaceTabs()
    qtbot.addWidget(widget)
    first, second = QWidget(), QWidget()
    assert widget.addWidget(first) == 0
    assert widget.addWidget(second) == 1
    assert widget.widget(0) is first
    assert widget.keys() == ("overview", "isolates")
    assert widget.tabText(0) == tab_title("overview")


def test_a_page_show_hook_runs_when_its_tab_becomes_current_without_recursing(tabs):
    calls = []
    tabs.set_show_hook("compare", lambda: calls.append(tabs.current_key()))
    shown = []
    tabs.pageShown.connect(shown.append)
    tabs.show_page("compare")
    assert calls == ["compare"]
    assert shown == ["compare"]
    tabs.show_page("compare")
    assert calls == ["compare"], "re-showing the same page must not run the hook again"
    tabs.show_page("reports")
    assert calls == ["compare"]
    assert shown == ["compare", "reports"]


def test_a_show_hook_that_navigates_again_cannot_loop(tabs):
    """A hook may redirect once; two hooks pointing at each other must still terminate."""
    seen = []

    def hook():
        seen.append(tabs.current_key())
        tabs.show_page("overview")

    tabs.set_show_hook("settings", hook)
    tabs.show_page("settings")
    assert seen == ["settings"]
    assert tabs.current_key() == "overview"

    bounces = []
    tabs.set_show_hook("overview", lambda: (bounces.append("overview"),
                                            tabs.show_page("settings")))
    tabs.set_show_hook("settings", lambda: (bounces.append("settings"),
                                            tabs.show_page("overview")))
    tabs.show_page("compare")
    tabs.show_page("settings")
    assert len(bounces) <= 3


def test_removing_a_show_hook_leaves_navigation_working(tabs):
    calls = []
    tabs.set_show_hook("compare", lambda: calls.append(1))
    tabs.set_show_hook("compare", None)
    tabs.show_page("compare")
    assert calls == []
    assert tabs.show_hook("compare") is None


def test_sub_tabs_are_registered_per_page_and_extendable_by_later_work(tabs):
    inner = QTabWidget()
    inner.addTab(QWidget(), "Current tree")
    tabs.register_subtabs("compare", inner)
    assert tabs.subtabs("compare") is inner
    assert tabs.subtab_keys() == ("compare",)
    index = tabs.add_subtab("compare", QWidget(), "Original tree")
    assert index == 1
    assert inner.tabText(1) == "Original tree"
    assert tabs.add_subtab("compare", QWidget(), "Side by side", index=1) == 1
    assert [inner.tabText(i) for i in range(inner.count())] == [
        "Current tree", "Side by side", "Original tree"]
    assert tabs.show_subtab("compare", "Original tree") is True
    assert inner.currentIndex() == 2


def test_a_page_without_sub_tabs_refuses_a_sub_tab_instead_of_crashing(tabs):
    assert tabs.subtabs("schemes") is None
    assert tabs.add_subtab("schemes", QWidget(), "Grouped view") == -1
    assert tabs.show_subtab("schemes", "Grouped view") is False


def test_seven_tabs_fit_the_narrowest_supported_window_without_scroll_buttons(tabs):
    """At 1080 px the tab bar is the navigation; it must not need arrows to reach a page."""
    assert tabs.tabBar().sizeHint().width() <= 1032


def test_purpose_lines_state_their_limits_where_a_claim_could_be_read_in(tabs):
    assert set(PAGE_PURPOSE) == set(PAGE_KEYS)
    assert set(NEXT_STEP) == set(PAGE_KEYS)
    assert "not proof of transmission" in PAGE_PURPOSE["compare"]
    assert "not measured susceptibility" in PAGE_PURPOSE["evidence"]
    for key, purpose in PAGE_PURPOSE.items():
        assert purpose.endswith("."), key
    for key, (text, method) in NEXT_STEP.items():
        assert text and method.isidentifier(), key
