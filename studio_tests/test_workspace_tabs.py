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


def test_a_page_folds_into_a_samples_sub_tab_without_moving_any_index(tabs):
    """Samples is the hub; a page that also lives there must not be offered twice."""
    inner = QTabWidget()
    inner.addTab(QWidget(), "Samples")
    inner.addTab(QWidget(), "Schemes")
    tabs.register_subtabs("isolates", inner)
    assert tabs.fold_page("schemes") is True
    assert tabs.folded_pages() == {"schemes": ("isolates", "Schemes")}
    assert tabs.tabBar().isTabVisible(PAGE_KEYS.index("schemes")) is False
    # Every existing call site keeps working, by key and by its old number.
    assert tabs.page_index() == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert tabs.index_of("schemes") == PAGE_KEYS.index("schemes")
    assert tabs.show_page("schemes") is True
    assert tabs.current_key() == "isolates"
    assert inner.currentIndex() == 1
    assert tabs.show_page(PAGE_KEYS.index("schemes")) is True
    assert tabs.current_key() == "isolates"
    assert tabs.unfold_page("schemes") is True
    assert tabs.tabBar().isTabVisible(PAGE_KEYS.index("schemes")) is True
    assert tabs.show_page("schemes") is True
    assert tabs.current_key() == "schemes"


def test_folding_refuses_a_page_or_a_host_it_does_not_have(tabs):
    assert tabs.fold_page("does-not-exist") is False
    assert tabs.fold_page("compare", "does-not-exist", "Anything") is False
    assert tabs.fold_page("compare", "compare", "Itself") is False
    assert tabs.unfold_page("compare") is False
    assert tabs.folded_pages() == {}


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


# ---------------------------------------------------------------------------
# The real window: the pages it builds, the tabs it shows and the guarded
# navigation every existing call site still goes through.
# ---------------------------------------------------------------------------


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    # A modal warning box in a headless run blocks the event loop for good, so the
    # window records its errors here instead of showing them.
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    if widget.worker is not None and widget.worker.isRunning():
        widget.worker.cancel()
    settled(widget, qtbot)
    widget.close()


def page_content(window, key):
    """The page's own widget, behind the scroll area every page is wrapped in."""
    from PySide6.QtWidgets import QScrollArea
    page = window.pages.widget(window.page_index[key])
    return page.widget() if isinstance(page, QScrollArea) else page


def test_the_window_builds_seven_keyed_tabs_and_the_old_indices_still_address_them(window):
    assert isinstance(window.pages, WorkspaceTabs)
    assert window.pages.count() == 7
    assert window.page_index == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert window.pages.keys() == PAGE_KEYS
    for index, key in enumerate(PAGE_KEYS):
        assert window.pages.tabBar().tabData(index) == key
        assert TAB_LABELS[key] in window.pages.tabText(index)
        assert NAV_SYMBOLS[index] in window.pages.tabText(index)
        # The tooltip keeps the long workspace name the breadcrumb and Alt+1…7 use.
        assert window.nav_names[index] in window.pages.tabToolTip(index)
        assert PAGE_PURPOSE[key] in window.pages.tabToolTip(index)


def test_clicking_a_tab_updates_the_breadcrumb_without_navigating_twice(window):
    calls = []
    original = type(window).navigate
    type(window).navigate = lambda self, index: (calls.append(index), original(self, index))[1]
    try:
        for index, key in enumerate(PAGE_KEYS):
            calls.clear()
            window.pages.setCurrentIndex(index)
            assert window.pages.currentIndex() == index
            assert window.breadcrumb.text() == "WORKSPACE  /  " + window.nav_names[index].upper()
            # currentChanged dispatches navigate, which must not re-enter itself.
            assert len(calls) <= 2, (key, calls)
    finally:
        type(window).navigate = original


def test_an_index_outside_the_workspace_is_refused_rather_than_guessed(window):
    window.navigate(2)
    before = window.breadcrumb.text()
    window.navigate(-1)
    window.navigate(99)
    assert window.pages.currentIndex() == 2
    assert window.breadcrumb.text() == before


def test_page_six_is_a_real_settings_tab_and_no_longer_opens_a_dialog(window, monkeypatch):
    opened = []
    monkeypatch.setattr(type(window), "open_interface_settings",
                        lambda self: opened.append(True))
    window.navigate(window.page_index["settings"])
    assert window.pages.currentIndex() == window.page_index["settings"]
    assert opened == []
    assert window.settings_tabs.count() == 3
    # QTabBar reads a single "&" as a mnemonic, so the titles escape it.
    shown = [window.settings_tabs.tabText(i).replace("&&", "&") for i in range(3)]
    assert shown == ["Window, text & display", "Data & references",
                     "What this version can and cannot do"]
    # The user could not find the display controls at all, so the tab is named
    # for what it does and its tooltip repeats the words they would search for.
    assert "screen resolution" in window.settings_tabs.tabToolTip(0)


def test_the_settings_tab_keeps_the_boundaries_paragraph_word_for_word(window):
    guide = window.settings_guide.toHtml()
    assert "not a validated diagnostic device" in guide
    assert "does not predict measured susceptibility or prove transmission" in guide


def test_the_compare_page_refreshes_itself_when_its_tab_becomes_current(window):
    refreshes = []
    window.page_shown["compare"] = lambda: refreshes.append(True)
    window.navigate(window.page_index["overview"])
    window.navigate(window.page_index["compare"])
    assert refreshes == [True]


def test_every_tab_explains_itself_and_offers_one_obvious_next_step(window):
    from PySide6.QtWidgets import QLabel, QPushButton
    for key in PAGE_KEYS:
        content = page_content(window, key)
        purposes = [child for child in content.findChildren(QLabel)
                    if child.text() == PAGE_PURPOSE[key]]
        assert purposes, key
        text, method = NEXT_STEP[key]
        buttons = [child for child in content.findChildren(QPushButton)
                   if child.text() == text]
        # A next step whose method does not exist yet must not be offered at all.
        assert bool(buttons) is callable(getattr(window, method, None)), key


def test_the_orientation_strip_repeats_the_limits_a_novice_could_read_past(window):
    from PySide6.QtWidgets import QLabel
    compare = page_content(window, "compare")
    lines = [child.text() for child in compare.findChildren(QLabel)]
    assert any("not proof of transmission" in line for line in lines)
    evidence = page_content(window, "evidence")
    lines = [child.text() for child in evidence.findChildren(QLabel)]
    assert any("not measured susceptibility" in line for line in lines)


def test_the_pages_other_controllers_rearrange_are_left_exactly_as_they_were(window):
    """The evidence and scheme pages move their own layout items by position."""
    # ui_reports.build_hydra re-homed the imported-source controls into their own tab.
    titles = [window.evidence_tabs.tabText(i) for i in range(window.evidence_tabs.count())]
    assert "Advanced · imported source" in titles
    assert window.hydra_table.parent() is not None
    # ui_workbench.build_schemes inserted its own action row.
    assert window.scheme_table.columnCount() == 3
    # The evidence scope label stayed on the evidence page, not inside the strip.
    from PySide6.QtWidgets import QLabel
    evidence = page_content(window, "evidence")
    assert window.feature_scope_label in evidence.findChildren(QLabel)


def test_the_tab_bar_reaches_every_page_on_the_narrowest_supported_window(window, qtbot):
    """Seven tabs must be clickable without scroll arrows at the minimum window size."""
    window.resize(1000, 680)
    # Hiding the sidebar happens in resizeEvent; the layout that gives its width
    # to the tabs runs on the next pass, so wait for the settled geometry.
    qtbot.waitUntil(lambda: not window.sidebar.isVisible()
                    and window.pages.width() >= 940, timeout=5000)
    assert window.pages.tabBar().sizeHint().width() <= window.pages.width()
    window.resize(1380, 940)
    qtbot.waitUntil(lambda: window.sidebar.isVisible(), timeout=5000)
    assert window.pages.tabBar().sizeHint().width() <= window.pages.width()


# ---------------------------------------------------------------------------
# The menu bar this window owns: a dedicated Update menu for everything that
# can be installed, and the display entries a user could not find before.
# ---------------------------------------------------------------------------


def test_the_menu_bar_carries_a_dedicated_update_menu_before_help(window):
    titles = [entry.text().replace("&", "") for entry in window.menuBar().actions()
              if entry.menu() is not None]
    assert "Update" in titles
    assert titles.index("Update") == titles.index("Help") - 1
    assert titles.index("Update") > titles.index("View")


def test_the_update_menu_names_every_installable_thing_once(window):
    from wmlstudio.update_center import MENU_ENTRIES
    entries = [action.text() for action in window.update_menu.actions()
               if not action.isSeparator()]
    assert entries[0].startswith("What is installed")
    for _key, title in MENU_ENTRIES:
        assert title in entries, title
    assert len(entries) == len(set(entries)), "nothing is offered twice"
    # A menu entry that cannot be used is worse than no entry at all.
    for action in window.update_menu.actions():
        if not action.isSeparator():
            assert action.text().strip() and action.isEnabled(), action.text()
    # The same commands are reachable from the Ctrl+K search, which is where a
    # user who cannot find a menu looks next.
    commands = [title for title, _ in window.command_actions]
    assert any(title.startswith("What is installed") for title in commands)


def settled(window, qtbot):
    """Wait until the background probe has finished and its result has been handled."""
    # worker_role is cleared by analysis_finished, so an empty role means the
    # queued finish has actually run — closing the project before it does raises.
    qtbot.waitUntil(lambda: window.worker is None or
                    (not window.worker.isRunning() and not window.worker_role), timeout=60000)


def update_centre(window, qtbot):
    """Open the Update Centre and wait for its background probe to land."""
    from wmlstudio.update_center import open_update_center
    centre = open_update_center(window)
    qtbot.addWidget(centre)
    qtbot.waitUntil(lambda: centre.report is not None, timeout=60000)
    settled(window, qtbot)
    assert window.test_errors == []
    return centre


def test_the_update_centre_lists_installed_state_and_never_checks_a_server(window, qtbot):
    from wmlstudio.update_center import UpdateCenter, open_update_center
    centre = update_centre(window, qtbot)
    assert isinstance(centre, UpdateCenter)
    assert open_update_center(window) is centre, "one page per window, not one per click"
    settled(window, qtbot)
    keys = [item["key"] for item in centre.report["items"]]
    from wmlstudio.provisioning import ORDER
    assert keys == list(ORDER)
    assert centre.table.rowCount() == len(ORDER)
    assert "until you press a button" in centre.status.text() or centre.report["summary"]
    row = keys.index("species_panel")
    # This workspace has no downloaded panel, and the page says so rather than
    # implying an update was checked for online.
    assert centre.table.item(row, 1).text() == "Not installed"
    assert centre.table.cellWidget(row, 4).text() == "Install…"
    assert centre.table.cellWidget(row, 4).isEnabled() is True
    centre.close()


def test_an_update_row_uses_the_installer_the_application_already_has(window, qtbot,
                                                                     monkeypatch):
    centre = update_centre(window, qtbot)
    called = []
    monkeypatch.setattr(type(window), "install_species_panel",
                        lambda self: called.append("species"))
    assert centre.start("species_panel") is True
    assert called == ["species"]
    assert "Check what is installed" in centre.status.text()
    centre.close()


def test_an_item_with_no_installer_here_explains_how_it_is_installed(window, qtbot,
                                                                     monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    centre = update_centre(window, qtbot)
    shown = []
    monkeypatch.setattr(QMessageBox, "information",
                        lambda parent, title, text, *args: shown.append((title, text)))
    # BLAST+ ships beside the application; nothing downloads it on a user's behalf.
    assert centre.start("blast_tools") is False
    assert shown and "BLAST" in shown[0][0]
    assert "stage_bio_tools" in shown[0][1] or "Tools/blast" in shown[0][1]
    centre.close()
