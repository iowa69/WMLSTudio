"""The workspace laid out as the pipeline: build-order numbers, pipeline tabs, Clear."""

import pytest
from PySide6.QtWidgets import QLabel, QMessageBox, QPushButton, QTabWidget, QWidget

from wmlstudio.ui_tabs import (
    CLEAR_ACTIONS,
    NEXT_STEP,
    PAGE_KEYS,
    PAGE_PURPOSE,
    PIPELINE,
    PLANNED,
    STATION_PAGES,
    TAB_LABELS,
    WorkspaceTabs,
    clear_promise,
    page_key,
    page_name,
    tab_title,
)


@pytest.fixture
def tabs(qtbot):
    # Tab widths are measured from the font in force, and other test files leave
    # the application font scaled, so a bar that fits at one point size does not
    # at another. Pin it on this widget: setting it on the application would
    # re-polish every open window, which is slow and was a real defect once.
    from PySide6.QtGui import QFont
    widget = WorkspaceTabs()
    widget.setFont(QFont(widget.font().family(), 9))
    qtbot.addWidget(widget)
    widget.resize(1032, 640)
    for key in PAGE_KEYS:
        page = QWidget()
        page.setObjectName(f"page-{key}")
        widget.add_page(page)
    widget.show()
    return widget


def test_pages_keep_their_build_order_number_and_are_shown_in_pipeline_order(tabs):
    """Two orders, on purpose: numbers address the build, the bar shows the workflow."""
    assert tabs.count() == len(PAGE_KEYS) == len(PIPELINE)
    assert tabs.keys() == PAGE_KEYS
    assert tabs.tab_order() == PIPELINE
    assert tabs.page_index() == {key: index for index, key in enumerate(PAGE_KEYS)}
    for slot, key in enumerate(PAGE_KEYS):
        assert tabs.widget(slot).objectName() == f"page-{key}"
        assert tabs.key_at(slot) == key
        assert tabs.slot_of(key) == slot
        assert tabs.position_of(key) == PIPELINE.index(key)
        assert tabs.key_at_position(PIPELINE.index(key)) == key


def test_the_bar_reads_as_the_users_own_workflow(tabs):
    """Samples, read QC, assembly, typing, its tree, cgMLST, its tree, SNP, HYDRA."""
    shown = [tabs.tabText(position) for position in range(tabs.count())]
    assert shown == [
        "Overview", "Samples", "Read QC", "Assembly", "MLST", "MLST tree", "cgMLST",
        "cgMLST tree", "SNP tree", "HYDRA", "Report", "Update", "Schemes", "Settings"]
    # A seven-locus tree and a cgMLST tree are different quantities, so they are
    # different tabs and neither label can be mistaken for the other.
    assert shown.index("MLST tree") < shown.index("cgMLST") < shown.index("cgMLST tree")
    # Raw reads come first: trim, then assemble, then type what was assembled.
    assert shown.index("Read QC") == shown.index("Samples") + 1
    assert shown.index("Assembly") == shown.index("Read QC") + 1
    assert shown.index("MLST") == shown.index("Assembly") + 1
    # The one slot a later round still fills sits where that work belongs, not at
    # the end: SNP distances are a separate line of evidence from allele typing.
    assert shown.index("SNP tree") == shown.index("cgMLST tree") + 1
    # The reported bug: the tab labelled Update was the scheme library, so pressing
    # it rescanned this computer's scheme folders instead of looking online. Update
    # is its own page now, in the place the bar already put it, and the library
    # keeps a tab named for what it actually is.
    assert shown.index("Update") == shown.index("Report") + 1
    assert tabs.key_at_position(shown.index("Update")) == "update"
    assert tabs.key_at_position(shown.index("Schemes")) == "schemes"
    assert "Update" not in TAB_LABELS["schemes"]


def test_every_tab_carries_a_purpose_tooltip_and_a_readable_label(tabs):
    for position, key in enumerate(PIPELINE):
        assert tabs.tabText(position) == TAB_LABELS[key] == tab_title(key)
        assert PAGE_PURPOSE[key] in tabs.tabToolTip(position)
        assert tabs.tabToolTip(position).strip()


def test_the_old_integer_call_sites_still_address_the_same_pages(tabs):
    """`navigate(2)` and `pages.widget(3)` predate the pipeline bar and must not move."""
    assert tabs.show_page(2) is True
    assert tabs.current_key() == "compare"
    assert tabs.currentIndex() == 2, "an integer is the page number, not the tab position"
    assert tabs.tabBar().currentIndex() == PIPELINE.index("compare")
    assert tabs.widget(3) is tabs.page_for("schemes")
    assert tabs.show_page("reports") is True
    assert tabs.currentIndex() == PAGE_KEYS.index("reports")
    assert tabs.index_of(5) == 5
    assert tabs.index_of("settings") == PAGE_KEYS.index("settings")
    assert tabs.page_for(1) is tabs.widget(1)
    assert tabs.page_for("compare") is tabs.widget(2)
    tabs.setCurrentIndex(4)
    assert tabs.current_key() == "evidence"


def test_unknown_keys_and_out_of_range_indices_are_refused_instead_of_guessed(tabs):
    tabs.show_page("overview")
    assert tabs.show_page("does-not-exist") is False
    assert tabs.show_page(99) is False
    assert tabs.show_page(-1) is False
    assert tabs.index_of("does-not-exist") == -1
    assert tabs.page_for(42) is None
    assert tabs.current_key() == "overview"
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


def test_a_page_folds_into_another_pages_sub_tab_without_moving_any_number(tabs):
    """A page that grows a second home must have one place to look, not two."""
    inner = QTabWidget()
    inner.addTab(QWidget(), "Samples")
    inner.addTab(QWidget(), "Schemes")
    tabs.register_subtabs("isolates", inner)
    assert tabs.fold_page("schemes", "isolates", "Schemes") is True
    assert tabs.folded_pages() == {"schemes": ("isolates", "Schemes")}
    assert tabs.tabBar().isTabVisible(PIPELINE.index("schemes")) is False
    # Every existing call site keeps working, by key and by its old number.
    assert tabs.page_index() == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert tabs.index_of("schemes") == PAGE_KEYS.index("schemes")
    assert tabs.show_page("schemes") is True
    assert tabs.current_key() == "isolates"
    assert inner.currentIndex() == 1
    assert tabs.show_page(PAGE_KEYS.index("schemes")) is True
    assert tabs.current_key() == "isolates"
    assert tabs.unfold_page("schemes") is True
    assert tabs.tabBar().isTabVisible(PIPELINE.index("schemes")) is True
    assert tabs.show_page("schemes") is True
    assert tabs.current_key() == "schemes"


def test_folding_refuses_a_page_or_a_host_it_does_not_have(tabs):
    assert tabs.fold_page("does-not-exist") is False
    assert tabs.fold_page("compare", "does-not-exist", "Anything") is False
    assert tabs.fold_page("compare", "compare", "Itself") is False
    assert tabs.unfold_page("compare") is False
    assert tabs.folded_pages() == {}


def test_the_bar_gives_up_padding_before_it_gives_up_words(tabs):
    """Fourteen tabs at the narrowest supported window: tighter, never chopped."""
    from PySide6.QtCore import Qt

    from wmlstudio.ui_tabs import COMFORTABLE_PADDING, PADDING_STEPS
    assert tabs.tabBar().elideMode() == Qt.TextElideMode.ElideNone
    # Widths are derived from the bar itself, not written in. Tab text is measured
    # in the font in force, and other test files leave the application font
    # scaled, so a pixel figure that fits in one run does not in the next.
    widest = PADDING_STEPS[0]
    tabs._set_padding(widest)
    roomy = tabs.tabBar().sizeHint().width() + 80
    tabs._set_padding(COMFORTABLE_PADDING)
    comfortable = tabs.tabBar().sizeHint().width()

    # A window with room for the widest padding keeps it.
    tabs.resize(roomy, 640)
    assert tabs._fit_tab_bar(roomy) == widest, "a roomy window keeps the roomiest padding"

    # A window too narrow for that gives up padding, step by step, and still fits.
    squeezed = comfortable - 40
    tabs.resize(squeezed, 640)
    chosen = tabs._fit_tab_bar(squeezed)
    assert chosen < widest and chosen in PADDING_STEPS
    assert tabs.tabBar().sizeHint().width() <= squeezed
    # Whatever it gave up, it never gave up a word.
    for position in range(tabs.count()):
        assert "…" not in tabs.tabText(position)
    # The tightest steps exist for a window that has nothing else left to give. The
    # sidebar decision is not allowed to spend them, so the width it is measured
    # against is the comfortable one, never the last step the bar could survive on.
    assert COMFORTABLE_PADDING in PADDING_STEPS and PADDING_STEPS[-1] < COMFORTABLE_PADDING
    tabs._set_padding(COMFORTABLE_PADDING)
    assert tabs.minimum_bar_width() == tabs.tabBar().sizeHint().width()


def test_purpose_lines_state_their_limits_where_a_claim_could_be_read_in(tabs):
    assert set(PAGE_PURPOSE) == set(PAGE_KEYS)
    assert set(NEXT_STEP) == set(PAGE_KEYS)
    assert set(CLEAR_ACTIONS) == set(PAGE_KEYS)
    assert "not proof of transmission" in PAGE_PURPOSE["compare"]
    assert "not measured susceptibility" in PAGE_PURPOSE["evidence"]
    assert "never share a scale" in PAGE_PURPOSE["cgmlst_tree"]
    # Update says that checking reaches the internet, and the scheme library says
    # that it does not: the two words a user could not tell apart before.
    assert "internet" in PAGE_PURPOSE["update"]
    assert "Update" in PAGE_PURPOSE["schemes"] and "installed on this computer" \
        in PAGE_PURPOSE["schemes"]
    assert NEXT_STEP["update"] == ("Check online for updates…", "check_for_updates")
    for key in PLANNED:
        assert "Planned for a later round" in PAGE_PURPOSE[key], key
    for key, purpose in PAGE_PURPOSE.items():
        assert purpose.endswith("."), key
    for key, (text, method) in NEXT_STEP.items():
        assert text and method.isidentifier(), key
    for key in PAGE_KEYS:
        promise = clear_promise(key)
        assert promise["clears"] and promise["keeps"], key


def test_a_station_page_says_what_it_does_not_do(tabs):
    """A tab a later round fills must not read as finished work."""
    for key, definition in STATION_PAGES.items():
        assert definition["title"] and definition["subtitle"].endswith(".")
        assert definition["body"], key
    for key in PLANNED:
        # The page names the round that builds it rather than implying it is ready.
        assert "a later round" in " ".join(STATION_PAGES[key]["body"]).casefold(), key
        assert tabs.mark_planned(key) is True
        assert "planned" in tabs.tabToolTip(PIPELINE.index(key))
    assert tabs.mark_planned("does-not-exist") is False


# ---------------------------------------------------------------------------
# The real window: the pages it builds, the tabs it shows, the Clear on every
# one of them, and the guarded navigation every existing call site goes through.
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


def typed_sample(window, tmp_path, name="isolate.fasta"):
    """One imported sample carrying a stored result, without running an analysis."""
    path = tmp_path / name
    path.write_text(">contig\nACGTACGTACGT\n", encoding="ascii")
    window.import_paths([path])
    sample_id = window.project.samples()[-1]["id"]
    window.project.set_result(sample_id, {
        "sample_name": path.stem, "status": "complete", "scheme": "practice_7",
        "scheme_digest": "abc", "st": "7", "alleles": {f"locus{i}": "1" for i in range(7)},
    })
    window.refresh()
    return sample_id


def test_the_window_builds_every_pipeline_page_and_the_old_numbers_still_address_them(window):
    assert isinstance(window.pages, WorkspaceTabs)
    assert window.pages.count() == len(PAGE_KEYS)
    assert window.page_index == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert window.pages.keys() == PAGE_KEYS
    assert window.pages.tab_order() == PIPELINE
    for position, key in enumerate(PIPELINE):
        assert window.pages.tabBar().tabData(position) == key
        assert window.pages.tabText(position) == TAB_LABELS[key]
        # The tooltip keeps the long workspace name the breadcrumb and Alt+N use.
        assert page_name(key) in window.pages.tabToolTip(position)
        assert PAGE_PURPOSE[key] in window.pages.tabToolTip(position)
    assert window.nav_names == [page_name(key) for key in PAGE_KEYS]


def test_clicking_a_tab_updates_the_breadcrumb_without_navigating_twice(window):
    calls = []
    original = type(window).navigate
    type(window).navigate = lambda self, target: (calls.append(target), original(self, target))[1]
    try:
        for position, key in enumerate(PIPELINE):
            calls.clear()
            window.pages.tabBar().setCurrentIndex(position)
            assert window.pages.current_key() == key
            assert window.breadcrumb.text() == "WORKSPACE  /  " + page_name(key).upper()
            # The bar's own signal dispatches navigate, which must not re-enter itself.
            assert len(calls) <= 2, (key, calls)
    finally:
        type(window).navigate = original


def test_an_index_outside_the_workspace_is_refused_rather_than_guessed(window):
    window.navigate(2)
    before = window.breadcrumb.text()
    window.navigate(-1)
    window.navigate(99)
    window.navigate("no-such-page")
    assert window.pages.current_key() == "compare"
    assert window.breadcrumb.text() == before


def test_page_six_is_a_real_settings_tab_and_no_longer_opens_a_dialog(window, monkeypatch):
    opened = []
    monkeypatch.setattr(type(window), "open_interface_settings",
                        lambda self: opened.append(True))
    window.navigate(window.page_index["settings"])
    assert window.pages.current_key() == "settings"
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


def test_every_tab_offers_a_clear_that_names_what_it_keeps(window):
    assert set(window.clear_buttons) == set(PAGE_KEYS)
    for key in PAGE_KEYS:
        content = page_content(window, key)
        # The tab's own Clear, named for the tab, on the same strip on every page.
        buttons = [child for child in content.findChildren(QPushButton)
                   if child.accessibleName() == f"Clear the {page_name(key)} tab"]
        assert len(buttons) == 1, key
        assert buttons[0] is window.clear_buttons[key]
        promise = clear_promise(key)
        assert promise["clears"] in buttons[0].toolTip()
        assert promise["keeps"] in buttons[0].toolTip()


def test_clearing_a_tab_restarts_the_view_and_deletes_no_evidence(window, tmp_path,
                                                                  monkeypatch):
    """The whole point of Clear: start again with new, past or mixed samples."""
    sample_id = typed_sample(window, tmp_path)
    window.cohort_ids = {sample_id}
    window.search.setText("isolate")
    asked = []
    monkeypatch.setattr(QMessageBox, "question",
                        lambda parent, title, text, *args: (asked.append((title, text)),
                                                            QMessageBox.StandardButton.Yes)[1])
    assert window.clear_page("compare") is True
    title, text = asked[-1]
    assert "Clear MLST tree?" == title
    assert clear_promise("compare")["clears"] in text
    assert clear_promise("compare")["keeps"] in text
    assert "Nothing you have imported or analysed is deleted." in text
    assert window.cohort_ids == set()
    assert window.distance_rows == []
    assert window.tree.nodes == {}
    assert window.clear_page("isolates") is True
    assert window.search.text() == ""
    # The evidence itself is untouched: the sample, its input and its result stay.
    stored = window.project.get_sample(sample_id)
    assert stored["result"]["st"] == "7"
    assert len(window.project.samples()) == 1
    assert window.test_errors == []


def test_clearing_can_be_refused_and_then_changes_nothing(window, tmp_path, monkeypatch):
    sample_id = typed_sample(window, tmp_path)
    window.cohort_ids = {sample_id}
    monkeypatch.setattr(QMessageBox, "question",
                        lambda *args: QMessageBox.StandardButton.No)
    assert window.clear_page("compare") is False
    assert window.cohort_ids == {sample_id}


def test_the_orientation_strip_repeats_the_limits_a_novice_could_read_past(window):
    compare = page_content(window, "compare")
    lines = [child.text() for child in compare.findChildren(QLabel)]
    assert any("not proof of transmission" in line for line in lines)
    evidence = page_content(window, "evidence")
    lines = [child.text() for child in evidence.findChildren(QLabel)]
    assert any("not measured susceptibility" in line for line in lines)
    cgmlst_tree = page_content(window, "cgmlst_tree")
    lines = [child.text() for child in cgmlst_tree.findChildren(QLabel)]
    assert any("never share a scale" in line for line in lines)


def test_a_planned_station_says_so_instead_of_pretending_to_work(window):
    for key in PLANNED:
        content = page_content(window, key)
        lines = [child.text() for child in content.findChildren(QLabel)]
        assert any("Nothing runs on this tab yet" in line for line in lines), key
        assert "planned" in window.pages.tabToolTip(window.pages.position_of(key))
        assert window.station_status[key].text().startswith("Nothing runs on this tab yet")


def test_a_typing_station_counts_what_this_project_actually_has(window, tmp_path):
    assert window.station_status["mlst"].text().startswith("0 of 0 samples")
    typed_sample(window, tmp_path)
    assert window.station_status["mlst"].text().startswith("1 of 1 samples have a stored "
                                                           "seven-locus result")
    # A cgMLST count is never inferred from a seven-locus result.
    assert window.station_status["cgmlst"].text().startswith("0 of 1 samples")
    assert window.station_status["cgmlst_tree"].text().startswith("0 samples carry a cgMLST")


def test_the_two_pipeline_pages_are_adopted_without_moving_a_single_number(window):
    """Read QC and Assembly stopped being placeholders; their numbers did not move."""
    for key in ("reads", "assembly"):
        assert key not in PLANNED
        assert window.stations[key]["adopted"] is not None, key
        assert window.stations[key]["placeholder"].isVisibleTo(window) is False
        assert window.station_status[key].isVisibleTo(window) is False
        assert "planned" not in window.pages.tabToolTip(window.pages.position_of(key))
    assert window.page_index == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert window.pages.tab_order() == PIPELINE


def test_a_station_hands_its_page_over_without_moving_a_single_number(window):
    """The seam a later round drops its finished page into."""
    page = QWidget()
    page.setObjectName("built-elsewhere")
    assert window.adopt_station("cgmlst_tree", page, title="cgMLST tree") is True
    assert page.parent() is not None
    assert window.stations["cgmlst_tree"]["placeholder"].isVisibleTo(window) is False
    assert window.page_index == {key: index for index, key in enumerate(PAGE_KEYS)}
    assert window.pages.tab_order() == PIPELINE
    assert window.adopt_station("no-such-station", QWidget()) is False


def test_a_station_only_offers_an_action_this_build_can_perform(window):
    """The cgMLST tab links to the table of calls only where this build has one."""
    labels = [child.text() for child in page_content(window, "cgmlst").findChildren(QPushButton)]
    assert ("Open the table of calls →" in labels) is (getattr(window, "cgmlst_calls", None)
                                                       is not None)
    tree_labels = [child.text() for child in
                   page_content(window, "cgmlst_tree").findChildren(QPushButton)]
    assert ("Draw it on the tree page →" in tree_labels) is callable(
        getattr(window, "show_cgmlst_tree", None))


def test_the_cgmlst_tree_station_asks_the_tree_page_for_a_cgmlst_graph(window, monkeypatch):
    """Until it draws its own, the station takes you to the page and says so."""
    drawn = []
    monkeypatch.setattr(type(window), "show_cgmlst_tree",
                        lambda self: (drawn.append(self.pages.current_key()), "cgmlst")[1])
    assert window.goto_cgmlst_tree_view() == "cgmlst"
    assert drawn == ["compare"], "the page is asked only once it is the page in front"


def test_an_adopted_station_takes_over_its_navigation_and_its_own_clear(window, monkeypatch):
    class Page(QWidget):
        cleared = 0

        def clear(self):
            type(self).cleared += 1

    page = Page()
    assert window.adopt_station("cgmlst", page) is True
    window.navigate("overview")
    window.goto_cgmlst_calls()
    assert window.pages.current_key() == "cgmlst", "the calls live on the cgMLST tab now"
    monkeypatch.setattr(QMessageBox, "question",
                        lambda *args: QMessageBox.StandardButton.Yes)
    assert window.clear_page("cgmlst") is True
    assert Page.cleared == 1


def test_a_station_action_takes_you_where_that_work_happens_today(window):
    window.navigate("assembly")
    assert window.goto_assembly_step() is True
    assert window.pages.current_key() == "isolates"
    assert window.sample_tabs.tabText(window.sample_tabs.currentIndex()) == "Assembly"
    window.navigate("snp")
    window.goto_mlst_tree()
    assert window.pages.current_key() == "compare"


def test_the_pages_other_controllers_rearrange_are_left_exactly_as_they_were(window):
    """The evidence and scheme pages move their own layout items by position."""
    # ui_reports.build_hydra re-homed the imported-source controls into their own tab.
    titles = [window.evidence_tabs.tabText(i) for i in range(window.evidence_tabs.count())]
    assert "Advanced · imported source" in titles
    assert window.hydra_table.parent() is not None
    # ui_workbench.build_schemes inserted its own action row.
    assert window.scheme_table.columnCount() == 4
    # The evidence scope label stayed on the evidence page, not inside the strip.
    evidence = page_content(window, "evidence")
    assert window.feature_scope_label in evidence.findChildren(QLabel)


def test_the_scheme_library_names_each_scheme_and_says_which_kind_it_is(window):
    """The reported bug: a downloaded cgMLST scheme could be installed and not found."""
    assert [window.scheme_table.horizontalHeaderItem(column).text()
            for column in range(4)] == ["SCHEME", "TYPE", "SOURCE", "LOCATION"]
    kinds = {window.scheme_table.item(row, 1).text()
             for row in range(window.scheme_table.rowCount())
             if window.scheme_table.item(row, 1) is not None}
    assert kinds <= {"MLST", "cgMLST", "Not classified"}
    for row in range(window.scheme_table.rowCount()):
        title, kind = (window.scheme_table.item(row, column).text() for column in (0, 1))
        # The name is the scheme's own, never the install folder read out loud:
        # "cgmlst org kpneumoniae abcdef0123456789" is what the user could not find.
        assert not title.startswith("cgmlst org "), title
        if kind in {"MLST", "cgMLST"}:
            assert kind in title or " · " in title, title


def test_the_tab_bar_reaches_every_page_on_the_narrowest_supported_window(window, qtbot):
    """Thirteen tabs must stay readable, and every page reachable, at 1000x680."""
    window.resize(1000, 680)
    # Hiding the sidebar happens in resizeEvent; the layout that gives its width
    # to the tabs runs on the next pass, so wait for the settled geometry.
    qtbot.waitUntil(lambda: not window.sidebar.isVisible()
                    and window.pages.width() >= 940, timeout=5000)
    assert window.pages.tabBar().sizeHint().width() <= window.pages.width()
    for position, key in enumerate(PIPELINE):
        assert "…" not in window.pages.tabText(position)
        window.pages.tabBar().setCurrentIndex(position)
        assert window.pages.current_key() == key
    window.resize(1380, 940)
    qtbot.waitUntil(lambda: window.sidebar.isVisible(), timeout=5000)
    assert window.pages.tabBar().sizeHint().width() <= window.pages.width()


def test_the_sidebar_stands_down_before_a_tab_goes_behind_a_scroll_arrow(window, qtbot):
    """Navigation first: the decorative column is what gives up its width."""
    window.resize(1380, 940)
    qtbot.waitUntil(lambda: window.sidebar.isVisible(), timeout=5000)
    needed = window.pages.minimum_bar_width()
    assert needed > 0
    # Just too narrow to carry both: the sidebar goes, the tabs stay whole.
    window.resize(needed + 48 + window.sidebar.width() - 30, 940)
    qtbot.waitUntil(lambda: not window.sidebar.isVisible()
                    and window.pages.width() >= needed, timeout=5000)
    assert window.pages.tabBar().sizeHint().width() <= window.pages.width()


def test_a_larger_interface_scale_re_measures_the_bar_instead_of_running_off(window, qtbot):
    """Text scaling has to keep working: the bar is measured again, never clipped."""
    window.resize(1380, 940)
    qtbot.waitUntil(lambda: window.pages.width() > 1000, timeout=5000)
    small = window.pages.minimum_bar_width()
    window.set_ui_scale(130)
    try:
        qtbot.waitUntil(lambda: window.pages.minimum_bar_width() > small, timeout=5000)
        qtbot.waitUntil(lambda: window.pages.tabBar().sizeHint().width()
                        <= window.pages.width(), timeout=5000)
        # Bigger letters buy their room from the padding and the sidebar, never
        # from the words: a tab that reads "cgMLST t…" has stopped being navigation.
        assert window.pages.tabBar().elideMode().name == "ElideNone"
        for position in range(window.pages.count()):
            assert "…" not in window.pages.tabText(position)
    finally:
        window.set_ui_scale(100)


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
    """Show the Update tab and wait for its own read of this computer to land."""
    window.navigate("update")
    centre = window.update_center
    qtbot.waitUntil(lambda: centre.report is not None, timeout=60000)
    qtbot.waitUntil(lambda: centre.worker is None or not centre.worker.isRunning(),
                    timeout=60000)
    settled(window, qtbot)
    assert window.test_errors == []
    return centre


def answer(centre, monkeypatch, *, yes=True):
    """Answer this page's own confirmations, since a modal box never returns headless."""
    asked, told = [], []
    monkeypatch.setattr(type(centre), "ask",
                        lambda self, title, text: (asked.append((title, text)), yes)[1])
    monkeypatch.setattr(type(centre), "tell",
                        lambda self, title, text: told.append((title, text)))
    return asked, told


def fake_plan(monkeypatch, **release):
    """A check that reaches a provider, or does not, without touching the network."""
    from wmlstudio import provisioning
    answer_release = {"installed": "2026-08-07.1", "latest": "2026-08-07.1",
                      "newer_available": False, "message": "", "error": "", **release}
    monkeypatch.setattr(provisioning, "hydra_update_available",
                        lambda *args, **options: dict(answer_release))
    return answer_release


def test_the_update_tab_is_the_update_centre_and_the_library_keeps_its_own_tab(window, qtbot):
    """The reported bug: pressing Update rescanned this computer instead of looking online."""
    from wmlstudio.update_center import UpdateCenter, open_update_center
    centre = update_centre(window, qtbot)
    assert isinstance(centre, UpdateCenter)
    assert window.pages.current_key() == "update"
    assert centre.isVisibleTo(window)
    # The menu, the command search and the tab are one page, not three copies.
    assert open_update_center(window) is centre
    assert window.pages.current_key() == "update"
    # The library is still there, under its own name, and is not this page.
    window.navigate("schemes")
    assert window.pages.tabText(window.pages.position_of("schemes")) == "Schemes"
    assert centre.isVisibleTo(window) is False
    assert window.scheme_table.rowCount() >= 0


def test_the_update_page_lists_every_installable_thing_with_its_size(window, qtbot):
    from wmlstudio.provisioning import ORDER
    centre = update_centre(window, qtbot)
    keys = [item["key"] for item in centre.report["items"]]
    assert keys == list(ORDER)
    listed = [row["key"] for row in centre.rows()]
    # Every requirement, and every reference set under the database that holds it:
    # "a database is installed" says nothing about which of fifteen sets is there.
    assert [key for key in listed if not key.startswith("database:")] == list(ORDER)
    sets = [key for key in listed if key.startswith("database:")]
    assert len(sets) == len(centre.catalogue["entries"]) >= 2
    assert listed.index("database:" + centre.catalogue["entries"][0]["name"]) == \
        listed.index("hydra_database") + 1
    assert centre.table.rowCount() == len(listed)
    for index, row in enumerate(centre.rows()):
        for column in range(5):
            assert centre.table.item(index, column).text().strip(), (row["key"], column)
        assert centre.table.cellWidget(index, 5) is not None
        # Size is measured where something is installed and stated as a download
        # where it is not; a row that has neither says so rather than guessing.
        assert ("here" in row["size"] or "download" in row["size"]
                or row["size"].startswith("Not published")), row
    panel = next(row for row in centre.rows() if row["key"] == "species_panel")
    assert panel["installed"] == "Not installed"
    assert panel["action"] == "Install…" and panel["enabled"]


def test_nothing_is_compared_with_a_provider_until_the_online_check_is_pressed(window, qtbot):
    """"Not checked" and "up to date" are different sentences, and stay different."""
    centre = update_centre(window, qtbot)
    assert centre.check is None
    assert "Not checked is not the same as up to date" in centre.check_line.text()
    assert "No check has ever reached a provider" in centre.check_line.text()
    installed = [row for row in centre.rows() if row["installed"].startswith(("Installed",
                                                                             "Bundled"))]
    assert installed, "this workspace ships something, or the assertion below proves nothing"
    for row in installed:
        assert "up to date" not in row["update"].casefold()
        assert ("Not yet checked" in row["update"]
                or "cannot be confirmed" in row["update"]
                or "not from here" in row["update"]
                or "nothing here compares" in row["update"]
                or "publishes no version" in row["update"]), row
    # Reading this computer and asking a provider are two buttons with two words.
    assert "Rescan what is installed here" == centre.rescan_button.text()
    assert "Check online for updates" == centre.check_button.text()
    assert "No server is contacted" in centre.rescan_button.toolTip()
    assert "over the internet" in centre.check_button.toolTip()


def test_a_failed_check_says_so_and_never_reads_as_up_to_date(window, qtbot, monkeypatch):
    centre = update_centre(window, qtbot)
    fake_plan(monkeypatch, error="Could not reach NCBI to check for a newer release: no network",
              message="Could not reach NCBI to check for a newer release: no network")
    assert centre.check_online() is True
    qtbot.waitUntil(lambda: centre.check is not None, timeout=60000)
    qtbot.waitUntil(lambda: not centre.worker.isRunning(), timeout=60000)
    assert centre.check["ok"] is False
    line = centre.check_line.text()
    assert "did not get through" in line and "no network" in line
    assert "A failed check is not an up-to-date result." in line
    assert "No check has ever reached a provider from this computer." in line
    for row in centre.rows():
        assert "up to date" not in row["update"].casefold()
    databases = [row for row in centre.rows()
                 if row["key"].startswith("database:") and row["installed"] != "Not installed"]
    assert databases
    for row in databases:
        assert "The check failed" in row["update"] and "unknown" in row["update"]
    assert window.test_errors == []


def test_a_check_that_reaches_the_provider_is_remembered_for_the_next_failure(window, qtbot,
                                                                             monkeypatch):
    """A failed check has one thing left to say: when the answer was last confirmed."""
    from wmlstudio.update_center import check_record
    centre = update_centre(window, qtbot)
    fake_plan(monkeypatch, latest="2026-09-01.1", newer_available=True,
              message="NCBI publishes release 2026-09-01.1 today.")
    assert centre.check_online() is True
    qtbot.waitUntil(lambda: centre.check is not None and centre.check.get("ok"), timeout=60000)
    qtbot.waitUntil(lambda: not centre.worker.isRunning(), timeout=60000)
    assert "Checked online" in centre.check_line.text()
    assert "2026-09-01.1" in centre.check_line.text()
    remembered = check_record(window.root)
    assert remembered["latest"] == "2026-09-01.1" and remembered["when"]
    # And now one that does not get through: it names the check that did.
    fake_plan(monkeypatch, error="connection refused", message="connection refused")
    assert centre.check_online() is True
    qtbot.waitUntil(lambda: centre.check.get("ok") is False, timeout=60000)
    qtbot.waitUntil(lambda: not centre.worker.isRunning(), timeout=60000)
    assert remembered["when"] in centre.check_line.text()
    assert "connection refused" in centre.check_line.text()


def test_a_set_whose_terms_must_be_read_asks_before_anything_is_downloaded(window, qtbot,
                                                                           monkeypatch):
    """No licence is ever accepted on the user's behalf, and no set is silently skipped."""
    centre = update_centre(window, qtbot)
    restricted = next((row for row in centre.catalogue["entries"]
                       if not row["open_licence"] and not row["installed"]
                       and row["download"] == "automatic"), None)
    if restricted is None:
        restricted = next(row for row in centre.catalogue["entries"]
                          if not row["open_licence"] and not row["installed"])
    shown = next(row for row in centre.rows()
                 if row["key"] == "database:" + restricted["name"])
    assert "Needs your decision" in shown["update"]
    assert restricted["licence"] in shown["update"]
    assert "will not accept those terms for you" in shown["update"]
    assert shown["action"] == "Read the terms…"
    asked, told = answer(centre, monkeypatch, yes=False)
    launched = []
    monkeypatch.setattr(type(window), "launch_task",
                        lambda self, operation, role, completed=None: launched.append(role))
    assert centre.start("database:" + restricted["name"]) is False
    assert launched == [], "nothing may be fetched before the terms are answered"
    if restricted["download"] == "automatic":
        assert asked and restricted["licence"] in asked[0][1]
        assert "cannot accept these terms for you" in asked[0][1]
        assert "needs its provider's terms accepted first" in centre.status.text()
    else:
        assert told and "installed by hand" in told[0][1]


def test_install_everything_lists_the_licences_it_refuses_to_accept(window, qtbot, monkeypatch):
    from wmlstudio import provisioning
    from wmlstudio.update_center import UpdateCenter, decisions_needed
    update_centre(window, qtbot)
    fake_plan(monkeypatch)
    plan = provisioning.installation_plan(window.root, "")
    waiting = decisions_needed(plan)
    assert waiting, "the catalogue lists sets whose terms a person has to read"
    text = UpdateCenter.plan_text(plan)
    assert "Needs your decision, and is not included here:" in text
    for row in waiting:
        assert row["title"] in text
        assert row["name"] not in plan["databases"], "never fetched without an answer"
    assert "will not accept those terms for you" in text


def test_a_plan_from_a_failed_check_says_it_is_not_a_clean_bill_of_health(window, qtbot,
                                                                          monkeypatch):
    from wmlstudio import provisioning
    from wmlstudio.update_center import UpdateCenter
    update_centre(window, qtbot)
    fake_plan(monkeypatch, error="no route to host", message="no route to host")
    text = UpdateCenter.plan_text(provisioning.installation_plan(window.root, ""))
    assert text.startswith("The online check did not get through")
    assert "not a statement that you are up to date" in text


def test_an_update_row_uses_the_installer_the_application_already_has(window, qtbot,
                                                                     monkeypatch):
    """And says which of the two things happened, rather than assuming the good one."""
    centre = update_centre(window, qtbot)
    called = []
    monkeypatch.setattr(type(window), "install_species_panel",
                        lambda self: called.append("species"))
    assert centre.start("species_panel") is True
    assert called == ["species"]
    # This installer asked and was answered no — nothing started — so the page
    # says so rather than reporting an install that is not running.
    assert "Nothing has been installed yet" in centre.status.text()
    assert centre.working == "", "no row is mid-install when nothing was started"

    # The same press, with the installer actually starting work on the window's
    # own worker: now the page says what is being installed and locks the rows.
    def launches(self):
        called.append("launched")
        self.launch_task(lambda cancelled, progress: {"species_count": 18}, "species_panel")

    monkeypatch.setattr(type(window), "install_species_panel", launches)
    assert centre.start("species_panel") is True
    assert "Installing Species reference panel" in centre.status.text()
    qtbot.waitUntil(lambda: not window.worker.isRunning(), timeout=60000)
    settled(window, qtbot)
    # Finishing re-reads this computer by itself. The reported bug was a species
    # panel that really did install while its row went on saying "Install…".
    qtbot.waitUntil(lambda: centre.working == "", timeout=60000)


def test_installing_a_reference_set_is_never_announced_as_an_analysis(window, qtbot,
                                                                     monkeypatch):
    """The reported bug, on the page it was reported from.

    Pressing Install raised a dialog headed "Your analysis is running" over a
    download, which is how "when I update it performs analysis instead of looking
    online" survived the Update tab getting its own page.
    """
    centre = update_centre(window, qtbot)
    row = next((entry for entry in centre.catalogue["entries"]
                if entry.get("open_licence") and not entry.get("installed")
                and entry.get("download") == "automatic"), None)
    assert row is not None, "the catalogue lists a set this engine may fetch without asking"
    launched = []
    monkeypatch.setattr(
        type(window), "launch_task",
        lambda self, operation, role, completed=None, *, caption=None:
            bool(launched.append((role, caption))) or True)
    assert centre.start("database:" + row["name"]) is True
    role, caption = launched[0]
    assert role == "update:database:" + row["name"]
    assert caption is not None, "a download must not borrow the analysis wording"
    assert caption.title == "Installing reference data"
    assert row["title"] in caption.heading
    assert "not an analysis" in caption.note
    # While it runs, its own row says so and no other row can be pressed: two
    # downloads into the same store would not both survive.
    assert centre.working == "database:" + row["name"]
    installing_row = next(index for index, entry in enumerate(centre.rows())
                          if entry["key"] == centre.working)
    assert centre.table.cellWidget(installing_row, 5).text() == "Installing…"
    assert not any(centre.table.cellWidget(index, 5).isEnabled()
                   for index in range(centre.table.rowCount()))
    assert "re-reads itself" in centre.status.text()


def test_an_install_that_fails_hands_its_own_button_back(window, qtbot):
    """Clearing the flag only on success disabled the button for the whole session.

    "Install and update everything" records each step's failure rather than
    raising, but the call around them can still fail — and that path never
    reached the handler that re-enables the button, so the one control a person
    would press again was the one thing the failure took away.
    """
    centre = update_centre(window, qtbot)
    centre.installing = True
    centre.working = "everything"
    centre.everything_button.setEnabled(False)
    assert window.refresh_update_center() is True
    # Wait on the button, not on the thread: the worker stops first and the
    # handler that hands the button back runs on the next pass of the event loop.
    qtbot.waitUntil(lambda: centre.everything_button.isEnabled(), timeout=60000)
    settled(window, qtbot)
    assert centre.installing is False
    assert centre.working == ""
    assert all(centre.table.cellWidget(index, 5).isEnabled()
               for index in range(centre.table.rowCount())
               if centre.rows()[index]["enabled"])


def test_an_item_with_no_installer_here_explains_how_it_is_installed(window, qtbot,
                                                                     monkeypatch):
    centre = update_centre(window, qtbot)
    _asked, told = answer(centre, monkeypatch)
    # BLAST+ ships beside the application; nothing downloads it on a user's behalf.
    assert centre.start("blast_tools") is False
    assert told and "BLAST" in told[0][0]
    assert "stage_bio_tools" in told[0][1] or "Tools/blast" in told[0][1]


def test_clearing_the_update_tab_forgets_the_view_and_keeps_the_check_record(window, qtbot,
                                                                            monkeypatch):
    from wmlstudio.update_center import check_record
    centre = update_centre(window, qtbot)
    fake_plan(monkeypatch, latest="2026-09-01.1", newer_available=True, message="newer release")
    centre.check_online()
    qtbot.waitUntil(lambda: centre.check is not None and centre.check.get("ok"), timeout=60000)
    qtbot.waitUntil(lambda: not centre.worker.isRunning(), timeout=60000)
    monkeypatch.setattr(QMessageBox, "question",
                        lambda *args: QMessageBox.StandardButton.Yes)
    assert window.clear_page("update") is True
    assert centre.check is None and centre.plan is None
    qtbot.waitUntil(lambda: centre.report is not None, timeout=60000)
    qtbot.waitUntil(lambda: not centre.worker.isRunning(), timeout=60000)
    # The one thing a failed check has left to say survives a cleared view.
    assert check_record(window.root)["latest"] == "2026-09-01.1"
    assert clear_promise("update")["keeps"].endswith("when a check last reached a provider")
    assert window.test_errors == []
