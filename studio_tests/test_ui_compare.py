import json
import platform
import threading
from pathlib import Path

import pytest
from PySide6.QtCore import QPoint, QTimer, qVersion
from PySide6.QtWidgets import QApplication, QSizePolicy, QSplitter, QWidget

from wmlstudio import ui_compare
from wmlstudio.app import MainWindow
from wmlstudio.investigation import InvestigationStore
from wmlstudio.project import Project
from wmlstudio.sequence import check_cancelled


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    monkeypatch.setattr(MainWindow, 'error', lambda self, message: pytest.fail(str(message)))
    widget = MainWindow(storage_root=tmp_path / 'application')
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


def profile(window, index):
    sid = window.project.add_profile(f'isolate-{index}', {
        'sample_name': f'isolate-{index}', 'scheme': 'Example cgMLST', 'scheme_digest': 'same-reference',
        'status': 'profile_imported', 'alleles': {'a': str(1 + index % 2), 'b': '1'},
        'calls': [],
    }, {'organism': {'genus': 'Enterococcus', 'species': 'faecium'}})
    window.cohort_ids = set(window.cohort_ids or ()) | {sid}
    return sid


def test_project_comparison_revision_changes_for_results_metadata_and_secondary(tmp_path):
    project = Project(tmp_path / 'project.sqlite')
    try:
        assert project.comparison_revision() == ()
        sid = project.add_profile('one', {'scheme': 's', 'scheme_digest': 'd', 'alleles': {'a': '1'}})
        first = project.comparison_revision()
        assert len(first) == 1 and len(first[0]) == 4
        project.update_metadata(sid, {'ward': 'ICU'})
        second = project.comparison_revision()
        assert second != first
        project.set_analysis(sid, {'scheme': 'other', 'scheme_digest': 'other', 'alleles': {'a': '2'}})
        third = project.comparison_revision()
        assert third != second
        project.set_setting('unrelated_graph_zoom', 2)
        assert project.comparison_revision() == third
    finally:
        project.close()


def test_small_cohort_stays_synchronous_and_reuses_distances_for_layout(window, monkeypatch):
    profile(window, 0)
    profile(window, 1)
    calls = []
    original = ui_compare.pairwise_distances
    def measured(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)
    monkeypatch.setattr(ui_compare, 'pairwise_distances', measured)
    window.refresh_comparison()
    assert window.comparison_worker is None
    assert window.distance_rows[0]['distance'] == 1
    window.cluster_threshold.setValue(4)
    window.refresh_comparison()
    assert len(calls) == 1
    window.overlap.setValue(.90)
    window.refresh_comparison()
    assert len(calls) == 2


def test_large_cohort_runs_off_gui_and_new_request_supersedes_old(window, qtbot, monkeypatch):
    ids = [profile(window, index) for index in range(32)]
    entered = threading.Event()
    requests = []
    original = ui_compare._calculate_comparison
    def gated(project, request, cancelled=None):
        requests.append(request['chosen'])
        if len(requests) == 1:
            entered.set()
            while not cancelled():
                threading.Event().wait(.005)
            check_cancelled(cancelled)
        return original(project, request, cancelled)
    monkeypatch.setattr(ui_compare, '_calculate_comparison', gated)
    window.cohort_ids = set(ids[:-1])
    window.refresh_comparison()
    qtbot.waitUntil(entered.is_set, timeout=3000)
    ticks = []
    timer = QTimer(window)
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start(5)
    qtbot.waitUntil(lambda: len(ticks) >= 3, timeout=1000)
    window.cohort_ids = set(ids)
    window.refresh_comparison()
    qtbot.waitUntil(lambda: len(window._last_comparison) == 32, timeout=5000)
    qtbot.waitUntil(lambda: window.comparison_worker is None, timeout=3000)
    timer.stop()
    assert len(window.distance_rows) == 32 * 31 // 2
    assert len(requests) == 2
    assert window.worker is None  # Analysis queue and comparison use separate workers.


def test_close_cancels_background_comparison_before_database_closes(window, qtbot, monkeypatch):
    for index in range(31):
        profile(window, index)
    entered = threading.Event()
    def gated(project, request, cancelled=None):
        entered.set()
        while not cancelled():
            threading.Event().wait(.005)
        check_cancelled(cancelled)
    monkeypatch.setattr(ui_compare, '_calculate_comparison', gated)
    window.refresh_comparison()
    qtbot.waitUntil(entered.is_set, timeout=3000)
    assert window.cancel_comparison_for_close() is False
    qtbot.waitUntil(lambda: window.comparison_worker is None, timeout=3000)
    assert window._last_comparison == []


def test_installed_schemes_filter_taxa_using_metadata_only(window):
    profile(window, 0)
    folders = []
    for name, organism in [('match', 'Enterococcus faecium'), ('other', 'Klebsiella pneumoniae'), ('unknown', '')]:
        folder = window.root / 'schemes' / name
        folder.mkdir(parents=True)
        (folder / 'scheme.json').write_text(json.dumps({'name': name, 'organism': organism}))
        (folder / 'locus.fasta').write_text('not a valid FASTA; selectors must never parse allele data')
        folders.append(folder)
    window.scheme_paths = folders
    window.refresh_cohort_table()
    window.compare_genus.setCurrentIndex(window.compare_genus.findData('Enterococcus'))
    available = {window.compare_scheme.itemData(index) for index in range(window.compare_scheme.count())}
    assert str(folders[0]) in available
    assert str(folders[1]) not in available
    assert str(folders[2]) in available
    window.show_all_comparison_schemes.setChecked(True)
    assert window.compare_scheme.findData(str(folders[1])) >= 0


@pytest.mark.parametrize('metadata_file', ['scheme.json', 'example_info.json', None])
def test_repeated_scheme_labels_never_stat_glob_or_read(tmp_path, monkeypatch, metadata_file):
    folder = tmp_path / 'example_scheme'
    folder.mkdir()
    if metadata_file:
        (folder / metadata_file).write_text(json.dumps({'name': 'Reference A', 'organism': 'Enterococcus faecium'}))
    workspace = ui_compare.ComparisonWorkspaceMixin()
    workspace.invalidate_comparison_scheme_labels()
    expected = ('Reference A', 'Enterococcus', 'faecium') if metadata_file else ('example scheme', '', '')
    assert workspace.comparison_scheme_label(folder) == expected

    def forbidden(*args, **kwargs):
        raise AssertionError('Cached scheme labels must not touch the filesystem.')

    with monkeypatch.context() as patch:
        for operation in ('stat', 'glob', 'read_text', 'resolve'):
            patch.setattr(Path, operation, forbidden)
        for _ in range(20):
            assert workspace.comparison_scheme_label(str(folder)) == expected


def test_explicit_reference_refresh_invalidates_replaced_label_metadata(window):
    folder = window.root / 'schemes' / 'replaceable'
    folder.mkdir(parents=True)
    metadata = folder / 'scheme.json'
    metadata.write_text(json.dumps({'name': 'Old reference', 'organism': 'Enterococcus faecium'}))
    assert window.comparison_scheme_label(folder) == ('Old reference', 'Enterococcus', 'faecium')
    metadata.write_text(json.dumps({'name': 'New reference', 'organism': 'Klebsiella pneumoniae'}))
    # Ordinary sample/progress refreshes use the immutable snapshot cache.
    assert window.comparison_scheme_label(folder)[0] == 'Old reference'
    window.populate_schemes()
    assert window.comparison_scheme_label(folder) == ('New reference', 'Klebsiella', 'pneumoniae')
    metadata.write_text(json.dumps({'name': 'Revised reference'}))
    window.invalidate_comparison_scheme_labels(folder)
    assert window.comparison_scheme_label(folder)[0] == 'Revised reference'


def test_graph_legend_uses_escaped_labels_and_safe_color_swatches(window):
    window.update_graph_legend({'<b>Ward A</b>': '#ef1234', 'other': 'invalid css; background:url(secret)'})
    text = window.graph_legend.text()
    assert 'color:#ef1234' in text
    assert '&lt;b&gt;Ward A&lt;/b&gt;' in text
    assert 'background' not in text
    assert '●' in text


@pytest.mark.parametrize('size', [(1380, 940), (1080, 720)])
def test_comparison_fits_normal_desktop_sizes_without_outer_scroll(window, qtbot, size):
    profile(window, 0)
    profile(window, 1)
    window.resize(*size)
    window.navigate(2)
    window.refresh_comparison()
    assert_comparison_geometry(window, qtbot, f'{size[0]}x{size[1]}-normal')


def assert_comparison_geometry(window, qtbot, case):
    """Keep exact geometry evidence in CI artifacts, including on assertion failure."""
    page = window.pages.widget(2)
    window.refresh_cohort_table()
    # Allow queued height-for-width layouts and the graph's 80 ms resize-fit
    # timer to settle before measuring; an early zero scroll range can lie.
    qtbot.wait(120)
    try:
        qtbot.waitUntil(lambda: page.verticalScrollBar().maximum() == 0
                       and page.horizontalScrollBar().maximum() == 0
                       and window.tree.height() >= 200
                       and window.tree_status.height() >= window.tree_status.heightForWidth(window.tree_status.width()),
                       timeout=3000)
    finally:
        destination = Path('artifacts') / 'compare-geometry'
        destination.mkdir(parents=True, exist_ok=True)
        stem = destination / f'{platform.system().lower()}-{case}'
        content = page.widget()
        def dimensions(size):
            return [size.width(), size.height()]
        widgets = []
        for widget in [content, *content.findChildren(QWidget)]:
            if not widget.isVisible():
                continue
            position = widget.mapTo(content, QPoint(0, 0))
            widgets.append({
                'class': type(widget).__name__, 'name': widget.objectName(),
                'text': (widget.text()[:180] if callable(getattr(widget, 'text', None)) else ''),
                'position': [position.x(), position.y()], 'size': dimensions(widget.size()),
                'minimum': dimensions(widget.minimumSize()), 'hint': dimensions(widget.sizeHint()),
                'minimum_hint': dimensions(widget.minimumSizeHint()),
                'font': widget.font().toString(),
                'splitter_sizes': widget.sizes() if isinstance(widget, QSplitter) else None,
            })
        evidence = {
            'qt': qVersion(), 'platform': QApplication.platformName(),
            'window': dimensions(window.size()), 'viewport': dimensions(page.viewport().size()),
            'content': dimensions(content.size()), 'content_minimum': dimensions(content.minimumSizeHint()),
            'scroll_maximum': {'vertical': page.verticalScrollBar().maximum(),
                               'horizontal': page.horizontalScrollBar().maximum()},
            'graph': dimensions(window.tree.size()), 'widgets': widgets,
        }
        stem.with_suffix('.json').write_text(json.dumps(evidence, indent=2), encoding='utf-8')
        window.grab().save(str(stem.with_suffix('.png')))
        print(f'Comparison geometry: {stem}.json; scroll={evidence["scroll_maximum"]}; '
              f'viewport={evidence["viewport"]}; content={evidence["content"]}; graph={evidence["graph"]}')
    assert window.tree.height() >= 200


def test_comparison_fits_larger_font_metrics_without_outer_scroll(window, qtbot):
    profile(window, 0)
    profile(window, 1)
    page = window.pages.widget(2)
    # Stress font metrics independently of whichever Windows/Linux font is installed.
    page.setStyleSheet(STRESS_FONTS)
    window.resize(1080, 720)
    window.navigate(2)
    window.refresh_comparison()
    assert_comparison_geometry(window, qtbot, '1080x720-larger-font')


STRESS_FONTS = ('QPushButton, QToolButton, QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox, QTabBar '
                '{ font-size: 15px; } QLabel#small { font-size: 13px; } '
                'QLabel#cardTitle { font-size: 17px; }')


def pinned_baseline(window):
    """A saved investigation whose first build becomes the pinned baseline tree."""
    plan = InvestigationStore(window.project).save(
        'Geometry check', sorted(window.cohort_ids or ()), scheme='Example cgMLST',
        scheme_digest='same-reference', threshold=1, min_overlap=0.95,
        protocol='Synthetic test protocol; not a clinical cutoff')
    window.select_investigation(plan['id'])
    assert window._current_snapshot, window.tree_status.text()
    return plan['id']


@pytest.mark.parametrize('size,stress', [((1380, 940), False), ((1080, 720), False), ((1080, 720), True)])
def test_both_trees_keep_a_usable_height_at_normal_desktop_sizes(window, qtbot, size, stress):
    profile(window, 0)
    profile(window, 1)
    page = window.pages.widget(2)
    if stress:
        page.setStyleSheet(STRESS_FONTS)
    pinned_baseline(window)
    window.resize(*size)
    window.navigate(2)
    window.dual_toggle.setChecked(True)
    window.refresh_comparison()
    case = f'{size[0]}x{size[1]}-dual' + ('-larger-font' if stress else '')
    assert_comparison_geometry(window, qtbot, case)
    assert window.baseline_pane.isVisible() and window.baseline_tree.isVisible()
    assert window.baseline_tree.height() >= 200, window.baseline_tree.height()
    assert window.baseline_caption.isVisible() and window.current_caption.isVisible()
    assert window.baseline_caption.text().startswith('Baseline · ')
    # A caption must never be what pushes the page into a horizontal scroll: the
    # size policy ignores its preferred width, so a long title clips instead.
    assert window.baseline_caption.sizePolicy().horizontalPolicy() == QSizePolicy.Policy.Ignored


def test_turning_the_baseline_off_restores_the_single_tree_layout(window, qtbot):
    profile(window, 0)
    profile(window, 1)
    pinned_baseline(window)
    window.resize(1080, 720)
    window.navigate(2)
    window.dual_toggle.setChecked(True)
    qtbot.wait(60)
    window.dual_toggle.setChecked(False)
    assert not window.baseline_pane.isVisible()
    assert not window.current_caption.isVisible()
    assert not window.export_tree_choice.isVisible()
    assert window.baseline_tree._results == {}
    assert_comparison_geometry(window, qtbot, '1080x720-dual-off')


# ---------------------------------------------------------------------------
# One tree in a window of its own, and starting the tab again
# ---------------------------------------------------------------------------

CG_TARGETS = [f'target{index:03d}' for index in range(1, 41)]


def cg_profile(window, index):
    """One core-genome profile: forty targets, so its typing kind is its own evidence."""
    alleles = {locus: '1' for locus in CG_TARGETS}
    alleles['target001'] = str(1 + index % 2)
    sid = window.project.add_profile(f'core-{index}', {
        'sample_name': f'core-{index}', 'scheme': 'Demo cgMLST', 'scheme_digest': 'core-reference',
        'status': 'profile_imported', 'alleles': alleles, 'calls': [], 'input_sha256': 'a' * 64,
        'analysis_kind': 'cgmlst'}, {'organism': {'genus': 'Klebsiella', 'species': 'pneumoniae'}})
    window.cohort_ids = set(window.cohort_ids or ()) | {sid}
    return sid


def built(window, kind='cgmlst'):
    window.refresh_cohort_table()
    window.show_typing_view(kind)
    window.refresh_comparison()
    assert window.tree.nodes, window.tree_status.text()
    return window.tree


def test_the_drawn_tree_opens_in_its_own_window_and_leaves_the_page_untouched(window, qtbot):
    cg_profile(window, 0)
    cg_profile(window, 1)
    built(window)
    detached = window.open_graph_window()
    qtbot.addWidget(detached)
    try:
        assert window._graph_windows == [detached]
        assert detached.view is not window.tree
        assert sorted(detached.view.nodes) == sorted(window.tree.nodes)
        positions = {key: (node.pos().x(), node.pos().y())
                     for key, node in window.tree.nodes.items()}
        moved = sorted(detached.view.nodes)[0]
        detached.view.nodes[moved].setPos(321, 123)
        detached.view.update_edges()
        detached.view.set_node_colors({moved: '#ff0055'})
        assert {key: (node.pos().x(), node.pos().y())
                for key, node in window.tree.nodes.items()} == positions
        assert window.tree.nodes[moved].slices != [('#ff0055', 1)]
    finally:
        detached.close()
    assert window._graph_windows == []
    assert sorted(window.tree.nodes) == sorted(positions)


def test_a_window_states_the_typing_kind_reference_and_target_count_it_is_showing(window, qtbot):
    cg_profile(window, 0)
    cg_profile(window, 1)
    built(window)
    window.cluster_threshold.setValue(5)
    detached = window.open_graph_window()
    qtbot.addWidget(detached)
    try:
        assert detached.identity.kind == 'cgmlst'
        assert detached.identity.targets == len(CG_TARGETS)
        assert detached.identity.caption() == 'cgMLST · Demo cgMLST · 40 targets'
        assert detached.identity.threshold_words() == 'link ≤ 5 of 40 targets'
        assert '2 isolates' in detached.identity.cohort_words()
        title = detached.windowTitle()
        assert title.startswith('cgMLST forest · link ≤ 5 of 40 targets · 2 isolates · Demo cgMLST')
        assert detached.headline.text() == 'cgMLST · Demo cgMLST · 40 targets'
    finally:
        detached.close()


def test_a_classical_window_and_a_core_genome_window_never_share_a_scale(window, qtbot):
    cg_profile(window, 0)
    cg_profile(window, 1)
    profile(window, 0)
    profile(window, 1)
    built(window, 'cgmlst')
    core = window.open_graph_window()
    qtbot.addWidget(core)
    built(window, 'mlst')
    classical = window.open_graph_window()
    qtbot.addWidget(classical)
    try:
        assert (core.identity.kind, classical.identity.kind) == ('cgmlst', 'mlst')
        assert core.identity.targets == 40 and classical.identity.targets == 2
        assert core.identity.target_word == 'targets' and classical.identity.target_word == 'loci'
        assert core.windowTitle() != classical.windowTitle()
        assert sorted(core.view.nodes) != sorted(classical.view.nodes)
        assert core.identity.separation == classical.identity.separation
        assert 'never share a scale' in core.identity.separation
    finally:
        core.close()
        classical.close()


def test_opening_a_window_before_anything_is_drawn_opens_nothing_and_says_so(window):
    assert window.open_graph_window() is None
    assert window._graph_windows == []
    assert 'nothing drawn' in window.progress_text.text()


def test_a_tree_that_is_not_on_screen_is_not_offered_as_a_second_opinion(window):
    cg_profile(window, 0)
    cg_profile(window, 1)
    built(window)
    assert window.open_graph_window('baseline') is None
    assert window._graph_windows == []
    assert 'not on screen' in window.progress_text.text()


def test_a_window_opened_from_the_baseline_tree_names_the_stored_snapshot(window, qtbot):
    profile(window, 0)
    profile(window, 1)
    pinned_baseline(window)
    window.dual_toggle.setChecked(True)
    qtbot.waitUntil(lambda: bool(window.baseline_tree.nodes), timeout=5000)
    detached = window.open_graph_window('baseline')
    qtbot.addWidget(detached)
    try:
        assert detached.identity.cohort == 'Geometry check'
        assert detached.identity.scheme == 'Example cgMLST'
        assert 'stored baseline snapshot' in detached.identity.note
        assert sorted(detached.view.nodes) == sorted(window.baseline_tree.nodes)
        assert detached.view is not window.baseline_tree
    finally:
        detached.close()


def test_clearing_the_compare_tab_empties_it_and_keeps_every_stored_profile(window, qtbot):
    first = cg_profile(window, 0)
    cg_profile(window, 1)
    built(window)
    window.cohort_search.setText('core')
    window.graph_search.setText('core-0')
    detached = window.open_graph_window()
    qtbot.addWidget(detached)
    window.clear_compare_tab()
    assert window.cohort_ids == set()
    assert window.project.get_setting('comparison_cohort') == []
    assert window._last_comparison == [] and window.tree.nodes == {}
    assert window._pending_graph_state is None
    assert window.cohort_search.text() == '' and window.graph_search.text() == ''
    assert window.active_investigation_id is None
    assert window._graph_windows == []
    assert window.cgmlst_calls.calls is None
    assert 'Nothing stored was changed' in window.tree_status.text()
    # Everything that was cleared was on screen. The evidence is untouched.
    assert len(window.project.samples()) == 2
    assert window.project.get_sample(first)['result']['alleles']['target001'] == '1'
    window.cohort_ids = {first}
    window.refresh_comparison()
    assert len(window.tree.nodes) == 1


def dual_typed(window, index):
    """One isolate carrying both a classical profile and a core-genome profile."""
    sid = cg_profile(window, index)
    window.project.set_analysis(sid, {
        'sample_name': f'core-{index}', 'scheme': 'Study MLST', 'scheme_digest': 'study-reference',
        'status': 'complete', 'alleles': dict(zip('abcdefg', '111111' + str(1 + index % 2))),
        'calls': [], 'st': '20', 'input_sha256': 'a' * 64, 'analysis_kind': 'mlst'})
    return sid


def test_both_typing_trees_open_in_their_own_windows_on_their_own_scales(window, qtbot):
    dual_typed(window, 0)
    dual_typed(window, 1)
    built(window, 'cgmlst')
    window.counterpart_toggle.setChecked(True)
    qtbot.waitUntil(lambda: bool(window.counterpart_tree.nodes), timeout=5000)
    core = window.open_graph_window('current')
    qtbot.addWidget(core)
    classical = window.open_graph_window('counterpart')
    qtbot.addWidget(classical)
    try:
        assert (core.identity.kind, classical.identity.kind) == ('cgmlst', 'mlst')
        assert (core.identity.targets, classical.identity.targets) == (40, 7)
        assert core.identity.scheme == 'Demo cgMLST' and classical.identity.scheme == 'Study MLST'
        assert 'of 40 targets' in core.identity.threshold_words()
        assert 'of 7 loci' in classical.identity.threshold_words()
        assert 'the other typing view' in classical.identity.note
        assert len(window._graph_windows) == 2
        # Two windows on the same isolates, and neither one can borrow the other's
        # arrangement or the other's numbers.
        assert core.view is not classical.view
        assert sorted(core.view.nodes) == sorted(classical.view.nodes)
    finally:
        core.close()
        classical.close()
    assert window._graph_windows == []
