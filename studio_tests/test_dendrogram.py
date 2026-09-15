"""The dendrogram on screen: what it draws, and what it must never leave unsaid."""

import pytest
from PySide6.QtCore import QRectF

from wmlstudio.app import MainWindow
from wmlstudio.dendrogram import DendrogramPanel, DendrogramView
from wmlstudio.hierarchical import build_tree

LOCI = [f'locus{index:02d}' for index in range(10)]


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    monkeypatch.setattr(MainWindow, 'error', lambda self, message: pytest.fail(str(message)))
    widget = MainWindow(storage_root=tmp_path / 'application')
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None
                    and widget.dendrogram_worker is None, timeout=5000)
    widget.close()


def rows(distances):
    """Pairwise rows in the shape the comparison page hands to the clustering core."""
    pairs = []
    for (left, right), value in sorted(distances.items()):
        pairs.append({'source': left, 'target': right, 'distance': value,
                      'comparable': value is not None,
                      'reason': '' if value is not None else 'Below the minimum overlap.'})
    return pairs


def add(window, name, alleles, *, organism=('Klebsiella', 'pneumoniae'), scheme='Demo cgMLST',
        digest='core-reference', kind='cgmlst'):
    """One stored profile with exactly the alleles a test wants to compare."""
    sid = window.project.add_profile(name, {
        'sample_name': name, 'scheme': scheme, 'scheme_digest': digest,
        'status': 'profile_imported', 'alleles': dict(alleles), 'calls': [],
        'input_sha256': 'a' * 64, 'analysis_kind': kind},
        {'organism': {'genus': organism[0], 'species': organism[1]}})
    window.cohort_ids = set(window.cohort_ids or ()) | {sid}
    return sid


def spread(prefix_changes):
    """A ten-locus profile whose first ``prefix_changes`` loci carry allele 2."""
    return {locus: ('2' if index < prefix_changes else '1') for index, locus in enumerate(LOCI)}


def four_isolates(window):
    """a-b 1, a-c 2, b-c 1, a-d 5, b-d 4, c-d 3: every linkage gives a different tree."""
    add(window, 'iso-a', spread(0))
    add(window, 'iso-b', spread(1))
    add(window, 'iso-c', spread(2))
    add(window, 'iso-d', spread(5))


def opened(window, kind='cgmlst'):
    """Build the comparison and put the dendrogram tab in front, as a user would."""
    window.refresh_cohort_table()
    window.show_typing_view(kind)
    window.refresh_comparison()
    window.graph_tabs.setCurrentWidget(window.dendrogram)
    return window.dendrogram


def heights(panel):
    return [merge['height'] for merge in panel.tree().get('merges', [])]


# ---------------------------------------------------------------------------
# The widget itself
# ---------------------------------------------------------------------------

def test_the_widget_draws_every_leaf_and_join_of_a_known_small_tree(qtbot):
    """A dendrogram that silently draws nothing would read as a cohort with no
    structure, so the geometry and the paint are both checked against a tree whose
    three leaves and two joins are known by hand."""
    view = DendrogramView()
    qtbot.addWidget(view)
    view.resize(640, 420)
    tree = build_tree(rows({('a', 'b'): 2, ('a', 'c'): 6, ('b', 'c'): 5}), linkage='single',
                      scale_caption='cgMLST · Demo cgMLST · 40 targets')
    view.set_tree(tree, names={'a': 'iso-a', 'b': 'iso-b', 'c': 'iso-c'}, cut=2)
    plan = view.plan(QRectF(0, 0, 640, 420))
    assert [leaf['name'] for leaf in plan['leaves']] == ['iso-a', 'iso-b', 'iso-c']
    assert len(plan['branches']) == 2
    # The join at 2 sits left of the join at 5, and both sit right of the leaves.
    joins = sorted(branch['x'] for branch in plan['branches'])
    assert plan['leaves'][0]['x'] < joins[0] < joins[1]
    assert plan['cut']['height'] == 2
    assert not view.grab().isNull()


def test_the_height_axis_is_captioned_with_the_quantity_the_tree_was_built_from(qtbot):
    """A height of 5 on seven loci and a height of 5 on 2 358 targets are different
    numbers. An axis without its quantity is how one is read as the other."""
    view = DendrogramView()
    qtbot.addWidget(view)
    view.set_tree(build_tree(rows({('a', 'b'): 2}), scale_caption='cgMLST · Demo · 2358 targets'))
    assert view.axis_caption() == 'cgMLST · Demo · 2358 targets'
    view.set_tree(build_tree(rows({('a', 'b'): 2})))
    assert 'not recorded' in view.axis_caption()


def test_dragging_the_cut_line_regroups_the_drawn_clusters(qtbot):
    """The cut is the control the reader actually uses; if it did not change the
    clusters under it the picture would disagree with the groups it claims to show."""
    panel = DendrogramPanel()
    qtbot.addWidget(panel)
    panel.resize(700, 500)
    panel.show_tree(build_tree(rows({('a', 'b'): 1, ('b', 'c'): 1, ('a', 'c'): 2}),
                               linkage='single', scale_caption='cgMLST · Demo · 40 targets'),
                    cut=0)
    assert panel.view.groups()['clusters'] == 0
    panel.view.set_cut(1)
    assert panel.view.groups()['clusters'] == 1
    assert panel.cut.value() == 1
    coloured = {leaf['colour'] for leaf in panel.view.plan(QRectF(0, 0, 700, 500))['leaves']}
    assert len(coloured) == 1 and coloured != {'#8b93a7'}


# ---------------------------------------------------------------------------
# The comparison page
# ---------------------------------------------------------------------------

def test_the_dendrogram_tab_sits_beside_the_graph_and_groups_tabs(window):
    """The user asked for hierarchical clustering where the trees already are; a
    tab that is not on the comparison page is a feature nobody finds."""
    titles = [window.graph_tabs.tabText(index) for index in range(window.graph_tabs.count())]
    assert titles[:3] == ['Graph', 'Groups', 'Dendrogram']
    assert [window.dendrogram_linkage.itemData(index)
            for index in range(window.dendrogram_linkage.count())] == ['single', 'complete', 'average']


def test_choosing_another_linkage_changes_the_heights_the_merges_are_drawn_at(window):
    """Single, complete and average linkage answer different questions about the same
    distances. A selector that did not change the tree would be a decoration."""
    four_isolates(window)
    panel = opened(window)
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('single'))
    assert heights(panel) == [1, 1, 3]
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('complete'))
    assert heights(panel) == [1, 2, 5]
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('average'))
    assert heights(panel) == [1, 1.5, 4]
    assert panel.tree()['linkage'] == 'average'


def test_a_cluster_held_together_by_a_chain_says_so_under_the_picture(window):
    """Single-linkage chaining is the classic way an outbreak is invented: a and c
    differ by 2 but land in one group cut at 1. The reader has to be told."""
    add(window, 'iso-a', spread(0))
    add(window, 'iso-b', spread(1))
    add(window, 'iso-c', spread(2))
    window.cluster_threshold.setValue(1)
    panel = opened(window)
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('single'))
    notes = panel.notes_text()
    assert 'Chaining' in notes
    assert 'furthest compared pair differs by 2' in notes
    assert panel.view.groups()['chained_groups']


def test_clusters_kept_apart_for_want_of_evidence_are_named_as_such(window):
    """A pair nobody could compare must never disappear into the picture. The tree
    refuses the join, and the page has to say it was refused for missing evidence."""
    add(window, 'iso-a', spread(0))
    add(window, 'iso-b', spread(1))
    # Its own extra targets halve the shared fraction, so every pair with it falls
    # below the minimum overlap and carries no distance at all.
    add(window, 'iso-d', {**spread(0), **{f'extra{index}': '1' for index in range(10)}})
    panel = opened(window)
    notes = panel.notes_text()
    assert 'could not be joined at any height' in notes
    assert 'never compared' in notes
    assert 'not because they were measured as far apart' in notes
    assert panel.tree()['blocked_merges']


def test_a_classical_tree_and_a_core_genome_tree_share_no_caption_and_no_cut(window):
    """The two typing views are different quantities. One caption or one cut line
    carried across them would put a seven-locus height on a core-genome axis."""
    for index in range(2):
        sid = add(window, f'both-{index}', spread(index))
        window.project.set_analysis(sid, {
            'sample_name': f'both-{index}', 'scheme': 'Study MLST', 'scheme_digest': 'study-reference',
            'status': 'complete', 'alleles': dict(zip('abcdefg', '111111' + str(1 + index))),
            'calls': [], 'st': '20', 'input_sha256': 'a' * 64, 'analysis_kind': 'mlst'})
    # Each view keeps its own threshold, so each is set after its own view is open.
    opened(window, 'cgmlst')
    window.cluster_threshold.setValue(3)
    core = window.dendrogram.tree()['scale_caption']
    core_cut = window.dendrogram.cut_height()
    opened(window, 'mlst')
    window.cluster_threshold.setValue(1)
    classical = window.dendrogram.tree()['scale_caption']
    assert 'cgMLST' in core and 'targets' in core
    assert 'MLST' in classical and 'loci' in classical
    assert core != classical
    assert (core_cut, window.dendrogram.cut_height()) == (3.0, 1.0)


def test_the_cut_and_the_group_threshold_are_one_number_until_someone_moves_it(window):
    """Two cluster lists with no explanation is worse than one. The cut opens on the
    page's own Group ≤, and once it is moved each list names the method behind it."""
    four_isolates(window)
    panel = opened(window)
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('single'))
    window.cluster_threshold.setValue(2)
    assert panel.cut_height() == 2
    assert 'Both lists place every isolate the same way.' in panel.notes_text()
    panel.view.set_cut(4)
    notes = panel.notes_text()
    assert 'single linkage, cut at 4' in notes
    assert 'The Groups tab: single linkage, link ≤ 2' in notes
    assert 'placed differently by the two lists' in notes


def test_an_average_linkage_cut_is_never_offered_as_a_group_threshold(window):
    """Group ≤ is a single-linkage rule. Handing it an average-linkage height would
    relabel one method's number as another's, which is the disagreement to avoid."""
    four_isolates(window)
    panel = opened(window)
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('average'))
    panel.view.set_cut(1.5)
    assert not panel.adopt.isEnabled()
    assert 'different method' in panel.adopt.toolTip()
    window.dendrogram_linkage.setCurrentIndex(window.dendrogram_linkage.findData('single'))
    panel.view.set_cut(3)
    assert panel.adopt.isEnabled()
    panel.adopt.click()
    assert window.cluster_threshold.value() == 3


def test_a_cohort_with_no_comparison_says_nothing_was_clustered(window):
    """An empty picture must never be readable as "nothing was found"; it has to
    say the analysis was not run."""
    panel = opened(window)
    assert panel.tree() == {}
    assert 'has not been run' in panel.notes_text()
    assert 'unrelated' in panel.notes_text()


def test_a_large_cohort_is_clustered_off_the_gui_thread(window, qtbot):
    """A 300-isolate cohort is about three seconds of clustering. Doing it on the GUI
    thread would freeze the window while it ran."""
    for index in range(34):
        add(window, f'iso-{index:02d}', spread(index % 6))
    window.refresh_cohort_table()
    window.show_typing_view('cgmlst')
    window.refresh_comparison()
    qtbot.waitUntil(lambda: len(window._last_comparison) == 34, timeout=10000)
    window.graph_tabs.setCurrentWidget(window.dendrogram)
    assert window.dendrogram_worker is not None
    assert 'in the background' in window.dendrogram.caption.text()
    qtbot.waitUntil(lambda: window.dendrogram_worker is None, timeout=15000)
    assert len(window.dendrogram.tree()['labels']) == 34


# ---------------------------------------------------------------------------
# The user's own operational cutoff
# ---------------------------------------------------------------------------

KP_TARGETS = [f'target{index:04d}' for index in range(2358)]


def klebsiella(window, index):
    """One profile on the 2 358-target K. pneumoniae scheme the local cutoff names."""
    alleles = {locus: '1' for locus in KP_TARGETS}
    alleles['target0000'] = str(1 + index)
    return add(window, f'kp-{index}', alleles, organism=('Klebsiella', 'pneumoniae'),
               scheme='cgmlst.org K. pneumoniae', digest='kp-2358')


def test_the_group_threshold_is_seeded_from_your_own_cutoff_and_never_from_a_publication(window):
    """The 5 this laboratory runs on has no citation behind it. It must reach the
    screen through the "your own setting" path, or a report would print it as
    published evidence that does not exist."""
    klebsiella(window, 0)
    klebsiella(window, 1)
    window.refresh_cohort_table()
    window.show_typing_view('cgmlst')
    window.refresh_comparison()
    assert window.cluster_threshold.value() == 5
    banner = window.operational_banner.text()
    assert not window.operational_row.isHidden()
    assert banner.startswith('Your own operational cutoff for Klebsiella pneumoniae: at most 5')
    assert 'no publication reviewed here establishes it' in banner
    assert 'seeded' in banner
    # Nothing here may travel by the published route: no citation was recorded,
    # and the published banner has nothing to say about this number.
    assert window.threshold_evidence == {}
    assert window.guidance_row.isHidden()
    assert 'operational' not in window.guidance_banner.text()
    assert 'locally declared operational cutoff' in window.cluster_threshold.toolTip()


def test_the_whole_provenance_of_your_own_cutoff_is_on_the_page(window, monkeypatch):
    """The banner has room for the first claim only. The rest — who declared it, that
    a targeted search found no publication for it, and what a stricter cutoff costs —
    must still be reachable without leaving the page."""
    klebsiella(window, 0)
    klebsiella(window, 1)
    window.refresh_cohort_table()
    window.show_typing_view('cgmlst')
    window.refresh_comparison()
    shown = {}
    monkeypatch.setattr('wmlstudio.ui_compare.QMessageBox.information',
                        lambda parent, title, text: shown.update(title=title, text=text))
    window.show_operational_cutoff()
    assert shown['title'] == 'Your own operational cutoff'
    assert 'No citation was supplied with it.' in shown['text']
    assert 'not a published one' in shown['text']
    assert 'A stricter cutoff is not automatically the safer one' in shown['text']
    assert 'not a published one' in window.operational_banner.toolTip()


def test_a_threshold_somebody_already_chose_is_never_overwritten_by_the_cutoff(window):
    """A cutoff that is offered is not a cutoff that is applied. Once this view has a
    threshold of its own, the local number is stated beside it and left alone."""
    window.show_typing_view('cgmlst')
    window.cluster_threshold.setValue(12)
    klebsiella(window, 0)
    klebsiella(window, 1)
    window.refresh_cohort_table()
    window.refresh_comparison()
    assert window.cluster_threshold.value() == 12
    assert 'Group ≤ is currently 12, so this cutoff is offered here and is not applied.' \
        in window.operational_banner.text()


def test_a_cohort_of_another_organism_is_offered_no_local_cutoff(window):
    """The cutoff is bound to one taxon and one target count. Offering it anywhere
    else would rescale a number that was never measured on that scheme."""
    add(window, 'ef-0', spread(0), organism=('Enterococcus', 'faecium'))
    add(window, 'ef-1', spread(1), organism=('Enterococcus', 'faecium'))
    window.refresh_cohort_table()
    window.show_typing_view('cgmlst')
    window.refresh_comparison()
    assert window.operational_row.isHidden()
    assert window.cluster_threshold.value() == 1
