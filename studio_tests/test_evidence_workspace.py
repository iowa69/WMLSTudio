import csv

import pytest
from PySide6.QtWidgets import QDialog

from wmlstudio.app import MainWindow
from wmlstudio.cohort_picker import CohortPickerDialog


@pytest.fixture
def window(qtbot, tmp_path):
    window = MainWindow(storage_root=tmp_path / 'app')
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(30)
    yield window
    window.close()


def sample(window, name, gene):
    digest = 'a' * 64
    return window.project.add_profile(name, {'sample_name': name, 'st': '20', 'input_sha256': digest,
        'scheme': 'MLST', 'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}},
        {'organism': {'genus': 'Staphylococcus', 'species': 'aureus'},
         'annotations': {'ward': 'ICU'}, 'hydra': {'evidence_input_sha256': digest,
         'source_sample': name, 'hits': [{'gene': gene, 'element_type': 'AMR', 'primary': True}]}})


def test_evidence_cohort_is_independent_and_exports_only_visible_columns(window, tmp_path, monkeypatch):
    first, second = sample(window, 'A', 'blaA'), sample(window, 'B', 'vanB')
    window.selection_ids = {second}
    window.report_ids = {second}
    window.refresh_features()
    assert window.feature_model.rows == [] and window.amr_model.rows == []
    def choose(dialog):
        assert dialog.selected_ids == set()
        dialog.selected_ids = {first}
        return QDialog.DialogCode.Accepted
    monkeypatch.setattr(CohortPickerDialog, 'exec', choose)
    window.choose_feature_cohort()
    assert window.feature_ids == {first}
    assert window.selection_ids == window.report_ids == {second}
    assert [row['_sample_id'] for row in window.feature_model.rows] == [first]
    assert 'blaA' in window.amr_model.headers and 'vanB' not in window.amr_model.headers
    opened = []
    monkeypatch.setattr(window, 'open_isolate_record', opened.append)
    window.feature_table.doubleClicked.emit(window.feature_model.index(0, 0))
    assert opened == [first]
    path = tmp_path / 'evidence.tsv'
    window.write_feature_table(path, 'amr')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1 and rows[0]['sample_id'] == first
    assert second not in path.read_text() and 'vanB' not in path.read_text()
    window.gene_filter.setText('not-present')
    window.write_feature_table(path, 'amr')
    assert 'blaA' not in path.read_text()


def test_right_click_in_the_evidence_tables_routes_isolates_without_moving_that_cohort(window):
    from wmlstudio.context_menus import SEPARATOR, Selection
    first, second = sample(window, 'A', 'blaA'), sample(window, 'B', 'vanB')
    window.feature_ids = {first, second}
    window.refresh_features()
    assert {'evidence.features', 'evidence.amr'} <= set(window._context_adapters)

    for view, model in (('evidence.features', window.feature_model), ('evidence.amr', window.amr_model)):
        adapter = window._context_adapters[view]
        # The ids come from the model row, never from the visual row number.
        assert [adapter.id_at(row) for row in range(len(model.rows))] == [row['_sample_id'] for row in model.rows]
        selection = Selection(view, (first,))
        titles = [entry.format_title(selection) for entry in window.context_menu_plan(selection)
                  if entry is not SEPARATOR]
        assert 'Add to report' in titles and 'Copy sample ID' in titles

    window.context_add_to_report(Selection('evidence.features', (first,)))
    assert window.report_ids == {first}
    assert window.feature_ids == {first, second}, 'routing one isolate must not edit the evidence cohort'
    assert {s['id'] for s in window.project.samples()} == {first, second}


def test_advanced_source_controls_are_not_in_main_cohort_view(window):
    tabs = window.evidence_tabs
    assert any('Advanced' in tabs.tabText(index) for index in range(tabs.count()))
    assert not window.hydra_table.isVisible()
    window.navigate(4)
    assert not window.hydra_table.isVisible()
    assert window.feature_scope_label.isVisible()


def characterized(window, name, *, replicons=(), links=(), elsewhere=(), stale=False):
    """A sample whose characterization carries plasmid evidence for this assembly."""
    digest = 'a' * 64
    characterization = {
        'input_sha256': 'b' * 64 if stale else digest,
        'plasmid_hypotheses': {
            'status': 'completed', 'reason': 'Replicon assay was supplied.',
            'replicons': [{'gene': gene, 'element_type': 'PLASMID'} for gene in replicons],
            'contig_associations': [{'replicon': replicon, 'marker': marker, 'marker_type': 'AMR',
                                     'contig': 'Contig_2'} for replicon, marker in links],
            'determinant_placement': [{'gene': marker, 'marker_type': 'AMR', 'contig': 'Contig_2',
                                       'placement': 'co_located_with_replicon'}
                                      for _, marker in links]
            + [{'gene': marker, 'marker_type': 'AMR', 'contig': 'Contig_9',
                'placement': 'no_replicon_on_this_contig'} for marker in elsewhere]}}
    return window.project.add_profile(name, {'sample_name': name, 'st': '20', 'input_sha256': digest,
        'scheme': 'MLST', 'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}},
        {'organism': {'genus': 'Klebsiella', 'species': 'pneumoniae'},
         'characterization': characterization})


def table_rows(table):
    return [[table.item(row, column).text() for column in range(table.columnCount())]
            for row in range(table.rowCount())]


def test_the_cohort_plasmid_table_counts_only_the_isolates_whose_assay_actually_ran(window):
    """An isolate with no plasmid assay is named, never counted as carrying nothing."""
    carrier = characterized(window, 'A', replicons=('IncFIB',), links=(('IncFIB', 'blaKPC-2'),))
    apart = characterized(window, 'B', replicons=('IncFIB',), elsewhere=('blaKPC-2',))
    outdated = characterized(window, 'C', replicons=('IncFIB',), stale=True)
    window.feature_ids = {carrier, apart, outdated}
    window.refresh_features()

    replicons = table_rows(window.plasmid_replicon_table)
    assert [row[:3] for row in replicons] == [['IncFIB', '2', '2']], 'the denominator travels'
    assert 'A' in replicons[0][3] and 'B' in replicons[0][3] and 'C' not in replicons[0][3]
    # Co-located and merely co-present are separate columns: collapsing them
    # would turn an assembly artefact into a shared finding.
    pairs = table_rows(window.plasmid_pair_table)
    assert pairs == [['IncFIB', 'blaKPC-2', 'AMR', '1', '1', '2']]
    assert '2 of 3 isolates' in window.plasmid_scope.text()
    assert 'not a distance and not a cluster' in window.plasmid_scope.text()
    notes = window.plasmid_notes.toPlainText()
    assert 'C' in notes and 'earlier or different assembly' in notes
    assert 'not MOB-suite' in notes and 'not transmission' in notes


def test_the_plasmid_cohort_tab_exports_its_two_tables_as_two_files(window, tmp_path):
    carrier = characterized(window, 'A', replicons=('IncFIB',), links=(('IncFIB', 'blaKPC-2'),))
    window.feature_ids = {carrier}
    window.refresh_features()
    tabs = window.evidence_tabs
    assert any('Plasmid evidence across the cohort' == tabs.tabText(index)
               for index in range(tabs.count()))

    path = tmp_path / 'replicons.tsv'
    window.write_feature_table(path, 'plasmid_replicons')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert rows[0]['Replicon marker'] == 'IncFIB' and rows[0]['Of isolates assayed'] == '1'

    path = tmp_path / 'pairs.tsv'
    window.write_feature_table(path, 'plasmid_pairs')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert rows[0]['Determinant'] == 'blaKPC-2'
    assert rows[0]['Co-located on one contig'] == '1'
    assert rows[0]['Both present, not co-located'] == '0'


def test_the_evidence_summary_names_the_reference_sets_and_searches_behind_each_row(window):
    """A gene list means nothing without the question it answered."""
    digest = 'a' * 64
    sample_id = window.project.add_profile('A', {'sample_name': 'A', 'st': '20', 'input_sha256': digest,
        'scheme': 'MLST', 'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}},
        {'hydra': {'evidence_input_sha256': digest, 'source_sample': 'A',
                   'databases': ['ncbi', 'protein'],
                   'hits': [{'gene': 'blaKPC-2', 'element_type': 'AMR', 'primary': True}],
                   'execution_provenance': {
                       'organism': {'resolved': 'Escherichia', 'point_mutation_level': 'protein_only'},
                       'virulence': {'enabled': True, 'organism_curated': True},
                       'reference_release': {'release': '2026-01-01.1'}}}})
    window.feature_ids = {sample_id}
    window.refresh_features()
    row = window.feature_model.rows[0]

    assert row['AMR reference sets'] == 'ncbi; protein'
    assert row['AMR reference release'] == '2026-01-01.1'
    assert 'no DNA-level target was assessed' in row['Point mutations searched']
    assert row['Virulence elements searched'] == 'yes, curated for Escherichia'


def test_an_isolate_with_no_current_amr_evidence_says_so_in_every_method_column(window):
    sample_id = sample(window, 'A', 'blaA')
    window.project.set_metadata(sample_id, {'hydra': {'evidence_input_sha256': 'c' * 64,
                                                      'source_sample': 'A', 'hits': []}})
    window.feature_ids = {sample_id}
    window.refresh_features()
    row = window.feature_model.rows[0]

    assert row['AMR reference sets'] == 'No current AMR evidence'
    assert row['Virulence elements searched'] == 'No current AMR evidence'
    assert row['AMR evidence state'] == 'stale'
