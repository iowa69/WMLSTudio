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


def screened(window, name, *, hits=(), databases=('ncbi', 'protein'), catalogue='Escherichia',
             level='dna_and_protein', linked=True, current=True):
    """An isolate carrying a linked AMR report that records what its run could search."""
    digest = 'a' * 64
    evidence = {'evidence_input_sha256': digest if current else 'b' * 64, 'source_sample': name,
                'databases': list(databases), 'hits': list(hits),
                'execution_provenance': {
                    'organism': {'requested': 'Escherichia coli', 'resolved': catalogue,
                                 'point_mutations': True, 'point_mutation_level': level,
                                 'reason': ''},
                    'virulence': {'enabled': True, 'organism_curated': True},
                    'reference_release': {'release': '2026-01-01.1'}}}
    metadata = {'organism': {'genus': 'Escherichia', 'species': 'coli'}}
    if linked:
        metadata['hydra'] = evidence
    return window.project.add_profile(name, {
        'sample_name': name, 'st': '20', 'input_sha256': digest, 'scheme': 'MLST',
        'scheme_digest': 'reference', 'status': 'profile_imported', 'alleles': {'a': '1'}}, metadata)


def acquired(gene='blaKPC-2', drug_class='BETA-LACTAM', **changes):
    return dict({'gene': gene, 'element_type': 'AMR', 'element_subtype': 'AMR', 'primary': True,
                 'database': 'ncbi', 'class': drug_class, 'subclass': drug_class,
                 'method': 'BLASTN', 'resolution': 'COMPLETE', 'identity_pct': 100.0,
                 'coverage_pct': 100.0}, **changes)


def mutation(gene='gyrA', note='gyrA_S83L (L)', **changes):
    return acquired(gene, 'QUINOLONE', element_subtype='POINT', resolution='POINT',
                    method='POINTX', database='protein', note=note, **changes)


def test_the_cohort_matrix_shows_determinants_by_class_with_each_cell_state_in_words(window):
    """Colour is never the only carrier: the cell text alone has to survive printing."""
    carrier = screened(window, 'A', hits=[acquired(), acquired('tet(A)', 'TETRACYCLINE')])
    partial = screened(window, 'B', hits=[acquired(resolution='PARTIAL', method='PARTIALN')])
    unscreened = screened(window, 'C', linked=False)
    window.feature_ids = {carrier, partial, unscreened}
    window.refresh_features()

    tabs = window.evidence_tabs
    assert any(tabs.tabText(index) == 'AMR determinants by class' for index in range(tabs.count()))
    table = window.determinant_table
    headers = [table.horizontalHeaderItem(column).text() for column in range(table.columnCount())]
    assert headers == ['Antimicrobial class', 'Determinant', 'Reference subclass', 'Detected in',
                       'A', 'B', 'C']
    assert table_rows(table) == [
        ['BETA-LACTAM', 'blaKPC-2', 'BETA-LACTAM', '2 of 2 screened · 1 unknown',
         'Detected', 'Detected · partial', 'No report'],
        ['TETRACYCLINE', 'tet(A)', 'TETRACYCLINE', '1 of 2 screened · 1 unknown',
         'Detected', 'Not detected', 'No report']]
    # The distinction is carried by the words in the cell, and the tooltip says why.
    assert 'never an absence' in table.item(0, 6).toolTip()
    assert 'No linked HYDRA report' in table.item(0, 6).toolTip()
    assert 'Not a susceptibility result' in window.determinant_scope.text()
    assert 'Not detected: 1' in window.determinant_legend.text()
    assert 'No report: 2' in window.determinant_legend.text()


def test_filtering_and_ungrouping_the_matrix_changes_no_count_on_the_page(window):
    window.feature_ids = {screened(window, 'A', hits=[acquired(), acquired('tet(A)',
                                                                          'TETRACYCLINE')])}
    window.refresh_features()
    before = window.determinant_scope.text()

    window.determinant_grouping.setChecked(False)
    assert [row[1] for row in table_rows(window.determinant_table)] == ['blaKPC-2', 'tet(A)']
    assert [row[0] for row in table_rows(window.determinant_table)] == ['BETA-LACTAM',
                                                                       'TETRACYCLINE']
    window.determinant_filter.setText('tetra')
    assert [row[1] for row in table_rows(window.determinant_table)] == ['tet(A)']
    # Only the "showing" count moves: the cohort, the rows and the cell denominator
    # beside them are what was screened, not what is on screen.
    after = window.determinant_scope.text()
    assert '2 determinant(s) across 1 isolate(s); showing 1.' in after
    assert before.split('showing')[0] == after.split('showing')[0]
    assert '2 row(s) × 1 isolate(s)' in after and '0 of 2 cells are unknown' in after


def test_the_point_mutation_tab_names_the_gene_the_substitution_and_the_catalogue(window):
    found = screened(window, 'A', hits=[acquired(), mutation()])
    without = screened(window, 'B', hits=[acquired()], catalogue='', level='none')
    window.feature_ids = {found, without}
    window.refresh_features()

    tabs = window.evidence_tabs
    assert any(tabs.tabText(index) == 'Resistance point mutations' for index in range(tabs.count()))
    table = window.mutation_table
    headers = [table.horizontalHeaderItem(column).text() for column in range(table.columnCount())]
    assert headers == ['Antimicrobial class', 'Gene', 'Substitution', 'Organism catalogue',
                       'Read from', 'Detected in', 'A', 'B']
    assert table_rows(table) == [['QUINOLONE', 'gyrA', 'S83L', 'Escherichia',
                                  'the protein catalogue', '1 of 1 screened · 1 unknown',
                                  'Detected', 'No catalogue for this organism']]
    # The acquired gene is not repeated here: the two are different evidence.
    assert 'blaKPC-2' not in [row[1] for row in table_rows(table)]
    assert 'No other organism' in table.item(0, 7).toolTip()


def test_an_isolate_with_no_curated_catalogue_is_listed_rather_than_left_looking_clean(window):
    window.feature_ids = {screened(window, 'A', hits=[mutation()]),
                          screened(window, 'B', catalogue='', level='none'),
                          screened(window, 'C', linked=False)}
    window.refresh_features()

    assert '2 of 3 isolate(s) were not screened for point mutations at all' in \
        window.mutation_catalogue_note.text()
    rows = table_rows(window.mutation_catalogue_table)
    assert [row[0] for row in rows] == ['A', 'B', 'C']
    assert rows[0][2] == 'Escherichia' and 'protein and DNA' in rows[0][3]
    assert rows[1][2] == 'none chosen' and 'no point-mutation catalogue' in rows[1][3]
    assert rows[2][4] == 'missing' and 'No linked HYDRA report' in rows[2][5]


def test_an_empty_evidence_cohort_says_unknown_rather_than_showing_an_empty_matrix(window):
    window.feature_ids = set()
    window.refresh_features()

    assert window.determinant_table.rowCount() == 0
    assert 'unknown evidence, not an isolate carrying nothing' in window.determinant_scope.text()
    assert 'no point-mutation catalogue was read' in window.mutation_scope.text()
    assert 'Choose evidence isolates' in window.mutation_catalogue_note.text()
    with pytest.raises(ValueError, match='would read as a cohort with nothing in it'):
        window.write_amr_matrix('unused.tsv', 'determinants')


def test_both_cohort_tables_export_with_their_distinctions_intact(window, tmp_path):
    window.feature_ids = {screened(window, 'A', hits=[acquired(), mutation()]),
                          screened(window, 'B', linked=False)}
    window.refresh_features()

    path = tmp_path / 'determinants.tsv'
    window.write_amr_matrix(path, 'determinants')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert rows[0][:2] == ['Antimicrobial class', 'Determinant'] and rows[0][-2:] == ['A', 'B']
    assert rows[1][:2] == ['BETA-LACTAM', 'blaKPC-2'] and rows[1][-2:] == ['Detected', 'No report']
    assert 'Not a susceptibility result.' in path.read_text(encoding='utf-8-sig')

    path = tmp_path / 'mutations.tsv'
    window.write_amr_matrix(path, 'point_mutations')
    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert rows[1][:3] == ['QUINOLONE', 'gyrA', 'S83L']
    assert rows[1][-2:] == ['Detected', 'No report']
    assert 'no other organism' in path.read_text(encoding='utf-8-sig').casefold()


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
