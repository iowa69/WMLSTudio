import csv

import pytest

from wmlstudio.amr_matrix import (
    NOT_A_PHENOTYPE,
    build_determinant_matrix,
    build_mutation_matrix,
    matrix_table,
    mutation_identity,
    write_matrix,
)

DIGEST = 'a' * 64


def gene_hit(gene='blaKPC-2', **changes):
    """One acquired-determinant hit in the shape the engine reports one."""
    return dict({'gene': gene, 'element_type': 'AMR', 'element_subtype': 'AMR', 'primary': True,
                 'database': 'ncbi', 'class': 'BETA-LACTAM', 'subclass': 'CARBAPENEM',
                 'method': 'BLASTN', 'resolution': 'COMPLETE', 'identity_pct': 100.0,
                 'coverage_pct': 100.0, 'product': 'carbapenem-hydrolysing class A'}, **changes)


def point_hit(gene='gyrA', note='gyrA_S83L (L)', method='POINTX', **changes):
    """One catalogued point-mutation hit, as the protein or DNA mutation engine emits it."""
    return gene_hit(gene, element_subtype='POINT', resolution='POINT', method=method, note=note,
                    database='protein', product='Escherichia quinolone resistant GyrA',
                    **{'class': 'QUINOLONE', 'subclass': 'QUINOLONE'}, **changes)


def record(identifier, name=None, hits=(), *, organism=('Escherichia', 'coli'),
           level='dna_and_protein', catalogue='Escherichia', databases=('ncbi', 'protein'),
           point_mutations=True, provenance=True, input_sha256=DIGEST, evidence_sha256=DIGEST):
    """A project sample record carrying a linked AMR report, in the stored shape."""
    execution = {'organism': {'requested': ' '.join(organism), 'resolved': catalogue,
                              'point_mutations': point_mutations, 'point_mutation_level': level,
                              'reason': ''},
                 'virulence': {'enabled': True, 'organism_curated': True},
                 'reference_release': {'release': '2026-01-01.1'}} if provenance else {}
    return {'id': identifier, 'name': name or identifier,
            'result': {'sample_name': name or identifier, 'input_sha256': input_sha256},
            'metadata': {'organism': {'genus': organism[0], 'species': organism[1]},
                         'hydra': {'evidence_input_sha256': evidence_sha256,
                                   'source_sample': name or identifier,
                                   'databases': list(databases), 'hits': list(hits),
                                   'execution_provenance': execution}}}


def unreported(identifier, name=None):
    """An isolate with no linked AMR report at all."""
    return {'id': identifier, 'name': name or identifier,
            'result': {'sample_name': name or identifier, 'input_sha256': DIGEST},
            'metadata': {'organism': {'genus': 'Escherichia', 'species': 'coli'}}}


def cells(table, row_label):
    """The state of every isolate's cell in one row, keyed by isolate name."""
    row = next(row for row in table['rows'] if row['label'] == row_label)
    names = {column['sample_id']: column['sample_name'] for column in table['samples']}
    return {names[sample_id]: cell['state'] for sample_id, cell in row['cells'].items()}


def test_the_matrix_keeps_a_shared_pattern_visible_and_grouped_by_antimicrobial_class():
    table = build_determinant_matrix([
        record('one', 'A', [gene_hit(), gene_hit('tet(A)', **{'class': 'TETRACYCLINE',
                                                              'subclass': 'TETRACYCLINE'})]),
        record('two', 'B', [gene_hit()]),
        record('three', 'C', [gene_hit(), gene_hit('aadA2', **{'class': 'AMINOGLYCOSIDE',
                                                               'subclass': 'STREPTOMYCIN'})])])

    assert table['isolate_count'] == 3 and table['row_count'] == 3
    # Rows arrive grouped by class, and every row keeps the class it was grouped under.
    assert [(row['group'], row['label']) for row in table['rows']] == [
        ('AMINOGLYCOSIDE', 'aadA2'), ('BETA-LACTAM', 'blaKPC-2'), ('TETRACYCLINE', 'tet(A)')]
    shared = next(row for row in table['rows'] if row['label'] == 'blaKPC-2')
    assert shared['detected_in'] == 3 and shared['answered'] == 3 and shared['unknown'] == 0
    assert shared['detected_summary'] == '3 of 3 screened'
    assert cells(table, 'tet(A)') == {'A': 'detected', 'B': 'not_detected', 'C': 'not_detected'}


def test_a_partial_and_a_disrupted_match_are_not_the_same_answer_as_a_complete_one():
    table = build_determinant_matrix([
        record('one', 'A', [gene_hit()]),
        record('two', 'B', [gene_hit(resolution='PARTIAL', method='PARTIALN', coverage_pct=41.0)]),
        record('three', 'C', [gene_hit(resolution='INTERNAL_STOP')]),
        record('four', 'D', [gene_hit(resolution='', method='BLASTN')])])

    assert cells(table, 'blaKPC-2') == {'A': 'detected', 'B': 'detected_partial',
                                        'C': 'detected_disrupted', 'D': 'detected_unclassified'}
    # All four are carriage of a sequence, so all four count toward the row's numerator.
    row = table['rows'][0]
    assert row['detected_in'] == 4 and row['answered'] == 4
    for state in ('detected', 'detected_partial', 'detected_disrupted', 'detected_unclassified'):
        assert table['states'][state]['evidence'] == 'detected'


def test_no_report_and_a_report_belonging_to_another_assembly_are_unknown_not_absence():
    stale = record('three', 'C', [gene_hit()], evidence_sha256='b' * 64)
    table = build_determinant_matrix([record('one', 'A', [gene_hit()]), unreported('two', 'B'),
                                      stale])

    states = cells(table, 'blaKPC-2')
    assert states == {'A': 'detected', 'B': 'no_report', 'C': 'report_not_current'}
    for state in ('no_report', 'report_not_current', 'not_searched'):
        assert table['states'][state]['evidence'] == 'unknown'
        assert 'never an absence' in table['states'][state]['meaning']
    # The denominator is the isolates that could answer, never the cohort.
    row = table['rows'][0]
    assert row['detected_in'] == 1 and row['answered'] == 1 and row['unknown'] == 2
    assert row['detected_summary'] == '1 of 1 screened · 2 unknown'
    assert table['unknown_cells'] == 2 and table['answered_cells'] == 1
    assert 'unknown rather than negative' in table['denominator_note']


def test_a_determinant_from_a_reference_set_an_isolate_never_searched_is_not_a_negative():
    """The report ran, but not against the catalogue that names this determinant."""
    table = build_determinant_matrix([
        record('one', 'A', [gene_hit('mdtK', database='card', **{'class': 'EFFLUX'})],
               databases=('ncbi', 'protein', 'card')),
        record('two', 'B', [gene_hit()], databases=('ncbi', 'protein'))])

    assert cells(table, 'mdtK') == {'A': 'detected', 'B': 'not_searched'}
    assert cells(table, 'blaKPC-2') == {'A': 'not_detected', 'B': 'detected'}
    row = next(row for row in table['rows'] if row['label'] == 'mdtK')
    assert row['cells']['two']['state'] == 'not_searched'
    assert 'card' in row['cells']['two']['reason'] and 'ncbi' in row['cells']['two']['reason']


def test_the_matrix_says_it_is_not_a_susceptibility_result_where_it_is_read():
    table = build_determinant_matrix([record('one', 'A', [gene_hit()])])

    assert table['headline'] == NOT_A_PHENOTYPE
    assert table['limitations'][0] == NOT_A_PHENOTYPE
    assert 'no determinant detected is not a susceptible isolate' in NOT_A_PHENOTYPE
    assert any('not AMRFinderPlus' in note for note in table['limitations'])
    assert any('not transmission' in note for note in table['limitations'])
    assert 'not a susceptibility result' in table['states']['not_detected']['meaning'].casefold()


def test_point_mutations_are_a_separate_table_naming_gene_substitution_and_catalogue():
    table = build_mutation_matrix([
        record('one', 'A', [gene_hit(), point_hit()]),
        record('two', 'B', [point_hit(note='gyrA_S83L')])])

    assert table['row_count'] == 1, 'the acquired gene belongs to the other table'
    row = table['rows'][0]
    assert (row['gene'], row['substitution'], row['label']) == ('gyrA', 'S83L', 'gyrA S83L')
    assert row['catalogue'] == 'Escherichia' and row['target'] == 'the protein catalogue'
    assert row['group'] == 'QUINOLONE'
    assert row['catalogue_entry'] == 'Escherichia quinolone resistant GyrA'
    assert cells(table, 'gyrA S83L') == {'A': 'detected', 'B': 'detected'}
    assert 'observed L' in row['cells']['one']['reason']
    assert not build_determinant_matrix([record('one', 'A', [point_hit()])])['rows']


def test_an_organism_with_no_curated_catalogue_is_said_plainly_and_never_left_blank():
    table = build_mutation_matrix([
        record('one', 'A', [point_hit()]),
        record('two', 'B', [], organism=('Staphylococcus', 'aureus'), catalogue='',
               level='none'),
        record('three', 'C', [], provenance=False)])

    assert cells(table, 'gyrA S83L') == {'A': 'detected', 'B': 'no_catalogue',
                                         'C': 'catalogue_not_recorded'}
    assert table['states']['no_catalogue']['evidence'] == 'unknown'
    assert 'No other organism' in table['states']['no_catalogue']['meaning']
    assert table['without_catalogue'] == ['B', 'C']
    assert '2 of 3 isolate(s) were not screened for point mutations at all' in \
        table['catalogue_note']
    coverage = {row['sample_name']: row for row in table['catalogue_coverage']}
    assert coverage['A']['searched'] is True and coverage['A']['catalogue'] == 'Escherichia'
    assert coverage['B']['searched'] is False and coverage['B']['catalogue'] == 'none chosen'
    assert 'no point-mutation catalogue' in coverage['B']['level_words']
    assert coverage['C']['level'] == 'unrecorded'


def test_a_catalogue_that_covers_only_one_level_did_not_assess_the_other_one():
    """A protein-only release never looked in 23S rRNA; that blank is not an absence."""
    table = build_mutation_matrix([
        record('one', 'A', [point_hit('23S', '23S_A2059G', method='POINTN')],
               level='dna_and_protein'),
        record('two', 'B', [point_hit()], level='protein_only'),
        record('three', 'C', [], level='dna_only')])

    assert cells(table, '23S A2059G') == {'A': 'detected', 'B': 'target_not_assessed',
                                          'C': 'not_detected'}
    assert cells(table, 'gyrA S83L') == {'A': 'not_detected', 'B': 'detected',
                                         'C': 'target_not_assessed'}
    row = next(row for row in table['rows'] if row['gene'] == '23S')
    assert row['target'] == 'the DNA catalogue'
    assert 'no DNA-level target' in row['cells']['two']['reason']


def test_point_mutations_switched_off_is_not_an_isolate_without_mutations():
    table = build_mutation_matrix([
        record('one', 'A', [point_hit()]),
        record('two', 'B', [], point_mutations=False, level='unknown')])

    assert cells(table, 'gyrA S83L') == {'A': 'detected', 'B': 'not_searched'}
    assert 'asked not to search' in next(row for row in table['rows'])['cells']['two']['reason']
    assert table['without_catalogue'] == ['B']


def test_a_finding_that_is_not_a_substitution_is_never_printed_as_one():
    divergent = point_hit('blaTEM-1', 'divergent from susceptible reference (<98% identity)',
                          method='SUSCEPTIBLEX')
    identity = mutation_identity(divergent)
    assert identity['recorded'] is False and identity['substitution'] == ''
    assert identity['finding'] == 'divergent from susceptible reference (<98% identity)'

    table = build_mutation_matrix([record('one', 'A', [divergent])])
    row = table['rows'][0]
    assert row['substitution'] == 'no substitution recorded — see the finding'
    assert row['finding'] == identity['finding']
    assert row['label'].startswith('blaTEM-1: divergent from susceptible reference')


@pytest.mark.parametrize('note, gene, substitution', [
    ('gyrA_S83L (L)', 'gyrA', 'S83L'),
    ('pbp4_T-266A', 'pbp4', 'T-266A'),
    ('mgrB_Q30STOP', 'mgrB', 'Q30STOP'),
    ('rplD_WR65del', 'rplD', 'WR65del'),
    ('cirA_S90YfsTer15ins2', 'cirA', 'S90YfsTer15ins2'),
    ('23S_A2059G', '23S', 'A2059G'),
])
def test_the_catalogues_own_symbol_is_read_and_never_rewritten(note, gene, substitution):
    identity = mutation_identity(point_hit(gene, note))
    assert identity['gene'] == gene and identity['substitution'] == substitution
    assert identity['recorded'] is True


def catalogued_symbols(root):
    """Every mutation symbol an installed reference release actually catalogues."""
    symbols = []
    protein = root / 'prot' / 'protein' / 'AMRProt-mutation.tsv'
    if protein.is_file():
        symbols += [line.split('\t')[3] for line in protein.open(encoding='utf-8')
                    if not line.startswith('#') and len(line.split('\t')) > 3]
    for table in sorted((root / 'mutation' / 'dna').glob('*.tsv')):
        symbols += [line.split('\t')[2] for line in table.open(encoding='utf-8')
                    if not line.startswith('#') and len(line.split('\t')) > 3]
    return [symbol for symbol in symbols if symbol]


def test_every_symbol_a_real_installed_catalogue_uses_survives_this_reader():
    """Read against the reference data itself, not against a hand-written example.

    A substitution this table prints must be the catalogue's own text, character for
    character: a symbol that were rewritten, or silently dropped because it is a
    frameshift rather than a single substitution, would be a different mutation from
    the one the release catalogues.
    """
    from wmlstudio.hydra_runtime import bundled_database_root

    root = bundled_database_root()
    if root is None:
        pytest.skip('No reference release is staged in this checkout.')
    symbols = catalogued_symbols(root)
    if not symbols:
        pytest.skip('The staged release carries no point-mutation catalogue.')
    rewritten, dropped = [], []
    for symbol in symbols:
        gene = symbol.rpartition('_')[0] or symbol
        identity = mutation_identity(point_hit(gene, symbol))
        if identity['finding'] != symbol:
            dropped.append(symbol)
        elif identity['recorded'] and f"{identity['gene']}_{identity['substitution']}" != symbol:
            rewritten.append((symbol, identity['substitution']))
    assert rewritten == [] and dropped == []
    # The forms this reader must not fold together are all present in a real release.
    assert any('-' in symbol.rpartition('_')[2] for symbol in symbols), 'promoter positions'
    assert any('del' in symbol or 'ins' in symbol for symbol in symbols), 'indels'


def test_every_isolate_appears_once_and_an_empty_cohort_is_refused():
    with pytest.raises(ValueError, match='Repeated identifiers'):
        build_determinant_matrix([record('one', 'A'), record('one', 'A again')])
    with pytest.raises(ValueError, match='at least one isolate'):
        build_determinant_matrix([])
    with pytest.raises(ValueError, match='at least one isolate'):
        build_mutation_matrix([])
    with pytest.raises(ValueError, match='needs a sample id'):
        build_determinant_matrix([{'name': 'nameless'}])


def test_two_isolates_sharing_a_name_keep_separate_columns_in_the_export():
    table = build_determinant_matrix([record('one', 'A', [gene_hit()]),
                                      record('two', 'A', [])])
    headers, _ = matrix_table(table)
    assert headers[-2:] == ['A · one', 'A · two']


def test_the_exported_file_keeps_every_distinction_and_the_sentences_that_frame_it(tmp_path):
    table = build_determinant_matrix([
        record('one', 'A', [gene_hit()]),
        record('two', 'B', [gene_hit(resolution='PARTIAL', method='PARTIALN')]),
        unreported('three', 'C')])
    path = tmp_path / 'matrix.tsv'
    write_matrix(table, path)

    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert rows[0] == ['Antimicrobial class', 'Determinant', 'Reference subclass',
                       'Detected in', 'A', 'B', 'C']
    assert rows[1] == ['BETA-LACTAM', 'blaKPC-2', 'CARBAPENEM', '2 of 2 screened · 1 unknown',
                       'Detected', 'Detected · partial', 'No report']
    text = path.read_text(encoding='utf-8-sig')
    # The four distinctions survive the file, and so do the sentences that frame them.
    assert 'Detected · partial' in text and 'No report' in text and 'Not detected' in text
    assert NOT_A_PHENOTYPE in text
    assert 'What each cell state means' in text
    assert 'Unknown, never an absence.' in text


def test_a_csv_export_is_still_a_matrix_and_guards_spreadsheet_formula_cells(tmp_path):
    table = build_mutation_matrix([record('one', '=cmd', [point_hit()])])
    path = tmp_path / 'mutations.csv'
    write_matrix(table, path)

    with path.open(encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle))
    # A sample name is a column heading here, so the heading row is guarded like a cell.
    assert rows[0] == ['Antimicrobial class', 'Gene', 'Substitution', 'Organism catalogue',
                       'Read from', 'Detected in', "'=cmd"]
    assert rows[1][:3] == ['QUINOLONE', 'gyrA', 'S83L']
    assert rows[1][-1] == 'Detected'
