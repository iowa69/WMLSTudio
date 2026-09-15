import pytest

from wmlstudio.plasmid_evidence import (
    BACKBONE_MIN_BP,
    MOB_SUITE_GAP,
    SUPPORT_STATES,
    backbone_depth,
    cohort_replicon_cooccurrence,
    contig_plasmid_evidence,
    contig_records,
    marker_location,
    parse_contig_header,
)
from wmlstudio.sequence import AnalysisCancelled


def replicon(gene='IncFIB', contig='Contig_2', start=101, end=400):
    """One accepted replicon marker, in the shape the AMR/plasmid assay reports."""
    return {'gene': gene, 'element_type': 'PLASMID', 'accession': f'{gene}_1',
            'sequence': contig, 'start': start, 'end': end, 'resolution': 'COMPLETE'}


def determinant(gene='blaKPC-2', contig='Contig_2', start=1000, end=1800, kind='AMR'):
    return {'gene': gene, 'element_type': kind, 'class': 'BETA-LACTAM', 'subclass': 'CARBAPENEM',
            'sequence': contig, 'start': start, 'end': end, 'resolution': 'COMPLETE'}, kind


# Header shapes taken from real assembler output. The SKESA ones are the exact
# conventions seen in this project's own validated Windows assemblies.
@pytest.mark.parametrize('header,convention,coverage,closure', [
    ('Contig_2_201.434_Circ', 'skesa', 201.434, 'declared_circular'),
    ('Contig_1_32.7795', 'skesa', 32.7795, 'not_declared'),
    ('NODE_3_length_4372_cov_201.434', 'spades', 201.434, 'unknown'),
    ('7 length=4372 depth=4.11x circular=true', 'key_value_header', 4.11, 'declared_circular'),
    ('7 length=4372 depth=1.00x circular=false', 'key_value_header', 1.0, 'not_declared'),
    ('tig00000001', 'unrecognized', None, 'unknown'),
])
def test_assembler_claims_are_read_only_from_structured_header_tokens(header, convention, coverage, closure):
    parsed = parse_contig_header(header)
    assert parsed['assembler_convention'] == convention
    assert parsed['declared_coverage'] == coverage
    assert parsed['closure'] == closure


def test_a_header_that_merely_says_plasmid_is_never_read_as_a_plasmid_claim():
    """Depositor prose is not evidence; only machine-readable tokens are parsed."""
    parsed = parse_contig_header('NZ_CP012345.1 Klebsiella pneumoniae plasmid pKPC-2, complete sequence')
    assert parsed['closure'] == 'unknown' and parsed['declared_coverage'] is None
    assert parsed['assembler_convention'] == 'unrecognized'
    bracketed = parse_contig_header('NZ_CP012345.1 plasmid pKPC-2 [topology=circular]')
    assert bracketed['closure'] == 'declared_circular'


def test_a_declared_length_never_overrides_the_counted_sequence_length():
    records = contig_records({'NODE_1': 4000}, {'NODE_1': 'NODE_1_length_9999_cov_50.0'})
    assert records[0]['length_bp'] == 4000
    assert records[0]['declared_length_disagrees'] is True


def test_unknown_closure_is_kept_distinct_from_a_declared_absence_of_circularity():
    """SKESA marks its circular contigs; SPAdes says nothing. Those differ."""
    skesa = parse_contig_header('Contig_9_52.9944')['closure']
    spades = parse_contig_header('NODE_9_length_500_cov_52.9944')['closure']
    assert skesa == 'not_declared' and spades == 'unknown'


def test_backbone_depth_carries_its_denominator_and_refuses_to_guess():
    lengths = {'Contig_1_50': 2_000_000, 'Contig_2_400_Circ': 2000, 'Contig_3_9': 400}
    headers = {'Contig_1_50': 'Contig_1_50', 'Contig_2_400_Circ': 'Contig_2_400_Circ',
               'Contig_3_9': 'Contig_3_9'}
    basis = backbone_depth(contig_records(lengths, headers))
    assert basis['status'] == 'available' and basis['median_declared_coverage'] == 50.0
    # The sub-floor contig is excluded from the estimate but stays in the total.
    assert basis['contigs'] == 2 and basis['bases'] == 2_002_000
    assert basis['assembly_bases'] == 2_002_400
    assert 0.99 < basis['assembly_fraction'] < 1.0
    unknown = backbone_depth(contig_records({'tig1': 5000}, {'tig1': 'tig1'}))
    assert unknown['status'] == 'unavailable' and unknown['median_declared_coverage'] is None
    assert unknown['assembly_fraction'] is None


def test_a_backbone_of_zero_coverage_never_becomes_a_depth_ratio():
    basis = backbone_depth(contig_records({'Contig_1_0': 5000}, {'Contig_1_0': 'Contig_1_0'}))
    assert basis['status'] == 'unavailable' and 'zero' in basis['reason']


@pytest.mark.parametrize('start,end', [(0, 100), (1, 100000), (True, 100), (None, 100)])
def test_a_marker_outside_its_contig_is_not_placed(start, end):
    assert marker_location({'sequence': 'c', 'start': start, 'end': end}, {'c': 2000}) is None
    assert marker_location({'sequence': 'absent', 'start': 1, 'end': 10}, {'c': 2000}) is None


def small_circular_assembly():
    """A 2.8 Mb backbone plus one 4.4 kb contig the assembler called circular.

    The shape and the numbers follow this project's real validated SKESA output
    for DRR428002, where the circular contig sat at 4.1x the backbone median.
    """
    lengths = {'Contig_1_48.6581': 2_800_000, 'Contig_2_201.434_Circ': 4372,
               'Contig_4_34.2467': 50_000}
    return lengths, {name: name for name in lengths}


def test_closure_and_depth_are_two_independent_signals_and_are_named_separately():
    lengths, headers = small_circular_assembly()
    result = contig_plasmid_evidence(lengths, [replicon(contig='Contig_2_201.434_Circ')],
                                     [determinant(contig='Contig_2_201.434_Circ')],
                                     contig_headers=headers)
    row = result['contigs'][0]
    assert row['support_state'] == 'replicon_marker_plus_closure_and_depth'
    assert sorted(row['signals']) == ['closure_claim', 'depth_departure']
    assert 4.0 < row['depth_ratio'] < 4.2 and row['length_bp'] == 4372
    assert row['closure'] == 'declared_circular'
    # Stripping the circularity claim must drop exactly one signal, not the row.
    plain = dict(headers, **{'Contig_2_201.434_Circ': 'Contig_2_201.434'})
    only_depth = contig_plasmid_evidence(lengths, [replicon(contig='Contig_2_201.434_Circ')],
                                         [determinant(contig='Contig_2_201.434_Circ')],
                                         contig_headers=plain)['contigs'][0]
    assert only_depth['support_state'] == 'replicon_marker_plus_depth_departure'
    assert only_depth['signals'] == ['depth_departure']


def test_without_any_assembler_metadata_support_is_not_assessable_rather_than_absent():
    lengths = {'tig1': 2_800_000, 'tig2': 4372}
    result = contig_plasmid_evidence(lengths, [replicon(contig='tig2')], [determinant(contig='tig2')],
                                     contig_headers={name: name for name in lengths})
    row = result['contigs'][0]
    assert row['support_state'] == 'not_assessable' and row['signals'] == []
    assert row['depth_ratio'] is None
    assert result['depth_basis']['status'] == 'unavailable'
    assert row['support_state'] in SUPPORT_STATES


def test_a_determinant_with_no_replicon_on_its_contig_is_reported_not_omitted():
    """An omitted row would read as absence of the determinant, which it is not."""
    lengths, headers = small_circular_assembly()
    result = contig_plasmid_evidence(lengths, [replicon(contig='Contig_2_201.434_Circ')],
                                     [determinant(gene='blaSHV-1', contig='Contig_1_48.6581',
                                                  start=10, end=800)],
                                     contig_headers=headers)
    placement = result['determinant_placement']
    assert [row['gene'] for row in placement] == ['blaSHV-1']
    assert placement[0]['placement'] == 'no_replicon_on_this_contig'
    assert 'does not place the determinant on the chromosome' in placement[0]['interpretation']
    chromosomal = [row for row in result['contigs'] if row['contig'] == 'Contig_1_48.6581'][0]
    assert chromosomal['support_state'] == 'no_replicon_marker'


def test_a_determinant_that_cannot_be_placed_is_unplaced_not_separated():
    lengths, headers = small_circular_assembly()
    nowhere, kind = determinant(contig='missing_contig')
    result = contig_plasmid_evidence(lengths, [replicon(contig='Contig_2_201.434_Circ')],
                                     [(nowhere, kind)], contig_headers=headers)
    assert result['unplaced_determinants'] == 1
    assert result['determinant_placement'][0]['placement'] == 'unplaced'
    assert result['determinant_placement'][0]['contig'] is None


def test_every_per_isolate_result_states_what_mob_suite_would_have_added():
    lengths, headers = small_circular_assembly()
    result = contig_plasmid_evidence(lengths, [], [], contig_headers=headers)
    assert result['mob_suite_gap'] == MOB_SUITE_GAP
    joined = ' '.join(MOB_SUITE_GAP)
    for missing in ('relaxase', 'MPF', 'oriT', 'reconstruction', 'cluster code'):
        assert missing in joined
    assert 'not equivalent' in joined


def cohort_entry(sample_id, replicons=(), links=(), placements=(), status='completed', reason=''):
    return {'sample_id': sample_id, 'sample_name': sample_id.upper(),
            'plasmid_hypotheses': {'status': status, 'reason': reason,
                                   'replicons': [{'gene': gene} for gene in replicons],
                                   'contig_associations': [{'replicon': rep, 'marker': marker,
                                                            'marker_type': 'AMR'} for rep, marker in links],
                                   'determinant_placement': list(placements)}}


def test_cohort_co_occurrence_separates_co_located_from_merely_both_present():
    cohort = [
        cohort_entry('a', ['IncFIB'], [('IncFIB', 'blaKPC-2')],
                     [{'gene': 'blaKPC-2', 'marker_type': 'AMR', 'placement': 'co_located_with_replicon'}]),
        cohort_entry('b', ['IncFIB'], [('IncFIB', 'blaKPC-2')],
                     [{'gene': 'blaKPC-2', 'marker_type': 'AMR', 'placement': 'co_located_with_replicon'}]),
        cohort_entry('c', ['IncFIB'], [],
                     [{'gene': 'blaKPC-2', 'marker_type': 'AMR', 'placement': 'no_replicon_on_this_contig'}]),
    ]
    result = cohort_replicon_cooccurrence(cohort)
    assert result['denominator'] == 3 and result['isolates_assayed'] == ['a', 'b', 'c']
    pair = result['co_occurrence'][0]
    assert pair['replicon'] == 'IncFIB' and pair['marker'] == 'blaKPC-2'
    assert pair['co_located_isolates'] == ['a', 'b'] and pair['co_located_count'] == 2
    assert pair['present_but_not_co_located'] == ['c'] and pair['not_co_located_count'] == 1
    assert pair['denominator'] == 3 and pair['status'] == 'hypothesis'
    assert 'not transmission' in pair['interpretation']
    assert result['replicons'][0] == {'replicon': 'IncFIB', 'isolates': ['a', 'b', 'c'],
                                      'isolate_count': 3, 'denominator': 3}


def test_an_isolate_whose_plasmid_assay_did_not_run_is_named_not_counted_as_negative():
    cohort = [cohort_entry('a', ['IncFIB'], [('IncFIB', 'blaKPC-2')],
                           [{'gene': 'blaKPC-2', 'marker_type': 'AMR',
                             'placement': 'co_located_with_replicon'}]),
              cohort_entry('b', status='not_run', reason='No plasmid-reference assay was run.'),
              {'sample_id': 'c'}]
    result = cohort_replicon_cooccurrence(cohort)
    assert result['denominator'] == 1 and result['isolates_assayed'] == ['a']
    assert [row['sample_id'] for row in result['isolates_excluded']] == ['b', 'c']
    assert result['isolates_excluded'][0]['reason'] == 'No plasmid-reference assay was run.'
    assert 'missing' in result['isolates_excluded'][1]['reason']
    assert result['co_occurrence'][0]['denominator'] == 1


def test_an_unplaced_determinant_is_never_counted_as_evidence_of_separation():
    cohort = [cohort_entry('a', ['IncFIB'], [],
                           [{'gene': 'blaKPC-2', 'marker_type': 'AMR', 'placement': 'unplaced'}])]
    result = cohort_replicon_cooccurrence(cohort)
    assert result['co_occurrence'] == []
    assert result['unplaced_determinants'] == [{'sample_id': 'a', 'gene': 'blaKPC-2', 'marker_type': 'AMR'}]


def test_the_cohort_view_produces_no_distance_and_no_threshold():
    """Plasmid co-occurrence must never acquire a scale MLST/cgMLST/SNP could share."""
    cohort = [cohort_entry('a', ['IncFIB'], [('IncFIB', 'blaKPC-2')],
                           [{'gene': 'blaKPC-2', 'marker_type': 'AMR',
                             'placement': 'co_located_with_replicon'}])]
    result = cohort_replicon_cooccurrence(cohort)
    forbidden = {'distance', 'threshold', 'similarity', 'cluster', 'score', 'snp', 'allele_differences'}
    flat = [result, *result['co_occurrence'], *result['replicons']]
    assert not any(key in forbidden for block in flat for key in block)


def test_a_cohort_entry_without_an_identifier_is_refused_rather_than_merged():
    with pytest.raises(ValueError, match='sample_id'):
        cohort_replicon_cooccurrence([{'plasmid_hypotheses': {'status': 'completed'}}])


def test_cohort_work_is_cancellable():
    with pytest.raises(AnalysisCancelled):
        cohort_replicon_cooccurrence([cohort_entry('a')], cancelled=lambda: True)


def test_the_backbone_floor_is_a_stated_constant_not_a_hidden_rule():
    assert BACKBONE_MIN_BP == 1000
    lengths = {'Contig_1_50': 5000, 'Contig_2_900_Circ': BACKBONE_MIN_BP - 1}
    basis = backbone_depth(contig_records(lengths, {name: name for name in lengths}))
    assert basis['contigs'] == 1 and basis['min_contig_bp'] == BACKBONE_MIN_BP
    assert f'{BACKBONE_MIN_BP} bp' in basis['reason']
