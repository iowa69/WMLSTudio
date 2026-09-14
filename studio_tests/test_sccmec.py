import json
import os
import random
from pathlib import Path

import pytest

from wmlstudio.characterization_refs import _split_fasta, reference_digest, sccmec_rules
from wmlstudio.marker_panel import union_coverage
from wmlstudio.sccmec_evidence import interpret_sccmec, sccmec_html, summarize_sccmec, type_sccmec
from wmlstudio.sequence import file_sha256

TARGETS = ('ccrA1', 'ccrB1', 'ccrA2', 'ccrB2', 'ccrC1', 'mecA', 'mecR1', 'mecI',
           'IS431', 'IS431_1', 'IS431_2', 'IS1272', 'IS12960D')
RULES = {
    'targets': list(TARGETS),
    'aliases': [{'name': 'ccr Type 1', 'targets': ['ccrA1', 'ccrB1']},
                {'name': 'ccr Type 2', 'targets': ['ccrA2', 'ccrB2']},
                {'name': 'ccr Type 5', 'targets': ['ccrC1']},
                {'name': 'mec Class A', 'targets': ['IS431', 'mecA', 'mecR1', 'mecI']},
                {'name': 'mec Class B', 'targets': ['IS431', 'mecA', 'mecR1', 'IS1272']},
                {'name': 'mec Class C', 'targets': ['IS431_1', 'mecA', 'mecR1', 'IS431_2']}],
    'types': [{'name': 'II', 'targets': ['ccr Type 2', 'mec Class A'], 'excludes': []},
              {'name': 'IV', 'targets': ['ccr Type 2', 'mec Class B'], 'excludes': []},
              {'name': 'V', 'targets': ['ccr Type 5', 'mec Class C'], 'excludes': ['mecI', 'IS12960D']}],
}
TYPE_IV = ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1', 'IS1272')
TYPE_V = ('ccrC1', 'IS431_1', 'IS431_2', 'mecA', 'mecR1')


def states(detected=(), ambiguous=()):
    calls = {target: 'not_detected' for target in TARGETS}
    calls.update({target: 'detected' for target in detected})
    calls.update({target: 'ambiguous' for target in ambiguous})
    return calls


def placements(detected, contig='contig_1', overrides=None):
    places = {target: [(contig, 100, 900, '+')] for target in detected}
    for target, other in (overrides or {}).items():
        places[target] = [(other, 100, 900, '+')]
    return places


def test_type_requires_a_single_candidate_with_both_complexes_resolved():
    result = interpret_sccmec(states(TYPE_IV), placements(TYPE_IV), RULES)
    assert result['status'] == 'completed' and result['type'] == 'IV'
    assert result['candidate_types'] == ['IV'] and result['official_type'] is None
    assert result['contig_fragmented'] is False and result['ccr_mec_same_contig'] is True
    assert {row['name']: row['state'] for row in result['ccr_complexes']}['ccr Type 2'] == 'present'
    assert {row['name']: row['state'] for row in result['mec_classes']}['mec Class B'] == 'present'
    assert summarize_sccmec(result) == 'IV (candidate) - ccr2 + mec B'


def test_two_candidate_types_are_ambiguous_not_a_call():
    detected = (*TYPE_IV, 'mecI')
    result = interpret_sccmec(states(detected), placements(detected), RULES)
    assert sorted(result['candidate_types']) == ['II', 'IV']
    assert result['status'] == 'ambiguous' and result['type'] is None
    assert 'more than one SCCmec element' in result['reason']
    assert summarize_sccmec(result) == 'ambiguous - 2 candidate types'


def test_excludes_rule_prevents_a_type_when_a_forbidden_target_is_present():
    allowed = interpret_sccmec(states(TYPE_V), placements(TYPE_V), RULES)
    assert allowed['type'] == 'V'
    detected = (*TYPE_V, 'IS12960D')
    blocked = interpret_sccmec(states(detected), placements(detected), RULES)
    assert blocked['candidate_types'] == [] and blocked['type'] is None
    assert blocked['status'] == 'ambiguous' and 'mecA is present' in blocked['reason']


def test_required_targets_on_different_contigs_withhold_the_type_and_flag_fragmentation():
    result = interpret_sccmec(states(TYPE_IV), placements(TYPE_IV, overrides={'mecA': 'contig_9'}), RULES)
    assert result['candidate_types'] == ['IV'] and result['type'] is None
    assert result['status'] == 'ambiguous' and result['contig_fragmented'] is True
    assert result['contigs'] == ['contig_1', 'contig_9']
    assert result['ccr_mec_same_contig'] is True  # ccr and part of the mec class still share a contig
    assert 'spans contig breaks' in result['reason']
    assert summarize_sccmec(result) == 'ambiguous - IV withheld'


def test_an_ambiguous_target_leaves_its_complex_unresolved_and_withholds_the_type():
    result = interpret_sccmec(states(TYPE_IV, ambiguous=('ccrA1',)), placements(TYPE_IV), RULES)
    assert result['candidate_types'] == ['IV'] and result['type'] is None
    assert result['status'] == 'ambiguous' and 'ccr Type 1' in result['reason']
    assert {row['name']: row['state'] for row in result['ccr_complexes']}['ccr Type 1'] == 'unresolved'


def test_partial_complexes_are_reported_as_partial_not_absent():
    detected = ('ccrA2', 'IS431', 'mecA')
    result = interpret_sccmec(states(detected), placements(detected), RULES)
    complexes = {row['name']: row['state'] for row in (*result['ccr_complexes'], *result['mec_classes'])}
    assert complexes['ccr Type 2'] == 'partial' and complexes['mec Class B'] == 'partial'
    assert complexes['ccr Type 1'] == 'absent'
    assert result['type'] is None and result['status'] == 'ambiguous'


def test_meca_present_with_no_matching_definition_is_ambiguous_not_negative():
    result = interpret_sccmec(states(('mecA',)), placements(('mecA',)), RULES)
    assert result['status'] == 'ambiguous' and result['type'] is None and result['candidate_types'] == []
    assert 'novel, composite or truncated element' in result['reason']
    assert summarize_sccmec(result) == 'ambiguous - untypeable by this panel'


def test_negative_is_withheld_when_assembly_qc_gates_fail():
    adequate = interpret_sccmec(states(), {}, RULES, qc_adequate=True)
    assert adequate['status'] == 'not_detected' and 'mecC is not in this panel' in adequate['reason']
    withheld = interpret_sccmec(states(), {}, RULES, qc_adequate=False)
    assert withheld['status'] == 'ambiguous' and withheld['qc_adequate'] is False
    assert 'withheld rather than reported as absence' in withheld['reason']


def test_mecc_is_reported_as_not_assayed_and_never_as_absent():
    for calls in (states(), states(TYPE_IV), states(('mecA',))):
        result = interpret_sccmec(calls, {}, RULES)
        assert result['mecC'] == 'not_assayed'
        assert 'not_detected' != result['mecC'] and result['mecC'] != 'absent'
    assert any('mecC is not in this reference panel' in text for text in sccmec_html({'status': 'not_detected'}).split('<p'))


def test_region_support_never_promotes_a_subtype_to_the_type_call():
    support = {'status': 'provisional_region_match', 'subtype': None,
               'best': {'subtype': 'IVa', 'accession': 'AB063172', 'coverage_pct': 96.5, 'identity_pct': 99.1,
                        'contig_count': 3, 'aligned_bp': 23_450, 'reference_bp': 24_304}}
    detected = ('mecA',)
    result = interpret_sccmec(states(detected), placements(detected), RULES, region_support=support)
    assert result['type'] is None and result['candidate_types'] == []
    assert result['region_support']['best']['subtype'] == 'IVa'
    assert result['region_support']['subtype'] is None
    assert 'IVa' not in summarize_sccmec(result)
    assert 'Never promoted to a type or subtype call' in sccmec_html(result)


def test_official_type_is_never_populated():
    for calls, places in ((states(TYPE_IV), placements(TYPE_IV)), (states(), {})):
        assert interpret_sccmec(calls, places, RULES)['official_type'] is None


def test_rendered_evidence_escapes_values_and_repeats_every_limitation():
    import html as html_module

    from wmlstudio.sccmec_evidence import LIMITATIONS
    result = interpret_sccmec(states(TYPE_IV), placements(TYPE_IV, contig='<script>x</script>'), RULES)
    rendered = sccmec_html(dict(result, limitations=list(LIMITATIONS)))
    assert '<script>' not in rendered and '&lt;script&gt;' in rendered
    unescaped = html_module.unescape(rendered)
    for limitation in LIMITATIONS:
        assert limitation in unescaped


def test_union_coverage_sums_intervals_across_contigs_without_double_counting():
    hits = [{'reference_start': 1, 'reference_end': 500, 'contig': 'c1', 'identity_pct': 100, 'bitscore': 900},
            {'reference_start': 400, 'reference_end': 700, 'contig': 'c1', 'identity_pct': 100, 'bitscore': 500},
            {'reference_start': 701, 'reference_end': 1000, 'contig': 'c2', 'identity_pct': 90, 'bitscore': 400}]
    coverage = union_coverage(hits, 1000)
    assert coverage['aligned_bp'] == 1000 and coverage['coverage_pct'] == 100
    assert coverage['contigs'] == ['c1', 'c2'] and coverage['contig_count'] == 2
    assert coverage['fragmented'] is True
    assert round(coverage['identity_pct'], 3) == round((700 * 100 + 300 * 90) / 1000, 3)


def test_union_coverage_ignores_a_fully_contained_duplicate_alignment():
    hits = [{'reference_start': 1, 'reference_end': 800, 'contig': 'c1', 'identity_pct': 99, 'bitscore': 900},
            {'reference_start': 100, 'reference_end': 400, 'contig': 'c7', 'identity_pct': 80, 'bitscore': 100}]
    coverage = union_coverage(hits, 1600)
    assert coverage['aligned_bp'] == 800 and coverage['coverage_pct'] == 50
    assert coverage['contigs'] == ['c1'] and coverage['fragmented'] is False
    assert coverage['identity_pct'] == 99


def test_union_coverage_rejects_coordinates_outside_the_reference():
    with pytest.raises(ValueError, match='outside the reference'):
        union_coverage([{'reference_start': 1, 'reference_end': 200, 'contig': 'c', 'identity_pct': 99,
                         'bitscore': 1}], 100)
    with pytest.raises(ValueError, match='positive reference length'):
        union_coverage([], 0)
    assert union_coverage([], 500) == {'aligned_bp': 0, 'coverage_pct': 0.0, 'identity_pct': 0.0,
                                       'contigs': [], 'contig_count': 0, 'fragmented': False}


def rule_document(**changes):
    document = {'metadata': {'version': '1.2.0'}, 'engine': {'params': {'min_pident': 90, 'min_coverage': 80}},
                'targets': list(TARGETS), 'aliases': json.loads(json.dumps(RULES['aliases'])),
                'types': json.loads(json.dumps(RULES['types']))}
    document.update(changes)
    return document


def rule_table(rows):
    lines = ['file\tformat\taccession\ttype\ttarget\tstrand\tstart\tstop']
    for sccmec_type, targets in rows.items():
        for target in targets:
            lines.append(f'type-{sccmec_type}.gb\tgenbank\tX1\t{sccmec_type}\t{target}\t+\t1\t100')
    return '\n'.join(lines) + '\n'


UPSTREAM_TABLE = rule_table({'II': ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1', 'mecI'),
                             'IV': TYPE_IV,
                             'V': ('ccrC1', 'IS431_1', 'IS431_2', 'mecA', 'mecR1')})


def test_sccmec_yaml_rules_cross_check_against_the_upstream_tsv():
    rules = sccmec_rules(rule_document(), UPSTREAM_TABLE, TARGETS)
    assert rules['schema_version'] == '1.2.0'
    assert [definition['name'] for definition in rules['types']] == ['II', 'IV', 'V']
    assert rules['aliases'][0] == {'name': 'ccr Type 1', 'targets': ['ccrA1', 'ccrB1']}
    disagreeing = rule_table({'II': ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1', 'mecI'),
                              'IV': ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1'),
                              'V': ('ccrC1', 'IS431_1', 'IS431_2', 'mecA', 'mecR1')})
    with pytest.raises(ValueError, match='does not match the upstream target table'):
        sccmec_rules(rule_document(), disagreeing, TARGETS)


def test_a_yaml_alias_naming_an_unknown_target_aborts_staging():
    document = rule_document(aliases=[{'name': 'ccr Type 2', 'targets': ['ccrA2', 'ccrB99']},
                                      *RULES['aliases'][3:]])
    with pytest.raises(ValueError, match='names targets absent from the panel'):
        sccmec_rules(document, UPSTREAM_TABLE, TARGETS)
    unknown_type = rule_document(types=[{'name': 'IV', 'targets': ['ccr Type 42', 'mec Class B']}])
    with pytest.raises(ValueError, match='neither a target nor an alias'):
        sccmec_rules(unknown_type, UPSTREAM_TABLE, TARGETS)


def test_rules_are_refused_when_the_panel_or_the_engine_floors_disagree():
    with pytest.raises(ValueError, match='exactly the staged target set'):
        sccmec_rules(rule_document(), UPSTREAM_TABLE, [*TARGETS, 'mecC'])
    with pytest.raises(ValueError, match='engine parameters differ'):
        sccmec_rules(rule_document(engine={'params': {'min_pident': 80, 'min_coverage': 80}}),
                     UPSTREAM_TABLE, TARGETS)
    without_rows = rule_table({'II': ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1', 'mecI')})
    with pytest.raises(ValueError, match='no rows in the upstream target table'):
        sccmec_rules(rule_document(), without_rows, TARGETS)


def test_an_excluded_target_recorded_as_present_upstream_aborts_staging():
    contradictory = rule_table({'II': ('ccrA2', 'ccrB2', 'IS431', 'mecA', 'mecR1', 'mecI'),
                                'IV': TYPE_IV,
                                'V': ('ccrC1', 'IS431_1', 'IS431_2', 'mecA', 'mecR1', 'mecI')})
    with pytest.raises(ValueError, match='does not match the upstream target table'):
        sccmec_rules(rule_document(), contradictory, TARGETS)


def test_target_fasta_split_preserves_every_sequence_byte_and_deduplicates_identifiers(tmp_path):
    from wmlstudio.characterization_refs import _sccmec_split_name
    from wmlstudio.sequence import SequenceReader
    stage = tmp_path / 'stage'
    (stage / 'sources/sccmec').mkdir(parents=True)
    records = [('ccrA1 AB033763.2|I', 'ACGTACGTAA'), ('ccrA1 AB505628.1|IX', 'TTGGCCAATT'),
               ('mecA D86934.2|II', 'ATGAAATAAG')]
    (stage / 'sources/sccmec/sccmec-targets.fasta').write_text(
        ''.join(f'>{header}\n{sequence}\n' for header, sequence in records))
    written = _split_fasta(stage, 'sources/sccmec/sccmec-targets.fasta', 'sccmec/targets', _sccmec_split_name, None)
    assert sorted(group for group, _, _ in written) == ['ccrA1', 'mecA']
    by_group = {group: relative for group, relative, _ in written}
    with SequenceReader(stage / by_group['ccrA1'], None) as reader:
        staged = {record.identifier: record.sequence for record in reader}
    assert staged == {'ccrA1__AB033763.2__I': 'ACGTACGTAA', 'ccrA1__AB505628.1__IX': 'TTGGCCAATT'}


def test_a_split_that_would_produce_an_unsafe_filename_is_refused(tmp_path):
    from wmlstudio.characterization_refs import _sccmec_split_name
    stage = tmp_path / 'stage'
    (stage / 'sources/sccmec').mkdir(parents=True)
    (stage / 'sources/sccmec/sccmec-targets.fasta').write_text('>../escape AB1|I\nACGT\n')
    with pytest.raises(ValueError, match='not a safe identifier'):
        _split_fasta(stage, 'sources/sccmec/sccmec-targets.fasta', 'sccmec/targets', _sccmec_split_name, None)
    (stage / 'sources/sccmec/duplicate.fasta').write_text('>ccrA1 AB1|I\nACGT\n>ccrA1 AB1|I\nTTTT\n')
    with pytest.raises(ValueError, match='is not unique within'):
        _split_fasta(stage, 'sources/sccmec/duplicate.fasta', 'sccmec/dup', _sccmec_split_name, None)


def sccmec_panel(tmp_path, targets):
    root = tmp_path / 'sccmec-references'
    (root / 'sccmec/targets').mkdir(parents=True)
    entries = []
    for gene, sequence in targets.items():
        relative = f'sccmec/targets/{gene}.fasta'
        (root / relative).write_text(f'>{gene}__X1__IV\n{sequence}\n')
        entries.append({'gene': gene, 'path': relative, 'reference_count': 1})
    manifest = {'format_version': 2, 'source_revision': 'synthetic-truth',
                'source_repository': 'synthetic-test-fixture',
                'sources': {'sccmec': {'repository': 'synthetic', 'revision': 'synthetic-truth',
                                       'license': 'MIT', 'license_file': 'source-LICENSE-sccmec'}},
                'species': [], 'virulence': {}, 'locus_profiles': {}, 'capsule': {},
                'sccmec': {'source': 'sccmec', 'schema_version': '1.2.0', 'targets': entries, 'regions': [],
                           'rules': {'targets': sorted(targets), 'aliases': RULES['aliases'], 'types': RULES['types']},
                           'engine_defaults': {'target_min_pident': 90, 'target_min_coverage': 80},
                           'not_assayed': ['mecC']},
                'files': []}
    for path in sorted((root / 'sccmec/targets').iterdir()):
        relative = f'sccmec/targets/{path.name}'
        manifest['files'].append({'path': relative, 'bytes': path.stat().st_size, 'sha256': file_sha256(path)})
    manifest['reference_digest'] = reference_digest(manifest)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    return root


def test_native_blast_sccmec_target_screen_calls_a_single_candidate_type(tmp_path):
    tools = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not tools:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for the real native BLAST assay.')
    source = random.Random(9021)
    sequences = {gene: ''.join(source.choices('ACGT', k=900)) for gene in TYPE_IV}
    references = sccmec_panel(tmp_path, sequences)
    filler = ''.join(source.choices('ACGT', k=2_500_000))
    assembly = tmp_path / 'isolate.fasta'
    assembly.write_text('>contig_1\n' + filler[:500_000] + ''.join(sequences.values()) + filler[500_000:] + '\n')
    suffix = '.exe' if os.name == 'nt' else ''
    result = type_sccmec(assembly, references,
                         blastn_path=Path(tools) / ('blastn' + suffix),
                         makeblastdb_path=Path(tools) / ('makeblastdb' + suffix))
    assert result['assay']['adequate_negative_assay'] is True
    assert result['targets']['mecA'] == 'detected' and result['mecC'] == 'not_assayed'
    assert result['candidate_types'] == ['IV'] and result['type'] == 'IV'
    assert result['status'] == 'completed' and result['official_type'] is None
    assert result['contig_fragmented'] is False and result['contigs'] == ['contig_1']
    assert result['region_support']['status'] == 'not_run'
