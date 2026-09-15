import os
import random
from pathlib import Path

import pytest

from wmlstudio.cgmlst_calls import (
    CALL_STATES,
    RECHECK_VERDICTS,
    build_calls_table,
    called_targets,
    classify_call,
    interpret_target_recheck,
    missing_targets,
    plan_target_recheck,
    recheck_missing_targets,
)
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled, file_sha256, file_signature
from wmlstudio.typing import Scheme, call_assembly, load_scheme, reverse_complement

DIGEST = 'a' * 64


def call_row(locus, status, allele=None, **extra):
    return {'locus': locus, 'status': status, 'allele': allele, 'reason': f'{status} evidence', **extra}


def hit(contig='contig1', start=1, end=300):
    return {'contig': contig, 'start': start, 'end': end, 'strand': '+'}


def profile(calls, *, sample_id='s1', name=None, scheme='cgmlst-demo', digest=DIGEST, status='incomplete'):
    return {'sample_id': sample_id, 'sample_name': name or sample_id, 'scheme': scheme,
            'scheme_digest': digest, 'status': status, 'input_sha256': 'b' * 64,
            'alleles': {call['locus']: call.get('allele') for call in calls}, 'calls': list(calls)}


def every_state_profile(sample_id='s1'):
    return profile([
        call_row('a', 'exact', '1'),
        call_row('b', 'novel_validated', 'NOVEL_' + 'b' * 64),
        call_row('c', 'missing'),
        call_row('d', 'partial'),
        call_row('e', 'low_similarity'),
        call_row('f', 'ambiguous', hits=[hit(), hit(start=1, end=300)]),
        call_row('g', 'ambiguous', hits=[hit(), hit('contig2', 5000, 5299)]),
        call_row('h', 'mixed', hits=[hit(), hit('contig2', 5000, 5299)]),
        call_row('i', 'a_state_from_the_future'),
    ], sample_id=sample_id)


def support_locus(locus, status, candidates=None, compatible=()):
    rows = candidates if candidates is not None else [
        {'locus': locus, 'allele': '1', 'status': 'no_support', 'breadth_min_depth': 0.0,
         'mean_depth': 0.0, 'supporting_reads': 0}]
    return {'locus': locus, 'status': status, 'candidates': rows,
            'compatible_candidates': list(compatible), 'assigned_allele': None,
            'interpretation': 'Support is restricted to tested references and aligned reads.'}


def support_payload(loci, *, complete=True, pairs=5000, digest=DIGEST):
    return {'format_version': 1, 'status': 'completed', 'scheme': 'cgmlst-demo', 'scheme_digest': digest,
            'loci': list(loci), 'sampling': {'complete_files': complete, 'sampled': not complete,
                                             'pairs_checked': pairs, 'reads_aligned': pairs * 2},
            'panel': {'loci': [row['locus'] for row in loci]}, 'parameters': {'min_depth': 3},
            'reads': [], 'input_sha256': 'b' * 64, 'provenance': {}}


def write_pairs(tmp_path, pairs, stem='reads'):
    paths = [tmp_path / (stem + suffix) for suffix in ('_R1.fastq', '_R2.fastq')]
    for mate, path in enumerate(paths, 1):
        path.write_text(''.join(f'@read{index}/{mate}\n{pair[mate - 1]}\n+\n{"I" * len(pair[mate - 1])}\n'
                                for index, pair in enumerate(pairs)))
    return paths


def tiled_pairs(sequence, step=50, depth=3):
    starts = sorted(set([*range(0, max(1, len(sequence) - 99), step), max(0, len(sequence) - 100)]))
    return [(sequence[start:start + 100], reverse_complement(sequence[start:start + 100]))
            for start in starts for _ in range(depth)]


# --- the table of calls -----------------------------------------------------

def test_every_reason_for_no_call_stays_its_own_state_and_never_merges_into_missing():
    table = build_calls_table([every_state_profile()])
    states = {row['locus']: row['cells']['s1']['state'] for row in table['loci']}
    assert states == {'a': 'called', 'b': 'called_novel', 'c': 'missing', 'd': 'partial',
                      'e': 'low_similarity', 'f': 'ambiguous', 'g': 'duplicated', 'h': 'mixed',
                      'i': 'unknown'}
    assert table['state_totals']['missing'] == 1
    assert sum(table['state_totals'].values()) == 9
    assert [row['cells']['s1']['display'] for row in table['loci']][:2] == ['1', 'local novel · bbbbbbbbbbbb']
    assert dict(states)['g'] == 'duplicated'
    assert table['loci'][6]['cells']['s1']['copies'] == 2
    assert table['loci'][6]['cells']['s1']['display'] == 'duplicated ×2'
    # The state vocabulary travels with the payload so the page never invents wording.
    assert set(table['states']) == set(CALL_STATES)
    assert 'not evidence that the target is absent' in table['states']['missing']['meaning']


def test_percentages_carry_the_schemes_full_target_count_as_their_denominator():
    table = build_calls_table([every_state_profile()])
    sample = table['samples'][0]
    assert table['target_count'] == 9
    assert sample['targets'] == 9 and sample['called'] == 2 and sample['without_call'] == 7
    assert sample['called_reference'] == 1 and sample['called_novel'] == 1
    assert sample['called_pct'] == pytest.approx(100 * 2 / 9)
    assert sample['summary'] == '2 of 9 targets called (22.2% of the scheme)'
    assert '9 targets' in table['denominator']
    assert 'not over the targets that happened to call' in table['denominator_note']


def test_locus_rows_and_isolate_summaries_count_the_same_cells():
    table = build_calls_table([every_state_profile('one'), every_state_profile('two')])
    assert table['isolate_count'] == 2
    assert [row['sample_id'] for row in table['samples']] == ['one', 'two']
    assert all(row['isolates'] == 2 and row['called_pct'] in (0.0, 100.0) for row in table['loci'])
    assert table['targets_called_in_every_isolate'] == 2
    assert table['targets_called_in_no_isolate'] == 7
    for state in CALL_STATES:
        assert table['state_totals'][state] == sum(row['states'][state] for row in table['loci'])
        assert table['state_totals'][state] == sum(row['states'][state] for row in table['samples'])


def test_a_stored_allele_and_its_evidence_must_agree_or_the_cell_is_unknown():
    assert classify_call(call_row('a', 'exact'), None)['state'] == 'unknown'
    assert 'no allele is stored' in classify_call(call_row('a', 'exact'), None)['reason']
    assert classify_call(call_row('a', 'missing'), '7')['state'] == 'unknown'
    assert 'assigns none' in classify_call(call_row('a', 'missing'), '7')['reason']
    # A placeholder under a no-allele state is the absence of an allele, not a conflict.
    assert classify_call(call_row('a', 'missing'), '-')['state'] == 'missing'


def test_an_imported_profile_without_call_evidence_is_read_without_inventing_evidence():
    imported = {'sample_id': 'i1', 'scheme': 'cgmlst-demo', 'scheme_digest': DIGEST,
                'alleles': {'a': '4', 'b': '-', 'c': None, 'd': 'novel'}}
    table = build_calls_table([imported])
    states = {row['locus']: row['cells']['i1'] for row in table['loci']}
    assert states['a']['state'] == 'called' and 'without per-target call evidence' in states['a']['reason']
    assert states['b']['state'] == 'missing' and states['c']['state'] == 'missing'
    assert states['d']['state'] == 'unknown' and 'without a sequence identity' in states['d']['reason']
    assert table['samples'][0]['called'] == 1


def test_duplicate_evidence_beyond_the_stored_cap_is_reported_as_a_lower_bound():
    cell = classify_call(call_row('g', 'ambiguous', hits=[hit(), hit('contig2', 10, 309)],
                                  evidence_truncated=True), None)
    assert cell['state'] == 'duplicated' and cell['copies_lower_bound'] is True
    assert classify_call(call_row('g', 'ambiguous', hits=[hit(), hit('contig2', 10, 309)],
                                  evidence_truncated=True), None)['copies'] == 2


def test_overlapping_hits_for_one_place_are_one_copy_not_a_duplication():
    cell = classify_call(call_row('f', 'ambiguous', hits=[hit(start=1, end=300), hit(start=250, end=549)]), None)
    assert cell['copies'] == 1 and cell['state'] == 'ambiguous'


def test_cross_locus_cds_candidates_count_as_separate_copies():
    cell = classify_call(call_row('g', 'ambiguous', cds_candidates=['cds1', 'cds2']), None)
    assert cell['copies'] == 2 and cell['state'] == 'duplicated'


def test_one_table_covers_one_scheme_and_one_target_set():
    left, right = every_state_profile('one'), every_state_profile('two')
    with pytest.raises(ValueError, match='different scheme fingerprints'):
        build_calls_table([left, dict(right, scheme_digest='c' * 64)])
    with pytest.raises(ValueError, match='different scheme names'):
        build_calls_table([left, dict(right, scheme='another')])
    with pytest.raises(ValueError, match='different targets'):
        build_calls_table([left, dict(right, alleles={**right['alleles'], 'zz': None})])
    with pytest.raises(ValueError, match='fingerprints'):
        build_calls_table([dict(left, scheme_digest='')])
    with pytest.raises(ValueError, match='appears once'):
        build_calls_table([left, every_state_profile('one')])
    with pytest.raises(ValueError, match='at least one'):
        build_calls_table([])


def test_a_loaded_scheme_supplies_the_target_order_and_is_checked_against_the_profiles(tmp_path):
    alleles = {locus: {'1': 'ACGT' * 25} for locus in 'abcdefghi'}
    scheme = Scheme('cgmlst-demo', tmp_path, ('i', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'), alleles, {}, {}, DIGEST)
    table = build_calls_table([every_state_profile()], scheme=scheme)
    assert [row['locus'] for row in table['loci']][0] == 'i'
    mismatched = Scheme('cgmlst-demo', tmp_path, tuple('abcdefghi'), alleles, {}, {}, 'c' * 64)
    with pytest.raises(ValueError, match='not the one these profiles were called against'):
        build_calls_table([every_state_profile()], scheme=mismatched)
    with pytest.raises(ValueError, match='loaded Scheme'):
        build_calls_table([every_state_profile()], scheme=str(tmp_path))


def test_the_typing_scale_is_stated_and_never_guessed_from_a_call_count():
    small = build_calls_table([every_state_profile()])
    assert small['typing']['kind'] == 'mlst'
    assert 'different quantities' in small['typing']['note']
    wide = profile([call_row(f'tag{index:03d}', 'exact', '1') for index in range(40)])
    assert build_calls_table([wide])['typing']['kind'] == 'cgmlst'
    assert build_calls_table([wide])['samples'][0]['analysis_kind']['kind'] == 'cgmlst'


def test_building_the_table_is_cancellable_and_reports_progress():
    with pytest.raises(AnalysisCancelled):
        build_calls_table([every_state_profile()], cancelled=lambda: True)
    messages = []
    build_calls_table([every_state_profile()], progress=lambda done, total, message: messages.append((done, total, message)))
    assert messages[0][1] == 9 and messages[-1][0] == 9
    assert 'Table of calls ready' in messages[-1][2]


def test_missing_and_called_targets_follow_the_tables_own_row_order():
    table = build_calls_table([every_state_profile()])
    assert missing_targets(table, 's1') == ['c', 'd', 'e', 'f', 'g', 'h', 'i']
    assert called_targets(table, 's1') == ['a']  # A local novel sequence is not in the scheme.
    with pytest.raises(ValueError, match='not one of the isolates'):
        missing_targets(table, 'absent')


# --- planning the read recheck ---------------------------------------------

def wide_table(uncalled=5, called=10, sample_id='s1'):
    calls = [call_row(f'called{index:03d}', 'exact', str(index + 1)) for index in range(called)]
    calls += [call_row(f'gap{index:03d}', 'missing') for index in range(uncalled)]
    return build_calls_table([profile(calls, sample_id=sample_id)])


def test_the_plan_bounds_the_panel_and_names_every_deferred_target():
    plan = plan_target_recheck(wide_table(uncalled=30), 's1')
    assert len(plan['panel_loci']) == 20
    assert len(plan['targets']) == 18 and len(plan['controls']) == 2
    assert len(plan['deferred_targets']) == 12
    assert plan['uncalled_total'] == 30
    assert set(plan['targets']) | set(plan['deferred_targets']) == set(missing_targets(wide_table(uncalled=30), 's1'))
    assert '12 are listed in this plan and are not checked until you run them' in plan['batch_note']


def test_controls_are_exact_calls_spread_across_the_scheme_not_its_first_corner():
    plan = plan_target_recheck(wide_table(), 's1', control_count=2)
    assert plan['controls'] == ['called000', 'called005']
    assert all(control in called_targets(wide_table(), 's1') for control in plan['controls'])
    assert 'which is what makes a negative result readable' in plan['control_basis']


def test_a_panel_without_a_possible_control_says_so_instead_of_pretending():
    table = build_calls_table([profile([call_row('gap0', 'missing'), call_row('gap1', 'partial')])])
    plan = plan_target_recheck(table, 's1')
    assert plan['controls'] == []
    assert 'stays inconclusive' in plan['control_basis']


def test_the_plan_refuses_targets_that_already_carry_a_call_or_do_not_exist():
    table = wide_table()
    with pytest.raises(ValueError, match='already carry a call'):
        plan_target_recheck(table, 's1', loci=['called000'])
    with pytest.raises(ValueError, match='not in this scheme'):
        plan_target_recheck(table, 's1', loci=['nowhere'])
    with pytest.raises(ValueError, match='appears once'):
        plan_target_recheck(table, 's1', loci=['gap000', 'gap000'])
    with pytest.raises(ValueError, match='must be one this isolate called exactly'):
        plan_target_recheck(table, 's1', controls=['gap000'])
    complete = build_calls_table([profile([call_row('a', 'exact', '1')])])
    with pytest.raises(ValueError, match='nothing to recheck'):
        plan_target_recheck(complete, 's1')


def test_the_plan_carries_the_isolates_verified_read_pair_or_refuses_to_plan(tmp_path):
    assembly = tmp_path / 'assembly.fasta'
    assembly.write_text('>contig\nACGTACGT\n')
    paths = write_pairs(tmp_path, [('ACGT' * 25, 'ACGT' * 25)])
    reads = [{'path': str(path), 'sha256': file_sha256(path), 'mate': index}
             for index, path in enumerate(paths, 1)]
    with Project(tmp_path / 'project.wmlstudio') as project:
        sample_id = project.add_sample(assembly)
        project.update_metadata(sample_id, {'reads': {'reads': reads}})
        record = project.get_sample(sample_id)
    table = wide_table(sample_id=record['id'])
    plan = plan_target_recheck(table, record['id'], sample=record)
    assert plan['read_paths'] == [str(path) for path in paths]
    assert plan['expected_read_sha256'] == [file_sha256(path) for path in paths]
    assert plan['assembly_path'] == record['input_path']
    with pytest.raises(ValueError, match='verified original FASTQ pair'):
        plan_target_recheck(table, record['id'], sample={**record, 'metadata': {}})


def test_a_recheck_without_attached_reads_never_reaches_a_tool(tmp_path, monkeypatch):
    from wmlstudio import cgmlst_calls

    monkeypatch.setattr(cgmlst_calls, 'investigate_read_support',
                        lambda *args, **kwargs: pytest.fail('No reads were verified for this isolate'))
    plan = plan_target_recheck(wide_table(), 's1')
    with pytest.raises(ValueError, match='verified original FASTQ pair'):
        recheck_missing_targets(plan, Scheme('cgmlst-demo', tmp_path, ('gap000',), {'gap000': {'1': 'ACGT'}}, {}, {}, DIGEST))


def test_a_recheck_is_cancelled_before_any_tool_runs(tmp_path, monkeypatch):
    from wmlstudio import cgmlst_calls

    monkeypatch.setattr(cgmlst_calls, 'investigate_read_support',
                        lambda *args, **kwargs: pytest.fail('Cancellation must precede the read assay'))
    plan = plan_target_recheck(wide_table(), 's1')
    plan.update(read_paths=['r1', 'r2'], expected_read_sha256=['a' * 64, 'b' * 64])
    with pytest.raises(AnalysisCancelled):
        recheck_missing_targets(plan, None, cancelled=lambda: True)


# --- reading the reads ------------------------------------------------------

def test_reads_that_recover_an_uncalled_target_report_presence_and_never_an_allele():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', 'supported', compatible=['1', '2']),
                               support_locus(plan['controls'][0], 'supported'),
                               support_locus(plan['controls'][1], 'supported')])
    result = interpret_target_recheck(support, plan, table=table)
    row = result['targets'][0]
    assert row['verdict'] == 'present_in_reads' and row['evidence'] == 'present'
    assert row['assigned_allele'] is None
    assert row['compatible_candidates'] == ['1', '2']
    assert 'Compatible is not called' in row['candidate_note']
    assert row['assembly_state'] == 'missing' and row['assembly_state_label'] == 'no call'
    assert 'coverage gap, an assembly break or a contig edge' in row['explanation']
    assert result['summary'] == {'present': 1, 'absent': 0, 'insufficient': 0, 'checked': 1,
                                 'uncalled_total': 1,
                                 'verdicts': {**{key: 0 for key in RECHECK_VERDICTS}, 'present_in_reads': 1}}


def test_absence_is_readable_only_with_a_recovered_control_and_complete_read_files():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    controls = plan['controls']
    passing = support_payload([support_locus('gap000', 'no_support'),
                               support_locus(controls[0], 'supported'),
                               support_locus(controls[1], 'no_support')])
    result = interpret_target_recheck(passing, plan)
    assert result['targets'][0]['verdict'] == 'consistent_with_absence'
    assert result['targets'][0]['evidence'] == 'absent'
    assert 'It is not proof' in result['targets'][0]['explanation']
    assert result['controls']['passed'] == 1 and result['controls']['absence_readable'] is True

    failed = support_payload([support_locus('gap000', 'no_support'),
                              support_locus(controls[0], 'no_support'),
                              support_locus(controls[1], 'no_support')])
    result = interpret_target_recheck(failed, plan)
    assert result['targets'][0]['verdict'] == 'no_evidence_either_way'
    assert result['targets'][0]['evidence'] == 'insufficient'
    assert 'not shown to work on this read set' in result['targets'][0]['explanation']
    assert result['controls']['absence_readable'] is False


def test_a_sampled_read_prefix_can_never_establish_absence():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', 'no_support'),
                               *(support_locus(locus, 'supported') for locus in plan['controls'])],
                              complete=False, pairs=1000)
    result = interpret_target_recheck(support, plan)
    assert result['targets'][0]['verdict'] == 'no_evidence_either_way'
    assert 'only the first 1,000 read pairs were examined' in result['targets'][0]['explanation']


def test_a_panel_with_no_control_at_all_leaves_every_negative_inconclusive():
    table = build_calls_table([profile([call_row('gap0', 'missing')])])
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap0', 'no_support')])
    result = interpret_target_recheck(support, plan)
    assert result['controls']['tested'] == 0
    assert result['targets'][0]['verdict'] == 'no_evidence_either_way'
    assert 'No control target was included' in result['targets'][0]['explanation']


@pytest.mark.parametrize('status,verdict,evidence', [
    ('supported', 'present_in_reads', 'present'),
    ('ambiguous', 'present_in_reads', 'present'),
    ('mixed_support', 'present_with_mixed_bases', 'present'),
    ('incomplete_or_discordant', 'partial_read_evidence', 'insufficient'),
    ('no_support', 'consistent_with_absence', 'absent'),
])
def test_each_read_support_state_maps_to_one_stated_verdict(status, verdict, evidence):
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', status),
                               *(support_locus(locus, 'supported') for locus in plan['controls'])])
    row = interpret_target_recheck(support, plan)['targets'][0]
    assert (row['verdict'], row['evidence']) == (verdict, evidence)
    assert row['assigned_allele'] is None


def test_a_read_support_state_this_version_does_not_know_is_never_read_as_absence():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', 'a_state_from_the_future'),
                               *(support_locus(locus, 'supported') for locus in plan['controls'])])
    row = interpret_target_recheck(support, plan)['targets'][0]
    assert row['verdict'] == 'no_evidence_either_way' and row['evidence'] == 'insufficient'
    assert 'does not recognise' in row['explanation']


def test_partial_read_evidence_is_neither_presence_nor_absence():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', 'incomplete_or_discordant'),
                               *(support_locus(locus, 'supported') for locus in plan['controls'])])
    row = interpret_target_recheck(support, plan)['targets'][0]
    assert 'neither the presence of a complete target nor evidence of absence' in row['explanation']


def test_a_recheck_refuses_a_different_scheme_or_an_incomplete_investigation():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    with pytest.raises(ValueError, match='different scheme'):
        interpret_target_recheck(support_payload([], digest='c' * 64), plan)
    with pytest.raises(ValueError, match='did not cover every target'):
        interpret_target_recheck(support_payload([support_locus('gap000', 'supported')]), plan)


def test_the_returned_payload_states_its_limits_in_the_words_the_page_shows():
    table = wide_table(uncalled=1)
    plan = plan_target_recheck(table, 's1')
    support = support_payload([support_locus('gap000', 'supported'),
                               *(support_locus(locus, 'supported') for locus in plan['controls'])])
    result = interpret_target_recheck(support, plan)
    joined = ' '.join(result['limitations'])
    assert 'It never assigns an allele' in joined
    assert 'is not evidence of absence in the isolate' in joined
    assert 'no distance and no tree' in joined
    # The read assay's own limits travel with the verdicts, not only with the assay.
    assert any('does not establish locus absence' in note for note in result['limitations'])
    assert set(result['verdict_meanings']) == set(RECHECK_VERDICTS)
    assert result['headline'].startswith('1 present in the reads · 0 consistent with absence')
    assert 'different quantities' in ' '.join(build_calls_table([profile([call_row('a', 'exact', '1')])])['limitations'])


# --- end to end against real reads -----------------------------------------

def test_real_blast_recheck_tells_an_assembly_gap_apart_from_an_absent_target(tmp_path):
    binary_dir = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not binary_dir:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for the native BLAST read recheck.')
    source = random.Random(4242)
    loci = {f'tag{index:02d}': ''.join(source.choices('ACGT', k=300)) for index in range(32)}
    loci['gap_target'] = ''.join(source.choices('ACGT', k=300))
    loci['absent_target'] = ''.join(source.choices('ACGT', k=300))
    directory = tmp_path / 'scheme'
    directory.mkdir()
    for locus, sequence in loci.items():
        (directory / f'{locus}.fasta').write_text(f'>{locus}_1\n{sequence}\n')
    scheme = load_scheme(directory)
    assert len(scheme.loci) == 34

    # The assembly holds every target except the two under test; the reads hold the
    # assembly plus gap_target, so an assembly gap and a true absence differ only in
    # what the reads show, which is exactly the distinction the button must make.
    assembled = ''.join(loci[locus] for locus in sorted(loci) if locus not in {'gap_target', 'absent_target'})
    assembly = tmp_path / 'isolate.fasta'
    assembly.write_text('>contig1\n' + assembled + '\n')
    paths = write_pairs(tmp_path, tiled_pairs(assembled + loci['gap_target']))
    result = call_assembly(assembly, scheme)
    table = build_calls_table([dict(result, sample_id='isolate')])
    assert table['typing']['kind'] == 'cgmlst'
    assert table['target_count'] == 34
    assert table['samples'][0]['called'] == 32
    assert missing_targets(table, 'isolate') == ['absent_target', 'gap_target']

    reads = [{'path': str(path), 'sha256': file_sha256(path), 'signature': list(file_signature(path)), 'mate': index}
             for index, path in enumerate(paths, 1)]
    with Project(tmp_path / 'project.wmlstudio') as project:
        sample_id = project.add_sample(assembly)
        project.update_metadata(sample_id, {'reads': {'reads': reads}})
        record = project.get_sample(sample_id)
    table = build_calls_table([dict(result, sample_id=record['id'])])
    plan = plan_target_recheck(table, record['id'], sample=record, control_count=2)
    assert plan['targets'] == ['absent_target', 'gap_target']
    assert len(plan['controls']) == 2

    suffix = '.exe' if os.name == 'nt' else ''
    verdicts = recheck_missing_targets(
        plan, scheme, table=table,
        blastn_path=str(Path(binary_dir) / ('blastn' + suffix)),
        makeblastdb_path=str(Path(binary_dir) / ('makeblastdb' + suffix)))
    by_locus = {row['locus']: row for row in verdicts['targets']}
    assert verdicts['controls']['passed'] == 2 and verdicts['controls']['absence_readable'] is True
    assert by_locus['gap_target']['verdict'] == 'present_in_reads'
    assert by_locus['gap_target']['compatible_candidates'] == ['1']
    assert by_locus['gap_target']['assigned_allele'] is None
    assert by_locus['absent_target']['verdict'] == 'consistent_with_absence'
    assert by_locus['absent_target']['supporting_reads'] == 0
    assert verdicts['summary']['present'] == 1 and verdicts['summary']['absent'] == 1
    # The read check leaves the table of calls exactly as it found it.
    assert missing_targets(table, record['id']) == ['absent_target', 'gap_target']
