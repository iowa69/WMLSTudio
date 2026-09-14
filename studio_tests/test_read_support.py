import copy
import os
import random
from pathlib import Path

import pytest

from wmlstudio import read_support
from wmlstudio.project import Project
from wmlstudio.read_support import (
    _aligned_bases,
    _summarize_candidate,
    current_read_support,
    investigate_read_support,
    persist_read_support,
    preflight_read_support,
)
from wmlstudio.sequence import AnalysisCancelled, SequenceRecord, file_sha256, file_signature
from wmlstudio.typing import Scheme, reverse_complement


def scheme_at(tmp_path, alleles=None):
    alleles = alleles or {'abc': {'1': 'ACGT' * 100}}
    return Scheme('truth', tmp_path, tuple(alleles), alleles, {}, {}, 'a' * 64)


def blast_tools():
    root = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not root:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for real native BLAST read-support controls.')
    suffix = '.exe' if os.name == 'nt' else ''
    return {name + '_path': str(Path(root) / (name + suffix)) for name in ('blastn', 'makeblastdb')}


def write_pairs(tmp_path, pairs, stem='reads', quality='I'):
    paths = [tmp_path / (stem + suffix) for suffix in ('_R1.fastq', '_R2.fastq')]
    for mate, path in enumerate(paths, 1):
        path.write_text(''.join(f'@read{index}/{mate}\n{pair[mate - 1]}\n+\n{quality * len(pair[mate - 1])}\n'
                                for index, pair in enumerate(pairs)))
    return paths


def tiled_pairs(sequence, depth=3):
    starts = sorted(set([*range(0, len(sequence) - 99, 25), len(sequence) - 100]))
    return [(sequence[start:start + 100], reverse_complement(sequence[start:start + 100]))
            for start in starts for _ in range(depth)]


def test_preflight_bounds_before_any_fastq_or_tool_access(tmp_path, monkeypatch):
    scheme = scheme_at(tmp_path, {'a': {'1': 'ACGT' * 50}, 'b': {'1': 'ACGT' * 50}})
    result = preflight_read_support(scheme, ['a', 'b'])
    assert result['allele_count'] == 2 and result['reference_bases'] == 400
    assert result['all_alleles_in_selected_loci_tested'] is True
    assert len(result['panel_sha256']) == 64
    monkeypatch.setattr(read_support, 'file_sha256', lambda *args: pytest.fail('Preflight must precede input hashing'))
    with pytest.raises(ValueError, match='exceeds'):
        investigate_read_support('absent', 'absent2', scheme, ['a', 'b'], max_loci=1)
    with pytest.raises(ValueError, match='exceeds'):
        preflight_read_support(scheme, ['a', 'b'], max_alleles=1)
    with pytest.raises(ValueError, match='exceeds'):
        preflight_read_support(scheme, ['a'], max_reference_bases=199)
    with pytest.raises(ValueError, match='unique'):
        preflight_read_support(scheme, ['a', 'a'])
    with pytest.raises(ValueError, match='unique'):
        preflight_read_support(scheme, ['unknown'])


def test_cancelled_does_not_touch_missing_input(tmp_path):
    with pytest.raises(AnalysisCancelled):
        investigate_read_support('missing', 'also-missing', scheme_at(tmp_path), ['abc'], cancelled=lambda: True)


def test_alignment_orientation_gaps_and_phred_base_filtering():
    record = SequenceRecord('truth', 'ACGT', 'I!II')
    row = {'qstart': '1', 'qend': '4', 'sstart': '4', 'send': '1', 'qseq': 'ACGT', 'sseq': 'ACGT'}
    observations, deletions, insertions = _aligned_bases(row, record, 'ACGT', 20)
    assert observations == {3: 'T', 1: 'C', 0: 'A'}
    assert not deletions and not insertions
    row = dict(row, sstart='1', send='5', qseq='AC-GT', sseq='ACAGT')
    observations, deletions, insertions = _aligned_bases(row, record, 'ACAGT', 20)
    assert observations == {0: 'A', 3: 'G', 4: 'T'} and not deletions  # Low-quality left flank.
    assert 2 not in observations  # A deleted reference base is not covered by a read base.
    assert _aligned_bases(row, SequenceRecord('truth', 'ACGT', 'IIII'), 'ACAGT', 20)[1] == {2}
    with pytest.raises(ValueError, match='actual queried sequence'):
        _aligned_bases(dict(row, qseq='AT-GT'), record, 'ACAGT', 20)


def test_summary_withholds_function_and_allele_assignment_for_disruption():
    from array import array

    entry = {'locus': 'truth', 'allele': '1', 'sequence': 'ACGT'}
    base = {'bases': [array('I', row) for row in ([6, 0, 0, 0], [0, 6, 0, 0], [0, 0, 6, 0], [0, 0, 0, 6])],
            'deletions': array('I', [0, 0, 0, 0]), 'insertions': {}, 'reads': 6, 'unique': 6, 'tied': 0}
    result = _summarize_candidate(entry, base, min_depth=3, min_breadth=.95)
    assert result['status'] == 'supported' and result['min_depth'] == 6
    assert result['breadth_min_depth'] == 1 and result['observed_identity_pct'] == 100
    mixed = copy.deepcopy(base)
    mixed['bases'][0][1] = 6
    result = _summarize_candidate(entry, mixed, min_depth=3, min_breadth=.95)
    assert result['status'] == 'mixed_support' and result['mixed_positions'] == 1
    disrupted = copy.deepcopy(base)
    disrupted['deletions'][2] = 6
    assert _summarize_candidate(entry, disrupted, min_depth=3, min_breadth=.95)['status'] == 'discordant_support'


def test_divergent_reference_review_does_not_override_a_consistent_supported_candidate():
    candidates = [{'locus': 'a', 'allele': '1', 'status': 'supported', 'uniquely_best_reads': 0},
                  {'locus': 'a', 'allele': '2', 'status': 'mixed_support', 'uniquely_best_reads': 0}]
    result = read_support._summarize_loci(candidates, ['a'])[0]
    assert result['status'] == 'supported' and result['compatible_candidates'] == ['1']
    assert result['candidates'][1]['status'] == 'mixed_support'


def test_real_blast_positive_alternative_disrupted_mixed_negative_and_sampled_controls(tmp_path):
    options = blast_tools()
    random_source = random.Random(7221)
    sequence = ''.join(random_source.choices('ACGT', k=350))
    sequence = sequence[:174] + 'ACG' + sequence[177:]
    alternative = sequence[:175] + ('A' if sequence[175] != 'A' else 'C') + sequence[176:]
    unrelated = ''.join(random_source.choices('ACGT', k=350))
    scheme = scheme_at(tmp_path, {'target': {'1': sequence, '2': alternative}, 'other': {'1': unrelated}})
    paths = write_pairs(tmp_path, tiled_pairs(sequence))
    result = investigate_read_support(*paths, scheme, ['target', 'other'], **options)
    first, second, negative = result['loci'][0]['candidates'] + result['loci'][1]['candidates']
    assert first['status'] == 'supported' and first['breadth_min_depth'] == 1
    assert first['min_depth'] >= 3 and first['observed_identity_pct'] == 100
    assert second['status'] == 'discordant_support' and second['discordant_positions'] == 1
    assert result['loci'][0]['compatible_candidates'] == ['1'] and result['loci'][0]['assigned_allele'] is None
    assert negative['status'] == 'no_support' and negative['breadth_1x'] == 0
    assert result['sampling']['complete_files'] and not result['sampling']['sampled']
    assert [read['sha256'] for read in result['reads']] == [file_sha256(path) for path in paths]
    assert len(result['provenance']['binary_sha256']) == 64
    assert first['tied_best_reads'] > 0  # Shared read stretches are not allele-unique.
    assert first['uniquely_best_reads'] > 0

    ambiguous = scheme_at(tmp_path, {'a': {'1': sequence}, 'b': {'1': sequence}})
    result = investigate_read_support(*paths, ambiguous, ['a', 'b'], **options)
    for locus in result['loci']:
        assert locus['status'] == 'ambiguous'
        assert locus['candidates'][0]['uniquely_best_reads'] == 0
        assert locus['candidates'][0]['tied_best_reads'] > 0

    mixed_paths = write_pairs(tmp_path, tiled_pairs(sequence) + tiled_pairs(alternative), 'mixed')
    result = investigate_read_support(*mixed_paths, scheme, ['target'], **options)
    assert result['loci'][0]['status'] == 'mixed_support'
    assert all(row['mixed_positions'] == 1 for row in result['loci'][0]['candidates'])

    deletion = sequence[:175] + sequence[176:]
    deletion_paths = write_pairs(tmp_path, tiled_pairs(deletion), 'deletion')
    result = investigate_read_support(*deletion_paths, scheme, ['target'], **options)
    assert result['loci'][0]['candidates'][0]['status'] == 'discordant_support'
    assert result['loci'][0]['candidates'][0]['deletion_positions'] == 1
    assert result['loci'][0]['candidates'][0]['breadth_1x'] < 1

    result = investigate_read_support(*paths, scheme, ['other'], max_pairs=1, **options)
    assert result['sampling']['sampled'] is True and result['sampling']['complete_files'] is False
    assert result['sampling']['unexamined_tail_structurally_validated'] is False
    assert result['loci'][0]['status'] == 'no_support'
    assert any('does not establish locus absence' in note for note in result['limitations'])


def test_real_blast_low_quality_and_prefix_pairing_and_explicit_output_limits(tmp_path):
    options = blast_tools()
    sequence = ''.join(random.Random(978).choices('ACGT', k=350))
    scheme = scheme_at(tmp_path, {'a': {'1': sequence}})
    paths = write_pairs(tmp_path, tiled_pairs(sequence), quality='!')
    result = investigate_read_support(*paths, scheme, ['a'], **options)
    assert result['loci'][0]['status'] == 'no_support'
    assert result['loci'][0]['candidates'][0]['breadth_1x'] == 0
    with pytest.raises(ValueError, match='safety bound'):
        investigate_read_support(*paths, scheme, ['a'], max_output_bytes=1, **options)
    with pytest.raises(ValueError, match='SHA-256'):
        investigate_read_support(*paths, scheme, ['a'], expected_read_sha256=['f' * 64, 'f' * 64], **options)
    paths[1].write_text(paths[1].read_text().replace('@read0/2', '@different/2'))
    with pytest.raises(ValueError, match='mate identity'):
        investigate_read_support(*paths, scheme, ['a'], **options)


def test_read_support_persistence_preserves_calls_and_excludes_changed_reads(tmp_path):
    assembly = tmp_path / 'assembly.fasta'
    assembly.write_text('>contig\nACGTACGT\n')
    paths = write_pairs(tmp_path, [('ACGT' * 20, 'ACGT' * 20)])
    with Project(tmp_path / 'project.wmlstudio') as project:
        sid = project.add_sample(assembly)
        project.set_result(sid, {'status': 'missing', 'alleles': {'a': None}, 'input_sha256': file_sha256(assembly)})
        reads = [{'path': str(path), 'sha256': file_sha256(path), 'signature': list(file_signature(path)), 'mate': index}
                 for index, path in enumerate(paths, 1)]
        project.update_metadata(sid, {'reads': {'reads': reads}})
        result = {'format_version': 1, 'input_path': str(assembly), 'input_sha256': file_sha256(assembly),
                  'assembly_signature': list(file_signature(assembly)), 'reads': reads, 'scheme_digest': 'a' * 64,
                  'panel': {'loci': ['a']}, 'loci': [{'locus': 'a', 'status': 'supported', 'assigned_allele': None}]}
        before = project.get_sample(sid)['result']
        record = persist_read_support(project, sid, result)
        assert record['result'] == before
        assert current_read_support(record)['status'] == 'current'
        record['metadata']['reads']['reads'][0]['sha256'] = 'b' * 64
        assert current_read_support(record)['status'] == 'stale'
        assert current_read_support(record)['evidence'] is None
        paths[0].write_text(paths[0].read_text() + '\n')
        with pytest.raises(ValueError, match='input changed'):
            persist_read_support(project, sid, result)
