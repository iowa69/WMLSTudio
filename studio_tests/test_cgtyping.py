import copy
import hashlib
import os
import random
import sys
from pathlib import Path

import pytest

from wmlstudio.cgtyping import (
    _apply_inference,
    _run,
    call_cgassembly,
    full_cds_qc,
    predict_cds,
)
from wmlstudio.comparison import pairwise_distances
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.typing import call_assembly, load_scheme, reverse_complement


def cds(seed=12, codons=120):
    rng = random.Random(seed)
    sense = ['GCT', 'GCA', 'AAA', 'GAA', 'GAT', 'CTG', 'ATT', 'GGT', 'TTC', 'TAT']
    return 'ATG' + ''.join(rng.choice(sense) for _ in range(codons)) + 'TAA'


def gene(sequence, identifier='cds0', start=20, **qc):
    return {'id': identifier, 'sequence': sequence,
            'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'contig': 'contig', 'start': start, 'end': start + len(sequence) - 1, 'strand': '+',
            'cds_qc': full_cds_qc(sequence, **qc)}


def alignment(locus='abc', identifier='cds0', **kwargs):
    return {'locus': locus, 'gene_id': identifier, 'reference_allele': '1',
            'identity': 0.99, 'query_coverage': 1.0, 'subject_coverage': 1.0,
            'bitscore': 500, **kwargs}


@pytest.fixture
def inference_setup(tmp_path):
    directory = tmp_path / 'schema'
    directory.mkdir()
    reference = cds()
    (directory / 'abc.fasta').write_text(f'>abc_1\n{reference}\n')
    (directory / 'profiles.tsv').write_text('ST\tabc\n42\t1\n')
    path = tmp_path / 'input.fa'
    path.write_text('>contig\n' + reference[:90] + 'C' + reference[91:] + '\n')
    scheme = load_scheme(directory)
    baseline = call_assembly(path, scheme)
    assert baseline['alleles'] == {'abc': None}
    return path, scheme, baseline


@pytest.mark.parametrize('sequence,kwargs,reason', [
    ('ATGAAATAA', {}, None),
    ('GTGAAATAG', {}, None),
    ('ATGAAATAA', {'partial_begin': True}, 'contig-edge'),
    ('ATGAAATAA', {'partial_end': True}, 'contig-edge'),
    ('ATGANA TAA'.replace(' ', ''), {}, 'ambiguous'),
    ('ATGAAAATAA', {}, 'multiple of three'),
    ('ACGAAATAA', {}, 'start codon'),
    ('ATGAAAGGG', {}, 'terminal stop'),
    ('ATGTAAAAATAA', {}, 'internal'),
    ('ATGTGATAA', {}, 'internal'),
    ('ATGTGATAA', {'genetic_code': 4}, None),
    ('ATGAAATAA', {'genetic_code': 2}, 'unsupported'),
])
def test_full_cds_quality_gates(sequence, kwargs, reason):
    result = full_cds_qc(sequence, **kwargs)
    assert result['valid'] is (reason is None)
    if reason:
        assert any(reason in message for message in result['reasons'])


def test_novel_identity_is_full_sequence_sha_never_nearest_st(inference_setup):
    _, scheme, baseline = inference_setup
    novel = cds(seed=13)
    result = _apply_inference(baseline, scheme, [gene(novel)], [alignment()], .9, .98)
    digest = hashlib.sha256(novel.encode()).hexdigest()
    assert result['alleles'] == {'abc': 'NOVEL_' + digest}
    assert result['st'] is None
    assert result['status'] == 'novel_alleles'
    assert result['calls'][0]['status'] == 'novel_validated'
    assert result['novel_sequences'][0]['sequence'] == novel
    assert result['cg_profile_complete'] is True
    assert result['cg_profile_id'].startswith('LOCAL_')


@pytest.mark.parametrize('kwargs', [{'partial_begin': True}, {'partial_end': True}])
def test_high_identity_partial_cds_is_never_an_allele(inference_setup, kwargs):
    _, scheme, baseline = inference_setup
    result = _apply_inference(baseline, scheme, [gene(cds(13), **kwargs)], [alignment()], .9, .98)
    assert result['alleles']['abc'] is None
    assert result['calls'][0]['status'] == 'partial'
    assert result['novel_sequences'] == []
    assert result['cg_profile_complete'] is False


@pytest.mark.parametrize('kwargs', [{'identity': .899}, {'query_coverage': .979}, {'subject_coverage': .979}])
def test_both_full_alignment_coverages_and_identity_required(inference_setup, kwargs):
    _, scheme, baseline = inference_setup
    result = _apply_inference(baseline, scheme, [gene(cds(13))], [alignment(**kwargs)], .9, .98)
    assert result['alleles']['abc'] is None
    assert result['calls'][0]['status'] == 'low_similarity'


def test_two_full_cds_copies_do_not_collapse_even_when_sequences_equal(inference_setup):
    _, scheme, baseline = inference_setup
    genes = [gene(cds(13)), gene(cds(13), 'cds1', start=1000)]
    result = _apply_inference(baseline, scheme, genes,
                              [alignment(), alignment(identifier='cds1', bitscore=490)], .9, .98)
    assert result['alleles']['abc'] is None
    assert result['status'] == 'ambiguous'


def test_cross_locus_homology_cannot_assign_same_gene_twice(inference_setup):
    _, scheme, baseline = inference_setup
    baseline['calls'].append({'locus': 'def', 'status': 'missing', 'allele': None})
    baseline['alleles']['def'] = None
    result = _apply_inference(baseline, scheme, [gene(cds(13))],
                              [alignment(), alignment(locus='def')], .9, .98)
    assert result['alleles'] == {'abc': None, 'def': None}
    assert result['status'] == 'ambiguous'


def test_exact_plus_novel_full_copy_is_mixed(inference_setup):
    path, scheme, _ = inference_setup
    path.write_text('>contig\n' + scheme.alleles['abc']['1'] + '\n')
    baseline = call_assembly(path, scheme)
    assert baseline['st'] == '42'
    result = _apply_inference(baseline, scheme, [gene(cds(13), start=1000)], [alignment()], .9, .98)
    assert result['st'] is None
    assert result['status'] == 'mixed'
    assert result['alleles']['abc'] is None


def test_novel_comparison_uses_sequence_evidence_not_placeholder(inference_setup):
    _, scheme, baseline = inference_setup
    first = _apply_inference(copy.deepcopy(baseline), scheme, [gene(cds(13))], [alignment()], .9, .98)
    first['sample_id'] = 'first'
    second = copy.deepcopy(first)
    second['sample_id'] = 'second'
    assert pairwise_distances([first, second])[0]['distance'] == 0
    second = _apply_inference(copy.deepcopy(baseline), scheme, [gene(cds(14))], [alignment()], .9, .98)
    second['sample_id'] = 'second'
    assert pairwise_distances([first, second])[0]['distance'] == 1
    second['calls'][0]['sequence_sha256'] = 'a' * 64
    assert pairwise_distances([first, second])[0]['comparable'] is False
    second['calls'][0]['cds_qc']['valid'] = False
    assert pairwise_distances([first, second])[0]['shared_loci'] == 0


def test_subprocess_cancellation_terminates_real_process():
    calls = 0
    def cancelled():
        nonlocal calls
        calls += 1
        return calls > 2
    with pytest.raises(AnalysisCancelled):
        _run([sys.executable, '-c', 'import time; time.sleep(20)'], cancelled=cancelled)


def test_real_pyrodigal_returns_complete_sense_cds(tmp_path):
    sequence = cds(codons=400)
    path = tmp_path / 'synthetic.fa'
    path.write_text('>contig\n' + 'TAA' * 100 + sequence + 'TAA' * 100 + '\n')
    genes, provenance = predict_cds(path)
    assert provenance['tool'] == 'pyrodigal'
    assert provenance['mode'] == 'meta-small-input'
    assert any(g['sequence'] == sequence and g['cds_qc']['valid'] for g in genes)


def test_native_blast_and_pyrodigal_novel_and_cache(tmp_path):
    binary_dir = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not binary_dir:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for native BLAST integration.')
    sequence = cds(codons=400)
    path = tmp_path / 'sample.fa'
    path.write_text('>contig\n' + 'TAA' * 100 + sequence + 'TAA' * 100 + '\n')
    # The input contains a complete novel gene; the reference differs by one
    # internal codon. Truth does not depend on the returned alignment.
    reference = sequence[:300] + ('GCT' if sequence[300:303] != 'GCT' else 'GCA') + sequence[303:]
    schema = tmp_path / 'scheme'
    schema.mkdir()
    (schema / 'abc.fasta').write_text('>abc_1\n' + reference + '\n')
    (schema / 'profiles.tsv').write_text('ST\tabc\n42\t1\n')
    suffix = '.exe' if os.name == 'nt' else ''
    kwargs = {'blastn_path': Path(binary_dir) / ('blastn' + suffix),
              'makeblastdb_path': Path(binary_dir) / ('makeblastdb' + suffix),
              'cache_dir': tmp_path / 'cache'}
    result = call_cgassembly(path, schema, **kwargs)
    assert result['alleles']['abc'] == 'NOVEL_' + hashlib.sha256(sequence.encode()).hexdigest()
    assert result['st'] is None
    assert result['calls'][0]['cds_qc']['valid'] is True
    assert result['novel_sequences'][0]['sequence'] == sequence
    copy_path = tmp_path / 'renamed.fa'
    copy_path.write_bytes(path.read_bytes())
    cached = call_cgassembly(copy_path, schema, **kwargs)
    assert cached['sample_name'] == 'renamed'
    assert cached['input_path'] == str(copy_path)
    assert cached['alleles'] == result['alleles']
    reverse_path = tmp_path / 'reverse.fa'
    reverse_path.write_text('>contig\n' + reverse_complement('TAA' * 100 + sequence + 'TAA' * 100) + '\n')
    reverse_result = call_cgassembly(reverse_path, schema, **kwargs)
    assert reverse_result['alleles'] == result['alleles']
    assert reverse_result['calls'][0]['hits'][0]['strand'] == '-'
