import hashlib
import os
from pathlib import Path

import pytest

from wmlstudio.adhoc import _select_cohort_loci, create_adhoc_scheme
from wmlstudio.cgtyping import full_cds_qc
from wmlstudio.sequence import AnalysisCancelled
from wmlstudio.typing import call_assembly, load_scheme


def gene(identifier, sequence='ATGAAATAA', valid=True):
    return {'id': identifier, 'sequence': sequence,
            'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'cds_qc': {**full_cds_qc(sequence), 'valid': valid}}


def test_cohort_selection_requires_prevalence_and_keeps_sequence_variants():
    seeds = [gene('a'), gene('b', 'ATGCCCTAA')]
    cohort = [seeds, [gene('x', 'ATGAAGTAA')], [gene('y')]]
    relations = [{'a': {'a'}, 'b': {'b'}}, {'a': {'x'}}, {'a': {'y'}}]
    retained, excluded = _select_cohort_loci(seeds, cohort, relations, 1.0)
    assert list(retained) == ['a']
    assert set(retained['a']['variants'].values()) == {'ATGAAATAA', 'ATGAAGTAA'}
    assert retained['a']['present_in'] == [0, 1, 2]
    assert '1/3' in excluded['b']


def test_paralog_in_one_sample_excludes_locus_even_at_low_prevalence():
    seed = gene('a')
    retained, excluded = _select_cohort_loci([seed], [[seed], [gene('x'), gene('y')]],
                                            [{'a': {'a'}}, {'a': {'x', 'y'}}], .5)
    assert retained == {}
    assert 'multiple CDS' in excluded['a']


def test_cross_locus_shared_cds_excludes_both_loci():
    seeds = [gene('a'), gene('b')]
    retained, excluded = _select_cohort_loci(seeds, [seeds], [{'a': {'a'}, 'b': {'a'}}], 1.)
    assert retained == {}
    assert set(excluded) == {'a', 'b'}


def test_partial_cds_does_not_meet_prevalence():
    seed = gene('a')
    retained, excluded = _select_cohort_loci([seed], [[seed], [gene('x', valid=False)]],
                                            [{'a': {'a'}}, {'a': {'x'}}], 1.)
    assert retained == {}
    assert '1/2' in excluded['a']


def test_builder_cancel_before_tools_or_io(tmp_path):
    with pytest.raises(AnalysisCancelled):
        create_adhoc_scheme([tmp_path / 'missing.fa'], tmp_path / 'output', 'name', cancelled=lambda: True)
    assert not (tmp_path / 'output').exists()


def test_native_adhoc_builder_is_atomic_reproducible_and_callable(tmp_path):
    binary_dir = os.environ.get('WMLSTUDIO_TEST_BLAST_DIR')
    if not binary_dir:
        pytest.skip('Set WMLSTUDIO_TEST_BLAST_DIR for native BLAST integration.')
    import random
    rng = random.Random(12)
    sequence = 'ATG' + ''.join(rng.choice(['GCT', 'GCA', 'AAA', 'GAA', 'GAT', 'CTG', 'ATT', 'GGT', 'TTC', 'TAT'])
                                for _ in range(400)) + 'TAA'
    first = tmp_path / 'first.fa'
    first.write_text('>contig\n' + 'TAA' * 100 + sequence + 'TAA' * 100 + '\n')
    variant = sequence[:300] + ('GCT' if sequence[300:303] != 'GCT' else 'GCA') + sequence[303:]
    second = tmp_path / 'second.fa'
    second.write_text('>contig\n' + 'TAA' * 100 + variant + 'TAA' * 100 + '\n')
    suffix = '.exe' if os.name == 'nt' else ''
    kwargs = {'blastn_path': Path(binary_dir) / ('blastn' + suffix),
              'makeblastdb_path': Path(binary_dir) / ('makeblastdb' + suffix)}
    output = tmp_path / 'library'
    result = create_adhoc_scheme([first, second], output, 'Local truth', **kwargs)
    assert result['locus_count'] == 1
    schema = load_scheme(result['path'])
    assert schema.metadata['provider'] == 'local-ad-hoc'
    assert schema.allele_count == 2
    assert call_assembly(first, schema)['alleles'][schema.loci[0]] == 'SHA256_' + hashlib.sha256(sequence.encode()).hexdigest()
    assert call_assembly(second, schema)['alleles'][schema.loci[0]] == 'SHA256_' + hashlib.sha256(variant.encode()).hexdigest()
    repeated = create_adhoc_scheme([first, second], output, 'Local truth', **kwargs)
    assert repeated['created'] is False
    assert repeated['scheme_digest'] == result['scheme_digest']
    assert len(list(output.iterdir())) == 1
    original = first.read_bytes()
    with pytest.raises(ValueError, match='Duplicate assembly contents'):
        clone = tmp_path / 'clone.fa'
        clone.write_bytes(original)
        create_adhoc_scheme([first, clone], output, 'Invalid duplicate cohort', **kwargs)
    assert first.read_bytes() == original
