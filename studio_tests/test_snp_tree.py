"""SKA2 SNP cohorts: their own scale, their own denominators, their own refusals.

The synthetic genomes here are written into the test's own directory and differ
at positions this file chooses, so every distance asserted below is a number
that can be counted by hand. They are a wiring control, not an epidemiological
cohort. Tests that need the native engine assert the honest behaviour whether or
not SKA2 is staged on the machine running them.
"""

import json
import random
from pathlib import Path

import pytest

from wmlstudio.sequence import file_sha256
from wmlstudio.ska_runtime import runtime_capabilities
from wmlstudio.snp_tree import (
    DRAWN_TREE,
    KIND,
    LIMITATIONS,
    NO_LINK_THRESHOLD,
    PUBLISHED_PROTOCOLS,
    SEPARATION,
    alignment_handoff,
    at_minimum_shared_fraction,
    build_snp_tree,
    protocol,
    snp_forest,
    snp_pairs,
    snp_payload,
    snp_scale,
    threshold_binding,
)
from wmlstudio.widgets import TreeView

SITES = [400 + 300 * index for index in range(15)]
SWAP = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}


def genomes():
    """One 20 kb sequence, two isolates differing from it at chosen sites, one fragment."""
    base = ''.join(random.Random(11).choices('ACGT', k=20000))

    def mutate(positions):
        bases = list(base)
        for position in positions:
            bases[position] = SWAP[bases[position]]
        return ''.join(bases)

    # iso-d is the first 8 kb only: a fragmentary assembly that carries no
    # difference at all, and must still never be called close to anything.
    return {'iso-a': base, 'iso-b': mutate(SITES[:3]), 'iso-c': mutate(SITES[3:]), 'iso-d': base[:8000]}


def cohort(tmp_path):
    samples = []
    for name, sequence in genomes().items():
        path = tmp_path / f'{name}.fasta'
        path.write_text(f'>contig1\n{sequence}\n', encoding='ascii')
        samples.append({'id': name, 'name': name.upper(), 'input_path': str(path),
                        'input_sha256': file_sha256(path)})
    return samples


def fake_run(rows=None, *, inputs=('iso-a', 'iso-b'), k=31, engine='SKA2',
             method='reference-free-assembly-split-kmer-SNPs'):
    """A SKA2 result of the exact shape run_ska returns, without running anything."""
    held = {name: 1000 + index for index, name in enumerate(inputs)}
    return {'format_version': 1, 'status': 'completed', 'engine': engine, 'version': '0.5.1',
            'run_id': 'ska2-test', 'method': method,
            'parameters': {'k': k, 'threads': 1, 'ambiguous_bases': 'excluded',
                           'min_frequency': 0.0, 'minimum_shared_fraction': 0.95},
            'binary_sha256': 'a' * 64,
            'inputs': [{'sample_id': name, 'sample_name': name.upper(), 'input_path': f'/read-only/{name}.fa',
                        'input_sha256': f'{index}' * 64} for index, name in enumerate(inputs)],
            'rows': list(rows if rows is not None else []),
            'split_kmers': {'k': k, 'cohort_split_kmers': sum(held.values()), 'per_sample': held},
            'alignment': {'status': 'completed', 'columns': 20, 'rows': [],
                          'limitations': ['Cohort alignment columns depend on this cohort.']},
            'mapped': None, 'output_directory': '/read-only/run', 'threshold': None,
            'limitations': ['SNP similarity is not direct transmission or direction of spread.']}


def pair(source, target, distance, *, shared=1000, fraction=0.99, comparable=True, observed=None):
    return {'source': source, 'target': target, 'source_name': source.upper(), 'target_name': target.upper(),
            'distance': distance if comparable else None,
            'observed_snp_count': distance if observed is None else observed,
            'shared_split_kmers': shared, 'unshared_split_kmers': 10, 'shared_fraction': fraction,
            'source_split_kmers': 1000, 'target_split_kmers': 1001,
            'shared_fraction_of_smaller': fraction, 'comparable': comparable,
            'reason': '' if comparable else 'Insufficient shared unambiguous split-kmers; no accepted distance.'}


def test_the_snp_scale_names_its_own_quantity_and_borrows_no_allele_words():
    scale = snp_scale(fake_run(), cohort='Ward B review')
    assert scale['kind'] == KIND == 'snp' and scale['difference_word'] == 'SNPs'
    assert scale['target_word'] == 'shared split k-mers' and scale['distance_phrase'] == 'SNP-distance'
    assert 'allele' not in scale['unit'] and 'allele' not in scale['caption']
    assert scale['scheme'] == 'ska2:assembly-split-kmer-k31'
    # No single cohort target count exists: every pair has its own denominator.
    assert scale['targets'] == 0 and 'own denominator' in scale['denominator_note']
    assert 'Ward B review' in scale['caption'] and 'k=31' in scale['caption']
    assert scale['separation'] == SEPARATION
    for word in ('MLST', 'cgMLST', 'no threshold'):
        assert word in SEPARATION


def test_a_reference_mapped_or_foreign_result_is_never_drawn_on_the_snp_tree_scale():
    with pytest.raises(ValueError, match='reference-free SKA2 assembly distances'):
        protocol(fake_run(method='reference-mapped-filtered-SNPs'))
    with pytest.raises(ValueError, match='reference-free SKA2 assembly distances'):
        snp_forest(fake_run(engine='some other engine'))


def test_pairs_sharing_too_little_sequence_have_no_distance_no_edge_and_stay_visible():
    rows = [pair('iso-a', 'iso-b', 3), pair('iso-a', 'iso-c', 12),
            pair('iso-b', 'iso-c', 15), pair('iso-a', 'iso-d', None, shared=40,
                                             fraction=0.39, comparable=False, observed=0),
            pair('iso-b', 'iso-d', None, shared=40, fraction=0.39, comparable=False, observed=3),
            pair('iso-c', 'iso-d', None, shared=40, fraction=0.39, comparable=False, observed=12)]
    payload = snp_payload(fake_run(rows, inputs=('iso-a', 'iso-b', 'iso-c', 'iso-d')))
    excluded = {(row['source'], row['target']) for row in payload['not_comparable']}
    assert excluded == {('iso-a', 'iso-d'), ('iso-b', 'iso-d'), ('iso-c', 'iso-d')}
    assert all(row['distance'] is None for row in payload['not_comparable'])
    # The zero that was seen is kept, and is never promoted to a distance.
    assert [row['observed_snp_count'] for row in payload['not_comparable']] == [0, 3, 12]
    assert all('iso-d' not in (edge['source'], edge['target']) for edge in payload['graph']['edges'])
    assert payload['summary']['isolates_without_a_comparable_pair'] == ['iso-d']
    assert payload['summary'] == {**payload['summary'], 'comparable_pairs': 3, 'excluded_pairs': 3,
                                  'minimum_snps': 3, 'median_snps': 12, 'maximum_snps': 15, 'edges': 2}
    matrix = payload['matrix']
    position = {row['sample_id']: index for index, row in enumerate(matrix['samples'])}
    assert matrix['distance'][position['iso-a']][position['iso-d']] is None
    assert matrix['observed_snp_count'][position['iso-a']][position['iso-d']] == 0
    assert matrix['comparable'][position['iso-a']][position['iso-d']] is False
    assert 'never means zero differences' in matrix['missing']


def test_every_distance_carries_the_sequence_it_was_measured_over():
    rows = [pair('iso-a', 'iso-b', 3, shared=19880, fraction=0.9906)]
    payload = snp_payload(fake_run(rows))
    row = payload['pairs'][0]
    assert '19,880 split k-mers shared' in row['denominator_label']
    assert "99.1% of the pair's combined set" in row['denominator_label']
    assert row['unit'] == 'SNPs' and row['shared_split_kmers'] == 19880
    edge = payload['graph']['edges'][0]
    assert edge['denominator_label'] == row['denominator_label'] and edge['distance'] == 3
    matrix = payload['matrix']
    assert matrix['shared_split_kmers'][0][1] == matrix['shared_split_kmers'][1][0] == 19880
    assert matrix['distance'][0][0] == 0 and 'by definition' in matrix['diagonal']
    assert matrix['shared_split_kmers'][0][0] == payload['split_kmers']['per_sample']['iso-a']
    # A pair that shared all but a few hundred k-mers is never rounded up to 100%.
    nearly = snp_payload(fake_run([pair('iso-a', 'iso-b', 1, shared=2_999_820, fraction=0.9999399)]))
    assert 'just under 100% of the pair' in nearly['pairs'][0]['denominator_label']
    whole = snp_payload(fake_run([pair('iso-a', 'iso-b', 1, shared=2_999_820, fraction=1.0)]))
    assert '100% of the pair' in whole['pairs'][0]['denominator_label']


def test_cohort_alignment_columns_are_kept_apart_from_the_pairwise_distance():
    run = fake_run([pair('iso-a', 'iso-b', 3), pair('iso-a', 'iso-c', None, fraction=0.4, comparable=False,
                                                    observed=1)], inputs=('iso-a', 'iso-b', 'iso-c'))
    run['alignment']['rows'] = [
        {'source': 'iso-a', 'target': 'iso-b', 'comparable_columns': 15, 'differing_columns': 3,
         'alignment_columns': 15, 'comparable_fraction': 1.0},
        {'source': 'iso-a', 'target': 'iso-c', 'comparable_columns': 15, 'differing_columns': 1,
         'alignment_columns': 15, 'comparable_fraction': 1.0}]
    rows = {(row['source'], row['target']): row for row in snp_pairs(run)}
    accepted, refused = rows[('iso-a', 'iso-b')], rows[('iso-a', 'iso-c')]
    # The alignment number never appears in a distance column, only under its own key.
    assert accepted['cohort_alignment']['differing_columns'] == 3 and accepted['distance'] == 3
    assert refused['distance'] is None and refused['cohort_alignment']['differing_columns'] == 1
    assert refused['cohort_alignment']['pairwise_distance_accepted'] is False
    assert 'must not be read as one' in refused['cohort_alignment']['note']
    assert all(row['cohort_alignment']['cohort_dependent'] for row in rows.values())
    limitations = snp_payload(run)['limitations']
    assert 'Cohort alignment columns depend on this cohort.' in limitations
    assert any('Recombination is not detected or removed' in line for line in limitations)


def test_a_cohort_where_nothing_is_comparable_says_so_instead_of_drawing_nothing():
    """Real E. faecium assemblies share 0.53-0.78 of their combined split k-mers."""
    rows = [pair('iso-a', 'iso-b', None, shared=2344102, fraction=0.662, comparable=False, observed=2356),
            pair('iso-a', 'iso-c', None, shared=1928027, fraction=0.567, comparable=False, observed=7510),
            pair('iso-b', 'iso-c', None, shared=1911287, fraction=0.526, comparable=False, observed=7922)]
    for row, value in zip(rows, (0.827, 0.772, 0.766), strict=True):
        row['shared_fraction_of_smaller'] = value
    payload = snp_payload(fake_run(rows, inputs=('iso-a', 'iso-b', 'iso-c')))
    state = payload['comparability']
    assert state['status'] == 'all_pairs_excluded' and state['accepted_pairs'] == 0
    assert state['minimum_shared_fraction'] == 0.95 and state['basis'] == 'combined'
    assert state['observed_combined_fraction'] == [0.526, 0.662]
    assert state['observed_fraction_of_smaller'] == [0.766, 0.827]
    assert 'no tree was drawn' in state['message'] and 'large accessory genome' in state['message']
    assert 'at_minimum_shared_fraction' in state['how_to_change_it']
    assert payload['graph']['edges'] == []
    relaxed = snp_payload(at_minimum_shared_fraction(fake_run(rows, inputs=('iso-a', 'iso-b', 'iso-c')), 0.5))
    assert relaxed['comparability']['status'] == 'all_pairs_comparable'
    assert relaxed['comparability_reread']['to'] == 0.5 and len(relaxed['graph']['edges']) == 2
    # The two denominators are different questions and the answer records which was asked.
    by_smaller = snp_payload(at_minimum_shared_fraction(
        fake_run(rows, inputs=('iso-a', 'iso-b', 'iso-c')), 0.8, basis='smaller'))
    assert by_smaller['comparability']['basis'] == 'smaller'
    assert by_smaller['comparability']['accepted_pairs'] == 1
    assert 'smaller isolate' in by_smaller['comparability']['basis_words']
    assert by_smaller['comparability_reread'] == {**by_smaller['comparability_reread'],
                                                 'from_basis': 'combined', 'to_basis': 'smaller'}
    with pytest.raises(ValueError, match="'combined' or 'smaller'"):
        at_minimum_shared_fraction(fake_run(rows), 0.8, basis='core genome')


def test_a_lower_comparability_floor_admits_pairs_without_moving_a_single_snp_count():
    run = fake_run([pair('iso-a', 'iso-b', 3, fraction=0.99),
                    pair('iso-a', 'iso-c', None, shared=700, fraction=0.91, comparable=False, observed=8)],
                   inputs=('iso-a', 'iso-b', 'iso-c'))
    assert snp_payload(run)['summary']['comparable_pairs'] == 1
    relaxed = snp_payload(at_minimum_shared_fraction(run, 0.9))
    admitted = next(row for row in relaxed['pairs'] if row['target'] == 'iso-c')
    assert admitted['comparable'] is True and admitted['distance'] == admitted['observed_snp_count'] == 8
    assert relaxed['summary']['comparable_pairs'] == 2 and relaxed['summary']['excluded_pairs'] == 0
    assert relaxed['protocol']['minimum_shared_fraction'] == 0.9
    reread = at_minimum_shared_fraction(run, 0.9)['comparability_reread']
    assert reread['from'] == 0.95 and reread['to'] == 0.9
    assert 'does not make its distance better evidence' in reread['note']
    # Raising the floor withdraws a distance again rather than keeping a stale one.
    strict = snp_payload(at_minimum_shared_fraction(run, 0.995))
    assert strict['summary']['comparable_pairs'] == 0 and strict['graph']['edges'] == []
    assert run['rows'][0]['comparable'] is True  # the original run is not mutated
    for value in (-0.1, 1.5, True, 'high'):
        with pytest.raises(ValueError, match='between 0 and 1'):
            at_minimum_shared_fraction(run, value)


def test_published_snp_cutoff_is_refused_for_a_protocol_it_was_not_measured_on():
    binding = threshold_binding(fake_run(), 'Enterococcus faecium')
    assert binding['status'] == 'refused_protocol_mismatch'
    assert binding['threshold'] is None and binding['applied'] is False
    assert binding['published_threshold'] == 7
    assert binding['protocol_differences'] == [
        'engine: published SKA, this run SKA2',
        'input material: published short_reads, this run assemblies',
        'split k-mer size k: published 15, this run 31']
    assert 'no number is offered and none was applied' in binding['message']
    # An organism the catalog holds no bound SNP entry for is an evidence gap,
    # never a licence to reuse another organism's number.
    other = threshold_binding(fake_run(), 'Klebsiella pneumoniae')
    assert other['status'] == 'no_curated_cutoff' and other['published_threshold'] is None
    assert 'evidence gap' in other['message'] and other['threshold'] is None
    assert threshold_binding(fake_run(), '')['status'] == 'no_curated_cutoff'


def test_a_matching_protocol_still_only_suggests_and_never_applies(monkeypatch):
    monkeypatch.setitem(PUBLISHED_PROTOCOLS, 'higgs2022:ska-short-reads-k15',
                        {'engine': 'SKA2', 'input': 'assemblies', 'k': 31,
                         'distance': 'reference-free split k-mer SNPs',
                         'description': 'a hypothetical protocol identical to this run'})
    binding = threshold_binding(fake_run(), 'Enterococcus faecium')
    assert binding['status'] == 'suggestion_requires_review' and binding['protocol_differences'] == []
    assert binding['threshold'] is None and binding['applied'] is False
    assert 'not as a setting in use' in binding['message'] and binding['notice']
    monkeypatch.setitem(PUBLISHED_PROTOCOLS, 'higgs2022:ska-short-reads-k15', None)
    monkeypatch.delitem(PUBLISHED_PROTOCOLS, 'higgs2022:ska-short-reads-k15')
    unknown = threshold_binding(fake_run(), 'Enterococcus faecium')
    assert unknown['status'] == 'refused_unknown_protocol' and unknown['threshold'] is None


def test_a_reader_chosen_link_threshold_is_recorded_as_a_choice_not_a_cutoff():
    payload = snp_payload(fake_run([pair('iso-a', 'iso-b', 3)]), organism='Enterococcus faecium',
                          link_threshold=5)
    binding = payload['threshold']
    assert binding['link_threshold'] == 5 and binding['link_threshold_source'] == 'reader_selected_unvalidated'
    assert 'single linkage' in binding['link_threshold_warning'] and binding['threshold'] is None
    assert payload['graph']['cluster_threshold'] == 5
    # Nothing chosen means nothing grouped: every edge is drawn as unlinked.
    default = snp_payload(fake_run([pair('iso-a', 'iso-b', 3)]))
    assert default['graph']['cluster_threshold'] == NO_LINK_THRESHOLD == -1
    assert default['threshold']['link_threshold_source'] == 'none'
    for value in (-1, 2.5, True, '3'):
        with pytest.raises(ValueError, match='nonnegative integer'):
            snp_forest(fake_run(), link_threshold=value)


def test_typing_evidence_never_rides_into_a_snp_forest_record():
    supplied = {'iso-a': {'sample_name': 'Ward B isolate', 'st': '80', 'metadata': {'ward': 'B'},
                          'alleles': {'locus1': '1'}, 'calls': [{'locus': 'locus1', 'status': 'exact'}],
                          'scheme': 'efaecium', 'scheme_digest': 'digest'},
                'not-in-this-cohort': {'sample_name': 'someone else'}}
    records = snp_forest(fake_run(), records=supplied)['results']
    first = next(row for row in records if row['sample_id'] == 'iso-a')
    assert first['sample_name'] == 'Ward B isolate' and first['st'] == '80'
    assert first['metadata'] == {'ward': 'B'}
    assert not {'alleles', 'calls', 'scheme', 'scheme_digest'} & set(first)
    assert [row['sample_id'] for row in records] == ['iso-a', 'iso-b']
    assert supplied['iso-a']['alleles'] == {'locus1': '1'}  # the caller's record is untouched


def test_forest_contents_have_the_shape_a_tree_view_already_draws(qtbot):
    view = TreeView()
    qtbot.addWidget(view)
    forest = snp_forest(fake_run([pair('iso-a', 'iso-b', 3)]), cohort='Ward B review')
    assert set(forest) == set(view.graph_contents())
    view.set_scale(forest['scale'])
    view.draw_results(forest['results'], [], NO_LINK_THRESHOLD)
    assert sorted(view._results) == ['iso-a', 'iso-b']
    assert view.scale_caption() == forest['scale']['caption']
    assert 'allele' not in view.scale_caption()


def test_native_ska2_cohort_measures_the_differences_that_were_written(tmp_path):
    # The organism is the one the catalog holds a SNP entry for, so the run also
    # exercises the refusal of a cutoff measured on another protocol.
    organism, samples = 'Enterococcus faecium', cohort(tmp_path)
    if not runtime_capabilities()['available']:
        with pytest.raises(ValueError, match='not staged'):
            build_snp_tree(samples, tmp_path / 'out', organism=organism, threads=1)
        return
    payload = build_snp_tree(samples, tmp_path / 'out', organism=organism, threads=1, k=31)
    distances = {(row['source'], row['target']): row['distance'] for row in payload['pairs']}
    assert distances[('iso-a', 'iso-b')] == 3 and distances[('iso-a', 'iso-c')] == 12
    assert distances[('iso-b', 'iso-c')] == 15
    # The fragment differs nowhere, yet is comparable to nothing: a short overlap
    # is an unknown distance, not a small one.
    assert [distances[key] for key in distances if 'iso-d' in key] == [None, None, None]
    fragment = next(row for row in payload['pairs'] if row['target'] == 'iso-d' and row['source'] == 'iso-a')
    assert fragment['observed_snp_count'] == 0 and fragment['shared_fraction'] < 0.5
    assert [(edge['source'], edge['target'], edge['distance']) for edge in payload['graph']['edges']] == [
        ('iso-a', 'iso-b', 3), ('iso-a', 'iso-c', 12)]
    assert payload['alignment']['status'] == 'completed' and payload['alignment']['columns'] == 15
    assert payload['split_kmers']['per_sample']['iso-d'] < payload['split_kmers']['per_sample']['iso-a']
    assert payload['threshold']['status'] == 'refused_protocol_mismatch'
    written = json.loads((tmp_path / 'out' / payload['run_id'] / 'snp-tree.json').read_text(encoding='utf-8'))
    assert written['summary'] == payload['summary'] and written['limitations'][:len(LIMITATIONS)] == list(LIMITATIONS)


def test_native_ska2_alignment_is_refused_rather_than_truncated(tmp_path):
    samples = cohort(tmp_path)
    if not runtime_capabilities()['available']:
        pytest.skip('Native SKA2 is not staged; the refusal path needs a real alignment to refuse.')
    payload = build_snp_tree(samples, tmp_path / 'out', threads=1, max_alignment_columns=2)
    assert payload['alignment']['status'] == 'refused' and payload['alignment']['columns'] == 15
    assert 'no truncated alignment' in payload['alignment']['reason'].lower()
    # Refusing the cohort alignment never touches the pairwise distances.
    assert payload['summary']['comparable_pairs'] == 3
    assert all(row['cohort_alignment'] is None for row in payload['pairs'])


def test_the_cohort_alignment_is_handed_over_by_name_and_is_never_called_a_phylogeny():
    """A request for "an ML tree" must not be answered by renaming the drawn picture.

    The alignment is the only file in a SKA2 run a maximum-likelihood tree can be
    inferred from, and before this it had no path anywhere in the interface. This
    checks that the path is offered, that the programs that would use it are
    named, and that nothing in the hand-off claims a tree was inferred here.
    """
    run = fake_run([pair('iso-a', 'iso-b', 3)])
    run['alignment'] = dict(run['alignment'], file='alignment.fasta', min_freq=0.9,
                            sha256='b' * 64,
                            method='cohort-variable-site-split-kmer-alignment')
    handoff = alignment_handoff(run)
    assert handoff['status'] == 'completed'
    assert handoff['path'] == str(Path('/read-only/run') / 'alignment.fasta')
    assert handoff['columns'] == 20 and handoff['minimum_kmer_frequency'] == 0.9
    assert handoff['tree_builders'] == ['IQ-TREE', 'FastTree', 'RAxML-NG']
    assert handoff['tree_built_here'] is False
    assert 'not a maximum-likelihood phylogeny' in handoff['drawn_tree']
    assert 'no branch lengths were estimated' in handoff['drawn_tree']
    assert 'none of them is part of this application' in handoff['maximum_likelihood_route']
    # The alignment column count is not the pairwise denominator and says so.
    assert 'not the pairwise shared split k-mer denominator' in handoff['column_meaning']
    payload = snp_payload(run)
    assert payload['alignment_handoff'] == handoff and payload['drawn_tree'] == DRAWN_TREE
    assert any('not a maximum-likelihood tree' in line for line in payload['limitations'])


def test_an_alignment_that_was_refused_or_never_run_offers_no_path_to_build_a_tree_from():
    """An unwritten alignment must read as absent, never as a file waiting on disk.

    A refused alignment is one this run declined to describe as a cohort
    alignment; offering its path would invite a tree to be built from it anyway.
    """
    refused = fake_run([pair('iso-a', 'iso-b', 3)])
    refused['alignment'] = dict(refused['alignment'], status='refused', file='alignment.fasta',
                                reason='The cohort alignment is 300,000,000 bytes, above the bound.')
    handoff = alignment_handoff(refused)
    assert handoff['path'] is None and handoff['status'] == 'refused'
    assert 'refused by this run' in handoff['message'] and 'above the bound' in handoff['message']
    absent = fake_run([pair('iso-a', 'iso-b', 3)])
    absent['alignment'] = {'status': 'not_run', 'rows': [],
                           'reason': 'No cohort alignment was requested.'}
    missing = alignment_handoff(absent)
    assert missing['path'] is None and missing['file'] is None
    assert 'nothing to take to a tree builder' in missing['message']
    assert 'No cohort alignment was requested.' in missing['message']
    # A run that recorded nothing at all is unrecorded, not empty.
    silent = fake_run([pair('iso-a', 'iso-b', 3)])
    silent['alignment'] = {}
    assert alignment_handoff(silent)['status'] == 'unknown'
    assert 'not an absent one' in alignment_handoff(silent)['message']
