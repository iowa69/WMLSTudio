"""A SNP cohort as a drawable forest, on a scale that belongs to nothing else.

SKA2 split k-mer SNPs are reference-free: a distance counts the unambiguous
middle-base differences between the split k-mers two assemblies actually share.
What is shared is the denominator, and it differs for every pair, so it travels
in the same row as the distance it qualifies. A pair compared over a small part
of its sequence has an *unknown* distance, never a small one, and gets no edge.

Three quantities are kept apart here and never merged, summed, ranked together
or drawn on one axis: a SNP distance, a classical seven-locus allele distance,
and a cgMLST target distance. Only the edge-selection algorithm is shared with
:mod:`wmlstudio.comparison`; no scale, column or threshold is.

The forest this module builds is a minimum spanning tree of pairwise distances.
No phylogeny is inferred here or anywhere else in this application, so nothing in
this module may be called a maximum-likelihood tree. What it does offer is the
file such an inference needs: :func:`alignment_handoff` names the cohort
alignment SKA2 wrote so a reader can take it to a tree builder themselves.

A published SNP cutoff is organism- and protocol-specific. The one SNP entry the
threshold catalog holds was measured with the original SKA on short reads at
k=15, so this module compares protocols before it offers a number, and refuses
when engine, input type or k differ. Refusal is the expected answer, not a bug.
"""

from __future__ import annotations

import json
import statistics
from datetime import datetime, timezone
from pathlib import Path

from .comparison import forest_from_distances
from .sequence import check_cancelled
from .threshold_guidance import SUGGESTION_NOTICE, short_citation, suggested_threshold

#: The typing kind this forest carries into a view, a window title and an export.
KIND = 'snp'
TITLE = 'SKA2 split k-mer SNPs'
#: The words a view must use for this forest's denominator and its differences.
#: "loci" and "allele differences" are another quantity's words and never apply.
TARGET_WORD = 'shared split k-mers'
DIFFERENCE_WORD = 'SNPs'
DISTANCE_PHRASE = 'SNP-distance'
QUESTION = 'How many single-nucleotide differences separate these isolates in the sequence they share?'

SEPARATION = (
    'A SNP distance, a classical MLST allele distance and a cgMLST target distance are three '
    'different quantities. They share no scale, no axis, no column and no threshold, and a "5" on '
    'one of them is not a "5" on another.')
INTERPRETATION = (
    'Small SNP distances are a research signal for epidemiological review, not proof of direct '
    'transmission, of its direction, or of a common source. Unshared sequence is unknown, not '
    'identical. A published SNP cutoff belongs to the organism and the protocol it was measured '
    'on and is not validated for this run merely because the species name matches.')

#: No link threshold is claimed until one is justified. Every SNP distance is a
#: nonnegative integer, so a negative threshold groups nothing and draws every
#: edge as unlinked, which is the honest picture when no cutoff has been bound.
NO_LINK_THRESHOLD = -1

LIMITATIONS = (
    'Split k-mer SNPs are reference-free: they count unambiguous middle-base differences between '
    'the split k-mers two assemblies actually share, not differences across whole genomes.',
    'Repetitive and accessory sequence is excluded, because a split k-mer occurring more than once '
    'cannot be placed. A SNP distance describes shared, unique sequence only.',
    'Recombination is not detected or removed by this analysis. A single imported fragment can '
    'contribute more SNPs than years of vertical evolution, and nothing here separates the two.',
    'A SNP threshold is organism- and protocol-specific. It does not survive a change of engine, of '
    'input type, of k, or of masking, and no cutoff is applied by this module.',
    'A minimum spanning forest is a layout of pairwise SNP distances. It is not a phylogeny, not a '
    'time line and not a transmission chain, and line length carries no meaning.',
    'A pair that shares too little sequence has an unknown distance, never zero. Such pairs carry '
    'no edge and are listed separately with the fraction they did share.',
    'Assembly choices propagate: sequence missing from an assembly is missing from its split k-mers, '
    'and an assembly artefact is indistinguishable here from a real difference.',
    'SNP distances are never merged with, added to, or plotted against classical MLST or cgMLST '
    'allele differences.',
    'No phylogeny is inferred anywhere in this application. The picture is a minimum spanning tree '
    'of the pairwise distances, not a maximum-likelihood tree: no substitution model was fitted, no '
    'branch lengths were estimated and no support values were computed.')

#: What the drawn picture is, and the thing readers most often mistake it for.
#: Printed in the view and carried in the payload, because a picture separates
#: from its caption and this is the sentence that must travel with it.
DRAWN_TREE = (
    'What is drawn here is a minimum spanning tree of the pairwise SNP distances in the table above. '
    'It is not a maximum-likelihood phylogeny and no phylogeny was inferred: no substitution model '
    'was fitted, no branch lengths were estimated, no support values were computed, and no node on '
    'the screen is an ancestor of any other. The length of a line carries no meaning at all.')

#: The programs that do build a maximum-likelihood tree. None of them is bundled
#: with this application, and naming one is not the same as having run it.
TREE_BUILDERS = ('IQ-TREE', 'FastTree', 'RAxML-NG')

ML_TREE_ROUTE = (
    'A maximum-likelihood tree is inferred from a sequence alignment by a phylogenetics program — '
    + ', '.join(TREE_BUILDERS[:-1]) + ' or ' + TREE_BUILDERS[-1] + ' — and none of them is part of '
    'this application. The cohort alignment SKA2 wrote is named here so it can be taken to one. The '
    'tree such a program returns is a different object from this picture, with branch lengths that '
    'do mean something, and nothing here will redraw it as this picture or this picture as it.')

#: What each published SNP cutoff in :mod:`wmlstudio.threshold_guidance` was
#: actually measured on, keyed by the catalog's own ``scheme_key``. A number is
#: only offered when every field here matches the run that is asking for it.
PUBLISHED_PROTOCOLS = {
    'higgs2022:ska-short-reads-k15': {
        'engine': 'SKA', 'input': 'short_reads', 'k': 15,
        'distance': 'reference-free split k-mer SNPs',
        'description': 'Original SKA (version 1, not SKA2) built directly from short reads at k=15, '
                       'with that study\'s own read quality and coverage filters.'},
}

_RUN_METHOD = 'reference-free-assembly-split-kmer-SNPs'
# A caller may hand in richer display records than the SKA2 run knows about, but
# typing evidence must not ride into a SNP forest: a view that merges nodes with
# identical allele profiles would merge two isolates whose SNP distance is not
# zero, which is exactly the conflation this module exists to prevent.
_TYPING_KEYS = ('alleles', 'calls', 'scheme', 'scheme_digest', 'st_evidence', 'profile')


def protocol(result):
    """The exact protocol one SKA2 run measured with, as a key a cutoff can bind to.

    Raises when handed anything but reference-free SKA2 assembly distances: a
    reference-mapped or coding-only run is separate evidence with its own
    denominators, and must not be drawn on this scale.
    """
    method = str(result.get('method') or '')
    engine = str(result.get('engine') or '')
    parameters = dict(result.get('parameters') or {})
    k = parameters.get('k')
    if engine != 'SKA2' or method != _RUN_METHOD:
        raise ValueError('A SNP tree is built from reference-free SKA2 assembly distances. '
                         f'This result is {engine or "an unnamed engine"} / {method or "an unnamed method"}; '
                         'reference-mapped and filtered distances are separate evidence with their own scale.')
    if isinstance(k, bool) or not isinstance(k, int):
        raise ValueError('The SKA2 result does not record the split k-mer size it used.')
    version = str(result.get('version') or 'unrecorded version')
    return {'key': f'ska2:assembly-split-kmer-k{k}', 'engine': 'SKA2', 'version': version,
            'input': 'assemblies', 'k': k, 'distance': 'reference-free split k-mer SNPs',
            'ambiguous_bases': str(parameters.get('ambiguous_bases') or 'unrecorded'),
            'minimum_shared_fraction': parameters.get('minimum_shared_fraction'),
            'repetitive_and_accessory_sequence': 'excluded by split k-mer uniqueness',
            'recombination_masking': 'not_run',
            'description': f'SKA2 {version} split k-mer SNP distances between assemblies at k={k}, '
                           'ambiguous bases excluded, recombination not masked.'}


def threshold_binding(result, organism='', *, link_threshold=None):
    """Offer a published SNP cutoff only when the protocol it was measured on matches.

    Nothing is ever applied. A number is offered only when organism, engine,
    input type and k all match a curated entry, and even then it is a suggestion
    that still requires local review. ``link_threshold`` is the number a reader
    chose for the picture; it is recorded as their unvalidated choice, kept apart
    from any published value, and it never becomes a finding.
    """
    facts = protocol(result)
    guidance = suggested_threshold(organism, method='snp') if str(organism).strip() else None
    suggestion = (guidance or {}).get('suggestion')
    published = PUBLISHED_PROTOCOLS.get(str((suggestion or {}).get('scheme_key') or ''))
    differences, status = [], 'no_curated_cutoff'
    if suggestion is None:
        message = ('No reviewed SNP cutoff for ' + (str(organism).strip() or 'an unnamed organism')
                   + ' is bound to a protocol in this catalog. This is an evidence gap, not proof that '
                     'no publications exist, and not permission to pick a number.')
    elif published is None:
        status = 'refused_unknown_protocol'
        message = ('The catalog publishes a SNP cutoff under ' + str(suggestion['scheme_key'])
                   + ', but this module holds no description of that protocol, so it cannot check whether '
                     'this run matches it. The number is not offered.')
    else:
        for field, word in (('engine', 'engine'), ('input', 'input material'), ('k', 'split k-mer size k')):
            if facts[field] != published[field]:
                differences.append(f'{word}: published {published[field]}, this run {facts[field]}')
        status = 'refused_protocol_mismatch' if differences else 'suggestion_requires_review'
        message = ((short_citation(suggestion['source']) + ' publishes at most '
                    + str(suggestion['published_threshold']) + ' ' + suggestion['unit'] + ' for '
                    + suggestion['organism'] + ', measured on ' + published['description'])
                   + (' This run differs in ' + '; '.join(differences)
                      + '. A SNP cutoff does not transfer across a changed protocol, so no number is '
                        'offered and none was applied.' if differences else
                      ' The protocol matches, so the number is offered as a suggestion to review, not as '
                      'a setting in use.'))
    link, source, warning = NO_LINK_THRESHOLD, 'none', ''
    if link_threshold is not None:
        if isinstance(link_threshold, bool) or not isinstance(link_threshold, int) or link_threshold < 0:
            raise ValueError('A SNP link threshold must be a nonnegative integer number of SNPs.')
        link, source = int(link_threshold), 'reader_selected_unvalidated'
        warning = (f'Grouping at {link} SNPs is a choice made in this view, not a validated cutoff. '
                   'Groups are single linkage, so a chain of small differences can join two isolates '
                   'whose own distance is larger than this number.')
    return {'method': 'snp', 'organism': str(organism).strip(), 'status': status, 'applied': False,
            'threshold': None, 'protocol': facts, 'published_protocol': published,
            'published_threshold': (suggestion or {}).get('published_threshold'),
            'published_scheme_key': (suggestion or {}).get('scheme_key'),
            'protocol_differences': differences, 'message': message,
            'link_threshold': link, 'link_threshold_source': source, 'link_threshold_warning': warning,
            'notice': SUGGESTION_NOTICE, 'interpretation': INTERPRETATION, 'guidance': guidance}


#: Which denominator a comparability floor is applied to, in the words a reader
#: has to be given before choosing one. Neither is a core-genome fraction.
FRACTION_BASES = {
    'combined': {'field': 'shared_fraction',
                 'words': 'the pair\'s combined split k-mer set',
                 'reading': 'This is SKA2\'s own shared proportion. In a species with a large accessory '
                            'genome two isolates share a modest fraction of their combined sets however '
                            'closely related they are, so a high floor here rejects close pairs for '
                            'carrying accessory sequence rather than for being poorly compared.'},
    'smaller': {'field': 'shared_fraction_of_smaller',
                'words': 'the smaller isolate\'s own split k-mer set',
                'reading': 'This asks how much of the smaller isolate was compared at all, so the other '
                           'isolate\'s accessory content cannot depress it. It is still not the fraction '
                           'of a core genome that was compared.'},
}


def at_minimum_shared_fraction(result, minimum, *, basis='combined'):
    """Re-read one run's pairs at a different comparability floor, recomputing nothing.

    Only which pairs count as comparable changes, and it changes from counts
    SKA2 already reported: no sequence is read again and no SNP count moves.
    Lowering the floor admits pairs compared over less shared sequence, and a
    distance is not more certain for having been admitted, so the floor and the
    denominator it was applied to are both recorded in the returned result and
    travel into every payload built from it.

    ``basis`` names that denominator: ``'combined'`` is the pair's combined split
    k-mer set, which is what :func:`wmlstudio.ska_runtime.run_ska` applied, and
    ``'smaller'`` is the smaller isolate's own set. They answer different
    questions and :data:`FRACTION_BASES` states both.
    """
    protocol(result)
    if basis not in FRACTION_BASES:
        raise ValueError("The comparability floor applies to 'combined' or 'smaller'; no other basis exists.")
    if isinstance(minimum, bool) or not isinstance(minimum, (int, float)) or not 0 <= minimum <= 1:
        raise ValueError('The minimum shared split k-mer fraction must be a number between 0 and 1.')
    minimum, field = float(minimum), FRACTION_BASES[basis]['field']
    rows = []
    for row in result.get('rows') or []:
        value = row.get(field)
        if value is None:
            raise ValueError(f'This run does not record {FRACTION_BASES[basis]["words"]}, so a floor cannot '
                             'be applied to it. Run the cohort again to record both denominators.')
        comparable = row['shared_split_kmers'] > 0 and value >= minimum
        rows.append(dict(row, comparable=comparable,
                         distance=row['observed_snp_count'] if comparable else None,
                         reason='' if comparable else
                                'Insufficient shared unambiguous split-kmers; no accepted distance.'))
    parameters = dict(result.get('parameters') or {})
    return dict(result, rows=rows,
                parameters={**parameters, 'minimum_shared_fraction': minimum,
                            'minimum_shared_fraction_basis': basis},
                comparability_reread={
                    'from': parameters.get('minimum_shared_fraction'),
                    'from_basis': parameters.get('minimum_shared_fraction_basis', 'combined'),
                    'to': minimum, 'to_basis': basis,
                    'source': 'SKA2 shared and unshared split k-mer counts already recorded by this run',
                    'note': 'A comparability floor decides which pairs may carry a distance. It does not '
                            'change any SNP count, and admitting a pair compared over less sequence does '
                            'not make its distance better evidence.'})


def snp_scale(result, *, cohort='', profiles=None):
    """Name the quantity this forest measures, in the mapping a graph view reads.

    ``targets`` is deliberately zero: there is no single cohort target count to
    put beside a SNP distance, because every pair is compared over the split
    k-mers that pair shares. The per-pair denominator travels on the edge.
    """
    facts = protocol(result)
    inventory = dict(result.get('split_kmers') or {})
    count = len(result.get('inputs') or ()) if profiles is None else int(profiles)
    caption = ' · '.join(filter(None, [
        TITLE, f'reference-free, k={facts["k"]}',
        f'{count} isolates' if count else '', cohort]))
    return {'kind': KIND, 'title': TITLE, 'target_word': TARGET_WORD,
            'difference_word': DIFFERENCE_WORD, 'distance_phrase': DISTANCE_PHRASE, 'question': QUESTION,
            'scheme': facts['key'], 'scheme_digest': str(result.get('binary_sha256') or ''),
            'targets': 0, 'profiles': count,
            'cohort_split_kmers': inventory.get('cohort_split_kmers'),
            'unit': 'SNPs over the split k-mers each pair shares',
            'caption': caption, 'separation': SEPARATION,
            'denominator_note': 'Every pair carries its own denominator: the split k-mers those two '
                                'isolates share. There is no single cohort target count to divide by.'}


def _percent(value):
    """One decimal, and never a partial overlap rounded up to a flat 100%."""
    if value >= 1:
        return '100%'
    text = f'{value:.1%}'
    return 'just under 100%' if text == '100.0%' else text


def _denominator_label(row):
    """The one line that must appear wherever this edge's distance appears."""
    parts = [f'{row["shared_split_kmers"]:,} split k-mers shared',
             f'{_percent(row["shared_fraction"])} of the pair\'s combined set']
    smaller = row.get('shared_fraction_of_smaller')
    if smaller is not None:
        parts.append(f'{_percent(smaller)} of the smaller isolate\'s')
    return ', '.join(parts)


def snp_pairs(result, *, cancelled=None):
    """Every unordered pair, its SNP distance and the denominators behind that distance.

    The cohort alignment's differing-column count is carried in a nested key of
    its own. It is a cohort-filtered quantity that changes when the cohort
    changes, so it is never flattened into the pairwise distance columns.
    """
    protocol(result)
    alignment = {(str(row['source']), str(row['target'])): dict(row)
                 for row in (dict(result.get('alignment') or {}).get('rows') or [])}
    pairs = []
    for index, row in enumerate(result.get('rows') or []):
        if index % 256 == 0:
            check_cancelled(cancelled)
        pair = dict(row)
        pair['unit'] = DIFFERENCE_WORD
        pair['denominator_label'] = _denominator_label(row)
        columns = alignment.get((str(row['source']), str(row['target'])))
        if columns is not None:
            columns['cohort_dependent'] = True
            columns['pairwise_distance_accepted'] = bool(row['comparable'])
            columns['note'] = ('Cohort-wide variable columns. ' + (
                'The pairwise distance beside them was accepted over the split k-mers this pair shares.'
                if row['comparable'] else
                'The pair did not share enough sequence for a pairwise distance to be accepted, so these '
                'columns are not that pair\'s distance and must not be read as one.'))
        pair['cohort_alignment'] = columns
        pairs.append(pair)
    pairs.sort(key=lambda pair: (str(pair['source']), str(pair['target'])))
    return pairs


def snp_matrix(result, *, cancelled=None):
    """A square SNP matrix whose every cell carries the sequence it was measured over.

    A cell is ``None`` when the pair did not share enough sequence to be
    compared. That is an unknown distance and must never be read, sorted or
    coloured as a zero. The observed count is kept in its own grid so nothing is
    hidden, and the diagonal is zero by definition rather than by measurement.
    """
    records = sorted(({'sample_id': str(entry['sample_id']), 'sample_name': str(entry['sample_name'])}
                      for entry in result.get('inputs') or []), key=lambda entry: entry['sample_id'])
    index = {entry['sample_id']: position for position, entry in enumerate(records)}
    size = len(records)
    held = dict(dict(result.get('split_kmers') or {}).get('per_sample') or {})
    grids = {name: [[None] * size for _ in range(size)]
             for name in ('distance', 'observed_snp_count', 'shared_split_kmers', 'shared_fraction', 'comparable')}
    for position, entry in enumerate(records):
        own = held.get(entry['sample_id'])
        grids['distance'][position][position] = 0
        grids['observed_snp_count'][position][position] = 0
        grids['shared_split_kmers'][position][position] = own
        grids['shared_fraction'][position][position] = 1.0 if own else None
        grids['comparable'][position][position] = True
    for count, row in enumerate(result.get('rows') or []):
        if count % 256 == 0:
            check_cancelled(cancelled)
        left, right = index.get(str(row['source'])), index.get(str(row['target']))
        if left is None or right is None:
            raise ValueError('A SKA2 distance row names an isolate that is not in the cohort.')
        for name, value in (('distance', row['distance']), ('observed_snp_count', row['observed_snp_count']),
                            ('shared_split_kmers', row['shared_split_kmers']),
                            ('shared_fraction', row['shared_fraction']), ('comparable', bool(row['comparable']))):
            grids[name][left][right] = grids[name][right][left] = value
    return {'samples': records, 'unit': DIFFERENCE_WORD, **grids,
            'missing': 'A null distance means the pair did not share enough sequence to be compared. '
                       'It never means zero differences.',
            'observed_meaning': 'The observed count is what SKA2 saw over whatever the pair shared. It is '
                                'reported for every pair, including the pairs whose distance was refused.',
            'diagonal': 'Zero by definition; an isolate is not compared with itself, and the shared count '
                        'on the diagonal is that isolate\'s own split k-mer set.',
            'separation': SEPARATION}


def alignment_handoff(result):
    """Where SKA2 left the cohort alignment, and what can honestly be done with it.

    The alignment is the input a tree-building program needs, so it is named
    rather than left inside a working directory nobody is told about. The file
    itself is not read, hashed or re-checked here: this reports what the run
    recorded, and says plainly when the run recorded nothing.
    """
    protocol(result)
    alignment = dict(result.get('alignment') or {})
    status = str(alignment.get('status') or 'unknown')
    directory, name = str(result.get('output_directory') or ''), str(alignment.get('file') or '')
    path = str(Path(directory) / name) if directory and name else ''
    reason = str(alignment.get('reason') or '')
    if status == 'completed' and path:
        message = ('SKA2 wrote the cohort variable-site alignment to this file. It is the input a '
                   'maximum-likelihood tree is built from, and it is the only file here that one '
                   'can be built from: the distances and the picture are not.')
    elif status == 'completed':
        message = ('This run reports a completed cohort alignment but records no file name for it, '
                   'so no path can be offered. Run the cohort again to have the alignment written '
                   'where it can be found.')
    elif status == 'refused':
        message = ('The cohort alignment was refused by this run, so it is not offered as a '
                   'hand-off. ' + reason)
    elif status == 'not_run':
        message = ('No cohort alignment was written by this run, so there is nothing to take to a '
                   'tree builder. ' + reason)
    else:
        message = ('This run does not record what became of the cohort alignment, so no file is '
                   'offered. An unrecorded alignment is not an absent one.')
    return {'status': status,
            # Only a completed alignment is handed on. A refused one exists on
            # disk, but offering its path would invite a tree to be built from a
            # file this run has already declined to describe as a cohort alignment.
            'path': (path or None) if status == 'completed' else None,
            'file': name or None, 'directory': directory or None,
            'columns': alignment.get('columns'), 'sha256': alignment.get('sha256'),
            'method': alignment.get('method'),
            'minimum_kmer_frequency': alignment.get('min_freq'),
            'reason': reason, 'message': message,
            'drawn_tree': DRAWN_TREE, 'maximum_likelihood_route': ML_TREE_ROUTE,
            'tree_builders': list(TREE_BUILDERS), 'tree_built_here': False,
            'column_meaning': 'Variable columns only, and only the columns this cohort at this '
                              'minimum k-mer frequency produced. The count is not a genome length, '
                              'and it is not the pairwise shared split k-mer denominator.'}


def _graph_records(result, records=None):
    """One display record per isolate: names and recorded detail, never typing evidence."""
    extra = {str(key): dict(value) for key, value in dict(records or {}).items()}
    rows = []
    for entry in sorted(result.get('inputs') or [], key=lambda item: str(item['sample_id'])):
        sample_id = str(entry['sample_id'])
        supplied = extra.get(sample_id, {})
        for key in _TYPING_KEYS:
            supplied.pop(key, None)
        rows.append({**supplied, 'sample_id': sample_id,
                     'sample_name': str(supplied.get('sample_name') or entry['sample_name']),
                     'input_sha256': str(entry.get('input_sha256') or '')})
    return rows


def snp_forest(result, *, records=None, link_threshold=None, cohort='', cancelled=None):
    """Graph contents in the shape a tree view already draws, with a SNP identity.

    The same deterministic Kruskal selection as the allele forests is reused, and
    nothing else: the scale, the unit words and the per-edge denominator all say
    SNPs. Pairs that were not comparable never reach the selection, so no edge is
    ever drawn across sequence the two isolates did not share.
    """
    pairs = snp_pairs(result, cancelled=cancelled)
    threshold = NO_LINK_THRESHOLD
    if link_threshold is not None:
        if isinstance(link_threshold, bool) or not isinstance(link_threshold, int) or link_threshold < 0:
            raise ValueError('A SNP link threshold must be a nonnegative integer number of SNPs.')
        threshold = int(link_threshold)
    return {'results': _graph_records(result, records),
            'edges': [dict(edge) for edge in forest_from_distances(pairs, cancelled=cancelled)],
            'cluster_threshold': threshold, 'groups': [],
            'scale': snp_scale(result, cohort=cohort)}


def _comparability(result, pairs):
    """Say plainly what the comparability floor let through, and what it did not.

    A cohort whose pairs all fall below the floor draws no tree at all, and an
    empty picture with no sentence beside it is the worst outcome this module
    could produce. The observed fractions are reported either way, because the
    floor is a decision a reader has to make with numbers in front of them.
    """
    parameters = dict(result.get('parameters') or {})
    floor = parameters.get('minimum_shared_fraction')
    basis = parameters.get('minimum_shared_fraction_basis', 'combined')
    words = FRACTION_BASES.get(basis, FRACTION_BASES['combined'])
    combined = [pair['shared_fraction'] for pair in pairs]
    smaller = [pair['shared_fraction_of_smaller'] for pair in pairs
               if pair.get('shared_fraction_of_smaller') is not None]
    accepted = sum(bool(pair['comparable']) for pair in pairs)
    status = ('no_pairs' if not pairs else 'all_pairs_excluded' if not accepted
              else 'all_pairs_comparable' if accepted == len(pairs) else 'some_pairs_excluded')
    applied = f'the minimum shared split k-mer fraction of {floor} over {words["words"]}'
    observed = (f'Observed fractions ran from {min(combined):.3f} to {max(combined):.3f} of each pair\'s '
                'combined set' + (f', and from {min(smaller):.3f} to {max(smaller):.3f} of the smaller '
                                  'isolate\'s own set' if smaller else '') + '. ') if combined else ''
    if status == 'all_pairs_excluded':
        message = (f'No pair reached {applied}, so this cohort has no accepted SNP distance and no tree was '
                   f'drawn. {observed}{words["reading"]} Choose the floor and the set it applies to '
                   'deliberately and record why: lowering it admits pairs compared over less shared '
                   'sequence, and a distance is not better evidence for having been admitted.')
    elif status == 'some_pairs_excluded':
        message = (f'{len(pairs) - accepted} of {len(pairs)} pairs fell below {applied} and carry no '
                   f'distance. {observed}{words["reading"]}')
    elif status == 'all_pairs_comparable':
        message = f'Every pair reached {applied}. {observed}{words["reading"]}'
    else:
        message = f'No pairs were compared against {applied}.'
    return {'minimum_shared_fraction': floor, 'basis': basis, 'basis_words': words['words'],
            'status': status, 'accepted_pairs': accepted,
            'observed_combined_fraction': [min(combined), max(combined)] if combined else None,
            'observed_fraction_of_smaller': [min(smaller), max(smaller)] if smaller else None,
            'message': message,
            'how_to_change_it': 'wmlstudio.snp_tree.at_minimum_shared_fraction re-reads the same run at a '
                                'different floor, and against either denominator, without recomputing a '
                                'single SNP count.'}


def _summary(records, pairs, edges):
    comparable = [pair for pair in pairs if pair['comparable']]
    distances = sorted(pair['distance'] for pair in comparable)
    # Read off the pairs themselves rather than off the drawn edges: an isolate
    # with no comparable pair is a fact about the sequence, not about the layout.
    compared = {str(key) for pair in comparable for key in (pair['source'], pair['target'])}
    return {'isolates': len(records), 'pairs': len(pairs), 'comparable_pairs': len(comparable),
            'excluded_pairs': len(pairs) - len(comparable), 'edges': len(edges),
            'minimum_snps': distances[0] if distances else None,
            'median_snps': statistics.median(distances) if distances else None,
            'maximum_snps': distances[-1] if distances else None,
            'minimum_shared_fraction_observed': min((pair['shared_fraction'] for pair in pairs), default=None),
            'isolates_without_a_comparable_pair': sorted(
                entry['sample_id'] for entry in records if entry['sample_id'] not in compared),
            'unit': DIFFERENCE_WORD,
            'note': 'Minimum, median and maximum are taken over comparable pairs only. Excluded pairs have '
                    'no distance to summarise and are not counted as large or as small.'}


def snp_payload(result, *, organism='', records=None, link_threshold=None, cohort='', cancelled=None):
    """Everything a SNP tab needs from one SKA2 run: matrix, forest, limits and refusals.

    Nothing here recomputes a distance. The run's own numbers are arranged so a
    reader always sees a distance beside the sequence it was measured over, and
    so a picture cannot be drawn without the sentences that say what it is not.
    """
    facts = protocol(result)
    pairs = snp_pairs(result, cancelled=cancelled)
    binding = threshold_binding(result, organism, link_threshold=link_threshold)
    forest = snp_forest(result, records=records, link_threshold=link_threshold, cohort=cohort,
                        cancelled=cancelled)
    alignment = dict(result.get('alignment') or {})
    return {'format_version': 1, 'kind': KIND, 'title': TITLE,
            'engine': f'{result.get("engine")} {result.get("version")}', 'method': result.get('method'),
            'protocol': facts, 'created_at': datetime.now(timezone.utc).isoformat(),
            'run_id': result.get('run_id'), 'output_directory': result.get('output_directory'),
            'cohort': forest['results'], 'pairs': pairs,
            'not_comparable': [dict(pair) for pair in pairs if not pair['comparable']],
            'matrix': snp_matrix(result, cancelled=cancelled), 'graph': forest, 'threshold': binding,
            'split_kmers': result.get('split_kmers'), 'alignment': alignment,
            'alignment_handoff': alignment_handoff(result), 'drawn_tree': DRAWN_TREE,
            'comparability': _comparability(result, pairs),
            'comparability_reread': result.get('comparability_reread'),
            'summary': _summary(forest['results'], pairs, forest['edges']),
            'separation': SEPARATION, 'interpretation': INTERPRETATION,
            'limitations': [*LIMITATIONS, *(result.get('limitations') or []),
                            *(alignment.get('limitations') or [])]}


def build_snp_tree(samples, output_root, *, organism='', records=None, link_threshold=None,
                   cohort='', cancelled=None, progress=None, **options):
    """Run SKA2 over a chosen cohort and return the payload a SNP tab draws.

    The payload is written beside the run's own outputs so the picture, the
    matrix and the refusal that produced them stay together on disk.
    """
    from .export import _atomic_text
    from .ska_runtime import run_ska
    result = run_ska(samples, output_root, cancelled=cancelled, progress=progress, **options)
    payload = snp_payload(result, organism=organism, records=records, link_threshold=link_threshold,
                          cohort=cohort, cancelled=cancelled)
    with _atomic_text(Path(result['output_directory']) / 'snp-tree.json') as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, allow_nan=False)
    return payload
