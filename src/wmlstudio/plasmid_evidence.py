"""Contig-level plasmid evidence from markers already in hand. Not MOB-suite.

Nothing here reconstructs a plasmid. A replicon marker is evidence about the
contig it sits on; it is not a plasmid, not a replicon copy number, not a
mobility prediction and not proof that a neighbouring resistance determinant
travels on a mobile element. Closure and coverage are read from the assembler's
own contig header: they are claims made by the assembler about its own output,
never measurements made here.

This module computes counts and set membership only. It defines no distance and
shares no scale, axis, column or threshold with MLST, cgMLST or SNP distances.

Why not MOB-suite (checked 2026-09-15). MOB-suite 3.1.9 is a Linux-only Python
package (`Operating System :: POSIX :: Linux`) pinning numpy<1.23.5 and
pandas<=1.5.3, neither of which builds or ships a wheel for this application's
Python 3.12 interpreter; its cluster assignment requires Mash, which publishes
Linux and macOS binaries only and whose source build needs Cap'n Proto plus GSL
(GPL-3.0) or Boost; and its database is a 472,851,578 byte archive, larger than
the entire current application ZIP. It is therefore not bundled and not run.
``MOB_SUITE_GAP`` states, in the words a user reads, exactly what a real
MOB-suite run would add that this screen does not provide.
"""

from __future__ import annotations

import re
from datetime import UTC, datetime

from .sequence import check_cancelled

# Naming conventions that carry structured, machine-readable claims. A header's
# free text is never mined: the word "plasmid" in a description is not evidence.
_SKESA = re.compile(r'^Contig_(\d+)_(\d+(?:\.\d+)?)(_Circ)?$')
_SPADES = re.compile(r'^NODE_\d+_length_(\d+)_cov_(\d+(?:\.\d+)?)$')
_KEY_VALUE = re.compile(r'\b(length|depth|cov|circular|topology)\s*=\s*([A-Za-z0-9.]+)', re.IGNORECASE)
_BRACKET = re.compile(r'\[\s*(topology|completeness)\s*=\s*([A-Za-z]+)\s*\]', re.IGNORECASE)

# Contigs at or above this length carry the backbone depth estimate. Below it a
# draft assembly is mostly short repeat fragments whose declared coverage is
# noise. Measured on three real SKESA assemblies (DRR428002, SRR26465518,
# SRR26465495): a 1 kb floor keeps 96-100% of assembled bases in the estimate,
# where a 10 kb floor kept only 43% of the most fragmented one.
BACKBONE_MIN_BP = 1000
# A reporting threshold for flagging a departure from backbone coverage, not a
# plasmid criterion and not a copy-number estimate. On the same three real
# assemblies the assembler-declared circular contigs sat at 4.1x, 7.8x, 4.1x and
# 20.0x the backbone median, and no non-circular contig above 1 kb reached 1.5x.
DEPTH_DEPARTURE_RATIO = 1.5

SUPPORT_STATES = ('no_replicon_marker', 'not_assessable', 'replicon_marker_only',
                  'replicon_marker_plus_depth_departure', 'replicon_marker_plus_closure_claim',
                  'replicon_marker_plus_closure_and_depth')

MOB_SUITE_GAP = [
    'This is a WMLSTudio contig screen over replicon markers already reported by the AMR/plasmid assay. '
    'It is not MOB-suite, not MOB-typer and not MOB-recon, and it is not equivalent to them.',
    'No relaxase (MOB) family is typed and no mate-pair formation (MPF) type is assigned, so no '
    'conjugative, mobilizable or non-mobilizable mobility prediction is made or implied.',
    'No origin of transfer (oriT) is searched for.',
    'No plasmid reconstruction or contig binning is performed: contigs are not grouped into predicted '
    'plasmids, so the number of plasmids in an isolate is not estimated here and cannot be read off this.',
    'No primary or secondary cluster code is assigned against a curated closed-plasmid database, so two '
    'isolates cannot be said to carry the same plasmid on this evidence.',
    'A full MOB-suite run supplies all of the above. Where those questions matter, run MOB-suite on a '
    'Linux host, or resolve the structure with long-read or hybrid assembly.',
]

LIMITATIONS = [
    'A replicon marker is evidence about the contig it was found on. It is not a plasmid, not a count of '
    'plasmids and not proof that the contig is extrachromosomal; replicon-family sequence also occurs on '
    'chromosomes and on integrated elements.',
    'Co-location of a replicon and a resistance or virulence determinant on one assembled contig is a '
    'hypothesis. Short-read assemblies split, collapse and chimerically join repeated plasmid sequence, so '
    'confirmation needs long-read or hybrid assembly, or a closed genome.',
    'Closure and coverage are copied from the assembler\'s own contig header. A circularity claim is the '
    'assembler\'s, not a verified closed molecule, and declared coverage is an assembly statistic, not a '
    'measured read depth. An assembler that states neither leaves both unknown, which is distinct from absent.',
    'The same replicon name in two isolates does not identify the same plasmid, and recurrent '
    'replicon/determinant co-occurrence across a cohort is not transmission, not a shared plasmid and not '
    'a clonal relationship. It is a description of what was seen in the isolates that were assayed.',
]


def parse_contig_header(header):
    """Read structured assembler claims from one FASTA header. No free-text mining.

    Returns declared length, coverage, closure and the convention they came from.
    ``closure`` is 'declared_circular', 'not_declared' when the convention marks
    circular contigs and this one is unmarked, or 'unknown' when the convention
    says nothing about topology at all. Those three are deliberately distinct.
    """
    header = str(header or '').strip()
    identifier = header.split()[0] if header else ''
    record = {'contig': identifier, 'declared_length': None, 'declared_coverage': None,
              'closure': 'unknown', 'assembler_convention': 'unrecognized'}
    skesa = _SKESA.match(identifier)
    spades = _SPADES.match(identifier)
    if skesa:
        record.update(assembler_convention='skesa', declared_coverage=float(skesa.group(2)),
                      closure='declared_circular' if skesa.group(3) else 'not_declared')
    elif spades:
        record.update(assembler_convention='spades', declared_length=int(spades.group(1)),
                      declared_coverage=float(spades.group(2)))
    fields = {key.lower(): value for key, value in _KEY_VALUE.findall(header)}
    fields.update({key.lower(): value for key, value in _BRACKET.findall(header)})
    if 'circular' in fields or 'topology' in fields:
        claim = (fields.get('circular') or fields.get('topology') or '').casefold()
        record['closure'] = ('declared_circular' if claim in {'true', 'yes', '1', 'circular'}
                             else 'not_declared')
        record['assembler_convention'] = 'key_value_header'
    for key in ('depth', 'cov'):
        if key in fields and record['declared_coverage'] is None:
            try:
                record['declared_coverage'] = float(fields[key].rstrip('xX'))
            except ValueError:
                record['declared_coverage'] = None
            record['assembler_convention'] = 'key_value_header'
    if 'length' in fields and record['declared_length'] is None:
        try:
            record['declared_length'] = int(fields['length'])
        except ValueError:
            record['declared_length'] = None
    return record


def contig_records(contig_lengths, contig_headers=None):
    """Join measured contig lengths with whatever the assembler declared about them.

    ``contig_lengths`` is authoritative: it was counted from the sequence bytes.
    A declared length that disagrees is reported, never substituted.
    """
    headers = dict(contig_headers or {})
    records = []
    for contig in sorted(contig_lengths):
        parsed = parse_contig_header(headers.get(contig, contig))
        declared = parsed.pop('declared_length')
        records.append(dict(parsed, contig=contig, length_bp=contig_lengths[contig],
                            declared_length_disagrees=declared is not None and declared != contig_lengths[contig]))
    return records


def backbone_depth(records):
    """Length-weighted median declared coverage of contigs at or above the floor.

    The denominator travels: how many contigs and how many assembled bases the
    estimate rests on, and what fraction of the assembly that is. A low fraction
    means the estimate is weak evidence about the chromosome, and the caller is
    told so rather than left to assume.
    """
    total_bp = sum(record['length_bp'] for record in records)
    usable = [record for record in records
              if record['length_bp'] >= BACKBONE_MIN_BP and isinstance(record['declared_coverage'], float)]
    basis = {'status': 'unavailable', 'median_declared_coverage': None, 'contigs': len(usable),
             'bases': sum(record['length_bp'] for record in usable), 'assembly_bases': total_bp,
             'assembly_fraction': None, 'min_contig_bp': BACKBONE_MIN_BP,
             'reason': 'No contig at or above the length floor carried an assembler-declared coverage.',
             'basis': 'Length-weighted median of assembler-declared contig coverage. This is an assembly '
                      'statistic reported by the assembler, not a read depth measured here.'}
    if not usable or not total_bp:
        return basis
    basis['assembly_fraction'] = basis['bases'] / total_bp
    running = 0
    for record in sorted(usable, key=lambda item: (item['declared_coverage'], item['contig'])):
        running += record['length_bp']
        if running * 2 >= basis['bases']:
            basis['median_declared_coverage'] = record['declared_coverage']
            break
    if not basis['median_declared_coverage']:
        basis['reason'] = 'The backbone median declared coverage is zero; no ratio can be formed from it.'
        return basis
    basis['status'] = 'available'
    basis['reason'] = (f'Median over {basis["contigs"]} contigs of at least {BACKBONE_MIN_BP} bp, '
                       f'covering {basis["assembly_fraction"]:.0%} of assembled bases.')
    return basis


def marker_location(hit, contig_lengths):
    """The one accepted placement rule: a real contig and coordinates inside it."""
    contig = hit.get('contig') or hit.get('sequence')
    start, end = hit.get('start'), hit.get('end')
    if (contig not in contig_lengths or isinstance(start, bool) or isinstance(end, bool)
            or not isinstance(start, int) or not isinstance(end, int)):
        return None
    start, end = sorted((start, end))
    return (contig, start, end) if 1 <= start <= end <= contig_lengths[contig] else None


def _signals(record, basis):
    """Independent signals that a contig behaves unlike the assembly backbone."""
    present, ratio = [], None
    if record['closure'] == 'declared_circular':
        present.append('closure_claim')
    median = basis['median_declared_coverage']
    if basis['status'] == 'available' and isinstance(record['declared_coverage'], float) and median:
        ratio = record['declared_coverage'] / median
        if ratio >= DEPTH_DEPARTURE_RATIO or ratio <= 1 / DEPTH_DEPARTURE_RATIO:
            present.append('depth_departure')
    return present, ratio


def _support_state(record, signals, basis):
    if record['closure'] == 'unknown' and (basis['status'] != 'available'
                                           or not isinstance(record['declared_coverage'], float)):
        return 'not_assessable'
    if 'closure_claim' in signals and 'depth_departure' in signals:
        return 'replicon_marker_plus_closure_and_depth'
    if 'closure_claim' in signals:
        return 'replicon_marker_plus_closure_claim'
    if 'depth_departure' in signals:
        return 'replicon_marker_plus_depth_departure'
    return 'replicon_marker_only'


def _marker_row(hit, position, kind):
    return {'gene': hit.get('gene'), 'accession': hit.get('accession'), 'marker_type': kind,
            'class': hit.get('class'), 'subclass': hit.get('subclass'),
            'start': position[1], 'end': position[2]}


def contig_plasmid_evidence(contig_lengths, replicon_hits, determinant_hits, *, contig_headers=None):
    """Per-contig evidence and per-determinant placement. Nothing is reconstructed.

    ``replicon_hits`` are the accepted replicon markers and ``determinant_hits``
    the accepted AMR/virulence markers, each already filtered by the caller's own
    resolution and identity gates. Every determinant is placed, including the
    ones that land on a contig with no replicon and the ones that cannot be
    placed at all: an absent row would read as absence of the determinant.
    """
    records = {record['contig']: record for record in contig_records(contig_lengths, contig_headers)}
    basis = backbone_depth(list(records.values()))
    # Every replicon is placed before any determinant, so a contig's replicon
    # list is already complete when that contig's determinants are interpreted.
    replicons_by_contig, determinants_by_contig, placement = {}, {}, []
    for hit in replicon_hits:
        position = marker_location(hit, contig_lengths)
        if position is not None:
            replicons_by_contig.setdefault(position[0], []).append(_marker_row(hit, position, 'PLASMID'))
    for hit, kind in determinant_hits:
        position = marker_location(hit, contig_lengths)
        if position is None:
            placement.append({'gene': hit.get('gene'), 'marker_type': kind, 'contig': None,
                              'placement': 'unplaced', 'replicons_on_contig': [],
                              'interpretation': 'No usable contig coordinates were reported for this '
                                                'determinant, so nothing is claimed about where it sits.'})
            continue
        row = _marker_row(hit, position, kind)
        determinants_by_contig.setdefault(position[0], []).append(row)
        co_located = sorted({entry['gene'] for entry in replicons_by_contig.get(position[0], []) if entry['gene']})
        reading = ('Shares one assembled contig with a replicon marker: a co-location hypothesis, not a '
                   'plasmid-borne determinant.') if co_located else (
                  'No replicon marker was found on this contig. That does not place the determinant on the '
                  'chromosome: a plasmid contig can assemble without its replicon, and the replicon may '
                  'have assembled separately.')
        placement.append({'gene': row['gene'], 'marker_type': kind, 'contig': position[0],
                          'placement': 'co_located_with_replicon' if co_located else 'no_replicon_on_this_contig',
                          'replicons_on_contig': co_located, 'interpretation': reading})
    contigs = []
    for contig in sorted(set(replicons_by_contig) | set(determinants_by_contig)):
        record = records[contig]
        signals, ratio = _signals(record, basis)
        replicons = replicons_by_contig.get(contig, [])
        contigs.append({
            'contig': contig, 'length_bp': record['length_bp'],
            'declared_coverage': record['declared_coverage'], 'closure': record['closure'],
            'assembler_convention': record['assembler_convention'],
            'declared_length_disagrees': record['declared_length_disagrees'],
            'depth_ratio': ratio, 'signals': signals,
            'support_state': _support_state(record, signals, basis) if replicons else 'no_replicon_marker',
            'replicons': replicons, 'determinants': determinants_by_contig.get(contig, []),
            'interpretation': 'Evidence about this assembled contig only. No plasmid is reconstructed, '
                              'no mobility is predicted and no plasmid count is implied.'})
    return {'contigs': contigs, 'determinant_placement': placement, 'depth_basis': basis,
            'depth_departure_ratio': DEPTH_DEPARTURE_RATIO,
            'unplaced_determinants': sum(row['placement'] == 'unplaced' for row in placement),
            'mob_suite_gap': list(MOB_SUITE_GAP), 'limitations': list(LIMITATIONS)}


def _assayed(entry):
    """An isolate counts towards the denominator only if its plasmid assay ran."""
    block = entry.get('plasmid_hypotheses') or {}
    status = block.get('status')
    if status == 'completed':
        return block, ''
    return None, block.get('reason') or f'Plasmid-marker assay status is {status or "missing"}.'


def cohort_replicon_cooccurrence(isolates, *, cancelled=None):
    """Which replicons, and which replicon/determinant pairs, recur across a cohort.

    Every count carries its denominator, and isolates whose plasmid assay did not
    run are named with the reason rather than counted as negative. Isolates where
    a replicon and a determinant are both present but on different contigs are
    reported separately from the co-located ones, because collapsing the two
    would turn an assembly artefact into a shared finding.

    This is co-occurrence across isolates. It is not transmission, not a shared
    plasmid and not a distance; it produces no threshold anyone could cluster on.
    """
    assayed, excluded, unplaced = [], [], []
    replicons, co_located, apart = {}, {}, {}
    for index, entry in enumerate(isolates):
        if index % 32 == 0:
            check_cancelled(cancelled)
        sample_id = str(entry.get('sample_id') or entry.get('id') or '')
        if not sample_id:
            raise ValueError('Every cohort entry needs a sample_id.')
        block, reason = _assayed(entry)
        if block is None:
            excluded.append({'sample_id': sample_id, 'sample_name': str(entry.get('sample_name') or sample_id),
                             'reason': reason})
            continue
        assayed.append(sample_id)
        present_replicons = {str(hit.get('gene')) for hit in block.get('replicons') or [] if hit.get('gene')}
        for replicon in present_replicons:
            replicons.setdefault(replicon, set()).add(sample_id)
        pairs = {(str(link.get('replicon')), str(link.get('marker')), str(link.get('marker_type') or ''))
                 for link in block.get('contig_associations') or []}
        for key in pairs:
            co_located.setdefault(key, set()).add(sample_id)
        for row in block.get('determinant_placement') or []:
            gene, kind = str(row.get('gene')), str(row.get('marker_type') or '')
            # A determinant with no usable coordinates is neither co-located nor
            # apart. Counting it as apart would read as evidence of separation.
            if row.get('placement') == 'unplaced':
                unplaced.append({'sample_id': sample_id, 'gene': gene, 'marker_type': kind})
                continue
            for replicon in present_replicons:
                if (replicon, gene, kind) not in pairs:
                    apart.setdefault((replicon, gene, kind), set()).add(sample_id)
    denominator = len(assayed)
    pair_keys = sorted(set(co_located) | set(apart))
    return {
        'isolates_assayed': sorted(assayed), 'denominator': denominator,
        'isolates_excluded': sorted(excluded, key=lambda row: row['sample_id']),
        'unplaced_determinants': sorted(unplaced, key=lambda row: (row['sample_id'], row['gene'])),
        'replicons': [{'replicon': name, 'isolates': sorted(members), 'isolate_count': len(members),
                       'denominator': denominator}
                      for name, members in sorted(replicons.items())],
        'co_occurrence': [{'replicon': replicon, 'marker': marker, 'marker_type': kind,
                           'co_located_isolates': sorted(co_located.get(key, ())),
                           'co_located_count': len(co_located.get(key, ())),
                           'present_but_not_co_located': sorted(apart.get(key, ())),
                           'not_co_located_count': len(apart.get(key, ())),
                           'denominator': denominator, 'status': 'hypothesis',
                           'interpretation': 'The same replicon marker and the same determinant shared one '
                                             'assembled contig in these isolates. Recurrence is worth '
                                             'investigating; it is not transmission, not one plasmid and '
                                             'not a relatedness measure.'}
                          for key in pair_keys for replicon, marker, kind in [key]],
        'limitations': list(LIMITATIONS), 'mob_suite_gap': list(MOB_SUITE_GAP),
    }


# --- was the question even asked ------------------------------------------
# A reference store holding no plasmid set reports no replicon for any isolate
# on earth. Rendering that as an empty replicon column is the false negative
# this half of the module exists to prevent, so "not installed", "installed but
# not read", "read and found nothing" and "no assay to read" are four separate
# answers that never collapse into one blank cell.

#: The screen states, in the order of how much was actually done. Every one of
#: them is a statement about the evidence, never about the isolate.
SCREEN_STATES = ('replicons_detected', 'none_detected', 'not_searched',
                 'plasmid_database_not_installed', 'not_assayed')

#: The same five states in the words a reader sees in a cell.
SCREEN_WORDS = {
    'replicons_detected': 'replicon marker(s) found',
    'none_detected': 'searched, none found',
    'not_searched': 'not searched — no plasmid reference set was read',
    'plasmid_database_not_installed': 'not searched — no plasmid database is installed',
    'not_assayed': 'no plasmid assay to read for this isolate',
}

#: What a matrix cell says. A replicon that was never looked for must not share a
#: word, a blank or a colour with a replicon that was looked for and not found.
MATRIX_WORDS = {'replicons_detected': ('present', 'not detected'),
                'none_detected': ('present', 'not detected'),
                'not_searched': ('present', 'not searched'),
                'plasmid_database_not_installed': ('present', 'not searched'),
                'not_assayed': ('present', 'not assayed')}

NO_DATABASE_REASON = (
    'No reference set that reports plasmid replicon types is installed, so no replicon could be '
    'reported for any isolate in this project. Install plasmidfinder from the AMR / plasmid '
    'reference databases page. An empty replicon table here means the question was never asked; '
    'it is not evidence that these isolates carry no plasmid.')

UNKNOWN_DATABASE_REASON = (
    'The reference store was not inspected, so whether a plasmid reference set is installed is '
    'unknown here. An empty replicon table is therefore not a negative result.')

#: The limit that belongs beside each table, printed with the numbers rather than
#: in a help page. Each is severe enough that a reader who misses it over-reads
#: the table it sits under.
TABLE_NOTES = {
    'replicons_by_isolate':
        'A replicon marker is a sequence match to a replicon family. It is not a plasmid, not a '
        'count of plasmids and not proof that anything in this isolate is extrachromosomal: '
        'replicon-family sequence also occurs on chromosomes and on integrated elements.',
    'replicon_matrix':
        'Each cell says what was seen in one isolate. The same replicon name in two isolates does '
        'not identify the same plasmid, and this grid is not a distance: it carries no scale and '
        'no threshold anybody could cluster on.',
    'determinant_colocation':
        'Sharing one assembled contig is a hypothesis, not a plasmid-borne determinant. A '
        'determinant on a contig with no replicon is not thereby chromosomal, and a determinant '
        'with no usable coordinates is neither co-located nor apart.',
    'cohort_cooccurrence':
        'Counts of what was seen in the isolates that were assayed, each row carrying its own '
        'denominator. Recurrence across isolates is not transmission, not one shared plasmid and '
        'not a relatedness measure.',
    'contigs':
        "Closure and coverage are the assembler's own claims about its own output, quoted back "
        'here and measured nowhere. A short-read assembly cannot resolve plasmid structure, so no '
        'row in this table is a plasmid and no plasmid is counted from it.',
    'isolates_not_assayed':
        'These isolates are named rather than counted as carrying nothing. Every count on this '
        'page is against the isolates that were assayed, and they are not among them.',
}


def plasmid_reference_sets():
    """Which reference sets can report a replicon at all, from the engine's own map.

    Read from :mod:`wmlstudio.hydra_runtime` rather than restated, so a set added
    there cannot silently stop counting as a plasmid database here.
    """
    from .hydra_runtime import ELEMENT_SOURCES
    return tuple(ELEMENT_SOURCES['PLASMID'])


def replicon_database_state(installed=None):
    """Whether any installed reference set can report a replicon marker at all.

    ``installed`` names what the reference store holds. ``None`` means nobody
    asked it, which is a third answer and not a negative one: only the store can
    say whether a set is absent, so without it nothing is called absent.
    """
    sources = list(plasmid_reference_sets())
    if installed is None:
        return {'status': 'unknown', 'known': False, 'installed': [], 'missing': sources,
                'reason': UNKNOWN_DATABASE_REASON}
    held = {str(name) for name in installed}
    present = [name for name in sources if name in held]
    missing = [name for name in sources if name not in held]
    if present:
        return {'status': 'installed', 'known': True, 'installed': present, 'missing': missing,
                'reason': f"{', '.join(present)} is installed, so a replicon marker can be "
                          'reported for an isolate that was screened against it.'}
    return {'status': 'not_installed', 'known': True, 'installed': [], 'missing': missing,
            'reason': NO_DATABASE_REASON}


def _entry_databases(entry, block):
    """Which reference sets this isolate's own screen recorded having read."""
    names = entry.get('databases')
    if names is None:
        names = (block.get('source') or {}).get('databases')
    return [str(name) for name in (names or [])]


def isolate_replicon_view(entry, *, database=None):
    """One isolate's replicon screen: what ran, what it read, and what it found.

    The state is decided from the isolate's own recorded run first and from the
    store second, because a run that did read a plasmid set really did look, even
    if that set has since been removed from the store.
    """
    database = replicon_database_state() if database is None else database
    sample_id = str(entry.get('sample_id') or entry.get('id') or '')
    if not sample_id:
        raise ValueError('Every cohort entry needs a sample_id.')
    block = entry.get('plasmid_hypotheses') or {}
    searched = [name for name in plasmid_reference_sets() if name in _entry_databases(entry, block)]
    replicons = sorted({str(hit.get('gene')) for hit in block.get('replicons') or [] if hit.get('gene')})
    if block.get('status') != 'completed':
        state = 'not_assayed'
        reason = block.get('reason') or f"Plasmid-marker assay status is {block.get('status') or 'missing'}."
    elif replicons:
        state, reason = 'replicons_detected', block.get('reason') or ''
    elif searched:
        state = 'none_detected'
        reason = (f"{', '.join(searched)} was searched and no replicon marker met this run's "
                  'thresholds. A screen over public reference data is not proof of absence.')
    elif database['status'] == 'not_installed':
        state, reason = 'plasmid_database_not_installed', NO_DATABASE_REASON
    else:
        state = 'not_searched'
        reason = ('No reference set that reports plasmid replicon types was read for this isolate, '
                  'so nothing was looked for. That is not a negative result.')
    return {'sample_id': sample_id, 'sample_name': str(entry.get('sample_name') or sample_id),
            'state': state, 'state_words': SCREEN_WORDS[state], 'reason': reason,
            'assayed': block.get('status') == 'completed', 'searched': searched,
            'replicons': replicons, 'contigs': list(block.get('contig_evidence') or []),
            'determinant_placement': list(block.get('determinant_placement') or []),
            'depth_basis': dict(block.get('depth_basis') or {})}


def _table(key, title, columns, rows, empty_reason):
    return {'key': key, 'title': title, 'columns': list(columns), 'rows': [list(row) for row in rows],
            'note': TABLE_NOTES[key], 'empty_reason': '' if rows else empty_reason,
            'row_count': len(rows)}


def _genes(row, field):
    return '; '.join(entry['gene'] for entry in row.get(field) or [] if entry.get('gene')) or '—'


def _replicons_by_isolate_table(views, database):
    rows = [[view['sample_name'], view['state_words'],
             '; '.join(view['replicons']) or '—', len(view['replicons']),
             ', '.join(view['searched']) or 'none', view['reason']]
            for view in views]
    empty = (NO_DATABASE_REASON if database['status'] == 'not_installed'
             else 'No isolate is in this cohort yet, so nothing was read either way.')
    return _table('replicons_by_isolate', 'Replicon markers, isolate by isolate',
                  ['Isolate', 'Replicon screen', 'Replicon markers', 'Markers found',
                   'Reference sets read', 'What this row means'], rows, empty)


def _matrix_table(views, database):
    names = sorted({name for view in views for name in view['replicons']})
    rows = []
    for view in views:
        present, absent = MATRIX_WORDS[view['state']]
        rows.append([view['sample_name'], view['state_words'],
                     *[present if name in view['replicons'] else absent for name in names]])
    empty = (NO_DATABASE_REASON if database['status'] == 'not_installed' else
             'No replicon marker was reported for any isolate here, so the grid has no columns to '
             'draw. Read the screen state beside each isolate before reading that as absence.')
    return _table('replicon_matrix', 'Replicon by isolate',
                  ['Isolate', 'Replicon screen', *names], rows if names else [], empty)


def _colocation_table(views):
    rows = []
    for view in views:
        for row in view['determinant_placement']:
            rows.append([view['sample_name'], row.get('gene') or '—', row.get('marker_type') or '—',
                         row.get('contig') or 'no usable contig coordinates',
                         '; '.join(row.get('replicons_on_contig') or []) or '—',
                         str(row.get('placement') or '—').replace('_', ' '),
                         row.get('interpretation') or ''])
    return _table('determinant_colocation',
                  'Determinants and the contigs they sit on',
                  ['Isolate', 'Determinant', 'Marker type', 'Contig', 'Replicons on that contig',
                   'Placement', 'What this means'], rows,
                  'No resistance or virulence determinant was placed on a contig for any assayed '
                  'isolate, so there is nothing to co-locate. That is not a statement that these '
                  'isolates carry no determinants.')


def _cooccurrence_table(cohort, names):
    rows = [[row['replicon'], row['marker'], row['marker_type'] or '—',
             row['co_located_count'], row['not_co_located_count'], row['denominator'],
             '; '.join(names.get(value, value) for value in row['co_located_isolates']) or '—',
             row['interpretation']]
            for row in cohort['co_occurrence']]
    return _table('cohort_cooccurrence', 'Replicon and determinant on one assembled contig',
                  ['Replicon marker', 'Determinant', 'Marker type', 'Co-located on one contig',
                   'Both present, not co-located', 'Of isolates assayed', 'Isolates',
                   'What this is'], rows,
                  'No replicon and determinant shared an assembled contig in the isolates that '
                  'were assayed. Contigs that carry a determinant and no replicon are in the '
                  'table above, and they place nothing on a chromosome.')


def _contig_table(views):
    rows = []
    for view in views:
        for row in view['contigs']:
            ratio = row.get('depth_ratio')
            coverage = row.get('declared_coverage')
            rows.append([view['sample_name'], row.get('contig'), row.get('length_bp'),
                         str(row.get('closure') or 'unknown').replace('_', ' '),
                         f'{coverage:g}' if isinstance(coverage, float) else 'not declared',
                         f'{ratio:.1f}×' if isinstance(ratio, float) else 'unknown',
                         _genes(row, 'replicons'), _genes(row, 'determinants'),
                         str(row.get('support_state') or '').replace('_', ' '),
                         row.get('assembler_convention') or 'unrecognized'])
    return _table('contigs', 'Every contig carrying a replicon or a determinant',
                  ['Isolate', 'Contig', 'Length (bp)', "Closure (assembler's claim)",
                   'Declared coverage', 'Coverage vs backbone', 'Replicon markers',
                   'Determinants on the same contig', 'Evidence present', 'Header convention'],
                  rows,
                  'No contig in this cohort carries a replicon marker or a placed determinant, so '
                  'there is no contig to describe. Whether anything was searched for is in the '
                  'screen state beside each isolate.')


def _not_assayed_table(views):
    rows = [[view['sample_name'], view['state_words'], view['reason']]
            for view in views if not view['assayed']]
    return _table('isolates_not_assayed', 'Isolates with no plasmid result to read',
                  ['Isolate', 'Replicon screen', 'Why'], rows,
                  'Every isolate in this cohort had a plasmid-marker assay to read.')


def plasmid_payload(entries, *, installed=None, cohort='', cancelled=None):
    """Every plasmid table this application can honestly draw, as plain data.

    One payload serves the page and the report, so a table a reader saw on screen
    and the same table in a document are the same numbers with the same words
    beside them. Nothing here is a distance: no row shares a scale, an axis, a
    column or a threshold with an MLST, cgMLST or SNP measurement.

    ``installed`` names what the reference store holds, so "no plasmid database
    is installed" never renders as "no replicon was found".
    """
    database = replicon_database_state(installed)
    views, assayable = [], []
    for index, entry in enumerate(entries):
        if index % 32 == 0:
            check_cancelled(cancelled)
        views.append(isolate_replicon_view(entry, database=database))
        assayable.append({'sample_id': views[-1]['sample_id'],
                          'sample_name': views[-1]['sample_name'],
                          'plasmid_hypotheses': entry.get('plasmid_hypotheses') or {}})
    views.sort(key=lambda view: (view['sample_name'].casefold(), view['sample_id']))
    cohort_counts = cohort_replicon_cooccurrence(assayable, cancelled=cancelled)
    names = {view['sample_id']: view['sample_name'] for view in views}
    tables = {'replicons_by_isolate': _replicons_by_isolate_table(views, database),
              'replicon_matrix': _matrix_table(views, database),
              'determinant_colocation': _colocation_table(views),
              'cohort_cooccurrence': _cooccurrence_table(cohort_counts, names),
              'contigs': _contig_table(views),
              'isolates_not_assayed': _not_assayed_table(views)}
    searched = sorted({name for view in views for name in view['searched']})
    return {
        'kind': 'plasmid_replicon_evidence',
        'title': 'Plasmid replicon evidence',
        'cohort': str(cohort or ''),
        'created_utc': datetime.now(UTC).isoformat(),
        'database': database,
        'databases_read': searched,
        'isolate_count': len(views),
        'denominator': cohort_counts['denominator'],
        'isolates_assayed': cohort_counts['isolates_assayed'],
        'isolates': [dict(view) for view in views],
        'unplaced_determinants': list(cohort_counts['unplaced_determinants']),
        'tables': tables,
        'table_order': list(TABLE_NOTES),
        'headline': headline(views, database, cohort_counts['denominator']),
        'limitations': list(LIMITATIONS),
        'mob_suite_gap': list(MOB_SUITE_GAP),
    }


def headline(views, database, denominator):
    """One sentence saying what was asked, of how many, before any table is read."""
    if not views:
        return ('No isolate is in this cohort yet. Choose isolates to see what their plasmid '
                'screens reported.')
    if database['status'] == 'not_installed' and not any(view['searched'] for view in views):
        return NO_DATABASE_REASON
    found = sum(1 for view in views if view['replicons'])
    return (f'{denominator} of {len(views)} isolate(s) had a plasmid-marker assay to read, and '
            f'{found} of those carried at least one replicon marker. The isolates without an '
            'assay are listed rather than counted as carrying nothing.')
