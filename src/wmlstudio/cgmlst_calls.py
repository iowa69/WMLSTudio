"""The cgMLST table of calls, and an honest read recheck of targets that did not call.

Two Qt-free things live here, so the dedicated cgMLST page can be tested without a screen.

THE TABLE is one row per target and one column per isolate. A cgMLST scheme is
thousands of targets and a cohort is tens of isolates, so this orientation keeps the
row count equal to the scheme's own size -- the denominator is the length of the table
-- and the column count equal to the cohort; the transpose would need thousands of
columns, and the recheck button acts on targets, not on isolates. The per-isolate
reading is not lost: every isolate carries its own summary beside the grid.

Every cell states what was actually found. An assigned reference allele, a validated
local full-CDS identity, and each separate reason there is no allele (missing, partial,
below similarity thresholds, ambiguous, duplicated, mixed) are distinct states and never
collapse into one another. Every percentage carries the scheme's full target count,
because a numerator without its denominator is the misreading this application exists to
prevent.

THE RECHECK reruns the existing bounded read-support assay (read_support.py) on targets
the assembly did not call, and reports what the reads show: the target is present, the
reads are consistent with it being absent, or there is not enough evidence to say. A
target missing from an assembly is not evidence of absence in the isolate -- a coverage
gap, an assembly break or a contig edge leave exactly the same empty cell -- so a
negative verdict is offered only when a control target that the assembly did call was
itself recovered from these same reads and the complete read files were examined. No
verdict here ever becomes an allele call: reads support a target being present, they do
not assign a validated allele.
"""

from __future__ import annotations

from collections.abc import Mapping
from datetime import UTC, datetime

from . import __version__
from .project import classify_typing
from .read_support import LIMITATIONS as READ_SUPPORT_LIMITATIONS
from .read_support import MAX_PANEL_LOCI, investigate_read_support
from .sample_workflow import _sha256
from .sequence import check_cancelled
from .threshold_guidance import typing_scale
from .typing import Scheme, _natural_key

# Every state a cell can hold, with the sentence the page prints beside it. The
# order is the order a legend and a summary bar read in: what was called first,
# then each distinct reason there is no allele. `has_allele` is the only thing
# that may be counted as a call; nothing else is ever folded into it.
CALL_STATES = {
    'called': {
        'label': 'allele called', 'has_allele': True,
        'meaning': 'An exact match to one reference allele of this target in the installed scheme.'},
    'called_novel': {
        'label': 'novel sequence (local)', 'has_allele': True,
        'meaning': 'A complete predicted coding sequence that matches no reference allele. It is identified by '
                   'the SHA-256 of its own sequence and is not a registered allele number in any nomenclature.'},
    'missing': {
        'label': 'no call', 'has_allele': False,
        'meaning': 'No reference allele matched and no homolog qualified. This is not evidence that the target is '
                   'absent from the isolate: a coverage gap, an assembly break or a contig edge look the same.'},
    'partial': {
        'label': 'partial (incomplete CDS)', 'has_allele': False,
        'meaning': 'A homolog was found but it is not a complete coding sequence, typically a contig edge or a '
                   'truncation. No allele is assigned.'},
    'low_similarity': {
        'label': 'below similarity thresholds', 'has_allele': False,
        'meaning': 'A homolog was found but fell below the identity or coverage thresholds in force for this run. '
                   'No allele is assigned.'},
    'ambiguous': {
        'label': 'ambiguous', 'has_allele': False,
        'meaning': 'More than one candidate, or sequence shared with another target, prevents a unique assignment.'},
    'duplicated': {
        'label': 'duplicated', 'has_allele': False,
        'meaning': 'The target was found at more than one separate place in this assembly, so no single copy can '
                   'be assigned.'},
    'mixed': {
        'label': 'mixed', 'has_allele': False,
        'meaning': 'Different allele candidates occur at distinct positions. Review the assembly for duplication or '
                   'a mixed culture; this is a property of the assembly, not a contamination diagnosis.'},
    'unknown': {
        'label': 'unknown', 'has_allele': False,
        'meaning': 'This profile does not carry a call state this version recognises, so nothing is claimed for '
                   'this cell. It is not a missing target and not a call.'},
}
STATE_ORDER = tuple(CALL_STATES)

# The caller statuses typing.call_assembly and cgtyping.call_cgassembly actually
# emit. Anything else becomes 'unknown' rather than being guessed into a state.
_STATUS_STATES = {
    'exact': 'called', 'novel_validated': 'called_novel', 'missing': 'missing',
    'partial': 'partial', 'low_similarity': 'low_similarity', 'ambiguous': 'ambiguous',
    'mixed': 'mixed',
}
# Placeholders an imported profile table uses for "no allele here". They are read
# as missing, never as an allele named '0' or '-'.
_NO_ALLELE_TOKENS = frozenset({'', '-', '?', '0', 'n', 'na', 'n/a', 'none', 'null', 'unknown', 'missing'})

TABLE_LIMITATIONS = [
    "Percentages are always over the scheme's full target count. A target with no call is unknown: it is never "
    'counted as a match and never as a difference.',
    'An assigned allele, a validated local novel sequence, and each reason for no call are distinct states here '
    'and are never merged into one another.',
    'A local novel identity is the SHA-256 of the observed sequence, not a registered allele number, and two '
    'isolates share one only when their sequences are identical.',
    'A duplicated or mixed target is a property of this assembly. It is not a species assignment, a contamination '
    'diagnosis, or any statement about antimicrobial susceptibility.',
    'A cgMLST distance over these targets and a classical seven-locus MLST distance are different quantities: they '
    'share no scale, no column and no threshold.',
]

RECHECK_LIMITATIONS = [
    'A read check reports what the reads show. It never assigns an allele: a target recovered from reads stays '
    'uncalled until an allele call is made from an assembly.',
    'A target missing from an assembly is not evidence of absence in the isolate. A coverage gap, an assembly '
    'break or a contig edge leave exactly the same empty cell.',
    'Absence is offered only when a control target that the assembly did call was itself recovered from these '
    'reads and the complete read files were examined. Even then it is consistent with absence, not proof.',
    'Reads are tested against the reference alleles of this scheme only. A divergent or unrepresented allele can '
    'look exactly like an absent target.',
    'Depth counts independently aligned reads, including overlapping mates and duplicates; it is not the number of '
    'DNA molecules observed.',
    'This check changes no allele call, no cgMLST profile, no distance and no tree.',
]

# What each verdict means, in the words the page prints. `evidence` is the
# three-way answer the button promises: present, absent, or not enough to say.
RECHECK_VERDICTS = {
    'present_in_reads': {
        'evidence': 'present', 'label': 'Present in the reads',
        'meaning': 'The reads carry this target even though the assembly produced no call. A coverage gap, an '
                   'assembly break or a contig edge can all leave a target uncalled in an assembly of an isolate '
                   'that has it. This is read evidence of presence, not an allele.'},
    'present_with_mixed_bases': {
        'evidence': 'present', 'label': 'Present, with disagreeing bases',
        'meaning': 'Reads cover this target but disagree at one or more positions. That is consistent with more '
                   'than one template in the sample, with a paralogous copy, or with sequencing error. The target '
                   'is present in the reads; which allele it carries is not established.'},
    'partial_read_evidence': {
        'evidence': 'insufficient', 'label': 'Partial read evidence only',
        'meaning': 'Reads align to part of this target only, or align with disruptions. That is neither the '
                   'presence of a complete target nor evidence of absence.'},
    'consistent_with_absence': {
        'evidence': 'absent', 'label': 'Consistent with absence',
        'meaning': 'No read in the complete pair aligned to any reference allele of this target, while a control '
                   'target the assembly did call was recovered from these same reads. That is consistent with the '
                   'target being absent from this isolate. It is not proof: a divergent allele, an allele absent '
                   'from this scheme, or a target outside the tested reference set would look the same.'},
    'no_evidence_either_way': {
        'evidence': 'insufficient', 'label': 'Not enough evidence to say',
        'meaning': 'No read aligned to this target, and the run does not license reading that as absence.'},
}


def _entry(record):
    """Accept either a stored analysis result or a project sample record holding one."""
    if not isinstance(record, Mapping):
        raise ValueError('Every profile in a calls table must be a stored result or a sample record.')
    result = record
    if 'alleles' not in record and isinstance(record.get('result'), Mapping):
        result = record['result']
    key = (record.get('sample_id') or record.get('id') or result.get('sample_id')
           or result.get('id') or record.get('name') or result.get('sample_name'))
    if key is None or not str(key).strip():
        raise ValueError('Every profile in a calls table needs a sample_id, id or sample_name.')
    name = record.get('name') or result.get('sample_name') or str(key)
    return str(key), str(name), result


def _distinct_copies(call):
    """How many separate places in the assembly this target was seen.

    Overlapping evidence for one place is one copy; the per-locus evidence list is
    capped at 50 hits upstream, so a truncated list makes this a lower bound and
    the cell says so rather than reporting a confident count.
    """
    intervals = set()
    for hit in call.get('hits') or ():
        if not isinstance(hit, Mapping):
            continue
        contig, start, end = hit.get('contig'), hit.get('start'), hit.get('end')
        if contig is None or isinstance(start, bool) or isinstance(end, bool) \
                or not isinstance(start, int) or not isinstance(end, int):
            continue
        intervals.add((str(contig), min(start, end), max(start, end)))
    copies, current = 0, None
    for contig, start, end in sorted(intervals):
        if current is None or current[0] != contig or current[2] < start:
            copies, current = copies + 1, (contig, start, end)
        else:
            current = (contig, current[1], max(current[2], end))
    candidates = call.get('cds_candidates')
    if isinstance(candidates, (list, tuple)):
        copies = max(copies, len({str(value) for value in candidates}))
    return copies


def classify_call(call, allele):
    """Return the cell state for one target, refusing to guess an unrecognised one.

    `call` is the per-locus evidence row a caller produced, or an empty mapping for
    an imported profile that carries alleles without evidence. The allele and the
    evidence must agree: a state that owns an allele without one, or an allele under
    a state that owns none, is reported as unknown rather than silently believed.
    """
    call = call if isinstance(call, Mapping) else {}
    token = '' if allele is None else str(allele).strip()
    copies = _distinct_copies(call)
    truncated = bool(call.get('evidence_truncated'))
    status = str(call.get('status') or '').strip().casefold()
    reason = str(call.get('reason') or '')
    if not call:
        # An imported profile states alleles without per-target evidence. The
        # allele is the call; a placeholder, or nothing at all, is the absence of
        # one, and says only that this profile records no allele here.
        if token.casefold() in _NO_ALLELE_TOKENS:
            state, reason = 'missing', 'This profile records no allele for this target and carries no call evidence.'
        elif token.casefold() == 'novel':
            state, reason = 'unknown', ('The profile says "novel" without a sequence identity, so no novel call '
                                        'can be shown or compared.')
        else:
            state, reason = 'called', 'Imported profile: an allele was supplied without per-target call evidence.'
        return {'state': state, 'allele': token or None, 'status': status or 'imported',
                'copies': copies, 'copies_lower_bound': truncated, 'reason': reason}
    state = _STATUS_STATES.get(status, 'unknown')
    if state == 'unknown' and not reason:
        reason = f'Call state {status or "(absent)"!r} is not one this version recognises.'
    # Ambiguity caused by more than one copy is its own answer to the user's
    # question, so it is named rather than left inside 'ambiguous'.
    if state == 'ambiguous' and copies > 1:
        state = 'duplicated'
    if CALL_STATES[state]['has_allele'] and not token:
        state, reason = 'unknown', (f'The evidence says {status!r} but no allele is stored for this target. ' + reason).strip()
    elif not CALL_STATES[state]['has_allele'] and token and token.casefold() not in _NO_ALLELE_TOKENS:
        state, reason = 'unknown', (f'An allele is stored under call state {status!r}, which assigns none. ' + reason).strip()
    return {'state': state, 'allele': token or None, 'status': status,
            'copies': copies, 'copies_lower_bound': truncated, 'reason': reason}


def _display(cell):
    """The short text a grid cell shows; the honest wording lives in one place."""
    state = cell['state']
    if state == 'called':
        return cell['allele']
    if state == 'called_novel':
        token = (cell['allele'] or '').removeprefix('NOVEL_')
        return 'local novel · ' + token[:12] if token else CALL_STATES[state]['label']
    if state == 'duplicated' and cell['copies'] > 1:
        return f'duplicated ×{cell["copies"]}' + ('+' if cell['copies_lower_bound'] else '')
    return CALL_STATES[state]['label']


def _scheme_identity(entries, scheme):
    """One table, one scheme. Different fingerprints are different quantities."""
    digests = sorted({str(result.get('scheme_digest') or '').strip() for _, _, result in entries})
    names = sorted({str(result.get('scheme') or '').strip() for _, _, result in entries})
    if len(digests) != 1 or not digests[0]:
        raise ValueError('A table of calls covers one scheme. These profiles carry '
                         f'{len(digests)} different scheme fingerprints ({", ".join(value[:12] or "(none)" for value in digests)}). '
                         'Build one table per scheme; profiles called against different schemes are not comparable.')
    if len(names) != 1:
        raise ValueError('These profiles share a scheme fingerprint but record different scheme names '
                         f'({", ".join(names)}). Resolve the naming before building one table from them.')
    locus_sets = {frozenset(str(locus) for locus in (result.get('alleles') or {})) for _, _, result in entries}
    if len(locus_sets) != 1:
        raise ValueError('These profiles claim the same scheme fingerprint but list different targets. '
                         'One of them was edited or imported against a different target set; no table is built from both.')
    targets = next(iter(locus_sets))
    if not targets:
        raise ValueError('None of these profiles carries an allelic profile, so there is no table of calls to show.')
    if scheme is not None:
        if not isinstance(scheme, Scheme):
            raise ValueError('Pass a loaded Scheme, or None to take the target set from the profiles themselves.')
        if scheme.digest != digests[0]:
            raise ValueError('The loaded scheme is not the one these profiles were called against '
                             f'({scheme.digest[:12]} against {digests[0][:12]}).')
        if set(scheme.loci) != targets:
            raise ValueError('The loaded scheme and these profiles list different targets under the same fingerprint.')
        return names[0], digests[0], tuple(scheme.loci)
    return names[0], digests[0], tuple(sorted(targets, key=_natural_key))


def build_calls_table(results, *, scheme=None, cancelled=None, progress=None):
    """One row per target, one column per isolate, with both summaries and the denominator.

    `results` are stored analyses (or project sample records holding one) that were
    all called against the same scheme. Pass the loaded `scheme` to take its own
    target order and to have the fingerprint checked against it.
    """
    check_cancelled(cancelled)
    entries = [_entry(record) for record in results]
    if not entries:
        raise ValueError('Select at least one typed isolate to show a table of calls.')
    keys = [key for key, _, _ in entries]
    duplicated = sorted({key for key in keys if keys.count(key) > 1})
    if duplicated:
        raise ValueError('Each isolate appears once in a table of calls. Repeated identifiers: '
                         + ', '.join(duplicated))
    name, digest, targets = _scheme_identity(entries, scheme)
    profiles = [(key, label, result, result.get('alleles') or {},
                 {str(call.get('locus')): call for call in (result.get('calls') or ())
                  if isinstance(call, Mapping)}) for key, label, result in entries]
    rows, state_totals = [], dict.fromkeys(STATE_ORDER, 0)
    per_sample = {key: dict.fromkeys(STATE_ORDER, 0) for key, _, _ in entries}
    for index, locus in enumerate(targets):
        if index % 256 == 0:
            check_cancelled(cancelled)
            if progress:
                progress(index, len(targets), f'Building the table of calls · {index:,} of {len(targets):,} targets')
        cells, states = {}, dict.fromkeys(STATE_ORDER, 0)
        for key, _, _, alleles, calls in profiles:
            cell = classify_call(calls.get(locus), alleles.get(locus))
            cell['display'] = _display(cell)
            cells[key] = cell
            states[cell['state']] += 1
            state_totals[cell['state']] += 1
            per_sample[key][cell['state']] += 1
        called = sum(states[state] for state in STATE_ORDER if CALL_STATES[state]['has_allele'])
        rows.append({
            'locus': locus, 'cells': cells, 'states': states,
            'called': called, 'without_call': len(profiles) - called,
            'isolates': len(profiles),
            'called_pct': 100 * called / len(profiles),
            'called_in_every_isolate': called == len(profiles),
            'called_in_no_isolate': called == 0,
        })
    check_cancelled(cancelled)
    samples = []
    for key, label, result, _, _ in profiles:
        states = per_sample[key]
        called = sum(states[state] for state in STATE_ORDER if CALL_STATES[state]['has_allele'])
        samples.append({
            'sample_id': key, 'sample_name': label, 'states': states,
            'targets': len(targets), 'called': called,
            'called_reference': states['called'], 'called_novel': states['called_novel'],
            'without_call': len(targets) - called,
            'called_pct': 100 * called / len(targets) if targets else 0.0,
            'status': str(result.get('status') or ''),
            'analysis_kind': classify_typing(result),
            'input_sha256': _sha256(result.get('input_sha256')),
            'cg_profile_id': result.get('cg_profile_id'),
            'summary': f'{called:,} of {len(targets):,} targets called '
                       f'({100 * called / len(targets) if targets else 0:.1f}% of the scheme)',
        })
    if progress:
        progress(len(targets), len(targets), f'Table of calls ready · {len(targets):,} targets × {len(samples):,} isolates')
    return {
        'format_version': 1,
        'scheme': name, 'scheme_digest': digest, 'target_count': len(targets),
        'typing': typing_scale(len(targets)),
        'isolate_count': len(samples),
        'samples': samples, 'loci': rows,
        'states': {state: dict(CALL_STATES[state]) for state in STATE_ORDER},
        'state_totals': state_totals,
        'targets_called_in_every_isolate': sum(row['called_in_every_isolate'] for row in rows),
        'targets_called_in_no_isolate': sum(row['called_in_no_isolate'] for row in rows),
        'denominator': f'{len(targets):,} targets in {name}',
        'denominator_note': f'Every percentage on this page is over the {len(targets):,} targets this scheme '
                            'defines, not over the targets that happened to call.',
        'provenance': {'software': 'WMLSTudio', 'version': __version__,
                       'created_utc': datetime.now(UTC).isoformat()},
        'limitations': list(TABLE_LIMITATIONS),
    }


def missing_targets(table, sample_id):
    """The targets with no allele for one isolate, in the table's own row order."""
    if not any(row['sample_id'] == sample_id for row in table['samples']):
        raise ValueError(f'{sample_id} is not one of the isolates in this table of calls.')
    return [row['locus'] for row in table['loci']
            if not CALL_STATES[row['cells'][sample_id]['state']]['has_allele']]


def called_targets(table, sample_id):
    """The targets that matched a reference allele exactly, in row order.

    Only exact reference calls qualify: a control target has to have its own
    sequence in the scheme for the read assay to be able to recover it.
    """
    if not any(row['sample_id'] == sample_id for row in table['samples']):
        raise ValueError(f'{sample_id} is not one of the isolates in this table of calls.')
    return [row['locus'] for row in table['loci'] if row['cells'][sample_id]['state'] == 'called']


def _spread(values, count):
    """Pick `count` entries spread evenly across `values`, deterministically.

    Controls taken from one alphabetical corner of a scheme can all sit in one
    region of the genome. Spreading them is a better test that the reads reach
    the genome at all.
    """
    if count <= 0 or not values:
        return []
    if count >= len(values):
        return list(values)
    step = len(values) / count
    return [values[min(len(values) - 1, int(index * step))] for index in range(count)]


def plan_target_recheck(table, sample_id, *, loci=None, controls=None, control_count=2,
                        max_panel=MAX_PANEL_LOCI, sample=None):
    """Choose one bounded panel: uncalled targets to ask about, plus positive controls.

    The panel obeys read_support's own locus bound. When more targets have no call
    than fit, the batch and the deferred remainder are both named in the plan;
    nothing is dropped quietly. Pass the project `sample` record to carry its
    attached FASTQ pair and assembly identity into the plan.
    """
    if isinstance(control_count, bool) or not isinstance(control_count, int) or control_count < 0:
        raise ValueError('The number of control targets must be a non-negative integer.')
    if isinstance(max_panel, bool) or not isinstance(max_panel, int) or not 1 <= max_panel <= MAX_PANEL_LOCI:
        raise ValueError(f'A read-check panel holds between 1 and {MAX_PANEL_LOCI} targets, the bound the read '
                         'assay enforces.')
    summary = next((row for row in table['samples'] if row['sample_id'] == sample_id), None)
    if summary is None:
        raise ValueError(f'{sample_id} is not one of the isolates in this table of calls.')
    without_call = missing_targets(table, sample_id)
    if loci is None:
        targets = list(without_call)
    else:
        requested = [str(locus) for locus in loci]
        if len(set(requested)) != len(requested):
            raise ValueError('Each target appears once in a read check. Remove the repeated selections.')
        unknown = [locus for locus in requested if locus not in {row['locus'] for row in table['loci']}]
        if unknown:
            raise ValueError('These targets are not in this scheme: ' + ', '.join(sorted(unknown)))
        already = [locus for locus in requested if locus not in set(without_call)]
        if already:
            raise ValueError('These targets already carry a call for this isolate, so there is nothing to recheck: '
                             + ', '.join(sorted(already)) + '. Use them as controls instead.')
        order = {row['locus']: index for index, row in enumerate(table['loci'])}
        targets = sorted(requested, key=lambda locus: order[locus])
    if not targets:
        raise ValueError(f'{summary["sample_name"]} called every one of the {table["target_count"]:,} targets in '
                         f'{table["scheme"]}, so there is nothing to recheck against the reads.')
    available = called_targets(table, sample_id)
    if controls is None:
        chosen = _spread(available, min(control_count, max_panel - 1))
    else:
        chosen = [str(locus) for locus in controls]
        if len(set(chosen)) != len(chosen):
            raise ValueError('Each control target appears once. Remove the repeated selections.')
        invalid = [locus for locus in chosen if locus not in set(available)]
        if invalid:
            raise ValueError('A control target must be one this isolate called exactly against the scheme, so the '
                             'reads can be asked to recover a sequence that is in the scheme. Not usable: '
                             + ', '.join(sorted(invalid)))
    chosen = [locus for locus in chosen if locus not in set(targets)][:max(0, max_panel - 1)]
    room = max_panel - len(chosen)
    batch, deferred = targets[:room], targets[room:]
    order = {row['locus']: index for index, row in enumerate(table['loci'])}
    panel = sorted(set(batch) | set(chosen), key=lambda locus: order[locus])
    reads, expected, assembly = [], [], None
    if sample is not None:
        attached = sorted(((sample.get('metadata') or {}).get('reads') or {}).get('reads') or [],
                          key=lambda row: row.get('mate', 0))
        if len(attached) != 2 or [row.get('mate') for row in attached] != [1, 2] \
                or not all(_sha256(row.get('sha256')) and row.get('path') for row in attached):
            raise ValueError('Attach a verified original FASTQ pair to this isolate before rechecking its targets '
                             'against the reads. The assembly alone cannot answer whether a target is present.')
        reads = [str(row['path']) for row in attached]
        expected = [_sha256(row['sha256']) for row in attached]
        assembly = sample.get('input_path')
    control_basis = ('Controls are targets this isolate called exactly. If the reads recover them, the assay works '
                     'on this read set, which is what makes a negative result readable at all.')
    if not chosen:
        control_basis = ('No control target could be included: this isolate has no exact reference call to test the '
                         'assay against. Every "no reads" result in this run stays inconclusive.')
    batch_note = ''
    if deferred:
        batch_note = (f'{len(batch):,} of {len(targets):,} uncalled targets are in this batch, with {len(chosen)} '
                      f'control target(s); the read assay is bounded at {max_panel} targets per run. The remaining '
                      f'{len(deferred):,} are listed in this plan and are not checked until you run them.')
    return {
        'format_version': 1, 'sample_id': sample_id, 'sample_name': summary['sample_name'],
        'scheme': table['scheme'], 'scheme_digest': table['scheme_digest'],
        'target_count': table['target_count'],
        'targets': batch, 'deferred_targets': deferred, 'controls': chosen, 'panel_loci': panel,
        'uncalled_total': len(without_call), 'control_basis': control_basis, 'batch_note': batch_note,
        'read_paths': reads, 'expected_read_sha256': expected, 'assembly_path': assembly,
        'limits': {'max_panel': max_panel, 'requested_control_count': control_count},
        'limitations': list(RECHECK_LIMITATIONS),
    }


def _control_verdicts(support, plan):
    """Did the assay recover targets this assembly did call, from these reads?"""
    by_locus = {row['locus']: row for row in support.get('loci') or ()}
    rows = []
    for locus in plan['controls']:
        status = str(by_locus.get(locus, {}).get('status') or 'not_tested')
        passed = status in {'supported', 'ambiguous', 'mixed_support'}
        rows.append({
            'locus': locus, 'support_status': status, 'passed': passed,
            'explanation': ('The reads recovered this called target, so the assay reached the genome on this read set.'
                            if status == 'supported' else
                            'The reads recovered this called target, but more than one reference allele is compatible '
                            'with them.' if status == 'ambiguous' else
                            'The reads recovered this called target with disagreeing bases; the assay reached the '
                            'genome, and this isolate may hold more than one template.' if status == 'mixed_support' else
                            'The reads did not recover a target this assembly called exactly. Until that is explained '
                            '- wrong read pair, too few pairs read, or coverage - no negative result can be read here.'),
        })
    passed = sum(row['passed'] for row in rows)
    complete = bool((support.get('sampling') or {}).get('complete_files'))
    if not rows:
        reason = ('No control target was included in this panel, so nothing shows that the assay works on this read '
                  'set. Every "no reads" result below stays inconclusive.')
    elif not passed:
        reason = (f'None of the {len(rows)} control target(s) this assembly called was recovered from these reads. '
                  'The assay is not shown to work on this read set, so no absence can be read from it.')
    elif not complete:
        reason = (f'{passed} of {len(rows)} control target(s) were recovered, but only the first '
                  f'{(support.get("sampling") or {}).get("pairs_checked", 0):,} read pairs were examined. Absence '
                  'cannot be established from part of a file; re-run over the complete pair to ask that question.')
    else:
        reason = (f'{passed} of {len(rows)} control target(s) this assembly called were recovered from the complete '
                  'read files, so a target with no reads at all can be read as consistent with absence.')
    return {'loci': rows, 'tested': len(rows), 'passed': passed, 'complete_files': complete,
            'absence_readable': bool(passed and complete), 'reason': reason}


def interpret_target_recheck(support, plan, *, table=None):
    """Turn one read-support payload into a present / absent / not-enough answer per target.

    This is the honest heart of the recheck and holds no file or tool access, so
    every verdict boundary is testable on its own.
    """
    if str(support.get('scheme_digest') or '') != str(plan.get('scheme_digest') or ''):
        raise ValueError('This read investigation was run against a different scheme than the plan it is being read '
                         'against; no verdict is produced from the two together.')
    by_locus = {row['locus']: row for row in support.get('loci') or ()}
    absent_from_run = [locus for locus in plan['panel_loci'] if locus not in by_locus]
    if absent_from_run:
        raise ValueError('The read investigation did not cover every target in the plan (' +
                         ', '.join(sorted(absent_from_run)[:5]) + '). No partial verdict is produced.')
    controls = _control_verdicts(support, plan)
    sampling = support.get('sampling') or {}
    states = {}
    if table is not None:
        states = {row['locus']: row['cells'][plan['sample_id']] for row in table['loci']
                  if plan['sample_id'] in row['cells']}
    rows, tally = [], {'present': 0, 'absent': 0, 'insufficient': 0}
    for locus in plan['targets']:
        evidence = by_locus[locus]
        status = str(evidence.get('status') or '')
        note = ''
        if status in {'supported', 'ambiguous'}:
            verdict = 'present_in_reads'
            if status == 'ambiguous':
                note = ('More than one reference allele is compatible with these reads, so no single allele is '
                        'even a candidate here.')
        elif status == 'mixed_support':
            verdict = 'present_with_mixed_bases'
        elif status == 'incomplete_or_discordant':
            verdict = 'partial_read_evidence'
        elif status != 'no_support':
            # Only an explicit no-support result can ever be read towards absence.
            # A state this version does not know is not evidence of anything.
            verdict = 'no_evidence_either_way'
            note = (f'The read assay reported state {status or "(absent)"!r}, which this version does not '
                    'recognise, so nothing is read from it in either direction.')
        elif controls['absence_readable']:
            verdict = 'consistent_with_absence'
        else:
            verdict, note = 'no_evidence_either_way', controls['reason']
        meaning = (RECHECK_VERDICTS[verdict]['meaning'] + ' ' + note).strip()
        cell = states.get(locus) or {}
        candidates = list(evidence.get('compatible_candidates') or ())
        best = max((row for row in evidence.get('candidates') or ()),
                   key=lambda row: (row.get('breadth_min_depth') or 0, row.get('supporting_reads') or 0),
                   default={})
        row = {
            'locus': locus, 'verdict': verdict, 'evidence': RECHECK_VERDICTS[verdict]['evidence'],
            'label': RECHECK_VERDICTS[verdict]['label'], 'explanation': meaning,
            'support_status': status, 'assembly_state': cell.get('state'),
            'assembly_state_label': CALL_STATES[cell['state']]['label'] if cell.get('state') else '',
            'assigned_allele': None,
            'compatible_candidates': candidates,
            'candidate_note': (f'{len(candidates):,} reference allele(s) are compatible with these reads. '
                               'Compatible is not called: no allele has been assigned to this target.'
                               if candidates else 'No reference allele of this target is compatible with these reads.'),
            'best_breadth': best.get('breadth_min_depth'), 'best_mean_depth': best.get('mean_depth'),
            'supporting_reads': best.get('supporting_reads', 0),
            'read_support': evidence,
        }
        tally[row['evidence']] += 1
        rows.append(row)
    return {
        'format_version': 1, 'sample_id': plan['sample_id'], 'sample_name': plan['sample_name'],
        'scheme': plan['scheme'], 'scheme_digest': plan['scheme_digest'],
        'target_count': plan['target_count'],
        'targets': rows, 'controls': controls,
        'deferred_targets': list(plan.get('deferred_targets') or ()),
        'summary': {**tally, 'checked': len(rows), 'uncalled_total': plan.get('uncalled_total', len(rows)),
                    'verdicts': {key: sum(row['verdict'] == key for row in rows) for key in RECHECK_VERDICTS}},
        'headline': (f'{tally["present"]:,} present in the reads · {tally["absent"]:,} consistent with absence · '
                     f'{tally["insufficient"]:,} not enough evidence, of {len(rows):,} rechecked targets '
                     f'({plan.get("uncalled_total", len(rows)):,} have no call in total).'),
        'sampling': sampling, 'panel': support.get('panel'), 'parameters': support.get('parameters'),
        'reads': support.get('reads'), 'input_sha256': support.get('input_sha256'),
        'provenance': support.get('provenance'),
        'verdict_meanings': {key: dict(value) for key, value in RECHECK_VERDICTS.items()},
        'limitations': list(RECHECK_LIMITATIONS) + list(READ_SUPPORT_LIMITATIONS),
    }


def recheck_missing_targets(plan, scheme, cancelled=None, progress=None, *, table=None,
                            max_pairs=100_000, **options):
    """Run the bounded read assay over a planned panel and return the read verdicts.

    Cancellation and progress pass straight through to the read investigation, which
    is the slow part; nothing is written and no call is changed.
    """
    check_cancelled(cancelled)
    if len(plan.get('read_paths') or ()) != 2 or len(plan.get('expected_read_sha256') or ()) != 2:
        raise ValueError('This read check needs the isolate\'s verified original FASTQ pair. Attach the reads to the '
                         'isolate and plan the check again.')
    if isinstance(scheme, Scheme) and scheme.digest != plan['scheme_digest']:
        raise ValueError('The loaded scheme is not the one these calls were made against, so a read check against it '
                         'would answer a different question.')
    if progress:
        progress(0, 0, f'Rechecking {len(plan["targets"]):,} uncalled target(s) and {len(plan["controls"]):,} '
                       f'control target(s) for {plan["sample_name"]} against the reads')
    support = investigate_read_support(
        plan['read_paths'][0], plan['read_paths'][1], scheme, plan['panel_loci'],
        cancelled=cancelled, progress=progress, assembly_path=plan.get('assembly_path'),
        expected_read_sha256=plan['expected_read_sha256'], max_pairs=max_pairs, **options)
    check_cancelled(cancelled)
    result = interpret_target_recheck(support, plan, table=table)
    if progress:
        progress(1, 1, result['headline'])
    return result
