"""Extended Klebsiella typing: virulence locus STs and wzi/wzc capsule markers.

Both assays run through the same exact-allele engine the MLST path uses, over
public reference data pinned by commit: Kleborate's allele and profile tables
and Kaptive v2.0.9's wzi/wzc marker database. Neither assay is Kleborate or
Kaptive, and neither is equivalent to them. No hypervirulence phenotype, no
aggregate virulence or resistance score, and no K locus or capsule serotype is
inferred: a wzi allele number is reported as a wzi allele number.
"""

from __future__ import annotations

import csv
import html
import threading
from pathlib import Path

from .characterization_refs import validate_characterization_references
from .organism_modules import OrganismMatch, OrganismModule, register
from .sequence import AnalysisCancelled, SequenceError, check_cancelled, file_sha256, file_signature
from .typing import SchemeError, call_assembly, load_scheme

LOCUS_ST_LIMITATIONS = [
    'Exact nucleotide allele matching only: an inexact or novel allele yields an incomplete call and no '
    'locus ST, never the nearest ST.',
    'The profile tables are a pinned snapshot. A locus ST registered upstream after the pinned revision is '
    'absent here, and a profile referencing an allele with no reference sequence cannot be assigned.',
    'The lineage value is a published Kleborate profile-table lookup reported verbatim, not an inference by '
    'this software.',
    'A locus ST is not a hypervirulence phenotype, not plasmid identity and not a Kleborate virulence or '
    'resistance score; those aggregate scores are not computed here.',
    'This is a WMLSTudio exact-allele call using Kleborate public allele and profile data. It is not '
    'Kleborate output and is not equivalent to Kleborate.',
    'rmpA2 has no upstream profile table and therefore never acquires a locus ST; it stays allele-only in '
    'the virulence screen.',
]
CAPSULE_LIMITATIONS = [
    'The result is a wzi or wzc allele number and nothing else. The published wzi-allele to K-type '
    'associations are neither shipped nor applied: no K locus, capsule type or serotype is inferred.',
    'This is a WMLSTudio exact-allele screen over the public Kaptive wzi/wzc marker database. It is not '
    'Kaptive output and is not equivalent to Kaptive: the whole K/O locus references, Kaptive match-'
    'confidence grading and the O-locus special logic are not implemented here.',
    'wzc entries are short (115-151 bp) variable-region fragments and are more prone to several exact '
    'matches; ambiguous and mixed calls are surfaced unchanged rather than tie-broken.',
    'Exact nucleotide matching only: a novel or inexact marker sequence is reported as missing, never as '
    'the nearest allele.',
    'No capsule marker allele establishes virulence, transmission or a phenotype.',
]
KLEBSIELLA_MATCH = OrganismMatch(genera=frozenset({'Klebsiella'}))
# The ybt allele set alone is about 11 MB, so exactly one locus scheme is held at
# a time; the automaton for the previous locus is released before the next build.
_CACHE = {}
_LOCK = threading.Lock()


def _scheme(directory, cancelled=None):
    directory = Path(directory).resolve()
    key = str(directory)
    with _LOCK:
        cached = _CACHE.get(key)
        if cached and all(file_signature(path) == signature for path, signature in cached[1]):
            return cached[0]
    scheme = load_scheme(directory, cancelled)
    signatures = [(path, file_signature(path)) for path in sorted(directory.iterdir()) if path.is_file()]
    with _LOCK:
        _CACHE.clear()  # One bounded scheme at a time; allele sets here are large.
        _CACHE[key] = (scheme, signatures)
    return scheme


def _lineage_table(path, entry):
    field = entry.get('lineage_field')
    if not field:
        return {}
    table = {}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            value = (row.get(field) or '').strip()
            key = (row.get(entry.get('st_field') or 'ST') or '').strip()
            if key and value:
                table[key] = value
    return table


def type_virulence_loci(path, reference_root, cancelled=None, progress=None, *, threads=2):
    """Assign Kleborate virulence locus STs through the exact-allele engine.

    ``threads`` is accepted for a uniform assay signature; exact allele calling is
    single-threaded and does not use it. A locus ST appears only when every locus
    gene produced one unambiguous exact match and the complete vector maps to
    exactly one ST in the pinned profile table.
    """
    root = Path(reference_root).resolve()
    check_cancelled(cancelled)
    manifest = validate_characterization_references(root, cancelled=cancelled)
    section = manifest.get('locus_profiles') or {}
    if not section:
        return {'status': 'not_run',
                'reason': 'This reference snapshot carries no Kleborate locus-ST profile tables. Install an '
                          'updated characterization snapshot.'}
    path = Path(path).resolve()
    signature = file_signature(path)
    input_hash = file_sha256(path, cancelled)
    loci, failed, unresolved = {}, False, False
    for index, (locus, entry) in enumerate(sorted(section.items()), 1):
        check_cancelled(cancelled)
        if progress:
            progress(index, len(section), f'Calling exact {locus} alleles against the pinned Kleborate profile table…')
        profile_path = root / entry['path']
        try:
            scheme = _scheme(profile_path.parent, cancelled)
            call = call_assembly(path, scheme, cancelled)
        except AnalysisCancelled:
            raise
        except (SchemeError, SequenceError, ValueError, OSError) as error:
            loci[locus] = {'locus_st': None, 'status': 'failed', 'reason': str(error)}
            failed = True
            continue
        locus_st = call['st'] if call['status'] == 'complete' else None
        unresolved = unresolved or call['status'] in {'ambiguous', 'mixed'}
        loci[locus] = {'locus_st': locus_st, 'status': call['status'], 'alleles': call['alleles'],
                       'calls': call['calls'], 'notes': call['notes'], 'scheme_digest': call['scheme_digest'],
                       'lineage': _lineage_table(profile_path, entry).get(locus_st) if locus_st else None,
                       'lineage_field': entry.get('lineage_field'),
                       'lineage_source': 'Kleborate profile table, reported verbatim; not inferred by WMLSTudio',
                       'module': entry.get('module'), 'phenotype': 'not_inferred'}
    if file_signature(path) != signature or file_sha256(path, cancelled) != input_hash:
        raise ValueError('Assembly changed during locus-ST calling; no result was accepted.')
    status = 'failed' if failed else ('ambiguous' if unresolved else 'completed')
    return {'status': status, 'input_sha256': input_hash, 'reference_digest': manifest['reference_digest'],
            'loci': loci, 'phenotype': 'not_inferred',
            'official_type': None,
            'assay': {'name': 'WMLSTudio Kleborate-profile virulence locus ST call',
                      'method': 'exact-nucleotide allele matching, complete vector to a single pinned ST',
                      'loci': sorted(section), 'profile_source': 'kleborate',
                      'not_typed': ['rmpA2']},
            'limitations': list(LOCUS_ST_LIMITATIONS)}


def locus_st_assignments(evidence):
    """Locus STs this assay actually assigned, in the shape the virulence screen accepts."""
    assignments = {}
    for locus, entry in ((evidence or {}).get('loci') or {}).items():
        if entry.get('locus_st'):
            assignments[locus] = {'locus_st': entry['locus_st'],
                                  'source': 'WMLSTudio exact-allele call against the pinned Kleborate '
                                            f'{locus} profile table'}
    return assignments


def type_capsule_markers(path, reference_root, cancelled=None, progress=None, *, threads=2):
    """Call wzi and wzc marker alleles exactly; no K locus is inferred from them.

    With no profile table in the capsule panel the engine correctly reports
    ``profile_unavailable`` and never invents a sequence type.
    """
    root = Path(reference_root).resolve()
    check_cancelled(cancelled)
    manifest = validate_characterization_references(root, cancelled=cancelled)
    section = manifest.get('capsule') or {}
    entries = section.get('loci') or []
    if not entries:
        return {'status': 'not_run',
                'reason': 'This reference snapshot carries no wzi/wzc capsule marker panel. Install an '
                          'updated characterization snapshot.'}
    directories = {(root / entry['path']).parent for entry in entries}
    if len(directories) != 1:
        raise ValueError('Capsule marker references must share one panel directory.')
    path = Path(path).resolve()
    if progress:
        progress(0, 0, 'Calling exact wzi and wzc capsule marker alleles…')
    scheme = _scheme(directories.pop(), cancelled)
    call = call_assembly(path, scheme, cancelled)
    markers = {}
    for row in call['calls']:
        markers[row['locus']] = {'allele': row['allele'], 'status': row['status'],
                                 'candidates': row['candidates'], 'reason': row['reason']}
    status = 'ambiguous' if call['status'] in {'ambiguous', 'mixed'} else 'completed'
    result = {'status': status, 'input_sha256': call['input_sha256'],
              'reference_digest': manifest['reference_digest'], 'markers': markers,
              'scheme_status': call['status'], 'scheme_digest': call['scheme_digest'],
              'st': None, 'k_locus': None, 'k_locus_status': 'not_assayed',
              'k_locus_reason': 'wzi-allele to K-type associations are not shipped and not applied by this '
                                'software; only the marker allele number is reported.',
              'phenotype': 'not_inferred', 'notes': call['notes'],
              'assay': {'name': 'WMLSTudio Kaptive-reference wzi/wzc capsule marker screen',
                        'method': 'exact-nucleotide allele matching; no profile table, so no sequence type',
                        'loci': list(scheme.loci), 'allele_counts': {entry['gene']: entry.get('allele_count')
                                                                     for entry in entries},
                        'marker_source': 'kaptive'},
              'limitations': list(CAPSULE_LIMITATIONS)}
    for gene in ('wzi', 'wzc'):
        result[gene] = markers.get(gene) or {'allele': None, 'status': 'not_assayed',
                                             'reason': 'This marker is not present in the staged panel.'}
    return result


def _escape(value):
    return html.escape(str(value if value is not None else ''))


def summarize_locus_sts(evidence):
    status = (evidence or {}).get('status')
    if status in (None, 'not_run'):
        return 'not_run'
    if status == 'failed':
        return 'failed'
    parts = [f'{locus} {entry["locus_st"]}' for locus, entry in sorted((evidence.get('loci') or {}).items())
             if entry.get('locus_st')]
    if parts:
        return '; '.join(parts)
    return 'no locus ST assigned'


def summarize_capsule(evidence):
    status = (evidence or {}).get('status')
    if status in (None, 'not_run'):
        return 'not_run'
    if status == 'failed':
        return 'failed'
    parts = []
    for gene in ('wzi', 'wzc'):
        marker = evidence.get(gene) or {}
        parts.append(f'{gene} {marker["allele"]}' if marker.get('allele') else f'{gene} {marker.get("status", "missing")}')
    return '; '.join(parts)


def locus_st_html(evidence):
    evidence = evidence or {}
    status = evidence.get('status', 'not_run')
    parts = ['<h4>Klebsiella virulence locus STs (Kleborate reference data)</h4>',
             '<p><b>Result:</b> ' + _escape(summarize_locus_sts(evidence)) + ' &middot; status ' + _escape(status) + '</p>']
    if status in {'not_run', 'failed'}:
        parts.append('<p>' + _escape(evidence.get('reason') or 'This assay was not run for this isolate.') + '</p>')
    else:
        parts.append('<p><b>Applicability:</b> ' + _escape(evidence.get('applicability') or 'unknown_organism') +
                     ' &middot; ' + _escape(evidence.get('applicability_reason') or '') + '</p>')
        rows = ''.join('<tr><td>' + _escape(locus) + '</td><td>' + _escape(entry.get('locus_st') or 'not assigned') +
                       '</td><td>' + _escape(entry.get('status')) + '</td><td>' +
                       _escape(entry.get('lineage') or 'not reported') + '</td></tr>'
                       for locus, entry in sorted((evidence.get('loci') or {}).items()))
        parts.append('<table><tr><th>Locus</th><th>Locus ST</th><th>Call status</th><th>Lineage (table lookup)</th></tr>'
                     + rows + '</table>')
    for limitation in evidence.get('limitations') or LOCUS_ST_LIMITATIONS:
        parts.append('<p class="muted">' + _escape(limitation) + '</p>')
    return ''.join(parts)


def capsule_html(evidence):
    evidence = evidence or {}
    status = evidence.get('status', 'not_run')
    parts = ['<h4>Klebsiella wzi / wzc capsule markers (Kaptive reference data)</h4>',
             '<p><b>Result:</b> ' + _escape(summarize_capsule(evidence)) + ' &middot; status ' + _escape(status) + '</p>']
    if status in {'not_run', 'failed'}:
        parts.append('<p>' + _escape(evidence.get('reason') or 'This assay was not run for this isolate.') + '</p>')
    else:
        parts.append('<p><b>Applicability:</b> ' + _escape(evidence.get('applicability') or 'unknown_organism') +
                     ' &middot; ' + _escape(evidence.get('applicability_reason') or '') + '</p>')
        rows = ''.join('<tr><td>' + _escape(gene) + '</td><td>' + _escape((evidence.get(gene) or {}).get('allele') or 'not assigned') +
                       '</td><td>' + _escape((evidence.get(gene) or {}).get('status')) + '</td></tr>'
                       for gene in ('wzi', 'wzc'))
        parts.append('<table><tr><th>Marker</th><th>Allele</th><th>Call status</th></tr>' + rows + '</table>')
        parts.append('<p><b>K locus:</b> not assayed &middot; ' + _escape(evidence.get('k_locus_reason') or '') + '</p>')
    for limitation in evidence.get('limitations') or CAPSULE_LIMITATIONS:
        parts.append('<p class="muted">' + _escape(limitation) + '</p>')
    return ''.join(parts)


LOCUS_ST_MODULE = register(OrganismModule(
    key='klebsiella_locus_st', title='Klebsiella virulence locus STs (ybt, clb, iuc, iro, rmp)',
    column_title='Locus STs', match=KLEBSIELLA_MATCH, runner=type_virulence_loci,
    option_keys=('threads',), manifest_sections=('locus_profiles',),
    summary=summarize_locus_sts, detail_html=locus_st_html, report_default=True,
    purpose='Matches the ybt, clb, iuc, iro and rmp virulence loci against published Kleborate profile tables '
            'and reports a locus sequence type only when every allele matches exactly. It does not establish '
            'hypervirulence, plasmid identity or a Kleborate virulence score, and the lineage shown is a table '
            'lookup rather than an inference.',
    locus_st_provider=locus_st_assignments, limitations=tuple(LOCUS_ST_LIMITATIONS)))

CAPSULE_MODULE = register(OrganismModule(
    key='klebsiella_capsule', title='Klebsiella wzi / wzc capsule markers',
    column_title='Capsule markers', match=KLEBSIELLA_MATCH, runner=type_capsule_markers,
    option_keys=('threads',), manifest_sections=('capsule.loci',),
    summary=summarize_capsule, detail_html=capsule_html, report_default=True,
    purpose='Matches the wzi and wzc capsule marker genes against the public Kaptive marker database and '
            'reports the exact allele number of each. It does not establish a K locus, a capsule type or a '
            'serotype: no wzi-allele to K-type mapping is shipped or applied.',
    limitations=tuple(CAPSULE_LIMITATIONS)))
