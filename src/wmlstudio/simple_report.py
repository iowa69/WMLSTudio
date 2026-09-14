"""A short, plain-language outbreak summary for readers who are not bioinformaticians.

The renderer is pure and Qt-free: the caller supplies the project records, an
already-built comparison snapshot and an already-rasterised picture, so the
document can be produced and checked without a display. Nothing here infers
transmission, a resistance gene is never presented as a measured susceptibility
result, and a pair that could not be compared is never reported as a distance of
zero. Measured on A4 at the 144-dpi report resolution, the page runs to two
sheets for a handful of isolates and three for ten; the caveats are not dropped
to reach a single sheet. ``wmlstudio.export`` imports this module lazily, so the
import below cannot close a cycle.
"""

from __future__ import annotations

import base64
from datetime import datetime, timezone
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.export import REPORT_PRESETS, _escape, _snapshot, investigation_document

# A different document shape, not a different set of facts: the same snapshot
# rows feed the existing presets. Only the layout and the wording change. The
# registry in export.py is the one definition; this is the name the renderer and
# its tests reach for, copied so a caller cannot edit the registry through it.
SIMPLE_REPORT_PRESET = dict(REPORT_PRESETS['simple'])

# Escape-stable: typographic quotes keep the constant identical in the document,
# so no option flag and no escaping pass can quietly reword it.
SUSCEPTIBILITY_CAVEAT = (
    'Genes found in a genome are not measured antibiotic susceptibility. An isolate carrying a '
    'resistance gene may still test susceptible in the laboratory, and an isolate with no gene found '
    'may still test resistant. Use validated laboratory susceptibility testing for any treatment or '
    'infection-control decision. “Not assessed” and “unknown” are never “susceptible”.'
)

LIMITATIONS = (
    'It does not prove transmission, or its direction. Genomic similarity is one line of evidence for '
    'epidemiological review.',
    'Resistance genes are not measured susceptibility. Use laboratory testing for treatment decisions.',
    'Loci that could not be called are unknown, not identical — the shared / total column shows how much '
    'was actually compared.',
    'Isolates are grouped by single linkage, so members of one group can differ by more than the threshold.',
    'Sections you switched off, and assays that were not run, are absent from this report — absence here '
    'is not a negative result.',
)

NO_COMPARISON = (
    'No comparison has been built for these isolates, so no distances were calculated and no picture '
    'could be drawn. This is not a statement that the isolates are unrelated.'
)

COMPARISON_OFF = (
    'The comparison section was switched off for this report. That is not a statement that the isolates '
    'are unrelated or that no distances exist.'
)

_STATE_WORDS = {
    'current': 'Checked against this exact sequence file',
    'unverified': 'Linked by name, not verified',
    'stale': 'Out of date',
    'missing': 'Not run',
}

_ORGANISM_SOURCE_WORDS = {
    'assigned': 'you assigned this organism',
    'local_scheme_detection': 'suggested by the local typing scheme',
    'unknown': 'organism not determined',
}

_GROUP_WORDS = {
    'cluster': 'isolates linked at this threshold',
    'singleton': 'no other isolate within the threshold',
    'not_comparable': 'no usable comparison',
}

# The same restricted rich-text subset the other reports print with: attribute
# table borders, class selectors, and a data: image source the CSP already allows.
_HEAD_OPEN = (
    '<!doctype html><html><head><meta charset="utf-8">'
    '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; style-src \'unsafe-inline\'; img-src data:">'
    '<title>'
)
_HEAD_CLOSE = (
    '</title>'
    '<style>body{font-family:Segoe UI,Arial,sans-serif;color:#203c45;max-width:1080px;margin:auto;padding:28px;line-height:1.3;font-size:12px}'
    'h1,h2{color:#156e65}h2{border-bottom:2px solid #d6e7e3;padding-bottom:8px}'
    'table{border-collapse:collapse;width:100%;font-size:11px}td,th{border:1px solid #d6e0e4;padding:5px;text-align:left;vertical-align:top}'
    'th{background:#edf5f3}.notice{padding:12px;background:#fff4dd;border-left:4px solid #c89031}'
    '.muted{color:#5c6f76}img{max-width:100%;height:auto}'
    '@media print{body{padding:0}tr{break-inside:avoid}h2{break-after:avoid}}</style></head><body>'
)

# Measured, not guessed: QTextDocument.print_ lays the page out in 96-dpi units,
# so an A4 sheet at the 144-dpi report resolution leaves about 660 usable units.
# A wider picture is silently cut off at the right edge of the printed sheet.
_IMAGE_WIDTH = 640


def _muted(text: str) -> str:
    return '<br><span class="muted">' + _escape(text) + '</span>'


def _differences(value) -> str:
    """Plain wording for a count a reader will read aloud in a meeting."""
    return _escape(value) + (' allele difference' if value == 1 else ' allele differences')


def _row(cells) -> str:
    """Cells are already-escaped markup; the callers escape every user value."""
    return '<tr>' + ''.join('<td>' + cell + '</td>' for cell in cells) + '</tr>'


def _image_source(image, mime) -> str:
    """Inline the caller's picture so the page stays offline and Qt-free here.

    ``image`` is raw bytes, a ``data:`` URI, or a readable file path: the
    rasterisation itself belongs to the Qt layer that drew the graph.
    """
    if not image:
        return ''
    if isinstance(image, (bytes, bytearray, memoryview)):
        data = bytes(image)
    else:
        text = str(image)
        if text.startswith('data:'):
            return _escape(text)
        try:
            data = Path(text).expanduser().read_bytes()
        except OSError:
            return ''
    if not data:
        return ''
    return 'data:' + _escape(mime or 'image/png') + ';base64,' + base64.b64encode(data).decode('ascii')


def _scheme_scale(snapshot) -> tuple[str, int]:
    """Name the reference that was actually used; never claim a core-genome scale."""
    profiles = snapshot.get('profiles') or []
    loci = max((profile.get('total_loci') or 0) for profile in profiles) if profiles else 0
    return str(snapshot.get('scheme') or 'Unnamed reference'), loci


def _catalog_lines(rows) -> list[str]:
    """Published cutoffs that exist for this organism, explicitly none of them applied."""
    from wmlstudio.threshold_guidance import guidance_for
    organisms = {str(row.get('organism') or '').strip() for row in rows}
    organisms.discard('')
    if len(organisms) != 1:
        return []
    organism = organisms.pop()
    guidance = guidance_for(organism)
    if not guidance['entries']:
        return ['<p class="muted">' + _escape(guidance['message']) + '</p>']
    cited = []
    for entry in guidance['entries'][:2]:
        value = entry.get('published_threshold')
        cited.append(('no numeric cutoff curated' if value is None else
                      'at most ' + _escape(value) + ' ' + _escape(entry.get('unit') or 'allele differences')) +
                     ' (' + _escape(entry['source']['citation']) + ')')
    return ['<p class="muted">Published cutoffs exist for ' + _escape(organism) + ': ' + '; '.join(cited) +
            '. <b>None of them is applied to this report.</b> ' + _escape(guidance['interpretation']) + '</p>']


def _threshold_lines(snapshot, rows) -> list[str]:
    """State where the threshold came from, or that nothing published backs it."""
    from wmlstudio.investigation import threshold_guidance_status
    binding = threshold_guidance_status(snapshot)
    evidence = snapshot.get('threshold_evidence') or {}
    citation = _escape(evidence.get('citation') or 'No citation recorded')
    doi = _escape(evidence.get('doi') or 'no DOI recorded')
    status = binding['status']
    if status == 'no_guidance':
        stated = ('no published cutoff is attached to it. At most ' + _differences(snapshot.get('threshold')) +
                  ' is a local exploratory setting, not a validated rule.')
    elif status == 'citation_only':
        stated = 'cited for context only — ' + citation + ' (' + doi + '). No numeric cutoff was adopted from it.'
    elif status == 'matches_reviewed_context':
        stated = (citation + ' (' + doi + ') — published cutoff at most ' + _escape(evidence.get('published_threshold')) +
                  '; applied here as at most ' + _escape(evidence.get('approved_threshold')) +
                  '. Matching an organism or a locus count is not clinical validation.')
    else:
        stated = _escape(binding['reason']) + ' Citation retained: ' + citation + ' (' + doi + ').'
    parts = ['<p><b>Where this threshold comes from:</b> ' + stated + '</p>']
    # A catalog entry is context that exists, never a rule that was applied here.
    return (parts + _catalog_lines(rows)) if status == 'no_guidance' else parts


def _comparison_section(snapshot, rows, options, graph_png, graph_mime) -> list[str]:
    parts = ['<h2>How close are these isolates?</h2>',
             '<p>This shows how many allele differences separate the isolates that were compared. '
             'It does not show who infected whom, in which direction, or when.</p>']
    scheme, loci = _scheme_scale(snapshot)
    if loci and loci < 100:
        parts.append('<p class="notice"><b>This comparison used a small, classical MLST-scale reference (' + str(loci) +
                     ' loci).</b> A reference this size cannot separate isolates within one outbreak: two unrelated '
                     'isolates of the same sequence type look identical here. Core-genome references compare hundreds '
                     'to thousands of loci and are installed separately.</p>')
    source = _image_source(graph_png, graph_mime) if options['graph'] else ''
    if source:
        parts.append('<img width="' + str(_IMAGE_WIDTH) + '" src="' + source +
                     '" alt="Allele-distance minimum spanning forest with the report isolates highlighted">')
        parts.append('<p class="muted">The picture is a minimum spanning forest of allele differences: each line joins '
                     'an isolate to a close comparable neighbour and the number on it is the count of differing loci. '
                     'It is not a phylogeny and not a transmission tree, and the positions carry no meaning.</p>')
    elif not options['graph']:
        parts.append('<p class="muted">The picture was switched off for this report. The numbers below still describe '
                     'the comparison.</p>')
    else:
        parts.append('<p class="muted">No picture is included here. The numbers below still describe the comparison; '
                     'a missing picture is not a statement about relatedness.</p>')
    parts.append('<p><b>Reference used:</b> ' + _escape(scheme) + ' · ' + _escape(loci) + ' loci per profile</p>')
    overlap = snapshot.get('min_overlap')
    share = f'{overlap:.0%}' if isinstance(overlap, (int, float)) and not isinstance(overlap, bool) else overlap
    parts.append('<p><b>Close</b> in this report means at most ' + _differences(snapshot.get('threshold')) +
                 ', counted only at loci both isolates could call, with at least ' +
                 _escape(share) + ' of the loci shared.</p>')
    parts.extend(_threshold_lines(snapshot, rows))
    if snapshot.get('preview'):
        parts.append('<p class="notice"><b>Unsaved threshold preview.</b> The saved protocol threshold was ' +
                     _escape(snapshot.get('saved_threshold')) + '; this report explicitly uses the threshold shown above.</p>')
    return parts


def _gene_cell(row) -> str:
    """Stale evidence is withheld, and nothing found is never reported as susceptible."""
    status = row.get('hydra_evidence_status') or 'missing'
    genes = [str(gene) for gene in (row.get('amr_genes') or [])]
    classes = [str(name) for name in (row.get('amr_classes') or [])]
    if status == 'stale':
        return _escape('Not shown — the saved AMR result belongs to a different sequence file')
    if status == 'missing':
        return _escape('Not assessed')
    if genes:
        cell = _escape(', '.join(genes))
        if classes:
            cell += _muted('Drug classes named by the reference database: ' + ', '.join(classes))
        if status == 'unverified':
            cell += _muted('The source report’s identity was not confirmed.')
        return cell
    if status == 'unverified':
        return _escape('None reported (source identity not confirmed)')
    return _escape('No resistance determinants reported by the AMR database used')


def _resistance_section(rows, options) -> list[str]:
    parts = ['<h2>Resistance genes found</h2>',
             '<p>This lists the resistance genes detected in each genome and how well that detection matches the '
             'current sequence file. It does not show whether an isolate is resistant or susceptible to any drug.</p>']
    if not options['amr']:
        parts.append('<p class="notice">The resistance section was switched off for this report. That is not a '
                     'statement that no resistance genes exist. ' + _escape(SUSCEPTIBILITY_CAVEAT) + '</p>')
        return parts
    parts.append('<table border="1" cellpadding="5" cellspacing="0"><tr><th>Isolate</th><th>Organism</th>'
                 '<th>Sequence type</th><th>Resistance genes detected</th><th>Evidence state</th></tr>')
    for row in rows:
        isolate = _escape(row.get('sample_name') or row.get('sample_id'))
        if row.get('missing_input'):
            isolate += _muted('The original sequence file is not available; stored results only.')
        if row.get('error'):
            isolate += _muted(str(row['error']))
        organism = _escape(row.get('organism') or 'Not assigned')
        organism += _muted(_ORGANISM_SOURCE_WORDS.get(row.get('organism_source'), 'not determined'))
        status = row.get('hydra_evidence_status') or 'missing'
        state = _escape(_STATE_WORDS.get(status, status))
        if status != 'current' and row.get('hydra_evidence_reason'):
            state += _muted(str(row['hydra_evidence_reason']))
        parts.append(_row([isolate, organism, _escape(row.get('st') or row.get('primary_st') or 'Not assigned'),
                           _gene_cell(row), state]))
    # The caveat closes the table in one string: no option can separate them.
    parts.append('</table><p class="notice"><b>' + _escape(SUSCEPTIBILITY_CAVEAT) + '</b></p>')
    return parts


def _proximity_cells(focal, threshold) -> list[str]:
    distance = focal.get('nearest_distance')
    nearest = focal.get('nearest') or []
    names = [str(pair.get('neighbour_name') or pair.get('neighbour_id')) for pair in nearest]
    closest = (_escape(', '.join(names[:3])) + (_escape(f' +{len(names) - 3} more') if len(names) > 3 else '')
               if names else _escape('No comparable isolate in this comparison'))
    if nearest:
        shown = nearest[:3]
        denominators = [f"{pair['shared_loci']}/{pair['total_loci']}" for pair in shown]
        if len(set(denominators)) == 1:
            shared = _escape(denominators[0]) + _muted(f"{shown[0]['overlap']:.0%} of loci shared")
        else:
            # Equally close isolates can have been compared over different locus
            # sets. One figure for all of them would hide that difference.
            shared = _escape(' · '.join(denominators)) + _muted(
                'One denominator per closest isolate, in the order listed: they were not all compared '
                'over the same loci.')
    else:
        shared = '—'
    group = _escape(focal.get('group_name') or 'Not grouped')
    group += _muted(_GROUP_WORDS.get(focal.get('group_status'), 'no group recorded'))
    if focal.get('chained'):
        group += _muted('Joined through intermediate isolates.')
    comparable = distance is not None and isinstance(threshold, (int, float)) and not isinstance(threshold, bool)
    return [_escape(focal.get('sample_name') or focal.get('sample_id')),
            _escape('Not comparable') if distance is None else _escape(distance), closest, shared,
            _escape('Yes' if distance <= threshold else 'No') if comparable else '—', group]


def _proximity_section(snapshot, rows, names) -> list[str]:
    doc = investigation_document(snapshot, [row.get('sample_id') for row in rows])
    threshold = snapshot.get('threshold')
    parts = ['<h2>Closest matches</h2>',
             '<p>This shows each isolate’s closest comparable isolate and how much of the reference the two shared. '
             'It does not show a transmission link, and a locus that could not be called is never counted as a match.</p>',
             '<table border="1" cellpadding="5" cellspacing="0"><tr><th>Isolate</th><th>Allele differences to closest</th>'
             '<th>Closest isolate(s)</th><th>Shared / total loci</th><th>Within threshold?</th><th>Group</th></tr>']
    for focal in doc['proximity']:
        parts.append(_row(_proximity_cells(focal, threshold)))
    for sample_id in doc['outside_snapshot_ids']:
        parts.append(_row([_escape(names.get(sample_id, sample_id)),
                           _escape('Not in this comparison — no allele profile for this reference'),
                           '—', '—', '—', '—']))
    parts.append('</table>')
    reasons = sorted({str(pair.get('reason')) for focal in doc['proximity']
                      for pair in focal.get('excluded') or [] if pair.get('reason')})
    if reasons:
        parts.append('<p class="muted">Some pairs could not be compared at all: ' + _escape('; '.join(reasons)) +
                     ' A pair that cannot be compared has no distance; it is not a distance of zero.</p>')
    parts.append('<p class="muted">Distances count allele differences only at loci both isolates could call. Missing '
                 'calls are unknown, not identical. Isolates are grouped by single linkage, so two members of one group '
                 'can differ by more than the threshold if an intermediate isolate links them.</p>')
    return parts


def simple_report_html(records, *, selected_ids, investigation=None, settings=None, graph_png=None,
                       graph_mime='image/png', scope_note=None, scope_implicit=False) -> str:
    """One short page in plain words: the picture, the genes, and each closest match.

    ``investigation`` is an already-built comparison snapshot and ``graph_png``
    the picture the caller rasterised (raw bytes, a ``data:`` URI, or a file
    path, described by ``graph_mime``). Missing evidence is printed, never
    dropped to keep the page tidy.
    """
    options = {**SIMPLE_REPORT_PRESET, **(settings or {})}
    rows = _snapshot(records, selected_ids, investigation=investigation)
    names = {row.get('sample_id'): row.get('sample_name') or row.get('sample_id') for row in rows}
    snapshot = investigation if options['investigation'] else None
    absent = _escape(COMPARISON_OFF if investigation else NO_COMPARISON)
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    parts = [_HEAD_OPEN, _escape(options['title']), _HEAD_CLOSE,
             '<h1>' + _escape(options['title']) + '</h1>',
             f'<p>WMLSTudio {_escape(__version__)} · {_escape(stamp)} · {len(rows)} isolate(s)</p>',
             ('<p class="notice">' if scope_implicit else '<p>') +
             _escape(scope_note or f'This report covers the {len(rows)} isolate(s) chosen for it.') + '</p>',
             '<p class="notice"><b>Research and review evidence, not a clinical diagnosis.</b> Genomic similarity does '
             'not prove transmission, or its direction. ' + _escape(SUSCEPTIBILITY_CAVEAT) + '</p>']
    if snapshot is None:
        parts.append('<h2>How close are these isolates?</h2><p class="notice">' + absent + '</p>')
    else:
        parts.extend(_comparison_section(snapshot, rows, options, graph_png, graph_mime))
    parts.extend(_resistance_section(rows, options))
    if snapshot is None:
        parts.append('<h2>Closest matches</h2><p class="notice">' + absent + '</p>')
    else:
        parts.extend(_proximity_section(snapshot, rows, names))
    parts.append('<h2>What this report does not tell you</h2><ul>')
    parts.extend('<li>' + _escape(limitation) + '</li>' for limitation in LIMITATIONS)
    parts.append('</ul>')
    if snapshot is None:
        parts.append('<p class="muted">No comparison snapshot: this report carries no distances and no picture.</p>')
    else:
        parts.append('<p class="muted">Snapshot ' + _escape(snapshot.get('snapshot_id')) + ' · reference SHA-256 ' +
                     _escape(str(snapshot.get('scheme_digest') or '')[:12]) + ' · distance method ' +
                     _escape(snapshot.get('metric_version')) + '</p>')
    parts.append('</body></html>')
    return ''.join(parts)
