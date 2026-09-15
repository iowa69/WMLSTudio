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
from wmlstudio.export import (
    REPORT_PRESETS,
    _escape,
    _snapshot,
    investigation_document,
    threshold_provenance,
)

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

# Escape-stable like the caveat above, and printed joined to the plasmid table so
# no option can put the table on one page and its boundary on another.
PLASMID_BOUNDARY = (
    'A replicon marker sitting on an assembled contig is evidence about that contig. No plasmid was '
    'reconstructed, no plasmid was counted, no relaxase or mate-pair-formation type was assigned and no '
    'mobility was predicted: this is not MOB-suite and is not equivalent to it. Two isolates carrying the '
    'same replicon name are not thereby carrying the same plasmid.'
)

LIMITATIONS = (
    'It does not prove transmission, or its direction. Genomic similarity is one line of evidence for '
    'epidemiological review.',
    'Resistance genes are not measured susceptibility. Use laboratory testing for treatment decisions.',
    'Loci that could not be called are unknown, not identical — the shared / total column shows how much '
    'was actually compared.',
    'Isolates are grouped by single linkage, so members of one group can differ by more than the threshold.',
    'A sequence type (7 loci), a core-genome comparison (hundreds to thousands of targets) and a SNP '
    'distance are three different quantities. They never share a scale, a column or a threshold, a cutoff '
    'published for one is not a cutoff for the other, and every distance in this report is the one named '
    'above it — SNP distances are not reported here at all.',
    'A replicon marker is evidence about the assembled contig it was found on. It is not a plasmid, not a '
    'count of plasmids and not proof that a resistance gene beside it can transfer.',
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


def _typing_sentence(provenance) -> str:
    """Which typing this summary is based on, in one sentence a reader can repeat."""
    scheme, loci = provenance['scheme'], provenance['total_loci']
    if provenance['typing']['kind'] == 'cgmlst':
        return ('This summary is based on core-genome typing: every isolate was compared at the ' + str(loci) +
                ' targets of “' + scheme + '”. Those distances are not sequence-type distances.')
    if provenance['typing']['kind'] == 'mlst':
        return ('This summary is based on classical MLST: every isolate was compared at the ' + str(loci) +
                ' loci of “' + scheme + '”. That is a sequence-type comparison, not a core-genome one.')
    return ('The reference “' + scheme + '” records no target count, so this summary cannot say which typing '
            'it is based on.')


def _suggestion_lines(provenance) -> list[str]:
    """What has been published for this organism — offered, never applied."""
    if provenance['suggestion_blocked']:
        return ['<p class="muted">' + _escape(provenance['suggestion_blocked']) + '</p>']
    suggestion = provenance['suggestion']
    if not suggestion:
        return []
    # The headline already leads with "Suggested, not applied"; a second label
    # in front of it would only push those words further from the number.
    parts = ['<p class="muted">' + _escape(suggestion['headline']) + '</p>']
    entry = suggestion.get('suggestion')
    if entry is None:
        return parts
    parts.append('<p class="muted">' + _escape(entry['source']['citation']) + ' (' +
                 _escape(entry['source']['doi']) + '). <b>The authors’ own caveat:</b> “' +
                 _escape(entry['caveat']) + '”</p>')
    if suggestion['scheme_match']['checked']:
        parts.append('<p class="muted">' + _escape(suggestion['scheme_match']['reason']) + '</p>')
    if suggestion['disagreement']:
        parts.append('<p class="muted">' + _escape(suggestion['disagreement']) + '</p>')
    parts.append('<p class="muted"><b>None of this is applied to this report.</b> ' +
                 _escape(suggestion['notice']) + '</p>')
    return parts


def _threshold_lines(provenance) -> list[str]:
    """Say which of the two the report is showing: a published cutoff, or your own.

    Only a cutoff whose reference, number and minimum overlap still match the
    reviewed context reads as published. A citation whose binding has changed is
    printed as history beside a local setting, never as approval of it.
    """
    citation = _escape(provenance['citation'] or 'No citation recorded')
    doi = _escape(provenance['doi'] or 'no DOI recorded')
    own = 'At most ' + _differences(provenance['threshold']) + ' over ' + _escape(provenance['total_loci']) + \
          ' targets is a local setting for this comparison, not a validated rule.'
    if provenance['source'] == 'adopted_publication':
        parts = ['<p><b>Where this threshold comes from:</b> it is a published cutoff that was reviewed and '
                 'adopted for this comparison.</p>',
                 '<p class="muted">' + citation + ' (' + doi + ') — published cutoff at most ' +
                 _escape(provenance['published_threshold']) + ' ' + _escape(provenance['published_unit']) +
                 (' on ' + _escape(provenance['published_scheme']) if provenance['published_scheme'] else '') +
                 (' over ' + _escape(provenance['published_locus_count']) + ' targets'
                  if provenance['published_locus_count'] else '') +
                 '; applied here as at most ' + _differences(provenance['threshold']) +
                 '. Matching an organism or a locus count is not clinical validation.</p>']
        if provenance['caveat']:
            parts.append('<p class="notice"><b>The authors’ own caveat:</b> “' +
                         _escape(provenance['caveat']) + '”</p>')
        return parts
    if provenance['status'] == 'citation_only':
        stated = ('it is your own setting. A publication is attached for context only — ' + citation +
                  ' (' + doi + ') — and no number was adopted from it. ' + own)
    elif provenance['status'] == 'no_guidance':
        stated = 'no published cutoff is attached to it. ' + own
    else:
        stated = (_escape(provenance['reason']) + ' Citation retained: ' + citation + ' (' + doi + '). ' + own)
    # A catalog entry is context that exists, never a rule that was applied here.
    return ['<p><b>Where this threshold comes from:</b> ' + stated + '</p>'] + _suggestion_lines(provenance)


def _comparison_section(snapshot, options, graph_png, graph_mime, provenance) -> list[str]:
    parts = ['<h2>How close are these isolates?</h2>',
             '<p>This shows how many allele differences separate the isolates that were compared. '
             'It does not show who infected whom, in which direction, or when.</p>',
             '<p>' + _escape(_typing_sentence(provenance)) + '</p>']
    scheme, loci = provenance['scheme'], provenance['total_loci']
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
    parts.append('<p><b>Reference used:</b> ' + _escape(scheme) + ' · ' + _escape(loci) + ' loci per profile · ' +
                 _escape(provenance['typing']['label']) + '</p>')
    overlap = snapshot.get('min_overlap')
    share = f'{overlap:.0%}' if isinstance(overlap, (int, float)) and not isinstance(overlap, bool) else overlap
    parts.append('<p><b>Close</b> in this report means at most ' + _differences(snapshot.get('threshold')) +
                 ' across the ' + _escape(loci) + ' loci of this reference, counted only at loci both isolates '
                 'could call, with at least ' + _escape(share) + ' of the loci shared.</p>')
    parts.extend(_threshold_lines(provenance))
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


def _plasmid_cells(row) -> list[str]:
    """One isolate's replicon markers and same-contig co-locations, each with its own gate.

    The markers come from the AMR/plasmid assay and the co-locations from the
    characterization, so the two are gated separately: evidence belonging to an
    earlier assembly is withheld and named, never printed as this one's.
    """
    from wmlstudio.characterization import current_characterization

    status = row.get('hydra_evidence_status') or 'missing'
    replicons = [str(name) for name in (row.get('plasmid_replicons') or [])]
    if status == 'stale':
        markers = _escape('Not shown — the saved plasmid result belongs to a different sequence file')
    elif status == 'missing':
        markers = _escape('Not assessed')
    elif replicons:
        markers = _escape(', '.join(replicons))
        if status == 'unverified':
            markers += _muted('The source report’s identity was not confirmed.')
    else:
        markers = _escape('No replicon marker reported by the reference database used')
    state = current_characterization(row)
    evidence = (state.get('evidence') or {}).get('plasmid_hypotheses') or {}
    links = evidence.get('contig_associations') or []
    if state['status'] != 'current':
        shared = _escape('Not assessed') + _muted(state['reason'])
    elif links:
        shared = _escape('; '.join(sorted({f"{link['marker']} with {link['replicon']} on "
                                           f"{link['contig']}" for link in links})))
        shared += _muted('Same assembled contig only — a hypothesis, not a plasmid-borne gene.')
    else:
        shared = _escape('No resistance or virulence gene shared a contig with a replicon marker')
        shared += _muted('A plasmid contig can assemble without its replicon, so this does not place '
                         'those genes on the chromosome.')
    return [_escape(row.get('sample_name') or row.get('sample_id')), markers, shared,
            _escape(_STATE_WORDS.get(status, status))]


def _plasmid_section(rows) -> list[str]:
    parts = ['<h2>Plasmid markers</h2>',
             '<p>This lists the plasmid replicon markers found in each genome, and any resistance or '
             'virulence gene that sat on the same assembled contig as one of them. It does not show '
             'which plasmids an isolate carries, or that any gene can move between isolates.</p>',
             '<table border="1" cellpadding="5" cellspacing="0"><tr><th>Isolate</th>'
             '<th>Replicon markers found</th><th>Genes on the same contig as a replicon</th>'
             '<th>Evidence state</th></tr>']
    for row in rows:
        parts.append(_row(_plasmid_cells(row)))
    # The boundary closes the table in one string, so no option can separate them.
    parts.append('</table><p class="notice"><b>' + _escape(PLASMID_BOUNDARY) + '</b></p>')
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


def _proximity_section(snapshot, rows, names, provenance) -> list[str]:
    doc = investigation_document(snapshot, [row.get('sample_id') for row in rows])
    threshold = snapshot.get('threshold')
    parts = ['<h2>Closest matches</h2>',
             '<p>This shows each isolate’s closest comparable isolate and how much of the reference the two shared. '
             'It does not show a transmission link, and a locus that could not be called is never counted as a match.</p>',
             '<p class="muted">Every distance below is measured on “' + _escape(provenance['scheme']) + '”, ' +
             _escape(provenance['typing']['label']) + '. ' + _escape(provenance['typing']['note']) + '</p>',
             '<table border="1" cellpadding="5" cellspacing="0"><tr><th>Isolate</th>'
             '<th>Allele differences to closest (of ' + _escape(provenance['total_loci']) + ' targets)</th>'
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
    # The organisms come from the report rows, not from the snapshot profiles: an
    # imported profile may carry no organism at all, and a suggestion must never
    # be attached to an organism this report cannot name.
    provenance = (threshold_provenance(snapshot, [row.get('organism') for row in rows])
                  if snapshot is not None else None)
    if snapshot is None:
        parts.append('<h2>How close are these isolates?</h2><p class="notice">' + absent + '</p>')
    else:
        parts.extend(_comparison_section(snapshot, options, graph_png, graph_mime, provenance))
    parts.extend(_resistance_section(rows, options))
    # Off by default in this preset and printed only when it is asked for: a
    # section nobody selected is covered by the limitation saying that what is
    # absent from this report is not a negative result.
    if options.get('plasmid_hypotheses'):
        parts.extend(_plasmid_section(rows))
    if snapshot is None:
        parts.append('<h2>Closest matches</h2><p class="notice">' + absent + '</p>')
    else:
        parts.extend(_proximity_section(snapshot, rows, names, provenance))
    parts.append('<h2>What this report does not tell you</h2><ul>')
    parts.extend('<li>' + _escape(limitation) + '</li>' for limitation in LIMITATIONS)
    parts.append('</ul>')
    if snapshot is None:
        parts.append('<p class="muted">No comparison snapshot: this report carries no distances and no picture.</p>')
    else:
        parts.append('<p class="muted">Snapshot ' + _escape(snapshot.get('snapshot_id')) + ' · reference SHA-256 ' +
                     _escape(str(snapshot.get('scheme_digest') or '')[:12]) + ' · distance method ' +
                     _escape(snapshot.get('metric_version')) + ' · ' + _escape(provenance['typing']['label']) +
                     ' · threshold ' + _escape('adopted from a publication'
                                               if provenance['source'] == 'adopted_publication'
                                               else 'set locally for this comparison') + '</p>')
    parts.append('</body></html>')
    return ''.join(parts)
