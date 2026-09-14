"""Atomic, offline exports from an explicit snapshot of project records.

The caller supplies the export scope. Passing ``project.samples()`` exports all
records, including queued/failed samples, regardless of any view filter.
"""

from __future__ import annotations

import csv
import html
import json
import math
import os
import re
import tempfile
from collections.abc import Iterable, Mapping
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

from wmlstudio import __version__
from wmlstudio.sample_workflow import feature_fields, select_records

_FIELDS = (
    "sample_id", "sample_name", "input_path", "kind", "job_status", "status", "scheme",
    "scheme_digest", "st", "input_sha256", "error", "notes", "metadata", "qc",
    "genus", "species", "organism", "organism_source", "typing_mode", "amr_genes",
    "amr_classes", "virulence_genes", "plasmid_replicons", "stress_genes",
    "hydra_source_sample", "hydra_report_sha256", "cluster_label", "cluster_highlight",
    "additional_profiles",
    "hydra_evidence_status", "hydra_evidence_reason",
    "investigation_cluster", "investigation_cluster_id", "investigation_cluster_status",
    "nearest_allele_distance", "nearest_isolates", "within_threshold_isolates",
    "comparison_scheme_digest", "comparison_threshold", "comparison_min_overlap",
    "organism_typing",
)


def ensure_separate_destination(destination: str | Path, protected_paths: Iterable[str | Path]) -> None:
    """Never replace a sequence input or live project with an exported result."""
    target = Path(destination).expanduser().resolve()
    for value in protected_paths:
        if value and target == Path(value).expanduser().resolve():
            raise ValueError("Choose a different output file. This path belongs to an input or open project.")


def _snapshot(
    records: Iterable[Mapping[str, Any]], selected_ids=None, highlight_clusters=None,
    investigation=None,
) -> list[dict[str, Any]]:
    proximity = {}
    if investigation:
        from wmlstudio.investigation import proximity_rows
        proximity = {p['sample_id']: p for p in proximity_rows(investigation)}
    rows = []
    for record in select_records(records, selected_ids):
        if "result" in record:
            result = record.get("result")
            row = dict(result) if isinstance(result, Mapping) else {}
            row.update({
                "sample_id": record.get("id", row.get("sample_id", "")),
                "sample_name": record.get("name", row.get("sample_name", "")),
                "input_path": record.get("input_path", row.get("input_path", "")),
                "job_status": record.get("status", ""),
                "error": record.get("error", ""),
                "metadata": record.get("metadata", {}),
                "created_at": record.get("created_at", ""),
                "updated_at": record.get("updated_at", ""),
                "missing_input": bool(record.get("missing_input", False)),
                "additional_profiles": record.get("analyses", []),
            })
        else:
            row = dict(record)
        if "result" in record or "sample_name" in record:
            sample_id = str(row.get("sample_id", ""))
            row.update(feature_fields(record, (highlight_clusters or {}).get(sample_id)))
            if sample_id in proximity:
                close = proximity[sample_id]
                row.update(investigation_cluster=close['group_name'], investigation_cluster_id=close['group_id'],
                           investigation_cluster_status=close['group_status'], nearest_allele_distance=close['nearest_distance'],
                           nearest_isolates=[p['neighbour_name'] for p in close['nearest']],
                           within_threshold_isolates=[p['neighbour_name'] for p in close['within_threshold']],
                           comparison_scheme_digest=investigation['scheme_digest'], comparison_threshold=investigation['threshold'],
                           comparison_min_overlap=investigation['min_overlap'])
        rows.append(row)
    # Validate and detach nested data before opening any destination.
    return json.loads(json.dumps(rows, ensure_ascii=False, allow_nan=False))


def _protected_paths(rows):
    for row in rows:
        yield row.get("input_path", "")
        metadata = row.get("metadata", {})
        if isinstance(metadata, dict):
            workflow = metadata.get("workflow", {})
            if isinstance(workflow, dict):
                yield workflow.get("source_path", "")
                yield workflow.get('managed_path', '')
            reads = metadata.get('reads') or {}
            if isinstance(reads, dict):
                for entry in reads.get('reads', []):
                    if isinstance(entry, dict):
                        yield entry.get('path', '')
            assembly = metadata.get('assembly') or {}
            if isinstance(assembly, dict):
                for entry in (assembly.get('provenance') or {}).get('inputs', []):
                    if isinstance(entry, dict):
                        yield entry.get('path', '')
            evidence = metadata.get('hydra') or {}
            if isinstance(evidence, dict):
                yield (evidence.get('provenance') or {}).get('source_path', '')


def protected_input_paths(records):
    """Original reads/assemblies and linked source evidence, not only active inputs."""
    return _protected_paths(records)


@contextmanager
def _atomic_text(path: str | Path, *, newline: str | None = None) -> Iterable[TextIO]:
    destination = Path(path).expanduser().resolve()
    temporary: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8-sig" if newline == "" else "utf-8",
            newline=newline, prefix=f".{destination.name}.", suffix=".tmp",
            dir=destination.parent, delete=False,
        ) as handle:
            temporary = handle.name
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)


def _display(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True)
    return str(value)


def _csv_cell(value: Any) -> str | int | float:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value
    text = _display(value)
    # Spreadsheet programs may ignore leading whitespace before interpreting a
    # formula; protect names, paths, scheme IDs, and allele IDs alike.
    if text and (text[0] in "\t\r\n" or text.lstrip().startswith(("=", "+", "-", "@"))):
        return "'" + text
    return text


def _write_delimited(
    records: Iterable[Mapping[str, Any]], path: str | Path, delimiter: str,
    selected_ids=None, highlight_clusters=None, investigation=None,
) -> Path:
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters, investigation)
    loci = sorted({str(locus) for row in rows for locus in row.get("alleles", {})})
    columns = [*_FIELDS, *(f"allele:{locus}" for locus in loci)]
    with _atomic_text(path, newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(columns)
        for row in rows:
            writer.writerow([
                *(_csv_cell(row.get(field)) for field in _FIELDS),
                *(_csv_cell(row.get("alleles", {}).get(locus)) for locus in loci),
            ])
    return Path(path).expanduser().resolve()


def write_csv(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None, investigation=None) -> Path:
    """Write Excel-friendly UTF-8 CSV with guarded textual formula cells."""
    return _write_delimited(records, path, ",", selected_ids, highlight_clusters, investigation)


def write_tsv(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None, investigation=None) -> Path:
    return _write_delimited(records, path, "\t", selected_ids, highlight_clusters, investigation)


def write_json(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None, investigation=None) -> Path:
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters, investigation)
    document = {
        "format_version": 1,
        "application": f"WMLSTudio {__version__}",
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "samples": rows,
        "report_scope": {"mode": "selected" if selected_ids is not None else "provided",
                         "sample_ids": [row.get("sample_id") for row in rows]},
    }
    if investigation:
        document['investigation'] = investigation_document(investigation, [r['sample_id'] for r in rows])
    with _atomic_text(path) as handle:
        json.dump(document, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write("\n")
    return Path(path).expanduser().resolve()


def write_distances(
    rows: Iterable[Mapping[str, Any]], path: str | Path, min_overlap: float = 0.95,
) -> Path:
    """Export distance evidence with its comparison policy and explicit pair scope."""
    if isinstance(min_overlap, bool) or not math.isfinite(min_overlap) or not 0 <= min_overlap <= 1:
        raise ValueError("Minimum overlap must be a finite number between 0 and 1.")
    document = {
        "format_version": 1,
        "application": f"WMLSTudio {__version__}",
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "min_overlap": min_overlap,
        "metric": "allele differences",
        "missing_policy": (
            "Count differences only among shared, unambiguous, known alleles. "
            "Overlap is shared known loci divided by the union of profile loci. "
            "Mixed samples and pairs without shared calls are excluded; "
            "non-comparable distances are null."
        ),
        "scheme_requirement": "Identical nonempty scheme names and scheme fingerprints.",
        "interpretation": "Allele similarity alone does not establish an outbreak or transmission.",
        "pairs": _snapshot(rows),
    }
    with _atomic_text(path) as handle:
        json.dump(document, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write("\n")
    return Path(path).expanduser().resolve()


def _escape(value: Any) -> str:
    return html.escape(_display(value), quote=True)


def investigation_document(snapshot, selected_ids=None):
    """Explicit focal isolates plus labelled cohort context; no scope ambiguity."""
    from wmlstudio.investigation import proximity_rows
    profiles = {p['sample_id']: p for p in snapshot.get('profiles', [])}
    selected = set(profiles) if selected_ids is None else set(selected_ids)
    focus = selected & profiles.keys()
    metadata = {key: snapshot.get(key) for key in ('snapshot_id', 'created_at', 'investigation_id',
        'investigation_name', 'protocol', 'scheme', 'scheme_digest', 'threshold', 'min_overlap',
        'missing_policy', 'metric_version', 'interpretation', 'preview', 'saved_threshold', 'saved_min_overlap', 'reuse', 'threshold_evidence')}
    metadata.update(focal_sample_ids=sorted(focus), outside_snapshot_ids=sorted(selected - profiles.keys()),
                    cohort_size=len(profiles), focal_profiles=[profiles[sid] for sid in sorted(focus)],
                    groups=[g for g in snapshot.get('groups', []) if focus.intersection(g['members'])],
                    review_groups=[g for g in snapshot.get('review_groups', []) if focus.intersection(g['sample_ids'])],
                    proximity=proximity_rows(snapshot, focus),
                    scope_note='Only focal isolates are sample records. Named neighbours and group members are explicitly labelled comparison-cohort context.')
    return json.loads(json.dumps(metadata, ensure_ascii=False, allow_nan=False))


def investigation_html(snapshot, selected_ids=None, sections=None):
    """Offline, escaped investigation and isolate-proximity section, also printable."""
    doc = investigation_document(snapshot, selected_ids)
    parts = ['<section><h2>Investigation · ' + _escape(doc.get('investigation_name') or 'Unsaved cohort') + '</h2>',
             '<p>' + _escape(doc['scope_note']) + '</p>',
             '<p><b>Reference:</b> ' + _escape(doc['scheme']) + '<br><b>Scheme SHA-256:</b> ' + _escape(doc['scheme_digest']) +
             '<br><b>Snapshot:</b> ' + _escape(doc['snapshot_id']) + ' · ' + _escape(doc['created_at']) + '</p>',
             f"<p><b>Single-link threshold:</b> ≤ {_escape(doc['threshold'])} allele differences. "
             f"<b>Minimum shared / total loci:</b> {_escape(doc['min_overlap'])}. "
             f"<b>Cohort:</b> {doc['cohort_size']} isolates; {len(doc['focal_sample_ids'])} focal isolates.</p>",
             '<p><b>Local protocol:</b> ' + _escape(doc.get('protocol') or 'Exploratory; not clinically calibrated') + '</p>',
             '<p class="notice">' + _escape(doc['interpretation']) + '</p>']
    if doc.get('preview'):
        parts.append('<p class="notice"><b>Unsaved threshold preview.</b> The saved protocol threshold was ' +
                     _escape(doc.get('saved_threshold')) + '; this report explicitly uses the threshold shown above.</p>')
    evidence = doc.get('threshold_evidence') or {}
    if evidence:
        from wmlstudio.investigation import threshold_guidance_status
        binding = threshold_guidance_status(snapshot)
        parts.append('<h3>Published threshold guidance and local application</h3>')
        parts.append('<p class="notice"><b>' + _escape(binding['status'].replace('_', ' ')) + ':</b> ' + _escape(binding['reason']) + '</p>')
        for key, value in evidence.items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                parts.append('<p><b>' + _escape(key.replace('_', ' ')) + ':</b> ' + _escape(value) + '</p>')
            elif key in {'citation', 'citations', 'scheme_scope', 'context'}:
                parts.append('<p><b>' + _escape(key.replace('_', ' ')) + ':</b> ' + _escape(value) + '</p>')
        parts.append('<p>A published cutoff is scoped to its study, scheme version and method; a local variant requires an explicit justification and is not automatic clinical validation.</p>')
    if doc['outside_snapshot_ids']:
        parts.append(f"<p>{len(doc['outside_snapshot_ids'])} selected sample(s) have no profile in this snapshot; no distance was assigned.</p>")
    parts.append('<h3>Threshold groups and membership changes</h3><table border="1" cellpadding="5" cellspacing="0"><tr><th>Group / stable ID</th><th>State</th><th>Members</th><th>Change</th><th>Within-group evidence</th></tr>')
    for group in doc['groups']:
        warning = ('Chaining: direct distance can exceed link threshold. ' if group.get('chained') else '')
        values = [group['name'] + ' · ' + group['id'], group['status'], len(group['members']),
                  group['change'], warning + f"Maximum direct distance: {group['max_direct_distance']}; "
                  f"unassessed pairs: {group['unassessed_within_pairs']}. " + group.get('reason', '')]
        parts.append('<tr>' + ''.join('<td>' + _escape(value) + '</td>' for value in values) + '</tr>')
    parts.append('</table>')
    for group in doc['review_groups']:
        parts.append('<p><b>Frozen review group:</b> ' + _escape(group['name']) + ' · ' + _escape(group['frozen_at']) +
                     ' · ' + str(len(group['sample_ids'])) + ' fixed isolate IDs. Automatic clustering does not rewrite these members.</p>')
    profile_map = {p['sample_id']: p for p in doc['focal_profiles']}
    for focal in doc['proximity']:
        profile = profile_map[focal['sample_id']]
        parts.append('<h3>Isolate proximity · ' + _escape(focal['sample_name']) + '</h3>')
        parts.append('<p><b>Primary MLST/ST:</b> ' + _escape(profile.get('primary_st') or profile.get('st') or 'Unassigned') +
                     ' · <b>Callable loci:</b> ' + str(profile['callable_loci']) + '/' + str(profile['total_loci']) +
                     (' · <b>AMR context:</b> ' + _escape(profile.get('amr_genes') or 'Not recorded') + ' · ' + _escape(profile.get('amr_evidence_status') or 'Identity unverified') if sections is None or sections.get('amr') else '') + '</p>')
        parts.append('<p>Nearest accepted distance: ' + _escape(focal['nearest_distance'] if focal['nearest_distance'] is not None else 'Not comparable') +
                     f". {focal['compared_count']} neighbours compared; {focal['excluded_count']} comparisons excluded.</p>")
        rows = focal['within_threshold'] or focal['nearest']
        if rows:
            parts.append('<table border="1" cellpadding="5" cellspacing="0"><tr><th>Comparison-cohort neighbour</th><th>Allele differences</th><th>Shared / total loci</th><th>Overlap</th></tr>')
            for pair in rows:
                parts.append('<tr>' + ''.join('<td>' + _escape(value) + '</td>' for value in (
                    pair['neighbour_name'], pair['distance'], f"{pair['shared_loci']}/{pair['total_loci']}",
                    f"{pair['overlap']:.3f}")) + '</tr>')
            parts.append('</table>')
        if focal['excluded']:
            reasons = sorted({row['reason'] for row in focal['excluded']})
            parts.append('<p><b>Exclusions:</b> ' + _escape('; '.join(reasons)) + '</p>')
    parts.append('<h3>Core and accessory evidence are separate</h3><p>The distances above use only the pinned allele-profile scheme. '
                 'AMR, virulence, plasmid markers and custom annotations are context, not extra distance loci. '
                 'Identical AMR or replicon names do not establish identical plasmids or transfer; genomic drug associations are not susceptibility testing.</p></section>')
    return ''.join(parts)


def characterization_html(record, sections=None):
    """Compact characterization evidence with an input-identity gate, not raw arrays."""
    from wmlstudio.characterization import current_characterization
    state = current_characterization(record)
    if state['status'] == 'missing':
        return ''
    evidence = state['evidence']
    if not evidence:
        return '<h3>Characterization evidence</h3><p>Archived or unverified input identity: concrete characterization calls are not presented as current.</p>'
    parts = ['<h3>Independent species and accessory evidence</h3>']
    species = evidence.get('species_evidence') or {}
    parts.append('<p><b>Species evidence:</b> ' + _escape(species.get('status', 'not_run')) + ' · ' +
                 _escape(' '.join(str(species.get(k) or '') for k in ('genus', 'species', 'subspecies')).strip()) +
                 ' · ' + _escape(species.get('method') or '') + '</p>')
    for key, title, rows_key, name_key in [('virulence', 'Virulence assay', 'hits', 'gene'),
                                          ('plasmid_hypotheses', 'Plasmid-marker assay', 'replicons', 'gene'),
                                          ('drug_associations', 'Reference-reported drug associations', 'associations', 'class')]:
        if sections is not None and not sections.get(key, True):
            continue
        module = evidence.get(key) or {}
        names = sorted({str(row.get(name_key)) for row in module.get(rows_key, []) if isinstance(row, dict) and row.get(name_key)}) if module.get('status') in {'completed', 'detected', 'not_detected'} else []
        parts.append('<p><b>' + title + ':</b> ' + _escape(module.get('status', 'not_run')) + ' · ' +
                     _escape('; '.join(names[:30]) or module.get('reason') or 'No determinants reported for this assay.') + '</p>')
        if len(names) > 30:
            parts.append(f'<p>{len(names)} annotations total; complete evidence retained in JSON.</p>')
        for limitation in module.get('limitations', [])[:3]:
            parts.append('<p class="muted">' + _escape(limitation) + '</p>')
    return ''.join(parts)


REPORT_PRESETS = {
    'isolate': {'title': 'Isolate evidence review', 'investigation': False, 'qc': True, 'amr': True,
                'virulence': True, 'plasmid_hypotheses': True, 'drug_associations': True, 'graph': False,
                'provenance': False, 'graph_jpeg': False},
    'cohort': {'title': 'Selected cohort review', 'investigation': True, 'qc': True, 'amr': True,
               'virulence': True, 'plasmid_hypotheses': True, 'drug_associations': True, 'graph': True,
               'provenance': False, 'graph_jpeg': False},
    'ipc': {'title': 'IPC cluster review', 'investigation': True, 'qc': True, 'amr': True,
            'virulence': True, 'plasmid_hypotheses': True, 'drug_associations': True, 'graph': True,
            'provenance': False, 'graph_jpeg': False},
    'proximity': {'title': 'Isolate proximity review', 'investigation': True, 'qc': True, 'amr': True,
                  'virulence': False, 'plasmid_hypotheses': False, 'drug_associations': False, 'graph': True,
                  'provenance': False, 'graph_jpeg': False},
    # A different document shape, not a different set of facts: one short page in
    # plain words for a reader who is not a bioinformatician. Kept last so no
    # preset that a saved template or a combo index already names can move.
    'simple': {'title': 'Simple outbreak summary', 'investigation': True, 'qc': False, 'amr': True,
               'virulence': False, 'plasmid_hypotheses': False, 'drug_associations': False, 'graph': True,
               'provenance': False, 'layout': 'one_page', 'graph_jpeg': True},
}


def review_report_html(records, *, selected_ids, investigation=None, options=None, graph_png=None,
                       graph_mime='image/png', scope_note=None, scope_implicit=False):
    """Readable research/IPC report with explicit scope and opt-in raw appendix.

    ``graph_png`` is raw image bytes whose type ``graph_mime`` states, so a
    caller may embed JPEG or PNG without changing this signature. ``scope_note``
    replaces the default header scope sentence, and ``scope_implicit`` marks a
    scope the user did not choose so the document says so itself.
    """
    import base64
    settings = {**REPORT_PRESETS['cohort'], **(options or {})}
    if settings.get('layout') == 'one_page':
        # Lazy on purpose: simple_report imports this module at module level.
        from wmlstudio.simple_report import simple_report_html
        return simple_report_html(records, selected_ids=selected_ids, investigation=investigation,
                                  settings=settings, graph_png=graph_png, graph_mime=graph_mime,
                                  scope_note=scope_note, scope_implicit=scope_implicit)
    rows = _snapshot(records, selected_ids, investigation=investigation)
    parts = ['<!doctype html><html><head><meta charset="utf-8">',
             '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; style-src \'unsafe-inline\'; img-src data:">',
             '<title>' + _escape(settings['title']) + '</title>',
             '<style>body{font-family:Segoe UI,Arial,sans-serif;color:#203c45;max-width:1080px;margin:auto;padding:28px;line-height:1.45}'
             'h1,h2,h3{color:#156e65}h2{border-bottom:2px solid #d6e7e3;padding-bottom:8px}'
             'table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #d6e0e4;padding:7px;text-align:left;vertical-align:top}'
             'th{background:#edf5f3}.notice{padding:12px;background:#fff4dd;border-left:4px solid #c89031}'
             '.muted{color:#5c6f76}section{margin:24px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:10px}'
             'img{max-width:100%;height:auto}@media print{body{padding:0}tr{break-inside:avoid}h2,h3{break-after:avoid}}</style></head><body>',
             '<h1>' + _escape(settings['title']) + '</h1>',
             f'<p>WMLSTudio {_escape(__version__)} · {_escape(datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"))} · '
             + _escape(scope_note or f'{len(rows)} explicitly selected isolate(s)') + '</p>',
             '<p class="notice"><b>Research / review evidence, not a clinical diagnosis.</b> Genomic similarity does not prove transmission. '
             'Resistance determinants are not measured susceptibility; use validated laboratory testing and epidemiological review.</p>']
    if scope_implicit:
        # A scope the user did not choose is stated in the document itself, not
        # only in a message that disappears once the file has been saved.
        parts.append('<p class="notice"><b>This scope was not chosen for this report.</b> ' +
                     _escape(scope_note or '') + '</p>')
    parts.append('<h2>At a glance</h2><table cellpadding="6" cellspacing="0" border="1"><tr><th>Isolate</th><th>Organism</th><th>Primary ST</th><th>Typing / QC state</th>' +
                 ('<th>AMR evidence</th>' if settings['amr'] else '') + '</tr>')
    for row in rows:
        values = [row.get('sample_name'), row.get('organism') or 'Unassigned', row.get('st') or 'Unassigned', row.get('status') or row.get('job_status') or 'Not analysed']
        if settings['amr']:
            genes = row.get('amr_genes') or []
            values.append('; '.join(genes) if genes else 'No determinants reported' if row.get('hydra_evidence_status') == 'current' else row.get('hydra_evidence_status') or 'Not assessed')
        parts.append('<tr>' + ''.join('<td>' + _escape(value) + '</td>' for value in values) + '</tr>')
    parts.append('</table>')
    if settings['graph']:
        if graph_png:
            encoded = base64.b64encode(graph_png).decode('ascii')
            parts.append('<h2>Comparison-cohort context</h2><p>Selected focal isolates are highlighted. Other nodes are context only, not additional reported sample records. Layout is not a transmission tree.</p>'
                         '<img width="920" src="data:' + _escape(graph_mime or 'image/png') + ';base64,' + encoded + '" alt="Allele-distance graph with focal isolates highlighted">')
        else:
            parts.append('<p class="muted">Graph omitted: no matching comparison snapshot is currently available.</p>')
    if settings['investigation']:
        parts.append(investigation_html(investigation, [r['sample_id'] for r in rows], settings) if investigation else
                     '<p class="notice">No comparison snapshot selected. No proximity or cluster inference was made.</p>')
    for row in rows:
        parts.append('<section><h2>' + _escape(row.get('sample_name')) + '</h2>')
        if row.get('cluster_label'):
            parts.append('<p class="notice"><b>User-defined review highlight:</b> ' + _escape(row['cluster_label']) + ' (not an inferred transmission assignment)</p>')
        parts.append('<p><b>Primary typing:</b> ' + _escape(row.get('scheme') or 'Not assessed') + ' · ST ' + _escape(row.get('st') or 'Unassigned') + '</p>')
        if row.get('error') or row.get('missing_input'):
            parts.append('<p class="notice">' + _escape(row.get('error') or 'Original sequence input unavailable; stored profile evidence only.') + '</p>')
        if settings['qc']:
            qc = row.get('qc') or {}
            summary = {key: value for key, value in qc.items() if isinstance(value, (int, float, str, bool)) or value is None}
            parts.append('<h3>Input quality evidence</h3><p>' + ('; '.join(_escape(k.replace('_', ' ')) + ': ' + _escape(v) for k, v in summary.items()) or 'No sequence QC recorded (profile-only imports do not imply QC passed).') + '</p>')
        if settings['amr']:
            parts.append('<h3>AMR determinants</h3><p><b>Evidence state:</b> ' + _escape(row.get('hydra_evidence_status') or 'Not assessed') + '. ' + _escape(row.get('hydra_evidence_reason') or '') + '</p><p>' +
                         _escape('; '.join(row.get('amr_genes') or []) or 'No current determinants to report. Not assessed is not a negative result.') + '</p>')
        parts.append(characterization_html(row, settings))
        profiles = [r for r in row.get('additional_profiles', []) if r.get('scheme_digest') != row.get('scheme_digest')]
        if profiles:
            parts.append('<h3>Additional allele profiles</h3><ul>')
            for profile in profiles:
                alleles = profile.get('alleles') or {}
                parts.append('<li>' + _escape(profile.get('scheme')) + f" · {sum(v is not None for v in alleles.values())}/{len(alleles)} called loci · " + _escape(profile.get('status')) + '</li>')
            parts.append('</ul>')
        custom = {key: value for key, value in (row.get('metadata') or {}).items()
                  if key not in {'workflow', 'assembly', 'hydra', 'characterization', 'input_identity'} and not isinstance(value, (dict, list))}
        custom.update({key: value for key, value in ((row.get('metadata') or {}).get('annotations') or {}).items()
                       if isinstance(value, (str, int, float, bool))})
        if custom:
            parts.append('<h3>Isolate annotations</h3><p>' + '<br>'.join(_escape(k) + ': ' + _escape(v) for k, v in custom.items()) + '</p>')
        if settings['provenance']:
            parts.append('<h3>Provenance appendix</h3><p>Full stored evidence below is archival; current/stale assessment above controls interpretation.</p><pre>' + _escape(row) + '</pre>')
        parts.append('</section>')
    parts.append('<footer>Core allele differences, AMR, virulence and plasmid hypotheses are separate evidence domains. Report sections are customizable; omitted assays must not be interpreted as absent.</footer></body></html>')
    return ''.join(parts)


def write_review_report(records, path, *, selected_ids, investigation=None, options=None, graph_png=None,
                        graph_mime='image/png', scope_note=None, scope_implicit=False):
    records = list(records)
    ensure_separate_destination(path, _protected_paths(records))
    markup = review_report_html(records, selected_ids=selected_ids, investigation=investigation, options=options,
                                graph_png=graph_png, graph_mime=graph_mime, scope_note=scope_note,
                                scope_implicit=scope_implicit)
    with _atomic_text(path) as handle:
        handle.write(markup)
    return Path(path).expanduser().resolve()


def write_html(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None, investigation=None) -> Path:
    """Write an escaped, self-contained report that also prints without scripts."""
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters, investigation)
    exported = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    output = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width,initial-scale=1">',
        '<meta http-equiv="Content-Security-Policy" '
        'content="default-src \'none\'; style-src \'unsafe-inline\'">',
        '<title>WMLSTudio · Sequence typing report</title>',
        '<style>:root{color-scheme:light;font-family:Segoe UI,Arial,sans-serif;color:#173042;'
        'background:#eef4f6}body{max-width:1160px;margin:auto;padding:40px}header{padding:28px;'
        'background:#143f4b;color:white;border-radius:18px}h1{margin:0 0 8px;font-size:30px}'
        'h2{font-size:20px;margin-top:0}.muted{color:#547080}.notice{padding:16px;'
        'border-left:4px solid #c58b2a;background:#fff7e6;margin:22px 0}'
        'section{background:white;padding:24px;border-radius:14px;margin:24px 0;'
        'overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px}'
        'th,td{text-align:left;padding:10px;border-bottom:1px solid #dce7ea;vertical-align:top;'
        'overflow-wrap:anywhere}th{background:#f1f7f8}dl{display:grid;'
        'grid-template-columns:130px minmax(0,1fr);gap:8px}dt{font-weight:600}dd{margin:0}'
        'pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f4f8f9;padding:14px;'
        'border-radius:8px}.scroll{overflow:auto}footer{font-size:12px;color:#547080}'
        '@media print{body{max-width:none;padding:0;background:white;font-size:11px}'
        'header{background:white;color:#173042;border-bottom:2px solid #143f4b;border-radius:0;'
        'padding:12px 0}section{padding:10px 0;margin:12px 0}.scroll{overflow:visible}'
        'tr{break-inside:avoid}thead{display:table-header-group}th,td{padding:5px}'
        'h2,h3{break-after:avoid}a{color:inherit}}'
        '</style></head><body><header><h1>WMLSTudio</h1>'
        '<div>Sequence typing report</div></header>',
        f'<p class="muted">{len(rows)} sample record(s) · Exported {_escape(exported)} · '
        f'Application {_escape(__version__)}</p>',
        f'<p>Report scope: {"selected samples" if selected_ids is not None else "provided project records"}. '
        'Highlighted groups are user-defined and do not establish transmission.</p>',
        '<p class="notice">Results describe the supplied sequence data and scheme snapshot. '
        'Missing, mixed, or ambiguous calls must not be treated as allele matches. '
        'Allele similarity alone does not establish an outbreak or transmission.</p>',
        investigation_html(investigation, [row['sample_id'] for row in rows]) if investigation else '',
        '<section><h2>Sample overview</h2><div class="scroll"><table><thead><tr>'
        '<th>Sample</th><th>Job</th><th>Typing result</th><th>Scheme</th><th>ST</th>'
        '<th>Known alleles</th><th>Organism</th><th>AMR genes</th><th>Group</th></tr></thead><tbody>',
    ]
    for row in rows:
        alleles = row.get("alleles", {})
        called = sum(value is not None for value in alleles.values())
        cells = (
            row.get("sample_name", ""), row.get("job_status", ""), row.get("status", ""),
            row.get("scheme", ""), row.get("st", ""), f"{called}/{len(alleles)}",
            row.get("organism", ""), ", ".join(row.get("amr_genes", [])), row.get("cluster_label", ""),
        )
        color = row.get("cluster_color", "#2F8A78")
        if not isinstance(color, str) or not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
            color = "#2F8A78"
        style = f' style="border-left:5px solid {color};background:#f1f7ed"' if row.get("cluster_highlight") else ""
        output.append(f"<tr{style}>" + "".join(f"<td>{_escape(cell)}</td>" for cell in cells) + "</tr>")
    output.append("</tbody></table></div></section>")
    for row in rows:
        output.append(f'<section><h2>{_escape(row.get("sample_name", "Unnamed sample"))}</h2><dl>')
        for title, key in (
            ("Sample ID", "sample_id"), ("Input", "input_path"), ("Input SHA-256", "input_sha256"),
            ("Scheme", "scheme"), ("Scheme SHA-256", "scheme_digest"),
            ("Job state", "job_status"), ("Typing result", "status"), ("Sequence type", "st"),
            ("Organism", "organism"), ("Organism source", "organism_source"),
            ("AMR genes", "amr_genes"), ("Virulence genes", "virulence_genes"),
            ("Plasmid replicons", "plasmid_replicons"), ("Group", "cluster_label"),
            ("HYDRA sample", "hydra_source_sample"), ("HYDRA report SHA-256", "hydra_report_sha256"),
            ("AMR evidence state", "hydra_evidence_status"), ("AMR evidence note", "hydra_evidence_reason"),
        ):
            output.append(f"<dt>{title}</dt><dd>{_escape(row.get(key)) or '—'}</dd>")
        output.append("</dl>")
        output.append(characterization_html(row))
        if row.get("missing_input"):
            output.append('<p class="notice">The original input is currently unavailable; '
                          'this report contains the stored result snapshot.</p>')
        for title, key in (("Error", "error"), ("Notes", "notes"), ("Metadata", "metadata"),
                           ("Quality summary", "qc")):
            if row.get(key):
                if key == "metadata" and row.get("hydra_evidence_status") == "stale":
                    title = "Archived metadata (old AMR evidence, not current calls)"
                shown = row[key]
                if key == 'metadata':
                    shown = {field: value for field, value in row[key].items() if field not in {'characterization', 'hydra'}}
                output.append(f"<h3>{title}</h3><pre>{_escape(shown)}</pre>")
        additional = [profile for profile in row.get("additional_profiles", [])
                      if profile.get("scheme_digest") != row.get("scheme_digest")]
        if additional:
            output.append("<h3>Additional typing snapshots</h3><p>These profiles do not replace the primary MLST result. Novel sequence identifiers are local evidence, not centrally registered allele numbers.</p>")
            for profile in additional:
                counts = profile.get("alleles", {})
                output.append(f"<h4>{_escape(profile.get('scheme'))}</h4><p>Called {sum(value is not None for value in counts.values())}/{len(counts)} loci · {_escape(profile.get('status'))}<br>Scheme SHA-256: {_escape(profile.get('scheme_digest'))}</p>")
                output.append("<details><summary>Complete stored profile and provenance</summary><pre>" + _escape(profile) + "</pre></details>")
        alleles = row.get("alleles", {})
        if alleles:
            statuses = {
                str(call.get("locus")): call.get("status", "")
                for call in row.get("calls", []) if isinstance(call, Mapping)
            }
            output.append("<h3>Allele calls</h3><table><thead><tr><th>Locus</th><th>Allele</th>"
                          "<th>Call status</th></tr></thead><tbody>")
            for locus, allele in sorted(alleles.items()):
                output.append(f"<tr><td>{_escape(locus)}</td><td>"
                              f"{_escape(allele) if allele is not None else 'Unresolved'}</td>"
                              f"<td>{_escape(statuses.get(locus, ''))}</td></tr>")
            output.append("</tbody></table>")
        output.append("</section>")
    output.append("<footer>Generated locally by WMLSTudio. This report contains no external "
                  "resources or tracking.</footer></body></html>\n")
    with _atomic_text(path) as handle:
        handle.write("".join(output))
    return Path(path).expanduser().resolve()


def export_results(
    records: Iterable[Mapping[str, Any]], destination: str | Path, format: str | None = None,
    *, selected_ids=None, highlight_clusters=None, investigation=None,
) -> Path:
    """Export CSV, TSV, JSON, or HTML, inferring format from the suffix by default."""
    selected = (format or Path(destination).suffix.lstrip(".")).lower()
    writers = {"csv": write_csv, "tsv": write_tsv, "json": write_json,
               "html": write_html, "htm": write_html}
    if selected not in writers:
        raise ValueError("Choose a CSV, TSV, JSON, or HTML export format.")
    return writers[selected](records, destination, selected_ids=selected_ids,
                             highlight_clusters=highlight_clusters, investigation=investigation)
