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
)


def ensure_separate_destination(destination: str | Path, protected_paths: Iterable[str | Path]) -> None:
    """Never replace a sequence input or live project with an exported result."""
    target = Path(destination).expanduser().resolve()
    for value in protected_paths:
        if value and target == Path(value).expanduser().resolve():
            raise ValueError("Choose a different output file. This path belongs to an input or open project.")


def _snapshot(
    records: Iterable[Mapping[str, Any]], selected_ids=None, highlight_clusters=None,
) -> list[dict[str, Any]]:
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
    selected_ids=None, highlight_clusters=None,
) -> Path:
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters)
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


def write_csv(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None) -> Path:
    """Write Excel-friendly UTF-8 CSV with guarded textual formula cells."""
    return _write_delimited(records, path, ",", selected_ids, highlight_clusters)


def write_tsv(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None) -> Path:
    return _write_delimited(records, path, "\t", selected_ids, highlight_clusters)


def write_json(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None) -> Path:
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters)
    document = {
        "format_version": 1,
        "application": f"WMLSTudio {__version__}",
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "samples": rows,
        "report_scope": {"mode": "selected" if selected_ids is not None else "provided",
                         "sample_ids": [row.get("sample_id") for row in rows]},
    }
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


def write_html(records: Iterable[Mapping[str, Any]], path: str | Path, *, selected_ids=None, highlight_clusters=None) -> Path:
    """Write an escaped, self-contained report that also prints without scripts."""
    all_records = list(records)
    ensure_separate_destination(path, _protected_paths(all_records))
    rows = _snapshot(all_records, selected_ids, highlight_clusters)
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
        if row.get("missing_input"):
            output.append('<p class="notice">The original input is currently unavailable; '
                          'this report contains the stored result snapshot.</p>')
        for title, key in (("Error", "error"), ("Notes", "notes"), ("Metadata", "metadata"),
                           ("Quality summary", "qc")):
            if row.get(key):
                if key == "metadata" and row.get("hydra_evidence_status") == "stale":
                    title = "Archived metadata (old AMR evidence, not current calls)"
                output.append(f"<h3>{title}</h3><pre>{_escape(row[key])}</pre>")
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
    *, selected_ids=None, highlight_clusters=None,
) -> Path:
    """Export CSV, TSV, JSON, or HTML, inferring format from the suffix by default."""
    selected = (format or Path(destination).suffix.lstrip(".")).lower()
    writers = {"csv": write_csv, "tsv": write_tsv, "json": write_json,
               "html": write_html, "htm": write_html}
    if selected not in writers:
        raise ValueError("Choose a CSV, TSV, JSON, or HTML export format.")
    return writers[selected](records, destination, selected_ids=selected_ids,
                             highlight_clusters=highlight_clusters)
