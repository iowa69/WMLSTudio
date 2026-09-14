"""Reviewed bulk organism assignment from a spreadsheet, with re-filing.

Organism assignment is a workflow action with result-invalidation and file-moving
consequences, so it deliberately does not travel through the epidemiology
annotation path: metadata_ingest.PROTECTED_FIELDS stays exactly as it is.

A spreadsheet is an operator statement, never genomic evidence. Rows applied from
here are recorded with basis 'csv_import' and confidence 'unresolved' so that a
typed-in label can never be displayed as a reference-supported identification.
Every row resolves to exactly one target or is reported as a problem; an
ambiguous name is never resolved by guessing.
"""

from __future__ import annotations

import csv
from pathlib import Path

from wmlstudio.export import _atomic_text, _csv_cell
from wmlstudio.project import Project
from wmlstudio.sequence import check_cancelled
from wmlstudio.storage import (
    QUARANTINE_BUCKETS,
    confirm_organism,
    quarantine_samples,
    quarantine_token,
    refile_samples,
)

IDENTITY_COLUMNS = ("sample_id", "sample_name", "file")
ASSIGNMENT_COLUMNS = ("genus", "species", "typing_mode", "scheme_path", "quarantine")
TEMPLATE_COLUMNS = ("sample_id", "sample_name", "file", *ASSIGNMENT_COLUMNS)
TYPING_MODES = {"auto", "manual", "unknown"}
MAX_ROWS = 20000
MAX_COLUMNS = 50
MAX_CELL = 4096
ASSIGNMENT_BASIS = "csv_import"
ASSIGNMENT_CONFIDENCE = "unresolved"


def _read_table(path) -> list[dict]:
    source = Path(path).expanduser()
    if source.stat().st_size > 20 * 1024 * 1024:
        raise ValueError("Organism assignment tables must be at most 20 MiB.")
    with source.open(encoding="utf-8-sig", newline="") as handle:
        head = handle.read(8192)
        handle.seek(0)
        delimiter = "\t" if source.suffix.casefold() in {".tsv", ".tab"} else ","
        if source.suffix.casefold() not in {".csv", ".tsv", ".tab"}:
            try:
                delimiter = csv.Sniffer().sniff(head, delimiters=",\t").delimiter
            except csv.Error as error:
                raise ValueError("Choose a CSV or TSV file with a header row.") from error
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        headers = [value.strip().casefold().replace(" ", "_") for value in next(reader, [])]
        if not headers or len(headers) > MAX_COLUMNS or len(headers) != len(set(headers)) \
                or any(not header for header in headers):
            raise ValueError("The header must contain unique, non-empty column names.")
        if not set(IDENTITY_COLUMNS).intersection(headers):
            raise ValueError("Include one of sample_id, sample_name or file to say which isolate each row is.")
        rows = []
        for number, values in enumerate(reader, 2):
            if not values or not any(value.strip() for value in values):
                continue
            if len(values) != len(headers):
                raise ValueError(f"Row {number}: expected {len(headers)} cells, found {len(values)}.")
            if len(rows) >= MAX_ROWS or any(len(value) > MAX_CELL for value in values):
                raise ValueError("The table exceeds the row or cell-size limit.")
            rows.append({"row_number": number, **dict(zip(headers, values, strict=True))})
    return rows


def _identity_column(rows, headers) -> str:
    present = [column for column in IDENTITY_COLUMNS if column in headers]
    if len(present) > 1:
        supplied = [column for column in present
                    if any(str(row.get(column, "")).strip() for row in rows)]
        present = supplied or present[:1]
    if len(present) != 1:
        raise ValueError("Use exactly one identity column: sample_id, sample_name or file.")
    return present[0]


def _targets(samples, paths):
    by_id = {}
    by_name: dict[str, list[str]] = {}
    for sample in samples or []:
        by_id[sample["id"]] = sample
        by_name.setdefault(str(sample.get("name") or ""), []).append(sample["id"])
    by_path: dict[str, list[str]] = {}
    for value in paths or []:
        resolved = Path(value).expanduser().resolve()
        by_path.setdefault(str(resolved), [str(resolved)])
        by_path.setdefault(resolved.name.casefold(), []).append(str(resolved))
    return by_id, by_name, by_path


def read_assignments(path, *, samples=None, paths=None) -> tuple[list[dict], list[str]]:
    """Parse a CSV/TSV of organism assignments; nothing is written by this function.

    Returns (rows, problems). A row that cannot be matched to exactly one target,
    or that asks for an organism the storage layer would refuse, becomes a problem
    and is left out of the returned rows.
    """
    table = _read_table(path)
    problems: list[str] = []
    if not table:
        return [], ["The table has a header but no rows."]
    headers = {key for row in table for key in row if key != "row_number"}
    identity = _identity_column(table, headers)
    by_id, by_name, by_path = _targets(samples, paths)
    rows, claimed = [], {}
    for row in table:
        number = row["row_number"]
        key = str(row.get(identity, "")).strip()
        entry = {"row_number": number, "identity": identity, "key": key, "sample_id": None,
                 "path": None, "genus": str(row.get("genus", "")).strip(),
                 "species": str(row.get("species", "")).strip(),
                 "typing_mode": str(row.get("typing_mode", "")).strip().casefold() or "manual",
                 "scheme_path": str(row.get("scheme_path", "")).strip() or None,
                 "quarantine": str(row.get("quarantine", "")).strip().casefold() or None}
        if not key:
            problems.append(f"Row {number}: {identity} is empty.")
            continue
        if identity == "sample_id":
            if key not in by_id:
                problems.append(f"Row {number}: unknown sample_id {key!r}; a supplied ID is never replaced by a guessed name.")
                continue
            entry["sample_id"] = key
        elif identity == "sample_name":
            matches = by_name.get(key, [])
            if len(matches) != 1:
                problems.append(f"Row {number}: sample name {key!r} is "
                                + ("ambiguous; use sample_id." if matches else "not in this project."))
                continue
            entry["sample_id"] = matches[0]
        else:
            candidates = by_path.get(str(Path(key).expanduser().resolve())) or by_path.get(key.casefold())
            candidates = sorted(set(candidates or []))
            if len(candidates) != 1:
                problems.append(f"Row {number}: file {key!r} is "
                                + ("ambiguous; use the full path." if candidates else "not in this import."))
                continue
            entry["path"] = candidates[0]
        target = entry["sample_id"] or entry["path"]
        if target in claimed:
            problems.append(f"Row {number}: {key!r} was already assigned on row {claimed[target]}.")
            continue
        if entry["typing_mode"] not in TYPING_MODES:
            problems.append(f"Row {number}: typing_mode must be one of {', '.join(sorted(TYPING_MODES))}.")
            continue
        try:
            entry["quarantine"] = quarantine_token(entry["quarantine"])
        except ValueError:
            problems.append(f"Row {number}: quarantine must be one of {', '.join(sorted(QUARANTINE_BUCKETS))}.")
            continue
        if not entry["quarantine"] and not entry["genus"]:
            problems.append(f"Row {number}: give a genus, or a quarantine reason to send this isolate to Needs review.")
            continue
        if entry["quarantine"] and entry["genus"]:
            problems.append(f"Row {number}: an isolate cannot be both assigned to {entry['genus']} and quarantined.")
            continue
        claimed[target] = number
        rows.append(entry)
    return rows, problems


def apply_assignments(project: Project, rows, *, storage_root=None, cancelled=None,
                      progress=None) -> dict:
    """Apply parsed rows to project samples, then move the managed copies to match.

    Every row is resolved before anything is written, so a table with one unusable
    row changes nothing. A corrected label that did not move the file would leave
    the file sitting in a folder that contradicts it, so re-filing is not optional.
    """
    rows = list(rows)
    unresolved = [row for row in rows if not row.get("sample_id")]
    if unresolved:
        raise ValueError("Rows matched to a file rather than to an isolate must be applied "
                         "through the import review, not to an existing project.")
    known = {sample["id"] for sample in project.samples()}
    missing = [row["key"] for row in rows if row["sample_id"] not in known]
    if missing:
        raise ValueError(f"These isolates are no longer in the project: {', '.join(sorted(missing)[:5])}")
    grouped: dict[tuple, list[str]] = {}
    quarantined: dict[str, list[str]] = {}
    for row in rows:
        if row["quarantine"]:
            quarantined.setdefault(row["quarantine"], []).append(row["sample_id"])
        else:
            key = (row["genus"], row["species"], row["typing_mode"], row["scheme_path"])
            grouped.setdefault(key, []).append(row["sample_id"])
    with project.transaction():
        for (genus, species, typing_mode, scheme_path), identifiers in grouped.items():
            check_cancelled(cancelled)
            confirm_organism(project, identifiers, genus, species, scheme_path=scheme_path,
                             typing_mode=typing_mode, basis=ASSIGNMENT_BASIS,
                             confidence=ASSIGNMENT_CONFIDENCE)
        for reason, identifiers in quarantined.items():
            quarantine_samples(project, identifiers, reason)
    identifiers = [row["sample_id"] for row in rows]
    filing = refile_samples(project, identifiers, storage_root=storage_root,
                            cancelled=cancelled, progress=progress)
    return {"assigned": sum(len(value) for value in grouped.values()),
            "quarantined": sum(len(value) for value in quarantined.values()),
            "sample_ids": identifiers, "filing": filing}


def write_template(path, targets) -> Path:
    """Write the current proposals as a spreadsheet the user can edit offline."""
    destination = Path(path).expanduser().resolve()
    with _atomic_text(destination, newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(TEMPLATE_COLUMNS)
        for target in targets:
            proposed = target.get("proposed") if isinstance(target.get("proposed"), dict) else {}
            row = {"sample_id": target.get("sample_id") or target.get("id") or "",
                   "sample_name": target.get("sample_name") or target.get("name") or "",
                   "file": target.get("file") or target.get("input_path") or target.get("path") or "",
                   "genus": target.get("genus", proposed.get("genus", "")),
                   "species": target.get("species", proposed.get("species", "")),
                   "typing_mode": target.get("typing_mode") or "manual",
                   "scheme_path": target.get("scheme_path") or "",
                   "quarantine": target.get("quarantine") or ""}
            writer.writerow([_csv_cell(row[column]) for column in TEMPLATE_COLUMNS])
    return destination
