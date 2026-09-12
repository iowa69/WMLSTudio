"""Profile interoperability, portable evidence bundles, and a persistent library.

The import rules preserve the correctness fixes documented in MLSTudio 2.0:
missing tokens and padded zeroes are not alleles, ragged rows are rejected, and
novel markers retain their identity. This implementation is independent and
keeps inferred external alleles out of exact-known-allele comparisons.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import sqlite3
import uuid
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.export import _atomic_text, _csv_cell, ensure_separate_destination
from wmlstudio.sample_workflow import feature_fields, select_records

_MISSING = frozenset({"", "0", "-", "--", "?", "N", "NA", "N/A", "NAN", "NULL", "NONE",
                      "LNF", "PLOT", "PLOT3", "PLOT5", "PLOT3'", "PLOT5'", "LOTSC",
                      "NIPH", "NIPHEM", "NIPL", "PAMA", "ASM", "ALM", "INF"})
_ZERO = re.compile(r"[+-]?0+(?:\.0+)?", re.ASCII)
_ALLELE = re.compile(r"(INF-)?([0-9]+)(~([0-9a-f]{6})?)?", re.IGNORECASE | re.ASCII)
_CONTENT_ALLELE = re.compile(r'(SHA256_|NOVEL_)([0-9a-f]{64})', re.IGNORECASE | re.ASCII)
MAX_BUNDLE_BYTES = 256 * 1024 * 1024


def _cell(raw: str) -> tuple[str | None, str | None, str]:
    value = raw.strip()
    if value.upper() in _MISSING or _ZERO.fullmatch(value):
        return None, None, "missing"
    content = _CONTENT_ALLELE.fullmatch(value)
    if content:
        # A TSV preserves an asserted sequence identity, but it carries no CDS
        # sequence/validation evidence. Keep it round-trippable, not callable.
        return None, content[1].upper() + content[2].lower(), 'external_sequence_unverified'
    match = _ALLELE.fullmatch(value)
    if match is None:
        return None, None, "unknown"
    number = match[2].lstrip("0")
    if not number:
        return None, None, "missing"
    if match[1] or match[3]:
        canonical = number + "~" + (match[4] or "").lower()
        return None, canonical, "inferred"
    return number, number, "external_exact"


def parse_profile_table(text: str, scheme_name: str, scheme_digest: str | None = None, expected_loci=None) -> dict:
    if not str(scheme_name).strip():
        raise ValueError("Choose a scheme name for the imported profile table.")
    rows = [row for row in csv.reader(io.StringIO(text.lstrip("\ufeff")), delimiter="\t") if row]
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError("A profile table needs a sample column, locus columns, and at least one sample.")
    header = [value.strip() for value in rows[0]]
    if any(not value for value in header):
        raise ValueError("Profile table column names cannot be empty.")
    st_indexes = [index for index, value in enumerate(header[1:], 1) if value.upper() in {"ST", "SEQUENCE TYPE"}]
    if len(st_indexes) > 1:
        raise ValueError("Profile table has multiple sequence-type columns.")
    positions = {}
    for index, locus in enumerate(header[1:], 1):
        if index not in st_indexes:
            positions.setdefault(locus, []).append(index)
    if not positions:
        raise ValueError("Profile table has no allele locus columns.")
    expected = list(dict.fromkeys(expected_loci or []))
    if expected and not set(expected).intersection(positions):
        raise ValueError("The profile table's loci do not match the selected scheme.")
    loci = list(dict.fromkeys([*expected, *positions]))
    digest = scheme_digest or "external-table:" + hashlib.sha256(text.encode("utf-8")).hexdigest()
    warnings = []
    if not scheme_digest:
        warnings.append("No verified scheme fingerprint supplied; comparisons are restricted to this table snapshot.")
    if len(positions) < len(header) - 1 - len(st_indexes):
        warnings.append("Repeated locus columns were merged only where their nonmissing calls agree.")
    if expected and set(expected) - positions.keys():
        warnings.append("Absent scheme columns are retained as missing in the full scheme denominator.")
    names = set()
    parsed = []
    unknown = Counter()
    inferred = 0
    for row_number, values in enumerate(rows[1:], 2):
        if len(values) != len(header):
            raise ValueError(f"Ragged profile table: row {row_number} has {len(values)} columns, expected {len(header)}.")
        name = values[0].strip()
        # Undo only the exact spreadsheet-escape prefixes produced by this tool.
        if name.startswith("'") and name[1:].lstrip().startswith(("=", "+", "-", "@")):
            name = name[1:]
        if not name or name in names:
            raise ValueError(f"Empty or duplicate profile sample name on row {row_number}: {name!r}")
        names.add(name)
        alleles, external, calls = {}, {}, []
        for locus in loci:
            candidates = [_cell(values[index]) for index in positions.get(locus, [])]
            identities = {candidate[1] for candidate in candidates if candidate[1] is not None}
            if len(identities) > 1:
                raise ValueError(f"Conflicting duplicate locus {locus!r} in sample {name!r}.")
            chosen = next((candidate for candidate in candidates if candidate[1] is not None), (None, None, "missing"))
            allele, original, status = chosen
            if status in {'inferred', 'external_sequence_unverified'}:
                inferred += 1
            for index in positions.get(locus, []):
                if _cell(values[index])[2] == "unknown":
                    unknown[values[index]] += 1
            alleles[locus], external[locus] = allele, original
            calls.append({"locus": locus, "allele": allele, "external_allele": original, "status": status})
        st = values[st_indexes[0]].strip() if st_indexes else None
        if st and not re.fullmatch(r"[0-9]+", st, re.ASCII):
            warnings.append(f"{name}: external ST {st!r} is an annotation, not a verified sequence type.")
        parsed.append({"sample_name": name, "kind": "profile", "input_path": "", "qc": {},
                       "status": "profile_imported", "scheme": scheme_name.strip(), "scheme_digest": digest,
                       "st": st or None, "alleles": alleles, "external_alleles": external, "calls": calls,
                       "notes": ["Imported external allele identifiers; sequence identity was not recalculated."],
                       "provenance": {"origin": "external_profile_table", "scheme_fingerprint_verified": bool(scheme_digest)}})
    if unknown:
        warnings.append("Unrecognised tokens imported as missing: " + ", ".join(f"{key!r} ({count})" for key, count in unknown.items()))
    if inferred:
        warnings.append(f"{inferred} inferred/novel or sequence-addressed external calls retained for round-trip export and excluded from callable-allele comparisons; a plain table does not carry full-CDS validation. Use a validated evidence bundle to retain native calling evidence.")
    for result in parsed:
        result["notes"].extend(warnings)
    return {"results": parsed, "loci": loci, "warnings": warnings, "scheme_digest": digest}


def import_profile_table(
    project, path, scheme_name: str, scheme_digest: str | None = None, expected_loci=None,
    organism: dict | None = None,
) -> dict:
    source = Path(path).resolve()
    raw = source.read_bytes()
    parsed = parse_profile_table(raw.decode("utf-8-sig"), scheme_name, scheme_digest, expected_loci)
    digest = hashlib.sha256(raw).hexdigest()
    identifiers = []
    with project.transaction():
        for result in parsed["results"]:
            result["provenance"].update({"source_path": str(source), "source_sha256": digest})
            identifiers.append(project.add_profile(result["sample_name"], result, {
                "organism": organism or {}, "provenance": result["provenance"],
                "workflow": {"typing_mode": "unknown", "source_kind": "profile"},
            }))
    return {"sample_ids": identifiers, "imported": len(identifiers), "loci": len(parsed["loci"]),
            "warnings": parsed["warnings"], "scheme_digest": parsed["scheme_digest"]}


def export_profile_table(project, path, sample_ids=None, scheme_digest=None) -> Path:
    samples = select_records(project.samples(), sample_ids)
    results = []
    for sample in samples:
        if scheme_digest:
            result = next((result for result in project.analysis_results(sample["id"])
                           if result.get("scheme_digest") == scheme_digest), None)
        else:
            result = sample.get("result")
        if not result or not result.get("alleles"):
            raise ValueError(f"{sample['name']} has no allele profile for the chosen scheme.")
        results.append((sample, result))
    if not results:
        raise ValueError("Select at least one profile to export.")
    fingerprints = {(result.get("scheme"), result.get("scheme_digest")) for _, result in results}
    if len(fingerprints) != 1 or not next(iter(fingerprints))[1]:
        raise ValueError("Profile table export requires one identical scheme fingerprint.")
    if len({sample["name"] for sample, _ in results}) != len(results):
        raise ValueError("Profile table sample names must be unique; use a portable bundle to retain duplicate names.")
    loci = list(dict.fromkeys(locus for _, result in results for locus in result["alleles"]))
    ensure_separate_destination(path, [project.path, *(sample["input_path"] for sample, _ in results)])
    with _atomic_text(path, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["FILE", *loci])
        for sample, result in results:
            values = result.get("external_alleles") or result["alleles"]
            writer.writerow([_csv_cell(sample["name"]), *(_csv_cell(values.get(locus) or "0") for locus in loci)])
    return Path(path).resolve()


def export_bundle(project, path, sample_ids=None) -> Path:
    samples = select_records(project.samples(), sample_ids)
    included = {sample["id"] for sample in samples}
    payload = {
        "format": "WMLSTudio-evidence-bundle", "format_version": 1,
        "application_version": __version__, "exported_at": datetime.now(UTC).isoformat(),
        "contains_sequences": False,
        "samples": [{"id": sample["id"], "name": sample["name"], "metadata": sample["metadata"],
                     "result": sample["result"], "analyses": project.analysis_results(sample["id"]),
                     "source_input_path": sample["input_path"]} for sample in samples],
        "collections": [{**entry, "sample_ids": [sid for sid in entry["sample_ids"] if sid in included]}
                        for entry in project.collections() if included.intersection(entry["sample_ids"])],
    }
    ensure_separate_destination(path, [project.path, *(sample["input_path"] for sample in samples)])
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, indent=2)
    with _atomic_text(path) as handle:
        handle.write(encoded + "\n")
    return Path(path).resolve()


def import_bundle(project, path) -> dict:
    source = Path(path).resolve()
    with source.open("rb") as handle:
        raw = handle.read(MAX_BUNDLE_BYTES + 1)
    if len(raw) > MAX_BUNDLE_BYTES:
        raise ValueError("Evidence bundle exceeds the 256 MiB import limit.")

    def unique(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"Duplicate bundle field: {key}")
            value[key] = item
        return value

    payload = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=unique)
    # Reject NaN/Infinity/overflow before any database mutation.
    json.dumps(payload, allow_nan=False)
    if not isinstance(payload, dict) or payload.get("format") != "WMLSTudio-evidence-bundle" or payload.get("format_version") != 1:
        raise ValueError("This file is not a supported WMLSTudio evidence bundle.")
    records = payload.get("samples")
    if not isinstance(records, list):
        raise ValueError("Bundle samples must be a list.")
    names = []
    for sample in records:
        if not isinstance(sample, dict) or not isinstance(sample.get("id"), str) or not sample["id"]:
            raise ValueError("Every bundled sample requires a stable identifier.")
        if not isinstance(sample.get("name"), str) or not sample["name"].strip():
            raise ValueError("Every bundled sample requires a name.")
        if not isinstance(sample.get("metadata", {}), dict) or not isinstance(sample.get("analyses", []), list):
            raise ValueError("Invalid sample metadata or scheme analyses in bundle.")
        for result in [sample.get("result"), *sample.get("analyses", [])]:
            if result is not None and (not isinstance(result, dict) or not isinstance(result.get("alleles", {}), dict)):
                raise ValueError("Invalid allele result in bundle.")
        names.append(sample["id"])
    if len(names) != len(set(names)):
        raise ValueError("Bundle contains duplicate sample identifiers.")
    existing = {sample["id"] for sample in project.samples()}
    mapping = {sid: uuid.uuid4().hex if sid in existing else sid for sid in names}
    digest = hashlib.sha256(raw).hexdigest()
    with project.transaction():
        for sample in records:
            metadata = sample.get("metadata", {})
            metadata["bundle"] = {"original_sample_id": sample["id"], "source_path": str(source),
                                   "source_sha256": digest, "original_input_path": sample.get("source_input_path", "")}
            result = sample.get("result") or {"sample_name": sample["name"], "kind": "profile",
                                              "status": "profile_imported", "alleles": {}, "notes": []}
            identifier = project.add_profile(sample["name"], result, metadata, sample_id=mapping[sample["id"]])
            for analysis in sample.get("analyses", []):
                project.set_analysis(identifier, analysis)
        for collection in payload.get("collections", []):
            if not isinstance(collection, dict) or not isinstance(collection.get("sample_ids"), list):
                raise ValueError("Invalid collection in evidence bundle.")
            if any(sid not in mapping for sid in collection["sample_ids"]):
                raise ValueError("Bundle collection references a sample outside the bundle.")
            collection_id = project.create_collection(collection["name"])
            project.set_collection_members(collection_id, [mapping[sid] for sid in collection["sample_ids"]])
    return {"sample_ids": list(mapping.values()), "id_mapping": mapping, "imported": len(mapping)}


class Library:
    """A cross-project search index of saved evidence, independent of input files."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.path, timeout=10)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS library_samples (project_path TEXT NOT NULL, sample_id TEXT NOT NULL, "
            "snapshot TEXT NOT NULL, indexed_at TEXT NOT NULL, PRIMARY KEY(project_path,sample_id))"
        )

    def index_project(self, project) -> int:
        samples = project.samples()
        collections = project.collections()
        with self.connection:
            self.connection.execute("DELETE FROM library_samples WHERE project_path = ?", (str(project.path),))
            for sample in samples:
                snapshot = {**sample, "analyses": project.analysis_results(sample["id"]),
                            "project_path": str(project.path), "features": feature_fields(sample),
                            "collections": [{"id": entry["id"], "name": entry["name"]}
                                            for entry in collections if sample["id"] in entry["sample_ids"]]}
                self.connection.execute(
                    "INSERT INTO library_samples VALUES (?, ?, ?, ?)",
                    (str(project.path), sample["id"], json.dumps(snapshot, ensure_ascii=False, allow_nan=False),
                     datetime.now(UTC).isoformat()),
                )
        return len(samples)

    def search(self, q: str = "", *, genus=None, species=None, st=None, amr_gene=None,
               scheme_digest=None, collection=None):
        matches = []
        for row in self.connection.execute("SELECT snapshot FROM library_samples ORDER BY project_path,sample_id"):
            sample = json.loads(row["snapshot"])
            feature = sample["features"]
            result = sample.get("result") or {}
            if q and q.casefold() not in json.dumps(sample, ensure_ascii=False).casefold():
                continue
            if genus and feature.get("genus", "").casefold() != str(genus).casefold():
                continue
            if species and feature.get("species", "").casefold() != str(species).casefold():
                continue
            if st is not None and str(result.get("st")) != str(st):
                continue
            if amr_gene and not any(str(amr_gene).casefold() in gene.casefold() for gene in feature.get("amr_genes", [])):
                continue
            if scheme_digest and not any(analysis.get("scheme_digest") == scheme_digest for analysis in sample["analyses"]):
                continue
            if collection and not any(str(collection).casefold() in {entry["id"].casefold(), entry["name"].casefold()}
                                      for entry in sample.get("collections", [])):
                continue
            matches.append(sample)
        return matches

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        self.close()
