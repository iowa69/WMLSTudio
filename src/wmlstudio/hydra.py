"""Import HYDRA v1 JSON reports without running engines or opening reported inputs.

Schema inspected from HYDRA 1.4.0, commit
6d36c109491c16544e8919fe6962b4b62e97d3d7, report/writer.py and records.py.
The report remains external evidence; importing it does not recalculate calls.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import UTC, datetime
from pathlib import Path

UPSTREAM_SCHEMA_REVISION = "6d36c109491c16544e8919fe6962b4b62e97d3d7"
MAX_REPORT_BYTES = 64 * 1024 * 1024
# Every element type the engine puts in a hit row, in the order a reader is
# shown them. An element type absent from a report is either an element type
# nothing searched for or one nothing was found for, and those are different
# answers; keeping the whole list here is what lets a caller say which.
ELEMENT_TYPES = ("AMR", "VIRULENCE", "STRESS", "PLASMID")
ELEMENT_TITLES = {
    "AMR": "acquired resistance genes",
    "VIRULENCE": "virulence genes",
    "STRESS": "biocide, metal and other stress-response genes",
    "PLASMID": "plasmid replicon types",
}


class HydraImportError(ValueError):
    """The selected file cannot be read as a supported HYDRA JSON report."""


def _object(value: object, label: str) -> dict:
    if not isinstance(value, dict):
        raise HydraImportError(f"{label} must be a JSON object.")
    return value


def _list(value: object, label: str) -> list:
    if not isinstance(value, list):
        raise HydraImportError(f"{label} must be a JSON list.")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str):
        raise HydraImportError(f"{label} must be text.")
    return value


def _number(value: object, label: str, minimum: float, maximum: float | None = None) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise HydraImportError(f"{label} must be a finite number.")
    if value < minimum or (maximum is not None and value > maximum):
        raise HydraImportError(f"{label} is outside its valid range.")


def _unique_keys(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise HydraImportError(f"Duplicate JSON key {key!r}; report is ambiguous.")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise HydraImportError(f"Invalid JSON numeric constant {value}.")


def _finite_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise HydraImportError("JSON number is outside the supported finite range.")
    return result


def is_point_mutation(hit: dict) -> bool:
    """A catalogued resistance point mutation, as the summary counts one.

    ``VARIANTR`` is the engine's label for a read-level variant it could not
    match to a catalogue entry, so it is not a catalogued mutation and is
    excluded here. One predicate, used by the summary and by every caller that
    shows mutations, so the count and the list can never disagree.
    """
    return (hit.get("element_type") == "AMR" and hit.get("resolution") == "POINT"
            and hit.get("method") != "VARIANTR")


def _validate_sample(value: object, index: int, warnings: list[str]) -> dict:
    sample = _object(value, f"Sample {index}")
    label = _text(sample.get("sample"), f"Sample {index} name")
    if not label.strip():
        raise HydraImportError(f"Sample {index} has an empty name.")
    hits = _list(sample.get("hits"), f"{label}: hits")
    for field in ("species", "mlst", "scores", "qc"):
        if field in sample:
            _object(sample[field], f"{label}: {field}")
    for field in ("name", "genus", "species", "confidence", "evidence"):
        if field in sample.get("species", {}):
            _text(sample["species"][field], f"{label}: species {field}")
    for field in ("scheme", "sequence_type", "source", "note"):
        if field in sample.get("mlst", {}):
            _text(sample["mlst"][field], f"{label}: MLST {field}")
    for field in ("inputs", "warnings"):
        for text in _list(sample.get(field, []), f"{label}: {field}"):
            _text(text, f"{label}: {field} entry")
    for item in _list(sample.get("typing", []), f"{label}: typing"):
        _object(item, f"{label}: typing entry")
    if "runtime_seconds" in sample:
        _number(sample["runtime_seconds"], f"{label}: runtime_seconds", 0)
    mlst = sample.get("mlst", {})
    if "alleles" in mlst:
        alleles = _object(mlst["alleles"], f"{label}: MLST alleles")
        if any(not isinstance(allele, (str, int)) or isinstance(allele, bool)
               for allele in alleles.values()):
            raise HydraImportError(f"{label}: MLST allele values must be text or integers.")
    for hit_index, hit_value in enumerate(hits, 1):
        hit = _object(hit_value, f"{label}: hit {hit_index}")
        for field in ("gene", "database", "element_type", "method"):
            _text(hit.get(field), f"{label}: hit {hit_index} {field}")
        for field in ("note", "class", "subclass", "resolution", "product", "accession"):
            if field in hit:
                _text(hit[field], f"{label}: hit {hit_index} {field}")
        if "sample" in hit and hit["sample"] != label:
            raise HydraImportError(f"{label}: hit {hit_index} belongs to a different sample.")
        if "primary" in hit and not isinstance(hit["primary"], bool):
            raise HydraImportError(f"{label}: hit {hit_index} primary must be true or false.")
        if "primary" not in hit:
            warnings.append(f"{label}: hit {hit_index} omits primary; excluded from summary counts.")
        for field in ("coverage_pct", "identity_pct"):
            if field in hit:
                _number(hit[field], f"{label}: {field}", 0, 100)
        if hit.get("allele_fraction") is not None:
            _number(hit["allele_fraction"], f"{label}: allele_fraction", 0, 1)
        if hit.get("depth") is not None:
            _number(hit["depth"], f"{label}: depth", 0)
    primary = [hit for hit in hits if hit.get("primary") is True]
    summary = {
        name: len({hit["gene"] for hit in primary if hit["element_type"] == kind})
        for name, kind in (("amr_genes", "AMR"), ("virulence_genes", "VIRULENCE"),
                           ("plasmid_replicons", "PLASMID"), ("stress_genes", "STRESS"))
    }
    summary["point_mutations"] = sum(is_point_mutation(hit) for hit in primary)
    summary["heteroresistant_sites"] = sum(
        hit["method"] == "POINTR" and hit["element_type"] == "AMR"
        and "HETERORESISTANT" in str(hit.get("note", "")).upper() for hit in hits)
    summary["total_hits"] = len(hits)
    # Reserved adapter fields are not silently allowed to overwrite upstream evidence.
    if "summary" in sample:
        sample["upstream_summary"] = sample["summary"]
    sample["summary"] = summary
    return sample


def load_hydra_report(path: str | Path) -> dict:
    """Return upstream report data plus import provenance and primary-aware counts.

    Unknown versions retaining the v1 shape are accepted with a visible warning.
    Duplicate sample IDs, conflicting hit owners and invalid quantities fail the
    whole import. File paths reported by HYDRA are retained as text only.
    """
    path = Path(path).resolve()
    try:
        with path.open("rb") as handle:
            raw = handle.read(MAX_REPORT_BYTES + 1)
        if len(raw) > MAX_REPORT_BYTES:
            raise HydraImportError("HYDRA report exceeds the 64 MiB import limit.")
        report = json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_unique_keys,
                            parse_constant=_invalid_constant, parse_float=_finite_float)
    except (OSError, UnicodeError, json.JSONDecodeError, RecursionError) as error:
        raise HydraImportError(f"Cannot read HYDRA JSON report: {error}") from error
    report = _object(report, "HYDRA report")
    if "hydra_version" not in report:
        raise HydraImportError("Not a HYDRA report: hydra_version is missing.")
    if any(key in report for key in ("import_provenance", "import_warnings")):
        raise HydraImportError("Import the original HYDRA JSON, not an already normalized report.")
    version = report["hydra_version"]
    if version is not None and not isinstance(version, str):
        raise HydraImportError("hydra_version must be text or null.")
    samples = _list(report.get("samples"), "HYDRA samples")
    _object(report.get("parameters", {}), "HYDRA parameters")
    _text(report.get("command", ""), "HYDRA command")
    for database in _list(report.get("databases", []), "HYDRA databases"):
        _text(database, "HYDRA database name")
    warnings = []
    if version != "1.4.0":
        warnings.append(f"HYDRA version {version or 'unspecified'} has not been validated; "
                        "the report matches the supported v1 structure.")
    names = set()
    for index, sample in enumerate(samples, 1):
        validated = _validate_sample(sample, index, warnings)
        if validated["sample"] in names:
            raise HydraImportError(f"Duplicate HYDRA sample name: {validated['sample']!r}.")
        names.add(validated["sample"])
    report["import_provenance"] = {
        "source_path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
        "imported_utc": datetime.now(UTC).isoformat(), "adapter_version": 1,
        "upstream_schema_revision": UPSTREAM_SCHEMA_REVISION,
        "analysis_origin": "external_hydra_report",
    }
    report["import_warnings"] = warnings
    return report
