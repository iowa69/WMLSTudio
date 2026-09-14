"""Conservative allele distances and deterministic minimum spanning forests.

Distances count unequal shared, unambiguous known or validated full-CDS novel
alleles. Novel identities must be content-qualified full SHA-256 tokens. Each edge carries
its denominator; a missing call is never silently counted as a match. These
distances are descriptive and do not establish an outbreak or transmission.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from itertools import combinations
from typing import Any

from .sequence import check_cancelled

_UNKNOWN = frozenset({"", "-", "?", "0", "na", "n/a", "none", "unknown", "missing", "novel"})
_EXCLUDED_CALLS = frozenset({"mixed", "ambiguous", "missing", "unknown", "novel", "partial"})


def _sample_key(result: Mapping[str, Any]) -> str:
    value = result.get("sample_id") or result.get("id") or result.get("sample_name")
    if value is None or not str(value).strip():
        raise ValueError("Every comparison result needs a sample_id, id, or sample_name.")
    return str(value)


def _known_calls(result: Mapping[str, Any]) -> dict[str, str]:
    excluded = {
        str(call.get("locus"))
        for call in result.get("calls", [])
        if isinstance(call, Mapping) and str(call.get("status", "")).lower() in _EXCLUDED_CALLS
    }
    known = {}
    validated_novel = {
        str(call.get("locus")): call
        for call in result.get("calls", [])
        if isinstance(call, Mapping) and call.get("status") == "novel_validated"
        and isinstance(call.get("cds_qc"), Mapping) and call["cds_qc"].get("valid") is True
    }
    for locus, allele in result.get("alleles", {}).items():
        if str(locus) in excluded or allele is None or isinstance(allele, (list, dict, tuple, bool)):
            continue
        value = str(allele).strip()
        if value.startswith("NOVEL_"):
            evidence = validated_novel.get(str(locus), {})
            if (not re.fullmatch(r"NOVEL_[0-9a-f]{64}", value)
                    or evidence.get("sequence_sha256") != value.removeprefix("NOVEL_")
                    or evidence.get("allele") != value):
                continue
        if value.lower() in _UNKNOWN or any(separator in value for separator in (",", ";", "|", "/")):
            continue
        known[str(locus)] = value
    return known


def pairwise_distances(
    results: Iterable[Mapping[str, Any]], min_overlap: float = 0.95, *, cancelled=None,
    pair_keys=None,
) -> list[dict[str, Any]]:
    """Return every unordered pair, including explicit non-comparable pairs.

    ``overlap = shared_loci / total_loci`` where total_loci is the union of
    profile locus keys, including missing calls. Both scheme name and nonempty
    scheme digest must agree. Rejected pairs have ``distance=None`` and a reason.
    A caller can request a lower overlap, but a pair with no shared calls is
    never comparable. Equal distances may use different shared loci.
    """
    if isinstance(min_overlap, bool) or not math.isfinite(min_overlap) or not 0 <= min_overlap <= 1:
        raise ValueError("Minimum overlap must be a finite number between 0 and 1.")
    prepared = []
    seen = set()
    for result in results:
        check_cancelled(cancelled)
        key = _sample_key(result)
        if key in seen:
            raise ValueError(f"Duplicate sample identifier in comparison: {key}")
        seen.add(key)
        alleles = result.get("alleles", {})
        if not isinstance(alleles, Mapping):
            raise ValueError(f"Alleles must be a mapping for sample {key}.")
        prepared.append((key, result, _known_calls(result)))
    prepared.sort(key=lambda entry: entry[0])
    pairs = []
    for index, ((left_id, left, left_calls), (right_id, right, right_calls)) in enumerate(combinations(prepared, 2)):
        if index % 32 == 0:
            check_cancelled(cancelled)
        if pair_keys is not None and (left_id, right_id) not in pair_keys:
            continue
        loci = set(map(str, left.get("alleles", {}))) | set(map(str, right.get("alleles", {})))
        shared = sorted(left_calls.keys() & right_calls.keys())
        overlap = len(shared) / len(loci) if loci else 0.0
        reason = ""
        if not left.get("scheme") or left.get("scheme") != right.get("scheme"):
            reason = "Different or missing typing schemes."
        elif not left.get("scheme_digest") or left.get("scheme_digest") != right.get("scheme_digest"):
            reason = "Different or missing scheme fingerprints."
        elif any(str(item.get("status", "")).lower() == "mixed" for item in (left, right)):
            reason = "A sample is flagged as mixed."
        elif not shared:
            reason = "No shared unambiguous known alleles."
        elif overlap < min_overlap:
            reason = "Shared callable loci are below the minimum overlap."
        distance = None if reason else sum(left_calls[locus] != right_calls[locus] for locus in shared)
        pairs.append({
            "source": left_id,
            "target": right_id,
            "source_name": str(left.get("sample_name", left_id)),
            "target_name": str(right.get("sample_name", right_id)),
            "distance": distance,
            "shared_loci": len(shared),
            "total_loci": len(loci),
            "overlap": overlap,
            "comparable": not reason,
            "reason": reason,
        })
    return pairs


def minimum_spanning_forest(
    results: Iterable[Mapping[str, Any]], min_overlap: float = 0.95,
) -> list[dict[str, Any]]:
    """Select Kruskal edges, leaving incompatible/isolated samples disconnected.

    Equal distance edges are ordered by stable source and target identifiers.
    Isolated samples remain in the caller's input; no fabricated zero-distance
    edges are added to connect them.
    """
    return forest_from_distances(pairwise_distances(results, min_overlap))


def forest_from_distances(pairs: Iterable[Mapping[str, Any]], *, cancelled=None) -> list[dict[str, Any]]:
    """Reuse precomputed pairwise rows without repeating O(samples² × loci) work."""
    pairs = list(pairs)
    check_cancelled(cancelled)
    parent = {item[key]: item[key] for item in pairs for key in ("source", "target")}

    def root(node: str) -> str:
        while node != parent[node]:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    selected = []
    candidates = sorted(
        (pair for pair in pairs if pair["comparable"]),
        key=lambda pair: (pair["distance"], pair["source"], pair["target"]),
    )
    for index, pair in enumerate(candidates):
        if index % 256 == 0:
            check_cancelled(cancelled)
        left, right = root(pair["source"]), root(pair["target"])
        if left != right:
            parent[right] = left
            selected.append(pair)
    return selected
