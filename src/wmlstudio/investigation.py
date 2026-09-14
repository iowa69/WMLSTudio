"""Versioned local investigations; core-profile similarity is not transmission.

All thresholds are explicit local protocol parameters. AMR, virulence and other
annotations are context, never additional cgMLST loci. No sequence file is read
by this module: immutable profile signatures allow incremental cohort updates.
"""

from __future__ import annotations

import hashlib
import json
import math
import uuid
from copy import deepcopy
from datetime import datetime, timezone
from functools import wraps
from itertools import combinations

from .comparison import _known_calls, _sample_key, pairwise_distances
from .sequence import check_cancelled

METRIC_VERSION = "shared-unambiguous-alleles-v1"
INTERPRETATION = (
    "Single-linkage clusters describe allele similarity, not proven transmission. "
    "A chain can join isolates whose direct distance exceeds the threshold. "
    "Thresholds require organism-, scheme- and protocol-specific justification. "
    "AMR/accessory features are reported separately and never added to the core distance."
)


def _now():
    return datetime.now(timezone.utc).isoformat()


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def profile_signature(result):
    """Fingerprint only evidence that can affect distances, not display metadata."""
    return _digest({"metric": METRIC_VERSION, "scheme": result.get("scheme"),
                    "digest": result.get("scheme_digest"), "status": result.get("status"),
                    "input_sha256": result.get("input_sha256"),
                    "loci": sorted(map(str, result.get("alleles", {}))),
                    "called": _known_calls(result)})


def incremental_distances(results, min_overlap=0.95, previous=None, *, cancelled=None):
    """Reuse unchanged pairs across sessions; recompute changed/new profiles only.

    Checksummed caches are an optimization, never authoritative input evidence.
    An invalid cache is discarded. Renames and annotations do not invalidate
    biological distances; input/profile changes do.
    """
    results = list(results)
    records = {_sample_key(r): r for r in results}
    if len(records) != len(results):
        raise ValueError("Duplicate sample identifiers in investigation.")
    signatures = {}
    for sid, result in records.items():
        check_cancelled(cancelled)
        signatures[sid] = profile_signature(result)
    prior = previous or {}
    payload = {key: prior.get(key) for key in ("metric_version", "min_overlap", "signatures", "rows")}
    valid = (prior.get("metric_version") == METRIC_VERSION and prior.get("min_overlap") == min_overlap
             and prior.get("checksum") == _digest(payload))
    old_signatures = prior.get("signatures", {}) if valid else {}
    old_rows = {(p["source"], p["target"]): p for p in prior.get("rows", [])} if valid else {}
    reused, required = {}, set()
    for index, pair in enumerate(combinations(sorted(records), 2)):
        if index % 256 == 0:
            check_cancelled(cancelled)
        if (pair in old_rows and all(old_signatures.get(sid) == signatures[sid] for sid in pair)):
            row = dict(old_rows[pair])
            row.update(source_name=str(records[pair[0]].get("sample_name", pair[0])),
                       target_name=str(records[pair[1]].get("sample_name", pair[1])))
            reused[pair] = row
        else:
            required.add(pair)
    fresh = pairwise_distances(results, min_overlap, cancelled=cancelled, pair_keys=required)
    combined = {**reused, **{(p["source"], p["target"]): p for p in fresh}}
    rows = [combined[key] for key in sorted(combined)]
    cache = {"metric_version": METRIC_VERSION, "min_overlap": min_overlap,
             "signatures": signatures, "rows": rows}
    cache["checksum"] = _digest(cache)
    return rows, cache, {"reused_pairs": len(reused), "computed_pairs": len(fresh),
                         "new_profiles": sorted(set(records) - set(old_signatures)),
                         "changed_profiles": sorted(sid for sid in records if sid in old_signatures
                                                    and signatures[sid] != old_signatures[sid])}


def threshold_clusters(results, rows, threshold, *, previous=None, cancelled=None):
    """All-pair single linkage with stable IDs and explicit split/merge lineage."""
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 0:
        raise ValueError("Cluster threshold must be a nonnegative integer.")
    results = list(results)
    records = {_sample_key(r): r for r in results}
    if len(records) != len(results):
        raise ValueError("Duplicate sample identifiers in investigation.")
    adjacency = {sid: set() for sid in records}
    comparable = {sid: set() for sid in records}
    lookup = {}
    for row in rows:
        check_cancelled(cancelled)
        a, b = str(row["source"]), str(row["target"])
        if a not in records or b not in records or a == b:
            raise ValueError("Pairwise evidence does not belong to this investigation cohort.")
        lookup[frozenset((a, b))] = row
        if row.get("comparable"):
            distance = row.get("distance")
            if isinstance(distance, bool) or not isinstance(distance, int) or distance < 0:
                raise ValueError("Comparable pairs require a nonnegative integer allele distance.")
            comparable[a].add(b)
            comparable[b].add(a)
            if distance <= threshold:
                adjacency[a].add(b)
                adjacency[b].add(a)
    components, remaining = [], set(records)
    while remaining:
        pending, members = [min(remaining)], set()
        while pending:
            sid = pending.pop()
            if sid not in members:
                members.add(sid)
                pending.extend(adjacency[sid] - members)
        remaining.difference_update(members)
        components.append(sorted(members))
    old = [g for g in (previous or {}).get("groups", []) if g.get("status") == "cluster"]
    next_number = max([int((previous or {}).get("next_cluster_number", 1)),
                       *[int(g.get("number", 0)) + 1 for g in old]])
    old_sets = {g["id"]: set(g["members"]) for g in old}
    overlaps = {gid: sum(bool(members.intersection(group)) for group in components)
                for gid, members in old_sets.items()}
    groups = []
    for members in components:
        check_cancelled(cancelled)
        parents = [g for g in old if old_sets[g["id"]].intersection(members)]
        member_set = set(members)
        if len(members) == 1:
            status = "singleton" if comparable[members[0]] else "not_comparable"
            group_id = f"{status}:{members[0]}"
            name = "Unlinked isolate" if status == "singleton" else "Not comparable"
            number = None
            change = "split" if parents else "unlinked"
        else:
            status = "cluster"
            exact = next((g for g in parents if old_sets[g["id"]] == member_set), None)
            growing = parents[0] if len(parents) == 1 and overlaps[parents[0]["id"]] == 1 else None
            retained = exact or growing
            if retained:
                group_id, name, number = retained["id"], retained["name"], retained["number"]
                change = "unchanged" if exact else "grown" if old_sets[group_id] <= member_set else "reduced"
            else:
                number = next_number
                next_number += 1
                group_id = "cluster:" + str(uuid.uuid5(uuid.NAMESPACE_URL,
                    "wmlstudio:" + _digest([sorted(records[r].get("scheme_digest", "") for r in members), members,
                                             threshold, sorted(g["id"] for g in parents)])))
                name = f"Cluster {number:03d}"
                change = "merged" if len(parents) > 1 else "split" if parents else "new"
        internal = [lookup.get(frozenset(pair)) for pair in combinations(members, 2)]
        distances = [p["distance"] for p in internal if p and p.get("comparable")]
        unassessed = sum(not p or not p.get("comparable") for p in internal)
        max_distance = max(distances) if distances else None
        old_members = set().union(*(old_sets[g["id"]] for g in parents)) if parents else set()
        groups.append({"id": group_id, "name": name, "number": number, "status": status,
                       "members": members, "change": change, "parents": [g["id"] for g in parents],
                       "added_members": sorted(member_set - old_members),
                       "removed_members": sorted(old_members - member_set),
                       "max_direct_distance": max_distance,
                       "chained": max_distance is not None and max_distance > threshold,
                       "unassessed_within_pairs": unassessed,
                       "reason": "No accepted pairwise comparison; review missing/mixed evidence or cohort size."
                       if status == "not_comparable" else ""})
    return {"groups": groups, "next_cluster_number": next_number}


def build_snapshot(results, rows, threshold, min_overlap, *, previous=None, reuse=None, cancelled=None):
    results = list(results)
    schemes = {(r.get("scheme"), r.get("scheme_digest")) for r in results}
    if len(schemes) > 1 or any(not name or not digest for name, digest in schemes):
        raise ValueError("An investigation snapshot must use one nonempty scheme fingerprint.")
    if isinstance(min_overlap, bool) or not math.isfinite(min_overlap) or not 0 <= min_overlap <= 1:
        raise ValueError("Minimum shared-locus fraction must be between zero and one.")
    # Cluster identities are reference-scoped. A deliberately changed scheme
    # starts a new lineage; identical sample names do not bridge protocols.
    lineage = previous if previous and results and previous.get('scheme_digest') == results[0].get('scheme_digest') else None
    grouping = threshold_clusters(results, rows, threshold, previous=lineage, cancelled=cancelled)
    profiles = []
    for result in results:
        check_cancelled(cancelled)
        known = _known_calls(result)
        profiles.append({key: deepcopy(result.get(key)) for key in (
            "sample_id", "sample_name", "scheme", "scheme_digest", "status", "st", "primary_st",
            "primary_scheme", "input_sha256", "metadata", "amr_genes", "amr_evidence_status", "organism")})
        profiles[-1].update(sample_id=_sample_key(result), known_alleles=known,
                            profile_loci=sorted(map(str, result.get("alleles", {}))),
                            signature=profile_signature(result), callable_loci=len(known),
                            total_loci=len(result.get("alleles", {})))
    return {"format_version": 1, "snapshot_id": uuid.uuid4().hex, "created_at": _now(),
            "scheme": results[0].get("scheme") if results else "",
            "scheme_digest": results[0].get("scheme_digest") if results else "",
            "threshold": threshold, "min_overlap": min_overlap, "metric_version": METRIC_VERSION,
            "missing_policy": "Shared callable loci / union of profile loci; excluded pairs have no distance.",
            "interpretation": INTERPRETATION, "profiles": profiles, "pairs": deepcopy(rows),
            "reuse": deepcopy(reuse or {}), **grouping}


def proximity_rows(snapshot, sample_ids=None):
    """One nearest-neighbour summary per explicitly selected isolate, ties retained."""
    profiles = {p["sample_id"]: p for p in snapshot.get("profiles", [])}
    selected = set(profiles) if sample_ids is None else set(sample_ids)
    if not selected <= profiles.keys():
        raise ValueError("Proximity selection contains samples outside this snapshot.")
    grouping = {sid: group for group in snapshot.get("groups", []) for sid in group["members"]}
    output = []
    for sid in sorted(selected):
        neighbours, excluded = [], []
        for pair in snapshot.get("pairs", []):
            if sid not in {pair["source"], pair["target"]}:
                continue
            other = pair["target"] if pair["source"] == sid else pair["source"]
            row = dict(pair, neighbour_id=other, neighbour_name=profiles[other].get("sample_name", other))
            (neighbours if pair.get("comparable") else excluded).append(row)
        neighbours.sort(key=lambda p: (p["distance"], p["neighbour_id"]))
        best = neighbours[0]["distance"] if neighbours else None
        group = grouping.get(sid, {})
        output.append({"sample_id": sid, "sample_name": profiles[sid].get("sample_name", sid),
                       "nearest_distance": best,
                       "nearest": [p for p in neighbours if p["distance"] == best],
                       "within_threshold": [p for p in neighbours if p["distance"] <= snapshot["threshold"]],
                       "compared_count": len(neighbours), "excluded_count": len(excluded),
                       "excluded": excluded, "group_id": group.get("id"),
                       "group_name": group.get("name"), "group_status": group.get("status"),
                       "chained": group.get("chained", False)})
    return output


def threshold_guidance_status(snapshot):
    """A citation persists, but approval never silently transfers to a new policy."""
    evidence = snapshot.get('threshold_evidence') or {}
    if not evidence:
        return {'status': 'no_guidance', 'reason': 'No published threshold guidance attached; this is an exploratory local setting.'}
    if evidence.get('approved_threshold') is None:
        return {'status': 'citation_only', 'reason': 'Citation retained for context only. No numeric cutoff was adopted from it.'}
    context = evidence.get('context') or {}
    if not context.get('scheme_digest') or context.get('min_overlap') is None:
        return {'status': 'unverified_binding', 'reason': 'The recorded citation lacks a complete reference/overlap binding; it is not approval of this snapshot threshold.'}
    changed = []
    for label, source, current in [('scheme fingerprint', context['scheme_digest'], snapshot.get('scheme_digest')),
                                    ('inclusive threshold', evidence['approved_threshold'], snapshot.get('threshold')),
                                    ('minimum shared-locus fraction', context['min_overlap'], snapshot.get('min_overlap'))]:
        if source != current:
            changed.append(label)
    if changed:
        return {'status': 'changed_since_review', 'reason': 'Changed since publication review: ' + ', '.join(changed) + '. The citation is historical context, not approval of the current comparison policy.'}
    return {'status': 'matches_reviewed_context', 'reason': 'Reference, cutoff and minimum overlap match the recorded local adaptation. This does not establish clinical validation or transmission.'}


def _transactional(method):
    @wraps(method)
    def synchronized(self, *args, **kwargs):
        with self.project.transaction():
            return method(self, *args, **kwargs)
    return synchronized


class InvestigationStore:
    """Small plan catalog plus independently persisted immutable snapshots."""

    def __init__(self, project):
        self.project = project

    def list(self):
        return self.project.get_setting("investigations.index", [])

    def get(self, investigation_id):
        if investigation_id not in {entry["id"] for entry in self.list()}:
            raise KeyError(investigation_id)
        return self.project.get_setting(f"investigation.{investigation_id}")

    @_transactional
    def save(self, name, sample_ids, *, investigation_id=None, scheme_digest="", scheme="",
             scheme_path="", threshold=1, min_overlap=0.95, protocol="User-defined exploratory threshold",
             include_new=False, filters=None, label_fields=None, threshold_evidence=None):
        if not str(name).strip() or len(str(name).strip()) > 160:
            raise ValueError("Investigation names must contain 1–160 characters.")
        if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 0:
            raise ValueError("Threshold must be a nonnegative integer.")
        if isinstance(min_overlap, bool) or not math.isfinite(min_overlap) or not 0 <= min_overlap <= 1:
            raise ValueError("Minimum overlap must be between zero and one.")
        ids = sorted(set(map(str, sample_ids)))
        known = {s["id"] for s in self.project.samples()}
        if not set(ids) <= known:
            raise ValueError("Investigation contains unknown project sample IDs.")
        old = self.get(investigation_id) if investigation_id else {}
        investigation_id = investigation_id or uuid.uuid4().hex
        plan = {**old, "id": investigation_id, "name": str(name).strip(), "sample_ids": ids,
                "scheme_digest": str(scheme_digest), "scheme": str(scheme), "scheme_path": str(scheme_path),
                "threshold": threshold, "min_overlap": min_overlap, "protocol": str(protocol),
                "threshold_evidence": deepcopy(old.get('threshold_evidence', {}) if threshold_evidence is None else threshold_evidence),
                "include_new": bool(include_new), "filters": deepcopy(filters or {}),
                "label_fields": list(label_fields or ["sample_name", "primary_st"]),
                "created_at": old.get("created_at", _now()), "updated_at": _now(),
                "review_groups": old.get("review_groups", []), "snapshots": old.get("snapshots", [])}
        with self.project.transaction():
            self.project.set_setting(f"investigation.{investigation_id}", plan)
            index = [entry for entry in self.list() if entry["id"] != investigation_id]
            index.append({"id": investigation_id, "name": plan["name"], "updated_at": plan["updated_at"],
                          "sample_count": len(ids), "scheme": plan["scheme"]})
            self.project.set_setting("investigations.index", sorted(index, key=lambda x: (x["name"].casefold(), x["id"])))
            self.project.record_history(None, "investigation_saved", {"id": investigation_id, "name": plan["name"],
                                                                      "sample_ids": ids, "threshold": threshold,
                                                                      "min_overlap": min_overlap})
        return deepcopy(plan)

    def snapshot(self, investigation_id, snapshot_id=None):
        plan = self.get(investigation_id)
        snapshot_id = snapshot_id or (plan["snapshots"][-1] if plan["snapshots"] else None)
        if snapshot_id and snapshot_id not in plan["snapshots"]:
            raise KeyError(snapshot_id)
        return self.project.get_setting(f"investigation.{investigation_id}.snapshot.{snapshot_id}") if snapshot_id else None

    def cache(self, investigation_id):
        self.get(investigation_id)
        return self.project.get_setting(f"investigation.{investigation_id}.cache", {})

    @_transactional
    def save_snapshot(self, investigation_id, snapshot, cache=None, *, expected_revision=None,
                      expected_project_revision=None, cancelled=None):
        plan = self.get(investigation_id)
        check_cancelled(cancelled)
        if expected_revision is not None and plan.get('updated_at') != expected_revision:
            raise ValueError('Investigation settings changed before snapshot persistence. Build again.')
        if expected_project_revision is not None and self.project.comparison_revision() != expected_project_revision:
            raise ValueError('Profile evidence changed before snapshot persistence. Build again.')
        snapshot = deepcopy(snapshot)
        sid = snapshot["snapshot_id"]
        if sid in plan["snapshots"]:
            raise ValueError("Saved investigation snapshots are immutable.")
        if plan.get('scheme_digest') and snapshot.get('scheme_digest') != plan['scheme_digest']:
            raise ValueError('The investigation is pinned to another scheme fingerprint. Edit its protocol or save a new investigation.')
        if not plan.get('scheme_digest') and snapshot.get('scheme_digest'):
            plan['scheme_digest'], plan['scheme'] = snapshot['scheme_digest'], snapshot['scheme']
        snapshot.update(investigation_id=investigation_id, investigation_name=plan["name"],
                        protocol=plan["protocol"], review_groups=deepcopy(plan["review_groups"]),
                        threshold_evidence=deepcopy(plan.get('threshold_evidence', {})))
        with self.project.transaction():
            self.project.set_setting(f"investigation.{investigation_id}.snapshot.{sid}", snapshot)
            plan["snapshots"].append(sid)
            plan["updated_at"] = _now()
            plan["summary"] = {"snapshot_id": sid, "created_at": snapshot["created_at"],
                               "profiles": len(snapshot["profiles"]),
                               "clusters": sum(g["status"] == "cluster" for g in snapshot["groups"]),
                               "not_comparable": sum(g["status"] == "not_comparable" for g in snapshot["groups"]),
                               "reuse": snapshot.get("reuse", {})}
            self.project.set_setting(f"investigation.{investigation_id}", plan)
            if cache is not None:
                self.project.set_setting(f"investigation.{investigation_id}.cache", cache)
            self.project.record_history(None, "investigation_snapshot", {"id": investigation_id,
                "snapshot_id": sid, "sample_ids": [p["sample_id"] for p in snapshot["profiles"]],
                "scheme_digest": snapshot["scheme_digest"], "threshold": snapshot["threshold"],
                "changes": [{"id": g["id"], "change": g["change"], "parents": g["parents"]}
                            for g in snapshot["groups"]], "reuse": snapshot.get("reuse", {})})
        check_cancelled(cancelled)
        return snapshot

    @_transactional
    def freeze_group(self, investigation_id, name, sample_ids, color="#E9AD66"):
        plan = self.get(investigation_id)
        ids = sorted(set(sample_ids))
        if not str(name).strip() or not ids or not set(ids) <= {s["id"] for s in self.project.samples()}:
            raise ValueError("A frozen review group requires a name and existing sample IDs.")
        group = {"id": uuid.uuid4().hex, "name": str(name).strip(), "sample_ids": ids,
                 "color": str(color), "frozen_at": _now(), "membership_policy": "fixed; never updated by clustering"}
        with self.project.transaction():
            plan["review_groups"].append(group)
            plan['updated_at'] = _now()
            self.project.set_setting(f"investigation.{investigation_id}", plan)
            self.project.record_history(None, "investigation_review_group_frozen",
                                        {"investigation_id": investigation_id, **group})
        return deepcopy(group)
