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

from .comparison import _known_calls, _sample_key, forest_from_distances, pairwise_distances
from .sequence import check_cancelled

METRIC_VERSION = "shared-unambiguous-alleles-v1"
INTERPRETATION = (
    "Single-linkage clusters describe allele similarity, not proven transmission. "
    "A chain can join isolates whose direct distance exceeds the threshold. "
    "Thresholds require organism-, scheme- and protocol-specific justification. "
    "AMR/accessory features are reported separately and never added to the core distance."
)
# Emitted with every snapshot comparison. Two forests side by side invite exactly
# these five misreadings, so the statements travel with the data, not with a view.
DIFF_CAVEATS = (
    "A minimum spanning forest is a layout of pairwise allele distances. It is not a phylogeny, "
    "not a time-ordered tree, and not a transmission chain.",
    "Two forests built from different cohorts can differ for cohort reasons alone: adding one "
    "isolate can re-route edges between isolates whose own pairwise distances did not change.",
    "Among equal distances the selected edge is deterministic but arbitrary. An edge present in "
    "one forest and absent from the other is not, on its own, evidence of change.",
    "A distance change accompanied by a changed shared-locus denominator is a measurement over a "
    "different locus set, not observed allele change.",
    "Cluster membership changes describe allele similarity under a local threshold. They are not "
    "proof that transmission started, stopped or occurred.",
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


def snapshot_graph(snapshot):
    """Replay the forest a stored snapshot displayed; no distance is recomputed.

    ``pairs`` are kept verbatim by :func:`build_snapshot` and
    ``forest_from_distances`` is deterministic (Kruskal ordered by distance,
    then source, then target), so the returned edges are the ones that were
    shown rather than a fresh measurement. ``known_alleles`` deliberately omits
    every uncallable locus, so ``len(result["alleles"])`` equals the profile's
    ``callable_loci`` and never its ``total_loci``: completeness must be read
    from those scalars, because a re-rendered snapshot read through the returned
    mapping would look fully called. Returned values are copies a renderer may
    annotate freely. What they draw is a layout of allele distances, not a
    phylogeny and not a transmission chain.
    """
    results = [dict(profile, alleles=dict(profile.get("known_alleles") or {}))
               for profile in snapshot.get("profiles", [])]
    edges = [dict(edge) for edge in forest_from_distances(snapshot.get("pairs", []))]
    return results, edges, deepcopy(list(snapshot.get("groups", [])))


def _snapshot_facts(snapshot):
    return {"snapshot_id": snapshot.get("snapshot_id", ""), "created_at": snapshot.get("created_at", ""),
            "investigation_id": snapshot.get("investigation_id"),
            "investigation_name": snapshot.get("investigation_name", ""),
            "scheme": snapshot.get("scheme", ""), "scheme_digest": snapshot.get("scheme_digest", ""),
            "threshold": snapshot.get("threshold"), "min_overlap": snapshot.get("min_overlap"),
            "metric_version": snapshot.get("metric_version", ""),
            "cohort_size": len(snapshot.get("profiles", []))}


def _diff_policy(baseline, current):
    """Decide whether two snapshots are on one measurement scale at all."""
    differences = [{"field": field, "baseline": baseline[field], "current": current[field]}
                   for field in ("scheme", "scheme_digest", "metric_version", "threshold", "min_overlap")
                   if baseline[field] != current[field]]
    blocking = sorted({field for field in ("scheme_digest", "metric_version")
                       if baseline[field] != current[field] or not baseline[field] or not current[field]})
    if blocking:
        reason = ("These snapshots were not measured on one scale (" + ", ".join(blocking) + " differs or is "
                  "unrecorded). Allele distances from different references or metrics are not comparable, so "
                  "no distance, edge or cluster comparison is offered here; only which isolates were added "
                  "or removed.")
    elif not differences:
        reason = ("Same reference, same link threshold and same minimum shared-locus fraction: the two "
                  "forests were built under one comparison policy.")
    else:
        reason = ("The comparison policy changed (" + "; ".join(
            f"{item['field']} {item['baseline']} → {item['current']}" for item in differences)
            + "). Differences between the two forests are partly a consequence of that change, not only "
              "of new evidence.")
    return {"comparable": not blocking, "identical": not blocking and not differences,
            "differences": differences,
            "distance_scale_attributable": baseline["min_overlap"] != current["min_overlap"],
            "grouping_attributable": baseline["threshold"] != current["threshold"], "reason": reason}


def _diff_isolates(baseline_profiles, current_profiles, names, *, comparable):
    """Cohort membership first: what is absent was never silently unchanged."""
    added = sorted(set(current_profiles) - set(baseline_profiles))
    removed = sorted(set(baseline_profiles) - set(current_profiles))
    retained = sorted(set(baseline_profiles) & set(current_profiles))
    evidence = None
    if comparable:
        evidence = [{"sample_id": sid, "sample_name": names[sid],
                     "baseline_callable": baseline_profiles[sid].get("callable_loci"),
                     "baseline_total": baseline_profiles[sid].get("total_loci"),
                     "current_callable": current_profiles[sid].get("callable_loci"),
                     "current_total": current_profiles[sid].get("total_loci")}
                    for sid in retained
                    if baseline_profiles[sid].get("signature") != current_profiles[sid].get("signature")]
    return {"added": [{"sample_id": sid, "sample_name": names[sid]} for sid in added],
            "removed": [{"sample_id": sid, "sample_name": names[sid]} for sid in removed],
            "retained": retained, "evidence_changed": evidence}


def _pair_key(row):
    return tuple(sorted((str(row["source"]), str(row["target"]))))


def _diff_pairs(baseline, current, retained, names):
    """Every row carries both denominators; a pair one snapshot never held is not 'unchanged'."""
    baseline_rows = {_pair_key(row): row for row in baseline.get("pairs", [])}
    current_rows = {_pair_key(row): row for row in current.get("pairs", [])}
    changed, became_comparable, became_excluded = [], [], []
    unchanged = denominator_only = still_excluded = not_assessed = 0
    for key in sorted(set(baseline_rows) | set(current_rows)):
        source, target = key
        before, after = baseline_rows.get(key), current_rows.get(key)
        if before is None or after is None or source not in retained or target not in retained:
            not_assessed += 1
            continue
        identity = {"source": source, "target": target,
                    "source_name": names.get(source, source), "target_name": names.get(target, target)}
        if before.get("comparable") and after.get("comparable"):
            denominator_changed = ((before.get("shared_loci"), before.get("total_loci"))
                                   != (after.get("shared_loci"), after.get("total_loci")))
            distance_changed = before.get("distance") != after.get("distance")
            if not denominator_changed and not distance_changed:
                unchanged += 1
                continue
            denominator_only += denominator_changed and not distance_changed
            changed.append({**identity, "baseline_distance": before.get("distance"),
                            "current_distance": after.get("distance"),
                            "baseline_shared": before.get("shared_loci"), "baseline_total": before.get("total_loci"),
                            "current_shared": after.get("shared_loci"), "current_total": after.get("total_loci"),
                            "distance_changed": distance_changed, "denominator_changed": denominator_changed})
        elif after.get("comparable"):
            became_comparable.append({**identity, "current_distance": after.get("distance"),
                                      "current_shared": after.get("shared_loci"),
                                      "current_total": after.get("total_loci"),
                                      "baseline_reason": before.get("reason", "")})
        elif before.get("comparable"):
            became_excluded.append({**identity, "baseline_distance": before.get("distance"),
                                    "baseline_shared": before.get("shared_loci"),
                                    "baseline_total": before.get("total_loci"),
                                    "current_reason": after.get("reason", "")})
        else:
            still_excluded += 1
    return {"distance_changed": changed, "became_comparable": became_comparable,
            "became_excluded": became_excluded, "unchanged_count": unchanged,
            "denominator_only_count": denominator_only, "still_excluded_count": still_excluded,
            "not_assessed_count": not_assessed}


def _diff_edges(baseline, current, names):
    """Presentation only: which lines were drawn, never which transmission happened."""
    baseline_edges = {_pair_key(edge): edge for edge in forest_from_distances(baseline.get("pairs", []))}
    current_edges = {_pair_key(edge): edge for edge in forest_from_distances(current.get("pairs", []))}

    def rows(edges, other):
        return [{"source": key[0], "target": key[1], "source_name": names.get(key[0], key[0]),
                 "target_name": names.get(key[1], key[1]), "distance": edge.get("distance")}
                for key, edge in sorted(edges.items()) if key not in other]

    return {"only_in_baseline": rows(baseline_edges, current_edges),
            "only_in_current": rows(current_edges, baseline_edges),
            "shared": len(set(baseline_edges) & set(current_edges)),
            "note": "Spanning-forest edges are a drawing choice. Among equal distances the edge is picked "
                    "deterministically but arbitrarily, and one added isolate can re-route edges between "
                    "isolates whose own distances did not change."}


def _diff_groups(baseline, current, retained, added_ids):
    """Merge and split are read from isolates present in both snapshots, nothing else."""
    baseline_groups, current_groups = list(baseline.get("groups", [])), list(current.get("groups", []))
    baseline_kept = {group["id"]: set(group["members"]) & retained for group in baseline_groups}
    current_kept = {group["id"]: set(group["members"]) & retained for group in current_groups}
    crosswalk, merged = [], []
    for group in current_groups:
        kept = current_kept[group["id"]]
        origins = [{"id": parent["id"], "name": parent.get("name", ""), "status": parent.get("status", ""),
                    "shared_members": sorted(baseline_kept[parent["id"]] & kept)}
                   for parent in baseline_groups if baseline_kept[parent["id"]] & kept]
        crosswalk.append({"current_group_id": group["id"], "current_group_name": group.get("name", ""),
                          "current_status": group.get("status", ""), "current_members": list(group["members"]),
                          "from_baseline_groups": origins,
                          "added_members": sorted(set(group["members"]) & added_ids),
                          "carried_members": sorted(kept)})
        if sum(origin["status"] == "cluster" for origin in origins) > 1:
            merged.append(group["id"])
    current_ids = {group["id"] for group in current_groups}
    dissolved, split = [], []
    for group in baseline_groups:
        kept = baseline_kept[group["id"]]
        went_to = [{"id": heir["id"], "name": heir.get("name", ""), "status": heir.get("status", ""),
                    "shared_members": sorted(current_kept[heir["id"]] & kept)}
                   for heir in current_groups if current_kept[heir["id"]] & kept]
        if group.get("status") == "cluster" and len(went_to) > 1:
            split.append(group["id"])
        if group["id"] not in current_ids:
            dissolved.append({"baseline_group_id": group["id"], "name": group.get("name", ""),
                              "status": group.get("status", ""), "members": list(group["members"]),
                              "went_to": went_to,
                              "lost_members": sorted(set(group["members"]) - retained),
                              "members_intact": bool(kept) and len(went_to) == 1
                              and current_kept[went_to[0]["id"]] == kept})
    return {"crosswalk": crosswalk, "dissolved": dissolved, "merged": merged, "split": split,
            "note": "A baseline group is listed as dissolved when no current group carries its identifier. "
                    "Read went_to before concluding anything: the same isolates can still be together "
                    "under a new cluster identity."}


def snapshot_diff(baseline, current):
    """Compare two stored snapshots, stating plainly what cannot be compared.

    Nothing is recomputed from sequence; both snapshots are read as the evidence
    they recorded. A different reference fingerprint or metric version makes the
    two sets of distances incommensurable, so ``pairs``, ``mst_edges`` and
    ``groups`` are then ``None`` rather than silently diffed and only cohort
    membership is reported. Every pair row carries both shared-locus
    denominators: an unchanged number measured over a different locus set is a
    different measurement rather than stability, and a changed number over a
    changed denominator is not observed allele change. Cluster merges and splits
    are derived only from isolates present in both snapshots, so an isolate that
    merely left the project never reads as a split; added and removed isolates
    are reported separately, and pairs touching them are counted as not assessed
    rather than unchanged. Threshold and minimum-overlap changes are flagged as
    policy-attributable, because they move groups and pair comparability by
    themselves. ``None`` anywhere in this structure means "not assessed", never
    "no change".
    """
    if not baseline or not current:
        raise ValueError("A snapshot comparison needs both a baseline and a current snapshot.")
    baseline_facts, current_facts = _snapshot_facts(baseline), _snapshot_facts(current)
    policy = _diff_policy(baseline_facts, current_facts)
    baseline_profiles = {str(p["sample_id"]): p for p in baseline.get("profiles", [])}
    current_profiles = {str(p["sample_id"]): p for p in current.get("profiles", [])}
    names = {sid: str(profile.get("sample_name", sid))
             for sid, profile in {**baseline_profiles, **current_profiles}.items()}
    isolates = _diff_isolates(baseline_profiles, current_profiles, names, comparable=policy["comparable"])
    retained = set(isolates["retained"])
    added_ids = {row["sample_id"] for row in isolates["added"]}
    caveats = list(DIFF_CAVEATS)
    pairs = edges = groups = None
    if policy["comparable"]:
        pairs = _diff_pairs(baseline, current, retained, names)
        edges = _diff_edges(baseline, current, names)
        groups = _diff_groups(baseline, current, retained, added_ids)
    else:
        caveats.append(policy["reason"])
    if policy["grouping_attributable"]:
        caveats.append(f"The link threshold changed from {baseline_facts['threshold']} to "
                       f"{current_facts['threshold']}. Group merges, splits and membership changes are at "
                       "least partly a consequence of that cutoff, not of new evidence.")
    if policy["distance_scale_attributable"]:
        caveats.append(f"The minimum shared-locus fraction changed from {baseline_facts['min_overlap']} to "
                       f"{current_facts['min_overlap']}. Pairs that became comparable or excluded are at "
                       "least partly a consequence of that policy, not of new evidence.")
    if isolates["added"] or isolates["removed"]:
        caveats.append(f"The cohort changed: {len(isolates['added'])} isolate(s) added and "
                       f"{len(isolates['removed'])} removed. Pairs touching them existed in only one "
                       "snapshot and were not assessed for change.")
    caveats.append(INTERPRETATION)
    summary = {"isolates_added": len(isolates["added"]), "isolates_removed": len(isolates["removed"]),
               "isolates_retained": len(retained),
               "evidence_changed": None if isolates["evidence_changed"] is None else len(isolates["evidence_changed"]),
               "distance_changed": None if pairs is None else len(pairs["distance_changed"]),
               "denominator_changed": None if pairs is None else sum(
                   row["denominator_changed"] for row in pairs["distance_changed"]),
               "became_comparable": None if pairs is None else len(pairs["became_comparable"]),
               "became_excluded": None if pairs is None else len(pairs["became_excluded"]),
               "not_assessed_pairs": None if pairs is None else pairs["not_assessed_count"],
               "edges_only_in_baseline": None if edges is None else len(edges["only_in_baseline"]),
               "edges_only_in_current": None if edges is None else len(edges["only_in_current"]),
               "groups_merged": None if groups is None else len(groups["merged"]),
               "groups_split": None if groups is None else len(groups["split"]),
               "groups_dissolved": None if groups is None else len(groups["dissolved"])}
    return {"format_version": 1, "baseline": baseline_facts, "current": current_facts, "policy": policy,
            "isolates": isolates, "pairs": pairs, "mst_edges": edges, "groups": groups,
            "summary": summary, "caveats": caveats}


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

    def baseline_snapshot_id(self, investigation_id):
        """The pinned baseline, else the first snapshot ever built for this plan.

        The baseline is a pointer into the append-only snapshot list, never a
        second copy: re-analysis, added isolates and reference updates append a
        new snapshot and cannot rewrite what the baseline recorded.
        """
        plan = self.get(investigation_id)
        pinned = plan.get("baseline_snapshot_id")
        if pinned and pinned in plan["snapshots"]:
            return pinned
        return plan["snapshots"][0] if plan["snapshots"] else None

    @_transactional
    def set_baseline(self, investigation_id, snapshot_id):
        """Repin the baseline pointer; stored snapshots stay immutable either way.

        ``updated_at`` is deliberately left alone: it is the optimistic-concurrency
        fence for in-flight comparisons, and choosing which stored snapshot to
        show beside the current one changes no comparison policy.
        """
        plan = self.get(investigation_id)
        if snapshot_id is not None and snapshot_id not in plan["snapshots"]:
            raise KeyError(snapshot_id)
        plan["baseline_snapshot_id"] = snapshot_id
        self.project.set_setting(f"investigation.{investigation_id}", plan)
        self.project.record_history(None, "investigation_baseline_pinned",
                                    {"id": investigation_id, "snapshot_id": snapshot_id})
        return snapshot_id

    def snapshot_catalog(self, investigation_id):
        """Picker labels read from history rows instead of the snapshot blobs.

        Every label is evidence the snapshot itself recorded (cohort, threshold,
        reference fingerprint), so two snapshots are never offered for comparison
        without the policy that produced them. Reads the whole history table:
        call it when a picker opens, never on every refresh.
        """
        plan = self.get(investigation_id)
        baseline = self.baseline_snapshot_id(investigation_id)
        rows = {entry["details"].get("snapshot_id"): entry for entry in self.project.history(None)
                if entry["action"] == "investigation_snapshot"
                and entry["details"].get("id") == investigation_id}
        catalog = []
        for number, sid in enumerate(plan["snapshots"], 1):
            details = rows.get(sid, {}).get("details", {})
            catalog.append({"snapshot_id": sid, "number": number,
                            "created_at": rows.get(sid, {}).get("created_at"),
                            "cohort_size": len(details.get("sample_ids", [])) or None,
                            "threshold": details.get("threshold"),
                            "scheme_digest": details.get("scheme_digest"),
                            "is_baseline": sid == baseline})
        return catalog

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
