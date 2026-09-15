"""Agglomerative hierarchical clustering over comparable pairwise rows.

This module builds a dendrogram from the same pairwise rows the rest of the
application already produces (:func:`wmlstudio.comparison.pairwise_distances`
for allele distances, :func:`wmlstudio.snp_tree.at_minimum_shared_fraction` for
SNP distances). It reads only ``source``, ``target``, ``comparable``,
``distance`` and ``reason``, so any row of that shape can be clustered, and it
holds no Qt objects: it is a computation, not a view.

**A pair that is not comparable is missing data, and missing data is not a
distance.** Two isolates that share too few loci, or too little sequence, have
no measured difference at all. Reading that absence as zero, or as "no
difference found", would invent an outbreak out of an unsequenced locus, so this
module never does it. Two policies are offered and the one in force is recorded
in the returned structure under ``missing_policy`` so a view can state it:

``refuse`` (the default)
    Two clusters may merge only when every pair of members across them was
    actually compared. A single incomparable pair blocks the merge for good, and
    the two clusters stay apart at every height. The tree is therefore a forest,
    ``complete`` is False, and ``blocked_merges`` names the clusters that were
    kept apart for want of evidence rather than for being far apart.

``provisional``
    The merge proceeds on the pairs that were compared, ignoring none of them
    and inventing none. Every affected merge, and every merge above it, carries
    ``provisional=True`` and the count of cross pairs that were never compared,
    so a view can mark the cluster as resting on incomplete evidence. A
    provisional average or complete-linkage height is computed over fewer pairs
    than it appears to span and is, in particular, a lower bound.

Heights are the distances of the rows themselves. A 7-locus MLST distance, a
cgMLST target distance and a SNP distance are different quantities and never
share a tree; ``scale_caption`` travels with the result so the axis of a
dendrogram can say which one is being drawn.

Nothing here decides that isolates are an outbreak. A dendrogram is a summary of
the distance rows and of nothing else.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from itertools import combinations
from typing import Any

from .sequence import check_cancelled

# Heights can be integer allele or SNP counts, or averages of them. The slack is
# only ever used to keep binary floating point noise from moving a cluster
# across a cut; it is far below any distance a user can enter.
TOLERANCE = 1e-9

LINKAGES = {
    "single": "Single linkage joins two clusters at the distance between their closest pair, so "
              "one near pair joins two groups however far apart the rest of their members are.",
    "complete": "Complete linkage joins two clusters at the distance between their furthest pair, "
                "so every member of the joined cluster lies within that height of every other.",
    "average": "Average linkage (UPGMA) joins two clusters at the mean distance over the pairs "
               "across them, which is neither their closest nor their furthest pair.",
}

MISSING_POLICIES = {
    "refuse": "Clusters were joined only where every pair across them had an accepted distance. "
              "Pairs that were not comparable blocked the join instead of counting as a "
              "difference of zero, so some groups are apart for lack of evidence, not for being "
              "far apart.",
    "provisional": "Clusters were joined on the pairs that had an accepted distance, and joins "
                   "that spanned pairs which were never compared are marked provisional. The "
                   "height of such a join was measured over fewer pairs than it spans.",
}

MISSING_DATA_RULE = (
    "A pair that was not comparable carries no distance. It is missing data: never zero, never "
    "'no difference found', and never evidence that two isolates are alike."
)

CHAINING_RULE = (
    "A cluster can be held together by a chain of short links while its two furthest members are "
    "much further apart than the height at which they were joined. Where that happened, the "
    "cluster reports the distance between its furthest compared pair beside its linkage height."
)

# Rows describing different quantities must not be clustered together. These are
# the keys the application's row builders use to name what was measured.
_SCALE_KEYS = ("typing_kind", "metric", "metric_version")


def _label(row: Mapping[str, Any], key: str) -> str:
    value = row.get(key)
    if value is None or not str(value).strip():
        raise ValueError(f"Every pairwise row needs a {key} identifier.")
    return str(value)


def _check_distance(value: Any, source: str, target: str) -> float | int:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Comparable pair {source}/{target} needs a finite numeric distance.")
    if value < 0:
        raise ValueError(f"Comparable pair {source}/{target} has a negative distance.")
    return value


def _read_rows(rows: Iterable[Mapping[str, Any]], labels: Sequence[str], cancelled):
    """Split the rows into accepted distances and explicitly uncompared pairs."""
    distances: dict[tuple[str, str], float | int] = {}
    uncompared: dict[tuple[str, str], str] = {}
    names: set[str] = {str(name) for name in labels}
    scales: dict[str, set[str]] = {key: set() for key in _SCALE_KEYS}
    seen: set[tuple[str, str]] = set()
    for index, row in enumerate(rows):
        if index % 256 == 0:
            check_cancelled(cancelled)
        source, target = _label(row, "source"), _label(row, "target")
        if source == target:
            raise ValueError(f"A pairwise row compares {source} with itself.")
        key = (source, target) if source <= target else (target, source)
        if key in seen:
            raise ValueError(f"Duplicate pairwise row for {key[0]} and {key[1]}.")
        seen.add(key)
        names.update(key)
        for scale_key in _SCALE_KEYS:
            value = row.get(scale_key)
            if value not in (None, ""):
                scales[scale_key].add(str(value))
        if row.get("comparable"):
            distances[key] = _check_distance(row.get("distance"), *key)
        else:
            uncompared[key] = str(row.get("reason") or "The pair was not comparable.")
    for scale_key, values in scales.items():
        if len(values) > 1:
            raise ValueError(
                f"These rows mix different quantities ({scale_key}: {', '.join(sorted(values))}). "
                "Distances of different kinds never share a tree.")
    # A pair nobody produced a row for was not compared either. Saying so here
    # keeps the clustering from silently treating an absent row as agreement.
    ordered = sorted(names)
    for pair in combinations(ordered, 2):
        if pair not in seen:
            uncompared[pair] = "No pairwise comparison was recorded for this pair."
    return ordered, distances, uncompared


def _cross_state(comparable: bool, distance: float | int | None) -> list:
    """[compared pairs, uncompared pairs, minimum, maximum, sum of distances]."""
    if comparable:
        return [1, 0, distance, distance, distance]
    return [0, 1, None, None, 0]


def _combine(left: list, right: list) -> list:
    minimums = [value for value in (left[2], right[2]) if value is not None]
    maximums = [value for value in (left[3], right[3]) if value is not None]
    return [left[0] + right[0], left[1] + right[1],
            min(minimums) if minimums else None,
            max(maximums) if maximums else None,
            left[4] + right[4]]


def _height(state: list, linkage: str) -> float | int:
    if linkage == "single":
        return state[2]
    if linkage == "complete":
        return state[3]
    return state[4] / state[0]


def build_tree(
    rows: Iterable[Mapping[str, Any]], *, linkage: str = "average", labels: Sequence[str] = (),
    missing: str = "refuse", scale_caption: str = "", cancelled=None,
) -> dict[str, Any]:
    """Cluster the rows agglomeratively and return the merge sequence with heights.

    ``labels`` adds isolates that appear in no row at all; every endpoint of a
    row is included automatically. ``missing`` selects the policy for pairs that
    were not comparable, which is the whole difficulty here and is described in
    the module docstring. The result is deterministic: it depends on the rows and
    on nothing else, ties included.
    """
    if linkage not in LINKAGES:
        raise ValueError(f"Linkage must be one of {', '.join(sorted(LINKAGES))}.")
    if missing not in MISSING_POLICIES:
        raise ValueError(f"The policy for uncompared pairs must be one of "
                         f"{', '.join(sorted(MISSING_POLICIES))}.")
    names, distances, uncompared = _read_rows(rows, labels, cancelled)

    clusters: dict[str, dict[str, Any]] = {}
    for name in names:
        clusters["leaf:" + name] = {
            "id": "leaf:" + name, "members": (name,), "key": (name,), "size": 1,
            "height": 0, "diameter": None, "uncompared_within": 0, "provisional": False,
        }

    def state_key(left: str, right: str) -> tuple[str, str]:
        return (left, right) if clusters[left]["key"] <= clusters[right]["key"] else (right, left)

    state: dict[tuple[str, str], list] = {}
    for one, two in combinations(names, 2):
        pair = (one, two)
        state[state_key("leaf:" + one, "leaf:" + two)] = _cross_state(
            pair in distances, distances.get(pair))

    merges: list[dict[str, Any]] = []
    while True:
        check_cancelled(cancelled)
        best = None
        for pair, entry in state.items():
            if not entry[0] or (entry[1] and missing == "refuse"):
                continue
            height = _height(entry, linkage)
            candidate = (height, clusters[pair[0]]["key"], clusters[pair[1]]["key"])
            if best is None or candidate < best[0]:
                best = (candidate, pair, entry, height)
        if best is None:
            break
        _, pair, entry, height = best
        left, right = clusters[pair[0]], clusters[pair[1]]
        compared, missed, _, furthest, _ = entry
        members = tuple(sorted(left["members"] + right["members"]))
        diameters = [value for value in (left["diameter"], right["diameter"], furthest)
                     if value is not None]
        diameter = max(diameters) if diameters else None
        child_height = max(left["height"], right["height"])
        node = {
            "id": f"node:{len(merges) + 1}", "members": members, "key": members,
            "size": len(members), "height": height, "diameter": diameter,
            "uncompared_within": left["uncompared_within"] + right["uncompared_within"] + missed,
            "provisional": bool(missed) or left["provisional"] or right["provisional"],
        }
        chained = diameter is not None and diameter > height + TOLERANCE
        merges.append({
            "step": len(merges), "id": node["id"], "height": height,
            "left": left["id"], "right": right["id"],
            "left_members": list(left["members"]), "right_members": list(right["members"]),
            "members": list(members), "size": node["size"],
            "compared_pairs_across": compared, "uncompared_pairs_across": missed,
            "provisional": node["provisional"],
            "max_within_distance": diameter,
            "uncompared_within_pairs": node["uncompared_within"],
            "chained": chained,
            "chain_note": (
                f"Joined at {height:g} by its closest compared pair, but its furthest compared "
                f"pair differs by {diameter:g}." if chained else ""),
            # A height below the height of a child would draw as an inverted
            # branch. It cannot happen under the refuse policy with these three
            # linkages; it can once a provisional merge measures fewer pairs.
            "inversion": height < child_height - TOLERANCE,
        })
        for other in list(clusters):
            if other in (left["id"], right["id"]):
                continue
            combined = _combine(state.pop(state_key(left["id"], other)),
                                state.pop(state_key(right["id"], other)))
            state[(other, node["id"])] = combined
        state.pop(pair)
        del clusters[left["id"]], clusters[right["id"]]
        clusters[node["id"]] = node
        # The new cluster's key decides the ordering of its state entries, so the
        # entries parked under the old ordering are rewritten once it exists.
        for other in list(clusters):
            if other == node["id"]:
                continue
            state[state_key(node["id"], other)] = state.pop((other, node["id"]))

    return _finish(clusters, merges, names, distances, uncompared, state,
                   linkage, missing, scale_caption)


def _blocked_merges(clusters, state, uncompared, missing):
    """Name the cluster pairs that could not be joined, and why."""
    blocked = []
    for pair, entry in sorted(state.items()):
        if entry[0] and not (entry[1] and missing == "refuse"):
            continue
        left, right = clusters[pair[0]], clusters[pair[1]]
        reasons = sorted({uncompared[(one, two) if one <= two else (two, one)]
                          for one in left["members"] for two in right["members"]
                          if ((one, two) if one <= two else (two, one)) in uncompared})
        blocked.append({
            "left": left["id"], "right": right["id"],
            "left_members": list(left["members"]), "right_members": list(right["members"]),
            "compared_pairs": entry[0], "uncompared_pairs": entry[1],
            "reasons": reasons,
            "note": "These were kept apart because pairs across them were never compared, "
                    "not because they were measured as far apart.",
        })
    return blocked


def _leaf_order(merges, roots):
    """Left subtree before right, roots in their own deterministic order."""
    children = {merge["id"]: (merge["left"], merge["right"]) for merge in merges}
    order = []
    for root in roots:
        pending = [root["id"]]
        while pending:
            node = pending.pop()
            if node in children:
                pending.extend(reversed(children[node]))
            else:
                order.append(node.removeprefix("leaf:"))
    return order


def _finish(clusters, merges, names, distances, uncompared, state,
            linkage, missing, scale_caption):
    roots = sorted(clusters.values(), key=lambda cluster: cluster["key"])
    compared_with_someone = {name for pair in distances for name in pair}
    blocked = _blocked_merges(clusters, state, uncompared, missing)
    chained = [merge["id"] for merge in merges if merge["chained"]]
    return {
        "format_version": 1,
        "linkage": linkage,
        "linkage_caption": LINKAGES[linkage],
        "missing_policy": missing,
        "missing_policy_note": MISSING_POLICIES[missing],
        "missing_data_rule": MISSING_DATA_RULE,
        "chaining_rule": CHAINING_RULE,
        "scale_caption": scale_caption,
        "labels": list(names),
        "leaf_order": _leaf_order(merges, roots),
        "merges": merges,
        "roots": [{"id": root["id"], "members": list(root["members"]), "size": root["size"],
                   "height": root["height"], "max_within_distance": root["diameter"],
                   "provisional": root["provisional"],
                   "uncompared_within_pairs": root["uncompared_within"]} for root in roots],
        "complete": len(roots) <= 1,
        "blocked_merges": blocked,
        "compared_pairs": len(distances),
        "uncompared_pairs": len(uncompared),
        "isolates_without_a_comparable_pair": sorted(set(names) - compared_with_someone),
        "chained_merges": chained,
        "inversions": [merge["id"] for merge in merges if merge["inversion"]],
    }


def merge_heights(tree: Mapping[str, Any]) -> list[float | int]:
    """The distinct heights at which this tree can be cut, lowest first."""
    return sorted({merge["height"] for merge in tree.get("merges", [])})


def cut_tree(tree: Mapping[str, Any], height: float | int) -> dict[str, Any]:
    """Cut the dendrogram at ``height`` and name the groups deterministically.

    A merge is applied only when its own height and the heights of everything
    below it are at or under the cut, which is what cutting a dendrogram with a
    horizontal line means and which stays correct if a provisional merge ever
    drew an inverted branch. Isolates kept apart by the policy for uncompared
    pairs stay apart here at every height; the result says so rather than letting
    the cut read as a measured separation.
    """
    if isinstance(height, bool) or not isinstance(height, (int, float)) \
            or not math.isfinite(height) or height < 0:
        raise ValueError("A cut height must be a finite nonnegative number.")
    below: dict[str, bool] = {}
    applied = []
    for merge in tree.get("merges", []):
        inside = (merge["height"] <= height + TOLERANCE
                  and below.get(merge["left"], True) and below.get(merge["right"], True))
        below[merge["id"]] = inside
        if inside:
            applied.append(merge)
    parent = {name: name for name in tree.get("labels", [])}

    def root(name: str) -> str:
        while name != parent[name]:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for merge in applied:
        anchor = root(merge["members"][0])
        for member in merge["members"]:
            parent[root(member)] = anchor
    components: dict[str, list[str]] = {}
    for name in tree.get("labels", []):
        components.setdefault(root(name), []).append(name)
    formed = {frozenset(merge["members"]): merge for merge in applied}
    links: dict[frozenset, list] = {}
    for merge in applied:
        links.setdefault(frozenset(components[root(merge["members"][0])]), []).append(merge)
    lonely = set(tree.get("isolates_without_a_comparable_pair", []))

    groups = []
    for members in sorted(components.values(), key=lambda group: (-len(group), group[0])):
        members = sorted(members)
        top = formed.get(frozenset(members))
        inner = sorted(merge["height"] for merge in links.get(frozenset(members), []))
        diameter = top["max_within_distance"] if top else None
        chained = diameter is not None and diameter > height + TOLERANCE
        status = "cluster" if len(members) > 1 else \
            "not_comparable" if members[0] in lonely else "singleton"
        groups.append({
            "id": "hc:" + members[0], "anchor": members[0], "number": None,
            "name": {"cluster": "", "singleton": "Unlinked isolate",
                     "not_comparable": "Not comparable"}[status],
            "status": status, "members": members, "size": len(members),
            "height": top["height"] if top else None,
            "link_heights": inner,
            "max_within_distance": diameter,
            "chained": chained,
            "chain_note": (
                f"Held together by a chain of links: the group was cut at {height:g}, but its "
                f"furthest compared pair differs by {diameter:g}." if chained else ""),
            "provisional": bool(top and top["provisional"]),
            "uncompared_within_pairs": top["uncompared_within_pairs"] if top else 0,
            "reason": "No pairwise comparison of this isolate was accepted; review missing or "
                      "mixed evidence." if status == "not_comparable" else "",
        })
    number = 0
    for group in groups:
        if group["status"] == "cluster":
            number += 1
            group["number"] = number
            group["name"] = f"HC group {number:03d}"

    blocked = tree.get("blocked_merges", [])
    return {
        "format_version": 1,
        "height": height,
        "threshold_applied": height,
        "linkage": tree.get("linkage", ""),
        "linkage_caption": tree.get("linkage_caption", ""),
        "missing_policy": tree.get("missing_policy", ""),
        "missing_policy_note": tree.get("missing_policy_note", ""),
        "missing_data_rule": MISSING_DATA_RULE,
        "scale_caption": tree.get("scale_caption", ""),
        "groups": groups,
        "clusters": sum(group["status"] == "cluster" for group in groups),
        "clustered_isolates": sum(group["size"] for group in groups
                                  if group["status"] == "cluster"),
        "chained_groups": [group["id"] for group in groups if group["chained"]],
        "provisional_groups": [group["id"] for group in groups if group["provisional"]],
        "groups_kept_apart_by_uncompared_pairs": len(blocked),
        "note": "This cut applied a height of "
                f"{height:g} to a {tree.get('linkage', '')}-linkage tree."
                + (" Some groups are separate because pairs across them were never compared, not "
                   "because they were measured as far apart." if blocked else ""),
    }


def dendrogram_layout(tree: Mapping[str, Any]) -> dict[str, Any]:
    """Coordinates for drawing the tree, with no drawing toolkit involved.

    Leaves sit at ``x = 0, 1, 2, ...`` in ``leaf_order`` and at ``y = 0``; a node
    sits at the mean ``x`` of its children and at ``y`` equal to its height. A
    view scales these to pixels; the arithmetic that has to be right belongs
    here, where it can be tested without a screen.
    """
    positions = {"leaf:" + name: (float(index), 0.0)
                 for index, name in enumerate(tree.get("leaf_order", []))}
    nodes = []
    for merge in tree.get("merges", []):
        left, right = positions[merge["left"]], positions[merge["right"]]
        point = ((left[0] + right[0]) / 2, float(merge["height"]))
        positions[merge["id"]] = point
        nodes.append({
            "id": merge["id"], "x": point[0], "y": point[1],
            "left": merge["left"], "left_x": left[0], "left_y": left[1],
            "right": merge["right"], "right_x": right[0], "right_y": right[1],
            "provisional": merge["provisional"], "chained": merge["chained"],
        })
    return {
        "leaves": [{"label": name, "x": float(index), "y": 0.0}
                   for index, name in enumerate(tree.get("leaf_order", []))],
        "nodes": nodes,
        "roots": [root["id"] for root in tree.get("roots", [])],
        "max_height": max((node["y"] for node in nodes), default=0.0),
        "height_caption": tree.get("scale_caption", ""),
        "linkage_caption": tree.get("linkage_caption", ""),
    }
