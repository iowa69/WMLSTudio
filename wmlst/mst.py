# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Minimum spanning tree over typed isolates -- the maths only, no drawing.

An MLST run answers "what is this isolate?" one file at a time. The question a
reader asks next is "how do these isolates relate to each other?", and the
standard answer is a **minimum spanning tree** (MST) over the allelic distance:
every isolate is a node, every edge is labelled with the number of loci at which
the two profiles disagree, and the tree keeps the cheapest set of edges that
still joins everything up. Cut the edges longer than one locus and what remains
are the classic single-locus-variant clonal complexes.

This module computes that tree and nothing else. It has **no GUI code, no
tkinter import and no I/O**: it takes plain sample records in and returns plain
dataclasses out, so it can be unit tested headless and reused by a report.

Three rules shape the whole file:

1. **Missing data never counts as a difference.** A locus is compared only when
   BOTH isolates called it; ``-``, ``0`` and an empty cell are absent data, not
   allele zero. Two profiles therefore carry both a distance and the number of
   loci that distance was measured over, and the second number travels with the
   first (:class:`AllelicDistance` is an ``int`` that remembers it) so a caller
   can show "2 loci differ, of 7 compared" and a reader can discount a distance
   measured over three loci.
2. **Only isolates typed with the SAME scheme may be joined.** An
   ``saureus`` profile and a ``kpneumoniae`` profile share neither loci nor
   allele numbering; comparing them position-by-position would be nonsense. The
   tree is therefore a forest: one component per scheme, plus a singleton for
   every isolate that has no scheme, no calls, or no comparable partner.
3. **The same input always gives the same tree, on every machine.** Prim's
   algorithm has ties whenever two edges weigh the same -- which, on allelic
   distances that are small integers, is most of the time. Every tie is broken
   by the (distance, label pair) key below, every iteration order is a sorted
   one, and the layout is a seeded radial construction: there is no
   ``random`` import in this module, and there must never be one.

Public API::

    allelic_distance(a, b)          -> AllelicDistance (an int, .n_compared)
    distance_matrix(samples)        -> tuple of rows of AllelicDistance
    build_mst(samples)              -> Mst
    clonal_groups(mst, threshold=1) -> tuple of ClonalGroup, biggest first
    layout(mst, width, height)      -> tuple of MstNode with x/y filled in
    MstSample / MstNode / MstEdge / Mst / ClonalGroup   dataclasses
    sample_from_result(r)           adapt one engine.SampleResult
    samples_from_results(rs)        adapt a whole run
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "MISSING_CODES",
    "AllelicDistance",
    "ClonalGroup",
    "Mst",
    "MstEdge",
    "MstNode",
    "MstSample",
    "allelic_distance",
    "build_mst",
    "clonal_groups",
    "codes_match",
    "distance_matrix",
    "edge_length",
    "layout",
    "normalise_code",
    "sample_from_result",
    "samples_from_results",
]

# ---------------------------------------------------------------------------
# 1. Allele codes
# ---------------------------------------------------------------------------

#: Codes that mean "this locus was not called", case-insensitively. ``0`` is
#: upstream's null allele and ``-`` its missing one; neither is an allele
#: number, so neither may ever make two isolates look different.
MISSING_CODES = frozenset({"", "-", "0", "?", "n/a", "na", "none", "null"})


def normalise_code(code: Any) -> str:
    """Reduce one allele code to a comparable form, or ``""`` when missing.

    Upstream decorates allele numbers: ``~16`` is an inexact (novel) match to
    allele 16, ``16?`` a partial one, ``1,2`` two equally good alleles at one
    locus. The decorations describe how confidently the allele was *seen*, not
    which allele it is, so they are stripped; a multiple call is reduced to its
    sorted set of numbers, kept as a comma-joined string.
    """
    if code is None:
        return ""
    text = str(code).strip()
    if not text:
        return ""
    parts: List[str] = []
    for raw in text.split(","):
        item = raw.strip().lstrip("~").rstrip("?").strip()
        if not item or item.lower() in MISSING_CODES:
            continue
        parts.append(item)
    if not parts:
        return ""
    return ",".join(sorted(set(parts), key=_code_sort_key))


def _code_sort_key(code: str) -> Tuple[int, float, str]:
    """Numeric-first ordering so ``2`` sorts before ``10`` and before ``2a``."""
    try:
        return (0, float(code), code)
    except ValueError:
        return (1, 0.0, code)


def codes_match(a: Any, b: Any) -> Optional[bool]:
    """``True`` same allele, ``False`` different, ``None`` not comparable.

    ``None`` -- one side or both not called -- is the case that must not be
    folded into ``False``: it is the reason :class:`AllelicDistance` reports
    how many loci it could actually look at.

    A multiple call such as ``1,2`` matches any profile that carries one of
    those alleles: the isolates are not demonstrably different at that locus.
    """
    ca = normalise_code(a)
    cb = normalise_code(b)
    if not ca or not cb:
        return None
    if ca == cb:
        return True
    return bool(set(ca.split(",")) & set(cb.split(",")))


class AllelicDistance(int):
    """The number of differing loci, carrying the size of the comparison.

    It IS an ``int`` -- ``d == 2``, ``sorted()``, ``sum()`` and ``"%d" % d`` all
    behave -- with two extra read-only numbers attached:

    ``n_compared``
        loci where both isolates made a call. ``0`` means the profiles share no
        called locus at all: the distance is then meaningless, ``comparable``
        is ``False``, and :func:`build_mst` refuses to draw such an edge.
    ``n_loci``
        loci lined up at all, i.e. the length of the shorter profile. The gap
        between the two numbers is how much data was missing.
    """

    def __new__(cls, distance: int, n_compared: int = 0, n_loci: int = 0):
        self = super().__new__(cls, distance)
        self.n_compared = int(n_compared)
        self.n_loci = int(n_loci)
        return self

    @property
    def comparable(self) -> bool:
        """Were there any loci both isolates called?"""
        return self.n_compared > 0

    @property
    def n_missing(self) -> int:
        """Loci skipped because at least one isolate had no call there."""
        return max(0, self.n_loci - self.n_compared)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return "AllelicDistance(%d, n_compared=%d, n_loci=%d)" % (
            int(self), self.n_compared, self.n_loci)


def allelic_distance(a: Sequence[Any], b: Sequence[Any]) -> AllelicDistance:
    """Pairwise loci that differ between two allele-code tuples.

    Profiles are compared position by position -- the engine emits allele codes
    in profile-header gene order, so position *is* locus -- over the length of
    the shorter one. A locus counts as different only when both isolates called
    it and the calls disagree; anything missing on either side is skipped and
    shows up in ``.n_compared`` instead.
    """
    codes_a = tuple(a or ())
    codes_b = tuple(b or ())
    n_loci = min(len(codes_a), len(codes_b))
    n_compared = 0
    distance = 0
    for i in range(n_loci):
        same = codes_match(codes_a[i], codes_b[i])
        if same is None:
            continue
        n_compared += 1
        if not same:
            distance += 1
    return AllelicDistance(distance, n_compared, n_loci)


# ---------------------------------------------------------------------------
# 2. Records
# ---------------------------------------------------------------------------

#: Scheme values that mean "no scheme was called". Such isolates are never
#: joined to anything: without a scheme there is no locus order to compare.
_NO_SCHEME = frozenset({"", "-", "none"})


@dataclass(frozen=True)
class MstSample:
    """One typed isolate, as the tree wants it.

    ``alleles`` is the profile in the scheme's gene order (``engine.
    SampleResult.alleles`` codes); ``st`` is carried for display only and is
    never used to compute a distance -- two isolates with the same ST are at
    distance 0 because their profiles agree, which is a property the tests
    check rather than assume.
    """

    label: str
    scheme: str = "-"
    st: str = "-"
    alleles: Tuple[str, ...] = ()

    @property
    def has_scheme(self) -> bool:
        return self.scheme.strip().lower() not in _NO_SCHEME

    @property
    def n_called(self) -> int:
        """Loci with a real allele code."""
        return sum(1 for code in self.alleles if normalise_code(code))


@dataclass
class MstNode:
    """One isolate placed in the tree. ``x``/``y`` are filled in by :func:`layout`."""

    label: str
    st: str = "-"
    scheme: str = "-"
    alleles: Tuple[str, ...] = ()
    group: int = 0          # clonal group id at the build threshold; 1-based
    x: float = 0.0
    y: float = 0.0


@dataclass(frozen=True)
class MstEdge:
    """A kept edge: two labels and the loci that differ between them."""

    a: str
    b: str
    distance: int

    @property
    def key(self) -> Tuple[str, str]:
        """The unordered label pair, ordered -- the deterministic tie-break."""
        return (self.a, self.b) if self.a <= self.b else (self.b, self.a)


@dataclass(frozen=True)
class ClonalGroup:
    """A connected component once long edges are cut. Ids are 1-based."""

    id: int
    labels: Tuple[str, ...]
    size: int
    scheme: str = "-"
    threshold: int = 1


@dataclass
class Mst:
    """The whole forest.

    ``components`` lists the labels of each connected component, biggest first,
    singletons included, so a caller never has to recompute connectivity.
    ``scheme`` is the one scheme every typed isolate shares, or ``""`` when the
    input mixes schemes. ``n_loci_compared`` is the SMALLEST number of loci any
    kept edge was measured over -- the honest headline number for "distances
    over N loci"; with no edges it is the smallest number of called loci.
    """

    nodes: Tuple[MstNode, ...] = ()
    edges: Tuple[MstEdge, ...] = ()
    components: Tuple[Tuple[str, ...], ...] = ()
    scheme: str = ""
    n_loci_compared: int = 0
    threshold: int = 1
    groups: Tuple[ClonalGroup, ...] = field(default=())

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def n_components(self) -> int:
        return len(self.components)

    def node(self, label: str) -> Optional[MstNode]:
        """The node carrying ``label``, or ``None``."""
        for item in self.nodes:
            if item.label == label:
                return item
        return None

    def index_of(self, label: str) -> int:
        """Position of ``label`` in :attr:`nodes`, or ``-1``."""
        for i, item in enumerate(self.nodes):
            if item.label == label:
                return i
        return -1

    def neighbours(self, label: str) -> Tuple[str, ...]:
        """Labels joined to ``label`` by a kept edge, in edge order."""
        out: List[str] = []
        for edge in self.edges:
            if edge.a == label:
                out.append(edge.b)
            elif edge.b == label:
                out.append(edge.a)
        return tuple(out)

    @property
    def max_distance(self) -> int:
        return max((edge.distance for edge in self.edges), default=0)


# ---------------------------------------------------------------------------
# 3. Adapting engine results
# ---------------------------------------------------------------------------

def sample_from_result(result: Any) -> MstSample:
    """Adapt one :class:`engine.SampleResult` (or a dict, or an MstSample).

    Anything with ``label``/``scheme``/``st``/``alleles`` is accepted, and an
    ``alleles`` tuple of :class:`engine.AlleleCall` objects is reduced to its
    codes -- so a caller can hand a run's samples straight over.
    """
    if isinstance(result, MstSample):
        return result
    if isinstance(result, dict):
        get = result.get
    else:
        def get(name, default=None):
            return getattr(result, name, default)
    label = str(get("label", "") or "")
    scheme = str(get("scheme", "-") or "-")
    st = str(get("st", "-") or "-")
    codes: List[str] = []
    for allele in get("alleles", ()) or ():
        code = getattr(allele, "code", None)
        if code is None and isinstance(allele, dict):
            code = allele.get("code")
        if code is None:
            code = allele
        codes.append(str(code))
    return MstSample(label=label, scheme=scheme, st=st, alleles=tuple(codes))


def samples_from_results(results: Iterable[Any]) -> Tuple[MstSample, ...]:
    """Adapt a whole run, skipping nothing: failures become singletons."""
    return tuple(sample_from_result(item) for item in results)


def _unique_labels(samples: Sequence[MstSample]) -> Tuple[MstSample, ...]:
    """Make labels unique, because labels are how edges name their ends.

    Two input files can share a basename, and the engine's label is that
    basename. The second and later collisions get a ``#2`` suffix rather than
    an exception: an MST of 40 genomes must not fail because two came from
    different directories.
    """
    seen: Dict[str, int] = {}
    out: List[MstSample] = []
    for sample in samples:
        label = sample.label or "(unnamed)"
        count = seen.get(label, 0) + 1
        seen[label] = count
        if count > 1:
            label = "%s #%d" % (label, count)
        out.append(sample if label == sample.label else replace(sample, label=label))
    return tuple(out)


# ---------------------------------------------------------------------------
# 4. The distance graph
# ---------------------------------------------------------------------------

def distance_matrix(samples: Sequence[Any]) -> Tuple[Tuple[AllelicDistance, ...], ...]:
    """Full symmetric matrix of :func:`allelic_distance`, diagonal zeroed.

    Computed once and reused by Prim's algorithm: the profiles are short but
    the graph is complete, so recomputing a distance inside the inner loop
    would cube the work for no reason.
    """
    items = [sample_from_result(item) for item in samples]
    n = len(items)
    rows: List[List[AllelicDistance]] = [
        [AllelicDistance(0, 0, 0)] * n for _ in range(n)
    ]
    for i in range(n):
        rows[i][i] = AllelicDistance(0, items[i].n_called, len(items[i].alleles))
        for j in range(i + 1, n):
            d = allelic_distance(items[i].alleles, items[j].alleles)
            rows[i][j] = d
            rows[j][i] = d
    return tuple(tuple(row) for row in rows)


def _scheme_key(sample: MstSample) -> str:
    """The grouping key: the scheme name, or ``""`` for "unknown, stand alone"."""
    return sample.scheme.strip() if sample.has_scheme else ""


def _prim(samples: Sequence[MstSample], idxs: Sequence[int],
          matrix: Sequence[Sequence[AllelicDistance]],
          min_shared: int) -> Tuple[List[Tuple[int, int, AllelicDistance]], List[List[int]]]:
    """Prim's algorithm over one scheme group, restarted per component.

    Returns ``(edges, components)`` with edges as ``(parent, child, distance)``
    index triples in the order the tree grew.

    Determinism, which is the whole point of the bookkeeping below: candidate
    edges are compared by ``(distance, label pair, parent index)``, never by
    ``distance`` alone and never by dictionary order, so two runs -- and two
    machines -- grow the identical tree. The restart loop exists because a
    group is not guaranteed connected: a profile that shares no called locus
    with anyone (``n_compared < min_shared``) has no finite edge and becomes
    its own component rather than being attached at a fabricated distance.
    """
    def node_key(i: int) -> Tuple[str, int]:
        return (samples[i].label, i)

    def edge_key(i: int, j: int) -> Tuple[Tuple[str, int], Tuple[str, int]]:
        ki, kj = node_key(i), node_key(j)
        return (ki, kj) if ki <= kj else (kj, ki)

    def grow(root: int, unvisited: List[int],
             edges: List[Tuple[int, int, AllelicDistance]]) -> List[int]:
        """Grow one component out from ``root``, consuming ``unvisited``."""
        component = [root]
        # best[j] = (distance, edge key, parent) -- the cheapest known edge
        # from the growing component to the outside node j.
        best: Dict[int, Tuple[int, Any, int]] = {}

        def relax(i: int) -> None:
            for j in unvisited:
                d = matrix[i][j]
                if d.n_compared < min_shared:
                    continue
                candidate = (int(d), edge_key(i, j), i)
                current = best.get(j)
                if current is None or candidate < current:
                    best[j] = candidate

        relax(root)
        while True:
            pick = None
            pick_key = None
            for j in unvisited:
                candidate = best.get(j)
                if candidate is None:
                    continue
                if pick_key is None or candidate < pick_key:
                    pick, pick_key = j, candidate
            if pick is None:
                return component
            parent = pick_key[2]
            edges.append((parent, pick, matrix[parent][pick]))
            component.append(pick)
            unvisited.remove(pick)
            best.pop(pick, None)
            relax(pick)

    unvisited = sorted(idxs, key=node_key)
    edges: List[Tuple[int, int, AllelicDistance]] = []
    components: List[List[int]] = []
    while unvisited:
        components.append(grow(unvisited.pop(0), unvisited, edges))
    return edges, components


def build_mst(samples: Iterable[Any], *, threshold: int = 1,
              min_shared: int = 1) -> Mst:
    """Minimum spanning forest over typed isolates (Prim's algorithm).

    ``samples`` are :class:`MstSample` records, ``engine.SampleResult`` objects
    or dicts -- anything :func:`sample_from_result` understands. Isolates are
    grouped by scheme and each group gets its own spanning tree; isolates with
    no scheme, with no called locus, or with no comparable partner come back as
    singleton components with no edges.

    ``min_shared`` is the number of loci two profiles must both have called
    before an edge between them is allowed to exist at all (1 = any overlap).
    ``threshold`` is the clonal-group cut recorded in :attr:`Mst.groups` and in
    each node's ``group``; it does not change the tree.

    Degenerate inputs are answers, not errors: no samples gives an empty tree,
    one sample gives one node and no edges, identical profiles give distance-0
    edges (they are the same clone and must stay joined), and a mix of schemes
    gives one component per scheme.
    """
    items = _unique_labels([sample_from_result(item) for item in samples])
    n = len(items)
    nodes = tuple(
        MstNode(label=s.label, st=s.st or "-", scheme=s.scheme or "-",
                alleles=tuple(s.alleles))
        for s in items
    )
    if n == 0:
        return Mst(nodes=(), edges=(), components=(), scheme="",
                   n_loci_compared=0, threshold=threshold, groups=())

    matrix = distance_matrix(items)

    groups_by_scheme: Dict[str, List[int]] = {}
    for i, sample in enumerate(items):
        groups_by_scheme.setdefault(_scheme_key(sample), []).append(i)

    raw_edges: List[Tuple[int, int, AllelicDistance]] = []
    raw_components: List[List[int]] = []
    for key in sorted(groups_by_scheme):
        idxs = groups_by_scheme[key]
        if not key:
            # No scheme, no comparison: every such isolate stands alone.
            for i in sorted(idxs, key=lambda i: (items[i].label, i)):
                raw_components.append([i])
            continue
        part_edges, part_components = _prim(items, idxs, matrix, min_shared)
        raw_edges.extend(part_edges)
        raw_components.extend(part_components)

    edges = tuple(
        MstEdge(items[i].label, items[j].label, int(d))
        for i, j, d in sorted(
            raw_edges,
            key=lambda e: (int(e[2]), items[e[0]].label, items[e[1]].label),
        )
    )
    components = tuple(
        tuple(items[i].label for i in sorted(comp, key=lambda i: (items[i].label, i)))
        for comp in sorted(
            raw_components,
            key=lambda comp: (-len(comp), min(items[i].label for i in comp)),
        )
    )

    schemes = {_scheme_key(s) for s in items if _scheme_key(s)}
    scheme = schemes.pop() if len(schemes) == 1 else ""

    if raw_edges:
        n_loci_compared = min(d.n_compared for _, _, d in raw_edges)
    else:
        n_loci_compared = min((s.n_called for s in items), default=0)

    mst = Mst(nodes=nodes, edges=edges, components=components, scheme=scheme,
              n_loci_compared=n_loci_compared, threshold=threshold)
    mst.groups = clonal_groups(mst, threshold)
    by_label = {}
    for group in mst.groups:
        for label in group.labels:
            by_label[label] = group.id
    for node in mst.nodes:
        node.group = by_label.get(node.label, 0)
    return mst


# ---------------------------------------------------------------------------
# 5. Clonal groups
# ---------------------------------------------------------------------------

def clonal_groups(mst: Mst, threshold: int = 1) -> Tuple[ClonalGroup, ...]:
    """Connected components once edges longer than ``threshold`` are cut.

    ``threshold=1`` is the classic single-locus-variant clonal complex: two
    isolates stay in one group while they differ at no more than one locus,
    transitively. ``threshold=0`` groups only identical profiles.

    Group ids are stable -- 1-based, assigned after sorting by size (largest
    first) and then by first label -- so the same tree always paints the same
    isolate the same colour, whatever dictionary order the interpreter chose.
    """
    labels = [node.label for node in mst.nodes]
    position = {label: i for i, label in enumerate(labels)}
    parent = list(range(len(labels)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    for edge in mst.edges:
        if edge.distance > threshold:
            continue
        i, j = position.get(edge.a), position.get(edge.b)
        if i is not None and j is not None:
            union(i, j)

    members: Dict[int, List[str]] = {}
    scheme_of: Dict[int, set] = {}
    for i, label in enumerate(labels):
        root = find(i)
        members.setdefault(root, []).append(label)
        scheme_of.setdefault(root, set()).add(mst.nodes[i].scheme)

    ordered = sorted(members.items(), key=lambda kv: (-len(kv[1]), kv[1][0]))
    out: List[ClonalGroup] = []
    for gid, (root, group_labels) in enumerate(ordered, start=1):
        schemes = scheme_of[root]
        out.append(ClonalGroup(
            id=gid,
            labels=tuple(group_labels),
            size=len(group_labels),
            scheme=schemes.pop() if len(schemes) == 1 else "-",
            threshold=threshold,
        ))
    return tuple(out)


# ---------------------------------------------------------------------------
# 6. Layout
# ---------------------------------------------------------------------------

#: Shortest edge drawn, in pixels: two isolates at distance 0 are still two
#: isolates and must not land on top of each other.
EDGE_BASE = 46.0
#: Extra pixels per differing locus.
EDGE_SCALE = 16.0
#: Longest edge drawn. Without a cap a single distant outlier flattens the rest
#: of the tree into an unreadable dot.
EDGE_MAX = 230.0
#: Closest two node centres may ever end up.
MIN_SEP = 34.0
#: How far the layout may be enlarged to fill a roomy canvas.
MAX_ZOOM = 2.0


def edge_length(distance: int, *, base: float = EDGE_BASE,
                scale: float = EDGE_SCALE, cap: float = EDGE_MAX) -> float:
    """Pixels for an edge of ``distance`` loci: ``base + scale * distance``, capped.

    Monotonic, so a longer edge is never drawn shorter than a shorter one, and
    bounded at both ends so the drawing stays readable for both a clonal
    cluster (all zeros) and an outgroup.
    """
    return min(cap, base + scale * max(0, int(distance)))


def _components_from(mst: Mst) -> List[List[int]]:
    """Connected components of the kept edges, recomputed from nodes + edges.

    Deliberately not read from :attr:`Mst.components`: :func:`layout` must also
    work on a tree a caller assembled by hand, and the edge list is the truth.
    """
    labels = [node.label for node in mst.nodes]
    position = {label: i for i, label in enumerate(labels)}
    adjacency: Dict[int, List[int]] = {i: [] for i in range(len(labels))}
    for edge in mst.edges:
        i, j = position.get(edge.a), position.get(edge.b)
        if i is None or j is None:
            continue
        adjacency[i].append(j)
        adjacency[j].append(i)
    seen = [False] * len(labels)
    out: List[List[int]] = []
    for start in range(len(labels)):
        if seen[start]:
            continue
        seen[start] = True
        stack = [start]
        comp = []
        while stack:
            i = stack.pop()
            comp.append(i)
            for j in sorted(adjacency[i]):
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        out.append(sorted(comp))
    return out


def _radial(mst: Mst, comp: Sequence[int], lengths: Dict[Tuple[int, int], float],
            adjacency: Dict[int, List[int]]) -> Dict[int, Tuple[float, float]]:
    """Seeded radial tree layout of one component.

    The root is the most connected node (ties: the smallest label), so a star
    cluster is drawn as a star. Each node owns an angular wedge, splits it
    between its children in proportion to how many leaves each child's subtree
    carries -- a big subtree gets room, a twig does not -- and sits at its
    wedge's bisector, one edge length further out than its parent. That is
    entirely determined by the tree, which is what makes it reproducible: no
    random seeds, no clock, no iteration over a dict.
    """
    if not comp:
        return {}
    degree = {i: len(adjacency[i]) for i in comp}
    root = min(comp, key=lambda i: (-degree[i], mst.nodes[i].label, i))

    # Breadth-first spanning order; children sorted by label for determinism.
    parent: Dict[int, Optional[int]] = {root: None}
    order = [root]
    children: Dict[int, List[int]] = {i: [] for i in comp}
    queue = [root]
    while queue:
        i = queue.pop(0)
        for j in sorted(adjacency[i], key=lambda j: (mst.nodes[j].label, j)):
            if j in parent:
                continue
            parent[j] = i
            children[i].append(j)
            order.append(j)
            queue.append(j)

    # Leaf counts, computed from the deepest node upwards.
    leaves: Dict[int, int] = {}
    for i in reversed(order):
        kids = children[i]
        leaves[i] = 1 if not kids else sum(leaves[k] for k in kids)

    pos: Dict[int, Tuple[float, float]] = {root: (0.0, 0.0)}
    radius: Dict[int, float] = {root: 0.0}
    angle: Dict[int, float] = {root: 0.0}
    wedge: Dict[int, Tuple[float, float]] = {root: (0.0, 2.0 * math.pi)}

    for i in order:
        kids = children[i]
        if not kids:
            continue
        lo, hi = wedge[i]
        span = hi - lo
        total = float(leaves[i])
        cursor = lo
        for k in kids:
            share = span * (leaves[k] / total)
            child_lo, child_hi = cursor, cursor + share
            cursor = child_hi
            a = 0.5 * (child_lo + child_hi)
            r = radius[i] + lengths[(min(i, k), max(i, k))]
            wedge[k] = (child_lo, child_hi)
            angle[k] = a
            radius[k] = r
            pos[k] = (r * math.cos(a), r * math.sin(a))
    return pos


def _relieve(points: List[List[float]], min_sep: float, passes: int = 240) -> None:
    """Push overlapping nodes apart, in place, deterministically.

    Pairs are visited in index order and each pass moves both members of an
    overlapping pair half the shortfall, so the result depends only on the
    starting coordinates. Two nodes sitting exactly on top of each other are
    separated along an angle derived from their indices -- the golden angle, a
    fixed number, not a random one.
    """
    n = len(points)
    if n < 2:
        return
    want = min_sep * min_sep
    for _ in range(passes):
        moved = False
        for i in range(n - 1):
            xi, yi = points[i]
            for j in range(i + 1, n):
                dx = points[j][0] - xi
                dy = points[j][1] - yi
                d2 = dx * dx + dy * dy
                if d2 >= want:
                    continue
                d = math.sqrt(d2)
                if d < 1e-9:
                    a = 2.39996322972865332 * i + 0.7 * j
                    dx, dy, d = math.cos(a), math.sin(a), 1.0
                push = 0.5 * (min_sep - d) + 0.01
                ux, uy = dx / d, dy / d
                points[i][0] -= ux * push
                points[i][1] -= uy * push
                points[j][0] += ux * push
                points[j][1] += uy * push
                xi, yi = points[i]
                moved = True
        if not moved:
            return


def _bbox(points: Sequence[Sequence[float]]) -> Tuple[float, float, float, float]:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return (min(xs), min(ys), max(xs), max(ys))


def layout(mst: Mst, width: float = 800.0, height: float = 600.0, *,
           min_sep: float = MIN_SEP, padding: Optional[float] = None,
           base: float = EDGE_BASE, scale: float = EDGE_SCALE,
           cap: float = EDGE_MAX, max_zoom: float = MAX_ZOOM) -> Tuple[MstNode, ...]:
    """2D coordinates for drawing. Pure, deterministic, no randomness.

    Returns COPIES of ``mst.nodes`` with ``x``/``y`` filled in, in the same
    order -- the input tree is not mutated, so a caller may lay the same tree
    out at several canvas sizes.

    Each component is laid out radially (see :func:`_radial`), overlaps are
    relieved, the components are packed into rows, and the whole drawing is
    centred in the canvas. Edge length grows with allelic distance, and no two
    node centres end up closer than ``min_sep``.

    The one promise the canvas size does NOT get is shrink-to-fit: 200 isolates
    cannot be drawn inside 400x300 pixels without piling nodes on top of each
    other, and an unreadable picture is worse than a scrollable one. The layout
    enlarges to fill a roomy canvas (up to ``max_zoom``) but never scales below
    1, so coordinates may exceed ``width``/``height`` -- the caller scrolls or
    zooms out.
    """
    nodes = tuple(replace(node) for node in mst.nodes)
    n = len(nodes)
    if n == 0:
        return ()
    pad = min_sep if padding is None else padding

    if n == 1:
        nodes[0].x = width / 2.0
        nodes[0].y = height / 2.0
        return nodes

    position = {node.label: i for i, node in enumerate(mst.nodes)}
    adjacency: Dict[int, List[int]] = {i: [] for i in range(n)}
    lengths: Dict[Tuple[int, int], float] = {}
    for edge in mst.edges:
        i, j = position.get(edge.a), position.get(edge.b)
        if i is None or j is None or i == j:
            continue
        adjacency[i].append(j)
        adjacency[j].append(i)
        lengths[(min(i, j), max(i, j))] = edge_length(
            edge.distance, base=base, scale=scale, cap=cap)

    components = _components_from(mst)
    placed: List[List[float]] = [[0.0, 0.0] for _ in range(n)]

    # 1. Lay out and de-overlap each component on its own.
    boxes: List[Tuple[float, float, float, float]] = []
    for comp in components:
        pos = _radial(mst, comp, lengths, adjacency)
        local = [[pos[i][0], pos[i][1]] for i in comp]
        _relieve(local, min_sep)
        for i, point in zip(comp, local):
            placed[i] = point
        boxes.append(_bbox(local))

    # 2. Pack the components into rows, biggest first, left to right.
    gap = min_sep * 1.5
    order = sorted(range(len(components)),
                   key=lambda c: (-len(components[c]),
                                  mst.nodes[components[c][0]].label))
    widths = [boxes[c][2] - boxes[c][0] for c in range(len(components))]
    heights = [boxes[c][3] - boxes[c][1] for c in range(len(components))]
    total_width = sum(widths) + gap * max(0, len(components) - 1)
    aspect = (width / height) if height > 0 else 1.3333
    row_limit = max(max(widths, default=0.0),
                    math.sqrt(max(total_width, 1.0)
                              * max(sum(heights), 1.0) * max(aspect, 0.2)))

    cursor_x = 0.0
    row_y = 0.0
    row_height = 0.0
    for c in order:
        w, h = widths[c], heights[c]
        if cursor_x > 0.0 and cursor_x + w > row_limit:
            cursor_x = 0.0
            row_y += row_height + gap
            row_height = 0.0
        dx = cursor_x - boxes[c][0]
        dy = row_y - boxes[c][1]
        for i in components[c]:
            placed[i][0] += dx
            placed[i][1] += dy
        cursor_x += w + gap
        row_height = max(row_height, h)

    # 3. One global relief pass: packing can bring two components close.
    _relieve(placed, min_sep)

    # 4. Centre, and enlarge (never shrink) towards the canvas size.
    min_x, min_y, max_x, max_y = _bbox(placed)
    span_x = max(max_x - min_x, 1e-9)
    span_y = max(max_y - min_y, 1e-9)
    avail_x = max(width - 2.0 * pad, 1.0)
    avail_y = max(height - 2.0 * pad, 1.0)
    zoom = min(avail_x / span_x, avail_y / span_y)
    zoom = min(max(zoom, 1.0), max(max_zoom, 1.0))

    draw_w = span_x * zoom
    draw_h = span_y * zoom
    off_x = max(pad, (width - draw_w) / 2.0)
    off_y = max(pad, (height - draw_h) / 2.0)
    for i, node in enumerate(nodes):
        node.x = (placed[i][0] - min_x) * zoom + off_x
        node.y = (placed[i][1] - min_y) * zoom + off_y
    return nodes
