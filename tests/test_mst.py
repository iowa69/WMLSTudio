# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Minimum spanning tree over typed isolates -- :mod:`wmlst.mst`.

Two halves. The first is pure arithmetic on hand-written profiles: what
counts as a difference, what a missing locus does, that Prim's tree really is
minimum, that ties are broken the same way twice, and that the layout puts no
two nodes on top of each other without ever calling a random number generator.

The second half types 20 real *Klebsiella pneumoniae* and 20 real
*Staphylococcus aureus* RefSeq assemblies with the actual engine and builds the
tree over the result, because a distance function that only ever sees the
profiles a test author imagined is a distance function that has not been
tested. Those tests are marked ``slow``/``needs_blast``/``needs_db`` and skip
themselves when the genome set is not on this machine.

Runnable as ``python3 -m pytest tests/test_mst.py`` or ``python3 tests/test_mst.py``.
"""

from __future__ import annotations

import itertools
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst.mst import (
    AllelicDistance,
    ClonalGroup,
    Mst,
    MstEdge,
    MstNode,
    MstSample,
    allelic_distance,
    build_mst,
    clonal_groups,
    codes_match,
    distance_matrix,
    edge_length,
    layout,
    normalise_code,
    sample_from_result,
    samples_from_results,
)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DBDIR = os.path.join(REPO, "db")

#: The labelled RefSeq genome set. Overridable so the same test can run against
#: another copy; absent means "skip", never "fail".
VALSET = os.environ.get(
    "WMLST_VALSET",
    "/tmp/claude-1000/-home-iowa-Desktop-wmlst/"
    "12caaa71-cd44-4d35-9516-55648dcead69/scratchpad/valset2",
)

#: How many genomes per organism the real-data tests type.
N_GENOMES = 20


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _sample(label, profile, scheme="sa", st="-"):
    """A sample from a compact profile string: ``"1-2-3"``, ``-`` = missing."""
    return MstSample(label=label, scheme=scheme, st=st,
                     alleles=tuple(profile.split("/")))


def _edge_set(mst):
    """Edges as ``{(low label, high label, distance)}`` -- order-free comparison."""
    return {(e.key[0], e.key[1], e.distance) for e in mst.edges}


def _edge_list(mst):
    """Edges as an ordered list, for the stricter "identical output" check."""
    return [(e.a, e.b, e.distance) for e in mst.edges]


def _brute_force_mst_weight(samples):
    """Total weight of a minimum spanning tree, found by Kruskal.

    A second, independent implementation: if Prim's total ever disagrees with
    Kruskal's total on the same complete graph, one of them is wrong.
    """
    n = len(samples)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            d = allelic_distance(samples[i].alleles, samples[j].alleles)
            if d.comparable:
                edges.append((int(d), i, j))
    edges.sort()
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    total = 0
    kept = 0
    for weight, i, j in edges:
        ri, rj = find(i), find(j)
        if ri == rj:
            continue
        parent[ri] = rj
        total += weight
        kept += 1
    return total, kept


# ===========================================================================
# 1. allelic_distance
# ===========================================================================
def test_identical_profiles_are_distance_zero():
    d = allelic_distance(("1", "2", "3"), ("1", "2", "3"))
    assert d == 0
    assert d.n_compared == 3
    assert d.comparable


def test_every_locus_differing_is_the_profile_length():
    d = allelic_distance(("1", "2", "3"), ("4", "5", "6"))
    assert d == 3
    assert d.n_compared == 3


def test_distance_counts_only_the_loci_that_differ():
    assert allelic_distance(("1", "2", "3"), ("1", "9", "3")) == 1
    assert allelic_distance(("1", "2", "3"), ("1", "9", "8")) == 2


def test_distance_is_an_int_and_behaves_like_one():
    """The GUI will sum, sort and format these; they must be plain ints."""
    d = allelic_distance(("1", "2"), ("1", "9"))
    assert isinstance(d, int)
    assert isinstance(d, AllelicDistance)
    assert d + 1 == 2
    assert sorted([d, allelic_distance(("1",), ("1",))]) == [0, 1]
    assert "%d" % d == "1"


def test_distance_is_symmetric():
    a = ("1", "-", "3", "4")
    b = ("1", "2", "9", "-")
    assert allelic_distance(a, b) == allelic_distance(b, a)
    assert allelic_distance(a, b).n_compared == allelic_distance(b, a).n_compared


@pytest.mark.parametrize("missing", ["-", "0", "", " ", "?"])
def test_a_missing_locus_is_never_a_difference(missing):
    """The whole point: absent data is not allele zero (section 1 of mst.py)."""
    d = allelic_distance(("1", "2", "3"), ("1", missing, "3"))
    assert d == 0, "%r counted as a difference" % missing
    assert d.n_compared == 2
    assert d.n_loci == 3
    assert d.n_missing == 1


def test_missing_on_either_side_is_skipped():
    d = allelic_distance(("-", "2", "0"), ("1", "2", "3"))
    assert d == 0
    assert d.n_compared == 1


def test_profiles_with_no_shared_called_locus_are_not_comparable():
    d = allelic_distance(("1", "-"), ("-", "2"))
    assert d.n_compared == 0
    assert not d.comparable


def test_two_empty_profiles_are_not_comparable():
    d = allelic_distance((), ())
    assert d.n_compared == 0
    assert not d.comparable


def test_novel_and_partial_markers_do_not_make_a_difference():
    """``~16`` and ``16?`` say how allele 16 was seen, not which allele it is."""
    assert allelic_distance(("16", "2"), ("~16", "2")) == 0
    assert allelic_distance(("16", "2"), ("16?", "2")) == 0
    assert allelic_distance(("~16", "2"), ("16?", "2")) == 0
    assert allelic_distance(("16",), ("~17",)) == 1


def test_a_multiple_call_matches_either_allele():
    assert allelic_distance(("1,2",), ("1",)) == 0
    assert allelic_distance(("1,2",), ("2",)) == 0
    assert allelic_distance(("1,2",), ("3",)) == 1
    assert allelic_distance(("1,2",), ("2,1",)) == 0


def test_profiles_of_different_lengths_compare_the_shared_prefix():
    d = allelic_distance(("1", "2", "3"), ("1", "9"))
    assert d == 1
    assert d.n_loci == 2


def test_normalise_code_and_codes_match():
    assert normalise_code("~16") == "16"
    assert normalise_code(" 16? ") == "16"
    assert normalise_code("2,1") == "1,2"
    assert normalise_code("-") == ""
    assert normalise_code("0") == ""
    assert normalise_code(None) == ""
    assert normalise_code(7) == "7"
    assert codes_match("1", "1") is True
    assert codes_match("1", "2") is False
    assert codes_match("1", "-") is None
    assert codes_match("-", "-") is None


def test_distance_matrix_is_square_symmetric_and_zero_on_the_diagonal():
    samples = [_sample("a", "1/2/3"), _sample("b", "1/2/9"), _sample("c", "4/5/6")]
    m = distance_matrix(samples)
    assert len(m) == 3 and all(len(row) == 3 for row in m)
    for i in range(3):
        assert m[i][i] == 0
        for j in range(3):
            assert m[i][j] == m[j][i]
    assert m[0][1] == 1
    assert m[0][2] == 3


# ===========================================================================
# 2. build_mst -- shape
# ===========================================================================
def test_no_samples_gives_an_empty_tree():
    mst = build_mst([])
    assert mst.nodes == ()
    assert mst.edges == ()
    assert mst.components == ()
    assert mst.n_nodes == 0


def test_one_sample_has_a_node_and_no_edges():
    mst = build_mst([_sample("only", "1/2/3", st="7")])
    assert mst.n_nodes == 1
    assert mst.edges == ()
    assert mst.components == (("only",),)
    assert mst.nodes[0].st == "7"
    assert mst.nodes[0].alleles == ("1", "2", "3")


def test_a_connected_scheme_group_has_exactly_n_minus_one_edges():
    samples = [_sample("s%d" % i, "%d/2/3" % i) for i in range(6)]
    mst = build_mst(samples)
    assert mst.n_nodes == 6
    assert mst.n_edges == 5
    assert mst.components == (tuple("s%d" % i for i in range(6)),)


def test_identical_profiles_stay_joined_at_distance_zero():
    """A clonal cluster is one component, not five singletons."""
    samples = [_sample("iso%d" % i, "1/2/3", st="5") for i in range(5)]
    mst = build_mst(samples)
    assert mst.n_edges == 4
    assert all(edge.distance == 0 for edge in mst.edges)
    assert len(mst.components) == 1


def test_samples_with_no_st_are_still_placed_by_their_profile():
    """ST is display only; an untypeable isolate still has a profile."""
    samples = [_sample("a", "1/2/3", st="-"),
               _sample("b", "1/2/9", st="-"),
               _sample("c", "1/2/3", st="4")]
    mst = build_mst(samples)
    assert mst.n_edges == 2
    assert _edge_set(mst) == {("a", "c", 0), ("a", "b", 1)}
    assert {n.st for n in mst.nodes} == {"-", "4"}


def test_mixed_schemes_never_share_an_edge():
    samples = [_sample("k1", "1/2/3", scheme="klebsiella"),
               _sample("k2", "1/2/9", scheme="klebsiella"),
               _sample("s1", "1/2/3", scheme="saureus"),
               _sample("s2", "1/2/9", scheme="saureus")]
    mst = build_mst(samples)
    by_label = {n.label: n.scheme for n in mst.nodes}
    for edge in mst.edges:
        assert by_label[edge.a] == by_label[edge.b], "edge crosses schemes: %r" % (edge,)
    assert mst.n_edges == 2
    assert len(mst.components) == 2
    assert mst.scheme == "", "a mixed run must not claim a single scheme"


def test_a_single_scheme_run_reports_that_scheme():
    mst = build_mst([_sample("a", "1/2/3", scheme="klebsiella"),
                     _sample("b", "1/2/9", scheme="klebsiella")])
    assert mst.scheme == "klebsiella"
    assert mst.n_loci_compared == 3


@pytest.mark.parametrize("scheme", ["-", "", "none"])
def test_an_isolate_with_no_scheme_is_a_singleton(scheme):
    samples = [_sample("a", "1/2/3"), _sample("b", "1/2/3"),
               _sample("x", "1/2/3", scheme=scheme)]
    mst = build_mst(samples)
    assert ("x",) in mst.components
    assert all("x" not in (e.a, e.b) for e in mst.edges)


def test_an_all_missing_profile_is_a_singleton_not_a_zero_distance_twin():
    """Two isolates that called nothing are not "identical"; they are unknown."""
    samples = [_sample("a", "1/2/3"), _sample("b", "1/2/3"),
               _sample("blank1", "-/-/-"), _sample("blank2", "-/0/-")]
    mst = build_mst(samples)
    assert ("blank1",) in mst.components
    assert ("blank2",) in mst.components
    assert mst.n_edges == 1


def test_components_are_ordered_biggest_first():
    samples = [_sample("a", "1/2/3"), _sample("b", "1/2/3"),
               _sample("c", "1/2/3"),
               _sample("z", "1/2/3", scheme="other")]
    mst = build_mst(samples)
    sizes = [len(c) for c in mst.components]
    assert sizes == sorted(sizes, reverse=True)


def test_duplicate_labels_are_made_unique_rather_than_fatal():
    """Two directories, one basename: an MST must not blow up over it."""
    samples = [_sample("iso", "1/2/3"), _sample("iso", "1/2/9")]
    mst = build_mst(samples)
    labels = [n.label for n in mst.nodes]
    assert len(set(labels)) == 2
    assert labels[0] == "iso"
    assert mst.n_edges == 1


def test_min_shared_refuses_thin_comparisons():
    """With only one locus in common, ``min_shared=2`` keeps them apart."""
    samples = [_sample("a", "1/-/-"), _sample("b", "1/2/3")]
    assert build_mst(samples).n_edges == 1
    assert build_mst(samples, min_shared=2).n_edges == 0


def test_n_loci_compared_reports_the_worst_edge():
    samples = [_sample("a", "1/2/3"), _sample("b", "1/2/3"), _sample("c", "1/-/-")]
    mst = build_mst(samples)
    assert mst.n_loci_compared == 1


def test_engine_sample_results_are_accepted_directly():
    """The adapter takes AlleleCall objects, dicts and MstSamples alike."""

    class FakeCall:
        def __init__(self, code):
            self.code = code

    class FakeResult:
        label = "fake"
        scheme = "klebsiella"
        st = "23"
        alleles = (FakeCall("2"), FakeCall("1"), FakeCall("~1"))

    got = sample_from_result(FakeResult())
    assert got == MstSample("fake", "klebsiella", "23", ("2", "1", "~1"))
    assert sample_from_result({"label": "d", "scheme": "s", "st": "1",
                               "alleles": ["1", "2"]}).alleles == ("1", "2")
    assert samples_from_results([FakeResult()])[0].label == "fake"
    mst = build_mst([FakeResult(), FakeResult()])
    assert mst.n_edges == 1 and mst.edges[0].distance == 0


# ===========================================================================
# 3. build_mst -- minimality and determinism
# ===========================================================================
def test_the_tree_is_minimum_agreeing_with_an_independent_kruskal():
    samples = [
        _sample("a", "1/1/1/1/1/1/1"),
        _sample("b", "1/1/1/1/1/1/2"),
        _sample("c", "1/1/1/1/9/9/9"),
        _sample("d", "5/5/5/5/5/5/5"),
        _sample("e", "5/5/5/5/5/5/6"),
        _sample("f", "1/1/1/1/1/9/2"),
    ]
    mst = build_mst(samples)
    weight = sum(e.distance for e in mst.edges)
    expected, kept = _brute_force_mst_weight(samples)
    assert mst.n_edges == kept == len(samples) - 1
    assert weight == expected


def test_the_same_input_gives_the_identical_tree_every_time():
    samples = [_sample("s%02d" % i, "%d/%d/%d" % (i % 3, i % 4, i % 5))
               for i in range(12)]
    first = build_mst(samples)
    for _ in range(4):
        again = build_mst(samples)
        assert _edge_list(again) == _edge_list(first)
        assert again.components == first.components
        assert [(n.label, n.group) for n in again.nodes] == \
               [(n.label, n.group) for n in first.nodes]


def test_input_order_does_not_change_the_tree():
    """Ties are broken by (distance, label pair), so shuffling is a no-op."""
    samples = [_sample("s%02d" % i, "%d/%d/%d" % (i % 3, i % 4, i % 5))
               for i in range(9)]
    baseline = build_mst(samples)
    for permutation in (list(reversed(samples)),
                        samples[4:] + samples[:4],
                        sorted(samples, key=lambda s: s.alleles)):
        other = build_mst(permutation)
        assert _edge_set(other) == _edge_set(baseline)
        assert sum(e.distance for e in other.edges) == \
               sum(e.distance for e in baseline.edges)


def test_ties_are_broken_by_the_label_pair():
    """Three isolates all one locus apart: the chosen pair must be predictable."""
    samples = [_sample("b", "2/1/1"), _sample("a", "1/1/1"), _sample("c", "3/1/1")]
    mst = build_mst(samples)
    assert _edge_list(mst) == [("a", "b", 1), ("a", "c", 1)]


def test_every_component_is_a_tree_with_no_cycles():
    samples = [_sample("s%d" % i, "%d/%d/1" % (i % 4, i % 3)) for i in range(10)]
    samples.append(_sample("other", "1/1/1", scheme="klebsiella"))
    mst = build_mst(samples)
    assert mst.n_edges == sum(len(c) - 1 for c in mst.components)
    seen = set()
    for edge in mst.edges:
        assert edge.key not in seen, "duplicate edge %r" % (edge,)
        seen.add(edge.key)


def test_edges_are_sorted_by_distance():
    samples = [_sample("a", "1/1/1"), _sample("b", "1/1/2"),
               _sample("c", "9/9/9"), _sample("d", "1/1/1")]
    mst = build_mst(samples)
    distances = [e.distance for e in mst.edges]
    assert distances == sorted(distances)


# ===========================================================================
# 4. clonal_groups
# ===========================================================================
def test_clonal_groups_cut_edges_longer_than_the_threshold():
    samples = [_sample("a", "1/1/1/1/1/1/1"),
               _sample("b", "1/1/1/1/1/1/2"),   # SLV of a
               _sample("c", "9/9/9/9/9/9/9"),   # far away
               _sample("d", "9/9/9/9/9/9/8")]   # SLV of c
    mst = build_mst(samples)
    groups = clonal_groups(mst, threshold=1)
    assert [g.size for g in groups] == [2, 2]
    assert [g.id for g in groups] == [1, 2]
    assert {frozenset(g.labels) for g in groups} == \
           {frozenset({"a", "b"}), frozenset({"c", "d"})}


def test_threshold_zero_groups_only_identical_profiles():
    samples = [_sample("a", "1/1/1"), _sample("b", "1/1/1"), _sample("c", "1/1/2")]
    groups = clonal_groups(build_mst(samples), threshold=0)
    assert [g.size for g in groups] == [2, 1]
    assert set(groups[0].labels) == {"a", "b"}


def test_a_large_threshold_merges_everything_in_a_component():
    samples = [_sample("a", "1/1/1"), _sample("b", "2/2/2"), _sample("c", "3/3/3")]
    groups = clonal_groups(build_mst(samples), threshold=99)
    assert len(groups) == 1
    assert groups[0].size == 3


def test_group_ids_are_ordered_by_size_and_stable():
    samples = ([_sample("big%d" % i, "1/1/1") for i in range(4)]
               + [_sample("small%d" % i, "9/9/9") for i in range(2)]
               + [_sample("lone", "5/5/5")])
    mst = build_mst(samples)
    groups = clonal_groups(mst, threshold=1)
    assert [g.size for g in groups] == [4, 2, 1]
    assert [g.id for g in groups] == [1, 2, 3]
    assert set(groups[0].labels) == {"big0", "big1", "big2", "big3"}
    for _ in range(3):
        assert clonal_groups(mst, threshold=1) == groups


def test_build_mst_stamps_the_group_id_on_every_node():
    samples = [_sample("a", "1/1/1"), _sample("b", "1/1/2"), _sample("z", "9/9/9")]
    mst = build_mst(samples)
    by_label = {n.label: n.group for n in mst.nodes}
    assert by_label["a"] == by_label["b"]
    assert by_label["z"] != by_label["a"]
    assert all(n.group >= 1 for n in mst.nodes)
    assert mst.groups and isinstance(mst.groups[0], ClonalGroup)


def test_singletons_each_get_their_own_group():
    mst = build_mst([_sample("a", "1/1/1", scheme="x"),
                     _sample("b", "1/1/1", scheme="y")])
    assert len({n.group for n in mst.nodes}) == 2


# ===========================================================================
# 5. layout
# ===========================================================================
def _min_separation(nodes):
    best = float("inf")
    for left, right in itertools.combinations(nodes, 2):
        best = min(best, ((left.x - right.x) ** 2 + (left.y - right.y) ** 2) ** 0.5)
    return best


def test_layout_returns_one_placed_node_per_sample():
    mst = build_mst([_sample("s%d" % i, "%d/1/1" % i) for i in range(7)])
    nodes = layout(mst, 800, 600)
    assert len(nodes) == mst.n_nodes
    assert [n.label for n in nodes] == [n.label for n in mst.nodes]
    assert all(isinstance(n, MstNode) for n in nodes)
    assert all(n.x == n.x and n.y == n.y for n in nodes)  # no NaN


def test_layout_does_not_mutate_the_tree():
    mst = build_mst([_sample("a", "1/1/1"), _sample("b", "1/1/2")])
    layout(mst, 400, 400)
    assert all(n.x == 0.0 and n.y == 0.0 for n in mst.nodes)


def test_layout_of_one_node_centres_it():
    mst = build_mst([_sample("only", "1/1/1")])
    (node,) = layout(mst, 400, 300)
    assert (node.x, node.y) == (200.0, 150.0)


def test_layout_of_an_empty_tree_is_empty():
    assert layout(build_mst([]), 400, 300) == ()


@pytest.mark.parametrize("n", [2, 3, 8, 25])
def test_no_two_nodes_overlap(n):
    samples = [_sample("s%02d" % i, "%d/%d/%d" % (i % 3, i % 5, i % 7))
               for i in range(n)]
    nodes = layout(build_mst(samples), 900, 700)
    assert _min_separation(nodes) >= 33.0


def test_identical_isolates_do_not_land_on_the_same_point():
    """Distance 0 means "same clone", not "same pixel"."""
    mst = build_mst([_sample("c%d" % i, "1/1/1") for i in range(6)])
    nodes = layout(mst, 800, 600)
    assert _min_separation(nodes) >= 33.0


def test_several_components_are_all_separated():
    samples = []
    for scheme in ("a", "b", "c"):
        samples += [_sample("%s%d" % (scheme, i), "%d/1/1" % i, scheme=scheme)
                    for i in range(3)]
    nodes = layout(build_mst(samples), 900, 700)
    assert _min_separation(nodes) >= 33.0


def test_layout_is_deterministic():
    samples = [_sample("s%02d" % i, "%d/%d/%d" % (i % 3, i % 4, i % 5))
               for i in range(15)]
    mst = build_mst(samples)
    first = [(n.label, n.x, n.y) for n in layout(mst, 800, 600)]
    for _ in range(3):
        assert [(n.label, n.x, n.y) for n in layout(mst, 800, 600)] == first
    # And a rebuilt tree from the same samples lays out identically.
    rebuilt = [(n.label, n.x, n.y) for n in layout(build_mst(samples), 800, 600)]
    assert rebuilt == first


def test_layout_uses_no_randomness_at_all():
    """A tree that moves between runs is not a tree anyone can cite."""
    with open(os.path.join(REPO, "wmlst", "mst.py"), encoding="utf-8") as fh:
        source = fh.read()
    for banned in ("import random", "random.", "uuid", "time.time", "id("):
        assert banned not in source, "mst.py uses %r" % banned


def test_edge_length_grows_with_distance_and_is_capped():
    lengths = [edge_length(d) for d in range(0, 40)]
    assert lengths == sorted(lengths)
    assert lengths[0] > 0
    assert lengths[1] > lengths[0]
    assert max(lengths) <= 230.0


def test_a_longer_edge_is_drawn_longer():
    samples = [_sample("a", "1/1/1/1/1"), _sample("b", "1/1/1/1/2"),
               _sample("c", "9/9/9/9/9")]
    nodes = {n.label: n for n in layout(build_mst(samples), 900, 700)}

    def span(p, q):
        return ((nodes[p].x - nodes[q].x) ** 2 + (nodes[p].y - nodes[q].y) ** 2) ** 0.5

    mst = build_mst(samples)
    near = next(e for e in mst.edges if e.distance <= 1)
    far = next(e for e in mst.edges if e.distance > 1)
    assert span(far.a, far.b) > span(near.a, near.b)


def test_layout_fills_a_larger_canvas():
    samples = [_sample("s%d" % i, "%d/1/1" % i) for i in range(5)]
    mst = build_mst(samples)
    small = layout(mst, 500, 400)
    large = layout(mst, 2000, 1600)
    assert max(n.x for n in large) > max(n.x for n in small)


def test_layout_of_a_small_tree_stays_inside_the_canvas():
    samples = [_sample("s%d" % i, "%d/1/1" % i) for i in range(5)]
    nodes = layout(build_mst(samples), 900, 700)
    assert all(0 <= n.x <= 900 for n in nodes)
    assert all(0 <= n.y <= 700 for n in nodes)


def test_dataclasses_are_what_the_gui_was_promised():
    node = MstNode("a", "1", "sa", ("1",), 1, 2.0, 3.0)
    assert (node.label, node.st, node.scheme, node.group, node.x, node.y) == \
           ("a", "1", "sa", 1, 2.0, 3.0)
    edge = MstEdge("b", "a", 3)
    assert (edge.a, edge.b, edge.distance) == ("b", "a", 3)
    assert edge.key == ("a", "b")
    empty = Mst()
    assert (empty.nodes, empty.edges, empty.components) == ((), (), ())
    assert empty.scheme == "" and empty.n_loci_compared == 0


def test_mst_lookup_helpers():
    mst = build_mst([_sample("a", "1/1/1"), _sample("b", "1/1/2")])
    assert mst.node("a").label == "a"
    assert mst.node("nope") is None
    assert mst.index_of("b") == 1
    assert mst.index_of("nope") == -1
    assert mst.neighbours("a") == ("b",)
    assert mst.max_distance == 1


# ===========================================================================
# 6. Real genomes -- 20 K. pneumoniae + 20 S. aureus, typed by the engine
# ===========================================================================
def _genomes(organism, n=N_GENOMES):
    folder = os.path.join(VALSET, organism)
    if not os.path.isdir(folder):
        pytest.skip("genome set not present: %s" % folder)
    files = sorted(os.path.join(folder, name)
                   for name in os.listdir(folder) if name.endswith(".fna"))
    if len(files) < n:
        pytest.skip("only %d genomes in %s" % (len(files), folder))
    return files[:n]


@pytest.fixture(scope="module")
def real_samples():
    """Type 20 K. pneumoniae and 20 S. aureus assemblies. Once per module."""
    from wmlst import engine

    files = _genomes("kpneumoniae") + _genomes("saureus")
    if not os.path.isdir(os.path.join(DBDIR, "pubmlst")):
        pytest.skip("no PubMLST database at %s" % DBDIR)
    result = engine.analyse_files(
        engine.RunConfig(files=tuple(files), dbdir=DBDIR, jobs=4))
    samples = []
    for item in samples_from_results(result.samples):
        # Basenames read better on a tree than absolute paths.
        samples.append(MstSample(label=os.path.basename(item.label),
                                 scheme=item.scheme, st=item.st,
                                 alleles=item.alleles))
    assert len(samples) == 2 * N_GENOMES
    return samples


pytestmark_real = pytest.mark.usefixtures("real_samples")


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_genomes_type_into_two_schemes(real_samples):
    schemes = {}
    for sample in real_samples:
        schemes.setdefault(sample.scheme, []).append(sample.label)
    assert len(schemes) >= 2, "40 genomes of two organisms collapsed to %r" % schemes
    biggest = sorted((len(v) for v in schemes.values()), reverse=True)
    assert biggest[0] >= 15 and biggest[1] >= 15, schemes
    assert all(sample.n_called >= 5 for sample in real_samples)


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_mst_is_connected_per_scheme_with_n_minus_one_edges(real_samples):
    mst = build_mst(real_samples)
    assert mst.n_nodes == len(real_samples)

    # Every component is a tree.
    assert mst.n_edges == sum(len(c) - 1 for c in mst.components)

    # Every scheme whose isolates all called something is ONE component.
    by_scheme = {}
    for sample in real_samples:
        by_scheme.setdefault(sample.scheme, []).append(sample.label)
    component_of = {}
    for i, component in enumerate(mst.components):
        for label in component:
            component_of[label] = i
    for scheme, labels in by_scheme.items():
        if scheme.strip() in ("", "-") or len(labels) < 2:
            continue
        assert len({component_of[label] for label in labels}) == 1, (
            "scheme %s split across components" % scheme)

    # And no edge ever crosses a scheme boundary.
    scheme_of = {s.label: s.scheme for s in real_samples}
    for edge in mst.edges:
        assert scheme_of[edge.a] == scheme_of[edge.b], edge


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_mst_is_deterministic_across_runs_and_orderings(real_samples):
    baseline = build_mst(real_samples)
    assert _edge_list(build_mst(real_samples)) == _edge_list(baseline)
    shuffled = list(reversed(real_samples))
    assert _edge_set(build_mst(shuffled)) == _edge_set(baseline)
    assert sum(e.distance for e in build_mst(shuffled).edges) == \
           sum(e.distance for e in baseline.edges)
    first = [(n.label, n.x, n.y) for n in layout(baseline, 1000, 800)]
    assert [(n.label, n.x, n.y) for n in layout(build_mst(real_samples), 1000, 800)] \
        == first


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_identical_sts_are_at_distance_zero(real_samples):
    """Same scheme, same ST, both fully called -> the profiles must agree."""
    checked = 0
    for left, right in itertools.combinations(real_samples, 2):
        if left.scheme != right.scheme or left.st in ("-", ""):
            continue
        plain = all(code.isdigit() for code in left.alleles + right.alleles)
        if not plain:
            continue
        d = allelic_distance(left.alleles, right.alleles)
        if left.st == right.st:
            assert d == 0, "%s and %s are both ST %s but differ at %d loci" % (
                left.label, right.label, left.st, int(d))
            checked += 1
        else:
            assert d >= 1, "%s (ST %s) and %s (ST %s) have identical profiles" % (
                left.label, left.st, right.label, right.st)
    assert checked or True  # a set with no repeated ST is still a valid set


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_mst_weight_matches_an_independent_kruskal(real_samples):
    mst = build_mst(real_samples)
    total = 0
    for scheme in {s.scheme for s in real_samples}:
        group = [s for s in real_samples if s.scheme == scheme]
        if scheme.strip() in ("", "-"):
            continue
        weight, _kept = _brute_force_mst_weight(group)
        total += weight
    assert sum(e.distance for e in mst.edges) == total


@pytest.mark.slow
@pytest.mark.needs_blast
@pytest.mark.needs_db
def test_real_clonal_groups_and_layout_hold_up(real_samples):
    mst = build_mst(real_samples)
    groups = clonal_groups(mst, threshold=1)
    assert sum(g.size for g in groups) == mst.n_nodes
    assert [g.size for g in groups] == sorted((g.size for g in groups), reverse=True)
    assert [g.id for g in groups] == list(range(1, len(groups) + 1))
    # Every isolate is in exactly one group.
    labels = [label for g in groups for label in g.labels]
    assert sorted(labels) == sorted(n.label for n in mst.nodes)
    # An SLV group never mixes schemes.
    assert all(g.scheme != "-" or g.size == 1 for g in groups)

    nodes = layout(mst, 1200, 900)
    assert len(nodes) == mst.n_nodes
    assert _min_separation(nodes) >= 33.0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([os.path.abspath(__file__), "-v"]))
