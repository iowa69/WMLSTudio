import itertools

import pytest

from wmlstudio.hierarchical import (
    build_tree,
    cut_tree,
    dendrogram_layout,
    merge_heights,
)
from wmlstudio.sequence import AnalysisCancelled


def row(source, target, distance=None, **extra):
    """One pairwise row in the shape comparison.pairwise_distances returns."""
    comparable = distance is not None
    return {"source": source, "target": target, "distance": distance,
            "comparable": comparable,
            "reason": "" if comparable else "Shared callable loci are below the minimum overlap.",
            **extra}


def square(distances):
    """Rows for every pair named in a {(left, right): distance} mapping."""
    return [row(left, right, value) for (left, right), value in sorted(distances.items())]


# Worked by hand for all three linkages in the tests below.
FOUR = {("a", "b"): 2, ("a", "c"): 6, ("a", "d"): 10,
        ("b", "c"): 5, ("b", "d"): 9, ("c", "d"): 4}

# Four isolates in a line, each one step from the next and three from the far end.
CHAIN = {("a", "b"): 1, ("b", "c"): 1, ("c", "d"): 1,
         ("a", "c"): 2, ("b", "d"): 2, ("a", "d"): 3}


def steps(tree):
    return [(merge["members"], merge["height"]) for merge in tree["merges"]]


def test_single_linkage_joins_clusters_at_their_closest_crossing_pair():
    """Guards the linkage arithmetic: single linkage must use the minimum, not the
    first or the mean, of the distances crossing between two clusters."""
    tree = build_tree(square(FOUR), linkage="single")
    assert steps(tree) == [(["a", "b"], 2), (["c", "d"], 4), (["a", "b", "c", "d"], 5)]
    assert tree["complete"] is True
    assert tree["roots"][0]["members"] == ["a", "b", "c", "d"]


def test_complete_linkage_joins_clusters_at_their_furthest_crossing_pair():
    """The same matrix must give a different top height under complete linkage;
    reusing a single-linkage height here would understate a cluster's spread."""
    tree = build_tree(square(FOUR), linkage="complete")
    assert steps(tree) == [(["a", "b"], 2), (["c", "d"], 4), (["a", "b", "c", "d"], 10)]


def test_average_linkage_takes_the_mean_over_every_crossing_pair():
    """UPGMA averages all four crossing distances (6+10+5+9)/4, not the average of
    the two cluster-level numbers, which would weight the smaller cluster wrongly."""
    tree = build_tree(square(FOUR), linkage="average")
    assert steps(tree) == [(["a", "b"], 2), (["c", "d"], 4), (["a", "b", "c", "d"], 7.5)]


def test_an_uncompared_pair_never_becomes_a_height_of_zero():
    """The central failure mode: isolates that share too few loci must not fall into
    one cluster at height zero, which would read on screen as an identical pair."""
    tree = build_tree([row("a", "b"), row("a", "c"), row("b", "c")], linkage="single")
    assert tree["merges"] == []
    assert len(tree["roots"]) == 3
    assert tree["isolates_without_a_comparable_pair"] == ["a", "b", "c"]
    groups = cut_tree(tree, 0)["groups"]
    assert [group["status"] for group in groups] == ["not_comparable"] * 3
    assert all(group["name"] == "Not comparable" for group in groups)


def test_refusing_to_merge_across_an_uncompared_pair_keeps_the_clusters_apart():
    """a and c were never compared, so the fact that both are 9 from b is no evidence
    about them; they must not be joined, and the reason must be visible."""
    tree = build_tree([row("a", "b", 9), row("b", "c", 9), row("a", "c")], linkage="single")
    assert steps(tree) == [(["a", "b"], 9)]
    assert tree["complete"] is False
    assert [root["members"] for root in tree["roots"]] == [["a", "b"], ["c"]]
    blocked = tree["blocked_merges"]
    assert len(blocked) == 1
    assert blocked[0]["uncompared_pairs"] == 1 and blocked[0]["compared_pairs"] == 1
    assert "below the minimum overlap" in blocked[0]["reasons"][0]


def test_no_cut_height_however_large_joins_clusters_across_an_uncompared_pair():
    """A user who raises the cut far above every measured distance must still not be
    shown a and c in one group, because nothing was ever measured between them."""
    tree = build_tree([row("a", "b", 9), row("b", "c", 9), row("a", "c")], linkage="single")
    cut = cut_tree(tree, 10_000)
    assert [group["members"] for group in cut["groups"]] == [["a", "b"], ["c"]]
    assert cut["groups_kept_apart_by_uncompared_pairs"] == 1
    assert "never compared" in cut["note"]
    assert cut["groups"][1]["status"] == "singleton"


def test_provisional_policy_merges_but_marks_what_it_could_not_measure():
    """The other allowed policy must still not hide the gap: the merge is labelled
    provisional and carries the count of pairs that were never compared."""
    rows = [row("a", "b", 9), row("b", "c", 9), row("a", "c")]
    tree = build_tree(rows, linkage="single", missing="provisional")
    assert steps(tree) == [(["a", "b"], 9), (["a", "b", "c"], 9)]
    top = tree["merges"][-1]
    assert top["provisional"] is True
    assert top["uncompared_pairs_across"] == 1
    assert top["compared_pairs_across"] == 1
    assert top["uncompared_within_pairs"] == 1
    group = cut_tree(tree, 9)["groups"][0]
    assert group["provisional"] is True and group["uncompared_within_pairs"] == 1
    assert cut_tree(tree, 9)["provisional_groups"] == ["hc:a"]


def test_the_policy_for_uncompared_pairs_is_stated_in_the_result():
    """A view cannot explain a missing branch it was never told about, so both the
    policy name and a sentence describing it travel with the tree and the cut."""
    tree = build_tree(square(FOUR))
    assert tree["missing_policy"] == "refuse"
    assert "never" in tree["missing_data_rule"]
    assert tree["missing_policy_note"] == cut_tree(tree, 3)["missing_policy_note"]
    assert "not comparable blocked the join" in tree["missing_policy_note"]


def test_merges_are_identical_whatever_order_tied_rows_arrive_in():
    """Four isolates at distance zero from each other tie at every step; the merge
    order must come from the isolate names, never from the order of the input."""
    rows = square({pair: 0 for pair in itertools.combinations("abcd", 2)})
    expected = [(["a", "b"], 0), (["a", "b", "c"], 0), (["a", "b", "c", "d"], 0)]
    for permutation in itertools.permutations(rows):
        assert steps(build_tree(permutation, linkage="average")) == expected


def test_single_linkage_reports_the_span_a_chain_of_short_links_hides():
    """Three one-step links join a and d, which are three steps apart. The cluster
    must report that span rather than let the height of 1 stand for the group."""
    tree = build_tree(square(CHAIN), linkage="single")
    assert [merge["height"] for merge in tree["merges"]] == [1, 1, 1]
    top = tree["merges"][-1]
    assert top["chained"] is True and top["max_within_distance"] == 3
    assert "furthest compared pair differs by 3" in top["chain_note"]
    assert tree["chained_merges"] == ["node:2", "node:3"]
    group = cut_tree(tree, 1)["groups"][0]
    assert group["size"] == 4 and group["chained"] is True
    assert group["max_within_distance"] == 3 and group["link_heights"] == [1, 1, 1]
    assert "cut at 1" in group["chain_note"]


def test_complete_linkage_groups_are_never_reported_as_chained():
    """Complete linkage bounds a cluster's furthest pair by its own height, so a
    chaining warning there would be a false alarm a user would learn to ignore."""
    tree = build_tree(square(CHAIN), linkage="complete")
    assert steps(tree) == [(["a", "b"], 1), (["c", "d"], 1), (["a", "b", "c", "d"], 3)]
    assert tree["chained_merges"] == []
    cut = cut_tree(tree, 3)
    assert cut["chained_groups"] == []
    assert cut["groups"][0]["max_within_distance"] == 3


def test_cutting_below_a_merge_keeps_its_two_halves_apart():
    """A cut must apply only the merges at or under it; applying the top merge of
    the hand-worked tree at height 4 would join clusters five apart."""
    tree = build_tree(square(FOUR), linkage="single")
    assert [group["members"] for group in cut_tree(tree, 4)["groups"]] == [["a", "b"], ["c", "d"]]
    assert [group["members"] for group in cut_tree(tree, 1)["groups"]] == [
        ["a"], ["b"], ["c"], ["d"]]
    assert cut_tree(tree, 5)["groups"][0]["members"] == ["a", "b", "c", "d"]


def test_groups_are_named_from_their_own_members_not_their_position():
    """Group identity has to survive a re-cut, so it is anchored on the smallest
    member; a purely positional identifier would silently point somewhere else."""
    distances = {("a", "b"): 1, ("a", "c"): 1, ("b", "c"): 1, ("d", "e"): 1}
    distances.update({pair: 20 for pair in itertools.product("abc", "de")})
    tree = build_tree(square(distances), linkage="single")
    cut = cut_tree(tree, 1)
    assert [(g["id"], g["name"], g["members"]) for g in cut["groups"]] == [
        ("hc:a", "HC group 001", ["a", "b", "c"]),
        ("hc:d", "HC group 002", ["d", "e"]),
    ]
    assert cut["clusters"] == 2 and cut["clustered_isolates"] == 5
    wider = cut_tree(tree, 20)
    assert wider["groups"][0]["id"] == "hc:a" and wider["groups"][0]["size"] == 5


def test_an_isolate_with_no_rows_at_all_still_appears_as_its_own_group():
    """A sample dropped from the pairwise stage must not vanish from the cohort; a
    tree that silently loses it would report fewer isolates than were examined."""
    tree = build_tree([row("a", "b", 1)], labels=["a", "b", "lonely"])
    assert tree["labels"] == ["a", "b", "lonely"]
    assert tree["isolates_without_a_comparable_pair"] == ["lonely"]
    groups = cut_tree(tree, 1)["groups"]
    assert [group["members"] for group in groups] == [["a", "b"], ["lonely"]]
    assert groups[1]["status"] == "not_comparable"


def test_rows_measuring_different_quantities_refuse_to_share_one_tree():
    """A 7-locus MLST distance and a cgMLST distance on one axis would be read as
    one scale. The tree refuses the mixture instead of drawing it."""
    rows = [row("a", "b", 1, typing_kind="mlst"), row("b", "c", 1, typing_kind="cgmlst"),
            row("a", "c", 1, typing_kind="mlst")]
    with pytest.raises(ValueError, match="mix different quantities"):
        build_tree(rows)


@pytest.mark.parametrize("distance", [None, -1, True, float("nan"), float("inf"), "3"])
def test_a_pair_marked_comparable_without_a_usable_distance_is_refused(distance):
    """A row claiming comparability but carrying no real number is a bug upstream;
    coercing it would put a fabricated height into the dendrogram."""
    rows = [{"source": "a", "target": "b", "comparable": True, "distance": distance}]
    with pytest.raises(ValueError):
        build_tree(rows)


@pytest.mark.parametrize("rows, message", [
    ([row("a", "a", 1)], "with itself"),
    ([row("a", "b", 1), row("b", "a", 2)], "Duplicate"),
    ([{"source": "", "target": "b", "comparable": True, "distance": 1}], "source identifier"),
])
def test_malformed_pairwise_input_is_rejected_by_name(rows, message):
    """Bad rows must fail loudly here rather than produce a tree that looks fine."""
    with pytest.raises(ValueError, match=message):
        build_tree(rows)


@pytest.mark.parametrize("bad", ["single-ish", "ward", ""])
def test_an_unknown_linkage_is_refused(bad):
    """Silently falling back to a default linkage would mislabel the whole figure."""
    with pytest.raises(ValueError, match="Linkage must be"):
        build_tree(square(FOUR), linkage=bad)


@pytest.mark.parametrize("bad", [-1, True, float("nan"), "3"])
def test_a_cut_height_must_be_a_finite_nonnegative_number(bad):
    """A nonsense cut would otherwise return groups that no threshold produced."""
    with pytest.raises(ValueError, match="cut height"):
        cut_tree(build_tree(square(FOUR)), bad)


def test_merge_heights_lists_each_distinct_height_once():
    """The heights a user may cut at are the heights the tree actually has."""
    assert merge_heights(build_tree(square(FOUR), linkage="average")) == [2, 4, 7.5]
    assert merge_heights(build_tree(square(CHAIN), linkage="single")) == [1]
    assert merge_heights(build_tree([])) == []


def test_an_empty_or_single_isolate_cohort_produces_no_merges():
    """The clustering has to survive the cohorts a user builds while still adding
    samples, rather than raising on the way to the second isolate."""
    assert build_tree([])["merges"] == []
    assert cut_tree(build_tree([]), 5)["groups"] == []
    alone = build_tree([], labels=["a"])
    assert alone["merges"] == [] and alone["leaf_order"] == ["a"]
    assert cut_tree(alone, 5)["groups"][0]["status"] == "not_comparable"


def test_dendrogram_layout_places_every_node_at_its_own_merge_height():
    """The drawing coordinates are arithmetic that has to be right before any view
    exists: a node drawn at the wrong height misstates the distance it marks."""
    tree = build_tree(square(FOUR), linkage="single")
    assert tree["leaf_order"] == ["a", "b", "c", "d"]
    layout = dendrogram_layout(tree)
    assert [leaf["x"] for leaf in layout["leaves"]] == [0.0, 1.0, 2.0, 3.0]
    assert [(node["x"], node["y"]) for node in layout["nodes"]] == [
        (0.5, 2.0), (2.5, 4.0), (1.5, 5.0)]
    assert layout["max_height"] == 5.0
    assert layout["roots"] == ["node:3"]


def test_the_measured_quantity_travels_with_the_tree_and_its_layout():
    """A dendrogram axis with no units invites reading a cgMLST height as a SNP
    count, so the caption the cohort was built with is carried through unchanged."""
    caption = "Allele differences over 2358 cgMLST targets"
    tree = build_tree(square(FOUR), scale_caption=caption)
    assert tree["scale_caption"] == caption
    assert cut_tree(tree, 3)["scale_caption"] == caption
    assert dendrogram_layout(tree)["height_caption"] == caption


def test_a_forest_lays_out_every_root_and_loses_no_leaf():
    """When uncompared pairs leave several trees, all of them must still be drawn;
    a layout covering only the first root would hide half the cohort."""
    rows = [row("a", "b", 1), row("c", "d", 1), row("a", "c"), row("a", "d"),
            row("b", "c"), row("b", "d")]
    tree = build_tree(rows, linkage="average")
    assert tree["complete"] is False
    layout = dendrogram_layout(tree)
    assert [leaf["label"] for leaf in layout["leaves"]] == ["a", "b", "c", "d"]
    assert layout["roots"] == ["node:1", "node:2"]
    assert [(node["x"], node["y"]) for node in layout["nodes"]] == [(0.5, 1.0), (2.5, 1.0)]


def test_clustering_stops_when_the_user_cancels():
    """A cohort of a few hundred isolates takes seconds, so the run happens off the
    interface thread and has to answer a cancel like every other analysis here."""
    with pytest.raises(AnalysisCancelled):
        build_tree(square(FOUR), cancelled=lambda: True)


def test_rows_from_a_partially_compared_cohort_keep_their_measured_heights():
    """Uncompared pairs must not drag an average down: the height of a provisional
    merge is the mean of the pairs that were compared, not of all of them."""
    rows = [row("a", "b", 2), row("a", "c", 4), row("b", "c"), row("a", "d", 10),
            row("b", "d", 10), row("c", "d", 10)]
    tree = build_tree(rows, linkage="average", missing="provisional")
    merge = next(m for m in tree["merges"] if m["members"] == ["a", "b", "c"])
    # (a,c) is the only compared pair crossing from {a,b} to {c}; (b,c) is absent.
    assert merge["height"] == 4
    assert merge["compared_pairs_across"] == 1 and merge["uncompared_pairs_across"] == 1
    assert tree["merges"][-1]["height"] == 10
