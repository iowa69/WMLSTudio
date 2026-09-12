import itertools

import pytest

from wmlstudio.comparison import minimum_spanning_forest, pairwise_distances


def result(name, alleles, **kwargs):
    return {"sample_name": name, "scheme": "demo", "scheme_digest": "fingerprint",
            "alleles": alleles, **kwargs}


def test_distance_is_shared_unambiguous_mismatch_count():
    a = result("a", {"x": "1", "y": "2", "z": None})
    b = result("b", {"x": "1", "y": "3", "z": "4"})
    pair = pairwise_distances([a, b], min_overlap=0.6)[0]
    assert pair["distance"] == 1
    assert pair["shared_loci"] == 2
    assert pair["total_loci"] == 3
    assert pair["overlap"] == pytest.approx(2 / 3)
    assert pairwise_distances([a, b])[0]["distance"] is None


@pytest.mark.parametrize("allele", [None, "", "?", "0", "unknown", "novel", "1,2", ["1", "2"]])
def test_no_overlap_never_becomes_zero_distance(allele):
    pair = pairwise_distances(
        [result("a", {"x": allele}), result("b", {"x": allele})], min_overlap=0,
    )[0]
    assert pair["distance"] is None
    assert pair["shared_loci"] == 0
    assert "No shared" in pair["reason"]


@pytest.mark.parametrize("override", [
    {"scheme": "other"}, {"scheme_digest": "new revision"}, {"scheme_digest": ""},
    {"scheme_digest": None}, {"status": "mixed"},
])
def test_incompatible_profiles_do_not_connect(override):
    samples = [result("a", {"x": "1"}), result("b", {"x": "1"}, **override)]
    assert pairwise_distances(samples)[0]["comparable"] is False
    assert minimum_spanning_forest(samples) == []


def test_mixed_and_ambiguous_calls_are_excluded_even_with_allele_ids():
    a = result("a", {"x": "1", "y": "2", "z": "3"}, calls=[
        {"locus": "x", "status": "mixed"}, {"locus": "y", "status": "ambiguous"},
    ])
    b = result("b", {"x": "1", "y": "2", "z": "4"})
    pair = pairwise_distances([a, b], min_overlap=0.3)[0]
    assert pair["distance"] == 1
    assert pair["shared_loci"] == 1
    assert pair["total_loci"] == 3


def test_locus_union_counts_missing_keys_in_overlap():
    samples = [result("a", {"x": "1"}), result("b", {"x": "1", "y": "2"})]
    pair = pairwise_distances(samples)[0]
    assert pair["overlap"] == 0.5
    assert pair["distance"] is None


def test_forest_is_deterministic_under_ties_and_disconnected_inputs():
    samples = [result(name, {"x": "1"}) for name in "cba"]
    samples.append(result("isolated", {"x": None}))
    expected = [("a", "b", 0), ("a", "c", 0)]
    for permutation in itertools.permutations(samples):
        edges = minimum_spanning_forest(permutation)
        assert [(edge["source"], edge["target"], edge["distance"]) for edge in edges] == expected


def test_kruskal_selects_minimum_total_weight():
    samples = [
        result("a", {"x": "1", "y": "1", "z": "1"}),
        result("b", {"x": "2", "y": "1", "z": "1"}),
        result("c", {"x": "2", "y": "2", "z": "1"}),
        result("d", {"x": "2", "y": "2", "z": "2"}),
    ]
    assert sum(edge["distance"] for edge in minimum_spanning_forest(samples)) == 3
    assert minimum_spanning_forest([]) == []
    assert minimum_spanning_forest(samples[:1]) == []


def test_stable_ids_allow_duplicate_display_names():
    samples = [result("same", {"x": "1"}, sample_id="id2"),
               result("same", {"x": "2"}, sample_id="id1")]
    assert pairwise_distances(samples)[0]["source"] == "id1"
    with pytest.raises(ValueError, match="Duplicate"):
        pairwise_distances([samples[0], samples[0]])


@pytest.mark.parametrize("threshold", [-1, 2, float("nan"), float("inf"), True])
def test_invalid_threshold_is_rejected(threshold):
    with pytest.raises(ValueError):
        pairwise_distances([], threshold)
