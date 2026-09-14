"""Baseline snapshot pointer, snapshot replay and the honest two-forest change summary."""

import json
from copy import deepcopy

import pytest

from wmlstudio.comparison import forest_from_distances, pairwise_distances
from wmlstudio.investigation import (
    DIFF_CAVEATS,
    INTERPRETATION,
    InvestigationStore,
    build_snapshot,
    snapshot_diff,
    snapshot_graph,
)
from wmlstudio.project import Project


def profile(sid, vector, **extra):
    return {"sample_id": sid, "sample_name": sid, "scheme": "cgMLST", "scheme_digest": "pinned",
            "status": "complete", "alleles": dict(zip("abcd", vector)), **extra}


def snapshot(results, threshold=1, overlap=0.95, previous=None):
    return build_snapshot(results, pairwise_distances(results, overlap), threshold, overlap,
                          previous=previous)


def stored_investigation(tmp_path, cohort):
    """A real project holding one saved investigation plus its imported profiles."""
    project = Project(tmp_path / "dual.wmlstudio")
    ids = [project.add_profile(sid, profile(sid, vector)) for sid, vector in cohort]
    store = InvestigationStore(project)
    plan = store.save("Ward A", ids, scheme_digest="pinned", scheme="cgMLST",
                      protocol="Local pilot protocol; threshold 1")
    return project, store, plan, ids


def test_baseline_points_at_the_first_snapshot_and_later_analysis_never_rewrites_it(tmp_path):
    project, store, plan, ids = stored_investigation(tmp_path, [("a", "1111"), ("b", "2111")])
    first_results = [profile(ids[0], "1111"), profile(ids[1], "2111")]
    first = store.save_snapshot(plan["id"], snapshot(first_results))
    added = project.add_profile("c", profile("c", "2211"))
    second_results = first_results + [profile(added, "2211")]
    second = store.save_snapshot(plan["id"], snapshot(second_results, previous=first))
    assert store.get(plan["id"])["snapshots"] == [first["snapshot_id"], second["snapshot_id"]]
    assert store.baseline_snapshot_id(plan["id"]) == first["snapshot_id"]
    frozen = store.snapshot(plan["id"], first["snapshot_id"])
    assert frozen == first
    assert [p["sample_id"] for p in frozen["profiles"]] == ids
    with pytest.raises(ValueError, match="immutable"):
        store.save_snapshot(plan["id"], first)
    project.close()


def test_pinning_a_baseline_records_history_without_touching_the_concurrency_fence(tmp_path):
    project, store, plan, ids = stored_investigation(tmp_path, [("a", "1111"), ("b", "2111")])
    results = [profile(ids[0], "1111"), profile(ids[1], "2111")]
    first = store.save_snapshot(plan["id"], snapshot(results))
    second = store.save_snapshot(plan["id"], snapshot(results, threshold=2, previous=first))
    fence = store.get(plan["id"])["updated_at"]
    store.set_baseline(plan["id"], second["snapshot_id"])
    reloaded = store.get(plan["id"])
    assert reloaded["updated_at"] == fence
    assert store.baseline_snapshot_id(plan["id"]) == second["snapshot_id"]
    pinned = [h for h in project.history() if h["action"] == "investigation_baseline_pinned"]
    assert [h["details"]["snapshot_id"] for h in pinned] == [second["snapshot_id"]]
    with pytest.raises(KeyError):
        store.set_baseline(plan["id"], "not-a-snapshot")
    store.set_baseline(plan["id"], None)
    assert store.baseline_snapshot_id(plan["id"]) == first["snapshot_id"]
    project.close()


def test_a_pinned_baseline_missing_from_the_snapshot_list_falls_back_to_the_first(tmp_path):
    project, store, plan, ids = stored_investigation(tmp_path, [("a", "1111"), ("b", "2111")])
    results = [profile(ids[0], "1111"), profile(ids[1], "2111")]
    first = store.save_snapshot(plan["id"], snapshot(results))
    project.set_setting(f"investigation.{plan['id']}",
                        {**store.get(plan["id"]), "baseline_snapshot_id": "vanished"})
    assert store.baseline_snapshot_id(plan["id"]) == first["snapshot_id"]
    project.close()


def test_the_baseline_pointer_survives_plan_edits_new_snapshots_and_reopening(tmp_path):
    project, store, plan, ids = stored_investigation(tmp_path, [("a", "1111"), ("b", "2111")])
    results = [profile(ids[0], "1111"), profile(ids[1], "2111")]
    first = store.save_snapshot(plan["id"], snapshot(results))
    second = store.save_snapshot(plan["id"], snapshot(results, threshold=2, previous=first))
    store.set_baseline(plan["id"], second["snapshot_id"])
    store.save("Ward A revised", ids, investigation_id=plan["id"], scheme_digest="pinned",
               scheme="cgMLST", threshold=3)
    assert store.baseline_snapshot_id(plan["id"]) == second["snapshot_id"]
    store.save_snapshot(plan["id"], snapshot(results, threshold=3, previous=second))
    assert store.baseline_snapshot_id(plan["id"]) == second["snapshot_id"]
    project.close()
    with Project(tmp_path / "dual.wmlstudio") as reopened:
        reloaded = InvestigationStore(reopened)
        assert reloaded.baseline_snapshot_id(plan["id"]) == second["snapshot_id"]
        assert reloaded.snapshot(plan["id"], second["snapshot_id"]) == second


def test_snapshot_catalog_labels_every_snapshot_in_order_with_its_policy(tmp_path):
    project, store, plan, ids = stored_investigation(tmp_path, [("a", "1111"), ("b", "2111")])
    results = [profile(ids[0], "1111"), profile(ids[1], "2111")]
    first = store.save_snapshot(plan["id"], snapshot(results))
    second = store.save_snapshot(plan["id"], snapshot(results, threshold=3, previous=first))
    catalog = store.snapshot_catalog(plan["id"])
    assert [entry["snapshot_id"] for entry in catalog] == [first["snapshot_id"], second["snapshot_id"]]
    assert [entry["number"] for entry in catalog] == [1, 2]
    assert [entry["threshold"] for entry in catalog] == [1, 3]
    assert [entry["cohort_size"] for entry in catalog] == [2, 2]
    assert all(entry["scheme_digest"] == "pinned" and entry["created_at"] for entry in catalog)
    assert [entry["is_baseline"] for entry in catalog] == [True, False]
    store.set_baseline(plan["id"], second["snapshot_id"])
    assert [entry["is_baseline"] for entry in store.snapshot_catalog(plan["id"])] == [False, True]
    project.close()


def test_snapshot_graph_replays_the_stored_edges_without_recomputing_distances():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211"), profile("d", "3333")]
    rows = pairwise_distances(records)
    stored = build_snapshot(records, rows, 1, 0.95)
    results, edges, groups = snapshot_graph(stored)
    expected = [(e["source"], e["target"], e["distance"]) for e in forest_from_distances(rows)]
    assert [(e["source"], e["target"], e["distance"]) for e in edges] == expected
    assert groups == stored["groups"]
    assert [r["sample_id"] for r in results] == ["a", "b", "c", "d"]
    for result, stored_profile in zip(results, stored["profiles"]):
        assert {"sample_id", "sample_name", "metadata", "primary_st"} <= set(result)
        assert len(result["alleles"]) == stored_profile["callable_loci"] == 4
        assert stored_profile["total_loci"] == 4


def test_snapshot_graph_leaves_an_unprofiled_isolate_isolated_instead_of_near():
    records = [profile("a", "1111"), profile("b", "2111"),
               profile("pending", [], status="not_profiled", comparison_unavailable=True)]
    stored = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    results, edges, groups = snapshot_graph(stored)
    assert {r["sample_id"] for r in results} == {"a", "b", "pending"}
    assert next(r for r in results if r["sample_id"] == "pending")["alleles"] == {}
    assert all("pending" not in (edge["source"], edge["target"]) for edge in edges)
    placeholder = next(p for p in stored["profiles"] if p["sample_id"] == "pending")
    assert (placeholder["callable_loci"], placeholder["total_loci"]) == (0, 0)
    assert next(g for g in groups if g["members"] == ["pending"])["status"] == "not_comparable"


def test_snapshot_graph_and_snapshot_diff_never_mutate_the_stored_snapshots():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    stored = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    baseline, current = deepcopy(stored), deepcopy(stored)
    results, edges, groups = snapshot_graph(current)
    results[0]["alleles"]["a"] = "999"
    edges[0]["distance"] = 99
    groups[0]["members"].append("ghost")
    snapshot_diff(baseline, current)
    assert current == stored and baseline == stored


def test_identical_snapshots_report_one_policy_and_no_change():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = build_snapshot(records, pairwise_distances(records), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert diff["policy"]["comparable"] is True and diff["policy"]["identical"] is True
    assert diff["policy"]["differences"] == []
    assert diff["isolates"]["added"] == [] and diff["isolates"]["removed"] == []
    assert diff["isolates"]["evidence_changed"] == []
    assert diff["pairs"]["distance_changed"] == [] and diff["pairs"]["unchanged_count"] == 3
    assert diff["pairs"]["not_assessed_count"] == 0
    assert diff["mst_edges"]["only_in_baseline"] == [] and diff["mst_edges"]["only_in_current"] == []
    assert diff["mst_edges"]["shared"] == 2
    assert diff["groups"]["merged"] == [] and diff["groups"]["split"] == []
    assert diff["groups"]["dissolved"] == []


def test_an_added_isolate_is_counted_as_not_assessed_never_as_an_unchanged_distance():
    original = [profile("a", "1111"), profile("b", "2111"), profile("d", "3333"), profile("e", "4333")]
    first = build_snapshot(original, pairwise_distances(original), 1, 0.95)
    grown = original + [profile("c", "2211")]
    second = build_snapshot(grown, pairwise_distances(grown), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert [row["sample_id"] for row in diff["isolates"]["added"]] == ["c"]
    assert diff["isolates"]["removed"] == [] and diff["isolates"]["retained"] == ["a", "b", "d", "e"]
    assert diff["pairs"]["distance_changed"] == []
    assert diff["pairs"]["not_assessed_count"] == 4
    assert diff["pairs"]["unchanged_count"] == 6
    assert diff["summary"]["isolates_added"] == 1
    assert any("cohort changed" in caveat for caveat in diff["caveats"])


def test_a_cluster_that_gains_a_new_isolate_reports_added_members_not_a_merge():
    original = [profile("a", "1111"), profile("b", "2111"), profile("d", "3333"), profile("e", "4333")]
    first = build_snapshot(original, pairwise_distances(original), 1, 0.95)
    grown = original + [profile("c", "2211")]
    second = build_snapshot(grown, pairwise_distances(grown), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    entry = next(row for row in diff["groups"]["crosswalk"] if "a" in row["current_members"])
    assert entry["added_members"] == ["c"] and entry["carried_members"] == ["a", "b"]
    assert [origin["shared_members"] for origin in entry["from_baseline_groups"]] == [["a", "b"]]
    assert diff["groups"]["merged"] == [] and diff["groups"]["split"] == []
    assert diff["groups"]["dissolved"] == []


def test_a_removed_isolate_never_reads_as_a_cluster_split():
    original = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    first = build_snapshot(original, pairwise_distances(original), 1, 0.95)
    smaller = original[:2]
    second = build_snapshot(smaller, pairwise_distances(smaller), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert [row["sample_id"] for row in diff["isolates"]["removed"]] == ["c"]
    assert diff["groups"]["split"] == [] and diff["groups"]["dissolved"] == []
    entry = next(row for row in diff["groups"]["crosswalk"] if row["current_status"] == "cluster")
    assert entry["carried_members"] == ["a", "b"] and entry["added_members"] == []
    assert diff["pairs"]["not_assessed_count"] == 2


def test_a_real_split_names_the_dissolved_baseline_group_and_both_heirs():
    records = [profile("a", "1111"), profile("b", "2111"), profile("d", "3311")]
    first = build_snapshot(records, pairwise_distances(records), 2, 0.95)
    second = build_snapshot(records, pairwise_distances(records), 1, 0.95, previous=first)
    parent = next(g for g in first["groups"] if g["status"] == "cluster")
    diff = snapshot_diff(first, second)
    assert diff["policy"]["comparable"] is True and diff["policy"]["grouping_attributable"] is True
    assert diff["groups"]["split"] == [parent["id"]]
    dissolved = next(row for row in diff["groups"]["dissolved"] if row["baseline_group_id"] == parent["id"])
    assert sorted(heir["shared_members"] for heir in dissolved["went_to"]) == [["a", "b"], ["d"]]
    assert dissolved["lost_members"] == [] and dissolved["members_intact"] is False
    assert any("link threshold changed from 2 to 1" in caveat for caveat in diff["caveats"])
    assert diff["summary"]["groups_split"] == 1


def test_a_threshold_change_merges_clusters_and_is_flagged_as_policy_attributable():
    records = [profile("a", "1111"), profile("b", "2111"), profile("d", "3333"), profile("e", "4333")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = build_snapshot(records, pairwise_distances(records), 4, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert diff["policy"]["identical"] is False and diff["policy"]["grouping_attributable"] is True
    assert {"field": "threshold", "baseline": 1, "current": 4} in diff["policy"]["differences"]
    assert diff["pairs"]["distance_changed"] == [] and diff["pairs"]["unchanged_count"] == 6
    merged = next(row for row in diff["groups"]["crosswalk"]
                  if row["current_group_id"] in diff["groups"]["merged"])
    assert len(merged["from_baseline_groups"]) == 2 and merged["added_members"] == []
    assert len(diff["groups"]["dissolved"]) == 2
    assert any("partly a consequence of that cutoff" in caveat for caveat in diff["caveats"])


def test_a_distance_over_fewer_shared_loci_is_flagged_as_a_changed_denominator():
    first = build_snapshot([profile("a", "1111"), profile("b", "2111")],
                           pairwise_distances([profile("a", "1111"), profile("b", "2111")], 0.5), 1, 0.5)
    retyped = [profile("a", "1111"), profile("b", ["2", "2", None, "1"])]
    second = build_snapshot(retyped, pairwise_distances(retyped, 0.5), 1, 0.5, previous=first)
    diff = snapshot_diff(first, second)
    row = diff["pairs"]["distance_changed"][0]
    assert (row["source"], row["target"]) == ("a", "b")
    assert (row["baseline_distance"], row["current_distance"]) == (1, 2)
    assert (row["baseline_shared"], row["baseline_total"]) == (4, 4)
    assert (row["current_shared"], row["current_total"]) == (3, 4)
    assert row["denominator_changed"] is True and row["distance_changed"] is True
    assert diff["summary"]["denominator_changed"] == 1
    assert [entry["sample_id"] for entry in diff["isolates"]["evidence_changed"]] == ["b"]
    assert diff["isolates"]["evidence_changed"][0]["baseline_callable"] == 4
    assert diff["isolates"]["evidence_changed"][0]["current_callable"] == 3
    assert any("different locus set" in caveat for caveat in diff["caveats"])


def test_an_unchanged_number_over_a_changed_locus_set_is_still_reported():
    first = build_snapshot([profile("a", "1111"), profile("b", "2111")],
                           pairwise_distances([profile("a", "1111"), profile("b", "2111")], 0.5), 1, 0.5)
    retyped = [profile("a", "1111"), profile("b", ["2", "1", None, "1"])]
    second = build_snapshot(retyped, pairwise_distances(retyped, 0.5), 1, 0.5, previous=first)
    diff = snapshot_diff(first, second)
    row = diff["pairs"]["distance_changed"][0]
    assert row["baseline_distance"] == row["current_distance"] == 1
    assert row["distance_changed"] is False and row["denominator_changed"] is True
    assert diff["pairs"]["unchanged_count"] == 0 and diff["pairs"]["denominator_only_count"] == 1


def test_comparability_gained_or_lost_is_reported_apart_from_distance_change():
    mixed = [profile("a", "1111"), profile("b", "2111", status="mixed"), profile("c", "2211")]
    clean = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    first = build_snapshot(mixed, pairwise_distances(mixed), 1, 0.95)
    second = build_snapshot(clean, pairwise_distances(clean), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert diff["pairs"]["distance_changed"] == []
    assert sorted((row["source"], row["target"]) for row in diff["pairs"]["became_comparable"]) == [
        ("a", "b"), ("b", "c")]
    assert all("mixed" in row["baseline_reason"] for row in diff["pairs"]["became_comparable"])
    assert diff["pairs"]["became_excluded"] == []
    reverse = snapshot_diff(second, first)
    assert sorted((row["source"], row["target"]) for row in reverse["pairs"]["became_excluded"]) == [
        ("a", "b"), ("b", "c")]
    assert reverse["pairs"]["became_comparable"] == []
    assert reverse["pairs"]["still_excluded_count"] == 0


def test_a_changed_reference_refuses_the_numeric_diff_instead_of_diffing_silently():
    records = [profile("a", "1111"), profile("b", "2111")]
    updated = [profile("a", "1111", scheme_digest="v2-panel"), profile("b", "2111", scheme_digest="v2-panel")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = build_snapshot(updated, pairwise_distances(updated), 1, 0.95)
    diff = snapshot_diff(first, second)
    assert diff["policy"]["comparable"] is False and diff["policy"]["identical"] is False
    assert diff["pairs"] is None and diff["mst_edges"] is None and diff["groups"] is None
    assert diff["isolates"]["retained"] == ["a", "b"]
    assert diff["isolates"]["added"] == [] and diff["isolates"]["removed"] == []
    assert diff["isolates"]["evidence_changed"] is None
    assert diff["summary"]["distance_changed"] is None and diff["summary"]["groups_split"] is None
    assert "scheme_digest" in diff["policy"]["reason"]
    assert any("not measured on one scale" in caveat for caveat in diff["caveats"])


def test_an_unrecorded_or_changed_metric_version_also_refuses_the_numeric_diff():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = {**deepcopy(first), "metric_version": "shared-unambiguous-alleles-v2"}
    assert snapshot_diff(first, second)["policy"]["comparable"] is False
    blank = {**deepcopy(first), "metric_version": ""}
    diff = snapshot_diff(blank, {**deepcopy(first), "metric_version": ""})
    assert diff["policy"]["comparable"] is False and diff["pairs"] is None


def test_spanning_forest_edge_differences_are_labelled_as_a_drawing_choice():
    original = [profile("a", "1111"), profile("b", "2111"), profile("d", "3333"), profile("e", "4333")]
    first = build_snapshot(original, pairwise_distances(original), 1, 0.95)
    grown = original + [profile("c", "2211")]
    second = build_snapshot(grown, pairwise_distances(grown), 1, 0.95, previous=first)
    diff = snapshot_diff(first, second)
    assert ("b", "c") in [(row["source"], row["target"]) for row in diff["mst_edges"]["only_in_current"]]
    assert diff["mst_edges"]["only_in_baseline"] == []
    assert "arbitrarily" in diff["mst_edges"]["note"]


def test_the_change_summary_is_json_safe_and_always_carries_its_caveats():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    grown = records + [profile("c", "2211")]
    second = build_snapshot(grown, pairwise_distances(grown), 5, 0.90, previous=first)
    diff = snapshot_diff(first, second)
    assert json.loads(json.dumps(diff, allow_nan=False))["format_version"] == 1
    assert all(caveat in diff["caveats"] for caveat in DIFF_CAVEATS)
    assert INTERPRETATION in diff["caveats"]
    assert any("not a phylogeny" in caveat for caveat in diff["caveats"])
    assert any("transmission" in caveat for caveat in diff["caveats"])
    assert diff["policy"]["distance_scale_attributable"] is True
    assert any("minimum shared-locus fraction changed" in caveat for caveat in diff["caveats"])
    assert diff["baseline"]["cohort_size"] == 2 and diff["current"]["cohort_size"] == 3


def test_a_change_summary_needs_both_halves():
    records = [profile("a", "1111"), profile("b", "2111")]
    stored = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    with pytest.raises(ValueError, match="baseline and a current"):
        snapshot_diff(stored, None)
    with pytest.raises(ValueError, match="baseline and a current"):
        snapshot_diff(None, stored)
