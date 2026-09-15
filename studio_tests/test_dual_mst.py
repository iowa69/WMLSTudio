"""Baseline snapshot pointer, snapshot replay and the honest two-forest change summary."""

import json
from copy import deepcopy

import pytest

from wmlstudio import ui_compare
from wmlstudio.comparison import forest_from_distances, pairwise_distances
from wmlstudio.investigation import (
    DIFF_CAVEATS,
    DIFF_ROW_FIELDS,
    INTERPRETATION,
    InvestigationStore,
    build_snapshot,
    snapshot_diff,
    snapshot_diff_caption,
    snapshot_diff_html,
    snapshot_diff_rows,
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


# ---------------------------------------------------------------------------
# The change summary as a reader sees it
# ---------------------------------------------------------------------------


def test_the_rendered_change_summary_escapes_names_and_states_its_limits():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    grown = records + [dict(profile("c", "2211"), sample_name="Ward <b>A</b>")]
    second = build_snapshot(grown, pairwise_distances(grown), 1, 0.95, previous=first)
    page = snapshot_diff_html(snapshot_diff(first, second))
    assert "Ward &lt;b&gt;A&lt;/b&gt;" in page and "<b>A</b>" not in page
    assert "not a phylogeny" in page and "transmission chain" in page
    assert "not assessed" in page.casefold()


def test_a_refused_comparison_renders_the_refusal_and_no_distance_table():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = build_snapshot([dict(row, scheme_digest="other") for row in records],
                            pairwise_distances([dict(row, scheme_digest="other") for row in records]), 1, 0.95)
    diff = snapshot_diff(first, second)
    page = snapshot_diff_html(diff)
    assert "not measured on one scale" in page
    assert "Baseline shared/total" not in page
    assert page.count("Not assessed") >= 3          # distances, groups and drawn lines
    assert snapshot_diff_caption(diff).endswith("distances, lines and groups not comparable")


def test_the_change_summary_table_names_what_was_not_assessed():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    grown = records + [profile("c", "2211")]
    second = build_snapshot(grown, pairwise_distances(grown), 1, 0.95, previous=first)
    rows = snapshot_diff_rows(snapshot_diff(first, second))
    assert all(tuple(row) == DIFF_ROW_FIELDS for row in rows)
    kinds = {row["change_type"]: row for row in rows}
    assert kinds["isolate_added"]["sample_id"] == "c"
    assert kinds["pairs_not_assessed"]["current_value"] == 2
    assert "not assessed is not" in kinds["pairs_not_assessed"]["note"].casefold()
    assert kinds["comparison_policy"]["note"].startswith("Same reference")


def test_a_caption_names_only_what_moved_and_says_so_when_nothing_did():
    records = [profile("a", "1111"), profile("b", "2111")]
    first = build_snapshot(records, pairwise_distances(records), 1, 0.95)
    second = build_snapshot(records, pairwise_distances(records), 1, 0.95, previous=first)
    assert snapshot_diff_caption(snapshot_diff(first, second)) == "no change since the baseline"
    changed = [profile("a", "1111"), profile("b", "2211")]
    third = build_snapshot(changed, pairwise_distances(changed), 1, 0.95, previous=first)
    caption = snapshot_diff_caption(snapshot_diff(first, third))
    assert "1 re-typed" in caption and "1 distance changed" in caption


# ---------------------------------------------------------------------------
# The two trees on screen
# ---------------------------------------------------------------------------


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "application")
    errors = []
    monkeypatch.setattr(widget, "error", lambda message: errors.append(str(message)))
    widget.test_errors = errors
    qtbot.addWidget(widget)
    # Laid out with resize, never shown: repainting a second TreeView after many
    # torn-down offscreen widgets crashes inside the graph label painter.
    widget.resize(1280, 860)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


def imported(window, name, vector="1111", ward="ICU"):
    return window.project.add_profile(name, {
        "sample_name": name, "scheme": "Study MLST", "scheme_digest": "study-reference",
        "status": "profile_imported", "alleles": dict(zip("abcd", vector)), "calls": [],
        "st": "20" if vector.startswith("1") else "44", "input_sha256": "a" * 64,
    }, {"organism": {"genus": "Staphylococcus", "species": "aureus"}, "annotations": {"ward": ward}})


def build(window, ids, threshold=1, investigation_id=None, name="Ward A"):
    """Save (or re-save) the investigation and build it, which stores a snapshot."""
    window.refresh_cohort_table()
    plan = InvestigationStore(window.project).save(
        name, list(ids), investigation_id=investigation_id, scheme="Study MLST",
        scheme_digest="study-reference", threshold=threshold, min_overlap=0.95,
        protocol="Synthetic test protocol; not a clinical cutoff")
    window.cohort_ids = set(ids)
    window.select_investigation(plan["id"])
    assert window._current_snapshot, window.tree_status.text()
    return plan["id"]


def grown_investigation(window):
    """A baseline of two isolates and a current comparison of three."""
    a, b = imported(window, "A"), imported(window, "B", "2111")
    investigation_id = build(window, [a, b])
    c = imported(window, "C", "2211", ward="Ward 2")
    build(window, [a, b, c], investigation_id=investigation_id)
    window.dual_toggle.setChecked(True)
    return investigation_id, a, b, c


def test_the_baseline_pane_is_absent_until_it_is_asked_for(window):
    a, b = imported(window, "A"), imported(window, "B", "2111")
    build(window, [a, b])
    assert window.baseline_pane.isVisibleTo(window.graph_split) is False
    assert window.current_caption.isVisibleTo(window.current_pane) is False
    assert window.export_tree_choice.isVisibleTo(window) is False
    assert window.graph_tabs.widget(0) is window.graph_split
    assert window.tree.parent() is window.current_pane
    assert window.baseline_tree._results == {}


def test_each_tree_holds_its_own_cohort_and_says_which_moment_it_is(window):
    _, a, b, c = grown_investigation(window)
    assert set(window.baseline_tree._results) == {a, b}
    assert set(window.tree._results) == {a, b, c}
    assert window.baseline_caption.text().startswith("Baseline · ")
    assert "2 isolates" in window.baseline_caption.text()
    assert "+1 isolate" in window.baseline_caption.text()
    assert window.current_caption.text() == "Current · Classical MLST · 4 loci · 3 isolates · link ≤ 1"
    assert "frozen" in window.baseline_caption.toolTip()
    assert "study-refer" in window.baseline_caption.toolTip()
    assert "nothing here is recalculated" in window.baseline_caption.toolTip()
    assert window.test_errors == []


def test_selecting_in_one_tree_selects_the_same_isolates_in_the_other_once(window):
    _, a, b, c = grown_investigation(window)
    emissions = []
    window.baseline_tree.selectionChanged.connect(lambda ids: emissions.append(tuple(ids)))
    window.tree.selectionChanged.connect(lambda ids: emissions.append(tuple(ids)))
    window.tree.select_ids([a, b])
    assert window.baseline_tree.selected_ids() == sorted([a, b])
    assert len(emissions) == 2                      # one per view, no ping-pong
    window.baseline_tree.select_ids([b])
    assert window.tree.selected_ids() == [b]
    assert len(emissions) == 4


def test_selecting_an_isolate_added_since_the_baseline_says_so_rather_than_nothing(window):
    _, a, b, c = grown_investigation(window)
    window.tree.select_ids([c])
    assert window.baseline_tree.selected_ids() == []
    assert "1 also in the baseline tree" not in window.graph_selection_label.text()
    assert "0 also in the baseline tree · 1 added since the baseline" in window.graph_selection_label.text()
    window.baseline_tree.select_ids([a])
    assert "1 still in the current comparison · 0 no longer in it" in window.graph_selection_label.text()


def test_linking_selection_can_be_switched_off_without_touching_either_tree(window):
    _, a, b, c = grown_investigation(window)
    window.tree.select_ids([a])
    assert window.baseline_tree.selected_ids() == [a]
    window.link_selection.setChecked(False)
    window.tree.select_ids([b])
    assert window.baseline_tree.selected_ids() == [a]


def test_each_tree_keeps_its_own_saved_arrangement(window):
    _, a, b, c = grown_investigation(window)
    window._pending_graph_state = None
    window._pending_baseline_graph_state = None
    window.tree.nodes[a].setPos(11, 12)
    window.tree.update_edges()
    window.baseline_tree.nodes[a].setPos(90, 91)
    window.baseline_tree.update_edges()
    current = window.project.get_setting("graph_style", {})
    baseline = window.project.get_setting("graph_style.baseline", {})
    assert current["positions"][a] == [11, 12]
    assert baseline["positions"][a] == [90, 91]
    assert current["positions"][a] != baseline["positions"][a]
    assert set(baseline["positions"]) == {a, b}


def test_highlighting_marks_matches_in_both_trees_and_stores_nothing(window):
    _, a, b, c = grown_investigation(window)
    stored = window.project.samples()
    window.graph_search.setText("Ward 2")
    assert window.tree.nodes[c].highlighted is True
    assert window.tree.nodes[a].highlighted is False
    assert all(not node.highlighted for node in window.baseline_tree.nodes.values())
    assert "1 matched here · 0 in the baseline tree" in window.graph_selection_label.text()
    window.graph_search.setText("A")
    assert window.baseline_tree.nodes[a].highlighted is True
    assert window.project.samples() == stored


def test_one_shared_colour_key_covers_both_trees_when_colouring_by_st(window):
    _, a, b, c = grown_investigation(window)
    window.color_by.setCurrentIndex(window.color_by.findData("st"))
    assert window.tree.color_by == "st" and window.baseline_tree.color_by == "st"
    shared = set(window.tree.color_categories()) & set(window.baseline_tree.color_categories())
    assert shared
    for category in shared:
        assert window.tree.legend()[category] == window.baseline_tree.legend()[category]
    assert all(category in window.graph_legend.toolTip()
               for category in window.tree.color_categories())


def test_colouring_by_cluster_pins_nothing_because_group_numbers_already_match(window):
    grown_investigation(window)
    window.color_by.setCurrentIndex(window.color_by.findData("cluster"))
    assert window.tree._pinned_legend == {}
    assert window.baseline_tree._pinned_legend == {}


def test_without_a_saved_investigation_the_pane_offers_to_save_one(window):
    imported(window, "A")
    imported(window, "B", "2111")
    window.refresh_comparison()
    window.dual_toggle.setChecked(True)
    assert window.baseline_tree.isVisibleTo(window.baseline_pane) is False
    assert "No baseline tree yet" in window.baseline_notice.text()
    assert window.baseline_action.text().startswith("Save this comparison as an investigation")
    assert "No baseline tree is pinned" in window.changes_view.toPlainText()
    assert window.baseline_caption.text() == "Baseline · none pinned"


def test_a_baseline_larger_than_the_autodraw_limit_waits_to_be_asked(window, monkeypatch):
    count = ui_compare.BASELINE_AUTODRAW_LIMIT + 10
    stored = {"snapshot_id": "big", "created_at": "2026-01-01T00:00:00+00:00", "threshold": 1,
              "min_overlap": 0.95, "scheme": "Study MLST", "scheme_digest": "study-reference",
              "groups": [], "pairs": [],
              "profiles": [{"sample_id": f"s{index}", "sample_name": f"Isolate {index}",
                            "known_alleles": {"a": "1"}, "callable_loci": 1, "total_loci": 4}
                           for index in range(count)]}
    monkeypatch.setattr(InvestigationStore, "baseline_snapshot_id", lambda self, iid: "big")
    monkeypatch.setattr(InvestigationStore, "snapshot", lambda self, iid, sid=None: stored)
    window.active_investigation_id = "pretend"
    window.dual_toggle.setChecked(True)
    assert window._baseline_drawn is False
    assert window.baseline_tree._results == {}
    assert f"({count} isolates; this can take several seconds)" in window.baseline_action.text()
    window.baseline_action.click()
    assert window._baseline_drawn is True
    assert len(window.baseline_tree._results) == count
    assert window.baseline_notice.isVisibleTo(window.baseline_pane) is False


def test_a_baseline_on_another_reference_refuses_the_comparison_on_screen(window, monkeypatch):
    _, a, b, c = grown_investigation(window)
    original = window._baseline_snapshot
    monkeypatch.setattr(InvestigationStore, "snapshot",
                        lambda self, iid, sid=None: dict(original, scheme_digest="another-reference"))
    window._reset_baseline_state()
    window.refresh_baseline_graph()
    assert window._baseline_diff["policy"]["comparable"] is False
    assert window._baseline_diff["pairs"] is None
    page = window.changes_view.toPlainText()
    assert "not measured on one scale" in page
    assert "Baseline shared/total" not in page
    assert "not comparable" in window.baseline_caption.text()
    assert "not measured on one scale" in window.baseline_caption.toolTip()


def test_reporting_from_the_baseline_tree_carries_the_baseline_snapshot(window):
    _, a, b, c = grown_investigation(window)
    window.report_graph_selection([a, b], role="baseline")
    assert window.report_sample_ids() == {a, b}
    assert window._report_investigation_snapshot is window._baseline_snapshot
    window.report_graph_selection([a, b, c], role="current")
    assert window._report_investigation_snapshot is window._current_snapshot


def test_switching_projects_forgets_the_other_projects_baseline(window, tmp_path, qtbot):
    investigation_id, a, b, c = grown_investigation(window)
    assert window._baseline_snapshot is not None
    window.switch_project(tmp_path / "second.wmlstudio")
    qtbot.waitUntil(lambda: window.comparison_worker is None, timeout=5000)
    assert window._baseline_snapshot is None
    assert window._baseline_diff is None
    assert window.baseline_tree._results == {}
    assert window.active_investigation_id is None


def test_pinning_a_later_snapshot_moves_the_left_tree_without_rebuilding(window, monkeypatch):
    investigation_id, a, b, c = grown_investigation(window)
    store = InvestigationStore(window.project)
    catalog = store.snapshot_catalog(investigation_id)
    assert [entry["is_baseline"] for entry in catalog] == [True, False]
    monkeypatch.setattr(ui_compare.QInputDialog, "getItem",
                        lambda *args, **kwargs: (args[3][1], True))
    before = window.project.get_setting(f"investigation.{investigation_id}")["updated_at"]
    window.choose_baseline_snapshot()
    assert store.baseline_snapshot_id(investigation_id) == catalog[1]["snapshot_id"]
    assert set(window.baseline_tree._results) == {a, b, c}
    assert window.project.get_setting(f"investigation.{investigation_id}")["updated_at"] == before
    assert window._baseline_diff["summary"]["isolates_added"] == 0


def choose_path(monkeypatch, path):
    monkeypatch.setattr(ui_compare.QFileDialog, "getSaveFileName",
                        lambda *args, **kwargs: (str(path), ""))


def test_both_trees_export_as_one_picture_carrying_the_interpretation_limit(window, monkeypatch, tmp_path):
    from PySide6.QtGui import QImage
    from PySide6.QtWidgets import QComboBox
    grown_investigation(window)
    combo = QComboBox()
    for index, suffix in ((10, "png"), (11, "jpg")):
        path = tmp_path / f"pair.{suffix}"
        choose_path(monkeypatch, path)
        window.export_graph_action(index, combo)
        assert path.stat().st_size > 0
        assert QImage(str(path)).width() == 3600
    assert window.test_errors == []


def test_the_change_summary_exports_as_a_spreadsheet_and_as_data(window, monkeypatch, tmp_path):
    import csv

    from PySide6.QtWidgets import QComboBox
    _, a, b, c = grown_investigation(window)
    combo = QComboBox()
    table = tmp_path / "changes.tsv"
    choose_path(monkeypatch, table)
    window.export_graph_action(12, combo)
    with table.open(encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("change_type")]
    added = [row for row in rows if row["change_type"] == "isolate_added"]
    assert [row["sample_id"] for row in added] == [c]
    assert any(row["change_type"] == "pairs_not_assessed" for row in rows)
    assert any(row["change_type"] == "how_to_read" and "not a phylogeny" in row["sample_id"]
               for row in rows)
    document = tmp_path / "changes.json"
    choose_path(monkeypatch, document)
    window.export_graph_action(13, combo)
    payload = json.loads(document.read_text(encoding="utf-8"))
    assert payload["format_version"] == 1 and payload["application"].startswith("WMLSTudio")
    assert all(caveat in payload["caveats"] for caveat in DIFF_CAVEATS)
    assert payload["summary"]["isolates_added"] == 1
    assert window.test_errors == []


def test_an_export_of_the_baseline_pane_writes_the_baseline_cohort_not_the_current_one(window, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QComboBox
    _, a, b, c = grown_investigation(window)
    combo = QComboBox()
    window.export_tree_choice.setCurrentIndex(window.export_tree_choice.findData("baseline"))
    profiles = tmp_path / "baseline-profiles.tsv"
    choose_path(monkeypatch, profiles)
    window.export_graph_action(8, combo)
    text = profiles.read_text(encoding="utf-8-sig")
    assert "A" in text and "B" in text and "\nC\t" not in text
    picture = tmp_path / "baseline.png"
    choose_path(monkeypatch, picture)
    window.export_graph_action(1, combo)
    assert picture.stat().st_size > 0
    distances = tmp_path / "baseline-distances.json"
    choose_path(monkeypatch, distances)
    window.export_graph_action(5, combo)
    payload = json.loads(distances.read_text(encoding="utf-8"))
    assert {row["source"] for row in payload["pairs"]} | {row["target"] for row in payload["pairs"]} == {a, b}
    assert window.test_errors == []


def test_exports_that_need_both_halves_say_which_half_is_missing(window, monkeypatch):
    from PySide6.QtWidgets import QComboBox
    imported(window, "A")
    imported(window, "B", "2111")
    window.refresh_comparison()
    window.dual_toggle.setChecked(True)
    said = []
    monkeypatch.setattr(window, "notify", lambda message: said.append(str(message)))
    monkeypatch.setattr(ui_compare.QFileDialog, "getSaveFileName",
                        lambda *args, **kwargs: pytest.fail("Nothing may be written without both trees."))
    for index in (10, 11, 12, 13):
        window.export_graph_action(index, QComboBox())
    assert said and all("no baseline tree yet" in message for message in said)


def test_the_comparison_views_offer_the_shared_right_click_actions(window):
    from PySide6.QtWidgets import QMenu
    _, a, b, c = grown_investigation(window)
    assert set(window._context_adapters) >= {"compare.cohort", "compare.groups", "compare.profiles"}
    window.tree.select_ids([a, b])
    menu = QMenu()
    dispatch = window.graph_context_entries(window.tree, menu, a)
    titles = [action.text() for action in menu.actions() if action.text()]
    assert "Add 2 isolates to comparison" in titles
    assert "Remove 2 isolates from comparison" not in titles   # not an action of this view
    assert any("Copy" in title for title in titles)
    assert dispatch
    menu.deleteLater()


def test_removing_from_the_comparison_leaves_the_isolate_and_its_evidence_in_place(window, monkeypatch):
    from wmlstudio.context_menus import Selection
    _, a, b, c = grown_investigation(window)
    said = []
    monkeypatch.setattr(window, "notify", lambda message: said.append(str(message)))
    window.context_remove_from_comparison(Selection("compare.cohort", (c,)))
    assert c not in window.cohort_ids
    assert window.project.get_setting("comparison_cohort") == sorted({a, b})
    assert {sample["id"] for sample in window.project.samples()} == {a, b, c}
    record = next(sample for sample in window.project.samples() if sample["id"] == c)
    assert record["result"]["alleles"] == {"a": "2", "b": "2", "c": "1", "d": "1"}
    assert "comparison cohort" in said[-1] and "removed" in said[-1]
    window.context_add_to_comparison(Selection("compare.cohort", (c,)))
    assert c in window.cohort_ids
    assert "added to the comparison cohort" in said[-1]


def test_a_right_click_on_a_group_row_acts_on_its_isolates_not_on_the_group_identifier(window):
    _, a, b, c = grown_investigation(window)
    group = next(g for g in window._current_snapshot["groups"] if len(g["members"]) > 1)
    adapter = window._context_adapters["compare.groups"]
    assert set(adapter.group_members({group["id"]})) == set(group["members"])
    assert group["id"] not in adapter.group_members({group["id"]})


# ---------------------------------------------------------------------------
# Separate trees for classical MLST and for cgMLST
# ---------------------------------------------------------------------------

CG_TARGETS = 1748                                 # the full Institut Pasteur BIGSdb-Lm target set
CG_SCHEME_KEY = "pasteur:lmonocytogenes-1748"


def core_profile(name, vector, targets=CG_TARGETS, key=CG_SCHEME_KEY):
    """A core-genome profile that declares its own kind and its own published key."""
    alleles = {f"LMO{index:05d}": "1" for index in range(targets)}
    for index, value in enumerate(vector):
        alleles[f"LMO{index:05d}"] = value
    return {"sample_name": name, "scheme": "Pasteur cgMLST", "scheme_digest": "cgmlst-reference",
            "status": "profile_imported", "alleles": alleles, "calls": [], "input_sha256": "a" * 64,
            "scheme_metadata": {"type": "cgmlst", "scheme_key": key}}


def core_typed(window, name, core="11", **kwargs):
    """An isolate with a core-genome profile and no classical ST at all."""
    return window.project.add_profile(name, core_profile(name, core, **kwargs),
                                      {"organism": {"genus": "Listeria", "species": "monocytogenes"}})


def both_typed(window, name, mlst="1111", core="11", ward="ICU", **kwargs):
    """One isolate carrying a classical ST and a core-genome profile, filed apart."""
    sid = imported(window, name, mlst, ward)
    window.project.set_analysis(sid, core_profile(name, core, **kwargs))
    return sid


def cohort(window, ids):
    window.cohort_ids = set(ids)
    window.refresh_cohort_table()
    window.refresh_comparison()


def test_the_two_typing_views_are_separate_trees_with_their_own_reference_and_targets(window):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    assert window.typing_kind == "mlst"
    assert set(window.tree._results) == {a, b}
    assert window._current_snapshot["typing_kind"] == "mlst"
    assert window._current_snapshot["target_loci"] == 4
    assert window._current_snapshot["scheme_digest"] == "study-reference"
    assert "Classical MLST · 4 loci" in window.current_caption.text()
    assert window.tree.scale["targets"] == 4
    window.show_typing_view("cgmlst")
    assert set(window.tree._results) == {a, b}
    assert window._current_snapshot["typing_kind"] == "cgmlst"
    assert window._current_snapshot["target_loci"] == CG_TARGETS
    assert window._current_snapshot["scheme_digest"] == "cgmlst-reference"
    assert f"cgMLST · {CG_TARGETS} targets" in window.current_caption.text()
    assert window.tree.scale["targets"] == CG_TARGETS
    assert window.test_errors == []


def test_a_threshold_set_on_one_tree_is_never_carried_into_the_other(window):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    window.cluster_threshold.setValue(1)
    assert " of 4 loci" in window.cluster_threshold.suffix()
    window.show_typing_view("cgmlst")
    window.cluster_threshold.setValue(12)
    assert window._current_snapshot["threshold"] == 12
    assert f" of {CG_TARGETS} targets" in window.cluster_threshold.suffix()
    window.show_typing_view("mlst")
    assert window.cluster_threshold.value() == 1
    assert window._current_snapshot["threshold"] == 1
    window.show_typing_view("cgmlst")
    assert window.cluster_threshold.value() == 12
    assert window.test_errors == []


def test_selecting_an_isolate_in_one_typing_tree_selects_it_in_the_other(window):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    window.counterpart_toggle.setChecked(True)
    assert set(window.counterpart_tree._results) == {a, b}
    assert window.counterpart_tree.scale["targets"] == CG_TARGETS
    window.tree.select_ids([a])
    assert window.counterpart_tree.selected_ids() == [a]
    assert "1 also have a cgMLST profile · 0 have no cgMLST profile" in window.graph_selection_label.text()
    window.counterpart_tree.select_ids([b])
    assert window.tree.selected_ids() == [b]
    assert window.test_errors == []


def test_an_isolate_with_no_core_genome_profile_is_absent_from_that_tree_not_near_in_it(window):
    a, b = both_typed(window, "A", "1111", "11"), imported(window, "B", "2111")
    cohort(window, [a, b])
    window.counterpart_toggle.setChecked(True)
    assert set(window.tree._results) == {a, b}
    assert set(window.counterpart_tree._results) == {a}
    window.tree.select_ids([a, b])
    assert "1 also have a cgMLST profile · 1 have no cgMLST profile" in window.graph_selection_label.text()
    assert window.counterpart_tree.edges == []
    assert window.test_errors == []


def test_each_tree_keeps_its_own_legend_instead_of_one_shared_key(window):
    # One allele apart on the classical loci, two apart on the core genome: the
    # same pair groups in one tree and not in the other, which is the whole point.
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "22")
    cohort(window, [a, b])
    window.counterpart_toggle.setChecked(True)
    assert list(window._legends["current"]) == ["Cluster 001"]
    assert list(window._legends["counterpart"]) == ["Unlinked isolate"]
    # The shared legend row belongs to the current and baseline trees, which are
    # one scale; the other kind's group names are never merged into it.
    assert "Cluster 001" in window.graph_legend.text()
    assert "Unlinked isolate" not in window.graph_legend.text()
    assert "Unlinked isolate" in window.counterpart_legend.text()
    assert "cgMLST groups only" in window.counterpart_legend.toolTip()
    assert window.test_errors == []


def test_a_baseline_of_another_typing_kind_refuses_the_numeric_comparison():
    records = [profile("a", "1111"), profile("b", "2111")]
    rows = pairwise_distances(records)
    classical = build_snapshot(records, rows, 1, 0.95, kind="mlst")
    core = build_snapshot(records, rows, 1, 0.95, kind="cgmlst")
    diff = snapshot_diff(classical, core)
    assert diff["policy"]["comparable"] is False
    assert diff["pairs"] is None and diff["groups"] is None and diff["mst_edges"] is None
    assert "Classical MLST comparison" in diff["policy"]["reason"]
    assert "cgMLST comparison" in diff["policy"]["reason"]
    assert any("different quantities" in caveat for caveat in diff["caveats"])
    assert diff["isolates"]["retained"] == ["a", "b"]
    assert diff["summary"]["distance_changed"] is None


def test_a_snapshot_with_no_recorded_typing_kind_is_not_refused_for_that_reason():
    records = [profile("a", "1111"), profile("b", "2111")]
    rows = pairwise_distances(records)
    older = {**build_snapshot(records, rows, 1, 0.95), "typing_kind": ""}
    current = build_snapshot(records, rows, 1, 0.95, kind="cgmlst")
    diff = snapshot_diff(older, current)
    assert diff["policy"]["comparable"] is True
    assert all(item["field"] != "typing_kind" for item in diff["policy"]["differences"])


def test_a_published_cutoff_for_this_exact_reference_is_offered_with_its_citation(window):
    a, b = core_typed(window, "A", "11"), core_typed(window, "B", "21")
    cohort(window, [a, b])
    assert window.typing_kind == "cgmlst", "a project with only core-genome profiles opens on that tree"
    suggestion = window._threshold_suggestion
    assert suggestion["status"] == "available"
    assert suggestion["auto_apply"] is False
    assert [entry["published_threshold"] for entry in suggestion["entries"]] == [7]
    assert window.cluster_threshold.value() == 1, "a published number is never applied on its own"
    banner = window.guidance_banner.text()
    assert "≤ 7" in banner and "1748 targets" in banner
    assert "10.1038/nmicrobiol.2016.185" in banner
    assert "Not applied" in banner
    assert window.guidance_row.isHidden() is False
    assert window.guidance_button.isEnabled() is True
    assert window.test_errors == []


def test_a_partial_target_set_is_refused_rather_than_scaled(window):
    a = core_typed(window, "A", "11", targets=40)
    b = core_typed(window, "B", "21", targets=40)
    cohort(window, [a, b])
    suggestion = window._threshold_suggestion
    assert suggestion["status"] == "target_count_mismatch"
    assert "no scaling" in suggestion["reason"]
    assert window.guidance_row.isHidden() is False
    assert window.guidance_button.isEnabled() is False
    assert window.cluster_threshold.value() == 1
    assert window.test_errors == []


def test_a_reference_with_no_curated_binding_reports_a_gap_not_a_number(window):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    assert window.typing_kind == "mlst"
    assert window._threshold_suggestion["status"] == "no_scheme_binding"
    assert window.guidance_row.isHidden() is True
    assert "gap in the evidence" in window.cluster_threshold.toolTip()
    assert window.test_errors == []


def test_the_crowded_graph_controls_start_behind_an_advanced_disclosure(window):
    imported(window, "A")
    assert window.advanced_panel.isHidden() is True
    assert window.advanced_toggle.text() == "Advanced ▸"
    for control in (window.dual_toggle, window.counterpart_toggle, window.graph_search,
                    window.color_by, window.export_tree_choice):
        assert control.parent() is window.advanced_panel
    assert window.cluster_threshold.isHidden() is False
    window.advanced_toggle.setChecked(True)
    assert window.advanced_panel.isHidden() is False
    assert window.advanced_toggle.text() == "Advanced ▾"
    assert window.project.get_setting("compare.advanced_open") is True


def test_each_tree_exports_under_its_own_typing_kind_and_target_count(window, monkeypatch, tmp_path):
    from pathlib import Path

    from PySide6.QtWidgets import QComboBox
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    window.counterpart_toggle.setChecked(True)
    offered = []

    def chosen(_parent, _title, suggested, _filters):
        offered.append(suggested)
        return str(tmp_path / Path(suggested).name), ""

    monkeypatch.setattr(ui_compare.QFileDialog, "getSaveFileName", chosen)
    combo = QComboBox()
    window.export_graph_action(4, combo)
    window.export_tree_choice.setCurrentIndex(window.export_tree_choice.findData("counterpart"))
    window.export_graph_action(4, combo)
    assert offered == ["mlst-comparison.nwk", "cgmlst-comparison.nwk"]
    assert "Classical MLST" in (tmp_path / "mlst-comparison.nwk").read_text(encoding="utf-8")
    core = (tmp_path / "cgmlst-comparison.nwk").read_text(encoding="utf-8")
    assert f"cgMLST · Pasteur cgMLST · {CG_TARGETS} targets" in core
    assert window.test_errors == []


def test_the_other_typing_view_groups_by_its_own_threshold_not_by_this_one(window):
    # Two allele differences on the core genome, one on the classical loci.
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "22")
    cohort(window, [a, b])
    window.show_typing_view("cgmlst")
    window.cluster_threshold.setValue(4)
    window.show_typing_view("mlst")
    window.cluster_threshold.setValue(1)
    window.counterpart_toggle.setChecked(True)
    assert window._current_snapshot["threshold"] == 1
    assert window._counterpart_snapshot["threshold"] == 4
    assert window.counterpart_tree.cluster_threshold == 4
    assert [g["status"] for g in window._counterpart_snapshot["groups"]] == ["cluster"]
    assert "link ≤ 4" in window.counterpart_caption.text()
    assert "of 4 loci" in window.cluster_threshold.suffix(), "the visible spinbox stays on this tree"
    assert window.test_errors == []


def test_opening_an_investigation_opens_the_typing_view_its_reference_belongs_to(window):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    assert window.typing_kind == "mlst"
    plan = InvestigationStore(window.project).save(
        "Core review", [a, b], scheme="Pasteur cgMLST", scheme_digest="cgmlst-reference",
        threshold=7, min_overlap=0.95, protocol="Synthetic test protocol; not a clinical cutoff")
    window.select_investigation(plan["id"])
    assert window.typing_kind == "cgmlst"
    assert window._current_snapshot["typing_kind"] == "cgmlst"
    assert window._current_snapshot["target_loci"] == CG_TARGETS
    assert window.cluster_threshold.value() == 7
    assert window.test_errors == []


def test_asking_for_a_tree_this_project_has_no_profiles_for_says_so(window):
    a, b = core_typed(window, "A", "11"), core_typed(window, "B", "21")
    cohort(window, [a, b])
    assert window.typing_kind == "cgmlst" and set(window.tree._results) == {a, b}
    window.show_typing_view("mlst")
    assert window.typing_kind == "mlst", "an explicit choice is not overridden by what happens to exist"
    assert window.tree._results == {}
    assert window._current_snapshot["profiles"] == []
    assert "Classical MLST" in window.tree_status.text()
    status = window.tree_status.text()
    assert "0 / 2 stored Classical MLST profiles" in status
    assert "of loci" not in status, "with no profiles there is no target set to count over"
    assert "link \u2264 1 allele differences" in status
    # The cohort table still shows what these isolates do have, so "not called"
    # in this view never reads as "not typed at all".
    stored = {window.cohort_table.item(row, 3).text() for row in range(window.cohort_table.rowCount())}
    assert stored == {"cgMLST: Pasteur cgMLST"}
    assert window.test_errors == []


def test_each_views_threshold_survives_leaving_and_reopening_the_project(window, tmp_path, qtbot):
    a, b = both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21")
    cohort(window, [a, b])
    window.cluster_threshold.setValue(2)
    window.show_typing_view("cgmlst")
    window.cluster_threshold.setValue(11)
    original = window.project.path
    assert window.switch_project(tmp_path / "elsewhere.wmlstudio") is True
    qtbot.waitUntil(lambda: window.comparison_worker is None, timeout=5000)
    assert window.switch_project(original) is True
    qtbot.waitUntil(lambda: window.comparison_worker is None, timeout=5000)
    assert window.typing_kind == "cgmlst"
    assert window.cluster_threshold.value() == 11
    window.show_typing_view("mlst")
    assert window.cluster_threshold.value() == 2
    assert window.test_errors == []


def test_the_cgmlst_tree_tab_really_draws_the_cgmlst_tree(window):
    """End to end on the tab the user could not draw on, with real profiles.

    "cgmlst tree i cannot draw" — because the tab was a signpost. This walks the
    route a person takes: isolates with core-genome profiles, choose the cohort,
    open the cgMLST tree tab, and check a graph is actually there, on the
    core-genome scale, with the seven-locus tab left on its own.
    """
    a, b, c = (both_typed(window, "A", "1111", "11"), both_typed(window, "B", "2111", "21"),
               both_typed(window, "C", "3111", "31"))
    cohort(window, [a, b, c])

    window.navigate("cgmlst_tree")
    assert window.pages.current_key() == "cgmlst_tree"
    assert window.typing_kind == "cgmlst"
    # A drawn graph: these are the isolates, and there are edges between them.
    assert set(window.tree._results) == {a, b, c}
    assert window._current_snapshot["typing_kind"] == "cgmlst"
    assert window._current_snapshot["target_loci"] == CG_TARGETS
    assert len(window.tree.edges) == 2, "three isolates make a two-edge spanning tree"
    assert "cgMLST" in window.current_caption.text()
    assert str(CG_TARGETS) in window.current_caption.text()
    cg_threshold = window.cluster_threshold.value()

    # The seven-locus tab is the seven-locus tab, and shares nothing with it.
    window.navigate("compare")
    assert window.typing_kind == "mlst"
    assert window._current_snapshot["typing_kind"] == "mlst"
    assert window._current_snapshot["target_loci"] == 4
    assert "Classical MLST" in window.current_caption.text()
    window.cluster_threshold.setValue(3)

    # Back again: this tab's own cutoff, not the one just set on the other scale.
    window.navigate("cgmlst_tree")
    assert window.typing_kind == "cgmlst"
    assert window.cluster_threshold.value() == cg_threshold
    assert window._current_snapshot["target_loci"] == CG_TARGETS
