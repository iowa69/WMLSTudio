from copy import deepcopy

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.investigation import (
    InvestigationStore,
    build_snapshot,
    incremental_distances,
    profile_signature,
    proximity_rows,
)
from wmlstudio.project import Project
from wmlstudio.sequence import AnalysisCancelled


def profile(sid, vector, **extra):
    return {"sample_id": sid, "sample_name": sid, "scheme": "cgMLST", "scheme_digest": "pinned",
            "status": "complete", "alleles": dict(zip("abcd", vector)), **extra}


def snapshot(results, threshold=1, previous=None):
    return build_snapshot(results, pairwise_distances(results), threshold, 0.95, previous=previous)


def test_all_clusters_chaining_and_not_comparable_are_distinct():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211"),
               profile("d", "3333"), profile("e", "4333"), profile("unlinked", "9999"),
               profile("missing", [None] * 4), profile("mixed", "1111", status="mixed")]
    data = snapshot(records)
    groups = {tuple(g["members"]): g for g in data["groups"]}
    assert groups[("a", "b", "c")]["chained"] is True
    assert groups[("a", "b", "c")]["max_direct_distance"] == 2
    assert groups[("d", "e")]["status"] == "cluster"
    assert groups[("unlinked",)]["status"] == "singleton"
    assert groups[("missing",)]["status"] == "not_comparable"
    assert groups[("mixed",)]["status"] == "not_comparable"
    assert data["min_overlap"] == 0.95


def test_clusters_keep_identity_when_growing_and_record_split_merge():
    original = [profile("a", "1111"), profile("b", "2111"), profile("d", "3333"), profile("e", "4333")]
    first = snapshot(original)
    grown = snapshot(original + [profile("c", "2211")], previous=first)
    group = next(g for g in grown["groups"] if "a" in g["members"])
    old = next(g for g in first["groups"] if "a" in g["members"])
    assert group["id"] == old["id"] and group["name"] == old["name"]
    assert group["change"] == "grown" and group["added_members"] == ["c"]
    merged = snapshot(original, threshold=4, previous=first)
    assert merged["groups"][0]["change"] == "merged"
    assert set(merged["groups"][0]["parents"]) == {g["id"] for g in first["groups"]}
    split = snapshot(original, threshold=1, previous=merged)
    assert all(g["change"] == "split" for g in split["groups"])
    assert len({g["id"] for g in split["groups"]}) == 2


def test_incremental_next_week_and_metadata_never_change_core_distances(monkeypatch):
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    rows, cache, stats = incremental_distances(records)
    assert stats["computed_pairs"] == 3
    from wmlstudio import investigation
    original = investigation.pairwise_distances
    requested = []

    def counted(*args, **kwargs):
        requested.append(kwargs["pair_keys"])
        return original(*args, **kwargs)

    monkeypatch.setattr(investigation, "pairwise_distances", counted)
    changed = deepcopy(records)
    changed[0]["metadata"] = {"hydra": {"hits": [{"gene": "mecA", "primary": True}]}}
    changed[0]["sample_name"] = "Renamed"
    new_rows, new_cache, stats = incremental_distances(changed + [profile("d", "3333")], previous=cache)
    assert requested == [{("a", "d"), ("b", "d"), ("c", "d")}]
    assert stats["reused_pairs"] == 3 and stats["computed_pairs"] == 3
    assert new_rows[0]["source_name"] == "Renamed"
    assert new_rows[0]["distance"] == rows[0]["distance"]
    assert profile_signature(records[0]) == profile_signature(changed[0])
    again, _, stats = incremental_distances(changed + [profile("d", "3333")], previous=new_cache)
    assert again == new_rows and stats["computed_pairs"] == 0


def test_changed_profile_overlap_and_corrupt_cache_recompute():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "2211")]
    _, cache, _ = incremental_distances(records)
    changed = deepcopy(records)
    changed[0]["alleles"]["a"] = "9"
    _, _, stats = incremental_distances(changed, previous=cache)
    assert stats["computed_pairs"] == 2 and stats["reused_pairs"] == 1
    _, _, stats = incremental_distances(records, min_overlap=1, previous=cache)
    assert stats["computed_pairs"] == 3
    cache["rows"][0]["distance"] = 999
    rows, _, stats = incremental_distances(records, previous=cache)
    assert stats["computed_pairs"] == 3 and rows[0]["distance"] == 1


def test_nearest_neighbours_ties_denominators_and_selection():
    records = [profile("a", "1111"), profile("b", "2111"), profile("c", "1211"),
               profile("missing", [None] * 4)]
    rows = proximity_rows(snapshot(records), {"a"})
    assert len(rows) == 1
    assert {p["neighbour_id"] for p in rows[0]["nearest"]} == {"b", "c"}
    assert all(p["shared_loci"] == p["total_loci"] == 4 for p in rows[0]["nearest"])
    assert rows[0]["excluded_count"] == 1
    with pytest.raises(ValueError, match="outside"):
        proximity_rows(snapshot(records), {"not-present"})


def test_persisted_snapshots_reopen_and_frozen_review_membership(tmp_path):
    project = Project(tmp_path / "investigation.wmlstudio")
    sid1 = project.add_profile("a", profile("a", "1111"))
    sid2 = project.add_profile("b", profile("b", "2111"))
    results = [profile(sid1, "1111"), profile(sid2, "2111")]
    store = InvestigationStore(project)
    plan = store.save("Ward A", [sid1, sid2], scheme_digest="pinned", scheme="cgMLST",
                      protocol="Local pilot protocol; threshold 1")
    frozen = store.freeze_group(plan["id"], "Review Monday", [sid1])
    rows, cache, _ = incremental_distances(results)
    saved = store.save_snapshot(plan["id"], build_snapshot(results, rows, 1, 0.95), cache)
    with pytest.raises(ValueError, match="immutable"):
        store.save_snapshot(plan["id"], saved)
    project.close()
    with Project(tmp_path / "investigation.wmlstudio") as reopened:
        store = InvestigationStore(reopened)
        assert store.snapshot(plan["id"])["snapshot_id"] == saved["snapshot_id"]
        assert store.get(plan["id"])["review_groups"] == [frozen]
        assert incremental_distances(results, previous=store.cache(plan["id"]))[2]["computed_pairs"] == 0
        assert any(h["action"] == "investigation_snapshot" for h in reopened.history())


def test_cancel_and_incompatible_schemes_do_not_invent_groups():
    with pytest.raises(AnalysisCancelled):
        incremental_distances([profile("a", "1111")], cancelled=lambda: True)
    with pytest.raises(ValueError, match="one nonempty"):
        snapshot([profile("a", "1111"), profile("b", "1111", scheme_digest="different")])
