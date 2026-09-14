from copy import deepcopy

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.investigation import (
    InvestigationStore,
    build_snapshot,
    incremental_distances,
    profile_signature,
    proximity_rows,
    scheme_key_from,
    scheme_threshold_suggestion,
    typing_scale,
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


# ---------------------------------------------------------------------------
# Which quantity a comparison measures, and which cutoff may be cited for it
# ---------------------------------------------------------------------------


def test_a_scale_counts_every_target_not_only_the_loci_that_were_called():
    records = [profile('a', '1111'), profile('b', ['2', '1', None, '1'])]
    scale = typing_scale(records, 'cgmlst')
    assert scale['kind'] == 'cgmlst' and scale['title'] == 'cgMLST'
    assert scale['targets'] == 4, 'a locus with no call is still a target'
    assert scale['unit'] == 'allele differences over 4 targets'
    assert '4 targets' in scale['caption']
    assert 'different quantities' in scale['separation']


def test_an_unstated_typing_kind_stays_unclassified_rather_than_being_guessed():
    scale = typing_scale([profile('a', '1111')])
    assert scale['kind'] == 'unclassified'
    assert scale['title'] == 'Unclassified typing'
    assert typing_scale([], 'mlst') == dict(typing_scale([], 'mlst'), targets=0, profiles=0)


def test_a_snapshot_records_the_kind_and_the_size_of_its_target_set():
    records = [profile('a', '1111'), profile('b', '2111')]
    stored = build_snapshot(records, pairwise_distances(records), 1, 0.95, kind='cgmlst')
    assert stored['typing_kind'] == 'cgmlst' and stored['target_loci'] == 4
    assert 'cgMLST' in stored['scale_caption']
    assert 'different quantities' in stored['scale_separation']


def test_a_publication_key_is_read_only_where_a_reference_states_it():
    assert scheme_key_from({'scheme_metadata': {'scheme_key': ' pasteur:lmonocytogenes-1748 '}}) == \
        'pasteur:lmonocytogenes-1748'
    assert scheme_key_from({'threshold_scheme_key': 'cgmlst.org:cdifficile-2270'}) == \
        'cgmlst.org:cdifficile-2270'
    assert scheme_key_from({'context': {'scheme_key': 'enterobase:senterica-cgmlst-3002'}}) == \
        'enterobase:senterica-cgmlst-3002'
    # An organism name is never a binding: only an explicit key counts.
    assert scheme_key_from({'scheme': 'Listeria monocytogenes cgMLST', 'organism': 'Listeria monocytogenes'}) == ''
    assert scheme_key_from(None, 'not a mapping', {}) == ''


def test_a_bound_publication_offers_its_number_with_its_citation_and_never_applies_it():
    offered = scheme_threshold_suggestion('pasteur:lmonocytogenes-1748', targets=1748)
    assert offered['status'] == 'available' and offered['auto_apply'] is False
    assert [entry['published_threshold'] for entry in offered['entries']] == [7]
    assert offered['entries'][0]['source']['doi'] == '10.1038/nmicrobiol.2016.185'
    assert 'Not applied' in offered['reason']
    assert 'Moura' in offered['headline'] and '1748 targets' in offered['headline']


def test_a_reduced_target_set_is_refused_rather_than_scaled():
    refused = scheme_threshold_suggestion('pasteur:lmonocytogenes-1748', targets=900)
    assert refused['status'] == 'target_count_mismatch'
    assert refused['entries'] == [] and 'no scaling' in refused['reason']
    assert '1748' in refused['reason'] and '900' in refused['reason']


def test_a_classical_mlst_tree_is_never_offered_a_core_genome_cutoff():
    classical = scheme_threshold_suggestion('pasteur:lmonocytogenes-1748', method='mlst', targets=1748)
    assert classical['status'] == 'no_curated_entry' and classical['entries'] == []
    assert 'mlst' in classical['reason']


def test_an_unbound_reference_reports_the_gap_instead_of_endorsing_the_current_number():
    unbound = scheme_threshold_suggestion('', targets=2358)
    assert unbound['status'] == 'no_scheme_binding' and unbound['entries'] == []
    assert 'gap in the evidence' in unbound['reason']
    assert 'not support for the threshold now in use' in unbound['reason']
    unknown = scheme_threshold_suggestion('vendor:invented-scheme-99', targets=99)
    assert unknown['status'] == 'no_curated_entry'
    assert 'vendor:invented-scheme-99' in unknown['reason']


def test_the_installed_scheme_catalogue_is_what_binds_a_cutoff_not_a_scheme_name():
    """The pinned-scheme lookup is keyed on the installed scheme, not on an organism."""
    offered = scheme_threshold_suggestion('', targets=1748,
                                          catalog_key='pasteur:lmonocytogenes-1748')
    assert offered['status'] == 'available' and offered['auto_apply'] is False
    assert offered['catalog_key'] == 'pasteur:lmonocytogenes-1748'
    assert offered['scheme_key'] == 'pasteur:lmonocytogenes-1748'
    assert [entry['published_threshold'] for entry in offered['entries']] == [7]
    assert 'Not applied' in offered['reason']
    truncated = scheme_threshold_suggestion('', targets=1200,
                                            catalog_key='pasteur:lmonocytogenes-1748')
    assert truncated['status'] == 'target_count_mismatch' and truncated['entries'] == []
    assert 'no scaling' in truncated['reason']
    # A classical MLST tree never reaches the core-genome catalogue at all.
    classical = scheme_threshold_suggestion('', method='mlst', targets=1748,
                                            catalog_key='pasteur:lmonocytogenes-1748')
    assert classical['status'] == 'no_scheme_binding'
    unknown = scheme_threshold_suggestion('', targets=10, catalog_key='vendor:not-in-the-catalogue')
    assert unknown['status'] == 'no_scheme_binding'
