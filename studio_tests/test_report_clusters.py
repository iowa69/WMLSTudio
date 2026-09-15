"""Clusters in a report, at the threshold that was actually chosen.

Every test here exists for one failure: a report that shows groups without
saying what produced them. A cutoff the user declared for their own laboratory
printed as though a journal published it, a core-genome number printed beside a
seven-locus picture, or a chain of near-threshold links printed as a set of
isolates that were all found close to each other — each of those is a document
asserting an outbreak that nobody measured.
"""

import hashlib

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.export import (
    cluster_linkage,
    cluster_section_html,
    cohort_picture_html,
    snp_section_html,
    threshold_provenance,
)
from wmlstudio.investigation import build_snapshot
from wmlstudio.simple_report import simple_report_html
from wmlstudio.threshold_guidance import record_decision

PANEL = tuple(f"locus{index:03d}" for index in range(12))


def record(sample_id, alleles, *, organism="Klebsiella pneumoniae", scheme="Ward core panel"):
    """A project record shaped exactly like Project.samples() returns one."""
    return {"id": sample_id, "name": sample_id, "status": "completed",
            "input_path": f"/local/{sample_id}.fasta", "metadata": {"organism": organism},
            "result": {"sample_name": sample_id, "status": "complete", "scheme": scheme,
                       "scheme_digest": "pinned", "st": "",
                       "input_sha256": hashlib.sha256(sample_id.encode()).hexdigest(),
                       "alleles": alleles}}


def panel(*differing, loci=PANEL):
    """A profile differing from the all-ones profile at exactly these positions."""
    return {locus: ("2" if index in differing else "1") for index, locus in enumerate(loci)}


def core(sample_id, differences=0, *, targets=2358, organism="Klebsiella pneumoniae"):
    """A core-genome-scale profile: 2,358 targets is the set the catalogue binds."""
    loci = tuple(f"target{index:05d}" for index in range(targets))
    return record(sample_id, panel(*range(differences), loci=loci), organism=organism,
                  scheme=f"{organism} cgMLST {targets}")


def snapshot_for(records, threshold=1, min_overlap=0.95):
    results = [{"sample_id": entry["id"], "sample_name": entry["name"],
                "organism": entry["metadata"]["organism"], **entry["result"]} for entry in records]
    return build_snapshot(results, pairwise_distances(results, min_overlap), threshold, min_overlap)


def klebsiella_five():
    """The user's own rule: five allele differences on the full 2,358-target scheme."""
    records = [core("kp-001"), core("kp-002", differences=3)]
    return records, snapshot_for(records, threshold=5)


def adoption(threshold=15):
    """The only route to a published number in force: every binding, all three reviews."""
    return record_decision(
        "klebsiella-pneumoniae-glasgow2025",
        {"method": "cgmlst", "organism": "Klebsiella pneumoniae",
         "scheme_key": "cgmlst.org:kpneumoniae-2358", "locus_count": 2358, "scheme_digest": "pinned",
         "min_overlap": 0.95, "missing_policy": "Pairwise ignore missing targets.",
         "caller": "WMLSTudio local caller"},
        selected_threshold=threshold, justification="Reviewed locally for this ward outbreak protocol.",
        protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)


# ---------------------------------------------------------------------------
# Where the number came from: published, or this laboratory's own rule
# ---------------------------------------------------------------------------


def test_the_operational_five_is_printed_as_the_users_own_cutoff_and_never_as_published():
    """Prevents an unsourced local rule acquiring the authority of a publication.

    Nobody published five allele differences for the 2,358-target Klebsiella
    scheme. A report that prints the clusters it produced must say whose number
    it is, or a reader will take the groups as a published cluster definition.
    """
    records, snapshot = klebsiella_five()
    block = cohort_picture_html(snapshot, None, sample_ids={entry["id"] for entry in records})

    assert "your own operational cutoff, never a published one" in block
    assert "found no publication establishing 5 allele differences" in block
    assert "cgmlst.org:kpneumoniae-2358" in block and "2358 targets" in block
    assert "at most 5 allele differences" in block
    # Not one of the words that would make the number sound published.
    assert "a published cutoff, reviewed and adopted here" not in block
    assert "published cutoff at most 5" not in block and "publishes at most 5" not in block


def test_the_operational_cutoff_is_shown_beside_the_published_number_it_departs_from():
    """Prevents a stricter local rule being read as the safer or the agreed one.

    Five is three times stricter than the only reviewed cutoff bound to this
    scheme. A reader has to see that, and has to see that stricter splits real
    chains as readily as it separates unrelated isolates.
    """
    records, snapshot = klebsiella_five()
    block = cohort_picture_html(snapshot, None, sample_ids={entry["id"] for entry in records})

    assert "at most 15 allele differences from Glasgow et al. (2025)" in block
    assert "This local cutoff is not one of them." in block
    assert "A stricter cutoff is not automatically the safer one." in block


def test_an_adopted_publication_prints_its_citation_beside_the_picture():
    """Prevents a published cutoff in force losing the citation that justifies it."""
    records = [core("kp-001"), core("kp-002", differences=3)]
    snapshot = snapshot_for(records, threshold=15)
    snapshot["threshold_evidence"] = adoption()
    block = cohort_picture_html(snapshot, None, sample_ids={entry["id"] for entry in records})

    assert "a published cutoff that was reviewed and adopted" in block
    assert "10.1128/jcm.00646-25" in block
    assert "The authors’ own caveat:" in block
    assert "your own operational cutoff" not in block


def test_a_local_setting_that_matches_no_declared_rule_says_it_has_nothing_behind_it():
    """Prevents an arbitrary number being read as either published or declared."""
    records = [core("kp-001"), core("kp-002", differences=3)]
    block = cohort_picture_html(snapshot_for(records, threshold=4), None)

    assert "your own setting for this comparison" in block
    assert "no declared operational rule matches it" in block
    assert "your own operational cutoff" not in block


# ---------------------------------------------------------------------------
# One threshold belongs to one quantity
# ---------------------------------------------------------------------------


def test_a_core_genome_cutoff_never_appears_beside_a_classical_mlst_picture():
    """Prevents the cgMLST cutoff for an organism being printed over its ST distances.

    The same organism has a core-genome catalogue entry and a declared local
    cgMLST rule. Neither may appear against a twelve-locus picture: a distance
    over twelve loci and a distance over 2,358 targets are different quantities.
    """
    records = [record("kp-001", panel()), record("kp-002", panel(0, 1, 2))]
    block = cohort_picture_html(snapshot_for(records, threshold=5), None,
                                sample_ids={entry["id"] for entry in records})

    assert "classical MLST over 12 loci" in block
    assert "cgmlst.org:kpneumoniae-2358" not in block
    assert "operational cutoff" not in block
    assert "Glasgow" not in block and "10.1128/" not in block
    assert "they share no scale, no column and no threshold" in block


def test_a_snapshot_that_calls_itself_mlst_is_offered_no_core_genome_number_either():
    """Prevents a target count overriding the typing the comparison recorded for itself.

    The catalogue refuses on the target count. This refuses again on the kind
    the snapshot declares, so a comparison that calls itself classical MLST is
    never handed a core-genome cutoff by arithmetic on its locus count.
    """
    records, snapshot = klebsiella_five()
    snapshot["typing_kind"] = "mlst"
    provenance = threshold_provenance(snapshot)

    assert provenance["operational"] == [] and provenance["operational_match"] is None
    assert provenance["operational_notice"] == ""
    assert provenance["suggestion"]["status"] == "scale_mismatch"
    assert provenance["suggestion"]["suggestion"] is None
    assert "suggests no cutoff for it" in provenance["suggestion"]["headline"]


def test_a_snp_distance_is_never_grouped_by_the_allele_cutoff():
    """Prevents SNP counts being read against the allele threshold printed above them."""
    payload = {"engine": "SKA2", "method": "reference-free split k-mers", "run_id": "ska2-cluster-test",
               "cohort": [{"sample_id": "kp-001", "sample_name": "kp-001"},
                          {"sample_id": "kp-002", "sample_name": "kp-002"}],
               "protocol": {"description": "test protocol", "minimum_shared_fraction": 0.95},
               "pairs": [{"source": "kp-001", "target": "kp-002", "distance": 4, "comparable": True,
                          "denominator_label": "1,000 split k-mers shared"}]}
    section = snp_section_html(payload, sample_ids={"kp-001", "kp-002"})

    assert "<b>4 SNPs</b>" in section
    assert "not applied to any number here" in section
    assert "A SNP distance is never grouped by an allele cutoff" in section


# ---------------------------------------------------------------------------
# The groups themselves: who is in them, how many, and what holds them together
# ---------------------------------------------------------------------------


def test_the_cluster_list_names_the_members_of_every_group_and_counts_the_groups():
    """Prevents a highlighted picture whose groups can only be guessed at by eye.

    Nobody can count groups off a forest or read which isolate sits in which.
    The report has to write them out or the picture is the only record.
    """
    records = [record("kp-001", panel()), record("kp-002", panel(0)),
               record("kp-003", panel(*range(9)))]
    block = cluster_section_html(snapshot_for(records, threshold=2),
                                 sample_ids={entry["id"] for entry in records})

    assert "1 group of two or more isolates" in block
    assert "1 isolate that linked to nothing else" in block
    assert "kp-001, kp-002" in block
    assert "2 members in the comparison cohort" in block
    assert "No other isolate came within the threshold" in block
    assert "border-left:5px solid" in block, "the groups in the picture are highlighted here too"


def test_a_group_held_together_by_a_chain_says_so_rather_than_reporting_one_close_set():
    """Prevents single linkage reading as "every isolate in this group is close".

    Two isolates six differences apart sit in one group at a cutoff of three
    because a third links them. A reader who is not told that will report a
    transmission cluster nobody measured.
    """
    records = [record("kp-001", panel()), record("kp-002", panel(0, 1, 2)),
               record("kp-003", panel(*range(6)))]
    block = cluster_section_html(snapshot_for(records, threshold=3))

    assert "the two most distant members differ by 6 allele differences" in block
    assert "more than the cutoff of 3 allele differences" in block
    assert "not because they were measured close to each other" in block


def test_a_group_held_together_only_by_links_at_the_cutoff_says_one_step_would_dissolve_it():
    """Prevents a group that exists only at this exact cutoff reading as a firm one.

    Every link that attaches kp-001 sits on the threshold itself. One step of
    the spinner and the group is gone, which is the opposite of the stability a
    reader assumes when a cluster is drawn.
    """
    records = [record("kp-001", panel()), record("kp-002", panel(0, 1, 2)),
               record("kp-003", panel(0, 1, 2))]
    snapshot = snapshot_for(records, threshold=3)
    linkage = cluster_linkage(snapshot, [entry["id"] for entry in records])
    block = cluster_section_html(snapshot)

    assert linkage["chained"] is False, "no direct distance exceeds the cutoff here"
    assert linkage["held_by_near_links"] is True
    assert "held together by links at the edge of the threshold" in block
    assert "a cutoff one step lower would not report it" in block


def test_a_group_whose_links_all_sit_well_inside_the_cutoff_says_that_plainly():
    """Prevents the chaining warning being printed where there is no chain to warn about."""
    records = [record("kp-001", panel()), record("kp-002", panel(0)), record("kp-003", panel(0, 1))]
    block = cluster_section_html(snapshot_for(records, threshold=5))

    assert "Every member is linked below the edge of the cutoff" in block
    assert "largest link inside the group is 2 allele differences" in block
    assert "held together by links at the edge of the threshold" not in block


def test_member_pairs_that_were_never_compared_are_counted_rather_than_read_as_close():
    """Prevents an unmeasured pair inside a group being read as a pair of distance zero.

    Single linkage puts two isolates in one group that share no callable locus
    at all, because a third was comparable with both. The group is real; the
    distance between those two was never measured, and must not read as nought.
    """
    whole = record("kp-001", {locus: "1" for locus in PANEL})
    first_half = record("kp-002", {locus: "1" for locus in PANEL[:6]})
    second_half = record("kp-003", {locus: "1" for locus in PANEL[6:]})
    snapshot = snapshot_for([whole, first_half, second_half], threshold=2, min_overlap=0.5)
    linkage = cluster_linkage(snapshot, ["kp-001", "kp-002", "kp-003"])
    block = cluster_section_html(snapshot)

    assert linkage["members"] == ["kp-001", "kp-002", "kp-003"], "all three landed in one group"
    assert linkage["unassessed_pairs"] == 1
    assert "1 member pair carries no accepted distance at all: unknown, never zero." in block


def test_the_cluster_section_names_the_snapshot_every_number_came_from():
    """Prevents a group that cannot be traced back to the comparison that produced it."""
    records = [record("kp-001", panel()), record("kp-002", panel(0))]
    snapshot = snapshot_for(records, threshold=2)
    block = cluster_section_html(snapshot)

    assert "Groups counted from snapshot " + snapshot["snapshot_id"] in block
    assert "reference SHA-256 pinned" in block
    assert "distance method " + snapshot["metric_version"] in block
    assert snapshot["created_at"] in block


def test_no_group_at_this_threshold_is_not_a_finding_that_the_isolates_are_unrelated():
    """Prevents an empty cluster list being read as evidence of no relationship."""
    records = [record("kp-001", panel()), record("kp-002", panel(0))]
    block = cluster_section_html(snapshot_for(records, threshold=2), sample_ids={"kp-404"})

    assert "not a finding that the isolates are unrelated" in block


# ---------------------------------------------------------------------------
# The same claims, in the plain-words summary
# ---------------------------------------------------------------------------


def test_the_plain_summary_lists_the_groups_and_labels_the_operational_cutoff_as_the_users_own():
    """Prevents the short report carrying the clusters without their provenance.

    The plain summary is the document a non-bioinformatician reads, so it is the
    one where an unsourced number is most likely to be taken as published.
    """
    records, snapshot = klebsiella_five()
    page = simple_report_html(records, selected_ids={entry["id"] for entry in records},
                              investigation=snapshot)

    assert "<h2>Which isolates group together at this threshold?</h2>" in page
    assert "kp-001, kp-002" in page
    assert "1 group of two or more isolates" in page
    assert "your own operational cutoff, never a published one" in page
    assert "found no publication establishing 5 allele differences" in page
    assert "a published cutoff, reviewed and adopted" not in page


@pytest.mark.parametrize("threshold,expected", [
    (5, "your own operational cutoff, never a published one"),
    (4, "Where this threshold comes from"),
])
def test_only_the_declared_number_itself_is_reported_as_the_declared_rule(threshold, expected):
    """Prevents a near miss on the local rule borrowing that rule's provenance."""
    records = [core("kp-001"), core("kp-002", differences=3)]
    page = simple_report_html(records, selected_ids={entry["id"] for entry in records},
                              investigation=snapshot_for(records, threshold=threshold))

    assert expected in page
    if threshold != 5:
        assert "your own operational cutoff" not in page
