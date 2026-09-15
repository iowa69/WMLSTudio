"""The report's evidence sections: SNP, point mutations, virulence, plasmids.

Every test here exists for one failure: an assay that nobody ran rendering as an
assay that found nothing. The point-mutation section is the sharpest case, so it
gets the most tests — a genome screened without an organism catalogue has had no
mutation search at all, and a reader who sees it in an empty table will read it
as clean.
"""

import hashlib

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.export import (
    REPORT_PRESETS,
    cohort_picture_html,
    mutation_section_html,
    plasmid_section_html,
    review_report_html,
    snp_section_html,
    virulence_section_html,
)
from wmlstudio.investigation import build_snapshot
from wmlstudio.snp_tree import snp_payload

LOCI = ("adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA")


def genome_hash(sample_id):
    return hashlib.sha256(sample_id.encode()).hexdigest()


def hydra(sample_id, *, genes=("blaKPC-2",), mutations=(), replicons=(),
          point_mutation_level="dna_and_protein", point_mutations=True, organism="Escherichia coli",
          databases=("ncbi", "plasmidfinder"), release="2026-03-01.1"):
    """One isolate's linked AMR evidence, with the run's own search scope recorded.

    ``point_mutation_level`` is the engine's own answer about what could be
    looked in. "none" is the case this module cares about most: it means no
    catalogue existed for the organism, so nothing was searched.
    """
    hits = [{"gene": gene, "element_type": "AMR", "element_subtype": "AMR", "class": "CARBAPENEM",
             "subclass": "MEROPENEM", "database": "ncbi", "resolution": "COMPLETE",
             "method": "EXACTX", "primary": True} for gene in genes]
    hits += [{"gene": gene, "element_type": "AMR", "element_subtype": "POINT", "class": "QUINOLONE",
              "subclass": "CIPROFLOXACIN", "database": "ncbi", "resolution": "POINT",
              "method": "POINTX", "note": note, "identity_pct": 100.0, "primary": True}
             for gene, note in mutations]
    hits += [{"gene": gene, "element_type": "PLASMID", "database": "plasmidfinder",
              "resolution": "COMPLETE", "method": "EXACTX", "primary": True} for gene in replicons]
    return {"report_sha256": hashlib.sha256(b"report").hexdigest(), "source_sample": sample_id,
            "evidence_input_sha256": genome_hash(sample_id), "linked_input_sha256": None,
            "databases": list(databases), "hits": hits,
            "provenance": {"sha256": hashlib.sha256(b"report").hexdigest()},
            "execution_provenance": {
                "version": "hydra-test",
                "reference_release": {"release": release},
                "organism": {"resolved": organism, "requested": organism,
                             "point_mutations": point_mutations,
                             "point_mutation_level": point_mutation_level,
                             "reason": ("No point-mutation catalogue is installed for "
                                        + organism if point_mutation_level == "none" else "")},
                "virulence": {"enabled": True, "organism_curated": True}}}


def record(sample_id, vector="1111111", *, evidence=None, characterization=None,
           organism="Escherichia coli"):
    """A project record shaped exactly like Project.samples() returns one."""
    metadata = {"organism": organism}
    if evidence is not None:
        metadata["hydra"] = evidence
    if characterization is not None:
        metadata["characterization"] = characterization
    return {"id": sample_id, "name": sample_id, "status": "completed",
            "input_path": f"/local/{sample_id}.fasta", "metadata": metadata,
            "result": {"sample_name": sample_id, "status": "complete", "scheme": "Demo MLST",
                       "scheme_digest": "pinned", "st": "17", "input_sha256": genome_hash(sample_id),
                       "alleles": dict(zip(LOCI, vector))}}


def virulence_block(sample_id, *, loci=(("ybt", "detected", 1, 1), ("clb", "not_detected", 0, 2)),
                    status="completed", reason="", stale=False):
    block = {"input_sha256": "b" * 64 if stale else genome_hash(sample_id)}
    if status != "completed":
        block["virulence"] = {"status": status, "reason": reason}
        return block
    block["virulence"] = {
        "status": "completed", "reference_digest": "c" * 64,
        "assay": {"name": "WMLSTudio KpSC-associated virulence allele screen",
                  "adequate_negative_assay": True, "gene_count": 6},
        "provenance": {"engine": "NCBI BLAST+", "version": "blastn: 2.15.0+"},
        "loci": [{"locus": name, "status": call, "genes_detected": found, "genes_total": total,
                  "intact_genes_detected": found, "phenotype": "not_inferred"}
                 for name, call, found, total in loci],
        "limitations": ["Not a Kleborate-equivalent result; this BLAST screen assigns no "
                        "virulence-locus lineage."]}
    return block


def snapshot_for(records, threshold=1, min_overlap=0.95):
    results = [{"sample_id": entry["id"], "sample_name": entry["name"],
                "organism": entry["metadata"].get("organism"), **entry["result"]}
               for entry in records]
    return build_snapshot(results, pairwise_distances(results, min_overlap), threshold, min_overlap)


def snp_run(pairs, *, inputs=("ward-A-001", "ward-A-002")):
    """A SKA2 payload of the shape snp_payload returns, without running the engine."""
    held = {name: 1000 + index for index, name in enumerate(inputs)}
    result = {
        "format_version": 1, "status": "completed", "engine": "SKA2", "version": "0.5.1",
        "run_id": "ska2-report-test", "method": "reference-free-assembly-split-kmer-SNPs",
        "parameters": {"k": 31, "threads": 1, "ambiguous_bases": "excluded",
                       "minimum_shared_fraction": 0.95},
        "binary_sha256": "a" * 64,
        "inputs": [{"sample_id": name, "sample_name": name, "input_path": f"/local/{name}.fa",
                    "input_sha256": f"{index}" * 64} for index, name in enumerate(inputs)],
        "rows": list(pairs),
        "split_kmers": {"k": 31, "cohort_split_kmers": sum(held.values()), "per_sample": held},
        "alignment": {"status": "completed", "columns": 20, "rows": [], "limitations": []},
        "output_directory": "/local/run", "limitations": ["Recombination is not masked."]}
    return snp_payload(result, organism="Escherichia coli")


def snp_pair(source, target, distance, *, shared=1000, fraction=0.99, comparable=True):
    return {"source": source, "target": target, "source_name": source, "target_name": target,
            "distance": distance if comparable else None, "observed_snp_count": distance,
            "shared_split_kmers": shared, "unshared_split_kmers": 10, "shared_fraction": fraction,
            "source_split_kmers": 1000, "target_split_kmers": 1001,
            "shared_fraction_of_smaller": fraction, "comparable": comparable,
            "reason": "" if comparable else
                      "Insufficient shared unambiguous split-kmers; no accepted distance."}


# ---------------------------------------------------------------------------
# SNP distances
# ---------------------------------------------------------------------------


def test_a_report_without_a_snp_run_says_the_snp_search_never_happened():
    """Prevents an empty SNP table reading as isolates with no SNP differences."""
    section = snp_section_html(None, sample_ids={"ward-A-001", "ward-A-002"})

    assert "<h2>SNP distances · SKA2 split k-mers</h2>" in section
    assert "No SNP analysis was run for these isolates" in section
    assert "not a statement that the isolates are identical, close or unrelated" in section
    assert ";base64," not in section and "<table" not in section


def test_the_snp_section_names_its_engine_run_and_per_pair_denominator():
    """Prevents a SNP count printed without the sequence it was actually measured over."""
    payload = snp_run([snp_pair("ward-A-001", "ward-A-002", 4)])
    section = snp_section_html(payload, sample_ids={"ward-A-001", "ward-A-002"})

    assert "SKA2 0.5.1" in section and "ska2-report-test" in section
    assert "k=31" in section
    assert "<b>4 SNPs</b>" in section
    assert "1,000 split k-mers shared" in section
    # No SNP threshold is ever in force in a report.
    assert "Applied to this report:</b> none" in section


def test_a_snp_pair_that_shared_too_little_sequence_is_unknown_and_never_zero():
    """Prevents a refused SNP comparison rendering as a distance of zero differences."""
    payload = snp_run([snp_pair("ward-A-001", "ward-A-002", 7, fraction=0.10, comparable=False)])
    section = snp_section_html(payload, sample_ids={"ward-A-001", "ward-A-002"})

    assert "<b>Not comparable</b>" in section
    assert "An unknown distance is never a distance of zero." in section
    assert "<b>0 SNPs</b>" not in section


def test_an_isolate_outside_the_snp_run_is_named_rather_than_left_out():
    """Prevents a reader assuming every isolate in the report was SNP-compared."""
    payload = snp_run([snp_pair("ward-A-001", "ward-A-002", 3)])
    section = snp_section_html(payload, sample_ids={"ward-A-001", "ward-A-002", "ward-A-009"})

    assert "ward-A-009" in section
    assert "<b>Not in this SNP run</b>" in section
    assert "no SNP distance to anything here" in section


# ---------------------------------------------------------------------------
# Resistance point mutations — the sharpest not-run case in the application
# ---------------------------------------------------------------------------


def test_a_genome_with_no_organism_catalogue_is_reported_as_never_searched():
    """Prevents "we had no catalogue to look in" rendering as "no mutations found"."""
    searched = record("ward-A-001", evidence=hydra(
        "ward-A-001", mutations=(("gyrA", "gyrA_S83L (L)"),)))
    unsearched = record("ward-A-002", "2111111", organism="Proteus mirabilis", evidence=hydra(
        "ward-A-002", organism="Proteus mirabilis", point_mutation_level="none"))
    section = mutation_section_html([searched, unsearched])

    assert "Was this genome searched for point mutations at all?" in section
    assert "<b>No — not searched</b>" in section
    assert "<b>Not assayed</b>" in section
    assert "so nothing below is an absence for it" in section
    # The one isolate that was searched still reports its finding, on its own row.
    assert "gyrA" in section and "S83L" in section
    assert "were not screened for point mutations at all" in section


def test_a_cohort_nobody_searched_prints_no_empty_mutation_table():
    """Prevents a mutation section with no rows reading as a cohort free of mutations."""
    records = [record("ward-A-001", evidence=hydra("ward-A-001", point_mutation_level="none")),
               record("ward-A-002", "2111111",
                      evidence=hydra("ward-A-002", point_mutation_level="none"))]
    section = mutation_section_html(records)

    assert section.count("<b>No — not searched</b>") == 2
    assert "No catalogued point mutation was reported for any isolate that was searched." in section
    assert "it says nothing at all about the isolate(s) whose search never ran" in section


def test_a_run_told_not_to_search_for_mutations_is_not_a_negative_result():
    """Prevents an opted-out mutation search being read as a clean mutation screen."""
    declined = record("ward-A-001", evidence=hydra("ward-A-001", point_mutations=False,
                                                   point_mutation_level="dna_and_protein"))
    section = mutation_section_html([declined])

    assert "<b>No — not searched</b>" in section
    assert "<b>Not assayed</b>" in section


def test_an_isolate_with_no_amr_report_is_unknown_rather_than_clear_of_mutations():
    """Prevents an isolate nobody screened at all sharing a word with one that was."""
    section = mutation_section_html([record("ward-A-001")])

    assert "<b>No — not searched</b>" in section
    assert "<b>Not assayed</b>" in section


def test_the_mutation_section_names_the_reference_release_behind_every_row():
    """Prevents a mutation finding that cannot be traced to the snapshot that made it."""
    section = mutation_section_html([record("ward-A-001", evidence=hydra(
        "ward-A-001", mutations=(("gyrA", "gyrA_S83L (L)"),), release="2026-03-01.1"))])

    assert "Reference release: 2026-03-01.1" in section
    assert "ncbi" in section
    assert "organism catalogue: Escherichia coli" in section


# ---------------------------------------------------------------------------
# Virulence
# ---------------------------------------------------------------------------


def test_an_isolate_without_a_virulence_screen_says_not_assayed():
    """Prevents a missing virulence screen reading as a genome free of virulence loci."""
    section = virulence_section_html([record("ward-A-001")])

    assert "<b>Not assayed</b>" in section
    assert "No characterization result is attached." in section
    assert "assigns no virulence-locus lineage" in section


def test_a_virulence_screen_belonging_to_another_assembly_is_withheld_not_shown():
    """Prevents an earlier assembly's virulence result being printed as this one's."""
    entry = record("ward-A-001", characterization=virulence_block("ward-A-001", stale=True))
    section = virulence_section_html([entry])

    assert "<b>Not assayed</b>" in section
    assert "ybt" not in section


def test_a_completed_virulence_screen_tabulates_each_locus_with_its_snapshot():
    """Prevents a virulence call printed without the assay and reference behind it."""
    entry = record("ward-A-001", characterization=virulence_block("ward-A-001"))
    section = virulence_section_html([entry])

    assert "<b>Detected</b>" in section and "ybt" in section
    assert "Not detected by this assay" in section
    assert "not proof of genomic absence" in section
    assert "WMLSTudio KpSC-associated virulence allele screen" in section
    assert "Reference snapshot cccccccccccc" in section


def test_a_virulence_screen_that_was_not_selected_names_its_own_reason():
    """Prevents an assay nobody selected being indistinguishable from one that ran."""
    entry = record("ward-A-001", characterization=virulence_block(
        "ward-A-001", status="not_run", reason="Assay not selected."))
    section = virulence_section_html([entry])

    assert "<b>Not assayed</b>" in section and "Assay not selected." in section


# ---------------------------------------------------------------------------
# Plasmid evidence
# ---------------------------------------------------------------------------


def test_the_plasmid_table_keeps_markers_and_co_locations_on_separate_gates():
    """Prevents a plasmid co-location surviving when its characterization is stale."""
    entry = record("ward-A-001", evidence=hydra("ward-A-001", replicons=("IncFIB",)))
    entry["metadata"]["characterization"] = {
        "input_sha256": "b" * 64,
        "plasmid_hypotheses": {"status": "completed", "contig_associations": [
            {"replicon": "IncFIB", "marker": "blaKPC-2", "contig": "Contig_2"}]}}
    rows = _rows([entry])
    section = plasmid_section_html(rows)

    assert "IncFIB" in section
    assert "blaKPC-2 with IncFIB" not in section
    assert "not MOB-suite" in section
    assert "Reference release: 2026-03-01.1" in section


def test_a_replicon_free_isolate_is_never_reported_as_having_no_plasmid():
    """Prevents "no replicon marker" being read as "this isolate carries no plasmid"."""
    entry = record("ward-A-001", evidence=hydra("ward-A-001"))
    section = plasmid_section_html(_rows([entry]))

    assert "No replicon marker reported by the reference database used" in section
    assert "no plasmid was counted" in section


def _rows(records):
    from wmlstudio.export import _snapshot
    return _snapshot(records, {entry["id"] for entry in records})


# ---------------------------------------------------------------------------
# The cohort picture, its typing and the threshold its clusters were drawn at
# ---------------------------------------------------------------------------


def test_the_cohort_picture_states_the_scheme_targets_and_threshold_beside_it():
    """Prevents a tree on a page with no sentence naming what it was drawn on."""
    records = [record("ward-A-001"), record("ward-A-002", "2111111")]
    block = cohort_picture_html(snapshot_for(records), b"\x89PNG pretend",
                                sample_ids={"ward-A-001", "ward-A-002"})

    assert ";base64," in block
    assert "Demo MLST" in block and "7 targets in the scheme" in block
    assert "classical MLST over 7 loci" in block
    assert "Clusters are highlighted at the threshold in force:</b> at most 1 allele difference" in block
    assert "your own setting, not a published cutoff" in block


@pytest.mark.parametrize("drawn", ["snp", "cgmlst"])
def test_a_tree_of_another_typing_is_refused_rather_than_printed(drawn):
    """Prevents an MLST forest being filed as the cgMLST picture of the same cohort."""
    records = [record("ward-A-001"), record("ward-A-002", "2111111")]
    block = cohort_picture_html(snapshot_for(records), b"\x89PNG pretend", graph_typing=drawn)

    assert ";base64," not in block
    assert "<b>No picture is shown.</b>" in block
    assert "never printed in a report about the other" in block


def test_a_view_that_recorded_no_typing_does_not_suppress_a_picture_the_snapshot_names():
    """Prevents an undeclared typing kind hiding a tree the report can describe in full."""
    records = [record("ward-A-001"), record("ward-A-002", "2111111")]
    block = cohort_picture_html(snapshot_for(records), b"\x89PNG pretend",
                                graph_typing="unclassified")

    assert ";base64," in block
    assert "No picture is shown." not in block
    assert "classical MLST over 7 loci" in block


def test_the_picture_block_names_the_clusters_that_exist_at_this_threshold():
    """Prevents a highlighted picture whose groups are nowhere stated in words."""
    records = [record("ward-A-001"), record("ward-A-002", "1111111")]
    snapshot = snapshot_for(records, threshold=1)
    block = cohort_picture_html(snapshot, None, sample_ids={"ward-A-001", "ward-A-002"})

    assert "<h3>Clusters at this threshold</h3>" in block
    assert "Isolates in this report" in block


def test_an_isolate_in_no_cluster_says_so_instead_of_showing_a_blank_row():
    """Prevents a singleton rendering as a cluster with an unstated, empty distance."""
    records = [record("ward-A-001"), record("ward-A-002", "2222222")]
    block = cohort_picture_html(snapshot_for(records, threshold=1), None,
                                sample_ids={"ward-A-001", "ward-A-002"})

    assert "No other isolate came within the threshold" in block
    assert "Not measured" in block
    assert "No pair inside this group carries an accepted distance." in block


def test_no_cluster_at_this_threshold_is_not_a_finding_that_isolates_are_unrelated():
    """Prevents an empty cluster table being read as evidence of no relationship."""
    records = [record("ward-A-001"), record("ward-A-002", "2222222")]
    block = cohort_picture_html(snapshot_for(records, threshold=1), None,
                                sample_ids={"ward-A-404"})

    assert "not a finding that the isolates are unrelated" in block


# ---------------------------------------------------------------------------
# The whole document
# ---------------------------------------------------------------------------


def test_the_cohort_report_carries_every_evidence_section_as_a_table():
    """Prevents a section quietly disappearing from the document the presets promise."""
    entry = record("ward-A-001", evidence=hydra("ward-A-001", replicons=("IncFIB",),
                                                mutations=(("gyrA", "gyrA_S83L (L)"),)),
                   characterization=virulence_block("ward-A-001"))
    other = record("ward-A-002", "2111111", evidence=hydra("ward-A-002"))
    records = [entry, other]
    page = review_report_html(records, selected_ids={"ward-A-001", "ward-A-002"},
                              investigation=snapshot_for(records),
                              options=REPORT_PRESETS["cohort"],
                              snp=snp_run([snp_pair("ward-A-001", "ward-A-002", 6)]))

    for heading in ("<h2>SNP distances · SKA2 split k-mers</h2>",
                    "<h2>Resistance point mutations</h2>",
                    "<h2>Virulence factors</h2>",
                    "<h2>Plasmid evidence · replicon markers on assembled contigs</h2>",
                    "<h2>Comparison-cohort picture</h2>"):
        assert heading in page, heading
    assert "<b>6 SNPs</b>" in page
    assert page.count("<table") >= 5


def test_one_isolate_report_and_one_cohort_report_carry_the_same_evidence_gates():
    """Prevents a per-isolate document silently dropping a section the cohort one keeps."""
    entry = record("ward-A-001", evidence=hydra("ward-A-001", point_mutation_level="none"))
    single = review_report_html([entry], selected_ids={"ward-A-001"},
                                options=REPORT_PRESETS["isolate"])

    assert "<h2>Resistance point mutations</h2>" in single
    assert "<b>No — not searched</b>" in single
    assert "<h2>SNP distances · SKA2 split k-mers</h2>" in single
    assert "No SNP analysis was run for these isolates" in single


def test_a_section_switched_off_is_still_not_a_negative_result():
    """Prevents an unticked section being mistaken for an assay that found nothing."""
    entry = record("ward-A-001", evidence=hydra("ward-A-001"))
    page = review_report_html([entry], selected_ids={"ward-A-001"},
                              options={**REPORT_PRESETS["cohort"], "snp": False,
                                       "point_mutations": False})

    assert "<h2>SNP distances" not in page
    assert "<h2>Resistance point mutations</h2>" not in page
    assert "an assay that was not run must never be interpreted as an absent finding" in page


# ---------------------------------------------------------------------------
# The two document shapes, driven through a real window
# ---------------------------------------------------------------------------


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "workspace")
    errors = []
    monkeypatch.setattr(widget, "error", lambda message: errors.append(str(message)))
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    qtbot.wait(40)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


def add_isolate(window, name, vector="1111"):
    return window.project.add_profile(name, {
        "sample_name": name, "scheme": "Ward panel", "scheme_digest": "ward-reference",
        "status": "profile_imported", "alleles": dict(zip("abcd", vector)), "calls": [],
        "st": "20", "input_sha256": "a" * 64,
    }, {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}})


def test_one_pdf_per_isolate_writes_one_document_for_each_chosen_isolate(window, tmp_path):
    """Prevents the per-isolate route silently folding several isolates into one file."""
    ids = [add_isolate(window, name) for name in ("ward-A-001", "ward-A-002", "ward-A-003")]
    window.report_ids = set(ids)
    window.refresh_report_table()

    folder = tmp_path / "per-isolate"
    folder.mkdir()
    written = window.export_report_per_isolate(str(folder))

    assert len(written) == 3
    assert {path.name for path in written} == {
        "wmlstudio-ward-A-001-report.pdf", "wmlstudio-ward-A-002-report.pdf",
        "wmlstudio-ward-A-003-report.pdf"}
    for path in written:
        assert path.read_bytes().startswith(b"%PDF")
    assert window.test_errors == []


def test_a_per_isolate_document_says_it_covers_one_isolate_out_of_the_selection(window, tmp_path):
    """Prevents a one-isolate PDF being read as the whole cohort's evidence."""
    ids = [add_isolate(window, name) for name in ("ward-A-001", "ward-A-002")]
    window.report_ids = set(ids)
    window.refresh_report_table()
    captured = {}
    original = window.write_pdf_report

    def spy(path, selected_ids=None, **kwargs):
        captured[path] = kwargs.get("scope_note", "")
        return original(path, selected_ids, **kwargs)

    window.write_pdf_report = spy
    folder = tmp_path / "per-isolate"
    folder.mkdir()
    window.export_report_per_isolate(str(folder))

    notes = list(captured.values())
    assert len(notes) == 2
    for note in notes:
        assert "covers one isolate" in note
        assert "out of the 2 isolate(s) in this report selection" in note
        assert "not summarised here" in note


def test_the_typing_guard_still_lets_a_real_comparison_print_its_own_tree(window, tmp_path):
    """Prevents the cgMLST/MLST picture guard suppressing the very tree it is about."""
    from wmlstudio.investigation import InvestigationStore
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    window.refresh_cohort_table()
    plan = InvestigationStore(window.project).save(
        "Ward review", ids, scheme="Ward panel", scheme_digest="ward-reference", threshold=1,
        include_new=True, protocol="Synthetic test protocol; not a clinical cutoff")
    window.select_investigation(plan["id"])
    assert window._current_snapshot
    window.report_ids = set(ids)
    window._report_investigation_snapshot = window._current_snapshot
    window.refresh_report_table()

    # The drawing view's own word for what it drew, which is what the guard reads.
    assert window.report_graph_typing() in {"mlst", "unclassified", ""}

    page = tmp_path / "cohort.html"
    from wmlstudio.export import write_review_report
    write_review_report(window.report_records(), page, selected_ids=set(ids),
                        investigation=window.report_context(), options=window.report_options(),
                        graph_png=window.report_graph_image(set(ids)),
                        graph_typing=window.report_graph_typing(),
                        snp=window.report_snp_payload())
    document = page.read_text()

    assert ";base64," in document, "the report refused its own comparison's picture"
    assert "No picture is shown." not in document
    assert "Ward panel" in document and "4 targets in the scheme" in document
    assert "Clusters are highlighted at the threshold in force" in document
    assert "No SNP analysis was run for these isolates" in document
    assert window.test_errors == []


def test_the_report_reads_the_snp_run_from_the_snp_page_and_nowhere_else(window):
    """Prevents a report inventing SNP distances, or reusing another cohort's run."""
    add_isolate(window, "ward-A-001")

    assert window.report_snp_payload() is None

    page = getattr(window, "snp_tree_page", None)
    if page is None:
        pytest.skip("this build did not mount the SNP tab")
    page.payload = snp_run([snp_pair("ward-A-001", "ward-A-002", 5)])
    assert window.report_snp_payload() is page.payload
    # An empty payload is no payload: the section must say the run never happened.
    page.payload = {}
    assert window.report_snp_payload() is None
