import hashlib
import re

import pytest

from wmlstudio.comparison import pairwise_distances
from wmlstudio.investigation import build_snapshot
from wmlstudio.simple_report import SUSCEPTIBILITY_CAVEAT, simple_report_html
from wmlstudio.threshold_guidance import record_decision

LOCI = ("adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA")


def genome_hash(sample_id):
    return hashlib.sha256(sample_id.encode()).hexdigest()


def hydra(sample_id, *, input_sha256, genes=("blaKPC-2",)):
    return {"report_sha256": hashlib.sha256(b"report").hexdigest(), "source_sample": sample_id,
            "summary": {"amr_genes": len(genes)}, "evidence_input_sha256": input_sha256,
            "linked_input_sha256": None, "provenance": {"sha256": hashlib.sha256(b"report").hexdigest()},
            "hits": [{"gene": gene, "element_type": "AMR", "class": "CARBAPENEM", "primary": True}
                     for gene in genes]}


def record(sample_id, vector, *, evidence=None, organism="Klebsiella pneumoniae", scheme="Demo MLST",
           loci=LOCI, st="17"):
    """A project record shaped exactly like Project.samples() returns one."""
    metadata = {"organism": organism} if organism else {}
    if evidence is not None:
        metadata["hydra"] = evidence
    return {"id": sample_id, "name": sample_id, "status": "completed", "input_path": f"/local/{sample_id}.fasta",
            "metadata": metadata,
            "result": {"sample_name": sample_id, "status": "complete", "scheme": scheme,
                       "scheme_digest": "pinned", "st": st, "input_sha256": genome_hash(sample_id),
                       "alleles": dict(zip(loci, vector))}}


def core_record(sample_id, *, targets=2358, differences=0, organism="Klebsiella pneumoniae", scheme=None):
    """A core-genome-scale profile: 2,358 targets is the set Glasgow 2025 binds."""
    loci = tuple(f"locus{index:05d}" for index in range(targets))
    vector = ["2" if index < differences else "1" for index in range(targets)]
    return record(sample_id, vector, organism=organism, loci=loci, st="",
                  scheme=scheme or f"{organism} cgMLST {targets}")


def profile(record_dict):
    return {"sample_id": record_dict["id"], "sample_name": record_dict["name"], **record_dict["result"]}


def snapshot_for(records, threshold=1, min_overlap=0.95):
    results = [profile(entry) for entry in records]
    return build_snapshot(results, pairwise_distances(results, min_overlap), threshold, min_overlap)


def render(records, snapshot=None, **kwargs):
    return simple_report_html(records, selected_ids={entry["id"] for entry in records},
                              investigation=snapshot, **kwargs)


def test_simple_report_never_separates_resistance_genes_from_the_susceptibility_caveat():
    with_genes = record("with-genes", "1111111")
    with_genes["metadata"]["hydra"] = hydra("with-genes", input_sha256=genome_hash("with-genes"))
    unverified = record("unverified", "1211111")
    unverified["metadata"]["hydra"] = hydra("unverified", input_sha256=None, genes=("blaOXA-48",))
    unverified_empty = record("unverified-empty", "1112111")
    unverified_empty["metadata"]["hydra"] = hydra("unverified-empty", input_sha256=None, genes=())
    nothing_found = record("nothing-found", "1121111")
    nothing_found["metadata"]["hydra"] = hydra("nothing-found", input_sha256=genome_hash("nothing-found"), genes=())
    never_run = record("never-run", "1111211")
    records = [with_genes, unverified, unverified_empty, nothing_found, never_run]
    report = render(records, snapshot_for(records))

    assert report.count(SUSCEPTIBILITY_CAVEAT) == 2  # the top notice and the resistance table
    assert "blaKPC-2" in report and "blaOXA-48" in report
    assert "None reported (source identity not confirmed)" in report
    assert "No resistance determinants reported by the AMR database used" in report
    assert "Not assessed" in report
    # No isolate is ever labelled susceptible: the word belongs to the caveats only.
    assert not [cell for cell in re.findall(r"<td>(.*?)</td>", report) if "susceptib" in cell.casefold()]
    # Switching the section off must not switch the caveat off with it.
    without = render(records, snapshot_for(records), settings={"amr": False})
    assert SUSCEPTIBILITY_CAVEAT in without
    assert "not a statement that no resistance genes exist" in without
    assert "blaKPC-2" not in without


def test_simple_report_withholds_stale_amr_and_labels_not_assessed_as_not_negative():
    stale = record("stale", "1111111")
    stale["metadata"]["hydra"] = hydra("stale", input_sha256=hashlib.sha256(b"an older assembly").hexdigest())
    never_run = record("never-run", "1211111")
    records = [stale, never_run]
    report = render(records, snapshot_for(records))

    assert "Not shown — the saved AMR result belongs to a different sequence file" in report
    assert "blaKPC-2" not in report
    assert "Out of date" in report and "Not run" in report
    assert "“Not assessed” and “unknown” are never “susceptible”." in report


def test_simple_report_reports_uncompared_isolates_instead_of_zero_distance():
    close = record("close-a", "1111111")
    neighbour = record("close-b", "2111111")
    uncallable = record("no-calls", [None] * 7)
    absent = record("not-compared", "1111111")
    compared = [close, neighbour, uncallable]
    report = render(compared + [absent], snapshot_for(compared))

    assert "Not comparable" in report
    assert "No comparable isolate in this comparison" in report
    assert "Not in this comparison — no allele profile for this reference" in report
    assert "No shared unambiguous known alleles." in report
    assert "it is not a distance of zero" in report
    assert "<td>0</td>" not in report
    assert re.search(r"<td>1</td>", report), "the comparable pair still reports its real distance"


def test_simple_report_gives_every_equally_close_isolate_its_own_denominator():
    """Two isolates can be equally close and still have been compared over different loci."""
    focal = record("ward-A-001", "1111111")
    complete = record("ward-A-002", "2111111")                       # 7 shared loci, 1 difference
    partial = record("ward-A-003", [None, None, "1", "1", "1", "1", "2"])  # 5 shared, 1 difference
    records = [focal, complete, partial]
    report = render(records, snapshot_for(records, min_overlap=0.5))

    assert "7/7 · 5/7" in report
    assert "they were not all compared over the same loci" in report
    # Where both closest isolates share the same loci, one denominator is enough.
    assert '<td>5/7<br><span class="muted">71% of loci shared</span></td>' in report


def adoption(scheme_digest="pinned", min_overlap=0.95):
    """The only route to an approved number: every binding field, all three reviews."""
    return record_decision(
        "klebsiella-pneumoniae-glasgow2025",
        {"method": "cgmlst", "organism": "Klebsiella pneumoniae", "scheme_key": "cgmlst.org:kpneumoniae-2358",
         "locus_count": 2358, "scheme_digest": scheme_digest, "min_overlap": min_overlap,
         "missing_policy": "Pairwise ignore missing targets.", "caller": "WMLSTudio local caller"},
        selected_threshold=15, justification="Reviewed locally against this ward outbreak protocol.",
        protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)


def test_simple_report_states_the_threshold_and_its_citation_or_that_none_is_bound():
    records = [core_record("a"), core_record("b", differences=3)]
    unbound = render(records, snapshot_for(records))
    assert "not a validated rule" in unbound
    assert "10.1128/jcm.01196-22" not in unbound
    # The most recent reviewed source is offered as a suggestion, never as a rule.
    assert "Suggested, not applied: Glasgow et al. (2025)" in unbound
    assert "cgmlst.org:kpneumoniae-2358" in unbound and "10.1128/jcm.00646-25" in unbound
    assert "Selection was based on prior clustering" in unbound  # the authors' own caveat
    assert "None of this is applied to this report." in unbound
    assert "The local reference matches the scheme binding recorded with this number" in unbound

    cited = snapshot_for(records)
    cited["threshold_evidence"] = record_decision("serratia-marcescens-kampmeier2022", {"method": "cgmlst"})
    citation_only = render(records, cited)
    assert "10.1128/jcm.01196-22" in citation_only
    assert "no number was adopted from it" in citation_only
    assert "it is your own setting" in citation_only

    bound = snapshot_for(records, threshold=15)
    bound["threshold_evidence"] = adoption()
    applied = render(records, bound)
    assert "it is a published cutoff that was reviewed and adopted for this comparison" in applied
    assert "published cutoff at most 15" in applied
    assert "applied here as at most 15 allele differences" in applied
    assert "Matching an organism or a locus count is not clinical validation." in applied
    assert "The authors’ own caveat:" in applied
    assert "Selection was based on prior clustering" in applied
    assert "Suggested, not applied" not in applied, "an adopted cutoff is not re-offered as a suggestion"


def test_an_adopted_cutoff_stops_being_in_force_when_the_target_count_is_not_the_published_one():
    """2,358 published targets and 7 compared loci are not the same measurement."""
    records = [record("a", "1111111"), record("b", "2111111")]
    snapshot = snapshot_for(records, threshold=15)
    snapshot["threshold_evidence"] = adoption()
    report = render(records, snapshot)

    assert "published over 2358 targets and this comparison measured 7" in report
    assert "not in force" in report
    assert "it is a published cutoff that was reviewed and adopted" not in report
    assert "threshold set locally for this comparison" in report


def test_simple_report_says_the_published_number_does_not_carry_to_a_different_target_set():
    records = [core_record("a", targets=120), core_record("b", targets=120, differences=2)]
    report = render(records, snapshot_for(records))

    assert "Suggested, not applied: Glasgow et al. (2025)" in report
    assert "the published number was measured over 2358 targets and this comparison used 120" in report
    assert "No scaling for a different target set exists." in report


def test_simple_report_shows_that_reviewed_sources_disagree_instead_of_choosing_for_you():
    records = [core_record("a", targets=1423, organism="Enterococcus faecium"),
               core_record("b", targets=1423, organism="Enterococcus faecium", differences=4)]
    report = render(records, snapshot_for(records))

    assert "Suggested, not applied: Glasgow et al. (2025)" in report  # the most recent bindable source
    assert "Reviewed sources for Enterococcus faecium do not agree on one number" in report
    assert "at most 25 allele differences (Higgs et al. (2022))" in report
    assert "They answer different questions and are not interchangeable." in report


def test_simple_report_names_the_newest_source_even_when_it_cannot_supply_a_number():
    records = [core_record("a", targets=2270, organism="Clostridioides difficile"),
               core_record("b", targets=2270, organism="Clostridioides difficile", differences=1)]
    report = render(records, snapshot_for(records))

    assert "Suggested, not applied: Bletz et al. (2018)" in report
    assert ("The most recent reviewed source for this organism, Siddall et al. (2025), is not bound to a "
            "scheme this catalog can bind, so it cannot supply a number.") in report
    assert "doubled it to 6 as a precaution" in report


def test_simple_report_suggests_nothing_when_the_isolates_are_not_one_organism():
    records = [core_record("a", organism="Klebsiella pneumoniae", scheme="Shared core scheme"),
               core_record("b", organism="Escherichia coli", scheme="Shared core scheme", differences=2)]
    report = render(records, snapshot_for(records))

    assert "not all the same organism" in report
    assert "Suggested, not applied" not in report


def test_simple_report_shows_failed_and_missing_inputs_instead_of_a_blank_row():
    broken = record("broth-mix", [None] * 7, organism="")
    broken["missing_input"] = True
    broken["error"] = "Assembly failed: coverage too low."
    records = [broken, record("ward-A-001", "1111111")]
    report = render(records, snapshot_for(records))

    assert "The original sequence file is not available; stored results only." in report
    assert "Assembly failed: coverage too low." in report
    assert "Not assigned" in report and "organism not determined" in report


def test_simple_report_says_an_unknown_organism_is_an_evidence_gap_not_an_absent_cutoff():
    records = [core_record("a", organism="Nocardia farcinica"),
               core_record("b", organism="Nocardia farcinica", differences=2)]
    report = render(records, snapshot_for(records))
    assert "This is an evidence gap, not proof that no publications exist." in report
    assert "10.1128/" not in report, "no other organism's citation may stand in for a missing one"


def test_no_core_genome_cutoff_is_offered_beside_a_seven_locus_distance():
    """The single most important separation: two scales, never one threshold."""
    records = [record("a", "1111111"), record("b", "2111111")]
    report = render(records, snapshot_for(records))

    assert "This summary is based on classical MLST" in report
    assert "That is a sequence-type comparison, not a core-genome one." in report
    assert "so this catalog suggests no cutoff for it" in report
    assert "they share no scale, no column and no threshold" in report
    # Not one published number, scheme key or DOI from the core-genome catalog.
    assert "cgmlst.org:kpneumoniae-2358" not in report
    assert "Glasgow" not in report and "10.1128/" not in report


def test_simple_report_names_the_scheme_and_the_typing_scale_it_measured_on():
    records = [record("a", "1111111"), record("b", "2111111")]
    report = render(records, snapshot_for(records))

    assert "<b>Reference used:</b> Demo MLST · 7 loci per profile · classical MLST over 7 loci" in report
    assert "classical MLST-scale reference (7 loci)" in report
    assert "Allele differences to closest (of 7 targets)" in report

    wide = [f"locus{index:04d}" for index in range(120)]
    large = [record("a", ["1"] * 120, scheme="Demo core scheme", loci=wide, organism="Nocardia farcinica"),
             record("b", ["2"] + ["1"] * 119, scheme="Demo core scheme", loci=wide, organism="Nocardia farcinica")]
    report = render(large, snapshot_for(large))
    assert "<b>Reference used:</b> Demo core scheme · 120 loci per profile · core-genome typing over 120 targets" in report
    assert "This summary is based on core-genome typing" in report
    assert "Those distances are not sequence-type distances." in report
    assert "classical MLST-scale reference" not in report


def test_simple_report_without_a_comparison_still_lists_the_genes_and_says_the_picture_is_missing():
    records = [record("a", "1111111")]
    records[0]["metadata"]["hydra"] = hydra("a", input_sha256=genome_hash("a"))
    report = render(records)

    assert "No comparison has been built for these isolates" in report
    assert "This is not a statement that the isolates are unrelated." in report
    assert "blaKPC-2" in report
    assert "no distances and no picture" in report
    assert "<img" not in report


def test_simple_report_inlines_the_callers_picture_from_bytes_a_data_uri_or_a_path(tmp_path):
    records = [record("a", "1111111"), record("b", "2111111")]
    snapshot = snapshot_for(records)
    picture = b"\xff\xd8\xff\xe0 not a real jpeg, only bytes"

    from_bytes = render(records, snapshot, graph_png=picture, graph_mime="image/jpeg")
    assert '<img width="640" src="data:image/jpeg;base64,' in from_bytes
    assert "minimum spanning forest" in from_bytes and "not a transmission tree" in from_bytes

    path = tmp_path / "graph.jpg"
    path.write_bytes(picture)
    assert from_bytes.split('src="')[1] == render(records, snapshot, graph_png=str(path),
                                                  graph_mime="image/jpeg").split('src="')[1]
    assert '<img width="640" src="data:image/png;base64,AAAA"' in render(
        records, snapshot, graph_png="data:image/png;base64,AAAA")
    # An unreadable picture degrades to an honest line, never to a broken image.
    missing = render(records, snapshot, graph_png=str(tmp_path / "absent.jpg"))
    assert "<img" not in missing
    assert "a missing picture is not a statement about relatedness" in missing


def test_simple_report_prints_an_implicit_whole_project_scope_as_a_notice():
    records = [record("a", "1111111")]
    note = "No isolates were chosen for this report, so it covers all 1 isolates in the project."
    chosen = render(records, scope_note="This report covers the 1 isolate(s) you chose for it.")
    implicit = render(records, scope_note=note, scope_implicit=True)

    assert '<p class="notice">' + note in implicit
    assert '<p>This report covers the 1 isolate(s) you chose for it.</p>' in chosen


def test_simple_report_escapes_isolate_names_and_keeps_the_page_offline():
    hostile = record('<script>alert("x")</script>', "1111111")
    report = render([hostile])

    assert "<script>" not in report
    assert "&lt;script&gt;" in report
    assert "http://" not in report and "https://" not in report.replace(
        "https://journals.asm.org", "")  # only cited publication URLs may appear, and none are printed
    assert "src=\"data:" not in report


@pytest.mark.parametrize("chained,expected", [(True, "Joined through intermediate isolates."),
                                              (False, "isolates linked at this threshold")])
def test_simple_report_marks_single_linkage_chains_in_the_group_column(chained, expected):
    records = [record("a", "1111111"), record("b", "2111111"), record("c", "2211111")]
    if not chained:
        records = records[:2]
    report = render(records, snapshot_for(records))
    assert expected in report
    assert "two members of one group can differ by more than the threshold" in report


# ---------------------------------------------------------------------------
# The one-click route from the Reports tab, driven through a real window.
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


def build_comparison(window, ids):
    """Take the proven investigation route to a snapshot the graph can draw."""
    from wmlstudio.investigation import InvestigationStore
    window.refresh_cohort_table()
    plan = InvestigationStore(window.project).save(
        "Ward review", ids, scheme="Ward panel", scheme_digest="ward-reference", threshold=1,
        include_new=True, protocol="Synthetic test protocol; not a clinical cutoff")
    window.select_investigation(plan["id"])
    assert window._current_snapshot, window.tree_status.text()
    return plan["id"]


def test_simple_report_on_a_fresh_project_covers_every_isolate_and_prints_that_choice(window, tmp_path):
    for name in ("ward-A-001", "ward-A-002"):
        add_isolate(window, name)
    window.refresh()
    assert window.report_sample_ids() == set(), "nothing is chosen for a report on a fresh project"

    pdf = tmp_path / "summary.pdf"
    window.simple_report(pdf, build_comparison=False)
    assert pdf.read_bytes().startswith(b"%PDF")

    page = tmp_path / "summary.html"
    window.simple_report(page, build_comparison=False)
    report = page.read_text()
    assert ('<p class="notice">No isolates were chosen for this report, so it covers all 2 '
            "isolates in the project.") in report
    assert "ward-A-001" in report and "ward-A-002" in report
    assert SUSCEPTIBILITY_CAVEAT in report
    assert window.test_errors == []


def test_declining_to_build_a_comparison_still_writes_the_summary_and_moves_no_cohort(window, tmp_path):
    for name in ("ward-A-001", "ward-A-002"):
        add_isolate(window, name)
    window.refresh()
    before = set(window.cohort_ids)
    page = tmp_path / "no-comparison.html"
    window.simple_report(page, build_comparison=False)
    report = page.read_text()

    assert "No comparison has been built for these isolates" in report
    assert "This is not a statement that the isolates are unrelated." in report
    assert "<img" not in report
    assert SUSCEPTIBILITY_CAVEAT in report
    assert set(window.cohort_ids) == before
    assert window.project.get_setting("comparison_cohort", None) is None
    assert window.test_errors == []


def test_agreeing_to_build_a_comparison_makes_it_cover_exactly_the_reported_isolates(window, tmp_path):
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    window.refresh()
    page = tmp_path / "built.html"
    window.simple_report(page, build_comparison=True)

    assert set(window.cohort_ids) == set(ids), "the picture must describe the reported isolates"
    assert window.project.get_setting("comparison_cohort", None) == sorted(ids)
    assert window.cohort_origins.entry("compare")["origin"] == "Simple summary report"
    report = page.read_text()
    assert "ward-A-001" in report and SUSCEPTIBILITY_CAVEAT in report
    assert window.test_errors == []


def test_the_simple_summary_embeds_a_jpeg_while_the_detailed_presets_still_embed_png(window, tmp_path, monkeypatch):
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    build_comparison(window, ids)
    window.report_ids = set(ids)
    window._report_investigation_snapshot = window._current_snapshot
    window.refresh_report_table()

    page = tmp_path / "simple.html"
    window.simple_report(page, build_comparison=False)
    simple = page.read_text()
    assert 'src="data:image/jpeg;base64,' in simple
    assert "Allele-distance minimum spanning forest" in simple
    assert "not a phylogeny and not a transmission tree" in simple
    assert "Reference used:" in simple and "Ward panel" in simple

    # End to end: the JPEG reaches the printed PDF as a JPEG stream, not as a
    # dropped picture or a re-encoded one.
    pdf = tmp_path / "simple.pdf"
    window.simple_report(pdf, build_comparison=False)
    printed = pdf.read_bytes()
    assert printed.startswith(b"%PDF") and b"DCTDecode" in printed

    detailed = tmp_path / "detailed.html"
    monkeypatch.setattr("wmlstudio.ui_reports.QFileDialog.getSaveFileName",
                        lambda *args: (str(detailed), ""))
    window.report_preset.setCurrentIndex(window.report_preset.findData("proximity"))
    window.export_report("html")
    assert "data:image/png;base64," in detailed.read_text()
    assert window.test_errors == []


def test_a_cohort_too_large_to_compare_in_line_still_finishes_the_summary(window, qtbot, tmp_path):
    """Above 30 profiles the comparison moves to a worker; the click must still land."""
    for index in range(35):
        vector = ["1", "1", "1", "1"]
        vector[index % 4] = str(1 + index % 7)
        add_isolate(window, f"iso-{index:03d}", vector)
    window.refresh()
    page = tmp_path / "large.html"
    window.simple_report(page, build_comparison=True)
    qtbot.waitUntil(page.exists, timeout=60000)

    report = page.read_text()
    assert 'src="data:image/jpeg;base64,' in report
    assert "No comparison has been built" not in report
    assert len(window.cohort_ids) == 35
    assert window.test_errors == []


def test_the_summary_falls_back_to_the_focused_isolates_and_names_where_they_came_from(window, tmp_path):
    focused = add_isolate(window, "ward-A-001")
    add_isolate(window, "ward-B-002")
    window.refresh()
    window.focus.set_focus([focused], "Graph selection")

    scope = window.resolve_report_scope()
    assert scope["ids"] == {focused} and scope["implicit"]
    page = tmp_path / "focused.html"
    window.simple_report(page, build_comparison=False)
    report = page.read_text()

    assert ('<p class="notice">No isolates were chosen for this report, so it covers the 1 '
            "isolate(s) you had selected — from Graph selection.") in report
    assert "ward-A-001" in report
    assert "ward-B-002" not in report
    assert window.test_errors == []


def test_the_reports_page_names_the_typing_the_report_will_be_about(window, qtbot, tmp_path):
    """A cgMLST report and an MLST report must never be mistakable for one another."""
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    build_comparison(window, ids)
    window.report_ids = set(ids)
    window._report_investigation_snapshot = window._current_snapshot
    window.refresh_report_table()

    line = window.report_typing_label.text()
    assert "classical MLST over 4 loci" in line
    assert "Ward panel · 4 targets in the reference" in line
    assert "The threshold in force is a local setting, not a published cutoff" in line
    assert "they share no scale, no column and no threshold" in line
    assert window.report_typing_kind() == "mlst"
    assert window.report_filename("report", "pdf") == "wmlstudio-mlst-report.pdf"
    assert window.pdf_title_suffix() == " — classical MLST, 4 targets"

    # The same page, a core-genome comparison: a different quantity, said so in
    # the sentence, in the offered filename and in the PDF's own title.
    core = [core_record("a", targets=40), core_record("b", targets=40, differences=3)]
    window._report_investigation_snapshot = snapshot_for(core)
    window.refresh_report_table()

    assert "core-genome typing over 40 targets" in window.report_typing_label.text()
    assert window.report_typing_kind() == "cgmlst"
    assert window.report_filename("report", "pdf") == "wmlstudio-cgmlst-report.pdf"
    assert window.pdf_title_suffix() == " — cgMLST, 40 targets"
    assert window.pdf_title_suffix(full_project=True) == "", "a whole-project export states no one scale"
    assert window.test_errors == []


def test_a_report_with_no_comparison_says_so_instead_of_naming_a_typing(window):
    add_isolate(window, "ward-A-001")
    window.refresh()
    window.report_ids = {s["id"] for s in window.project.samples()}
    window.refresh_report_table()

    assert window.report_typing() is None
    assert "No comparison is attached to this report" in window.report_typing_label.text()
    assert "not a statement that the isolates are unrelated" in window.report_typing_label.text()
    assert window.report_filename("summary", "pdf") == "wmlstudio-summary.pdf"
    assert window.pdf_title_suffix() == ""
    assert window.test_errors == []


def test_the_simple_summary_offers_a_filename_that_names_its_typing(window, qtbot, monkeypatch):
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    build_comparison(window, ids)
    window.report_ids = set(ids)
    window.refresh_report_table()
    offered = []
    monkeypatch.setattr("wmlstudio.ui_reports.QFileDialog.getSaveFileName",
                        lambda *args: offered.append(args[2]) or ("", ""))
    window.simple_report(None, build_comparison=False)

    assert offered == ["wmlstudio-mlst-summary.pdf"]
    assert window.test_errors == []


def test_right_click_on_the_report_table_changes_only_this_reports_scope(window, qtbot):
    from wmlstudio.context_menus import SEPARATOR, Selection
    first = add_isolate(window, "ward-A-001")
    second = add_isolate(window, "ward-B-002")
    window.report_ids = {first, second}
    window.refresh()
    window.navigate(5)
    qtbot.waitUntil(window.report_table.isVisible, timeout=3000)

    selection = Selection("reports", (second,))
    titles = [entry.format_title(selection) for entry in window.context_menu_plan(selection)
              if entry is not SEPARATOR]
    assert "Remove from report" in titles and "Add to report" in titles
    assert "context_remove_from_report" not in window.context_menu_report()

    window.context_remove_from_report(selection)
    assert window.report_ids == {first}
    assert {s["id"] for s in window.project.samples()} == {first, second}, "no record was deleted"
    window.context_add_to_report(selection)
    assert window.report_ids == {first, second}
    assert window.test_errors == []


def test_simple_report_picture_is_never_wider_than_the_printed_page(qapp, tmp_path):
    """The one Qt check: a picture wider than its text block is silently cut off."""
    from PySide6.QtCore import QBuffer, QByteArray, QIODevice
    from PySide6.QtGui import QColor, QImage, QPageSize, QPdfWriter, QTextDocument

    picture = QImage(1800, 1200, QImage.Format.Format_RGB32)  # the size ui_reports rasterises
    picture.fill(QColor("#12222b"))
    data = QByteArray()
    buffer = QBuffer(data)
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    assert picture.save(buffer, "JPEG", 92)
    records = [record("a", "1111111"), record("b", "2111111")]
    html = render(records, snapshot_for(records), graph_png=bytes(data), graph_mime="image/jpeg")
    writer = QPdfWriter(str(tmp_path / "geometry.pdf"))
    writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
    writer.setResolution(144)
    document = QTextDocument()
    document.setHtml(html)
    # QTextDocument.print_ lays the page out in 96-dpi units, whatever the printer resolution.
    document.setTextWidth(writer.width() * 96 / writer.resolution())

    pictures = []
    block = document.begin()
    while block.isValid():
        fragment = block.begin()
        while not fragment.atEnd():
            fmt = fragment.fragment().charFormat()
            if fmt.isImageFormat():
                pictures.append((fmt.toImageFormat().width(),
                                 document.documentLayout().blockBoundingRect(block).width()))
            fragment += 1
        block = block.next()
    assert pictures, "the report should carry exactly one embedded picture"
    for width, available in pictures:
        assert 0 < width <= available, f"a {width}px picture is cut off by a {available}px page body"


def plasmid_record(sample_id, vector="1111", *, replicons=("IncFIB",), links=(("IncFIB", "blaKPC-2"),),
                   stale=False):
    """An isolate with replicon markers and a characterization that placed them."""
    digest = genome_hash(sample_id)
    evidence = hydra(sample_id, input_sha256=digest)
    evidence["hits"] += [{"gene": gene, "element_type": "PLASMID", "primary": True}
                         for gene in replicons]
    entry = record(sample_id, vector, evidence=evidence)
    entry["metadata"]["characterization"] = {
        "input_sha256": "b" * 64 if stale else digest,
        "plasmid_hypotheses": {
            "status": "completed",
            "contig_associations": [{"replicon": replicon, "marker": marker, "contig": "Contig_2",
                                     "marker_type": "AMR", "gap_bp": 100}
                                    for replicon, marker in links]}}
    return entry


def test_the_summary_prints_plasmid_markers_only_when_that_section_was_asked_for():
    records = [plasmid_record("ward-A-001"), plasmid_record("ward-A-002", "2111")]
    snapshot = snapshot_for(records)

    without = render(records, snapshot)
    assert "Plasmid markers" not in without

    page = render(records, snapshot, settings={"plasmid_hypotheses": True})
    assert "<h2>Plasmid markers</h2>" in page
    assert "IncFIB" in page and "blaKPC-2 with IncFIB on Contig_2" in page
    # The boundary is joined to the table in one string, so no option can part them.
    assert "</table><p class=\"notice\"><b>A replicon marker sitting on an assembled contig" in page
    assert "not MOB-suite" in page
    assert "same replicon name are not thereby carrying the same plasmid" in page


def test_the_summary_withholds_plasmid_co_location_that_belongs_to_another_assembly():
    """Stale characterization is named, never printed as this assembly's evidence."""
    records = [plasmid_record("ward-A-001", stale=True)]
    page = render(records, snapshot_for(records + [plasmid_record("ward-A-002", "2111")]),
                  settings={"plasmid_hypotheses": True})

    assert "IncFIB" in page, "the replicon marker itself is still input-verified"
    assert "earlier or different assembly" in page
    assert "Contig_2" not in page


def test_a_replicon_free_isolate_is_not_reported_as_having_a_chromosomal_gene():
    records = [plasmid_record("ward-A-001", replicons=(), links=())]
    page = render(records, snapshot_for(records + [plasmid_record("ward-A-002", "2111")]),
                  settings={"plasmid_hypotheses": True})

    assert "No replicon marker reported by the reference database used" in page
    assert "does not place those genes on the chromosome" in page


def test_the_summary_keeps_snp_distance_on_its_own_scale_in_its_own_section():
    """Prevents a SNP count being read against the allele threshold printed above it."""
    records = [record("ward-A-001", "1111"), record("ward-A-002", "2111")]
    page = render(records, snapshot_for(records))

    assert "three different quantities" in page
    assert "SNP distances appear only in the SNP section, on their own scale" in page
    assert "SNP distances · SKA2 split k-mers" in page
    assert "evidence about the assembled contig it was found on" in page


def test_the_reports_page_says_where_the_third_quantity_is(window, qtbot, tmp_path):
    ids = [add_isolate(window, name, vector) for name, vector in
           [("ward-A-001", "1111"), ("ward-A-002", "2111")]]
    build_comparison(window, ids)
    window.report_ids = set(ids)
    window._report_investigation_snapshot = window._current_snapshot
    window.refresh_report_table()

    line = window.report_typing_label.text()
    assert "classical MLST over 4 loci" in line
    assert "SNP distances from the SNP tree are a third quantity" in line
    assert "never comparable with the allele figures above" in line
    # The sentence says where SNP distances are, not that they are absent: the
    # report now carries them, and a page promising otherwise would be a lie.
    assert "carries them in their own section, on their own scale" in line
    assert window.test_errors == []
