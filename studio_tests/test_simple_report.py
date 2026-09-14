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


def test_simple_report_states_the_threshold_and_its_citation_or_that_none_is_bound():
    records = [record("a", "1111111"), record("b", "2111111")]
    unbound = render(records, snapshot_for(records))
    assert "not a validated rule" in unbound
    assert "10.1128/jcm.01196-22" not in unbound
    # The catalog is shown as context that exists, never as a rule that was applied.
    assert "Published cutoffs exist for Klebsiella pneumoniae" in unbound
    assert "None of them is applied to this report." in unbound

    cited = snapshot_for(records)
    cited["threshold_evidence"] = record_decision("serratia-marcescens-kampmeier2022", {"method": "cgmlst"})
    citation_only = render(records, cited)
    assert "10.1128/jcm.01196-22" in citation_only
    assert "No numeric cutoff was adopted from it." in citation_only
    assert "Published cutoffs exist for" not in citation_only

    bound = snapshot_for(records, threshold=15)
    bound["threshold_evidence"] = record_decision(
        "klebsiella-pneumoniae-glasgow2025",
        {"method": "cgmlst", "organism": "Klebsiella pneumoniae", "scheme_key": "cgmlst.org:kpneumoniae-2358",
         "locus_count": 2358, "scheme_digest": "pinned", "min_overlap": 0.95,
         "missing_policy": "Pairwise ignore missing targets.", "caller": "WMLSTudio local caller"},
        selected_threshold=15, justification="Reviewed locally against this ward outbreak protocol.",
        protocol_reviewed=True, schema_reviewed=True, epi_reviewed=True)
    applied = render(records, bound)
    assert "published cutoff at most 15; applied here as at most 15" in applied
    assert "Matching an organism or a locus count is not clinical validation." in applied


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
    records = [record("a", "1111111", organism="Nocardia farcinica"),
               record("b", "2111111", organism="Nocardia farcinica")]
    report = render(records, snapshot_for(records))
    assert "This is an evidence gap, not proof that no publications exist." in report
    assert "Published cutoffs exist for" not in report


def test_simple_report_names_the_scheme_and_never_claims_cgmlst_for_a_seven_locus_comparison():
    records = [record("a", "1111111"), record("b", "2111111")]
    report = render(records, snapshot_for(records))

    assert "<b>Reference used:</b> Demo MLST · 7 loci per profile" in report
    assert "classical MLST-scale reference (7 loci)" in report
    assert "cgmlst" not in report.casefold()

    wide = [f"locus{index:04d}" for index in range(120)]
    large = [record("a", ["1"] * 120, scheme="Demo core scheme", loci=wide),
             record("b", ["2"] + ["1"] * 119, scheme="Demo core scheme", loci=wide)]
    report = render(large, snapshot_for(large))
    assert "<b>Reference used:</b> Demo core scheme · 120 loci per profile" in report
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
