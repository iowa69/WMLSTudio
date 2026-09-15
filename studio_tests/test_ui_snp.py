"""The SNP tree tab: its own cohort, its own denominators, its own words, its own refusals.

The synthetic genomes here are written into the test's own directory and differ at
positions this file chooses, so the distances the native engine measures can be
counted by hand. The fabricated SKA2 results have the exact shape
:func:`wmlstudio.ska_runtime.run_ska` returns, so the page is exercised end to end
without a subprocess on machines where the engine is not staged.
"""

import json
import random
from xml.etree import ElementTree

import pytest

from wmlstudio.app import MainWindow
from wmlstudio.graph_window import GraphWindow
from wmlstudio.sequence import file_sha256
from wmlstudio.ska_runtime import runtime_capabilities
from wmlstudio.snp_tree import NO_LINK_THRESHOLD
from wmlstudio.ui_snp import SnpTreePanel

LOCI = ["gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB"]
SITES = [400 + 300 * index for index in range(15)]
SWAP = {"A": "C", "C": "G", "G": "T", "T": "A"}


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    yield widget
    widget.cancel_comparison_for_close()
    qtbot.waitUntil(lambda: widget.comparison_worker is None, timeout=5000)
    widget.close()


@pytest.fixture
def panel(window, qtbot):
    page = SnpTreePanel(window)
    qtbot.addWidget(page)
    page.resize(1100, 760)
    return page


def genomes():
    """One 20 kb sequence, two isolates differing from it at chosen sites, one fragment."""
    base = "".join(random.Random(11).choices("ACGT", k=20000))

    def mutate(positions):
        bases = list(base)
        for position in positions:
            bases[position] = SWAP[bases[position]]
        return "".join(bases)

    # iso-d is the first 8 kb only: a fragmentary assembly that differs nowhere and
    # must still never be drawn close to anything.
    return {"iso-a": base, "iso-b": mutate(SITES[:3]), "iso-c": mutate(SITES[3:]),
            "iso-d": base[:8000]}


def assembled(window, tmp_path, name, sequence=">contig1\nACGTACGTAC\n", *, organism="Enterococcus",
              species="faecium", st=None):
    """An isolate whose current input is a reviewed assembly: the case this tab compares."""
    path = tmp_path / f"{name}.fasta"
    path.write_text(sequence if sequence.startswith(">") else f">contig1\n{sequence}\n",
                    encoding="ascii")
    sample_id = window.project.add_sample(path, name, sample_id=name)
    result = {"kind": "fasta", "input_sha256": file_sha256(path), "scheme": "efaecium 7-locus",
              "scheme_digest": "d" * 64, "st": st,
              "alleles": {locus: "1" for locus in LOCI},
              "calls": [{"locus": locus, "allele": "1", "status": "exact"} for locus in LOCI]}
    window.project.set_result(sample_id, result)
    window.project.set_status(sample_id, "completed")
    window.project.update_metadata(sample_id, {"organism": {"genus": organism, "species": species}})
    return sample_id


def fake_run(rows, *, inputs=("iso-a", "iso-b"), k=31, floor=0.95, output="/read-only/run"):
    """A SKA2 result of the exact shape run_ska returns, without running anything."""
    held = {name: 20000 - 500 * index for index, name in enumerate(inputs)}
    return {"format_version": 1, "status": "completed", "engine": "SKA2", "version": "0.5.1",
            "run_id": "ska2-test", "method": "reference-free-assembly-split-kmer-SNPs",
            "parameters": {"k": k, "threads": 1, "ambiguous_bases": "excluded",
                           "min_frequency": 0.0, "minimum_shared_fraction": floor},
            "binary_sha256": "a" * 64,
            "inputs": [{"sample_id": name, "sample_name": name.upper(),
                        "input_path": f"/read-only/{name}.fa", "input_sha256": f"{index}" * 64}
                       for index, name in enumerate(inputs)],
            "rows": list(rows),
            "split_kmers": {"k": k, "cohort_split_kmers": sum(held.values()), "per_sample": held},
            "alignment": {"status": "completed", "columns": 20, "rows": [],
                          "limitations": ["Cohort alignment columns depend on this cohort."]},
            "mapped": None, "output_directory": output, "threshold": None,
            "limitations": ["SNP similarity is not direct transmission or direction of spread."]}


def pair(source, target, distance, *, shared=19800, fraction=0.99, comparable=True, observed=None,
         smaller=None):
    return {"source": source, "target": target, "source_name": source.upper(),
            "target_name": target.upper(),
            "distance": distance if comparable else None,
            "observed_snp_count": distance if observed is None else observed,
            "shared_split_kmers": shared, "unshared_split_kmers": 200,
            "shared_fraction": fraction, "source_split_kmers": 20000, "target_split_kmers": 19500,
            "shared_fraction_of_smaller": fraction if smaller is None else smaller,
            "comparable": comparable,
            "reason": "" if comparable else
                      "Insufficient shared unambiguous split-kmers; no accepted distance."}


def measurement_words(panel):
    """Everything this page says about the quantity it is showing.

    Deliberately not the whole page: the separation statement and the limitations
    name classical MLST and cgMLST allele differences on purpose, to say that this
    is not those. What must never borrow an allele word is the description of this
    measurement itself.
    """
    tables = [panel.pairs_table.horizontalHeaderItem(column).text() for column in range(6)]
    for table in (panel.pairs_table, panel.matrix_table):
        for row in range(table.rowCount()):
            for column in range(table.columnCount()):
                item = table.item(row, column)
                tables.extend([item.text(), item.toolTip()] if item is not None else [])
    return "\n".join([panel.engine.text(), panel.matrix_note.text(), panel.tree_caption.text(),
                      panel.status.text(), panel.tree.scale_caption(), *tables])


def test_the_cohort_table_says_which_isolates_can_be_compared_and_why_the_others_cannot(
        panel, window, tmp_path):
    first = assembled(window, tmp_path, "iso-a")
    second = assembled(window, tmp_path, "iso-b")
    reads = window.project.add_sample(_fastq(tmp_path, "reads"), "reads-only")
    mate = window.project.add_sample(_fastq(tmp_path, "mate"), "read-mate")
    window.project.update_metadata(mate, {"workflow": {"source_kind": "read_mate",
                                                       "paired_with": reads}})
    unverified = window.project.add_sample(tmp_path / "iso-a.fasta", "never analysed")
    panel.refresh_cohort()
    rows = {row["id"]: row for row in panel.cohort}
    assert sorted(rows) == sorted([first, second, reads, mate, unverified])
    assert rows[first]["eligible"] is True and rows[second]["eligible"] is True
    assert rows[reads]["eligible"] is False and "Assemble this isolate first" in rows[reads]["reason"]
    assert rows[mate]["eligible"] is False and "second mate" in rows[mate]["reason"]
    assert rows[unverified]["eligible"] is False
    assert "reviewed input fingerprint" in rows[unverified]["reason"]
    reasons = {panel.cohort_table.item(row, 0).text(): panel.cohort_table.item(row, 3).text()
               for row in range(panel.cohort_table.rowCount())}
    assert reasons["iso-a"] == "yes"
    assert reasons["read-mate"].startswith("no — ")
    # Two comparable isolates are enough to offer the run, and the others are
    # listed rather than silently dropped.
    assert panel.run_button.isEnabled() is runtime_capabilities()["available"]
    assert "2 of 5 isolate(s)" in panel.status.text()
    assert window.test_errors == []


def test_a_measured_cohort_is_drawn_in_its_own_words_and_borrows_none_from_allele_typing(
        panel, window, tmp_path):
    assembled(window, tmp_path, "iso-a")
    assembled(window, tmp_path, "iso-b")
    panel.refresh_cohort()
    payload = panel.show_result(fake_run([pair("iso-a", "iso-b", 3, shared=19880, fraction=0.9906)]))
    scale = panel.tree.scale
    assert scale["kind"] == "snp" and scale["difference_word"] == "SNPs"
    assert scale["target_word"] == "shared split k-mers"
    assert "SKA2 split k-mer SNPs" in panel.tree.scale_caption()
    tooltip = panel.tree._edge_tooltip(payload["graph"]["edges"][0])
    assert tooltip.startswith("3 SNPs / 19,880 split k-mers shared")
    assert "99.1% of the pair's combined set" in tooltip
    assert "not evolutionary time or transmission" in tooltip
    # Nothing describing this measurement calls a SNP an allele.
    assert "allele" not in measurement_words(panel).casefold()
    assert "allele" not in tooltip.casefold()
    # The separation statement is the one place another quantity is named, and it
    # is named to be excluded.
    assert "share no scale, no axis, no column and no threshold" in panel.tree._view_tooltip()
    assert panel.matrix_table.item(0, 1).text() == "3"
    assert panel.matrix_table.item(0, 0).text() == "0"
    assert "not compared with itself" in panel.matrix_table.item(0, 0).toolTip()
    row = {panel.pairs_table.horizontalHeaderItem(column).text():
           panel.pairs_table.item(0, column).text() for column in range(6)}
    assert row["SNP distance"] == "3" and row["Status"] == "compared"
    assert "19,880 split k-mers shared" in row["Measured over"]
    assert window.test_errors == []


def test_a_pair_that_shared_too_little_sequence_gets_no_edge_and_never_a_zero(
        panel, window, tmp_path):
    for name in ("iso-a", "iso-b", "iso-d"):
        assembled(window, tmp_path, name)
    panel.refresh_cohort()
    rows = [pair("iso-a", "iso-b", 3),
            pair("iso-a", "iso-d", None, shared=7900, fraction=0.39, comparable=False, observed=0),
            pair("iso-b", "iso-d", None, shared=7900, fraction=0.39, comparable=False, observed=3)]
    payload = panel.show_result(fake_run(rows, inputs=("iso-a", "iso-b", "iso-d")))
    assert [(edge["source"], edge["target"]) for edge in payload["graph"]["edges"]] == [
        ("iso-a", "iso-b")]
    position = {row["sample_id"]: index for index, row in enumerate(payload["matrix"]["samples"])}
    refused = panel.matrix_table.item(position["iso-a"], position["iso-d"])
    assert refused.text() == "not comparable"
    assert "unknown, not zero" in refused.toolTip()
    texts = {(panel.pairs_table.item(row, 0).text(), panel.pairs_table.item(row, 1).text()):
             (panel.pairs_table.item(row, 2).text(), panel.pairs_table.item(row, 3).text())
             for row in range(panel.pairs_table.rowCount())}
    # The zero SKA2 actually saw is kept and is visibly not a distance.
    assert texts[("ISO-A", "ISO-D")] == ("no accepted distance", "0")
    assert texts[("ISO-A", "ISO-B")] == ("3", "3")
    assert "never means zero differences" in panel.matrix_legend.text()
    assert "2 of 3 pairs fell below" in panel.comparability.text()
    assert "2 pair(s) too little compared to be drawn" in panel.tree_caption.text()
    assert panel.cohort_table.rowCount() == 3
    assert window.test_errors == []


def test_a_published_snp_cutoff_is_refused_for_a_protocol_it_was_not_measured_on(
        panel, window, tmp_path):
    assembled(window, tmp_path, "iso-a")
    assembled(window, tmp_path, "iso-b")
    panel.refresh_cohort()
    assert panel.cohort_organism() == "Enterococcus faecium"
    payload = panel.show_result(fake_run([pair("iso-a", "iso-b", 3)]))
    assert payload["threshold"]["status"] == "refused_protocol_mismatch"
    assert payload["threshold"]["threshold"] is None and payload["threshold"]["applied"] is False
    said = panel.threshold_note.text()
    assert "published SKA, this run SKA2" in said and "published 15, this run 31" in said
    assert "no number is offered and none was applied" in said
    # Nothing is grouped until a reader chooses a number, and the picture says so
    # rather than outlining every isolate as a group of one.
    assert payload["graph"]["cluster_threshold"] == NO_LINK_THRESHOLD == -1
    assert panel.tree._halos == []
    assert [group["name"] for group in panel.tree.groups()] == ["Not grouped", "Not grouped"]
    panel.link.setValue(5)
    panel.apply_settings()
    assert panel.payload["threshold"]["link_threshold"] == 5
    assert panel.payload["threshold"]["link_threshold_source"] == "reader_selected_unvalidated"
    assert panel.payload["threshold"]["threshold"] is None
    assert "a choice made in this view, not a validated cutoff" in panel.threshold_note.text()
    assert window.test_errors == []


def test_the_comparability_floor_is_re_read_from_the_run_without_measuring_anything_again(
        panel, window, tmp_path):
    """Real E. faecium assemblies share 0.53-0.78 of their combined split k-mers."""
    for name in ("iso-a", "iso-b", "iso-c"):
        assembled(window, tmp_path, name)
    panel.refresh_cohort()
    rows = [pair("iso-a", "iso-b", None, shared=2344102, fraction=0.662, comparable=False,
                 observed=2356, smaller=0.827),
            pair("iso-a", "iso-c", None, shared=1928027, fraction=0.567, comparable=False,
                 observed=7510, smaller=0.772),
            pair("iso-b", "iso-c", None, shared=1911287, fraction=0.526, comparable=False,
                 observed=7922, smaller=0.766)]
    panel.show_result(fake_run(rows, inputs=("iso-a", "iso-b", "iso-c")))
    assert panel.payload["graph"]["edges"] == []
    assert "no tree was drawn" in panel.comparability.text()
    assert "large accessory genome" in panel.comparability.text()
    panel.floor.setValue(0.5)
    panel.apply_settings()
    assert panel.payload["run_id"] == "ska2-test", "the same run, re-read"
    assert len(panel.payload["graph"]["edges"]) == 2
    assert panel.payload["comparability_reread"]["from"] == 0.95
    assert panel.payload["comparability_reread"]["to"] == 0.5
    assert "already recorded by this run" in panel.payload["comparability_reread"]["source"]
    # Not one SNP count moved; only which pairs are allowed to carry one.
    assert [row["observed_snp_count"] for row in panel.payload["pairs"]] == [2356, 7510, 7922]
    assert "does not make its distance better evidence" in panel.status.text()
    # The other denominator answers a different question, and the answer says which.
    panel.basis.setCurrentIndex(panel.basis.findData("smaller"))
    assert "smaller isolate" in panel.basis_note.text()
    panel.floor.setValue(0.8)
    panel.apply_settings()
    assert panel.payload["comparability"]["basis"] == "smaller"
    assert panel.payload["comparability"]["accepted_pairs"] == 1
    assert window.test_errors == []


def test_the_forest_opens_in_its_own_window_stating_the_method_and_never_alleles(
        panel, window, tmp_path, qtbot):
    assembled(window, tmp_path, "iso-a", st="80")
    assembled(window, tmp_path, "iso-b", st="80")
    panel.refresh_cohort()
    panel.show_result(fake_run([pair("iso-a", "iso-b", 3, shared=19880, fraction=0.9906)]))
    detached = panel.open_window()
    qtbot.addWidget(detached)
    assert isinstance(detached, GraphWindow)
    assert detached.identity.kind == "snp"
    assert "SKA2 split k-mer SNPs forest" in detached.windowTitle()
    assert "no link threshold set" in detached.windowTitle()
    assert "Every isolate in this project" in detached.windowTitle()
    assert "ska2:assembly-split-kmer-k31" in detached.windowTitle()
    assert detached.identity.caption() == panel.tree.scale_caption()
    assert "target count not recorded" not in detached.identity.caption()
    assert detached.toggles["edges"].text() == "SNPs on edges"
    assert "layout of SNPs" in detached.honesty.text()
    assert "never a zero distance" in detached.honesty.text()
    assert "own denominator" in detached.honesty.text()
    assert "allele" not in detached.honesty.text().casefold()
    assert "allele" not in detached.identity.export_subtitle().casefold()
    assert "allele" not in detached.identity.export_title().casefold()
    # Arranging the window leaves the page's own forest exactly where it was.
    before = (panel.tree.nodes["iso-a"].pos().x(), panel.tree.nodes["iso-a"].pos().y())
    detached.view.nodes["iso-a"].setPos(321, 123)
    assert (panel.tree.nodes["iso-a"].pos().x(), panel.tree.nodes["iso-a"].pos().y()) == before
    picture = tmp_path / "snp-forest.svg"
    detached.export(picture)
    written = picture.read_text(encoding="utf-8")
    assert "SNP-distance minimum spanning forest" in written
    assert "Every isolate in this project" in written
    assert "no link threshold set" in written and "single-link threshold: none set" in written
    assert "not a phylogeny or transmission chain" in written
    assert "allele" not in written.casefold()
    # The page's own export names the quantity too, without a window's identity.
    panel.tree.save_svg(tmp_path / "page.svg")
    assert "edge labels are SNPs" in (tmp_path / "page.svg").read_text(encoding="utf-8")
    graph = tmp_path / "snp-forest.graphml"
    detached.export(graph)
    document = ElementTree.parse(graph)
    namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
    keys = [key.get("id") for key in document.findall(".//g:key", namespace)]
    assert "allele_differences" not in keys and "shared_loci" not in keys
    edge = {data.get("key"): data.text
            for data in document.findall(".//g:edge/g:data", namespace)}
    assert edge == {"distance": "3", "unit": "SNPs",
                    "shared_denominator": "19,880 split k-mers shared, 99.1% of the pair's "
                                          "combined set, 99.1% of the smaller isolate's"}
    detached.close()
    assert panel._windows == []
    assert window.test_errors == []


def test_clear_empties_this_page_and_keeps_every_sample_result_and_written_output(
        panel, window, tmp_path):
    first = assembled(window, tmp_path, "iso-a")
    assembled(window, tmp_path, "iso-b")
    panel.refresh_cohort()
    panel.show_result(fake_run([pair("iso-a", "iso-b", 3)]))
    assert panel.tree.nodes and panel.window_button.isEnabled()
    before = window.project.get_sample(first)
    panel.clear()
    assert panel.result is None and panel.payload is None
    assert panel.tree.nodes == {} and panel.pairs_table.rowCount() == 0
    assert panel.matrix_table.rowCount() == 0
    assert panel.link.value() == NO_LINK_THRESHOLD
    assert panel.window_button.isEnabled() is False
    assert panel.reread_button.isEnabled() is False
    assert "every sample, assembly and stored result in this project is untouched" in \
        panel.status.text()
    assert window.project.get_sample(first) == before
    assert len(window.project.samples()) == 2
    # The cohort itself is not a result, so it is listed again immediately.
    assert len(panel.cohort) == 2
    assert window.test_errors == []


def test_the_page_fits_the_snp_station_and_the_windows_own_clear_empties_it(
        panel, window, tmp_path, monkeypatch):
    """The wiring the workspace has to add: adopt the station, and Clear reaches here.

    `MainWindow.clear_page` already calls an adopted station's own `clear`, so this
    asserts the contract from the window's side rather than from the panel's.
    """
    from PySide6.QtWidgets import QMessageBox
    assembled(window, tmp_path, "iso-a")
    assembled(window, tmp_path, "iso-b")
    panel.refresh_cohort()
    panel.show_result(fake_run([pair("iso-a", "iso-b", 3)]))
    assert window.adopt_station("snp", panel, title="SNP tree") is True
    assert window.navigate("snp") is not False
    assert window.pages.current_key() == "snp"
    assert "planned" not in window.pages.tabToolTip(window.pages.position_of("snp"))
    monkeypatch.setattr(QMessageBox, "question",
                        lambda *args, **kwargs: QMessageBox.StandardButton.Yes)
    assert window.clear_page("snp") is True
    assert panel.payload is None and panel.tree.nodes == {}
    assert len(window.project.samples()) == 2
    assert window.test_errors == []


def test_the_tab_measures_real_assemblies_with_the_native_engine_or_says_it_has_none(
        panel, window, tmp_path, qtbot):
    for name, sequence in genomes().items():
        assembled(window, tmp_path, name, sequence)
    panel.refresh_cohort()
    assert len([row for row in panel.cohort if row["eligible"]]) == 4
    if not runtime_capabilities()["available"]:
        assert panel.run_button.isEnabled() is False
        assert "not staged" in panel.engine.text() or "not available" in panel.status.text()
        panel.run_tree()
        assert panel.payload is None
        return
    panel.run_tree()
    qtbot.waitUntil(lambda: panel.payload is not None and not window.worker_role, timeout=60000)
    assert window.test_errors == []
    distances = {(row["source"], row["target"]): row["distance"] for row in panel.payload["pairs"]}
    assert distances[("iso-a", "iso-b")] == 3 and distances[("iso-a", "iso-c")] == 12
    assert distances[("iso-b", "iso-c")] == 15
    # The fragment differs nowhere and is still comparable to nothing.
    assert [distances[key] for key in distances if "iso-d" in key] == [None, None, None]
    assert [(edge["source"], edge["target"], edge["distance"])
            for edge in panel.payload["graph"]["edges"]] == [("iso-a", "iso-b", 3),
                                                             ("iso-a", "iso-c", 12)]
    assert sorted(panel.tree.nodes) == ["iso-a", "iso-b", "iso-c", "iso-d"]
    written = json.loads((panel.output_root() / panel.payload["run_id"] / "snp-tree.json")
                         .read_text(encoding="utf-8"))
    assert written["summary"] == panel.payload["summary"]
    assert "3 of 6 pairs carry a SNP distance" in panel.status.text()
    assert "3 shared too little sequence" in panel.status.text()


def _fastq(directory, name):
    path = directory / f"{name}.fastq"
    path.write_text("@read1\nACGTACGTAC\n+\nIIIIIIIIII\n", encoding="ascii")
    return path
