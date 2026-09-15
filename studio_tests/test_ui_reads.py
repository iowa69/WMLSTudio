"""The two pipeline tabs between Samples and MLST: read trimming, and assembly.

These tests drive the real pages in the real window. The native engines are stood
in for, because a stub proves the wiring and never proves that fastp trims or that
an assembler assembles; the numbers on screen, though, are produced by the engines'
own readers (``read_tools.parse_report`` and ``assembly.assembly_metrics``) from
files on disk, so a change in what an engine reports is caught here too.
"""

import json
from pathlib import Path

import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QDialog, QLabel, QMessageBox, QPushButton

from wmlstudio import read_tools, ui_reads
from wmlstudio.assembly import assembly_metrics
from wmlstudio.ui_tabs import PAGE_PURPOSE, PIPELINE, PLANNED, clear_promise


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from wmlstudio.app import MainWindow
    # A modal warning box in a headless run blocks the event loop for good, so the
    # window records its errors here instead of showing them.
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    if widget.worker is not None and widget.worker.isRunning():
        widget.worker.cancel()
        qtbot.waitUntil(lambda: not widget.worker.isRunning(), timeout=15000)
    widget.close()


def settled(window, qtbot):
    """Wait until the background task has finished and its callbacks have run."""
    qtbot.waitUntil(lambda: window.worker is None or
                    (not window.worker.isRunning() and not window.worker_role), timeout=30000)


def read_pair(window, tmp_path, stem="isolate"):
    """Two imported read records, the way a user's paired FASTQ set arrives."""
    paths = []
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    for mate in (1, 2):
        path = tmp_path / f"{stem}_R{mate}.fastq"
        path.write_text("".join(f"@read{index}/{mate}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
                                for index in range(4)), encoding="ascii")
        paths.append(path)
    window.import_paths(paths)
    samples = window.project.samples()[-2:]
    return paths, [sample["id"] for sample in samples]


def fastp_report(tmp_path, *, before_reads=1000, after_reads=900):
    """A fastp JSON report, read back by the engine's own parser, never hand-built."""
    path = tmp_path / "fastp.json"
    path.write_text(json.dumps({
        "summary": {
            "fastp_version": read_tools.VERSION, "sequencing": "paired end (150 + 150 cycles)",
            "before_filtering": {"total_reads": before_reads, "total_bases": before_reads * 150,
                                 "q20_bases": before_reads * 140, "q30_bases": before_reads * 120,
                                 "q20_rate": 0.93, "q30_rate": 0.80, "gc_content": 0.33,
                                 "read1_mean_length": 150, "read2_mean_length": 150},
            "after_filtering": {"total_reads": after_reads, "total_bases": after_reads * 146,
                                "q20_bases": after_reads * 143, "q30_bases": after_reads * 131,
                                "q20_rate": 0.98, "q30_rate": 0.90, "gc_content": 0.33,
                                "read1_mean_length": 146, "read2_mean_length": 146}},
        "filtering_result": {"passed_filter_reads": after_reads,
                             "low_quality_reads": before_reads - after_reads,
                             "too_many_N_reads": 0, "adapter_dimer_reads": 0,
                             "too_short_reads": 0, "too_long_reads": 0},
        "duplication": {"rate": 0.021},
        "insert_size": {"peak": 308, "unknown": 27},
        "adapter_cutting": {"adapter_trimmed_reads": 121, "adapter_trimmed_bases": 5745,
                            "read1_adapter_sequence": "AGATCGGAAGAGC",
                            "read2_adapter_sequence": "AGATCGGAAGAGC"},
    }), encoding="utf-8")
    return read_tools.parse_report(path)


def record_trimming(window, identifiers, paths, tmp_path, *, report=None):
    """The record a completed fastp run leaves on both read records."""
    report = report or fastp_report(tmp_path)
    output = tmp_path / "trimmed"
    output.mkdir(exist_ok=True)
    html = output / "fastp.html"
    html.write_text("<h1>TEST ONLY fastp report fixture</h1>", encoding="utf-8")
    for index, (identifier, original) in enumerate(zip(identifiers, paths), 1):
        trimmed = output / f"trimmed_R{index}.fastq.gz"
        trimmed.write_bytes(b"TEST ONLY trimmed output fixture")
        window.project.update_metadata(identifier, {"read_trimming": {
            "status": "completed", "engine": "fastp", "version": read_tools.VERSION,
            "mate": index, "paired_with": identifiers[index % 2],
            "original_path": str(original), "original_sha256": "a" * 64,
            "trimmed_path": str(trimmed), "trimmed_sha256": "b" * 64,
            "output_directory": str(output), "report": report,
            "report_json_path": str(output / "fastp.json"), "report_html_path": str(html),
            "parameters": {"min_length": 15, "quality_threshold": 15, "detect_adapter": True},
            "completed_at": "2026-01-02T03:04:05+00:00",
            "notes": [read_tools.TRIMMING_DISCLAIMER],
        }})
    window.refresh()
    return report


def record_assembly(window, identifier, tmp_path, *, read_source, name="contigs"):
    """An assembly record whose metrics are measured from a real FASTA on disk."""
    path = tmp_path / f"{name}.fasta"
    path.write_text(">Contig_1_42.5\n" + "ACGTACGGGC" * 30 + "\n"
                    ">Contig_2_38.1_Circ\n" + "ACGTTTACGC" * 12 + "\n", encoding="ascii")
    metrics = assembly_metrics(path)
    window.project.set_input_path(identifier, path)
    window.project.update_metadata(identifier, {
        "assembly": {
            "assembly_path": str(path), "metrics": metrics,
            "coverage": {"input_read_bases": 5000, "read_bases_per_assembled_base": 12.5,
                         "denominator": "Bases in the reads given to the assembler, over the "
                                        "assembled length.",
                         "basis": "Not genome coverage: the genome size is unknown."},
            "pairing": {"explicit_mates": True}, "read_source": read_source,
            "provenance": {"engine": {"name": "SKESA", "version": "2.5.1"},
                           "completed_at": "2026-01-02T05:06:07+00:00",
                           "read_source": read_source,
                           "read_preprocessing": "TEST ONLY provenance fixture"},
        },
        "workflow": {"source_kind": "assembly"},
    })
    window.refresh()
    return metrics


def unavailable(monkeypatch, reason="Native fastp is not staged in this package, so read "
                                    "trimming is unavailable here."):
    monkeypatch.setattr(read_tools, "runtime_capabilities", lambda root=None: {
        "available": False, "binary": None, "version": read_tools.VERSION,
        "platform": "test-only", "reason": reason})
    return reason


def buttons(page):
    return [child.text() for child in page.findChildren(QPushButton)]


# ---------------------------------------------------------------------------
# The pipeline: two real tabs where two reserved slots were.
# ---------------------------------------------------------------------------


def test_the_read_and_assembly_tabs_are_real_pages_in_the_workflows_own_order(window):
    """Samples → Read QC → Assembly → MLST, with neither slot still reserved."""
    shown = [window.pages.tabText(position) for position in range(window.pages.count())]
    assert shown[:5] == ["Overview", "Samples", "Read QC", "Assembly", "MLST"]
    assert shown.index("SNP tree") == shown.index("cgMLST tree") + 1
    assert "reads" not in PLANNED and "assembly" not in PLANNED
    for key in ("reads", "assembly"):
        assert "planned" not in window.pages.tabToolTip(window.pages.position_of(key))
        assert window.stations[key]["adopted"] is not None
        assert window.stations[key]["placeholder"].isVisibleTo(window) is False
        assert "Planned for a later round" not in PAGE_PURPOSE[key]
    assert isinstance(window.read_trimming_page, ui_reads.ReadTrimmingPanel)
    assert isinstance(window.assembly_page, ui_reads.AssemblyPanel)
    # The pages are addressed exactly as before; adopting one moved no number.
    assert window.pages.tab_order() == PIPELINE
    assert window.page_index["reads"] == 7 and window.page_index["assembly"] == 8


def test_each_new_tab_states_its_limit_beside_the_work_it_offers(window):
    for key, limit in (("reads", "does not validate an isolate"),
                       ("assembly", "not a finished genome")):
        window.navigate(key)
        page = window.pages.currentWidget()
        lines = [child.text() for child in page.findChildren(QLabel)]
        assert any(limit in line for line in lines), key
    trimming = [child.text() for child in
                window.read_trimming_page.findChildren(QLabel)]
    assert any(read_tools.TRIMMING_DISCLAIMER in line for line in trimming)
    assembling = [child.text() for child in window.assembly_page.findChildren(QLabel)]
    assert any(ui_reads.ASSEMBLY_DISCLAIMER in line for line in assembling)


def test_the_next_step_button_on_each_tab_runs_that_tabs_own_work(window, monkeypatch):
    calls = []
    monkeypatch.setattr(ui_reads.ReadTrimmingPanel, "trim_selected",
                        lambda self: calls.append(("trim", self.host.pages.current_key())))
    monkeypatch.setattr(ui_reads.AssemblyPanel, "assemble_selected",
                        lambda self: calls.append(("assemble", self.host.pages.current_key())))
    window.navigate("overview")
    window.run_read_trimming()
    window.run_assembly_station()
    # Each one takes you to the tab that owns the work before doing it, so the
    # result of a run is never reported on a page the user cannot see.
    assert calls == [("trim", "reads"), ("assemble", "assembly")]


# ---------------------------------------------------------------------------
# Read QC: fastp's own report, and an honest refusal when fastp is not there.
# ---------------------------------------------------------------------------


def test_an_absent_fastp_is_reported_in_the_engines_own_words_with_no_way_around_it(
        window, monkeypatch):
    reason = unavailable(monkeypatch)
    panel = window.read_trimming_page
    panel.refresh()
    assert panel.capability.text() == reason, "the engine's sentence, not a paraphrase"
    assert panel.trim_button.isEnabled() is False
    # No affordance offers to trim anyway, or to substitute something for fastp.
    assert not [text for text in buttons(panel)
                if any(word in text.casefold() for word in ("anyway", "force", "without"))]


def test_trimming_is_refused_rather_than_substituted_when_the_tool_is_absent(
        window, tmp_path, monkeypatch):
    reason = unavailable(monkeypatch)
    read_pair(window, tmp_path)
    notices = []
    monkeypatch.setattr(type(window), "notify", lambda self, message: notices.append(message))
    panel = window.read_trimming_page
    panel.refresh()
    assert panel.trim_selected() == reason
    assert notices == [reason]
    assert window.worker is None or not window.worker.isRunning()
    assert window.test_errors == []


def test_an_untrimmed_read_file_is_listed_without_a_number_being_invented(window, tmp_path):
    read_pair(window, tmp_path)
    panel = window.read_trimming_page
    panel.refresh()
    assert panel.table.rowCount() == 2
    names = [panel.table.item(row, 0).text() for row in range(2)]
    assert names == ["isolate_R1", "isolate_R2"], "rows land in sample-name order"
    assert panel.table.item(0, 2).text() == "Not trimmed"
    # Every measured column is an em dash, never a zero: nothing was measured.
    assert [panel.table.item(0, column).text() for column in range(3, 9)] == ["—"] * 6
    assert "0 of 2 read files carry a trimming record" in panel.status.text()
    panel.table.selectRow(0)
    detail = panel.detail.toPlainText()
    assert "were not trimmed" in detail
    assert read_tools.TRIMMING_DISCLAIMER in detail


def test_a_recorded_run_is_shown_with_the_denominator_it_was_measured_over(window, tmp_path):
    paths, identifiers = read_pair(window, tmp_path)
    report = record_trimming(window, identifiers, paths, tmp_path)
    panel = window.read_trimming_page
    panel.refresh()
    row = [panel.table.item(0, column).text() for column in range(9)]
    assert row[2] == "Trimmed"
    assert row[3] == "1,000" and row[4] == "900"
    assert row[5] == "10.0%" and row[7] == "80.0% → 90.0%"
    # A percentage is never shown without saying what it is a percentage of.
    assert panel.table.item(0, 5).toolTip() == report["denominators"]["reads_removed_percent"]
    assert "not of reads" in panel.table.item(0, 7).toolTip()
    panel.table.selectRow(0)
    detail = panel.detail.toPlainText()
    assert "1,000 before → 900 after (100 removed, 10.0% of the reads present before filtering)" \
        in detail
    assert "fastp's estimate from read content; not a library-preparation measurement" in detail
    assert "Adapter-trimmed   121 reads, 5,745 bases · detected AGATCGGAAGAGC" in detail
    assert str(paths[0]) in detail, "the original is named, and named as unchanged"
    assert "The original FASTQ is unchanged and remains this record's input." in detail
    assert "2 of 2 samples carry a trimming record" in window.station_status["reads"].text()


def test_the_report_button_opens_the_file_fastp_wrote_or_says_there_is_none(
        window, tmp_path, monkeypatch):
    paths, identifiers = read_pair(window, tmp_path)
    record_trimming(window, identifiers, paths, tmp_path)
    panel = window.read_trimming_page
    panel.refresh()
    opened = []
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    panel.table.selectRow(0)
    panel.open_report()
    assert opened == [str(tmp_path / "trimmed/fastp.html")]
    # A report that is no longer on this computer is said to be missing, not opened.
    (tmp_path / "trimmed/fastp.html").unlink()
    message = panel.open_report()
    assert "has no fastp report on this computer" in message
    assert len(opened) == 1


def test_an_unconfirmed_pairing_runs_nothing_at_all(window, tmp_path, monkeypatch):
    read_pair(window, tmp_path)

    class Refused(QDialog):
        def __init__(self, samples, parent=None):
            super().__init__(parent)
            self.assignments = []

        def exec(self):
            return QDialog.DialogCode.Rejected

    monkeypatch.setattr(ui_reads, "PairReadsDialog", Refused)
    panel = window.read_trimming_page
    panel.refresh()
    assert panel.trim_selected() is False
    assert "nothing ran and nothing changed" in panel.status.text()
    assert window.worker is None or not window.worker.isRunning()


def test_one_read_file_on_its_own_is_not_treated_as_a_pair(window, tmp_path):
    paths, _identifiers = read_pair(window, tmp_path)
    panel = window.read_trimming_page
    panel.refresh()
    panel.table.selectRow(0)
    # Only the forward file is selected, so there is no mate to confirm and the
    # page says which second file it needs rather than guessing at one.
    assert panel.candidates() == [] or len(panel.candidates()) == 1
    panel.trim_selected()
    assert "at least two read files" in panel.status.text()
    assert "never evidence" in panel.status.text()
    assert window.worker is None or not window.worker.isRunning()


def test_trimming_records_every_confirmed_pair_and_then_shows_it(window, tmp_path,
                                                                 qtbot, monkeypatch):
    """The wiring, with a stub engine: a stub never proves that fastp trims."""
    paths, identifiers = read_pair(window, tmp_path)
    report = fastp_report(tmp_path)
    recorded = []

    class Confirmed(QDialog):
        def __init__(self, samples, parent=None):
            super().__init__(parent)
            self.assignments = [{"primary_id": samples[0]["id"], "mate_id": samples[1]["id"]}]

        def exec(self):
            return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ui_reads, "PairReadsDialog", Confirmed)
    monkeypatch.setattr(read_tools, "run_fastp",
                        lambda read1, read2, destination, **options: {
                            "read1_path": str(read1), "read2_path": str(read2),
                            "report": report, "provenance": {"parameters": options},
                            "notes": [read_tools.TRIMMING_DISCLAIMER]})
    monkeypatch.setattr(read_tools, "record_trimming",
                        lambda project, first, second, result: recorded.append(
                            (first, second, result["provenance"]["parameters"])))
    panel = window.read_trimming_page
    panel.refresh()
    panel.min_length.setValue(30)
    panel.deduplicate.setChecked(True)
    assert "Trimming 1 read pair(s)" in panel.trim_selected()
    settled(window, qtbot)
    assert [item[:2] for item in recorded] == [(identifiers[0], identifiers[1])]
    options = recorded[0][2]
    # The settings on screen are the settings the run was given, and are recorded.
    assert options["min_length"] == 30 and options["deduplicate"] is True
    assert options["quality_threshold"] == 15 and options["detect_adapter"] is True
    assert "fastp removed 100 reads in total" in panel.status.text()
    assert read_tools.TRIMMING_DISCLAIMER in panel.status.text()
    # A run is followed by several refreshes; the result must survive all of them,
    # because it is the one line the user stayed on this tab to read.
    window.refresh()
    panel.refresh()
    assert "fastp removed 100 reads in total" in panel.status.text()
    assert window.test_errors == []


# ---------------------------------------------------------------------------
# Assembly: the contigs, and which reads actually went into them.
# ---------------------------------------------------------------------------


def test_an_assembly_is_listed_with_the_metrics_its_own_run_measured(window, tmp_path):
    paths, identifiers = read_pair(window, tmp_path)
    metrics = record_assembly(window, identifiers[0], tmp_path, read_source={
        "kind": "fastp_trimmed", "tool": "fastp", "version": read_tools.VERSION,
        "description": "fastp 1.3.7 trimmed and quality-filtered reads; the original FASTQs "
                       "are unchanged and still on record.",
        "reads_removed": 100, "bases_removed": 14_600,
        "originals": [{"mate": 1, "path": str(paths[0])}, {"mate": 2, "path": str(paths[1])}]})
    panel = window.assembly_page
    panel.refresh()
    row = [panel.table.item(0, column).text() for column in range(9)]
    assert row[1] == "Assembled here"
    assert row[2] == f"{metrics['contigs']:,}" == "2"
    assert row[3] == f"{metrics['total_length']:,} bp" == "420 bp"
    assert row[4] == f"{metrics['n50']:,} bp" == "300 bp"
    assert row[5] == f"{metrics['largest_contig']:,} bp" == "300 bp"
    assert row[6] == "fastp-trimmed reads"
    assert row[7] == "SKESA 2.5.1"
    panel.table.selectRow(0)
    detail = panel.detail.toPlainText()
    assert "Contigs        2" in detail and "N50            300 bp" in detail
    assert "100 reads and 14,600 bases were removed before assembly." in detail
    # The assembler's own depth estimate is labelled as exactly that.
    assert "SKESA's own k-mer depth estimate, read from each contig name." in detail
    assert "1 contig(s) reported circular" in detail
    assert ui_reads.ASSEMBLY_DISCLAIMER in [child.text() for child in
                                            panel.findChildren(QLabel)]
    assert "Not genome coverage" in detail
    assert "1 of 2 samples carry an assembly made here" in window.station_status["assembly"].text()


@pytest.mark.parametrize("source,shown", [
    ({"kind": "original_reads", "tool": None, "version": None,
      "description": "The reads exactly as supplied; no trimming or quality filtering was "
                     "performed."}, "Original reads, untrimmed"),
    ({"kind": "unrecorded", "description": "This assembly predates read-source recording."},
     "Not recorded"),
])
def test_an_assembly_never_claims_a_preprocessing_step_it_cannot_show(window, tmp_path,
                                                                      source, shown):
    _paths, identifiers = read_pair(window, tmp_path)
    record_assembly(window, identifiers[0], tmp_path, read_source=source)
    panel = window.assembly_page
    panel.refresh()
    assert panel.table.item(0, 6).text() == shown
    assert source["description"] in panel.table.item(0, 6).toolTip()


def test_reads_that_were_never_assembled_say_so_rather_than_showing_a_zero(window, tmp_path):
    read_pair(window, tmp_path)
    panel = window.assembly_page
    panel.refresh()
    assert [panel.table.item(row, 1).text() for row in range(2)] == ["Reads, not assembled"] * 2
    assert [panel.table.item(0, column).text() for column in range(2, 9)] == ["—"] * 7
    assert "0 of 2 samples carry an assembly made here" in panel.status.text()


def test_a_number_column_sorts_by_its_quantity_and_keeps_the_unmeasured_apart(window, tmp_path):
    _paths, identifiers = read_pair(window, tmp_path)
    _third, more = read_pair(window, tmp_path / "second", stem="other")
    record_assembly(window, identifiers[0], tmp_path, name="small",
                    read_source={"kind": "original_reads", "description": "as supplied"})
    big = tmp_path / "big.fasta"
    big.write_text("".join(f">Contig_{index}_20.0\n{'ACGT' * 25}\n" for index in range(10)),
                   encoding="ascii")
    window.project.set_input_path(more[0], big)
    window.project.update_metadata(more[0], {
        "assembly": {"assembly_path": str(big), "metrics": assembly_metrics(big),
                     "provenance": {"engine": {"name": "SKESA", "version": "2.5.1"},
                                    "completed_at": "2026-01-03T00:00:00+00:00"}},
        "workflow": {"source_kind": "assembly"}})
    window.refresh()
    panel = window.assembly_page
    panel.refresh()
    panel.table.sortItems(2, Qt.SortOrder.AscendingOrder)
    contigs = [panel.table.item(row, 2).text() for row in range(panel.table.rowCount())]
    measured = [text for text in contigs if text != "—"]
    # "10" must not sort before "2" as text would have it, and a row with nothing
    # measured is not ranked as though it held a zero.
    assert measured == ["2", "10"]
    assert contigs[:2] == ["—", "—"]


def test_the_assembly_tab_assembles_without_typing_anything(window, tmp_path, qtbot,
                                                            monkeypatch):
    """The wiring, with a stub engine: a stub never proves that anything assembled."""
    _paths, identifiers = read_pair(window, tmp_path)
    from wmlstudio import assembly as assembly_module
    typed, attached = [], []

    class Confirmed(QDialog):
        def __init__(self, samples, parent=None):
            super().__init__(parent)
            self.assignments = [{"primary_id": samples[0]["id"], "mate_id": samples[1]["id"]}]

        def exec(self):
            return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ui_reads, "PairReadsDialog", Confirmed)
    monkeypatch.setattr(assembly_module, "run_skesa",
                        lambda read1, read2, destination, **options: {
                            "assembly_path": str(destination / "contigs.fasta"),
                            "read_source": options.get("read_source")
                            or assembly_module.RAW_READ_SOURCE})
    monkeypatch.setattr(assembly_module, "associate_assembly",
                        lambda project, primary, mate, result: attached.append(
                            (primary, mate, result["read_source"]["kind"])))
    monkeypatch.setattr(type(window), "begin_typing",
                        lambda self, samples, scheme, **kwargs: typed.append(samples))
    panel = window.assembly_page
    panel.refresh()
    assert "Assembling 1 read pair(s)" in panel.assemble_selected()
    settled(window, qtbot)
    assert attached == [(identifiers[0], identifiers[1], "original_reads")]
    # Nothing was typed: this tab assembles, and typing is a separate decision.
    assert typed == []
    assert "0 used trimmed reads and 1 used the reads as supplied" in panel.status.text()
    assert "type these isolates from the MLST or cgMLST tab" in panel.status.text()
    assert window.test_errors == []


def test_a_trimmed_pair_is_what_the_assembler_is_handed(window, tmp_path, qtbot, monkeypatch):
    paths, identifiers = read_pair(window, tmp_path)
    record_trimming(window, identifiers, paths, tmp_path)
    from wmlstudio import assembly as assembly_module
    given = []

    class Confirmed(QDialog):
        def __init__(self, samples, parent=None):
            super().__init__(parent)
            self.assignments = [{"primary_id": samples[0]["id"], "mate_id": samples[1]["id"]}]

        def exec(self):
            return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ui_reads, "PairReadsDialog", Confirmed)
    monkeypatch.setattr(assembly_module, "run_skesa",
                        lambda read1, read2, destination, **options: given.append(
                            (Path(read1).name, options["read_source"]["kind"])) or {
                            "assembly_path": str(destination / "contigs.fasta"),
                            "read_source": options["read_source"]})
    monkeypatch.setattr(assembly_module, "associate_assembly",
                        lambda project, primary, mate, result: None)
    panel = window.assembly_page
    panel.refresh()
    panel.assemble_selected()
    settled(window, qtbot)
    assert given == [("trimmed_R1.fastq.gz", "fastp_trimmed")]
    assert "1 used trimmed reads and 0 used the reads as supplied" in panel.status.text()


def test_an_assembly_that_cannot_be_attached_is_kept_and_named_rather_than_lost(
        window, tmp_path, qtbot, monkeypatch):
    """A finished assembly is minutes of work; a refused record must not discard it."""
    _paths, _identifiers = read_pair(window, tmp_path)
    from wmlstudio import assembly as assembly_module

    class Confirmed(QDialog):
        def __init__(self, samples, parent=None):
            super().__init__(parent)
            self.assignments = [{"primary_id": samples[0]["id"], "mate_id": samples[1]["id"]}]

        def exec(self):
            return QDialog.DialogCode.Accepted

    monkeypatch.setattr(ui_reads, "PairReadsDialog", Confirmed)
    monkeypatch.setattr(assembly_module, "run_skesa",
                        lambda read1, read2, destination, **options: {
                            "assembly_path": str(destination / "contigs.fasta"),
                            "read_source": assembly_module.RAW_READ_SOURCE})

    def refuse(project, primary, mate, result):
        raise ValueError("Assembly read provenance does not match the selected sample records.")

    monkeypatch.setattr(assembly_module, "associate_assembly", refuse)
    panel = window.assembly_page
    panel.refresh()
    panel.assemble_selected()
    settled(window, qtbot)
    message = panel.status.text()
    assert "0 read pair(s) were assembled and attached to their samples" in message
    assert "1 assembly(ies) finished but could not be attached" in message
    assert "does not match the selected sample records" in message
    assert "The contigs are kept at" in message
    # The refusal is reported, not raised: the rest of a queue still runs.
    assert window.test_errors == []


def test_metrics_can_be_read_from_the_file_without_being_written_back(window, tmp_path, qtbot):
    _paths, identifiers = read_pair(window, tmp_path)
    record_assembly(window, identifiers[0], tmp_path,
                    read_source={"kind": "original_reads", "description": "as supplied"})
    # An assembly recorded before per-contig metrics were kept.
    stored = window.project.get_sample(identifiers[0])["metadata"]["assembly"]
    window.project.update_metadata(identifiers[0], {"assembly": {**stored, "metrics": None}})
    window.refresh()
    panel = window.assembly_page
    panel.refresh()
    assert panel.table.item(0, 2).text() == "—"
    panel.table.selectRow(0)
    panel.measure_selected()
    settled(window, qtbot)
    detail = panel.detail.toPlainText()
    assert "Metrics read from the FASTA just now; nothing was written back" in detail
    assert "Contigs        2" in detail
    # The project still holds no metrics: reading is not recording.
    assert window.project.get_sample(identifiers[0])["metadata"]["assembly"]["metrics"] is None
    assert window.test_errors == []


# ---------------------------------------------------------------------------
# Clear, on both tabs: a view is emptied, evidence never is.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key,attribute", [("reads", "read_trimming_page"),
                                           ("assembly", "assembly_page")])
def test_clearing_one_of_these_tabs_empties_the_view_and_keeps_every_record(
        window, tmp_path, monkeypatch, key, attribute):
    paths, identifiers = read_pair(window, tmp_path)
    record_trimming(window, identifiers, paths, tmp_path)
    panel = getattr(window, attribute)
    panel.refresh()
    assert panel.table.rowCount() == 2
    asked = []
    monkeypatch.setattr(QMessageBox, "question",
                        lambda parent, title, text, *args: (asked.append(text),
                                                            QMessageBox.StandardButton.Yes)[1])
    assert window.clear_page(key) is True
    assert clear_promise(key)["keeps"] in asked[-1]
    assert panel.table.rowCount() == 0
    assert panel.detail.toPlainText() == ""
    # Nothing stored moved: the trimming record and both original files are intact.
    stored = window.project.get_sample(identifiers[0])["metadata"]["read_trimming"]
    assert stored["status"] == "completed"
    assert [path.is_file() for path in paths] == [True, True]
    panel.refresh()
    assert panel.table.rowCount() == 2
