"""The two pipeline pages between Samples and MLST: read trimming, and assembly.

Both pages render what the engines recorded and never recompute, round or relabel
a number of their own. :mod:`wmlstudio.read_tools` owns everything fastp says, and
:mod:`wmlstudio.assembly` owns everything the assembler says; a cell with no
measurement behind it prints an em dash rather than a zero.

The claim boundaries these pages have to hold, because they are the ones a
non-bioinformatician could read past:

* Trimming improves reads. It does not validate an isolate, prove purity,
  establish a species, or make a typing result correct.
* An assembly is a reconstruction. Contiguity is not completeness, not purity and
  not a species assignment, and no threshold on this page decides that an assembly
  is good enough -- the numbers are shown and the judgement stays the user's.
* The original FASTQ files are read-only inputs throughout. Trimming writes new
  files beside them; nothing is renamed, replaced or deleted.

Neither page offers a way to trim without the staged native tool: where fastp is
not in this package the tab says so in the engine's own words and the reads stay
usable exactly as supplied.
"""

from __future__ import annotations

import uuid
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QSpinBox,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import read_tools
from wmlstudio.assembly import assembly_metrics, assembly_view
from wmlstudio.pairing_dialog import PairReadsDialog
from wmlstudio.ui_common import FlowLayout, make_table
from wmlstudio.widgets import button, label

#: Reading one header tells reads from an assembly for a record that has never been
#: analysed. That is one small read per row, so a cohort larger than this is listed
#: from what is stored and the page says which rows it could not classify.
INSPECT_LIMIT = 250

#: What each engine status is called on screen. The words are the page's, the states
#: are the engines'; none of them is ever folded into another.
TRIM_STATUS = {
    "trimmed": "Trimmed",
    "not_trimmed": "Not trimmed",
    "not_reads": "Not a read file",
}

ASSEMBLY_STATUS = {
    "assembled": "Assembled here",
    "read_mate": "Second mate of a pair",
    "reads_not_assembled": "Reads, not assembled",
    "assembly_supplied": "Assembly supplied",
    "unknown": "Input not read",
}

#: Which reads an assembly actually consumed. "Not recorded" is its own answer and
#: never becomes "original": an assembly made before read sources were recorded
#: cannot say which files it was given.
READ_SOURCE_WORDS = {
    "fastp_trimmed": "fastp-trimmed reads",
    "original_reads": "Original reads, untrimmed",
    "unrecorded": "Not recorded",
}

ASSEMBLY_DISCLAIMER = (
    "An assembly is a reconstruction of a genome, not a finished genome. Contiguity is "
    "not completeness, not a purity check and not a species assignment, and no threshold "
    "on this page decides whether an assembly is good enough to type."
)


class SortableCell(QTableWidgetItem):
    """A cell that sorts by the quantity behind it, not by its printed text.

    "1,204" and "987" sort the wrong way round as strings, and a column of contig
    counts that puts 10 before 9 is a column nobody trusts twice. A row with
    nothing measured carries no key: those rows group together at one end instead
    of being ranked as though they held a zero.
    """

    def __init__(self, text, key=None, sample_id=None, tooltip=None):
        super().__init__(str(text))
        self.sort_key = key
        self.setData(Qt.ItemDataRole.UserRole, sample_id)
        self.setToolTip(str(tooltip) if tooltip else str(text))

    def __lt__(self, other):
        mine, theirs = self.sort_key, getattr(other, "sort_key", None)
        if mine is None or theirs is None:
            return theirs is not None
        return mine < theirs


def _measured(value):
    """True only for a real number. False is not a count and None is not a zero."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def number_cell(value, sample_id, *, unit="", tooltip=""):
    """A counted quantity, grouped for reading and kept exact for sorting."""
    if not _measured(value):
        return SortableCell("—", None, sample_id, tooltip or "Not measured for this record.")
    text = f"{value:,}" if isinstance(value, int) else f"{value:,.1f}"
    return SortableCell(f"{text} {unit}".strip(), value, sample_id, tooltip)


def percent_cell(value, sample_id, *, tooltip=""):
    """A percentage that always travels with the denominator it was measured over."""
    if not _measured(value):
        return SortableCell("—", None, sample_id, tooltip or "Not measured for this record.")
    return SortableCell(f"{value:.1f}%", value, sample_id, tooltip)


def text_cell(text, sample_id, *, tooltip=""):
    value = "—" if text in (None, "") else str(text)
    return SortableCell(value, value.casefold(), sample_id, tooltip)


def _rate(value):
    """fastp reports Q20/Q30 as a fraction of bases; this is the same number, shown."""
    return 100 * value if _measured(value) else None


def _short(moment):
    """An ISO timestamp as a date and a time, without inventing a timezone for it."""
    return str(moment).replace("T", " ")[:19] if moment else "—"


def is_read_record(sample):
    """True when this record's own input is still a read file.

    A record whose reads were assembled here now holds the assembly as its input,
    so it is no longer something trimming or assembly can be run on; its linked
    mate still is, and is excluded separately by whoever pairs the rows.
    """
    metadata = sample.get("metadata") or {}
    if (metadata.get("workflow") or {}).get("source_kind") == "assembly":
        return False
    if (sample.get("result") or {}).get("kind") == "fastq":
        return True
    name = str(sample.get("input_path") or "").lower()
    return name.removesuffix(".gz").removesuffix(".bz2").endswith((".fq", ".fastq"))


class PipelinePanel(QWidget):
    """Shared scaffolding for the read and assembly pages.

    The host window supplies ``project``, ``launch_task``, ``notify`` and
    ``project_path`` -- the same contract every other workspace page uses, so these
    panels can be adopted onto their stations without either side knowing more
    about the other.
    """

    REQUIRED = ("project", "launch_task", "notify", "project_path")
    COLUMNS: tuple[str, ...] = ()
    READY_MESSAGE = ""

    def __init__(self, window, parent=None):
        # Checked before the widget exists, so a wiring mistake is named here rather
        # than surfacing later as a missing attribute in the middle of a run.
        missing = [name for name in self.REQUIRED if not hasattr(window, name)]
        if missing:
            raise TypeError(f"{type(self).__name__} needs a workspace window providing "
                            + ", ".join(missing) + ".")
        super().__init__(parent if parent is not None else window)
        self.host = window
        self.views = {}     # sample id -> the engine view this page last drew
        self.outcome = ""   # what the last run on this page actually did
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self.build(layout)
        self.status = label(self.READY_MESSAGE, "small", True)
        layout.addWidget(self.status)

    # --- construction --------------------------------------------------------
    def build(self, layout):
        raise NotImplementedError

    def build_table(self, layout, hint):
        self.table = make_table(list(self.COLUMNS))
        self.table.setMinimumHeight(180)
        # Sample name, ascending, until the reader clicks a column of their own.
        # Qt's default indicator is neither of those, so a table left to it lands
        # in an order nobody chose.
        self.table.horizontalHeader().setSortIndicator(0, Qt.SortOrder.AscendingOrder)
        self.table.itemSelectionChanged.connect(self.show_detail)
        layout.addWidget(self.table, 1)
        self.detail = QTextBrowser()
        self.detail.setMinimumHeight(150)
        self.detail.setPlainText(hint)
        layout.addWidget(self.detail, 1)

    # --- the project's own records -------------------------------------------
    def say(self, message):
        self.status.setText(str(message))
        return message

    def said_with_outcome(self, summary):
        """The current summary, with the last run's own sentence still beside it.

        A run is followed immediately by a refresh -- its own, and the window's --
        and a summary that replaced the result would take away the one line the
        user was waiting to read. The sentence stays until the next run or Clear.
        """
        return self.say(f"{summary} {self.outcome}".strip() if self.outcome else summary)

    def samples(self):
        return list(self.host.project.samples())

    def selected_ids(self):
        """The sample ids of the selected rows, in the order the table shows them."""
        chosen, seen = [], set()
        for index in self.table.selectedIndexes():
            item = self.table.item(index.row(), 0)
            identifier = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
            if identifier and identifier not in seen:
                seen.add(identifier)
                chosen.append(identifier)
        return chosen

    def selected_view(self):
        """The engine view behind the first selected row, or None."""
        chosen = self.selected_ids()
        return self.views.get(chosen[0]) if chosen else None

    def fill(self, rows):
        """Draw one row per record. Sorting is suspended so keys are not reshuffled."""
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for index, cells in enumerate(rows):
            for column, item in enumerate(cells):
                self.table.setItem(index, column, item)
        self.table.setSortingEnabled(True)
        return len(rows)

    def output_root(self, name):
        return Path(self.host.project_path).with_suffix(".files") / name

    def pairs_from(self, records, *, verb):
        """Confirm the biological pairing explicitly; a filename match is only a hint."""
        if len(records) < 2:
            self.say(f"Select at least two read files to {verb}: one forward file and its "
                     "reverse mate. A filename that looks like a pair is a hint, never "
                     "evidence that two files belong together.")
            return []
        dialog = PairReadsDialog(records, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            self.say("No pairing was confirmed, so nothing ran and nothing changed.")
            return []
        return list(dialog.assignments)

    def candidates(self):
        """Read records this page can act on, in the order the project holds them.

        The selected rows when any of them is still a read file, every read record
        otherwise: a selection of rows this page cannot act on is not a refusal to
        work, and the pairing dialog is where the choice is confirmed anyway.
        """
        chosen = set(self.selected_ids())
        records = [sample for sample in self.samples() if is_read_record(sample)]
        picked = [sample for sample in records if sample["id"] in chosen]
        return picked or records

    def clear(self):
        """Empty what is on screen. Every stored file, report and record is kept."""
        self.table.clearSelection()
        self.table.setRowCount(0)
        self.views, self.outcome = {}, ""
        self.detail.setPlainText("")
        return self.say(self.READY_MESSAGE)

    def show_detail(self):
        raise NotImplementedError


class ReadTrimmingPanel(PipelinePanel):
    """fastp adapter trimming and quality filtering, with the tool's own report.

    Nothing on this page is a WMLSTudio calculation: every count, rate and curve
    comes from the fastp JSON report that run wrote, and is shown with the
    denominator fastp measured it over.
    """

    COLUMNS = ("Sample", "Read file", "Trimming", "Reads before", "Reads after",
               "Reads removed", "Bases removed", "Q30 before → after", "Trimmed at")
    READY_MESSAGE = ("Select the read files of one or more pairs, then trim them. The originals "
                     "stay exactly as they are; trimmed reads are new files beside them.")

    def build(self, layout):
        self.capability = label("", "small", True)
        layout.addWidget(self.capability)
        strip = FlowLayout()
        self.trim_button = button("Trim read pairs…", self.trim_selected, True)
        self.trim_button.setToolTip("Run fastp on the confirmed pairs. The original FASTQ files "
                                    "are opened read-only and are never modified.")
        strip.addWidget(self.trim_button)
        self.report_button = button("Open the fastp report", self.open_report)
        self.report_button.setToolTip("Open the HTML report fastp itself wrote for the selected "
                                      "sample, in your browser.")
        strip.addWidget(self.report_button)
        strip.addWidget(button("Refresh", self.refresh))
        layout.addLayout(strip)
        settings = FlowLayout()
        settings.addWidget(label("Minimum read length", "small"))
        self.min_length = self._spin(1, 10_000, 15, "Reads shorter than this after trimming are "
                                                    "discarded. This is fastp's own default.")
        settings.addWidget(self.min_length)
        settings.addWidget(label("Qualified base (Phred)", "small"))
        self.quality = self._spin(0, 60, 15, "A base at or above this Phred score counts as "
                                             "qualified. This is fastp's own default.")
        settings.addWidget(self.quality)
        self.detect_adapter = QCheckBox("Detect adapters for paired reads")
        self.detect_adapter.setChecked(True)
        self.detect_adapter.setToolTip("Let fastp infer the adapter sequence from the overlap of "
                                       "each pair, rather than assuming one.")
        settings.addWidget(self.detect_adapter)
        self.deduplicate = QCheckBox("Drop duplicate reads")
        self.deduplicate.setToolTip("fastp's own duplicate removal. Duplicates are estimated from "
                                    "read content; this is not a library-preparation measurement.")
        settings.addWidget(self.deduplicate)
        layout.addLayout(settings)
        layout.addWidget(label("Every value used is recorded with the run. " +
                               read_tools.TRIMMING_DISCLAIMER, "small", True))
        self.build_table(layout, "Select a row to read the fastp report behind it.")

    def _spin(self, low, high, value, tip):
        widget = QSpinBox()
        widget.setRange(low, high)
        widget.setValue(value)
        widget.setToolTip(tip)
        return widget

    # --- what this package can actually do -----------------------------------
    def capabilities(self):
        """Ask the engine, every time the tab is shown. Nothing is cached or assumed."""
        capability = read_tools.runtime_capabilities()
        if capability["available"]:
            self.capability.setText(
                f"fastp {capability['version']} is staged in this package "
                f"({capability['platform']}). Trimming runs on this computer; nothing is uploaded "
                "and nothing is downloaded.")
        else:
            # The engine's own sentence, word for word: it already says fastp is
            # missing, that nothing is downloaded, and that the reads stay usable.
            self.capability.setText(capability["reason"])
        self.trim_button.setEnabled(capability["available"])
        return capability

    # --- the table ------------------------------------------------------------
    def refresh(self):
        """One row per read record, from what the project stored. No file is read."""
        capability = self.capabilities()
        self.views, rows = {}, []
        for sample in self.samples():
            view = read_tools.trimming_view(sample)
            if view["status"] == "not_reads":
                continue
            self.views[sample["id"]] = view
            rows.append(self._row(sample, view))
        self.fill(rows)
        trimmed = sum(view["status"] == "trimmed" for view in self.views.values())
        if not rows:
            return self.said_with_outcome(
                "No read files are in this project yet. Import a paired FASTQ set in Samples; "
                "an assembly has no reads to trim.")
        summary = (f"{trimmed} of {len(rows)} read files carry a trimming record · "
                   f"{len(rows) - trimmed} do not.")
        if not capability["available"]:
            return self.said_with_outcome(
                summary + " Trimming cannot run in this package, and untrimmed reads remain "
                          "usable exactly as they were supplied.")
        return self.said_with_outcome(summary + " " + self.READY_MESSAGE)

    def _row(self, sample, view):
        identifier = sample["id"]
        report = view["report"] or {}
        before, after = report.get("before") or {}, report.get("after") or {}
        denominators = report.get("denominators") or {}
        return [
            text_cell(sample.get("name"), identifier,
                      tooltip="\n".join(view["notes"])),
            text_cell(Path(view["original_path"] or "").name, identifier,
                      tooltip=str(view["original_path"] or "This record has no input path.")),
            text_cell(TRIM_STATUS.get(view["status"], view["status"]), identifier,
                      tooltip=("fastp " + str(view["version"]) if view["status"] == "trimmed"
                               else "\n".join(view["notes"]))),
            number_cell(before.get("total_reads"), identifier,
                        tooltip="Reads present before filtering, as fastp counted them."),
            number_cell(after.get("total_reads"), identifier,
                        tooltip="Reads kept by the filters, as fastp counted them."),
            percent_cell(report.get("reads_removed_percent"), identifier,
                         tooltip=denominators.get("reads_removed_percent", "")),
            percent_cell(report.get("bases_removed_percent"), identifier,
                         tooltip="Of the bases present before filtering."),
            text_cell(self._q30(before, after), identifier,
                      tooltip=denominators.get("q20_rate/q30_rate", "")),
            text_cell(_short(view["completed_at"]), identifier,
                      tooltip="When this trimming run finished, in UTC."),
        ]

    def _q30(self, before, after):
        first, second = _rate(before.get("q30_rate")), _rate(after.get("q30_rate"))
        if first is None or second is None:
            return ""
        return f"{first:.1f}% → {second:.1f}%"

    # --- the report behind one row --------------------------------------------
    def show_detail(self):
        view = self.selected_view()
        if view is None:
            self.detail.setPlainText("Select a row to read the fastp report behind it.")
            return None
        lines = [f"{view['sample_name']} · {TRIM_STATUS.get(view['status'], view['status'])}"]
        report = view["report"] or {}
        if view["status"] == "trimmed" and report:
            lines += self._report_lines(view, report)
        lines += ["", "What this does and does not mean:"]
        lines += [f"  · {note}" for note in view["notes"]]
        self.detail.setPlainText("\n".join(lines))
        return view

    def _report_lines(self, view, report):
        before, after = report["before"], report["after"]
        denominators = report.get("denominators") or {}
        adapter, filtering = report.get("adapter") or {}, report.get("filtering") or {}
        insert = report.get("insert_size") or {}
        lines = [
            f"fastp {view['version']} · finished {_short(view['completed_at'])} · "
            f"{report.get('sequencing') or 'sequencing layout not reported'}",
            "",
            f"Reads   {before['total_reads']:,} before → {after['total_reads']:,} after "
            f"({report['reads_removed']:,} removed"
            + (f", {report['reads_removed_percent']:.1f}% "
               f"{denominators.get('reads_removed_percent', '')}"
               if _measured(report.get("reads_removed_percent")) else "") + ")",
            f"Bases   {before['total_bases']:,} before → {after['total_bases']:,} after "
            f"({report['bases_removed']:,} removed)",
        ]
        for name, key in (("Q20", "q20_rate"), ("Q30", "q30_rate")):
            first, second = _rate(before.get(key)), _rate(after.get(key))
            if first is not None and second is not None:
                lines.append(f"{name}     {first:.1f}% → {second:.1f}% "
                             f"({denominators.get('q20_rate/q30_rate', '')})")
        lengths = [f"R{mate} {before.get(f'read{mate}_mean_length')} → "
                   f"{after.get(f'read{mate}_mean_length')}" for mate in (1, 2)
                   if before.get(f"read{mate}_mean_length") is not None]
        if lengths:
            lines.append("Mean read length  " + " · ".join(lengths))
        if _measured(adapter.get("trimmed_reads")):
            sequences = " · ".join(filter(None, [adapter.get("read1_sequence"),
                                                 adapter.get("read2_sequence")]))
            trimmed_bases = adapter.get("trimmed_bases")
            lines.append(f"Adapter-trimmed   {adapter['trimmed_reads']:,} reads"
                         + (f", {trimmed_bases:,} bases" if _measured(trimmed_bases) else "")
                         + (f" · detected {sequences}" if sequences else ""))
        counted = " · ".join(f"{key.replace('_', ' ')}: {value:,}"
                             for key, value in filtering.items() if _measured(value))
        if counted:
            lines.append("Filters   " + counted)
        if _measured(report.get("duplication_rate")):
            lines.append(f"Duplication rate  {100 * report['duplication_rate']:.2f}% "
                         f"({denominators.get('duplication_rate', '')})")
        undetermined = insert.get("undetermined_pairs")
        if _measured(insert.get("peak")):
            lines.append(f"Insert size peak  {insert['peak']:,} bp"
                         + (f" · undetermined pairs {undetermined:,}"
                            if _measured(undetermined) else "")
                         + f" ({denominators.get('insert_size', '')})")
        parameters = view.get("parameters") or {}
        if parameters:
            lines += ["", "Settings this run used: " + " · ".join(
                f"{key}={value}" for key, value in sorted(parameters.items())
                if not isinstance(value, str) or len(value) < 40)]
        lines += ["", f"Original (unchanged)  {view['original_path']}",
                  f"Trimmed               {view['trimmed_path']}"]
        return lines

    def open_report(self):
        """Open fastp's own HTML report for the selected row, if that run wrote one."""
        view = self.selected_view()
        if view is None:
            return self.say("Select a trimmed sample first; the report belongs to one run.")
        path = view.get("report_html_path")
        if view["status"] != "trimmed" or not path or not Path(path).is_file():
            return self.say(f"{view['sample_name']} has no fastp report on this computer. A report "
                            "exists only where that run's output directory is still in place.")
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(path))))
        return self.say(f"Opened fastp's own report for {view['sample_name']}: {path}")

    # --- running fastp ---------------------------------------------------------
    def trim_selected(self):
        """Trim the confirmed pairs. Refused outright where fastp is not staged."""
        capability = self.capabilities()
        if not capability["available"]:
            self.host.notify(capability["reason"])
            return self.say(capability["reason"])
        records = self.candidates()
        pairs = self.pairs_from(records, verb="trim")
        if not pairs:
            return False
        already = [pair for pair in pairs
                   if (self.views.get(pair["primary_id"]) or {}).get("status") == "trimmed"]
        project = self.host.project
        tasks = [{"primary": project.get_sample(pair["primary_id"]),
                  "mate": project.get_sample(pair["mate_id"])} for pair in pairs]
        options = {"min_length": self.min_length.value(),
                   "quality_threshold": self.quality.value(),
                   "detect_adapter": self.detect_adapter.isChecked(),
                   "deduplicate": self.deduplicate.isChecked()}
        output_root = self.output_root("trimmed_reads")
        self.outcome = ""   # the previous run's sentence is not this run's result
        if not self._launch_trimming(tasks, options, output_root):
            return self.say("A background task is already running. Let it finish or cancel it, "
                            "then trim again.")
        note = (f" {len(already)} of them already carried a trimming record: a new run replaces "
                "that record and leaves the earlier files on disk."
                if already else "")
        return self.say(f"Trimming {len(tasks)} read pair(s) with fastp {capability['version']}."
                        + note)

    def _launch_trimming(self, tasks, options, output_root):
        from wmlstudio.scheduler import resources_for_run, run_bounded
        project = self.host.project
        try:
            allocation = resources_for_run({}, memory_gb=1)
        except ValueError as error:
            self.host.notify(str(error))
            return False

        def operation(cancelled, progress):
            trimmed = []

            def trim(task, resources, stopped, report):
                primary, mate = task["primary"], task["mate"]
                destination = output_root / primary["id"] / uuid.uuid4().hex
                return read_tools.run_fastp(
                    primary["input_path"], mate["input_path"], destination,
                    threads=resources.threads_per_sample, cancelled=stopped, **options,
                    progress=lambda done, total, message:
                        report(done, total, f"{primary['name']} · {message}"))

            def record(task, result):
                read_tools.record_trimming(project, task["primary"]["id"], task["mate"]["id"],
                                           result)
                trimmed.append({"name": task["primary"]["name"],
                                "reads_removed": result["report"]["reads_removed"]})

            run_bounded(tasks, trim, allocation, cancelled=cancelled, on_result=record,
                        progress=progress)
            return trimmed

        return self.host.launch_task(operation, "read_trimming", self.trimming_finished)

    def trimming_finished(self, trimmed):
        removed = sum(item["reads_removed"] for item in trimmed)
        self.outcome = (f"Trimmed {len(trimmed)} read pair(s); fastp removed {removed:,} reads in "
                        "total. The original FASTQ files are unchanged and still each record's "
                        "input. " + read_tools.TRIMMING_DISCLAIMER)
        self.host.notify(self.outcome)
        self.refresh()
        return self.outcome


class AssemblyPanel(PipelinePanel):
    """Assemble confirmed read pairs, and show what each assembly actually is.

    Every metric is the one the run recorded from the contigs it wrote, and each
    row says which reads went in -- the trimmed pair when one was recorded, the
    originals otherwise -- so an assembly can never imply a preprocessing step
    that did not produce its input.
    """

    COLUMNS = ("Sample", "Assembly", "Contigs", "Total length", "N50", "Largest contig",
               "Reads used", "Assembler", "Assembled at")
    READY_MESSAGE = ("Select the read files of one or more pairs, then assemble them. Nothing is "
                     "typed here: typing runs from the MLST and cgMLST tabs.")

    def build(self, layout):
        strip = FlowLayout()
        self.assemble_button = button("Assemble read pairs…", self.assemble_selected, True)
        self.assemble_button.setToolTip("Assemble the confirmed pairs with the bundled assembler. "
                                        "Where a pair was trimmed, the trimmed reads are used and "
                                        "the run records that it used them.")
        strip.addWidget(self.assemble_button)
        self.measure_button = button("Read metrics from the assembly file", self.measure_selected)
        self.measure_button.setToolTip("Compute contiguity from the FASTA on disk, for an assembly "
                                       "recorded before per-contig metrics were kept. Nothing is "
                                       "written back to the project.")
        strip.addWidget(self.measure_button)
        strip.addWidget(button("Refresh", self.refresh))
        layout.addLayout(strip)
        layout.addWidget(label(ASSEMBLY_DISCLAIMER, "small", True))
        self.build_table(layout, "Select a row to read the metrics behind it.")

    # --- the table ------------------------------------------------------------
    def refresh(self):
        """One row per sample, from what each run recorded."""
        samples = self.samples()
        # One header per unanalysed record separates reads from an assembly. On a
        # large cohort that is a read per row, so the view stays off disk and the
        # rows it cannot classify say "Input not read" instead of guessing.
        inspect = len(samples) <= INSPECT_LIMIT
        self.views, rows = {}, []
        for sample in samples:
            view = assembly_view(sample, inspect_input=inspect)
            self.views[sample["id"]] = view
            rows.append(self._row(sample, view))
        self.fill(rows)
        assembled = sum(view["status"] == "assembled" for view in self.views.values())
        if not rows:
            return self.said_with_outcome(
                "No samples in this project yet. Import a paired FASTQ set in Samples, then "
                "assemble the pairs here.")
        summary = f"{assembled} of {len(rows)} samples carry an assembly made here."
        if not inspect:
            summary += (f" More than {INSPECT_LIMIT} records, so inputs that have never been "
                        "analysed were not read from disk and are listed as not read.")
        return self.said_with_outcome(summary + " " + self.READY_MESSAGE)

    def _row(self, sample, view):
        identifier = sample["id"]
        metrics = view["metrics"] or {}
        source = view["read_source"] or {}
        notes = "\n".join(view["notes"])
        return [
            text_cell(sample.get("name"), identifier, tooltip=notes),
            text_cell(ASSEMBLY_STATUS.get(view["status"], view["status"]), identifier,
                      tooltip=str(view["assembly_path"] or notes)),
            number_cell(metrics.get("contigs"), identifier,
                        tooltip="Contigs in this assembly, counted from the FASTA it wrote."),
            number_cell(metrics.get("total_length"), identifier, unit="bp",
                        tooltip="Assembled bases. This is not a genome size."),
            number_cell(metrics.get("n50"), identifier, unit="bp",
                        tooltip="Half of the assembled length lies in contigs at least this long."),
            number_cell(metrics.get("largest_contig"), identifier, unit="bp",
                        tooltip="The longest single contig in this assembly."),
            text_cell(READ_SOURCE_WORDS.get(source.get("kind"), source.get("kind")), identifier,
                      tooltip=str(source.get("description") or
                                  "No read source was recorded for this record.")),
            text_cell(view["engine"], identifier,
                      tooltip="The assembler this run used, as it reported itself."),
            text_cell(_short(view["completed_at"]), identifier,
                      tooltip="When this assembly finished, in UTC."),
        ]

    # --- the metrics behind one row -------------------------------------------
    def show_detail(self, measured=None):
        view = self.selected_view()
        if view is None:
            self.detail.setPlainText("Select a row to read the metrics behind it.")
            return None
        lines = [f"{view['sample_name']} · "
                 f"{ASSEMBLY_STATUS.get(view['status'], view['status'])}"]
        if view["assembly_path"]:
            lines.append(f"Assembly  {view['assembly_path']}")
        if view["engine"]:
            lines.append(f"Assembler {view['engine']} · finished {_short(view['completed_at'])}")
        metrics = measured or view["metrics"]
        if metrics:
            if measured:
                lines.append("Metrics read from the FASTA just now; nothing was written back to "
                             "this project.")
            lines += self._metric_lines(metrics)
        lines += self._source_lines(view)
        coverage = view["coverage"] or {}
        ratio = coverage.get("read_bases_per_assembled_base")
        if _measured(ratio):
            lines += ["", f"Read bases per assembled base  {ratio:.1f}",
                      "  " + str(coverage.get("denominator", "")),
                      "  " + str(coverage.get("basis", ""))]
        pairing = view["pairing"] or {}
        if pairing:
            lines.append("")
            lines.append("Pairing   " + ("explicit mate markers in the read identifiers"
                                         if pairing.get("explicit_mates")
                                         else "user-assigned; the reads carried no mate markers"))
        lines += ["", "What this does and does not mean:"]
        lines += [f"  · {note}" for note in view["notes"]]
        self.detail.setPlainText("\n".join(lines))
        return view

    def _metric_lines(self, metrics):
        """Exactly what the run measured. A metric it did not keep is named, not filled in."""
        depth = metrics.get("assembler_depth") or {}
        lines = [""]
        for title, key, unit in (("Contigs", "contigs", ""), ("Total length", "total_length", " bp"),
                                 ("N50", "n50", " bp"), ("Largest contig", "largest_contig", " bp"),
                                 ("Smallest contig", "smallest_contig", " bp")):
            value = metrics.get(key)
            lines.append(f"{title:<15}" + (f"{value:,}{unit}" if _measured(value)
                                           else "not recorded for this assembly"))
        if _measured(metrics.get("gc_percent")):
            lines.append(f"GC             {metrics['gc_percent']:.1f}% of "
                         f"{metrics.get('gc_denominator', 'ACGT bases')}")
        ambiguous, unknown_bases = metrics.get("ambiguous_bases"), metrics.get("n_bases")
        if _measured(ambiguous) and _measured(unknown_bases):
            lines.append(f"Ambiguous      {ambiguous:,} bases, of which {unknown_bases:,} are N")
        if depth.get("reported"):
            lines.append(f"Assembler depth  mean {depth['length_weighted_mean']:.1f}× "
                         f"(length-weighted) · min {depth['minimum']:.1f}× · "
                         f"max {depth['maximum']:.1f}× · "
                         f"{depth.get('circular_contigs') or 0} contig(s) reported circular")
        lines.append("  " + str(depth.get("basis", "")))
        lines.append("  " + str(metrics.get("interpretation", "")))
        return lines

    def _source_lines(self, view):
        source = view["read_source"] or {}
        if not source:
            return []
        lines = ["", "Reads used     "
                 + READ_SOURCE_WORDS.get(source.get("kind"), str(source.get("kind"))),
                 "  " + str(source.get("description", ""))]
        removed_bases = source.get("bases_removed")
        if _measured(source.get("reads_removed")):
            lines.append(f"  {source['reads_removed']:,} reads"
                         + (f" and {removed_bases:,} bases" if _measured(removed_bases) else "")
                         + " were removed before assembly.")
        for original in source.get("originals") or ():
            lines.append(f"  original mate {original.get('mate')}  {original.get('path')}")
        return lines

    def measure_selected(self):
        """Read contiguity from the FASTA on disk, without writing anything back."""
        view = self.selected_view()
        if view is None:
            return self.say("Select an assembled sample first.")
        path = view.get("assembly_path")
        if not path or not Path(path).is_file():
            return self.say(f"{view['sample_name']} has no assembly file on this computer, so "
                            "there is nothing to measure.")

        def operation(cancelled, progress):
            progress(10, 100, f"Reading {Path(path).name}")
            return assembly_metrics(path, cancelled=cancelled)

        if not self.host.launch_task(operation, "assembly_metrics",
                                     lambda metrics: self.show_detail(measured=metrics)):
            return self.say("A background task is already running. Let it finish or cancel it, "
                            "then measure again.")
        return self.say(f"Reading contiguity metrics from {Path(path).name}. Nothing is written "
                        "back to this project.")

    # --- running the assembler -------------------------------------------------
    def assemble_selected(self):
        """Assemble the confirmed pairs, each from the reads it is actually given."""
        records = self.candidates()
        pairs = self.pairs_from(records, verb="assemble")
        if not pairs:
            return False
        project = self.host.project
        tasks = [{"primary": project.get_sample(pair["primary_id"]),
                  "mate": project.get_sample(pair["mate_id"])} for pair in pairs]
        self.outcome = ""   # the previous run's sentence is not this run's result
        if not self._launch_assembly(tasks, self.output_root("assemblies")):
            return self.say("A background task is already running. Let it finish or cancel it, "
                            "then assemble again.")
        return self.say(f"Assembling {len(tasks)} read pair(s). Each run records whether it used "
                        "trimmed or original reads.")

    def _launch_assembly(self, tasks, output_root):
        from wmlstudio.assembly import associate_assembly, run_skesa
        from wmlstudio.read_tools import reads_for_assembly
        from wmlstudio.scheduler import resources_for_run, run_bounded
        project = self.host.project
        try:
            allocation = resources_for_run({}, memory_gb=8)
        except ValueError as error:
            self.host.notify(str(error))
            return False

        def operation(cancelled, progress):
            assembled, unattached = [], []

            def assemble(task, resources, stopped, report):
                primary, mate = task["primary"], task["mate"]
                destination = output_root / primary["id"] / uuid.uuid4().hex
                # The trimmed pair when one was recorded, the originals otherwise;
                # either way the run records which reads it actually used.
                read1, read2, source = reads_for_assembly(primary, mate)
                return run_skesa(
                    read1, read2, destination, threads=resources.threads_per_sample,
                    memory_gb=resources.memory_gb, read_source=source, cancelled=stopped,
                    progress=lambda done, total, message:
                        report(done, total, f"{primary['name']} · {message}"))

            def attach(task, result):
                primary = task["primary"]
                entry = {"name": primary["name"], "path": result.get("assembly_path"),
                         "from_trimmed": (result.get("read_source") or {}).get("kind")
                         == "fastp_trimmed"}
                try:
                    associate_assembly(project, primary["id"], task["mate"]["id"], result)
                except ValueError as refusal:
                    # The contigs are written and on disk. Dropping them because the
                    # record would not take them costs the user the whole run, so the
                    # assembly is kept, named, and reported as unattached -- which is
                    # a different outcome from an assembly that was never made.
                    unattached.append({**entry, "reason": str(refusal)})
                else:
                    assembled.append(entry)

            run_bounded(tasks, assemble, allocation, cancelled=cancelled, on_result=attach,
                        progress=progress)
            return {"assembled": assembled, "unattached": unattached}

        return self.host.launch_task(operation, "tab_assembly", self.assembly_finished)

    def assembly_finished(self, result):
        assembled, unattached = result["assembled"], result["unattached"]
        trimmed = sum(item["from_trimmed"] for item in assembled)
        # "Assembled" and "attached to its sample" are two different things, and a
        # run that did the first but not the second must not be counted as either.
        message = (f"{len(assembled)} read pair(s) were assembled and attached to their samples; "
                   f"{trimmed} used trimmed reads and {len(assembled) - trimmed} used the reads "
                   "as supplied. Nothing was typed: type these isolates from the MLST or cgMLST "
                   "tab. " + ASSEMBLY_DISCLAIMER)
        if unattached:
            message += (f" {len(unattached)} assembly(ies) finished but could not be attached to "
                        "their samples, so they are not in this project: "
                        + " · ".join(f"{item['name']} — {item['reason']} The contigs are kept at "
                                     f"{item['path']}." for item in unattached))
        self.outcome = message
        self.host.notify(message)
        self.refresh()
        return message


def install_pipeline_panels(window):
    """Give the Read QC and Assembly stations their real pages.

    Returns the panels that were adopted. A station whose panel cannot be built
    keeps its own page, which says plainly that the work does not run there.
    """
    installed = {}
    for key, factory, attribute in (("reads", ReadTrimmingPanel, "read_trimming_page"),
                                    ("assembly", AssemblyPanel, "assembly_page")):
        try:
            panel = factory(window)
        except TypeError:
            continue
        setattr(window, attribute, panel)
        if window.adopt_station(key, panel):
            installed[key] = panel
    return installed
