"""The SNP tree tab: reference-free split k-mer SNP distances, on a scale of their own.

This page runs SKA2 over a chosen cohort of assemblies and shows three things a
reader needs together: who was compared, how many single-nucleotide differences
separate each pair, and how much sequence that count was measured over. The
denominator is never left off. Every pair is compared over the split k-mers those
two isolates happen to share, so each row carries its own, and a pair that shared
too little has an *unknown* distance -- no edge, no cell, no number -- rather than
a short one.

A SNP distance, a classical seven-locus allele distance and a cgMLST target
distance are three different quantities. Nothing here shares a scale, an axis, a
column or a threshold with an allele tree, the words on this page and in every
window and export opened from it come from :mod:`wmlstudio.snp_tree`, and no
published SNP cutoff is offered unless the protocol it was measured on matches
the run in front of the reader.

The picture this page draws is a minimum spanning tree of those distances. It is
not a maximum-likelihood phylogeny, nothing here infers one, and the page says so
beside the picture rather than in a help file. What it does hand over is the
cohort alignment SKA2 wrote, named with its full path, because that is the file a
tree-building program needs and this application ships none.

All of the science is in :mod:`wmlstudio.snp_tree` and :mod:`wmlstudio.ska_runtime`.
This module arranges what they return and re-words none of it.
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt
from PySide6.QtGui import QColor, QGuiApplication, QImage, QPainter
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QHeaderView,
    QSizePolicy,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import snp_tree
from wmlstudio.assembly import assembly_view
from wmlstudio.graph_window import EXPORT_FORMATS, GraphIdentity, GraphWindow
from wmlstudio.sample_workflow import current_input_sha256
from wmlstudio.ska_runtime import runtime_capabilities
from wmlstudio.theme import BACKGROUND
from wmlstudio.ui_common import FlowLayout, cell, make_table, organism_for
from wmlstudio.widgets import TreeView, button, card, label

#: The run bounds SKA2 itself enforces, stated here so the page can say why an
#: isolate is not in the cohort before a worker starts rather than after it fails.
MINIMUM_COHORT, MAXIMUM_COHORT = 2, 200

#: Beyond this many isolates the square matrix is not drawn: 200 x 200 cells is
#: 40,000 widgets nobody reads. The pair table carries every pair whatever the
#: cohort size, and no number is withheld -- only the grid is.
MATRIX_LIMIT = 60

#: Why an isolate cannot take part, in the words the cohort table prints. Each is
#: a statement about evidence, never a judgement about the isolate.
NOT_COMPARABLE = {
    "read_mate": "The second mate of a read pair; its assembly belongs to the sample it is paired "
                 "with, which is the record that takes part.",
    "no_input": "No input file is recorded for this sample.",
    "unverified": "No reviewed input fingerprint. Analyse this isolate first: SKA2 compares the "
                  "exact file an analysis was reviewed against, never an unverified one.",
    "reads": "Reads, not an assembly. Assemble this isolate first; split k-mers are built from "
             "assembled sequence here.",
}

PURPOSE = ("Single-nucleotide differences between assemblies, counted by SKA2 over the split "
           "k-mers each pair of isolates shares. Reference-free: no reference genome is chosen, "
           "and nothing is mapped.")

#: Printed where the tree is, because a picture separates from its caption and
#: this is the sentence that must never be separated from this one.
NOT_A_PHYLOGENY = ("A minimum spanning forest is a layout of pairwise SNP distances. It is not a "
                   "phylogeny, not a time line and not a transmission chain, and the length of a "
                   "line on the screen carries no meaning at all.")

#: The project setting holding this page's arrangement. It follows the key shape
#: the comparison page uses for its own trees, and it is a key of its own: a SNP
#: forest and an allele forest hold different isolates in different places, and
#: one must never be laid out on the other's saved coordinates.
GRAPH_STYLE_SETTING = "graph_style.snp"

#: Said where the tree is drawn, beside the picture rather than in a help page.
#: A reader who wants a maximum-likelihood tree is told here what this is not and
#: exactly which file to take elsewhere to get one.
NOT_A_MAXIMUM_LIKELIHOOD_TREE = snp_tree.DRAWN_TREE + " " + snp_tree.ML_TREE_ROUTE

#: What the alignment line says before any run has produced one. It states the
#: absence as an absence: no run, therefore no file, never "no alignment exists".
NO_ALIGNMENT_YET = ("No cohort alignment has been written yet. Run the SNP cohort and the file SKA2 "
                    "writes will be named here, with its full path, so it can be taken to a "
                    "tree-building program.")

CLEARED = ("Cleared. The SKA2 run's own output files are still on disk where they were written; "
           "every sample, assembly and stored result in this project is untouched.")

#: What the empty view says it is. A tab with nothing drawn still has to name its
#: quantity: the words belong to the method, not to a particular measurement, and
#: an empty SNP view must never ask a reader for allele profiles.
EMPTY_SCALE = {"kind": snp_tree.KIND, "title": snp_tree.TITLE,
               "target_word": snp_tree.TARGET_WORD, "difference_word": snp_tree.DIFFERENCE_WORD,
               "distance_phrase": snp_tree.DISTANCE_PHRASE, "targets": 0,
               "separation": snp_tree.SEPARATION}


class SnpTreePanel(QWidget):
    """One cohort's SKA2 SNP distances: the isolates, the matrix, the pairs, the tree.

    The host window supplies the project and the background task runner:
    ``project``, ``launch_task`` and ``notify``. ``project_path``, ``cohort_ids``,
    ``check_output`` and ``root`` are used when the host has them, so this panel
    can be mounted as a pipeline station or built on its own without changing.
    """

    REQUIRED = ("project", "launch_task", "notify")

    def __init__(self, window, parent=None):
        missing = [name for name in self.REQUIRED if not hasattr(window, name)]
        if missing:
            raise TypeError("A SNP tree panel needs a workspace window providing "
                            + ", ".join(missing) + ".")
        super().__init__(parent if parent is not None else window)
        self.host = window
        # The SKA2 run this page is showing, and the payload built from it. The run
        # is kept so a different comparability floor is a re-read of numbers SKA2
        # already reported, never a second pass over the sequences.
        self.result = None
        self.payload = None
        # True only while a kept arrangement is being put back on the tree.
        self._restoring = False
        self.cohort = []
        # The cohort alignment this run wrote, when it wrote one: the file a
        # reader takes to a tree builder, since nothing here infers a phylogeny.
        self.alignment_path = None
        self._windows = []
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self._build_header(layout)
        self._build_controls(layout)
        self._build_tabs(layout)
        self._build_footer(layout)
        self.refresh_cohort()

    # --- construction --------------------------------------------------------
    def _build_header(self, layout):
        frame, inner = card()
        inner.setContentsMargins(16, 12, 16, 12)
        inner.setSpacing(6)
        inner.addWidget(label(snp_tree.TITLE, "badge"))
        inner.addWidget(label(snp_tree.QUESTION, "cardTitle", wrap=True))
        inner.addWidget(label(PURPOSE, "muted", wrap=True))
        inner.addWidget(label(snp_tree.SEPARATION, "small", wrap=True))
        self.engine = label("", "small", wrap=True)
        inner.addWidget(self.engine)
        layout.addWidget(frame)

    def _build_controls(self, layout):
        strip = FlowLayout()
        self.scope = QComboBox()
        # Every isolate first: a SNP cohort is chosen here, and a tab that opened
        # empty because nothing happened to be selected elsewhere reads as broken.
        self.scope.addItem("Every isolate in this project", "all")
        self.scope.addItem("The isolates selected elsewhere", "cohort")
        self.scope.setToolTip("Which isolates this page offers to compare. Only assemblies with a "
                              "reviewed input fingerprint take part; the rest are listed with the "
                              "reason they cannot.")
        self.scope.currentIndexChanged.connect(self.refresh_cohort)
        strip.addWidget(self.scope)
        self.run_button = button("Run the SNP cohort…", self.run_tree, True)
        strip.addWidget(self.run_button)
        self.window_button = button("Open the tree in a window", self.open_window)
        self.window_button.setToolTip("A movable, editable copy of this forest, stating its own "
                                      "method, cohort and link threshold.")
        self.window_button.setEnabled(False)
        strip.addWidget(self.window_button)
        # A figure for a report should not require a second window to be opened
        # and arranged first: the picture on this page is already the evidence.
        self.picture_button = button("Export picture…", self.export_picture)
        self.picture_button.setToolTip("Save this forest as a picture or a graph file, with the "
                                       "method, cohort and grouping printed on it.")
        self.picture_button.setEnabled(False)
        strip.addWidget(self.picture_button)
        strip.addWidget(button("Clear", self.clear))
        layout.addLayout(strip)
        settings = FlowLayout()
        settings.addWidget(label("Comparable when the pair shares at least", "small"))
        self.floor = QDoubleSpinBox()
        self.floor.setRange(0.0, 1.0)
        self.floor.setSingleStep(0.05)
        self.floor.setDecimals(3)
        self.floor.setValue(0.95)
        self.floor.setToolTip("The comparability floor. Lowering it admits pairs compared over less "
                              "shared sequence; a distance is not better evidence for having been "
                              "admitted, and no SNP count changes either way.")
        settings.addWidget(self.floor)
        self.basis = QComboBox()
        for key, words in snp_tree.FRACTION_BASES.items():
            self.basis.addItem("of " + words["words"], key)
            self.basis.setItemData(self.basis.count() - 1, words["reading"], Qt.ItemDataRole.ToolTipRole)
        self.basis.currentIndexChanged.connect(self._basis_changed)
        settings.addWidget(self.basis)
        settings.addWidget(label("· group at", "small"))
        self.link = QSpinBox()
        self.link.setRange(snp_tree.NO_LINK_THRESHOLD, 100000)
        self.link.setValue(snp_tree.NO_LINK_THRESHOLD)
        self.link.setSpecialValueText("no grouping")
        self.link.setSuffix(" SNPs")
        self.link.setToolTip("Single-link grouping for this picture only. It is a choice made in "
                             "this view, never a validated cutoff, and nothing is grouped until "
                             "one is chosen.")
        settings.addWidget(self.link)
        self.reread_button = button("Re-read with these settings", self.apply_settings)
        self.reread_button.setToolTip("Re-reads the same SKA2 run at this floor and this grouping. "
                                      "No sequence is read again and no SNP count moves.")
        self.reread_button.setEnabled(False)
        settings.addWidget(self.reread_button)
        layout.addLayout(settings)
        self.basis_note = label(snp_tree.FRACTION_BASES["combined"]["reading"], "small", wrap=True)
        layout.addWidget(self.basis_note)

    def _build_tabs(self, layout):
        self.tabs = QTabWidget()
        self.tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.tabs.addTab(self._build_cohort_tab(), "Cohort")
        self.tabs.addTab(self._build_matrix_tab(), "SNP distances")
        self.tabs.addTab(self._build_pairs_tab(), "Pairs and denominators")
        self.tabs.addTab(self._build_tree_tab(), "Tree")
        layout.addWidget(self.tabs, 1)

    def _build_cohort_tab(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        column.addWidget(label("Which isolates this comparison can include, and why an isolate "
                               "cannot. An isolate that is not compared is unknown evidence here; "
                               "it is never counted as identical to anything.", "small", True))
        self.cohort_table = make_table(["Isolate", "Organism", "Assembly", "In this comparison",
                                        "Split k-mers held", "Compared with"])
        column.addWidget(self.cohort_table, 1)
        return page

    def _build_matrix_tab(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        self.matrix_note = label("Run a cohort to see its SNP distances.", "small", True)
        column.addWidget(self.matrix_note)
        self.matrix_table = QTableWidget(0, 0)
        self.matrix_table.setAlternatingRowColors(True)
        self.matrix_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        # Deliberately unsortable: a SNP matrix is square and symmetric, and sorting
        # one axis of it would put each row against the wrong column.
        self.matrix_table.setSortingEnabled(False)
        self.matrix_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.matrix_table.horizontalHeader().setDefaultSectionSize(120)
        column.addWidget(self.matrix_table, 1)
        self.matrix_legend = label("", "small", True)
        column.addWidget(self.matrix_legend)
        return page

    def _build_pairs_tab(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        column.addWidget(label("Every pair, with the sequence its distance was measured over. A "
                               "pair that did not share enough has no distance; what SKA2 observed "
                               "over the little it did share is shown beside it and is not a "
                               "distance.", "small", True))
        self.pairs_table = make_table(["Isolate", "Compared with", "SNP distance",
                                       "Observed by SKA2", "Measured over", "Status"])
        column.addWidget(self.pairs_table, 1)
        self.comparability = label("", "small", True)
        column.addWidget(self.comparability)
        return page

    def _build_tree_tab(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        self.tree_caption = label("No SNP forest is drawn yet.", "small", True)
        column.addWidget(self.tree_caption)
        self.tree = TreeView()
        self.tree.setMinimumHeight(240)
        self.tree.show_contents({"scale": dict(EMPTY_SCALE)})
        # Connected after the empty view is drawn, not before: showing contents is
        # itself a layout change, and saving that one would have overwritten a
        # reader's kept arrangement with an empty forest every time the page was
        # built. Without these three the page's own tree emitted its changes into
        # nothing, so dragging, recolouring or renaming a node survived only until
        # the forest was next drawn.
        for signal in ("layoutChanged", "colorsChanged", "labelsChanged"):
            getattr(self.tree, signal).connect(lambda *_ignored: self.store_arrangement())
        column.addWidget(self.tree, 1)
        column.addWidget(label(NOT_A_PHYLOGENY, "small", True))
        column.addWidget(label(NOT_A_MAXIMUM_LIKELIHOOD_TREE, "small", True))
        self.alignment_note = label(NO_ALIGNMENT_YET, "small", True)
        # The path is worth selecting with the mouse: it is meant to be pasted
        # into the command line of a tree builder this application does not run.
        self.alignment_note.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        column.addWidget(self.alignment_note)
        self.alignment_button = button("Copy the alignment file path", self.copy_alignment_path)
        self.alignment_button.setToolTip("Copies the path of the cohort alignment SKA2 wrote. That "
                                         "file is the input a maximum-likelihood tree is built "
                                         "from; this page draws no such tree.")
        self.alignment_button.setEnabled(False)
        column.addWidget(self.alignment_button)
        return page

    def _build_footer(self, layout):
        self.threshold_note = label("", "small", True)
        layout.addWidget(self.threshold_note)
        self.limitations = label("", "small", True)
        layout.addWidget(self.limitations)
        self.status = label("", "small", True)
        layout.addWidget(self.status)

    # --- the cohort ----------------------------------------------------------
    def say(self, message):
        self.status.setText(str(message))
        return message

    def cohort_ids(self):
        """The isolates this page offers, or None for every isolate in the project."""
        if self.scope.currentData() == "all":
            return None
        chosen = getattr(self.host, "cohort_ids", None)
        return None if chosen is None else set(chosen)

    @staticmethod
    def _eligibility(sample):
        """Whether one isolate can take part, and the reason when it cannot."""
        workflow = (sample.get("metadata") or {}).get("workflow") or {}
        if workflow.get("source_kind") == "read_mate":
            return False, "read_mate", "read mate"
        if not sample.get("input_path"):
            return False, "no_input", "no input file"
        status = assembly_view(sample, inspect_input=True)["status"]
        if status == "reads_not_assembled":
            return False, "reads", "reads, not assembled"
        if not current_input_sha256(sample):
            return False, "unverified", "fingerprint not reviewed"
        return True, "", {"assembled": "assembled here", "assembly_supplied": "supplied assembly",
                          "read_mate": "read mate"}.get(status, "input kind not established")

    def refresh_cohort(self):
        """List every isolate this page can offer, with what it is and whether it counts."""
        chosen = self.cohort_ids()
        rows = []
        for sample in self.host.project.samples():
            if chosen is not None and sample["id"] not in chosen:
                continue
            eligible, reason, assembly = self._eligibility(sample)
            genus, species, source = organism_for(sample)
            rows.append({"id": sample["id"], "name": sample["name"],
                         "input_path": sample.get("input_path"),
                         "input_sha256": current_input_sha256(sample),
                         "organism": " ".join(filter(None, [genus, species])) or "not established",
                         "organism_source": source, "assembly": assembly,
                         "eligible": eligible, "reason": NOT_COMPARABLE.get(reason, "")})
        rows.sort(key=lambda row: (not row["eligible"], str(row["name"]).casefold(), row["id"]))
        self.cohort = rows
        self._fill_cohort_table()
        count = sum(1 for row in rows if row["eligible"])
        capability = runtime_capabilities()
        self.engine.setText(
            f"SKA2 {capability['version']} ({capability['platform']}), bundled with this "
            "application. Reference-free split k-mer SNPs between assemblies."
            if capability["available"] else capability["reason"])
        self.run_button.setEnabled(capability["available"] and MINIMUM_COHORT <= count)
        if not capability["available"]:
            self.say("The native SNP engine is not available in this build, so nothing on this tab "
                     "can run. Nothing was downloaded and no distance was estimated another way.")
        elif count < MINIMUM_COHORT:
            self.say(f"{count} of {len(rows)} isolate(s) in this selection can be compared; at "
                     f"least {MINIMUM_COHORT} are needed. The table says why each of the others "
                     "cannot, and none of them is counted as similar to anything.")
        elif self.payload is None:
            self.say(f"{count} of {len(rows)} isolate(s) in this selection can be compared. Run the "
                     "SNP cohort to measure them.")
        return rows

    def _fill_cohort_table(self):
        held = dict((self.payload or {}).get("split_kmers", {}).get("per_sample") or {})
        summary = (self.payload or {}).get("summary") or {}
        alone = set(summary.get("isolates_without_a_comparable_pair") or ())
        compared = {row["sample_id"] for row in (self.payload or {}).get("cohort") or ()}
        table = self.cohort_table
        table.setSortingEnabled(False)
        table.setRowCount(len(self.cohort))
        for index, row in enumerate(self.cohort):
            if row["id"] in alone:
                partners = "no pair shared enough sequence"
            elif row["id"] in compared:
                partners = "compared in this run"
            else:
                partners = "not in this run"
            values = [row["name"], f"{row['organism']} · {row['organism_source'].casefold()}",
                      row["assembly"],
                      "yes" if row["eligible"] else "no — " + row["reason"],
                      f"{held[row['id']]:,}" if row["id"] in held else "—", partners]
            for column, value in enumerate(values):
                table.setItem(index, column, cell(value, row["id"]))
        table.setSortingEnabled(True)

    def cohort_name(self):
        """What this cohort is, for a caption and a window title.

        The isolate count is the scale's to state, so this says which selection was
        compared rather than repeating a number beside itself.
        """
        return self.scope.currentText()

    def cohort_organism(self):
        """The one organism this cohort is, or '' when it is not one organism.

        A published cutoff belongs to an organism, so a cohort of two species has
        no organism to ask about and says so instead of borrowing one of them.
        """
        names = {row["organism"] for row in self.cohort if row["eligible"]} - {"not established"}
        return next(iter(names)) if len(names) == 1 else ""

    def display_records(self):
        """Names, sequence types and recorded detail for the drawn nodes.

        The classical sequence type is display detail here and nothing more: it is
        a different quantity from anything this page measures, so it labels a node
        and never enters a distance, an edge or a group.
        """
        from wmlstudio.ui_compare import _available_profiles
        records = {}
        for sample in self.host.project.samples():
            profiles = _available_profiles(self.host.project, sample, "mlst")
            st = next((str(row["st"]) for row in profiles if row.get("st")), None)
            records[sample["id"]] = {"sample_name": sample["name"], "st": st, "primary_st": st,
                                     "metadata": sample.get("metadata") or {}}
        return records

    # --- running -------------------------------------------------------------
    def output_root(self):
        """Where SKA2's own output goes: beside the project's files, never in an input folder."""
        path = getattr(self.host, "project_path", None)
        base = Path(path).with_suffix(".files") if path else Path(self.host.project.path).parent
        return base / "snp"

    def run_tree(self):
        """Measure the cohort with SKA2, in the background, and draw what it returns."""
        capability = runtime_capabilities()
        if not capability["available"]:
            return self.say(capability["reason"])
        chosen = [row for row in self.cohort if row["eligible"]]
        if len(chosen) < MINIMUM_COHORT:
            return self.say(f"Choose at least {MINIMUM_COHORT} isolates with a reviewed assembly.")
        if len(chosen) > MAXIMUM_COHORT:
            return self.say(f"{len(chosen)} isolates were selected; this analysis is bounded at "
                            f"{MAXIMUM_COHORT}. Narrow the selection rather than comparing a "
                            "cohort nobody can read.")
        from wmlstudio.scheduler import plan_resources
        from wmlstudio.ska_runtime import run_ska
        samples = [{"id": row["id"], "name": row["name"], "input_path": row["input_path"],
                    "input_sha256": row["input_sha256"]} for row in chosen]
        organism, records = self.cohort_organism(), self.display_records()
        cohort_name = self.cohort_name()
        output_root, floor, basis = self.output_root(), self.floor.value(), self.basis.currentData()
        link = self.link.value()
        resources = plan_resources(threads_per_sample=4, memory_gb=2)

        def analyse(cancelled, report):
            result = run_ska(samples, output_root, threads=resources.threads_per_sample,
                             min_shared_fraction=floor, cancelled=cancelled, progress=report)
            if basis != "combined":
                result = snp_tree.at_minimum_shared_fraction(result, floor, basis=basis)
            payload = snp_tree.snp_payload(result, organism=organism, records=records,
                                           link_threshold=None if link < 0 else link,
                                           cohort=cohort_name, cancelled=cancelled)
            _write_payload(result, payload)
            return {"result": result, "payload": payload}

        if self.host.launch_task(analyse, "snp", self._run_finished) is False:
            return self.say("A background task is already running; finish or cancel it first.")
        worker = getattr(self.host, "worker", None)
        if worker is not None:
            # The same sentence the message box shows, kept on this page so the
            # status line never reads "running" after a run has stopped.
            worker.failed.connect(self.say)
        return self.say(f"SKA2 is comparing {len(samples)} isolate(s) at a floor of {floor:g} "
                        + str(snp_tree.FRACTION_BASES[basis]["words"])
                        + ". Nothing is written back onto any sample.")

    def _run_finished(self, bundle):
        bundle = dict(bundle or {})
        self.result, self.payload = bundle.get("result"), bundle.get("payload")
        if self.payload is None:
            return self.say("The SNP run returned nothing to show.")
        self.show_payload()
        summary = self.payload["summary"]
        return self.say(
            f"{summary['isolates']} isolate(s) · {summary['comparable_pairs']} of "
            f"{summary['pairs']} pairs carry a SNP distance · {summary['excluded_pairs']} shared "
            "too little sequence to be compared · written beside the run's own output in "
            f"{self.payload['output_directory']}.")

    def show_result(self, result):
        """Draw one SKA2 result: the same arrangement a finished run produces.

        Separate from :meth:`run_tree` so a result already on disk, or one a caller
        holds, reaches the page without the native engine being run again.
        """
        self.result = result
        self.payload = snp_tree.snp_payload(
            result, organism=self.cohort_organism(), records=self.display_records(),
            link_threshold=None if self.link.value() < 0 else self.link.value(),
            cohort=self.cohort_name())
        self.show_payload()
        return self.payload

    def apply_settings(self):
        """Re-read the run this page is showing at the chosen floor and grouping.

        SKA2 is not run again: which pairs count as comparable is decided from the
        shared and unshared counts the run already reported, and no SNP count moves.
        """
        if self.result is None:
            return self.say("Run a cohort first; there is nothing to re-read.")
        floor, basis = self.floor.value(), self.basis.currentData()
        try:
            self.result = snp_tree.at_minimum_shared_fraction(self.result, floor, basis=basis)
            self.payload = snp_tree.snp_payload(
                self.result, organism=self.cohort_organism(), records=self.display_records(),
                link_threshold=None if self.link.value() < 0 else self.link.value(),
                cohort=self.cohort_name())
        except ValueError as error:
            return self.say(str(error))
        self.show_payload()
        reread = self.payload.get("comparability_reread") or {}
        return self.say(f"Re-read at {reread.get('to')} {snp_tree.FRACTION_BASES[basis]['words']} "
                        f"from the counts SKA2 already reported. {reread.get('note', '')}")

    def _basis_changed(self):
        words = snp_tree.FRACTION_BASES.get(self.basis.currentData() or "combined")
        self.basis_note.setText(words["reading"])

    # --- what one run says ---------------------------------------------------
    def show_payload(self):
        """Put one payload on the page: cohort, matrix, pairs, forest and refusals."""
        payload = self.payload or {}
        self._fill_cohort_table()
        self._fill_matrix(payload.get("matrix") or {})
        self._fill_pairs(payload.get("pairs") or [])
        self._draw_forest(payload)
        self.comparability.setText((payload.get("comparability") or {}).get("message", ""))
        binding = payload.get("threshold") or {}
        self.threshold_note.setText(" ".join(filter(None, [
            binding.get("message", ""), binding.get("link_threshold_warning", ""),
            binding.get("notice", "")])))
        self.limitations.setText("\n".join(payload.get("limitations") or ()))
        self._show_alignment(payload.get("alignment_handoff") or {})
        self.reread_button.setEnabled(self.result is not None)
        self.window_button.setEnabled(bool(self.tree.nodes))
        self.picture_button.setEnabled(bool(self.tree.nodes))
        return payload

    def _show_alignment(self, handoff):
        """Name the alignment file SKA2 wrote, or say plainly that there is none.

        This is the only file in the run a phylogenetics program can be fed, and
        it is the whole of this application's answer to "where is the ML tree":
        the input is handed over, and no tree is drawn from it here.
        """
        path = handoff.get("path")
        message = str(handoff.get("message") or "")
        self.alignment_path = path
        self.alignment_note.setText(
            f"Cohort alignment: {path}\n{message}" if path else message or NO_ALIGNMENT_YET)
        self.alignment_button.setEnabled(bool(path))

    def _fill_matrix(self, matrix):
        """The square matrix, each cell carrying what it was measured over.

        A refused pair is left with no distance at all and says so. It is never a
        zero, never an empty cell a reader can take for one, and never sorted or
        coloured as though it were a small number.
        """
        samples = list(matrix.get("samples") or ())
        table = self.matrix_table
        table.clear()
        self.matrix_legend.setText(" ".join(filter(None, [
            matrix.get("missing", ""), matrix.get("observed_meaning", ""),
            matrix.get("diagonal", ""), matrix.get("separation", "")])))
        if not samples:
            table.setRowCount(0)
            table.setColumnCount(0)
            self.matrix_note.setText("Run a cohort to see its SNP distances.")
            return
        if len(samples) > MATRIX_LIMIT:
            table.setRowCount(0)
            table.setColumnCount(0)
            self.matrix_note.setText(
                f"{len(samples)} isolates is more than this grid shows ({MATRIX_LIMIT}). Every "
                "pair, its distance and the sequence it was measured over are in Pairs and "
                "denominators; nothing is withheld, only the square grid.")
            return
        names = [row["sample_name"] for row in samples]
        table.setRowCount(len(samples))
        table.setColumnCount(len(samples))
        table.setHorizontalHeaderLabels(names)
        table.setVerticalHeaderLabels(names)
        for row, source in enumerate(samples):
            for column, target in enumerate(samples):
                distance = matrix["distance"][row][column]
                shared = matrix["shared_split_kmers"][row][column]
                item = QTableWidgetItem("0" if row == column else
                                        str(distance) if distance is not None else "not comparable")
                item.setToolTip("\n".join(filter(None, [
                    f"{source['sample_name']} / {target['sample_name']}",
                    "An isolate is not compared with itself; this zero is by definition."
                    if row == column else
                    f"{distance} {matrix.get('unit', 'SNPs')} over {shared:,} shared split k-mers"
                    if distance is not None else
                    "This pair did not share enough sequence to be compared. Its distance is "
                    "unknown, not zero." + (f" It shared {shared:,} split k-mers." if shared else ""),
                ])))
                table.setItem(row, column, item)
        self.matrix_note.setText(
            f"SNP distances between {len(samples)} isolates, each counted over the split k-mers "
            "that pair shares. Rows and columns are in the same order.")

    def _fill_pairs(self, pairs):
        table = self.pairs_table
        table.setSortingEnabled(False)
        table.setRowCount(len(pairs))
        for index, row in enumerate(pairs):
            distance = (str(row["distance"]) if row["comparable"]
                        else "no accepted distance")
            status = "compared" if row["comparable"] else row.get("reason") or "not comparable"
            values = [row.get("source_name") or row["source"], row.get("target_name") or row["target"],
                      distance, row["observed_snp_count"], row.get("denominator_label", ""), status]
            for column, value in enumerate(values):
                table.setItem(index, column, cell(value, row["source"]))
            table.item(index, 3).setToolTip(
                "What SKA2 saw over whatever this pair shared. It is reported for every pair, "
                "including the pairs whose distance was refused, and is not a distance.")
        table.setSortingEnabled(True)

    def _draw_forest(self, payload):
        graph = payload.get("graph") or {}
        records = graph.get("results") or []
        # Drawing lays the forest out automatically, which the tree reports as a
        # layout change like any other. Saving that would overwrite the reader's
        # own arrangement with the automatic one a moment before it is read back,
        # so nothing is saved until the kept arrangement is in place.
        self._restoring = True
        try:
            self.tree.show_contents(graph)
        finally:
            self._restoring = False
        self._restore_arrangement()
        # "ST unassigned" under every node of a cohort nobody has typed classically
        # reads as a failed typing run. The line is shown only where an ST exists.
        # Applied after the stored arrangement, because a saved state carries its
        # own label fields and they may have been saved on a cohort that had STs.
        self.tree.set_label_fields(["sample_name", "primary_st"]
                                   if any(row.get("st") for row in records) else ["sample_name"])
        summary = payload.get("summary") or {}
        self.tree_caption.setText(" · ".join(filter(None, [
            (graph.get("scale") or {}).get("caption", ""),
            f"{summary.get('edges', 0)} edge(s) between {summary.get('isolates', 0)} isolate(s)",
            f"{summary.get('excluded_pairs', 0)} pair(s) too little compared to be drawn"
            if summary.get("excluded_pairs") else ""])) or "No SNP forest is drawn yet.")

    # --- the arrangement this reader made ------------------------------------
    def graph_identity(self):
        """What this forest is, in the words every window and every export uses.

        One identity for the page, the window and the saved picture: a figure
        that reached a report must carry the same method, cohort and grouping as
        the window it was arranged in, or the two are separate claims.
        """
        contents = self.tree.graph_contents()
        return GraphIdentity.from_scale(
            contents.get("scale"), threshold=contents.get("cluster_threshold"),
            cohort=self.cohort_name(), created=(self.payload or {}).get("created_at"),
            note="SKA2 split k-mer SNPs; not allele differences and not a transmission chain")

    def _window_arranged(self, state):
        """Take an arrangement made in a detached window and keep it.

        Only presentation crosses: positions, colors and display labels. No SNP
        count, denominator, edge or group comes back this way, so a rearranged
        picture is the same measurement in a different place on the screen.
        """
        self._restoring = True
        try:
            self.tree.restore_state(state)
        except ValueError:
            # Presentation state this build does not understand is not evidence;
            # the drawn forest stays exactly as it is rather than half-applying it.
            return None
        finally:
            self._restoring = False
        return self.store_arrangement(state)

    def store_arrangement(self, state=None):
        """Remember this page's arrangement in the project, under its own key.

        Silent while an arrangement is being put back, because restoring one
        makes the tree emit the very signals that ask for it to be saved: the
        arrangement a reader made was read back, re-exported mid-restore, and
        overwritten with the automatic layout it was replacing.
        """
        if self._restoring:
            return None
        state = self.tree.export_state() if state is None else state
        project = getattr(self.host, "project", None)
        if project is not None and hasattr(project, "set_setting"):
            project.set_setting(GRAPH_STYLE_SETTING, state)
        return state

    def _restore_arrangement(self):
        """Put a kept arrangement back onto a freshly drawn forest.

        Isolates the stored state does not name keep the automatic layout, and a
        state from another version is dropped whole: an arrangement is
        presentation, and a broken one must never stop a measurement being drawn.
        """
        project = getattr(self.host, "project", None)
        stored = project.get_setting(GRAPH_STYLE_SETTING, None) if project is not None else None
        if not isinstance(stored, dict) or not stored:
            return None
        self._restoring = True
        try:
            self.tree.restore_state(stored)
        except ValueError:
            return None
        finally:
            self._restoring = False
        return stored

    # --- a window of its own -------------------------------------------------
    def open_window(self):
        """Open this forest in the same movable, editable window the other trees use.

        The window draws a copy and states its own identity, so a SNP forest and an
        allele forest can stand side by side without either borrowing the other's
        numbers or its words.
        """
        if not self.tree.nodes:
            return self.say("Run a cohort first: there is nothing drawn to open in a window.")
        from wmlstudio.interface_settings import interface_preferences
        root = getattr(self.host, "root", None)
        window = GraphWindow.from_view(
            self.tree, self.graph_identity(), parent=self,
            preferences=interface_preferences(root) if root else None)
        window.check_output = getattr(self.host, "check_output", None)
        # An arrangement made in the window is the arrangement the reader wants:
        # it comes back to the page and is kept, so closing the window does not
        # throw the layout away and re-running the cohort does not undo it.
        window.stateChanged.connect(self._window_arranged)
        self._windows.append(window)
        window.closed.connect(lambda: self._windows.remove(window)
                              if window in self._windows else None)
        window.show()
        return window

    # --- a picture of this forest --------------------------------------------
    def copy_alignment_path(self):
        """Put the cohort alignment's path on the clipboard, for a tree builder."""
        path = getattr(self, "alignment_path", None)
        if not path:
            return self.say("This run wrote no cohort alignment, so there is no path to copy. "
                            "Nothing was inferred in its place.")
        clipboard = QGuiApplication.clipboard()
        if clipboard is not None:
            clipboard.setText(str(path))
        return self.say(f"Copied {path}. It is a cohort variable-site alignment, the input a "
                        "maximum-likelihood tree is built from by IQ-TREE, FastTree or RAxML-NG. "
                        "No tree was built here and the picture on this page is not one.")

    def export_picture(self):
        """Ask where to save a figure of this forest, and write it with its own words."""
        if not self.tree.nodes:
            return self.say("Run a cohort first: there is nothing drawn to save as a picture.")
        filters = ";;".join(f"{name} (*.{suffix})" for name, suffix in EXPORT_FORMATS)
        path, chosen = QFileDialog.getSaveFileName(self, "Save this SNP forest",
                                                   f"{snp_tree.KIND}-forest.png", filters)
        if not path:
            return None
        # A typed name without an extension gets the format the dialog was on,
        # rather than a file the computer cannot open again.
        picked = next((suffix for name, suffix in EXPORT_FORMATS
                       if chosen and chosen.startswith(name)), "png")
        target = Path(path)
        target = target if target.suffix else target.with_suffix("." + picked)
        try:
            return self.save_picture(target)
        except Exception as error:  # a full disk or a read-only folder, said plainly
            return self.say(f"Could not save: {error}")

    def save_picture(self, path, *, suffix=None):
        """Write the drawn forest to a file, carrying what it measures with it.

        The same formats and the same title and subtitle as a detached window
        writes, so a figure saved from the page and one saved from a window of it
        are the same picture of the same measurement, not two different claims.
        """
        path = Path(path)
        suffix = str(suffix or path.suffix.lstrip(".")).lower()
        check = getattr(self.host, "check_output", None)
        if callable(check):
            check(path)
        identity = self.graph_identity()
        if suffix in {"png", "jpg", "jpeg"}:
            self.tree.save_image(path, title=identity.export_title(),
                                 subtitle=identity.export_subtitle())
        elif suffix == "svg":
            self.tree.save_svg(path, title=identity.export_title(),
                               subtitle=identity.export_subtitle())
        elif suffix == "graphml":
            self.tree.save_graphml(path)
        elif suffix == "nwk":
            # A Newick file of a minimum spanning tree is a topology, not a
            # phylogram: it carries no inferred branch lengths, because none
            # were estimated.
            self.tree.save_newick(path)
        else:
            raise ValueError(f"Unsupported graph export format: {suffix or 'none given'}")
        self.say(f"Saved {path.name} · {identity.caption()} · {identity.threshold_words()}. "
                 "It is a minimum spanning tree of SNP distances, not a phylogeny.")
        return path

    def snp_graph_image(self, fmt="PNG"):
        """This forest as image bytes, for a report that prints the picture.

        Returns ``None`` when nothing is drawn, so a report can say that no
        picture was produced rather than printing an empty frame that a reader
        could take for a cohort with nothing in it.
        """
        if not self.tree.nodes:
            return None
        identity = self.graph_identity()
        # JPEG carries no alpha channel, so the canvas is opaque from the start;
        # the renderer fills the same background over it either way.
        picture = QImage(1800, 1200, QImage.Format.Format_RGB32)
        picture.fill(QColor(BACKGROUND))
        painter = QPainter(picture)
        try:
            self.tree._render(painter, 1800, 1200, title=identity.export_title(),
                              subtitle=identity.export_subtitle())
        finally:
            painter.end()
        data = QByteArray()
        buffer = QBuffer(data)
        buffer.open(QIODevice.OpenModeFlag.WriteOnly)
        if not picture.save(buffer, fmt, 92 if fmt == "JPEG" else -1):
            raise OSError("Could not render the SNP forest picture.")
        return bytes(data)

    def snp_graph_png(self):
        return self.snp_graph_image("PNG")

    # --- clearing ------------------------------------------------------------
    def clear(self):
        """Empty this page's view. No sample, assembly, result or output file is touched."""
        self.result = self.payload = None
        self.tree.show_contents({"scale": dict(EMPTY_SCALE)})
        self.matrix_table.clear()
        self.matrix_table.setRowCount(0)
        self.matrix_table.setColumnCount(0)
        self.pairs_table.setRowCount(0)
        for widget in (self.matrix_legend, self.comparability, self.threshold_note,
                       self.limitations):
            widget.setText("")
        self.matrix_note.setText("Run a cohort to see its SNP distances.")
        self.tree_caption.setText("No SNP forest is drawn yet.")
        self.alignment_path = None
        self.alignment_note.setText(NO_ALIGNMENT_YET)
        self.alignment_button.setEnabled(False)
        self.link.setValue(snp_tree.NO_LINK_THRESHOLD)
        self.reread_button.setEnabled(False)
        self.window_button.setEnabled(False)
        self.picture_button.setEnabled(False)
        self.refresh_cohort()
        return self.say(CLEARED)


def _write_payload(result, payload):
    """Keep the picture, the matrix and the refusal that produced them together on disk."""
    from wmlstudio.export import _atomic_text
    with _atomic_text(Path(result["output_directory"]) / "snp-tree.json") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, allow_nan=False)
