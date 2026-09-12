"""WMLSTudio native desktop entry point, project workflow and evidence views."""

import argparse
import html
import json
import sqlite3
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QLockFile, Qt, QTimer
from PySide6.QtGui import QAction, QColor, QKeySequence, QPageSize, QPdfWriter, QTextDocument
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import __version__
from wmlstudio.comparison import minimum_spanning_forest, pairwise_distances
from wmlstudio.demo import create_demo
from wmlstudio.export import ensure_separate_destination, export_results
from wmlstudio.jobs import AnalysisWorker, SchemeImportWorker
from wmlstudio.paths import data_root, scheme_locations
from wmlstudio.project import Project
from wmlstudio.theme import STYLE
from wmlstudio.widgets import DropZone, Helix, Metric, TreeView, button, card, label

FILE_FILTER = "Sequence files (*.fasta *.fa *.fna *.fastq *.fq *.gz *.bz2);;All files (*)"
STATUS_TEXT = {
    "queued": "Ready", "running": "Analysing", "completed": "Complete",
    "complete": "Exact match", "novel_profile": "New combination", "mixed": "Mixed alleles",
    "incomplete": "Missing loci", "ambiguous": "Ambiguous", "profile_unavailable": "Profile only",
    "qc_only": "Quality checked", "failed": "Needs attention", "interrupted": "Interrupted",
}


def table(headers):
    widget = QTableWidget(0, len(headers))
    widget.setHorizontalHeaderLabels(headers)
    widget.verticalHeader().hide()
    widget.verticalHeader().setDefaultSectionSize(49)
    widget.setShowGrid(False)
    widget.setAlternatingRowColors(True)
    widget.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    widget.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
    widget.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
    return widget


class BaseWindow(QMainWindow):
    def __init__(self, project_path=None, storage_root=None):
        super().__init__()
        self.root = Path(storage_root) if storage_root else data_root()
        self.root.mkdir(parents=True, exist_ok=True)
        self.project_path = Path(project_path) if project_path else self.root / "My workspace.wmlstudio"
        self.project_lock = self.lock_project(self.project_path)
        self.project = Project(self.project_path)
        self.worker = None
        self.closing_after_cancel = False
        self.scheme_paths = []
        self.current_samples = []
        self.setWindowTitle("WMLSTudio · Microbial genomics workspace")
        self.resize(1380, 940)
        self.setMinimumSize(1000, 680)
        if QApplication.platformName() != "offscreen" and self.screen():
            available = self.screen().availableGeometry()
            self.resize(min(1380, available.width() - 40), min(940, available.height() - 60))
        self.setAcceptDrops(True)
        main = QWidget()
        main.setObjectName("main")
        self.setCentralWidget(main)
        outer = QHBoxLayout(main)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        self.build_sidebar(outer)
        body = QVBoxLayout()
        body.setContentsMargins(33, 25, 33, 18)
        body.setSpacing(18)
        outer.addLayout(body, 1)
        top = QHBoxLayout()
        self.breadcrumb = label("WORKSPACE  /  OVERVIEW", "eyebrow")
        top.addWidget(self.breadcrumb)
        top.addStretch()
        top.addWidget(label("●  Local & private", "badge"))
        top.addWidget(button("Open project", self.open_project_dialog))
        body.addLayout(top)
        self.pages = QStackedWidget()
        body.addWidget(self.pages, 1)
        self.build_overview()
        self.build_samples()
        self.build_compare()
        self.build_schemes()
        self.build_hydra()
        self.build_reports()
        self.build_settings()
        bottom = QHBoxLayout()
        self.progress_text = label("Ready when you are. Your files stay on this computer.", "small")
        bottom.addWidget(self.progress_text, 1)
        self.cancel_button = button("Cancel analysis", self.cancel_analysis)
        self.cancel_button.hide()
        bottom.addWidget(self.cancel_button)
        body.addLayout(bottom)
        self.progress_bar = QProgressBar()
        self.progress_bar.hide()
        body.addWidget(self.progress_bar)
        self.statusBar().showMessage(f"WMLSTudio {__version__} · Research workbench")
        self.statusBar().addPermanentWidget(label("IOWA-BioTech", "small"))
        # Native child viewports must paint directly. A full-stack opacity effect
        # caches child surfaces and caused stale/blank pages on Windows.
        self.pages.setAutoFillBackground(True)
        self.motion_enabled = bool(self.project.get_setting("motion", True))
        self.set_motion(self.motion_enabled)
        self.populate_schemes()
        self.refresh()
        self.navigate(0)
        self.add_shortcuts()

    def build_sidebar(self, outer):
        sidebar = QWidget()
        sidebar.setObjectName("sidebar")
        sidebar.setFixedWidth(218)
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(18, 28, 18, 18)
        layout.setSpacing(6)
        brand = QHBoxLayout()
        monogram = label("W")
        monogram.setAlignment(Qt.AlignmentFlag.AlignCenter)
        monogram.setFixedSize(34, 36)
        monogram.setStyleSheet("background: #197D6C; color: white; border-radius: 10px; font-size: 24px; font-weight: 700;")
        brand.addWidget(monogram)
        brand.addWidget(label("WMLSTudio", "brand"))
        layout.addLayout(brand)
        layout.addWidget(label("MICROBIAL GENOMICS", "eyebrow"))
        layout.addSpacing(30)
        layout.addWidget(label("YOUR WORKSPACE", "eyebrow"))
        self.project_label = label(self.project_path.stem, "cardTitle", True)
        layout.addWidget(self.project_label)
        layout.addSpacing(18)
        self.nav_buttons = []
        self.nav_names = ["Overview", "Samples", "Compare", "Scheme library", "HYDRA insights", "Reports", "Settings and help"]
        symbols = ["◫", "▤", "⌘", "▥", "◈", "↗", "⚙"]
        for index, (name, symbol) in enumerate(zip(self.nav_names, symbols, strict=True)):
            item = button(f"{symbol}   {name}", lambda checked=False, i=index: self.navigate(i))
            item.setObjectName("nav")
            item.setCheckable(True)
            item.setAccessibleName(name)
            layout.addWidget(item)
            self.nav_buttons.append(item)
        layout.addStretch()
        tip, content = card()
        self.sidebar_tip = tip
        content.setContentsMargins(13, 15, 13, 15)
        content.addWidget(label("First time here?", "cardTitle"))
        content.addWidget(label("Learn with a synthetic dataset.", "small", True))
        self.demo_button = button("Practice project", self.load_demo)
        content.addWidget(self.demo_button)
        tip.setMinimumHeight(148)
        layout.addWidget(tip)
        layout.addSpacing(12)
        layout.addWidget(label("Built for discovery.\nDesigned for you.", "small"))
        outer.addWidget(sidebar)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "sidebar_tip"):
            self.sidebar_tip.setVisible(self.height() >= 810)

    def page(self):
        widget = QWidget()
        widget.setObjectName("page")
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(18)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(widget)
        self.pages.addWidget(scroll)
        return widget, layout

    def heading(self, layout, title, subtitle):
        layout.addWidget(label(title, "title"))
        layout.addWidget(label(subtitle, "muted", True))

    def build_overview(self):
        _, layout = self.page()
        hero = QFrame()
        hero.setObjectName("hero")
        content = QHBoxLayout(hero)
        content.setContentsMargins(28, 23, 18, 20)
        words = QVBoxLayout()
        words.addWidget(label("LESS FRICTION. MORE DISCOVERY.", "eyebrow"))
        words.addSpacing(8)
        words.addWidget(label("A clearer view of your isolates.", "heroTitle", True))
        words.addWidget(label("From sequences to meaningful comparisons.\nOne calm workspace for your microbial genomics.", "muted", True))
        words.addSpacing(10)
        row = QHBoxLayout()
        row.addWidget(button("＋  Import sequences", self.browse_files, True))
        row.addWidget(button("Explore your samples  →", lambda: self.navigate(1)))
        row.addStretch()
        words.addLayout(row)
        content.addLayout(words, 3)
        self.helix = Helix()
        content.addWidget(self.helix, 1)
        layout.addWidget(hero)
        metrics = QHBoxLayout()
        metrics.setSpacing(14)
        self.metrics = [Metric("Sequence files", "in this project"), Metric("Analysed", "results saved locally"), Metric("Exact profiles", "all loci assigned"), Metric("Review needed", "missing, mixed, or failed")]
        for item in self.metrics:
            metrics.addWidget(item)
        layout.addLayout(metrics)
        self.drop_zone = DropZone()
        self.drop_zone.filesDropped.connect(self.import_paths)
        self.drop_zone.browseRequested.connect(self.browse_files)
        layout.addWidget(self.drop_zone)
        recent, content = card()
        title = QHBoxLayout()
        title.addWidget(label("Your latest samples", "cardTitle"))
        title.addStretch()
        title.addWidget(button("View all  →", lambda: self.navigate(1)))
        content.addLayout(title)
        self.recent_table = table(["SAMPLE", "INPUT", "STATUS", "SEQUENCE TYPE", "LOCI CALLED"])
        self.recent_table.cellDoubleClicked.connect(lambda *_: self.navigate(1))
        content.addWidget(self.recent_table)
        layout.addWidget(recent, 1)
        self.practice_notice = label("Practice project · synthetic sequences for learning; not biological isolates.", "small")
        self.practice_notice.hide()
        layout.addWidget(self.practice_notice)

    def build_samples(self):
        _, layout = self.page()
        self.heading(layout, "Your samples", "Import sequences, choose a scheme, and start. Select a row to inspect its evidence.")
        controls = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search sample names…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_tables)
        controls.addWidget(self.search, 1)
        controls.addWidget(button("Add files", self.browse_files))
        controls.addWidget(button("Add folder", self.browse_folder))
        controls.addWidget(button("Edit metadata", self.edit_metadata))
        controls.addWidget(button("Remove row", self.remove_sample))
        layout.addLayout(controls)
        run, line = card()
        strip = QHBoxLayout()
        strip.addWidget(label("Typing scheme", "cardTitle"))
        self.scheme_combo = QComboBox()
        self.scheme_combo.setMinimumWidth(300)
        self.scheme_combo.setAccessibleName("Typing scheme")
        strip.addWidget(self.scheme_combo, 1)
        self.run_button = button("▶  Analyse pending", self.start_analysis, True)
        strip.addWidget(self.run_button)
        self.rerun_button = button("Reanalyse all", lambda: self.start_analysis(all_samples=True))
        strip.addWidget(self.rerun_button)
        line.addLayout(strip)
        line.addWidget(label("Exact known-allele matching on assemblies. FASTQ inputs receive read QC only.", "small", True))
        layout.addWidget(run)
        splitter = QSplitter(Qt.Orientation.Vertical)
        sample_card, content = card()
        content.setContentsMargins(0, 6, 0, 0)
        self.sample_table = table(["SAMPLE", "INPUT", "STATUS", "SEQUENCE TYPE", "LOCI CALLED"])
        self.sample_table.itemSelectionChanged.connect(self.show_sample_detail)
        content.addWidget(self.sample_table)
        splitter.addWidget(sample_card)
        self.detail = QTextBrowser()
        self.detail.setOpenExternalLinks(False)
        self.detail.setMinimumHeight(155)
        splitter.addWidget(self.detail)
        splitter.setSizes([350, 220])
        layout.addWidget(splitter, 1)

    def build_compare(self):
        _, layout = self.page()
        self.heading(layout, "See how your isolates connect", "Allele differences across comparable profiles. A tree describes genetic similarity; it does not establish transmission.")
        controls = QHBoxLayout()
        controls.addWidget(label("Minimum shared loci", "muted"))
        self.overlap = QDoubleSpinBox()
        self.overlap.setRange(0.01, 1.0)
        self.overlap.setDecimals(2)
        self.overlap.setSingleStep(0.05)
        self.overlap.setValue(0.95)
        self.overlap.setToolTip("Shared exact loci divided by all loci. 0.95 means at least 95%.")
        controls.addWidget(self.overlap)
        controls.addWidget(button("Build comparison", self.refresh_comparison, True))
        controls.addStretch()
        controls.addWidget(button("Fit view", lambda: self.tree.fit_tree()))
        controls.addWidget(button("Save PNG", self.save_tree))
        controls.addWidget(button("Export distances", self.export_distances))
        layout.addLayout(controls)
        grouping = QHBoxLayout()
        grouping.addWidget(label("Colour groups linked by at most", "small"))
        self.cluster_threshold = QSpinBox()
        self.cluster_threshold.setRange(0, 100000)
        self.cluster_threshold.setValue(1)
        self.cluster_threshold.setSuffix(" differences")
        self.cluster_threshold.setToolTip("Visual groups use single-linkage at this threshold. They are not outbreak assignments.")
        self.cluster_threshold.valueChanged.connect(self.refresh_comparison)
        grouping.addWidget(self.cluster_threshold)
        grouping.addWidget(label("Single-linkage groups · line lengths are for layout only", "small"))
        grouping.addStretch()
        layout.addLayout(grouping)
        self.tree = TreeView()
        layout.addWidget(self.tree, 1)
        self.tree_status = label("Identical scheme fingerprints are required. Missing calls never become zero-distance evidence.", "small", True)
        layout.addWidget(self.tree_status)
        self.distance_rows = []

    def build_schemes(self):
        _, layout = self.page()
        self.heading(layout, "Your scheme library", "Local, versioned allele collections for MLST and cgMLST. Each analysis records the database fingerprint.")
        row = QHBoxLayout()
        row.addWidget(button("Import a scheme folder", self.import_scheme, True))
        row.addStretch()
        row.addWidget(label("No automatic database updates", "badge"))
        layout.addLayout(row)
        self.scheme_table = table(["SCHEME", "SOURCE", "LOCATION"])
        layout.addWidget(self.scheme_table, 1)
        note, content = card()
        content.addWidget(label("What belongs in a scheme folder?", "cardTitle"))
        content.addWidget(label("One allele FASTA file per locus (.tfa, .fasta, .fa or .fna), with headers such as adk_1. Add a TSV profile table with an ST column to assign sequence types. cgMLST collections without an ST table can still produce allele profiles.", "muted", True))
        content.addWidget(label("Importing copies the scheme into your local library. Keep its source citation and confirm that its data terms allow your intended use.", "small", True))
        layout.addWidget(note)

    def build_hydra(self):
        _, layout = self.page()
        self.heading(layout, "HYDRA insights", "Bring resistance, virulence and lineage evidence alongside your typing work.")
        row = QHBoxLayout()
        row.addWidget(button("Import HYDRA JSON", self.import_hydra, True))
        row.addWidget(button("Export imported evidence", self.export_hydra))
        row.addStretch()
        row.addWidget(label("Result integration", "badge"))
        layout.addLayout(row)
        self.hydra_summary = label("No HYDRA report imported in this project.", "small", True)
        layout.addWidget(self.hydra_summary)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.hydra_table = table(["SAMPLE", "SPECIES", "SEQUENCE TYPE", "AMR GENES", "VIRULENCE", "PLASMIDS"])
        self.hydra_table.itemSelectionChanged.connect(self.show_hydra_sample)
        splitter.addWidget(self.hydra_table)
        self.hydra_view = QTextBrowser()
        self.hydra_view.setHtml("<h2>More context for each genome</h2><p>Import a HYDRA JSON report to inspect its elements, sample summaries and provenance in this native application.</p><p>HYDRA analysis execution is not bundled in this preview. Its Windows dependencies need a separate port and validation.</p>")
        splitter.addWidget(self.hydra_view)
        splitter.setSizes([230, 320])
        layout.addWidget(splitter, 1)

    def build_reports(self):
        _, layout = self.page()
        self.heading(layout, "Share the evidence", "Export every project sample, including failures and pending work. Search filters do not change these exports.")
        for title, description, formats in [
            ("A report you can read", "A self-contained report with sample-level evidence, quality measures and method notes.", [("Save PDF report", "pdf"), ("Save HTML report", "html")]),
            ("Ready for your next analysis", "Open tables in Excel, or retain the full structured results and provenance as JSON.", [("Save CSV", "csv"), ("Save TSV", "tsv"), ("Save JSON", "json")]),
        ]:
            frame, content = card()
            content.addWidget(label(title, "cardTitle"))
            content.addWidget(label(description, "muted", True))
            row = QHBoxLayout()
            for text, fmt in formats:
                row.addWidget(button(text, lambda checked=False, f=fmt: self.export_project(f)))
            row.addStretch()
            content.addLayout(row)
            layout.addWidget(frame)
        frame, content = card()
        content.addWidget(label("Take your workspace with you", "cardTitle"))
        content.addWidget(label("Save a project copy with all results and metadata. Sequence inputs remain at their original paths; copy them separately if you need to run the analyses on another computer.", "muted", True))
        content.addWidget(button("Save project copy…", self.save_project_copy))
        layout.addWidget(frame)
        layout.addStretch()

    def build_settings(self):
        _, layout = self.page()
        self.heading(layout, "Make yourself at home", "A few settings, and a straightforward guide to what this version can do.")
        frame, content = card()
        self.motion = QCheckBox("Gentle interface animations")
        self.motion.setChecked(bool(self.project.get_setting("motion", True)))
        self.motion.toggled.connect(self.set_motion)
        content.addWidget(self.motion)
        content.addWidget(label("Data location: " + str(self.root), "small", True))
        content.addWidget(button("Create a new project…", self.new_project))
        layout.addWidget(frame)
        guide = QTextBrowser()
        guide.setHtml("""<h2>From files to a comparison</h2>
        <p><b>1. Import and assign</b> FASTA or FASTQ files. Choose automatic MLST evidence, a manual organism/scheme, or unknown for QC/AMR-only work. Managed storage and ST filename suffixes affect copies only.</p>
        <p><b>2. Review the run plan</b> before analysis. Paired short reads can be validated and assembled with the installed native SKESA runtime. Unassembled reads receive clearly labelled sampled QC.</p>
        <p><b>3. Analyse</b> using the per-sample workflow. Automatic MLST provides provisional organism evidence, not independent species confirmation. Optional HYDRA runs against an installed reference snapshot and links evidence to stable sample IDs.</p>
        <p><b>4. Compare</b> an explicit cohort with one MLST/cgMLST/wgMLST snapshot. Additional profiles do not overwrite classical MLST. Change node colors, metadata groups, labels and layout without changing genetic distances.</p>
        <p><b>5. Report and reuse</b> selected isolates, user-highlighted groups and linked features. Saved libraries and portable bundles let you reuse profiles without reanalysing genomes. History retains prior evidence and input provenance.</p>
        <h3>Understanding results</h3><p><b>Exact match:</b> every locus has one known allele and the profile has a registered ST. <b>New combination:</b> the known alleles form an unregistered profile; this is not a new allele. <b>Missing loci:</b> exact sequence evidence is absent. <b>Mixed alleles:</b> more than one allele was found.</p>
        <h3>cgMLST and reference control</h3><p>Large schemes support exact-first, complete-CDS guarded novel calling. A local SHA-256 sequence identifier is not a registered allele number. Missing, duplicated, frame-disrupted or ambiguous loci require review. Local ad-hoc schemes are cohort-defined and cannot be treated as validated public nomenclature.</p>
        <p>Use Data to inspect/install versioned schemes and AMR reference snapshots. Downloads never send sequences. Provider rights and authentication can limit availability; a public snapshot is not necessarily the complete current database.</p>
        <h3>Current boundaries</h3><p>This is research software, not a validated diagnostic device or established SeqSphere+ equivalent. It does not predict measured susceptibility or prove transmission. Independent species confirmation, contamination quantification, direct-read AMR/pileup, long-read assembly and the full Kleborate/AMRFinderPlus/staphylococcal module stack are not validated here. Review biological quality, thresholds and database versions before interpreting a cluster.</p>
        <p>Keyboard: Ctrl+O import · Ctrl+Shift+O open project · Ctrl+S save project copy · Ctrl+R analyse selected · Ctrl+K commands · Alt+1…7 workspace pages.</p>""")
        layout.addWidget(guide, 1)

    def add_shortcuts(self):
        for shortcut, callback in [("Ctrl+O", self.browse_files), ("Ctrl+Shift+O", self.open_project_dialog), ("Ctrl+S", self.save_project_copy)]:
            action = QAction(self)
            action.setShortcut(QKeySequence(shortcut))
            action.triggered.connect(callback)
            self.addAction(action)

    def navigate(self, index):
        self.pages.setCurrentIndex(index)
        self.breadcrumb.setText("WORKSPACE  /  " + self.nav_names[index].upper())
        for i, item in enumerate(self.nav_buttons):
            item.setChecked(i == index)
        page = self.pages.currentWidget()
        if page:
            page.update()
            if isinstance(page, QScrollArea):
                page.viewport().update()
                page.widget().update()
        if index == 2:
            self.refresh_comparison()

    def set_motion(self, enabled):
        self.motion_enabled = enabled
        self.project.set_setting("motion", enabled)
        self.helix.timer.start(60) if enabled else self.helix.timer.stop()

    def notify(self, message):
        self.progress_text.setText(message)
        self.statusBar().showMessage(message, 15000)

    def error(self, message):
        QMessageBox.warning(self, "Let's resolve this", str(message))

    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Import sequence files", "", FILE_FILTER)
        if paths:
            self.import_paths(paths)

    def browse_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Import a sequence folder")
        if path:
            self.import_paths([path])

    def import_paths(self, paths):
        if self.worker and self.worker.isRunning():
            self.notify("Please finish or cancel the analysis before importing more files.")
            return
        supported = (".fa", ".fasta", ".fna", ".fq", ".fastq")
        def acceptable(path):
            name = path.name.lower()
            for compression in (".gz", ".bz2"):
                if name.endswith(compression):
                    name = name[:-len(compression)]
            return name.endswith(supported)
        count, skipped, errors = 0, 0, []
        existing = {str(Path(s["input_path"]).resolve()) for s in self.project.samples()}
        for source in paths:
            path = Path(source)
            candidates = path.rglob("*") if path.is_dir() else [path]
            for candidate in candidates:
                if not candidate.is_file() or not acceptable(candidate):
                    skipped += 1
                    continue
                try:
                    resolved = str(candidate.resolve())
                    if resolved in existing:
                        skipped += 1
                        continue
                    from wmlstudio.sequence import sample_name
                    self.project.add_sample(candidate, name=sample_name(candidate))
                    existing.add(resolved)
                    count += 1
                except Exception as exc:
                    errors.append(f"{candidate.name}: {exc}")
        self.refresh()
        self.navigate(1)
        self.notify(f"Added {count} sequence files. {skipped} unsupported entries skipped." + (f" {len(errors)} could not be imported." if errors else ""))
        if errors:
            self.error("\n".join(errors[:8]))

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        self.import_paths([url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()])
        event.acceptProposedAction()

    def populate_schemes(self):
        selected = self.project.get_setting("scheme_path", "")
        self.scheme_paths = scheme_locations(self.root)
        self.scheme_combo.clear()
        self.scheme_combo.addItem("Quality checks only", None)
        self.scheme_table.setRowCount(len(self.scheme_paths))
        for index, path in enumerate(self.scheme_paths):
            name = "Practice scheme · 7 loci (synthetic)" if path.name == "practice_7" else path.name.replace("_", " ")
            self.scheme_combo.addItem(name, str(path))
            if str(path) == selected:
                self.scheme_combo.setCurrentIndex(index + 1)
            for col, text in enumerate([name, "Local import" if self.root in path.parents else "Bundled snapshot", str(path)]):
                item = QTableWidgetItem(text)
                item.setToolTip(text)
                self.scheme_table.setItem(index, col, item)

    def import_scheme(self):
        if self.worker and self.worker.isRunning():
            self.notify("Please finish or cancel the current analysis first.")
            return
        source = QFileDialog.getExistingDirectory(self, "Select a folder containing one allele FASTA per locus")
        if not source:
            return
        self.worker = SchemeImportWorker(source, self.root, self)
        self.worker.imported.connect(self.scheme_imported)
        self.worker.failed.connect(self.error)
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.run_button.setEnabled(False)
        self.rerun_button.setEnabled(False)
        self.scheme_combo.setEnabled(False)
        self.demo_button.setEnabled(False)
        self.cancel_button.show()
        self.progress_bar.show()
        self.worker.start()

    def scheme_imported(self, destination, loci):
        self.project.set_setting("scheme_path", destination)
        self.populate_schemes()
        self.notify(f"Imported {loci} loci. Select Samples to start typing.")

    def refresh(self):
        self.current_samples = self.project.samples()
        results = [s["result"] for s in self.current_samples if s["status"] == "completed" and s.get("result")]
        values = [len(self.current_samples), len(results), sum(r.get("status") == "complete" for r in results), sum(s["status"] in {"failed", "interrupted"} or (s.get("result") or {}).get("status") in {"mixed", "incomplete", "ambiguous"} for s in self.current_samples)]
        for metric, value in zip(self.metrics, values, strict=True):
            metric.value.setText(str(value))
        self.practice_notice.setVisible(bool(self.project.get_setting("practice", False)))
        self.refresh_tables()
        self.restore_hydra()

    def refresh_tables(self):
        text = self.search.text().casefold() if hasattr(self, "search") else ""
        samples = [s for s in self.current_samples if text in s["name"].casefold()]
        self.fill_sample_table(self.sample_table, samples)
        self.fill_sample_table(self.recent_table, self.current_samples[-5:])
        self.show_sample_detail()

    def fill_sample_table(self, widget, samples):
        selected = widget.currentItem()
        selected_id = selected.data(Qt.ItemDataRole.UserRole) if selected else None
        widget.blockSignals(True)
        widget.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            result = sample.get("result") or {}
            status = result.get("status", sample["status"]) if sample["status"] == "completed" else sample["status"]
            alleles = result.get("alleles", {})
            called = sum(v is not None for v in alleles.values())
            name = sample["name"] + (" · file missing" if sample.get("missing_input") else "")
            kind = result.get("kind") or ("fastq" if any(s in Path(sample["input_path"]).suffixes for s in (".fq", ".fastq")) else "fasta")
            values = [name, "Reads" if kind == "fastq" else "Assembly", STATUS_TEXT.get(status, status), f"ST {result['st']}" if result.get("st") else "—", f"{called} / {len(alleles)}" if alleles else "—"]
            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, sample["id"])
                item.setToolTip(sample.get("error") or sample["input_path"])
                if col == 2:
                    item.setForeground(QColor("#197866" if status in {"complete", "qc_only"} else "#AF7737" if status in {"incomplete", "mixed", "failed", "interrupted"} else "#617A71"))
                widget.setItem(row, col, item)
            if sample["id"] == selected_id:
                widget.selectRow(row)
        widget.blockSignals(False)

    def selected_sample(self):
        item = self.sample_table.currentItem()
        if item:
            return next((s for s in self.current_samples if s["id"] == item.data(Qt.ItemDataRole.UserRole)), None)
        return None

    def show_sample_detail(self):
        sample = self.selected_sample()
        if not sample:
            self.detail.setHtml("<h3>Your evidence, one sample at a time</h3><p>Select a row above to see sequence quality, allele calls and the source of each result.</p>")
            return
        def e(value):
            return html.escape(str(value))
        result = sample.get("result") or {}
        qc = result.get("qc", {})
        body = f"<h2>{e(sample['name'])}</h2><p>{e(sample['input_path'])}</p>"
        if sample.get("error"):
            body += f"<p><b>Needs attention:</b> {e(sample['error'])}</p>"
        if sample.get("metadata"):
            body += "<p>" + " · ".join(f"{e(k)}: {e(v)}" for k, v in sample["metadata"].items()) + "</p>"
        if qc:
            shown = [("Records checked", qc.get("records")), ("Bases checked", qc.get("total_bases")), ("GC %", qc.get("gc_percent")), ("N50", qc.get("n50")), ("Q30 %", qc.get("q30_percent"))]
            body += "<p>" + " &nbsp; | &nbsp; ".join(f"<b>{key}:</b> {e(round(value, 2) if isinstance(value, float) else value)}" for key, value in shown if value is not None) + "</p>"
            if qc.get("sampled"):
                body += "<p><b>Sampled QC:</b> statistics describe the first records only; the rest of the read file has not been quality-validated.</p>"
        if result.get("calls"):
            body += "<table width='100%' cellspacing='0' cellpadding='7'><tr bgcolor='#EAF3ED'><th align='left'>Locus</th><th align='left'>Allele</th><th align='left'>Evidence</th></tr>"
            for call in result["calls"]:
                body += f"<tr><td>{e(call['locus'])}</td><td>{e(call.get('allele') or '—')}</td><td>{e(call['status'])}</td></tr>"
            body += "</table>"
        for note in result.get("notes", []):
            body += f"<p>{e(note)}</p>"
        if result.get("input_sha256"):
            body += f"<p><b>Input SHA-256:</b> {e(result['input_sha256'])}<br><b>Scheme fingerprint:</b> {e(result.get('scheme_digest') or 'Not used')}</p>"
        self.detail.setHtml(body)

    def start_analysis(self, checked=False, all_samples=False):
        if self.worker and self.worker.isRunning():
            return
        samples = [s for s in self.project.samples() if all_samples or s["status"] in {"queued", "failed", "interrupted"}]
        if not samples:
            self.notify("No pending files. Import sequences or choose Reanalyse all.")
            return
        scheme = self.scheme_combo.currentData()
        self.project.set_setting("scheme_path", scheme or "")
        self.worker = AnalysisWorker(samples, scheme, parent=self)
        self.worker.sample_started.connect(self.sample_started)
        self.worker.sample_finished.connect(self.sample_finished)
        self.worker.sample_failed.connect(self.sample_failed)
        self.worker.sample_cancelled.connect(lambda sid: self.project.set_status(sid, "interrupted", "Cancelled by user"))
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.run_button.setEnabled(False)
        self.rerun_button.setEnabled(False)
        self.scheme_combo.setEnabled(False)
        self.demo_button.setEnabled(False)
        self.cancel_button.show()
        self.progress_bar.show()
        self.progress_bar.setValue(0)
        self.worker.start()

    def sample_started(self, sample_id):
        self.project.set_status(sample_id, "running")
        self.refresh()

    def sample_finished(self, sample_id, result):
        self.project.set_result(sample_id, result)
        self.refresh()

    def sample_failed(self, sample_id, message):
        self.project.set_status(sample_id, "failed", message)
        self.refresh()

    def job_progress(self, percent, message):
        self.progress_bar.setValue(percent)
        self.notify(message)

    def cancel_analysis(self):
        if self.worker:
            self.worker.cancel()
            self.notify("Cancelling safely; completed results are already saved.")

    def analysis_finished(self):
        self.run_button.setEnabled(True)
        self.rerun_button.setEnabled(True)
        self.scheme_combo.setEnabled(True)
        self.demo_button.setEnabled(True)
        self.cancel_button.hide()
        self.progress_bar.hide()
        self.refresh()
        if self.closing_after_cancel:
            self.close()

    def comparison_results(self):
        return [dict(s["result"], sample_id=s["id"], sample_name=s["name"]) for s in self.project.samples() if s["status"] == "completed" and s.get("result", {}).get("alleles")]

    def refresh_comparison(self):
        results = self.comparison_results()
        if len(results) > 250:
            self.tree_status.setText("This preview displays up to 250 profiles per comparison. Use a smaller project for an interactive view.")
            self.tree.draw_results([], [])
            self.distance_rows = []
            return
        try:
            self.distance_rows = pairwise_distances(results, self.overlap.value())
            edges = minimum_spanning_forest(results, self.overlap.value())
            self.tree.draw_results(results, edges, self.cluster_threshold.value())
            excluded = sum(not row["comparable"] for row in self.distance_rows)
            self.tree_status.setText(f"{len(results)} profiles · {len(edges)} edges · {excluded} pairs excluded for insufficient or incompatible evidence. Edge labels are allele differences. Drag nodes, pan, or scroll to zoom.")
        except Exception as exc:
            self.tree_status.setText(str(exc))

    def save_tree(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save comparison image", "comparison.png", "PNG image (*.png)")
        if path:
            try:
                self.check_output(path)
                self.tree.save_png(path)
                self.notify("Comparison image saved.")
            except Exception as exc:
                self.error(exc)

    def export_distances(self):
        self.refresh_comparison()
        path, _ = QFileDialog.getSaveFileName(self, "Save distance evidence", "distances.json", "JSON (*.json)")
        if path:
            try:
                self.check_output(path)
                from wmlstudio.export import write_distances
                write_distances(self.distance_rows, path, self.overlap.value())
                self.notify("Distance evidence saved, including shared-locus denominators and exclusion reasons.")
            except Exception as exc:
                self.error(exc)

    def export_project(self, fmt):
        path, _ = QFileDialog.getSaveFileName(self, "Export project", self.project_path.stem + "." + fmt, f"{fmt.upper()} (*.{fmt})")
        if path:
            try:
                self.check_output(path)
                if fmt == "pdf":
                    self.write_pdf_report(path)
                else:
                    export_results(self.project.samples(), path, fmt)
                self.notify(f"Exported all {len(self.project.samples())} project samples.")
            except Exception as exc:
                self.error(exc)

    def write_pdf_report(self, path, samples=None):
        def escaped(value):
            return html.escape(str(value))

        samples = self.project.samples() if samples is None else samples
        sections = ["<html><body style='font-family: sans-serif; color: #203736'>",
                    "<h1 style='color: #187D6D'>WMLSTudio</h1><h2>Sequence typing report</h2>",
                    f"<p>{len(samples)} sample records · Application {__version__}</p>",
                    "<p>Results describe the supplied sequence data and scheme snapshot. Allele similarity alone does not establish transmission. Missing, mixed and ambiguous calls require review.</p>",
                    "<table width='100%' cellpadding='6' cellspacing='0' border='1'><tr><th>Sample</th><th>Analysis</th><th>Result</th><th>ST</th><th>Called loci</th></tr>"]
        for sample in samples:
            result = sample.get("result") or {}
            alleles = result.get("alleles", {})
            values = [sample["name"], STATUS_TEXT.get(sample["status"], sample["status"]), STATUS_TEXT.get(result.get("status"), result.get("status", "Not analysed")), result.get("st") or "Unassigned", f"{sum(v is not None for v in alleles.values())}/{len(alleles)}" if alleles else "Not typed"]
            highlight = sample.get("metadata", {}).get("cluster", {}).get("highlight")
            sections.append(("<tr bgcolor='#FFF0CD'>" if highlight else "<tr>") + "".join(f"<td>{escaped(v)}</td>" for v in values) + "</tr>")
        sections.append("</table>")
        for sample in samples:
            result = sample.get("result") or {}
            sections.append(f"<h2 style='page-break-before: always; color: #187D6D'>{escaped(sample['name'])}</h2>")
            group = sample.get("metadata", {}).get("cluster", {})
            if group.get("highlight"):
                sections.append(f"<p style='background-color: #FFF0CD'><b>Highlighted group: {escaped(group.get('label') or 'Selected group')}</b></p>")
            sections.append(f"<p><b>Input:</b> {escaped(sample['input_path'])}<br><b>Scheme:</b> {escaped(result.get('scheme') or 'Not used')}<br><b>Sequence type:</b> {escaped(result.get('st') or 'Unassigned')}<br><b>Job state:</b> {escaped(sample['status'])}</p>")
            if sample.get("error"):
                sections.append(f"<p><b>Needs attention:</b> {escaped(sample['error'])}</p>")
            from wmlstudio.sample_workflow import current_hydra_evidence, hydra_evidence_status
            from wmlstudio.ui_common import flattened_metadata, gene_names, organism_for
            genus, species, organism_evidence = organism_for(sample)
            amr_state = hydra_evidence_status(sample)
            gene_summary = "Archived evidence — no current AMR result" if amr_state["status"] == "stale" else "No linked report" if amr_state["status"] == "missing" else "; ".join(gene_names(sample)) or "No primary AMR genes reported"
            sections.append(f"<p><b>Organism:</b> {escaped(' '.join([genus, species]).strip() or 'Unknown')} ({escaped(organism_evidence)})<br><b>AMR genes:</b> {escaped(gene_summary)}</p>")
            metadata = flattened_metadata(sample)
            if metadata:
                title = "Sample metadata (includes archived AMR fields, not current calls)" if amr_state["status"] == "stale" else "Sample annotations and workflow"
                sections.append(f"<h3>{title}</h3><p>" + "<br>".join(f"<b>{escaped(k)}:</b> {escaped(v)}" for k, v in metadata.items()) + "</p>")
            if amr_state["status"] in {"stale", "unverified"}:
                sections.append(f"<p><b>AMR evidence: {escaped(amr_state['status'])}.</b> {escaped(amr_state['reason'])}</p>")
            hydra = current_hydra_evidence(sample)
            if hydra:
                sections.append("<h3>Linked HYDRA evidence</h3><p>Sequence detections are not measured susceptibility. Secondary hits are retained in the exported JSON evidence, not counted twice as primary genes.</p><table width='100%' cellpadding='5' cellspacing='0' border='1'><tr><th>Gene</th><th>Element</th><th>Method</th><th>Identity / coverage</th></tr>")
                for hit in hydra.get("hits", []):
                    if hit.get("primary") is True:
                        values = [hit.get("gene"), hit.get("element_type"), hit.get("method"), f"{hit.get('identity_pct', '—')} / {hit.get('coverage_pct', '—')}"]
                        sections.append("<tr>" + "".join(f"<td>{escaped(value)}</td>" for value in values) + "</tr>")
                sections.append("</table><p><b>HYDRA report SHA-256:</b><br>" + escaped(hydra.get("report_sha256", "Unknown")) + "</p>")
            for profile in sample.get("analyses", []):
                if profile.get("scheme_digest") != result.get("scheme_digest"):
                    alleles = profile.get("alleles", {})
                    sections.append(f"<h3>Additional profile · {escaped(profile.get('scheme'))}</h3><p>{sum(value is not None for value in alleles.values())}/{len(alleles)} loci called · {escaped(profile.get('status'))}<br>Scheme SHA-256: {escaped(profile.get('scheme_digest'))}</p><p>The complete additional locus matrix is retained in the JSON report and portable profile bundle.</p>")
            qc = result.get("qc", {})
            if qc:
                sections.append("<h3>Sequence quality</h3><table width='100%' cellpadding='5' cellspacing='0'>")
                metrics = [("records", "Records inspected"), ("total_bases", "Bases inspected"), ("n50", "N50 (bases)"), ("gc_percent", "GC (%)"), ("n_percent", "N bases (%)"), ("ambiguous_percent", "Ambiguous bases (%)"), ("mean_quality", "Mean Phred+33 quality"), ("q30_percent", "Q30 bases (%)")]
                for key, title in metrics:
                    value = qc.get(key)
                    if value is not None:
                        shown = f"{value:,.2f}" if isinstance(value, float) else f"{value:,}"
                        sections.append(f"<tr><td>{title}</td><td>{shown}</td></tr>")
                sections.append("</table>")
                if qc.get("sampled"):
                    sections.append("<p><b>Sampled read QC:</b> these values describe a prefix of the reads, not quality validation of the entire file.</p>")
            if result.get("notes"):
                sections.append("<h3>Method and review notes</h3><ul>" + "".join(f"<li>{escaped(note)}</li>" for note in result["notes"]) + "</ul>")
            if result.get("calls"):
                sections.append("<h3>Allele evidence</h3><table width='100%' cellpadding='5' cellspacing='0' border='1'><tr><th>Locus</th><th>Allele</th><th>Evidence</th></tr>")
                for call in result["calls"]:
                    sections.append(f"<tr><td>{escaped(call['locus'])}</td><td>{escaped(call.get('allele') or 'Unassigned')}</td><td>{escaped(call.get('reason') or call['status'])}</td></tr>")
                sections.append("</table>")
            sections.append("<h3>Provenance</h3>")
            for key, title in [("input_sha256", "Input SHA-256"), ("scheme_digest", "Scheme SHA-256")]:
                value = str(result.get(key) or "Not used")
                sections.append(f"<p><b>{title}</b><br><span style='font-size: 8pt'>{escaped(value[:32])}<br>{escaped(value[32:])}</span></p>")
        sections.append("</body></html>")
        with tempfile.TemporaryDirectory(prefix="wmlstudio-report-", dir=Path(path).parent) as folder:
            output = Path(folder) / "report.pdf"
            document = QTextDocument()
            document.setHtml("".join(sections))
            writer = QPdfWriter(str(output))
            writer.setTitle("WMLSTudio · Sequence typing report")
            writer.setCreator(f"WMLSTudio {__version__}")
            writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            writer.setResolution(144)
            document.print_(writer)
            del writer
            if not output.is_file() or output.stat().st_size < 100:
                raise OSError("PDF generation failed. Please choose another output location.")
            output.replace(path)

    def edit_metadata(self):
        sample = self.selected_sample()
        if not sample:
            self.notify("Select a sample first.")
            return
        text, accepted = QInputDialog.getMultiLineText(self, "Sample metadata", "One field per line, for example: location = Laboratory A", "\n".join(f"{key} = {value}" for key, value in sample.get("metadata", {}).items()))
        if accepted:
            try:
                metadata = {}
                for line in text.splitlines():
                    if not line.strip():
                        continue
                    key, sep, value = line.partition("=")
                    if not sep or not key.strip():
                        raise ValueError("Use field = value on each line.")
                    metadata[key.strip()] = value.strip()
                self.project.set_metadata(sample["id"], metadata)
                self.refresh()
            except Exception as exc:
                self.error(exc)

    def remove_sample(self):
        if self.worker and self.worker.isRunning():
            self.notify("Please finish or cancel the current analysis first.")
            return
        sample = self.selected_sample()
        if sample and QMessageBox.question(self, "Remove from project?", f"Remove {sample['name']} and its saved analysis from this project? The original sequence file stays untouched.") == QMessageBox.StandardButton.Yes:
            self.project.remove_sample(sample["id"])
            self.refresh()

    def switch_project(self, path):
        if self.worker and self.worker.isRunning():
            self.notify("Please finish or cancel the current analysis before changing projects.")
            return False
        if Path(path).resolve() == self.project_path.resolve():
            return True
        new_lock = None
        try:
            new_lock = self.lock_project(path)
            new = Project(path)
            self.project.close()
            self.project_lock.unlock()
            self.project_lock = new_lock
            self.project = new
            self.project_path = Path(path)
            self.project_label.setText(self.project_path.stem)
            enabled = bool(self.project.get_setting("motion", True))
            self.motion.blockSignals(True)
            self.motion.setChecked(enabled)
            self.motion.blockSignals(False)
            self.set_motion(enabled)
            self.populate_schemes()
            self.refresh()
            self.navigate(0)
            return True
        except Exception as exc:
            if new_lock is not None and new_lock is not self.project_lock:
                new_lock.unlock()
            self.error(exc)
            return False

    @staticmethod
    def lock_project(path):
        lock = QLockFile(str(Path(path).resolve()) + ".lock")
        if not lock.tryLock(0):
            raise ValueError("This project is already open in another WMLSTudio window, or its folder is not writable. Close the other window or choose a writable location.")
        return lock

    def open_project_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open WMLSTudio project", str(self.root), "WMLSTudio project (*.wmlstudio)")
        if path:
            self.switch_project(path)

    def new_project(self):
        path, _ = QFileDialog.getSaveFileName(self, "Create project", str(self.root / "New project.wmlstudio"), "WMLSTudio project (*.wmlstudio)")
        if path:
            if Path(path).exists():
                self.error("A project already exists there. Choose a new filename, or open the existing project.")
            else:
                self.switch_project(path)

    def save_project_copy(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save project copy", str(self.root / (self.project_path.stem + " copy.wmlstudio")), "WMLSTudio project (*.wmlstudio)")
        if not path:
            return
        if Path(path).resolve() == self.project_path.resolve():
            self.notify("Your project is already saved automatically.")
            return
        try:
            self.check_output(path)
            copy_lock = self.lock_project(path)
            source = sqlite3.connect(self.project_path)
            temporary = tempfile.NamedTemporaryFile(prefix=".wmlstudio-copy-", suffix=".tmp", dir=Path(path).parent, delete=False)
            temporary.close()
            temporary_path = Path(temporary.name)
            target = sqlite3.connect(temporary_path)
            try:
                source.backup(target)
            finally:
                source.close()
                target.close()
            temporary_path.replace(path)
            self.notify("Project copy saved. Input sequences remain at their original paths.")
        except Exception as exc:
            self.error(exc)
        finally:
            if "temporary_path" in locals():
                temporary_path.unlink(missing_ok=True)
            if "copy_lock" in locals():
                copy_lock.unlock()

    def check_output(self, path):
        protected = [self.project_path, *(sample["input_path"] for sample in self.project.samples())]
        report = self.project.get_setting("hydra_report", None)
        if report:
            protected.append(report.get("import_provenance", {}).get("source_path", ""))
        ensure_separate_destination(path, protected)

    def load_demo(self):
        if self.worker and self.worker.isRunning():
            return
        if not self.switch_project(self.root / "Practice project.wmlstudio"):
            return
        scheme, samples = create_demo(self.root)
        self.project.set_setting("practice", True)
        self.project.set_setting("scheme_path", str(scheme))
        existing = {str(Path(s["input_path"]).resolve()) for s in self.project.samples()}
        for path in samples:
            if str(path.resolve()) not in existing:
                self.project.add_sample(path, name=path.stem)
        self.populate_schemes()
        self.refresh()
        self.start_analysis()
        self.navigate(0)

    def import_hydra(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import HYDRA evidence", "", "HYDRA JSON (*.json)")
        if not path:
            return
        try:
            from wmlstudio.hydra import load_hydra_report
            report = load_hydra_report(path)
            self.project.set_setting("hydra_report", report)
            self.restore_hydra()
            self.notify("HYDRA report imported. It retains its own sample identities and provenance.")
        except Exception as exc:
            self.error(exc)

    def restore_hydra(self):
        report = self.project.get_setting("hydra_report", None)
        if not report:
            self.hydra_table.setRowCount(0)
            self.hydra_summary.setText("No HYDRA report imported in this project.")
            self.hydra_view.setHtml("<h2>AMR evidence connected to each isolate</h2><p>Select assemblies in Samples, then use Run HYDRA to analyse them against an installed reference snapshot. You can also import an existing HYDRA JSON report and explicitly map its sample identities.</p><p>Linked genes and mutation evidence are available alongside typing, quality metrics, and annotations. Genotypic detections are not measured susceptibility.</p>")
            return
        samples = report["samples"]
        provenance = report.get("import_provenance", {})
        self.hydra_summary.setText(f"HYDRA {report.get('hydra_version', 'unknown')} · {len(samples)} samples · Imported from {Path(provenance.get('source_path', 'report')).name}. Imported evidence remains separate from locally computed calls.")
        selected = max(0, self.hydra_table.currentRow())
        self.hydra_table.blockSignals(True)
        self.hydra_table.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            summary = sample.get("summary", {})
            values = [sample["sample"], sample.get("species", {}).get("name") or "Unknown", sample.get("mlst", {}).get("sequence_type") or "—", summary.get("amr_genes", 0), summary.get("virulence_genes", 0), summary.get("plasmid_replicons", 0)]
            for col, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                item.setToolTip(str(value))
                self.hydra_table.setItem(row, col, item)
        self.hydra_table.blockSignals(False)
        if samples:
            self.hydra_table.selectRow(min(selected, len(samples) - 1))
            self.show_hydra_sample()

    def show_hydra_sample(self):
        report = self.project.get_setting("hydra_report", None)
        row = self.hydra_table.currentRow()
        if not report or row < 0 or row >= len(report["samples"]):
            return
        sample = report["samples"][row]
        def e(value):
            return html.escape(str(value))
        species = sample.get("species", {})
        body = f"<h2>{e(sample['sample'])}</h2><p><b>{e(species.get('name') or 'Unknown species')}</b> · Evidence: {e(species.get('evidence') or 'Not reported')}</p>"
        body += "<p>Imported findings describe sequence evidence. They are not a measured antimicrobial susceptibility result.</p>"
        warnings = list(report.get("import_warnings", [])) + list(sample.get("warnings", []))
        for warning in warnings:
            body += f"<p><b>Review:</b> {e(warning)}</p>"
        typing = sample.get("typing", [])
        if typing:
            body += "<h3>Lineage and typing</h3><p>" + e(json.dumps(typing, ensure_ascii=False)) + "</p>"
        hits = sample.get("hits", [])
        body += f"<h3>Detected elements · {len(hits)}</h3>"
        if hits:
            body += "<table width='100%' cellpadding='7' cellspacing='0'><tr bgcolor='#EAF3ED'><th align='left'>Gene / marker</th><th align='left'>Category</th><th align='left'>Database</th><th align='left'>Identity %</th><th align='left'>Coverage %</th><th align='left'>Evidence</th></tr>"
            for hit in hits[:1000]:
                evidence = str(hit.get("method") or "")
                if hit.get("allele_fraction") is not None:
                    evidence += f"; fraction {hit['allele_fraction']}"
                if hit.get("depth") is not None:
                    evidence += f"; depth {hit['depth']}"
                values = [hit.get("gene", "—"), hit.get("element_subtype") or hit.get("element_type", "—"), hit.get("database", "—"), hit.get("identity_pct", "—"), hit.get("coverage_pct", "—"), evidence]
                body += "<tr>" + "".join(f"<td>{e(value)}</td>" for value in values) + "</tr>"
            body += "</table>"
            if len(hits) > 1000:
                body += f"<p>Showing the first 1,000 of {len(hits)} elements. Export imported evidence to retain all elements.</p>"
        else:
            body += "<p>No elements are present in this sample's imported report.</p>"
        body += "<h3>Run provenance</h3>"
        body += f"<p><b>HYDRA version:</b> {e(report.get('hydra_version', 'unknown'))}<br><b>Command:</b> {e(report.get('command', 'Not recorded'))}<br><b>Imported report SHA-256:</b> {e(report.get('import_provenance', {}).get('sha256', 'Not recorded'))}</p>"
        body += "<p><b>Parameters:</b> " + e(json.dumps(report.get("parameters", {}), ensure_ascii=False)) + "</p>"
        body += "<p><b>Databases:</b> " + e(json.dumps(report.get("databases", []), ensure_ascii=False)) + "</p>"
        self.hydra_view.setHtml(body)

    def export_hydra(self):
        report = self.project.get_setting("hydra_report", None)
        if not report:
            self.notify("Import a HYDRA report first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export complete HYDRA evidence", "hydra-evidence.json", "JSON (*.json)")
        if path:
            try:
                self.check_output(path)
                from wmlstudio.export import _atomic_text
                with _atomic_text(path) as handle:
                    json.dump(report, handle, indent=2, ensure_ascii=False, allow_nan=False)
                self.notify("All imported HYDRA evidence and its provenance saved.")
            except Exception as exc:
                self.error(exc)

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.closing_after_cancel = True
            self.cancel_analysis()
            event.ignore()
            return
        self.project.close()
        self.project_lock.unlock()
        event.accept()


from wmlstudio.ui_compare import ComparisonWorkspaceMixin  # noqa: E402
from wmlstudio.ui_reports import ReportWorkspaceMixin  # noqa: E402
from wmlstudio.ui_workbench import WorkbenchMixin  # noqa: E402


class MainWindow(WorkbenchMixin, ComparisonWorkspaceMixin, ReportWorkspaceMixin, BaseWindow):
    """Native workbench composed from focused workflow controllers."""


def main(argv=None):
    parser = argparse.ArgumentParser(description="WMLSTudio native desktop")
    parser.add_argument("--project", type=Path)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--screenshot-page", type=int, choices=range(7), default=0,
                        help="Workspace page index for reproducible screenshots (0–6)")
    parser.add_argument("--window-size", default="1380x940", help="Desktop size WIDTHxHEIGHT")
    parser.add_argument("--native-screenshot", action="store_true", help="Capture the native window surface, not an offscreen render")
    args = parser.parse_args(argv)
    try:
        width, height = map(int, args.window_size.lower().split("x"))
        if not 1000 <= width <= 7680 or not 680 <= height <= 4320:
            raise ValueError
    except ValueError:
        parser.error("--window-size must be WIDTHxHEIGHT, at least 1000x680 and at most 7680x4320")
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("WMLSTudio")
    app.setOrganizationName("IOWA-BioTech")
    app.setStyle("Fusion")
    app.setStyleSheet(STYLE)
    temp = tempfile.TemporaryDirectory(prefix="wmlstudio-smoke-") if args.smoke_test else None
    try:
        window = MainWindow(args.project, Path(temp.name) if temp else None)
    except (OSError, ValueError, sqlite3.Error) as exc:
        QMessageBox.critical(None, "WMLSTudio could not open the workspace", str(exc))
        if temp:
            temp.cleanup()
        return 1
    window.resize(width, height)
    window.show()
    if args.demo:
        window.load_demo()
    if args.smoke_test or args.screenshot:
        capture_prepared = False
        def finish():
            nonlocal capture_prepared
            if window.worker_role or (window.worker and window.worker.isRunning()):
                QTimer.singleShot(100, finish)
                return
            if not capture_prepared:
                capture_prepared = True
                window.navigate(args.screenshot_page)
                window.set_motion(False)
                QTimer.singleShot(300, finish)
                return
            if window._comparison_timer.isActive() or (window.comparison_worker and window.comparison_worker.isRunning()):
                QTimer.singleShot(100, finish)
                return
            if args.screenshot:
                args.screenshot.parent.mkdir(parents=True, exist_ok=True)
                # Capture the settled page even when an emulated/native animation
                # clock delivers its first frame after the screenshot timer.
                window.helix.timer.stop()
                app.processEvents()
                if args.screenshot_page == 2:
                    window.tree.fit_tree()
                    app.processEvents()
                capture = window.screen().grabWindow(window.winId()) if args.native_screenshot else window.grab()
                if not capture.save(str(args.screenshot)):
                    print("WMLSTudio: could not save the requested screenshot", file=sys.stderr)
                    window.close()
                    app.exit(1)
                    return
            if args.smoke_test:
                window.close()
                app.quit()
        QTimer.singleShot(700, finish)
    code = app.exec()
    if temp:
        temp.cleanup()
    return code


if __name__ == "__main__":
    raise SystemExit(main())
