"""WMLSTudio native desktop entry point, project workflow and evidence views."""

import argparse
import html
import json
import sqlite3
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QEvent, QLockFile, Qt, QTimer
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
    QMenu,
    QMessageBox,
    QProgressBar,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import __version__, display
from wmlstudio.comparison import minimum_spanning_forest, pairwise_distances
from wmlstudio.context_menus import ContextMenuMixin
from wmlstudio.demo import create_demo
from wmlstudio.display import SCALE_CHOICES, apply_display_settings
from wmlstudio.export import ensure_separate_destination, export_results
from wmlstudio.jobs import AnalysisWorker, SchemeImportWorker
from wmlstudio.paths import data_root, scheme_locations
from wmlstudio.project import Project
from wmlstudio.theme import STYLE
from wmlstudio.ui_characterization import CharacterizationWorkspaceMixin
from wmlstudio.ui_common import workspace_header
from wmlstudio.ui_tabs import (
    PAGE_KEYS,
    PAGE_PURPOSE,
    PLANNED,
    STATION_KEYS,
    STATION_PAGES,
    WorkspaceTabs,
    clear_promise,
    page_name,
)
from wmlstudio.widgets import DropZone, Helix, Metric, TreeView, button, card, label
from wmlstudio.workspace_focus import CohortLedger, FocusBar, FocusBus

#: The gap between the workspace and each edge of the window. Read by the sidebar
#: decision as well, so the two can never disagree about how much room the tab bar
#: actually has.
WORKSPACE_GUTTER = 16

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
        # Focus is the shared, temporary highlight; the ledger records where each
        # tab's own cohort came from. A mixin may already have created them.
        self.focus = getattr(self, "focus", None) or FocusBus(self)
        self.cohort_origins = getattr(self, "cohort_origins", None) or CohortLedger(self)
        self._navigating = False
        self.page_shown = {}
        self.setWindowTitle("WMLSTudio · Microbial genomics workspace")
        self.resize(1380, 940)
        self.setMinimumSize(1000, 680)
        if QApplication.platformName() != "offscreen" and self.screen():
            available = self.screen().availableGeometry()
            self.resize(min(1380, available.width() - 40), min(940, available.height() - 60))
        # A window size the user chose in Settings outlives the session, but never
        # a size this screen cannot show: that is how a window becomes unreachable.
        chosen = display.read_window_size()
        if chosen:
            self.resize(*display.clamp_window_size(*chosen, self.available_screen_size()))
        self.setAcceptDrops(True)
        main = QWidget()
        main.setObjectName("main")
        self.setCentralWidget(main)
        outer = QHBoxLayout(main)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        self.build_sidebar(outer)
        body = QVBoxLayout()
        # The tab bar and the per-tab orientation strip both cost vertical space;
        # the margins pay for them so the comparison graph keeps its usable height.
        # The side gutter is WORKSPACE_GUTTER rather than the old 24: fourteen tabs
        # have to be drawn whole in a 1000 px window, and 16 px of empty margin is
        # a cheaper thing to give the bar than a tab behind a scroll arrow.
        body.setContentsMargins(WORKSPACE_GUTTER, 12, WORKSPACE_GUTTER, 10)
        body.setSpacing(12)
        outer.addLayout(body, 1)
        top = QHBoxLayout()
        self.breadcrumb = label("WORKSPACE  /  OVERVIEW", "eyebrow")
        top.addWidget(self.breadcrumb)
        top.addStretch()
        self.focus_bar = FocusBar(self.focus)
        self.focus_bar.showRequested.connect(self.show_focus_in_library)
        self.focus_bar.clearRequested.connect(self.focus.clear)
        # Kept under its original name: refresh_journey still writes the project
        # count here, and the bar shows it whenever nothing is focused.
        self.scope_label = self.focus_bar
        top.addWidget(self.focus_bar)
        top.addWidget(label("●  Local & private", "badge"))
        top.addWidget(button("Open project", self.open_project_dialog))
        top.addWidget(button("Settings", lambda: self.navigate("settings")))
        body.addLayout(top)
        self.page_index = {}
        self.stations = {}
        self.clear_buttons = {}
        self.pages = WorkspaceTabs()
        body.addWidget(self.pages, 1)
        # Built in the order these seven pages have always been built, because an
        # integer anywhere in this application means "the Nth page built". The tab
        # bar shows them in the pipeline order instead (ui_tabs.PIPELINE).
        self.build_overview()
        self.build_samples()
        self.build_compare()
        self.build_schemes()
        self.build_hydra()
        self.build_reports()
        self.build_settings()
        self.build_stations()
        # Appended after the stations, so the update page takes the next build
        # number rather than displacing one: every `navigate(<int>)` in this
        # application still means the page it has always meant.
        self.build_update()
        self.install_pipeline_pages()
        self.install_tree_station()
        self.install_typing_stations()
        self.install_workspace_headers()
        self.fold_duplicate_pages()
        # A page that reads the project every time it is opened says so here rather
        # than recomputing on every keystroke somewhere else in the window.
        self.page_shown = {key: hook for key, hook in (
            # The two tree tabs are the same workspace bound to two scales, so
            # each one mounts it before redrawing. Neither inherits the other's
            # threshold: show_typing_view stores and restores per kind.
            ("compare", getattr(self, "show_mlst_tree_tab", None)),
            ("cgmlst_tree", getattr(self, "show_cgmlst_tree_tab", None)),
            ("cgmlst", getattr(getattr(self, "cgmlst_calls", None), "refresh", None)),
            ("snp", getattr(getattr(self, "snp_tree_page", None), "refresh_cohort", None)),
            ("reads", getattr(getattr(self, "read_trimming_page", None), "refresh", None)),
            ("assembly", getattr(getattr(self, "assembly_page", None), "refresh", None)),
            # The first time only, and it reads this computer's folders: opening
            # the Update tab must never contact a server or look like an analysis.
            ("update", getattr(getattr(self, "update_center", None), "first_look", None)),
        ) if callable(hook)}
        # The tab widget reports a tab position; navigate speaks keys and
        # build-order numbers, so the bar hands it the key it just showed.
        self.pages.pageShown.connect(self.navigate)
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
        self.apply_saved_resources()
        self.populate_schemes()
        self.refresh()
        self.navigate(0)
        self.add_shortcuts()

    def build_sidebar(self, outer):
        sidebar = QWidget()
        sidebar.setObjectName("sidebar")
        sidebar.setFixedWidth(205)
        self.sidebar = sidebar
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(14, 28, 14, 18)
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
        # The tab bar is the navigation now. nav_buttons stays as an empty list so
        # navigate's checked-state loop keeps working untouched.
        self.nav_buttons = []
        # One name per page, in build order, because the View menu binds Alt+N to
        # the Nth page built. The tab bar shows the same pages in pipeline order.
        self.nav_names = [page_name(key) for key in PAGE_KEYS]
        layout.addWidget(label("Your isolates, your evidence and your report each keep their own cohort. Nothing is included automatically.", "small", True))
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
        self.update_workspace_room()

    def changeEvent(self, event):
        super().changeEvent(event)
        # A larger text size costs the tab bar room without changing the window,
        # so the sidebar decision is made again once the new metrics are polished.
        if event.type() in (QEvent.Type.StyleChange, QEvent.Type.FontChange,
                            QEvent.Type.ApplicationFontChange):
            QTimer.singleShot(0, self.update_workspace_room)

    def update_workspace_room(self):
        """Decide whether this window can afford the sidebar beside the tab bar.

        The tab bar is the navigation, and a tab hidden behind a scroll arrow is a
        tab nobody finds. The sidebar is only afforded when every label can still be
        drawn whole beside it; below that it stands down, as it always did below
        1180 px. Nothing here changes a page, a cohort or a preference.
        """
        try:
            if not hasattr(self, "sidebar") or not hasattr(self, "pages"):
                return None
            room = self.width() - 2 * WORKSPACE_GUTTER - self.sidebar.width()
            afford = self.width() >= 1180 and room >= self.pages.minimum_bar_width()
            self.sidebar.setVisible(afford)
            return afford
        except RuntimeError:  # The window was closed before the queued call ran.
            return None

    def page(self):
        widget = QWidget()
        widget.setObjectName("page")
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(widget)
        slot = self.pages.add_page(scroll)
        key = self.pages.key_at(slot)
        if key:
            position = self.pages.position_of(key)
            self.pages.setTabToolTip(position, f"{page_name(key)} — {PAGE_PURPOSE.get(key, '')}")
            self.page_index[key] = slot
        return widget, layout

    def install_workspace_headers(self):
        """Give every tab the same orientation line once all the pages exist.

        This runs after `build_*` on purpose: the evidence and scheme pages move
        their own layout items by position, and inserting the strip beforehand
        would silently rearrange somebody else's page.
        """
        self.page_headers = {}
        for key, slot in self.page_index.items():
            page = self.pages.widget(slot)
            content = page.widget() if isinstance(page, QScrollArea) else page
            layout = content.layout() if content is not None else None
            if layout is None:
                continue
            strip = workspace_header(self, layout, key, index=0)
            self.page_headers[key] = strip
            self.install_clear_action(key, strip)

    def install_clear_action(self, key, strip):
        """Every tab gets the same Clear, in the same place, saying the same thing.

        A person works with new samples, past samples or a mixture of both, so each
        tab has to be restartable on its own. Clear resets that tab's view, its
        selection and its cohort; it never removes a sample, a result, a scheme or
        an imported report, and the confirmation says which is which.
        """
        row = strip.layout() if strip is not None else None
        if row is None:
            return None
        promise = clear_promise(key)
        action = button("Clear", lambda checked=False, k=key: self.clear_page(k))
        action.setObjectName("nextStep")
        action.setToolTip(f"Start this tab again. Clears {promise['clears']}; "
                          f"keeps {promise['keeps']}.")
        action.setAccessibleName(f"Clear the {page_name(key)} tab")
        # Before the "?" guide button, which is always last on the strip.
        row.insertWidget(max(0, row.count() - 1), action)
        self.clear_buttons[key] = action
        return action

    def fold_duplicate_pages(self):
        """Stop offering a top-level tab whose content is already a sub-tab elsewhere.

        Nothing folds today: each task has its own tab now, which is the point of
        this layout. The mechanism stays because a page that later grows a second
        home must have one place to look, not two. Nothing is folded unless the
        sub-tab really exists, and the page itself is kept: every `navigate(<int>)`
        call site keeps its meaning (ui_tabs.WorkspaceTabs.fold_page).
        """
        from wmlstudio.ui_tabs import FOLDABLE
        folded = []
        for key, (host, title) in FOLDABLE.items():
            tabs = self.pages.subtabs(host)
            if tabs is None or not any(tabs.tabText(index) == title
                                       for index in range(tabs.count())):
                continue
            if self.pages.fold_page(key, host, title):
                folded.append(key)
        return folded

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
        self.heading(layout, "Scheme library",
                     "The local, versioned allele collections this workspace types against, for "
                     "MLST and for cgMLST. Every analysis records the scheme fingerprint it used. "
                     "What is installed, what is published and what an update costs is the Update "
                     "tab; nothing here contacts a server.")
        row = QHBoxLayout()
        row.addWidget(button("Import a scheme folder", self.import_scheme, True))
        row.addWidget(button("What is installed, and what can be updated  →",
                             self.open_update_center))
        row.addStretch()
        row.addWidget(label("No automatic database updates", "badge"))
        layout.addLayout(row)
        # TYPE is its own column because a classical seven-locus scheme and a
        # cgMLST scheme are different quantities: selecting one where the other is
        # expected is wrong, and the two must never read as one alphabetical list.
        self.scheme_table = table(["SCHEME", "TYPE", "SOURCE", "LOCATION"])
        self.install_view_menu("schemes", self.scheme_table)
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
        self.install_view_menu("evidence.hydra", self.hydra_table)
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
        self.heading(layout, "Make yourself at home",
                     "Window size, text size, screen resolution, where your data lives, and a "
                     "straightforward guide to what this version can do.")
        # Page 6 used to be unreachable: the workbench intercepted it and opened a
        # modal dialog instead. It is a real tab now, so the dialog's controls and
        # this page share one implementation (InterfaceSettingsPanel).
        from wmlstudio.interface_settings import InterfaceSettingsPanel
        self.settings_tabs = QTabWidget()
        self.settings_tabs.setObjectName("settingsTabs")
        # "&&" because QTabBar reads a single "&" as a keyboard mnemonic. The tab is
        # named for what a person searches for — the size of the window and of the
        # text — because this is the page they reported as missing.
        self.interface_panel = InterfaceSettingsPanel(self, references=False)
        self.settings_tabs.addTab(self.interface_panel, "Window, text && display")
        self.settings_tabs.setTabToolTip(0, "Window size, text size, screen resolution and "
                                            "how much of this computer analyses may use.")
        data_page = QWidget()
        content = QVBoxLayout(data_page)
        self.motion = QCheckBox("Gentle interface animations")
        self.motion.setChecked(bool(self.project.get_setting("motion", True)))
        self.motion.toggled.connect(self.set_motion)
        content.addWidget(self.motion)
        content.addWidget(label("Data location: " + str(self.root), "small", True))
        content.addWidget(label(f"Active project: {self.project_path}", "small", True))
        for text, method in (("Scheme references…", "open_reference_manager"),
                             ("AMR / plasmid reference databases…", "open_amr_databases"),
                             ("Characterization reference panel…", "install_characterization_references"),
                             ("Open application data folder", "open_data_folder"),
                             ("Save project backup copy…", "save_project_copy"),
                             ("Create a new project…", "new_project")):
            handler = getattr(self, method, None)
            if callable(handler):
                content.addWidget(button(text, handler))
        content.addWidget(label("Updates create immutable reference snapshots; results you already have keep their original provenance.", "small", True))
        # One list of everything installable, so a person never has to know which
        # menu owns which download. The Update tab and the Update menu are the same
        # page; this takes you there rather than opening a second one.
        content.addWidget(button("What is installed, and what can be updated  →",
                                 self.open_update_center))
        content.addWidget(button("Problem → solution guide", self.open_workflow_guide, True))
        content.addStretch()
        self.settings_tabs.addTab(data_page, "Data && references")
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
        self.settings_guide = guide
        self.settings_tabs.addTab(guide, "What this version can and cannot do")
        layout.addWidget(self.settings_tabs, 1)

    def build_update(self):
        """The Update tab: the update centre itself, not a second copy of it.

        For several releases this tab was the scheme library, so the one button
        that looked like an update rescanned this computer's scheme folders and no
        provider was ever asked anything. The page here is the update centre that
        already existed behind a menu (update_center.UpdateCenter): what is
        installed, what is published, what an update costs, and one button each.
        """
        _, layout = self.page()
        self.heading(layout, "Updates and installed reference data",
                     "What is installed here, what each provider publishes today, how much room "
                     "it takes, and what an update would cost.")
        from wmlstudio.update_center import UpdateCenter
        self.update_center = UpdateCenter(self)
        layout.addWidget(self.update_center, 1)
        return self.update_center

    def check_for_updates(self):
        """The Update tab's next step: ask the providers, from wherever it was pressed."""
        self.navigate("update")
        centre = getattr(self, "update_center", None)
        return centre.check_online() if centre is not None else False

    # --- the pipeline stations this layout owns ------------------------------
    def build_stations(self):
        """One tab per remaining step of the workflow, in the order it is worked.

        These are the pages no controller has claimed: read quality, assembly, the
        two typing steps, the cgMLST tree and the SNP tree. Each one is a real page
        with a purpose, a live summary, its own actions and a Clear. The first two
        are handed their finished pages immediately afterwards by
        install_pipeline_pages; a station that cannot do its work yet says so
        instead of looking finished.
        """
        self.station_status = {}
        for key in STATION_KEYS:
            definition = STATION_PAGES.get(key)
            if definition is None:
                continue
            self.build_station(key, definition)
        for key in PLANNED:
            self.pages.mark_planned(key)
        return list(self.stations)

    def build_station(self, key, definition):
        _, layout = self.page()
        self.heading(layout, definition["title"], definition["subtitle"])
        if key in PLANNED:
            # Said in the bar's tooltip and here, where a reader who opened the tab
            # expecting to work cannot miss it. The stretch keeps the badge the width
            # of its own words rather than the width of the page.
            banner = QHBoxLayout()
            banner.addWidget(label("Planned — nothing on this tab runs yet", "badge"))
            banner.addStretch()
            layout.addLayout(banner)
        status = label("", "muted", True)
        status.setObjectName(f"stationStatus.{key}")
        self.station_status[key] = status
        layout.addWidget(status)
        strip = QHBoxLayout()
        for entry in definition.get("actions", ()):
            text, method = entry[0], entry[1]
            needs = entry[2] if len(entry) > 2 else ""
            handler = getattr(self, method, None)
            # An action this build cannot perform is not offered at all: a button
            # that leads nowhere costs more trust than a missing one.
            if callable(handler) and (not needs or getattr(self, needs, None) is not None):
                strip.addWidget(button(text, handler, text.endswith("…")))
        strip.addStretch()
        layout.addLayout(strip)
        card_frame, content = card()
        for paragraph in definition.get("body", ()):
            content.addWidget(label(paragraph, "muted", True))
        layout.addWidget(card_frame)
        # Where a finished page for this station is dropped in, without moving a
        # tab or changing a single page number (see adopt_station).
        holder = QWidget()
        holder_layout = QVBoxLayout(holder)
        holder_layout.setContentsMargins(0, 0, 0, 0)
        holder.hide()
        layout.addWidget(holder, 1)
        layout.addStretch()
        self.stations[key] = {"holder": holder, "layout": holder_layout,
                              "placeholder": card_frame, "status": status}
        return holder

    def install_pipeline_pages(self):
        """Give the two stations after Samples their real pages.

        They are built here rather than inside build_station because they need the
        finished window -- the project, the background task runner and the status
        line -- and because a build that cannot provide one must keep the station's
        own page, which says plainly that the work does not run there.
        """
        self.read_trimming_page = self.assembly_page = None
        try:
            from wmlstudio.ui_reads import install_pipeline_panels
        except ImportError:
            return {}
        return install_pipeline_panels(self)

    def install_typing_stations(self):
        """Give the cgMLST and SNP tree tabs the finished pages they already had.

        Both pages exist and are tested. The SNP tree panel had simply never been
        mounted anywhere, so its tab said "Nothing runs on this tab yet" about a
        page that was complete; the cgMLST calls table was reachable only as a
        sub-tab of the tree page, so the tab named cgMLST showed a signpost to it.
        A tab that names a task and does not carry it is worth less than no tab.
        """
        installed = {}
        calls = getattr(self, "cgmlst_calls", None)
        if calls is not None and self.adopt_station("cgmlst", calls):
            # It was a sub-tab of the comparison page; a widget lives in one place,
            # and its place is the tab that carries its name.
            subtabs = self.pages.subtabs("compare")
            if subtabs is not None:
                index = subtabs.indexOf(calls)
                if index >= 0:
                    subtabs.removeTab(index)
            installed["cgmlst"] = calls
        self.snp_tree_page = None
        try:
            from wmlstudio.ui_snp import SnpTreePanel
            panel = SnpTreePanel(self)
        except (ImportError, TypeError):
            return installed
        self.snp_tree_page = panel
        if self.adopt_station("snp", panel):
            installed["snp"] = panel
        return installed

    def install_tree_station(self):
        """Give the cgMLST tree tab somewhere the comparison workspace can live.

        The reported failure, in the user's words: "cgmlst tree i cannot draw".
        The tab said "Nothing is drawn on this tab yet" and offered a button to
        the other tree page, because there is exactly one comparison workspace —
        its widgets are attributes of this window, so a second one cannot exist —
        and two tabs that each need a tree.

        Rather than leave the core of this application behind a signpost, the
        workspace itself moves into whichever of the two tabs is in front and is
        bound to that tab's typing kind. Each kind already keeps its own
        reference, threshold, minimum overlap, investigation and legend, so the
        two tabs still never share a scale, an axis, a column or a threshold.
        """
        station = self.stations.get("cgmlst_tree")
        if station is None or getattr(self, "comparison_content", None) is None:
            return None
        holder = QScrollArea()
        holder.setWidgetResizable(True)
        self.cgmlst_tree_scroll = holder
        self.adopt_station("cgmlst_tree", holder)
        return holder

    def tree_holders(self):
        """The two scroll areas the one comparison workspace moves between."""
        return {"compare": self.pages.page_for("compare"),
                "cgmlst_tree": getattr(self, "cgmlst_tree_scroll", None)}

    def mount_comparison(self, key):
        """Move the comparison workspace into this tab, bound to this tab's kind.

        QScrollArea.takeWidget hands the content back without destroying it,
        which is what makes one workspace serving two tabs safe rather than a
        trick: the widget is never reparented behind Qt's back and never exists
        in two places.
        """
        content = getattr(self, "comparison_content", None)
        holders = self.tree_holders()
        holder = holders.get(key)
        if content is None or holder is None:
            return False
        if holder.widget() is not content:
            for other in holders.values():
                if other is not None and other.widget() is content:
                    other.takeWidget()
            holder.setWidget(content)
            content.show()
        kind = "cgmlst" if key == "cgmlst_tree" else "mlst"
        show = getattr(self, "show_typing_view", None)
        if callable(show):
            # chosen=False: which tree this tab shows is what the tab is for, and
            # must not overwrite the kind the user last chose for themselves.
            show(kind, chosen=False)
        return True

    def show_mlst_tree_tab(self):
        """The MLST tree tab: the workspace, on the seven-locus scale."""
        self.mount_comparison("compare")
        return self.refresh_comparison()

    def show_cgmlst_tree_tab(self):
        """The cgMLST tree tab: the same workspace, on the core-genome scale."""
        self.mount_comparison("cgmlst_tree")
        return self.refresh_comparison()

    def adopt_station(self, key, widget, *, title=None):
        """Give a pipeline station its real page, built elsewhere.

        The station's own explanation is hidden, the widget takes the page, and the
        tab stops being marked as planned. Nothing about the tab bar, the page
        numbering or any navigation call site changes.
        """
        station = self.stations.get(str(key))
        if station is None or widget is None:
            return False
        station["layout"].addWidget(widget)
        station["holder"].show()
        station["placeholder"].hide()
        station["status"].hide()
        if title:
            position = self.pages.position_of(key)
            if position >= 0:
                self.pages.setTabText(position, str(title))
        self.pages.mark_planned(key, False)
        station["adopted"] = widget
        return True

    def refresh_stations(self):
        """One honest line per station: what this project actually has, counted."""
        status = getattr(self, "station_status", None)
        if not status:
            return {}
        samples = self.current_samples or []
        try:
            index = self.project.typing_kind_index()
        except (sqlite3.Error, ValueError):
            index = {}
        counted = {}
        for kind in ("mlst", "cgmlst"):
            counted[kind] = sum(1 for sample in samples
                                if index.get(sample["id"], {}).get(kind))
        total = len(samples)
        trimmed = sum(1 for sample in samples
                      if (sample.get("metadata", {}).get("read_trimming") or {}).get("status")
                      == "completed")
        assembled = sum(1 for sample in samples if sample.get("metadata", {}).get("assembly", {})
                        .get("assembly_path"))
        lines = {
            "reads": f"{trimmed} of {total} samples carry a trimming record · "
                     f"{total - trimmed} do not. Untrimmed reads stay usable as supplied.",
            "assembly": f"{assembled} of {total} samples carry an assembly made here · "
                        f"{total - assembled} do not.",
            "mlst": f"{counted['mlst']} of {total} samples have a stored seven-locus result · "
                    f"{total - counted['mlst']} have none.",
            "cgmlst": f"{counted['cgmlst']} of {total} samples have a stored cgMLST result · "
                      f"{total - counted['cgmlst']} have none.",
            "cgmlst_tree": f"{counted['cgmlst']} samples carry a cgMLST profile. A tree needs at "
                           "least two profiles from the same scheme.",
            "snp": "Nothing runs on this tab yet.",
        }
        for key, widget in status.items():
            widget.setText(lines.get(key, ""))
        # A station that owns a real page redraws it, but only while it is the page
        # in front: listing a cohort costs a pass over every record, and no tab the
        # user cannot see should pay for one.
        current = self.pages.current_key()
        adopted = (self.stations.get(current) or {}).get("adopted")
        refresh = getattr(adopted, "refresh", None)
        if callable(refresh):
            refresh()
        return counted

    # --- the two pipeline pages between Samples and MLST ----------------------
    def run_read_trimming(self):
        """Trim read pairs on the Read QC tab, wherever the button was pressed."""
        return self.run_pipeline_page("reads", getattr(self, "read_trimming_page", None),
                                      "trim_selected",
                                      "Read trimming needs the Read QC page, which this build "
                                      "could not open.")

    def run_assembly_station(self):
        """Assemble read pairs on the Assembly tab, wherever the button was pressed."""
        return self.run_pipeline_page("assembly", getattr(self, "assembly_page", None),
                                      "assemble_selected",
                                      "Assembly on its own tab needs the Assembly page, which "
                                      "this build could not open. Samples → Assembly still runs "
                                      "the assembler.")

    def run_pipeline_page(self, key, page, action, refusal):
        """Show the page that owns this work, then ask it to do it. Never guesses."""
        runner = getattr(page, action, None)
        if not callable(runner):
            self.notify(refusal)
            return False
        self.navigate(key)
        return runner()

    # --- where each station's work happens today -----------------------------
    def goto_samples(self):
        """The hub, without changing any cohort."""
        self.navigate("isolates")

    def goto_step(self, title, step=""):
        """Samples, with one of its steps selected. Falls back to Samples itself."""
        self.navigate("isolates")
        show = getattr(self, "show_step", None)
        if step and callable(show):
            return bool(show(step))
        return self.pages.show_subtab("isolates", title)

    def goto_assembly_step(self):
        return self.goto_step("Assembly", "assembly")

    def goto_mlst_step(self):
        return self.goto_step("ST", "st")

    def goto_cgmlst_step(self):
        return self.goto_step("cgMLST", "cgmlst")

    def goto_cgmlst_station(self):
        self.navigate("cgmlst")

    def goto_mlst_tree(self):
        self.navigate("compare")

    def goto_cgmlst_calls(self):
        """The cgMLST table of calls, wherever this build keeps it.

        Once the calls page is adopted onto the cgMLST tab this is simply that tab;
        until then it is a sub-tab of the tree page, and the button says so.
        """
        if (self.stations.get("cgmlst") or {}).get("adopted") is not None:
            return self.navigate("cgmlst")
        if getattr(self, "cgmlst_calls", None) is None:
            self.notify("This build has no table of calls yet.")
            return False
        self.navigate("compare")
        return self.pages.show_subtab("compare", "cgMLST calls")

    def goto_cgmlst_tree_view(self):
        """Draw the cgMLST tree, on its own tab, with its own reference and threshold.

        The tab carries the comparison workspace bound to the core-genome scale and
        labels every graph with the scheme and the target count it was built from,
        so a cgMLST graph is never presented on a seven-locus scale. The fallback
        below is for a build where that tab could not be given the workspace.
        """
        if (self.stations.get("cgmlst_tree") or {}).get("adopted") is not None:
            self.navigate("cgmlst_tree")
            return self.typing_kind
        draw = getattr(self, "show_cgmlst_tree", None)
        if not callable(draw):
            self.notify("This build draws cgMLST trees from the tree page.")
            return False
        self.navigate("compare")
        return draw()

    def run_station_step(self, step, what):
        """Run one typing step on the isolates already selected in Samples."""
        runner = getattr(self, "run_step", None)
        if not callable(runner):
            self.notify(f"{what} runs from Samples in this build.")
            return self.goto_samples()
        return runner(step)

    def run_mlst_station(self):
        return self.run_station_step("st", "Seven-locus MLST")

    def run_cgmlst_station(self):
        return self.run_station_step("cgmlst", "cgMLST")

    # --- restarting one tab without losing any evidence ----------------------
    def clear_page(self, key=None, *, confirm=True):
        """Start one tab again: its view, its selection and its cohort, nothing else.

        The user asked for this on every tab so they can work with new samples, past
        samples or a mixture whenever they want. Clearing a view must never delete
        evidence, so the confirmation names what goes and what stays, and this
        method only ever touches interface state.
        """
        key = self.pages.current_key() if key is None else str(key)
        if not key:
            return False
        promise = clear_promise(key)
        question = (f"Start the {page_name(key)} tab again?\n\n"
                    f"This clears {promise['clears']}.\n"
                    f"This keeps {promise['keeps']}.\n\n"
                    "Nothing you have imported or analysed is deleted.")
        if confirm and QMessageBox.question(self, f"Clear {page_name(key)}?", question) \
                != QMessageBox.StandardButton.Yes:
            return False
        self.clear_page_state(key)
        self.notify(f"{page_name(key)} cleared: {promise['clears']}. "
                    f"Kept: {promise['keeps']}.")
        return True

    def clear_page_state(self, key):
        """The interface state one tab owns. Never project data."""
        from wmlstudio.workspace_focus import COHORT_ATTRIBUTES, REFRESH_METHODS
        focus = getattr(self, "focus", None)
        if key in {"overview", "isolates"} and focus is not None and hasattr(focus, "clear"):
            focus.clear()
        if key == "isolates":
            search = getattr(self, "search", None)
            if search is not None:
                search.clear()
            for table_name in ("sample_table", "recent_table"):
                widget = getattr(self, table_name, None)
                if widget is not None:
                    widget.clearSelection()
            for widget in (getattr(self, "step_tables", None) or {}).values():
                widget.clearSelection()
        attribute = COHORT_ATTRIBUTES.get(key)
        if attribute is not None and hasattr(self, attribute):
            setattr(self, attribute, set())
            ledger = getattr(self, "cohort_origins", None)
            if ledger is not None and hasattr(ledger, "forget"):
                ledger.forget(key)
        if key == "compare":
            self.distance_rows = []
            tree = getattr(self, "tree", None)
            if tree is not None:
                tree.draw_results([], [])
            status = getattr(self, "tree_status", None)
            if status is not None:
                status.setText("Cleared. Choose a cohort to build a comparison; every stored "
                               "profile is still in this project.")
        if key == "evidence":
            table_widget = getattr(self, "hydra_table", None)
            if table_widget is not None:
                table_widget.clearSelection()
        if key == "schemes":
            self.populate_schemes()
        if key == "update":
            # Forgets the list and the last check, then reads this computer again.
            # The record of when a check last reached a provider is kept: it is the
            # only thing a failed check has left to tell anybody.
            centre = getattr(self, "update_center", None)
            if centre is not None:
                centre.clear()
        if key == "settings":
            tabs = getattr(self, "settings_tabs", None)
            if tabs is not None:
                tabs.setCurrentIndex(0)
        # A typing station's selection is the step table it runs on, nothing wider.
        step = {"mlst": "st", "cgmlst": "cgmlst"}.get(key)
        table_widget = (getattr(self, "step_tables", None) or {}).get(step)
        if table_widget is not None:
            table_widget.clearSelection()
        refresh = getattr(self, REFRESH_METHODS.get(key, ""), None)
        if callable(refresh):
            refresh()
        self.refresh_stations()
        # A station that has been given a real page clears that page its own way,
        # and does it last: the station refresh above would otherwise redraw the
        # very table the user just asked to be emptied.
        adopted = (self.stations.get(key) or {}).get("adopted")
        if adopted is not None and callable(getattr(adopted, "clear", None)):
            adopted.clear()
        return True

    # --- menus this window owns ----------------------------------------------
    def menu_named(self, title):
        """The menu bar action carrying a menu with this title, or None."""
        for entry in self.menuBar().actions():
            if entry.menu() is not None and entry.text().replace("&", "") == title:
                return entry
        return None

    def register_command(self, action):
        """Make a menu action findable in the Ctrl+K command search, when there is one."""
        commands = getattr(self, "command_actions", None)
        if isinstance(commands, list) and action is not None:
            commands.append((action.text(), action))
        return action

    def extend_display_menu(self):
        """Put window size and the display settings where somebody looks for them.

        The controls existed in Settings and were reported as "not present", so the
        View menu now names them in the words a person searches with, and offers the
        window sizes directly.
        """
        entry = self.menu_named("View")
        menu = entry.menu() if entry is not None else self.menuBar().addMenu("&View")
        menu.addSeparator()
        sizes = menu.addMenu("Window size")
        for width, height in display.available_window_sizes(self.available_screen_size()):
            sizes.addAction(display.describe_window_size(width, height),
                            lambda checked=False, w=width, h=height: self.set_window_size((w, h)))
        sizes.addSeparator()
        sizes.addAction("Fit this screen", lambda: self.set_window_size(None))
        full = sizes.addAction("Full screen", self.toggle_full_screen)
        full.setShortcut(QKeySequence("F11"))
        self.register_command(full)
        self.register_command(menu.addAction(
            "Text too small? Window size, text size and screen resolution…",
            self.open_display_settings))
        return menu

    def open_update_center(self):
        """The Update tab, from a menu, the command search or another page.

        There is one update centre and it is a tab, so a menu entry takes you to
        it rather than opening a second copy in a dialog that could disagree with
        what the tab says.
        """
        from wmlstudio.update_center import open_update_center
        return open_update_center(self)

    def refresh_update_center(self):
        """Tell the Update page that something which installs has stopped.

        Called for every install, wherever it was started from. An install that
        finishes and leaves the row still reading "Not installed" with an
        "Install…" button beside it is indistinguishable from one that never ran,
        which is what a working 24-second download was reported as.
        """
        centre = getattr(self, "update_center", None)
        if centre is None or centre.report is None:
            # Nothing has read this computer yet, so there is no stale list to
            # correct and no reason to start a probe nobody asked for.
            return False
        return centre.install_stopped()

    def build_update_menu(self):
        """A menu for everything that can be installed or updated, and nothing else."""
        from wmlstudio.update_center import MENU_ENTRIES, ROUTES
        bar = self.menuBar()
        menu = QMenu("&Update", self)
        self.register_command(menu.addAction(
            "What is installed, and what can be updated (Update tab)", self.open_update_center))
        menu.addSeparator()
        for key, title in MENU_ENTRIES:
            handler = next((getattr(self, name) for name in ROUTES.get(key, ())
                            if callable(getattr(self, name, None))), None)
            # A thing with no installer on this window is opened in the Update
            # Centre, which states what it is and how it is installed instead.
            self.register_command(menu.addAction(
                title, handler or (lambda checked=False, k=key:
                                   self.open_update_center().start(k))))
        menu.addSeparator()
        for title, method in (("Open the data folder", "open_data_folder"),
                              ("Download practice data…", "download_practice_cohort")):
            handler = getattr(self, method, None)
            if callable(handler):
                menu.addAction(title, handler)
        help_entry = self.menu_named("Help")
        if help_entry is not None:
            bar.insertMenu(help_entry, menu)
        else:
            bar.addMenu(menu)
        self.update_menu = menu
        return menu

    def add_shortcuts(self):
        for shortcut, callback in [("Ctrl+O", self.browse_files), ("Ctrl+Shift+O", self.open_project_dialog), ("Ctrl+S", self.save_project_copy)]:
            action = QAction(self)
            action.setShortcut(QKeySequence(shortcut))
            action.triggered.connect(callback)
            self.addAction(action)

    def navigate(self, target):
        """The single guarded choke point for changing workspace page.

        Both a tab click and a `navigate(<int>)` call arrive here. An integer is a
        build-order page number, the address this application has always used, so
        `navigate(2)` still means the tree page wherever its tab now sits; a string
        is the page key. The guard stops the tab widget's own currentChanged from
        re-entering, and the page's show hook runs exactly once, outside the guard.
        """
        if self._navigating:
            return
        self._navigating = True
        key = ""
        try:
            key = self.pages.key_for(target)
            if not key or not self.pages.show_page(key):
                return
            key = self.pages.current_key() or key
            self.breadcrumb.setText("WORKSPACE  /  " + page_name(key).upper())
            current = self.pages.slot_of(key)
            for i, item in enumerate(self.nav_buttons):
                item.setChecked(i == current)
            page = self.pages.currentWidget()
            if page:
                page.update()
                if isinstance(page, QScrollArea):
                    page.viewport().update()
                    # A tree tab's scroll area is empty until its show hook, below,
                    # moves the one comparison workspace into it.
                    content = page.widget()
                    if content is not None:
                        content.update()
        finally:
            self._navigating = False
        hook = self.page_shown.get(key)
        if callable(hook):
            hook()

    def show_focus_in_library(self):
        """Take the user to the isolates the shared focus names. Changes no cohort."""
        self.navigate("isolates")
        select = getattr(self, "select_samples_in_table", None)
        if callable(select):
            select(sorted(self.focus.ids))

    def preview_interface_scale(self):
        """Show the display controls, where the sample strip previews the size."""
        return self.open_display_settings()

    def open_display_settings(self):
        """Window size, text size and screen resolution, in one obvious place."""
        self.navigate("settings")
        tabs = getattr(self, "settings_tabs", None)
        if tabs is not None:
            tabs.setCurrentIndex(0)
        panel = getattr(self, "interface_panel", None)
        if panel is not None and hasattr(panel, "refresh_preview"):
            panel.refresh_preview()
        return panel

    # --- window size ---------------------------------------------------------
    def available_screen_size(self):
        """The usable area of the screen this window is on, for clamping a size."""
        screen = self.screen()
        if screen is None:
            return (display.MINIMUM_WINDOW_WIDTH, display.MINIMUM_WINDOW_HEIGHT)
        return screen.availableGeometry().size()

    def set_window_size(self, size=None):
        """Resize the window itself. None means 'fit this screen' again.

        This is the window, not the screen: the display resolution belongs to the
        operating system and nothing here changes it.
        """
        width, height = display.clamp_window_size(*(size or (1380, 940)),
                                                  available_size=self.available_screen_size())
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()
        self.resize(width, height)
        self.notify(f"Window size {display.describe_window_size(width, height)}. "
                    + display.WINDOW_SIZE_NOTICE)
        return width, height

    def toggle_full_screen(self):
        """Full screen and back, for a small laptop screen."""
        self.showNormal() if self.isFullScreen() else self.showFullScreen()
        return self.isFullScreen()

    # --- how much of this computer to use ------------------------------------
    def apply_saved_resources(self):
        """Make this workspace's saved jobs × threads the default for new runs."""
        from wmlstudio.interface_settings import interface_preferences, saved_worksize
        from wmlstudio.scheduler import set_default_worksize
        try:
            return set_default_worksize(saved_worksize(interface_preferences(self.root)))
        except (OSError, ValueError):
            return None

    def analysis_allocation(self, memory_gb=1):
        """The admission plan a run should use when it states none of its own."""
        from wmlstudio.scheduler import resources_for_run
        try:
            return resources_for_run({}, memory_gb=memory_gb)
        except ValueError as error:
            self.notify(str(error))
            return None

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
        """Every installed scheme, under a name a person can read, labelled by kind.

        A downloaded cgMLST scheme used to appear as its folder name with the
        underscores turned into spaces, sorted in among 162 classical schemes: the
        user could install one and then not find it. The row title now comes from
        the scheme's own metadata, and TYPE says which library it belongs to.
        """
        selected = self.project.get_setting("scheme_path", "")
        try:
            from wmlstudio import cgmlst_schemes
            cgmlst_schemes.prepare_library(self.root)
            # A scheme moved into the cgMLST library keeps working for a project
            # that stored where it used to be.
            if selected:
                moved = str(cgmlst_schemes.resolve_migrated_path(self.root, selected))
                if moved != selected:
                    self.project.set_setting("scheme_path", moved)
                    selected = moved
        except (OSError, ValueError):
            pass  # A read-only or full data folder must not stop the workspace opening.
        self.scheme_paths = scheme_locations(self.root)
        from wmlstudio.reference_index import scheme_entries as installed_rows
        try:
            titles = {row["path"]: row for row in installed_rows(self.scheme_paths)}
        except (OSError, ValueError):
            titles = {}
        self.scheme_combo.clear()
        self.scheme_combo.addItem("Quality checks only", None)
        self.scheme_table.setRowCount(len(self.scheme_paths))
        for index, path in enumerate(self.scheme_paths):
            row = titles.get(str(path), {})
            name = ("Practice scheme · 7 loci (synthetic)" if path.name == "practice_7"
                    else row.get("title") or path.name)
            kind = {"mlst": "MLST", "cgmlst": "cgMLST"}.get(row.get("kind"), "Not classified")
            self.scheme_combo.addItem(name, str(path))
            if str(path) == selected:
                self.scheme_combo.setCurrentIndex(index + 1)
            for col, text in enumerate([name, kind,
                                        "Local import" if self.root in path.parents else "Bundled snapshot",
                                        str(path)]):
                item = QTableWidgetItem(text)
                item.setToolTip(row.get("kind_conflict") or row.get("kind_basis") or text)
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
        self.refresh_stations()
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
        self.worker = AnalysisWorker(samples, scheme, parent=self,
                                     resource_plan=self.analysis_allocation())
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
        # The organism this run established becomes the sample's own, so it fills
        # in on every other tab without anybody retyping it.
        from wmlstudio.storage import adopt_analysis_organism
        try:
            adopt_analysis_organism(self.project, sample_id, result)
        except (KeyError, ValueError):
            pass  # Recording a label must never lose the result that produced it.
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
            self.recover_moved_storage()
            self.refresh()
            self.navigate(0)
            return True
        except Exception as exc:
            if new_lock is not None and new_lock is not self.project_lock:
                new_lock.unlock()
            self.error(exc)
            return False

    def recover_moved_storage(self):
        """Re-point managed copies after a project folder moves drive or machine.

        Every stored path is absolute, so a project carried on a USB stick shows
        every sample as missing although its copies sit right beside it. A
        candidate is adopted only when its SHA-256 matches the recorded copy: a
        file that merely has the right name is reported, never accepted.
        """
        from wmlstudio import storage
        try:
            missing = storage.missing_managed_copies(self.project)
            if not missing:
                return
            result = storage.relocate_managed_storage(self.project)
        except (OSError, ValueError) as exc:
            self.notify(f"Some sequence copies could not be found: {exc}")
            return
        recovered, unresolved = len(result.get("relocated", ())), len(result.get("unresolved", ()))
        if recovered:
            self.notify(f"Found {recovered} sequence {'copy' if recovered == 1 else 'copies'} "
                        "beside the project after it moved. Nothing was copied or changed."
                        + (f" {unresolved} could not be matched and still need attention." if unresolved else ""))
        elif unresolved:
            self.notify(f"{unresolved} sequence {'copy is' if unresolved == 1 else 'copies are'} "
                        "missing from this project's folder. Re-link them from the isolate record.")

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
        from wmlstudio.export import protected_input_paths
        protected = [self.project_path, *protected_input_paths(self.project.samples())]
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
        # The Update tab's own probe only reads folders, so it never blocks a
        # close; it is asked to stop and waited for, so no thread outlives the
        # window. Anything that downloads runs on the window's worker, below.
        centre = getattr(self, "update_center", None)
        if centre is not None:
            centre.stop()
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


class MainWindow(ContextMenuMixin, WorkbenchMixin, ComparisonWorkspaceMixin, CharacterizationWorkspaceMixin, ReportWorkspaceMixin, BaseWindow):
    """Native workbench composed from focused workflow controllers."""

    def build_menus(self):
        """The workbench's menus, then the two this window owns.

        Update is inserted before Help so the bar reads File … View, Update, Help;
        the display entries extend the View menu the workbench already built.
        """
        super().build_menus()
        self.extend_display_menu()
        self.build_update_menu()


def main(argv=None):
    parser = argparse.ArgumentParser(description="WMLSTudio native desktop")
    parser.add_argument("--project", type=Path)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--screenshot-page", type=int, choices=range(len(PAGE_KEYS)), default=0,
                        help="Workspace page, by the build-order number navigate() uses "
                             f"(0–{len(PAGE_KEYS) - 1}: {', '.join(PAGE_KEYS)})")
    parser.add_argument("--window-size", help="Override screen-aware desktop size with WIDTHxHEIGHT")
    parser.add_argument("--display-scale", type=int, choices=SCALE_CHOICES,
                        help="Magnify the whole interface by this percentage for this run. "
                             "Start once with --display-scale 100 if a saved size made the window unusable.")
    parser.add_argument("--native-screenshot", action="store_true", help="Capture the native window surface, not an offscreen render")
    args = parser.parse_args(argv)
    if args.window_size:
        try:
            width, height = map(int, args.window_size.lower().split("x"))
            if not 1000 <= width <= 7680 or not 680 <= height <= 4320:
                raise ValueError
        except ValueError:
            parser.error("--window-size must be WIDTHxHEIGHT, at least 1000x680 and at most 7680x4320")
    # Both names must be set before the display preferences are read: data_root()
    # resolves the per-user application directory from them, and Qt reads its
    # scaling policy when the QApplication is constructed, not later.
    QCoreApplication.setOrganizationName("IOWA-BioTech")
    QCoreApplication.setApplicationName("WMLSTudio")
    apply_display_settings(args.display_scale)
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
    if args.window_size:
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
