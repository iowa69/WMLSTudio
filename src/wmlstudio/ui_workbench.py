"""Sample-centric native workspace: library, assignments, jobs and application menus."""

import html
import json
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, QUrl
from PySide6.QtGui import QAction, QColor, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLineEdit,
    QSplitter,
    QTabWidget,
    QTextBrowser,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.archive import active_samples, archived_samples, is_archived
from wmlstudio.background import FunctionWorker
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.storage import QUARANTINE_BUCKETS
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    flattened_metadata,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.widgets import DropZone, Helix, Metric, button, card, label


def quarantine_bucket(sample):
    """The needs-review folder this isolate sits in, or None when it is filed.

    Quarantine means WMLSTudio declined to decide what the organism is. It is not
    a claim that the organism is unusual, and it is not the user's own "unknown".
    """
    evidence = (sample.get("metadata") or {}).get("organism_evidence")
    if not isinstance(evidence, dict) or evidence.get("status") != "quarantined":
        return None
    token = str(evidence.get("quarantine_reason") or "awaiting_identification")
    return QUARANTINE_BUCKETS.get(token, QUARANTINE_BUCKETS["awaiting_identification"])


def organism_evidence_note(sample):
    """Where this genus and species came from, in the words the engine used."""
    evidence = (sample.get("metadata") or {}).get("organism_evidence")
    if not isinstance(evidence, dict):
        return ""
    from wmlstudio.organism_id import confidence_label, evidence_sentences
    from wmlstudio.workflow_dialogs import BASIS_LABELS
    word, explanation = confidence_label(evidence)
    basis = BASIS_LABELS.get(evidence.get("basis"), str(evidence.get("basis") or "Not identified"))
    lines = [f"{basis} · {word}", explanation,
             "A folder name is where the copy is stored; it is not a laboratory identification."]
    lines.extend(evidence_sentences(evidence))
    return "\n".join(line for line in lines if line)


class WorkbenchMixin:
    def __init__(self, *args, **kwargs):
        self.selection_ids = set()
        self.cohort_ids = set()
        self.feature_ids = set()
        self.report_ids = set()
        self.library_filter = None
        self.worker_role = ""
        self._filling = False
        self._refreshing = False
        self.reference_dialog = None
        self.amr_database_dialog = None
        self.library = None
        self._pending_graph_state = None
        self._run_plan = {}
        self._run_ids = set()
        self._run_cancelled = False
        self._task_succeeded = False
        self._typing_override = None
        self._progress_dialog = None
        self._pending_import = None
        self._import_notes = []
        self._practice_cohort = None
        self.ui_scale = 100
        super().__init__(*args, **kwargs)
        from wmlstudio import theme
        if hasattr(theme, "apply_dark_palette"):
            theme.apply_dark_palette(QApplication.instance())
        QApplication.instance().setStyleSheet(theme.STYLE)
        from wmlstudio.interface_settings import interface_preferences
        preferences = interface_preferences(self.root)
        try:
            initial_scale = int(preferences.value("scale", 100))
        except (TypeError, ValueError):
            initial_scale = 100
        self.set_ui_scale(initial_scale)
        for combo in self.findChildren(QComboBox):
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(int(combo.property("compactCharacters") or 18))
            combo.setMaximumWidth(360)
        self.setWindowTitle("WMLSTudio · Genomics workbench")
        self.build_menus()
        from wmlstudio.library import Library
        self.library = Library(self.root / "library-index.sqlite")
        self.library.index_project(self.project)
        saved_cohort = self.project.get_setting("comparison_cohort", None)
        if saved_cohort is not None:
            self.cohort_ids = set(saved_cohort)
        if hasattr(self.tree, "restore_state"):
            self._pending_graph_state = self.project.get_setting("graph_style", {"version": 1})
        self.refresh()
        self.statusBar().showMessage("Local genomics workbench · Import isolates, then choose an analysis")

    def build_overview(self):
        _, layout = self.page()
        hero, content = card()
        row = QHBoxLayout()
        words = QVBoxLayout()
        words.addWidget(label("FROM A QUESTION TO REVIEWABLE EVIDENCE", "eyebrow"))
        words.addWidget(label("Your investigation, connected.", "title", True))
        words.addWidget(label("Identify, characterize, compare and report the same isolates. Add new samples without losing earlier work.", "muted", True))
        actions = FlowLayout()
        actions.addWidget(button("Import sequences…", self.browse_files, True))
        actions.addWidget(button("New project…", self.new_project))
        actions.addWidget(button("Research saved library…", self.open_library_research))
        actions.addStretch()
        words.addLayout(actions)
        row.addLayout(words, 4)
        self.helix = Helix()
        self.helix.setMinimumSize(150, 130)
        self.helix.setMaximumWidth(220)
        self.helix.setMaximumHeight(160)
        row.addWidget(self.helix, 1)
        content.addLayout(row)
        layout.addWidget(hero)
        from wmlstudio.journey_widgets import InvestigationMap

        self.overview_tabs = QTabWidget()
        self.investigation_map = InvestigationMap()
        self.investigation_map.actionRequested.connect(self.journey_action)
        self.overview_tabs.addTab(self.investigation_map, "Investigation map")
        library_page = QWidget()
        library_layout = QVBoxLayout(library_page)
        library_layout.setContentsMargins(0, 14, 0, 0)
        self.overview_tabs.addTab(library_page, "Stored library")
        layout.addWidget(self.overview_tabs, 1)
        metrics = QHBoxLayout()
        self.metrics = [Metric("Isolates", "in the active library"), Metric("Analysed", "saved evidence"),
                        Metric("Exact MLST", "registered profiles"), Metric("Review", "incomplete or failed")]
        for metric in self.metrics:
            metrics.addWidget(metric)
        library_layout.addLayout(metrics)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        explorer, content = card()
        content.addWidget(label("Library navigator", "cardTitle"))
        self.library_tree = QTreeWidget()
        self.library_tree.setHeaderHidden(True)
        self.library_tree.setMinimumWidth(250)
        self.library_tree.itemActivated.connect(self.open_library_group)
        self.library_tree.itemClicked.connect(self.open_library_group)
        self.install_view_menu("library.tree", self.library_tree)
        content.addWidget(self.library_tree)
        splitter.addWidget(explorer)
        recent, content = card()
        content.addWidget(label("Recent samples · double-click to inspect", "cardTitle"))
        self.recent_table = make_table(["Sample", "Input", "Status", "ST", "Loci"])
        self.recent_table.cellDoubleClicked.connect(self.open_recent_sample)
        self.install_view_menu("overview.recent", self.recent_table)
        content.addWidget(self.recent_table)
        self.drop_zone = DropZone()
        self.drop_zone.filesDropped.connect(lambda paths: self.import_paths(paths, configure=True))
        self.drop_zone.browseRequested.connect(self.browse_files)
        content.addWidget(self.drop_zone)
        splitter.addWidget(recent)
        splitter.setSizes([320, 650])
        library_layout.addWidget(splitter, 1)
        self.practice_notice = label("Practice project · synthetic sequences, not biological isolates.", "small")
        self.practice_notice.hide()
        layout.addWidget(self.practice_notice)

    def build_samples(self):
        _, layout = self.page()
        self.heading(layout, "Isolate library", "Browse your stored isolates. Double-click a record for complete evidence; each analysis chooses its own cohort.")
        row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search names, organisms, STs, AMR genes, collections or metadata…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_tables)
        row.addWidget(self.search, 1)
        row.addWidget(button("Import…", self.browse_files, True))
        row.addWidget(button("Assign organism…", self.assign_selected))
        row.addWidget(button("Epidemiology grid…", self.open_metadata_grid))
        layout.addLayout(row)
        filters = FlowLayout()
        self.genus_filter = QComboBox()
        self.species_filter = QComboBox()
        self.collection_filter = QComboBox()
        for combo, title in [(self.genus_filter, "All genera"), (self.species_filter, "All species"),
                             (self.collection_filter, "All collections")]:
            combo.addItem(title, "")
            combo.currentIndexChanged.connect(self.refresh_tables)
            filters.addWidget(combo, 1)
        filters.addWidget(button("Clear filters", self.clear_filters))
        layout.addLayout(filters)
        run, content = card()
        strip = FlowLayout()
        self.scheme_combo = QComboBox()
        self.scheme_combo.setMinimumWidth(180)
        self.scheme_combo.setAccessibleName("Typing scheme override")
        self.scheme_combo.hide()  # Per-isolate schemes are reviewed in the analysis dialog.
        self.run_button = button("Analyse…", self.choose_and_analyse, True)
        self.rerun_button = button("Review pending…", lambda: self.choose_and_analyse(pending_only=True))
        strip.addWidget(self.run_button)
        strip.addWidget(self.rerun_button)
        strip.addWidget(button("Attach reads…", self.attach_reads_selected))
        content.addLayout(strip)
        self.selection_label = label("No samples selected · each sample keeps its own organism and typing workflow", "small", True)
        content.addWidget(self.selection_label)
        layout.addWidget(run)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.sample_table = make_table(["Sample", "Input", "Status", "ST", "Called loci", "Genus", "Species",
                                        "Evidence", "Scheme", "AMR genes", "Collection", "Storage"])
        self.sample_table.itemSelectionChanged.connect(self.sample_selection_changed)
        self.sample_table.cellDoubleClicked.connect(lambda row, column: self.open_isolate_record(self.sample_table.item(row, 0).data(Qt.ItemDataRole.UserRole)))
        for column in (7, 8, 10, 11):
            self.sample_table.setColumnHidden(column, True)
        self.install_view_menu("library", self.sample_table)
        splitter.addWidget(self.sample_table)
        self.detail = QTextBrowser()
        self.detail.setOpenExternalLinks(False)
        self.detail.setMinimumHeight(120)
        inspector = QTabWidget()
        inspector.addTab(self.detail, "Sample evidence")
        self.history_view = QTextBrowser()
        inspector.addTab(self.history_view, "History / provenance")
        inspector.setParent(self)
        inspector.hide()  # Complete evidence is available in the isolate popup.
        layout.addWidget(splitter, 1)

    def refresh(self):
        if self._refreshing:
            return
        self._refreshing = True
        try:
            super().refresh()
            # Archived isolates keep every record they carry; they simply leave
            # the working views until the user restores them.
            working = active_samples(self.current_samples)
            archived = len(self.current_samples) - len(working)
            isolates = [sample for sample in working if sample.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
            mates = len(working) - len(isolates)
            self.metrics[0].value.setText(str(len(isolates)))
            hint = "in the active library"
            if mates:
                hint = f"{len(isolates)} isolates · {mates} linked mates"
            if archived:
                hint += f" · {archived} archived"
            self.metrics[0].hint.setText(hint)
            if self.library is not None and not (self.worker and self.worker.isRunning()):
                self.library.index_project(self.project)
            valid = {s["id"] for s in working}
            self.selection_ids.intersection_update(valid)
            self.report_ids.intersection_update(valid)
            self.focus.prune(valid)
            self.refresh_library()
            if hasattr(self, "cohort_table"):
                self.refresh_cohort_table()
            if hasattr(self, "report_table"):
                self.refresh_report_table()
            if hasattr(self, "feature_table"):
                self.refresh_features()
            self.refresh_journey()
        finally:
            self._refreshing = False

    def refresh_journey(self):
        if not hasattr(self, "investigation_map"):
            return
        from wmlstudio.journey import journey_summary
        summary = journey_summary(active_samples(self.current_samples), selected_ids=self.selection_ids,
                                  comparison_ids=self.cohort_ids)
        investigation = self.investigation_summary() if hasattr(self, "investigation_summary") else {}
        description = None
        if investigation.get("saved"):
            description = (f"Investigation: {investigation.get('name')} · {investigation.get('snapshots', 0)} saved snapshots · "
                           "New samples are included only after you review the cohort.")
        self.investigation_map.set_summary(summary, description)
        if hasattr(self, "scope_label"):
            counts = summary["counts"]
            self.scope_label.setText(f"{counts['total']} isolates in project")

    def journey_action(self, action):
        """Route problem-oriented actions without silently broadening a cohort."""
        if action == "import":
            self.browse_files()
        elif action == "guide":
            self.open_workflow_guide()
        elif action == "samples":
            self.navigate(1)
        elif action == "review":
            from wmlstudio.journey import journey_summary
            identifiers = journey_summary(active_samples(self.current_samples))["review_ids"]
            self.clear_filters()
            self.selection_ids = set(identifiers)
            self.library_filter = ("ids", set(identifiers))
            self.refresh_tables()
            self.navigate(1)
            self.notify(f"{len(identifiers)} flagged isolates shown. Clear filters to return to all samples.")
        elif action == "analyse":
            self.choose_and_analyse()
        elif action == "characterize":
            self.run_characterization_selected()
        elif action == "features":
            self.choose_feature_cohort()
        elif action == "compare":
            self.navigate(2)
        elif action == "reports":
            self.navigate(5)

    def refresh_library(self):
        tree = self.library_tree
        tree.blockSignals(True)
        tree.clear()
        working = active_samples(self.current_samples)
        all_item = QTreeWidgetItem([f"All samples  ({len(working)})"])
        all_item.setData(0, Qt.ItemDataRole.UserRole, None)
        tree.addTopLevelItem(all_item)
        groups = {}
        collections = {}
        review = {}
        for sample in working:
            if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate":
                continue
            bucket = quarantine_bucket(sample)
            if bucket:
                # Needs review is not an organism, so it never joins the genus tree.
                review.setdefault(bucket, []).append(sample)
                continue
            genus, species, _ = organism_for(sample)
            genus, species = genus or "Unknown", species or "Unspecified"
            st = (sample.get("result") or {}).get("st")
            st_label = f"ST {st}" if st else "ST unassigned"
            groups.setdefault(genus, {}).setdefault(species, {}).setdefault(st_label, []).append(sample)
            for name in (sample.get("metadata") or {}).get("collections", []):
                collections.setdefault(name, []).append(sample)
        for genus, species_groups in sorted(groups.items()):
            genus_node = QTreeWidgetItem([genus])
            genus_node.setData(0, Qt.ItemDataRole.UserRole, ("genus", genus))
            all_item.addChild(genus_node)
            for species, st_groups in sorted(species_groups.items()):
                species_node = QTreeWidgetItem([species])
                species_node.setData(0, Qt.ItemDataRole.UserRole, ("organism", genus, species))
                genus_node.addChild(species_node)
                for st, samples in sorted(st_groups.items()):
                    node = QTreeWidgetItem([f"{st}  ({len(samples)})"])
                    node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                    species_node.addChild(node)
                    for sample in samples:
                        child = QTreeWidgetItem([sample["name"]])
                        child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                        node.addChild(child)
        if review:
            total = sum(len(samples) for samples in review.values())
            branch = QTreeWidgetItem([f"Needs review  ({total})"])
            branch.setData(0, Qt.ItemDataRole.UserRole,
                           ("ids", [s["id"] for samples in review.values() for s in samples]))
            branch.setToolTip(0, "WMLSTudio declined to decide what these organisms are, so nothing "
                                 "was filed into a genus folder on a guess. Right-click to assign an "
                                 "organism; the managed copy moves with the label.")
            tree.addTopLevelItem(branch)
            for bucket, samples in sorted(review.items()):
                node = QTreeWidgetItem([f"{bucket.replace('_', ' ')}  ({len(samples)})"])
                node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                branch.addChild(node)
                for sample in samples:
                    child = QTreeWidgetItem([sample["name"]])
                    child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                    child.setToolTip(0, organism_evidence_note(sample))
                    node.addChild(child)
            branch.setExpanded(True)
        if collections:
            branch = QTreeWidgetItem(["Collections"])
            tree.addTopLevelItem(branch)
            for name, samples in sorted(collections.items()):
                node = QTreeWidgetItem([f"{name}  ({len(samples)})"])
                node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                branch.addChild(node)
            branch.setExpanded(True)
        if self.library is not None:
            project_counts = {}
            for sample in self.library.search():
                path = sample["project_path"]
                if Path(path).resolve() != self.project_path.resolve():
                    project_counts[path] = project_counts.get(path, 0) + 1
            if project_counts:
                projects = QTreeWidgetItem(["Other saved projects"])
                tree.addTopLevelItem(projects)
                for path, count in sorted(project_counts.items()):
                    node = QTreeWidgetItem([f"{Path(path).stem}  ({count})"])
                    node.setData(0, Qt.ItemDataRole.UserRole, ("project", path))
                    projects.addChild(node)
                projects.setExpanded(True)
        all_item.setExpanded(True)
        stored = archived_samples(self.current_samples)
        if stored:
            branch = QTreeWidgetItem([f"Archived  ({len(stored)})"])
            branch.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in stored]))
            branch.setToolTip(0, "Archived isolates keep every result, allele call, file and history "
                                 "entry. They are hidden from the working views until you restore them.")
            tree.addTopLevelItem(branch)
            for sample in stored:
                child = QTreeWidgetItem([sample["name"]])
                child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                branch.addChild(child)
        mates = [sample for sample in active_samples(self.current_samples) if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate"]
        if mates:
            mate_branch = QTreeWidgetItem([f"Linked read mates  ({len(mates)})"])
            mate_branch.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"] for sample in mates]))
            tree.addTopLevelItem(mate_branch)
        tree.blockSignals(False)

    def open_library_group(self, item, column=0):
        value = item.data(0, Qt.ItemDataRole.UserRole)
        if value and value[0] == "project":
            self.switch_project(value[1])
            return
        self.library_filter = value
        self.refresh_tables()
        self.navigate(1)

    def open_recent_sample(self, row, column=0):
        item = self.recent_table.item(row, 0)
        if item:
            self.open_isolate_record(item.data(Qt.ItemDataRole.UserRole))

    @staticmethod
    def fill_filter(combo, values, title):
        selected = combo.currentData()
        combo.blockSignals(True)
        combo.clear()
        combo.addItem(title, "")
        for value in sorted(set(v for v in values if v)):
            combo.addItem(value, value)
        index = combo.findData(selected)
        combo.setCurrentIndex(max(0, index))
        combo.blockSignals(False)

    def refresh_tables(self):
        if not hasattr(self, "sample_table"):
            return
        working = active_samples(self.current_samples)
        if self.library_filter and self.library_filter[0] == "ids":
            # An explicitly named group — the Archived branch, for instance — shows
            # the isolates it names, so a hidden isolate is never simply missing.
            named = set(self.library_filter[1])
            working = working + [s for s in archived_samples(self.current_samples) if s["id"] in named]
        taxa = [organism_for(s) for s in working]
        self.fill_filter(self.genus_filter, [t[0] for t in taxa], "All genera")
        self.fill_filter(self.species_filter, [t[1] for t in taxa], "All species")
        self.fill_filter(self.collection_filter, [c for s in working for c in
                         (s.get("metadata") or {}).get("collections", [])], "All collections")
        query = self.search.text().strip().casefold()
        visible = []
        for sample in working:
            genus, species, evidence = organism_for(sample)
            metadata = sample.get("metadata") or {}
            result = sample.get("result") or {}
            searchable = " ".join([sample["name"], genus, species, str(result.get("st") or ""),
                                   " ".join(gene_names(sample)), json.dumps(flattened_metadata(sample), ensure_ascii=False)])
            if query and query not in searchable.casefold():
                continue
            if self.genus_filter.currentData() and genus != self.genus_filter.currentData():
                continue
            if self.species_filter.currentData() and species != self.species_filter.currentData():
                continue
            if self.collection_filter.currentData() and self.collection_filter.currentData() not in metadata.get("collections", []):
                continue
            if self.library_filter:
                kind, *values = self.library_filter
                if kind == "ids" and sample["id"] not in values[0]:
                    continue
                if kind == "genus" and (genus or "Unknown") != values[0]:
                    continue
                if kind == "organism" and ((genus or "Unknown"), (species or "Unspecified")) != tuple(values):
                    continue
            visible.append(sample)
        self.fill_sample_table(self.sample_table, visible)
        self.fill_sample_table(self.recent_table, active_samples(self.current_samples)[-12:])
        self.selection_label.setText(self.selection_summary(len(visible)))
        self.show_sample_detail()

    def selection_summary(self, visible):
        working = active_samples(self.current_samples)
        archived = len(self.current_samples) - len(working)
        review = sum(1 for sample in working if quarantine_bucket(sample))
        text = f"{len(self.selection_ids)} selected · {visible} visible · {len(working)} in project"
        if review:
            text += f" · {review} in Needs review"
        if archived:
            text += f" · {archived} archived"
        return text

    def fill_sample_table(self, widget, samples):
        self._filling = True
        widget.blockSignals(True)
        widget.setSortingEnabled(False)
        widget.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            result = sample.get("result") or {}
            metadata = sample.get("metadata") or {}
            status = result.get("status", sample["status"]) if sample["status"] == "completed" else sample["status"]
            alleles = result.get("alleles", {})
            kind = result.get("kind") or ("fastq" if ".fq" in sample["input_path"] or ".fastq" in sample["input_path"] else "fasta")
            if metadata.get("workflow", {}).get("source_kind") == "profile":
                kind = "profile"
            if metadata.get("workflow", {}).get("source_kind") == "read_mate":
                kind, status = "read_mate", "linked_mate"
            genus, species, evidence = organism_for(sample)
            bucket = quarantine_bucket(sample)
            note = organism_evidence_note(sample)
            storage = "Managed copy" if metadata.get("workflow", {}).get("managed") else "Linked original"
            if sample.get("missing_input") and kind != "profile":
                storage = "Input unavailable"
            if bucket:
                genus, species = genus or "Needs review", species or "—"
                evidence = f"Needs review · {bucket.replace('_', ' ')}"
            if is_archived(sample):
                evidence = f"Archived · {evidence}"
            values = [sample["name"], {"fastq": "Reads", "profile": "Profile only", "read_mate": "Linked reverse mate"}.get(kind, "Assembly"),
                      status.replace("_", " ").title(), f"ST {result['st']}" if result.get("st") else "—",
                      f"{sum(v is not None for v in alleles.values())} / {len(alleles)}" if alleles else "—",
                      genus or "Unknown", species or "—", evidence, result.get("scheme") or "—",
                      "; ".join(gene_names(sample)) or "—", "; ".join(metadata.get("collections", [])) or "—", storage]
            for column, value in enumerate(values[:widget.columnCount()]):
                item = cell(value, sample["id"])
                if column == 2:
                    item.setForeground(QColor("#5AE0BD" if status in {"complete", "qc_only", "completed"}
                                               else "#F5BE73" if status in {"incomplete", "mixed", "failed", "interrupted"}
                                               else "#A6B7D0"))
                if metadata.get("cluster", {}).get("highlight"):
                    item.setBackground(QColor("#343348"))
                if note and column in (5, 6, 7):
                    item.setToolTip(note)
                widget.setItem(row, column, item)
                if widget is self.sample_table and sample["id"] in self.selection_ids:
                    item.setSelected(True)
        widget.setSortingEnabled(True)
        widget.blockSignals(False)
        self._filling = False

    def sample_selection_changed(self):
        if self._filling:
            return
        visible = {self.sample_table.item(r, 0).data(Qt.ItemDataRole.UserRole)
                   for r in range(self.sample_table.rowCount())}
        selected = {item.data(Qt.ItemDataRole.UserRole) for item in self.sample_table.selectedItems()}
        self.selection_ids.difference_update(visible)
        self.selection_ids.update(selected)
        self.focus.set_focus(self.selection_ids, "Isolate table selection")
        self.selection_label.setText(self.selection_summary(len(visible)))
        self.show_sample_detail()

        self.refresh_journey()

    def select_visible_samples(self):
        self.sample_table.selectAll()
        self.sample_selection_changed()

    def clear_sample_selection(self):
        self.selection_ids.clear()
        self.sample_table.clearSelection()
        self.sample_selection_changed()

    def clear_filters(self):
        self.library_filter = None
        self.search.clear()
        for combo in (self.genus_filter, self.species_filter, self.collection_filter):
            combo.setCurrentIndex(0)
        self.refresh_tables()

    def selected_samples(self):
        return [s for s in self.project.samples() if s["id"] in self.selection_ids]

    def show_sample_detail(self, sample_id=None):
        sample = self.project.get_sample(sample_id) if sample_id else self.selected_sample()
        if not sample:
            self.detail.setHtml("<h3>Select an isolate to inspect its evidence</h3><p>Use Ctrl/Shift to select a cohort. Organism, typing, AMR and annotations stay linked by sample ID.</p>")
            if hasattr(self, "history_view"):
                self.history_view.clear()
            return
        def e(value):
            return html.escape(str(value))
        result = sample.get("result") or {}
        metadata = sample.get("metadata") or {}
        genus, species, evidence = organism_for(sample)
        body = f"<h2>{e(sample['name'])}</h2><p><b>{e(' '.join([genus, species]).strip() or 'Unknown organism')}</b> · {e(evidence)} · ST {e(result.get('st') or 'unassigned')}</p>"
        body += f"<p><b>Input:</b> {e(sample['input_path'] or 'Imported profile; no sequence attached')}<br><b>Sample ID:</b> {e(sample['id'])}</p>"
        if sample.get("error"):
            body += f"<p><b>Needs attention:</b> {e(sample['error'])}</p>"
        workflow = metadata.get("workflow", {})
        if workflow.get("paired_with"):
            try:
                paired = self.project.get_sample(workflow["paired_with"])
                body += f"<p><b>Paired record:</b> {e(paired['name'])} · {e(paired['id'])}<br><b>Read role:</b> {e(workflow.get('source_kind'))}</p>"
            except KeyError:
                body += "<p>The paired record is not present in this project; its provenance is retained.</p>"
        body += f"<p><b>Workflow:</b> {e(workflow.get('typing_mode', 'manual'))} · <b>Scheme:</b> {e(result.get('scheme') or workflow.get('scheme_path') or 'Not assigned')}</p>"
        if workflow.get("managed"):
            body += f"<p><b>Managed copy:</b> original remains at {e(workflow.get('source_path', 'recorded in provenance'))}</p>"
        qc = result.get("qc") or {}
        metrics = [("Records", "records"), ("Bases", "total_bases"), ("N50", "n50"), ("GC %", "gc_percent"), ("Q30 %", "q30_percent")]
        body += "<p>" + " · ".join(f"<b>{title}:</b> {e(round(qc[key], 2) if isinstance(qc[key], float) else qc[key])}" for title, key in metrics if qc.get(key) is not None) + "</p>"
        if qc.get("sampled"):
            body += "<p><b>Sampled read QC:</b> statistics describe an inspected prefix, not the entire read file.</p>"
        if metadata.get("hydra"):
            from wmlstudio.sample_workflow import hydra_evidence_status
            amr_state = hydra_evidence_status(sample)
            body += f"<h3>Linked HYDRA evidence · {e(amr_state['status'])}</h3><p>{e(amr_state['reason'])}</p>"
            if amr_state["status"] != "stale":
                body += "<p>" + e("; ".join(gene_names(sample)) or "No primary AMR genes reported in linked evidence") + "</p>"
        annotations = {k: v for k, v in metadata.items() if k not in {"workflow", "organism", "hydra", "assembly"}}
        if annotations:
            body += "<h3>Annotations</h3><p>" + "<br>".join(f"<b>{e(k)}:</b> {e(v)}" for k, v in annotations.items()) + "</p>"
        calls = result.get("calls") or []
        if calls:
            body += "<h3>Allele calls</h3><table width='100%' cellpadding='6'><tr bgcolor='#253650'><th>Locus</th><th>Allele</th><th>Evidence</th></tr>"
            for call in calls[:100]:
                body += f"<tr><td>{e(call['locus'])}</td><td>{e(call.get('allele') or '—')}</td><td>{e(call.get('status'))}</td></tr>"
            body += "</table>"
            if len(calls) > 100:
                body += f"<p>Showing 100 of {len(calls)} loci. The profile table and exports retain every locus.</p>"
        for note in result.get("notes", []):
            body += f"<p>{e(note)}</p>"
        body += f"<p><b>Input SHA-256:</b> {e(result.get('input_sha256') or workflow.get('source_sha256') or 'Not analysed')}<br><b>Scheme fingerprint:</b> {e(result.get('scheme_digest') or 'Not used')}</p>"
        self.detail.setHtml(body)
        if hasattr(self, "history_view"):
            history = self.project.history(sample["id"])
            rows = []
            for entry in reversed(history[-100:]):
                details = entry["details"]
                brief = {key: value for key, value in details.items() if not isinstance(value, (dict, list))}
                archived = details.get("result")
                if isinstance(archived, dict):
                    brief.update(archived_scheme=archived.get("scheme"), archived_st=archived.get("st"), archived_input_sha256=archived.get("input_sha256"))
                rows.append(f"<h3>{e(entry['action'].replace('_', ' '))}</h3><p>{e(entry['created_at'])}</p><p>" + "<br>".join(f"<b>{e(key)}:</b> {e(value)}" for key, value in brief.items()) + "</p>")
            self.history_view.setHtml(f"<h2>Sample audit trail · {e(sample['name'])}</h2><p>Showing {min(100, len(history))} of {len(history)} saved events. Original inputs and archived analysis evidence are retained.</p>" + "".join(rows))

    def scheme_entries(self):
        return [(p.name.replace("_", " "), str(p)) for p in self.scheme_paths]

    def populate_schemes(self):
        super().populate_schemes()
        self.scheme_combo.setItemText(0, "Use each sample's workflow")

    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Import sequences", "", "FASTA / FASTQ (*.fa *.fasta *.fna *.fq *.fastq *.gz *.bz2);;All files (*)")
        if paths:
            self.import_paths(paths, configure=True)

    def browse_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Import sequence folder")
        if path:
            self.import_paths([path], configure=True)

    def dropEvent(self, event):
        self.import_paths([u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()], configure=True)
        event.acceptProposedAction()

    def import_paths(self, paths, configure=False):
        if not configure:
            # Direct, non-interactive API retained for scripted integration tests.
            return super().import_paths(paths)
        if self.busy():
            return
        files, seen = [], set()
        for source in paths:
            path = Path(source)
            # Sorted, so a dropped folder is reviewed in the order the user sees it
            # in their own file manager rather than in filesystem order.
            candidates = sorted(path.rglob("*")) if path.is_dir() else [path]
            for candidate in candidates:
                name = candidate.name.lower().removesuffix(".gz").removesuffix(".bz2")
                if candidate.is_file() and name.endswith((".fa", ".fasta", ".fna", ".fq", ".fastq")):
                    resolved = str(candidate.resolve())
                    if resolved not in seen:
                        seen.add(resolved)
                        files.append(resolved)
        if not files:
            self.error("No supported FASTA or FASTQ files were found.")
            return
        from wmlstudio.workflow_dialogs import ImportSamplesDialog
        dialog = ImportSamplesDialog(files, self.scheme_entries(), self)
        dialog.storage_root.setText(str(self.project_path.with_suffix(".files")))
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if any(assignment.get("typing_mode") == "auto" for assignment in dialog.assignments):
            self.identify_before_import(dialog.assignments, dialog.options)
            return
        self.import_assignments(dialog.assignments, dialog.options)

    def installed_species_panel(self):
        """The broad ANI panel the user installed, or None when only the starter is present."""
        from wmlstudio import organism_panel, paths
        try:
            return organism_panel.installed_species_panel(paths.data_root())
        except OSError:
            return None

    def organism_suggestions(self):
        """Organism names worth offering: installed references first, then this project."""
        pairs = {}
        try:
            from wmlstudio.organism_panel import installed_species_panel
            from wmlstudio.paths import data_root
            from wmlstudio.reference_index import panel_entries
            from wmlstudio.reference_index import scheme_entries as installed_scheme_entries
            panel = installed_species_panel(data_root())
            entries = list(installed_scheme_entries(self.scheme_paths))
            entries += list(panel_entries(panel)) if panel else []
            for entry in entries:
                if entry.get("genus"):
                    pairs.setdefault((str(entry["genus"]), str(entry.get("species") or "")), None)
        except (OSError, ValueError):
            pass  # Suggestions are a convenience; a missing reference never blocks import.
        for sample in self.current_samples:
            genus, species, _ = organism_for(sample)
            if genus:
                pairs.setdefault((genus, species or ""), None)
        return sorted(pairs)

    def identify_before_import(self, assignments, options):
        """Identify the user's own files first, so nothing is copied to a folder it must leave.

        Identification runs on the originals, before a single byte is copied, so a
        cancelled or rejected review leaves the project exactly as it was.
        """
        from wmlstudio.characterization_refs import bundled_reference_root
        from wmlstudio.scheduler import resources_for_run
        originals = [assignment["path"] for assignment in assignments
                     if assignment.get("typing_mode") == "auto"]
        panel = self.installed_species_panel()
        scheme_paths = list(self.scheme_paths)
        starter = bundled_reference_root()

        def operation(cancelled, progress):
            from wmlstudio.organism_id import identify_batch
            return identify_batch(originals, allocation=resources_for_run({"threads": 2, "memory_gb": 2}),
                                  species_panel_root=panel, kpsc_panel_root=starter,
                                  scheme_paths=scheme_paths, cancelled=cancelled, progress=progress)

        self._pending_import = {"assignments": list(assignments), "options": dict(options),
                                "verdicts": []}
        if not self.launch_task(operation, "identify", self.identification_completed):
            self._pending_import = None

    def identification_completed(self, verdicts):
        if self._pending_import is not None:
            self._pending_import["verdicts"] = list(verdicts)

    def review_identification(self):
        """Show what was proposed, and copy only what the user accepts."""
        pending, self._pending_import = self._pending_import, None
        if not pending:
            return
        if not self._task_succeeded:
            self.notify("Identification stopped. Nothing was imported and your files are unchanged.")
            return
        from wmlstudio.workflow_dialogs import IdentificationReviewDialog
        options = pending["options"]
        root = options.get("storage_root") or str(self.project_path.with_suffix(".files"))
        names = {str(Path(a["path"]).expanduser().resolve()): a.get("name")
                 for a in pending["assignments"] if a.get("name")}
        dialog = IdentificationReviewDialog(pending["verdicts"], root, self,
                                            organisms=self.organism_suggestions(), names=names,
                                            policy=self.project)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            self.notify("Nothing was imported. Your files are where you left them and no folders "
                        "were created.")
            return
        reviewed = {str(Path(entry["path"]).expanduser().resolve()): entry
                    for entry in dialog.assignments}
        assignments = []
        for assignment in pending["assignments"]:
            if assignment.get("typing_mode") != "auto":
                assignments.append(assignment)
                continue
            chosen = reviewed.get(str(Path(assignment["path"]).expanduser().resolve()))
            if chosen is not None:
                # The review decides the organism, which is what files the copy. How
                # the typing scheme is chosen is a separate question the user already
                # answered in the import dialog, so their answer is kept.
                assignments.append({**assignment, **chosen,
                                    "typing_mode": assignment.get("typing_mode", "auto")})
        if not assignments:
            self.notify("No files were selected for import. Nothing was copied.")
            return
        self.import_assignments(assignments, options, duplicates="skip")

    def import_assignments(self, assignments, options, *, duplicates="allow"):
        from wmlstudio.storage import import_samples
        root = options.get("storage_root") or str(self.project_path.with_suffix(".files"))
        notes = []
        self._import_notes = notes
        def operation(cancelled, progress):
            return import_samples(
                self.project, assignments, storage_root=root, managed=options.get("managed", True),
                append_st=options.get("append_st", False), cancelled=cancelled, progress=progress,
                duplicates=duplicates, notes=notes)
        self.launch_task(operation, "import", lambda ids: self.import_completed(ids))

    def import_completed(self, ids):
        skipped = list(getattr(self, "_import_notes", ()) or ())
        self._import_notes = []
        self.selection_ids = set(ids)
        self.clear_filters()
        self.refresh()
        self.navigate(1)
        chosen = set(ids)
        review = [sample for sample in self.current_samples
                  if sample["id"] in chosen and quarantine_bucket(sample)]
        message = f"Imported {len(ids)} samples."
        if skipped:
            message += (f" {len(skipped)} file(s) already in this project were not copied again.")
        if review:
            message += (f" {len(review)} are in Needs review: no installed reference supported an "
                        "organism, so nothing was filed into a genus folder on a guess.")
        else:
            message += " Review their assignments, then analyse when ready."
        self.notify(message)

    def refile_completed(self, report):
        """Say what moved, what did not, and why — never a bare success."""
        self.refresh()
        moved, skipped = len(report.get("moved", ())), report.get("skipped", ())
        message = (f"{moved} managed copies were filed to match their organism. "
                   "Your original files were not moved.")
        if skipped:
            message += f" {len(skipped)} were left where they are: {skipped[0][1]}"
        self.notify(message)

    def apply_organism_assignments(self, assignments, *, notice=""):
        """Record corrected organisms, then move the managed copies to match them.

        A corrected label that did not move the file would leave the file sitting
        in a folder that contradicts it, so the two always travel together.
        """
        from wmlstudio.storage import confirm_organism, refile_samples
        assignments = [dict(assignment) for assignment in assignments if assignment.get("sample_id")]
        if not assignments:
            return False

        def operation(cancelled, progress):
            report = {"moved": [], "unchanged": [], "skipped": []}
            for index, assignment in enumerate(assignments):
                progress(index, len(assignments), f"Filing {index + 1} of {len(assignments)}…")
                sample_id = assignment["sample_id"]
                confirm_organism(self.project, [sample_id], assignment.get("genus", ""),
                                 assignment.get("species", ""),
                                 scheme_path=assignment.get("scheme_path"),
                                 typing_mode=assignment.get("typing_mode", "manual"))
                # Each sample keeps its own storage location; correcting a label
                # never moves a managed copy into a different root.
                outcome = refile_samples(self.project, [sample_id], cancelled=cancelled)
                for key in report:
                    report[key].extend(outcome[key])
            return report

        started = self.launch_task(operation, "refile", self.refile_completed)
        if started and notice:
            self.notify(notice)
        return started

    def assign_selected(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if not samples:
            self.notify("Select one or more sample rows first.")
            return
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        dialog = BatchAssignmentDialog(samples, self.scheme_entries(), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.apply_organism_assignments(dialog.assignments)

    def refile_selected(self):
        return self.context_refile(self.library_selection())

    def context_refile(self, selection):
        """Show where each managed copy would go, then move only what the user confirms."""
        if self.busy():
            return
        from PySide6.QtWidgets import QMessageBox

        from wmlstudio.storage import plan_filing, refile_samples
        identifiers = list(selection.sample_ids)
        if not identifiers:
            self.notify("Select the isolates whose managed copies should be re-filed.")
            return
        plans = [plan_filing(self.project, sample_id) for sample_id in identifiers]
        moving = [plan for plan in plans if plan["eligible"] and plan["changed"]]
        if not moving:
            reasons = {plan["reason"] for plan in plans if plan["reason"]}
            self.notify("Every selected managed copy is already in the folder its organism says. "
                        + (" ".join(sorted(reasons)) if reasons else ""))
            return
        lines = [f"{plan['name']}  →  {plan['relative']}" for plan in moving[:10]]
        if len(moving) > 10:
            lines.append(f"…and {len(moving) - 10} more.")
        if QMessageBox.question(self, "Re-file managed copies?",
                f"{len(moving)} managed copies will move to match the organism recorded for them. "
                "Your original files are not moved, and a folder name is a filing decision rather "
                "than a laboratory identification.\n\n" + "\n".join(lines)) != QMessageBox.StandardButton.Yes:
            return
        ids = [plan["sample_id"] for plan in moving]

        def operation(cancelled, progress):
            return refile_samples(self.project, ids, cancelled=cancelled, progress=progress)

        self.launch_task(operation, "refile", self.refile_completed)

    def download_practice_cohort(self):
        """Fetch a pinned teaching cohort after its caveats have been read."""
        if self.busy():
            return
        from wmlstudio import paths, practice_cohorts
        from wmlstudio.workflow_dialogs import PracticeCohortDialog
        cohorts = practice_cohorts.describe_cohorts()
        destinations = {entry["name"]: practice_cohorts.default_destination(paths.data_root(), entry["name"])
                        for entry in cohorts}
        dialog = PracticeCohortDialog(cohorts, destinations, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        name = dialog.chosen
        destination = destinations[name]

        def operation(cancelled, progress):
            return practice_cohorts.download_cohort(name, destination, cancelled=cancelled,
                                                    progress=progress)

        self._practice_cohort = None
        self.launch_task(operation, "practice_cohort",
                         lambda result: setattr(self, "_practice_cohort", result))

    def practice_cohort_ready(self):
        """Point the user at the downloaded files and start the ordinary import."""
        result, self._practice_cohort = self._practice_cohort, None
        if not result or not self._task_succeeded:
            return
        path = Path(result["path"])
        self.notify(f"{result['genomes']} practice genomes are ready in {path}. Every file was "
                    "checked against the checksum NCBI publishes for it.")
        self.import_paths([str(path)], configure=True)

    def install_species_panel(self):
        """Download the broader reference panel identification needs, on request only."""
        if self.busy():
            return
        from PySide6.QtWidgets import QMessageBox

        from wmlstudio import organism_panel, paths
        data_root = paths.data_root()
        installed = organism_panel.installed_species_panel(data_root)
        megabytes = sum(entry[5] for entry in organism_panel.SPECIES_PANEL) / (1024 * 1024)
        question = (f"Download {len(organism_panel.SPECIES_PANEL)} reference genomes "
                    f"({megabytes:.1f} MB) from NCBI RefSeq to {data_root}?\n\n"
                    "One reference per organism is a triage panel, not a representation of "
                    "within-species diversity. Nothing is uploaded and no genome leaves this "
                    "computer.")
        if installed:
            question = f"A panel is already installed at {installed}.\n\n" + question
        if QMessageBox.question(self, "Install broader species panel", question) != QMessageBox.StandardButton.Yes:
            return
        root = data_root / organism_panel.PANEL_DIRECTORY

        def operation(cancelled, progress):
            return organism_panel.provision_species_panel(root, cancelled=cancelled, progress=progress)

        self.launch_task(operation, "species_panel", lambda result: self.notify(
            f"Species panel installed: {result['species_count']} references at {result['path']}. "
            "New imports will be compared against it."))

    def busy(self):
        if self.worker and self.worker.isRunning():
            self.notify("A background task is active. You can navigate and inspect results; cancel or finish before changing inputs.")
            return True
        dialog = self.amr_database_dialog
        if dialog and dialog.worker and dialog.worker.isRunning():
            self.notify("A reference snapshot is being downloaded. Finish or cancel that task before starting another analysis.")
            return True
        return False

    def launch_task(self, operation, role, completed=None):
        if self.busy():
            return False
        self.worker_role = role
        self._task_succeeded = False
        self.worker = FunctionWorker(operation, self)
        self.worker.completed.connect(lambda result: setattr(self, "_task_succeeded", True))
        if completed:
            self.worker.completed.connect(completed)
        self.worker.failed.connect(self.error)
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.set_running(True)
        self.worker.start()
        return True

    def set_running(self, running):
        for widget in (self.run_button, self.rerun_button, self.scheme_combo, self.demo_button):
            widget.setEnabled(not running)
        self.cancel_button.setVisible(running)
        self.progress_bar.setVisible(running)
        if running:
            self.progress_bar.setValue(0)
            from wmlstudio.analysis_progress import AnalysisProgressDialog
            if self._progress_dialog is None or self._progress_dialog.finished_safely:
                self._progress_dialog = AnalysisProgressDialog(self)
            self._progress_dialog.show()
        elif self._progress_dialog is not None:
            self._progress_dialog.finish()

    def job_progress(self, percent, message):
        super().job_progress(percent, message)
        if self._progress_dialog is not None and not self._progress_dialog.finished_safely:
            self._progress_dialog.update_progress(percent, message)

    def start_analysis(self, checked=False, all_samples=False, selected_only=False, confirm=False, assemble=False, sample_ids=None):
        if self.busy():
            return
        chosen_ids = set(sample_ids) if sample_ids is not None else None
        samples = [s for s in self.project.samples() if
                   (s["id"] in chosen_ids if chosen_ids is not None else
                    s["id"] in self.selection_ids if selected_only else
                    all_samples or s["status"] in {"queued", "failed", "interrupted"})]
        samples = [s for s in samples if s.get("input_path") and s.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
        if not samples:
            self.notify("No sequence inputs in the requested selection. Import files, select rows or use saved profiles in Compare.")
            return
        scheme = self.scheme_combo.currentData()
        if scheme:
            samples = [dict(s, metadata={**s.get("metadata", {}), "workflow": {
                **s.get("metadata", {}).get("workflow", {}), "typing_mode": "manual", "scheme_path": scheme}})
                for s in samples]
        plan = self.review_run_plan(samples, scheme, assemble=assemble) if confirm else {}
        if plan is None:
            return
        self._run_plan = plan
        self._run_ids = {sample["id"] for sample in samples}
        self._run_cancelled = False
        self._typing_override = scheme
        pairs = None
        if plan.get("assemble"):
            from wmlstudio.pairing_dialog import PairReadsDialog
            reads = [sample for sample in samples if Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fq", ".fastq"))]
            if len(reads) < 2:
                self.notify("Select at least two paired FASTQ inputs to assemble, or disable assembly and run QC only.")
                self._run_plan = {}
                return
            pairing = PairReadsDialog(reads, self)
            if pairing.exec() != QDialog.DialogCode.Accepted:
                self._run_plan = {}
                return
            pairs = pairing.assignments
        if plan.get("fastqc"):
            self.begin_fastqc_plan(samples, scheme, pairs)
        elif pairs:
            self.assemble_pairs(pairs, plan)
        else:
            self.begin_typing(samples, scheme)

    def begin_fastqc_plan(self, samples, scheme, pairs=None):
        from wmlstudio.fastqc_dialog import reads_for_sample, run_project_fastqc
        from wmlstudio.scheduler import resources_for_run
        reads = [sample for sample in samples if reads_for_sample(sample)]
        if not reads:
            self.notify("No original reads are linked to the reviewed inputs; FastQC was not run.")
            self._run_plan = {}
            return
        try:
            allocation = resources_for_run(self._run_plan, memory_gb=1)
        except ValueError as exc:
            self.error(str(exc))
            self._run_plan = {}
            return
        self._fastqc_pending = {"ids": [sample["id"] for sample in samples], "scheme": scheme, "pairs": pairs}
        self._fastqc_flagged = 0
        project = self.project
        output = self.project_path.with_suffix(".files") / "fastqc_reports"
        def operation(cancelled, progress):
            return run_project_fastqc(project, reads, output, allocation, cancelled=cancelled, progress=progress)
        def completed(results):
            self._fastqc_flagged = sum(report["qc_status"] == "FAIL" for entry in results
                                      for report in entry["result"]["reports"])
            self.notify(f"Original FastQC completed for {len(results)} inputs; {self._fastqc_flagged} read reports have FAIL flags. No reads were trimmed.")
        self.launch_task(operation, "fastqc_pipeline", completed)

    def continue_after_fastqc(self):
        pending = getattr(self, "_fastqc_pending", None)
        self._fastqc_pending = None
        if not pending or self._run_cancelled or not self._task_succeeded:
            self._run_plan = {}
            return
        if self._fastqc_flagged:
            from PySide6.QtWidgets import QMessageBox
            answer = QMessageBox.question(self, "Review FastQC flags before continuing",
                f"{self._fastqc_flagged} read reports contain FastQC FAIL flags. Reports are saved with the isolates. "
                "These flags are not an automatic instruction to trim, but may affect downstream interpretation. "
                "Continue the reviewed analysis on the unchanged reads?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                self._run_plan = {}
                self.notify("Stopped after FastQC for review. Saved reports and original reads are retained.")
                return
        if pending["pairs"]:
            self.assemble_pairs(pending["pairs"], self._run_plan)
            return
        samples = [self.project.get_sample(identifier) for identifier in pending["ids"]]
        scheme = pending["scheme"]
        self.begin_typing(samples, scheme)

    def begin_typing(self, samples, scheme):
        from wmlstudio.scheduler import resources_for_run
        try:
            allocation = resources_for_run(self._run_plan)
        except ValueError as exc:
            self.error(str(exc))
            return
        if scheme:
            samples = [dict(sample, metadata={**sample.get("metadata", {}), "workflow": {
                **sample.get("metadata", {}).get("workflow", {}), "typing_mode": "manual", "scheme_path": scheme}})
                for sample in samples]
        self.project.set_setting("scheme_path", scheme or "")
        self.worker_role = "analysis"
        self.worker = AnalysisWorker(samples, scheme, parent=self, installed_scheme_paths=self.scheme_paths,
                                     resource_plan=allocation)
        self.worker.sample_started.connect(self.sample_started)
        self.worker.sample_finished.connect(self.sample_finished)
        self.worker.sample_failed.connect(self.sample_failed)
        self.worker.sample_cancelled.connect(lambda sid: self.project.set_status(sid, "interrupted", "Cancelled by user"))
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.set_running(True)
        self.worker.start()

    def assemble_pairs(self, pairs, plan):
        import uuid

        from wmlstudio.assembly import associate_assembly, run_skesa
        from wmlstudio.scheduler import resources_for_run, run_bounded
        project = self.project
        try:
            allocation = resources_for_run(plan, memory_gb=8)
        except ValueError as exc:
            self.error(str(exc))
            return
        tasks = [{"primary": project.get_sample(pair["primary_id"]),
                  "mate": project.get_sample(pair["mate_id"])} for pair in pairs]
        output_root = self.project_path.with_suffix(".files") / "assemblies"
        def operation(cancelled, progress):
            identifiers = []
            def assemble(task, resources, stopped, report):
                primary, mate = task["primary"], task["mate"]
                destination = output_root / primary["id"] / uuid.uuid4().hex
                return run_skesa(primary["input_path"], mate["input_path"], destination,
                    threads=resources.threads_per_sample, memory_gb=resources.memory_gb, cancelled=stopped,
                    progress=lambda done, total, message: report(done, total, f"{primary['name']} · {message}"))
            def attach(task, result):
                primary, mate = task["primary"], task["mate"]
                associate_assembly(project, primary["id"], mate["id"], result)
                identifiers.append(primary["id"])
            run_bounded(tasks, assemble, allocation, cancelled=cancelled, on_result=attach, progress=progress)
            return identifiers
        self.launch_task(operation, "assembly", lambda ids: self.notify(f"Assembled {len(ids)} read pairs. Starting the assigned typing workflow…"))

    def analysis_finished(self):
        role = self.worker_role
        self.worker_role = ""
        self.set_running(False)
        self.refresh()
        if self.closing_after_cancel:
            self.close()
            return
        if role == "fastqc_pipeline":
            self.continue_after_fastqc()
            return
        if role == "identify":
            self.review_identification()
            return
        if role == "practice_cohort":
            self.practice_cohort_ready()
            return
        if role == "assembly":
            if self._task_succeeded and not self._run_cancelled:
                samples = [sample for sample in self.project.samples() if sample["id"] in self._run_ids
                           and sample.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
                self.begin_typing(samples, self._typing_override)
            else:
                self._run_plan = {}
            return
        if role == "analysis" and not self._run_cancelled:
            managed = [s for s in self.project.samples() if s["id"] in self._run_ids and s["status"] == "completed"
                       and s.get("metadata", {}).get("workflow", {}).get("managed")]
            if managed:
                from wmlstudio.storage import refile_samples
                identifiers = [sample["id"] for sample in managed]

                def organize(cancelled, progress):
                    return refile_samples(self.project, identifiers, cancelled=cancelled,
                                          progress=progress)

                self.launch_task(organize, "organize", self.refile_completed)
                return
        if role in {"analysis", "organize"}:
            if self._run_plan.get("hydra") and not self._run_cancelled and (role != "organize" or self._task_succeeded):
                self.run_hydra_plan(self._run_ids, self._run_plan)
            self._run_plan = {}

    def cancel_analysis(self):
        self._run_cancelled = True
        self._run_plan = {}
        super().cancel_analysis()

    def active_amr_database(self):
        selected = self.project.get_setting("hydra_database_root", "")
        if selected:
            return selected
        from wmlstudio import hydra_runtime
        bundled = getattr(hydra_runtime, "bundled_database_root", lambda: None)()
        return str(bundled or self.root / "references" / "hydra")

    def review_run_plan(self, samples, scheme=None, hydra=False, assemble=False):
        from wmlstudio.analysis_plan import RunPlanDialog
        dialog = RunPlanDialog(samples, scheme, self.active_amr_database(), self)
        dialog.hydra.setChecked(hydra)
        dialog.assemble.setChecked(assemble)
        if hydra:
            dialog.assemble.setEnabled(False)
        dialog.manageDatabases.connect(lambda: self.manage_run_databases(dialog))
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        if dialog.plan.get("hydra"):
            self.project.set_setting("hydra_database_root", dialog.plan["db_root"])
        return dialog.plan

    def manage_run_databases(self, plan_dialog):
        from wmlstudio.amr_databases import AMRDatabaseDialog
        database_dialog = AMRDatabaseDialog(self.active_amr_database(), plan_dialog,
            update_root=self.root / "references" / "hydra")
        self.amr_database_dialog = database_dialog
        database_dialog.snapshotInstalled.connect(plan_dialog.database.setText)
        database_dialog.snapshotInstalled.connect(lambda path: self.project.set_setting("hydra_database_root", path))
        database_dialog.exec()

    def run_hydra_selected(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if not samples:
            self.notify("Select assembly samples in the Samples menu first.")
            return
        plan = self.review_run_plan(samples, hydra=True)
        if plan and plan.get("hydra"):
            self._run_cancelled = False
            self.run_hydra_plan({sample["id"] for sample in samples}, plan, require_completed=False)

    def run_hydra_plan(self, identifiers, plan, require_completed=True):
        from wmlstudio.hydra_runtime import run_assemblies
        from wmlstudio.sample_workflow import link_hydra
        from wmlstudio.scheduler import resources_for_run, run_bounded
        samples = [sample for sample in self.project.samples() if sample["id"] in identifiers
                   and sample.get("input_path") and (not require_completed or sample["status"] == "completed")]
        assemblies = [sample for sample in samples if (sample.get("result") or {}).get("kind", sample.get("kind")) == "fasta"
                      or Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fa", ".fasta", ".fna"))]
        if not assemblies:
            self.notify("No FASTA assemblies are ready for HYDRA in this selection. Raw reads must be assembled first.")
            return
        try:
            allocation = resources_for_run(plan, memory_gb=3)
        except ValueError as exc:
            self.error(str(exc))
            return

        def operation(cancelled, progress):
            import hashlib

            reports = []
            def analyse(sample, resources, stopped, report_progress):
                genus, species, _ = organism_for(sample)
                assigned = sample.get("metadata", {}).get("organism", {})
                # A provisional MLST lineage alone is not a verified mutation-catalog assignment.
                organism = " ".join([genus, species]) if assigned.get("genus") and assigned.get("species") else None
                return run_assemblies([sample["input_path"]], plan["db_root"], plan.get("databases"),
                    sample_names=[sample["id"]], organism=organism, threads=resources.threads_per_sample,
                    protein=plan.get("protein", True), cancelled=stopped,
                    point_mutations=plan.get("point_mutations", True),
                    **plan.get("thresholds", {}),
                    progress=lambda done, total, message: report_progress(done, total, f"{sample['name']} · {message}"))
            def attach(sample, report):
                link_hydra(self.project, report, {sample["id"]: sample["id"]})
                reports.append(report)
                combined = dict(report, samples=[entry for result in reports for entry in result["samples"]])
                components = [result["import_provenance"] for result in reports]
                combined["import_provenance"] = {"source_path": "Native HYDRA run", "components": components,
                    "sha256": hashlib.sha256(json.dumps(components, sort_keys=True).encode()).hexdigest()}
                combined["execution_provenance"] = {"sample_runs": [result.get("execution_provenance", {}) for result in reports]}
                self.project.set_setting("hydra_report", combined)
            run_bounded(assemblies, analyse, allocation, cancelled=cancelled, on_result=attach, progress=progress)
            return len(reports)

        self.launch_task(operation, "hydra", lambda count: self.notify(f"HYDRA completed for {count} assemblies. AMR evidence is linked to the sample library and reports."))

    def open_amr_databases(self):
        if self.busy():
            return
        from wmlstudio.amr_databases import AMRDatabaseDialog
        if self.amr_database_dialog is None:
            self.amr_database_dialog = AMRDatabaseDialog(self.active_amr_database(), self,
                update_root=self.root / "references" / "hydra")
            self.amr_database_dialog.snapshotInstalled.connect(lambda path: self.project.set_setting("hydra_database_root", path))
        self.amr_database_dialog.show()
        self.amr_database_dialog.raise_()
        self.amr_database_dialog.activateWindow()

    def switch_project(self, path):
        if self.worker and self.worker.isRunning():
            return super().switch_project(path)
        old_path = self.project_path.resolve()
        if self.library is not None:
            self.library.index_project(self.project)
        changed = super().switch_project(path)
        if changed and old_path != self.project_path.resolve():
            self.selection_ids.clear()
            self.report_ids.clear()
            cohort = self.project.get_setting("comparison_cohort", None)
            self.cohort_ids = set(cohort) if cohort is not None else set()
            self.feature_ids = set()
            for dialog in getattr(self, "isolate_dialogs", {}).values():
                dialog.close()
            self.isolate_dialogs = {}
            self.library_filter = None
            if hasattr(self, "restore_investigations"):
                self.restore_investigations()
            if hasattr(self, "tree") and hasattr(self.tree, "restore_state"):
                self._pending_graph_state = self.project.get_setting("graph_style", {"version": 1})
            self.refresh()
        return changed

    def edit_metadata(self):
        samples = self.selected_samples()
        if not samples:
            sample = self.selected_sample()
            samples = [sample] if sample else []
        if not samples:
            self.notify("Select sample rows first.")
            return
        annotations = samples[0].get("metadata", {}).get("annotations", {}) if len(samples) == 1 else {}
        text, accepted = QInputDialog.getMultiLineText(self, "Sample annotations", "field = value, one per line. These annotations apply to every selected sample; workflow and AMR evidence are preserved.", "\n".join(f"{key} = {value}" for key, value in annotations.items()))
        if not accepted:
            return
        try:
            values = {}
            for line in text.splitlines():
                if line.strip():
                    key, separator, value = line.partition("=")
                    if not separator or not key.strip():
                        raise ValueError("Use field = value on each line.")
                    values[key.strip()] = value.strip()
            for sample in samples:
                metadata = dict(sample.get("metadata") or {})
                metadata["annotations"] = {**metadata.get("annotations", {}), **values}
                self.project.set_metadata(sample["id"], metadata)
            self.refresh()
        except Exception as exc:
            self.error(exc)

    def add_collection(self):
        if not self.selection_ids:
            self.notify("Select the samples to add to a collection first.")
            return
        name, accepted = QInputDialog.getText(self, "Add to collection", "Collection name (for example, Ward 5 · September):")
        if accepted and name.strip():
            collection_id = self.project.create_collection(name.strip())
            self.project.set_collection_members(collection_id, self.selection_ids, add=True)
            for sample in self.selected_samples():
                metadata = dict(sample.get("metadata") or {})
                metadata["collections"] = sorted(set(metadata.get("collections", []) + [name.strip()]))
                self.project.set_metadata(sample["id"], metadata)
            self.refresh()

    def open_data_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project_path.parent)))

    def open_library_research(self):
        if self.busy():
            return
        from wmlstudio.library_dialog import LibraryDialog
        self.library.index_project(self.project)
        dialog = LibraryDialog(self.library, self.project, self)
        dialog.profilesImported.connect(self.import_completed)
        dialog.exec()

    def open_sample_folder(self):
        sample = self.selected_sample()
        if sample and sample.get("input_path"):
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(sample["input_path"]).parent)))

    def relink_selected_input(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if len(samples) != 1:
            self.notify("Select exactly one sample to relink. The replacement must match its recorded input SHA-256.")
            return
        sample = samples[0]
        path, _ = QFileDialog.getOpenFileName(self, f"Relink {sample['name']} · identical input bytes required", "",
            "Sequence files (*.fasta *.fa *.fna *.fastq *.fq *.gz *.bz2);;All files (*)")
        if not path:
            return
        from wmlstudio.storage import relink_input
        self.launch_task(lambda cancelled, progress: relink_input(self.project, sample["id"], path, cancelled=cancelled),
            "relink", lambda result: self.notify("Input relinked after SHA-256 verification. Saved profiles and original files are unchanged."))

    # --- right-click handlers -----------------------------------------------
    # One method per entry in context_menus.ACTIONS that acts on isolates. A view
    # is wired with a single install_view_menu call; an action whose handler is
    # absent is left out of the menu rather than offered and then failing.

    def context_records(self, selection):
        wanted = set(selection.sample_ids)
        return [sample for sample in self.project.samples() if sample["id"] in wanted]

    def context_open_record(self, selection):
        if selection.single:
            self.open_isolate_record(selection.single)

    def context_select_in_library(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.clear_filters()
        self.selection_ids = set(ids)
        self.library_filter = ("ids", ids)
        self.focus.set_focus(ids, "Right-click selection")
        self.refresh_tables()
        self.navigate(1)
        self.notify(f"{len(ids)} isolates shown. Clear filters to return to the whole library.")

    def update_comparison_cohort(self, ids, add=True):
        """Change the comparison cohort explicitly, and say where the change came from."""
        ids = {value for value in ids if value}
        if not ids:
            return
        self.cohort_ids = (set(self.cohort_ids) | ids) if add else (set(self.cohort_ids) - ids)
        self.project.set_setting("comparison_cohort", sorted(self.cohort_ids))
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("compare", self.cohort_ids, "Right-click in the isolate library")
        if hasattr(self, "cohort_table"):
            self.refresh_cohort_table()
        self.refresh_journey()
        verb = "added to" if add else "removed from"
        self.notify(f"{len(ids)} isolates {verb} the comparison cohort · {len(self.cohort_ids)} in it now. "
                    "Nothing else was included automatically.")

    def context_add_to_comparison(self, selection):
        self.update_comparison_cohort(selection.sample_ids, add=True)

    def context_remove_from_comparison(self, selection):
        self.update_comparison_cohort(selection.sample_ids, add=False)

    def context_add_to_characterization(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.feature_ids = set(self.feature_ids) | ids
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("evidence", self.feature_ids, "Right-click in the isolate library")
        if hasattr(self, "feature_table"):
            self.refresh_features()
        self.notify(f"{len(ids)} isolates added to the evidence cohort · {len(self.feature_ids)} in it now.")

    def context_add_to_report(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.report_ids = set(self.report_ids) | ids
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("reports", self.report_ids, "Right-click in the isolate library")
        if hasattr(self, "report_table"):
            self.refresh_report_table()
        self.notify(f"{len(ids)} isolates added to the report · {len(self.report_ids)} in it now.")

    def context_remove_from_report(self, selection):
        ids = set(selection.sample_ids)
        self.report_ids = set(self.report_ids) - ids
        if hasattr(self, "report_table"):
            self.refresh_report_table()
        self.notify(f"{len(ids)} isolates removed from the report · {len(self.report_ids)} in it now.")

    @staticmethod
    def label_assignment(sample, genus, species):
        """An organism correction that keeps the sample's own typing workflow."""
        workflow = (sample.get("metadata") or {}).get("workflow", {})
        mode = workflow.get("typing_mode") or "manual"
        if mode == "unknown" and genus:
            mode = "manual"
        return {"sample_id": sample["id"], "genus": genus, "species": species,
                "typing_mode": mode, "scheme_path": workflow.get("scheme_path")}

    def context_assign_organism(self, selection):
        if self.busy():
            return
        samples = self.context_records(selection)
        if not samples:
            return
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        dialog = BatchAssignmentDialog(samples, self.scheme_entries(), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.apply_organism_assignments(dialog.assignments)

    def context_assign_organism_quick(self, selection, organism):
        if self.busy():
            return
        genus, species = organism
        samples = self.context_records(selection)
        self.apply_organism_assignments([self.label_assignment(sample, genus, species)
                                         for sample in samples])

    def context_assign_scheme(self, selection):
        if self.busy():
            return
        entries = self.scheme_entries()
        samples = self.context_records(selection)
        if not entries or not samples:
            self.notify("No typing schemes are installed. Use Data ▸ Online scheme catalog first.")
            return
        names = [name for name, _ in entries]
        title, accepted = QInputDialog.getItem(self, "Assign typing scheme",
            f"Type these {len(samples)} isolates with:", names, 0, False)
        if not accepted:
            return
        path = entries[names.index(title)][1]
        assignments = []
        for sample in samples:
            genus, species, _ = organism_for(sample)
            if not genus:
                self.notify(f"{sample['name']} has no organism yet. Assign a genus first, or use "
                            "Unknown organism in the assignment dialog.")
                return
            assignments.append({"sample_id": sample["id"], "genus": genus, "species": species,
                                "typing_mode": "manual", "scheme_path": path})
        self.apply_organism_assignments(assignments)

    def context_rename_sample(self, selection):
        sample_id = selection.single
        if not sample_id or self.busy():
            return
        sample = self.project.get_sample(sample_id)
        name, accepted = QInputDialog.getText(self, "Rename isolate",
            "Display name. The sample identifier, its input file and every stored result stay "
            "exactly as they are:", text=sample["name"])
        if not accepted:
            return
        try:
            self.project.rename_sample(sample_id, name)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.notify("Renamed. No result, allele call or file was changed.")

    def context_rename_folder(self, selection):
        """Re-label every isolate in an organism folder, and move their copies with it."""
        if self.busy() or not selection.folder:
            return
        genus, species = selection.folder
        samples = self.context_records(selection)
        if not samples:
            self.notify("That folder holds no isolates to re-label.")
            return
        new_genus, accepted = QInputDialog.getText(self, "Rename this organism folder",
            f"Genus recorded for these {len(samples)} isolates. A folder name is a filing decision, "
            "not a laboratory identification:", text="" if genus == "Unknown" else genus)
        if not accepted:
            return
        new_species, accepted = QInputDialog.getText(self, "Rename this organism folder",
            "Species (leave empty if you only know the genus):",
            text="" if species in {"Unspecified", ""} else species)
        if not accepted:
            return
        self.apply_organism_assignments([self.label_assignment(sample, new_genus.strip(), new_species.strip())
                                         for sample in samples])

    def context_edit_annotations(self, selection):
        self.with_selection(selection, self.edit_metadata)

    def context_add_collection(self, selection):
        self.with_selection(selection, self.add_collection)

    def context_highlight(self, selection):
        self.with_selection(selection, self.highlight_selected)

    def with_selection(self, selection, action):
        """Run an existing selection-driven action on exactly what was right-clicked."""
        previous = self.selection_ids
        self.selection_ids = set(selection.sample_ids)
        try:
            action()
        finally:
            self.selection_ids = previous
        self.refresh_tables()

    def context_unhighlight(self, selection):
        from wmlstudio.sample_workflow import set_cluster
        identifiers = list(selection.sample_ids)
        if not identifiers:
            return
        for sample in self.context_records(selection):
            group = (sample.get("metadata") or {}).get("cluster", {})
            set_cluster(self.project, [sample["id"]], group.get("label") or "Highlighted",
                        group.get("color") or "#2F8A78", False)
        self.refresh()
        self.notify(f"Highlight removed from {len(identifiers)} isolates. The grouping you recorded "
                    "stays in each isolate's history.")

    def context_add_files(self, selection):
        self.browse_files()

    def context_add_folder(self, selection):
        self.browse_folder()

    def context_import_scheme(self, selection):
        self.import_scheme()

    def context_open_scheme_folder(self, selection):
        folder = next((Path(path) for path in selection.paths), None)
        if folder is None or not folder.is_dir():
            self.notify("That scheme folder is not available on this computer.")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    def context_remove_scheme(self, selection):
        import shutil

        from PySide6.QtWidgets import QMessageBox
        if self.busy():
            return
        root = Path(self.root).resolve()
        paths = [Path(path).resolve() for path in selection.paths]
        removable = [path for path in paths if root in path.parents and path.is_dir()]
        if not removable or len(removable) != len(paths):
            self.notify("Bundled reference snapshots are read-only; nothing was removed.")
            return
        names = ", ".join(path.name for path in removable)
        if QMessageBox.question(self, "Remove imported scheme?",
                f"Delete {names} from this computer? Results already produced with it keep their "
                "scheme fingerprint in the project, but you cannot rerun them until it is installed "
                "again.") != QMessageBox.StandardButton.Yes:
            return
        for path in removable:
            try:
                shutil.rmtree(path)
            except OSError as error:
                self.error(error)
                break
        self.populate_schemes()
        self.notify(f"Removed {len(removable)} imported scheme folders. Saved results are unchanged.")

    def context_archive(self, selection):
        from wmlstudio import archive
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        reason, accepted = QInputDialog.getText(self, f"Archive {len(identifiers)} isolates",
            "Archiving hides these isolates from the working views and destroys nothing: every "
            "result, allele call, file and history entry is kept, and you can restore them at any "
            "time.\n\nWhy are you archiving them? (optional)")
        if not accepted:
            return
        try:
            changed = archive.archive_samples(self.project, identifiers, reason.strip())
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        already = len(identifiers) - len(changed)
        self.notify(f"{len(changed)} isolates archived and hidden from the working views."
                    + (f" {already} were already archived." if already else "")
                    + " Nothing was deleted; find them under Archived in the library navigator.")

    def context_restore(self, selection):
        from wmlstudio import archive
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        try:
            changed = archive.restore_samples(self.project, identifiers)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.notify(f"{len(changed)} isolates restored to the working views with every record they "
                    "carried.")

    def context_remove(self, selection):
        from PySide6.QtWidgets import QCheckBox, QMessageBox

        from wmlstudio.storage import remove_samples
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        box = QMessageBox(self)
        box.setWindowTitle("Remove from this project?")
        box.setText(f"Remove {len(identifiers)} isolates and their saved analyses from this project?")
        box.setInformativeText(
            "Your original sequence files are never deleted. If you only want them out of the way, "
            "cancel and choose Archive instead — archiving keeps every result and can be undone.")
        box.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(QMessageBox.StandardButton.Cancel)
        check = QCheckBox("Also delete the managed copies WMLSTudio made in its own folders")
        box.setCheckBox(check)
        if box.exec() != QMessageBox.StandardButton.Yes:
            return
        try:
            report = remove_samples(self.project, identifiers,
                                    delete_managed_copy=check.isChecked())
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        message = (f"{len(report['removed'])} isolates removed. Their evidence is recorded in the "
                   "project history and can be restored from Samples ▸ Recently removed.")
        if report["deleted"]:
            message += f" {len(report['deleted'])} managed copies were deleted; your originals were not."
        if report["retained"]:
            message += f" {len(report['retained'])} copies were kept: {report['retained'][0][1]}"
        self.notify(message)

    def remove_sample(self):
        """Menu twin of the right-click removal, so both routes explain themselves."""
        from wmlstudio.context_menus import Selection
        identifiers = sorted(self.selection_ids)
        if not identifiers:
            sample = self.selected_sample()
            identifiers = [sample["id"]] if sample else []
        if not identifiers:
            self.notify("Select the isolates to remove first.")
            return
        self.context_remove(Selection("library", tuple(identifiers)))

    def fill_recently_removed(self, menu):
        """Offer the last removals back, saying plainly which ones cannot be restored."""
        menu.clear()
        entries = [entry for entry in self.project.history()
                   if entry["action"] == "sample_removed"][-20:]
        if not entries:
            menu.addAction("Nothing has been removed from this project").setEnabled(False)
            return
        for entry in reversed(entries):
            details = entry.get("details") or {}
            sample = details.get("sample") or {}
            action = menu.addAction(f"{sample.get('name') or 'Unnamed isolate'} · removed {entry['created_at']}")
            if details.get("format_version") != 2:
                action.setEnabled(False)
                action.setToolTip("This removal predates restorable records.")
                continue
            action.triggered.connect(lambda checked=False, i=entry["id"]: self.restore_removed(i))

    def restore_removed(self, history_id):
        try:
            sample_id = self.project.restore_removed_sample(history_id)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.selection_ids = {sample_id}
        self.refresh_tables()
        self.notify("Isolate restored with every saved analysis it had. Its original file was never "
                    "deleted.")

    def build_menus(self):
        self.command_actions = []

        def action(menu, title, callback, shortcut=None):
            item = QAction(title, self)
            if shortcut:
                item.setShortcut(QKeySequence(shortcut))
            item.triggered.connect(callback)
            menu.addAction(item)
            self.command_actions.append((title, item))
            return item

        bar = self.menuBar()
        file = bar.addMenu("&File")
        action(file, "New project…", self.new_project, "Ctrl+N")
        action(file, "Open project…", self.open_project_dialog)
        action(file, "Save project copy…", self.save_project_copy)
        file.addSeparator()
        imports = file.addMenu("Import")
        action(imports, "Sequence files…", self.browse_files)
        action(imports, "Sequence folder…", self.browse_folder)
        action(imports, "HYDRA report…", self.import_hydra)
        action(imports, "Allelic profile table…", self.import_profile_table)
        action(imports, "Sample bundle…", self.import_sample_bundle)
        exports = file.addMenu("Export")
        for title, fmt in [("Selected samples as CSV…", "csv"), ("Selected samples as TSV…", "tsv"), ("Selected samples as JSON…", "json"), ("HTML report…", "html"), ("PDF report…", "pdf")]:
            action(exports, title, lambda checked=False, f=fmt: self.export_project(f))
        action(exports, "Allelic profile table…", self.export_profile_table)
        action(exports, "Sample bundle…", self.export_sample_bundle)
        file.addSeparator()
        action(file, "Exit", self.close, "Alt+F4")
        samples = bar.addMenu("&Samples")
        action(samples, "Assign organism / workflow…", self.assign_selected)
        action(samples, "Edit annotations…", self.edit_metadata)
        action(samples, "Add to collection…", self.add_collection)
        action(samples, "Highlight cluster…", self.highlight_selected)
        action(samples, "Select visible", self.select_visible_samples)
        action(samples, "Clear selection", self.clear_sample_selection)
        action(samples, "Relink input (same bytes)…", self.relink_selected_input)
        action(samples, "Re-file managed copies now…", self.refile_selected)
        action(samples, "Attach original reads to assemblies…", self.attach_reads_selected)
        samples.addSeparator()
        action(samples, "Archive selected isolates…",
               lambda: self.context_archive(self.library_selection()))
        action(samples, "Restore selected isolates from the archive",
               lambda: self.context_restore(self.library_selection()))
        action(samples, "Remove sample…", self.remove_sample)
        recently_removed = samples.addMenu("Recently removed")
        recently_removed.aboutToShow.connect(lambda menu=recently_removed: self.fill_recently_removed(menu))
        analysis = bar.addMenu("&Analysis")
        action(analysis, "Choose isolates and analyse…", self.choose_and_analyse, "Ctrl+R")
        action(analysis, "Review pending isolates…", lambda: self.choose_and_analyse(pending_only=True))
        action(analysis, "Choose assemblies for HYDRA…", self.choose_hydra_cohort)
        action(analysis, "Characterize identity / virulence / accessory evidence…", self.run_characterization_selected)
        action(analysis, "Assemble and analyse read pairs…", lambda: self.choose_and_analyse(assemble=True))
        action(analysis, "Check missing loci with original reads…", self.open_read_support)
        action(analysis, "Cancel current task", self.cancel_analysis)
        analysis.addSeparator()
        action(analysis, "Build comparison from selected", self.compare_selected)
        action(analysis, "Compute selected comparison scheme…", self.type_comparison_scheme)
        action(analysis, "Create local scheme from selected assemblies…", self.create_adhoc_scheme)
        data = bar.addMenu("&Data")
        action(data, "Online scheme catalog…", self.open_reference_manager)
        action(data, "Research saved library…", self.open_library_research)
        action(data, "AMR databases / updates…", self.open_amr_databases)
        action(data, "Characterization references / updates…", self.install_characterization_references)
        action(data, "Install broader species panel…", self.install_species_panel)
        action(data, "Download practice data…", self.download_practice_cohort)
        action(data, "Import local scheme…", self.import_scheme)
        action(data, "Epidemiology grid / bulk import…", self.open_metadata_grid)
        action(data, "Published cluster threshold guidance…", self.open_threshold_guidance)
        action(data, "Open data folder", self.open_data_folder)
        action(data, "Open selected input folder", self.open_sample_folder)
        view = bar.addMenu("&View")
        for index, title in enumerate(self.nav_names):
            action(view, title, lambda checked=False, i=index: self.navigate(i), f"Alt+{index + 1}")
        action(view, "Refresh current workspace", self.refresh, "F5")
        action(view, "Fit comparison", lambda: self.tree.fit_tree())
        action(view, "Increase interface scale", lambda: self.set_ui_scale(self.ui_scale + 10), "Ctrl++")
        action(view, "Decrease interface scale", lambda: self.set_ui_scale(self.ui_scale - 10), "Ctrl+-")
        action(view, "Reset interface scale", lambda: self.set_ui_scale(100), "Ctrl+0")
        help_menu = bar.addMenu("&Help")
        action(help_menu, "Practice project", self.load_demo)
        action(help_menu, "Problem → solution guide", self.open_workflow_guide, "F1")
        action(help_menu, "Command search…", self.command_palette, "Ctrl+K")
        action(help_menu, "Workflow and limitations", self.open_workflow_guide)

    def command_palette(self):
        names = [title.replace("&", "") for title, _ in self.command_actions]
        title, accepted = QInputDialog.getItem(self, "Command search", "Choose a command (type to search):", names, 0, True)
        if accepted and title in names:
            self.command_actions[names.index(title)][1].trigger()

    def open_workflow_guide(self, topic=None):
        from wmlstudio.workflow_guide import WorkflowGuide
        try:
            if getattr(self, "workflow_guide", None) is None:
                self.workflow_guide = WorkflowGuide(self)
                self.workflow_guide.actionRequested.connect(self.journey_action)
            if isinstance(topic, str) and topic:
                self.workflow_guide.search.setText(topic)
            self.workflow_guide.show()
            self.workflow_guide.raise_()
        except (OSError, ValueError) as error:
            self.error(str(error))

    def attach_reads_selected(self):
        from wmlstudio.read_attachment_dialog import launch_read_attachment
        launch_read_attachment(self)

    def set_ui_scale(self, percent):
        from wmlstudio.interface_settings import interface_preferences, scaled_style
        percent = max(80, min(150, int(percent)))
        self.ui_scale = percent
        QApplication.instance().setStyleSheet(scaled_style(percent))
        interface_preferences(self.root).setValue("scale", percent)
        self.updateGeometry()

    def set_graph_text_scale(self, percent):
        """Resize the lettering on both comparison trees; the trees themselves do not move."""
        from wmlstudio.widgets import set_graph_text_scale
        percent = set_graph_text_scale(percent)
        from wmlstudio.interface_settings import interface_preferences
        interface_preferences(self.root).setValue("graph_scale", percent)
        for name in ("tree", "baseline_tree"):
            view = getattr(self, name, None)
            if view is not None and hasattr(view, "_redraw"):
                view._redraw()

    def open_interface_settings(self):
        from wmlstudio.interface_settings import InterfaceSettingsDialog
        if getattr(self, "interface_dialog", None) is None:
            self.interface_dialog = InterfaceSettingsDialog(self)
            self.interface_dialog.setWindowModality(Qt.WindowModality.WindowModal)
        self.interface_dialog.show()
        self.interface_dialog.raise_()

    def show_sample_detail_for_id(self, sample_id):
        self.show_sample_detail(sample_id)

    def open_isolate_record(self, sample_id):
        from wmlstudio.isolate_dialog import IsolateRecordDialog
        if not isinstance(sample_id, str):
            return
        try:
            self.project.get_sample(sample_id)
        except KeyError:
            self.notify("That isolate is no longer in this project.")
            return
        if not hasattr(self, "isolate_dialogs"):
            self.isolate_dialogs = {}
        if sample_id not in self.isolate_dialogs:
            self.isolate_dialogs[sample_id] = IsolateRecordDialog(self, sample_id)
        dialog = self.isolate_dialogs[sample_id]
        dialog.refresh_record()
        dialog.show()
        dialog.raise_()

    def open_metadata_grid(self, checked=False):
        from wmlstudio.metadata_grid import launch_metadata_grid
        launch_metadata_grid(self)

    def open_read_support(self, sample_id=None):
        from wmlstudio.read_support_dialog import launch_read_support
        launch_read_support(self, sample_id=sample_id if isinstance(sample_id, str) else None)

    def open_threshold_guidance(self):
        from wmlstudio.threshold_dialog import ThresholdGuideDialog
        dialog = ThresholdGuideDialog(self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.apply_threshold_guidance(dialog.evidence)

    def choose_feature_cohort(self):
        from wmlstudio.cohort_picker import CohortPickerDialog
        dialog = CohortPickerDialog(self.project.samples(), self.project, "Which isolates should the evidence workspace show?", self.feature_ids, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.feature_ids = set(dialog.selected_ids)
        self.refresh_features()
        self.navigate(4)

    def library_selection(self):
        """The isolate-library selection as a right-click Selection, for menu twins."""
        from wmlstudio.context_menus import Selection
        return Selection("library", tuple(sorted(self.selection_ids)))

    def choose_and_analyse(self, checked=False, assemble=False, pending_only=False):
        if self.busy():
            return
        from wmlstudio.cohort_picker import CohortPickerDialog
        samples = self.project.samples()
        initial = {sample["id"] for sample in samples if sample["status"] in {"queued", "interrupted"}} if pending_only else None
        dialog = CohortPickerDialog(samples, self.project, "Which isolates should be analysed?", initial, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.start_analysis(confirm=True, all_samples=True, sample_ids=dialog.selected_ids, assemble=assemble)

    def choose_hydra_cohort(self):
        if self.busy():
            return
        from wmlstudio.cohort_picker import CohortPickerDialog
        samples = self.project.samples()
        dialog = CohortPickerDialog(samples, self.project, "Which assemblies should HYDRA analyse?", parent=self, include_reads=False)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        chosen = [sample for sample in samples if sample["id"] in dialog.selected_ids]
        plan = self.review_run_plan(chosen, hydra=True)
        if plan and plan.get("hydra"):
            self._run_cancelled = False
            self.run_hydra_plan(dialog.selected_ids, plan, require_completed=False)

    def build_schemes(self):
        super().build_schemes()
        layout = self.pages.widget(3).widget().layout()
        row = FlowLayout()
        row.addWidget(button("Browse online / install updates…", self.open_reference_manager, True))
        row.addWidget(button("Refresh installed schemes", self.populate_schemes))
        row.addWidget(button("Create local scheme…", self.create_adhoc_scheme))
        layout.insertLayout(2, row)
        layout.insertWidget(3, label("Updates create new snapshots; existing results retain their original fingerprint.", "small", True))

    def open_reference_manager(self):
        from wmlstudio.reference_dialog import ReferenceManagerDialog
        if self.reference_dialog is None:
            self.reference_dialog = ReferenceManagerDialog(self.root, self)
            self.reference_dialog.schemeInstalled.connect(self.reference_installed)
        self.reference_dialog.show()
        self.reference_dialog.raise_()
        self.reference_dialog.activateWindow()

    def reference_installed(self, path):
        self.populate_schemes()
        self.refresh_cohort_table()
        self.notify(f"Reference snapshot installed: {Path(path).name}. Existing results retain their original scheme fingerprint.")

    def create_adhoc_scheme(self):
        if self.busy():
            return
        samples = [sample for sample in self.selected_samples() if sample.get("input_path")]
        if not samples:
            self.notify("Select the assembly cohort in Samples first. Profile-only imports cannot define new sequence loci.")
            return
        from wmlstudio.adhoc import create_adhoc_scheme
        from wmlstudio.adhoc_dialog import AdhocSchemeDialog
        dialog = AdhocSchemeDialog(samples, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        plan = dict(dialog.plan)
        self.launch_task(lambda cancelled, progress: create_adhoc_scheme(library_root=self.root / "schemes",
                         cancelled=cancelled, progress=progress, **plan), "adhoc_scheme",
                         lambda result: self.reference_installed(result["path"]))

    def closeEvent(self, event):
        if not self.cancel_comparison_for_close():
            event.ignore()
            return
        if self.worker and self.worker.isRunning():
            return super().closeEvent(event)
        if self.amr_database_dialog is not None:
            database_worker = self.amr_database_dialog.worker
            if database_worker and database_worker.isRunning():
                database_worker.cancel()
                event.ignore()
                QTimer.singleShot(100, self.close)
                return
            self.amr_database_dialog.close()
        if self.reference_dialog is not None:
            reference_worker = self.reference_dialog.worker
            if reference_worker and reference_worker.isRunning():
                reference_worker.cancel()
                event.ignore()
                QTimer.singleShot(100, self.close)
                return
            self.reference_dialog.close()
        if self.library is not None:
            self.library.index_project(self.project)
            self.library.close()
            self.library = None
        super().closeEvent(event)
