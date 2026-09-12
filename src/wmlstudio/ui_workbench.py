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
)

from wmlstudio.background import FunctionWorker
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    flattened_metadata,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.widgets import DropZone, Helix, Metric, button, card, label


class WorkbenchMixin:
    def __init__(self, *args, **kwargs):
        self.selection_ids = set()
        self.cohort_ids = None
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
        super().__init__(*args, **kwargs)
        from wmlstudio import theme
        if hasattr(theme, "apply_dark_palette"):
            theme.apply_dark_palette(QApplication.instance())
        QApplication.instance().setStyleSheet(theme.STYLE)
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
        self.statusBar().showMessage("Local genomics workbench · Select samples to begin")

    def build_overview(self):
        _, layout = self.page()
        hero, content = card()
        row = QHBoxLayout()
        words = QVBoxLayout()
        words.addWidget(label("YOUR GENOMICS LIBRARY", "eyebrow"))
        words.addWidget(label("Every isolate. One connected workspace.", "title", True))
        words.addWidget(label("Browse by organism and ST, keep your input copies safe, and return to any saved analysis.", "muted", True))
        actions = FlowLayout()
        actions.addWidget(button("Import sequences…", self.browse_files, True))
        actions.addWidget(button("New project…", self.new_project))
        actions.addWidget(button("Open data folder", self.open_data_folder))
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
        metrics = QHBoxLayout()
        self.metrics = [Metric("Isolates", "in the active library"), Metric("Analysed", "saved evidence"),
                        Metric("Exact MLST", "registered profiles"), Metric("Review", "incomplete or failed")]
        for metric in self.metrics:
            metrics.addWidget(metric)
        layout.addLayout(metrics)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        explorer, content = card()
        content.addWidget(label("Library navigator", "cardTitle"))
        self.library_tree = QTreeWidget()
        self.library_tree.setHeaderHidden(True)
        self.library_tree.setMinimumWidth(250)
        self.library_tree.itemActivated.connect(self.open_library_group)
        self.library_tree.itemClicked.connect(self.open_library_group)
        content.addWidget(self.library_tree)
        splitter.addWidget(explorer)
        recent, content = card()
        content.addWidget(label("Recent samples · double-click to inspect", "cardTitle"))
        self.recent_table = make_table(["Sample", "Input", "Status", "ST", "Loci"])
        self.recent_table.cellDoubleClicked.connect(self.open_recent_sample)
        content.addWidget(self.recent_table)
        self.drop_zone = DropZone()
        self.drop_zone.filesDropped.connect(lambda paths: self.import_paths(paths, configure=True))
        self.drop_zone.browseRequested.connect(self.browse_files)
        content.addWidget(self.drop_zone)
        splitter.addWidget(recent)
        splitter.setSizes([320, 650])
        layout.addWidget(splitter, 1)
        self.practice_notice = label("Practice project · synthetic sequences, not biological isolates.", "small")
        self.practice_notice.hide()
        layout.addWidget(self.practice_notice)

    def build_samples(self):
        _, layout = self.page()
        self.heading(layout, "Sample library", "Select isolates, assign their workflow, and inspect the connected evidence. Ctrl/Shift selects a cohort.")
        row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search names, organisms, STs, AMR genes, collections or metadata…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_tables)
        row.addWidget(self.search, 1)
        row.addWidget(button("Import…", self.browse_files, True))
        row.addWidget(button("Assign organism…", self.assign_selected))
        row.addWidget(button("Metadata…", self.edit_metadata))
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
        filters.addWidget(button("Select visible", self.select_visible_samples))
        filters.addWidget(button("Clear selection", self.clear_sample_selection))
        layout.addLayout(filters)
        run, content = card()
        strip = FlowLayout()
        self.scheme_combo = QComboBox()
        self.scheme_combo.setMinimumWidth(180)
        self.scheme_combo.setAccessibleName("Typing scheme override")
        strip.addWidget(self.scheme_combo, 1)
        self.run_button = button("Analyse pending…", lambda: self.start_analysis(confirm=True), True)
        self.rerun_button = button("Analyse selected…", lambda: self.start_analysis(selected_only=True, confirm=True))
        strip.addWidget(self.run_button)
        strip.addWidget(self.rerun_button)
        strip.addWidget(button("Compare selected", self.compare_selected))
        strip.addWidget(button("Report selected", self.report_selected))
        content.addLayout(strip)
        self.selection_label = label("No samples selected · each sample keeps its own organism and typing workflow", "small", True)
        content.addWidget(self.selection_label)
        layout.addWidget(run)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.sample_table = make_table(["Sample", "Input", "Status", "ST", "Called loci", "Genus", "Species",
                                        "Evidence", "Scheme", "AMR genes", "Collection", "Storage"])
        self.sample_table.itemSelectionChanged.connect(self.sample_selection_changed)
        self.sample_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.sample_table.customContextMenuRequested.connect(self.sample_context_menu)
        splitter.addWidget(self.sample_table)
        self.detail = QTextBrowser()
        self.detail.setOpenExternalLinks(False)
        self.detail.setMinimumHeight(120)
        inspector = QTabWidget()
        inspector.addTab(self.detail, "Sample evidence")
        self.history_view = QTextBrowser()
        inspector.addTab(self.history_view, "History / provenance")
        splitter.addWidget(inspector)
        splitter.setSizes([420, 190])
        layout.addWidget(splitter, 1)

    def refresh(self):
        if self._refreshing:
            return
        self._refreshing = True
        try:
            super().refresh()
            isolates = [sample for sample in self.current_samples if sample.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
            mates = len(self.current_samples) - len(isolates)
            self.metrics[0].value.setText(str(len(isolates)))
            self.metrics[0].hint.setText(f"{len(isolates)} isolates · {mates} linked mates" if mates else "in the active library")
            if self.library is not None and not (self.worker and self.worker.isRunning()):
                self.library.index_project(self.project)
            valid = {s["id"] for s in self.current_samples}
            self.selection_ids.intersection_update(valid)
            self.report_ids.intersection_update(valid)
            self.refresh_library()
            if hasattr(self, "cohort_table"):
                self.refresh_cohort_table()
            if hasattr(self, "report_table"):
                self.refresh_report_table()
            if hasattr(self, "feature_table"):
                self.refresh_features()
        finally:
            self._refreshing = False

    def refresh_library(self):
        tree = self.library_tree
        tree.blockSignals(True)
        tree.clear()
        all_item = QTreeWidgetItem([f"All samples  ({len(self.current_samples)})"])
        all_item.setData(0, Qt.ItemDataRole.UserRole, None)
        tree.addTopLevelItem(all_item)
        groups = {}
        collections = {}
        for sample in self.current_samples:
            if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate":
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
        mates = [sample for sample in self.current_samples if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate"]
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
            self.library_filter = ("ids", [item.data(Qt.ItemDataRole.UserRole)])
            self.selection_ids = set(self.library_filter[1])
            self.refresh_tables()
            self.navigate(1)

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
        taxa = [organism_for(s) for s in self.current_samples]
        self.fill_filter(self.genus_filter, [t[0] for t in taxa], "All genera")
        self.fill_filter(self.species_filter, [t[1] for t in taxa], "All species")
        self.fill_filter(self.collection_filter, [c for s in self.current_samples for c in
                         (s.get("metadata") or {}).get("collections", [])], "All collections")
        query = self.search.text().strip().casefold()
        visible = []
        for sample in self.current_samples:
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
        self.fill_sample_table(self.recent_table, self.current_samples[-12:])
        self.selection_label.setText(f"{len(self.selection_ids)} selected · {len(visible)} visible · {len(self.current_samples)} in project")
        self.show_sample_detail()

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
            storage = "Managed copy" if metadata.get("workflow", {}).get("managed") else "Linked original"
            if sample.get("missing_input") and kind != "profile":
                storage = "Input unavailable"
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
        self.selection_label.setText(f"{len(self.selection_ids)} selected · {len(visible)} visible · {len(self.current_samples)} in project")
        self.show_sample_detail()

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

    def show_sample_detail(self):
        sample = self.selected_sample()
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
            candidates = path.rglob("*") if path.is_dir() else [path]
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
        self.import_assignments(dialog.assignments, dialog.options)

    def import_assignments(self, assignments, options):
        from wmlstudio.storage import import_samples
        root = options.get("storage_root") or str(self.project_path.with_suffix(".files"))
        def operation(cancelled, progress):
            return import_samples(
                self.project, assignments, storage_root=root, managed=options.get("managed", True),
                append_st=options.get("append_st", False), cancelled=cancelled, progress=progress)
        self.launch_task(operation, "import", lambda ids: self.import_completed(ids))

    def import_completed(self, ids):
        self.selection_ids = set(ids)
        self.clear_filters()
        self.refresh()
        self.navigate(1)
        self.notify(f"Imported {len(ids)} samples. Review their assignments, then analyse when ready.")

    def assign_selected(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if not samples:
            self.notify("Select one or more sample rows first.")
            return
        from wmlstudio.storage import assign_organism
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        dialog = BatchAssignmentDialog(samples, self.scheme_entries(), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        for assignment in dialog.assignments:
            assign_organism(self.project, [assignment["sample_id"]], assignment["genus"], assignment["species"],
                            scheme_path=assignment.get("scheme_path"), typing_mode=assignment["typing_mode"])
        self.refresh()
        self.notify("Assignments saved. Previous analysis evidence remains in the sample history.")

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

    def start_analysis(self, checked=False, all_samples=False, selected_only=False, confirm=False, assemble=False):
        if self.busy():
            return
        samples = [s for s in self.project.samples() if (s["id"] in self.selection_ids if selected_only
                   else all_samples or s["status"] in {"queued", "failed", "interrupted"})]
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
            self.assemble_pairs(pairing.assignments, plan)
            return
        self.begin_typing(samples, scheme)

    def begin_typing(self, samples, scheme):
        if scheme:
            samples = [dict(sample, metadata={**sample.get("metadata", {}), "workflow": {
                **sample.get("metadata", {}).get("workflow", {}), "typing_mode": "manual", "scheme_path": scheme}})
                for sample in samples]
        self.project.set_setting("scheme_path", scheme or "")
        self.worker_role = "analysis"
        self.worker = AnalysisWorker(samples, scheme, parent=self, installed_scheme_paths=self.scheme_paths)
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
        project = self.project
        def operation(cancelled, progress):
            identifiers = []
            for index, pair in enumerate(pairs):
                primary = project.get_sample(pair["primary_id"])
                mate = project.get_sample(pair["mate_id"])
                destination = self.project_path.with_suffix(".files") / "assemblies" / primary["id"] / uuid.uuid4().hex
                result = run_skesa(primary["input_path"], mate["input_path"], destination,
                    threads=plan.get("threads", 4), memory_gb=plan.get("memory_gb", 8), cancelled=cancelled,
                    progress=lambda done, total, message: progress(index * 100 + done / max(1, total) * 100, len(pairs) * 100, f"{primary['name']} · {message}"))
                associate_assembly(project, primary["id"], mate["id"], result)
                identifiers.append(primary["id"])
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
                from wmlstudio.storage import organize_sample

                def organize(cancelled, progress):
                    for i, sample in enumerate(managed):
                        organize_sample(self.project, sample["id"], cancelled=cancelled)
                        progress(i + 1, len(managed), "Organising managed input copies by organism and ST…")
                    return len(managed)

                self.launch_task(organize, "organize", lambda count: self.notify(f"{count} managed sample copies organised. Originals unchanged."))
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
        samples = [sample for sample in self.project.samples() if sample["id"] in identifiers
                   and sample.get("input_path") and (not require_completed or sample["status"] == "completed")]
        assemblies = [sample for sample in samples if (sample.get("result") or {}).get("kind", sample.get("kind")) == "fasta"
                      or Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fa", ".fasta", ".fna"))]
        if not assemblies:
            self.notify("No FASTA assemblies are ready for HYDRA in this selection. Raw reads must be assembled first.")
            return

        def operation(cancelled, progress):
            import hashlib

            from wmlstudio.sequence import AnalysisCancelled
            reports = []
            for index, sample in enumerate(assemblies):
                if cancelled():
                    raise AnalysisCancelled("HYDRA cancelled; previously completed sample evidence is saved.")
                genus, species, _ = organism_for(sample)
                assigned = sample.get("metadata", {}).get("organism", {})
                # A provisional MLST lineage alone is not a verified mutation-catalog assignment.
                organism = " ".join([genus, species]) if assigned.get("genus") and assigned.get("species") else None
                progress(index, len(assemblies), f"HYDRA · {sample['name']}")
                report = run_assemblies([sample["input_path"]], plan["db_root"], plan.get("databases"),
                    sample_names=[sample["id"]], organism=organism, threads=plan.get("threads", 4),
                    protein=plan.get("protein", True), cancelled=cancelled,
                    point_mutations=plan.get("point_mutations", True),
                    **plan.get("thresholds", {}),
                    progress=lambda done, total, message: progress(index, len(assemblies), f"{sample['name']} · {message}"))
                link_hydra(self.project, report, {sample["id"]: sample["id"]})
                reports.append(report)
                combined = dict(report, samples=[entry for result in reports for entry in result["samples"]])
                components = [result["import_provenance"] for result in reports]
                combined["import_provenance"] = {"source_path": "Native HYDRA run", "components": components,
                    "sha256": hashlib.sha256(json.dumps(components, sort_keys=True).encode()).hexdigest()}
                combined["execution_provenance"] = {"sample_runs": [result.get("execution_provenance", {}) for result in reports]}
                self.project.set_setting("hydra_report", combined)
                progress(index + 1, len(assemblies), f"HYDRA evidence saved · {sample['name']}")
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
            self.cohort_ids = set(cohort) if cohort is not None else None
            self.library_filter = None
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

    def sample_context_menu(self, point):
        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        for title, callback in [("Assign organism / workflow…", self.assign_selected), ("Edit annotations…", self.edit_metadata),
                                ("Add to collection…", self.add_collection), ("Compare selected", self.compare_selected),
                                ("Report selected", self.report_selected), ("Highlight cluster…", self.highlight_selected),
                                ("Relink input (same bytes)…", self.relink_selected_input),
                                ("Open input folder", self.open_sample_folder), ("Remove from project…", self.remove_sample)]:
            menu.addAction(title, callback)
        menu.exec(self.sample_table.viewport().mapToGlobal(point))

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
        action(samples, "Remove sample…", self.remove_sample)
        analysis = bar.addMenu("&Analysis")
        action(analysis, "Analyse selected…", lambda: self.start_analysis(selected_only=True, confirm=True), "Ctrl+R")
        action(analysis, "Analyse pending…", lambda: self.start_analysis(confirm=True))
        action(analysis, "Run HYDRA on selected assemblies…", self.run_hydra_selected)
        action(analysis, "Assemble and analyse selected read pairs…", lambda: self.start_analysis(selected_only=True, confirm=True, assemble=True))
        action(analysis, "Cancel current task", self.cancel_analysis)
        analysis.addSeparator()
        action(analysis, "Build comparison from selected", self.compare_selected)
        action(analysis, "Compute selected comparison scheme…", self.type_comparison_scheme)
        action(analysis, "Create local scheme from selected assemblies…", self.create_adhoc_scheme)
        data = bar.addMenu("&Data")
        action(data, "Online scheme catalog…", self.open_reference_manager)
        action(data, "Research saved library…", self.open_library_research)
        action(data, "AMR databases / updates…", self.open_amr_databases)
        action(data, "Import local scheme…", self.import_scheme)
        action(data, "Open data folder", self.open_data_folder)
        action(data, "Open selected input folder", self.open_sample_folder)
        view = bar.addMenu("&View")
        for index, title in enumerate(self.nav_names):
            action(view, title, lambda checked=False, i=index: self.navigate(i), f"Alt+{index + 1}")
        action(view, "Refresh current workspace", self.refresh, "F5")
        action(view, "Fit comparison", lambda: self.tree.fit_tree())
        help_menu = bar.addMenu("&Help")
        action(help_menu, "Practice project", self.load_demo)
        action(help_menu, "Command search…", self.command_palette, "Ctrl+K")
        action(help_menu, "Workflow and limitations", lambda: self.navigate(6))

    def command_palette(self):
        names = [title.replace("&", "") for title, _ in self.command_actions]
        title, accepted = QInputDialog.getItem(self, "Command search", "Choose a command (type to search):", names, 0, True)
        if accepted and title in names:
            self.command_actions[names.index(title)][1].trigger()

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
