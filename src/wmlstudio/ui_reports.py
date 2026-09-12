"""Selection-scoped reports, linked AMR matrices and interoperable sample exchange."""


from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QInputDialog,
    QLineEdit,
    QTableView,
    QTabWidget,
    QVBoxLayout,
)

from wmlstudio.export import export_results
from wmlstudio.sample_workflow import (
    hydra_evidence_status,
    link_hydra,
    set_cluster,
    suggest_hydra_links,
)
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    flattened_metadata,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.ui_compare import EvidenceMatrixModel
from wmlstudio.widgets import button, card, label


class ReportWorkspaceMixin:
    def build_reports(self):
        _, layout = self.page()
        self.heading(layout, "Report a defined cohort", "Choose exactly which samples belong in the report. Highlighted isolates retain their group label, color and linked AMR evidence.")
        controls = FlowLayout()
        self.report_scope = QComboBox()
        self.report_scope.addItems(["Explicit report selection", "Current comparison cohort", "All project samples"])
        self.report_scope.setCurrentIndex(2)
        self.report_scope.currentIndexChanged.connect(self.refresh_report_table)
        controls.addWidget(self.report_scope, 1)
        controls.addWidget(button("Use Samples selection", self.report_selected))
        controls.addWidget(button("Select all", self.select_all_reports))
        controls.addWidget(button("Select none", self.clear_report_selection))
        controls.addWidget(button("Highlight selected…", self.highlight_selected))
        layout.addLayout(controls)
        self.report_count = label("Report scope: all project samples", "cardTitle")
        layout.addWidget(self.report_count)
        self.report_table = make_table(["Include", "Highlight", "Sample", "Organism", "ST", "AMR genes", "Group"])
        self.report_table.setColumnWidth(0, 65)
        self.report_table.setColumnWidth(1, 75)
        self.report_table.itemChanged.connect(self.report_item_changed)
        layout.addWidget(self.report_table, 1)
        self._report_filling = False
        panel, content = card()
        row = FlowLayout()
        for title, fmt in [("PDF report…", "pdf"), ("HTML report…", "html"), ("CSV…", "csv"), ("TSV…", "tsv"), ("JSON…", "json")]:
            row.addWidget(button(title, lambda checked=False, f=fmt: self.export_project(f), fmt == "pdf"))
        content.addLayout(row)
        row = FlowLayout()
        row.addWidget(button("Export profile table…", self.export_profile_table))
        row.addWidget(button("Export portable bundle…", self.export_sample_bundle))
        row.addWidget(button("Save project copy…", self.save_project_copy))
        content.addLayout(row)
        content.addWidget(label("Reports describe sequence evidence, not measured susceptibility or proof of transmission. Profile bundles carry results and metadata without sequence files.", "small", True))
        layout.addWidget(panel)

    def report_sample_ids(self):
        index = self.report_scope.currentIndex()
        if index == 2:
            return {s["id"] for s in self.project.samples() if s.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"}
        if index == 1:
            return {r["sample_id"] for r in self._last_comparison}
        return set(self.report_ids)

    def refresh_report_table(self):
        if not hasattr(self, "report_table") or getattr(self, "_report_filling", False):
            return
        self._report_filling = True
        try:
            chosen = self.report_sample_ids()
            samples = self.project.samples()
            self.report_table.blockSignals(True)
            self.report_table.setSortingEnabled(False)
            self.report_table.setRowCount(len(samples))
            for row, sample in enumerate(samples):
                genus, species, _ = organism_for(sample)
                group = sample.get("metadata", {}).get("cluster", {})
                values = ["", "", sample["name"], " ".join([genus, species]).strip() or "Unknown",
                          (sample.get("result") or {}).get("st") or "—", "; ".join(gene_names(sample)) or "—", group.get("label", "")]
                for column, value in enumerate(values):
                    item = cell(value, sample["id"])
                    if column < 2:
                        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                        checked = sample["id"] in chosen if column == 0 else bool(group.get("highlight"))
                        item.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
                    if group.get("highlight"):
                        item.setBackground(QColor("#343348"))
                    self.report_table.setItem(row, column, item)
            self.report_table.setSortingEnabled(True)
            self.report_table.blockSignals(False)
            highlights = sum(s["id"] in chosen and s.get("metadata", {}).get("cluster", {}).get("highlight", False) for s in samples)
            self.report_count.setText(f"{len(chosen)} included · {highlights} highlighted · {len(samples)} in project")
        finally:
            self._report_filling = False

    def report_item_changed(self, item):
        if self._report_filling or item.column() > 1:
            return
        sid = item.data(Qt.ItemDataRole.UserRole)
        checked = item.checkState() == Qt.CheckState.Checked
        if item.column() == 0:
            chosen = self.report_sample_ids()
            chosen.add(sid) if checked else chosen.discard(sid)
            self.report_ids = chosen
            self.report_scope.blockSignals(True)
            self.report_scope.setCurrentIndex(0)
            self.report_scope.blockSignals(False)
        else:
            sample = self.project.get_sample(sid)
            group = sample.get("metadata", {}).get("cluster", {})
            set_cluster(self.project, [sid], group.get("label") or "Highlighted", group.get("color") or "#E9AD66", checked)
        self.refresh_report_table()

    def report_selected(self):
        if not self.selection_ids:
            self.notify("Select sample rows first, or choose an explicit scope on Reports.")
            self.navigate(5)
            return
        self.report_ids = set(self.selection_ids)
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()
        self.navigate(5)

    def select_all_reports(self):
        self.report_ids = {s["id"] for s in self.project.samples()}
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()

    def clear_report_selection(self):
        self.report_ids.clear()
        self.report_scope.setCurrentIndex(0)
        self.refresh_report_table()

    def highlight_selected(self):
        ids = self.selection_ids or self.report_ids
        if not ids:
            self.notify("Select isolates in Samples or the graph before assigning a highlight.")
            return
        name, accepted = QInputDialog.getText(self, "Highlight sample group", "Group label (for example, Investigation A):")
        if not accepted or not name.strip():
            return
        color = QColorDialog.getColor(QColor("#E9AD66"), self, "Group highlight colour")
        if color.isValid():
            set_cluster(self.project, ids, name.strip(), color.name(), True)
            self.refresh()
            self.refresh_comparison()

    def export_project(self, fmt):
        ids = self.report_sample_ids()
        if not ids:
            self.notify("The report selection is empty. Select samples or choose All project samples.")
            return
        path, _ = QFileDialog.getSaveFileName(self, f"Export {len(ids)} samples", f"wmlstudio-report.{fmt}", f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        try:
            self.check_output(path)
            if fmt == "pdf":
                self.write_pdf_report(path, selected_ids=ids)
            else:
                export_results(self.report_records(), path, fmt, selected_ids=ids)
            scope = f"all {len(ids)} project samples" if self.report_scope.currentIndex() == 2 else f"{len(ids)} selected samples"
            self.notify(f"Exported {scope} with their linked metadata and highlights.")
        except Exception as exc:
            self.error(exc)

    def write_pdf_report(self, path, selected_ids=None):
        ids = self.report_sample_ids() if selected_ids is None else set(selected_ids)
        samples = [s for s in self.report_records() if s["id"] in ids]
        if not samples:
            raise ValueError("No samples selected for the PDF report.")
        return super().write_pdf_report(path, samples=samples)

    def report_records(self):
        return [dict(sample, analyses=self.available_profiles(sample)) for sample in self.project.samples()]

    def build_hydra(self):
        super().build_hydra()
        layout = self.pages.widget(4).widget().layout()
        actions = FlowLayout()
        actions.addWidget(button("Run HYDRA on selected…", self.run_hydra_selected, True))
        actions.addWidget(button("AMR databases…", self.open_amr_databases))
        actions.addWidget(button("Link report to project samples…", self.link_hydra_samples, True))
        actions.addWidget(button("Refresh linked features", self.refresh_features))
        actions.addWidget(label("Stable sample IDs join AMR, typing and annotations; unresolved names are not merged.", "small", True), 1)
        layout.insertLayout(3, actions)
        old_splitter = layout.takeAt(layout.count() - 1).widget()
        tabs = QTabWidget()
        tabs.addTab(old_splitter, "Source report evidence")
        self.feature_table = QTableView()
        self.feature_model = EvidenceMatrixModel(self.feature_table)
        self.feature_table.setModel(self.feature_model)
        self.feature_table.setSortingEnabled(True)
        self.feature_table.setAlternatingRowColors(True)
        tabs.addTab(self.feature_table, "All sample features")
        from PySide6.QtWidgets import QWidget
        matrix_page = QWidget()
        matrix_layout = QVBoxLayout(matrix_page)
        self.gene_filter = QLineEdit()
        self.gene_filter.setPlaceholderText("Filter AMR gene columns (e.g. vanA, blaKPC)…")
        self.gene_filter.textChanged.connect(self.refresh_features)
        matrix_layout.addWidget(self.gene_filter)
        self.amr_matrix = QTableView()
        self.amr_model = EvidenceMatrixModel(self.amr_matrix)
        self.amr_matrix.setModel(self.amr_model)
        self.amr_matrix.setAlternatingRowColors(True)
        matrix_layout.addWidget(self.amr_matrix)
        matrix_layout.addWidget(label("Present = detected in the linked report; not detected is not a susceptibility claim. No report = unknown, never absence.", "small", True))
        tabs.addTab(matrix_page, "AMR gene matrix")
        layout.addWidget(tabs, 1)

    def restore_hydra(self):
        super().restore_hydra()
        if hasattr(self, "feature_table"):
            self.refresh_features()

    def refresh_features(self):
        if not hasattr(self, "feature_model"):
            return
        rows = []
        samples = self.project.samples()
        for sample in samples:
            genus, species, status = organism_for(sample)
            result = sample.get("result") or {}
            amr_state = hydra_evidence_status(sample)
            row = {"Sample": sample["name"], "Genus": genus or "Unknown", "Species": species or "—", "Organism evidence": status,
                   "ST": result.get("st"), "Scheme": result.get("scheme"), "AMR genes": "; ".join(gene_names(sample)),
                   "HYDRA linked": bool(sample.get("metadata", {}).get("hydra")), "AMR evidence state": amr_state["status"],
                   "AMR evidence note": amr_state["reason"], **flattened_metadata(sample)}
            row.update({"QC." + key: value for key, value in result.get("qc", {}).items() if not isinstance(value, (dict, list))})
            rows.append(row)
        main = ["Sample", "Genus", "Species", "Organism evidence", "ST", "Scheme", "AMR genes", "HYDRA linked", "AMR evidence state", "AMR evidence note"]
        self.feature_model.replace(main + sorted({key for row in rows for key in row} - set(main)), rows)
        query = self.gene_filter.text().casefold() if hasattr(self, "gene_filter") else ""
        genes = sorted({gene for sample in samples for gene in gene_names(sample) if query in gene.casefold()})
        matrix = []
        for sample in samples:
            state = hydra_evidence_status(sample)["status"]
            present = set(gene_names(sample))
            def state_cell(gene):
                if state == "stale":
                    return "Stale evidence — rerun / review"
                if state == "missing":
                    return "No report"
                value = "Present" if gene in present else "Not detected"
                return value
            matrix.append({"Sample": sample["name"], "Evidence state": state, **{gene: state_cell(gene) for gene in genes}})
        self.amr_model.replace(["Sample", "Evidence state", *genes], matrix)

    def import_hydra(self):
        before = self.project.get_setting("hydra_report", None)
        super().import_hydra()
        after = self.project.get_setting("hydra_report", None)
        if after and after != before and self.project.samples():
            self.link_hydra_samples()

    def link_hydra_samples(self):
        report = self.project.get_setting("hydra_report", None)
        if not report:
            self.notify("Import a HYDRA JSON report first.")
            return
        samples = self.project.samples()
        if not samples:
            self.notify("Import the associated sequence files or a sample profile bundle before linking report evidence.")
            return
        proposed = suggest_hydra_links(self.project, report)
        dialog = QDialog(self)
        dialog.setWindowTitle("Map HYDRA evidence to project samples")
        dialog.resize(780, 480)
        layout = QVBoxLayout(dialog)
        layout.addWidget(label("Review every mapping. Duplicate or unmatched names stay unlinked until you choose a sample. Local MLST calls are not overwritten.", "muted", True))
        table = make_table(["HYDRA sample", "Project sample"])
        table.setSortingEnabled(False)
        table.setRowCount(len(report["samples"]))
        mappings = []
        for row, source in enumerate(report["samples"]):
            name = source["sample"]
            table.setItem(row, 0, cell(name))
            combo = QComboBox()
            combo.addItem("Leave unlinked", None)
            for sample in samples:
                combo.addItem(f"{sample['name']} · {sample['id'][:8]}", sample["id"])
            combo.setCurrentIndex(max(0, combo.findData(proposed["mapping"].get(name))))
            table.setCellWidget(row, 1, combo)
            mappings.append((name, combo))
        layout.addWidget(table)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        controls.accepted.connect(dialog.accept)
        controls.rejected.connect(dialog.reject)
        layout.addWidget(controls)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            try:
                linked = link_hydra(self.project, report, {name: combo.currentData() for name, combo in mappings if combo.currentData()})
                self.refresh()
                self.notify(f"Linked HYDRA evidence to {len(linked)} project samples. AMR metadata is now available in tables and reports.")
            except Exception as exc:
                self.error(exc)

    def import_profile_table(self):
        if self.busy():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import allelic profile table", "", "TSV / CSV (*.tsv *.tab *.txt *.csv)")
        if not path:
            return
        name, accepted = QInputDialog.getText(self, "Profile scheme identity", "Scheme name for this imported table. Without a verified fingerprint, it remains a distinct external snapshot:")
        if not accepted or not name.strip():
            return
        from wmlstudio.library import import_profile_table
        self.launch_task(lambda cancelled, progress: import_profile_table(self.project, path, name.strip()),
                         "profile_import", lambda result: self.import_completed(result["sample_ids"]))

    def export_profile_table(self):
        results = self.comparison_results()
        if not results:
            self.notify("Build a comparison with one scheme snapshot before exporting its profile table.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export full allele matrix", "profiles.tsv", "TSV (*.tsv)")
        if path:
            try:
                self.check_output(path)
                from wmlstudio.library import export_profile_table
                export_profile_table(self.project, path, [r["sample_id"] for r in results], results[0]["scheme_digest"])
                self.notify(f"Exported all loci for {len(results)} profiles.")
            except Exception as exc:
                self.error(exc)

    def export_sample_bundle(self):
        ids = self.report_sample_ids()
        if not ids:
            self.notify("Choose a non-empty report selection before exporting a bundle.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export portable sample bundle", "sample-bundle.json", "JSON (*.json)")
        if path:
            try:
                self.check_output(path)
                from wmlstudio.library import export_bundle
                export_bundle(self.project, path, ids)
                self.notify("Bundle exported: profiles, metadata and collections; no sequence files.")
            except Exception as exc:
                self.error(exc)

    def import_sample_bundle(self):
        if self.busy():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import sample bundle", "", "JSON (*.json)")
        if not path:
            return
        from wmlstudio.library import import_bundle
        self.launch_task(lambda cancelled, progress: import_bundle(self.project, path), "bundle_import",
                         lambda result: self.import_completed(result["sample_ids"]))
