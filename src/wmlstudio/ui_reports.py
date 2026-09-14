"""Selection-scoped reports, linked AMR matrices and interoperable sample exchange."""


import os
import tempfile
from pathlib import Path

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt
from PySide6.QtGui import QColor, QImage, QPageSize, QPainter, QPdfWriter, QTextDocument
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QInputDialog,
    QLineEdit,
    QMenu,
    QTableView,
    QTabWidget,
    QTextBrowser,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.export import REPORT_PRESETS, export_results, review_report_html, write_review_report
from wmlstudio.sample_workflow import (
    hydra_evidence_status,
    link_hydra,
    set_cluster,
    suggest_hydra_links,
)
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.ui_compare import EvidenceMatrixModel
from wmlstudio.widgets import button, card, label


class ReportWorkspaceMixin:
    def build_reports(self):
        _, layout = self.page()
        self.heading(layout, "Prepare a review report", "Choose a report template and its own isolates. Clear evidence summaries come first; technical provenance is optional.")
        controls = FlowLayout()
        self.report_scope = QComboBox(self)
        self.report_scope.hide()  # Kept as the explicit scope model, not a global selector.
        self.report_scope.addItems(["Explicit report selection", "Current comparison cohort", "All project samples"])
        self.report_scope.setCurrentIndex(0)
        self.report_scope.currentIndexChanged.connect(self.refresh_report_table)
        self.report_preset = QComboBox()
        for key, value in REPORT_PRESETS.items():
            self.report_preset.addItem(value['title'], key)
        self.report_preset.setCurrentIndex(1)
        self.report_preset.currentIndexChanged.connect(self.report_preset_changed)
        controls.addWidget(self.report_preset, 1)
        controls.addWidget(button('Choose report isolates…', self.choose_report_cohort, True))
        controls.addWidget(button('Customize sections…', self.customize_report))
        controls.addWidget(button('Save template', self.save_report_template))
        layout.addLayout(controls)
        self.report_count = label("Report scope: all project samples", "cardTitle")
        layout.addWidget(self.report_count)
        self.report_table = make_table(["Include", "Highlight", "Sample", "Organism", "ST", "AMR genes", "Group"])
        self.report_table.setParent(self)
        self.report_table.hide()
        self.report_table.setColumnWidth(0, 65)
        self.report_table.setColumnWidth(1, 75)
        self.report_table.itemChanged.connect(self.report_item_changed)
        self.report_preview = QTextBrowser()
        self.report_preview.setOpenExternalLinks(False)
        layout.addWidget(self.report_preview, 1)
        self._report_filling = False
        self._report_options_override = None
        panel, content = card()
        row = FlowLayout()
        for title, fmt in [("PDF report…", "pdf"), ("HTML report…", "html"), ("CSV…", "csv"), ("TSV…", "tsv"), ("JSON…", "json")]:
            row.addWidget(button(title, lambda checked=False, f=fmt: self.export_report(f), fmt == "pdf"))
        content.addLayout(row)
        row = FlowLayout()
        row.addWidget(button("Export profile table…", self.export_profile_table))
        row.addWidget(button("Export portable bundle…", self.export_sample_bundle))
        row.addWidget(button("Save project copy…", self.save_project_copy))
        row.addWidget(button('Graph PNG / JPEG…', self.export_report_graph))
        content.addLayout(row)
        content.addWidget(label("Reports describe sequence evidence, not measured susceptibility or proof of transmission. Profile bundles carry results and metadata without sequence files.", "small", True))
        layout.addWidget(panel)

    def report_options(self):
        key = self.report_preset.currentData()
        saved = self.project.get_setting('reports.template.' + key, {})
        return {**REPORT_PRESETS[key], **saved, **(self._report_options_override or {})}

    def report_preset_changed(self):
        self._report_options_override = None
        self.refresh_report_table()

    def choose_report_cohort(self):
        from wmlstudio.cohort_picker import CohortPickerDialog
        dialog = CohortPickerDialog(self.project.samples(), self.project, 'Choose isolates for this report',
                                    selected_ids=self.report_sample_ids(), parent=self, include_reads=True)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.report_ids = set(dialog.selected_ids)
            self.report_scope.setCurrentIndex(0)
            self._report_investigation_snapshot = getattr(self, '_current_snapshot', None)
            self.refresh_report_table()

    def customize_report(self):
        options = self.report_options()
        dialog = QDialog(self)
        dialog.setWindowTitle('Report sections')
        layout = QVBoxLayout(dialog)
        title = QLineEdit(options['title'])
        layout.addWidget(title)
        checks = {}
        for key, text in [('investigation', 'Groups and nearest-neighbour comparisons'), ('qc', 'Input quality evidence'),
                          ('amr', 'AMR determinants'), ('virulence', 'Virulence assay'),
                          ('plasmid_hypotheses', 'Plasmid-marker hypotheses'), ('drug_associations', 'Reference-reported drug classes'),
                          ('graph', 'Graph with focal isolates highlighted'), ('provenance', 'Full technical provenance appendix')]:
            check = QCheckBox(text)
            check.setChecked(options[key])
            checks[key] = check
            layout.addWidget(check)
        layout.addWidget(label('Omitted or unassessed assays are not negative results. Genomic drug annotations do not replace measured susceptibility.', 'small', True))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._report_options_override = {'title': title.text().strip() or options['title'], **{key: check.isChecked() for key, check in checks.items()}}
            self.refresh_report_table()

    def save_report_template(self):
        self.project.set_setting('reports.template.' + self.report_preset.currentData(), self.report_options())
        self.notify('Saved this report template in the project. Its isolate selection remains separate for each report.')

    def report_context(self):
        context = getattr(self, '_report_investigation_snapshot', None)
        if self.report_scope.currentIndex() == 1:
            context = getattr(self, '_current_snapshot', None)
        return context if self.report_options()['investigation'] else None

    def report_graph_png(self, ids=None):
        context = self.report_context()
        current = getattr(self, '_current_snapshot', None)
        if not self.report_options()['graph'] or not context or not current or context['snapshot_id'] != current['snapshot_id']:
            return None
        tree = self.tree
        selected = tree.selected_ids()
        tree.blockSignals(True)
        try:
            focus = set(ids if ids is not None else self.report_sample_ids()) & tree._results.keys()
            tree.select_ids(focus)
            picture = QImage(1800, 1200, QImage.Format.Format_ARGB32)
            painter = QPainter(picture)
            tree._render(painter, 1800, 1200)
            painter.end()
            data = QByteArray()
            buffer = QBuffer(data)
            buffer.open(QIODevice.OpenModeFlag.WriteOnly)
            if not picture.save(buffer, 'PNG'):
                raise OSError('Could not render the report graph.')
            return bytes(data)
        finally:
            tree.select_ids(selected)
            tree.blockSignals(False)

    def export_report_graph(self):
        data = self.report_graph_png()
        if not data:
            self.notify('Choose a comparison-aware report with a matching current graph first.')
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Export focal-isolate graph', 'investigation.png', 'PNG (*.png);;JPEG (*.jpg)')
        if not path:
            return
        try:
            self.check_output(path)
            image = QImage.fromData(data, 'PNG')
            with tempfile.TemporaryDirectory(prefix='.wmlstudio-image-', dir=Path(path).resolve().parent) as directory:
                temporary = Path(directory) / Path(path).name
                if not image.save(str(temporary)):
                    raise OSError('The image could not be written.')
                os.replace(temporary, path)
            self.notify('Exported the graph with this report’s focal isolates highlighted.')
        except Exception as error:
            self.error(error)

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
            if hasattr(self, 'report_preview'):
                if chosen:
                    self.report_preview.setHtml(review_report_html(self.report_records(), selected_ids=chosen,
                        investigation=self.report_context(), options=self.report_options(), graph_png=self.report_graph_png(chosen)))
                else:
                    self.report_preview.setHtml('<h2>Start a focused report</h2><p>Choose a template above, then choose this report’s isolates. No sample is silently included from a library filter.</p>')
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
        self._report_investigation_snapshot = None
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
        """Global project action: never inherit a report/filter/graph selection."""
        records = self.report_records()
        path, _ = QFileDialog.getSaveFileName(self, 'Export complete project', f'wmlstudio-project.{fmt}', f'{fmt.upper()} (*.{fmt})')
        if not path:
            return
        try:
            self.check_output(path)
            if fmt == 'pdf':
                self.write_pdf_report(path, selected_ids={sample['id'] for sample in records}, full_project=True)
            else:
                export_results(records, path, fmt)
            self.notify(f'Exported all {len(records)} project samples. Report selections and visible filters did not change this scope.')
        except Exception as error:
            self.error(error)

    def export_report(self, fmt):
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
            elif fmt == 'html':
                write_review_report(self.report_records(), path, selected_ids=ids, investigation=self.report_context(),
                                    options=self.report_options(), graph_png=self.report_graph_png(ids))
            else:
                export_results(self.report_records(), path, fmt, selected_ids=ids, investigation=self.report_context())
            scope = f"all {len(ids)} project samples" if self.report_scope.currentIndex() == 2 else f"{len(ids)} selected samples"
            self.notify(f"Exported {scope} with their linked metadata and highlights.")
        except Exception as exc:
            self.error(exc)

    def write_pdf_report(self, path, selected_ids=None, *, full_project=False):
        self.check_output(path)
        ids = self.report_sample_ids() if selected_ids is None else set(selected_ids)
        samples = [s for s in self.report_records() if s["id"] in ids]
        if not samples:
            raise ValueError("No samples selected for the PDF report.")
        document = QTextDocument()
        document.setHtml(review_report_html(samples, selected_ids=ids,
                                           investigation=None if full_project else self.report_context(),
                                           options=REPORT_PRESETS['cohort'] if full_project else self.report_options(),
                                           graph_png=None if full_project else self.report_graph_png(ids)))
        with tempfile.TemporaryDirectory(prefix='.wmlstudio-pdf-', dir=Path(path).resolve().parent) as directory:
            output = Path(directory) / 'report.pdf'
            writer = QPdfWriter(str(output))
            writer.setTitle(self.report_options()['title'])
            writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            writer.setResolution(144)
            document.print_(writer)
            del writer
            if not output.is_file() or output.stat().st_size < 100:
                raise OSError('PDF generation failed; the previous destination was preserved.')
            os.replace(output, path)

    def report_records(self):
        return [dict(sample, analyses=self.available_profiles(sample)) for sample in self.project.samples()]

    def build_hydra(self):
        super().build_hydra()
        layout = self.pages.widget(4).widget().layout()
        layout.itemAt(0).widget().setText('Review isolate evidence')
        layout.itemAt(1).widget().setText('Choose this workspace’s cohort. Identity, AMR and accessory evidence remain separate from core allele distances.')
        # Preserve source-import functionality, but move its independent scope
        # to an explicitly advanced tab instead of mixing it with cohort data.
        source_controls = layout.takeAt(2).layout()
        source_summary = layout.takeAt(2).widget()
        old_splitter = layout.takeAt(2).widget()
        source_page = QWidget()
        source_layout = QVBoxLayout(source_page)
        source_layout.addWidget(label('Advanced: the imported source report can contain unlinked samples outside this workspace cohort. Review mapping before associating evidence.', 'small', True))
        source_layout.addLayout(source_controls)
        source_layout.addWidget(source_summary)
        source_layout.addWidget(old_splitter, 1)
        self.hydra_view.setHtml('<h2>Imported source evidence</h2><p>Import a HYDRA report to review its original records and provenance. Its names do not automatically become project-isolate identities.</p>')
        actions = FlowLayout()
        actions.addWidget(button('Choose evidence isolates…', self.choose_feature_cohort, True))
        actions.addWidget(button('Run HYDRA…', self.choose_hydra_cohort))
        tools = QToolButton()
        tools.setText('Evidence actions ▾')
        tools.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(tools)
        menu.addAction('AMR database snapshots…', self.open_amr_databases)
        menu.addAction('Link imported source report…', self.link_hydra_samples)
        menu.addAction('Refresh scoped evidence', self.refresh_features)
        for title, kind in [('Export feature table…', 'features'), ('Export AMR matrix…', 'amr'), ('Export characterization table…', 'characterization')]:
            menu.addAction(title, lambda checked=False, kind=kind: self.export_feature_table(kind))
        tools.setMenu(menu)
        actions.addWidget(tools)
        layout.addLayout(actions)
        self.feature_scope_label = label('No evidence isolates selected. This scope is independent of the library, comparison and reports.', 'small', True)
        layout.addWidget(self.feature_scope_label)
        tabs = QTabWidget()
        self.evidence_tabs = tabs
        self.feature_table = QTableView()
        self.feature_model = EvidenceMatrixModel(self.feature_table)
        self.feature_table.setModel(self.feature_model)
        self.feature_table.setSortingEnabled(True)
        self.feature_table.setAlternatingRowColors(True)
        self.feature_table.doubleClicked.connect(lambda index: self.open_isolate_record(self.feature_model.rows[index.row()]['_sample_id']))
        tabs.addTab(self.feature_table, 'Isolate feature summary')
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
        self.amr_matrix.setSortingEnabled(True)
        self.amr_matrix.doubleClicked.connect(lambda index: self.open_isolate_record(self.amr_model.rows[index.row()]['_sample_id']))
        matrix_layout.addWidget(self.amr_matrix)
        matrix_layout.addWidget(label("Present = detected in the linked report; not detected is not a susceptibility claim. No report = unknown, never absence.", "small", True))
        tabs.addTab(matrix_page, "AMR gene matrix")
        tabs.addTab(source_page, 'Advanced · imported source')
        layout.addWidget(tabs, 1)

    def restore_hydra(self):
        super().restore_hydra()
        if hasattr(self, "feature_table"):
            self.refresh_features()

    def refresh_features(self):
        if not hasattr(self, "feature_model"):
            return
        rows = []
        from wmlstudio.metadata_ingest import annotation_values
        samples = [sample for sample in self.project.samples() if sample['id'] in getattr(self, 'feature_ids', set())
                   and sample.get('metadata', {}).get('workflow', {}).get('source_kind') != 'read_mate']
        for sample in samples:
            genus, species, status = organism_for(sample)
            result = sample.get("result") or {}
            amr_state = hydra_evidence_status(sample)
            row = {'_sample_id': sample['id'], "Sample": sample["name"], "Genus": genus or "Unknown", "Species": species or "—", "Organism evidence": status,
                   "ST": result.get("st"), "Scheme": result.get("scheme"), "AMR genes": "; ".join(gene_names(sample)),
                   "HYDRA linked": bool(sample.get("metadata", {}).get("hydra")), "AMR evidence state": amr_state["status"],
                   "AMR evidence note": amr_state["reason"], **{'Annotation.' + key: value for key, value in annotation_values(sample).items()}}
            row.update({"QC." + key: value for key, value in result.get("qc", {}).items() if not isinstance(value, (dict, list))})
            rows.append(row)
        main = ["Sample", "Genus", "Species", "Organism evidence", "ST", "Scheme", "AMR genes", "HYDRA linked", "AMR evidence state", "AMR evidence note"]
        self.feature_model.replace(main + sorted({key for row in rows for key in row} - set(main) - {'_sample_id'}), rows)
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
            matrix.append({'_sample_id': sample['id'], "Sample": sample["name"], "Evidence state": state, **{gene: state_cell(gene) for gene in genes}})
        self.amr_model.replace(["Sample", "Evidence state", *genes], matrix)
        self.feature_scope_label.setText(f'{len(samples)} explicitly chosen isolates · {len(genes)} visible AMR-gene columns · double-click an isolate for its complete record. No report is not a negative result.')

    def export_feature_table(self, kind='features'):
        path, _ = QFileDialog.getSaveFileName(self, 'Export this evidence cohort table', f'evidence-{kind}.tsv', 'TSV (*.tsv);;CSV (*.csv)')
        if not path:
            return
        try:
            self.check_output(path)
            self.write_feature_table(path, kind)
            self.notify('Exported only the evidence workspace cohort and visible table columns; core distances were not changed.')
        except Exception as error:
            self.error(error)

    def write_feature_table(self, path, kind='features'):
        import csv

        from wmlstudio.export import _atomic_text, _csv_cell
        if kind == 'characterization':
            table = self.characterization_table
            headers = ['sample_id', *[table.horizontalHeaderItem(column).text() for column in range(table.columnCount())]]
            rows = [[table.item(row, 0).data(Qt.ItemDataRole.UserRole), *[table.item(row, column).text() for column in range(table.columnCount())]] for row in range(table.rowCount())]
        else:
            model = self.feature_model if kind == 'features' else self.amr_model
            headers = ['sample_id', *model.headers]
            rows = [[row['_sample_id'], *[row.get(key) for key in model.headers]] for row in model.rows]
        with _atomic_text(path, newline='') as handle:
            writer = csv.writer(handle, delimiter=',' if Path(path).suffix.casefold() == '.csv' else '\t')
            writer.writerow(headers)
            writer.writerows([_csv_cell(value) for value in row] for row in rows)

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
