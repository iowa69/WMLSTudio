"""Selection-scoped reports, linked AMR matrices and interoperable sample exchange."""


import os
import tempfile
from pathlib import Path

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt, QTimer
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
    QMessageBox,
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
from wmlstudio.theme import BACKGROUND
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.ui_compare import EvidenceMatrixModel
from wmlstudio.widgets import button, card, label

# How long the one-click summary waits for a background comparison before it
# writes the document anyway, saying plainly that no comparison was available.
_COMPARISON_WAIT_STEP_MS = 250
_COMPARISON_WAIT_LIMIT_MS = 120_000

# Three quantities, three methods, never one axis. A report built on allele
# distances says where the third one is rather than leaving a reader to assume
# every number in the workspace belongs on the same scale.
METHOD_SEPARATION = (
    "SNP distances from the SNP tree are a third quantity, produced by a different method over the "
    "split k-mers two isolates share rather than over loci. They are not in this report and are "
    "never comparable with the figures above.")

PLASMID_SCOPE = (
    "Plasmid evidence in a report is a screen over replicon markers on assembled contigs: a marker "
    "is evidence about the contig it sits on, not a reconstructed plasmid, a plasmid count or a "
    "mobility prediction, and it is not MOB-suite.")

# How the engine's point-mutation levels read to somebody deciding what a blank
# mutation column means. "none" is not "no mutations": it is "nothing to look in".
MUTATION_WORDS = {
    "dna_and_protein": "DNA and protein catalogues",
    "protein_only": "protein catalogue only — no DNA-level target was assessed",
    "dna_only": "DNA catalogue only — no curated protein mutation was assessed",
    "none": "no catalogue for this organism — genes only, which is not an absence of mutations",
    "unknown": "not recorded by this run",
}


def amr_method_fields(sample):
    """Which reference sets, release and searches produced this isolate's AMR evidence.

    A determinant list means nothing without the question it answers: a run
    against two reference sets and a run against ten look identical in a table of
    gene names. The run's own provenance is read, never inferred, and an older
    result that recorded none of this says so rather than borrowing today's.
    """
    from wmlstudio.sample_workflow import current_hydra_evidence

    evidence = current_hydra_evidence(sample)
    if not evidence:
        return dict.fromkeys(("AMR reference sets", "AMR reference release",
                              "Point mutations searched", "Virulence elements searched"),
                             "No current AMR evidence")
    execution = evidence.get("execution_provenance") or {}
    organism = execution.get("organism") or {}
    virulence = execution.get("virulence") or {}
    release = (execution.get("reference_release") or {}).get("release")
    databases = evidence.get("databases") or list(
        ((execution.get("reference_snapshot") or {}).get("databases") or {}))
    return {"AMR reference sets": "; ".join(sorted(str(name) for name in databases)) or "not recorded",
            "AMR reference release": str(release or "not recorded"),
            "Point mutations searched": MUTATION_WORDS.get(
                str(organism.get("point_mutation_level") or "unknown"), "not recorded by this run"),
            "Virulence elements searched": (
                "yes, curated for " + str(organism.get("resolved") or "this isolate's organism")
                if virulence.get("enabled") and virulence.get("organism_curated") else
                "yes, with no organism curation" if virulence.get("enabled") else
                "no — the translated search was acquired resistance only"
                if virulence else "not recorded by this run")}


class ReportWorkspaceMixin:
    def build_reports(self):
        _, layout = self.page()
        self.heading(layout, "Prepare a review report", "Choose a report template and its own isolates. Clear evidence summaries come first; technical provenance is optional.")
        # The one-click route, above the template controls, because a reader who
        # is not a bioinformatician should not have to choose a template first.
        simple, simple_layout = card()
        simple_layout.addWidget(label('Just give me a summary', 'cardTitle'))
        simple_layout.addWidget(label(
            'One short document in plain words: a picture of how close the isolates are, the resistance '
            'genes that were found, and each isolate’s closest match. It prints which isolates it covered '
            'and where that choice came from, says which typing it is based on and how many targets the '
            'reference holds, and it never calls a gene a susceptibility result.',
            'small', True))
        simple_row = FlowLayout()
        simple_row.addWidget(button('Make a simple summary (PDF)…', lambda: self.simple_report(), True))
        simple_layout.addLayout(simple_row)
        layout.addWidget(simple)
        controls = FlowLayout()
        self.report_scope = QComboBox(self)
        self.report_scope.hide()  # Kept as the explicit scope model, not a global selector.
        self.report_scope.addItems(["Explicit report selection", "Current comparison cohort", "All project samples"])
        self.report_scope.setCurrentIndex(0)
        self.report_scope.currentIndexChanged.connect(self.refresh_report_table)
        self.report_preset = QComboBox()
        for key, value in REPORT_PRESETS.items():
            self.report_preset.addItem(value['title'], key)
        # By key, not by position: appending a preset must never repoint the default.
        self.report_preset.setCurrentIndex(self.report_preset.findData('cohort'))
        self.report_preset.currentIndexChanged.connect(self.report_preset_changed)
        controls.addWidget(self.report_preset, 1)
        controls.addWidget(button('Choose report isolates…', self.choose_report_cohort, True))
        controls.addWidget(button('Customize sections…', self.customize_report))
        controls.addWidget(button('Save template', self.save_report_template))
        layout.addLayout(controls)
        self.report_count = label("Report scope: all project samples", "cardTitle")
        layout.addWidget(self.report_count)
        # Which typing this report is about, on the page and in the document, from
        # the one call the document itself prints: a cgMLST report and an MLST
        # report must never be mistakable for one another, and the page is where a
        # reader finds out before they export.
        self.report_typing_label = label("", "small", True)
        layout.addWidget(self.report_typing_label)
        self.report_table = make_table(["Include", "Highlight", "Sample", "Organism", "ST", "AMR genes", "Group"])
        self.report_table.setColumnWidth(0, 65)
        self.report_table.setColumnWidth(1, 75)
        # Visible and compact: the reader can see exactly which isolates the
        # report covers, tick one in or out, and right-click to act on them.
        self.report_table.setMaximumHeight(214)
        self.report_table.itemChanged.connect(self.report_item_changed)
        layout.addWidget(self.report_table)
        self.install_view_menu('reports', self.report_table)
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
        # Stated on the page as well as in the document, so a reader knows what
        # to look for before they open the PDF.
        content.addWidget(label("Every report that shows a tree, a cluster or a distance names the reference and its "
                                "full target count, the threshold in force, and whether that threshold is a published "
                                "cutoff you adopted or your own setting — with the citation, its DOI and the authors’ "
                                "own caveat when it is published. A sequence type and a core-genome profile are "
                                "reported as separate quantities and never share a threshold. " + METHOD_SEPARATION,
                                "small", True))
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
                          ('plasmid_hypotheses', 'Plasmid evidence: replicon markers on assembled contigs'),
                          ('drug_associations', 'Reference-reported drug classes'),
                          ('graph', 'Graph with focal isolates highlighted'),
                          ('graph_jpeg', 'Embed the picture as JPEG (smaller file, slightly softer text)'),
                          ('provenance', 'Full technical provenance appendix')]:
            check = QCheckBox(text)
            check.setChecked(options[key])
            checks[key] = check
            layout.addWidget(check)
        layout.addWidget(label('Omitted or unassessed assays are not negative results. Genomic drug annotations do not replace measured susceptibility.', 'small', True))
        layout.addWidget(label(PLASMID_SCOPE, 'small', True))
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
        options = self.report_options()
        if not options['investigation']:
            return None
        context = getattr(self, '_report_investigation_snapshot', None)
        current = getattr(self, '_current_snapshot', None)
        if self.report_scope.currentIndex() == 1:
            context = current
        if options.get('layout') == 'one_page' and (
                context is None or current is None
                or context.get('snapshot_id') != current.get('snapshot_id')):
            # The simple summary always describes the comparison that is on
            # screen, so its picture, its threshold and its distances cannot
            # disagree. Every other preset keeps the strict stored context.
            context = current
        return context

    def report_typing(self):
        """What this report will say it measured, or None when it measures nothing.

        Built by export.threshold_provenance — the same call the document itself
        prints — so the page and the PDF can never state two different schemes,
        two different target counts or two different thresholds.
        """
        context = self.report_context()
        if not context:
            return None
        from wmlstudio.export import threshold_provenance
        return threshold_provenance(context)

    def report_typing_line(self):
        """One sentence naming the typing, the reference, its size and the threshold."""
        provenance = self.report_typing()
        if provenance is None:
            return ("No comparison is attached to this report, so it will state no typing scale, no "
                    "reference and no threshold. That is not a statement that the isolates are "
                    "unrelated.")
        typing = provenance["typing"]
        return (f"This report is {typing['label']}: {provenance['scheme']} · "
                f"{provenance['total_loci']} targets in the reference. {provenance['statement']} "
                f"{typing['note']} {METHOD_SEPARATION}")

    def report_typing_kind(self):
        """'mlst', 'cgmlst' or '' — the word that keeps two reports apart on a desk."""
        provenance = self.report_typing()
        kind = ((provenance or {}).get("typing") or {}).get("kind", "")
        return kind if kind in {"mlst", "cgmlst"} else ""

    def report_filename(self, stem, suffix):
        """A default filename that names the typing, so the two cannot be swapped."""
        kind = self.report_typing_kind()
        return f"wmlstudio-{kind}-{stem}.{suffix}" if kind else f"wmlstudio-{stem}.{suffix}"

    def report_graph_available(self):
        """Whether a picture would describe the same snapshot as the numbers."""
        context = self.report_context()
        current = getattr(self, '_current_snapshot', None)
        return bool(self.report_options()['graph'] and context and current
                    and context['snapshot_id'] == current['snapshot_id'])

    def report_graph_image(self, ids=None, *, fmt='PNG'):
        """Focal isolates highlighted on the graph that is actually on screen."""
        if not self.report_graph_available():
            return None
        tree = self.tree
        selected = tree.selected_ids()
        tree.blockSignals(True)
        try:
            focus = set(ids if ids is not None else self.report_sample_ids()) & tree._results.keys()
            tree.select_ids(focus)
            # JPEG carries no alpha channel, so the canvas is opaque from the
            # start; _render fills the same background over it.
            picture = QImage(1800, 1200, QImage.Format.Format_RGB32)
            picture.fill(QColor(BACKGROUND))
            painter = QPainter(picture)
            tree._render(painter, 1800, 1200)
            painter.end()
            data = QByteArray()
            buffer = QBuffer(data)
            buffer.open(QIODevice.OpenModeFlag.WriteOnly)
            if not picture.save(buffer, fmt, 92 if fmt == 'JPEG' else -1):
                raise OSError('Could not render the report graph.')
            return bytes(data)
        finally:
            tree.select_ids(selected)
            tree.blockSignals(False)

    def report_graph_png(self, ids=None):
        return self.report_graph_image(ids, fmt='PNG')

    def export_report_graph(self):
        if not self.report_graph_available():
            self.notify('Choose a comparison-aware report with a matching current graph first.')
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Export focal-isolate graph',
                                              self.report_filename('graph', 'png'),
                                              'PNG (*.png);;JPEG (*.jpg *.jpeg)')
        if not path:
            return
        fmt = 'JPEG' if Path(path).suffix.casefold() in {'.jpg', '.jpeg'} else 'PNG'
        try:
            self.check_output(path)
            data = self.report_graph_image(fmt=fmt)
            if not data:
                raise OSError('The image could not be written.')
            with tempfile.TemporaryDirectory(prefix='.wmlstudio-image-', dir=Path(path).resolve().parent) as directory:
                temporary = Path(directory) / Path(path).name
                temporary.write_bytes(data)
                os.replace(temporary, path)
            # A picture carries no scale once it leaves the application, so the
            # typing it was drawn on is named here and in its filename.
            self.notify(f'Exported the graph as {fmt} with this report’s focal isolates highlighted. '
                        + self.report_typing_line())
        except Exception as error:
            self.error(error)

    def resolve_report_scope(self, selected_ids=None):
        """Decide a report's scope once, and state that decision inside the report.

        A fallback is announced in the document, never inherited silently from a
        view filter. An explicit caller wins; then the isolates chosen for this
        report; then whatever the workspace has focused, named by where it came
        from; and only then every isolate in the project.
        """
        records = self.report_records()
        known = {sample['id'] for sample in records}
        eligible = {sample['id'] for sample in records
                    if sample.get('metadata', {}).get('workflow', {}).get('source_kind') != 'read_mate'}
        if selected_ids is not None:
            ids = {str(value) for value in selected_ids} & known
            if not ids:
                raise ValueError('None of the isolates given for this report are in the project. '
                                 'Choose this report’s isolates, or import sequence files first.')
            return {'ids': ids, 'implicit': False,
                    'note': f'This report covers {len(ids)} isolate(s) chosen by the caller.'}
        ids = self.report_sample_ids() & eligible
        if ids:
            return {'ids': ids, 'implicit': False,
                    'note': f'This report covers the {len(ids)} isolate(s) you chose for it.'}
        focus = getattr(self, 'focus', None)
        focused = {str(value) for value in getattr(focus, 'ids', ()) or ()} & eligible
        if focused:
            origin = getattr(focus, 'origin', '') or 'your current selection'
            return {'ids': focused, 'implicit': True,
                    'note': (f'No isolates were chosen for this report, so it covers the {len(focused)} '
                             f'isolate(s) you had selected — from {origin}. '
                             'Use “Choose report isolates…” to set the scope yourself.')}
        if not eligible:
            raise ValueError('This project has no isolates yet. Import sequence files before making a report.')
        return {'ids': eligible, 'implicit': True,
                'note': (f'No isolates were chosen for this report, so it covers all {len(eligible)} '
                         'isolates in the project. Use “Choose report isolates…” to narrow it.')}

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
            self.report_typing_label.setText(self.report_typing_line())
            if hasattr(self, 'report_preview'):
                if chosen:
                    options = self.report_options()
                    jpeg = bool(options.get('graph_jpeg'))
                    self.report_preview.setHtml(review_report_html(self.report_records(), selected_ids=chosen,
                        investigation=self.report_context(), options=options,
                        graph_png=self.report_graph_image(chosen, fmt='JPEG' if jpeg else 'PNG'),
                        graph_mime='image/jpeg' if jpeg else 'image/png'))
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

    # The right-click handlers these three views need (open a record, route a
    # selection to another tab, archive, export) all live on WorkbenchMixin,
    # which sits earlier in the MainWindow MRO. Wiring the views is what is
    # needed here; a second copy of those handlers would only shadow-compete.

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
        try:
            scope = self.resolve_report_scope()
        except ValueError as error:
            self.notify(str(error))
            return
        ids = scope['ids']
        path, _ = QFileDialog.getSaveFileName(self, f"Export {len(ids)} samples",
                                              self.report_filename("report", fmt),
                                              f"{fmt.upper()} (*.{fmt})")
        if not path:
            return
        try:
            self.check_output(path)
            options = self.report_options()
            jpeg = bool(options.get('graph_jpeg'))
            if fmt == "pdf":
                self.write_pdf_report(path, selected_ids=ids, scope_note=scope['note'],
                                      scope_implicit=scope['implicit'])
            elif fmt == 'html':
                write_review_report(self.report_records(), path, selected_ids=ids, investigation=self.report_context(),
                                    options=options, graph_png=self.report_graph_image(ids, fmt='JPEG' if jpeg else 'PNG'),
                                    graph_mime='image/jpeg' if jpeg else 'image/png',
                                    scope_note=scope['note'], scope_implicit=scope['implicit'])
            else:
                export_results(self.report_records(), path, fmt, selected_ids=ids, investigation=self.report_context())
            self.notify(f"Exported {len(ids)} samples with their linked metadata and highlights. " + scope['note'])
        except Exception as exc:
            self.error(exc)

    def write_pdf_report(self, path, selected_ids=None, *, full_project=False, options=None, graph=None,
                         scope_note=None, scope_implicit=False):
        """Print the chosen report to PDF, stating in the document which isolates it covered.

        ``options`` and ``graph`` (an ``(image bytes, mime)`` pair) let the
        one-click summary pin its own preset and its own already-rendered
        picture; both default to today's behaviour.
        """
        self.check_output(path)
        scope = self.resolve_report_scope(selected_ids)
        ids = scope['ids']
        samples = [s for s in self.report_records() if s["id"] in ids]
        settings = options if options is not None else (REPORT_PRESETS['cohort'] if full_project else self.report_options())
        if graph is not None:
            image, mime = graph
        elif full_project:
            image, mime = None, 'image/png'
        else:
            jpeg = bool(settings.get('graph_jpeg'))
            image = self.report_graph_image(ids, fmt='JPEG' if jpeg else 'PNG')
            mime = 'image/jpeg' if jpeg else 'image/png'
        document = QTextDocument()
        document.setHtml(review_report_html(samples, selected_ids=ids,
                                           investigation=None if full_project else self.report_context(),
                                           options=settings, graph_png=image, graph_mime=mime,
                                           scope_note=scope_note or scope['note'],
                                           scope_implicit=scope_implicit or (scope['implicit'] and scope_note is None)))
        with tempfile.TemporaryDirectory(prefix='.wmlstudio-pdf-', dir=Path(path).resolve().parent) as directory:
            output = Path(directory) / 'report.pdf'
            writer = QPdfWriter(str(output))
            # The document's own title carries the typing too: a PDF filed away
            # for a week must still say which quantity it measured.
            writer.setTitle(settings['title'] + self.pdf_title_suffix(full_project))
            writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            writer.setResolution(144)
            document.print_(writer)
            del writer
            if not output.is_file() or output.stat().st_size < 100:
                raise OSError('PDF generation failed; the previous destination was preserved.')
            os.replace(output, path)

    def pdf_title_suffix(self, full_project=False):
        """' — cgMLST, 2358 targets' for the PDF's own title, or '' when unknown."""
        provenance = None if full_project else self.report_typing()
        typing = (provenance or {}).get('typing') or {}
        if typing.get('kind') not in {'mlst', 'cgmlst'}:
            return ''
        word = 'cgMLST' if typing['kind'] == 'cgmlst' else 'classical MLST'
        return f" — {word}, {provenance['total_loci']} targets"

    def report_records(self):
        return [dict(sample, analyses=self.available_profiles(sample)) for sample in self.project.samples()]

    # ------------------------------------------------------------------
    # The one-click plain-language summary
    # ------------------------------------------------------------------
    def simple_report(self, path=None, *, build_comparison=None):
        """One click: the picture, the genes found, and each isolate's closest match.

        ``build_comparison`` decides whether to build a comparison when none
        exists; ``None`` asks the user, so a test never meets a modal dialog.
        """
        try:
            scope = self.resolve_report_scope()
        except ValueError as error:
            self.notify(str(error))
            return
        if getattr(self, '_current_snapshot', None) is None:
            if build_comparison is None:
                build_comparison = self._ask_to_build_comparison(len(scope['ids']))
            if build_comparison:
                from wmlstudio import workspace_focus
                # An explicit, announced cohort change: the picture then shows
                # exactly the isolates this report is about, and the Compare tab
                # records where its cohort came from.
                workspace_focus.send_selection(self, 'compare', sorted(scope['ids']),
                                               origin='Simple summary report', navigate=False)
                if getattr(self, '_current_snapshot', None) is None:
                    self._await_comparison_for_summary(path, getattr(self, '_comparison_generation', 0))
                    return
        self._write_simple_report(path, scope)

    def _ask_to_build_comparison(self, count):
        answer = QMessageBox.question(
            self, 'Build a comparison first?',
            f'No comparison has been built yet, so the summary would carry no picture and no '
            f'closest-match table.\n\nBuild one now from these {count} isolate(s)? This replaces the '
            f'current comparison cohort.\n\nChoosing No still writes the summary; it will say plainly '
            f'that no comparison was available.',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes)
        return answer == QMessageBox.StandardButton.Yes

    def _await_comparison_for_summary(self, path, generation, waited=0):
        """Finish the summary once a background comparison lands, or say it did not.

        The comparison belongs to the Compare workspace; this watches the state
        it publishes rather than reaching into its worker, abandons quietly when
        the user starts a different comparison or closes the project, and never
        leaves the click without a document.
        """
        if getattr(self, '_comparison_closing', False) or generation != getattr(self, '_comparison_generation', generation):
            return
        idle = (getattr(self, 'comparison_worker', None) is None
                and getattr(self, '_comparison_pending', None) is None)
        ready = getattr(self, '_current_snapshot', None) is not None
        if not (ready or waited >= _COMPARISON_WAIT_LIMIT_MS or (idle and waited >= 3 * _COMPARISON_WAIT_STEP_MS)):
            QTimer.singleShot(_COMPARISON_WAIT_STEP_MS, lambda: self._await_comparison_for_summary(
                path, generation, waited + _COMPARISON_WAIT_STEP_MS))
            return
        if not ready:
            self.notify('The comparison did not produce a snapshot for these isolates. '
                        'The summary says so; it does not report them as unrelated.')
        try:
            scope = self.resolve_report_scope()
        except ValueError as error:
            self.notify(str(error))
            return
        self._write_simple_report(path, scope)

    def _write_simple_report(self, path, scope):
        options = {**REPORT_PRESETS['simple']}
        # One click must produce the same document whatever preset the combo was
        # left on, and report_context()/report_graph_image() both read
        # report_options(). Evaluate the whole write — the offered filename
        # included, since that is where the typing is named — against the simple
        # preset.
        previous = self._report_options_override
        self._report_options_override = options
        try:
            if path is None:
                path, _ = QFileDialog.getSaveFileName(self, 'Save simple summary',
                                                      self.report_filename('summary', 'pdf'),
                                                      'PDF (*.pdf);;HTML (*.html)')
                if not path:
                    return
            self.check_output(path)
            investigation = self.report_context()
            image = self.report_graph_image(scope['ids'], fmt='JPEG')
            if Path(path).suffix.casefold() == '.pdf':
                self.write_pdf_report(path, selected_ids=scope['ids'], options=options,
                                      graph=(image, 'image/jpeg'), scope_note=scope['note'],
                                      scope_implicit=scope['implicit'])
            else:
                records = [s for s in self.report_records() if s['id'] in scope['ids']]
                write_review_report(records, path, selected_ids=scope['ids'], investigation=investigation,
                                    options=options, graph_png=image, graph_mime='image/jpeg',
                                    scope_note=scope['note'], scope_implicit=scope['implicit'])
        except Exception as error:
            self.error(error)
            return
        finally:
            self._report_options_override = previous
        self.notify('Simple summary saved. ' + scope['note'])

    def build_hydra(self):
        super().build_hydra()
        layout = self.pages.page_for('evidence').widget().layout()
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
        for title, kind in [('Export feature table…', 'features'), ('Export AMR matrix…', 'amr'),
                            ('Export characterization table…', 'characterization'),
                            ('Export replicon markers across the cohort…', 'plasmid_replicons'),
                            ('Export replicon / determinant co-location…', 'plasmid_pairs')]:
            menu.addAction(title, lambda checked=False, kind=kind: self.export_feature_table(kind))
        tools.setMenu(menu)
        actions.addWidget(tools)
        layout.addLayout(actions)
        self.feature_scope_label = label('No evidence isolates selected. This scope is independent of the library, comparison and reports.', 'small', True)
        layout.addWidget(self.feature_scope_label)
        tabs = QTabWidget()
        self.evidence_tabs = tabs
        self.pages.register_subtabs('evidence', self.evidence_tabs)
        self.feature_table = QTableView()
        self.feature_model = EvidenceMatrixModel(self.feature_table)
        self.feature_table.setModel(self.feature_model)
        self.feature_table.setSortingEnabled(True)
        self.feature_table.setAlternatingRowColors(True)
        self.feature_table.doubleClicked.connect(lambda index: self.open_isolate_record(self.feature_model.rows[index.row()]['_sample_id']))
        self.install_view_menu('evidence.features', self.feature_table)
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
        self.install_view_menu('evidence.amr', self.amr_matrix)
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
                   "AMR evidence note": amr_state["reason"], **amr_method_fields(sample),
                   **{'Annotation.' + key: value for key, value in annotation_values(sample).items()}}
            row.update({"QC." + key: value for key, value in result.get("qc", {}).items() if not isinstance(value, (dict, list))})
            rows.append(row)
        main = ["Sample", "Genus", "Species", "Organism evidence", "ST", "Scheme", "AMR genes",
                "HYDRA linked", "AMR evidence state", "AMR evidence note", "AMR reference sets",
                "AMR reference release", "Point mutations searched", "Virulence elements searched"]
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
        # The widget tables are written exactly as they are read on screen, each
        # count beside the denominator it was counted against. The two plasmid
        # tables stay two files: merging replicon counts with co-location counts
        # would put two different observations under one heading.
        widgets = {'characterization': getattr(self, 'characterization_table', None),
                   'plasmid_replicons': getattr(self, 'plasmid_replicon_table', None),
                   'plasmid_pairs': getattr(self, 'plasmid_pair_table', None)}
        if kind in widgets:
            table = widgets[kind]
            if table is None:
                raise ValueError('That evidence table has not been built in this window yet.')
            headers = [*(['sample_id'] if kind == 'characterization' else []),
                       *[table.horizontalHeaderItem(column).text() for column in range(table.columnCount())]]
            rows = [[*([table.item(row, 0).data(Qt.ItemDataRole.UserRole)] if kind == 'characterization' else []),
                     *[table.item(row, column).text() if table.item(row, column) is not None else ''
                       for column in range(table.columnCount())]] for row in range(table.rowCount())]
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
