"""Publication context first; explicit, auditable local adaptation second."""

import html
import json
from pathlib import Path
from urllib.parse import urlencode

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QFormLayout, QHBoxLayout, QLineEdit, QListWidget, QListWidgetItem,
    QPlainTextEdit, QSpinBox, QSplitter, QTextBrowser, QVBoxLayout, QWidget)

from .threshold_guidance import CATALOG_VERSION, ORGANISMS, REVIEWED_ON, guidance_for, record_decision
from .ui_common import organism_for
from .widgets import button, label


class ThresholdGuideDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window, self.evidence = window, None
        self.setWindowTitle("Published cluster guidance · research review")
        self.resize(1060, 740)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Which published context fits your question?", "title", True))
        layout.addWidget(label(f"Catalog {CATALOG_VERSION} · reviewed {REVIEWED_ON}. Dated curation, not a live clinical guideline. Numbers never apply automatically.", "muted", True))
        top = QHBoxLayout()
        self.organism = QComboBox()
        self.organism.addItems(sorted(ORGANISMS))
        top.addWidget(self.organism, 1)
        top.addWidget(button("Check newer publications online", self.search_newer))
        layout.addLayout(top)
        split = QSplitter(Qt.Orientation.Horizontal)
        self.entries = QListWidget()
        self.entries.setMinimumWidth(260)
        self.text = QTextBrowser()
        self.text.setOpenLinks(False)
        self.text.anchorClicked.connect(self.open_source)
        split.addWidget(self.entries)
        split.addWidget(self.text)
        split.setSizes([325, 680])
        layout.addWidget(split, 1)
        self.apply_value = QCheckBox("Adapt this published number to the current cgMLST comparison…")
        layout.addWidget(self.apply_value)
        self.advanced = QWidget()
        form = QFormLayout(self.advanced)
        self.threshold = QSpinBox()
        self.threshold.setRange(0, 100000)
        form.addRow("Local inclusive cutoff ≤", self.threshold)
        self.bound_scheme = QLineEdit()
        self.bound_scheme.setPlaceholderText("Enter the exact scheme key shown above after verifying its target set/source")
        form.addRow("Verified scheme binding", self.bound_scheme)
        self.schema_review = QCheckBox("I verified the exact published target set and the local reference snapshot.")
        self.protocol_review = QCheckBox("I reviewed caller, missing-data policy, minimum overlap and clustering differences.")
        self.epi_review = QCheckBox("I reviewed the question, time window, epidemiological setting and QC limitations.")
        form.addRow(self.schema_review)
        form.addRow(self.protocol_review)
        form.addRow(self.epi_review)
        self.justification = QPlainTextEdit()
        self.justification.setMaximumHeight(85)
        self.justification.setPlaceholderText("Why this threshold? Record local validation or why an exploratory adaptation is appropriate; describe departures from the source.")
        form.addRow("Protocol justification", self.justification)
        self.advanced.hide()
        self.apply_value.toggled.connect(self.advanced.setVisible)
        layout.addWidget(self.advanced)
        self.feedback = label("Saving a citation alone does not change the threshold. Local adaptations remain labelled unvalidated for this application.", "small", True)
        layout.addWidget(self.feedback)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        controls.accepted.connect(self.accept)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)
        self.organism.currentTextChanged.connect(self.refresh_entries)
        self.entries.currentRowChanged.connect(self.show_entry)
        cohort = [s for s in window.project.samples() if s['id'] in (window.cohort_ids or set())]
        taxa = {" ".join(organism_for(s)[:2]).strip() for s in cohort}
        if len(taxa) == 1 and next(iter(taxa)) in ORGANISMS:
            self.organism.setCurrentText(next(iter(taxa)))
        self.refresh_entries()

    def context(self):
        results = getattr(self.window, '_last_comparison', [])
        cache = getattr(self.window, '_comparison_cache', None)
        current_key = self.window._comparison_key(self.window._comparison_request(), self.window.project.comparison_revision())
        if not cache or cache[0] != current_key:
            results = []  # A changed/pending selector cannot inherit the old displayed reference.
        digests = {r.get('scheme_digest') for r in results}
        counts = {len(r.get('alleles') or {}) for r in results}
        cohort = [s for s in self.window.project.samples() if s['id'] in (self.window.cohort_ids or set())]
        taxa = {" ".join(organism_for(s)[:2]).strip() for s in cohort}
        return {'method': getattr(self.window, 'comparison_mode', 'unspecified'),
                'organism': next(iter(taxa)) if len(taxa) == 1 else 'Mixed / unknown cohort',
                'scheme_digest': next(iter(digests)) if len(digests) == 1 else None,
                'locus_count': next(iter(counts)) if len(counts) == 1 else None,
                'scheme_key': self.bound_scheme.text().strip(),
                'caller': 'WMLSTudio saved profile caller; inspect each result provenance',
                'missing_policy': 'Shared unambiguous calls / union of profile loci; low-overlap pairs excluded',
                'min_overlap': self.window.overlap.value() if hasattr(self.window, 'overlap') else None,
                'clustering': 'all-comparable-pair single linkage',
                'sample_ids': sorted(s['id'] for s in cohort)}

    def refresh_entries(self):
        self.entries.clear()
        guidance = guidance_for(self.organism.currentText())
        for entry in guidance['entries']:
            value = entry['published_threshold']
            title = f"{entry['method']} · {'No numeric rule' if value is None else '≤ ' + str(value)} · {entry['source']['published']}"
            item = QListWidgetItem(title)
            item.setData(Qt.ItemDataRole.UserRole, entry)
            item.setToolTip(entry['scope'])
            self.entries.addItem(item)
        self.apply_value.setChecked(False)
        if self.entries.count():
            self.entries.setCurrentRow(0)
        else:
            self.text.setPlainText(guidance['message'] + '\n\n' + guidance['interpretation'])
            self.apply_value.setEnabled(False)

    def show_entry(self, row):
        item = self.entries.item(row)
        if item is None:
            return
        entry = item.data(Qt.ItemDataRole.UserRole)
        e = lambda value: html.escape(str(value))
        source = entry['source']
        value = entry['published_threshold']
        body = f"<h2>{e(entry['organism'])}</h2><h3>{e(entry['method'])}: {e('No numeric rule curated' if value is None else '≤ ' + str(value) + ' ' + entry['unit'])}</h3>"
        body += f"<p>{e(entry['scope'])}</p><p><b>Scheme key:</b> {e(entry['scheme_key'] or 'Not curated; citation only')}<br><b>Target count:</b> {e(entry['locus_count'] or 'Not established here')}<br><b>Missing data:</b> {e(entry['missing_policy'])}</p>"
        body += f"<p><b>Limitations:</b> {e(entry['limitations'])}</p><p>{e(source['citation'])}<br><a href='{e(source['url'])}'>{e(source['doi'])}</a><br>{e(source['locator'])}</p>"
        body += "<p><b>Your protocol:</b></p><pre>" + e(json.dumps(self.context(), indent=2)) + "</pre>"
        self.text.setHtml(body)
        self.apply_value.setChecked(False)
        self.apply_value.setEnabled(value is not None and entry['scheme_key'] is not None and entry['method'] == 'cgmlst')
        self.threshold.setValue(value or 0)
        self.bound_scheme.clear()
        for check in (self.schema_review, self.protocol_review, self.epi_review):
            check.setChecked(False)

    def search_newer(self):
        query = self.organism.currentText() + ' (cgMLST OR "single nucleotide polymorphism") (outbreak OR surveillance)'
        QDesktopServices.openUrl(QUrl('https://pubmed.ncbi.nlm.nih.gov/?' + urlencode({'term': query, 'sort': 'date'})))

    @staticmethod
    def open_source(url):
        if url.scheme() == 'https' and url.host() in {'journals.asm.org', 'www.nature.com'}:
            QDesktopServices.openUrl(url)

    def accept(self):
        item = self.entries.currentItem()
        if item is None:
            self.feedback.setText('No curated source selected. No cutoff or citation has been saved.')
            return
        try:
            self.evidence = record_decision(item.data(Qt.ItemDataRole.UserRole)['id'], self.context(),
                selected_threshold=self.threshold.value() if self.apply_value.isChecked() else None,
                justification=self.justification.toPlainText(), protocol_reviewed=self.protocol_review.isChecked(),
                schema_reviewed=self.schema_review.isChecked(), epi_reviewed=self.epi_review.isChecked())
        except ValueError as error:
            self.feedback.setText(str(error))
            return
        super().accept()
