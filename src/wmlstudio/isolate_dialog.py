"""A complete isolate record reached by double-click, not a cramped side panel."""

import html
import json

from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QFileDialog, QTabWidget, QTextBrowser, QVBoxLayout,
)

from wmlstudio.ui_characterization import characterization_html
from wmlstudio.ui_common import FlowLayout
from wmlstudio.widgets import button, label


class IsolateRecordDialog(QDialog):
    def __init__(self, window, sample_id):
        super().__init__(window)
        self.window, self.sample_id = window, sample_id
        self.setWindowTitle("Complete isolate record")
        self.resize(1030, 750)
        self.setMinimumSize(760, 480)
        layout = QVBoxLayout(self)
        self.heading = label("", "title", True)
        layout.addWidget(self.heading)
        controls = FlowLayout()
        controls.addWidget(button("Edit organism / workflow…", self.edit_organism, True))
        controls.addWidget(button("Edit metadata…", self.edit_metadata))
        controls.addWidget(button("Analyse this isolate…", self.analyse))
        controls.addWidget(button("Export complete record…", self.export_record))
        layout.addLayout(controls)
        self.tabs = QTabWidget()
        self.views = {}
        for title in ("Overview / typing", "Identity / genes / plasmids", "Original reads", "History"):
            browser = QTextBrowser()
            browser.setOpenLinks(False)
            self.views[title] = browser
            self.tabs.addTab(browser, title)
        layout.addWidget(self.tabs, 1)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)
        self.refresh_record()

    def refresh_record(self):
        sample = self.window.project.get_sample(self.sample_id)
        self.heading.setText(sample["name"])
        escape = lambda value: html.escape(str(value))
        self.window.show_sample_detail_for_id(self.sample_id)
        self.views["Overview / typing"].setHtml(self.window.detail.toHtml())
        self.views["History"].setHtml(self.window.history_view.toHtml())
        self.views["Identity / genes / plasmids"].setHtml(characterization_html(sample))
        reads = (sample.get("metadata") or {}).get("reads") or {}
        entries = reads.get("reads") or []
        if entries:
            content = "<h2>Attached original reads</h2><p>Association is user-confirmed; read-ID agreement does not prove these reads generated the assembly.</p>"
            for entry in entries:
                content += (f"<h3>Mate {escape(entry.get('mate'))}</h3><p>{escape(entry.get('path'))}<br>"
                            f"SHA-256: {escape(entry.get('sha256'))}</p>")
            content += f"<p>Assembly SHA-256 at attachment: {escape(reads.get('assembly_sha256_at_link'))}</p>"
        else:
            content = "<h2>Original reads</h2><p>No explicit FASTQ attachment is recorded. Use Attach reads from the library to review and attach a pair.</p>"
            assembly = (sample.get("metadata") or {}).get("assembly") or {}
            if assembly:
                content += "<h3>Assembly provenance</h3><pre>" + escape(json.dumps(assembly, indent=2, ensure_ascii=False)) + "</pre>"
        self.views["Original reads"].setHtml(content)

    def edit_organism(self):
        if self.window.busy():
            return
        from wmlstudio.storage import assign_organism
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        sample = self.window.project.get_sample(self.sample_id)
        dialog = BatchAssignmentDialog([sample], self.window.scheme_entries(), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        for assignment in dialog.assignments:
            assign_organism(self.window.project, [self.sample_id], assignment["genus"], assignment["species"],
                            scheme_path=assignment.get("scheme_path"), typing_mode=assignment["typing_mode"])
        self.window.refresh()
        self.refresh_record()
        self.window.notify("Organism assignment saved across the project. Previous results are archived; review analysis to produce evidence for the changed workflow.")

    def edit_metadata(self):
        previous = self.window.selection_ids
        self.window.selection_ids = {self.sample_id}
        try:
            self.window.edit_metadata()
        finally:
            self.window.selection_ids = previous
        self.refresh_record()

    def analyse(self):
        self.accept()
        self.window.start_analysis(confirm=True, all_samples=True, sample_ids={self.sample_id})

    def export_record(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export complete isolate evidence", "isolate-evidence.json", "JSON (*.json)")
        if not path:
            return
        try:
            from wmlstudio.export import _atomic_text
            self.window.check_output(path)
            sample = self.window.project.get_sample(self.sample_id)
            payload = {"sample": sample, "additional_profiles": self.window.project.analysis_results(self.sample_id),
                       "history": self.window.project.history(self.sample_id)}
            with _atomic_text(path) as handle:
                json.dump(payload, handle, ensure_ascii=False, indent=2, allow_nan=False)
        except (OSError, ValueError, ImportError) as error:
            self.window.error(str(error))
