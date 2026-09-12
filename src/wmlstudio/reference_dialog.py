"""Native, cancellable public reference catalog and snapshot download dialog."""

import threading
from pathlib import Path

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
)

from .reference_catalog import CGMLSTOrgCatalog, PasteurCatalog, PubMLSTCatalog
from .sequence import AnalysisCancelled


class _ReferenceWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int, str)

    def __init__(self, operation, parent=None):
        super().__init__(parent)
        self.operation = operation
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        try:
            result = self.operation(self.cancel_event.is_set, self.progress.emit)
            self.succeeded.emit(result)
        except AnalysisCancelled:
            self.progress.emit(0, 1, "Reference operation cancelled")
        except Exception as error:
            self.failed.emit(str(error))


class ReferenceManagerDialog(QDialog):
    schemeInstalled = Signal(str)
    installed = Signal(str)

    def __init__(self, root, parent=None, *, catalog=None):
        super().__init__(parent)
        self.root = Path(root)
        self.catalog = catalog or PubMLSTCatalog(self.root / "reference_cache")
        self.catalogs = {"PubMLST": self.catalog, "BIGSdb-Pasteur": PasteurCatalog(),
                         "cgMLST.org": CGMLSTOrgCatalog()}
        self.worker = None
        self.entries = []
        self.catalog_errors = []
        self._close_pending = False
        self.setWindowTitle("Online reference library")
        self.resize(980, 650)
        layout = QVBoxLayout(self)
        intro = QLabel("Find public MLST and cgMLST references by organism. Downloads create a new, validated snapshot in your local library.")
        intro.setWordWrap(True)
        layout.addWidget(intro)
        row = QHBoxLayout()
        self.provider = QComboBox()
        self.provider.addItems(list(self.catalogs))
        row.addWidget(self.provider)
        self.query = QLineEdit()
        self.query.setPlaceholderText("Organism, for example Staphylococcus epidermidis")
        self.query.returnPressed.connect(self.search)
        self.search_button = QPushButton("Search catalog")
        self.search_button.clicked.connect(self.search)
        row.addWidget(self.query, 1)
        row.addWidget(self.search_button)
        layout.addLayout(row)
        filters = QHBoxLayout()
        self.kind = QComboBox()
        self.kind.addItems(["All scheme types", "MLST", "cgMLST", "wgMLST", "Other"])
        self.min_loci, self.max_loci = QSpinBox(), QSpinBox()
        for spin in (self.min_loci, self.max_loci):
            spin.setRange(0, 100000)
        self.max_loci.setValue(100000)
        filters.addWidget(self.kind)
        filters.addWidget(QLabel("Loci from"))
        filters.addWidget(self.min_loci)
        filters.addWidget(QLabel("to"))
        filters.addWidget(self.max_loci)
        filters.addStretch()
        layout.addLayout(filters)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Organism", "Scheme", "Loci", "Type", "Metadata updated"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self.show_entry)
        layout.addWidget(self.table, 1)
        self.notice = QTextBrowser()
        self.notice.setMaximumHeight(140)
        self.notice.setPlainText("Your sample sequences are never sent to the reference service. Public access may exclude newer alleles; any access restriction will appear here and in the downloaded snapshot.")
        layout.addWidget(self.notice)
        self.status = QLabel("Ready")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        layout.addWidget(self.progress)
        buttons = QHBoxLayout()
        self.download_button = QPushButton("Download selected snapshot")
        self.download_button.setEnabled(False)
        self.download_button.clicked.connect(self.download)
        self.cancel_button = QPushButton("Cancel operation")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        buttons.addWidget(self.download_button)
        buttons.addStretch()
        buttons.addWidget(self.cancel_button)
        buttons.addWidget(close)
        layout.addLayout(buttons)

    def _start(self, operation, completed):
        if self.worker and self.worker.isRunning():
            return
        self.worker = _ReferenceWorker(operation, self)
        self.worker.succeeded.connect(completed)
        self.worker.failed.connect(lambda text: self.status.setText(f"Could not complete: {text}"))
        self.worker.progress.connect(self._progress)
        self.worker.finished.connect(self._finished)
        self.search_button.setEnabled(False)
        self.provider.setEnabled(False)
        self.download_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress.setRange(0, 0)
        self.worker.start()

    def _progress(self, current, total, text):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(current)
        self.status.setText(text)

    def _finished(self):
        self.search_button.setEnabled(True)
        self.provider.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        self.show_entry()
        if self._close_pending:
            super().reject()

    def search(self):
        query = self.query.text().strip()
        if not query:
            self.status.setText("Enter an organism name to search.")
            return
        if self.min_loci.value() > self.max_loci.value():
            self.status.setText("The minimum locus count must not exceed the maximum.")
            return
        filters = {"min_loci": self.min_loci.value(), "max_loci": self.max_loci.value(),
                   "scheme_type": self.kind.currentText() if self.kind.currentIndex() else None}
        self.catalog = self.catalogs[self.provider.currentText()]
        self.status.setText("Searching the public reference catalog…")
        self._start(lambda cancelled, progress: self.catalog.search_schemes(
            query, cancelled=cancelled, progress=progress, **filters), self._searched)

    def _searched(self, document):
        self.entries = document["schemes"]
        self.catalog_errors = document['errors']
        self.table.setRowCount(len(self.entries))
        for row, entry in enumerate(self.entries):
            values = [entry["organism"], entry["name"], entry["locus_count"], entry["type"], entry["last_updated"] or "Unknown"]
            for column, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                if column == 2:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(row, column, item)
        self.table.resizeColumnsToContents()
        self.status.setText(f"Found {len(self.entries)} schemes across {document['organisms_searched']} matching organism catalogs. {len(document['errors'])} catalog errors.")
        if document["errors"]:
            self.notice.setPlainText("Some catalogs could not be read:\n" + "\n".join(
                f"{item['organism']}: {item['error']}" for item in document["errors"]))
        if self.entries:
            self.table.selectRow(0)

    def show_entry(self):
        row = self.table.currentRow()
        selected = 0 <= row < len(self.entries)
        running = self.worker is not None and self.worker.isRunning()
        self.download_button.setEnabled(selected and not running)
        if selected:
            entry = self.entries[row]
            text = (f"{entry['organism']} · {entry['name']} · {entry['locus_count']} loci\n"
                    f"Source: {entry['url']}\n")
            text += entry["access_notice"] or "No anonymous-access restriction was reported in the scheme metadata."
            if entry.get('terms_notice'):
                text += '\n\n' + entry['terms_notice'] + '\n' + entry['terms_url']
            if self.catalog_errors:
                text += '\n\nCatalog errors:\n' + '\n'.join(
                    f"{item['organism']}: {item['error']}" for item in self.catalog_errors)
            self.notice.setPlainText(text)

    def download(self):
        row = self.table.currentRow()
        if not 0 <= row < len(self.entries):
            return
        entry = dict(self.entries[row])
        kwargs = {}
        if isinstance(self.catalog, CGMLSTOrgCatalog):
            choice = QMessageBox.question(self, 'Reference data usage restrictions',
                entry['terms_notice'] + '\n\n' + entry['terms_url']
                + '\n\nHave you reviewed the policy and confirmed your use is permitted?',
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if choice != QMessageBox.StandardButton.Yes:
                return
            kwargs['terms_acknowledged'] = True
        self._start(lambda cancelled, progress: self.catalog.download_scheme(
            entry, self.root / "schemes", cancelled=cancelled, progress=progress, **kwargs), self._downloaded)

    def _downloaded(self, result):
        self.status.setText(f"{'Installed' if result['created'] else 'Already installed'} {result['locus_count']} loci in a verified public snapshot.")
        self.notice.setPlainText("\n".join(result["notes"]) or "Snapshot validated.")
        self.schemeInstalled.emit(result["path"])
        self.installed.emit(result["path"])

    def cancel(self):
        if self.worker and self.worker.isRunning():
            self.status.setText("Cancelling after the current network read…")
            self.worker.cancel()

    def reject(self):
        if self.worker and self.worker.isRunning():
            self._close_pending = True
            self.cancel()
        else:
            super().reject()

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self._close_pending = True
            self.cancel()
            event.ignore()
        else:
            event.accept()
