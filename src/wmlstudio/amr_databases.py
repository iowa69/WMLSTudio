"""User-initiated AMR reference management; downloads never run with an analysis."""

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QVBoxLayout,
)

from wmlstudio.background import FunctionWorker
from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import button, label


class AMRDatabaseDialog(QDialog):
    snapshotInstalled = Signal(str)

    def __init__(self, root, parent=None, update_root=None):
        super().__init__(parent)
        self.setWindowTitle("HYDRA reference databases")
        self.resize(820, 530)
        self.worker = None
        self.close_pending = False
        self.update_root = str(update_root or root)
        layout = QVBoxLayout(self)
        layout.addWidget(label("AMR reference snapshots", "title"))
        layout.addWidget(
            label(
                "Install a versioned local snapshot. Updates retain previous snapshots and never send sequence data. Review each provider's terms before sharing its data.",
                "muted",
                True,
            )
        )
        row = QHBoxLayout()
        self.path = QLineEdit(str(root))
        row.addWidget(self.path, 1)
        row.addWidget(button("Choose existing…", self.browse))
        layout.addLayout(row)
        self.table = make_table(["Database", "Version / date", "Source path"])
        layout.addWidget(self.table, 1)
        row = QHBoxLayout()
        self.providers = QLineEdit("ncbi")
        self.providers.setPlaceholderText("Explicit provider names, comma-separated")
        row.addWidget(label("Install / update", "small"))
        row.addWidget(self.providers, 1)
        self.install = button("Download snapshot…", self.download, True)
        row.addWidget(self.install)
        row.addWidget(button("Use selected snapshot", self.choose))
        layout.addLayout(row)
        self.status = label(
            "NCBI is the starter provider. Other HYDRA providers can be selected by their upstream names.",
            "small",
            True,
        )
        layout.addWidget(self.status)
        self.progress = QProgressBar()
        layout.addWidget(self.progress)
        layout.addWidget(button("Cancel download / close", self.reject))
        self.path.editingFinished.connect(self.refresh)
        self.refresh()

    def refresh(self):
        from wmlstudio.hydra_runtime import installed_databases

        try:
            entries = installed_databases(self.path.text())
            self.table.setSortingEnabled(False)
            self.table.setRowCount(len(entries))
            for row, (name, entry) in enumerate(entries.items()):
                values = [
                    name,
                    entry.get("version") or entry.get("downloaded_at") or "Recorded in manifest",
                    entry.get("path"),
                ]
                for column, value in enumerate(values):
                    self.table.setItem(row, column, cell(value))
            self.table.setSortingEnabled(True)
        except ValueError as exc:
            self.status.setText(str(exc))

    def browse(self):
        path = QFileDialog.getExistingDirectory(self, "Choose HYDRA snapshot")
        if path:
            self.path.setText(path)
            self.refresh()

    def choose(self):
        from wmlstudio.hydra_runtime import installed_databases

        try:
            if not installed_databases(self.path.text()):
                raise ValueError("This folder has no usable HYDRA database manifest.")
            self.snapshotInstalled.emit(str(Path(self.path.text()).resolve()))
            self.accept()
        except ValueError as exc:
            self.status.setText(str(exc))

    def download(self):
        if self.worker and self.worker.isRunning():
            return
        names = [value.strip() for value in self.providers.text().split(",") if value.strip()]
        if not names:
            self.status.setText("Specify at least one upstream database name.")
            return
        if (
            QMessageBox.question(
                self,
                "Download reference data",
                "Download these HYDRA providers?\n\n"
                + ", ".join(names)
                + "\n\nThis may be large and take several minutes. Existing snapshots are retained. Your samples are not uploaded.",
            )
            != QMessageBox.StandardButton.Yes
        ):
            return
        from wmlstudio.hydra_runtime import update_databases

        root = self.update_root
        self.worker = FunctionWorker(
            lambda cancelled, progress: update_databases(
                root, names, cancelled=cancelled, progress=progress
            ),
            self,
        )
        self.worker.completed.connect(self.installed)
        self.worker.failed.connect(self.status.setText)
        self.worker.progress.connect(
            lambda percent, message: (self.progress.setValue(percent), self.status.setText(message))
        )
        self.worker.finished.connect(self.finished_task)
        self.install.setEnabled(False)
        self.worker.start()

    def installed(self, result):
        path = result.get("database_root") or result["reference_snapshot"]["root"]
        self.path.setText(path)
        self.refresh()
        self.snapshotInstalled.emit(path)
        self.status.setText(
            "New snapshot installed. Existing analyses keep their original database provenance."
        )

    def finished_task(self):
        self.install.setEnabled(True)
        if self.close_pending:
            super().reject()

    def reject(self):
        if self.worker and self.worker.isRunning():
            self.close_pending = True
            self.worker.cancel()
            self.status.setText("Cancelling after the current network read…")
        else:
            super().reject()

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.reject()
            event.ignore()
        else:
            event.accept()
