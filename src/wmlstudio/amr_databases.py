"""User-initiated AMR reference management; downloads never run with an analysis.

This is the surface a person reaches from the Updates menu, and the one they are
sent to when a HYDRA run is refused. It answers three questions before anything
else: what reference data is installed, which release is it and how old is that,
and what can this screen not look for until something else is installed. A store
that exists but is empty is the failure that made HYDRA look broken, so an empty
store is stated in words here rather than shown as a table with no rows.
"""

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

BOUNDARY = ("Gene and mutation evidence is not a measured susceptibility phenotype, and this is "
            "not AMRFinderPlus or Kleborate: it is a screen over NCBI's published reference data.")


class AMRDatabaseDialog(QDialog):
    snapshotInstalled = Signal(str)

    def __init__(self, root, parent=None, update_root=None, project=None):
        super().__init__(parent)
        self.setWindowTitle("HYDRA reference databases")
        self.resize(880, 600)
        self.worker = None
        self.close_pending = False
        self.project = project
        self.update_root = str(update_root or root)
        layout = QVBoxLayout(self)
        layout.addWidget(label("AMR reference databases", "title"))
        layout.addWidget(
            label(
                "HYDRA searches these reference sets; it downloads nothing while an analysis runs. "
                "A release only contains determinants named before its date, so updating changes "
                "what the next run can find and never what a recorded run reported. " + BOUNDARY,
                "muted",
                True,
            )
        )
        row = QHBoxLayout()
        self.path = QLineEdit(str(root))
        row.addWidget(self.path, 1)
        row.addWidget(button("Choose existing…", self.browse))
        layout.addLayout(row)
        self.summary = label("", "small", True)
        layout.addWidget(self.summary)
        self.table = make_table(["Reference set", "Release", "Installed", "What it is searched for"])
        layout.addWidget(self.table, 1)
        self.gaps = label("", "small", True)
        layout.addWidget(self.gaps)
        row = QHBoxLayout()
        self.providers = QLineEdit("ncbi, protein")
        self.providers.setPlaceholderText("Explicit provider names, comma-separated")
        row.addWidget(label("Install / update", "small"))
        row.addWidget(self.providers, 1)
        self.install = button("Update to the current NCBI release", self.download, True)
        row.addWidget(self.install)
        row.addWidget(button("Use selected snapshot", self.choose))
        layout.addLayout(row)
        self.status = label(
            "NCBI is the starter provider. Other HYDRA providers can be selected by their upstream "
            "names; review each provider's terms before sharing its data.",
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
        """Read the chosen store and say plainly what it holds and what it does not."""
        from wmlstudio.hydra_runtime import database_status

        status = database_status(self.path.text() or None)
        self.summary.setText(status["label"])
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(status["installed"]))
        for row, (name, entry) in enumerate(sorted(status["installed"].items())):
            values = [name, entry["version"] or "unrecorded",
                      (entry["installed"] or "unrecorded")[:10], entry["purpose"]]
            for column, value in enumerate(values):
                self.table.setItem(row, column, cell(value))
        self.table.setSortingEnabled(True)
        self.gaps.setText(self.gap_sentence(status))
        if status["error"]:
            self.status.setText(status["error"])

    @staticmethod
    def gap_sentence(status):
        """What a run against this store could not report, before anyone starts one."""
        if status["error"]:
            return ("This folder's database manifest cannot be read, so no run may use it until "
                    "another snapshot is selected or downloaded.")
        if not status["installed"]:
            return ("No reference set is installed here, so a HYDRA run would have nothing to "
                    "search and would report no determinant for every isolate. That is not a "
                    "negative result. Download the reference sets, or select a folder that "
                    "already holds them.")
        parts = [f"Missing: {entry['name']} — without it, {entry['purpose']} does not run and "
                 "nothing is reported for it." for entry in status["missing"]]
        if status["organisms"]:
            parts.append(f"{len(status['organisms'])} organisms are accepted for point-mutation "
                         f"catalogues, {len(status['point_mutation_organisms'])} of them with a "
                         "DNA catalogue; an isolate outside that list is screened for genes only.")
        if status["stale"]:
            parts.append(f"This release is {status['age_days']} days old. Determinants named "
                         "after it are not in it.")
        return " ".join(parts)

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
                + "\n\nThis may be large and take several minutes. Existing snapshots are retained, "
                "and analyses already recorded keep the reference snapshot they were run against. "
                "Your samples are not uploaded.",
            )
            != QMessageBox.StandardButton.Yes
        ):
            return
        from wmlstudio.provisioning import update_hydra_databases

        root, project = self.update_root, self.project
        self.worker = FunctionWorker(
            lambda cancelled, progress: update_hydra_databases(
                root, names, project=project, cancelled=cancelled, progress=progress
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
        # The owner of this dialog records the choice from snapshotInstalled, so the
        # sentence here describes the snapshot rather than claiming the selection.
        self.snapshotInstalled.emit(path)
        self.status.setText(
            (result.get("database_status") or {}).get("label", "New snapshot installed.")
            + " Analyses already recorded keep the reference snapshot they were run against."
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
