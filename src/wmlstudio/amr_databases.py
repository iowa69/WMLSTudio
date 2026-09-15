"""User-initiated AMR reference management; downloads never run with an analysis.

This is the surface a person reaches from the Updates menu, and the one they are
sent to when a HYDRA run is refused. It answers three questions before anything
else: what reference data is installed, which release is it and how old is that,
and what can this screen not look for until something else is installed. A store
that exists but is empty is the failure that made HYDRA look broken, so an empty
store is stated in words here rather than shown as a table with no rows.

The list is the whole catalogue, not only what is installed. HYDRA is powerful
and silent about its inputs: it runs against whatever happens to be there and
reports nothing for everything else, and in a result those two are
indistinguishable. So every reference set the engine can search is named here
whether or not this computer has it, with what it is for, who publishes it, under
what licence and how large it is, and the user ticks the ones they want.
"""

from pathlib import Path

from PySide6.QtCore import Qt, Signal
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

COLUMNS = ["Reference set", "What it is searched for", "Size", "Published by", "Licence",
           "Release", "Installed"]


class AMRDatabaseDialog(QDialog):
    snapshotInstalled = Signal(str)

    def __init__(self, root, parent=None, update_root=None, project=None):
        super().__init__(parent)
        self.setWindowTitle("HYDRA reference databases")
        self.resize(1080, 660)
        self.worker = None
        self.close_pending = False
        self.project = project
        self.update_root = str(update_root or root)
        # The sets the user has ticked, kept across refreshes so re-reading the
        # store does not silently undo a choice somebody just made.
        self.chosen = set()
        self.catalogue = {"entries": []}
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
        self.table = make_table(COLUMNS)
        self.table.itemChanged.connect(self.remember)
        layout.addWidget(self.table, 1)
        self.gaps = label("", "small", True)
        layout.addWidget(self.gaps)
        self.organisms = label("", "small", True)
        layout.addWidget(self.organisms)
        row = QHBoxLayout()
        self.install = button("Install / update the ticked sets", self.download, True)
        row.addWidget(self.install)
        row.addWidget(button("Tick the core sets", self.choose_core))
        row.addWidget(button("Use selected snapshot", self.choose))
        row.addStretch()
        layout.addLayout(row)
        self.status = label(
            "The two core sets ship with this application, so a first run works with no download. "
            "Everything else is listed above and fetched only when you tick it; review each "
            "provider's own terms before installing or sharing its data.",
            "small",
            True,
        )
        layout.addWidget(self.status)
        self.progress = QProgressBar()
        layout.addWidget(self.progress)
        layout.addWidget(button("Cancel download / close", self.reject))
        self.path.editingFinished.connect(self.refresh)
        self.refresh()

    # --- reading the store --------------------------------------------------
    def refresh(self):
        """Read the chosen store and say plainly what it holds and what it does not."""
        from wmlstudio.hydra_runtime import database_catalogue, database_status

        chosen = self.path.text() or None
        status = database_status(chosen)
        self.catalogue = database_catalogue(chosen)
        if not self.chosen:
            self.chosen = set(self.catalogue["installed"]) or set(self.catalogue["core"])
        self.summary.setText(status["label"] + " " + self.catalogue["summary"])
        self.fill(self.catalogue["entries"])
        self.gaps.setText(self.gap_sentence(status))
        self.organisms.setText(self.organism_sentence(status))
        if status["error"] or self.catalogue["error"]:
            self.status.setText(status["error"] or self.catalogue["error"])

    def fill(self, entries):
        """One row per reference set, installed or not, with a tick where one is possible."""
        self.table.setSortingEnabled(False)
        self.table.blockSignals(True)
        self.table.setRowCount(len(entries))
        for row, entry in enumerate(entries):
            name = cell(entry["name"])
            name.setData(Qt.ItemDataRole.UserRole, entry["name"])
            manual = entry["download"] != "automatic"
            name.setFlags(name.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            if manual:
                name.setFlags(name.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
            name.setCheckState(Qt.CheckState.Checked
                               if entry["name"] in self.chosen and not manual
                               else Qt.CheckState.Unchecked)
            name.setToolTip("\n\n".join(filter(None, [
                entry["title"], entry["purpose"], entry["notes"], entry["citation"], entry["url"],
                entry["licence_note"],
                "This engine cannot fetch this set; it must be installed by hand."
                if manual else ""])))
            self.table.setItem(row, 0, name)
            values = [entry["purpose"],
                      entry["size"] or "not published; measured once installed",
                      entry["provider"] or "not recorded", entry["licence"],
                      entry["version"] or ("—" if not entry["installed"] else "unrecorded"),
                      self.installed_word(entry)]
            for column, value in enumerate(values, start=1):
                item = cell(value)
                item.setToolTip(entry["size_basis"] if column == 2 else str(value))
                self.table.setItem(row, column, item)
        self.table.blockSignals(False)
        self.table.setSortingEnabled(True)
        self.table.resizeRowsToContents()

    @staticmethod
    def installed_word(entry):
        """Installed, bundled, or the honest reason nothing is here for it."""
        if entry["bundled"]:
            return "Bundled with this release"
        if entry["installed"]:
            return "Installed" + (f" on {entry['staged'][:10]}" if entry["staged"] else "")
        return "Not installed — nothing is reported for it"

    def remember(self, item):
        """Keep a tick when the table is rebuilt, and never tick what cannot be fetched."""
        if item.column() != 0:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if not name:
            return
        if item.checkState() == Qt.CheckState.Checked:
            self.chosen.add(name)
        else:
            self.chosen.discard(name)

    def choose_core(self):
        self.chosen = set(self.catalogue["core"])
        self.fill(self.catalogue["entries"])

    @staticmethod
    def organism_sentence(status):
        """Which isolates get point mutations, and which get a gene screen only.

        Three different things are counted, because they answer three different
        questions: how many organisms the release will accept at all, how many
        have a curated protein mutation catalogue, and how many also have a DNA
        one. An isolate outside the last two is screened for genes, and that is
        not evidence that it carries no resistance mutation.
        """
        if not status["organisms"]:
            return ""
        return (f"{len(status['organisms'])} organisms are accepted by this release. "
                f"{len(status['protein_point_mutation_organisms'])} have a curated protein "
                f"point-mutation catalogue and {len(status['dna_point_mutation_organisms'])} also "
                "have a DNA catalogue, which is what 23S rRNA and the other non-coding targets are "
                "read from. An isolate whose genus and species fall outside those lists is "
                "screened for genes only; that is not evidence that it carries no mutation, and no "
                "other organism's catalogue is ever substituted for it.")

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

    # --- installing ---------------------------------------------------------
    def selected_names(self):
        """The ticked sets the engine can actually fetch, in a stable order."""
        fetchable = {entry["name"] for entry in self.catalogue["entries"]
                     if entry["download"] == "automatic"}
        return sorted(self.chosen & fetchable)

    def confirmation(self, names):
        """Name the sets, their sizes and whose terms apply, before anything is fetched."""
        rows = {entry["name"]: entry for entry in self.catalogue["entries"]}
        lines = [f"  • {name} — {rows[name]['provider'] or 'provider not recorded'}"
                 + (f", {rows[name]['size']}" if rows[name]["size"] else ", size not published")
                 + f", licence: {rows[name]['licence']}" for name in names]
        terms = sorted(name for name in names if rows[name]["licence_note"])
        return ("Download these HYDRA reference sets?\n\n" + "\n".join(lines)
                + ("\n\nThe provider's own terms apply to "
                   + ", ".join(terms)
                   + "; read them before installing or sharing that data." if terms else "")
                + "\n\nThis may be large and take several minutes. Existing snapshots are "
                  "retained, and analyses already recorded keep the reference snapshot they were "
                  "run against. Your samples are not uploaded.")

    def download(self):
        if self.worker and self.worker.isRunning():
            return
        names = self.selected_names()
        if not names:
            self.status.setText("Tick at least one reference set this engine can fetch. A set "
                                "marked as installed by hand is listed with where to get it.")
            return
        if (
            QMessageBox.question(self, "Download reference data", self.confirmation(names))
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
