"""One explicit, reusable native isolate selection dialog per analysis/report."""

from __future__ import annotations

from collections import defaultdict

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLineEdit,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
)

from wmlstudio.journey import isolate_records
from wmlstudio.ui_common import cell, make_table, organism_for
from wmlstudio.widgets import button, label


class CohortPickerDialog(QDialog):
    """Picking a cohort never modifies samples, profiles, or another view's scope."""

    def __init__(self, samples, project, title="Choose isolates", selected_ids=None,
                 parent=None, *, include_reads=True):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1060, 690)
        self.setMinimumSize(780, 480)
        self.samples = isolate_records(samples)
        if not include_reads:
            self.samples = [sample for sample in self.samples
                            if (sample.get("result") or {}).get("kind") != "fastq"
                            and not sample.get("input_path", "").lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fq", ".fastq"))]
        valid = {sample["id"] for sample in self.samples}
        self.selected_ids = set(selected_ids or ()) & valid
        self._folder = None
        self._filling = False
        self._profiles = {}
        for sample in self.samples:
            primary = sample.get("result") or {}
            secondary = project.analysis_summaries(sample["id"]) if project else []
            self._profiles[sample["id"]] = sorted({str(profile.get("scheme") or "unnamed")
                for profile in secondary if profile.get("scheme_digest") != primary.get("scheme_digest")
                and profile.get("input_sha256") and primary.get("input_sha256")
                and sample.get("status") == "completed"
                and primary["input_sha256"] == profile["input_sha256"]})
        layout = QVBoxLayout(self)
        layout.addWidget(label(title, "title", True))
        layout.addWidget(label("Choose an organism folder, then check all or individual isolates. This selection belongs only to the action you are starting.", "muted", True))
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find an isolate, date, organism, ST or saved cgMLST scheme…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_rows)
        layout.addWidget(self.search)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.folders = QTreeWidget()
        self.folders.setHeaderLabel("Organism folders")
        self.folders.setMinimumWidth(180)
        everything = QTreeWidgetItem([f"All isolates ({len(self.samples)})"])
        everything.setData(0, Qt.ItemDataRole.UserRole, None)
        self.folders.addTopLevelItem(everything)
        groups = defaultdict(lambda: defaultdict(int))
        for sample in self.samples:
            genus, species, _ = organism_for(sample)
            groups[genus or "Unknown"][species or "Unspecified"] += 1
        for genus in sorted(groups):
            parent_item = QTreeWidgetItem([f"{genus} ({sum(groups[genus].values())})"])
            parent_item.setData(0, Qt.ItemDataRole.UserRole, (genus,))
            everything.addChild(parent_item)
            for species, count in sorted(groups[genus].items()):
                item = QTreeWidgetItem([f"{species} ({count})"])
                item.setData(0, Qt.ItemDataRole.UserRole, (genus, species))
                parent_item.addChild(item)
        self.folders.expandAll()
        self.folders.itemClicked.connect(self.choose_folder)
        splitter.addWidget(self.folders)
        self.table = make_table(["Include", "Isolate", "Collection date", "Organism", "ST", "Additional saved profiles", "Quality / state"])
        self.table.setColumnWidth(0, 60)
        self.table.itemChanged.connect(self.item_changed)
        self.table.cellDoubleClicked.connect(self.toggle_row)
        splitter.addWidget(self.table)
        splitter.setSizes([225, 810])
        layout.addWidget(splitter, 1)
        controls = QHBoxLayout()
        controls.addWidget(button("Select shown", self.select_shown))
        controls.addWidget(button("Deselect shown", self.deselect_shown))
        controls.addWidget(button("Clear all", self.clear_all))
        self.count = label("", "small", True)
        controls.addWidget(self.count, 1)
        layout.addLayout(controls)
        self.feedback = label("Organism folders show assigned or provisional labels; inspect identity evidence before interpretation.", "small", True)
        layout.addWidget(self.feedback)
        self.controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.controls.button(QDialogButtonBox.StandardButton.Ok).setText("Use selected isolates")
        self.controls.accepted.connect(self.accept)
        self.controls.rejected.connect(self.reject)
        layout.addWidget(self.controls)
        self.refresh_rows()

    def choose_folder(self, item, column=0):
        self._folder = item.data(0, Qt.ItemDataRole.UserRole)
        self.refresh_rows()

    def refresh_rows(self):
        self._filling = True
        self.table.setSortingEnabled(False)
        self.table.blockSignals(True)
        query = self.search.text().strip().casefold()
        rows = []
        for sample in self.samples:
            genus, species, _ = organism_for(sample)
            taxon = (genus or "Unknown", species or "Unspecified")
            if self._folder and taxon[:len(self._folder)] != tuple(self._folder):
                continue
            metadata = sample.get("metadata") or {}
            result = sample.get("result") or {}
            date = (metadata.get("annotations") or {}).get("collection_date") or metadata.get("collection_date") or metadata.get("isolation_date") or metadata.get("date") or "Not recorded"
            values = ["", sample["name"], date, " ".join(taxon), result.get("st") if result.get("st") is not None else "Unassigned",
                      "; ".join(self._profiles[sample["id"]]) or "Not called", result.get("status") or sample.get("status")]
            if query and query not in " ".join(map(str, values)).casefold():
                continue
            rows.append((sample["id"], values))
        self.table.setRowCount(len(rows))
        for row, (identifier, values) in enumerate(rows):
            for column, value in enumerate(values):
                item = cell(value, identifier)
                if column == 0:
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked if identifier in self.selected_ids else Qt.CheckState.Unchecked)
                self.table.setItem(row, column, item)
        self.table.blockSignals(False)
        self.table.setSortingEnabled(True)
        self._filling = False
        self.update_count()

    def visible_ids(self):
        return {self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) for row in range(self.table.rowCount())}

    def update_count(self):
        hidden = len(self.selected_ids - self.visible_ids())
        self.count.setText(f"{len(self.selected_ids)} selected · {self.table.rowCount()} shown" +
                           (f" · {hidden} selected outside this filter" if hidden else ""))

    def item_changed(self, item):
        if self._filling or item.column() != 0:
            return
        identifier = item.data(Qt.ItemDataRole.UserRole)
        if item.checkState() == Qt.CheckState.Checked:
            self.selected_ids.add(identifier)
        else:
            self.selected_ids.discard(identifier)
        self.update_count()

    def toggle_row(self, row, column):
        if column != 0:
            item = self.table.item(row, 0)
            item.setCheckState(Qt.CheckState.Unchecked if item.checkState() == Qt.CheckState.Checked else Qt.CheckState.Checked)

    def select_shown(self):
        self.selected_ids.update(self.visible_ids())
        self.refresh_rows()

    def deselect_shown(self):
        self.selected_ids.difference_update(self.visible_ids())
        self.refresh_rows()

    def clear_all(self):
        self.selected_ids.clear()
        self.refresh_rows()

    def accept(self):
        if not self.selected_ids:
            self.feedback.setText("Choose at least one isolate, or cancel. No project records have been changed.")
            return
        super().accept()
