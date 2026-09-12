"""Native assignment forms; filesystem work is left to the caller's worker."""

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)


class ImportSamplesDialog(QDialog):
    """Return per-file assignments and copy options after the user accepts.

    ``scheme_entries`` is a list of (display name, path) tuples or paths.
    ``assignments`` and ``options`` contain plain Python data, safe to give to a
    worker. This dialog does not read sequences, load schemes, or copy files.
    """

    def __init__(self, paths, scheme_entries=(), parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import samples and assign organisms")
        self.resize(900, 680)
        self.assignments = [{"path": str(path), "typing_mode": "auto", "genus": "",
                             "species": "", "scheme_path": None} for path in paths]
        self.options = {"managed": True, "storage_root": "", "append_st": False}
        layout = QVBoxLayout(self)
        description = QLabel(
            "Assign an organism to each sample, or identify it from your local schemes. "
            "Select several rows to apply the same choice together."
        )
        description.setWordWrap(True)
        layout.addWidget(description)
        self.table = QTableWidget(len(self.assignments), 4)
        self.table.setHorizontalHeaderLabels(["Sample / file", "Assignment", "Organism", "Scheme"])
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table, 1)
        form = QFormLayout()
        self.mode = QComboBox()
        self.mode.addItem("Identify from local schemes", "auto")
        self.mode.addItem("Assign an organism", "manual")
        self.mode.addItem("Unknown organism / AMR only", "unknown")
        form.addRow("Organism choice", self.mode)
        self.genus = QLineEdit()
        self.genus.setPlaceholderText("For example, Klebsiella")
        self.species = QLineEdit()
        self.species.setPlaceholderText("For example, pneumoniae; optional if genus only")
        form.addRow("Genus", self.genus)
        form.addRow("Species", self.species)
        self.scheme = QComboBox()
        self.scheme.addItem("Choose from organism / local detection", None)
        for entry in scheme_entries:
            if isinstance(entry, (tuple, list)) and len(entry) == 2:
                name, path = entry
            else:
                path = Path(entry)
                name = path.name
            self.scheme.addItem(str(name), str(path))
        form.addRow("Typing scheme", self.scheme)
        layout.addLayout(form)
        buttons = QHBoxLayout()
        apply_selected = QPushButton("Apply to selected rows")
        apply_selected.clicked.connect(self.apply_selected)
        apply_all = QPushButton("Apply to all rows")
        apply_all.clicked.connect(self.apply_all)
        buttons.addWidget(apply_selected)
        buttons.addWidget(apply_all)
        buttons.addStretch()
        layout.addLayout(buttons)
        self.managed = QCheckBox("Keep managed copies organised by genus, species and ST")
        self.managed.setChecked(True)
        layout.addWidget(self.managed)
        storage = QHBoxLayout()
        self.storage_root = QLineEdit()
        self.storage_root.setPlaceholderText("Default: Sequences folder beside the project")
        storage.addWidget(self.storage_root, 1)
        browse = QPushButton("Choose storage folder…")
        browse.clicked.connect(self.browse_storage)
        storage.addWidget(browse)
        layout.addLayout(storage)
        self.append_st = QCheckBox("Append the assigned ST to managed filenames")
        layout.addWidget(self.append_st)
        self.feedback = QLabel("Original sequence files stay unchanged.")
        self.feedback.setWordWrap(True)
        layout.addWidget(self.feedback)
        actions = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        actions.accepted.connect(self.accept)
        actions.rejected.connect(self.reject)
        layout.addWidget(actions)
        self.mode.currentIndexChanged.connect(self.update_form)
        self.update_form()
        self.refresh_table()

    def update_form(self):
        manual = self.mode.currentData() == "manual"
        self.genus.setEnabled(manual)
        self.species.setEnabled(manual)
        self.scheme.setEnabled(self.mode.currentData() != "unknown")

    def browse_storage(self):
        path = QFileDialog.getExistingDirectory(self, "Choose managed sequence storage")
        if path:
            self.storage_root.setText(path)

    def apply_selected(self):
        rows = sorted({index.row() for index in self.table.selectionModel().selectedRows()})
        if not rows:
            self.feedback.setText("Select one or more sample rows first.")
            return
        self._apply(rows)

    def apply_all(self):
        self._apply(range(len(self.assignments)))

    def _apply(self, rows):
        mode = self.mode.currentData()
        genus = self.genus.text().strip() if mode == "manual" else ""
        species = self.species.text().strip() if mode == "manual" else ""
        if mode == "manual" and not genus:
            self.feedback.setText("Enter a genus for manual organism assignment.")
            return
        for row in rows:
            self.assignments[row].update({"typing_mode": mode, "genus": genus, "species": species,
                                          "scheme_path": self.scheme.currentData() if mode != "unknown" else None})
        self.feedback.setText("Assignment applied. Original sequence files stay unchanged.")
        self.refresh_table()

    def refresh_table(self):
        labels = {"auto": "Identify locally", "manual": "Assigned", "unknown": "Unknown / AMR only"}
        for row, assignment in enumerate(self.assignments):
            cells = [assignment.get("name") or Path(assignment["path"]).name,
                     labels.get(assignment["typing_mode"], assignment["typing_mode"]),
                     " ".join(filter(None, [assignment.get("genus"), assignment.get("species")])),
                     Path(assignment["scheme_path"]).name if assignment.get("scheme_path") else "Automatic"]
            for column, value in enumerate(cells):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, row)
                item.setToolTip(assignment["path"] if column == 0 else value)
                self.table.setItem(row, column, item)

    def accept(self):
        if not self.assignments:
            self.feedback.setText("No sequence files were selected.")
            return
        self.options = {"managed": self.managed.isChecked(), "storage_root": self.storage_root.text().strip(),
                        "append_st": self.append_st.isChecked()}
        super().accept()


class BatchAssignmentDialog(ImportSamplesDialog):
    def __init__(self, samples, scheme_entries=(), parent=None):
        samples = list(samples)
        super().__init__([sample.get("input_path", "") for sample in samples], scheme_entries, parent)
        self.setWindowTitle("Assign organisms to selected samples")
        self.managed.setChecked(False)
        for assignment, sample in zip(self.assignments, samples, strict=True):
            metadata = sample.get("metadata", {})
            workflow = metadata.get("workflow", {})
            organism = metadata.get("organism", {})
            assignment.update({"sample_id": sample["id"], "name": sample["name"],
                               "typing_mode": workflow.get("typing_mode", "auto"),
                               "genus": organism.get("genus", ""), "species": organism.get("species", ""),
                               "scheme_path": workflow.get("scheme_path")})
        self.refresh_table()
