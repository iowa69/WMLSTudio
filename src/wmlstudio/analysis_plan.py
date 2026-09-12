"""Human-readable launch review: explicit analysis scope, resources and AMR options."""

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
)

from wmlstudio.ui_common import cell, make_table, organism_for
from wmlstudio.widgets import label


class RunPlanDialog(QDialog):
    manageDatabases = Signal()

    def __init__(self, samples, scheme=None, db_root=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Review analysis plan")
        self.resize(850, 670)
        self.plan = {}
        layout = QVBoxLayout(self)
        layout.addWidget(label(f"Run plan · {len(samples)} sample inputs", "title"))
        layout.addWidget(
            label(
                "Review the selected samples and optional analyses before starting. All computation stays on this computer.",
                "muted",
                True,
            )
        )
        table = make_table(["Sample", "Organism", "Workflow", "Scheme"])
        table.setSortingEnabled(False)
        table.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            genus, species, _ = organism_for(sample)
            workflow = sample.get("metadata", {}).get("workflow", {})
            values = [
                sample["name"],
                " ".join([genus, species]).strip() or "Unknown",
                "Scheme override" if scheme else workflow.get("typing_mode", "QC only"),
                Path(scheme or workflow.get("scheme_path") or "Automatic / not assigned").name,
            ]
            for column, value in enumerate(values):
                table.setItem(row, column, cell(value, sample["id"]))
        table.setSortingEnabled(True)
        layout.addWidget(table, 1)
        layout.addWidget(label("Typing / quality", "cardTitle"))
        self.assemble = QCheckBox("Assemble paired short reads with SKESA before typing")
        layout.addWidget(self.assemble)
        layout.addWidget(
            label(
                "Assemblies follow the assigned MLST workflow; additional cgMLST/wgMLST schemes are selected in Compare. FASTQ inputs receive quality checks unless assembled first.",
                "small",
                True,
            )
        )
        self.hydra = QCheckBox("Also run HYDRA AMR analysis on the selected assemblies")
        layout.addWidget(self.hydra)
        form = QFormLayout()
        self.database = QLineEdit(str(db_root or ""))
        database_row = QHBoxLayout()
        database_row.addWidget(self.database)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self.browse_database)
        database_row.addWidget(browse)
        manage = QPushButton("Install / update…")
        manage.clicked.connect(self.manageDatabases.emit)
        database_row.addWidget(manage)
        form.addRow("HYDRA reference snapshot", database_row)
        self.database_choice = QComboBox()
        self.database_choice.addItem("All installed nucleotide / protein databases", None)
        form.addRow("Database selection", self.database_choice)
        self.protein = QCheckBox("Include translated protein search")
        self.protein.setChecked(True)
        form.addRow("AMR methods", self.protein)
        self.point_mutations = QCheckBox("Organism-specific mutation catalog (confirmed labels only)")
        self.point_mutations.setChecked(True)
        form.addRow("Mutation evidence", self.point_mutations)
        self.threshold_controls = {}
        thresholds = QHBoxLayout()
        for key, title, default in [("min_identity", "Nucleotide identity", 80),
                                     ("min_coverage", "Nucleotide coverage", 60),
                                     ("protein_min_identity", "Protein identity", 90),
                                     ("protein_min_coverage", "Protein complete boundary", 90)]:
            field = QDoubleSpinBox()
            field.setRange(0, 100)
            field.setValue(default)
            field.setSuffix(" %")
            field.setToolTip(title + "; HYDRA's protein complete boundary does not exclude partial hits.")
            group = QVBoxLayout()
            group.addWidget(label(title, "small", True))
            group.addWidget(field)
            thresholds.addLayout(group)
            self.threshold_controls[key] = field
        form.addRow("AMR thresholds", thresholds)
        self.threads = QSpinBox()
        self.threads.setRange(1, 64)
        self.threads.setValue(4)
        form.addRow("CPU threads", self.threads)
        self.memory = QSpinBox()
        self.memory.setRange(3, 512)
        self.memory.setValue(8)
        self.memory.setSuffix(" GB")
        form.addRow("Assembly memory limit", self.memory)
        layout.addLayout(form)
        self.feedback = label(
            "No databases are downloaded during analysis. AMR detections are sequence evidence, not measured susceptibility.",
            "small",
            True,
        )
        layout.addWidget(self.feedback)
        self.database.textChanged.connect(self.refresh_databases)
        self.refresh_databases()
        actions = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        actions.button(QDialogButtonBox.StandardButton.Ok).setText("Start analysis")
        actions.accepted.connect(self.accept)
        actions.rejected.connect(self.reject)
        layout.addWidget(actions)

    def browse_database(self):
        path = QFileDialog.getExistingDirectory(self, "Select HYDRA reference snapshot")
        if path:
            self.database.setText(path)

    def refresh_databases(self):
        from wmlstudio.hydra_runtime import installed_databases

        selected = self.database_choice.currentData()
        self.database_choice.clear()
        self.database_choice.addItem("All installed nucleotide / protein databases", None)
        try:
            for name in installed_databases(self.database.text() or "."):
                self.database_choice.addItem(name, name)
            self.database_choice.setCurrentIndex(max(0, self.database_choice.findData(selected)))
        except ValueError as exc:
            self.feedback.setText(str(exc))

    def accept(self):
        if self.assemble.isChecked():
            from wmlstudio.assembly import resolve_skesa
            try:
                resolve_skesa()
            except RuntimeError as exc:
                self.feedback.setText(str(exc))
                return
        if self.hydra.isChecked():
            from wmlstudio.hydra_runtime import runtime_capabilities

            try:
                capabilities = runtime_capabilities(self.database.text())
                if not capabilities["available"]:
                    self.feedback.setText(capabilities["message"])
                    return
                if not capabilities["databases"]:
                    self.feedback.setText(
                        "Install or select a HYDRA database snapshot before requesting AMR analysis."
                    )
                    return
            except ValueError as exc:
                self.feedback.setText(str(exc))
                return
        choice = self.database_choice.currentData()
        self.plan = {
            "hydra": self.hydra.isChecked(),
            "assemble": self.assemble.isChecked(),
            "memory_gb": self.memory.value(),
            "db_root": self.database.text(),
            "databases": [choice] if choice else None,
            "protein": self.protein.isChecked(),
            "point_mutations": self.point_mutations.isChecked(),
            "thresholds": {key: field.value() for key, field in self.threshold_controls.items()},
            "threads": self.threads.value(),
        }
        super().accept()
