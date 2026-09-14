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
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.scheduler import GIB, detect_hardware, plan_resources
from wmlstudio.ui_common import cell, make_table, organism_for
from wmlstudio.widgets import label


class RunPlanDialog(QDialog):
    manageDatabases = Signal()

    def __init__(self, samples, scheme=None, db_root=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Review analysis plan")
        self.resize(850, 670)
        self.plan = {}
        self.hardware = detect_hardware()
        outer = QVBoxLayout(self)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)
        scroll.setWidget(content)
        outer.addWidget(scroll)
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
        from wmlstudio.fastqc import runtime_capabilities as fastqc_capabilities
        fastqc = fastqc_capabilities()
        self.fastqc = QCheckBox("Generate complete FastQC read reports (no trimming)")
        self.fastqc.setEnabled(fastqc["available"])
        self.fastqc.setToolTip(fastqc["message"])
        layout.addWidget(self.fastqc)
        layout.addWidget(
            label(
                "Assemblies follow the assigned MLST workflow; additional cgMLST/wgMLST schemes are selected in Compare. FASTQ inputs receive quality checks unless assembled first.",
                "small",
                True,
            )
        )
        self.hydra = QCheckBox("Also run HYDRA AMR analysis on the selected assemblies")
        layout.addWidget(self.hydra)
        self.resource_summary = label("", "small", True)
        layout.addWidget(self.resource_summary)
        self.advanced_button = QPushButton("Advanced settings…")
        self.advanced_button.setCheckable(True)
        layout.addWidget(self.advanced_button)
        self.advanced = QWidget()
        form = QFormLayout(self.advanced)
        self.advanced.setVisible(False)
        self.advanced_button.toggled.connect(self.advanced.setVisible)
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
        self.threads.setRange(1, self.hardware.cpus)
        self.threads.setValue(min(4, self.hardware.cpus))
        form.addRow("Threads per sample", self.threads)
        self.memory = QSpinBox()
        self.memory.setRange(3, 512)
        self.memory.setValue(3)
        self.memory.setSuffix(" GB")
        form.addRow("RAM reservation per sample", self.memory)
        self.resource_policy = QComboBox()
        for title, value in [("Balanced", "balanced"), ("Fast · smaller RAM reserve", "fast"),
                             ("Low memory · one sample at a time", "low_memory")]:
            self.resource_policy.addItem(title, value)
        form.addRow("Scheduling policy", self.resource_policy)
        resource_row = QHBoxLayout()
        self.cpu_budget = QSpinBox()
        self.cpu_budget.setRange(1, self.hardware.cpus)
        self.cpu_budget.setValue(self.hardware.cpus)
        resource_row.addWidget(label("Total CPU budget", "small"))
        resource_row.addWidget(self.cpu_budget)
        self.parallel = QSpinBox()
        self.parallel.setRange(0, self.hardware.cpus)
        self.parallel.setSpecialValueText("Automatic")
        resource_row.addWidget(label("Concurrent samples", "small"))
        resource_row.addWidget(self.parallel)
        form.addRow("Queue limits", resource_row)
        layout.addWidget(self.advanced)
        self.resource_feedback = label("", "small", True)
        form.addRow(self.resource_feedback)
        for field in (self.threads, self.memory, self.cpu_budget, self.parallel):
            field.valueChanged.connect(self.refresh_resources)
        self.resource_policy.currentIndexChanged.connect(self.refresh_resources)
        self.refresh_resources()
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
        outer.addWidget(actions)

    def browse_database(self):
        path = QFileDialog.getExistingDirectory(self, "Select HYDRA reference snapshot")
        if path:
            self.database.setText(path)

    def calculate_resources(self):
        return plan_resources(threads_per_sample=self.threads.value(), memory_gb=self.memory.value(),
                              cpu_budget=self.cpu_budget.value(), max_parallel=self.parallel.value() or None,
                              policy=self.resource_policy.currentData(), hardware=self.hardware)

    def refresh_resources(self):
        try:
            plan = self.calculate_resources()
            self.resource_summary.setText(
                f"Automatic queue: up to {plan.max_parallel} samples × {plan.threads_per_sample} threads. "
                "CPU and available memory are checked again when work starts.")
            available = f"{self.hardware.available_memory / GIB:.1f} GiB available" if self.hardware.available_memory is not None else "RAM unavailable: single-job fallback"
            self.resource_feedback.setText(
                f"Detected {self.hardware.cpus} available CPUs · {available}. "
                f"Up to {plan.max_parallel} samples × {plan.threads_per_sample} threads; "
                f"{plan.memory_gb} GiB reserved per sample, {plan.reserve_gb:.1f} GiB left for other work. "
                "RAM is an admission estimate, not an OS-enforced process limit. Cancellation remains available.")
        except ValueError as exc:
            self.resource_summary.setText(str(exc))
            self.resource_feedback.setText(str(exc))

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
        try:
            allocation = self.calculate_resources()
        except ValueError as exc:
            self.feedback.setText(str(exc))
            return
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
            "fastqc": self.fastqc.isChecked(),
            "hydra": self.hydra.isChecked(),
            "assemble": self.assemble.isChecked(),
            "memory_gb": self.memory.value(),
            "db_root": self.database.text(),
            "databases": [choice] if choice else None,
            "protein": self.protein.isChecked(),
            "point_mutations": self.point_mutations.isChecked(),
            "thresholds": {key: field.value() for key, field in self.threshold_controls.items()},
            "threads": self.threads.value(),
            "resource_plan": allocation.to_dict(),
        }
        super().accept()
