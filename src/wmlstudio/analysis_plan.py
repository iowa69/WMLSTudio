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

from wmlstudio.scheduler import GIB, detect_hardware
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
        # Two numbers, not five. A threads box, a RAM box, a policy list, a CPU
        # budget and a concurrency box overlapped and contradicted one another;
        # how much of a computer to use is "how many at once" and "how big each".
        from wmlstudio.scheduler import auto_worksize, thread_ceiling, worksize_constraint
        ceiling = thread_ceiling(self.hardware.cpus)
        suggested = auto_worksize(self.hardware)
        self.jobs = QSpinBox()
        self.jobs.setRange(1, ceiling)
        self.jobs.setValue(suggested.jobs)
        form.addRow("Samples at a time", self.jobs)
        self.threads = QSpinBox()
        self.threads.setRange(1, ceiling)
        self.threads.setValue(suggested.threads)
        form.addRow("CPU threads for each sample", self.threads)
        self.memory = QSpinBox()
        self.memory.setRange(3, 512)
        self.memory.setValue(3)
        self.memory.setSuffix(" GB")
        form.addRow("Memory for each sample", self.memory)
        form.addRow(label(worksize_constraint(self.hardware.cpus), "small", True))
        layout.addWidget(self.advanced)
        # resource_summary already exists above, outside the Advanced disclosure,
        # so the plain sentence about this computer is readable without opening it.
        self.resource_feedback = label("", "small", True)
        form.addRow(self.resource_feedback)
        for field in (self.jobs, self.threads, self.memory):
            field.valueChanged.connect(self.refresh_resources)
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
        from wmlstudio.scheduler import plan_for
        return plan_for(self.jobs.value(), self.threads.value(),
                        memory_gb=self.memory.value(), hardware=self.hardware)

    def refresh_resources(self):
        from wmlstudio.scheduler import clamp_worksize, describe_worksize
        try:
            plan = self.calculate_resources()
            size = clamp_worksize(self.jobs.value(), self.threads.value(), self.hardware.cpus)
            self.resource_summary.setText(describe_worksize(size))
            available = (f"{self.hardware.available_memory / GIB:.1f} GiB free"
                         if self.hardware.available_memory is not None
                         else "free memory unknown, so one sample at a time")
            trimmed = (" The request was reduced to fit this computer."
                       if getattr(plan, "reduced_for_memory", False) else "")
            self.resource_feedback.setText(
                f"{self.hardware.cpus} CPU threads · {available}. "
                f"{plan.memory_gb} GiB is set aside for each sample and "
                f"{plan.reserve_gb:.1f} GiB is left for everything else.{trimmed} "
                "That figure reserves a place in the queue rather than capping the program, "
                "and free memory is checked again when the run starts. You can always cancel.")
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
            # Refuse here, naming every missing piece, rather than starting a run
            # that dies partway with a message about a tool the user never chose.
            from wmlstudio import provisioning

            try:
                check = provisioning.hydra_prerequisites(selected=self.database.text() or None)
                if not check["ready"]:
                    self.feedback.setText(check["message"])
                    return
                if check.get("warnings"):
                    self.feedback.setText(check["database"]["label"] + " " + " ".join(check["warnings"]))
            except (OSError, ValueError) as exc:
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
