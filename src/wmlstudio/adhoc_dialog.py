"""Explicit cohort and scientific thresholds for a local, non-public scheme."""

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
)

from wmlstudio.widgets import label


class AdhocSchemeDialog(QDialog):
    def __init__(self, samples, parent=None):
        super().__init__(parent)
        self.samples, self.plan = samples, {}
        self.setWindowTitle("Create a local cohort scheme")
        self.resize(690, 570)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Build an ad-hoc scheme", "title"))
        layout.addWidget(label(f"Use {len(samples)} selected assemblies to find unique, complete coding loci. The reference anchors the locus set. This is a local research scheme, not a public cgMLST nomenclature or registered ST system.", "muted", True))
        form = QFormLayout()
        self.name = QLineEdit()
        self.name.setPlaceholderText("For example: E. faecium ward cohort · 2026-09")
        form.addRow("Scheme name", self.name)
        self.organism = QLineEdit()
        self.organism.setPlaceholderText("Optional user-assigned organism label")
        form.addRow("Organism", self.organism)
        self.anchor = QComboBox()
        for sample in samples:
            self.anchor.addItem(sample["name"], sample["id"])
        form.addRow("Reference assembly", self.anchor)
        self.thresholds = {}
        for key, title, default in [("min_prevalence", "Minimum cohort prevalence", 1.0),
                                     ("min_identity", "Minimum nucleotide identity", 0.95),
                                     ("min_coverage", "Minimum reciprocal coverage", 0.98)]:
            control = QDoubleSpinBox()
            control.setRange(0.5, 1.0)
            control.setSingleStep(0.01)
            control.setValue(default)
            self.thresholds[key] = control
            form.addRow(title, control)
        self.threads = QSpinBox()
        self.threads.setRange(1, 64)
        self.threads.setValue(4)
        form.addRow("CPU threads", self.threads)
        layout.addLayout(form)
        layout.addWidget(label("Full-CDS prediction uses bacterial genetic code 11. Duplicate copies and cross-locus homologs are excluded at the stated thresholds. Quality and reference selection matter: cohort-dependent distances are not interchangeable with a validated public scheme.", "small", True))
        self.feedback = label("No genomes are uploaded. A new fingerprinted scheme is saved locally; source assemblies stay unchanged.", "small", True)
        layout.addWidget(self.feedback)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        controls.button(QDialogButtonBox.StandardButton.Ok).setText("Build local scheme")
        controls.accepted.connect(self.accept)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)

    def accept(self):
        if not self.name.text().strip():
            self.feedback.setText("Give the local scheme a descriptive name before starting.")
            return
        anchor = self.anchor.currentData()
        ordered = sorted(self.samples, key=lambda sample: sample["id"] != anchor)
        self.plan = {"assemblies": [sample["input_path"] for sample in ordered],
            "name": self.name.text().strip(), "organism": self.organism.text().strip(),
            "threads": self.threads.value(), **{key: control.value() for key, control in self.thresholds.items()}}
        super().accept()
