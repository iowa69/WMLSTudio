"""Review conservative read candidates or browse explicit files for each assembly."""

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
)

from wmlstudio.read_attachments import (
    attach_read_pair,
    suggest_read_attachments,
    validate_read_attachment,
)
from wmlstudio.scheduler import plan_resources, run_bounded
from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import label


class AttachReadsDialog(QDialog):
    def __init__(self, assemblies, reads=(), parent=None):
        super().__init__(parent)
        self.assemblies, self.assignments = list(assemblies), []
        reads = list(reads)
        self.setWindowTitle("Attach original reads to assemblies")
        self.resize(920, 590)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Link reads to existing isolates", "title"))
        layout.addWidget(label("Review filename suggestions or browse R1/R2 for any row. The assembly, existing typing and AMR evidence keep their stable IDs. Every paired record and both file hashes are checked in the background; no trimming occurs.", "muted", True))
        candidates = suggest_read_attachments(self.assemblies, reads)
        suggestions = {item["sample_id"]: item for item in candidates["suggestions"]}
        self.table = make_table(["Existing assembly", "Forward reads (R1)", "Reverse reads (R2)"])
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.assemblies))
        self.rows = []
        for index, sample in enumerate(self.assemblies):
            self.table.setItem(index, 0, cell(sample["name"], sample["id"]))
            combos = []
            for mate in (1, 2):
                combo = QComboBox()
                combo.addItem("Not attached", None)
                for read in reads:
                    combo.addItem(read["name"], {"path": read["input_path"], "sample_id": read["id"]})
                    if suggestions.get(sample["id"], {}).get(f"read{mate}_id") == read["id"]:
                        combo.setCurrentIndex(combo.count() - 1)
                self.table.setCellWidget(index, mate, combo)
                combos.append(combo)
            self.rows.append(combos)
        layout.addWidget(self.table, 1)
        browse_row = QHBoxLayout()
        for mate in (1, 2):
            button = QPushButton(f"Browse R{mate} for selected assembly…")
            button.clicked.connect(lambda checked=False, value=mate: self.browse(value))
            browse_row.addWidget(button)
        layout.addLayout(browse_row)
        self.feedback = label(f"{len(suggestions)} unique suggestions; {len(candidates['ambiguous'])} ambiguous rows withheld. Confirm that these reads belong to the biological isolate. Filename/read-ID agreement does not prove assembly provenance.", "small", True)
        layout.addWidget(self.feedback)
        actions = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        actions.button(QDialogButtonBox.StandardButton.Ok).setText("Validate and attach reads")
        actions.accepted.connect(self.accept)
        actions.rejected.connect(self.reject)
        layout.addWidget(actions)
        if self.assemblies:
            self.table.selectRow(0)

    def browse(self, mate):
        row = self.table.currentRow()
        if row < 0:
            self.feedback.setText("Select an assembly row before choosing its reads.")
            return
        path, _ = QFileDialog.getOpenFileName(self, f"Choose original R{mate} FASTQ", "", "FASTQ (*.fastq *.fq *.fastq.gz *.fq.gz *.fastq.bz2 *.fq.bz2);;All files (*)")
        if path:
            combo = self.rows[row][mate - 1]
            combo.addItem(Path(path).name, {"path": path, "sample_id": None})
            combo.setCurrentIndex(combo.count() - 1)

    def accept(self):
        assignments, used = [], set()
        for sample, combos in zip(self.assemblies, self.rows, strict=True):
            first, second = [combo.currentData() for combo in combos]
            if first is None and second is None:
                continue
            if first is None or second is None:
                self.feedback.setText(f"Choose both R1 and R2 for {sample['name']}, or leave both unattached.")
                return
            paths = [str(Path(item["path"]).expanduser().resolve()) for item in (first, second)]
            if len(set(paths)) != 2 or used.intersection(paths):
                self.feedback.setText("Each original read file may be assigned to only one isolate in this review.")
                return
            used.update(paths)
            assignments.append({"sample_id": sample["id"], "read1": paths[0], "read2": paths[1],
                                "read_sample_ids": [first["sample_id"], second["sample_id"]]})
        if not assignments:
            self.feedback.setText("Choose an R1/R2 pair for at least one assembly, or cancel.")
            return
        self.assignments = assignments
        super().accept()


def launch_read_attachment(window):
    if window.busy():
        return
    records = window.project.samples()
    def kind(record):
        known = (record.get("result") or {}).get("kind") or record.get("kind")
        name = Path(record.get("input_path", "")).name.lower().removesuffix(".gz").removesuffix(".bz2")
        return known or ("fastq" if name.endswith((".fq", ".fastq")) else "fasta" if name.endswith((".fa", ".fas", ".fasta", ".fna")) else None)
    assemblies = [record for record in records if record.get("input_path") and kind(record) == "fasta"]
    reads = [record for record in records if record.get("input_path") and kind(record) == "fastq"]
    if not assemblies:
        window.notify("Import or select an existing FASTA assembly before attaching its original reads.")
        return
    from wmlstudio.cohort_picker import CohortPickerDialog
    picker = CohortPickerDialog(assemblies, window.project, "Which assemblies need original reads?", parent=window, include_reads=False)
    if picker.exec() != QDialog.DialogCode.Accepted:
        return
    assemblies = [record for record in assemblies if record["id"] in picker.selected_ids]
    dialog = AttachReadsDialog(assemblies, reads, window)
    if dialog.exec() != QDialog.DialogCode.Accepted:
        return
    snapshots = {sample["id"]: sample for sample in assemblies}
    project = window.project
    try:
        allocation = plan_resources(threads_per_sample=1, memory_gb=1, max_parallel=2)
    except ValueError as exc:
        window.error(str(exc))
        return
    def operation(cancelled, progress):
        identifiers = []
        def validate(task, resources, stopped, report):
            return validate_read_attachment(snapshots[task["sample_id"]], task["read1"], task["read2"],
                read_sample_ids=task["read_sample_ids"], cancelled=stopped, progress=report)
        def attach(task, evidence):
            attach_read_pair(project, task["sample_id"], evidence)
            identifiers.append(task["sample_id"])
        run_bounded(dialog.assignments, validate, allocation, cancelled=cancelled, on_result=attach, progress=progress)
        return identifiers
    window.launch_task(operation, "read_attachment", lambda ids: window.notify(
        f"Attached fully validated original read pairs to {len(ids)} existing isolates. Assembly IDs and prior evidence are unchanged."))
