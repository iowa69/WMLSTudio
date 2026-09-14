"""Optional native FastQC review, bounded execution, and evidence attachment."""

from __future__ import annotations

import copy
import uuid
from pathlib import Path

from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout

from wmlstudio.fastqc import run_fastqc, runtime_capabilities
from wmlstudio.scheduler import plan_resources, resources_for_run, run_bounded
from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import label


def reads_for_sample(sample):
    metadata = sample.get("metadata") or {}
    attached = metadata.get("reads") or {}
    if attached.get("reads"):
        return [item["path"] for item in sorted(attached["reads"], key=lambda item: item["mate"])]
    workflow = metadata.get("workflow") or {}
    if workflow.get("read1_path") and workflow.get("read2_path"):
        return [workflow["read1_path"], workflow["read2_path"]]
    path = sample.get("input_path") or ""
    name = path.lower().removesuffix(".gz").removesuffix(".bz2")
    if (sample.get("result") or {}).get("kind") == "fastq" or name.endswith((".fq", ".fastq")):
        return [path]
    return []


class FastQCDialog(QDialog):
    def __init__(self, samples, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Review complete FastQC checks")
        self.resize(790, 530)
        self.samples = [sample for sample in samples if reads_for_sample(sample)]
        self.allocation = None
        layout = QVBoxLayout(self)
        layout.addWidget(label("FastQC · original read reports", "title"))
        layout.addWidget(label("Every read is analyzed by the original FastQC engine. No trimming, filtering, assembly replacement, or external data upload occurs. Module warnings and failures remain visible for review.", "muted", True))
        table = make_table(["Isolate", "Original read files"])
        table.setSortingEnabled(False)
        table.setRowCount(len(self.samples))
        for index, sample in enumerate(self.samples):
            table.setItem(index, 0, cell(sample["name"], sample["id"]))
            table.setItem(index, 1, cell("; ".join(Path(path).name for path in reads_for_sample(sample))))
        table.setSortingEnabled(True)
        layout.addWidget(table, 1)
        capabilities = runtime_capabilities()
        try:
            self.allocation = plan_resources(threads_per_sample=2, memory_gb=1)
            summary = f"Up to {self.allocation.max_parallel} isolates concurrently; 2 CPU threads and 1 GiB RAM reservation per isolate."
        except ValueError as exc:
            summary = str(exc)
        skipped = len(samples) - len(self.samples)
        self.feedback = label(capabilities["message"] + "\n" + summary +
                              (f"\n{skipped} selected isolates have no linked FASTQs and will be skipped." if skipped else ""), "small", True)
        layout.addWidget(self.feedback)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Run FastQC")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            capabilities["available"] and self.allocation is not None and bool(self.samples))
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)


def run_project_fastqc(project, samples, output_root, allocation, *, cancelled=None, progress=None, root=None):
    """Background worker API; metadata writes are serial, completed evidence survives cancellation."""
    snapshots = copy.deepcopy([sample for sample in samples if reads_for_sample(sample)])
    completed = []
    def analyse(sample, resources, stopped, report):
        return run_fastqc(reads_for_sample(sample), Path(output_root) / f"{sample['id']}-{uuid.uuid4().hex[:10]}",
                          threads=resources.threads_per_sample, memory_gb=resources.memory_gb,
                          root=root, cancelled=stopped, progress=report)
    def attach(sample, result):
        current = project.get_sample(sample["id"])
        if current["input_path"] != sample["input_path"] or reads_for_sample(current) != reads_for_sample(sample):
            raise ValueError("Isolate inputs or attached reads changed while FastQC was running; report not linked")
        with project.transaction():
            metadata = copy.deepcopy(current.get("metadata") or {})
            if metadata.get("fastqc"):
                project.record_history(sample["id"], "fastqc_superseded", {"evidence": metadata["fastqc"]})
            metadata["fastqc"] = result
            project.set_metadata(sample["id"], metadata)
            project.record_history(sample["id"], "fastqc_completed", {
                "engine": result["engine"], "version": result["version"],
                "input_sha256": [entry["sha256"] for entry in result["inputs"]]})
        completed.append({"sample_id": sample["id"], "result": result})
    run_bounded(snapshots, analyse, allocation, cancelled=cancelled, on_result=attach, progress=progress)
    return completed


def launch_fastqc(window, sample_ids=None, plan=None, on_complete=None):
    """UI entrypoint. Explicit plan/IDs can be used by a reviewed workflow chain."""
    if window.busy():
        return
    samples = window.project.samples()
    if sample_ids is None:
        from wmlstudio.cohort_picker import CohortPickerDialog
        picker = CohortPickerDialog(samples, window.project, "Choose isolates for FastQC", parent=window)
        if picker.exec() != QDialog.DialogCode.Accepted:
            return
        sample_ids = picker.selected_ids
    samples = [sample for sample in samples if sample["id"] in set(sample_ids)]
    if plan is None:
        review = FastQCDialog(samples, window)
        if review.exec() != QDialog.DialogCode.Accepted:
            return
        samples, allocation = review.samples, review.allocation
    else:
        try:
            allocation = resources_for_run(plan, memory_gb=1)
        except ValueError as exc:
            window.error(str(exc))
            return
    project = window.project
    output = window.project_path.with_suffix(".files") / "fastqc_reports"
    def operation(cancelled, progress):
        return run_project_fastqc(project, samples, output, allocation, cancelled=cancelled, progress=progress)
    def finish(results):
        window.notify(f"FastQC reports saved for {len(results)} isolates. Module flags are available in each original HTML report; reads were not trimmed.")
        if on_complete:
            on_complete(results)
    window.launch_task(operation, "fastqc", finish)
