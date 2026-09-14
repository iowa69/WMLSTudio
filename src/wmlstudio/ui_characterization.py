"""Native, selected-isolate characterization review and linked evidence inspector."""

from __future__ import annotations

import html
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QScrollArea,
    QSplitter,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.journey import characterization_state, isolate_records
from wmlstudio.paths import resource_root
from wmlstudio.scheduler import plan_resources, run_bounded
from wmlstudio.ui_common import FlowLayout, cell, make_table
from wmlstudio.widgets import button, label


def characterization_html(sample):
    """Escaped, compact drill-down; no unverified imported values become HTML."""
    def escape(value):
        return html.escape(str(value))

    state, evidence = characterization_state(sample)
    body = f"<h2>{escape(sample['name'])}</h2><p><b>Characterization:</b> {escape(state)}</p>"
    if not evidence:
        return body + "<p>No characterization was run. An empty result is not evidence of absence.</p>"
    if state in {"stale", "unverified"}:
        body += "<p><b>Review input identity:</b> this stored evidence is not confirmed current for this isolate.</p>"
    species = evidence.get("species_evidence") or {}
    body += f"<h3>Independent species evidence · {escape(species.get('status', 'not_run'))}</h3>"
    name = " ".join(str(species.get(key) or "") for key in ("genus", "species")).strip()
    if name:
        body += f"<p><b>Reference-supported candidate:</b> {escape(name)}</p>"
    body += f"<p>{escape(species.get('reason', ''))}<br>Subspecies: {escape(species.get('subspecies_status', 'unresolved'))}</p>"
    nearest = species.get("nearest") or {}
    if nearest:
        body += "<p>" + " · ".join(f"<b>{title}:</b> {escape(nearest.get(key, 'not available'))}"
                                    for title, key in (("ANI %", "ani"), ("Query aligned fraction", "query_fraction"),
                                                       ("Reference aligned fraction", "reference_fraction"), ("Reference", "reference_id"))) + "</p>"
    virulence = evidence.get("virulence") or {}
    body += f"<h3>Defined virulence panel · {escape(virulence.get('status', 'not_run'))}</h3>"
    groups = virulence.get("loci") or virulence.get("groups") or []
    if groups:
        body += "<table cellpadding='5'><tr><th>Locus</th><th>Evidence</th><th>Genes detected</th><th>Intact CDS</th></tr>"
        for group in groups:
            body += (f"<tr><td>{escape(group.get('locus'))}</td><td>{escape(group.get('status'))}</td>"
                     f"<td>{escape(group.get('genes_detected'))}/{escape(group.get('genes_total'))}</td>"
                     f"<td>{escape(group.get('intact_genes_detected'))}</td></tr>")
        body += "</table><p>Defined-locus screening, not all virulence factors or a virulence phenotype.</p>"
    elif virulence.get("reason"):
        body += f"<p>{escape(virulence['reason'])}</p>"
    drugs = evidence.get("drug_associations") or {}
    body += f"<h3>Reference-backed drug associations · {escape(drugs.get('status', 'not_run'))}</h3>"
    if drugs.get("associations"):
        body += "<table cellpadding='5'><tr><th>Determinant</th><th>Class</th><th>Reference subclass</th><th>Method</th></tr>"
        for entry in drugs["associations"]:
            body += "<tr>" + "".join(f"<td>{escape(entry.get(key) or '—')}</td>" for key in ("gene", "class", "subclass", "method")) + "</tr>"
        body += "</table>"
    body += "<p>Not an AST result. No determinant detected does not mean susceptible.</p>"
    plasmids = evidence.get("plasmid_hypotheses") or {}
    body += f"<h3>Plasmid hypotheses · {escape(plasmids.get('status', 'not_run'))}</h3>"
    body += f"<p>{escape(plasmids.get('reason', ''))}</p>"
    for association in plasmids.get("contig_associations", []):
        body += (f"<p><b>{escape(association.get('replicon'))}</b> + {escape(association.get('marker'))}"
                 f" on {escape(association.get('contig'))}; {escape(association.get('gap_bp'))} bp apart. "
                 "Same-contig evidence only; complete plasmid identity and transmission are unproven.</p>")
    body += f"<h3>Provenance</h3><p>Assembly SHA-256: {escape(evidence.get('input_sha256'))}</p>"
    for section in (species, virulence):
        digest = section.get("reference_digest") or section.get("provenance", {}).get("reference_digest")
        if digest:
            body += f"<p>Reference fingerprint: {escape(digest)}</p>"
    for note in evidence.get("limitations", []):
        body += f"<p>{escape(note)}</p>"
    return body


class CharacterizationPlanDialog(QDialog):
    installRequested = Signal()

    def __init__(self, samples, reference_path=None, parent=None):
        super().__init__(parent)
        self.plan = None
        self.setWindowTitle("Review identity, virulence and accessory evidence")
        self.resize(920, 740)
        outer = QVBoxLayout(self)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QVBoxLayout(panel)
        scroll.setWidget(panel)
        outer.addWidget(scroll, 1)
        layout.addWidget(label(f"Characterize {len(samples)} selected assemblies", "title", True))
        layout.addWidget(label("Review the exact cohort and assays. Classical MLST and additional profiles are preserved. Downloads and analyses require explicit confirmation.", "muted", True))
        self.table = make_table(["Include", "Isolate", "Existing characterization"])
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(samples))
        self.table.setMinimumHeight(160)
        for row, sample in enumerate(samples):
            item = cell("", sample["id"])
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked)
            self.table.setItem(row, 0, item)
            self.table.setItem(row, 1, cell(sample["name"]))
            self.table.setItem(row, 2, cell(characterization_state(sample)[0]))
        self.table.setColumnWidth(0, 65)
        layout.addWidget(self.table)
        self.species = QCheckBox("Independent species evidence: genome ANI, aligned fractions and alternatives")
        self.species.setChecked(True)
        self.virulence = QCheckBox("Defined KpSC-associated virulence loci: ybt, clb, iuc, iro, rmp, rmpA2")
        self.virulence.setChecked(True)
        self.hydra = QCheckBox("Run HYDRA AMR if no input-verified report is available; reuse verified existing evidence")
        self.hydra.setChecked(True)
        for control in (self.species, self.virulence, self.hydra):
            layout.addWidget(control)
        layout.addWidget(label("Species/virulence reference panel: focused Klebsiella representatives and E. coli outgroup; not a comprehensive taxonomy database. Existing HYDRA plasmid assays may support co-location hypotheses; the AMR starter alone does not assay plasmids.", "small", True))
        form = QFormLayout()
        self.reference = QLineEdit(str(reference_path or ""))
        reference_row = QHBoxLayout()
        reference_row.addWidget(self.reference)
        reference_row.addWidget(button("Browse…", self.browse_reference))
        reference_row.addWidget(button("Install / update…", self.request_install))
        form.addRow("Characterization snapshot", reference_row)
        self.policy = QComboBox()
        for title, value in (("Balanced", "balanced"), ("Faster throughput", "fast"), ("Low memory", "low_memory")):
            self.policy.addItem(title, value)
        form.addRow("Automatic resource policy", self.policy)
        layout.addLayout(form)
        self.resource_text = label("", "small", True)
        layout.addWidget(self.resource_text)
        self.feedback = label("No susceptible/resistant phenotype or transmission verdict will be inferred.", "muted", True)
        layout.addWidget(self.feedback)
        self.policy.currentIndexChanged.connect(self.update_resources)
        self.update_resources()
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        controls.button(QDialogButtonBox.StandardButton.Ok).setText("Run reviewed characterization")
        controls.accepted.connect(self.accept)
        controls.rejected.connect(self.reject)
        outer.addWidget(controls)

    def update_resources(self):
        try:
            allocation = plan_resources(threads_per_sample=4, memory_gb=3, policy=self.policy.currentData())
            self.resource_text.setText(f"Up to {allocation.max_parallel} simultaneous isolates × {allocation.threads_per_sample} threads; "
                                      f"{allocation.memory_gb} GiB reserved per job, {allocation.reserve_gb:g} GiB retained for the system. RAM is an estimate, not an enforced process cap.")
        except ValueError as error:
            self.resource_text.setText(str(error))

    def browse_reference(self):
        path = QFileDialog.getExistingDirectory(self, "Choose a verified characterization snapshot")
        if path:
            self.reference.setText(path)

    def request_install(self):
        self.reject()
        self.installRequested.emit()

    def accept(self):
        identifiers = [self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) for row in range(self.table.rowCount())
                       if self.table.item(row, 0).checkState() == Qt.CheckState.Checked]
        if not identifiers:
            self.feedback.setText("Include at least one assembly, or cancel.")
            return
        reference = self.reference.text().strip()
        if (self.species.isChecked() or self.virulence.isChecked()) and not (Path(reference) / "manifest.json").is_file():
            self.feedback.setText("Install or select a characterization snapshot. No reference download occurs during analysis.")
            return
        try:
            allocation = plan_resources(threads_per_sample=4, memory_gb=3, policy=self.policy.currentData())
        except ValueError as error:
            self.feedback.setText(str(error))
            return
        self.plan = {"sample_ids": identifiers, "reference_path": reference or None,
                     "species": self.species.isChecked(), "virulence": self.virulence.isChecked(),
                     "hydra": self.hydra.isChecked(), "resources": allocation.to_dict()}
        super().accept()


class CharacterizationWorkspaceMixin:
    def build_hydra(self):
        super().build_hydra()
        tabs = getattr(self, "evidence_tabs", None)
        if tabs is None:
            return
        panel = QWidget()
        layout = QVBoxLayout(panel)
        controls = FlowLayout()
        controls.addWidget(button("Choose isolates…", self.choose_feature_cohort, True))
        controls.addWidget(button("Run characterization…", self.run_characterization_selected))
        controls.addWidget(button("Explain these results", lambda: self.open_workflow_guide("resistance")))
        layout.addLayout(controls)
        self.characterization_table = make_table(["Sample", "Evidence state", "Species candidate", "ANI %", "Classical ST", "Virulence loci", "Drug classes", "Replicons"])
        self.characterization_table.itemSelectionChanged.connect(self.show_characterization_detail)
        self.characterization_table.cellDoubleClicked.connect(lambda row, column: self.open_isolate_record(self.characterization_table.item(row, 0).data(Qt.ItemDataRole.UserRole)))
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self.characterization_table)
        self.characterization_detail = QTextBrowser()
        self.characterization_detail.setOpenExternalLinks(False)
        splitter.addWidget(self.characterization_detail)
        splitter.setSizes([320, 230])
        layout.addWidget(splitter, 1)
        tabs.insertTab(0, panel, "Identity / virulence / plasmid hypotheses")

    def refresh_features(self):
        super().refresh_features()
        if not hasattr(self, "characterization_table"):
            return
        table = self.characterization_table
        selected = {item.data(Qt.ItemDataRole.UserRole) for item in table.selectedItems()}
        table.blockSignals(True)
        table.setSortingEnabled(False)
        samples = [sample for sample in isolate_records(self.project.samples()) if sample["id"] in getattr(self, "feature_ids", set())]
        table.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            state, evidence = characterization_state(sample)
            usable = evidence if state not in {"stale", "unverified"} else {}
            species = usable.get("species_evidence") or {}
            virulence = usable.get("virulence") or {}
            groups = virulence.get("loci") or virulence.get("groups") or []
            values = [sample["name"], state, " ".join(str(species.get(key) or "") for key in ("genus", "species")).strip() or "Unresolved",
                      (species.get("nearest") or {}).get("ani"), (sample.get("result") or {}).get("st"),
                      "; ".join(group["locus"] for group in groups if group.get("status") == "detected") or virulence.get("status", "not_run"),
                      "; ".join(sorted({entry.get("class") or entry.get("subclass") or "" for entry in usable.get("drug_associations", {}).get("associations", [])})) or usable.get("drug_associations", {}).get("status", "not_run"),
                      "; ".join(entry["gene"] for entry in usable.get("plasmid_hypotheses", {}).get("replicons", [])) or usable.get("plasmid_hypotheses", {}).get("status", "not_run")]
            for column, value in enumerate(values):
                item = cell(value, sample["id"])
                table.setItem(row, column, item)
                item.setSelected(sample["id"] in selected)
        table.setSortingEnabled(True)
        table.blockSignals(False)
        self.show_characterization_detail()

    def show_characterization_detail(self):
        if not hasattr(self, "characterization_table"):
            return
        item = self.characterization_table.item(self.characterization_table.currentRow(), 0)
        if item is None:
            self.characterization_detail.setHtml("<h3>Select an isolate to inspect the methods and evidence</h3><p>Not run is not negative. Input hashes link these results to the original isolate.</p>")
            return
        try:
            sample = self.project.get_sample(item.data(Qt.ItemDataRole.UserRole))
            self.characterization_detail.setHtml(characterization_html(sample))
        except KeyError:
            self.characterization_detail.clear()

    def active_characterization_reference(self):
        path = self.project.get_setting("characterization_reference", "")
        if path:
            return Path(path)
        bundled = resource_root() / "characterization" / "starter"
        return bundled if (bundled / "manifest.json").is_file() else None

    def install_characterization_references(self):
        if self.busy():
            return
        answer = QMessageBox.question(self, "Install immutable characterization references",
            "Download the pinned public Kleborate-derived reference panel?\n\n"
            "It covers defined Klebsiella representatives plus an E. coli outgroup and six defined virulence loci. "
            "This is not the complete Kleborate pipeline or a comprehensive species database. "
            "Source revision, hashes, accessions and GPL notices are preserved; no sample sequences are uploaded. "
            "Existing snapshots remain unchanged.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
        if answer != QMessageBox.StandardButton.Yes:
            return
        from wmlstudio.characterization_refs import provision_characterization_references
        root = self.root / "characterization-references"
        project = self.project
        def finished(result):
            project.set_setting("characterization_reference", result["path"])
            self.notify(f"Characterization snapshot installed: {result['species_count']} references. Review characterization to select and run assays.")
        self.launch_task(lambda cancelled, progress: provision_characterization_references(root, cancelled=cancelled, progress=progress), "characterization_references", finished)

    def run_characterization_selected(self):
        if self.busy():
            return
        records = isolate_records(self.project.samples())
        from wmlstudio.cohort_picker import CohortPickerDialog
        picker = CohortPickerDialog(records, self.project, "Which assemblies should be characterized?", parent=self, include_reads=False)
        if picker.exec() != QDialog.DialogCode.Accepted:
            return
        records = [sample for sample in records if sample["id"] in picker.selected_ids]
        assemblies = [sample for sample in records if sample.get("input_path") and not sample.get("missing_input")
                      and ((sample.get("result") or {}).get("kind") == "fasta" or
                           Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fa", ".fna", ".fasta", ".fas")))]
        if not assemblies:
            self.notify("Select/import FASTA assemblies for characterization. Review paired FASTQ assembly first if you only have reads.")
            self.navigate(1)
            return
        dialog = CharacterizationPlanDialog(assemblies, self.active_characterization_reference(), self)
        dialog.installRequested.connect(self.install_characterization_references)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        plan = dialog.plan
        self.feature_ids = set(plan["sample_ids"])
        samples = [sample for sample in assemblies if sample["id"] in plan["sample_ids"]]
        project = self.project
        database_root = self.active_amr_database()
        def operation(cancelled, progress):
            from wmlstudio.characterization import (
                characterize_assembly,
                hydra_report_for_record,
                persist_characterization,
            )
            from wmlstudio.hydra_runtime import run_assemblies
            from wmlstudio.sample_workflow import link_hydra
            completed, failures = [], []
            def characterize(sample, resources, stopped, report):
                hydra = hydra_report_for_record(sample)
                generated = None
                if plan["hydra"] and hydra is None:
                    generated = run_assemblies([sample["input_path"]], database_root,
                                              sample_names=[sample["id"]], threads=resources.threads_per_sample,
                                              organism=None, point_mutations=False, cancelled=stopped, progress=report)
                    hydra = generated
                result = characterize_assembly(sample["input_path"], plan["reference_path"], hydra_report=hydra,
                                               cancelled=stopped, progress=report, threads=resources.threads_per_sample,
                                               species=plan["species"], virulence=plan["virulence"])
                result["provenance"]["resource_plan"] = resources.to_dict()
                return result, generated
            def attach(sample, result):
                evidence, hydra = result
                persist_characterization(project, sample["id"], evidence)
                if hydra:
                    link_hydra(project, hydra, {sample["id"]: sample["id"]})
                completed.append(sample["id"])
            def failed(sample, error):
                failures.append({"sample_id": sample["id"], "message": str(error)})
                project.record_history(sample["id"], "characterization_failed", {"error": str(error)})
            run_bounded(samples, characterize, plan["resources"], cancelled=cancelled,
                        on_result=attach, on_error=failed, progress=progress)
            return {"completed": completed, "failures": failures}
        def finished(result):
            self.notify(f"Characterization saved for {len(result['completed'])} isolates; {len(result['failures'])} failed. Inspect individual assay states and provenance.")
            self.navigate(4)
            if hasattr(self, "evidence_tabs"):
                self.evidence_tabs.setCurrentIndex(0)
        self.launch_task(operation, "characterization", finished)
