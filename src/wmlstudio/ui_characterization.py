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
    QGroupBox,
    QHBoxLayout,
    QLineEdit,
    QMenu,
    QMessageBox,
    QScrollArea,
    QSplitter,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.context_menus import install_context_menu
from wmlstudio.journey import characterization_state, isolate_records
from wmlstudio.organism_modules import (
    BOUNDARY,
    cohort_applicability,
    cohort_sentence,
    record_skipped,
    registered_modules,
    selection_for_record,
)
from wmlstudio.paths import resource_root
from wmlstudio.scheduler import plan_resources, run_bounded
from wmlstudio.ui_common import FlowLayout, cell, make_table, organism_for
from wmlstudio.ui_workbench import CLASSICAL_KINDS
from wmlstudio.widgets import button, label

HIT_LIMIT = 40


def organism_line(sample):
    """'Klebsiella pneumoniae (Assigned)' — the source is never dropped."""
    genus, species, source = organism_for(sample)
    name = f"{genus} {species}".strip()
    return f"{name} ({source})" if name else "Unknown organism"


def classical_st(sample):
    """The classical ST, or a plain statement that this result is another quantity.

    A core-genome profile is not a sequence type: a cgST printed under a column
    headed 'Classical ST' is exactly the confusion this table must not create.
    """
    from wmlstudio.threshold_guidance import typing_scale
    result = sample.get("result") or {}
    count = len(result.get("alleles") or {}) or int(result.get("total_loci") or 0)
    scale = typing_scale(count)
    if scale["kind"] == "cgmlst":
        return f"Not a classical ST — {scale['label']}"
    return result.get("st") or ""


def hits_html(evidence, limit=HIT_LIMIT):
    """The individual reference matches behind a screen, with identity and coverage.

    A user who is asked to trust a call should be able to see what it was made
    from. Hits are shown highest identity first and the total is always stated,
    so a truncated table never looks like the whole evidence.
    """
    def escape(value):
        return html.escape(str(value if value is not None else "—"))

    def percent(hit, key):
        return f"{float(hit.get(key) or 0):.1f}"

    hits = (evidence or {}).get("hits") or []
    if not hits:
        return ""
    ordered = sorted(hits, key=lambda hit: (-float(hit.get("identity_pct") or 0),
                                            -float(hit.get("coverage_pct") or 0), str(hit.get("gene") or "")))
    body = ("<p><b>Individual reference matches</b></p><table cellpadding='5'><tr><th>Marker</th><th>Contig</th>"
            "<th>Identity %</th><th>Reference covered %</th><th>Position</th><th>Complete CDS</th></tr>")
    for hit in ordered[:limit]:
        qc = hit.get("cds_qc") or {}
        reasons = "; ".join(qc.get("reasons") or [])
        verdict = "yes" if qc.get("valid") else (f"no: {reasons}" if reasons else "not assessed")
        body += (f"<tr><td>{escape(hit.get('gene'))}</td><td>{escape(hit.get('contig'))}</td>"
                 f"<td>{escape(percent(hit, 'identity_pct'))}</td>"
                 f"<td>{escape(percent(hit, 'coverage_pct'))}</td>"
                 f"<td>{escape(hit.get('start'))}–{escape(hit.get('end'))} ({escape(hit.get('strand'))})</td>"
                 f"<td>{escape(verdict)}</td></tr>")
    body += "</table>"
    if len(ordered) > limit:
        body += f"<p>Showing the {limit} highest-identity matches of {len(ordered)}.</p>"
    return body + ("<p>A match is sequence similarity to a reference allele at this assay's threshold. It is not "
                   "proof that the gene is intact, expressed, or part of the element being typed.</p>")


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
    body += hits_html(virulence)
    # Each module renders its own table and repeats its own limitations; the shell
    # adds only the organism the tools were offered on, and the claim boundary.
    modules = registered_modules()
    if modules:
        body += (f"<h3>Organism-specific tools</h3><p>Offered for this isolate as <b>{escape(organism_line(sample))}</b>."
                 f" {escape(BOUNDARY)}</p>")
        for key, module in modules.items():
            block = evidence.get(key) or {}
            body += module.detail_html(block) + hits_html(block)
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

    def __init__(self, samples, reference_path=None, parent=None, *, project=None, scheme_entries=(),
                 database_root=None):
        super().__init__(parent)
        self.plan = None
        self.samples = list(samples)
        self.project = project
        self.database_root = str(database_root or "")
        self.hydra_check = None
        self.scheme_entries = list(scheme_entries)
        self.module_boxes = {}
        self.module_notes = {}
        self.module_touched = set()
        self.setWindowTitle("Review identity, virulence and accessory evidence")
        self.resize(920, 740)
        outer = QVBoxLayout(self)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QVBoxLayout(panel)
        scroll.setWidget(panel)
        outer.addWidget(scroll, 1)
        layout.addWidget(label(f"Characterize {len(self.samples)} selected assemblies", "title", True))
        layout.addWidget(label("Review the exact cohort and assays. Classical MLST and additional profiles are preserved. Downloads and analyses require explicit confirmation.", "muted", True))
        self.table = make_table(["Include", "Isolate", "Organism", "Existing characterization"])
        self.table.setSortingEnabled(False)
        self.table.setMinimumHeight(160)
        self.table.setColumnWidth(0, 65)
        self.refresh_plan_table()
        install_context_menu(self.table, "characterization.plan", self.show_plan_menu)
        layout.addWidget(self.table)
        layout.addWidget(label("Right-click a row to correct the organism. Which organism-specific tools are offered follows that assignment.", "small", True))
        self.species = QCheckBox("Independent species evidence: genome ANI, aligned fractions and alternatives")
        self.species.setChecked(True)
        self.virulence = QCheckBox("Defined KpSC-associated virulence loci: ybt, clb, iuc, iro, rmp, rmpA2")
        self.virulence.setChecked(True)
        self.hydra = QCheckBox("Run HYDRA AMR if no input-verified report is available; reuse verified existing evidence")
        self.hydra.setChecked(True)
        for control in (self.species, self.virulence, self.hydra):
            layout.addWidget(control)
        # What HYDRA is about to search, and what it cannot search for, stated
        # before the run rather than discovered in an empty result afterwards.
        self.hydra_state = label("", "small", True)
        layout.addWidget(self.hydra_state)
        self.refresh_hydra_state()
        layout.addWidget(self.build_module_group())
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

    def refresh_hydra_state(self):
        """Probe the AMR engine, tools and reference release once, when the plan opens.

        A run that cannot finish must not be offered as if it could, so an
        unavailable HYDRA is switched off and disabled with the reason beside it,
        rather than accepted here and failed per isolate after every other assay
        has already been done. The ready case still says which release will be
        searched and which isolates it holds no point-mutation catalogue for.
        """
        from wmlstudio import provisioning

        try:
            self.hydra_check = provisioning.hydra_prerequisites(
                selected=self.database_root or None)
        except (OSError, ValueError) as error:
            self.hydra_check = None
            self.hydra.setChecked(False)
            self.hydra.setEnabled(False)
            self.hydra_state.setText(f"The AMR reference store could not be read: {error}")
            return
        check = self.hydra_check
        self.hydra.setEnabled(check["ready"])
        if not check["ready"]:
            self.hydra.setChecked(False)
            self.hydra.setToolTip(check["message"])
            self.hydra_state.setText(check["message"])
            return
        notes = [check["database"]["label"]]
        covered = set(check["database"]["organisms"])
        uncovered = sorted({organism_line(sample).split(" (")[0] for sample in self.samples
                            if organism_line(sample) != "Unknown organism"
                            and not self.organism_covered(sample, covered)})
        if uncovered:
            notes.append("No point-mutation catalogue is installed for " + "; ".join(uncovered)
                         + ". Those isolates are screened for genes only, which is not evidence "
                         "that they carry no resistance mutation.")
        notes.extend(check["warnings"])
        self.hydra_state.setText(" ".join(notes))

    @staticmethod
    def organism_covered(sample, accepted):
        """True when the installed release has a catalogue this isolate can be run against."""
        from wmlstudio.hydra_runtime import match_organism

        genus, species, _ = organism_for(sample)
        name = f"{genus} {species}".strip()
        return bool(name) and match_organism(name, accepted=accepted) is not None

    def refresh_plan_table(self):
        """Fill the cohort rows, keeping whatever the user has already excluded."""
        excluded = {self.table.item(row, 0).data(Qt.ItemDataRole.UserRole)
                    for row in range(self.table.rowCount())
                    if self.table.item(row, 0) is not None
                    and self.table.item(row, 0).checkState() != Qt.CheckState.Checked}
        self.table.setRowCount(len(self.samples))
        for row, sample in enumerate(self.samples):
            item = cell("", sample["id"])
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked if sample["id"] in excluded else Qt.CheckState.Checked)
            self.table.setItem(row, 0, item)
            self.table.setItem(row, 1, cell(sample["name"]))
            self.table.setItem(row, 2, cell(organism_line(sample)))
            self.table.setItem(row, 3, cell(characterization_state(sample)[0]))

    def build_module_group(self):
        """One checkbox per registered module, with what it does and who it runs on."""
        group = QGroupBox("Organism-specific tools")
        box = QVBoxLayout(group)
        self.include_off_panel = QCheckBox("Also run these tools on isolates outside their reference taxa")
        self.module_warning = label("", "small", True)
        modules = registered_modules()
        if not modules:
            box.addWidget(label("No organism-specific tools are registered in this build.", "small", True))
            return group
        for key, module in modules.items():
            control = QCheckBox(module.title)
            control.setToolTip(module.purpose)
            # `clicked` fires only for the user, so a deliberate choice is never
            # overwritten by the next recount of the cohort.
            control.clicked.connect(lambda checked=False, name=key: self.module_toggled(name))
            self.module_boxes[key] = control
            note = label("", "small", True)
            self.module_notes[key] = note
            box.addWidget(control)
            if module.purpose:
                box.addWidget(label(module.purpose, "small", True))
            box.addWidget(note)
        self.include_off_panel.setToolTip("Off by default: a Staphylococcus panel run on a Klebsiella isolate "
                                          "produces a blank that is easily misread as a negative.")
        self.include_off_panel.stateChanged.connect(self.refresh_module_badges)
        box.addWidget(self.include_off_panel)
        box.addWidget(self.module_warning)
        box.addWidget(label(BOUNDARY, "small", True))
        self.refresh_module_badges()
        return group

    def module_toggled(self, key):
        self.module_touched.add(key)
        self.refresh_module_badges()

    def refresh_module_badges(self):
        """Recount applicability after any change to the cohort or the assignments."""
        included = self.included_samples()
        include_off_panel = self.include_off_panel.isChecked()
        warnings = []
        for key, module in registered_modules().items():
            counts = cohort_applicability(module, included)
            control, note = self.module_boxes[key], self.module_notes[key]
            if key not in self.module_touched:
                control.setChecked(counts["recommended"] > 0)
            note.setText(cohort_sentence(module, counts, off_panel_included=include_off_panel)
                         if included else "Include at least one isolate above.")
            if control.isChecked() and included and not (counts["recommended"] or counts["possible"]
                                                         or counts["unknown_organism"]):
                warnings.append(module.title)
        if warnings:
            self.module_warning.setText(
                "No selected isolate is within the reference taxa of: " + "; ".join(warnings) +
                ". Nothing will be recorded for them unless you also tick the box above, and an off-panel "
                "result is weak evidence in either direction.")
        else:
            self.module_warning.setText("")

    def included_samples(self):
        chosen = {self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) for row in range(self.table.rowCount())
                  if self.table.item(row, 0) is not None
                  and self.table.item(row, 0).checkState() == Qt.CheckState.Checked}
        return [sample for sample in self.samples if sample["id"] in chosen]

    def show_plan_menu(self, selection, position):
        """Right-click on the cohort: correct the organism, or change what is included."""
        identifiers = [value for value in selection.sample_ids if value]
        menu = QMenu(self)
        menu.setToolTipsVisible(True)
        assign = menu.addAction("Assign genus / species…" if len(identifiers) < 2
                                else f"Assign genus / species to {len(identifiers)} isolates…")
        assign.setEnabled(bool(identifiers) and self.project is not None)
        if self.project is None:
            assign.setToolTip("This plan was opened without the project, so assignments cannot be saved here.")
        menu.addSeparator()
        include = menu.addAction("Include in this run")
        exclude = menu.addAction("Exclude from this run")
        for action in (include, exclude):
            action.setEnabled(bool(identifiers))
        chosen = menu.exec(position)
        menu.deleteLater()
        if chosen is assign and identifiers:
            self.assign_organisms(identifiers)
        elif chosen in (include, exclude) and identifiers:
            self.set_included(identifiers, chosen is include)

    def set_included(self, sample_ids, included):
        wanted = set(sample_ids)
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item is not None and item.data(Qt.ItemDataRole.UserRole) in wanted:
                item.setCheckState(Qt.CheckState.Checked if included else Qt.CheckState.Unchecked)
        self.refresh_module_badges()

    def assign_organisms(self, sample_ids):
        """Manual genus/species override, applied through the one storage entry point."""
        from wmlstudio.storage import assign_organism
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        wanted = set(sample_ids)
        chosen = [sample for sample in self.samples if sample["id"] in wanted]
        if not chosen or self.project is None:
            return
        dialog = BatchAssignmentDialog(chosen, self.scheme_entries, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            for assignment in dialog.assignments:
                assign_organism(self.project, [assignment["sample_id"]], assignment["genus"],
                                assignment["species"], scheme_path=assignment.get("scheme_path"),
                                typing_mode=assignment["typing_mode"])
        except (ValueError, KeyError) as error:
            self.feedback.setText(str(error))
            return
        records = {sample["id"]: sample for sample in self.project.samples()}
        self.samples = [records.get(sample["id"], sample) for sample in self.samples]
        self.refresh_plan_table()
        self.refresh_module_badges()
        # Which isolates have a point-mutation catalogue follows the assignment too.
        self.refresh_hydra_state()
        self.feedback.setText("Assignment saved. The organism-specific tools offered now follow it.")

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
        # Refuse here rather than halfway through the cohort: without the engine,
        # the tools or a populated reference store every isolate would fail after
        # its species and virulence work had already been done.
        if self.hydra.isChecked() and self.hydra_check is not None and not self.hydra_check["ready"]:
            self.feedback.setText(self.hydra_check["message"])
            return
        try:
            allocation = plan_resources(threads_per_sample=4, memory_gb=3, policy=self.policy.currentData())
        except ValueError as error:
            self.feedback.setText(str(error))
            return
        modules = {key: control.isChecked() for key, control in self.module_boxes.items()}
        if any(modules.values()) and not (Path(reference) / "manifest.json").is_file():
            self.feedback.setText("Organism-specific tools need a characterization snapshot. Install or select one, "
                                  "or clear those tools.")
            return
        self.plan = {"sample_ids": identifiers, "reference_path": reference or None,
                     "species": self.species.isChecked(), "virulence": self.virulence.isChecked(),
                     "hydra": self.hydra.isChecked(), "modules": modules,
                     "modules_off_panel": self.include_off_panel.isChecked(),
                     "resources": allocation.to_dict()}
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
        # One extra column per registered organism-specific tool, so a tool that is
        # added later needs no change here.
        self.characterization_modules = list(registered_modules().values())
        self.characterization_table = make_table(
            ["Sample", "Evidence state", "Species candidate", "ANI %", "Classical ST", "Virulence loci",
             "Drug classes", "Replicons"] + [module.column_title for module in self.characterization_modules])
        for offset, module in enumerate(self.characterization_modules):
            header = self.characterization_table.horizontalHeaderItem(8 + offset)
            if header is not None:
                header.setToolTip(module.purpose or module.title)
        self.install_view_menu("characterization", self.characterization_table)
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
                      (species.get("nearest") or {}).get("ani"), classical_st(sample),
                      "; ".join(group["locus"] for group in groups if group.get("status") == "detected") or virulence.get("status", "not_run"),
                      "; ".join(sorted({entry.get("class") or entry.get("subclass") or "" for entry in usable.get("drug_associations", {}).get("associations", [])})) or usable.get("drug_associations", {}).get("status", "not_run"),
                      "; ".join(entry["gene"] for entry in usable.get("plasmid_hypotheses", {}).get("replicons", [])) or usable.get("plasmid_hypotheses", {}).get("status", "not_run")]
            # The stale/unverified gate above applies here too: `usable` is empty and
            # every module reads not_run rather than showing an earlier assembly's call.
            values += [module.summary(usable.get(module.key) or {})
                       for module in getattr(self, "characterization_modules", [])]
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
        except KeyError:
            self.characterization_detail.clear()
            return
        self.characterization_detail.setHtml(characterization_html(sample))
        if self.characterization_table.selectedItems():
            self.focus.set_focus([sample["id"]], "Evidence table selection")

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

    def context_run_characterization(self, selection):
        """Right-click route into the same reviewed plan, for the clicked isolates only."""
        if self.busy():
            return
        wanted = {value for value in selection.sample_ids if value}
        records = [sample for sample in isolate_records(self.project.samples()) if sample["id"] in wanted]
        if not records:
            self.notify("Select one or more isolates first.")
            return
        self.characterize_records(records)

    def run_characterization_selected(self):
        if self.busy():
            return
        records = isolate_records(self.project.samples())
        from wmlstudio.cohort_picker import CohortPickerDialog
        picker = CohortPickerDialog(records, self.project, "Which assemblies should be characterized?", parent=self, include_reads=False)
        if picker.exec() != QDialog.DialogCode.Accepted:
            return
        self.characterize_records([sample for sample in records if sample["id"] in picker.selected_ids])

    def characterize_records(self, records):
        """Review, then run: the one plan every entry point goes through."""
        assemblies = [sample for sample in records if sample.get("input_path") and not sample.get("missing_input")
                      and ((sample.get("result") or {}).get("kind") == "fasta" or
                           Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fa", ".fna", ".fasta", ".fas")))]
        if not assemblies:
            self.notify("Select/import FASTA assemblies for characterization. Review paired FASTQ assembly first if you only have reads.")
            self.navigate(1)
            return
        # The assignment this plan can make sets the scheme a sample is TYPED
        # against, so it offers classical schemes only: a 2,000-target cgMLST
        # scheme chosen where a seven-locus one is expected is a different
        # quantity under the same label.
        dialog = CharacterizationPlanDialog(
            assemblies, self.active_characterization_reference(), self, project=self.project,
            scheme_entries=self.scheme_entries(CLASSICAL_KINDS),
            database_root=self.active_amr_database())
        dialog.installRequested.connect(self.install_characterization_references)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        plan = dialog.plan
        self.feature_ids = set(plan["sample_ids"])
        # Re-read: the plan dialog may have corrected an organism assignment, and
        # which modules apply follows that assignment.
        current = {sample["id"]: sample for sample in self.project.samples()}
        samples = [current.get(sample["id"], sample) for sample in assemblies if sample["id"] in plan["sample_ids"]]
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
                genus, species, _ = organism_for(sample)
                assigned = (sample.get("metadata") or {}).get("organism") or {}
                # A provisional detection is not a verified mutation-catalog assignment,
                # so only an assigned genus AND species chooses a point-mutation catalog.
                # The engine's taxgroups are underscore-joined and genus-level for some
                # organisms, so run_assemblies resolves this name against the installed
                # catalogue: an unknown one stops the whole run upstream, and a name it
                # cannot resolve is reported rather than replaced by a near neighbour.
                organism = f"{genus} {species}" if assigned.get("genus") and assigned.get("species") else None
                if plan["hydra"] and hydra is None:
                    generated = run_assemblies([sample["input_path"]], database_root,
                                              sample_names=[sample["id"]], threads=resources.threads_per_sample,
                                              organism=organism, point_mutations=True, cancelled=stopped, progress=report)
                    hydra = generated
                modules, skipped = selection_for_record(sample, plan.get("modules"),
                                                        include_off_panel=plan.get("modules_off_panel", False))
                result = characterize_assembly(sample["input_path"], plan["reference_path"], hydra_report=hydra,
                                               cancelled=stopped, progress=report, threads=resources.threads_per_sample,
                                               species=plan["species"], virulence=plan["virulence"],
                                               modules=modules, organism=(genus, species))
                record_skipped(result, skipped)
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
            chosen = [key for key, value in (plan.get("modules") or {}).items() if value]
            tools = f" Organism-specific tools requested: {', '.join(chosen)}." if chosen else ""
            self.notify(f"Characterization saved for {len(result['completed'])} isolates; {len(result['failures'])} failed."
                        f"{tools} Inspect individual assay states and provenance.")
            self.navigate(4)
            if hasattr(self, "evidence_tabs"):
                self.evidence_tabs.setCurrentIndex(0)
        self.launch_task(operation, "characterization", finished)
