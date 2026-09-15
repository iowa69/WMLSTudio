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


def assigned_organism_name(sample):
    """The isolate's assigned 'Genus species', or '' when one was never assigned.

    A provisional detection is not a verified assignment, so only an assigned
    genus AND species chooses a point-mutation catalogue or an organism-curated
    virulence search. The plan dialog and the run itself read this one rule, so
    what the dialog promises is what the run does.
    """
    genus, species, source = organism_for(sample)
    return f"{genus} {species}".strip() if source == "Assigned" and genus and species else ""


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


def plasmid_html(plasmids):
    """Contig-level plasmid evidence, and the sentences saying what it is not.

    Every row is about one assembled contig: which replicon markers sit on it,
    which determinants share it, what the assembler itself claimed about closure
    and coverage, and which of those signals are present. Nothing here is a
    reconstructed plasmid, a plasmid count or a mobility prediction, so the
    engine's own gap notes are printed word for word rather than summarised into
    something that would read as a smaller caveat than it is.
    """
    def escape(value):
        return html.escape(str(value if value is not None else "—"))

    def genes(row, key):
        return "; ".join(entry["gene"] for entry in row.get(key) or [] if entry.get("gene")) or "—"

    body = (f"<h3>Plasmid evidence · {escape(plasmids.get('status', 'not_run'))}</h3>"
            f"<p>{escape(plasmids.get('reason', ''))}</p>")
    rows = plasmids.get("contig_evidence") or []
    if rows:
        body += ("<table cellpadding='5'><tr><th>Contig</th><th>Length (bp)</th>"
                 "<th>Closure (assembler's claim)</th><th>Coverage vs backbone</th>"
                 "<th>Replicon markers</th><th>Determinants on the same contig</th>"
                 "<th>Evidence present</th></tr>")
        for row in rows:
            ratio = row.get("depth_ratio")
            body += ("<tr>" + f"<td>{escape(row.get('contig'))}</td>"
                     f"<td>{escape(row.get('length_bp'))}</td>"
                     f"<td>{escape(row.get('closure'))}</td>"
                     f"<td>{escape('%.1f×' % ratio if ratio else 'unknown')}</td>"
                     f"<td>{escape(genes(row, 'replicons'))}</td>"
                     f"<td>{escape(genes(row, 'determinants'))}</td>"
                     f"<td>{escape(row.get('support_state'))}</td></tr>")
        basis = plasmids.get("depth_basis") or {}
        body += f"</table><p>{escape(basis.get('reason', ''))} {escape(basis.get('basis', ''))}</p>"
    for association in plasmids.get("contig_associations", []):
        body += (f"<p><b>{escape(association.get('replicon'))}</b> + {escape(association.get('marker'))}"
                 f" on {escape(association.get('contig'))}; {escape(association.get('gap_bp'))} bp apart. "
                 "Same-contig evidence only; complete plasmid identity and transmission are unproven.</p>")
    # Every determinant is listed, the unplaced ones included: a determinant that
    # is absent from this list would read as a determinant that was not found.
    for row in plasmids.get("determinant_placement", []):
        where = f" on {escape(row.get('contig'))}" if row.get("contig") else ""
        body += (f"<p><b>{escape(row.get('gene'))}</b> — {escape(row.get('placement'))}{where}. "
                 f"{escape(row.get('interpretation'))}</p>")
    if plasmids.get("mobility") or plasmids.get("plasmid_count"):
        body += (f"<p>Mobility: {escape(plasmids.get('mobility'))} · plasmid count: "
                 f"{escape(plasmids.get('plasmid_count'))} · reconstruction: "
                 f"{escape(plasmids.get('reconstruction'))}.</p>")
    for note in plasmids.get("mob_suite_gap", []):
        body += f"<p>{escape(note)}</p>"
    return body


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
    body += plasmid_html(evidence.get("plasmid_hypotheses") or {})
    body += f"<h3>Provenance</h3><p>Assembly SHA-256: {escape(evidence.get('input_sha256'))}</p>"
    for section in (species, virulence):
        digest = section.get("reference_digest") or section.get("provenance", {}).get("reference_digest")
        if digest:
            body += f"<p>Reference fingerprint: {escape(digest)}</p>"
    for note in evidence.get("limitations", []):
        body += f"<p>{escape(note)}</p>"
    return body


def plasmid_cohort(records, *, cancelled=None):
    """Replicon and replicon/determinant co-occurrence across the chosen isolates.

    An isolate whose characterization is missing or stale is named with that
    reason rather than counted as carrying nothing: an empty row and an unasked
    question are the same picture otherwise. Recurrence across isolates is not
    transmission, not one plasmid and not a relatedness measure, and this
    produces no distance and no threshold anybody could cluster on.
    """
    from wmlstudio.characterization import current_characterization
    from wmlstudio.plasmid_evidence import cohort_replicon_cooccurrence

    entries = []
    for record in records:
        state = current_characterization(record)
        evidence = state.get("evidence") or {}
        entries.append({"sample_id": record["id"], "sample_name": record.get("name") or record["id"],
                        "plasmid_hypotheses": evidence.get("plasmid_hypotheses")
                        or {"status": "not_run", "reason": state["reason"]}})
    return cohort_replicon_cooccurrence(entries, cancelled=cancelled)


class CharacterizationPlanDialog(QDialog):
    installRequested = Signal()
    databasesRequested = Signal()

    def __init__(self, samples, reference_path=None, parent=None, *, project=None, scheme_entries=(),
                 database_root=None):
        super().__init__(parent)
        self.plan = None
        self.samples = list(samples)
        self.project = project
        self.database_root = str(database_root or "")
        self.hydra_check = None
        self.catalogue = {"entries": [], "summary": ""}
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
        # HYDRA is powerful and silent about its inputs. What it is about to
        # search, which reference sets this computer holds, which isolates have a
        # mutation catalogue and whether virulence elements are read at all are
        # all answered here, before the run, rather than discovered afterwards in
        # a result where "not searched for" and "not found" look the same.
        hydra_options = QFormLayout()
        self.point_mutations = QCheckBox(
            "Search the isolate's own point-mutation catalogue where the release has one")
        self.point_mutations.setChecked(True)
        self.point_mutations.stateChanged.connect(self.render_hydra_state)
        hydra_options.addRow("Mutation evidence", self.point_mutations)
        self.hydra_virulence = QComboBox()
        self.hydra_virulence.addItem("Where the isolate's organism is established", None)
        self.hydra_virulence.addItem("Always, even with no organism (uncurated)", True)
        self.hydra_virulence.addItem("Never; acquired resistance only", False)
        self.hydra_virulence.setToolTip(
            "The engine's virulence and stress curation is keyed to the organism, so a run with no "
            "assigned organism is a different, uncurated search and is reported as one.")
        self.hydra_virulence.currentIndexChanged.connect(self.refresh_hydra_state)
        hydra_options.addRow("Virulence and stress elements", self.hydra_virulence)
        layout.addLayout(hydra_options)
        self.hydra_state = label("", "small", True)
        layout.addWidget(self.hydra_state)
        self.hydra_databases = label("", "small", True)
        layout.addWidget(self.hydra_databases)
        databases_row = QHBoxLayout()
        databases_row.addWidget(button("Choose / download reference databases…", self.request_databases))
        databases_row.addStretch()
        layout.addLayout(databases_row)
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
                selected=self.database_root or None,
                virulence=self.hydra_virulence.currentData())
            self.catalogue = provisioning.hydra_database_catalogue(
                selected=self.database_root or None, measure=False)
        except (OSError, ValueError) as error:
            self.hydra_check = None
            self.hydra.setChecked(False)
            self.hydra.setEnabled(False)
            self.hydra_state.setText(f"The AMR reference store could not be read: {error}")
            self.hydra_databases.setText("")
            return
        self.render_hydra_state()

    def render_hydra_state(self):
        """Re-word what the last probe found; ticking a box re-reads no reference table."""
        check = self.hydra_check
        if check is None:
            return
        self.hydra_databases.setText(self.database_sentence(check, self.catalogue))
        self.hydra.setEnabled(check["ready"])
        if not check["ready"]:
            self.hydra.setChecked(False)
            self.hydra.setToolTip(check["message"])
            self.hydra_state.setText(check["message"])
            return
        notes = [check["database"]["label"], self.mutation_sentence(check),
                 self.virulence_sentence(check)]
        # The engine's own organism and virulence sentences are answered per
        # isolate above, so the cohort-level copies of them are not repeated:
        # this check was made with no organism, and saying "no organism was
        # given" beside a cohort that has them would be false.
        answered = {check["organism"]["reason"], check["virulence"]["reason"]}
        notes.extend(note for note in check["warnings"] if note not in answered)
        self.hydra_state.setText(" ".join(filter(None, notes)))

    def mutation_sentence(self, check):
        """Which of these isolates the release actually holds a mutation catalogue for."""
        if not self.point_mutations.isChecked():
            return ("Point mutations are switched off for this run, so none is reported for any "
                    "isolate here. That is not evidence that none is present.")
        accepted, catalogues = (check["database"]["organisms"],
                                check["database"]["point_mutation_organisms"])
        uncovered = sorted({organism_line(sample).split(" (")[0] for sample in self.samples
                            if organism_line(sample) != "Unknown organism"
                            and not self.organism_covered(sample, accepted, catalogues)})
        unassigned = sum(1 for sample in self.samples if not assigned_organism_name(sample))
        parts = []
        if uncovered:
            parts.append("No point-mutation catalogue is installed for " + "; ".join(uncovered)
                         + ". Those isolates are screened for genes only, which is not evidence "
                         "that they carry no resistance mutation.")
        if unassigned:
            parts.append(f"{unassigned} of {len(self.samples)} isolates have no assigned genus and "
                         "species, so no catalogue is chosen for them; assign one to have their "
                         "mutations assessed.")
        return " ".join(parts)

    def virulence_sentence(self, check):
        """What the translated search will do about virulence and stress elements.

        The engine's curation is keyed to the organism, so an uncurated screen and
        an organism's curated one are two different results and are never reported
        as the same one. Whatever this says, the nucleotide catalogues report their
        own virulence-typed genes regardless: "off" never means none was sought.
        """
        block = check["virulence"]
        if block["requested"] is False or not block["available"]:
            return block["reason"]
        counts = (f"{block['virulence_elements']} virulence and {block['stress_elements']} stress "
                  "elements in the protein reference")
        if block["requested"] is True:
            return (f"{counts} will be searched for every isolate here, including those with no "
                    "established organism: that search is uncurated, so nothing an organism's "
                    "curation would have suppressed is suppressed. A gene is not a demonstrated "
                    "virulence phenotype.")
        established = sum(1 for sample in self.samples if assigned_organism_name(sample))
        return (f"{counts} will be searched for the {established} of {len(self.samples)} isolates "
                "whose genus and species are assigned, curated for that organism. The rest are "
                "searched for acquired resistance only, which is not evidence that they carry no "
                "virulence gene. A gene is not a demonstrated virulence phenotype.")

    @staticmethod
    def database_sentence(check, catalogue):
        """Which reference sets this run reads, and which are simply not here."""
        searched = ", ".join(check["databases"]) or "none"
        return f"This run will search: {searched}. " + (catalogue.get("summary") or "")

    @staticmethod
    def organism_covered(sample, accepted, catalogues=None):
        """True when the installed release has a catalogue this isolate can be run against.

        ``accepted`` is every taxgroup the release will take at all; ``catalogues``,
        when given, is the smaller set that actually holds point mutations. Being
        accepted is not being covered, and both questions are decided by the
        engine's own parent/child taxgroup rule rather than by exact membership.
        """
        from wmlstudio.hydra_runtime import catalogue_covers, match_organism

        genus, species, _ = organism_for(sample)
        name = f"{genus} {species}".strip()
        resolved = match_organism(name, accepted=accepted) if name else None
        if resolved is None:
            return False
        return True if catalogues is None else catalogue_covers(resolved, catalogues)

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

    def request_databases(self):
        """Ask for the reference-database list. Nothing is downloaded from this dialog."""
        self.databasesRequested.emit()

    def set_database_root(self, path):
        """Follow a snapshot chosen or installed from that list, and re-read it."""
        self.database_root = str(path or "")
        self.refresh_hydra_state()

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
                     # HYDRA's own options, kept apart from the defined virulence
                     # panel above: these two search different references and
                     # answer different questions.
                     "point_mutations": self.point_mutations.isChecked(),
                     "hydra_virulence": self.hydra_virulence.currentData(),
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
        tabs.insertTab(0, panel, "Identity / virulence / plasmid evidence")
        tabs.insertTab(1, self.build_plasmid_cohort(), "Plasmid evidence across the cohort")

    def build_plasmid_cohort(self):
        """Which replicon markers, and which replicon/determinant pairs, recur here.

        Two tables, never one: a replicon and a determinant that share an
        assembled contig and a pair that is merely present in the same isolate
        are different observations, and collapsing them would turn an assembly
        artefact into a shared finding. Neither table is sorted or coloured by
        count, because a ranked list of co-occurrence counts reads as clustering
        and nothing here is a distance.
        """
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.plasmid_scope = label("", "small", True)
        layout.addWidget(self.plasmid_scope)
        layout.addWidget(label("Replicon markers detected, by isolate count", "cardTitle"))
        self.plasmid_replicon_table = make_table(
            ["Replicon marker", "Isolates carrying it", "Of isolates assayed", "Isolates"])
        self.plasmid_replicon_table.setSortingEnabled(False)
        layout.addWidget(self.plasmid_replicon_table, 1)
        layout.addWidget(label("Replicon and determinant on one assembled contig", "cardTitle"))
        self.plasmid_pair_table = make_table(
            ["Replicon marker", "Determinant", "Marker type", "Co-located on one contig",
             "Both present, not co-located", "Of isolates assayed"])
        self.plasmid_pair_table.setSortingEnabled(False)
        layout.addWidget(self.plasmid_pair_table, 1)
        self.plasmid_notes = QTextBrowser()
        self.plasmid_notes.setOpenExternalLinks(False)
        layout.addWidget(self.plasmid_notes, 1)
        return panel

    def refresh_plasmid_cohort(self, samples):
        """Recount co-occurrence over the evidence cohort, denominators and all."""
        if not hasattr(self, "plasmid_replicon_table"):
            return
        result = plasmid_cohort(samples)
        names = {sample["id"]: sample["name"] for sample in samples}
        denominator = result["denominator"]
        self.plasmid_replicon_table.setRowCount(len(result["replicons"]))
        for row, entry in enumerate(result["replicons"]):
            members = "; ".join(sorted(names.get(value, value) for value in entry["isolates"]))
            for column, value in enumerate([entry["replicon"], entry["isolate_count"],
                                            denominator, members]):
                self.plasmid_replicon_table.setItem(row, column, cell(value))
        self.plasmid_pair_table.setRowCount(len(result["co_occurrence"]))
        for row, entry in enumerate(result["co_occurrence"]):
            values = [entry["replicon"], entry["marker"], entry["marker_type"] or "—",
                      entry["co_located_count"], entry["not_co_located_count"], denominator]
            for column, value in enumerate(values):
                item = cell(value)
                item.setToolTip(entry["interpretation"])
                self.plasmid_pair_table.setItem(row, column, item)
        excluded = result["isolates_excluded"]
        self.plasmid_scope.setText(
            f"{denominator} of {len(samples)} isolates in this evidence cohort had a plasmid-marker "
            f"assay to read; {len(excluded)} did not and are listed below rather than counted as "
            "carrying nothing. Rows are in alphabetical order, not ranked: these are counts of what "
            "was seen, not a distance and not a cluster.")
        self.plasmid_notes.setHtml(self.plasmid_notes_html(result, names))

    @staticmethod
    def plasmid_notes_html(result, names):
        """The exclusions, the unplaced determinants and the claim boundary, verbatim."""
        def escape(value):
            return html.escape(str(value if value is not None else "—"))

        body = ""
        if result["isolates_excluded"]:
            body += "<h3>Not assayed, so not counted either way</h3>"
            for entry in result["isolates_excluded"]:
                body += (f"<p><b>{escape(names.get(entry['sample_id'], entry['sample_name']))}</b> — "
                         f"{escape(entry['reason'])}</p>")
        if result["unplaced_determinants"]:
            body += ("<h3>Determinants with no usable contig coordinates</h3><p>These are neither "
                     "co-located nor apart: counting them as apart would read as evidence of "
                     "separation.</p><p>"
                     + escape("; ".join(sorted({f"{names.get(row['sample_id'], row['sample_id'])}: "
                                                f"{row['gene']}"
                                                for row in result["unplaced_determinants"]})))
                     + "</p>")
        body += "<h3>What this table is not</h3>"
        for note in [*result["limitations"], *result["mob_suite_gap"]]:
            body += f"<p>{escape(note)}</p>"
        return body

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
        self.refresh_plasmid_cohort(samples)

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

    def manage_characterization_databases(self, plan_dialog):
        """Show the whole reference-database list over the plan, then re-read the store.

        The plan stays open: a person who discovers here that a reference set is
        missing should be able to install it and carry on, not lose the cohort
        they had just reviewed. Nothing is downloaded without their confirmation.
        """
        from wmlstudio.amr_databases import AMRDatabaseDialog
        dialog = AMRDatabaseDialog(self.active_amr_database(), plan_dialog,
                                   update_root=self.root / "references" / "hydra",
                                   project=self.project)
        self.amr_database_dialog = dialog
        dialog.snapshotInstalled.connect(
            lambda path: self.project.set_setting("hydra_database_root", path))
        dialog.snapshotInstalled.connect(plan_dialog.set_database_root)
        dialog.exec()
        plan_dialog.refresh_hydra_state()

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
        dialog.databasesRequested.connect(lambda: self.manage_characterization_databases(dialog))
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
                # A provisional detection is not a verified mutation-catalog assignment,
                # so only an assigned genus AND species chooses a point-mutation catalog.
                # The engine's taxgroups are underscore-joined and genus-level for some
                # organisms, so run_assemblies resolves this name against the installed
                # catalogue: an unknown one stops the whole run upstream, and a name it
                # cannot resolve is reported rather than replaced by a near neighbour.
                organism = assigned_organism_name(sample) or None
                if plan["hydra"] and hydra is None:
                    generated = run_assemblies([sample["input_path"]], database_root,
                                              sample_names=[sample["id"]], threads=resources.threads_per_sample,
                                              organism=organism,
                                              point_mutations=plan.get("point_mutations", True),
                                              virulence=plan.get("hydra_virulence"),
                                              cancelled=stopped, progress=report)
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
