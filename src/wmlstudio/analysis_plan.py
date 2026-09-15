"""Human-readable launch review: explicit analysis scope, resources and AMR options."""

from pathlib import Path

from PySide6.QtCore import Qt, Signal
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

from wmlstudio.sample_workflow import cgmlst_scheme_applies, sample_configuration
from wmlstudio.scheduler import GIB, detect_hardware
from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import label

# A cgMLST run and a classical MLST run are two different quantities produced in
# one launch, never one analysis with two labels. The chosen cgMLST scheme travels
# with the flag so the plan says which target set it would call against.
CGMLST_SCOPE = (
    "cgMLST is run against the installed target set you choose here, alongside the MLST scheme "
    "above. The two distances are different quantities: they never share a scale, an axis or a "
    "threshold, and neither is a phylogeny.")


def _cgmlst_choice(row) -> dict:
    """One offerable cgMLST scheme, flattened from a cgMLST library or catalogue row.

    Both spellings the library uses are accepted -- a catalogue row names itself
    'key' and 'scheme_name', an installed row 'catalog_key' and 'name' -- so this
    dialog does not have to care which listing it was handed.
    """
    row = row if isinstance(row, dict) else {}
    installed = row.get("installed") if isinstance(row.get("installed"), dict) else {}
    return {"key": row.get("key") or row.get("catalog_key"),
            "organism": str(row.get("organism") or ""),
            "genus": str(row.get("genus") or ""), "species": str(row.get("species") or ""),
            "scheme_name": str(row.get("scheme_name") or row.get("name") or ""),
            "locus_count": installed.get("locus_count") or row.get("locus_count"),
            "path": str(installed.get("path") or row.get("path") or ""),
            "ready": bool(row.get("ready", bool(installed)))}


def installed_cgmlst_choices(rows=None, data_root=None) -> list[dict]:
    """The cgMLST schemes this computer can actually run, in the order they were listed.

    A scheme that is catalogued but not installed is not offered: choosing it would
    promise a run that has no target set behind it. Reading the library is optional,
    so a caller that already holds the rows -- the cgMLST schemes tab does -- never
    pays for a second folder scan.
    """
    if rows is None:
        if not data_root:
            return []
        from wmlstudio.cgmlst_schemes import library_status
        try:
            rows = library_status(data_root)
        except (OSError, ValueError):
            return []
    choices = [_cgmlst_choice(row) for row in rows]
    return [choice for choice in choices if choice["path"] and choice["key"]]


class RunPlanDialog(QDialog):
    manageDatabases = Signal()

    def __init__(self, samples, scheme=None, db_root=None, parent=None, *,
                 cgmlst_schemes=None, data_root=None):
        super().__init__(parent)
        self.samples = list(samples)
        self.cgmlst_choices = installed_cgmlst_choices(cgmlst_schemes, data_root)
        self.setWindowTitle("Review analysis plan")
        self.resize(850, 670)
        self.plan = {}
        self.catalogue = {"entries": [], "summary": "", "error": ""}
        # The release and its organism tables are read from disk, so they are read
        # once per store rather than again on every tick of a checkbox.
        self._status_cache = (None, None)
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
            # Read from the one sample configuration, so this table shows what the
            # other menus show rather than a second copy of the same choices.
            configuration = sample_configuration(sample)
            values = [
                sample["name"],
                configuration["organism"] or "Unknown",
                "Scheme override" if scheme else configuration["typing_mode"] or "QC only",
                Path(scheme or configuration["scheme_path"] or "Automatic / not assigned").name,
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
                "Assemblies follow the assigned MLST workflow, and a cgMLST run can be added to the same launch below. FASTQ inputs receive quality checks unless assembled first.",
                "small",
                True,
            )
        )
        self.hydra = QCheckBox("Also run HYDRA AMR analysis on the selected assemblies")
        layout.addWidget(self.hydra)
        # HYDRA runs against whatever reference data happens to be installed and
        # reports nothing for everything else; in a result those two are
        # indistinguishable. So which sets are here, which the run will read, and
        # which of these isolates have a mutation catalogue at all are stated
        # before the run rather than discovered in an empty report afterwards.
        self.hydra_state = label("", "small", True)
        layout.addWidget(self.hydra_state)
        self.hydra_gaps = label("", "small", True)
        layout.addWidget(self.hydra_gaps)
        self.cgmlst = QCheckBox("Also run cgMLST on the selected assemblies")
        layout.addWidget(self.cgmlst)
        self.cgmlst_choice = QComboBox()
        for choice in self.cgmlst_choices:
            targets = f"{choice['locus_count']} targets" if choice["locus_count"] else "installed"
            self.cgmlst_choice.addItem(
                f"{choice['organism'] or 'Unknown organism'} · {choice['scheme_name']} ({targets})",
                choice)
        self.cgmlst.setEnabled(bool(self.cgmlst_choices))
        if not self.cgmlst_choices:
            self.cgmlst.setToolTip("No cgMLST scheme is installed yet. Install or download one in "
                                   "the cgMLST schemes tab, then run it from here.")
        cgmlst_row = QHBoxLayout()
        cgmlst_row.addWidget(label("cgMLST scheme", "small", True))
        cgmlst_row.addWidget(self.cgmlst_choice, 1)
        layout.addLayout(cgmlst_row)
        self.cgmlst_choice.setEnabled(False)
        self.cgmlst.toggled.connect(self.cgmlst_choice.setEnabled)
        layout.addWidget(label(CGMLST_SCOPE, "small", True))
        self.restore_sample_flags()
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
        self.point_mutations.toggled.connect(self.refresh_databases)
        form.addRow("Mutation evidence", self.point_mutations)
        # The engine's virulence and stress curation is keyed to the organism, so
        # a run with no established organism is a different, uncurated search. It
        # is offered as three explicit choices and reported for what each one is.
        self.virulence = QComboBox()
        self.virulence.addItem("Where the isolate's organism is established", None)
        self.virulence.addItem("Always, even with no organism (uncurated)", True)
        self.virulence.addItem("Never; acquired resistance only", False)
        self.virulence.currentIndexChanged.connect(self.refresh_databases)
        form.addRow("Virulence and stress elements", self.virulence)
        form.addRow(label("Virulence and stress elements come from the protein reference. The "
                          "nucleotide catalogues report their own virulence-typed genes whatever "
                          "this is set to, and a gene found is never a demonstrated phenotype.",
                          "small", True))
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

    def restore_sample_flags(self):
        """Open with what the samples themselves already say they want run.

        The flags are part of the one sample configuration, so a cohort that was
        configured to run cgMLST arrives here already ticked instead of asking again.
        A flag is honoured only when every selected sample carries it; a mixed
        selection is left for the person to decide.
        """
        configurations = [sample_configuration(sample) for sample in self.samples]
        if configurations and all(config["run_hydra"] for config in configurations):
            self.hydra.setChecked(True)
        keys = {config["cgmlst_scheme_key"] for config in configurations}
        if not configurations or len(keys) != 1 or not all(config["run_cgmlst"] for config in configurations):
            return
        wanted = next(iter(keys))
        for row in range(self.cgmlst_choice.count()):
            if (self.cgmlst_choice.itemData(row) or {}).get("key") == wanted:
                self.cgmlst_choice.setCurrentIndex(row)
                self.cgmlst.setChecked(self.cgmlst.isEnabled())
                return

    def selected_cgmlst_scheme(self):
        """The cgMLST scheme this plan would call against, or None when none is chosen."""
        return self.cgmlst_choice.currentData() if self.cgmlst.isChecked() else None

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
        """Read the chosen store and say what it holds, what it does not, and what runs.

        The whole catalogue, not only the installed part: a reference set that is
        absent reports nothing, and nothing reported is indistinguishable from a
        clean isolate unless the absence is named before the run. Each installed
        set is offered by name with what it is searched for, so choosing one is a
        decision about a question rather than about a word.
        """
        from wmlstudio.hydra_runtime import database_catalogue

        selected = self.database_choice.currentData()
        self.database_choice.blockSignals(True)
        self.database_choice.clear()
        self.database_choice.addItem("All installed nucleotide / protein databases", None)
        try:
            self.catalogue = database_catalogue(self.database.text() or None, measure=False)
        except (OSError, ValueError) as exc:
            self.database_choice.blockSignals(False)
            self.feedback.setText(str(exc))
            self.hydra_state.setText(str(exc))
            return
        for entry in self.catalogue["entries"]:
            if not entry["installed"]:
                continue
            self.database_choice.addItem(f"{entry['name']} — {entry['purpose']}", entry["name"])
            self.database_choice.setItemData(
                self.database_choice.count() - 1,
                "\n\n".join(filter(None, [entry["title"], entry["purpose"], entry["provider"],
                                          f"Licence: {entry['licence']}"])),
                Qt.ItemDataRole.ToolTipRole)
        self.database_choice.setCurrentIndex(max(0, self.database_choice.findData(selected)))
        self.database_choice.blockSignals(False)
        self.hydra_state.setText(self.catalogue["error"] or self.catalogue["summary"])
        self.hydra_gaps.setText(self.coverage_sentence())

    def database_state(self):
        """The chosen store's release and organism tables, read once per store."""
        from wmlstudio.hydra_runtime import database_status

        path = self.database.text() or ""
        if self._status_cache[0] != path:
            self._status_cache = (path, database_status(path or None))
        return self._status_cache[1]

    def coverage_sentence(self):
        """What this store could not report for these isolates, before anyone starts.

        Being accepted by a release and having a point-mutation catalogue in it
        are two different things, and an isolate outside the second is screened
        for genes only — which is not evidence that it carries no mutation.
        """
        from wmlstudio.hydra_runtime import catalogue_covers, match_organism

        status = self.database_state()
        if status["error"] or not status["installed"]:
            return (status["error"] or "No reference set is installed in this store, so a HYDRA "
                    "run would have nothing to search and would report no determinant for any "
                    "isolate. That is not a negative result.")
        parts = []
        if self.point_mutations.isChecked():
            configurations = [sample_configuration(sample) for sample in self.samples]
            named = [config for config in configurations
                     if config["organism_source"] == "assigned" and config["genus"]
                     and config["species"]]
            uncovered = sorted({config["organism"] for config in named
                                if not catalogue_covers(
                                    match_organism(config["organism"],
                                                   accepted=status["organisms"]) or "",
                                    status["point_mutation_organisms"])})
            if len(named) < len(configurations):
                parts.append(f"{len(configurations) - len(named)} of {len(configurations)} inputs "
                             "have no assigned genus and species, so no point-mutation catalogue "
                             "is chosen for them and none is reported; that is not evidence that "
                             "they carry no resistance mutation.")
            if uncovered:
                parts.append("This release holds no point-mutation catalogue for "
                             + "; ".join(uncovered) + ". Those isolates are screened for genes "
                             "only, and no other organism's catalogue is substituted for them.")
        else:
            parts.append("Point mutations are switched off, so none is reported for any isolate in "
                         "this run. That is not evidence that none is present.")
        if self.virulence.currentData() is False:
            parts.append("The translated search is limited to acquired resistance, so no virulence "
                         "or stress element is read from the protein reference.")
        elif self.virulence.currentData() is True:
            parts.append("Virulence and stress elements are searched for every isolate here, "
                         "including those with no established organism: that search is uncurated.")
        if status["stale"]:
            parts.append(f"This release is {status['age_days']} days old ({status['release']}); "
                         "determinants named after it are not in it.")
        return " ".join(parts)

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
                check = provisioning.hydra_prerequisites(
                    selected=self.database.text() or None,
                    virulence=self.virulence.currentData())
                if not check["ready"]:
                    self.feedback.setText(check["message"])
                    return
                if check.get("warnings"):
                    self.feedback.setText(check["database"]["label"] + " " + " ".join(check["warnings"]))
            except (OSError, ValueError) as exc:
                self.feedback.setText(str(exc))
                return
        cgmlst_scheme = self.selected_cgmlst_scheme()
        if self.cgmlst.isChecked():
            # Refuse here, naming the organism the scheme is for, rather than calling
            # a genome against a target set that was never defined for it.
            if cgmlst_scheme is None:
                self.feedback.setText("Choose an installed cgMLST scheme, or clear the cgMLST "
                                      "option. Nothing was run.")
                return
            verdict = cgmlst_scheme_applies(self.samples, cgmlst_scheme)
            if not verdict["applies"]:
                self.feedback.setText(verdict["message"] + " Nothing was run.")
                return
            if verdict["warnings"]:
                self.feedback.setText(" ".join([verdict["message"], *verdict["warnings"]]))
        choice = self.database_choice.currentData()
        self.plan = {
            "fastqc": self.fastqc.isChecked(),
            "hydra": self.hydra.isChecked(),
            "cgmlst": self.cgmlst.isChecked(),
            "cgmlst_scheme": cgmlst_scheme,
            "assemble": self.assemble.isChecked(),
            "memory_gb": self.memory.value(),
            "db_root": self.database.text(),
            "databases": [choice] if choice else None,
            "protein": self.protein.isChecked(),
            "point_mutations": self.point_mutations.isChecked(),
            "virulence": self.virulence.currentData(),
            "thresholds": {key: field.value() for key, field in self.threshold_controls.items()},
            "threads": self.threads.value(),
            "resource_plan": allocation.to_dict(),
        }
        super().accept()
