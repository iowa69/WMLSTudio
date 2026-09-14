"""Native assignment forms; filesystem work is left to the caller's worker."""

from datetime import UTC, datetime
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QMenu,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
)

from wmlstudio.ui_common import FlowLayout

# What produced a label, in words a microbiologist reads rather than a token.
# Every one of them is a statement about an installed reference panel, never a
# statement that the organism was identified in a laboratory.
BASIS_LABELS = {
    "genomic_ani": "Whole-genome comparison",
    "genomic_ani_kpsc": "Whole-genome comparison (focused panel)",
    "mlst_panel": "Typing-panel match",
    "user_assigned": "You assigned it",
    "csv_import": "From your spreadsheet",
    "none": "Not identified",
}


class ImportSamplesDialog(QDialog):
    """Return per-file assignments and copy options after the user accepts.

    ``scheme_entries`` is a list of (display name, path) tuples or paths.
    ``assignments`` and ``options`` contain plain Python data, safe to give to a
    worker. This dialog does not read sequences, load schemes, or copy files.
    """

    def __init__(self, paths, scheme_entries=(), parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import samples and assign organisms")
        self.resize(900, 680)
        self.assignments = [{"path": str(path), "typing_mode": "auto", "genus": "",
                             "species": "", "scheme_path": None} for path in paths]
        self.options = {"managed": True, "storage_root": "", "append_st": False}
        layout = QVBoxLayout(self)
        description = QLabel(
            "Assign an organism to each sample, or identify it from your local schemes. "
            "Select several rows to apply the same choice together."
        )
        description.setWordWrap(True)
        layout.addWidget(description)
        self.table = QTableWidget(len(self.assignments), 4)
        self.table.setHorizontalHeaderLabels(["Sample / file", "Assignment", "Organism", "Scheme"])
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table, 1)
        form = QFormLayout()
        self.mode = QComboBox()
        self.mode.addItem("Identify from local schemes", "auto")
        self.mode.addItem("Assign an organism", "manual")
        self.mode.addItem("Unknown organism / AMR only", "unknown")
        form.addRow("Organism choice", self.mode)
        self.genus = QLineEdit()
        self.genus.setPlaceholderText("For example, Klebsiella")
        self.species = QLineEdit()
        self.species.setPlaceholderText("For example, pneumoniae; optional if genus only")
        form.addRow("Genus", self.genus)
        form.addRow("Species", self.species)
        self.scheme = QComboBox()
        self.scheme.addItem("Choose from organism / local detection", None)
        for entry in scheme_entries:
            if isinstance(entry, (tuple, list)) and len(entry) == 2:
                name, path = entry
            else:
                path = Path(entry)
                name = path.name
            self.scheme.addItem(str(name), str(path))
        form.addRow("Typing scheme", self.scheme)
        layout.addLayout(form)
        buttons = QHBoxLayout()
        apply_selected = QPushButton("Apply to selected rows")
        apply_selected.clicked.connect(self.apply_selected)
        apply_all = QPushButton("Apply to all rows")
        apply_all.clicked.connect(self.apply_all)
        buttons.addWidget(apply_selected)
        buttons.addWidget(apply_all)
        buttons.addStretch()
        layout.addLayout(buttons)
        self.managed = QCheckBox("Keep managed copies organised by genus, species and ST")
        self.managed.setChecked(True)
        layout.addWidget(self.managed)
        storage = QHBoxLayout()
        self.storage_root = QLineEdit()
        self.storage_root.setPlaceholderText("Default: Sequences folder beside the project")
        storage.addWidget(self.storage_root, 1)
        browse = QPushButton("Choose storage folder…")
        browse.clicked.connect(self.browse_storage)
        storage.addWidget(browse)
        layout.addLayout(storage)
        self.append_st = QCheckBox("Append the assigned ST to managed filenames")
        layout.addWidget(self.append_st)
        self.feedback = QLabel("Original sequence files stay unchanged.")
        self.feedback.setWordWrap(True)
        layout.addWidget(self.feedback)
        actions = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        actions.accepted.connect(self.accept)
        actions.rejected.connect(self.reject)
        layout.addWidget(actions)
        self.mode.currentIndexChanged.connect(self.update_form)
        self.update_form()
        self.refresh_table()

    def update_form(self):
        manual = self.mode.currentData() == "manual"
        self.genus.setEnabled(manual)
        self.species.setEnabled(manual)
        self.scheme.setEnabled(self.mode.currentData() != "unknown")

    def browse_storage(self):
        path = QFileDialog.getExistingDirectory(self, "Choose managed sequence storage")
        if path:
            self.storage_root.setText(path)

    def apply_selected(self):
        rows = sorted({index.row() for index in self.table.selectionModel().selectedRows()})
        if not rows:
            self.feedback.setText("Select one or more sample rows first.")
            return
        self._apply(rows)

    def apply_all(self):
        self._apply(range(len(self.assignments)))

    def _apply(self, rows):
        mode = self.mode.currentData()
        genus = self.genus.text().strip() if mode == "manual" else ""
        species = self.species.text().strip() if mode == "manual" else ""
        if mode == "manual" and not genus:
            self.feedback.setText("Enter a genus for manual organism assignment.")
            return
        for row in rows:
            self.assignments[row].update({"typing_mode": mode, "genus": genus, "species": species,
                                          "scheme_path": self.scheme.currentData() if mode != "unknown" else None})
        self.feedback.setText("Assignment applied. Original sequence files stay unchanged.")
        self.refresh_table()

    def refresh_table(self):
        labels = {"auto": "Identify locally", "manual": "Assigned", "unknown": "Unknown / AMR only"}
        for row, assignment in enumerate(self.assignments):
            cells = [assignment.get("name") or Path(assignment["path"]).name,
                     labels.get(assignment["typing_mode"], assignment["typing_mode"]),
                     " ".join(filter(None, [assignment.get("genus"), assignment.get("species")])),
                     Path(assignment["scheme_path"]).name if assignment.get("scheme_path") else "Automatic"]
            for column, value in enumerate(cells):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, row)
                item.setToolTip(assignment["path"] if column == 0 else value)
                self.table.setItem(row, column, item)

    def accept(self):
        if not self.assignments:
            self.feedback.setText("No sequence files were selected.")
            return
        self.options = {"managed": self.managed.isChecked(), "storage_root": self.storage_root.text().strip(),
                        "append_st": self.append_st.isChecked()}
        super().accept()


class IdentificationReviewDialog(QDialog):
    """Review every proposed organism before one byte is copied anywhere.

    The dialog does no filesystem work: it edits verdicts and hands back the
    assignments ``storage.import_samples`` consumes. Its "Goes to" column is
    rendered by ``storage.preview_target`` — the same function the import itself
    calls — so the folder shown here and the folder created later cannot differ.

    Nothing here is an identification. A folder name is where a copy is stored,
    and a proposal nobody accepted is filed under Needs review rather than pushed
    into the genus somebody expected.
    """

    FOOTER = ("A folder name is where the file is stored. It is not a laboratory "
              "identification. Confirm the organism before clinical interpretation.")
    ACTION_CHOICES = (("Accept this organism", "accept"),
                      ("Send to Needs review", "quarantine"),
                      ("Do not import this file", "skip"))
    COLUMNS = ("File", "Genus", "Species", "Basis", "Confidence", "Goes to", "What happens")
    PLACEHOLDER = "new-isolate"

    def __init__(self, verdicts, storage_root, parent=None, *, organisms=(), names=None,
                 policy=None):
        super().__init__(parent)
        from wmlstudio.organism_id import AUTO_FILE_MIN, confidence_rank, meets_policy
        self.setWindowTitle("Review automatic organism identification")
        self.resize(1020, 660)
        self.setMinimumSize(840, 520)
        self.verdicts = [dict(verdict) for verdict in verdicts]
        self.root = Path(storage_root)
        self.policy = policy
        self.names = dict(names or {})
        self.assignments = []
        self.rows = []
        self._filling = True
        layout = QVBoxLayout(self)
        heading = QLabel(f"{len(self.verdicts)} files were examined on disk, where they are. "
                         "Nothing has been copied yet. Choose what to do with each file; you can "
                         "change any organism now or later.")
        heading.setWordWrap(True)
        layout.addWidget(heading)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.table = QTableWidget(len(self.verdicts), len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(list(self.COLUMNS))
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.itemSelectionChanged.connect(self.show_evidence)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_row_menu)
        splitter.addWidget(self.table)
        self.evidence = QTextBrowser()
        self.evidence.setMinimumHeight(110)
        splitter.addWidget(self.evidence)
        splitter.setSizes([420, 160])
        layout.addWidget(splitter, 1)
        genera, species = self._suggestions(organisms)
        for row, verdict in enumerate(self.verdicts):
            proposed = verdict.get("proposed") or {}
            genus_box = self._combo(genera, proposed.get("genus", ""))
            species_box = self._combo(species, proposed.get("species", ""))
            action_box = QComboBox()
            for title, value in self.ACTION_CHOICES:
                action_box.addItem(title, value)
            strong = (confidence_rank(verdict.get("confidence")) >= confidence_rank(AUTO_FILE_MIN)
                      and bool(proposed.get("genus")))
            chosen = "accept" if strong or meets_policy(verdict, policy) else "quarantine"
            action_box.setCurrentIndex(action_box.findData(chosen))
            for box in (genus_box, species_box, action_box):
                box.setMinimumWidth(112)
            genus_box.currentTextChanged.connect(lambda _text, r=row: self.update_row(r))
            species_box.currentTextChanged.connect(lambda _text, r=row: self.update_row(r))
            action_box.currentIndexChanged.connect(lambda _index, r=row: self.update_row(r))
            self.rows.append((genus_box, species_box, action_box))
            self.table.setCellWidget(row, 1, genus_box)
            self.table.setCellWidget(row, 2, species_box)
            self.table.setCellWidget(row, 6, action_box)
        bulk = FlowLayout()
        for title, tip, handler in (
                ("Accept all strong", "Accept only the proposals that met every configured gate.",
                 self.accept_strong),
                ("Accept all proposals", "Accept every proposed organism, including the ones the "
                 "evidence does not support strongly. Each label keeps its own basis.",
                 self.accept_proposed),
                ("Send the rest to Needs review", "Import the undecided files without filing them "
                 "under an organism.", self.quarantine_rest),
                ("Save as spreadsheet…", "Write this list out to edit the organisms elsewhere.",
                 self.save_spreadsheet),
                ("Load from spreadsheet…", "Read organisms back in. A spreadsheet is your "
                 "statement, never genomic evidence.", self.load_spreadsheet)):
            item = QPushButton(title)
            item.setToolTip(tip)
            item.clicked.connect(handler)
            bulk.addWidget(item)
        layout.addLayout(bulk)
        self.feedback = QLabel("")
        self.feedback.setWordWrap(True)
        layout.addWidget(self.feedback)
        footer = QLabel(self.FOOTER)
        footer.setWordWrap(True)
        layout.addWidget(footer)
        self.controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                         | QDialogButtonBox.StandardButton.Cancel)
        self.controls.button(QDialogButtonBox.StandardButton.Ok).setText("Import the accepted files")
        self.controls.button(QDialogButtonBox.StandardButton.Cancel).setText("Cancel · import nothing")
        self.controls.accepted.connect(self.accept)
        self.controls.rejected.connect(self.reject)
        layout.addWidget(self.controls)
        self._filling = False
        self.refresh_table()
        if self.verdicts:
            self.table.selectRow(0)

    # --- row model ----------------------------------------------------------
    @staticmethod
    def _combo(values, current):
        box = QComboBox()
        box.setEditable(True)
        box.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        box.addItem("")
        for value in values:
            box.addItem(value)
        box.setCurrentText(current or "")
        return box

    def _suggestions(self, organisms):
        genera, species = set(), set()
        for entry in organisms or ():
            genus = str(entry[0] or "").strip() if isinstance(entry, (tuple, list)) else str(entry).strip()
            child = str(entry[1] or "").strip() if isinstance(entry, (tuple, list)) and len(entry) > 1 else ""
            if genus:
                genera.add(genus)
            if child:
                species.add(child)
        for verdict in self.verdicts:
            proposed = verdict.get("proposed") or {}
            if proposed.get("genus"):
                genera.add(str(proposed["genus"]))
            if proposed.get("species"):
                species.add(str(proposed["species"]))
        return sorted(genera), sorted(species)

    def row_action(self, row):
        return self.rows[row][2].currentData()

    def row_verdict(self, row):
        """The verdict as the user has it now: their organism, their decision."""
        verdict = dict(self.verdicts[row])
        genus = self.rows[row][0].currentText().strip()
        species = self.rows[row][1].currentText().strip()
        proposed = verdict.get("proposed") or {}
        verdict["proposed"] = {"genus": genus, "species": species}
        if (genus, species) != (proposed.get("genus", ""), proposed.get("species", "")):
            # A label typed over the proposal is the operator's statement, so it
            # does not inherit the confidence of the evidence it replaced.
            verdict.update(basis="user_assigned", confidence="unresolved")
        if self.row_action(row) == "accept":
            verdict.update(status="confirmed", accepted={"genus": genus, "species": species},
                           quarantine_reason=None, confirmed_by="user",
                           confirmed_utc=datetime.now(UTC).isoformat())
        else:
            # Why a file could not be identified is a fact about the evidence. It
            # is only "you deferred it" when there was a proposal to defer.
            reason = None if genus else verdict.get("quarantine_reason")
            verdict.update(status="quarantined", quarantine_reason=reason or "user_deferred",
                           accepted={"genus": "", "species": ""}, confirmed_by=None,
                           confirmed_utc=None)
        return verdict

    def destination_path(self, row):
        """The folder this row's file would be copied into, computed by storage itself."""
        from wmlstudio.organism_id import proposed_destination
        from wmlstudio.storage import preview_target
        verdict = self.row_verdict(row)
        quarantine, genus, species = proposed_destination(verdict, policy=self.policy)
        target = preview_target(self.root, self.PLACEHOLDER, Path(verdict["input_path"]),
                                genus, species, quarantine=quarantine)
        return target.parent.parent

    def destination_text(self, row):
        return self.destination_path(row).relative_to(self.root.expanduser().resolve()).as_posix() + "/"

    # --- rendering ----------------------------------------------------------
    def refresh_table(self):
        for row in range(len(self.verdicts)):
            self.update_row(row)

    def update_row(self, row):
        if self._filling:
            return
        from wmlstudio.organism_id import confidence_label
        verdict = self.verdicts[row]
        path = Path(verdict.get("input_path") or "")
        word, explanation = confidence_label(verdict)
        basis = BASIS_LABELS.get(verdict.get("basis"), str(verdict.get("basis") or "Not identified"))
        skipped = self.row_action(row) == "skip"
        goes_to = "Not imported" if skipped else self.destination_text(row)
        values = {0: self.names.get(str(path)) or path.name, 3: basis, 4: word, 5: goes_to}
        for column, value in values.items():
            item = QTableWidgetItem(value)
            item.setData(Qt.ItemDataRole.UserRole, row)
            item.setToolTip({3: str(verdict.get("reason") or basis), 4: explanation,
                             5: "A folder for this isolate is created inside this one. "
                                "Your original file stays where it is.",
                             0: str(path)}.get(column, value))
            self.table.setItem(row, column, item)
        self.update_feedback()

    def update_feedback(self):
        actions = [self.row_action(row) for row in range(len(self.verdicts))]
        accepted = actions.count("accept")
        review = actions.count("quarantine")
        skipped = actions.count("skip")
        self.feedback.setText(
            f"{accepted} will be filed under the organism shown · {review} will be imported into "
            f"Needs review, where you can assign an organism later · {skipped} will not be imported.")

    def show_evidence(self):
        rows = {index.row() for index in self.table.selectionModel().selectedRows()}
        if not rows:
            self.evidence.setPlainText("Select a row to see why this organism was proposed.")
            return
        from wmlstudio.organism_id import confidence_label, evidence_sentences
        row = min(rows)
        verdict = self.verdicts[row]
        word, explanation = confidence_label(verdict)
        lines = [f"Why this organism? · {word}", explanation, ""]
        for hit in (verdict.get("detail") or {}).get("ani_top") or ():
            lines.append(f"  {hit.get('genus', '')} {hit.get('species', '')} · ANI {hit.get('ani')} · "
                         f"aligned {hit.get('query_fraction')} of this assembly, "
                         f"{hit.get('reference_fraction')} of the reference")
        for candidate in (verdict.get("detail") or {}).get("mlst_top") or ():
            lines.append(f"  {candidate.get('scheme')} · {candidate.get('matched_loci')} loci matched · "
                         f"coverage {candidate.get('coverage')} · {candidate.get('organism_label') or 'no organism recorded'}")
        if verdict.get("margin_ani") is not None:
            lines.append(f"  Separation from the next species: {verdict['margin_ani']}")
        lines.append("")
        lines.extend(evidence_sentences(verdict))
        for error in verdict.get("errors") or ():
            lines.append(f"{error.get('tier')}: {error.get('message')}")
        self.evidence.setPlainText("\n".join(line for line in lines if line is not None))

    def selected_rows(self):
        return sorted({index.row() for index in self.table.selectionModel().selectedRows()})

    def show_row_menu(self, point):
        """Decide several rows at once without leaving the table."""
        item = self.table.itemAt(point)
        if item is not None and item.row() not in self.selected_rows():
            self.table.selectRow(item.row())
        rows = self.selected_rows()
        if not rows:
            return
        menu = QMenu(self)
        for title, value in self.ACTION_CHOICES:
            entry = menu.addAction(f"{title} · {len(rows)} rows" if len(rows) > 1 else title)
            entry.setData(value)
        menu.addSeparator()
        menu.addAction("Copy file path").setData("copy")
        menu.addAction("Why this organism?").setData("why")
        chosen = menu.exec(self.table.viewport().mapToGlobal(point))
        value = chosen.data() if chosen is not None else None
        menu.deleteLater()
        if value in {"accept", "quarantine", "skip"}:
            self._set_action(rows, value)
        elif value == "copy":
            QApplication.clipboard().setText(
                "\n".join(str(self.verdicts[row].get("input_path") or "") for row in rows))
            self.feedback.setText(f"Copied {len(rows)} file paths. Your files were not moved.")
        elif value == "why":
            self.show_evidence()

    # --- bulk decisions -----------------------------------------------------
    def _set_action(self, rows, value):
        for row in rows:
            box = self.rows[row][2]
            box.setCurrentIndex(box.findData(value))

    def accept_strong(self):
        from wmlstudio.organism_id import AUTO_FILE_MIN, confidence_rank
        floor = confidence_rank(AUTO_FILE_MIN)
        rows = [row for row, verdict in enumerate(self.verdicts)
                if confidence_rank(verdict.get("confidence")) >= floor
                and self.rows[row][0].currentText().strip()]
        self._set_action(rows, "accept")
        self.feedback.setText(f"{len(rows)} strong proposals accepted. Weaker rows were left for you "
                              "to decide; nothing was strengthened.")

    def accept_proposed(self):
        rows = [row for row in range(len(self.verdicts)) if self.rows[row][0].currentText().strip()]
        self._set_action(rows, "accept")
        self.feedback.setText(f"{len(rows)} proposals accepted, including proposals the evidence does "
                              "not support strongly. The basis of each label is kept with the isolate.")

    def quarantine_rest(self):
        rows = [row for row in range(len(self.verdicts)) if self.row_action(row) != "accept"]
        self._set_action(rows, "quarantine")
        self.feedback.setText(f"{len(rows)} files will be imported into Needs review. Nothing is lost; "
                              "assign an organism whenever you are ready.")

    # --- spreadsheet override ----------------------------------------------
    def save_spreadsheet(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save this list as a spreadsheet",
                                              "organism-assignments.csv", "CSV (*.csv)")
        if not path:
            return
        from wmlstudio.organism_assignments import write_template
        # The file path is the only identity a not-yet-imported row has, so the
        # name column is left empty: two identity columns would be ambiguous.
        targets = [{"file": verdict.get("input_path", ""),
                    "genus": self.rows[row][0].currentText().strip(),
                    "species": self.rows[row][1].currentText().strip(),
                    "quarantine": "" if self.rows[row][0].currentText().strip() else "user_deferred"}
                   for row, verdict in enumerate(self.verdicts)]
        try:
            write_template(path, targets)
        except (OSError, ValueError) as error:
            self.feedback.setText(str(error))
            return
        self.feedback.setText(f"Saved to {path}. Edit the genus and species columns, then load it back.")

    def load_spreadsheet(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load organisms from a spreadsheet", "",
                                              "Tables (*.csv *.tsv *.txt);;All files (*)")
        if not path:
            return
        from wmlstudio.organism_assignments import read_assignments
        try:
            rows, problems = read_assignments(path, paths=[v.get("input_path") for v in self.verdicts])
        except (OSError, ValueError) as error:
            self.feedback.setText(str(error))
            return
        applied = 0
        by_path = {str(Path(verdict.get("input_path") or "").expanduser().resolve()): row
                   for row, verdict in enumerate(self.verdicts)}
        for entry in rows:
            row = by_path.get(str(Path(entry.get("path") or "").expanduser().resolve()))
            if row is None:
                continue
            self.rows[row][0].setCurrentText(entry.get("genus", ""))
            self.rows[row][1].setCurrentText(entry.get("species", ""))
            self._set_action([row], "accept" if entry.get("genus") else "quarantine")
            applied += 1
        self.feedback.setText(f"{applied} rows updated from the spreadsheet. "
                              + (f"{len(problems)} rows could not be matched: " + "; ".join(problems[:3])
                                 if problems else "A spreadsheet is your statement, not genomic evidence."))

    def accept(self):
        from wmlstudio.organism_id import assignment_for
        assignments = []
        for row, verdict in enumerate(self.verdicts):
            if self.row_action(row) == "skip":
                continue
            assignments.append(assignment_for(self.row_verdict(row), policy=self.policy,
                                              name=self.names.get(str(verdict.get("input_path")))))
        if not assignments:
            self.feedback.setText("Every file is set to 'Do not import'. Choose at least one file, "
                                  "or cancel; nothing has been copied.")
            return
        self.assignments = assignments
        super().accept()


class PracticeCohortDialog(QDialog):
    """Choose a pinned teaching cohort, having read what it is and is not.

    The caveats are rendered exactly as ``practice_cohorts`` states them: these
    are public reference assemblies for practice, never a validation set and
    never an expected answer key.
    """

    def __init__(self, cohorts, destinations, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Download practice data")
        self.resize(820, 560)
        self.setMinimumSize(700, 480)
        self.cohorts = list(cohorts)
        self.destinations = dict(destinations or {})
        self.chosen = self.cohorts[0]["name"] if self.cohorts else ""
        layout = QVBoxLayout(self)
        heading = QLabel("Practice datasets are downloaded from NCBI RefSeq on this computer and "
                         "verified against pinned checksums. They are teaching material: work "
                         "through them exactly as you would your own isolates.")
        heading.setWordWrap(True)
        layout.addWidget(heading)
        self.list = QListWidget()
        for entry in self.cohorts:
            megabytes = entry["download_bytes"] / (1024 * 1024)
            self.list.addItem(f"{entry['title']} · {entry['genomes']} genomes · {megabytes:.1f} MB")
        self.list.currentRowChanged.connect(self.show_cohort)
        layout.addWidget(self.list)
        self.detail = QTextBrowser()
        layout.addWidget(self.detail, 1)
        self.controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                         | QDialogButtonBox.StandardButton.Cancel)
        self.controls.button(QDialogButtonBox.StandardButton.Ok).setText("Download")
        self.controls.accepted.connect(self.accept)
        self.controls.rejected.connect(self.reject)
        layout.addWidget(self.controls)
        if self.cohorts:
            self.list.setCurrentRow(0)

    def show_cohort(self, row):
        if not 0 <= row < len(self.cohorts):
            return
        entry = self.cohorts[row]
        self.chosen = entry["name"]
        destination = self.destinations.get(entry["name"], "")
        lines = [entry["title"], "", entry["purpose"], "", entry.get("expect", ""), "",
                 f"{entry['genomes']} genomes across {len(entry['genera'])} genera: "
                 + ", ".join(entry["genera"]), "",
                 f"Downloads to: {destination}", "", "Before you use these files:"]
        lines.extend(f"  • {caveat}" for caveat in entry["caveats"])
        self.detail.setPlainText("\n".join(str(line) for line in lines))

    def accept(self):
        if not self.chosen:
            return
        super().accept()


class BatchAssignmentDialog(ImportSamplesDialog):
    def __init__(self, samples, scheme_entries=(), parent=None):
        samples = list(samples)
        super().__init__([sample.get("input_path", "") for sample in samples], scheme_entries, parent)
        self.setWindowTitle("Assign organisms to selected samples")
        self.managed.setChecked(False)
        for assignment, sample in zip(self.assignments, samples, strict=True):
            metadata = sample.get("metadata", {})
            workflow = metadata.get("workflow", {})
            organism = metadata.get("organism", {})
            assignment.update({"sample_id": sample["id"], "name": sample["name"],
                               "typing_mode": workflow.get("typing_mode", "auto"),
                               "genus": organism.get("genus", ""), "species": organism.get("species", ""),
                               "scheme_path": workflow.get("scheme_path")})
        self.refresh_table()
