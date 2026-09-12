"""Cross-project evidence research and explicit, sequence-free cohort reuse."""

import copy
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLineEdit, QMessageBox, QVBoxLayout

from wmlstudio.sample_workflow import current_input_sha256
from wmlstudio.ui_common import FlowLayout, cell, gene_names, make_table, organism_for
from wmlstudio.widgets import button, label


def reuse_library_profiles(project, snapshots):
    """Copy frozen evidence with origin IDs; never open or mutate a source project."""
    current = project.samples()
    origins = {(str(project.path), sample["id"]): sample["id"] for sample in current}
    for sample in current:
        origin = sample.get("metadata", {}).get("library_origin", {})
        if origin.get("project_path") and origin.get("sample_id"):
            origins[(origin["project_path"], origin["sample_id"])] = sample["id"]
    identifiers = []
    with project.transaction():
        for snapshot in snapshots:
            if snapshot.get('status') != 'completed':
                raise ValueError(f"{snapshot['name']} is not a completed current sample ({snapshot.get('status') or 'unknown state'}). Reanalyse or explicitly export reviewed historical evidence before reuse.")
            origin = (snapshot["project_path"], snapshot["id"])
            if origin in origins:
                identifiers.append(origins[origin])
                continue
            primary = copy.deepcopy(snapshot.get("result") or {})
            analyses = [copy.deepcopy(result) for result in snapshot.get("analyses", [])
                        if result.get("scheme_digest") and result.get("alleles")]
            current_hash = current_input_sha256(snapshot) or primary.get('input_sha256')
            analyses = [result for result in analyses if not current_hash or not result.get("input_sha256")
                        or result["input_sha256"] == current_hash]
            if current_hash and primary.get('input_sha256') and primary['input_sha256'] != current_hash:
                primary = {}
            if not primary.get("alleles"):
                primary = analyses[0] if analyses else {}
            if not primary:
                raise ValueError(f"{snapshot['name']} has no saved allelic profile to reuse.")
            metadata = copy.deepcopy(snapshot.get("metadata", {}))
            metadata["library_origin"] = {"project_path": origin[0], "sample_id": origin[1],
                "updated_at": snapshot.get("updated_at"), "source_input_sha256": current_hash,
                "meaning": "Frozen library evidence; source project and input files were not modified."}
            identifier = project.add_profile(snapshot["name"], primary, metadata)
            for analysis in analyses:
                project.set_analysis(identifier, analysis)
            for collection in snapshot.get("collections", []):
                collection_id = project.create_collection(collection["name"])
                project.set_collection_members(collection_id, [identifier], add=True)
            project.record_history(identifier, "library_evidence_reused", metadata["library_origin"])
            origins[origin] = identifier
            identifiers.append(identifier)
    return identifiers


class LibraryDialog(QDialog):
    profilesImported = Signal(list)

    def __init__(self, library, project, parent=None):
        super().__init__(parent)
        self.library, self.project = library, project
        self.matches = []
        self.setWindowTitle("Research the saved library")
        self.resize(1060, 650)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Research across saved projects", "title"))
        layout.addWidget(label("Search indexed sample evidence without reopening genomes. Select isolates to reuse their saved profiles in the active project; original projects stay unchanged.", "muted", True))
        self.query = QLineEdit()
        self.query.setPlaceholderText("Search any indexed name, date, QC metric, gene or metadata value…")
        self.query.returnPressed.connect(self.search)
        layout.addWidget(self.query)
        filters = FlowLayout()
        self.filters = {}
        for key, title in [("genus", "Genus"), ("species", "Species"), ("st", "Exact ST"),
                           ("amr_gene", "AMR gene"), ("collection", "Collection")]:
            entry = QLineEdit()
            entry.setPlaceholderText(title)
            entry.setMaximumWidth(180)
            entry.returnPressed.connect(self.search)
            self.filters[key] = entry
            filters.addWidget(entry)
        filters.addWidget(button("Search library", self.search, True))
        layout.addLayout(filters)
        self.table = make_table(["Sample", "Genus", "Species", "ST", "AMR genes", "Project", "Saved / updated"])
        layout.addWidget(self.table, 1)
        self.status = label("", "small", True)
        layout.addWidget(self.status)
        actions = QHBoxLayout()
        actions.addWidget(button("Reuse selected profiles…", self.reuse, True))
        actions.addWidget(button("Select visible", self.table.selectAll))
        actions.addStretch()
        actions.addWidget(button("Close", self.accept))
        layout.addLayout(actions)
        self.search()

    def search(self):
        self.matches = self.library.search(self.query.text().strip(),
            **{key: entry.text().strip() or None for key, entry in self.filters.items()})
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.matches))
        for row, sample in enumerate(self.matches):
            genus, species, _ = organism_for(sample)
            values = [sample["name"], genus or "Unknown", species or "—", (sample.get("result") or {}).get("st"),
                      "; ".join(gene_names(sample)), Path(sample["project_path"]).stem, sample.get("updated_at")]
            for column, value in enumerate(values):
                self.table.setItem(row, column, cell(value, row))
        self.table.setSortingEnabled(True)
        self.status.setText(f"{len(self.matches)} saved samples matched. This index reflects each project's last saved snapshot, not a live reanalysis.")

    def reuse(self):
        indexes = sorted({item.data(Qt.ItemDataRole.UserRole) for item in self.table.selectedItems()})
        if not indexes:
            self.status.setText("Select one or more rows first.")
            return
        if QMessageBox.question(self, "Reuse saved evidence", f"Bring {len(indexes)} saved profiles and their metadata into this project?\n\nNo sequence files are copied or reanalysed. Existing copies from the same origin are selected, not duplicated.") != QMessageBox.StandardButton.Yes:
            return
        try:
            identifiers = reuse_library_profiles(self.project, [self.matches[index] for index in indexes])
            self.profilesImported.emit(identifiers)
            self.accept()
        except Exception as exc:
            self.status.setText(str(exc))
