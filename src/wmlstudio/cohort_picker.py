"""One explicit, reusable native isolate selection dialog per analysis/report."""

from __future__ import annotations

from collections import defaultdict

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLineEdit,
    QMenu,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
)

from wmlstudio.archive import active_samples
from wmlstudio.context_menus import (
    SEPARATOR,
    action_state,
    install_context_menu,
    plan_for,
    tidy_plan,
)
from wmlstudio.journey import isolate_records
from wmlstudio.ui_common import cell, make_table, organism_for
from wmlstudio.widgets import button, label

# A chooser must not become a second place to edit the project. It offers the two
# actions a user genuinely needs while choosing — correct an organism folder, copy
# the identifiers — and leaves everything else to the workspace that owns it.
PICKER_HANDLERS = frozenset({"context_assign_organism", "context_copy_id"})


class CohortPickerDialog(QDialog):
    """Picking a cohort never modifies samples, profiles, or another view's scope."""

    def __init__(self, samples, project, title="Choose isolates", selected_ids=None,
                 parent=None, *, include_reads=True):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1060, 690)
        self.setMinimumSize(780, 480)
        # Archived isolates stay in the project with all their evidence; they are
        # simply not offered as a choice until the user restores them.
        self.samples = isolate_records(active_samples(samples))
        self._host = parent
        if not include_reads:
            self.samples = [sample for sample in self.samples
                            if (sample.get("result") or {}).get("kind") != "fastq"
                            and not sample.get("input_path", "").lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fq", ".fastq"))]
        valid = {sample["id"] for sample in self.samples}
        self.selected_ids = set(selected_ids or ()) & valid
        self._folder = None
        self._filling = False
        self._profiles = {}
        # The classical ST is read from the analysis stored as MLST, not from the
        # headline result: the headline is whichever analysis ran last, and a
        # core-genome run must not empty the seven-locus column.
        self._mlst = {}
        for sample in self.samples:
            primary = sample.get("result") or {}
            secondary = project.analysis_summaries(sample["id"]) if project else []
            for row in secondary:
                if row.get("typing_kind") == "mlst":
                    self._mlst[sample["id"]] = {"scheme": row.get("scheme") or "", "st": row.get("st")}
            self._profiles[sample["id"]] = sorted({str(profile.get("scheme") or "unnamed")
                for profile in secondary if profile.get("scheme_digest") != primary.get("scheme_digest")
                and profile.get("input_sha256") and primary.get("input_sha256")
                and sample.get("status") == "completed"
                and primary["input_sha256"] == profile["input_sha256"]})
        layout = QVBoxLayout(self)
        layout.addWidget(label(title, "title", True))
        layout.addWidget(label("Choose an organism folder, then check all or individual isolates. This selection belongs only to the action you are starting.", "muted", True))
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find an isolate, date, organism, ST or saved cgMLST scheme…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_rows)
        layout.addWidget(self.search)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.folders = QTreeWidget()
        self.folders.setHeaderLabel("Organism folders")
        self.folders.setMinimumWidth(180)
        everything = QTreeWidgetItem([f"All isolates ({len(self.samples)})"])
        everything.setData(0, Qt.ItemDataRole.UserRole, None)
        self.folders.addTopLevelItem(everything)
        groups = defaultdict(lambda: defaultdict(int))
        for sample in self.samples:
            genus, species, _ = organism_for(sample)
            groups[genus or "Unknown"][species or "Unspecified"] += 1
        for genus in sorted(groups):
            parent_item = QTreeWidgetItem([f"{genus} ({sum(groups[genus].values())})"])
            parent_item.setData(0, Qt.ItemDataRole.UserRole, (genus,))
            everything.addChild(parent_item)
            for species, count in sorted(groups[genus].items()):
                item = QTreeWidgetItem([f"{species} ({count})"])
                item.setData(0, Qt.ItemDataRole.UserRole, (genus, species))
                parent_item.addChild(item)
        self.folders.expandAll()
        self.folders.itemClicked.connect(self.choose_folder)
        install_context_menu(self.folders, "picker.folders", self.show_context_menu,
                             dialect="folders", folder_ids=self.folder_members)
        splitter.addWidget(self.folders)
        # "ST (7-locus)" is read from the classical profile alone: a core-genome run
        # is a different measurement and must never appear in this column.
        self.table = make_table(["Include", "Isolate", "Collection date", "Organism",
                                 "Organism review", "ST (7-locus)", "Additional saved profiles",
                                 "Quality / state"])
        self.table.setColumnWidth(0, 60)
        self.table.itemChanged.connect(self.item_changed)
        self.table.cellDoubleClicked.connect(self.toggle_row)
        install_context_menu(self.table, "picker.table", self.show_context_menu, id_column=0)
        splitter.addWidget(self.table)
        splitter.setSizes([225, 810])
        layout.addWidget(splitter, 1)
        controls = QHBoxLayout()
        controls.addWidget(button("Select shown", self.select_shown))
        controls.addWidget(button("Deselect shown", self.deselect_shown))
        controls.addWidget(button("Clear all", self.clear_all))
        self.count = label("", "small", True)
        controls.addWidget(self.count, 1)
        layout.addLayout(controls)
        self.feedback = label("Organism folders show assigned or provisional labels; inspect identity evidence before interpretation.", "small", True)
        layout.addWidget(self.feedback)
        self.controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.controls.button(QDialogButtonBox.StandardButton.Ok).setText("Use selected isolates")
        self.controls.accepted.connect(self.accept)
        self.controls.rejected.connect(self.reject)
        layout.addWidget(self.controls)
        self.refresh_rows()

    def choose_folder(self, item, column=0):
        self._folder = item.data(0, Qt.ItemDataRole.UserRole)
        self.refresh_rows()

    def folder_members(self, genus, species=""):
        """The isolates an organism folder node stands for; folder nodes carry no ids."""
        wanted = (genus, species) if species else (genus,)
        members = []
        for sample in self.samples:
            found, child, _ = organism_for(sample)
            taxon = (found or "Unknown", child or "Unspecified")
            if taxon[:len(wanted)] == wanted:
                members.append(sample["id"])
        return members

    def show_context_menu(self, selection, position):
        """Offer the owning window's own handlers, from a menu this dialog owns.

        The menu is parented to the dialog so it works while the chooser is modal,
        and the actions are the window's, so there is one implementation of each.
        """
        host = self._host
        if host is None:
            return None
        entries = tidy_plan([entry for entry in plan_for(selection.view_id, selection)
                             if entry is SEPARATOR or (entry.handler in PICKER_HANDLERS
                                                       and callable(getattr(host, entry.handler, None)))])
        if not entries:
            return None
        worker = getattr(host, "worker", None)
        busy = bool(worker is not None and worker.isRunning())
        chosen_ids = set(selection.sample_ids)
        samples = [sample for sample in self.samples if sample["id"] in chosen_ids]
        menu = QMenu(self)
        menu.setToolTipsVisible(True)
        for entry in entries:
            if entry is SEPARATOR:
                menu.addSeparator()
                continue
            enabled, reason = action_state(entry, selection, window=host, busy=busy, samples=samples)
            action = menu.addAction(entry.format_title(selection))
            action.setEnabled(enabled)
            if reason:
                action.setToolTip(reason)
            action.setData((entry.handler, None))
        chosen = menu.exec(position)
        data = chosen.data() if chosen is not None else None
        menu.deleteLater()
        if not data:
            return None
        result = host.run_context_action(data, selection)
        if data[0] != "context_copy_id":
            self.feedback.setText("Organism changes are saved to the project and the managed copies "
                                  "move with them. Reopen this chooser to see the new folders.")
        return result

    def refresh_rows(self):
        from wmlstudio.ui_workbench import organism_review, paint_review
        self._filling = True
        self.table.setSortingEnabled(False)
        self.table.blockSignals(True)
        query = self.search.text().strip().casefold()
        rows = []
        for sample in self.samples:
            genus, species, _ = organism_for(sample)
            taxon = (genus or "Unknown", species or "Unspecified")
            if self._folder and taxon[:len(self._folder)] != tuple(self._folder):
                continue
            metadata = sample.get("metadata") or {}
            result = sample.get("result") or {}
            mlst = self._mlst.get(sample["id"], {})
            date = (metadata.get("annotations") or {}).get("collection_date") or metadata.get("collection_date") or metadata.get("isolation_date") or metadata.get("date") or "Not recorded"
            # An organism nobody has settled travels with the isolate into whatever
            # is started from here, so the chooser says so before it is chosen.
            review = organism_review(sample, None, mlst)
            values = ["", sample["name"], date, " ".join(taxon),
                      review["label"] if review["needs_review"] else "Ready",
                      mlst.get("st") if mlst.get("st") is not None else "Unassigned",
                      "; ".join(self._profiles[sample["id"]]) or "Not called",
                      result.get("status") or sample.get("status")]
            if query and query not in " ".join(map(str, values)).casefold():
                continue
            rows.append((sample["id"], values, review))
        self.table.setRowCount(len(rows))
        for row, (identifier, values, review) in enumerate(rows):
            for column, value in enumerate(values):
                item = cell(value, identifier)
                if column == 0:
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked if identifier in self.selected_ids else Qt.CheckState.Unchecked)
                if column in (3, 4):
                    paint_review(item, review)
                    if not review["needs_review"]:
                        item.setToolTip(review["detail"])
                self.table.setItem(row, column, item)
        self.table.blockSignals(False)
        self.table.setSortingEnabled(True)
        self._filling = False
        self.update_count()

    def visible_ids(self):
        return {self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) for row in range(self.table.rowCount())}

    def update_count(self):
        hidden = len(self.selected_ids - self.visible_ids())
        self.count.setText(f"{len(self.selected_ids)} selected · {self.table.rowCount()} shown" +
                           (f" · {hidden} selected outside this filter" if hidden else ""))

    def item_changed(self, item):
        if self._filling or item.column() != 0:
            return
        identifier = item.data(Qt.ItemDataRole.UserRole)
        if item.checkState() == Qt.CheckState.Checked:
            self.selected_ids.add(identifier)
        else:
            self.selected_ids.discard(identifier)
        self.update_count()

    def toggle_row(self, row, column):
        if column != 0:
            item = self.table.item(row, 0)
            item.setCheckState(Qt.CheckState.Unchecked if item.checkState() == Qt.CheckState.Checked else Qt.CheckState.Checked)

    def select_shown(self):
        self.selected_ids.update(self.visible_ids())
        self.refresh_rows()

    def deselect_shown(self):
        self.selected_ids.difference_update(self.visible_ids())
        self.refresh_rows()

    def clear_all(self):
        self.selected_ids.clear()
        self.refresh_rows()

    def accept(self):
        if not self.selected_ids:
            self.feedback.setText("Choose at least one isolate, or cancel. No project records have been changed.")
            return
        super().accept()
