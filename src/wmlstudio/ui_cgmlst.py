"""The dedicated cgMLST page: the table of calls, and a read recheck of uncalled targets.

One row per target, one column per isolate, and the scheme's own target count beside
every number on the page. A count of called targets without the size of the scheme it
was counted over is the misreading this page exists to prevent, so the denominator is
never left to the reader to supply.

Every state a call can be in keeps its own cell text, its own colour and its own
sentence: an assigned reference allele, a validated local novel sequence, and each
separate reason there is no allele (missing, partial, below the similarity thresholds,
ambiguous, duplicated, mixed, unrecognised). None of them is folded into another, and
nothing without an allele is ever counted as a call.

The recheck button asks the isolate's own reads about the targets its assembly did not
call, and reports what they show: the target is present, the reads are consistent with
it being absent, or there is not enough evidence to say. That is evidence about this
assembly. It assigns no allele and writes nothing back onto a profile, a distance or a
tree.

All of the science is in :mod:`wmlstudio.cgmlst_calls`; this module renders what that
module returns and never re-words its statements.
"""

from __future__ import annotations

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QComboBox,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.cgmlst_calls import (
    STATE_ORDER,
    build_calls_table,
    plan_target_recheck,
    recheck_missing_targets,
)
from wmlstudio.identification import cached_scheme
from wmlstudio.investigation import TYPING_SCALES
from wmlstudio.ui_common import FlowLayout, cell, make_table
from wmlstudio.widgets import button, label

# One background per call state, dark enough for this theme and separable enough to
# scan a scheme-length column. Colour is only an aid: every cell also carries the
# state's own words in its text and its tooltip, so nothing here depends on seeing it.
# 'called' deliberately has none — the ordinary case must not be the loudest.
STATE_TINTS = {
    "called": "",
    "called_novel": "#1D3346",
    "missing": "#3A2230",
    "partial": "#3A2E1C",
    "low_similarity": "#332A1E",
    "ambiguous": "#2C2440",
    "duplicated": "#262046",
    "mixed": "#3B2334",
    "unknown": "#25303F",
}

# Which rows of the scheme the grid shows. The filter never changes a number: the
# denominator printed beside the grid stays the scheme's full target count whatever
# is on screen, and the row count says how much of it is being shown.
ROW_FILTERS = (
    ("all", "Every target in the scheme"),
    ("incomplete", "Targets without a call in at least one isolate"),
    ("absent", "Targets called in no isolate"),
)

# What the recheck button promises, in one line above its results. The full wording
# travels with each result from cgmlst_calls and is printed unchanged beneath them.
RECHECK_PURPOSE = ("Investigative read evidence about this assembly. No allele is assigned, and no "
                   "profile, distance or tree is changed by anything on this tab.")


def cgmlst_references(project, sample_ids=None):
    """The cgMLST references a cohort actually has current stored profiles against.

    One table of calls covers one reference, because two references are two different
    target sets and their percentages are not the same quantity. Grouping here is by
    the stored scheme fingerprint, never by a scheme's name.
    """
    from wmlstudio.ui_compare import _available_profiles
    grouped = {}
    for sample in project.samples():
        if sample_ids is not None and sample["id"] not in sample_ids:
            continue
        for result in _available_profiles(project, sample, "cgmlst"):
            key = (str(result.get("scheme") or "Unnamed reference"),
                   str(result.get("scheme_digest") or ""))
            grouped.setdefault(key, []).append(
                dict(result, sample_id=sample["id"], sample_name=sample["name"]))
    references = []
    for (name, digest), records in grouped.items():
        if not digest:
            continue  # A profile with no reference fingerprint is not comparable with any other.
        references.append({"scheme": name, "scheme_digest": digest, "records": records,
                           "isolates": len(records),
                           "targets": max(len(row.get("alleles") or {}) for row in records)})
    references.sort(key=lambda row: (-row["isolates"], row["scheme"].casefold(), row["scheme_digest"]))
    return references


def installed_cgmlst_folders(scheme_paths):
    """Installed references a read check may be run against, cgMLST ones first.

    A folder whose typing kind was not recorded is still offered, labelled as such:
    the read check verifies the scheme's own fingerprint against the one these calls
    were made with, so a wrong choice is refused rather than silently answered.
    """
    from wmlstudio.reference_index import scheme_entries
    try:
        rows = scheme_entries(list(scheme_paths or ()))
    except (OSError, ValueError):
        return []
    offered = []
    for kind, suffix in (("cgmlst", ""), ("unknown", " · typing kind not recorded")):
        for row in rows:
            if row.get("kind") == kind:
                offered.append({"title": (row.get("title") or row.get("name") or row["id"]) + suffix,
                                "name": str(row.get("name") or ""), "path": row["path"]})
    return offered


class CallsTableModel(QAbstractTableModel):
    """The grid itself: rows are the scheme's targets, columns the cohort's isolates.

    Cells are rendered on demand rather than pre-built, because a cgMLST scheme is
    thousands of rows. Nothing here decides what a cell means: the text, the reason
    and the state all come from :func:`wmlstudio.cgmlst_calls.build_calls_table`.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.columns, self.rows, self.states = [], [], {}

    def set_table(self, table):
        """Show a whole table of calls; an empty mapping empties the grid."""
        self.beginResetModel()
        table = dict(table or {})
        self.columns = list(table.get("samples") or [])
        self.states = dict(table.get("states") or {})
        self.rows = list(table.get("loci") or [])
        self.endResetModel()

    def set_rows(self, rows):
        """Show a subset of the same table's rows; the columns and states are kept."""
        self.beginResetModel()
        self.rows = list(rows)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else (len(self.columns) + 1 if self.columns else 0)

    def cell_at(self, row, column):
        """The call behind one grid position, or None for the target column itself."""
        if not 0 <= row < len(self.rows) or not 1 <= column <= len(self.columns):
            return None
        return self.rows[row]["cells"][self.columns[column - 1]["sample_id"]]

    def _state_words(self, state):
        return self.states.get(state, {})

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row = self.rows[index.row()]
        if index.column() == 0:
            if role == Qt.ItemDataRole.DisplayRole:
                return row["locus"]
            if role == Qt.ItemDataRole.ToolTipRole:
                counts = " · ".join(f"{self._state_words(state).get('label', state)}: {count}"
                                    for state, count in row["states"].items() if count)
                return (f"{row['locus']} · {row['called']} of {row['isolates']} isolate(s) carry an "
                        f"allele for this target.\n{counts}")
            return None
        call = self.cell_at(index.row(), index.column())
        if call is None:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return call["display"]
        if role == Qt.ItemDataRole.ToolTipRole:
            words = self._state_words(call["state"])
            return "\n".join(filter(None, [
                f"{row['locus']} · {words.get('label', call['state'])}",
                call["reason"], words.get("meaning", "")]))
        if role == Qt.ItemDataRole.BackgroundRole:
            tint = STATE_TINTS.get(call["state"], "")
            return QBrush(QColor(tint)) if tint else None
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Vertical:
            return section + 1 if role == Qt.ItemDataRole.DisplayRole else None
        if section == 0:
            if role == Qt.ItemDataRole.DisplayRole:
                return "Target"
            if role == Qt.ItemDataRole.ToolTipRole:
                return "One row per target this scheme defines, in the scheme's own order."
            return None
        if not 0 <= section - 1 < len(self.columns):
            return None
        column = self.columns[section - 1]
        if role == Qt.ItemDataRole.DisplayRole:
            return column["sample_name"]
        if role == Qt.ItemDataRole.ToolTipRole:
            return f"{column['sample_name']}\n{column['summary']}"
        return None


class CgmlstCallsPanel(QWidget):
    """The cgMLST table of calls for one cohort, with the read recheck beside it.

    The host window supplies the project, the background task runner and the status
    line: ``project``, ``launch_task``, ``notify`` and ``scheme_paths``. It is the same
    contract every other workspace page uses, so this panel can be mounted as a tab of
    the comparison page or as a page of its own without changing.
    """

    REQUIRED = ("project", "launch_task", "notify", "scheme_paths")

    def __init__(self, window, parent=None):
        # Checked before the widget exists, so a wiring mistake is named rather than
        # surfacing later as a missing attribute in the middle of a background task.
        missing = [name for name in self.REQUIRED if not hasattr(window, name)]
        if missing:
            raise TypeError("A cgMLST calls panel needs a workspace window providing "
                            + ", ".join(missing) + ".")
        super().__init__(parent if parent is not None else window)
        self.host = window
        self.calls = None       # the built table of calls, or None before one is built
        self.recheck = None     # the last read recheck of one isolate, or None
        self.plan = None        # the plan the last recheck was run from, or None
        # Assigned by a host that also shows a cgMLST tree; it is handed this page's
        # isolates. Left None, the button that offers it is simply not shown.
        self.show_tree_requested = None
        self.references = []
        self._records = {}      # sample id -> the project's own sample record
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self._build_controls(layout)
        self._build_tabs(layout)
        self._build_footer(layout)
        self.refresh_references()

    # --- construction --------------------------------------------------------
    def _build_controls(self, layout):
        strip = FlowLayout()
        self.scope = QComboBox()
        self.scope.addItem("The comparison cohort", "cohort")
        self.scope.addItem("Every isolate in this project", "all")
        self.scope.setToolTip("Which isolates this page reads stored cgMLST profiles for.")
        self.scope.currentIndexChanged.connect(self.refresh_references)
        strip.addWidget(self.scope)
        self.reference = QComboBox()
        self.reference.setMinimumWidth(240)
        self.reference.setToolTip("One table of calls covers one reference. Two references are two "
                                  "different target sets, and their percentages are not the same "
                                  "quantity.")
        strip.addWidget(self.reference)
        self.build_button = button("Build the table of calls", self.build_table, True)
        strip.addWidget(self.build_button)
        self.tree_button = button("Show these isolates in the cgMLST tree", self.show_tree)
        self.tree_button.setToolTip(
            "Open the cgMLST tree for these isolates. It is a separate comparison on its own "
            "scale: a distance over this scheme's targets is never read against a seven-locus "
            "classical MLST distance.")
        self.tree_button.hide()
        strip.addWidget(self.tree_button)
        strip.addWidget(button("Clear", self.clear))
        layout.addLayout(strip)

    def _build_tabs(self, layout):
        self.tabs = QTabWidget()
        self.tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.tabs.addTab(self._build_grid(), "Targets")
        self.tabs.addTab(self._build_isolates(), "Isolates")
        self.tabs.addTab(self._build_recheck(), "Read recheck")
        layout.addWidget(self.tabs, 1)

    def _build_grid(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        filters = FlowLayout()
        self.row_filter = QComboBox()
        for key, title in ROW_FILTERS:
            self.row_filter.addItem(title, key)
        self.row_filter.setToolTip("Which rows are shown. Filtering changes no number: every "
                                   "percentage on this page stays over the scheme's full target count.")
        self.row_filter.currentIndexChanged.connect(self.refresh_rows)
        filters.addWidget(self.row_filter)
        self.locus_search = QLineEdit()
        self.locus_search.setPlaceholderText("Find a target…")
        self.locus_search.setMinimumWidth(140)
        self.locus_search.textChanged.connect(self.refresh_rows)
        filters.addWidget(self.locus_search)
        column.addLayout(filters)
        self.calls_view = QTableView()
        self.model = CallsTableModel(self.calls_view)
        self.calls_view.setModel(self.model)
        self.calls_view.setAlternatingRowColors(True)
        # Deliberately unsortable: the rows are the scheme's own target order, and
        # sorting a column of allele numbers would order targets by a value that
        # means nothing across rows.
        self.calls_view.setSortingEnabled(False)
        self.calls_view.setSelectionBehavior(QTableView.SelectionBehavior.SelectItems)
        self.calls_view.setMinimumHeight(200)
        # Fixed column widths rather than resizeColumnsToContents: measuring every
        # cell of a scheme-length table costs seconds for nothing a reader gains.
        self.calls_view.horizontalHeader().setDefaultSectionSize(150)
        column.addWidget(self.calls_view, 1)
        self.legend = QLabel()
        self.legend.setWordWrap(True)
        self.legend.setTextFormat(Qt.TextFormat.RichText)
        column.addWidget(self.legend)
        return page

    def _build_isolates(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        column.addWidget(label("How much of this scheme each isolate called. A target with no call is "
                               "unknown evidence: it is never counted as a match and never as a "
                               "difference.", "small", True))
        self.isolate_table = make_table(["Isolate", "Targets called", "Reference alleles",
                                         "Local novel", "No call", "Typing", "Assembly SHA-256"])
        column.addWidget(self.isolate_table, 1)
        return page

    def _build_recheck(self):
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        column.addWidget(label(RECHECK_PURPOSE, "small", True))
        strip = FlowLayout()
        self.recheck_isolate = QComboBox()
        self.recheck_isolate.setMinimumWidth(180)
        self.recheck_isolate.setToolTip("The isolate whose uncalled targets are asked about. Its own "
                                        "verified original FASTQ pair is what answers.")
        strip.addWidget(self.recheck_isolate)
        self.scheme_folder = QComboBox()
        self.scheme_folder.setMinimumWidth(220)
        self.scheme_folder.setToolTip("The installed reference these calls were made against. Its "
                                      "fingerprint is checked against the stored profiles before any "
                                      "read is examined.")
        strip.addWidget(self.scheme_folder)
        self.recheck_button = button("Recheck missing targets against the reads…", self.recheck_missing)
        strip.addWidget(self.recheck_button)
        column.addLayout(strip)
        self.recheck_headline = label("", "cardTitle", True)
        column.addWidget(self.recheck_headline)
        self.recheck_controls = label("", "small", True)
        column.addWidget(self.recheck_controls)
        self.recheck_table = make_table(["Target", "What the reads show", "Evidence",
                                         "In this assembly", "Compatible reference alleles",
                                         "What that means"])
        column.addWidget(self.recheck_table, 1)
        self.recheck_limits = label("", "small", True)
        column.addWidget(self.recheck_limits)
        return page

    def _build_footer(self, layout):
        self.denominator = label("", "small", True)
        layout.addWidget(self.denominator)
        self.limitations = label("", "small", True)
        layout.addWidget(self.limitations)
        self.status = label("Choose a reference and build the table of calls.", "small", True)
        layout.addWidget(self.status)

    # --- the cohort ----------------------------------------------------------
    def say(self, message):
        self.status.setText(str(message))
        return message

    def cohort_ids(self):
        """The isolates this page reads, or None for every isolate in the project."""
        if self.scope.currentData() == "all":
            return None
        chosen = getattr(self.host, "cohort_ids", None)
        return None if chosen is None else set(chosen)

    def refresh_references(self):
        """List the references this cohort has stored cgMLST profiles against."""
        self.references = cgmlst_references(self.host.project, self.cohort_ids())
        selected = self.reference.currentData()
        self.reference.blockSignals(True)
        self.reference.clear()
        for entry in self.references:
            self.reference.addItem(
                f"{entry['scheme']} · {entry['targets']:,} targets · "
                f"{entry['isolates']} isolate(s) · {entry['scheme_digest'][:8]}",
                entry["scheme_digest"])
        self.reference.setCurrentIndex(max(0, self.reference.findData(selected)))
        self.reference.blockSignals(False)
        self.build_button.setEnabled(bool(self.references))
        self.tree_button.setVisible(callable(self.show_tree_requested))
        self.tree_button.setEnabled(bool(self.references))
        if not self.references:
            self.say("No isolate in this selection has a stored cgMLST profile yet. Type isolates "
                     "against a cgMLST reference first; an untyped isolate is unknown evidence here, "
                     "not an isolate with no targets.")
        elif self.calls is None:
            entry = self.references[0]
            self.say(f"{entry['isolates']} isolate(s) have a stored profile against "
                     f"{entry['scheme']} ({entry['targets']:,} targets). Build the table of calls.")
        return self.references

    def show_tree(self):
        """Hand this page's isolates to the cgMLST tree, which is its own comparison.

        The tree and this table answer two different questions over the same target
        set. Nothing is carried across but which isolates are being looked at: the
        tree keeps its own reference, threshold and legend.
        """
        if not callable(self.show_tree_requested):
            return self.say("This window has no cgMLST tree to show these isolates in.")
        ids = [row["sample_id"] for row in (self.calls or {}).get("samples") or ()]
        self.show_tree_requested(ids or None)
        if not ids:
            return self.say("The cgMLST tree is showing this cohort on its own scale. Build the "
                            "table of calls first to send it exactly these isolates.")
        return self.say(
            f"The cgMLST tree is showing {len(ids):,} isolate(s). Its distances are counted over "
            "this scheme's targets and share no scale, axis or threshold with a classical "
            "seven-locus MLST distance.")

    def selected_reference(self):
        digest = self.reference.currentData()
        return next((entry for entry in self.references
                     if entry["scheme_digest"] == digest), None) or (self.references[0]
                                                                     if self.references else None)

    # --- the table of calls --------------------------------------------------
    def build_table(self):
        """Build the table off the interface thread; it is cancellable and reports progress.

        The loaded scheme is deliberately not required here: the target set, its order
        and its size all come from the stored profiles themselves, so showing the table
        never waits for thousands of allele files to be read. The read recheck, which
        does need the scheme's own sequences, loads and verifies it then.
        """
        entry = self.selected_reference()
        if entry is None:
            return self.say("There is no stored cgMLST profile in this selection to tabulate.")
        records = [dict(row) for row in entry["records"]]
        self._records = {sample["id"]: sample for sample in self.host.project.samples()}

        def operation(cancelled, progress):
            return build_calls_table(records, cancelled=cancelled, progress=progress)

        self.say(f"Building the table of calls for {len(records)} isolate(s) against "
                 f"{entry['scheme']}…")
        return self.host.launch_task(operation, "cgmlst_calls", self.table_ready)

    def table_ready(self, table):
        """Show a built table of calls, with its own denominator and its own legend."""
        self.calls = table
        self.model.set_table(table)
        self.fill_isolates(table)
        self.fill_legend(table)
        self.fill_recheck_choices(table)
        self.denominator.setText(table["denominator_note"])
        self.limitations.setText("\n".join(table["limitations"]))
        self.refresh_rows()
        self.tabs.setCurrentIndex(0)
        return self.say(
            f"{table['scheme']} · {table['target_count']:,} targets × {table['isolate_count']:,} "
            f"isolate(s) · {table['targets_called_in_every_isolate']:,} targets called in every "
            f"isolate · {table['targets_called_in_no_isolate']:,} called in none.")

    def visible_rows(self):
        """The rows the current filter and search leave on screen, in the scheme's order."""
        if self.calls is None:
            return []
        key = self.row_filter.currentData() or "all"
        query = self.locus_search.text().strip().casefold()
        rows = []
        for row in self.calls["loci"]:
            if key == "incomplete" and row["called_in_every_isolate"]:
                continue
            if key == "absent" and not row["called_in_no_isolate"]:
                continue
            if query and query not in row["locus"].casefold():
                continue
            rows.append(row)
        return rows

    def refresh_rows(self):
        """Re-apply the row filter. The denominator beside the grid never moves with it."""
        if self.calls is None:
            self.model.set_rows([])
            return 0
        rows = self.visible_rows()
        self.model.set_rows(rows)
        self.denominator.setText(
            f"{self.calls['denominator_note']} Showing {len(rows):,} of "
            f"{self.calls['target_count']:,} targets.")
        return len(rows)

    def fill_isolates(self, table):
        self.isolate_table.setSortingEnabled(False)
        self.isolate_table.setRowCount(len(table["samples"]))
        for index, row in enumerate(table["samples"]):
            kind = row["analysis_kind"] if isinstance(row["analysis_kind"], dict) else {}
            states = " · ".join(f"{table['states'][state]['label']}: {row['states'][state]:,}"
                                for state in STATE_ORDER if row["states"][state])
            values = [row["sample_name"], row["summary"], f"{row['called_reference']:,}",
                      f"{row['called_novel']:,}", f"{row['without_call']:,}",
                      TYPING_SCALES.get(str(kind.get("kind") or ""),
                                        TYPING_SCALES["unclassified"])["title"],
                      (row["input_sha256"] or "not recorded")[:12]]
            for column, value in enumerate(values):
                item = cell(value, row["sample_id"])
                item.setToolTip("\n".join(filter(None, [states, str(kind.get("basis") or "")]))
                                or "No call state was recorded for this isolate.")
                self.isolate_table.setItem(index, column, item)
        self.isolate_table.setSortingEnabled(True)

    def fill_legend(self, table):
        """The legend is the module's own wording; this page never invents a label."""
        entries = []
        for state in STATE_ORDER:
            words = table["states"][state]
            tint = STATE_TINTS.get(state) or "#7b8496"
            entries.append(f'<span style="color:{tint}">■</span> {words["label"]} '
                           f'({table["state_totals"][state]:,})')
        self.legend.setText(" &nbsp; ".join(entries))
        self.legend.setToolTip("\n\n".join(f"{table['states'][state]['label']}: "
                                           f"{table['states'][state]['meaning']}"
                                           for state in STATE_ORDER))

    # --- the read recheck ----------------------------------------------------
    def fill_recheck_choices(self, table):
        """Offer the isolates of this table, and the installed references to check against."""
        self.recheck_isolate.clear()
        for row in table["samples"]:
            self.recheck_isolate.addItem(
                f"{row['sample_name']} · {row['without_call']:,} of {row['targets']:,} without a call",
                row["sample_id"])
        self.scheme_folder.clear()
        folders = installed_cgmlst_folders(getattr(self.host, "scheme_paths", ()))
        for folder in folders:
            self.scheme_folder.addItem(folder["title"], folder["path"])
        matching = next((index for index, folder in enumerate(folders)
                         if folder["name"] and folder["name"] == table["scheme"]), None)
        if matching is not None:
            self.scheme_folder.setCurrentIndex(matching)
        self.recheck_button.setEnabled(bool(folders))
        return folders

    def recheck_missing(self):
        """Plan and run the bounded read assay over one isolate's uncalled targets."""
        if self.calls is None:
            return self.say("Build the table of calls first: the recheck asks about the targets it "
                            "shows as uncalled for one isolate.")
        sample_id = self.recheck_isolate.currentData()
        record = self._records.get(sample_id)
        if record is None:
            return self.say("Choose an isolate from this table to recheck.")
        # The isolate's own evidence is checked first: an isolate with no verified
        # read pair cannot be answered by any reference, and that is the refusal the
        # user has to act on.
        try:
            plan = plan_target_recheck(self.calls, sample_id, sample=record, control_count=2)
        except ValueError as error:
            return self.say(str(error))
        path = self.scheme_folder.currentData()
        if not path:
            return self.say("Choose the installed reference these calls were made against: the read "
                            "check needs the scheme's own allele sequences to ask the reads about.")
        if plan["batch_note"]:
            self.host.notify(plan["batch_note"])
        self.plan, table = plan, self.calls

        def operation(cancelled, progress):
            scheme = cached_scheme(path, cancelled)
            return recheck_missing_targets(plan, scheme, cancelled, progress, table=table)

        self.tabs.setCurrentIndex(2)
        self.say(f"Asking the reads of {plan['sample_name']} about {len(plan['targets']):,} uncalled "
                 f"target(s), with {len(plan['controls'])} control target(s) it did call. "
                 + plan["control_basis"])
        return self.host.launch_task(operation, "cgmlst_recheck", self.recheck_ready)

    def recheck_ready(self, result):
        """Show one read verdict per target. Nothing here becomes an allele."""
        self.recheck = result
        self.recheck_headline.setText(result["headline"])
        deferred = result.get("deferred_targets") or []
        self.recheck_controls.setText("\n".join(filter(None, [
            result["controls"]["reason"],
            (f"{len(deferred):,} further uncalled target(s) were not in this batch and have not been "
             "checked: " + ", ".join(deferred[:8]) + ("…" if len(deferred) > 8 else ".")) if deferred else "",
        ])))
        rows = result["targets"]
        self.recheck_table.setSortingEnabled(False)
        self.recheck_table.setRowCount(len(rows))
        for index, row in enumerate(rows):
            candidates = row["compatible_candidates"]
            values = [row["locus"], row["label"], row["evidence"],
                      row["assembly_state_label"] or "not recorded",
                      f"{len(candidates):,} compatible, none assigned",
                      row["explanation"]]
            for column, value in enumerate(values):
                item = cell(value, result["sample_id"])
                item.setToolTip(row["candidate_note"] if column == 4 else row["explanation"])
                self.recheck_table.setItem(index, column, item)
        self.recheck_table.setSortingEnabled(True)
        self.recheck_limits.setText("\n".join(result["limitations"]))
        self.tabs.setCurrentIndex(2)
        return self.say(result["headline"] + " " + RECHECK_PURPOSE)

    # --- clearing ------------------------------------------------------------
    def clear(self):
        """Empty this page so a new, a past or a mixed cohort can be started here.

        Only what is on screen is cleared. Every stored profile, read attachment and
        analysis is left exactly as it was.
        """
        self.calls = self.recheck = self.plan = None
        self._records = {}
        self.model.set_table({})
        self.isolate_table.setRowCount(0)
        self.recheck_table.setRowCount(0)
        self.recheck_isolate.clear()
        self.scheme_folder.clear()
        self.recheck_headline.setText("")
        self.recheck_controls.setText("")
        self.recheck_limits.setText("")
        self.legend.clear()
        self.legend.setToolTip("")
        self.denominator.setText("")
        self.limitations.setText("")
        self.locus_search.clear()
        self.row_filter.setCurrentIndex(0)
        self.tabs.setCurrentIndex(0)
        self.refresh_references()
        return self.say("Cleared. Nothing stored changed: choose a reference and build the table again.")
