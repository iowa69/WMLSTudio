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

Any cell opens the evidence the analysis stored behind it: where in the assembly the
target was found, how well it aligned, what the complete-CDS check said, and the full
sequence identity of a local novel allele. That window displays stored fields and
computes nothing, so what it shows is what was recorded at analysis time; a field the
analysis never wrote is named as not recorded rather than shown as a zero.

The review filter narrows the grid to the targets that carry no allele, for one chosen
isolate or for all of them. It hides rows and changes no number: a hidden target is not
a called target, and every percentage stays over the scheme's full target count.

The recheck button asks the isolate's own reads about the targets its assembly did not
call, and reports what they show: the target is present, the reads are consistent with
it being absent, or there is not enough evidence to say. That is evidence about this
assembly. It assigns no allele and writes nothing back onto a profile, a distance or a
tree.

All of the science is in :mod:`wmlstudio.cgmlst_calls`; this module renders what that
module returns and never re-words its statements.
"""

from __future__ import annotations

import json
from collections.abc import Mapping

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMenu,
    QSizePolicy,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.cgmlst_calls import (
    CALL_STATES,
    STATE_ORDER,
    build_calls_table,
    nomenclature_summary,
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
# is on screen, and the row count says how much of it is being shown. The order of
# these entries is the order of the control, so new ones are appended.
ROW_FILTERS = (
    ("all", "Every target in the scheme"),
    ("incomplete", "Targets without a call in at least one chosen isolate"),
    ("absent", "Targets called in none of the chosen isolates"),
)

# The review filter. 'no_allele' is the review task itself — missing, partial, below
# thresholds, ambiguous, duplicated, mixed and unrecognised are seven different
# answers, so the list also offers each of them on its own rather than only their
# union. Every label is the calling module's own word for the state.
STATE_FILTERS = (("any", "Any call state"),
                 ("no_allele", "Only targets with no allele — the review list"),
                 *((state, f"Only {CALL_STATES[state]['label']}") for state in STATE_ORDER))

# What the detail window is, in one line above the evidence. It displays what was
# stored at analysis time and recomputes nothing, so a field that is not there was
# never written rather than measured as nothing.
DETAIL_PURPOSE = ("Everything this analysis stored about this one isolate and this one target, "
                  "shown as it was recorded. Nothing here is recalculated, and nothing here "
                  "assigns or changes an allele.")

# The per-target evidence a call carries, in reading order: what was decided, what was
# matched, how well it aligned, and where. The note beside each says what the value is,
# never what it proves. Fields this version does not know are still shown, verbatim,
# under the name the analysis wrote them with.
EVIDENCE_FIELDS = (
    ("status", "Recorded call state",
     "The word the caller itself wrote for this target. The answer above is that word in this page's legend."),
    ("allele", "Allele recorded",
     "The reference allele of this target, or the local novel identity. Only the two states that assign an "
     "allele ever carry one."),
    ("sequence_sha256", "Novel sequence SHA-256",
     "The SHA-256 of the observed coding sequence itself. It is a local sequence identity, not a registered "
     "allele number in any nomenclature."),
    ("reason", "Why the caller said so",
     "The caller's own sentence for this target, unchanged."),
    ("candidates", "Exact reference alleles matched",
     "Every reference allele of this target whose complete sequence was found in this assembly."),
    ("cds_candidates", "Predicted coding sequences in play",
     "The predicted CDS this target's reference alleles aligned to. More than one is why no single copy could "
     "be assigned."),
    ("hit_count", "Places matched in the assembly",
     "How many matches were found, before the stored list of places was capped."),
    ("evidence_truncated", "Stored list of places capped",
     "Yes when more matches were found than the stored evidence holds, so the places below are part of them."),
    ("identity_percent", "Identity of the alignment (%)",
     "Percent identity between the predicted coding sequence and the reference allele it aligned to."),
    ("query_coverage", "Reference allele covered (fraction of 1)",
     "How much of the reference allele the alignment spans."),
    ("subject_coverage", "Predicted CDS covered (fraction of 1)",
     "How much of the predicted coding sequence the alignment spans."),
    ("nearest_reference_allele", "Nearest reference allele",
     "The reference allele this sequence aligned to best. It is not the allele of this target: a novel sequence "
     "is given no allele number here."),
    ("cds_qc", "Complete-CDS check",
     "The full coding-sequence validation of the predicted gene: a complete start, a terminal stop, no internal "
     "stop, a length in whole codons, and no contig edge."),
    ("best_alignment", "Best alignment that did not qualify",
     "The strongest homolog found for this target that failed the identity, coverage or complete-CDS criteria, "
     "which is why no allele was assigned from it."),
)
# Rendered elsewhere in the window, so they are not repeated as fields.
_DETAIL_SKIPPED = frozenset({"locus", "hits", "display", "state", "copies_lower_bound"})
# The run's own settings, shown because a partial or below-threshold call can only be
# read against the criteria that were in force when it was made.
PARAMETER_FIELDS = (
    ("min_identity", "Run setting · minimum identity (fraction of 1)"),
    ("min_coverage", "Run setting · minimum coverage (fraction of 1)"),
    ("genetic_code", "Run setting · genetic code"),
    ("method", "Run setting · calling method"),
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


def value_text(value):
    """One stored value as text, without rounding it into something it is not.

    Numbers are printed as they were stored rather than rescaled, because a
    coverage fraction and a percentage are read differently and the field's own
    name carries the unit. A long list says how much of it is on screen.
    """
    if value is None:
        return "not recorded"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, (list, tuple)):
        if not value:
            return "none recorded"
        shown = [value_text(item) if isinstance(item, (dict, list, tuple)) else str(item)
                 for item in value[:40]]
        return ", ".join(shown) + (f" … and {len(value) - 40:,} more" if len(value) > 40 else "")
    if isinstance(value, Mapping):
        return json.dumps(value, sort_keys=True, default=str)
    return str(value)


def evidence_rows(call):
    """Every field of one stored call, as (rows, fields the analysis did not write).

    Known fields come first with a plain sentence saying what the value is; anything
    this version does not recognise is still shown under the name it was stored with,
    so a field the caller writes can never be silently dropped from a review.
    """
    call = dict(call) if isinstance(call, Mapping) else {}
    rows, absent, seen = [], [], set()
    for key, title, note in EVIDENCE_FIELDS:
        seen.add(key)
        value = call.get(key)
        if key not in call or value is None:
            absent.append(title)
            continue
        if isinstance(value, Mapping):
            # A nested block is unfolded rather than printed as one blob, because its
            # parts are what a reviewer actually reads: which QC reason failed, which
            # allele the rejected alignment was against.
            for name, item in value.items():
                rows.append({"key": f"{key}.{name}", "name": f"{title} · {name}",
                             "value": value_text(item), "note": note})
            continue
        rows.append({"key": key, "name": title, "value": value_text(value), "note": note})
    for key in call:
        if key in seen or key in _DETAIL_SKIPPED:
            continue
        rows.append({"key": key, "name": key, "value": value_text(call[key]),
                     "note": "Recorded by the analysis. This version has no plain-language note for this "
                             "field, so it is shown exactly as it was stored."})
    return rows, absent


class TargetCallDialog(QDialog):
    """The stored evidence behind one isolate × one target, and nothing beyond it.

    Every number here was written by the analysis that made the call. This window
    reads them out; it computes no identity, no coverage and no verdict of its own,
    and it never turns read or alignment evidence into an allele.
    """

    def __init__(self, detail, parent=None):
        super().__init__(parent)
        self.detail = detail
        self.setWindowTitle(f"Stored evidence · {detail['locus']} · {detail['sample_name']}")
        self.resize(940, 700)
        layout = QVBoxLayout(self)
        layout.addWidget(label(f"{detail['locus']} · {detail['sample_name']}", "cardTitle", True))
        self.answer = label(detail["answer"], "title", True)
        layout.addWidget(self.answer)
        # The state's own meaning is printed in full and unedited. For a target with no
        # allele that sentence is the whole point of the window: it says what was not
        # found, not that the target is absent from the isolate.
        self.meaning = label(detail["meaning"], "small", True)
        layout.addWidget(self.meaning)
        if detail["reason"]:
            layout.addWidget(label(detail["reason"], "small", True))
        if detail["next_step"]:
            layout.addWidget(label(detail["next_step"], "small", True))
        layout.addWidget(label(DETAIL_PURPOSE, "muted", True))
        self.fields = make_table(["What the analysis recorded", "Value", "What that value is"])
        self.fields.setSortingEnabled(False)
        self.fields.setRowCount(len(detail["fields"]))
        for index, row in enumerate(detail["fields"]):
            for column, value in enumerate((row["name"], row["value"], row["note"])):
                item = cell(value, detail["sample_id"])
                item.setToolTip(row["note"])
                self.fields.setItem(index, column, item)
        layout.addWidget(self.fields, 2)
        self.absent_note = label(detail["absent_note"], "small", True)
        layout.addWidget(self.absent_note)
        self.places_note = label(detail["places_note"], "small", True)
        layout.addWidget(self.places_note)
        self.places = make_table(["Contig", "First base", "Last base", "Strand",
                                  "Allele matched here", "Sequence shared with another target"])
        self.places.setSortingEnabled(False)
        self.places.setRowCount(len(detail["places"]))
        for index, place in enumerate(detail["places"]):
            values = [place.get("contig"), place.get("start"), place.get("end"), place.get("strand"),
                      place.get("allele"), place.get("shared_reference")]
            for column, value in enumerate(values):
                item = cell(value_text(value), detail["sample_id"])
                item.setToolTip("One place in this assembly where the stored evidence puts this target.")
                self.places.setItem(index, column, item)
        layout.addWidget(self.places, 1)
        provenance = label(detail["provenance"], "small", True)
        provenance.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(provenance)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.addButton(button("Copy this evidence as JSON", self.copy_evidence),
                          QDialogButtonBox.ButtonRole.ActionRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def copy_evidence(self):
        """Put the stored call on the clipboard exactly as it is held, for a report."""
        payload = json.dumps(self.detail["stored"], indent=2, sort_keys=True, default=str)
        clipboard = QApplication.clipboard()
        if clipboard is not None:
            clipboard.setText(payload)
        return payload


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
        self.detail = None      # the open per-target evidence window, or None
        self.references = []
        self._records = {}      # sample id -> the project's own sample record
        # sample id -> {locus: the call evidence the analysis stored}. The table of
        # calls carries a cell's state and reason; the alignment, the QC and the
        # places behind it stay in the stored result, and the detail window reads
        # them from here rather than having the table carry thousands of copies.
        self._evidence = {}
        self._profiles = {}     # sample id -> the stored profile the calls came from
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
        self.state_filter = QComboBox()
        self.state_filter.setMinimumWidth(220)
        for key, title in STATE_FILTERS:
            self.state_filter.addItem(title, key)
        self.state_filter.setToolTip(
            "The review filter. Reviewing a scheme means looking at the targets that carry no allele, "
            "one call state at a time. A target hidden by this filter is not a called target.")
        self.state_filter.currentIndexChanged.connect(self.refresh_rows)
        filters.addWidget(self.state_filter)
        self.isolate_filter = QComboBox()
        self.isolate_filter.setMinimumWidth(180)
        self.isolate_filter.addItem("In any isolate of this table", None)
        self.isolate_filter.setToolTip("Which isolates the filters above are read over. Every isolate "
                                       "keeps its own column whatever is chosen here.")
        self.isolate_filter.currentIndexChanged.connect(self.refresh_rows)
        filters.addWidget(self.isolate_filter)
        self.locus_search = QLineEdit()
        self.locus_search.setPlaceholderText("Find a target…")
        self.locus_search.setMinimumWidth(140)
        self.locus_search.textChanged.connect(self.refresh_rows)
        filters.addWidget(self.locus_search)
        self.detail_button = button("Show the full evidence for this call…", self.open_call_detail)
        self.detail_button.setToolTip("Open everything the analysis stored about the selected cell. "
                                      "Double-clicking a cell, or right-clicking it, does the same.")
        filters.addWidget(self.detail_button)
        column.addLayout(filters)
        self.review_note = label("", "small", True)
        column.addWidget(self.review_note)
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
        # A scheme is thousands of rows, so the grid is built for scanning a column:
        # short fixed rows, narrow fixed columns, no grid lines and no wrapping. Fixed
        # sizes rather than resizeColumnsToContents, because measuring every cell of a
        # scheme-length table costs seconds for nothing a reader gains.
        self.calls_view.horizontalHeader().setDefaultSectionSize(132)
        self.calls_view.horizontalHeader().setMinimumSectionSize(64)
        self.calls_view.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        self.calls_view.verticalHeader().setDefaultSectionSize(22)
        self.calls_view.setShowGrid(False)
        self.calls_view.setWordWrap(False)
        self.calls_view.setTextElideMode(Qt.TextElideMode.ElideRight)
        self.calls_view.doubleClicked.connect(self.open_call_detail)
        self.calls_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.calls_view.customContextMenuRequested.connect(self.show_grid_menu)
        self.calls_view.selectionModel().currentChanged.connect(self.describe_selection)
        column.addWidget(self.calls_view, 1)
        self.selected_note = label("", "small", True)
        column.addWidget(self.selected_note)
        self.legend = QLabel()
        self.legend.setWordWrap(True)
        self.legend.setTextFormat(Qt.TextFormat.RichText)
        column.addWidget(self.legend)
        self.describe_selection()
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
                                         "Local novel", "No call", "Typing",
                                         "Published nomenclature", "Assembly SHA-256"])
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
        self._profiles = {str(row.get("sample_id")): row for row in records}
        self._evidence = {
            str(row.get("sample_id")): {str(call.get("locus")): call
                                        for call in (row.get("calls") or ())
                                        if isinstance(call, Mapping) and call.get("locus")}
            for row in records}

        def operation(cancelled, progress):
            return build_calls_table(records, cancelled=cancelled, progress=progress)

        self.say(f"Building the table of calls for {len(records)} isolate(s) against "
                 f"{entry['scheme']}…")
        return self.host.launch_task(operation, "cgmlst_calls", self.table_ready)

    def table_ready(self, table):
        """Show a built table of calls, with its own denominator and its own legend."""
        self.calls = table
        self.model.set_table(table)
        self.fill_isolate_filter(table)
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

    def chosen_isolates(self):
        """The isolates the filters are read over: one chosen column, or all of them."""
        if self.calls is None:
            return []
        chosen = self.isolate_filter.currentData()
        ids = [row["sample_id"] for row in self.calls["samples"]]
        return [chosen] if chosen in ids else ids

    def visible_rows(self):
        """The rows the current filters and search leave on screen, in the scheme's order.

        A row survives when at least one of the chosen isolates matches, because the
        review question is asked of isolates one at a time and a target uncalled in one
        of them is a target to review even when another called it.
        """
        if self.calls is None:
            return []
        key = self.row_filter.currentData() or "all"
        state = self.state_filter.currentData() or "any"
        ids = self.chosen_isolates()
        query = self.locus_search.text().strip().casefold()
        rows = []
        for row in self.calls["loci"]:
            cells = [row["cells"][sample_id] for sample_id in ids if sample_id in row["cells"]]
            called = [cell_ for cell_ in cells if CALL_STATES[cell_["state"]]["has_allele"]]
            if key == "incomplete" and len(called) == len(cells):
                continue
            if key == "absent" and called:
                continue
            if state == "no_allele" and len(called) == len(cells):
                continue
            if state not in {"any", "no_allele"} and not any(cell_["state"] == state for cell_ in cells):
                continue
            if query and query not in row["locus"].casefold():
                continue
            rows.append(row)
        return rows

    def filter_description(self):
        """What the filters are currently asking, in the same words the controls use."""
        chosen = self.isolate_filter.currentData()
        name = next((row["sample_name"] for row in (self.calls or {}).get("samples") or ()
                     if row["sample_id"] == chosen), None)
        where = f"in {name}" if name else f"in any of {len(self.chosen_isolates()):,} isolate(s)"
        parts = [self.row_filter.currentText()]
        if (self.state_filter.currentData() or "any") != "any":
            parts.append(self.state_filter.currentText().removeprefix("Only ").strip())
        query = self.locus_search.text().strip()
        if query:
            parts.append(f"name contains {query!r}")
        return " · ".join(parts) + " · " + where

    def refresh_rows(self):
        """Re-apply the filters. The denominator beside the grid never moves with them."""
        if self.calls is None:
            self.model.set_rows([])
            self.review_note.setText("")
            return 0
        rows = self.visible_rows()
        self.model.set_rows(rows)
        self.denominator.setText(
            f"{self.calls['denominator_note']} Showing {len(rows):,} of "
            f"{self.calls['target_count']:,} targets.")
        self.review_note.setText(
            f"{len(rows):,} of {self.calls['target_count']:,} targets shown · {self.filter_description()}. "
            "Filtering hides rows and changes nothing else: a target that is not on screen is not a "
            "called target, and every count on this page stays over the scheme's full target set.")
        self.describe_selection()
        return len(rows)

    def fill_isolate_filter(self, table):
        """Offer this table's isolates to the review filter, keeping the current choice.

        Every isolate keeps its column in the grid whatever is chosen here: this scopes
        the question being asked of the rows, not which isolates the page covers.
        """
        chosen = self.isolate_filter.currentData()
        self.isolate_filter.blockSignals(True)
        self.isolate_filter.clear()
        self.isolate_filter.addItem("In any isolate of this table", None)
        for row in table["samples"]:
            self.isolate_filter.addItem(f"In {row['sample_name']} only", row["sample_id"])
        self.isolate_filter.setCurrentIndex(max(0, self.isolate_filter.findData(chosen)))
        self.isolate_filter.blockSignals(False)
        return self.isolate_filter.count()

    def fill_isolates(self, table):
        self.isolate_table.setSortingEnabled(False)
        self.isolate_table.setRowCount(len(table["samples"]))
        for index, row in enumerate(table["samples"]):
            kind = row["analysis_kind"] if isinstance(row["analysis_kind"], dict) else {}
            states = " · ".join(f"{table['states'][state]['label']}: {row['states'][state]:,}"
                                for state in STATE_ORDER if row["states"][state])
            naming = nomenclature_summary(table, row["sample_id"])
            values = [row["sample_name"], row["summary"], f"{row['called_reference']:,}",
                      f"{row['called_novel']:,}", f"{row['without_call']:,}",
                      TYPING_SCALES.get(str(kind.get("kind") or ""),
                                        TYPING_SCALES["unclassified"])["title"],
                      naming["headline"],
                      (row["input_sha256"] or "not recorded")[:12]]
            for column, value in enumerate(values):
                item = cell(value, row["sample_id"])
                # Whether this profile can be set beside a published one is the
                # question somebody comparing with a paper is actually asking, so
                # the rule and the denominator travel with the cell that answers it.
                item.setToolTip("\n\n".join(filter(None, [
                    states, str(kind.get("basis") or ""), naming["comparability"],
                    *naming["notes"]]))
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

    # --- one call, in full ---------------------------------------------------
    def cell_position(self, index):
        """The (target, isolate) a grid position names, or None outside the cells.

        The first column is the target's own name and belongs to every isolate in the
        row, so it names no single call and opens nothing.
        """
        if self.calls is None or not isinstance(index, QModelIndex) or not index.isValid():
            return None
        row, column = index.row(), index.column()
        if not 0 <= row < len(self.model.rows) or not 1 <= column <= len(self.model.columns):
            return None
        return self.model.rows[row]["locus"], self.model.columns[column - 1]["sample_id"]

    def call_detail(self, locus, sample_id):
        """Gather everything stored about one isolate × one target, for display only.

        The state, its label and its meaning come from the table of calls; the
        alignment, the complete-CDS check and the places come from the stored result
        exactly as the analysis wrote them. Nothing is recomputed here, and a field the
        analysis never wrote is named as not recorded rather than shown as a zero.
        """
        if self.calls is None:
            return None
        row = next((entry for entry in self.calls["loci"] if entry["locus"] == locus), None)
        sample = next((entry for entry in self.calls["samples"]
                       if entry["sample_id"] == sample_id), None)
        if row is None or sample is None:
            return None
        call = row["cells"][sample_id]
        words = self.calls["states"].get(call["state"], {})
        stored = (self._evidence.get(sample_id) or {}).get(locus) or {}
        profile = self._profiles.get(sample_id) or {}
        fields, absent = evidence_rows(stored)
        for key, title in PARAMETER_FIELDS:
            value = (profile.get("parameters") or {}).get(key)
            if value is not None:
                fields.append({"key": key, "name": title, "value": value_text(value),
                               "note": "The criterion in force when this isolate was called. It is what "
                                       "this call was judged against, not a threshold this page applies."})
        novel = next((entry for entry in (profile.get("novel_sequences") or ())
                      if isinstance(entry, Mapping) and entry.get("locus") == locus), None)
        if novel and novel.get("sequence"):
            fields.append({"key": "novel_sequence", "name": "Novel sequence stored with this profile",
                           "value": f"{len(str(novel['sequence'])):,} bases · SHA-256 "
                                    f"{novel.get('sequence_sha256') or 'not recorded'}",
                           "note": "The complete coding sequence this local novel identity names. Two "
                                   "isolates share a local novel allele only when this sequence is identical."})
        places = [hit for hit in (stored.get("hits") or ()) if isinstance(hit, Mapping)]
        label_text = words.get("label", call["state"])
        answer = f"{label_text}: {call['display']}" if words.get("has_allele") else label_text
        if places:
            places_note = (f"{len(places):,} place(s) in this assembly carry the stored evidence for this "
                           "target.")
            if stored.get("evidence_truncated"):
                places_note += (f" The stored list is capped, so these are part of the "
                                f"{value_text(stored.get('hit_count'))} match(es) that were found.")
        elif words.get("has_allele"):
            places_note = "This analysis recorded no position for this target beside its allele."
        else:
            # The empty case is exactly the one a reader can misread, so the sentence
            # says what the emptiness is: nothing was recorded, which is not a finding
            # that the target is absent from the isolate.
            places_note = ("No place in this assembly was recorded for this target. That is the absence of "
                           "a recorded match, not a finding that the target is absent from the isolate.")
        next_step = ("" if words.get("has_allele") else
                     "The Read recheck tab asks this isolate's own reads about targets like this one. It "
                     "reports what the reads show and assigns no allele.")
        return {
            "locus": locus, "sample_id": sample_id, "sample_name": sample["sample_name"],
            "state": call["state"], "state_label": label_text, "answer": answer,
            "meaning": words.get("meaning", ""), "reason": call["reason"],
            "display": call["display"], "allele": call["allele"],
            "has_allele": bool(words.get("has_allele")), "next_step": next_step,
            "fields": fields, "places": places, "places_note": places_note,
            "absent": absent,
            "absent_note": ("This analysis recorded no value for: " + ", ".join(absent)
                            + ". A field with no value was not measured for this target; it is not a zero "
                              "and not a negative result.") if absent else
                           "Every field this version knows of was recorded for this target.",
            "provenance": (f"Reference: {self.calls['scheme']} · fingerprint "
                           f"{self.calls['scheme_digest']}\nAssembly SHA-256: "
                           f"{profile.get('input_sha256') or 'not recorded'}"),
            "stored": stored,
        }

    def open_call_detail(self, index=None):
        """Open the stored evidence behind the selected cell in its own window."""
        if self.calls is None:
            return self.say("Build the table of calls first: the evidence window shows what one isolate's "
                            "analysis recorded for one target.")
        if not isinstance(index, QModelIndex) or not index.isValid():
            index = self.calls_view.currentIndex()
        position = self.cell_position(index)
        if position is None:
            return self.say("Choose a cell under an isolate's column. The target column names the row; the "
                            "evidence belongs to one isolate's call for it.")
        detail = self.call_detail(*position)
        if detail is None:
            return self.say("That call is no longer in this table. Build it again to open its evidence.")
        self.detail = TargetCallDialog(detail, self)
        self.detail.show()
        self.say(f"{detail['locus']} · {detail['sample_name']} · {detail['answer']}. "
                 + (detail["meaning"] or ""))
        return self.detail

    def describe_selection(self, current=None, previous=None):
        """Say in full what the selected cell is, so no state is read off colour alone."""
        del previous
        position = self.cell_position(current if isinstance(current, QModelIndex)
                                      else self.calls_view.currentIndex())
        if position is None:
            self.selected_note.setText("Select a cell to read its call state here, or double-click it for "
                                       "everything the analysis stored about it.")
            return ""
        locus, sample_id = position
        row = next((entry for entry in self.calls["loci"] if entry["locus"] == locus), None)
        if row is None:
            return ""
        call = row["cells"][sample_id]
        words = self.calls["states"].get(call["state"], {})
        name = next((entry["sample_name"] for entry in self.calls["samples"]
                     if entry["sample_id"] == sample_id), sample_id)
        text = " · ".join(filter(None, [f"{locus} · {name}", words.get("label", call["state"]),
                                        call["reason"], words.get("meaning", "")]))
        self.selected_note.setText(text)
        return text

    def apply_state_filter(self, state, sample_id=None):
        """Point the review filter at one call state, optionally in one isolate."""
        if self.calls is None:
            return self.say("Build the table of calls first: there are no rows to filter yet.")
        self.state_filter.blockSignals(True)
        self.state_filter.setCurrentIndex(max(0, self.state_filter.findData(state)))
        self.state_filter.blockSignals(False)
        self.isolate_filter.blockSignals(True)
        self.isolate_filter.setCurrentIndex(max(0, self.isolate_filter.findData(sample_id)))
        self.isolate_filter.blockSignals(False)
        shown = self.refresh_rows()
        return self.say(f"Showing {shown:,} of {self.calls['target_count']:,} targets · "
                        f"{self.filter_description()}. Filtering hides rows and changes no number.")

    def grid_menu(self, index):
        """The right-click menu for one grid position, built so it can be tested."""
        menu = QMenu(self.calls_view)
        position = self.cell_position(index)
        if position is not None:
            locus, sample_id = position
            row = next(entry for entry in self.calls["loci"] if entry["locus"] == locus)
            call = row["cells"][sample_id]
            words = self.calls["states"].get(call["state"], {})
            name = next((entry["sample_name"] for entry in self.calls["samples"]
                         if entry["sample_id"] == sample_id), sample_id)
            evidence = menu.addAction(f"Show the full stored evidence for {locus} in {name}")
            evidence.triggered.connect(lambda checked=False, chosen=index: self.open_call_detail(chosen))
            same = menu.addAction(f"Show only {words.get('label', call['state'])} in {name}")
            same.triggered.connect(lambda checked=False, state=call["state"], sid=sample_id:
                                   self.apply_state_filter(state, sid))
            review = menu.addAction(f"Show only targets with no allele in {name}")
            review.triggered.connect(lambda checked=False, sid=sample_id:
                                     self.apply_state_filter("no_allele", sid))
        reset = menu.addAction("Show every target again")
        reset.triggered.connect(lambda checked=False: self.apply_state_filter("any", None))
        return menu

    def show_grid_menu(self, point):
        if self.calls is None:
            return None
        menu = self.grid_menu(self.calls_view.indexAt(point))
        menu.exec(self.calls_view.viewport().mapToGlobal(point))
        return menu

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
        self._evidence = {}
        self._profiles = {}
        if self.detail is not None:
            # An evidence window outlives nothing: the table it read from is gone, so
            # leaving it open would invite reading it as the state of this page.
            self.detail.close()
            self.detail = None
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
        self.state_filter.setCurrentIndex(0)
        self.isolate_filter.clear()
        self.isolate_filter.addItem("In any isolate of this table", None)
        self.review_note.setText("")
        self.selected_note.setText("")
        self.tabs.setCurrentIndex(0)
        self.refresh_references()
        return self.say("Cleared. Nothing stored changed: choose a reference and build the table again.")
