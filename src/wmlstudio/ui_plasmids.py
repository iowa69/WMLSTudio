"""The plasmid page: replicon markers, the contigs they sit on, and what none of it proves.

This page exists because "there is no plasmid analysis" was true of the window
even though the evidence was there: the replicon tables sat in the fourth sub-tab
of a strip that ran off the edge, and the per-isolate contig evidence existed only
as HTML inside a detail pane, where it could not be sorted, widened or copied. It
is all tabulated here, in one place, in the tables every other page uses.

Everything on this page is contig-level evidence from markers already in hand. It
is not MOB-suite: no relaxase or mate-pair-formation type is assigned, no origin
of transfer is searched for, no contig is binned into a predicted plasmid and no
plasmid is counted. Those sentences are printed beside the numbers rather than in
a help page, because a reader who misses them reads a replicon marker as a
plasmid, and a shared contig as proof that a carbapenemase is mobile.

The distinction this page guards hardest is between a question that was answered
and a question that was never asked. A reference store holding no plasmid set
reports no replicon for any isolate on earth, and an empty replicon column is then
the same picture as a cohort of clean isolates. So "no plasmid database is
installed", "installed but not read by this isolate's run", "read and nothing met
the thresholds" and "no plasmid assay to read" are four different cells here and
never collapse into one blank.

Nothing here is a distance. No table shares a scale, an axis, a column or a
threshold with a seven-locus MLST distance, a cgMLST target distance or a SNP
distance, and no count on this page can be clustered on.

All of the science is in :mod:`wmlstudio.plasmid_evidence`. This module arranges
what it returns and re-words none of it.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import plasmid_evidence
from wmlstudio.ui_common import FlowLayout, cell, make_table
from wmlstudio.widgets import button, card, label

TITLE = "Plasmid replicon evidence"

QUESTION = "Which plasmid replicon markers are in these isolates, and what sits on the same contig?"

PURPOSE = ("Replicon markers reported by the AMR / plasmid screen, the assembled contigs they were "
           "found on, and the resistance or virulence determinants that share those contigs. "
           "Contig-level evidence only: nothing here is reconstructed, predicted or counted as a "
           "plasmid.")

#: Printed under the header, where every reader passes it. These are the three
#: sentences that must never be separated from the numbers on this page.
BOUNDARY = ("A replicon marker is not a plasmid. Co-location on one assembled contig is not proof "
            "that a determinant is plasmid-borne. A short-read assembly cannot resolve plasmid "
            "structure, so neither a plasmid's identity nor its transferability can be read here.")

#: The tab each table gets, in the order a reader works through them: what was
#: found, who shares it, what sits beside it, what recurs, the raw contigs, and
#: finally who could not be asked at all.
TAB_TITLES = {
    "replicons_by_isolate": "Replicons by isolate",
    "replicon_matrix": "Replicon × isolate",
    "determinant_colocation": "Determinants on contigs",
    "cohort_cooccurrence": "Recurring co-locations",
    "contigs": "Contig evidence",
    "isolates_not_assayed": "Not assayed",
}

#: Beyond this many isolates the payload is built on the background runner. Below
#: it the work is a few hundred stored dictionaries and a progress dialog for that
#: is noise. No number changes either way: it is the same function on both paths.
BACKGROUND_COHORT = 40

#: Tables whose rows must never be re-ordered by a column. A grid whose columns are
#: replicons sorts by a word, not by evidence, and a co-occurrence table ranked by
#: count reads as a cluster list when nothing here is a distance.
UNSORTABLE = ("replicon_matrix", "cohort_cooccurrence")

NOTHING_YET = ("No isolate is in this cohort yet. Choose isolates to see what their plasmid "
               "screens reported, and what they were never asked.")


class PlasmidPanel(QWidget):
    """One cohort's plasmid evidence: replicons, contigs, co-locations and refusals.

    The host window supplies the project and the background task runner:
    ``project``, ``launch_task`` and ``notify``. ``cohort_ids`` and
    ``active_amr_database`` are used when the host has them, so this panel can be
    mounted as a page or built on its own without changing.
    """

    REQUIRED = ("project", "launch_task", "notify")

    def __init__(self, window, parent=None):
        missing = [name for name in self.REQUIRED if not hasattr(window, name)]
        if missing:
            raise TypeError("A plasmid panel needs a workspace window providing "
                            + ", ".join(missing) + ".")
        super().__init__(parent if parent is not None else window)
        self.host = window
        # The payload this page is showing. It is the same dictionary the report
        # consumes, so a table on screen and the same table in a document are one
        # set of numbers with one set of words beside them.
        self.payload = None
        self.tables = {}
        self.empties = {}
        # The host's status bar may not exist while the window is still building
        # its pages, so the first fill speaks only on this page's own status line.
        self._ready = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self._build_header(layout)
        self._build_controls(layout)
        self._build_tabs(layout)
        self._build_footer(layout)
        self.refresh()
        self._ready = True

    # --- construction --------------------------------------------------------
    def _build_header(self, layout):
        frame, inner = card()
        inner.setContentsMargins(16, 12, 16, 12)
        inner.setSpacing(6)
        inner.addWidget(label(TITLE, "badge"))
        inner.addWidget(label(QUESTION, "cardTitle", wrap=True))
        inner.addWidget(label(PURPOSE, "muted", wrap=True))
        inner.addWidget(label(BOUNDARY, "small", wrap=True))
        layout.addWidget(frame)

    def _build_controls(self, layout):
        strip = FlowLayout()
        self.scope = QComboBox()
        # Every isolate first: a page that opened empty because nothing happened
        # to be selected elsewhere reads as a page where nothing was found.
        self.scope.addItem("Every isolate in this project", "all")
        self.scope.addItem("The isolates selected elsewhere", "cohort")
        self.scope.setToolTip("Which isolates these tables cover. Isolates with no plasmid result "
                              "to read are listed with the reason, never counted as carrying "
                              "nothing.")
        self.scope.currentIndexChanged.connect(self.refresh)
        strip.addWidget(self.scope)
        refresh = button("Re-read the stored evidence", self.refresh, True)
        refresh.setToolTip("Reads the plasmid evidence already stored on these isolates. No "
                           "assembly is read again and no screen is run.")
        strip.addWidget(refresh)
        layout.addLayout(strip)
        # The database state is the first thing to read on this page and the last
        # thing that may be hidden in a tab: an empty table means nothing until a
        # reader knows whether anything could have filled it.
        self.database_state = label("", "badge")
        self.database_state.setWordWrap(True)
        layout.addWidget(self.database_state)
        self.headline = label(NOTHING_YET, "small", True)
        layout.addWidget(self.headline)

    def _build_tabs(self, layout):
        self.tabs = QTabWidget()
        self.tabs.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        for key, title in TAB_TITLES.items():
            self.tabs.addTab(self._build_table_tab(key), title)
        layout.addWidget(self.tabs, 1)

    def _build_table_tab(self, key):
        """One table, with the limit that belongs to it printed above the numbers."""
        page = QWidget()
        column = QVBoxLayout(page)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        note = label(plasmid_evidence.TABLE_NOTES[key], "small", True)
        column.addWidget(note)
        table = make_table(["Isolate"])
        # Qt's untouched sort indicator is column 0 descending, so re-enabling
        # sorting after a fill would silently reverse the isolate order. The
        # indicator is set once here; a reader who clicks a heading still gets
        # their own order back on the next fill.
        table.sortByColumn(0, Qt.SortOrder.AscendingOrder)
        if key in UNSORTABLE:
            table.setSortingEnabled(False)
        self.tables[key] = table
        column.addWidget(table, 1)
        # Said where the missing rows would have been. An empty table with no
        # sentence under it is the false negative this whole page guards against.
        self.empties[key] = label("", "small", True)
        column.addWidget(self.empties[key])
        return page

    def _build_footer(self, layout):
        self.limitations = label("", "small", True)
        layout.addWidget(self.limitations)
        self.gap = label("", "small", True)
        self.gap.setToolTip("What a MOB-suite run would have added that this screen does not.")
        layout.addWidget(self.gap)
        self.status = label("", "small", True)
        layout.addWidget(self.status)

    # --- what this page is looking at ----------------------------------------
    def say(self, message):
        self.status.setText(str(message))
        notify = getattr(self.host, "notify", None)
        if self._ready and callable(notify):
            notify(str(message))
        return message

    def cohort_ids(self):
        """The isolates this page covers, or None for every isolate in the project."""
        if self.scope.currentData() == "all":
            return None
        chosen = getattr(self.host, "cohort_ids", None)
        if callable(chosen):
            chosen = chosen()
        return None if chosen is None else {str(value) for value in chosen}

    def cohort_name(self):
        return self.scope.currentText()

    def installed_sets(self):
        """What the reference store holds, or ``None`` when it cannot be asked.

        ``None`` is not an empty store. Only the store can say that a reference
        set is absent, so a page that never reached one reports the question as
        unanswered rather than reporting the set as missing.
        """
        resolve = getattr(self.host, "active_amr_database", None)
        try:
            root = resolve() if callable(resolve) else None
        except (OSError, ValueError):
            return None
        if not root:
            return None
        from wmlstudio.hydra_runtime import database_status
        try:
            status = database_status(root)
        except (OSError, ValueError):
            return None
        # A manifest that could not be read says nothing about what is installed,
        # and must never be reported as a store holding no plasmid set.
        return None if status.get("error") else sorted(status.get("installed") or ())

    def entries(self):
        """One plain record per isolate: its identity and its stored plasmid block.

        An isolate whose characterization is missing, stale or unverified carries
        that reason into the table rather than an empty replicon list, because an
        empty list and an unasked question are the same picture otherwise.
        """
        from wmlstudio.characterization import current_characterization
        chosen = self.cohort_ids()
        rows = []
        for sample in self.host.project.samples():
            if chosen is not None and str(sample["id"]) not in chosen:
                continue
            state = current_characterization(sample)
            evidence = state.get("evidence") or {}
            rows.append({"sample_id": str(sample["id"]),
                         "sample_name": sample.get("name") or str(sample["id"]),
                         "plasmid_hypotheses": evidence.get("plasmid_hypotheses")
                         or {"status": "not_run", "reason": state["reason"]}})
        return rows

    # --- building the payload ------------------------------------------------
    def refresh(self):
        """Rebuild every table from the evidence already stored on these isolates."""
        entries, installed = self.entries(), self.installed_sets()
        cohort = self.cohort_name()
        if len(entries) <= BACKGROUND_COHORT:
            return self.show_payload(plasmid_evidence.plasmid_payload(
                entries, installed=installed, cohort=cohort))

        def build(cancelled, report):
            return plasmid_evidence.plasmid_payload(entries, installed=installed, cohort=cohort,
                                                    cancelled=cancelled)

        if self.host.launch_task(build, "plasmids", self.show_payload) is False:
            return self.say("A background task is already running; finish or cancel it first. "
                            "Nothing on this page changed.")
        return self.say(f"Reading the stored plasmid evidence of {len(entries)} isolate(s). "
                        "No assembly is read again and no screen is run.")

    def show_payload(self, payload):
        """Put one payload on the page: the state of the question, then the tables."""
        payload = payload if isinstance(payload, dict) else {}
        self.payload = payload or None
        database = payload.get("database") or {}
        self.database_state.setText(_database_words(database))
        self.database_state.setToolTip(database.get("reason", ""))
        self.headline.setText(payload.get("headline") or NOTHING_YET)
        tables = payload.get("tables") or {}
        for key in TAB_TITLES:
            self._fill(key, tables.get(key))
        self.limitations.setText("\n".join(payload.get("limitations") or ()))
        self.gap.setText("\n".join(payload.get("mob_suite_gap") or ()))
        return self.say(payload.get("headline") or NOTHING_YET)

    def _fill(self, key, spec):
        """One table's columns, rows and empty-state sentence, in that order."""
        table = self.tables[key]
        spec = spec if isinstance(spec, dict) else {}
        columns = list(spec.get("columns") or ["Isolate"])
        rows = list(spec.get("rows") or [])
        sorting = table.isSortingEnabled()
        table.setSortingEnabled(False)
        table.clearContents()
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels([str(name) for name in columns])
        table.setRowCount(len(rows))
        for index, values in enumerate(rows):
            for column, value in enumerate(values):
                item = cell(value)
                # Every cell can be asked what the table it sits in is not, because
                # a cell is what gets copied into a message and read on its own.
                item.setToolTip("\n\n".join(filter(None, [
                    str(value) if value is not None else "Not available",
                    plasmid_evidence.TABLE_NOTES[key]])))
                table.setItem(index, column, item)
        table.setSortingEnabled(sorting)
        self.empties[key].setText(spec.get("empty_reason") or "")
        self.empties[key].setVisible(bool(spec.get("empty_reason")))
        return table

    # --- what the report consumes --------------------------------------------
    def report_payload(self):
        """The tables this page is showing, as plain data, for a report to print.

        Returns ``None`` when nothing has been built, so a report says that the
        plasmid page was never read rather than printing empty tables a reader
        would take for a cohort with no plasmid evidence in it.
        """
        return self.payload


def _database_words(database):
    """The database state as the one line a reader must not be able to miss."""
    status = str(database.get("status") or "unknown")
    if status == "installed":
        return "Plasmid reference set installed: " + ", ".join(database.get("installed") or ())
    if status == "not_installed":
        return "No plasmid database is installed — nothing below was searched for"
    return "Whether a plasmid database is installed is unknown on this page"
