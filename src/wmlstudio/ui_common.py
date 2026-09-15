"""Shared native table and sample-presentation helpers, independent of controllers."""

import re

from PySide6.QtCore import QPoint, QRect, QSize, Qt, QTimer
from PySide6.QtGui import QFontMetrics
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLayout,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)

from wmlstudio import theme


class FlowLayout(QLayout):
    """Wrap action strips at the available width instead of clipping controls."""

    def __init__(self, parent=None, spacing=8):
        super().__init__(parent)
        self.items = []
        self.setContentsMargins(0, 0, 0, 0)
        self.setSpacing(spacing)

    def addItem(self, item):
        self.items.append(item)

    def addWidget(self, widget, stretch=0, alignment=Qt.AlignmentFlag(0)):
        super().addWidget(widget)

    def addStretch(self, stretch=0):
        pass

    def addLayout(self, layout, stretch=0):
        wrapper = QWidget()
        wrapper.setLayout(layout)
        self.addWidget(wrapper)

    def count(self):
        return len(self.items)

    def itemAt(self, index):
        return self.items[index] if 0 <= index < len(self.items) else None

    def takeAt(self, index):
        return self.items.pop(index) if 0 <= index < len(self.items) else None

    def expandingDirections(self):
        return Qt.Orientation(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self.arrange(QRect(0, 0, width, 0), True)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self.arrange(rect, False)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        size = QSize()
        for item in self.items:
            size = size.expandedTo(item.minimumSize())
        return size

    def arrange(self, rect, measure):
        x, y, height = rect.x(), rect.y(), 0
        for item in self.items:
            if item.isEmpty():
                continue
            hint = item.sizeHint()
            width = min(hint.width(), max(item.minimumSize().width(), rect.width()))
            if x > rect.x() and x + width > rect.right() + 1:
                x, y, height = rect.x(), y + height + self.spacing(), 0
            if not measure:
                item.setGeometry(QRect(QPoint(x, y), QSize(width, hint.height())))
            x += width + self.spacing()
            height = max(height, hint.height())
        return y + height - rect.y()


def workspace_header(window, layout, key, *, index=None):
    """The orientation strip for one tab: what it is for, and one obvious next step.

    Someone who is not a bioinformatician should be able to read a tab and know what
    it answers without opening the guide, so every page carries the same three things
    in the same place: the plain-language purpose, the next action, and the guide.

    The strip is inserted rather than built into `build_*`, because the evidence and
    scheme pages rearrange their own layout items by position; adding this after the
    pages are built keeps that surgery working untouched.
    """
    from wmlstudio.ui_tabs import NEXT_STEP, PAGE_PURPOSE
    from wmlstudio.widgets import button, label
    from wmlstudio.workspace_focus import COHORT_ATTRIBUTES, cohort_bar

    purpose_text = PAGE_PURPOSE.get(key, "")
    strip = QFrame()
    strip.setObjectName("purposeStrip")
    strip.setProperty("pageKey", key)
    row = QHBoxLayout(strip)
    row.setContentsMargins(*theme.density_metrics(table_density())["strip_margins"])
    row.setSpacing(8)
    purpose = label(purpose_text, "purpose", True)
    purpose.setToolTip(purpose_text)
    row.addWidget(purpose, 1)
    if key in COHORT_ATTRIBUTES and getattr(window, "focus", None) is not None:
        # A tab that owns a cohort says here what it is reviewing and where that
        # came from. The focus is only ever copied in by the button beside it.
        strip.cohort_bar = cohort_bar(window, key, compact=True)
        row.addWidget(strip.cohort_bar)
    text, method = NEXT_STEP.get(key, ("", ""))
    # Several next-step methods are added by later work; a missing one must not
    # raise while the window is being built, so the button simply does not appear.
    handler = getattr(window, method, None) if method else None
    if text and callable(handler):
        step = button(text, handler)
        step.setObjectName("nextStep")
        step.setToolTip(f"The usual next step on this tab: {text}")
        row.addWidget(step)
    guide = getattr(window, "open_workflow_guide", None)
    if callable(guide):
        help_button = button("?", lambda checked=False, topic=key: guide(topic))
        help_button.setObjectName("pageGuide")
        help_button.setAccessibleName(f"Guide for {key}")
        help_button.setToolTip("Open the problem → solution guide for this tab")
        row.addWidget(help_button)
    if index is None:
        layout.addWidget(strip)
    else:
        layout.insertWidget(index, strip)
    return strip


# --- table density ---------------------------------------------------------
# One chosen density for the whole application. Tables are built by about thirty
# call sites that must not each grow a settings lookup, so the choice is held
# here and `make_table` reads it. The window sets it once from the saved
# preference before it builds any page.
_density = {"name": theme.DEFAULT_DENSITY}

# What a header section costs beyond its text: the section padding on both sides
# plus room for the sort indicator, which every one of these tables can show.
HEADER_CHROME = 30
# No column is opened wider than this. A very long cell would otherwise push
# every other column off the screen, which is the truncation problem again.
MAXIMUM_COLUMN = 320
# Wider tables are left alone: measuring contents across a cgMLST-sized matrix
# costs more than the tidier columns are worth.
AUTOFIT_COLUMN_LIMIT = 30


def table_density() -> str:
    return _density["name"]


def set_table_density(name) -> str:
    """Choose the density new tables are built at; returns the name that took."""
    _density["name"] = theme.DEFAULT_DENSITY if name is None else str(name).strip().lower()
    if _density["name"] not in theme.TABLE_DENSITY:
        _density["name"] = theme.DEFAULT_DENSITY
    return _density["name"]


def header_width(table, column, *, maximum=MAXIMUM_COLUMN) -> int:
    """How wide a column must be for its own heading to be readable in full.

    A truncated heading ("R evidence st" for "AMR evidence state") is not a
    cosmetic problem: the reader cannot tell which column they are looking at,
    which is the one thing a table of evidence must never leave in doubt.
    """
    heading = table.model().headerData(column, Qt.Orientation.Horizontal)
    metrics = QFontMetrics(table.horizontalHeader().font())
    return min(int(maximum), metrics.horizontalAdvance(str(heading or "")) + HEADER_CHROME)


def _columns(table) -> int:
    """How many columns a view has, whether or not it owns its own items.

    `restyle_tables` reaches plain QTableViews as well as QTableWidgets, and only
    the latter has `columnCount`; the model answers for both.
    """
    model = table.model()
    return 0 if model is None else model.columnCount()


def apply_table_density(table, density=None) -> str:
    """Row height, grid and column widths for one table at the chosen density."""
    name = table_density() if density is None else str(density).strip().lower()
    if name not in theme.TABLE_DENSITY:
        name = theme.DEFAULT_DENSITY
    metrics = theme.density_metrics(name)
    table.setProperty("tableDensity", name)
    rows = table.verticalHeader()
    # Qt's own minimum is derived from the font and is taller than a compact row,
    # so without lowering it first the chosen height is silently ignored.
    rows.setMinimumSectionSize(min(16, metrics["row_height"]))
    rows.setDefaultSectionSize(metrics["row_height"])
    table.setShowGrid(bool(metrics["grid"]))
    header = table.horizontalHeader()
    header.setDefaultSectionSize(metrics["column_width"])
    for column in range(_columns(table)):
        header.resizeSection(column, max(metrics["column_width"],
                                         header_width(table, column)))
    return name


def fit_table_columns(table, *, maximum=MAXIMUM_COLUMN, grow_only=True) -> None:
    """Widen columns to what they actually hold, without ever hiding a heading.

    Only ever widens by default, so a caller that deliberately narrowed a
    checkbox column keeps it. Measurement is capped at a sample of rows, so this
    stays cheap on a table that has just been filled with hundreds of isolates.
    """
    if _columns(table) > AUTOFIT_COLUMN_LIMIT:
        return
    header = table.horizontalHeader()
    header.setResizeContentsPrecision(50)
    for column in range(_columns(table)):
        wanted = max(header_width(table, column, maximum=maximum),
                     min(int(maximum), table.sizeHintForColumn(column) + 12))
        if grow_only:
            wanted = max(wanted, header.sectionSize(column))
        header.resizeSection(column, wanted)


def _fit_later(table) -> None:
    """Coalesce the refit: a fill inserts rows one at a time, but fits once."""
    try:
        if table.property("awaitingColumnFit"):
            return
        table.setProperty("awaitingColumnFit", True)
    except RuntimeError:
        return

    def run():
        try:
            table.setProperty("awaitingColumnFit", False)
            fit_table_columns(table)
        except RuntimeError:
            # The table was destroyed between the fill and this callback, which
            # happens whenever a dialog is closed while it is still filling.
            return

    QTimer.singleShot(0, run)


def make_table(headers, multiple=True, *, density=None, autofit=True):
    widget = QTableWidget(0, len(headers))
    widget.setHorizontalHeaderLabels(headers)
    # The row numbers down the left say nothing: the isolate is named in the
    # first column, and the strip costs width on every table in the application.
    widget.verticalHeader().hide()
    widget.setAlternatingRowColors(True)
    widget.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    widget.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection if multiple
                            else QTableWidget.SelectionMode.SingleSelection)
    widget.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    widget.setWordWrap(False)
    widget.setTextElideMode(Qt.TextElideMode.ElideRight)
    widget.setCornerButtonEnabled(False)
    header = widget.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    header.setStretchLastSection(True)
    # Bold-on-click headers reflow the heading and can elide it mid-reading.
    header.setHighlightSections(False)
    widget.setSortingEnabled(True)
    apply_table_density(widget, density)
    if autofit:
        model = widget.model()
        model.rowsInserted.connect(lambda *_, table=widget: _fit_later(table))
        model.modelReset.connect(lambda table=widget: _fit_later(table))
    return widget


def restyle_tables(root, density=None) -> int:
    """Re-apply the density to tables that are already on screen. Returns how many."""
    tables = root.findChildren(QTableView)
    for table in tables:
        apply_table_density(table, density)
    return len(tables)


_NUMBER = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\s*%?$")


def is_numeric_text(value) -> bool:
    """Whether a cell holds a quantity, so the column can line its digits up."""
    if isinstance(value, bool) or value is None:
        return False
    if isinstance(value, (int, float)):
        return True
    return bool(_NUMBER.match(str(value).strip().replace(",", "")))


def cell(value, sample_id=None, *, align=None):
    item = QTableWidgetItem(str(value) if value is not None else "—")
    item.setData(Qt.ItemDataRole.UserRole, sample_id)
    item.setToolTip(str(value) if value is not None else "Not available")
    if align is None and is_numeric_text(value):
        # Quantities are read by comparing them down the column, which only works
        # when the digits line up. Text stays where the reader expects it.
        align = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
    if align is not None:
        item.setTextAlignment(align)
    return item


def organism_for(sample):
    metadata = sample.get("metadata") or {}
    assigned = metadata.get("organism") or {}
    if isinstance(assigned, str):
        parts = assigned.split(maxsplit=1)
        assigned = {"genus": parts[0] if parts else "", "species": parts[1] if len(parts) > 1 else ""}
    if assigned.get("genus") or assigned.get("species"):
        return assigned.get("genus", ""), assigned.get("species", ""), "Assigned"
    result = sample.get("result") or {}
    identification = result.get("identification") or {}
    detected = identification.get("organism") or {}
    genus = detected.get("genus") or identification.get("genus") or ""
    species = detected.get("species") or identification.get("species") or ""
    return genus, species, "Provisional" if genus else "Unknown"


def gene_names(sample):
    from wmlstudio.sample_workflow import current_hydra_evidence
    evidence = current_hydra_evidence(sample)
    return sorted({str(hit.get("gene")) for hit in evidence.get("hits", [])
                   if hit.get("gene") and hit.get("element_type") == "AMR"
                   and hit.get("primary") is True})


def flattened_metadata(sample):
    output = {}

    def flatten(value, prefix=""):
        if isinstance(value, dict):
            for key, child in value.items():
                # Evidence arrays belong in the drill-down, not a multi-megabyte cell.
                if key not in {"hits", "provenance", "execution_provenance", "analyses", "upstream", "mate_record"}:
                    flatten(child, f"{prefix}.{key}" if prefix else str(key))
        elif isinstance(value, list):
            output[prefix] = "; ".join(str(v) for v in value if not isinstance(v, dict))
        else:
            output[prefix] = value

    flatten(sample.get("metadata") or {})
    return output
