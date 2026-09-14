"""Shared native table and sample-presentation helpers, independent of controllers."""

from PySide6.QtCore import QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLayout,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)


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
    row.setContentsMargins(12, 4, 12, 4)
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


def make_table(headers, multiple=True):
    widget = QTableWidget(0, len(headers))
    widget.setHorizontalHeaderLabels(headers)
    widget.verticalHeader().hide()
    widget.verticalHeader().setDefaultSectionSize(39)
    widget.setShowGrid(False)
    widget.setAlternatingRowColors(True)
    widget.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    widget.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection if multiple
                            else QTableWidget.SelectionMode.SingleSelection)
    widget.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    widget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    widget.horizontalHeader().setStretchLastSection(True)
    widget.horizontalHeader().setDefaultSectionSize(145)
    widget.setSortingEnabled(True)
    return widget


def cell(value, sample_id=None):
    item = QTableWidgetItem(str(value) if value is not None else "—")
    item.setData(Qt.ItemDataRole.UserRole, sample_id)
    item.setToolTip(str(value) if value is not None else "Not available")
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
