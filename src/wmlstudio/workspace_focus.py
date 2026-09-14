"""Shared, temporary focus across workspace tabs — never an automatic cohort.

Selecting isolates in one tab should help you find them in another. It must never
quietly widen an analysis, a comparison or a report, because those cohorts are the
record of what was actually reviewed.

So the model is deliberately two-layered:

* **Focus** is transient and presentational. It is shared by every tab, it is never
  saved to the project, never exported, and it decides nothing.
* **A cohort** belongs to one tab, is written only when the user asks for it, and
  keeps a visible note of where it came from.

`adopt_focus` and `send_selection` are the only functions here that write a cohort,
and both are user-triggered. Nothing in this module observes a selection and copies
it into a cohort on its own.
"""

from __future__ import annotations

from PySide6.QtCore import QObject, QTime, Signal
from PySide6.QtWidgets import QHBoxLayout, QWidget

from wmlstudio.widgets import button, label

# The four tabs that own an explicit cohort, and the window attribute each uses.
COHORT_ATTRIBUTES = {
    "isolates": "selection_ids",
    "compare": "cohort_ids",
    "evidence": "feature_ids",
    "reports": "report_ids",
}

# The dialog each tab opens to choose its own cohort, by method name on the window.
CHOOSE_METHODS = {
    "compare": "choose_comparison_cohort",
    "evidence": "choose_feature_cohort",
    "reports": "choose_report_cohort",
}

# What to refresh once a cohort changed, by method name on the window.
REFRESH_METHODS = {
    "isolates": "refresh_tables",
    "compare": "refresh_comparison",
    "evidence": "refresh_features",
    "reports": "refresh_report_table",
}

COHORT_NOUNS = {
    "isolates": "selected",
    "compare": "in this comparison",
    "evidence": "in this evidence review",
    "reports": "in this report",
}

# A focused row is a temporary view state. It must look different from the saved
# metadata.cluster highlight (#343348), which is evidence the user recorded.
FOCUS_ACCENT = "#153B3A"
SAVED_HIGHLIGHT_ACCENT = "#343348"
FOCUS_LEGEND = ("Teal edge = isolates you are focused on right now (temporary, shared between "
                "tabs).\nPurple = a highlight group you saved with these isolates.")

EMPTY_COHORT_TEXT = ("No isolates chosen yet for this tab — nothing here is inherited from "
                     "another tab.")


class FocusBus(QObject):
    """Transient, presentational cross-tab focus. Never a cohort, never persisted."""

    focusChanged = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.ids: set[str] = set()
        self.origin = ""
        self.at = ""

    @property
    def count(self) -> int:
        return len(self.ids)

    def set_focus(self, ids, origin) -> None:
        """Record what the user just selected, and where. Changes no cohort."""
        chosen = {str(value) for value in (ids or ()) if value}
        origin = str(origin or "")
        if chosen == self.ids and origin == self.origin:
            return
        self.ids = chosen
        self.origin = origin
        self.at = QTime.currentTime().toString("HH:mm") if chosen else ""
        self.focusChanged.emit()

    def clear(self) -> None:
        if not (self.ids or self.origin):
            return
        self.ids, self.origin, self.at = set(), "", ""
        self.focusChanged.emit()

    def prune(self, valid_ids) -> None:
        """Drop isolates that left the project, so focus cannot name a missing row."""
        valid = {str(value) for value in (valid_ids or ())}
        remaining = self.ids & valid
        if remaining == self.ids:
            return
        self.ids = remaining
        if not remaining:
            self.origin, self.at = "", ""
        self.focusChanged.emit()

    def describe(self) -> str:
        if not self.ids:
            return "No isolates focused."
        where = f" · from {self.origin}" if self.origin else ""
        when = f" at {self.at}" if self.at else ""
        return f"Focus · {len(self.ids)} isolate{'s' if len(self.ids) != 1 else ''}{where}{when}"


class CohortLedger(QObject):
    """Where each tab's cohort came from, so a reader can see how it was chosen."""

    changed = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.entries: dict[str, dict] = {}

    def record(self, key, ids, origin) -> dict:
        entry = {"count": len(set(ids or ())), "origin": str(origin or ""),
                 "at": QTime.currentTime().toString("HH:mm")}
        self.entries[str(key)] = entry
        self.changed.emit(str(key))
        return entry

    def forget(self, key) -> None:
        if self.entries.pop(str(key), None) is not None:
            self.changed.emit(str(key))

    def entry(self, key) -> dict:
        return self.entries.get(str(key), {})

    def describe(self, key, count=None) -> str:
        """One line naming what this tab is reviewing and where that came from."""
        entry = self.entry(key)
        total = entry.get("count", 0) if count is None else int(count)
        if not total:
            return EMPTY_COHORT_TEXT
        noun = COHORT_NOUNS.get(str(key), "in this tab")
        line = f"Reviewing {total} isolate{'s' if total != 1 else ''} {noun}"
        origin = entry.get("origin")
        if origin and entry.get("count") == total:
            line += f" · chosen from {origin}"
            if entry.get("at"):
                line += f" at {entry['at']}"
        return line + " · this cohort belongs to this tab only."


def cohort_ids(window, key) -> set[str]:
    """Read one tab's own cohort without assuming which attribute holds it."""
    attribute = COHORT_ATTRIBUTES.get(str(key))
    return set(getattr(window, attribute, set()) or ()) if attribute else set()


def _call(window, name, *args):
    method = getattr(window, name, None) if name else None
    if callable(method):
        method(*args)
        return True
    return False


def _notify(window, message):
    return _call(window, "notify", message)


def adopt_focus(window, key, focus=None, *, ledger=None, ids=None, origin=None) -> set[str]:
    """Copy the shared focus into one tab's own cohort, on an explicit user action.

    Returns the ids that were adopted, or an empty set when there was nothing to
    adopt. The cohort of every other tab is left exactly as it was.
    """
    key = str(key)
    attribute = COHORT_ATTRIBUTES.get(key)
    if attribute is None:
        raise ValueError(f"'{key}' does not own a cohort.")
    focus = getattr(window, "focus", None) if focus is None else focus
    chosen = {str(value) for value in (ids if ids is not None else getattr(focus, "ids", ()) or ())}
    if not chosen:
        _notify(window, "Select isolates first. Nothing is included automatically.")
        return set()
    if origin is None:
        origin = getattr(focus, "origin", "") or "Current focus"
    setattr(window, attribute, set(chosen))
    _remember(window, key, chosen, origin, ledger)
    _persist(window, key, chosen)
    _call(window, REFRESH_METHODS.get(key))
    _notify(window, f"{len(chosen)} isolates adopted from {origin}. "
                    "This cohort belongs to this tab only.")
    return set(chosen)


def send_selection(window, target, ids=None, origin="Selection", *, ledger=None,
                   navigate=True) -> set[str]:
    """Give one tab a cohort and go there. The origin is recorded and displayed."""
    target = str(target)
    if target not in COHORT_ATTRIBUTES:
        raise ValueError(f"'{target}' does not own a cohort.")
    focus = getattr(window, "focus", None)
    chosen = {str(value) for value in (ids if ids is not None else getattr(focus, "ids", ()) or ())}
    if not chosen:
        _notify(window, "Select isolates first. Nothing is included automatically.")
        return set()
    setattr(window, COHORT_ATTRIBUTES[target], set(chosen))
    _remember(window, target, chosen, origin, ledger)
    _persist(window, target, chosen)
    if target == "reports":
        snapshot = getattr(window, "_current_snapshot", None)
        if snapshot is not None:
            window._report_investigation_snapshot = snapshot
    _call(window, REFRESH_METHODS.get(target))
    if navigate:
        _call(window, "navigate", _page_index(window, target))
    _notify(window, f"{len(chosen)} isolates sent to {target}. "
                    "This cohort belongs to that tab only.")
    return set(chosen)


def _remember(window, key, ids, origin, ledger=None):
    ledger = getattr(window, "cohort_origins", None) if ledger is None else ledger
    if isinstance(ledger, CohortLedger):
        ledger.record(key, ids, origin)


def _persist(window, key, ids):
    """Only the comparison cohort is a saved project setting today; keep that true."""
    if key != "compare":
        return
    project = getattr(window, "project", None)
    if project is not None and hasattr(project, "set_setting"):
        project.set_setting("comparison_cohort", sorted(ids))


def _page_index(window, key):
    mapping = getattr(window, "page_index", None)
    if isinstance(mapping, dict) and key in mapping:
        return mapping[key]
    pages = getattr(window, "pages", None)
    index = pages.index_of(key) if hasattr(pages, "index_of") else -1
    return index


class FocusBar(QWidget):
    """Top-strip line: how many isolates are focused, and where they came from."""

    showRequested = Signal()
    clearRequested = Signal()

    def __init__(self, focus, parent=None):
        super().__init__(parent)
        self.focus = focus
        self.project_count = 0
        self.project_text = "0 isolates in project"
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)
        self.text = label("0 isolates in project", "small")
        self.text.setToolTip("Focus is a temporary highlight shared by the tabs. "
                             "Each analysis, comparison and report keeps its own cohort.")
        row.addWidget(self.text)
        self.show_button = button("Show in Isolates", lambda: self.showRequested.emit())
        self.clear_button = button("Clear focus", lambda: self.clearRequested.emit())
        for control in (self.show_button, self.clear_button):
            control.hide()
            row.addWidget(control)
        if focus is not None:
            focus.focusChanged.connect(self.refresh)
        self.refresh()

    def set_project_count(self, count) -> None:
        self.project_count = int(count)
        self.project_text = f"{self.project_count} isolates in project"
        self.refresh()

    def setText(self, text) -> None:      # noqa: N802 - drop-in for the label it replaced
        """Compatibility with the plain scope label this bar replaced in the top strip.

        The window still writes the project count here; the bar shows it whenever
        nothing is focused, and shows the focus line when something is.
        """
        self.project_text = str(text)
        self.refresh()

    def refresh(self) -> None:
        focused = bool(getattr(self.focus, "ids", set()))
        self.text.setText(self.focus.describe() if focused else self.project_text)
        for control in (self.show_button, self.clear_button):
            control.setVisible(focused)


class CohortBar(QWidget):
    """One row under a tab header: what it reviews, where that came from, how to change it.

    `compact` is for hosts that already carry the tab's own chooser — the orientation
    strip does — so the bar shows only the scope line and the one button that turns
    focus into this tab's cohort. The full sentence stays in the tooltip; it is the
    honest one, and it must be readable wherever the short form is shown.
    """

    def __init__(self, window, key, focus=None, *, ledger=None, extra=None, compact=False,
                 parent=None):
        super().__init__(parent)
        self.window_ref = window
        self.key = str(key)
        self.compact = bool(compact)
        self.setObjectName("cohortBar")
        self.focus = getattr(window, "focus", None) if focus is None else focus
        self.ledger = getattr(window, "cohort_origins", None) if ledger is None else ledger
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)
        self.text = label(EMPTY_COHORT_TEXT, "cohortScope" if self.compact else "small",
                          not self.compact)
        row.addWidget(self.text, 0 if self.compact else 1)
        if extra is not None:
            row.addWidget(extra)
        chooser = CHOOSE_METHODS.get(self.key)
        if chooser and not self.compact and hasattr(window, chooser):
            row.addWidget(button("Choose isolates…", getattr(window, chooser)))
        self.adopt_button = button("Use current focus", self.adopt)
        if self.compact:
            self.adopt_button.setObjectName("cohortAdopt")
        self.adopt_button.setToolTip("Copy the isolates you are focused on into this tab's "
                                     "own cohort. Nothing is copied until you click.")
        row.addWidget(self.adopt_button)
        if self.focus is not None:
            self.focus.focusChanged.connect(self.refresh)
        if isinstance(self.ledger, CohortLedger):
            self.ledger.changed.connect(lambda key: self.refresh())
        self.refresh()

    def adopt(self) -> None:
        adopt_focus(self.window_ref, self.key, self.focus, ledger=self.ledger)
        self.refresh()

    def describe(self) -> str:
        count = len(cohort_ids(self.window_ref, self.key))
        if isinstance(self.ledger, CohortLedger):
            return self.ledger.describe(self.key, count)
        if count:
            noun = COHORT_NOUNS.get(self.key, "in this tab")
            return (f"Reviewing {count} isolate{'s' if count != 1 else ''} {noun} · "
                    "this cohort belongs to this tab only.")
        return EMPTY_COHORT_TEXT

    def short_text(self) -> str:
        count = len(cohort_ids(self.window_ref, self.key))
        if not count:
            return "Nothing chosen here yet"
        entry = self.ledger.entry(self.key) if isinstance(self.ledger, CohortLedger) else {}
        # Name the origin only while it still describes the cohort on screen; a
        # cohort edited since was not chosen that way.
        origin = entry.get("origin", "") if entry.get("count") == count else ""
        return f"Reviewing {count} · from {origin}" if origin else f"Reviewing {count}"

    def refresh(self) -> None:
        full = self.describe()
        self.text.setText(self.short_text() if self.compact else full)
        self.text.setToolTip(full)
        focused = len(getattr(self.focus, "ids", set()) or ())
        self.adopt_button.setText(f"Use current focus ({focused})" if focused
                                  else "Use current focus")
        self.adopt_button.setEnabled(bool(focused))


def cohort_bar(window, key, focus=None, *, ledger=None, extra=None, compact=False) -> CohortBar:
    """Build the explicit focus-to-cohort bridge for one tab."""
    return CohortBar(window, key, focus, ledger=ledger, extra=extra, compact=compact)
