"""Top-level workspace tabs. One page per workflow question, addressed by a stable key.

The window used to hold its pages in a QStackedWidget addressed by integer index,
so every navigation site hard-coded a number. WorkspaceTabs keeps those indices
working while giving each page a stable key, a visible tab, a purpose line and an
optional "about to be shown" hook, so a later renumbering cannot silently send the
user to the wrong workspace.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTabWidget, QWidget

PAGE_KEYS = ("overview", "isolates", "compare", "schemes", "evidence", "reports", "settings")

# Kept byte-identical to the sidebar glyphs the pages were built with, so the tab
# bar and the View menu name the same seven places.
NAV_SYMBOLS = ("◫", "▤", "⌘", "▥", "◈", "↗", "⚙")

TAB_LABELS = {
    "overview": "Overview",
    "isolates": "Isolates",
    "compare": "Compare",
    "schemes": "Schemes",
    "evidence": "Evidence",
    "reports": "Reports",
    "settings": "Settings",
}

# One plain sentence per tab. A tab that carries a scientific claim states its
# limit in the same sentence; the interface must never imply more than the data.
PAGE_PURPOSE = {
    "overview": "See where this investigation stands and what to do next.",
    "isolates": "Keep every isolate, its organism and its files in one place.",
    "compare": "See which isolates have similar allele profiles. Similarity is not proof of transmission.",
    "schemes": "Manage the local allele databases your typing runs against.",
    "evidence": "Review identity, resistance and virulence evidence for isolates you choose. Genotype is not measured susceptibility.",
    "reports": "Turn reviewed evidence into a document you can share.",
    "settings": "Make text size, screen scaling and data locations comfortable.",
}

# (button text, name of the MainWindow method the button calls). Resolved with
# getattr at build time so a missing method degrades to a hidden button.
NEXT_STEP = {
    "overview": ("Import sequences…", "browse_files"),
    "isolates": ("Import…", "browse_files"),
    "compare": ("Choose cohort…", "choose_comparison_cohort"),
    "schemes": ("Browse online / install updates…", "open_reference_manager"),
    "evidence": ("Choose isolates…", "choose_feature_cohort"),
    "reports": ("Choose report isolates…", "choose_report_cohort"),
    "settings": ("Try this size on a sample row", "preview_interface_scale"),
}


def page_key(index: int) -> str:
    """The stable key of a build-order index, or '' when the index is unknown."""
    return PAGE_KEYS[index] if 0 <= index < len(PAGE_KEYS) else ""


def tab_title(key: str) -> str:
    """Tab bar text: the page glyph plus its short label."""
    index = PAGE_KEYS.index(key) if key in PAGE_KEYS else -1
    symbol = NAV_SYMBOLS[index] if index >= 0 else ""
    return f"{symbol}  {TAB_LABELS.get(key, key.title())}".strip()


class WorkspaceTabs(QTabWidget):
    """Page container addressed by key, still addressable by the old integer index."""

    pageShown = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("workspaceTabs")
        self.setDocumentMode(False)
        self.setMovable(False)
        self.setUsesScrollButtons(True)
        self.tabBar().setObjectName("workspaceTabBar")
        self.tabBar().setExpanding(False)
        self.tabBar().setElideMode(Qt.TextElideMode.ElideRight)
        self.tabBar().setAccessibleName("Workspace sections")
        self.setAccessibleName("Workspace sections")
        # Native child viewports must paint directly; a transparent stack caused
        # stale or blank pages on Windows.
        self.setAutoFillBackground(True)
        self._keys: list[str] = []
        self._hooks: dict[str, object] = {}
        self._subtabs: dict[str, QTabWidget] = {}
        self._depth = 0
        self.currentChanged.connect(self._current_changed)

    # --- building -----------------------------------------------------------
    def add_page(self, widget, key=None, *, title=None, tooltip=None, purpose=None) -> int:
        """Add a page and return its index. The key defaults to the build order."""
        if key is None:
            key = page_key(self.count()) or f"page{self.count()}"
        index = self.addTab(widget, title if title is not None else tab_title(key))
        self._keys.insert(index, str(key))
        self.tabBar().setTabData(index, str(key))
        if tooltip is None:
            tooltip = purpose if purpose is not None else PAGE_PURPOSE.get(key, "")
        self.setTabToolTip(index, str(tooltip))
        self.setTabWhatsThis(index, str(tooltip))
        return index

    def addWidget(self, widget) -> int:
        """QStackedWidget compatibility: keep `pages.addWidget(page)` call sites working."""
        return self.add_page(widget)

    # --- addressing ---------------------------------------------------------
    def keys(self) -> tuple[str, ...]:
        return tuple(self._keys)

    def key_at(self, index) -> str:
        return self._keys[index] if isinstance(index, int) and 0 <= index < len(self._keys) else ""

    def index_of(self, target) -> int:
        """Resolve a key or an integer index to a page index, or -1 when unknown."""
        if isinstance(target, bool):
            return -1
        if isinstance(target, int):
            return target if 0 <= target < self.count() else -1
        return self._keys.index(target) if target in self._keys else -1

    def page_index(self) -> dict[str, int]:
        """{key: index} for the sites that still need a number."""
        return {key: index for index, key in enumerate(self._keys)}

    def current_key(self) -> str:
        return self.key_at(self.currentIndex())

    def page_for(self, target) -> QWidget | None:
        index = self.index_of(target)
        return self.widget(index) if index >= 0 else None

    def show_page(self, target) -> bool:
        """Show a page by key or index. Returns False when the page does not exist."""
        index = self.index_of(target)
        if index < 0:
            return False
        self.setCurrentIndex(index)
        return True

    # --- per-page show hook -------------------------------------------------
    def set_show_hook(self, key, callback) -> None:
        """Run `callback` whenever that page becomes visible; None removes the hook."""
        if callback is None:
            self._hooks.pop(str(key), None)
        else:
            self._hooks[str(key)] = callback

    def show_hook(self, key):
        return self._hooks.get(str(key))

    def _current_changed(self, index):
        # A hook may legitimately redirect to another page, so nesting is allowed;
        # the depth limit only stops two hooks from bouncing off each other forever.
        key = self.key_at(index)
        if not key or self._depth >= 3:
            return
        self._depth += 1
        try:
            self.pageShown.emit(key)
            hook = self._hooks.get(key)
            if hook is not None:
                hook()
        finally:
            self._depth -= 1

    # --- sub-tab registry ---------------------------------------------------
    def register_subtabs(self, key, tabs) -> None:
        """Record the inner QTabWidget a page owns, so later work can extend it."""
        self._subtabs[str(key)] = tabs

    def subtabs(self, key):
        return self._subtabs.get(str(key))

    def subtab_keys(self) -> tuple[str, ...]:
        return tuple(self._subtabs)

    def add_subtab(self, key, widget, title, *, index=None) -> int:
        """Insert a sub-tab into a page's inner tabs. Returns -1 when it has none."""
        tabs = self._subtabs.get(str(key))
        if tabs is None:
            return -1
        if index is None:
            return tabs.addTab(widget, title)
        return tabs.insertTab(int(index), widget, title)

    def show_subtab(self, key, title) -> bool:
        tabs = self._subtabs.get(str(key))
        if tabs is None:
            return False
        for index in range(tabs.count()):
            if tabs.tabText(index) == title:
                tabs.setCurrentIndex(index)
                return True
        return False
