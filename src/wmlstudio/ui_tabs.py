"""Top-level workspace tabs: the pipeline, in the order the work is actually done.

The window used to hold its pages in a QStackedWidget addressed by integer index,
so every navigation site hard-coded a number. WorkspaceTabs keeps those numbers
working while giving each page a stable key, a visible tab, a purpose line, a
Clear action and an optional "about to be shown" hook.

Two orders live here and they are deliberately different:

* ``PAGE_KEYS`` is the **build order**. It is the compatibility address space: an
  integer anywhere in this application -- ``navigate(2)``, ``pages.widget(3)``,
  ``page_index[key]`` -- means "the Nth page that was built", never "the Nth tab
  in the bar". Those numbers must keep their meaning, so nothing is ever inserted
  into the middle of this tuple; new pages are appended.
* ``PIPELINE`` is the **display order**: the tab bar reads as the user's own
  workflow -- load the samples, look at the reads, assemble, type, draw the tree,
  do it again for cgMLST, screen for resistance, report, update, settle the
  settings. Moving a tab in the bar therefore costs nothing elsewhere.

A station that a later round builds is listed in ``PLANNED``. Its tab exists so
the pipeline has no gap where a step belongs, and its page says plainly that
nothing runs there yet and where that work happens today. A tab must never look
finished when it is not.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTabWidget, QWidget

#: Build order. The integer address space every existing call site uses.
#: Indices 0-6 are the seven pages this application shipped with and must not move.
PAGE_KEYS = (
    "overview", "isolates", "compare", "schemes", "evidence", "reports", "settings",
    # Appended by the pipeline layout. Their numbers are new, so nothing depends
    # on them; they exist so `navigate(<int>)` can still reach every page.
    "reads", "assembly", "mlst", "cgmlst", "cgmlst_tree", "snp",
)

#: Display order: one dedicated tab per task, in the order the tasks happen.
PIPELINE = (
    "overview",     # where this investigation stands
    "isolates",     # Samples -- the hub everything else works from
    "reads",        # planned: read trimming and quality filtering (fastp)
    "assembly",     # planned: read pairs to assemblies
    "mlst",         # classical seven-locus typing
    "compare",      # the MLST minimum spanning tree
    "cgmlst",       # core-genome typing and its table of calls
    "cgmlst_tree",  # the cgMLST minimum spanning tree -- its own scale, always
    "snp",          # planned: a SNP tree from read alignment
    "evidence",     # HYDRA: resistance, virulence and lineage evidence
    "reports",      # the report
    "schemes",      # Update: scheme libraries and reference databases
    "settings",
)

#: Stations a later round builds. Their pages say so; they never pretend to work.
PLANNED = ("reads", "assembly", "snp")

#: Pages this layout adds on top of the seven the window has always built.
STATION_KEYS = PAGE_KEYS[7:]

TAB_LABELS = {
    "overview": "Overview",
    # Samples is the hub: every other page acts on a cohort chosen here, and
    # "samples" is the word the people using this application use.
    "isolates": "Samples",
    "reads": "Read QC",
    "assembly": "Assembly",
    "mlst": "MLST",
    # Seven-locus differences only. The cgMLST tree is a separate tab on purpose.
    "compare": "MLST tree",
    "cgmlst": "cgMLST",
    "cgmlst_tree": "cgMLST tree",
    "snp": "SNP tree",
    "evidence": "HYDRA",
    "reports": "Report",
    # Schemes are reference data you install and update, so they live here with
    # the AMR databases rather than in a second place of their own.
    "schemes": "Update",
    "settings": "Settings",
}

#: The page each key folds into once its content also exists as a sub-tab there,
#: with the sub-tab title to select. `WorkspaceTabs.fold_page` applies one: the
#: top-level tab stops being shown, every index keeps its meaning, and navigating
#: to the key still lands on the same content. Nothing folds today -- the point of
#: this layout is that each task has its own tab -- but a page that later grows a
#: second home must have one place, not two.
FOLDABLE: dict[str, tuple[str, str]] = {}

# One plain sentence per tab. A tab that carries a scientific claim states its
# limit in the same sentence; the interface must never imply more than the data.
PAGE_PURPOSE = {
    "overview": "See where this investigation stands and what to do next.",
    "isolates": "Every sample, its organism and its files in one place — the hub every other "
                "tab works from.",
    "reads": "Read trimming and quality filtering. Planned for a later round; nothing runs "
             "on this tab yet.",
    "assembly": "Read pairs into assemblies, with their own metrics. Planned for a later round; "
                "today assemblies are made in Samples.",
    "mlst": "Classical seven-locus MLST: one sequence type per isolate from a curated allele panel.",
    "compare": "A minimum spanning tree of seven-locus allele differences. An MST is not a "
               "phylogeny, and similarity is not proof of transmission.",
    "cgmlst": "Core-genome MLST: hundreds to thousands of targets, with its own scheme and its "
              "own missing-target count.",
    "cgmlst_tree": "A minimum spanning tree of cgMLST target differences. A cgMLST distance and a "
                   "seven-locus distance are different quantities and never share a scale.",
    "snp": "A SNP tree from read alignment. Planned for a later round; nothing runs on this "
           "tab yet.",
    "evidence": "Review identity, resistance and virulence evidence for isolates you choose. "
                "Genotype is not measured susceptibility.",
    "reports": "Turn reviewed evidence into a document you can share.",
    "schemes": "Install and update the MLST and cgMLST scheme libraries and the reference "
               "databases this workspace types against.",
    "settings": "Window size, text size, screen resolution and where your data lives.",
}

# (button text, name of the MainWindow method the button calls). Resolved with
# getattr at build time so a missing method degrades to a hidden button.
NEXT_STEP = {
    "overview": ("Import sequences…", "browse_files"),
    "isolates": ("Add samples…", "browse_files"),
    "reads": ("Check reads in Samples →", "goto_assembly_step"),
    "assembly": ("Assemble in Samples →", "goto_assembly_step"),
    "mlst": ("Type selected isolates…", "run_mlst_station"),
    "compare": ("Choose cohort…", "choose_comparison_cohort"),
    "cgmlst": ("Call cgMLST on selected…", "run_cgmlst_station"),
    "cgmlst_tree": ("Call cgMLST first →", "goto_cgmlst_station"),
    "snp": ("Compare allele differences instead →", "goto_mlst_tree"),
    "evidence": ("Choose isolates…", "choose_feature_cohort"),
    "reports": ("Choose report isolates…", "choose_report_cohort"),
    "schemes": ("Browse online / install updates…", "open_reference_manager"),
    "settings": ("Window size, text size and screen resolution", "open_display_settings"),
}

#: What Clear does on each tab, in the words the confirmation uses. "clears" is
#: always a view, a selection or a cohort; "keeps" is always the evidence. Clearing
#: a tab must never delete a sample, a result, a scheme or an imported report.
CLEAR_ACTIONS = {
    "overview": {
        "clears": "the highlight shared between tabs",
        "keeps": "every sample, result and imported report in this project",
    },
    "isolates": {
        "clears": "the search box, the row selection and the shared highlight",
        "keeps": "every sample, its files, its organism and its saved results",
    },
    "reads": {
        "clears": "nothing — no read-quality work runs on this tab yet",
        "keeps": "every sample and every file exactly as it is",
    },
    "assembly": {
        "clears": "nothing — no assembly runs on this tab yet",
        "keeps": "every assembly already made in Samples",
    },
    "mlst": {
        "clears": "this tab's selection and its summary line",
        "keeps": "every stored seven-locus result and the scheme library",
    },
    "compare": {
        "clears": "the drawn tree, the distance table and this tab's cohort",
        "keeps": "every profile, saved investigation and exported file",
    },
    "cgmlst": {
        "clears": "this tab's selection and its summary line",
        "keeps": "every stored cgMLST result and the cgMLST scheme library",
    },
    "cgmlst_tree": {
        "clears": "the drawn cgMLST tree and this tab's cohort",
        "keeps": "every cgMLST profile and saved investigation",
    },
    "snp": {
        "clears": "nothing — no SNP work runs on this tab yet",
        "keeps": "every sample and every file exactly as it is",
    },
    "evidence": {
        "clears": "this tab's cohort and its table selection",
        "keeps": "every imported HYDRA report and every stored evidence record",
    },
    "reports": {
        "clears": "this tab's report cohort",
        "keeps": "every report you have already saved, and every stored result",
    },
    "schemes": {
        "clears": "the library view, which is then read from disk again",
        "keeps": "every installed scheme and reference database",
    },
    "settings": {
        "clears": "this page's view, returning it to the first section",
        "keeps": "every preference you have saved, including window and text size",
    },
}

#: The pages this layout builds itself: the pipeline stations no controller owns
#: yet. Each one says what it is for, what it does not do, and where that work
#: happens today. `MainWindow.adopt_station(key, widget)` replaces the placeholder
#: with a real page without moving a tab or renumbering anything.
STATION_PAGES = {
    "reads": {
        "title": "Read quality",
        "subtitle": "Trimming and quality filtering of raw read pairs, before anything is "
                    "assembled or typed.",
        "body": [
            "A later round puts read trimming and filtering here, with the filtered read "
            "files kept beside the originals and the original files never modified.",
            "Until then, the read statistics this application does compute are shown per "
            "sample in Samples. They describe the reads as supplied: sampled read statistics "
            "are a prefix of the file, not quality validation of the whole file.",
        ],
        "actions": [("Open Samples", "goto_samples")],
    },
    "assembly": {
        "title": "Assembly",
        "subtitle": "Read pairs into assemblies, with the metrics that say whether an assembly "
                    "is worth typing.",
        "body": [
            "A tab of its own is planned for a later round: assembler choice, contig metrics "
            "and a per-isolate record of how each assembly was produced.",
            "Assembling already works today from Samples, where the Assembly step runs the "
            "bundled assembler on the read pairs you select. An assembly is a reconstruction, "
            "not a finished genome, and a typing result inherits its limits.",
        ],
        "actions": [("Open Samples → Assembly", "goto_assembly_step")],
    },
    "mlst": {
        "title": "Classical MLST",
        "subtitle": "Seven loci, one sequence type per isolate, against a curated allele panel.",
        "body": [
            "Runs on the isolates selected in Samples, so one selection carries through the "
            "whole pipeline. The organism established by the first analysis fills in here; "
            "typing is not repeated unless the scheme or the isolate's configuration changes.",
            "A seven-locus sequence type is a lineage label from a curated panel. It is not a "
            "species confirmation, and its allele differences are a different quantity from "
            "cgMLST target differences — the two never share a scale, a column or a threshold.",
        ],
        "actions": [("Type selected isolates…", "run_mlst_station"),
                    ("Open the MLST table in Samples →", "goto_mlst_step"),
                    ("MLST schemes…", "open_reference_manager")],
    },
    "cgmlst": {
        "title": "Core-genome MLST",
        "subtitle": "Hundreds to thousands of targets, against a cgMLST scheme chosen for the "
                    "organism.",
        "body": [
            "cgMLST has its own scheme library, its own target count and its own missing-target "
            "count. Missing, ambiguous, mixed and unknown targets stay distinct from assigned "
            "alleles wherever they are counted or drawn.",
            "The table of calls and the read-based recheck of missing targets belong on this "
            "tab. A recheck is investigative evidence about whether a target's sequence is "
            "present in the reads; it never becomes an allele call, a distance or a tree edge.",
        ],
        # A third entry in an action is the attribute that has to exist before the
        # button is offered: a button that leads nowhere is worse than no button.
        "actions": [("Call cgMLST on selected…", "run_cgmlst_station"),
                    ("Open the table of calls →", "goto_cgmlst_calls", "cgmlst_calls"),
                    ("Open the cgMLST table in Samples →", "goto_cgmlst_step"),
                    ("cgMLST schemes…", "open_reference_manager")],
    },
    "cgmlst_tree": {
        "title": "cgMLST tree",
        "subtitle": "A minimum spanning tree of cgMLST target differences, on its own scale.",
        "body": [
            "Nothing is drawn on this tab yet. Until it draws its own graph, the button "
            "above takes you to the tree page and switches it to cgMLST: that page states "
            "the scheme and the target count of whatever it has drawn, so a 2,000-target "
            "distance is never presented as a seven-locus one.",
            "A minimum spanning tree is a layout of similarity, not a phylogeny, and "
            "similarity is not proof of transmission. Distances carry their shared-target "
            "denominator: too little overlap is reported as insufficient evidence, never as "
            "zero distance.",
        ],
        "actions": [("Draw it on the tree page →", "goto_cgmlst_tree_view", "show_cgmlst_tree"),
                    ("Call cgMLST on selected…", "run_cgmlst_station")],
    },
    "snp": {
        "title": "SNP tree",
        "subtitle": "Single-nucleotide differences from read alignment, as a separate line of "
                    "evidence from allele typing.",
        "body": [
            "Nothing runs on this tab yet. A later round puts reference-free SNP distances "
            "here, with their own alignment fraction and their own scale.",
            "SNP distances and allele differences answer different questions and are never "
            "merged, summed or plotted on one axis.",
        ],
        "actions": [("Open the MLST tree →", "goto_mlst_tree")],
    },
}

#: The horizontal tab padding to try, widest first. Thirteen tabs have to stay
#: readable in a 1000 px window, and a chopped-off label ("cgMLST t…") is worse
#: than a tighter one, so the bar gives up padding before it gives up words.
PADDING_STEPS = (12, 10, 9, 8, 7, 6, 5)


def page_key(index: int) -> str:
    """The stable key of a build-order index, or '' when the index is unknown."""
    return PAGE_KEYS[index] if 0 <= index < len(PAGE_KEYS) else ""


def tab_title(key: str) -> str:
    """Tab bar text: the page's short label.

    Thirteen tabs have to fit a 1000 px window, so the label carries no glyph and
    no decoration; the long name and the purpose line live in the tooltip.
    """
    return TAB_LABELS.get(key, str(key).replace("_", " ").title())


def page_name(key: str) -> str:
    """The name the breadcrumb, the View menu and Alt+N use for a page.

    Deliberately the same word as its tab: a page called one thing in the bar and
    another in a menu is a page somebody looks for twice.
    """
    return tab_title(key)


def clear_promise(key: str) -> dict:
    """What Clear removes and what it keeps on one tab, for the confirmation."""
    return dict(CLEAR_ACTIONS.get(key, {"clears": "this tab's view",
                                        "keeps": "everything this project has stored"}))


class WorkspaceTabs(QTabWidget):
    """Pages addressed by key or by their build-order number, shown in pipeline order.

    Integers here are **page numbers** -- the order the pages were built -- not tab
    positions: `widget(3)`, `currentIndex()`, `index_of(3)` and
    `page_index()["schemes"]` all mean the fourth page built, wherever its tab now
    sits in the bar. Qt's own tab API (`tabText`, `tabToolTip`, `setTabVisible`,
    `tabBar()`) still counts tab positions; ask `position_of(key)` for one.
    """

    pageShown = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("workspaceTabs")
        self.setDocumentMode(False)
        self.setMovable(False)
        self.setUsesScrollButtons(True)
        self.tabBar().setObjectName("workspaceTabBar")
        self.tabBar().setExpanding(False)
        # Never elide: a tab whose words are chopped cannot be read, and this bar
        # is the navigation. If the labels cannot all fit, the bar scrolls instead.
        self.tabBar().setElideMode(Qt.TextElideMode.ElideNone)
        self.tabBar().setAccessibleName("Workspace sections")
        self.setAccessibleName("Workspace sections")
        # Native child viewports must paint directly; a transparent stack caused
        # stale or blank pages on Windows.
        self.setAutoFillBackground(True)
        self._slots: list[str] = []
        self._hooks: dict[str, object] = {}
        self._subtabs: dict[str, QTabWidget] = {}
        self._folded: dict[str, tuple[str, str]] = {}
        self._depth = 0
        self._padding = None
        self._minimum_bar = None
        self._fitting = False
        self.currentChanged.connect(self._current_changed)

    # --- building -----------------------------------------------------------
    def add_page(self, widget, key=None, *, title=None, tooltip=None, purpose=None) -> int:
        """Add a page and return its build-order slot, not its tab position.

        The tab itself is inserted where PIPELINE says it belongs, so the bar
        reads as the workflow while the returned number keeps meaning "the Nth
        page built" for every existing call site.
        """
        if key is None:
            key = page_key(len(self._slots)) or f"page{len(self._slots)}"
        key = str(key)
        position = self._pipeline_position(key)
        self.insertTab(position, widget, title if title is not None else tab_title(key))
        self.tabBar().setTabData(position, key)
        slot = len(self._slots)
        self._slots.append(key)
        self._minimum_bar = None
        if tooltip is None:
            tooltip = purpose if purpose is not None else PAGE_PURPOSE.get(key, "")
        self.setTabToolTip(position, str(tooltip))
        self.setTabWhatsThis(position, str(tooltip))
        self._fit_tab_bar()
        return slot

    def addWidget(self, widget) -> int:
        """QStackedWidget compatibility: keep `pages.addWidget(page)` call sites working."""
        return self.add_page(widget)

    def _pipeline_position(self, key) -> int:
        """Where this key's tab belongs in the bar, given the tabs already there."""
        rank = PIPELINE.index(key) if key in PIPELINE else len(PIPELINE)
        for position in range(self.count()):
            other = self.tabBar().tabData(position)
            other_rank = PIPELINE.index(other) if other in PIPELINE else len(PIPELINE)
            if other_rank > rank:
                return position
        return self.count()

    # --- addressing ---------------------------------------------------------
    def keys(self) -> tuple[str, ...]:
        """Every page key in build order: the integer address space."""
        return tuple(self._slots)

    def tab_order(self) -> tuple[str, ...]:
        """Every page key in the order its tab is shown."""
        return tuple(str(self.tabBar().tabData(position) or "")
                     for position in range(self.count()))

    def key_at(self, slot) -> str:
        """The key of a build-order slot."""
        return self._slots[slot] if isinstance(slot, int) and not isinstance(slot, bool) \
            and 0 <= slot < len(self._slots) else ""

    def key_at_position(self, position) -> str:
        """The key of a tab position in the bar."""
        if not isinstance(position, int) or isinstance(position, bool):
            return ""
        return str(self.tabBar().tabData(position) or "") if 0 <= position < self.count() else ""

    def key_for(self, target) -> str:
        """Resolve a key or a build-order number to a key, or '' when unknown."""
        if isinstance(target, bool):
            return ""
        if isinstance(target, int):
            return self.key_at(target)
        return str(target) if str(target) in self._slots else ""

    def slot_of(self, target) -> int:
        """The build-order number of a page, or -1 when it does not exist."""
        key = self.key_for(target)
        return self._slots.index(key) if key else -1

    def position_of(self, target) -> int:
        """The tab position of a page, or -1 when it does not exist."""
        key = self.key_for(target)
        if not key:
            return -1
        for position in range(self.count()):
            if self.tabBar().tabData(position) == key:
                return position
        return -1

    def index_of(self, target) -> int:
        """Resolve a key or a number to the build-order number every call site uses."""
        return self.slot_of(target)

    def page_index(self) -> dict[str, int]:
        """{key: build-order number} for the sites that still need a number."""
        return {key: slot for slot, key in enumerate(self._slots)}

    def current_key(self) -> str:
        return self.key_at_position(self.tabBar().currentIndex())

    # Page numbers, not tab positions. QTabWidget numbers its pages by where their
    # tab sits; this window has always numbered them by the order they were built,
    # and `pages.widget(3)` means "the fourth page built" in code that predates the
    # pipeline layout. These three keep every such call site pointing at the page it
    # named, whichever tab now holds it. Qt's own tab API (tabText, tabToolTip,
    # setTabVisible, the tab bar itself) still takes a position: use position_of.
    def widget(self, index):
        """The page with this build-order number."""
        if isinstance(index, int) and not isinstance(index, bool) \
                and 0 <= index < len(self._slots):
            return super().widget(self.position_of(self._slots[index]))
        return super().widget(index)

    def currentIndex(self) -> int:  # noqa: N802 - Qt's own name
        """The build-order number of the page being shown."""
        return self.slot_of(self.current_key())

    def setCurrentIndex(self, index) -> None:  # noqa: N802 - Qt's own name
        """Show the page with this build-order number."""
        position = self.position_of(index) if self.key_for(index) else int(index)
        super().setCurrentIndex(position)

    def page_for(self, target) -> QWidget | None:
        position = self.position_of(target)
        return super().widget(position) if position >= 0 else None

    def show_page(self, target) -> bool:
        """Show a page by key or build-order number. False when it does not exist."""
        key = self.key_for(target)
        if not key:
            return False
        host = self._folded.get(key)
        if host is not None:
            # This page now lives as a sub-tab of another one. Every call site that
            # asks for it by key or by its old number still arrives at the content.
            host_position = self.position_of(host[0])
            if host_position >= 0:
                QTabWidget.setCurrentIndex(self, host_position)
                self.show_subtab(host[0], host[1])
                return True
        position = self.position_of(key)
        if position < 0:
            return False
        QTabWidget.setCurrentIndex(self, position)
        return True

    # --- the tab bar has to stay readable -----------------------------------
    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._fit_tab_bar()

    def changeEvent(self, event):
        super().changeEvent(event)
        # A larger interface scale re-polishes every widget with new text sizes,
        # so the bar has to be measured again or the words run off the window.
        self._minimum_bar = None
        self._fit_tab_bar()

    def minimum_bar_width(self) -> int:
        """The narrowest the bar can be drawn with every label still whole.

        The window asks before deciding whether it can afford the sidebar: losing
        the sidebar is a smaller loss than pushing tabs behind a scroll arrow.
        """
        if self._minimum_bar is None:
            if not self.count():
                return 0
            current = self._padding
            self._set_padding(PADDING_STEPS[-1])
            self._minimum_bar = self.tabBar().sizeHint().width()
            if current is not None:
                self._set_padding(current)
        return self._minimum_bar

    def _fit_tab_bar(self, available=None) -> int:
        """Trade tab padding for room until every label fits, widest padding first.

        Returns the padding in use. The words themselves are never shortened: a
        reader who cannot tell "cgMLST tree" from "cgMLST" has lost the navigation.
        """
        if self._fitting or not self.count():
            return self._padding or PADDING_STEPS[0]
        room = self.width() if available is None else int(available)
        if room <= 0:
            return self._padding or PADDING_STEPS[0]
        self._fitting = True
        try:
            for padding in PADDING_STEPS:
                self._set_padding(padding)
                if self.tabBar().sizeHint().width() <= room:
                    return padding
            return PADDING_STEPS[-1]
        finally:
            self._fitting = False

    def _set_padding(self, padding) -> None:
        """Only the sides: the height of a tab belongs to the interface scale."""
        if padding == self._padding:
            return
        self._padding = padding
        self.tabBar().setStyleSheet(
            "QTabBar#workspaceTabBar::tab { padding-left: %dpx; padding-right: %dpx; }"
            % (padding, padding))

    def mark_planned(self, key, planned=True) -> bool:
        """Say on the tab that a station is a place a later round fills.

        The tab stays reachable and keeps its full name -- its page explains what
        will live there and where that work happens today -- but it must not read as
        finished work. The message goes in the tooltip rather than in a colour: the
        stylesheet gives every tab its colour, so a per-tab one would not be painted.
        """
        position = self.position_of(key)
        if position < 0:
            return False
        tip = self.tabToolTip(position)
        suffix = " (planned — nothing runs here yet)"
        if planned and not tip.endswith(suffix):
            self.setTabToolTip(position, tip + suffix)
        elif not planned and tip.endswith(suffix):
            self.setTabToolTip(position, tip[:-len(suffix)])
        return True

    # --- folding a page into another page's sub-tabs ------------------------
    def fold_page(self, key, host_key=None, subtab_title=None) -> bool:
        """Stop showing a top-level tab whose content is now a sub-tab elsewhere.

        The page itself is kept, so `widget(slot)`, `page_index()` and every
        `navigate(<int>)` call site keep meaning what they meant; only the tab stops
        being offered twice. Returns False when the page or its host is unknown.
        """
        if host_key is None or subtab_title is None:
            host_key, subtab_title = FOLDABLE.get(str(key), (None, None))
        position, host_position = self.position_of(key), self.position_of(host_key)
        if position < 0 or host_position < 0 or position == host_position:
            return False
        self._folded[str(key)] = (str(host_key), str(subtab_title))
        self._minimum_bar = None
        self.setTabVisible(position, False)
        if self.tabBar().currentIndex() == position:
            self.show_page(key)
        self._fit_tab_bar()
        return True

    def unfold_page(self, key) -> bool:
        """Offer a folded page as its own tab again."""
        position = self.position_of(key)
        if position < 0 or self._folded.pop(str(key), None) is None:
            return False
        self._minimum_bar = None
        self.setTabVisible(position, True)
        self._fit_tab_bar()
        return True

    def folded_pages(self) -> dict[str, tuple[str, str]]:
        return dict(self._folded)

    # --- per-page show hook -------------------------------------------------
    def set_show_hook(self, key, callback) -> None:
        """Run `callback` whenever that page becomes visible; None removes the hook."""
        if callback is None:
            self._hooks.pop(str(key), None)
        else:
            self._hooks[str(key)] = callback

    def show_hook(self, key):
        return self._hooks.get(str(key))

    def _current_changed(self, position):
        # A hook may legitimately redirect to another page, so nesting is allowed;
        # the depth limit only stops two hooks from bouncing off each other forever.
        key = self.key_at_position(position)
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
