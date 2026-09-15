"""Sample-centric native workspace: library, assignments, jobs and application menus."""

import html
import json
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, QUrl
from PySide6.QtGui import QAction, QColor, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLineEdit,
    QSplitter,
    QTabWidget,
    QTextBrowser,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.archive import active_samples, archived_samples, is_archived
from wmlstudio.background import FunctionWorker
from wmlstudio.jobs import AnalysisWorker
from wmlstudio.storage import QUARANTINE_BUCKETS
from wmlstudio.ui_common import (
    FlowLayout,
    cell,
    flattened_metadata,
    gene_names,
    make_table,
    organism_for,
)
from wmlstudio.widgets import DropZone, Helix, Metric, button, card, label


def quarantine_bucket(sample):
    """The needs-review folder this isolate sits in, or None when it is filed.

    Quarantine means WMLSTudio declined to decide what the organism is. It is not
    a claim that the organism is unusual, and it is not the user's own "unknown".

    storage.review_bucket is the one place that answers this, so the count shown
    here and the folder the file actually sits in cannot drift apart: a proposal
    nobody accepted also waits under _Unresolved, and it is not "filed".
    """
    from wmlstudio.storage import review_bucket
    token = review_bucket((sample.get("metadata") or {}))
    return QUARANTINE_BUCKETS.get(token) if token else None


# Two kinds of trouble, two colours — and the words are always in the cell too.
# Colour alone is unreadable to a colour-blind reader and to a printed screenshot,
# so every amber row also says "Needs organism review" in text.
REVIEW_STYLES = {"conflict": ("#E08A2E", "#20160A"), "undecided": ("#F2C94C", "#241D07")}
REVIEW_PREFIX = "Needs organism review"
STEP_KEYS = ("assembly", "st", "cgmlst", "hydra")
#: Task roles that put reference data on this computer rather than computing
#: anything from a sample. When one of these stops, the Update page re-reads the
#: disk: a row that still says "Not installed" after a download that worked is
#: indistinguishable from a download that never ran.
INSTALL_ROLES = ("species_panel", "characterization_references", "scheme_import")
# The heading each step is known by. "ST" is the classical seven-locus scheme and
# cgMLST is a different quantity against a different scheme; they never share a
# column, a tab or a threshold.
STEP_TITLES = {"assembly": "Assembly", "st": "ST", "cgmlst": "cgMLST", "hydra": "HYDRA"}
STEP_PURPOSE = {
    "assembly": "Turn read pairs into assemblies. Read quality is not typing, and an "
                "assembly is not a finished genome.",
    "st": "Classical seven-locus MLST: one ST per isolate, from a curated allele panel.",
    "cgmlst": "Core-genome MLST: hundreds to thousands of targets. A cgMLST distance and a "
              "seven-locus distance are different quantities and never share a scale.",
    "hydra": "Screen assemblies for resistance, virulence and plasmid markers. A genotype is "
             "not measured susceptibility.",
}
STEP_ACTIONS = {"assembly": "Assemble selected read pairs…", "st": "Type selected · 7-locus MLST…",
                "cgmlst": "Call cgMLST on selected…", "hydra": "Run HYDRA on selected…"}
# The two typing steps are the ones that can already be answered by a stored
# result, so they are the two that offer to run again anyway.
RERUN_STEPS = ("st", "cgmlst")
# What a picker for the classical seven-locus workflow may offer. A scheme whose
# kind could not be read is still offered there — hiding a scheme somebody
# installed helps nobody — but a cgMLST target set never is: the workflow scheme
# is what an ST is called against, and an ST is not a core-genome profile.
CLASSICAL_KINDS = ("mlst", "unknown")
RERUN_TOGGLE = "Re-run even when nothing changed"
STANDING_NOTE = ("Isolates whose stored result nothing has invalidated are left alone, and this "
                 "tab names what changed for the ones that do run.")
RERUN_TOGGLE_TIP = (
    "Off, a stored result that nothing has invalidated is left alone and this tab says so "
    "instead of repeating the work. A result counts as standing only when the sequence input "
    "and the scheme's own fingerprint both still match what it was produced from; a result "
    "that never recorded them is re-run, not assumed current.")

# What a Clear button does and, just as importantly, what it does not. Clearing is
# about what you are working on at this moment; it never removes an isolate, a
# result, an allele call, a file or a history entry.
CLEAR_KEEPS = ("Nothing was removed — every isolate, result, allele call, file and history entry "
               "is untouched, so the next thing you do can use new, past or a mix of samples.")
CLEAR_NOTE = ("Cleared: selection, search, organism filters and the library branch. " + CLEAR_KEEPS)

# Dropping files is the fast route the user asked to have back. It is fast because
# it understands what was dropped and where it landed, not because it skips
# anything: dropped sequences still go through identification-before-copy and the
# same review dialog the menus use.
SEQUENCE_SUFFIXES = (".fa", ".fasta", ".fna", ".fq", ".fastq")
READ_SUFFIXES = (".fq", ".fastq")
TABLE_SUFFIXES = (".csv", ".tsv", ".tab")
SCHEME_MARKERS = ("profiles.tsv", "profiles.txt", "profiles.tab", "scheme.json")
DROP_HINT = ("Drop sequences to import them, a scheme folder on the Schemes tab to install it, "
             "a CSV/TSV to review epidemiology, or read files on one assembly row to attach them.")


def plain_suffix(path):
    """A file's own extension, with any compression wrapper taken off."""
    name = Path(path).name.casefold().removesuffix(".gz").removesuffix(".bz2")
    return Path(name).suffix


def looks_like_scheme(path):
    """Whether a dropped folder is a typing scheme rather than a folder of isolates.

    Only a scheme's own unambiguous markers count — a profile table, or allele files
    under the .tfa extension curated schemes use. A folder of .fasta files is a
    cohort of assemblies everywhere except the scheme library, where a folder can
    mean nothing else and the caller says so.
    """
    folder = Path(path)
    if not folder.is_dir():
        return False
    if any((folder / marker).is_file() for marker in SCHEME_MARKERS):
        return True
    try:
        return any(entry.is_file() and plain_suffix(entry) == ".tfa" for entry in folder.iterdir())
    except OSError:
        return False  # An unreadable folder is not claimed to be anything.


def sort_dropped(paths, *, schemes_expected=False):
    """Split a drop into the kinds of thing this workspace knows what to do with.

    Nothing here touches a file: it reads names and asks whether a folder carries a
    scheme's markers, so a drop can be routed before anything is opened or copied.
    """
    sorted_paths = {"sequences": [], "schemes": [], "tables": [], "unusable": []}
    for value in paths:
        path = Path(value)
        if path.is_dir():
            kind = "schemes" if (schemes_expected or looks_like_scheme(path)) else "sequences"
        elif plain_suffix(path) in SEQUENCE_SUFFIXES:
            kind = "sequences"
        elif plain_suffix(path) in TABLE_SUFFIXES:
            kind = "tables"
        else:
            kind = "unusable"
        sorted_paths[kind].append(str(path))
    return sorted_paths


def input_kind(sample):
    """What this isolate actually is on disk: reads, an assembly, or a profile only."""
    metadata = sample.get("metadata") or {}
    workflow = metadata.get("workflow") or {}
    if workflow.get("source_kind") == "profile" or sample.get("profile_only") or not sample.get("input_path"):
        return "profile"
    if workflow.get("source_kind") == "read_mate":
        return "read_mate"
    if metadata.get("assembly"):
        return "assembly"  # Assembled here from its own reads.
    name = Path(sample["input_path"]).name.casefold().removesuffix(".gz").removesuffix(".bz2")
    return "reads" if name.endswith((".fq", ".fastq")) else "assembly"


def accepted_organism(sample):
    """The organism somebody or something actually decided on, or None.

    A proposal nobody accepted is not a decision, so the pipeline may not act on
    it. This is the fact the ST, cgMLST and HYDRA gates ask for.
    """
    if quarantine_bucket(sample):
        return None
    genus, species, _ = organism_for(sample)
    return (genus, species) if genus else None


def organism_text(organism):
    return " ".join(part for part in (organism or ()) if part) or "No organism set"


def proposed_organism(sample):
    """The organism the evidence put forward but nobody has accepted, or None.

    Worth showing — it is often the classical panel's suggestion, which is what a
    user is asking for when they drop files in — but never worth acting on: it is
    always shown beside the words "not accepted".
    """
    evidence = (sample.get("metadata") or {}).get("organism_evidence") or {}
    proposed = evidence.get("proposed") or {}
    if not proposed.get("genus"):
        return None
    return (str(proposed["genus"]), str(proposed.get("species") or ""))


def organism_cell_text(sample):
    """What to write in an Organism column: the decision, or the unaccepted proposal."""
    accepted = accepted_organism(sample)
    if accepted:
        return organism_text(accepted)
    suggestion = proposed_organism(sample)
    return f"Proposed: {organism_text(suggestion)} · not accepted" if suggestion else "No organism set"


def scheme_organism_index(entries):
    """{scheme name or folder id → its entry} for the installed references."""
    index = {}
    for entry in entries or ():
        for key in (entry.get("name"), entry.get("id")):
            if key:
                index.setdefault(str(key).casefold(), entry)
    return index


def mlst_corroboration(sample, schemes=None, mlst=None):
    """The organism label carried by the classical scheme this isolate matched.

    The user asked for MLST to count as identification evidence, and it does —
    as corroboration. A seven-locus profile is compatibility with a curated panel,
    so the genus belongs to the reference that panel was built from, not to this
    genome. It is reported beside the genomic evidence and never in place of it.
    """
    if mlst is None:
        from wmlstudio.sample_workflow import typing_profiles
        mlst = typing_profiles(sample)["mlst"] or {}
    name = str(mlst.get("scheme") or "")
    entry = (schemes or {}).get(name.casefold()) if name else None
    if entry and entry.get("genus"):
        return {"genus": str(entry["genus"]), "species": str(entry.get("species") or ""),
                "scheme": str(entry.get("name") or name), "st": mlst.get("st"),
                "source": "the seven-locus scheme this isolate typed against"}
    evidence = (sample.get("metadata") or {}).get("organism_evidence") or {}
    for candidate in (evidence.get("detail") or {}).get("mlst_top") or ():
        parts = str(candidate.get("organism_label") or "").split()
        if parts:
            return {"genus": parts[0], "species": parts[1] if len(parts) > 1 else "",
                    "scheme": str(candidate.get("scheme") or ""), "st": candidate.get("st"),
                    "source": "the typing panel that matched during identification"}
    return None


def identification_support(sample, schemes=None, mlst=None):
    """What this isolate's organism rests on, and whether the two sources agree.

    There are only two sources here and they are not the same kind of thing: the
    claim on the record — a genome comparison, or a person's own assignment — and
    the classical typing panel, which can corroborate a claim but never make one.
    They are compared at genus level only: a curated panel routinely covers a whole
    species complex, so a species difference inside one genus is not a conflict.
    """
    from wmlstudio.workflow_dialogs import BASIS_LABELS
    evidence = (sample.get("metadata") or {}).get("organism_evidence") or {}
    basis = str(evidence.get("basis") or "")
    label = evidence.get("accepted") or {}
    if not label.get("genus"):
        label = evidence.get("proposed") or {}
    claim = None
    if label.get("genus") and basis not in {"", "none", "mlst_panel"}:
        claim = {"genus": str(label["genus"]), "species": str(label.get("species") or ""),
                 "source": BASIS_LABELS.get(basis, basis),
                 "genomic": basis.startswith("genomic_")}
    panel = mlst_corroboration(sample, schemes, mlst)
    if claim and panel:
        agreement = "agree" if claim["genus"].casefold() == panel["genus"].casefold() else "conflict"
    else:
        agreement = "claim_only" if claim else "panel_only" if panel else "none"
    sentences = {
        "agree": lambda: (f"{claim['source']} and {panel['source']} both point to "
                          f"{panel['genus']}. The panel corroborates the call; a panel match is "
                          "not an independent identification."),
        "conflict": lambda: (f"{claim['source']}: {claim['genus']} {claim['species']}".rstrip()
                             + f". {panel['source']}: labelled {panel['genus']} "
                               f"{panel['species']}".rstrip()
                             + ". These disagree, so nothing is assumed: review the organism "
                               "before analysing this isolate."),
        "claim_only": lambda: (f"{claim['source']} supported this organism. "
                               + ("No classical typing panel has corroborated it yet."
                                  if claim["genomic"] else
                                  "No installed reference has corroborated it.")),
        "panel_only": lambda: (f"Only {panel['source']} supports this organism. Panel "
                               "compatibility is not a species identification."),
        "none": lambda: "Nothing installed has supported an organism for this isolate.",
    }
    return {"claim": claim, "genomic": claim if claim and claim["genomic"] else None,
            "panel": panel, "agreement": agreement, "sentence": sentences[agreement]()}


def organism_review(sample, schemes=None, mlst=None):
    """Whether this isolate needs a person to look at its organism, and why.

    Returned as words first: `label` is written into the row and `level` only
    chooses which amber the row is painted.
    """
    evidence = (sample.get("metadata") or {}).get("organism_evidence") or {}
    support = identification_support(sample, schemes, mlst)
    bucket = quarantine_bucket(sample)
    confidence = str(evidence.get("confidence") or "")
    settled = {"needs_review": False, "level": "", "label": "", "detail": support["sentence"],
               "support": support}
    if support["agreement"] == "conflict":
        return {**settled, "needs_review": True, "level": "conflict",
                "label": f"{REVIEW_PREFIX} · genome comparison and typing panel disagree"}
    suggestion = proposed_organism(sample)
    # The proposal is shown because it is often what the user wants to accept, and
    # it always carries "not accepted" with it: it decided nothing.
    offer = f" · proposed {organism_text(suggestion)}, not accepted" if suggestion else ""
    if bucket:
        words = bucket.replace("_", " ").casefold()
        level = "conflict" if "conflicting" in words else "undecided"
        return {**settled, "needs_review": True, "level": level,
                "label": f"{REVIEW_PREFIX} · {words}{offer}",
                "detail": organism_evidence_note(sample) or support["sentence"]}
    if confidence == "complex_only":
        return {**settled, "needs_review": True, "level": "conflict",
                "label": f"{REVIEW_PREFIX} · species complex, not separated",
                "detail": organism_evidence_note(sample) or support["sentence"]}
    if accepted_organism(sample) is None:
        return {**settled, "needs_review": True, "level": "undecided",
                "label": f"{REVIEW_PREFIX} · no organism set{offer}"}
    return settled


def paint_review(item, review):
    """Repeat a review state as colour on a cell whose text already says it."""
    if not review.get("needs_review"):
        return item
    background, foreground = REVIEW_STYLES.get(review["level"], REVIEW_STYLES["undecided"])
    item.setBackground(QColor(background))
    item.setForeground(QColor(foreground))
    item.setToolTip(f"{review['label']}\n\n{review['detail']}")
    return item


def organism_evidence_note(sample):
    """Where this genus and species came from, in the words the engine used."""
    evidence = (sample.get("metadata") or {}).get("organism_evidence")
    if not isinstance(evidence, dict):
        return ""
    from wmlstudio.organism_id import confidence_label, evidence_sentences
    from wmlstudio.workflow_dialogs import BASIS_LABELS
    word, explanation = confidence_label(evidence)
    basis = BASIS_LABELS.get(evidence.get("basis"), str(evidence.get("basis") or "Not identified"))
    lines = [f"{basis} · {word}", explanation,
             "A folder name is where the copy is stored; it is not a laboratory identification."]
    lines.extend(evidence_sentences(evidence))
    return "\n".join(line for line in lines if line)


def apply_application_style(style):
    """Set the application stylesheet only when it would actually change.

    Qt re-polishes every widget of every open window on each assignment, so
    reassigning the same sheet per window costs more the more windows are open.
    Opening projects one after another slowed down for the same reason.
    """
    application = QApplication.instance()
    if application is None or application.styleSheet() == style:
        return False
    application.setStyleSheet(style)
    return True


class WorkbenchMixin:
    def __init__(self, *args, **kwargs):
        self.selection_ids = set()
        self.cohort_ids = set()
        self.feature_ids = set()
        self.report_ids = set()
        self.library_filter = None
        self.worker_role = ""
        self._filling = False
        self._refreshing = False
        self.reference_dialog = None
        self.amr_database_dialog = None
        self.library = None
        self._pending_graph_state = None
        self._run_plan = {}
        self._run_ids = set()
        self._run_cancelled = False
        self._task_succeeded = False
        self._typing_override = None
        self._progress_dialog = None
        # What the progress dialog should call the running task. Set per launch;
        # None means the analysis wording, which is what most callers are.
        self._task_caption = None
        self._pending_import = None
        self._import_notes = []
        self._practice_cohort = None
        self._intake_report = None
        self._scheme_rows = []
        self._scheme_rows_key = None
        self._amr_state = None
        self._amr_state_key = None
        self._visible_samples = []
        self._typing_overview = {}
        self._called_cache = {}
        self.step_tables = {}
        self.step_buttons = {}
        self.step_notes = {}
        self.step_toggles = {}
        self._drop_step = ""
        self._pending_typing = None
        self._typing_currency = None
        self._installed_scheme = None
        # What a stored result already answers, per (step, isolate). Kept after a
        # check so the line under the button keeps saying which isolates it left
        # alone, and dropped for an isolate the moment it is typed again.
        self._typing_verdicts = {}
        self.last_typing_report = ""
        self.ui_scale = 100
        super().__init__(*args, **kwargs)
        from wmlstudio import theme
        application = QApplication.instance()
        if hasattr(theme, "apply_dark_palette"):
            theme.apply_dark_palette(application)
        # Setting the application stylesheet re-polishes every widget of every
        # open window, so doing it unconditionally per window costs more the more
        # windows exist. Apply it only when it would actually change.
        apply_application_style(theme.STYLE)
        from wmlstudio.interface_settings import interface_preferences
        preferences = interface_preferences(self.root)
        try:
            initial_scale = int(preferences.value("scale", 100))
        except (TypeError, ValueError):
            initial_scale = 100
        self.set_ui_scale(initial_scale)
        for combo in self.findChildren(QComboBox):
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(int(combo.property("compactCharacters") or 18))
            combo.setMaximumWidth(360)
        self.setWindowTitle("WMLSTudio · Genomics workbench")
        self.build_menus()
        from wmlstudio.library import Library
        self.library = Library(self.root / "library-index.sqlite")
        self.library.index_project(self.project)
        saved_cohort = self.project.get_setting("comparison_cohort", None)
        if saved_cohort is not None:
            self.cohort_ids = set(saved_cohort)
        if hasattr(self.tree, "restore_state"):
            self._pending_graph_state = self.project.get_setting("graph_style", {"version": 1})
        self.refresh()
        self.statusBar().showMessage("Local genomics workbench · Import isolates, then choose an analysis")

    def build_overview(self):
        _, layout = self.page()
        hero, content = card()
        row = QHBoxLayout()
        words = QVBoxLayout()
        words.addWidget(label("FROM A QUESTION TO REVIEWABLE EVIDENCE", "eyebrow"))
        words.addWidget(label("Your investigation, connected.", "title", True))
        words.addWidget(label("Identify, characterize, compare and report the same isolates. Add new samples without losing earlier work.", "muted", True))
        actions = FlowLayout()
        actions.addWidget(button("Import sequences…", self.browse_files, True))
        actions.addWidget(button("New project…", self.new_project))
        actions.addWidget(button("Research saved library…", self.open_library_research))
        actions.addStretch()
        words.addLayout(actions)
        row.addLayout(words, 4)
        self.helix = Helix()
        self.helix.setMinimumSize(150, 130)
        self.helix.setMaximumWidth(220)
        self.helix.setMaximumHeight(160)
        row.addWidget(self.helix, 1)
        content.addLayout(row)
        layout.addWidget(hero)
        from wmlstudio.journey_widgets import InvestigationMap

        self.overview_tabs = QTabWidget()
        self.investigation_map = InvestigationMap()
        self.investigation_map.actionRequested.connect(self.journey_action)
        self.overview_tabs.addTab(self.investigation_map, "Investigation map")
        library_page = QWidget()
        library_layout = QVBoxLayout(library_page)
        library_layout.setContentsMargins(0, 14, 0, 0)
        self.overview_tabs.addTab(library_page, "Stored library")
        layout.addWidget(self.overview_tabs, 1)
        metrics = QHBoxLayout()
        self.metrics = [Metric("Isolates", "in the active library"), Metric("Analysed", "saved evidence"),
                        Metric("Exact MLST", "registered profiles"), Metric("Review", "incomplete or failed")]
        for metric in self.metrics:
            metrics.addWidget(metric)
        library_layout.addLayout(metrics)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        explorer, content = card()
        content.addWidget(label("Library navigator", "cardTitle"))
        self.library_tree = QTreeWidget()
        self.library_tree.setHeaderHidden(True)
        self.library_tree.setMinimumWidth(250)
        self.library_tree.itemActivated.connect(self.open_library_group)
        self.library_tree.itemClicked.connect(self.open_library_group)
        self.install_view_menu("library.tree", self.library_tree)
        content.addWidget(self.library_tree)
        splitter.addWidget(explorer)
        recent, content = card()
        content.addWidget(label("Recent samples · double-click to inspect", "cardTitle"))
        self.recent_table = make_table(["Sample", "Input", "Status", "ST", "Loci"])
        self.recent_table.cellDoubleClicked.connect(self.open_recent_sample)
        self.install_view_menu("overview.recent", self.recent_table)
        content.addWidget(self.recent_table)
        self.drop_zone = DropZone()
        # Routed rather than sent straight to intake: a scheme folder or an
        # epidemiology table dropped here means what it means anywhere else.
        self.drop_zone.filesDropped.connect(self.handle_drop)
        self.drop_zone.browseRequested.connect(self.browse_files)
        content.addWidget(self.drop_zone)
        splitter.addWidget(recent)
        splitter.setSizes([320, 650])
        library_layout.addWidget(splitter, 1)
        self.practice_notice = label("Practice project · synthetic sequences, not biological isolates.", "small")
        self.practice_notice.hide()
        layout.addWidget(self.practice_notice)

    def build_samples(self):
        """Samples is the hub: loaded once here, then seen from four angles.

        Assembly, ST, cgMLST and HYDRA are the same isolates, each tab showing the
        facts that step needs and running only that step. ST is the classical
        seven-locus scheme; cgMLST has its own scheme, its own target count and its
        own called/missing loci, and the two are never merged into one column.
        """
        _, layout = self.page()
        self.heading(layout, "Samples",
                     "Load your sequences once. Each tab below is the same isolates seen from one "
                     "angle and runs that one step for the isolates you select.")
        row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search names, organisms, STs, AMR genes, collections or metadata…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.refresh_tables)
        row.addWidget(self.search, 1)
        row.addWidget(button("Add samples…", self.browse_files, True))
        row.addWidget(button("Assign organism…", self.assign_selected))
        self.samples_advanced_button = button("Advanced ▾", None)
        self.samples_advanced_button.setCheckable(True)
        self.samples_advanced_button.setToolTip(
            "Filters, storage options and the older import dialog. Nothing here is needed for an "
            "ordinary investigation.")
        row.addWidget(self.samples_advanced_button)
        layout.addLayout(row)
        self.selection_label = label("No samples selected · each sample keeps its own organism and typing workflow", "small", True)
        layout.addWidget(self.selection_label)
        layout.addWidget(self.build_samples_advanced())
        self.sample_tabs = QTabWidget()
        self.sample_tabs.setAccessibleName("Sample steps")
        for key in STEP_KEYS:
            self.sample_tabs.addTab(self.build_step_page(key), STEP_TITLES[key])
            self.sample_tabs.setTabToolTip(self.sample_tabs.count() - 1, STEP_PURPOSE[key])
        self.sample_tabs.currentChanged.connect(lambda _index: self.update_step_gates())
        if hasattr(self.pages, "register_subtabs"):
            self.pages.register_subtabs("isolates", self.sample_tabs)
        layout.addWidget(self.sample_tabs, 1)
        self.apply_column_visibility()
        self.detail = QTextBrowser()
        self.detail.setOpenExternalLinks(False)
        self.detail.setMinimumHeight(120)
        inspector = QTabWidget()
        inspector.addTab(self.detail, "Sample evidence")
        self.history_view = QTextBrowser()
        inspector.addTab(self.history_view, "History / provenance")
        inspector.setParent(self)
        inspector.hide()  # Complete evidence is available in the isolate popup.

    def build_samples_advanced(self):
        """Rarely-used controls, hidden until asked for.

        The default view shows what somebody importing their first isolates needs.
        Filters, the managed-storage options and the older per-file import dialog
        live here, because reaching for them is the exception.
        """
        panel, content = card()
        filters = FlowLayout()
        self.genus_filter = QComboBox()
        self.species_filter = QComboBox()
        self.collection_filter = QComboBox()
        for combo, title in [(self.genus_filter, "All genera"), (self.species_filter, "All species"),
                             (self.collection_filter, "All collections")]:
            combo.addItem(title, "")
            combo.currentIndexChanged.connect(self.refresh_tables)
            filters.addWidget(combo, 1)
        filters.addWidget(button("Clear filters", self.clear_filters))
        content.addLayout(filters)
        strip = FlowLayout()
        self.scheme_combo = QComboBox()
        self.scheme_combo.setMinimumWidth(180)
        self.scheme_combo.setAccessibleName("Typing scheme override")
        self.scheme_combo.hide()  # Per-isolate schemes are reviewed in the analysis dialog.
        self.run_button = button("Analyse…", self.choose_and_analyse, True)
        self.run_button.setToolTip("Choose any isolates and any analysis, in one review dialog.")
        self.rerun_button = button("Review pending…", lambda: self.choose_and_analyse(pending_only=True))
        for widget in (self.run_button, self.rerun_button,
                       button("Import with options…", self.browse_files_with_options),
                       button("Attach reads…", self.attach_reads_selected),
                       button("Re-file managed copies…", self.refile_selected),
                       button("Epidemiology grid…", self.open_metadata_grid)):
            strip.addWidget(widget)
        content.addLayout(strip)
        self.sample_columns_button = button("Show every column", None)
        self.sample_columns_button.setCheckable(True)
        self.sample_columns_button.toggled.connect(self.apply_column_visibility)
        content.addWidget(self.sample_columns_button)
        content.addWidget(label("Import with options… opens the per-file dialog: it asks for an "
                                "organism and a scheme before anything is copied. Dropping files "
                                "on the window does the same identification without the dialog. "
                                + DROP_HINT, "small", True))
        panel.setVisible(False)
        self.samples_advanced = panel
        self.samples_advanced_button.toggled.connect(panel.setVisible)
        return panel

    def build_step_page(self, key):
        """One sub-tab: the isolates from this step's angle, and only this step's action."""
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 10, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(label(STEP_PURPOSE[key], "small", True))
        strip = FlowLayout()
        action = button(STEP_ACTIONS[key], lambda checked=False, k=key: self.run_step(k), True)
        self.step_buttons[key] = action
        strip.addWidget(action)
        if key == "hydra":
            strip.addWidget(button("AMR databases / updates…", self.open_amr_databases))
        if key in {"st", "cgmlst"}:
            strip.addWidget(button("Install schemes…", self.open_reference_manager))
        if key == "assembly":
            strip.addWidget(button("Identify waiting samples again", self.reidentify_waiting))
        clear = button("Clear", lambda checked=False, k=key: self.clear_step(k))
        clear.setToolTip(CLEAR_NOTE)
        strip.addWidget(clear)
        layout.addLayout(strip)
        if key in RERUN_STEPS:
            from PySide6.QtWidgets import QCheckBox
            toggle = QCheckBox(RERUN_TOGGLE)
            toggle.setToolTip(RERUN_TOGGLE_TIP)
            self.step_toggles[key] = toggle
            layout.addWidget(toggle)
        note = label("", "small", True)
        note.setObjectName("stepGate")
        self.step_notes[key] = note
        layout.addWidget(note)
        if key == "hydra":
            layout.addWidget(self.build_hydra_organism_control())
        table = self.build_step_table(key)
        self.step_tables[key] = table
        layout.addWidget(table, 1)
        return page

    # Column layouts. ST and cgMLST deliberately have different column names and
    # different denominators: seven loci and several thousand targets are not the
    # same measurement and must never be read off the same column.
    STEP_COLUMNS = {
        "assembly": ["Sample", "Input", "Status", "ST (7-locus)", "MLST loci", "Genus", "Species",
                     "Organism evidence", "Scheme", "AMR genes", "Collection", "Storage"],
        "st": ["Sample", "Organism", "Organism review", "7-locus scheme", "ST", "Loci called",
               "Status"],
        "cgmlst": ["Sample", "Organism", "Organism review", "cgMLST scheme", "Targets in scheme",
                   "Loci called", "Missing", "Status"],
        "hydra": ["Sample", "Organism", "Organism review", "HYDRA organism", "AMR genes",
                  "Evidence state"],
    }
    # What a first-time user needs on the roster. The rest is a click away.
    ASSEMBLY_HIDDEN = (3, 4, 8, 9, 10)

    def build_step_table(self, key):
        if key == "assembly":
            # The roster keeps its long-standing name: every other view, export and
            # test addresses the isolate library through it.
            table = self.sample_table = make_table(self.STEP_COLUMNS["assembly"])
            self.install_view_menu("library", table)
        else:
            table = make_table(self.STEP_COLUMNS[key])
            self.install_view_menu(f"library.{key}", table)
        table.itemSelectionChanged.connect(lambda k=key: self.step_selection_changed(k))
        table.cellDoubleClicked.connect(
            lambda row, column, t=table: self.open_isolate_record(t.item(row, 0).data(Qt.ItemDataRole.UserRole)))
        return table

    def apply_column_visibility(self, show_all=None):
        """Hide the roster columns that belong to another tab, unless asked for them."""
        if not hasattr(self, "sample_table"):
            return
        show_all = self.sample_columns_button.isChecked() if show_all is None else bool(show_all)
        for column in range(self.sample_table.columnCount()):
            self.sample_table.setColumnHidden(column, not show_all and column in self.ASSEMBLY_HIDDEN)

    def refresh(self):
        if self._refreshing:
            return
        self._refreshing = True
        try:
            super().refresh()
            # Archived isolates keep every record they carry; they simply leave
            # the working views until the user restores them.
            working = active_samples(self.current_samples)
            archived = len(self.current_samples) - len(working)
            isolates = [sample for sample in working if sample.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
            mates = len(working) - len(isolates)
            self.metrics[0].value.setText(str(len(isolates)))
            hint = "in the active library"
            if mates:
                hint = f"{len(isolates)} isolates · {mates} linked mates"
            if archived:
                hint += f" · {archived} archived"
            self.metrics[0].hint.setText(hint)
            if self.library is not None and not (self.worker and self.worker.isRunning()):
                self.library.index_project(self.project)
            valid = {s["id"] for s in working}
            self.selection_ids.intersection_update(valid)
            self.report_ids.intersection_update(valid)
            self.focus.prune(valid)
            self.refresh_library()
            if hasattr(self, "cohort_table"):
                self.refresh_cohort_table()
            if hasattr(self, "report_table"):
                self.refresh_report_table()
            if hasattr(self, "feature_table"):
                self.refresh_features()
            self.refresh_journey()
        finally:
            self._refreshing = False

    def refresh_journey(self):
        if not hasattr(self, "investigation_map"):
            return
        from wmlstudio.journey import journey_summary
        summary = journey_summary(active_samples(self.current_samples), selected_ids=self.selection_ids,
                                  comparison_ids=self.cohort_ids)
        investigation = self.investigation_summary() if hasattr(self, "investigation_summary") else {}
        description = None
        if investigation.get("saved"):
            description = (f"Investigation: {investigation.get('name')} · {investigation.get('snapshots', 0)} saved snapshots · "
                           "New samples are included only after you review the cohort.")
        self.investigation_map.set_summary(summary, description)
        if hasattr(self, "scope_label"):
            counts = summary["counts"]
            self.scope_label.setText(f"{counts['total']} isolates in project")

    def journey_action(self, action):
        """Route problem-oriented actions without silently broadening a cohort."""
        if action == "import":
            self.browse_files()
        elif action == "guide":
            self.open_workflow_guide()
        elif action == "samples":
            self.navigate(1)
        elif action == "review":
            from wmlstudio.journey import journey_summary
            identifiers = journey_summary(active_samples(self.current_samples))["review_ids"]
            self.clear_filters()
            self.selection_ids = set(identifiers)
            self.library_filter = ("ids", set(identifiers))
            self.refresh_tables()
            self.navigate(1)
            self.notify(f"{len(identifiers)} flagged isolates shown. Clear filters to return to all samples.")
        elif action == "analyse":
            self.choose_and_analyse()
        elif action == "characterize":
            self.run_characterization_selected()
        elif action == "features":
            self.choose_feature_cohort()
        elif action == "compare":
            self.navigate(2)
        elif action == "reports":
            self.navigate(5)

    def refresh_library(self):
        tree = self.library_tree
        tree.blockSignals(True)
        tree.clear()
        working = active_samples(self.current_samples)
        all_item = QTreeWidgetItem([f"All samples  ({len(working)})"])
        all_item.setData(0, Qt.ItemDataRole.UserRole, None)
        tree.addTopLevelItem(all_item)
        groups = {}
        collections = {}
        review = {}
        for sample in working:
            if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate":
                continue
            bucket = quarantine_bucket(sample)
            if bucket:
                # Needs review is not an organism, so it never joins the genus tree.
                review.setdefault(bucket, []).append(sample)
                continue
            genus, species, _ = organism_for(sample)
            genus, species = genus or "Unknown", species or "Unspecified"
            st = (sample.get("result") or {}).get("st")
            st_label = f"ST {st}" if st else "ST unassigned"
            groups.setdefault(genus, {}).setdefault(species, {}).setdefault(st_label, []).append(sample)
            for name in (sample.get("metadata") or {}).get("collections", []):
                collections.setdefault(name, []).append(sample)
        for genus, species_groups in sorted(groups.items()):
            genus_node = QTreeWidgetItem([genus])
            genus_node.setData(0, Qt.ItemDataRole.UserRole, ("genus", genus))
            all_item.addChild(genus_node)
            for species, st_groups in sorted(species_groups.items()):
                species_node = QTreeWidgetItem([species])
                species_node.setData(0, Qt.ItemDataRole.UserRole, ("organism", genus, species))
                genus_node.addChild(species_node)
                for st, samples in sorted(st_groups.items()):
                    node = QTreeWidgetItem([f"{st}  ({len(samples)})"])
                    node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                    species_node.addChild(node)
                    for sample in samples:
                        child = QTreeWidgetItem([sample["name"]])
                        child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                        node.addChild(child)
        if review:
            total = sum(len(samples) for samples in review.values())
            branch = QTreeWidgetItem([f"Needs review  ({total})"])
            branch.setData(0, Qt.ItemDataRole.UserRole,
                           ("ids", [s["id"] for samples in review.values() for s in samples]))
            branch.setToolTip(0, "WMLSTudio declined to decide what these organisms are, so nothing "
                                 "was filed into a genus folder on a guess. Right-click to assign an "
                                 "organism; the managed copy moves with the label.")
            tree.addTopLevelItem(branch)
            for bucket, samples in sorted(review.items()):
                node = QTreeWidgetItem([f"{bucket.replace('_', ' ')}  ({len(samples)})"])
                node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                branch.addChild(node)
                for sample in samples:
                    child = QTreeWidgetItem([sample["name"]])
                    child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                    child.setToolTip(0, organism_evidence_note(sample))
                    node.addChild(child)
            branch.setExpanded(True)
        if collections:
            branch = QTreeWidgetItem(["Collections"])
            tree.addTopLevelItem(branch)
            for name, samples in sorted(collections.items()):
                node = QTreeWidgetItem([f"{name}  ({len(samples)})"])
                node.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in samples]))
                branch.addChild(node)
            branch.setExpanded(True)
        if self.library is not None:
            project_counts = {}
            for sample in self.library.search():
                path = sample["project_path"]
                if Path(path).resolve() != self.project_path.resolve():
                    project_counts[path] = project_counts.get(path, 0) + 1
            if project_counts:
                projects = QTreeWidgetItem(["Other saved projects"])
                tree.addTopLevelItem(projects)
                for path, count in sorted(project_counts.items()):
                    node = QTreeWidgetItem([f"{Path(path).stem}  ({count})"])
                    node.setData(0, Qt.ItemDataRole.UserRole, ("project", path))
                    projects.addChild(node)
                projects.setExpanded(True)
        all_item.setExpanded(True)
        stored = archived_samples(self.current_samples)
        if stored:
            branch = QTreeWidgetItem([f"Archived  ({len(stored)})"])
            branch.setData(0, Qt.ItemDataRole.UserRole, ("ids", [s["id"] for s in stored]))
            branch.setToolTip(0, "Archived isolates keep every result, allele call, file and history "
                                 "entry. They are hidden from the working views until you restore them.")
            tree.addTopLevelItem(branch)
            for sample in stored:
                child = QTreeWidgetItem([sample["name"]])
                child.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"]]))
                branch.addChild(child)
        mates = [sample for sample in active_samples(self.current_samples) if sample.get("metadata", {}).get("workflow", {}).get("source_kind") == "read_mate"]
        if mates:
            mate_branch = QTreeWidgetItem([f"Linked read mates  ({len(mates)})"])
            mate_branch.setData(0, Qt.ItemDataRole.UserRole, ("ids", [sample["id"] for sample in mates]))
            tree.addTopLevelItem(mate_branch)
        tree.blockSignals(False)

    def open_library_group(self, item, column=0):
        value = item.data(0, Qt.ItemDataRole.UserRole)
        if value and value[0] == "project":
            self.switch_project(value[1])
            return
        self.library_filter = value
        self.refresh_tables()
        self.navigate(1)

    def open_recent_sample(self, row, column=0):
        item = self.recent_table.item(row, 0)
        if item:
            self.open_isolate_record(item.data(Qt.ItemDataRole.UserRole))

    @staticmethod
    def fill_filter(combo, values, title):
        selected = combo.currentData()
        combo.blockSignals(True)
        combo.clear()
        combo.addItem(title, "")
        for value in sorted(set(v for v in values if v)):
            combo.addItem(value, value)
        index = combo.findData(selected)
        combo.setCurrentIndex(max(0, index))
        combo.blockSignals(False)

    def refresh_tables(self):
        if not hasattr(self, "sample_table"):
            return
        working = active_samples(self.current_samples)
        if self.library_filter and self.library_filter[0] == "ids":
            # An explicitly named group — the Archived branch, for instance — shows
            # the isolates it names, so a hidden isolate is never simply missing.
            named = set(self.library_filter[1])
            working = working + [s for s in archived_samples(self.current_samples) if s["id"] in named]
        taxa = [organism_for(s) for s in working]
        self.fill_filter(self.genus_filter, [t[0] for t in taxa], "All genera")
        self.fill_filter(self.species_filter, [t[1] for t in taxa], "All species")
        self.fill_filter(self.collection_filter, [c for s in working for c in
                         (s.get("metadata") or {}).get("collections", [])], "All collections")
        query = self.search.text().strip().casefold()
        visible = []
        for sample in working:
            genus, species, evidence = organism_for(sample)
            metadata = sample.get("metadata") or {}
            result = sample.get("result") or {}
            searchable = " ".join([sample["name"], genus, species, str(result.get("st") or ""),
                                   " ".join(gene_names(sample)), json.dumps(flattened_metadata(sample), ensure_ascii=False)])
            if query and query not in searchable.casefold():
                continue
            if self.genus_filter.currentData() and genus != self.genus_filter.currentData():
                continue
            if self.species_filter.currentData() and species != self.species_filter.currentData():
                continue
            if self.collection_filter.currentData() and self.collection_filter.currentData() not in metadata.get("collections", []):
                continue
            if self.library_filter:
                kind, *values = self.library_filter
                if kind == "ids" and sample["id"] not in values[0]:
                    continue
                if kind == "genus" and (genus or "Unknown") != values[0]:
                    continue
                if kind == "organism" and ((genus or "Unknown"), (species or "Unspecified")) != tuple(values):
                    continue
            visible.append(sample)
        self._visible_samples = visible
        recent = active_samples(self.current_samples)[-12:]
        self._typing_overview = self.typing_overview(
            list({sample["id"]: sample for sample in visible + recent}.values()))
        self.fill_sample_table(self.sample_table, visible)
        self.fill_sample_table(self.recent_table, recent)
        for key in STEP_KEYS:
            if key != "assembly" and key in self.step_tables:
                self.fill_step_table(key, visible)
        self.selection_label.setText(self.selection_summary(len(visible)))
        self.populate_hydra_organisms()
        self.update_step_gates()
        self.show_sample_detail()

    def selection_summary(self, visible):
        working = active_samples(self.current_samples)
        archived = len(self.current_samples) - len(working)
        review = sum(1 for sample in working if quarantine_bucket(sample))
        text = f"{len(self.selection_ids)} selected · {visible} visible · {len(working)} in project"
        if review:
            text += f" · {review} in Needs review"
        if archived:
            text += f" · {archived} archived"
        return text

    def scheme_organisms(self):
        """Installed references keyed by the name a stored result records for them."""
        return scheme_organism_index(self.installed_scheme_rows())

    def called_loci(self, sample, kind, row):
        """How many loci of one stored profile were actually called.

        The scheme, the ST and the number of targets come straight out of SQLite
        without decoding anything. Only this number needs the profile itself, so it
        is taken from the copy the sample already carries when that is the profile
        in question, and cached on the scheme fingerprint otherwise: a stored
        profile never changes under a fingerprint.
        """
        digest = row.get("scheme_digest")
        key = (sample["id"], kind, digest)
        if key in self._called_cache:
            return self._called_cache[key]
        headline = sample.get("result") or {}
        alleles = None
        if digest and headline.get("scheme_digest") == digest:
            alleles = headline.get("alleles")
        if not isinstance(alleles, dict):
            try:
                profile = self.project.latest_analysis(sample["id"], kind) or {}
            except (KeyError, ValueError):
                profile = {}
            alleles = profile.get("alleles") if isinstance(profile.get("alleles"), dict) else {}
        called = sum(1 for value in alleles.values() if value)
        if len(self._called_cache) > 8192:
            self._called_cache.clear()
        self._called_cache[key] = called
        return called

    def typing_overview(self, samples):
        """Each isolate's MLST and cgMLST profile summaries, kept strictly apart.

        The headline result is whichever analysis ran last, so reading a seven-locus
        ST off it would let a core-genome run erase it. These come from the stored
        analyses, addressed by typing kind, and a kind with no profile is absent —
        never filled in from the other one.
        """
        overview = {}
        for sample in samples:
            try:
                rows = self.project.analysis_summaries(sample["id"])
            except KeyError:
                rows = []  # Removed between the refresh and this read.
            latest = {}
            for row in rows:
                if row.get("typing_kind") in {"mlst", "cgmlst"}:
                    latest[row["typing_kind"]] = row  # Ordered oldest first; the last wins.
            entries = {}
            for kind, row in latest.items():
                loci = int(row.get("locus_count") or 0)
                called = self.called_loci(sample, kind, row)
                entries[kind] = {"scheme": row.get("scheme") or "", "st": row.get("st"),
                                 "status": row.get("status") or "", "loci": loci, "called": called,
                                 "scheme_digest": row.get("scheme_digest")}
            overview[sample["id"]] = entries
        return overview

    def typing_entry(self, sample, kind):
        return (self._typing_overview.get(sample["id"]) or {}).get(kind) or {}

    def fill_sample_table(self, widget, samples):
        """The roster: what each isolate is, what its organism is, and where it lives.

        The ST and loci columns are read from the isolate's classical MLST profile
        only. A cgMLST run is a different measurement against a different scheme, so
        it never writes an ST or a seven-locus denominator into these cells.
        """
        schemes = self.scheme_organisms()
        self._filling = True
        widget.blockSignals(True)
        widget.setSortingEnabled(False)
        widget.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            result = sample.get("result") or {}
            metadata = sample.get("metadata") or {}
            status = result.get("status", sample["status"]) if sample["status"] == "completed" else sample["status"]
            mlst = self.typing_entry(sample, "mlst")
            kind = result.get("kind") or ("fastq" if ".fq" in sample["input_path"] or ".fastq" in sample["input_path"] else "fasta")
            if metadata.get("workflow", {}).get("source_kind") == "profile":
                kind = "profile"
            if metadata.get("workflow", {}).get("source_kind") == "read_mate":
                kind, status = "read_mate", "linked_mate"
            genus, species, evidence = organism_for(sample)
            review = organism_review(sample, schemes, mlst)
            note = review["detail"] or organism_evidence_note(sample)
            storage = "Managed copy" if metadata.get("workflow", {}).get("managed") else "Linked original"
            if sample.get("missing_input") and kind != "profile":
                storage = "Input unavailable"
            if review["needs_review"]:
                genus, species = genus or "Needs review", species or "—"
                evidence = review["label"]
            if is_archived(sample):
                evidence = f"Archived · {evidence}"
            values = [sample["name"], {"fastq": "Reads", "profile": "Profile only", "read_mate": "Linked reverse mate"}.get(kind, "Assembly"),
                      status.replace("_", " ").title(), f"ST {mlst['st']}" if mlst.get("st") else "—",
                      f"{mlst['called']} / {mlst['loci']}" if mlst.get("loci") else "—",
                      genus or "Unknown", species or "—", evidence, mlst.get("scheme") or "—",
                      "; ".join(gene_names(sample)) or "—", "; ".join(metadata.get("collections", [])) or "—", storage]
            for column, value in enumerate(values[:widget.columnCount()]):
                item = cell(value, sample["id"])
                if column == 2:
                    item.setForeground(QColor("#5AE0BD" if status in {"complete", "qc_only", "completed"}
                                               else "#F5BE73" if status in {"incomplete", "mixed", "failed", "interrupted"}
                                               else "#A6B7D0"))
                if metadata.get("cluster", {}).get("highlight"):
                    item.setBackground(QColor("#343348"))
                if note and column in (5, 6, 7):
                    item.setToolTip(note)
                if column in (5, 6, 7):
                    paint_review(item, review)
                widget.setItem(row, column, item)
                if widget is self.sample_table and sample["id"] in self.selection_ids:
                    item.setSelected(True)
        widget.setSortingEnabled(True)
        widget.blockSignals(False)
        self._filling = False

    def step_row_values(self, sample, key, review):
        """One row of a step view, in that step's own units."""
        from wmlstudio.sample_workflow import hydra_evidence_status
        organism = organism_cell_text(sample)
        marker = review["label"] if review["needs_review"] else "Ready"
        if key == "st":
            mlst = self.typing_entry(sample, "mlst")
            return [sample["name"], organism, marker, mlst.get("scheme") or "—",
                    f"ST {mlst['st']}" if mlst.get("st") else "—",
                    f"{mlst['called']} / {mlst['loci']}" if mlst.get("loci") else "—",
                    (mlst.get("status") or "Not typed").replace("_", " ").title()]
        if key == "cgmlst":
            cgmlst = self.typing_entry(sample, "cgmlst")
            loci, called = cgmlst.get("loci") or 0, cgmlst.get("called") or 0
            return [sample["name"], organism, marker, cgmlst.get("scheme") or "—",
                    str(loci) if loci else "—", str(called) if loci else "—",
                    str(loci - called) if loci else "—",
                    (cgmlst.get("status") or "Not called").replace("_", " ").title()]
        state = hydra_evidence_status(sample)
        return [sample["name"], organism, marker, self.hydra_organism_for(sample) or "No organism flags",
                "; ".join(gene_names(sample)) or "—", state["status"].replace("_", " ").title()]

    def fill_step_table(self, key, samples):
        schemes = self.scheme_organisms()
        widget = self.step_tables[key]
        self._filling = True
        widget.blockSignals(True)
        widget.setSortingEnabled(False)
        widget.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            review = organism_review(sample, schemes, self.typing_entry(sample, "mlst"))
            for column, value in enumerate(self.step_row_values(sample, key, review)):
                item = cell(value, sample["id"])
                if column in (1, 2):
                    paint_review(item, review)
                if column == 2 and not review["needs_review"]:
                    item.setToolTip(review["detail"])
                widget.setItem(row, column, item)
                if sample["id"] in self.selection_ids:
                    item.setSelected(True)
        widget.setSortingEnabled(True)
        widget.blockSignals(False)
        self._filling = False

    def step_selection_changed(self, key):
        """One selection, whichever view the user clicked in."""
        if self._filling:
            return
        widget = self.step_tables.get(key)
        if widget is None:
            return
        visible = {widget.item(r, 0).data(Qt.ItemDataRole.UserRole) for r in range(widget.rowCount())}
        selected = {item.data(Qt.ItemDataRole.UserRole) for item in widget.selectedItems()}
        self.selection_ids.difference_update(visible)
        self.selection_ids.update(selected)
        self.focus.set_focus(self.selection_ids, "Isolate table selection")
        self.selection_label.setText(self.selection_summary(len(visible)))
        self.mirror_selection(key)
        self.update_step_gates()
        self.show_sample_detail()
        self.refresh_journey()

    def mirror_selection(self, source):
        """Show the same isolates as selected in the other step views."""
        self._filling = True
        try:
            for key, widget in self.step_tables.items():
                if key == source:
                    continue
                widget.blockSignals(True)
                widget.clearSelection()
                for row in range(widget.rowCount()):
                    item = widget.item(row, 0)
                    if item is not None and item.data(Qt.ItemDataRole.UserRole) in self.selection_ids:
                        widget.selectRow(row)
                widget.blockSignals(False)
        finally:
            self._filling = False

    def sample_selection_changed(self):
        self.step_selection_changed("assembly")

    def select_visible_samples(self):
        self.sample_table.selectAll()
        self.sample_selection_changed()

    def clear_sample_selection(self):
        self.selection_ids.clear()
        self.sample_table.clearSelection()
        self.sample_selection_changed()

    def clear_filters(self):
        self.library_filter = None
        self.search.clear()
        for combo in (self.genus_filter, self.species_filter, self.collection_filter):
            combo.setCurrentIndex(0)
        self.refresh_tables()

    def selected_samples(self):
        return [s for s in self.project.samples() if s["id"] in self.selection_ids]

    # --- what a step needs before it may run --------------------------------
    # The user's own words: "only when basic information for the pipeline to
    # continue are set then the fasta can be processed". So each step states the
    # facts it needs, the action stays disabled until they are set, and the note
    # under the button names the missing fact and where to set it — instead of
    # letting a run start and fail with something only a bioinformatician reads.

    def installed_scheme_path(self, path):
        """A recorded scheme path, followed to wherever the library moved it.

        A scheme the cgMLST library reorganised has to keep working for a sample
        that recorded where it used to be. Anything else is a configuration that
        silently stops running, which is exactly what the user reported.
        """
        if not path:
            return ""
        from wmlstudio import cgmlst_schemes
        try:
            return str(cgmlst_schemes.resolve_migrated_path(self.root, str(path)))
        except (OSError, ValueError):
            return str(path)

    def installed_scheme_rows(self):
        """Installed references with the organism each one records for itself.

        Cached on the list of scheme folders: reference_index already caches on
        each folder's file signature, so a refresh costs nothing until something
        is installed or removed.
        """
        key = tuple(str(path) for path in self.scheme_paths)
        if self._scheme_rows_key != key:
            from wmlstudio.reference_index import scheme_entries
            try:
                self._scheme_rows = scheme_entries(self.scheme_paths)
            except (OSError, ValueError):
                self._scheme_rows = []  # A reference that cannot be read blocks nothing.
            self._scheme_rows_key = key
        return self._scheme_rows

    def step_schemes(self, key, organism):
        """Installed schemes of this step's kind that name this organism.

        Exact species matches come first; a scheme labelled with the genus alone
        follows, because many curated panels cover a whole species complex.
        """
        if organism is None or key not in {"st", "cgmlst"}:
            return []
        genus, species = organism
        wanted = "mlst" if key == "st" else "cgmlst"
        matches = [entry for entry in self.installed_scheme_rows()
                   if entry.get("kind") == wanted
                   and str(entry.get("genus") or "").casefold() == str(genus).casefold()]
        exact = [e for e in matches if species and str(e.get("species") or "").casefold() == str(species).casefold()]
        return exact + [e for e in matches if e not in exact]

    def preferred_schemes(self, key, organism):
        """The best-matching tier only: exact species matches when there are any.

        One reference named for this exact species is not a choice worth
        interrupting for; two are, and the run says which one it used either way.
        """
        matches = self.step_schemes(key, organism)
        species = (organism or ("", ""))[1]
        exact = [entry for entry in matches
                 if species and str(entry.get("species") or "").casefold() == str(species).casefold()]
        return exact or matches

    def amr_database_state(self):
        """Whether HYDRA has a database to search, asked of the installed store itself."""
        key = self.active_amr_database()
        if self._amr_state_key != key:
            from wmlstudio import paths, provisioning
            try:
                requirement = provisioning.hydra_database_requirement(paths.data_root(), key)
                detail = requirement.detail or {}
                self._amr_state = {
                    # A partial store still runs and says what it could not search;
                    # only an empty or unusable one blocks the step.
                    "ready": requirement.state in {"ready", "partial"},
                    "reason": requirement.reason,
                    "organisms": list(detail.get("organisms") or []),
                    "route": requirement.action.ui_route if requirement.action else "",
                }
            except (OSError, ValueError) as error:
                self._amr_state = {"ready": False, "reason": str(error), "organisms": [],
                                   "route": "Data ▸ AMR databases / updates…"}
            self._amr_state_key = key
        return self._amr_state

    def step_gate(self, sample, key):
        """What this isolate still needs before this step may run, in plain words."""
        metadata = sample.get("metadata") or {}
        workflow = metadata.get("workflow") or {}
        kind = input_kind(sample)
        missing = []
        if kind == "profile":
            missing.append(("no sequence file", "This isolate was imported as a saved allele "
                            "profile. Profiles can be compared, but not assembled, typed or "
                            "screened."))
        elif kind == "read_mate":
            missing.append(("nothing of its own to run", "This record is the reverse mate of "
                            "another isolate; run the step on its forward partner."))
        elif sample.get("missing_input"):
            missing.append(("its sequence file", "The file recorded for it is not on this "
                            "computer. Use Samples ▸ Relink input (same bytes)…"))
        elif key == "assembly" and kind != "reads":
            missing.append(("nothing to assemble", "This isolate is already an assembly. Type it "
                            "on the ST or cgMLST tab."))
        elif key != "assembly" and kind == "reads":
            missing.append(("an assembly", "These are raw reads. Assemble them on the Assembly "
                            "tab first; reads are never typed directly."))
        notes = []
        if key in {"st", "cgmlst", "hydra"} and accepted_organism(sample) is None:
            if key == "hydra":
                # HYDRA itself runs without one; only the organism-specific tools
                # need it, so this is stated rather than used to block the run.
                notes.append(f"{sample['name']} has no accepted organism, so AMRFinderPlus "
                             "organism rules, Kleborate and SCCmec typing stay unavailable for "
                             "it. Choose one in the organism box above to enable them.")
            else:
                missing.append(("an accepted organism", "Use Assign organism…, or accept the "
                                "proposal in the Needs review branch of the library tree."))
        if key in {"st", "cgmlst"}:
            organism = accepted_organism(sample)
            pinned = workflow.get("scheme_path")
            entry = self.scheme_organisms().get(str(Path(pinned).name).casefold()) if pinned else None
            wanted = "mlst" if key == "st" else "cgmlst"
            if entry is not None and entry.get("kind") == wanted:
                pass  # An explicitly pinned scheme of the right kind settles it.
            elif organism is not None and not self.step_schemes(key, organism):
                title = "seven-locus MLST" if key == "st" else "cgMLST"
                missing.append((f"an installed {title} scheme for {organism_text(organism)}",
                                "Install one from Schemes ▸ Browse online / install updates…"))
        if key == "hydra" and not self.amr_database_state()["ready"]:
            state = self.amr_database_state()
            missing.append(("an installed AMR database",
                            state["reason"] or "Use Data ▸ AMR databases / updates…"))
        return {"ready": not missing, "missing": missing, "notes": notes}

    def step_partition(self, key, samples=None):
        """Split a selection into what this step can run and what it cannot, with reasons."""
        if samples is None:
            samples = self.selected_samples()
        ready, blocked = [], []
        for sample in samples:
            gate = self.step_gate(sample, key)
            (ready if gate["ready"] else blocked).append((sample, gate))
        return [sample for sample, _ in ready], blocked

    def step_gate_sentence(self, key):
        """The line under a step's button: what it would run, or exactly what is missing."""
        selected = self.selected_samples()
        if not selected:
            return (False, f"Select isolates in any tab, then {STEP_ACTIONS[key][:-1].casefold()}. "
                           "Nothing runs on isolates you have not chosen.")
        ready, blocked = self.step_partition(key, selected)
        caveats = " ".join(dict.fromkeys(
            note for sample in ready for note in self.step_gate(sample, key)["notes"]))
        toggle = self.step_toggles.get(key)
        if ready and toggle is not None and not toggle.isChecked():
            settled = sum(1 for sample in ready
                          if (self._typing_verdicts.get((key, sample["id"])) or {}).get("status")
                          == "current")
            standing = (f"{settled} of them already have a {STEP_TITLES[key]} result nothing has "
                        "invalidated, so they would be left alone. " if settled else "")
            caveats = f"{standing}{STANDING_NOTE} {caveats}".strip()
        if ready and not blocked:
            return (True, f"{len(ready)} of {len(selected)} selected isolates are ready for "
                          f"{STEP_TITLES[key]}. {caveats}".strip())
        reasons = []
        for sample, gate in blocked[:3]:
            fact, fix = gate["missing"][0]
            reasons.append(f"{sample['name']} needs {fact} — {fix}")
        more = f" …and {len(blocked) - 3} more." if len(blocked) > 3 else ""
        if ready:
            return (True, f"{len(ready)} of {len(selected)} are ready. "
                          f"{len(blocked)} cannot run yet: " + " ".join(reasons) + more
                          + (f" {caveats}" if caveats else ""))
        return (False, f"None of the {len(selected)} selected isolates can run {STEP_TITLES[key]} "
                       "yet: " + " ".join(reasons) + more)

    def update_step_gates(self):
        """Enable each step's action only when it would actually work."""
        if not self.step_buttons:
            return
        # The window's own idea of "a task is in progress": the thread object can
        # still report itself running inside its own finished handler, and a step
        # that is ready must not stay greyed until the user clicks something.
        busy = bool(self.worker_role and self.worker and self.worker.isRunning())
        for key, action in self.step_buttons.items():
            enabled, sentence = self.step_gate_sentence(key)
            action.setEnabled(enabled and not busy)
            action.setToolTip(sentence)
            note = self.step_notes.get(key)
            if note is not None:
                note.setText(sentence)

    def show_sample_detail(self, sample_id=None):
        sample = self.project.get_sample(sample_id) if sample_id else self.selected_sample()
        if not sample:
            self.detail.setHtml("<h3>Select an isolate to inspect its evidence</h3><p>Use Ctrl/Shift to select a cohort. Organism, typing, AMR and annotations stay linked by sample ID.</p>")
            if hasattr(self, "history_view"):
                self.history_view.clear()
            return
        def e(value):
            return html.escape(str(value))
        result = sample.get("result") or {}
        metadata = sample.get("metadata") or {}
        genus, species, evidence = organism_for(sample)
        # Addressed by typing kind: the headline result is only whichever analysis
        # ran last, and a core-genome run must not hide the seven-locus ST.
        mlst = self.project.latest_analysis(sample["id"], "mlst") or {}
        cgmlst = self.project.latest_analysis(sample["id"], "cgmlst") or {}
        review = organism_review(sample, self.scheme_organisms(), mlst)
        body = f"<h2>{e(sample['name'])}</h2><p><b>{e(' '.join([genus, species]).strip() or 'Unknown organism')}</b> · {e(evidence)} · ST {e(mlst.get('st') or 'unassigned')}</p>"
        if review["needs_review"]:
            body += f"<p><b>{e(review['label'])}</b><br>{e(review['detail'])}</p>"
        else:
            body += f"<p>{e(review['detail'])}</p>"
        # The two typings are reported separately and labelled with their own
        # denominators: seven loci and several thousand targets are not comparable,
        # and one must never stand in for the other.
        body += "<h3>Typing</h3><p>"
        if mlst:
            alleles = mlst.get("alleles") or {}
            called = sum(1 for value in alleles.values() if value)
            body += (f"<b>Classical MLST (7-locus):</b> {e(mlst.get('scheme') or 'scheme not recorded')} · "
                     f"ST {e(mlst.get('st') or 'unassigned')} · {called} of {len(alleles)} loci called<br>")
        else:
            body += "<b>Classical MLST (7-locus):</b> no profile. An ST is not implied by any other result.<br>"
        if cgmlst:
            targets = cgmlst.get("alleles") or {}
            called = sum(1 for value in targets.values() if value)
            body += (f"<b>cgMLST:</b> {e(cgmlst.get('scheme') or 'scheme not recorded')} · "
                     f"{called} of {len(targets)} targets called · {len(targets) - called} missing")
        else:
            body += "<b>cgMLST:</b> no profile."
        body += ("</p><p>A seven-locus distance and a core-genome distance are different "
                 "quantities against different schemes; they share no scale and no threshold.</p>")
        body += f"<p><b>Input:</b> {e(sample['input_path'] or 'Imported profile; no sequence attached')}<br><b>Sample ID:</b> {e(sample['id'])}</p>"
        if sample.get("error"):
            body += f"<p><b>Needs attention:</b> {e(sample['error'])}</p>"
        workflow = metadata.get("workflow", {})
        if workflow.get("paired_with"):
            try:
                paired = self.project.get_sample(workflow["paired_with"])
                body += f"<p><b>Paired record:</b> {e(paired['name'])} · {e(paired['id'])}<br><b>Read role:</b> {e(workflow.get('source_kind'))}</p>"
            except KeyError:
                body += "<p>The paired record is not present in this project; its provenance is retained.</p>"
        # One configuration, read the one way, with who set each field and where.
        # That is what makes a choice made in any menu visible in every other one.
        from wmlstudio.sample_workflow import sample_configuration
        configuration = sample_configuration(sample)
        body += (f"<p><b>Workflow:</b> {e(configuration['typing_mode'] or 'manual')} · "
                 f"<b>MLST scheme:</b> {e(result.get('scheme') or configuration['scheme_path'] or 'Not assigned')}"
                 f" · <b>cgMLST scheme:</b> {e(configuration['cgmlst_scheme_path'] or 'Not chosen')}"
                 f" · <b>Run flags:</b> HYDRA {'on' if configuration['run_hydra'] else 'off'}, "
                 f"cgMLST {'on' if configuration['run_cgmlst'] else 'off'}</p>")
        stamps = "; ".join(
            f"{field} by {entry.get('by') or 'unknown'}"
            + (f" in {entry['surface']}" if entry.get("surface") else "")
            + (" (automatic)" if entry.get("automatic") else "")
            for field, entry in sorted(configuration["set_by"].items()))
        if stamps:
            body += f"<p><b>Configuration set:</b> {e(stamps)}</p>"
        if configuration["organism_note"]:
            body += f"<p>{e(configuration['organism_note'])}</p>"
        if workflow.get("managed"):
            body += f"<p><b>Managed copy:</b> original remains at {e(workflow.get('source_path', 'recorded in provenance'))}</p>"
        qc = result.get("qc") or {}
        metrics = [("Records", "records"), ("Bases", "total_bases"), ("N50", "n50"), ("GC %", "gc_percent"), ("Q30 %", "q30_percent")]
        body += "<p>" + " · ".join(f"<b>{title}:</b> {e(round(qc[key], 2) if isinstance(qc[key], float) else qc[key])}" for title, key in metrics if qc.get(key) is not None) + "</p>"
        if qc.get("sampled"):
            body += "<p><b>Sampled read QC:</b> statistics describe an inspected prefix, not the entire read file.</p>"
        if metadata.get("hydra"):
            from wmlstudio.sample_workflow import hydra_evidence_status
            amr_state = hydra_evidence_status(sample)
            body += f"<h3>Linked HYDRA evidence · {e(amr_state['status'])}</h3><p>{e(amr_state['reason'])}</p>"
            if amr_state["status"] != "stale":
                body += "<p>" + e("; ".join(gene_names(sample)) or "No primary AMR genes reported in linked evidence") + "</p>"
        annotations = {k: v for k, v in metadata.items() if k not in {"workflow", "organism", "hydra", "assembly"}}
        if annotations:
            body += "<h3>Annotations</h3><p>" + "<br>".join(f"<b>{e(k)}:</b> {e(v)}" for k, v in annotations.items()) + "</p>"
        calls = result.get("calls") or []
        if calls:
            # Named, because these calls belong to one scheme — the most recent run —
            # and not to whichever typing the reader has in mind.
            body += (f"<h3>Allele calls · {e(result.get('scheme') or 'most recent analysis')}</h3>"
                     "<table width='100%' cellpadding='6'><tr bgcolor='#253650'><th>Locus</th>"
                     "<th>Allele</th><th>Evidence</th></tr>")
            for call in calls[:100]:
                body += f"<tr><td>{e(call['locus'])}</td><td>{e(call.get('allele') or '—')}</td><td>{e(call.get('status'))}</td></tr>"
            body += "</table>"
            if len(calls) > 100:
                body += f"<p>Showing 100 of {len(calls)} loci. The profile table and exports retain every locus.</p>"
        for note in result.get("notes", []):
            body += f"<p>{e(note)}</p>"
        body += f"<p><b>Input SHA-256:</b> {e(result.get('input_sha256') or workflow.get('source_sha256') or 'Not analysed')}<br><b>Scheme fingerprint:</b> {e(result.get('scheme_digest') or 'Not used')}</p>"
        self.detail.setHtml(body)
        if hasattr(self, "history_view"):
            history = self.project.history(sample["id"])
            rows = []
            for entry in reversed(history[-100:]):
                details = entry["details"]
                brief = {key: value for key, value in details.items() if not isinstance(value, (dict, list))}
                archived = details.get("result")
                if isinstance(archived, dict):
                    brief.update(archived_scheme=archived.get("scheme"), archived_st=archived.get("st"), archived_input_sha256=archived.get("input_sha256"))
                rows.append(f"<h3>{e(entry['action'].replace('_', ' '))}</h3><p>{e(entry['created_at'])}</p><p>" + "<br>".join(f"<b>{e(key)}:</b> {e(value)}" for key, value in brief.items()) + "</p>")
            self.history_view.setHtml(f"<h2>Sample audit trail · {e(sample['name'])}</h2><p>Showing {min(100, len(history))} of {len(history)} saved events. Original inputs and archived analysis evidence are retained.</p>" + "".join(rows))

    def scheme_entries(self, kind=None):
        """(title, path) for installed schemes, optionally of one kind only.

        The title is reference_index's row title — organism, scheme, target count,
        provider, version — so a downloaded cgMLST scheme reads as "Klebsiella
        pneumoniae sensu lato · cgMLST · 2358 targets · cgMLST.org" rather than as
        its folder name with the underscores rubbed out. Pass kind="mlst" or
        kind="cgmlst" wherever a picker offers a scheme to type against: a
        seven-locus distance and a 2,000-target distance are different quantities
        and must not be offered from one list. Several kinds may be named together,
        which is how a classical picker keeps offering a local scheme whose kind
        could not be read: leaving it out would hide a scheme somebody installed.
        """
        rows = self.installed_scheme_rows()
        if not rows:  # A reference that cannot be read must still be selectable.
            return [(p.name, str(p)) for p in self.scheme_paths]
        wanted = {kind} if isinstance(kind, str) else set(kind) if kind is not None else None
        chosen = [row for row in rows if wanted is None or row.get("kind") in wanted]
        # Two snapshots of one scheme read alike, and a picker that offers the same
        # words twice cannot say which one a result was produced against.
        repeated = {row["title"] for row in chosen
                    if sum(other["title"] == row["title"] for other in chosen) > 1}
        return [(f"{row['title']} · {Path(row['path']).name}" if row["title"] in repeated
                 else row["title"], row["path"]) for row in chosen]

    def populate_schemes(self):
        super().populate_schemes()
        self.scheme_combo.setItemText(0, "Use each sample's workflow")
        # Installing or removing a reference changes which steps can run at all.
        self._scheme_rows_key = None
        if self.step_buttons:
            self.update_step_gates()

    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Add samples", "", "FASTA / FASTQ (*.fa *.fasta *.fna *.fq *.fastq *.gz *.bz2);;All files (*)")
        if paths:
            self.intake_paths(paths)

    def browse_files_with_options(self):
        """The per-file import dialog, for the cases the automatic route cannot answer."""
        paths, _ = QFileDialog.getOpenFileNames(self, "Import sequences", "", "FASTA / FASTQ (*.fa *.fasta *.fna *.fq *.fastq *.gz *.bz2);;All files (*)")
        if paths:
            self.import_paths(paths, configure=True)

    def browse_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Add a folder of samples")
        if path:
            self.intake_paths([path])

    # --- drag and drop ------------------------------------------------------
    # The fastest way to start work, and the reason it is fast is that a drop
    # means what the thing dropped means, where it landed. It never means less
    # care: sequences are still identified where they lie before a byte is copied,
    # and the same review dialog still decides what is filed under which organism.

    def dragEnterEvent(self, event):
        if not (event.mimeData().hasUrls()
                and any(url.isLocalFile() for url in event.mimeData().urls())):
            return
        event.acceptProposedAction()
        self.statusBar().showMessage(DROP_HINT, 8000)

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls() and any(url.isLocalFile() for url in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        position = event.position().toPoint() if hasattr(event, "position") else None
        paths = [url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()]
        event.acceptProposedAction()
        self.handle_drop(paths, self.drop_target(position))

    def show_step(self, key):
        """Bring one Samples step sub-tab forward, found by the title it carries."""
        tabs = getattr(self, "sample_tabs", None)
        if tabs is None or key not in STEP_TITLES:
            return False
        for index in range(tabs.count()):
            if tabs.tabText(index) == STEP_TITLES[key]:
                tabs.setCurrentIndex(index)
                return True
        return False

    def drop_target(self, position=None):
        """Where a drop landed: the page, the step sub-tab, and the row under it.

        The pointer decides when it is over a row of isolates, because that is the
        one case where a drop is about one named isolate. Everywhere else the
        visible page decides, which is also what happens on a platform that hands
        us no position at all.
        """
        target = {"page": "", "step": "", "sample_id": ""}
        if hasattr(self, "pages"):
            target["page"] = self.pages.current_key()
        if target["page"] == "isolates" and hasattr(self, "sample_tabs"):
            title = self.sample_tabs.tabText(self.sample_tabs.currentIndex())
            target["step"] = next((key for key, value in STEP_TITLES.items() if value == title), "")
        if position is None:
            return target
        widget, table, step = self.childAt(position), None, ""
        while widget is not None and table is None:
            for key, candidate in self.step_tables.items():
                if widget is candidate or widget is candidate.viewport():
                    table, step = candidate, key
            widget = widget.parentWidget()
        if table is None:
            return target
        target["step"] = step
        item = table.itemAt(table.viewport().mapFrom(self, position))
        if item is not None:
            target["sample_id"] = item.data(Qt.ItemDataRole.UserRole) or ""
        return target

    def handle_drop(self, paths, target=None):
        """Do the obvious thing with what was dropped where. Returns the route taken."""
        target = dict(target or {})
        if self.busy():
            return ""
        dropped = sort_dropped(paths, schemes_expected=target.get("page") == "schemes")
        if not any(dropped[kind] for kind in ("sequences", "schemes", "tables")):
            self.error("Nothing here can be used: drop FASTA or FASTQ sequences, a folder of them, "
                       "a scheme folder, or a CSV/TSV of epidemiology metadata.")
            return ""
        if dropped["schemes"] and not dropped["sequences"]:
            return self.install_scheme_folders(dropped["schemes"])
        if dropped["tables"] and not dropped["sequences"]:
            return self.review_dropped_table(dropped["tables"][0])
        reads = [path for path in dropped["sequences"] if plain_suffix(path) in READ_SUFFIXES]
        if target.get("sample_id") and reads and len(reads) == len(dropped["sequences"]):
            return self.attach_dropped_reads(target["sample_id"], reads)
        # Everything else is a cohort arriving. The step the drop landed on is
        # where the user is working, so that is the tab they are returned to.
        self._drop_step = target.get("step") or ""
        self.intake_paths(dropped["sequences"])
        return "import"

    def install_scheme_folders(self, folders):
        """Install dropped folders of allele FASTA as local typing schemes.

        One at a time, because each is validated and copied by the same worker the
        Import local scheme… menu uses: a folder that is not a readable scheme is
        refused by name rather than half-installed.
        """
        from wmlstudio.jobs import SchemeImportWorker
        folder = Path(folders[0])
        if len(folders) > 1:
            self.notify(f"{len(folders)} folders were dropped; installing {folder.name} first. "
                        "Drop the rest once this one is in.")
        self.worker_role = "scheme_import"
        self._task_succeeded = False
        self._installed_scheme = None
        self.worker = SchemeImportWorker(folder, self.root, self)
        self.worker.imported.connect(
            lambda path, loci: setattr(self, "_installed_scheme", (path, loci)))
        self.worker.failed.connect(self.error)
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.set_running(True)
        self.worker.start()
        return "scheme"

    def dropped_scheme_installed(self):
        """Say what was installed and which library it joined; choose nothing for the user."""
        installed, self._installed_scheme = getattr(self, "_installed_scheme", None), None
        if installed is None:
            self.notify("That folder was not installed as a scheme; nothing was copied.")
            return
        destination, loci = installed
        self.populate_schemes()
        # Followed to where it now is: preparing the library moves a cgMLST scheme
        # out of the classical folder, and the row is found by where it landed.
        landed = Path(self.installed_scheme_path(destination)).expanduser().resolve()
        kind = {"mlst": "classical MLST", "cgmlst": "cgMLST"}.get(
            next((row.get("kind") for row in self.installed_scheme_rows()
                  if Path(row.get("path") or "") == landed), ""), "unclassified")
        self.notify(f"Installed {Path(destination).name} · {loci} loci · into the {kind} library. "
                    "Nothing was selected for you: choose it on the step that uses it. Results "
                    "already produced keep their own scheme fingerprint.")

    def attach_dropped_reads(self, sample_id, paths):
        """Offer dropped FASTQ files as the original reads of the isolate they landed on.

        The same review and the same full validation as the menu route: the pair is
        shown for confirmation and both files are hashed before anything is
        recorded. Filename or read-ID agreement is not proof that these reads
        produced that assembly, which is why a person still confirms it.
        """
        from wmlstudio.read_attachment_dialog import AttachReadsDialog
        try:
            sample = self.project.get_sample(sample_id)
        except KeyError:
            self.notify("That isolate is no longer in this project.")
            return ""
        if input_kind(sample) != "assembly":
            self.notify(f"{sample['name']} is not an assembly, so original reads cannot be attached "
                        "to it. Drop reads on an assembled isolate, or import them as their own "
                        "samples.")
            return ""
        dialog = AttachReadsDialog([sample], parent=self)
        for mate, path in zip((1, 2), sorted(paths)[:2], strict=False):
            combo = dialog.rows[0][mate - 1]
            combo.addItem(Path(path).name, {"path": path, "sample_id": None})
            combo.setCurrentIndex(combo.count() - 1)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return ""
        self.launch_read_validation(dialog.assignments, {sample_id: sample})
        return "attach_reads"

    def launch_read_validation(self, assignments, snapshots):
        """Hash and check each proposed read pair in a worker, then record what passed."""
        from wmlstudio.read_attachments import attach_read_pair, validate_read_attachment
        from wmlstudio.scheduler import plan_resources, run_bounded
        project = self.project
        try:
            allocation = plan_resources(threads_per_sample=1, memory_gb=1, max_parallel=2)
        except ValueError as exc:
            self.error(str(exc))
            return False

        def operation(cancelled, progress):
            identifiers = []

            def validate(task, resources, stopped, report):
                return validate_read_attachment(snapshots[task["sample_id"]], task["read1"],
                                                task["read2"], read_sample_ids=task["read_sample_ids"],
                                                cancelled=stopped, progress=report)

            def attach(task, evidence):
                attach_read_pair(project, task["sample_id"], evidence)
                identifiers.append(task["sample_id"])

            run_bounded(assignments, validate, allocation, cancelled=cancelled, on_result=attach,
                        progress=progress)
            return identifiers

        return self.launch_task(operation, "read_attachment", lambda ids: self.notify(
            f"Attached fully validated original read pairs to {len(ids)} isolates. Their assembly "
            "IDs, typing and AMR evidence are unchanged, and the read files were not moved."))

    def review_dropped_table(self, path):
        """Preview a dropped epidemiology table against this project before anything is saved."""
        from wmlstudio.metadata_grid import MetadataPreviewDialog
        from wmlstudio.metadata_ingest import read_metadata_table
        try:
            rows = read_metadata_table(path)
        except (OSError, ValueError) as error:
            self.error(str(error))
            return ""
        dialog = MetadataPreviewDialog(self.project, rows, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.refresh()
            self.notify(f"{len(dialog.changed_ids)} isolates updated from {Path(path).name}. "
                        "A spreadsheet is your statement about these isolates, not genomic evidence.")
        return "metadata"

    # --- clear --------------------------------------------------------------

    #: The Samples step each pipeline station's Clear stands for. A station clear
    #: and this hub's own Clear must leave the same thing behind, or one tab would
    #: keep a selection another says it just dropped.
    STEP_FOR_PAGE = {"isolates": None, "mlst": "st", "cgmlst": "cgmlst"}

    def clear_step(self, key=None, *, announce=True):
        """Put a tab back to 'nothing chosen' so the next thing can be anything.

        Clearing is about what you are working on at this moment. No isolate,
        result, allele call, file or history entry is touched, so after a clear the
        same tab is ready for new samples, old ones, or a mix of both. A step clears
        the isolates chosen for it and its own options; the hub also clears the
        search box, the organism filters and the library branch.
        """
        key = key if isinstance(key, str) and key in STEP_KEYS else None
        self.selection_ids.clear()
        for widget in self.step_tables.values():
            widget.clearSelection()
        for name, toggle in self.step_toggles.items():
            if key in (None, name):
                toggle.setChecked(False)
        if key is None:
            self.library_filter = None
            self.search.clear()
            for combo in (self.genus_filter, self.species_filter, self.collection_filter):
                combo.setCurrentIndex(0)
        self.focus.clear()
        self.refresh_tables()
        if announce:
            self.notify(CLEAR_NOTE if key is None else
                        f"{STEP_TITLES[key]} cleared: the isolates chosen for it and this tab's "
                        "own options. " + CLEAR_KEEPS)
        return {"page": key or "isolates", "selected": 0, "removed": 0}

    def clear_isolates(self):
        """The Samples tab's own Clear, under the name a page-level Clear looks for."""
        return self.clear_step()

    def clear_page_state(self, key):
        """Extend the window's page-level Clear with the chooser state this hub owns.

        The window clears what it can see. The selection set, the library branch and
        the organism filters live here, so both halves go together and a tab that
        says it started again really has.
        """
        parent = getattr(super(), "clear_page_state", None)
        if callable(parent):
            parent(key)
        if key in self.STEP_FOR_PAGE and hasattr(self, "sample_table"):
            self.clear_step(self.STEP_FOR_PAGE[key], announce=False)

    @staticmethod
    def sequence_files(paths):
        """Every readable sequence file in a selection or a dropped folder, in reading order."""
        files, seen = [], set()
        for source in paths:
            path = Path(source)
            # Sorted, so a dropped folder is reviewed in the order the user sees it
            # in their own file manager rather than in filesystem order. The key is
            # explicit because sorting Path objects is case-insensitive on Windows
            # and case-sensitive elsewhere: the same folder would otherwise import
            # in a different order on each platform, and when two files hold the
            # same bytes it decides which one becomes the isolate.
            candidates = (sorted(path.rglob("*"), key=lambda entry: (str(entry).casefold(), str(entry)))
                          if path.is_dir() else [path])
            for candidate in candidates:
                name = candidate.name.lower().removesuffix(".gz").removesuffix(".bz2")
                if candidate.is_file() and name.endswith((".fa", ".fasta", ".fna", ".fq", ".fastq")):
                    resolved = str(candidate.resolve())
                    if resolved not in seen:
                        seen.add(resolved)
                        files.append(resolved)
        return files

    def intake_paths(self, paths):
        """Drag and drop: identify each original, then store it under its organism.

        One action, one pass over the user's own files. Identification runs on the
        originals before a byte is copied; a strong genome-comparison call files
        itself into its genus and species folder, and everything weaker is still
        imported and waits under a named reason instead of being guessed at.
        """
        if self.busy():
            return
        files = self.sequence_files(paths)
        if not files:
            self.error("No supported FASTA or FASTQ files were found.")
            return
        from wmlstudio.characterization_refs import bundled_reference_root
        from wmlstudio.storage import intake_samples
        project = self.project
        root = str(self.project_path.with_suffix(".files"))
        panel = self.installed_species_panel()
        starter = bundled_reference_root()
        scheme_paths = list(self.scheme_paths)

        def operation(cancelled, progress):
            return intake_samples(project, files, storage_root=root, species_panel_root=panel,
                                  kpsc_panel_root=starter, scheme_paths=scheme_paths,
                                  cancelled=cancelled, progress=progress)

        self._intake_report = None
        self.launch_task(operation, "intake",
                         lambda report: setattr(self, "_intake_report", report))

    def intake_completed(self):
        """Say what was loaded, what was already here, and what still needs a person."""
        report, self._intake_report = self._intake_report, None
        step, self._drop_step = self._drop_step, ""
        if not report or not self._task_succeeded:
            self.notify("Loading stopped. Nothing was imported and your files are unchanged.")
            return
        imported = list(report.get("imported") or ())
        self.selection_ids = set(imported)
        self.clear_filters()
        self.refresh()
        self.navigate(1)
        if step in STEP_KEYS and hasattr(self, "sample_tabs"):
            # Dropped on a step, so the isolates arrive selected on that step, with
            # its own button and its own gate already saying what they still need.
            # Found by its title, because a later round adds sub-tabs of its own.
            self.show_step(step)
            self.mirror_selection("")
            self.update_step_gates()
        waiting = list(report.get("needs_review") or ())
        parts = [f"{len(imported)} isolates loaded."]
        # Everything imported is "filed" somewhere; only what left the review tree
        # was filed under an organism, and that is what this sentence may claim.
        by_organism = len(imported) - len(waiting)
        if by_organism:
            parts.append(f"{by_organism} were filed under the organism the genome comparison "
                         "supported.")
        skipped = list(report.get("skipped") or ())
        if skipped:
            parts.append(f"{len(skipped)} files were already in this project and were not copied "
                         "again.")
        if waiting:
            parts.append(f"{len(waiting)} are marked '{REVIEW_PREFIX}' and shown in amber: no "
                         "organism was decided for them, so nothing was filed on a guess.")
        notice = ((report.get("identification") or {}).get("notice") or "").strip()
        if notice:
            parts.append(notice)
        self.notify(" ".join(parts))

    def reidentify_waiting(self):
        """Run identification again on the isolates still waiting for an organism.

        Installing a reference panel does not change anything on its own: this is
        the second half of that sentence, and it re-identifies only what is waiting.
        A decision a person already made is never overwritten.
        """
        if self.busy():
            return
        waiting = [sample for sample in active_samples(self.current_samples)
                   if quarantine_bucket(sample)]
        if not waiting:
            self.notify("No isolates are waiting for an organism.")
            return
        from wmlstudio.characterization_refs import bundled_reference_root
        from wmlstudio.storage import reidentify_samples
        project = self.project
        root = str(self.project_path.with_suffix(".files"))
        panel = self.installed_species_panel()
        starter = bundled_reference_root()
        scheme_paths = list(self.scheme_paths)

        def operation(cancelled, progress):
            return reidentify_samples(project, species_panel_root=panel, kpsc_panel_root=starter,
                                      scheme_paths=scheme_paths, storage_root=root,
                                      cancelled=cancelled, progress=progress)

        self.launch_task(operation, "reidentify", self.reidentification_completed)

    def reidentification_completed(self, report):
        self.refresh()
        accepted = report.get("accepted") or ()
        waiting = report.get("waiting") or ()
        message = (f"{len(accepted)} isolates were identified and filed under their organism; "
                   f"{len(waiting)} are still waiting.")
        if waiting:
            buckets = sorted({bucket.replace('_', ' ') for _, bucket in waiting})
            message += " Reasons: " + ", ".join(buckets) + "."
        notice = ((report.get("identification") or {}).get("notice") or "").strip()
        if notice:
            message += " " + notice
        self.notify(message)

    def import_paths(self, paths, configure=False):
        if not configure:
            # Direct, non-interactive API retained for scripted integration tests.
            return super().import_paths(paths)
        if self.busy():
            return
        files = self.sequence_files(paths)
        if not files:
            self.error("No supported FASTA or FASTQ files were found.")
            return
        from wmlstudio.workflow_dialogs import ImportSamplesDialog
        dialog = ImportSamplesDialog(files, self.scheme_entries(CLASSICAL_KINDS), self)
        dialog.storage_root.setText(str(self.project_path.with_suffix(".files")))
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if any(assignment.get("typing_mode") == "auto" for assignment in dialog.assignments):
            self.identify_before_import(dialog.assignments, dialog.options)
            return
        self.import_assignments(dialog.assignments, dialog.options)

    def installed_species_panel(self):
        """The broad ANI panel the user installed, or None when only the starter is present."""
        from wmlstudio import organism_panel, paths
        try:
            return organism_panel.installed_species_panel(paths.data_root())
        except OSError:
            return None

    def organism_suggestions(self):
        """Organism names worth offering: installed references first, then this project."""
        pairs = {}
        try:
            from wmlstudio.organism_panel import installed_species_panel
            from wmlstudio.paths import data_root
            from wmlstudio.reference_index import panel_entries
            from wmlstudio.reference_index import scheme_entries as installed_scheme_entries
            panel = installed_species_panel(data_root())
            entries = list(installed_scheme_entries(self.scheme_paths))
            entries += list(panel_entries(panel)) if panel else []
            for entry in entries:
                if entry.get("genus"):
                    pairs.setdefault((str(entry["genus"]), str(entry.get("species") or "")), None)
        except (OSError, ValueError):
            pass  # Suggestions are a convenience; a missing reference never blocks import.
        for sample in self.current_samples:
            genus, species, _ = organism_for(sample)
            if genus:
                pairs.setdefault((genus, species or ""), None)
        return sorted(pairs)

    def identify_before_import(self, assignments, options):
        """Identify the user's own files first, so nothing is copied to a folder it must leave.

        Identification runs on the originals, before a single byte is copied, so a
        cancelled or rejected review leaves the project exactly as it was.
        """
        from wmlstudio.characterization_refs import bundled_reference_root
        from wmlstudio.scheduler import resources_for_run
        originals = [assignment["path"] for assignment in assignments
                     if assignment.get("typing_mode") == "auto"]
        panel = self.installed_species_panel()
        scheme_paths = list(self.scheme_paths)
        starter = bundled_reference_root()

        def operation(cancelled, progress):
            from wmlstudio.organism_id import identify_batch
            return identify_batch(originals, allocation=resources_for_run({"threads": 2, "memory_gb": 2}),
                                  species_panel_root=panel, kpsc_panel_root=starter,
                                  scheme_paths=scheme_paths, cancelled=cancelled, progress=progress)

        self._pending_import = {"assignments": list(assignments), "options": dict(options),
                                "verdicts": []}
        if not self.launch_task(operation, "identify", self.identification_completed):
            self._pending_import = None

    def identification_completed(self, verdicts):
        if self._pending_import is not None:
            self._pending_import["verdicts"] = list(verdicts)

    def review_identification(self):
        """Show what was proposed, and copy only what the user accepts."""
        pending, self._pending_import = self._pending_import, None
        if not pending:
            return
        if not self._task_succeeded:
            self.notify("Identification stopped. Nothing was imported and your files are unchanged.")
            return
        from wmlstudio.workflow_dialogs import IdentificationReviewDialog
        options = pending["options"]
        root = options.get("storage_root") or str(self.project_path.with_suffix(".files"))
        names = {str(Path(a["path"]).expanduser().resolve()): a.get("name")
                 for a in pending["assignments"] if a.get("name")}
        dialog = IdentificationReviewDialog(pending["verdicts"], root, self,
                                            organisms=self.organism_suggestions(), names=names,
                                            policy=self.project)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            self.notify("Nothing was imported. Your files are where you left them and no folders "
                        "were created.")
            return
        reviewed = {str(Path(entry["path"]).expanduser().resolve()): entry
                    for entry in dialog.assignments}
        assignments = []
        for assignment in pending["assignments"]:
            if assignment.get("typing_mode") != "auto":
                assignments.append(assignment)
                continue
            chosen = reviewed.get(str(Path(assignment["path"]).expanduser().resolve()))
            if chosen is not None:
                # The review decides the organism, which is what files the copy. How
                # the typing scheme is chosen is a separate question the user already
                # answered in the import dialog, so their answer is kept.
                assignments.append({**assignment, **chosen,
                                    "typing_mode": assignment.get("typing_mode", "auto")})
        if not assignments:
            self.notify("No files were selected for import. Nothing was copied.")
            return
        self.import_assignments(assignments, options, duplicates="skip")

    def import_assignments(self, assignments, options, *, duplicates="skip"):
        # Offering the same file twice is a mistake, not an instruction: the same
        # bytes stay one isolate whichever route imported them.
        from wmlstudio.storage import import_samples
        root = options.get("storage_root") or str(self.project_path.with_suffix(".files"))
        notes = []
        self._import_notes = notes
        def operation(cancelled, progress):
            return import_samples(
                self.project, assignments, storage_root=root, managed=options.get("managed", True),
                append_st=options.get("append_st", False), cancelled=cancelled, progress=progress,
                duplicates=duplicates, notes=notes)
        self.launch_task(operation, "import", lambda ids: self.import_completed(ids))

    def import_completed(self, ids):
        skipped = list(getattr(self, "_import_notes", ()) or ())
        self._import_notes = []
        self.selection_ids = set(ids)
        self.clear_filters()
        self.refresh()
        self.navigate(1)
        chosen = set(ids)
        review = [sample for sample in self.current_samples
                  if sample["id"] in chosen and quarantine_bucket(sample)]
        message = f"Imported {len(ids)} samples."
        if skipped:
            message += (f" {len(skipped)} file(s) already in this project were not copied again.")
        if review:
            message += (f" {len(review)} are in Needs review: no installed reference supported an "
                        "organism, so nothing was filed into a genus folder on a guess.")
        else:
            message += " Review their assignments, then analyse when ready."
        self.notify(message)

    def refile_completed(self, report):
        """Say what moved, what did not, and why — never a bare success."""
        self.refresh()
        moved, skipped = len(report.get("moved", ())), report.get("skipped", ())
        message = (f"{moved} managed copies were filed to match their organism. "
                   "Your original files were not moved.")
        if skipped:
            message += f" {len(skipped)} were left where they are: {skipped[0][1]}"
        self.notify(message)

    def apply_organism_assignments(self, assignments, *, notice="", surface="Assign organism"):
        """Record corrected organisms, then move the managed copies to match them.

        The write goes through the one sample configuration, so an organism or a
        scheme set here is what every other menu and submenu reads next, and the
        surface it was set from is stamped on the field. A corrected label that did
        not move the file would leave the file sitting in a folder that contradicts
        it, so the two always travel together.
        """
        from wmlstudio.storage import confirm_organism, refile_samples, set_sample_configuration
        assignments = [dict(assignment) for assignment in assignments if assignment.get("sample_id")]
        if not assignments:
            return False

        def operation(cancelled, progress):
            report = {"moved": [], "unchanged": [], "skipped": []}
            for index, assignment in enumerate(assignments):
                progress(index, len(assignments), f"Filing {index + 1} of {len(assignments)}…")
                sample_id = assignment["sample_id"]
                genus, species = assignment.get("genus", ""), assignment.get("species", "")
                mode = assignment.get("typing_mode", "manual")
                changes = {"genus": genus, "species": species, "typing_mode": mode,
                           "scheme_path": assignment.get("scheme_path")}
                try:
                    written = set_sample_configuration(self.project, [sample_id], changes,
                                                       surface=surface)
                    if genus and not written["applied"].get(sample_id):
                        # Re-affirming the label a waiting isolate already carries
                        # changes no field, and is still the decision that takes it
                        # out of the review tree. Accepting it is that decision.
                        confirm_organism(self.project, [sample_id], genus, species,
                                         scheme_path=assignment.get("scheme_path"),
                                         typing_mode=mode)
                except ValueError as error:
                    report["skipped"].append((sample_id, str(error)))
                    continue
                # Each sample keeps its own storage location; correcting a label
                # never moves a managed copy into a different root.
                outcome = refile_samples(self.project, [sample_id], cancelled=cancelled)
                for key in report:
                    report[key].extend(outcome[key])
            return report

        started = self.launch_task(operation, "refile", self.refile_completed)
        if started and notice:
            self.notify(notice)
        return started

    def assign_selected(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if not samples:
            self.notify("Select one or more sample rows first.")
            return
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        dialog = BatchAssignmentDialog(samples, self.scheme_entries(CLASSICAL_KINDS), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.apply_organism_assignments(dialog.assignments)

    def refile_selected(self):
        return self.context_refile(self.library_selection())

    def context_refile(self, selection):
        """Show where each managed copy would go, then move only what the user confirms."""
        if self.busy():
            return
        from PySide6.QtWidgets import QMessageBox

        from wmlstudio.storage import plan_filing, refile_samples
        identifiers = list(selection.sample_ids)
        if not identifiers:
            self.notify("Select the isolates whose managed copies should be re-filed.")
            return
        plans = [plan_filing(self.project, sample_id) for sample_id in identifiers]
        moving = [plan for plan in plans if plan["eligible"] and plan["changed"]]
        if not moving:
            reasons = {plan["reason"] for plan in plans if plan["reason"]}
            self.notify("Every selected managed copy is already in the folder its organism says. "
                        + (" ".join(sorted(reasons)) if reasons else ""))
            return
        lines = [f"{plan['name']}  →  {plan['relative']}" for plan in moving[:10]]
        if len(moving) > 10:
            lines.append(f"…and {len(moving) - 10} more.")
        if QMessageBox.question(self, "Re-file managed copies?",
                f"{len(moving)} managed copies will move to match the organism recorded for them. "
                "Your original files are not moved, and a folder name is a filing decision rather "
                "than a laboratory identification.\n\n" + "\n".join(lines)) != QMessageBox.StandardButton.Yes:
            return
        ids = [plan["sample_id"] for plan in moving]

        def operation(cancelled, progress):
            return refile_samples(self.project, ids, cancelled=cancelled, progress=progress)

        self.launch_task(operation, "refile", self.refile_completed)

    def download_practice_cohort(self):
        """Fetch a pinned teaching cohort after its caveats have been read."""
        if self.busy():
            return
        from wmlstudio import paths, practice_cohorts
        from wmlstudio.workflow_dialogs import PracticeCohortDialog
        cohorts = practice_cohorts.describe_cohorts()
        destinations = {entry["name"]: practice_cohorts.default_destination(paths.data_root(), entry["name"])
                        for entry in cohorts}
        dialog = PracticeCohortDialog(cohorts, destinations, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        name = dialog.chosen
        destination = destinations[name]

        def operation(cancelled, progress):
            return practice_cohorts.download_cohort(name, destination, cancelled=cancelled,
                                                    progress=progress)

        self._practice_cohort = None
        self.launch_task(operation, "practice_cohort",
                         lambda result: setattr(self, "_practice_cohort", result))

    def practice_cohort_ready(self):
        """Point the user at the downloaded files and start the ordinary import."""
        result, self._practice_cohort = self._practice_cohort, None
        if not result or not self._task_succeeded:
            return
        path = Path(result["path"])
        self.notify(f"{result['genomes']} practice genomes are ready in {path}. Every file was "
                    "checked against the checksum NCBI publishes for it.")
        self.import_paths([str(path)], configure=True)

    def install_species_panel(self):
        """Download the broader reference panel identification needs, on request only."""
        if self.busy():
            return
        from PySide6.QtWidgets import QMessageBox

        from wmlstudio import organism_panel, paths
        data_root = paths.data_root()
        installed = organism_panel.installed_species_panel(data_root)
        megabytes = sum(entry[5] for entry in organism_panel.SPECIES_PANEL) / (1024 * 1024)
        question = (f"Download {len(organism_panel.SPECIES_PANEL)} reference genomes "
                    f"({megabytes:.1f} MB) from NCBI RefSeq to {data_root}?\n\n"
                    "One reference per organism is a triage panel, not a representation of "
                    "within-species diversity. Nothing is uploaded and no genome leaves this "
                    "computer.")
        if installed:
            question = f"A panel is already installed at {installed}.\n\n" + question
        if QMessageBox.question(self, "Install broader species panel", question) != QMessageBox.StandardButton.Yes:
            return
        root = data_root / organism_panel.PANEL_DIRECTORY

        def operation(cancelled, progress):
            return organism_panel.provision_species_panel(root, cancelled=cancelled, progress=progress)

        from wmlstudio.analysis_progress import installing

        def finished(result):
            self.notify(
                f"Species panel installed: {result['species_count']} references at "
                f"{result['path']}. New imports will be compared against it. Use 'Identify "
                "waiting samples again' on the Assembly tab to re-run identification on the "
                "isolates already waiting.")
            self.refresh_update_center()

        self.launch_task(operation, "species_panel", finished,
                         caption=installing("the species reference panel"))

    # --- one step at a time -------------------------------------------------

    def run_step(self, key):
        """Run exactly this step, on the isolates in the selection that are ready for it."""
        if self.busy():
            return
        selected = self.selected_samples()
        if not selected:
            self.notify(f"Select the isolates to run {STEP_TITLES[key]} on first.")
            return
        ready, blocked = self.step_partition(key, selected)
        if not ready:
            _, sentence = self.step_gate_sentence(key)
            self.notify(sentence)
            return
        if blocked:
            self.notify(f"{len(blocked)} selected isolates are not ready for {STEP_TITLES[key]} "
                        f"and were left out. {blocked[0][0]['name']} needs "
                        f"{blocked[0][1]['missing'][0][0]}.")
        if key == "assembly":
            self.start_analysis(confirm=True, assemble=True,
                                sample_ids=[sample["id"] for sample in ready])
        elif key == "hydra":
            self.run_hydra_step(ready)
        else:
            self.run_typing_step(key, ready)

    def choose_step_scheme(self, key, organism):
        """Settle which reference this run uses, asking only when it is a real choice."""
        candidates = self.preferred_schemes(key, organism)
        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0]["path"]
        names = [f"{entry['name']} · {entry['locus_count']} loci" for entry in candidates]
        title = "seven-locus MLST" if key == "st" else "cgMLST"
        chosen, accepted = QInputDialog.getItem(
            self, f"Which {title} scheme?",
            f"{len(candidates)} installed schemes are labelled {organism_text(organism)}. "
            "The choice is recorded with every result it produces:", names, 0, False)
        return candidates[names.index(chosen)]["path"] if accepted and chosen in names else None

    def run_typing_step(self, key, samples):
        """Type the ready isolates against a scheme of this step's kind, and only that kind.

        A seven-locus ST and a core-genome profile are different measurements; the
        scheme is chosen for the step the user asked for, never inherited from the
        other tab's run.
        """
        from wmlstudio.sample_workflow import sample_configuration
        wanted = "mlst" if key == "st" else "cgmlst"
        prepared, chosen = [], {}
        for sample in samples:
            workflow = (sample.get("metadata") or {}).get("workflow") or {}
            configuration = sample_configuration(sample)
            # The sample's own configuration answers first, and each step reads its
            # own field: the cgMLST target set the user chose in any menu is what
            # the cgMLST tab runs, and it is never the classical workflow scheme.
            pinned = self.installed_scheme_path(
                (configuration["cgmlst_scheme_path"] if key == "cgmlst" else None)
                or configuration["scheme_path"])
            entry = self.scheme_organisms().get(str(Path(pinned).name).casefold()) if pinned else None
            if entry is not None and entry.get("kind") == wanted:
                path = pinned
            else:
                path = self.choose_step_scheme(key, accepted_organism(sample))
                if path is None:
                    self.notify("No scheme was chosen; nothing was run.")
                    return
            chosen[Path(path).name] = chosen.get(Path(path).name, 0) + 1
            prepared.append(dict(sample, metadata={**(sample.get("metadata") or {}), "workflow": {
                **workflow, "typing_mode": "manual", "scheme_path": str(path)}}))
        plan = self.review_run_plan(prepared)
        if plan is None:
            return
        if key == "cgmlst":
            # This tab is the cgMLST run. Chaining the launch flag on top of it
            # would call the same isolates against the same targets twice.
            plan = {**plan, "cgmlst": False}
        # Which reference produced a result is part of the result: say it before the
        # run as well, so nobody has to open a record to find out what was used.
        self.notify(f"{STEP_TITLES[key]} on {len(prepared)} isolates using "
                    + ", ".join(f"{name} ({count})" for name, count in sorted(chosen.items())) + ".")
        self._run_plan = plan
        self._run_ids = {sample["id"] for sample in prepared}
        self._run_cancelled = False
        self._typing_override = None
        toggle = self.step_toggles.get(key)
        if toggle is not None and not toggle.isChecked():
            self.check_typing_currency(key, prepared)
            return
        self.begin_typing(prepared, None, remember=False)

    # --- what a stored result already answers --------------------------------
    # "Already done, nothing changed" is only honest when the two things an allele
    # call is made against — the sequence input and the scheme's own fingerprint —
    # are both still what the stored result was produced from. The fingerprint is
    # the hash of the scheme's files, so answering the question means reading them;
    # that happens in a worker, through the same cache the run itself uses, so it
    # costs the run nothing.

    def check_typing_currency(self, key, prepared):
        """Ask which of these isolates a stored result already answers, then run the rest."""
        from wmlstudio.identification import cached_scheme
        from wmlstudio.sample_workflow import current_input_sha256, rerun_decision
        kind = "mlst" if key == "st" else "cgmlst"
        project = self.project
        samples = list(prepared)

        def operation(cancelled, progress):
            digests, decisions = {}, {}
            for index, sample in enumerate(samples):
                progress(index, len(samples), f"Checking what already stands for {sample['name']}…")
                path = ((sample.get("metadata") or {}).get("workflow") or {}).get("scheme_path") or ""
                if path not in digests:
                    digests[path] = cached_scheme(path, cancelled).digest
                try:
                    stored = project.latest_analysis(sample["id"], kind)
                except (KeyError, ValueError):
                    stored = None
                decisions[sample["id"]] = rerun_decision(
                    stored, {"input_sha256": current_input_sha256(sample),
                             "scheme_digest": digests[path]})
            return decisions

        self._pending_typing = {"key": key, "samples": samples}
        self._typing_currency = None
        if not self.launch_task(operation, "typing_currency",
                                lambda result: setattr(self, "_typing_currency", result)):
            self._pending_typing = None
            self._run_plan = {}

    def continue_after_currency(self):
        """Run only what a stored result does not already answer, naming what changed."""
        pending, self._pending_typing = self._pending_typing, None
        decisions, self._typing_currency = self._typing_currency or {}, None
        if not pending or self._run_cancelled or not self._task_succeeded:
            self._run_plan = {}
            return
        key, samples = pending["key"], pending["samples"]
        for sample_id, decision in decisions.items():
            self._typing_verdicts[(key, sample_id)] = decision
        settled = {s["id"] for s in samples
                   if (decisions.get(s["id"]) or {}).get("status") == "current"}
        standing = [s for s in samples if s["id"] in settled]
        changed = [s for s in samples if s["id"] not in settled]
        title = STEP_TITLES[key]
        if not changed:
            plan, self._run_plan = self._run_plan, {}
            self.last_typing_report = (
                f"Nothing was run. All {len(standing)} selected isolates already have a {title} "
                "result and nothing it was made against has changed — same sequence input, same "
                f"scheme fingerprint. Tick '{RERUN_TOGGLE}' to produce it again.")
            self.refresh()
            self.notify(self.last_typing_report)
            # The stored typing standing does not cancel the rest of the launch: a
            # cgMLST or HYDRA run this plan also asked for is a different question.
            self.start_reviewed_plan(samples, plan)
            return
        reasons = "; ".join(f"{sample['name']} — "
                            + ((decisions.get(sample["id"]) or {}).get("summary") or "never typed")
                            for sample in changed[:3])
        if len(changed) > 3:
            reasons += f"; and {len(changed) - 3} more."
        message = f"Running {title} on {len(changed)} of {len(samples)} isolates. {reasons}"
        if standing:
            message += (f" {len(standing)} were left alone: their stored {title} result stands, so "
                        "repeating it would produce the same answer.")
        self.last_typing_report = message
        self.notify(message)
        self._run_ids = {sample["id"] for sample in changed}
        self.begin_typing(changed, None, remember=False)

    # --- HYDRA organism -----------------------------------------------------

    def build_hydra_organism_control(self):
        """Choose the organism HYDRA is told about, and say what that unlocks.

        The organism-specific tools — AMRFinderPlus organism flags, Kleborate and
        SCCmec typing — only run when an organism is set, and the only values that
        mean anything are the ones the installed point-mutation catalogues cover.
        Choosing one here is your statement about the isolate, not new evidence.
        """
        panel, content = card()
        content.addWidget(label("Organism for the organism-specific tools", "cardTitle"))
        row = FlowLayout()
        self.hydra_organism = QComboBox()
        self.hydra_organism.setAccessibleName("HYDRA organism")
        self.hydra_organism.setMinimumWidth(220)
        row.addWidget(self.hydra_organism)
        row.addWidget(button("Apply to selected samples",
                             lambda: self.apply_hydra_organism("selected")))
        row.addWidget(button("Apply to every sample in this project",
                             lambda: self.apply_hydra_organism("project")))
        content.addLayout(row)
        self.hydra_organism_note = label("", "small", True)
        content.addWidget(self.hydra_organism_note)
        return panel

    def populate_hydra_organisms(self):
        """Offer only organisms the installed database actually has catalogues for."""
        if not hasattr(self, "hydra_organism"):
            return
        state = self.amr_database_state()
        names = list(state["organisms"])
        for sample in active_samples(self.current_samples):
            organism = accepted_organism(sample)
            if organism and organism_text(organism) not in names:
                names.append(organism_text(organism))
        current = self.hydra_organism.currentData()
        self.hydra_organism.blockSignals(True)
        self.hydra_organism.clear()
        self.hydra_organism.addItem("Unresolved · no organism flags", "")
        for name in sorted(dict.fromkeys(names)):
            self.hydra_organism.addItem(name, name)
        index = self.hydra_organism.findData(current)
        self.hydra_organism.setCurrentIndex(max(0, index))
        self.hydra_organism.blockSignals(False)
        if state["organisms"]:
            self.hydra_organism_note.setText(
                f"{len(state['organisms'])} organisms have point-mutation catalogues in the "
                "installed database. An organism is what enables the AMRFinderPlus organism "
                "rules, Kleborate and SCCmec typing; setting one here records it as your own "
                "assignment — it is not genomic evidence, and it never changes an organism the "
                "genome comparison supported. A detected gene is not measured susceptibility.")
        else:
            self.hydra_organism_note.setText(
                "The installed AMR database reports no organism catalogues, so the AMRFinderPlus "
                "organism rules are unavailable and no point mutation can be called. "
                + (state["reason"] or ""))

    def hydra_organism_for(self, sample):
        """The organism HYDRA would be told about for this isolate, or an empty string."""
        override = str(((sample.get("metadata") or {}).get("workflow") or {}).get("hydra_organism") or "")
        if override:
            return override
        organism = accepted_organism(sample)
        # HYDRA's organism flags need a species; a genus alone cannot select a catalogue.
        return organism_text(organism) if organism and organism[1] else ""

    def apply_hydra_organism(self, scope):
        """Record the chosen organism for HYDRA, on the selection or the whole project."""
        if self.busy():
            return
        name = self.hydra_organism.currentData() or ""
        samples = (self.selected_samples() if scope == "selected"
                   else active_samples(self.current_samples))
        if not samples:
            self.notify("Select the isolates this organism applies to first.")
            return
        for sample in samples:
            self.project.update_metadata(sample["id"], {"workflow": {"hydra_organism": name}})
        self.refresh()
        if name:
            self.notify(f"{len(samples)} isolates will be screened as {name}. This is your "
                        "assignment for the organism-specific tools; the organism recorded from "
                        "genomic evidence is unchanged.")
        else:
            self.notify(f"{len(samples)} isolates will be screened with no organism flags. "
                        "AMRFinderPlus organism rules, Kleborate and SCCmec typing stay "
                        "unavailable for them.")

    def run_hydra_step(self, samples):
        return self.start_reviewed_plan(samples, self.review_run_plan(samples, hydra=True))

    def start_reviewed_plan(self, samples, plan):
        """Start what one HYDRA-oriented launch asked for, in order, or nothing.

        cgMLST runs first when the plan carries that flag, because HYDRA screens
        whatever assemblies the launch ends with; the chain is then driven from
        analysis_finished exactly as an ordinary run is. Returns False when the plan
        asked for nothing, so the caller never claims a run that did not start.
        """
        if not plan:
            return False
        identifiers = {sample["id"] for sample in samples}
        self._run_cancelled = False
        self._run_ids = set(identifiers)
        if plan.get("cgmlst"):
            self._run_plan = {**plan, "cgmlst": False}
            if self.run_cgmlst_plan(identifiers, plan):
                return True
        self._run_plan = {}
        if plan.get("hydra"):
            self.run_hydra_plan(identifiers, plan, require_completed=False)
            return True
        return False

    def busy(self):
        if self.worker and self.worker.isRunning():
            self.notify("A background task is active. You can navigate and inspect results; cancel or finish before changing inputs.")
            return True
        dialog = self.amr_database_dialog
        if dialog and dialog.worker and dialog.worker.isRunning():
            self.notify("A reference snapshot is being downloaded. Finish or cancel that task before starting another analysis.")
            return True
        return False

    def launch_task(self, operation, role, completed=None, *, caption=None):
        """Run one background task, with the progress dialog captioned by its caller.

        ``caption`` is an analysis_progress.TaskCaption naming what is running.
        Without one the dialog says an analysis is running, which is true of most
        callers and was badly untrue of the reference-data downloads: a person who
        pressed Install was told "your analysis is running" and reported that the
        update button runs analyses.
        """
        if self.busy():
            return False
        self.worker_role = role
        self._task_succeeded = False
        self._task_caption = caption
        self.worker = FunctionWorker(operation, self)
        self.worker.completed.connect(lambda result: setattr(self, "_task_succeeded", True))
        if completed:
            self.worker.completed.connect(completed)
        self.worker.failed.connect(self.error)
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.set_running(True)
        self.worker.start()
        return True

    def set_running(self, running):
        for widget in (self.run_button, self.rerun_button, self.scheme_combo, self.demo_button):
            widget.setEnabled(not running)
        for action in self.step_buttons.values():
            # Re-enabled by the gate, not by the end of the task: a step whose facts
            # are still missing must not become clickable just because nothing is busy.
            action.setEnabled(False)
        if not running:
            self.update_step_gates()
        self.cancel_button.setVisible(running)
        self.progress_bar.setVisible(running)
        if running:
            self.progress_bar.setValue(0)
            from wmlstudio.analysis_progress import AnalysisProgressDialog
            caption = getattr(self, "_task_caption", None)
            if self._progress_dialog is None or self._progress_dialog.finished_safely:
                self._progress_dialog = AnalysisProgressDialog(self, caption)
            else:
                # A dialog being shown again must be re-labelled: it carries the
                # previous task's words until it is told what this one is.
                self._progress_dialog.set_caption(caption)
            self._progress_dialog.show()
        elif self._progress_dialog is not None:
            self._progress_dialog.finish()

    def job_progress(self, percent, message):
        super().job_progress(percent, message)
        if self._progress_dialog is not None and not self._progress_dialog.finished_safely:
            self._progress_dialog.update_progress(percent, message)

    def start_analysis(self, checked=False, all_samples=False, selected_only=False, confirm=False, assemble=False, sample_ids=None):
        if self.busy():
            return
        chosen_ids = set(sample_ids) if sample_ids is not None else None
        samples = [s for s in self.project.samples() if
                   (s["id"] in chosen_ids if chosen_ids is not None else
                    s["id"] in self.selection_ids if selected_only else
                    all_samples or s["status"] in {"queued", "failed", "interrupted"})]
        samples = [s for s in samples if s.get("input_path") and s.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
        if not samples:
            self.notify("No sequence inputs in the requested selection. Import files, select rows or use saved profiles in Compare.")
            return
        scheme = self.scheme_combo.currentData()
        if scheme:
            samples = [dict(s, metadata={**s.get("metadata", {}), "workflow": {
                **s.get("metadata", {}).get("workflow", {}), "typing_mode": "manual", "scheme_path": scheme}})
                for s in samples]
        plan = self.review_run_plan(samples, scheme, assemble=assemble) if confirm else {}
        if plan is None:
            return
        self._run_plan = plan
        self._run_ids = {sample["id"] for sample in samples}
        self._run_cancelled = False
        self._typing_override = scheme
        pairs = None
        if plan.get("assemble"):
            from wmlstudio.pairing_dialog import PairReadsDialog
            reads = [sample for sample in samples if Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fq", ".fastq"))]
            if len(reads) < 2:
                self.notify("Select at least two paired FASTQ inputs to assemble, or disable assembly and run QC only.")
                self._run_plan = {}
                return
            pairing = PairReadsDialog(reads, self)
            if pairing.exec() != QDialog.DialogCode.Accepted:
                self._run_plan = {}
                return
            pairs = pairing.assignments
        if plan.get("fastqc"):
            self.begin_fastqc_plan(samples, scheme, pairs)
        elif pairs:
            self.assemble_pairs(pairs, plan)
        else:
            self.begin_typing(samples, scheme)

    def begin_fastqc_plan(self, samples, scheme, pairs=None):
        from wmlstudio.fastqc_dialog import reads_for_sample, run_project_fastqc
        from wmlstudio.scheduler import resources_for_run
        reads = [sample for sample in samples if reads_for_sample(sample)]
        if not reads:
            self.notify("No original reads are linked to the reviewed inputs; FastQC was not run.")
            self._run_plan = {}
            return
        try:
            allocation = resources_for_run(self._run_plan, memory_gb=1)
        except ValueError as exc:
            self.error(str(exc))
            self._run_plan = {}
            return
        self._fastqc_pending = {"ids": [sample["id"] for sample in samples], "scheme": scheme, "pairs": pairs}
        self._fastqc_flagged = 0
        project = self.project
        output = self.project_path.with_suffix(".files") / "fastqc_reports"
        def operation(cancelled, progress):
            return run_project_fastqc(project, reads, output, allocation, cancelled=cancelled, progress=progress)
        def completed(results):
            self._fastqc_flagged = sum(report["qc_status"] == "FAIL" for entry in results
                                      for report in entry["result"]["reports"])
            self.notify(f"Original FastQC completed for {len(results)} inputs; {self._fastqc_flagged} read reports have FAIL flags. No reads were trimmed.")
        self.launch_task(operation, "fastqc_pipeline", completed)

    def continue_after_fastqc(self):
        pending = getattr(self, "_fastqc_pending", None)
        self._fastqc_pending = None
        if not pending or self._run_cancelled or not self._task_succeeded:
            self._run_plan = {}
            return
        if self._fastqc_flagged:
            from PySide6.QtWidgets import QMessageBox
            answer = QMessageBox.question(self, "Review FastQC flags before continuing",
                f"{self._fastqc_flagged} read reports contain FastQC FAIL flags. Reports are saved with the isolates. "
                "These flags are not an automatic instruction to trim, but may affect downstream interpretation. "
                "Continue the reviewed analysis on the unchanged reads?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                self._run_plan = {}
                self.notify("Stopped after FastQC for review. Saved reports and original reads are retained.")
                return
        if pending["pairs"]:
            self.assemble_pairs(pending["pairs"], self._run_plan)
            return
        samples = [self.project.get_sample(identifier) for identifier in pending["ids"]]
        scheme = pending["scheme"]
        self.begin_typing(samples, scheme)

    def begin_typing(self, samples, scheme, *, remember=True):
        """Type these samples. `remember` records the global scheme the user picked.

        A step run chooses a scheme per isolate for that step alone, so it passes
        remember=False: it must not overwrite the workspace-wide default.
        """
        from wmlstudio.scheduler import resources_for_run
        try:
            allocation = resources_for_run(self._run_plan)
        except ValueError as exc:
            self.error(str(exc))
            return
        if scheme:
            samples = [dict(sample, metadata={**sample.get("metadata", {}), "workflow": {
                **sample.get("metadata", {}).get("workflow", {}), "typing_mode": "manual", "scheme_path": scheme}})
                for sample in samples]
        if remember:
            self.project.set_setting("scheme_path", scheme or "")
        self.worker_role = "analysis"
        self.worker = AnalysisWorker(samples, scheme, parent=self, installed_scheme_paths=self.scheme_paths,
                                     resource_plan=allocation)
        self.worker.sample_started.connect(self.sample_started)
        self.worker.sample_finished.connect(self.sample_finished)
        self.worker.sample_failed.connect(self.sample_failed)
        self.worker.sample_cancelled.connect(lambda sid: self.project.set_status(sid, "interrupted", "Cancelled by user"))
        self.worker.progress.connect(self.job_progress)
        self.worker.finished.connect(self.analysis_finished)
        self.set_running(True)
        self.worker.start()

    def assemble_pairs(self, pairs, plan):
        import uuid

        from wmlstudio.assembly import associate_assembly, run_skesa
        from wmlstudio.scheduler import resources_for_run, run_bounded
        project = self.project
        try:
            allocation = resources_for_run(plan, memory_gb=8)
        except ValueError as exc:
            self.error(str(exc))
            return
        tasks = [{"primary": project.get_sample(pair["primary_id"]),
                  "mate": project.get_sample(pair["mate_id"])} for pair in pairs]
        output_root = self.project_path.with_suffix(".files") / "assemblies"
        def operation(cancelled, progress):
            identifiers = []
            def assemble(task, resources, stopped, report):
                primary, mate = task["primary"], task["mate"]
                destination = output_root / primary["id"] / uuid.uuid4().hex
                return run_skesa(primary["input_path"], mate["input_path"], destination,
                    threads=resources.threads_per_sample, memory_gb=resources.memory_gb, cancelled=stopped,
                    progress=lambda done, total, message: report(done, total, f"{primary['name']} · {message}"))
            def attach(task, result):
                primary, mate = task["primary"], task["mate"]
                associate_assembly(project, primary["id"], mate["id"], result)
                identifiers.append(primary["id"])
            run_bounded(tasks, assemble, allocation, cancelled=cancelled, on_result=attach, progress=progress)
            return identifiers
        self.launch_task(operation, "assembly", lambda ids: self.notify(f"Assembled {len(ids)} read pairs. Starting the assigned typing workflow…"))

    def sample_finished(self, sample_id, result):
        """Save the result, then let the organism it established fill itself in everywhere.

        This is the user's "the info for the first analysis should autofill
        everywhere": the moment an analysis establishes an organism it lands on the
        sample's one configuration, so no later menu asks again. It never overrules
        a person — their label is kept and the proposal is held — and a scheme match
        is adopted as panel compatibility, never as an independent identification.
        """
        self.project.set_result(sample_id, result)
        from wmlstudio.storage import adopt_analysis_organism
        try:
            adopt_analysis_organism(self.project, sample_id, result)
        except (KeyError, ValueError):
            pass  # Recording a label must never lose the result that produced it.
        self.refresh()

    def analysis_finished(self):
        role = self.worker_role
        self.worker_role = ""
        self.set_running(False)
        self.refresh()
        if self.closing_after_cancel:
            self.close()
            return
        if role in INSTALL_ROLES or role.startswith("update:"):
            # Whatever installed something has stopped — finished, failed or
            # cancelled — so the page that lists what is installed re-reads this
            # computer rather than leaving a row that a person has just filled
            # still reading "Not installed" beside an "Install…" button. After the
            # close check: a window on its way out starts no new probe.
            self.refresh_update_center()
        if role == "fastqc_pipeline":
            self.continue_after_fastqc()
            return
        if role == "identify":
            self.review_identification()
            return
        if role == "intake":
            self.intake_completed()
            return
        if role == "typing_currency":
            self.continue_after_currency()
            return
        if role == "scheme_import":
            self.dropped_scheme_installed()
            return
        if role == "practice_cohort":
            self.practice_cohort_ready()
            return
        if role == "assembly":
            if self._task_succeeded and not self._run_cancelled:
                samples = [sample for sample in self.project.samples() if sample["id"] in self._run_ids
                           and sample.get("metadata", {}).get("workflow", {}).get("source_kind") != "read_mate"]
                self.begin_typing(samples, self._typing_override)
            else:
                self._run_plan = {}
            return
        if role == "analysis":
            # An isolate that has just been typed again is no longer answered by the
            # verdict that sent it into this run, so that verdict is dropped.
            self._typing_verdicts = {key: decision for key, decision in self._typing_verdicts.items()
                                     if key[1] not in self._run_ids}
        if role == "analysis" and not self._run_cancelled:
            managed = [s for s in self.project.samples() if s["id"] in self._run_ids and s["status"] == "completed"
                       and s.get("metadata", {}).get("workflow", {}).get("managed")]
            if managed:
                from wmlstudio.storage import refile_samples
                identifiers = [sample["id"] for sample in managed]

                def organize(cancelled, progress):
                    return refile_samples(self.project, identifiers, cancelled=cancelled,
                                          progress=progress)

                self.launch_task(organize, "organize", self.refile_completed)
                return
        if role in {"analysis", "organize"}:
            usable = not self._run_cancelled and (role != "organize" or self._task_succeeded)
            # cgMLST first, because it is a typing run this launch also asked for and
            # HYDRA screens whatever assemblies the launch ends with. The flag is
            # taken off the plan before launching so the chain runs once.
            if usable and self._run_plan.get("cgmlst"):
                plan = dict(self._run_plan)
                self._run_plan = {**plan, "cgmlst": False}
                if self.run_cgmlst_plan(self._run_ids, plan):
                    return
            if self._run_plan.get("hydra") and usable:
                self.run_hydra_plan(self._run_ids, self._run_plan)
            self._run_plan = {}

    def cancel_analysis(self):
        self._run_cancelled = True
        self._run_plan = {}
        super().cancel_analysis()

    def active_amr_database(self):
        selected = self.project.get_setting("hydra_database_root", "")
        if selected:
            return selected
        from wmlstudio import hydra_runtime
        bundled = getattr(hydra_runtime, "bundled_database_root", lambda: None)()
        return str(bundled or self.root / "references" / "hydra")

    def review_run_plan(self, samples, scheme=None, hydra=False, assemble=False, cgmlst=False):
        """The one launch review: what runs, on which isolates, against which references.

        cgMLST sits here beside HYDRA because the user asked for it in the same
        place: one launch can produce a seven-locus ST and a core-genome profile,
        and they stay two separate measurements against two separate schemes. The
        dialog refuses a cgMLST scheme the cohort's organism contradicts, naming
        the organism the target set was defined on.
        """
        from wmlstudio.analysis_plan import RunPlanDialog
        dialog = RunPlanDialog(samples, scheme, self.active_amr_database(), self,
                               data_root=self.root)
        dialog.hydra.setChecked(hydra)
        dialog.assemble.setChecked(assemble)
        if cgmlst and dialog.cgmlst.isEnabled():
            dialog.cgmlst.setChecked(True)
        if hydra:
            dialog.assemble.setEnabled(False)
        dialog.manageDatabases.connect(lambda: self.manage_run_databases(dialog))
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        if dialog.plan.get("hydra"):
            self.project.set_setting("hydra_database_root", dialog.plan["db_root"])
        self.remember_run_flags(samples, dialog.plan)
        return dialog.plan

    def remember_run_flags(self, samples, plan):
        """Keep the run flags and the chosen cgMLST scheme on the samples themselves.

        These are configuration, not a one-off: a cohort the user decided should be
        screened for AMR and called against one cgMLST target set arrives at the
        next launch already saying so, instead of being asked again. They are
        workflow-only fields, so recording them never invalidates a stored result.
        """
        from wmlstudio.storage import set_sample_configuration
        chosen = plan.get("cgmlst_scheme") or {}
        changes = {"run_hydra": bool(plan.get("hydra")), "run_cgmlst": bool(plan.get("cgmlst")),
                   "cgmlst_scheme_key": chosen.get("key") or None,
                   "cgmlst_scheme_path": chosen.get("path") or None}
        identifiers = [sample["id"] for sample in samples if sample.get("id")]
        if not identifiers:
            return {}
        try:
            return set_sample_configuration(self.project, identifiers, changes,
                                            surface="Run plan")
        except (KeyError, ValueError) as error:
            # Remembering a preference must never stop a run the user just approved.
            self.notify(f"The run starts, but this choice could not be saved on the samples: {error}")
            return {}

    def manage_run_databases(self, plan_dialog):
        from wmlstudio.amr_databases import AMRDatabaseDialog
        database_dialog = AMRDatabaseDialog(self.active_amr_database(), plan_dialog,
            update_root=self.root / "references" / "hydra")
        self.amr_database_dialog = database_dialog
        database_dialog.snapshotInstalled.connect(plan_dialog.database.setText)
        database_dialog.snapshotInstalled.connect(lambda path: self.project.set_setting("hydra_database_root", path))
        database_dialog.exec()

    def run_hydra_selected(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if not samples:
            self.notify("Select assembly samples in the Samples menu first.")
            return
        self.start_reviewed_plan(samples, self.review_run_plan(samples, hydra=True))

    def run_cgmlst_plan(self, identifiers, plan):
        """Call cgMLST on this run's isolates, against the target set the plan named.

        This is the second half of the run-cgMLST flag the user asked for beside
        run-HYDRA: one launch, two measurements. The profile is stored under its own
        typing kind against its own scheme, so a core-genome distance and the
        seven-locus ST this run also produced never share a scale or a threshold.
        Returns False when there was nothing it could honestly call.
        """
        chosen = plan.get("cgmlst_scheme") or {}
        path = self.installed_scheme_path(chosen.get("path") or "")
        if not path:
            self.notify("No installed cgMLST scheme was chosen, so no cgMLST profile was called.")
            return False
        samples = [sample for sample in self.project.samples()
                   if sample["id"] in identifiers and sample.get("input_path")
                   and input_kind(sample) == "assembly" and not sample.get("missing_input")]
        if not samples:
            self.notify("No assembled isolate in this run could be called against the cgMLST "
                        "scheme. Reads are never typed directly; assemble them first.")
            return False
        prepared = [dict(sample, metadata={**(sample.get("metadata") or {}), "workflow": {
            **((sample.get("metadata") or {}).get("workflow") or {}),
            "typing_mode": "manual", "scheme_path": path}}) for sample in samples]
        targets = chosen.get("locus_count")
        self._run_ids = {sample["id"] for sample in prepared}
        self.notify(f"cgMLST on {len(prepared)} isolates against "
                    f"{chosen.get('scheme_name') or Path(path).name}"
                    + (f" · {targets} targets" if targets else "")
                    + ". A cgMLST distance and a seven-locus distance are different quantities "
                      "against different schemes.")
        self.begin_typing(prepared, None, remember=False)
        return True

    def run_hydra_plan(self, identifiers, plan, require_completed=True):
        from wmlstudio.hydra_runtime import run_assemblies
        from wmlstudio.sample_workflow import link_hydra
        from wmlstudio.scheduler import resources_for_run, run_bounded
        samples = [sample for sample in self.project.samples() if sample["id"] in identifiers
                   and sample.get("input_path") and (not require_completed or sample["status"] == "completed")]
        assemblies = [sample for sample in samples if (sample.get("result") or {}).get("kind", sample.get("kind")) == "fasta"
                      or Path(sample["input_path"]).name.lower().removesuffix(".gz").removesuffix(".bz2").endswith((".fa", ".fasta", ".fna"))]
        if not assemblies:
            self.notify("No FASTA assemblies are ready for HYDRA in this selection. Raw reads must be assembled first.")
            return
        try:
            allocation = resources_for_run(plan, memory_gb=3)
        except ValueError as exc:
            self.error(str(exc))
            return

        def operation(cancelled, progress):
            import hashlib

            reports = []
            def analyse(sample, resources, stopped, report_progress):
                # A provisional MLST lineage alone is not a verified mutation-catalog
                # assignment: only an accepted organism, or one the user chose on the
                # HYDRA tab, selects an organism catalogue.
                organism = self.hydra_organism_for(sample) or None
                return run_assemblies([sample["input_path"]], plan["db_root"], plan.get("databases"),
                    sample_names=[sample["id"]], organism=organism, threads=resources.threads_per_sample,
                    protein=plan.get("protein", True), cancelled=stopped,
                    point_mutations=plan.get("point_mutations", True),
                    **plan.get("thresholds", {}),
                    progress=lambda done, total, message: report_progress(done, total, f"{sample['name']} · {message}"))
            def attach(sample, report):
                link_hydra(self.project, report, {sample["id"]: sample["id"]})
                reports.append(report)
                combined = dict(report, samples=[entry for result in reports for entry in result["samples"]])
                components = [result["import_provenance"] for result in reports]
                combined["import_provenance"] = {"source_path": "Native HYDRA run", "components": components,
                    "sha256": hashlib.sha256(json.dumps(components, sort_keys=True).encode()).hexdigest()}
                combined["execution_provenance"] = {"sample_runs": [result.get("execution_provenance", {}) for result in reports]}
                self.project.set_setting("hydra_report", combined)
            run_bounded(assemblies, analyse, allocation, cancelled=cancelled, on_result=attach, progress=progress)
            return len(reports)

        self.launch_task(operation, "hydra", lambda count: self.notify(f"HYDRA completed for {count} assemblies. AMR evidence is linked to the sample library and reports."))

    def open_amr_databases(self):
        if self.busy():
            return
        from wmlstudio.amr_databases import AMRDatabaseDialog
        if self.amr_database_dialog is None:
            self.amr_database_dialog = AMRDatabaseDialog(self.active_amr_database(), self,
                update_root=self.root / "references" / "hydra")
            self.amr_database_dialog.snapshotInstalled.connect(lambda path: self.project.set_setting("hydra_database_root", path))
        self.amr_database_dialog.show()
        self.amr_database_dialog.raise_()
        self.amr_database_dialog.activateWindow()

    def switch_project(self, path):
        if self.worker and self.worker.isRunning():
            return super().switch_project(path)
        old_path = self.project_path.resolve()
        if self.library is not None:
            self.library.index_project(self.project)
        changed = super().switch_project(path)
        if changed and old_path != self.project_path.resolve():
            self.selection_ids.clear()
            self.report_ids.clear()
            cohort = self.project.get_setting("comparison_cohort", None)
            self.cohort_ids = set(cohort) if cohort is not None else set()
            self.feature_ids = set()
            for dialog in getattr(self, "isolate_dialogs", {}).values():
                dialog.close()
            self.isolate_dialogs = {}
            self.library_filter = None
            if hasattr(self, "restore_investigations"):
                self.restore_investigations()
            if hasattr(self, "tree") and hasattr(self.tree, "restore_state"):
                self._pending_graph_state = self.project.get_setting("graph_style", {"version": 1})
            self.refresh()
        return changed

    def edit_metadata(self):
        samples = self.selected_samples()
        if not samples:
            sample = self.selected_sample()
            samples = [sample] if sample else []
        if not samples:
            self.notify("Select sample rows first.")
            return
        annotations = samples[0].get("metadata", {}).get("annotations", {}) if len(samples) == 1 else {}
        text, accepted = QInputDialog.getMultiLineText(self, "Sample annotations", "field = value, one per line. These annotations apply to every selected sample; workflow and AMR evidence are preserved.", "\n".join(f"{key} = {value}" for key, value in annotations.items()))
        if not accepted:
            return
        try:
            values = {}
            for line in text.splitlines():
                if line.strip():
                    key, separator, value = line.partition("=")
                    if not separator or not key.strip():
                        raise ValueError("Use field = value on each line.")
                    values[key.strip()] = value.strip()
            for sample in samples:
                metadata = dict(sample.get("metadata") or {})
                metadata["annotations"] = {**metadata.get("annotations", {}), **values}
                self.project.set_metadata(sample["id"], metadata)
            self.refresh()
        except Exception as exc:
            self.error(exc)

    def add_collection(self):
        if not self.selection_ids:
            self.notify("Select the samples to add to a collection first.")
            return
        name, accepted = QInputDialog.getText(self, "Add to collection", "Collection name (for example, Ward 5 · September):")
        if accepted and name.strip():
            collection_id = self.project.create_collection(name.strip())
            self.project.set_collection_members(collection_id, self.selection_ids, add=True)
            for sample in self.selected_samples():
                metadata = dict(sample.get("metadata") or {})
                metadata["collections"] = sorted(set(metadata.get("collections", []) + [name.strip()]))
                self.project.set_metadata(sample["id"], metadata)
            self.refresh()

    def open_data_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project_path.parent)))

    def open_library_research(self):
        if self.busy():
            return
        from wmlstudio.library_dialog import LibraryDialog
        self.library.index_project(self.project)
        dialog = LibraryDialog(self.library, self.project, self)
        dialog.profilesImported.connect(self.import_completed)
        dialog.exec()

    def open_sample_folder(self):
        sample = self.selected_sample()
        if sample and sample.get("input_path"):
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(sample["input_path"]).parent)))

    def relink_selected_input(self):
        if self.busy():
            return
        samples = self.selected_samples()
        if len(samples) != 1:
            self.notify("Select exactly one sample to relink. The replacement must match its recorded input SHA-256.")
            return
        sample = samples[0]
        path, _ = QFileDialog.getOpenFileName(self, f"Relink {sample['name']} · identical input bytes required", "",
            "Sequence files (*.fasta *.fa *.fna *.fastq *.fq *.gz *.bz2);;All files (*)")
        if not path:
            return
        from wmlstudio.storage import relink_input
        self.launch_task(lambda cancelled, progress: relink_input(self.project, sample["id"], path, cancelled=cancelled),
            "relink", lambda result: self.notify("Input relinked after SHA-256 verification. Saved profiles and original files are unchanged."))

    # --- right-click handlers -----------------------------------------------
    # One method per entry in context_menus.ACTIONS that acts on isolates. A view
    # is wired with a single install_view_menu call; an action whose handler is
    # absent is left out of the menu rather than offered and then failing.

    def context_records(self, selection):
        wanted = set(selection.sample_ids)
        return [sample for sample in self.project.samples() if sample["id"] in wanted]

    def context_open_record(self, selection):
        if selection.single:
            self.open_isolate_record(selection.single)

    def context_select_in_library(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.clear_filters()
        self.selection_ids = set(ids)
        self.library_filter = ("ids", ids)
        self.focus.set_focus(ids, "Right-click selection")
        self.refresh_tables()
        self.navigate(1)
        self.notify(f"{len(ids)} isolates shown. Clear filters to return to the whole library.")

    def update_comparison_cohort(self, ids, add=True):
        """Change the comparison cohort explicitly, and say where the change came from."""
        ids = {value for value in ids if value}
        if not ids:
            return
        self.cohort_ids = (set(self.cohort_ids) | ids) if add else (set(self.cohort_ids) - ids)
        self.project.set_setting("comparison_cohort", sorted(self.cohort_ids))
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("compare", self.cohort_ids, "Right-click in the isolate library")
        if hasattr(self, "cohort_table"):
            self.refresh_cohort_table()
        self.refresh_journey()
        verb = "added to" if add else "removed from"
        self.notify(f"{len(ids)} isolates {verb} the comparison cohort · {len(self.cohort_ids)} in it now. "
                    "Nothing else was included automatically.")

    def context_add_to_comparison(self, selection):
        self.update_comparison_cohort(selection.sample_ids, add=True)

    def context_remove_from_comparison(self, selection):
        self.update_comparison_cohort(selection.sample_ids, add=False)

    def context_add_to_characterization(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.feature_ids = set(self.feature_ids) | ids
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("evidence", self.feature_ids, "Right-click in the isolate library")
        if hasattr(self, "feature_table"):
            self.refresh_features()
        self.notify(f"{len(ids)} isolates added to the evidence cohort · {len(self.feature_ids)} in it now.")

    def context_add_to_report(self, selection):
        ids = set(selection.sample_ids)
        if not ids:
            return
        self.report_ids = set(self.report_ids) | ids
        ledger = getattr(self, "cohort_origins", None)
        if ledger is not None:
            ledger.record("reports", self.report_ids, "Right-click in the isolate library")
        if hasattr(self, "report_table"):
            self.refresh_report_table()
        self.notify(f"{len(ids)} isolates added to the report · {len(self.report_ids)} in it now.")

    def context_remove_from_report(self, selection):
        ids = set(selection.sample_ids)
        self.report_ids = set(self.report_ids) - ids
        if hasattr(self, "report_table"):
            self.refresh_report_table()
        self.notify(f"{len(ids)} isolates removed from the report · {len(self.report_ids)} in it now.")

    @staticmethod
    def label_assignment(sample, genus, species):
        """An organism correction that keeps the sample's own typing workflow."""
        workflow = (sample.get("metadata") or {}).get("workflow", {})
        mode = workflow.get("typing_mode") or "manual"
        if mode == "unknown" and genus:
            mode = "manual"
        return {"sample_id": sample["id"], "genus": genus, "species": species,
                "typing_mode": mode, "scheme_path": workflow.get("scheme_path")}

    def context_assign_organism(self, selection):
        if self.busy():
            return
        samples = self.context_records(selection)
        if not samples:
            return
        from wmlstudio.workflow_dialogs import BatchAssignmentDialog
        dialog = BatchAssignmentDialog(samples, self.scheme_entries(CLASSICAL_KINDS), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.apply_organism_assignments(dialog.assignments)

    def context_assign_organism_quick(self, selection, organism):
        if self.busy():
            return
        genus, species = organism
        samples = self.context_records(selection)
        self.apply_organism_assignments([self.label_assignment(sample, genus, species)
                                         for sample in samples])

    def context_assign_scheme(self, selection):
        if self.busy():
            return
        entries = self.scheme_entries(CLASSICAL_KINDS)
        samples = self.context_records(selection)
        if not entries or not samples:
            self.notify("No typing schemes are installed. Use Data ▸ Online scheme catalog first.")
            return
        names = [name for name, _ in entries]
        title, accepted = QInputDialog.getItem(self, "Assign typing scheme",
            f"Type these {len(samples)} isolates with:", names, 0, False)
        if not accepted:
            return
        path = entries[names.index(title)][1]
        assignments = []
        for sample in samples:
            genus, species, _ = organism_for(sample)
            if not genus:
                self.notify(f"{sample['name']} has no organism yet. Assign a genus first, or use "
                            "Unknown organism in the assignment dialog.")
                return
            assignments.append({"sample_id": sample["id"], "genus": genus, "species": species,
                                "typing_mode": "manual", "scheme_path": path})
        self.apply_organism_assignments(assignments)

    def context_rename_sample(self, selection):
        sample_id = selection.single
        if not sample_id or self.busy():
            return
        sample = self.project.get_sample(sample_id)
        name, accepted = QInputDialog.getText(self, "Rename isolate",
            "Display name. The sample identifier, its input file and every stored result stay "
            "exactly as they are:", text=sample["name"])
        if not accepted:
            return
        try:
            self.project.rename_sample(sample_id, name)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.notify("Renamed. No result, allele call or file was changed.")

    def context_rename_folder(self, selection):
        """Re-label every isolate in an organism folder, and move their copies with it."""
        if self.busy() or not selection.folder:
            return
        genus, species = selection.folder
        samples = self.context_records(selection)
        if not samples:
            self.notify("That folder holds no isolates to re-label.")
            return
        new_genus, accepted = QInputDialog.getText(self, "Rename this organism folder",
            f"Genus recorded for these {len(samples)} isolates. A folder name is a filing decision, "
            "not a laboratory identification:", text="" if genus == "Unknown" else genus)
        if not accepted:
            return
        new_species, accepted = QInputDialog.getText(self, "Rename this organism folder",
            "Species (leave empty if you only know the genus):",
            text="" if species in {"Unspecified", ""} else species)
        if not accepted:
            return
        self.apply_organism_assignments([self.label_assignment(sample, new_genus.strip(), new_species.strip())
                                         for sample in samples])

    def context_edit_annotations(self, selection):
        self.with_selection(selection, self.edit_metadata)

    def context_add_collection(self, selection):
        self.with_selection(selection, self.add_collection)

    def context_highlight(self, selection):
        self.with_selection(selection, self.highlight_selected)

    def with_selection(self, selection, action):
        """Run an existing selection-driven action on exactly what was right-clicked."""
        previous = self.selection_ids
        self.selection_ids = set(selection.sample_ids)
        try:
            action()
        finally:
            self.selection_ids = previous
        self.refresh_tables()

    def context_unhighlight(self, selection):
        from wmlstudio.sample_workflow import set_cluster
        identifiers = list(selection.sample_ids)
        if not identifiers:
            return
        for sample in self.context_records(selection):
            group = (sample.get("metadata") or {}).get("cluster", {})
            set_cluster(self.project, [sample["id"]], group.get("label") or "Highlighted",
                        group.get("color") or "#2F8A78", False)
        self.refresh()
        self.notify(f"Highlight removed from {len(identifiers)} isolates. The grouping you recorded "
                    "stays in each isolate's history.")

    def context_add_files(self, selection):
        self.browse_files()

    def context_add_folder(self, selection):
        self.browse_folder()

    def context_import_scheme(self, selection):
        self.import_scheme()

    def context_open_scheme_folder(self, selection):
        folder = next((Path(path) for path in selection.paths), None)
        if folder is None or not folder.is_dir():
            self.notify("That scheme folder is not available on this computer.")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    def context_remove_scheme(self, selection):
        import shutil

        from PySide6.QtWidgets import QMessageBox
        if self.busy():
            return
        root = Path(self.root).resolve()
        paths = [Path(path).resolve() for path in selection.paths]
        removable = [path for path in paths if root in path.parents and path.is_dir()]
        if not removable or len(removable) != len(paths):
            self.notify("Bundled reference snapshots are read-only; nothing was removed.")
            return
        names = ", ".join(path.name for path in removable)
        if QMessageBox.question(self, "Remove imported scheme?",
                f"Delete {names} from this computer? Results already produced with it keep their "
                "scheme fingerprint in the project, but you cannot rerun them until it is installed "
                "again.") != QMessageBox.StandardButton.Yes:
            return
        for path in removable:
            try:
                shutil.rmtree(path)
            except OSError as error:
                self.error(error)
                break
        self.populate_schemes()
        self.notify(f"Removed {len(removable)} imported scheme folders. Saved results are unchanged.")

    def context_archive(self, selection):
        from wmlstudio import archive
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        reason, accepted = QInputDialog.getText(self, f"Archive {len(identifiers)} isolates",
            "Archiving hides these isolates from the working views and destroys nothing: every "
            "result, allele call, file and history entry is kept, and you can restore them at any "
            "time.\n\nWhy are you archiving them? (optional)")
        if not accepted:
            return
        try:
            changed = archive.archive_samples(self.project, identifiers, reason.strip())
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        already = len(identifiers) - len(changed)
        self.notify(f"{len(changed)} isolates archived and hidden from the working views."
                    + (f" {already} were already archived." if already else "")
                    + " Nothing was deleted; find them under Archived in the library navigator.")

    def context_restore(self, selection):
        from wmlstudio import archive
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        try:
            changed = archive.restore_samples(self.project, identifiers)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.notify(f"{len(changed)} isolates restored to the working views with every record they "
                    "carried.")

    def context_remove(self, selection):
        from PySide6.QtWidgets import QCheckBox, QMessageBox

        from wmlstudio.storage import remove_samples
        identifiers = list(selection.sample_ids)
        if not identifiers or self.busy():
            return
        box = QMessageBox(self)
        box.setWindowTitle("Remove from this project?")
        box.setText(f"Remove {len(identifiers)} isolates and their saved analyses from this project?")
        box.setInformativeText(
            "Your original sequence files are never deleted. If you only want them out of the way, "
            "cancel and choose Archive instead — archiving keeps every result and can be undone.")
        box.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(QMessageBox.StandardButton.Cancel)
        check = QCheckBox("Also delete the managed copies WMLSTudio made in its own folders")
        box.setCheckBox(check)
        if box.exec() != QMessageBox.StandardButton.Yes:
            return
        try:
            report = remove_samples(self.project, identifiers,
                                    delete_managed_copy=check.isChecked())
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        message = (f"{len(report['removed'])} isolates removed. Their evidence is recorded in the "
                   "project history and can be restored from Samples ▸ Recently removed.")
        if report["deleted"]:
            message += f" {len(report['deleted'])} managed copies were deleted; your originals were not."
        if report["retained"]:
            message += f" {len(report['retained'])} copies were kept: {report['retained'][0][1]}"
        self.notify(message)

    def remove_sample(self):
        """Menu twin of the right-click removal, so both routes explain themselves."""
        from wmlstudio.context_menus import Selection
        identifiers = sorted(self.selection_ids)
        if not identifiers:
            sample = self.selected_sample()
            identifiers = [sample["id"]] if sample else []
        if not identifiers:
            self.notify("Select the isolates to remove first.")
            return
        self.context_remove(Selection("library", tuple(identifiers)))

    def fill_recently_removed(self, menu):
        """Offer the last removals back, saying plainly which ones cannot be restored."""
        menu.clear()
        entries = [entry for entry in self.project.history()
                   if entry["action"] == "sample_removed"][-20:]
        if not entries:
            menu.addAction("Nothing has been removed from this project").setEnabled(False)
            return
        for entry in reversed(entries):
            details = entry.get("details") or {}
            sample = details.get("sample") or {}
            action = menu.addAction(f"{sample.get('name') or 'Unnamed isolate'} · removed {entry['created_at']}")
            if details.get("format_version") != 2:
                action.setEnabled(False)
                action.setToolTip("This removal predates restorable records.")
                continue
            action.triggered.connect(lambda checked=False, i=entry["id"]: self.restore_removed(i))

    def restore_removed(self, history_id):
        try:
            sample_id = self.project.restore_removed_sample(history_id)
        except (KeyError, ValueError) as error:
            self.error(error)
            return
        self.refresh()
        self.selection_ids = {sample_id}
        self.refresh_tables()
        self.notify("Isolate restored with every saved analysis it had. Its original file was never "
                    "deleted.")

    def build_menus(self):
        self.command_actions = []

        def action(menu, title, callback, shortcut=None):
            item = QAction(title, self)
            if shortcut:
                item.setShortcut(QKeySequence(shortcut))
            item.triggered.connect(callback)
            menu.addAction(item)
            self.command_actions.append((title, item))
            return item

        bar = self.menuBar()
        file = bar.addMenu("&File")
        action(file, "New project…", self.new_project, "Ctrl+N")
        action(file, "Open project…", self.open_project_dialog)
        action(file, "Save project copy…", self.save_project_copy)
        file.addSeparator()
        imports = file.addMenu("Import")
        action(imports, "Sequence files…", self.browse_files)
        action(imports, "Sequence folder…", self.browse_folder)
        action(imports, "HYDRA report…", self.import_hydra)
        action(imports, "Allelic profile table…", self.import_profile_table)
        action(imports, "Sample bundle…", self.import_sample_bundle)
        exports = file.addMenu("Export")
        for title, fmt in [("Selected samples as CSV…", "csv"), ("Selected samples as TSV…", "tsv"), ("Selected samples as JSON…", "json"), ("HTML report…", "html"), ("PDF report…", "pdf")]:
            action(exports, title, lambda checked=False, f=fmt: self.export_project(f))
        action(exports, "Allelic profile table…", self.export_profile_table)
        action(exports, "Sample bundle…", self.export_sample_bundle)
        file.addSeparator()
        action(file, "Exit", self.close, "Alt+F4")
        samples = bar.addMenu("&Samples")
        action(samples, "Add samples…", self.browse_files)
        action(samples, "Add a folder of samples…", self.browse_folder)
        action(samples, "Identify waiting samples again", self.reidentify_waiting)
        for key in STEP_KEYS:
            action(samples, STEP_ACTIONS[key], lambda checked=False, k=key: self.run_step(k))
        samples.addSeparator()
        action(samples, "Assign organism / workflow…", self.assign_selected)
        action(samples, "Edit annotations…", self.edit_metadata)
        action(samples, "Add to collection…", self.add_collection)
        action(samples, "Highlight cluster…", self.highlight_selected)
        action(samples, "Select visible", self.select_visible_samples)
        action(samples, "Clear selection", self.clear_sample_selection)
        action(samples, "Clear this tab (selection, search and filters)", self.clear_step)
        action(samples, "Relink input (same bytes)…", self.relink_selected_input)
        action(samples, "Re-file managed copies now…", self.refile_selected)
        action(samples, "Attach original reads to assemblies…", self.attach_reads_selected)
        samples.addSeparator()
        action(samples, "Archive selected isolates…",
               lambda: self.context_archive(self.library_selection()))
        action(samples, "Restore selected isolates from the archive",
               lambda: self.context_restore(self.library_selection()))
        action(samples, "Remove sample…", self.remove_sample)
        recently_removed = samples.addMenu("Recently removed")
        recently_removed.aboutToShow.connect(lambda menu=recently_removed: self.fill_recently_removed(menu))
        analysis = bar.addMenu("&Analysis")
        action(analysis, "Choose isolates and analyse…", self.choose_and_analyse, "Ctrl+R")
        action(analysis, "Review pending isolates…", lambda: self.choose_and_analyse(pending_only=True))
        action(analysis, "Choose assemblies for HYDRA…", self.choose_hydra_cohort)
        action(analysis, "Characterize identity / virulence / accessory evidence…", self.run_characterization_selected)
        action(analysis, "Assemble and analyse read pairs…", lambda: self.choose_and_analyse(assemble=True))
        action(analysis, "Check missing loci with original reads…", self.open_read_support)
        action(analysis, "Cancel current task", self.cancel_analysis)
        analysis.addSeparator()
        action(analysis, "Build comparison from selected", self.compare_selected)
        action(analysis, "Compute selected comparison scheme…", self.type_comparison_scheme)
        action(analysis, "Create local scheme from selected assemblies…", self.create_adhoc_scheme)
        data = bar.addMenu("&Data")
        action(data, "Online scheme catalog…", self.open_reference_manager)
        action(data, "Research saved library…", self.open_library_research)
        action(data, "AMR databases / updates…", self.open_amr_databases)
        action(data, "Characterization references / updates…", self.install_characterization_references)
        action(data, "Install broader species panel…", self.install_species_panel)
        action(data, "Download practice data…", self.download_practice_cohort)
        action(data, "Import local scheme…", self.import_scheme)
        action(data, "Epidemiology grid / bulk import…", self.open_metadata_grid)
        action(data, "Published cluster threshold guidance…", self.open_threshold_guidance)
        action(data, "Open data folder", self.open_data_folder)
        action(data, "Open selected input folder", self.open_sample_folder)
        view = bar.addMenu("&View")
        for index, title in enumerate(self.nav_names):
            action(view, title, lambda checked=False, i=index: self.navigate(i), f"Alt+{index + 1}")
        action(view, "Refresh current workspace", self.refresh, "F5")
        action(view, "Fit comparison", lambda: self.tree.fit_tree())
        action(view, "Increase interface scale", lambda: self.set_ui_scale(self.ui_scale + 10), "Ctrl++")
        action(view, "Decrease interface scale", lambda: self.set_ui_scale(self.ui_scale - 10), "Ctrl+-")
        action(view, "Reset interface scale", lambda: self.set_ui_scale(100), "Ctrl+0")
        help_menu = bar.addMenu("&Help")
        action(help_menu, "Practice project", self.load_demo)
        action(help_menu, "Problem → solution guide", self.open_workflow_guide, "F1")
        action(help_menu, "Command search…", self.command_palette, "Ctrl+K")
        action(help_menu, "Workflow and limitations", self.open_workflow_guide)

    def command_palette(self):
        names = [title.replace("&", "") for title, _ in self.command_actions]
        title, accepted = QInputDialog.getItem(self, "Command search", "Choose a command (type to search):", names, 0, True)
        if accepted and title in names:
            self.command_actions[names.index(title)][1].trigger()

    def open_workflow_guide(self, topic=None):
        from wmlstudio.workflow_guide import WorkflowGuide
        try:
            if getattr(self, "workflow_guide", None) is None:
                self.workflow_guide = WorkflowGuide(self)
                self.workflow_guide.actionRequested.connect(self.journey_action)
            if isinstance(topic, str) and topic:
                self.workflow_guide.search.setText(topic)
            self.workflow_guide.show()
            self.workflow_guide.raise_()
        except (OSError, ValueError) as error:
            self.error(str(error))

    def attach_reads_selected(self):
        from wmlstudio.read_attachment_dialog import launch_read_attachment
        launch_read_attachment(self)

    def set_ui_scale(self, percent):
        from wmlstudio.interface_settings import interface_preferences, scaled_style
        percent = max(80, min(150, int(percent)))
        self.ui_scale = percent
        apply_application_style(scaled_style(percent))
        interface_preferences(self.root).setValue("scale", percent)
        self.updateGeometry()

    def set_graph_text_scale(self, percent):
        """Resize the lettering on both comparison trees; the trees themselves do not move."""
        from wmlstudio.widgets import set_graph_text_scale
        percent = set_graph_text_scale(percent)
        from wmlstudio.interface_settings import interface_preferences
        interface_preferences(self.root).setValue("graph_scale", percent)
        for name in ("tree", "baseline_tree"):
            view = getattr(self, name, None)
            if view is not None and hasattr(view, "_redraw"):
                view._redraw()

    def open_interface_settings(self):
        from wmlstudio.interface_settings import InterfaceSettingsDialog
        if getattr(self, "interface_dialog", None) is None:
            self.interface_dialog = InterfaceSettingsDialog(self)
            self.interface_dialog.setWindowModality(Qt.WindowModality.WindowModal)
        self.interface_dialog.show()
        self.interface_dialog.raise_()

    def show_sample_detail_for_id(self, sample_id):
        self.show_sample_detail(sample_id)

    def open_isolate_record(self, sample_id):
        from wmlstudio.isolate_dialog import IsolateRecordDialog
        if not isinstance(sample_id, str):
            return
        try:
            self.project.get_sample(sample_id)
        except KeyError:
            self.notify("That isolate is no longer in this project.")
            return
        if not hasattr(self, "isolate_dialogs"):
            self.isolate_dialogs = {}
        if sample_id not in self.isolate_dialogs:
            self.isolate_dialogs[sample_id] = IsolateRecordDialog(self, sample_id)
        dialog = self.isolate_dialogs[sample_id]
        dialog.refresh_record()
        dialog.show()
        dialog.raise_()

    def open_metadata_grid(self, checked=False):
        from wmlstudio.metadata_grid import launch_metadata_grid
        launch_metadata_grid(self)

    def open_read_support(self, sample_id=None):
        from wmlstudio.read_support_dialog import launch_read_support
        launch_read_support(self, sample_id=sample_id if isinstance(sample_id, str) else None)

    def open_threshold_guidance(self):
        from wmlstudio.threshold_dialog import ThresholdGuideDialog
        dialog = ThresholdGuideDialog(self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.apply_threshold_guidance(dialog.evidence)

    def choose_feature_cohort(self):
        from wmlstudio.cohort_picker import CohortPickerDialog
        dialog = CohortPickerDialog(self.project.samples(), self.project, "Which isolates should the evidence workspace show?", self.feature_ids, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.feature_ids = set(dialog.selected_ids)
        self.refresh_features()
        self.navigate(4)

    def library_selection(self):
        """The isolate-library selection as a right-click Selection, for menu twins."""
        from wmlstudio.context_menus import Selection
        return Selection("library", tuple(sorted(self.selection_ids)))

    def choose_and_analyse(self, checked=False, assemble=False, pending_only=False):
        if self.busy():
            return
        from wmlstudio.cohort_picker import CohortPickerDialog
        samples = self.project.samples()
        initial = {sample["id"] for sample in samples if sample["status"] in {"queued", "interrupted"}} if pending_only else None
        dialog = CohortPickerDialog(samples, self.project, "Which isolates should be analysed?", initial, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.start_analysis(confirm=True, all_samples=True, sample_ids=dialog.selected_ids, assemble=assemble)

    def choose_hydra_cohort(self):
        if self.busy():
            return
        from wmlstudio.cohort_picker import CohortPickerDialog
        samples = self.project.samples()
        dialog = CohortPickerDialog(samples, self.project, "Which assemblies should HYDRA analyse?", parent=self, include_reads=False)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        chosen = [sample for sample in samples if sample["id"] in dialog.selected_ids]
        self.start_reviewed_plan(chosen, self.review_run_plan(chosen, hydra=True))

    def build_schemes(self):
        super().build_schemes()
        layout = self.pages.widget(3).widget().layout()
        row = FlowLayout()
        row.addWidget(button("Browse online / install updates…", self.open_reference_manager, True))
        row.addWidget(button("Refresh installed schemes", self.populate_schemes))
        row.addWidget(button("Create local scheme…", self.create_adhoc_scheme))
        layout.insertLayout(2, row)
        layout.insertWidget(3, label("Updates create new snapshots; existing results retain their original fingerprint.", "small", True))

    def open_reference_manager(self):
        from wmlstudio.reference_dialog import ReferenceManagerDialog
        if self.reference_dialog is None:
            self.reference_dialog = ReferenceManagerDialog(self.root, self)
            self.reference_dialog.schemeInstalled.connect(self.reference_installed)
        self.reference_dialog.show()
        self.reference_dialog.raise_()
        self.reference_dialog.activateWindow()

    def reference_installed(self, path):
        self.populate_schemes()
        self.refresh_cohort_table()
        self.notify(f"Reference snapshot installed: {Path(path).name}. Existing results retain their original scheme fingerprint.")

    def create_adhoc_scheme(self):
        if self.busy():
            return
        samples = [sample for sample in self.selected_samples() if sample.get("input_path")]
        if not samples:
            self.notify("Select the assembly cohort in Samples first. Profile-only imports cannot define new sequence loci.")
            return
        from wmlstudio.adhoc import create_adhoc_scheme
        from wmlstudio.adhoc_dialog import AdhocSchemeDialog
        dialog = AdhocSchemeDialog(samples, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        plan = dict(dialog.plan)
        self.launch_task(lambda cancelled, progress: create_adhoc_scheme(library_root=self.root / "schemes",
                         cancelled=cancelled, progress=progress, **plan), "adhoc_scheme",
                         lambda result: self.reference_installed(result["path"]))

    def closeEvent(self, event):
        if not self.cancel_comparison_for_close():
            event.ignore()
            return
        if self.worker and self.worker.isRunning():
            return super().closeEvent(event)
        if self.amr_database_dialog is not None:
            database_worker = self.amr_database_dialog.worker
            if database_worker and database_worker.isRunning():
                database_worker.cancel()
                event.ignore()
                QTimer.singleShot(100, self.close)
                return
            self.amr_database_dialog.close()
        if self.reference_dialog is not None:
            reference_worker = self.reference_dialog.worker
            if reference_worker and reference_worker.isRunning():
                reference_worker.cancel()
                event.ignore()
                QTimer.singleShot(100, self.close)
                return
            self.reference_dialog.close()
        if self.library is not None:
            self.library.index_project(self.project)
            self.library.close()
            self.library = None
        super().closeEvent(event)
