"""Two scheme libraries — classical MLST and cgMLST — and the public catalogs behind them.

One tab per library, because a seven-locus MLST distance and a two-thousand-target
cgMLST distance are different quantities. Each tab lists what is installed, what is
catalogued and can be downloaded, and what a provider only permits under its own
terms; every row names the organism, the scheme, the size in that library's own
unit, the provider and the version. Nothing is downloaded until a row is selected
and its terms are shown.
"""

import threading
import time
from pathlib import Path

from PySide6.QtCore import Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from . import cgmlst_schemes
from .paths import cgmlst_library_roots, resource_root, scheme_locations
from .progress_words import describe_progress  # noqa: F401  (this dialog's own wording)
from .reference_catalog import CGMLSTOrgCatalog, PasteurCatalog, PubMLSTCatalog
from .reference_index import scheme_library
from .sequence import AnalysisCancelled

KIND_TITLES = {"cgmlst": "cgMLST schemes", "mlst": "Classical MLST schemes"}
LIBRARY_NAMES = {"cgmlst": "cgMLST", "mlst": "classical MLST"}
UNIT = {"cgmlst": "targets", "mlst": "loci"}
SEPARATION = ("A classical seven-locus MLST distance and a 2,000-target cgMLST distance are different "
              "quantities: they never share a scale, an axis, a column or a threshold. That is why these "
              "are two libraries and two tabs, not one list.")
# The three target sets a provider may publish for one organism, named where the
# choice between them is made. The keys are cgmlst_schemes' own target_set_detail
# values, so a set this dialog has no words for is labelled as unrecorded rather
# than quietly described as a core genome.
TARGET_SET_PREFIXES = {"core": "Core target set",
                       "accessory": "Accessory target set",
                       "whole_genome": "Whole-genome set (core + accessory)"}
TARGET_SET_PURPOSE = {
    "core": ("Core: the targets expected in every isolate of this organism. This is what cgMLST "
             "means, and a published outbreak cutoff is derived on a set of this kind."),
    "accessory": ("Accessory: targets that some isolates of this organism carry and others do not, "
                  "by design. A target not called here is biology and not a failed call, so it is "
                  "reported as not assayed and never counted as a difference."),
    "whole_genome": ("Whole genome: core and accessory targets in ONE set, which is sometimes "
                     "called cgMLST plus accessory genes. Run it instead of a core scheme, never "
                     "beside one as something added to a core distance."),
}
# The point of offering the choice at all, and what the catalogue's own sentence
# about scales and cutoffs does not say: the consequence for the isolates. Shown
# wherever more than one target set is on offer, because choosing between them is
# exactly the moment at which two different quantities can be taken for two
# versions of one.
TARGET_SET_SEPARATION = (
    "These are DIFFERENT QUANTITIES, not versions of one measurement: two isolates typed against "
    "different target sets are not comparable at all. Type every isolate of one comparison against "
    "the set you choose here.")
TARGET_SET_NOT_RECORDED = (
    "Target set: not recorded. This folder does not say whether it holds a core target set, an "
    "accessory one or both, so a distance measured on it cannot be read as a cgMLST core distance "
    "and no published cutoff is offered for it.")
# Core first, then accessory, then the pan-genome set: the order in which a person
# decides, from the set a cgMLST scheme IS outwards.
_TARGET_SET_ORDER = {"core": 0, "accessory": 1, "whole_genome": 2}
# cgMLST.org publishes allele nomenclature but no central profile table, so a
# snapshot taken from it assigns no cgST. The same sentence the online catalog
# records on a searched entry is recorded on a catalogued one.
CGMLST_ORG_ACCESS = ("Allele nomenclature is public; no central ST/profile table is supplied by this "
                     "download.")
PROVIDER_CLIENTS = {"pubmlst": "PubMLST", "pasteur": "BIGSdb-Pasteur", "cgmlst_org": "cgMLST.org"}
# Installed first: after a download the thing a person is looking for is what they
# now have, not what they could get next.
STATE_ORDER = {"installed": 0, "catalogued": 1, "online": 2}


class _ReferenceWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int, str)

    def __init__(self, operation, parent=None):
        super().__init__(parent)
        self.operation = operation
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        try:
            result = self.operation(self.cancel_event.is_set, self.progress.emit)
            self.succeeded.emit(result)
        except AnalysisCancelled:
            self.progress.emit(0, 1, "Reference operation cancelled")
        except Exception as error:
            self.failed.emit(str(error))


def _count_words(count, kind):
    """'2358 targets' or '7 loci' — the unit is the point and is never dropped."""
    return f"{count} {UNIT[kind]}" if count else f"no {UNIT[kind]} installed"


def _row_target_set(row) -> str:
    """Which target set a row records, preferring the finer of the two fields.

    ``target_set`` says core or not-core, which is what most of the application
    stores; ``target_set_detail`` separates an accessory set from a pan-genome one,
    because 251 accessory targets and 1,907 pan-genome targets are themselves two
    different quantities. A row that records neither returns '' and is labelled as
    unrecorded — never as core, which is what an empty value would otherwise be
    read as.
    """
    return str(row.get("target_set_detail") or row.get("target_set") or "")


def _target_set_words(row) -> str:
    """'Core', 'Accessory', 'Whole genome' or 'Not recorded' for one row."""
    detail = _row_target_set(row)
    return cgmlst_schemes.target_set_label(detail) if detail else "Not recorded"


def download_entry(entry) -> dict:
    """The catalog-client entry that installs one pinned cgMLST scheme.

    Built from the pin, so a catalogued scheme can be installed straight from the
    library list without searching online first. The client still re-reads the
    provider's own definition before writing anything, so a scheme that has been
    redefined upstream is reported rather than silently accepted.
    """
    token = cgmlst_schemes.provider_token(entry["provider"])
    common = {"id": entry["key"], "organism": entry["organism"], "name": entry["scheme_name"],
              "locus_count": entry["locus_count"], "type": "cgMLST", "url": entry["source_url"],
              "provider": entry["provider_name"], "access_notice": "",
              # Carried on the entry so the confirmation a person is shown is this
              # provider's own restriction, never a nearby provider's.
              "terms_url": entry["terms_url"], "terms_notice": entry["licence_restriction"]}
    if token == "cgmlst_org":
        return {**common, "slug": entry["scheme_id"], "access_notice": CGMLST_ORG_ACCESS}
    return {**common, "database": entry["database"], "scheme_id": entry["scheme_id"]}


def _scheme_group(key) -> str:
    """The provider-and-organism family of a catalogue key, or '' when it has none.

    A key this catalogue does not pin resolves to no family rather than to a
    guessed one, because the family is what decides which target sets are offered
    beside each other as a choice.
    """
    if not key:
        return ""
    try:
        return cgmlst_schemes.entry_for(key)["scheme_group"]
    except cgmlst_schemes.SchemeCatalogError:
        return ""


def multi_set_organisms() -> list[str]:
    """Organisms whose catalogued schemes are not all the same kind of target set.

    Named in the library summary because the choice between running cgMLST and
    running cgMLST plus accessory genes is otherwise invisible until a row that
    has one happens to be selected, and a choice nobody can find is a choice
    nobody makes. Two core schemes from one provider are not such a pair: they are
    two schemes, and choosing between them is a different question.
    """
    families: dict[str, list[dict]] = {}
    try:
        for entry in cgmlst_schemes.catalog_entries():
            families.setdefault(entry["scheme_group"], []).append(entry)
    except (OSError, ValueError):
        return []
    return sorted({rows[0]["organism"] for rows in families.values()
                   if len({cgmlst_schemes.target_set_detail(row) for row in rows}) > 1})


def _installed_detail(row, source) -> str:
    """Why this folder is in this library, where it is, and what is bound to it."""
    lines = [source["title"], f"Installed at {source['path']}",
             f"Classified as {KIND_TITLES[source['kind']]}: {source['kind_basis']}."]
    if source.get("kind_conflict"):
        lines.append("Conflict: " + source["kind_conflict"])
    guidance = source.get("threshold") or {}
    if guidance.get("reason"):
        lines.append(guidance["reason"])
    if guidance.get("status") == "threshold_offerable" and guidance.get("threshold") is not None:
        lines.append(f"A reviewed cutoff of {guidance['threshold']} is bound to this scheme. It is "
                     "evidence to review, never a setting this application applies for you.")
    if source.get("catalogued") is False:
        lines.append("This scheme is not in the pinned catalogue, so no published cutoff is bound to "
                     "it. Distances within it are still perfectly usable.")
    if source.get("expected_locus_count") and not source.get("count_matched", True):
        lines.append(f"The catalogue pins {source['expected_locus_count']} {UNIT[source['kind']]} and "
                     f"this folder holds {source['locus_count']}. A cutoff published for the complete "
                     "set is not offered for a partial one.")
    if row["kind"] == "cgmlst":
        # Which target set is in this folder decides what a distance from it means,
        # so it is stated on the folder rather than left to the catalogue row a
        # person may never open.
        detail = _row_target_set(source)
        if detail:
            lines.append(f"Target set: {cgmlst_schemes.target_set_label(detail)}. "
                         + TARGET_SET_PURPOSE.get(detail, ""))
        else:
            lines.append(TARGET_SET_NOT_RECORDED)
        lines.append(cgmlst_schemes.INTERPRETATION)
    return "\n\n".join(lines)


def _catalogued_detail(entry) -> str:
    """Everything the download button must show before a byte is transferred."""
    plan = cgmlst_schemes.download_plan(entry["key"])
    # What this scheme measures comes before how it is fetched: a person about to
    # spend an hour downloading a 1,907-target pan-genome set must read that it is
    # not a core cgMLST scheme before the first byte moves, not afterwards.
    lines = [plan["title"],
             f"Target set: {plan['target_set_label']} · {plan['locus_count']} targets. "
             + " ".join(plan["target_set_notice"].split()),
             f"Not installed. It will install into your cgMLST library, in the folder "
             f"{entry['folder']}.",
             "How it downloads: " + plan["method"],
             f"Provider: {plan['provider_name']} · {plan['source_url']}",
             "Terms: " + plan["terms_notice"] + " " + plan["terms_url"]]
    if plan["licence_quote"]:
        lines.append("The provider's own words: “" + plan["licence_quote"] + "”")
    if plan["requires_acknowledgement"]:
        lines.append("This provider does not grant redistribution, so this scheme is never packed into "
                     "WMLSTudio. You will be asked to confirm that your intended use is permitted "
                     "before anything is downloaded.")
    guidance = entry.get("threshold") or {}
    if guidance.get("reason"):
        lines.append(guidance["reason"])
    lines.extend(plan["notes"])
    lines.append(cgmlst_schemes.INTERPRETATION)
    return "\n\n".join(lines)


def _online_detail(entry, kind) -> str:
    lines = [f"{entry['organism']} · {entry['name']} · "
             f"{_count_words(entry.get('locus_count') or 0, kind)}",
             f"Source: {entry['url']}",
             entry.get("access_notice")
             or "No anonymous-access restriction was reported in the scheme metadata."]
    if entry.get("terms_notice"):
        lines.append(entry["terms_notice"] + "\n" + str(entry.get("terms_url") or ""))
    lines.append("Found in the online catalog. Nothing is downloaded until you choose it.")
    return "\n\n".join(lines)


def installed_rows(root, *, cancelled=None) -> dict:
    """Every installed scheme folder, split into the two libraries plus the unreadable.

    The split is reference_index's classification, read from what each folder
    records and holds — never from which library folder it happens to sit in, so a
    cgMLST scheme a user dropped into the classical folder by hand is still listed
    as cgMLST here.
    """
    root = Path(root)
    grouped = scheme_library(scheme_locations(root), cancelled=cancelled)
    bundled = resource_root().resolve()
    catalogued = {str(Path(item["path"]).resolve()): item
                  for item in cgmlst_schemes.installed_entries(root, cancelled=cancelled)}
    rows = {"mlst": [], "cgmlst": [], "unknown": []}
    for kind, entries in grouped.items():
        for entry in entries:
            resolved = Path(entry["path"]).resolve()
            source = {**entry, **catalogued.get(str(resolved), {})}
            packed = bundled in resolved.parents
            key = source.get("catalog_key") or ""
            row = {"row_id": "path:" + str(resolved), "kind": kind, "state": "installed",
                   "status": "Installed · bundled" if packed else "Installed · your library",
                   "organism": entry["organism_label"], "scheme": entry["scheme_label"],
                   "count": entry["locus_count"], "target_set": source.get("target_set") or "",
                   "target_set_detail": source.get("target_set_detail") or "",
                   "provider": entry["provider"], "version": entry["version"],
                   "title": entry["title"], "location": entry["path"], "path": entry["path"],
                   "key": key, "entry": None,
                   # An installed scheme this catalogue recognises carries its
                   # provider-and-organism family, so the target-set menu offers the
                   # core and accessory sets of THIS scheme rather than a different
                   # provider's scheme for the same organism, which is not the
                   # accessory half of anything.
                   "scheme_group": _scheme_group(key),
                   "requires_terms": False, "downloadable": False,
                   "provider_key": ""}
            row["detail"] = _installed_detail(row, source)
            rows[kind].append(row)
    return rows


def catalogued_rows(root, *, installed_paths=()) -> list[dict]:
    """One row per pinned cgMLST scheme that is not already installed.

    The pinned catalogue is what makes a scheme findable before it exists on disk:
    the row names the organism, the target count, the provider and the folder the
    download will fill, so the place a person browses and the place a scheme lands
    are the same place.
    """
    known = {str(Path(path).resolve()) for path in installed_paths}
    rows = []
    for entry in cgmlst_schemes.library_status(root):
        installed = entry.get("installed") or {}
        if installed and str(Path(installed["path"]).resolve()) in known:
            continue  # Already listed from disk, with its real installed target count.
        rows.append({
            "row_id": "key:" + entry["key"], "kind": "cgmlst", "state": "catalogued",
            "status": ("Not installed · download under the provider's terms"
                       if entry["requires_terms_acknowledgement"] else "Not installed · download"),
            "organism": entry["organism"], "scheme": entry["scheme_name"],
            "count": entry["locus_count"], "target_set": entry["target_set"],
            "target_set_detail": entry["target_set_detail"],
            "provider": entry["provider_name"], "version": entry["version"],
            "title": entry["title"], "location": entry["folder"], "path": "",
            "key": entry["key"], "entry": download_entry(entry),
            "scheme_group": entry["scheme_group"],
            "requires_terms": bool(entry["requires_terms_acknowledgement"]),
            "downloadable": True,
            "provider_key": PROVIDER_CLIENTS.get(cgmlst_schemes.provider_token(entry["provider"]), ""),
            "detail": _catalogued_detail(entry)})
    return rows


def online_row(entry, kind, provider_key) -> dict:
    """A scheme found in a public catalog, in the same row shape as an installed one."""
    count = int(entry.get("locus_count") or 0)
    organism = str(entry.get("organism") or "Organism not recorded")
    scheme = str(entry.get("name") or "Unnamed scheme")
    return {"row_id": "online:" + str(entry.get("id") or entry.get("url") or scheme),
            "kind": kind, "state": "online", "status": "Found online · not installed",
            "organism": str(entry.get("organism") or ""), "scheme": scheme,
            # An online search result says nothing about which target set it is, and
            # an unlabelled set is an unknown quantity rather than a core genome.
            "count": count, "target_set": "", "target_set_detail": "",
            "provider": str(entry.get("provider") or provider_key),
            "version": str(entry.get("last_updated") or ""),
            "title": f"{organism} · {scheme} · {_count_words(count, kind)}",
            "location": "Not installed", "path": "", "key": "", "entry": dict(entry),
            "scheme_group": "", "requires_terms": bool(entry.get("terms_notice")),
            "downloadable": True, "provider_key": provider_key,
            "detail": _online_detail(entry, kind)}


class SchemeLibraryPage(QWidget):
    """One library surface: searchable rows, one detail pane, one unit of measure."""

    selectionChanged = Signal()

    COLUMNS = {"cgmlst": ["Organism", "Scheme", "Targets", "Target set", "Provider", "Version",
                          "Status", "Where it is"],
               "mlst": ["Organism", "Scheme", "Loci", "Provider", "Version", "Status",
                        "Where it is"]}

    def __init__(self, root, kind, parent=None):
        super().__init__(parent)
        self.root = Path(root)
        self.kind = kind
        self.rows = []
        self.visible = []
        self.online = []
        self.unreadable = 0
        self.entries = []  # the download descriptors of the rows currently shown
        # Read once: the pinned catalogue does not change while the dialog is open.
        self.multi_set = multi_set_organisms() if kind == "cgmlst" else []
        layout = QVBoxLayout(self)
        layout.addWidget(self._wrapped(
            f"Every {LIBRARY_NAMES[kind]} scheme this computer has, plus the ones WMLSTudio knows how "
            "to fetch. " + SEPARATION))
        row = QHBoxLayout()
        self.query = QLineEdit()
        self.query.setPlaceholderText("Filter by organism, scheme or provider — for example Klebsiella")
        self.query.textChanged.connect(self.apply_filter)
        self.provider = QComboBox()
        self.provider.addItems(["cgMLST.org", "PubMLST", "BIGSdb-Pasteur"] if kind == "cgmlst"
                               else ["PubMLST", "BIGSdb-Pasteur"])
        self.search_button = QPushButton("Search this provider online")
        row.addWidget(self.query, 1)
        row.addWidget(self.provider)
        row.addWidget(self.search_button)
        layout.addLayout(row)
        filters = QHBoxLayout()
        self.min_loci, self.max_loci = QSpinBox(), QSpinBox()
        for spin in (self.min_loci, self.max_loci):
            spin.setRange(0, 100000)
        self.max_loci.setValue(100000)
        filters.addWidget(QLabel(f"Online search: {UNIT[kind]} from"))
        filters.addWidget(self.min_loci)
        filters.addWidget(QLabel("to"))
        filters.addWidget(self.max_loci)
        filters.addStretch()
        self.folder_button = QPushButton("Open this library folder")
        self.folder_button.clicked.connect(self.open_library_folder)
        filters.addWidget(self.folder_button)
        layout.addLayout(filters)
        self.table = QTableWidget(0, len(self.COLUMNS[kind]))
        self.table.setHorizontalHeaderLabels(self.COLUMNS[kind])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self.selectionChanged.emit)
        layout.addWidget(self.table, 1)
        if kind == "cgmlst":
            # Requirement of the science, not of the layout: a core target set, an
            # accessory set and a whole-genome set are different quantities, so where
            # a provider publishes more than one the choice is explicit, and where it
            # publishes one that is said in words instead of shown as an empty menu.
            # It is a titled block of its own because choosing the target set is the
            # cgMLST decision, and it was reported as impossible to find.
            heading = self._wrapped("Target set to run — core, accessory, or core plus accessory")
            heading.setStyleSheet("font-weight: 600;")
            layout.addWidget(heading)
            variants = QHBoxLayout()
            variants.addWidget(QLabel("Target set:"))
            self.variant = QComboBox()
            self.variant.setMinimumWidth(460)
            self.variant.activated.connect(self.choose_variant)
            variants.addWidget(self.variant, 1)
            variants.addStretch()
            layout.addLayout(variants)
            self.variant_note = self._wrapped("")
            layout.addWidget(self.variant_note)
        else:
            self.variant, self.variant_note = None, None
        self.notice = QTextBrowser()
        self.notice.setMaximumHeight(160)
        self.notice.setPlainText(
            "Your sample sequences are never sent to a reference service. Public access may exclude "
            "newer alleles; any access restriction appears here and in the downloaded snapshot.")
        layout.addWidget(self.notice)
        self.summary = self._wrapped("")
        layout.addWidget(self.summary)

    @staticmethod
    def _wrapped(text):
        widget = QLabel(text)
        widget.setWordWrap(True)
        return widget

    # -- the library -----------------------------------------------------------

    def library_folder(self) -> Path:
        return (cgmlst_library_roots(self.root)[-1] if self.kind == "cgmlst"
                else self.root / "schemes")

    def open_library_folder(self):
        folder = self.library_folder()
        folder.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    def set_library(self, installed, catalogued=(), unreadable=0):
        """Replace the local rows, keeping whatever was found online this session."""
        self.rows = list(installed) + list(catalogued) + list(self.online)
        self.rows.sort(key=lambda row: (STATE_ORDER[row["state"]], row["organism"].casefold(),
                                        row["scheme"].casefold(), row["count"]))
        self.unreadable = unreadable
        self.apply_filter()

    def add_online(self, rows):
        """Add online results, dropping any that duplicate a row already listed."""
        known = {row["row_id"] for row in self.rows}
        added = [row for row in rows if row["row_id"] not in known]
        self.online.extend(added)
        self.rows.extend(added)
        self.rows.sort(key=lambda row: (STATE_ORDER[row["state"]], row["organism"].casefold(),
                                        row["scheme"].casefold(), row["count"]))
        self.apply_filter()
        if added:  # A search that found something shows what it found.
            self.select_row_id(added[0]["row_id"])
        return added

    def apply_filter(self):
        query = self.query.text().strip().casefold()
        selected = self.selected_row()
        self.visible = [row for row in self.rows if not query or query in " ".join(
            [row["organism"], row["scheme"], row["provider"], row["title"], row["status"]]).casefold()]
        self.entries = [row["entry"] for row in self.visible if row["entry"] is not None]
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.visible))
        for index, row in enumerate(self.visible):
            values = [row["organism"] or "Organism not recorded", row["scheme"], row["count"]]
            if self.kind == "cgmlst":
                values.append(_target_set_words(row))
            values += [row["provider"] or "Not recorded", row["version"] or "Not recorded",
                       row["status"], row["location"]]
            for column, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                item.setToolTip(row["title"])
                if column == 2:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(index, column, item)
        self.table.resizeColumnsToContents()
        self.table.blockSignals(False)
        if selected is not None:
            self.select_row_id(selected["row_id"])  # Keep the reader where they were.
        elif self.visible:
            self.table.selectRow(0)
        self.show_entry()
        self.refresh_summary()

    def refresh_summary(self):
        counts = {state: sum(row["state"] == state for row in self.rows)
                  for state in ("installed", "catalogued", "online")}
        text = (f"{counts['installed']} installed · {counts['catalogued']} catalogued and ready to "
                f"download · {counts['online']} found online · {len(self.visible)} shown. "
                f"Your {LIBRARY_NAMES[self.kind]} library: {self.library_folder()}")
        if self.unreadable:
            text += (f" · {self.unreadable} installed folder(s) record no readable typing kind and are "
                     "offered by neither library.")
        if self.multi_set:
            text += (" · More than one target set is catalogued for "
                     + ", ".join(self.multi_set)
                     + ": select one of their rows to choose between the core set and the "
                       "accessory or whole-genome set beside it.")
        self.summary.setText(text)

    # -- selection -------------------------------------------------------------

    def selected_row(self):
        index = self.table.currentRow()
        return self.visible[index] if 0 <= index < len(self.visible) else None

    def select_row_id(self, row_id) -> bool:
        for index, row in enumerate(self.visible):
            if row["row_id"] == row_id:
                self.table.selectRow(index)
                return True
        return False

    def select_path(self, path) -> bool:
        return self.select_row_id("path:" + str(Path(path).resolve()))

    def show_entry(self):
        row = self.selected_row()
        self.notice.setPlainText(row["detail"] if row else
                                 "Select a scheme to see its provider, its terms and where it "
                                 "installs. Nothing is downloaded until you choose it.")
        self.refresh_variants(row)

    # -- which target set to run -----------------------------------------------

    def refresh_variants(self, row):
        """Offer every target set catalogued for this scheme; say plainly when one is.

        Core, accessory and whole-genome sets are offered together, each with the
        number of targets it holds and what it is for, because "run cgMLST" and
        "run cgMLST plus accessory genes" are two different runs producing two
        quantities that are never compared.

        The choice is scoped to ONE provider's scheme family, because those sets
        are what one provider publishes for one organism. Another provider's
        scheme for the same organism is a different scheme altogether, not the
        accessory half of this one, so it is named in words rather than dropped
        into the same menu.
        """
        if self.variant is None:
            return
        self.variant.blockSignals(True)
        self.variant.clear()
        if row is None or not (row["scheme_group"] or row["organism"]):
            self.variant.setEnabled(False)
            self.variant_note.setText(
                "Select a scheme to see which target sets are catalogued for it: a core set, an "
                "accessory set, a whole-genome set, or only one of them." if row is None else
                "This scheme records no organism, so no core or accessory pairing can be looked up "
                "for it. That is a gap in what the scheme records, not proof that none exists.")
            self.variant.blockSignals(False)
            return
        group = row["scheme_group"] or None
        variants = cgmlst_schemes.scheme_variants(row["organism"] or None, group=group)
        # Built from the union rather than from the three split lists, so a target
        # set this dialog has no name for is still offered — labelled as unrecorded
        # — instead of disappearing from a menu that claims to hold every set.
        choices = sorted(((cgmlst_schemes.target_set_detail(entry), entry)
                          for entry in variants["core"] + variants["accessory"]),
                         key=lambda item: (_TARGET_SET_ORDER.get(item[0], 3),
                                           item[1]["locus_count"]))
        for position, (detail, entry) in enumerate(choices):
            self.variant.addItem(self._variant_text(detail, entry), entry["key"])
            self.variant.setItemData(position, f"{entry['title']}\n\n"
                                     + TARGET_SET_PURPOSE.get(detail, ""),
                                     Qt.ItemDataRole.ToolTipRole)
        # Enabled the moment there is something to choose between. Several organisms
        # now have two or three catalogued sets, and a menu that stays greyed out on
        # all of them is a choice a person cannot reach.
        self.variant.setEnabled(len(choices) > 1)
        index = self.variant.findData(row["key"])
        self.variant.setCurrentIndex(index if index >= 0 else 0)
        self.variant.blockSignals(False)
        note = variants["message"]
        if not choices:
            note = ("This scheme is not in the pinned catalogue, so no core or accessory pairing is "
                    "known for it. " + note)
        elif len(choices) == 1:
            note = f"{self._set_prefix(choices[0][0])} only. {note}"
        # What each set on offer is for, once each, in the order they are listed;
        # and, where there is a choice at all, why the numbers they produce cannot
        # be read against one another.
        said = [TARGET_SET_PURPOSE[detail] for detail in dict.fromkeys(item[0] for item in choices)
                if detail in TARGET_SET_PURPOSE]
        if len(choices) > 1:
            said.append(TARGET_SET_SEPARATION)
        self.variant_note.setText(" ".join([note, *said]) + self._other_schemes(row, choices))

    @staticmethod
    def _set_prefix(detail) -> str:
        return TARGET_SET_PREFIXES.get(detail, "Target set not recorded")

    def _variant_text(self, detail, entry) -> str:
        """One menu line: which set it is, how many targets it holds, whose it is.

        The target count is on the line because it is the difference between the
        sets as well as the reason their distances are not comparable: 1,649 core
        targets and 1,907 pan-genome targets are two counts of two different things.
        """
        return (f"{self._set_prefix(detail)} — {entry['locus_count']} targets — "
                f"{entry['scheme_name']} · {entry['provider_name']}")

    def _other_schemes(self, row, choices) -> str:
        """Name the other catalogued schemes for this organism, without offering them."""
        if not row["organism"]:
            return ""
        chosen = {entry["key"] for _, entry in choices}
        wide = cgmlst_schemes.scheme_variants(row["organism"])
        others = [entry for entry in wide["core"] + wide["accessory"] if entry["key"] not in chosen]
        if not others:
            return ""
        return (f" {len(others)} other cgMLST scheme(s) are catalogued for this organism: "
                + "; ".join(f"{entry['scheme_name']} · {entry['locus_count']} targets · "
                            f"{cgmlst_schemes.target_set_label(entry)} · {entry['provider_name']}"
                            for entry in others)
                + ". They are different target sets, not the accessory half of this one, and they "
                  "never share a distance, an axis or a cutoff with it.")

    def choose_variant(self, index):
        key = self.variant.itemData(index)
        if key and not self.select_row_id("key:" + str(key)):
            # The chosen variant is installed, so it is listed by its folder rather
            # than by its catalogue key.
            for row in self.visible:
                if row["key"] == key:
                    self.select_row_id(row["row_id"])
                    return


class ReferenceManagerDialog(QDialog):
    """The scheme manager: two libraries, one background worker, one download button."""

    schemeInstalled = Signal(str)
    installed = Signal(str)

    def __init__(self, root, parent=None, *, catalog=None):
        super().__init__(parent)
        self.root = Path(root)
        self.catalog = catalog or PubMLSTCatalog(self.root / "reference_cache")
        # An injected catalog stands in for every provider, so a test drives the
        # whole surface — including a cgMLST download — through one fake client.
        self.catalogs = ({name: catalog for name in PROVIDER_CLIENTS.values()} if catalog else
                         {"PubMLST": self.catalog, "BIGSdb-Pasteur": PasteurCatalog(),
                          "cgMLST.org": CGMLSTOrgCatalog()})
        self.worker = None
        self.catalog_errors = []
        self._close_pending = False
        self._started = time.monotonic()
        # Where the stage now counting has got to, so a stage that restarts its
        # count is recognised as a new one and timed from its own start.
        self._progress_at = 0
        self._progress_total = 0
        self.setWindowTitle("Scheme libraries — classical MLST and cgMLST")
        self.setSizeGripEnabled(True)
        # Everything above the buttons scrolls. 1080x760 was asked for
        # unconditionally, so on a laptop -- or at 125% Windows scaling, which
        # makes every widget bigger without making the screen bigger -- the
        # dialog was taller than the work area and the button row sat below the
        # bottom edge, unreachable. The buttons now live outside the scroll area
        # and are always on screen; the size is clamped to the screen that will
        # actually show it.
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        body = QWidget()
        self._scroll.setWidget(body)
        outer.addWidget(self._scroll, 1)
        layout = QVBoxLayout(body)
        intro = QLabel("Find, download and check typing schemes. cgMLST schemes install into their own "
                       "library under a name that says the organism, the provider and the target count.")
        intro.setWordWrap(True)
        layout.addWidget(intro)
        self.tabs = QTabWidget()
        self.pages = {}
        for kind in ("cgmlst", "mlst"):  # cgMLST first: it is the one that matters most here.
            page = SchemeLibraryPage(self.root, kind, self)
            page.search_button.clicked.connect(self.search)
            page.query.returnPressed.connect(self.search)
            page.selectionChanged.connect(self.show_entry)
            self.pages[kind] = page
            self.tabs.addTab(page, KIND_TITLES[kind])
        self.tabs.currentChanged.connect(self.show_entry)
        layout.addWidget(self.tabs, 1)
        self.status = QLabel("Ready")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        # The application bar is a 7px sliver with transparent text, which is right
        # for a background task and wrong for a download that runs for many
        # minutes: there was nothing on screen that visibly moved, so a working
        # download read as a frozen one. Here it is tall enough to show its own
        # percentage.
        self.progress.setTextVisible(True)
        self.progress.setStyleSheet(
            "QProgressBar { max-height: 18px; min-height: 18px; color: #E7EFF9; }")
        layout.addWidget(self.progress)
        buttons = QHBoxLayout()
        self.download_button = QPushButton("Download selected scheme")
        self.download_button.setEnabled(False)
        self.download_button.clicked.connect(self.download)
        self.refresh_button = QPushButton("Rescan installed schemes")
        self.refresh_button.clicked.connect(self.refresh_library)
        self.cancel_button = QPushButton("Cancel operation")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        buttons.addWidget(self.download_button)
        buttons.addWidget(self.refresh_button)
        buttons.addStretch()
        buttons.addWidget(self.cancel_button)
        buttons.addWidget(close)
        buttons.setContentsMargins(9, 6, 9, 9)
        outer.addLayout(buttons)
        self._fit_to_screen()
        self.refresh_library()

    def _fit_to_screen(self):
        """Open at the intended size, or the largest that fits, whichever is smaller."""
        preferred = (1080, 760)
        parent = self.parentWidget()
        window = parent.window() if parent is not None else None
        screen = ((window.screen() if window is not None else None)
                  or self.screen() or QGuiApplication.primaryScreen())
        if screen is None:
            self.resize(*preferred)
            return
        room = screen.availableGeometry()
        width = min(preferred[0], int(room.width() * .94))
        height = min(preferred[1], int(room.height() * .92))
        # A minimum larger than the screen would defeat the scroll area, so the
        # dialog is allowed to be small and scroll rather than refuse to shrink.
        self.setMinimumSize(min(560, width), min(420, height))
        self.resize(min(self.width() or width, width), min(self.height() or height, height))
        # A dialog that opened before the screen was known can already be off the
        # edge, so it is brought back rather than merely resized.
        centre = room.center()
        self.move(max(room.left(), centre.x() - self.width() // 2),
                  max(room.top(), centre.y() - self.height() // 2))

    # -- the current page, so one set of controls serves both libraries ---------

    @property
    def page(self):
        return self.tabs.currentWidget()

    @property
    def query(self):
        return self.page.query

    @property
    def table(self):
        return self.page.table

    @property
    def notice(self):
        return self.page.notice

    @property
    def provider(self):
        return self.page.provider

    @property
    def min_loci(self):
        return self.page.min_loci

    @property
    def max_loci(self):
        return self.page.max_loci

    @property
    def entries(self):
        return self.page.entries

    def show_page(self, kind):
        self.tabs.setCurrentWidget(self.pages[kind])

    # -- the library -----------------------------------------------------------

    def refresh_library(self):
        """Re-read both libraries from disk, and put a misfiled cgMLST scheme right.

        An earlier release installed every download into the classical schemes
        folder, where a 2,358-target scheme sat between two seven-locus ones under a
        digest-shaped name. Opening this dialog moves such a scheme into the cgMLST
        library and says so, so the scheme a person downloaded yesterday is where
        they are now told to look.
        """
        migrated = []
        try:
            report = cgmlst_schemes.prepare_library(self.root)
            migrated = report["migrated"]
        except (OSError, ValueError) as error:
            self.status.setText(f"The cgMLST library could not be prepared: {error}")
        try:
            grouped = installed_rows(self.root)
        except OSError as error:
            self.status.setText(f"Installed schemes could not be read: {error}")
            return
        unreadable = len(grouped["unknown"])
        paths = [row["path"] for row in grouped["cgmlst"]]
        try:
            catalogued = catalogued_rows(self.root, installed_paths=paths)
        except (OSError, ValueError):
            catalogued = []  # A pinned catalogue that cannot be read never hides what is installed.
        self.pages["cgmlst"].set_library(grouped["cgmlst"], catalogued, unreadable)
        self.pages["mlst"].set_library(grouped["mlst"], (), unreadable)
        if migrated:
            moved = "; ".join(f"{Path(row['from']).name} → {row['to']}" for row in migrated)
            self.status.setText(f"Moved {len(migrated)} cgMLST scheme(s) out of the classical schemes "
                                f"folder and into your cgMLST library: {moved}")
        self.show_entry()

    def showEvent(self, event):
        # Reopened after a download, an import or a hand-copied folder: the list a
        # person comes back to must be the library as it is now. A rescan during a
        # download would only describe a folder that is still being written.
        super().showEvent(event)
        # Fitted here as well as in the constructor: a window that has not been
        # shown yet reports the primary screen, which on a two-monitor desk is
        # often not the one the application is on. Reported as a window too big
        # to fit, with the controls along its bottom edge off the screen.
        self._fit_to_screen()
        if not (self.worker and self.worker.isRunning()):
            self.refresh_library()

    # -- background work -------------------------------------------------------

    def _start(self, operation, completed):
        if self.worker and self.worker.isRunning():
            return
        self.worker = _ReferenceWorker(operation, self)
        self.worker.succeeded.connect(completed)
        self.worker.failed.connect(lambda text: self.status.setText(f"Could not complete: {text}"))
        self.worker.progress.connect(self._progress)
        self.worker.finished.connect(self._finished)
        for page in self.pages.values():
            page.search_button.setEnabled(False)
            page.provider.setEnabled(False)
        self.download_button.setEnabled(False)
        self.refresh_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self._started = time.monotonic()
        self._progress_at, self._progress_total = 0, 0
        self.progress.setRange(0, 0)
        self.worker.start()

    def _progress(self, current, total, text):
        # A cgMLST scheme is thousands of loci and gigabytes on disk, so a step
        # naming only the locus it just finished gives no sense of position. Say
        # how far through it is and how long is left, because the difference
        # between "working" and "stuck" is the whole question a user has here.
        if total != self._progress_total or current < self._progress_at:
            # A new counted stage. A download that has just spent twenty minutes
            # extracting must not divide that time by this stage's first locus
            # and promise hours: each stage is timed from its own start.
            self._started = time.monotonic()
            self._progress_total = total
        self._progress_at = current
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(current)
        self.status.setText(describe_progress(current, total, text,
                                              time.monotonic() - self._started))

    def _finished(self):
        for page in self.pages.values():
            page.search_button.setEnabled(True)
            page.provider.setEnabled(True)
        self.refresh_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        self.show_entry()
        if self._close_pending:
            super().reject()

    # -- searching -------------------------------------------------------------

    def search(self):
        query = self.query.text().strip()
        if not query:
            self.status.setText("Enter an organism name to search the online catalog.")
            return
        if self.min_loci.value() > self.max_loci.value():
            self.status.setText(f"The minimum {UNIT[self.page.kind]} count must not exceed the maximum.")
            return
        filters = {"min_loci": self.min_loci.value(), "max_loci": self.max_loci.value()}
        for page in self.pages.values():
            # One search, both libraries: the words a person typed filter what they
            # already have as well as what comes back, in whichever tab it lands in.
            if page.query.text() != query:
                page.query.setText(query)
        provider_key = self.provider.currentText()
        self.catalog = self.catalogs[provider_key]
        self._provider_key = provider_key
        self.status.setText(f"Searching {provider_key} for {query}…")
        self._start(lambda cancelled, progress: self.catalog.search_schemes(
            query, cancelled=cancelled, progress=progress, **filters), self._searched)

    def _searched(self, document):
        """File every result into the library it would actually install into.

        Which library a scheme belongs to is decided by cgmlst_schemes.is_gene_by_gene
        — the same call the installer makes when it chooses a destination — so the
        tab a scheme is offered in is the tab it will appear in once installed.
        """
        self.catalog_errors = document["errors"]
        provider_key = getattr(self, "_provider_key", "")
        added = {}
        for kind in ("cgmlst", "mlst"):
            rows = [online_row(entry, kind, provider_key) for entry in document["schemes"]
                    if (cgmlst_schemes.is_gene_by_gene(entry)) == (kind == "cgmlst")]
            added[kind] = len(self.pages[kind].add_online(rows))
        found = len(document["schemes"])
        self.status.setText(
            f"Found {found} scheme(s) across {document['organisms_searched']} matching organism "
            f"catalog(s): {added['cgmlst']} listed under cgMLST schemes and {added['mlst']} under "
            f"classical MLST schemes. {len(document['errors'])} catalog error(s).")
        if found and not added[self.page.kind]:
            # Everything found belongs to the other library; open the tab holding it
            # rather than leaving a search that looks as if it returned nothing.
            self.show_page("cgmlst" if added["cgmlst"] else "mlst")

    def show_entry(self):
        for page in self.pages.values():
            page.show_entry()
        row = self.page.selected_row()
        running = self.worker is not None and self.worker.isRunning()
        self.download_button.setEnabled(bool(row and row["downloadable"]) and not running)
        self.download_button.setText("Download under the provider's terms…" if row
                                     and row["requires_terms"] else "Download selected scheme")
        if self.catalog_errors:
            # A catalog that could not be read is not an empty catalog; the last
            # search's failures stay visible until the next one replaces them.
            self.notice.setPlainText(self.notice.toPlainText() + "\n\nCatalog errors:\n" + "\n".join(
                f"{item['organism']}: {item['error']}" for item in self.catalog_errors))

    # -- downloading -----------------------------------------------------------

    def download(self):
        row = self.page.selected_row()
        if row is None or not row["downloadable"]:
            return
        entry = dict(row["entry"])
        client = self.catalogs.get(row["provider_key"]) or self.catalog
        kwargs = {}
        if row["requires_terms"]:
            notice = entry.get("terms_notice") or getattr(client, "TERMS_NOTICE", "")
            url = entry.get("terms_url") or getattr(client, "TERMS_URL", "")
            choice = QMessageBox.question(
                self, "Reference data usage restrictions",
                f"{notice}\n\n{url}\n\nHave you reviewed the policy and confirmed your use is "
                "permitted?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if choice != QMessageBox.StandardButton.Yes:
                self.status.setText("Nothing was downloaded. The provider's terms were not confirmed.")
                return
            if isinstance(client, CGMLSTOrgCatalog):
                kwargs["terms_acknowledged"] = True
        self._pending_kind = row["kind"]
        self.status.setText(f"Downloading {row['title']}…")
        # The classical library is named for both: a gene-by-gene scheme is
        # redirected into the cgMLST library by the installer itself.
        self._start(lambda cancelled, progress: client.download_scheme(
            entry, self.root / "schemes", cancelled=cancelled, progress=progress, **kwargs),
            self._downloaded)

    def _downloaded(self, result):
        """Say what was installed, where it now is, and show it in its own library.

        The reported bug was not a failed download: it was a scheme that could not
        be found afterwards. So the folder is named, the library is rescanned, the
        matching tab is opened and the row is selected — and if the rescan does not
        list it, that is said plainly rather than left to be discovered later.
        """
        path = Path(result["path"])
        kind = getattr(self, "_pending_kind", "cgmlst")
        self.refresh_library()
        page = self.pages["cgmlst"] if self._in_cgmlst_library(path) else self.pages[kind]
        self.tabs.setCurrentWidget(page)
        found = page.select_path(path)
        if not found and page.query.text():
            page.query.clear()  # A leftover filter must never hide what was just installed.
            found = page.select_path(path)
        row = page.selected_row() if found else None
        verb = "Installed" if result["created"] else "Already installed"
        where = ("your cgMLST library" if self._in_cgmlst_library(path)
                 else "your classical MLST scheme library")
        title = row["title"] if row else Path(result["path"]).name
        self.status.setText(
            f"{verb}: {title} — {result['locus_count']} {UNIT[page.kind]} in {where}, at {path}."
            + ("" if found else " It is not listed above; open the folder to check it."))
        page.notice.setPlainText("\n".join(result["notes"]) or "Snapshot validated.")
        self.schemeInstalled.emit(result["path"])
        self.installed.emit(result["path"])
        # The status line at the foot of a long dialog is easy to miss after a
        # download that ran for several minutes, and "did that work?" is the
        # question the person is left with. Say it plainly, once.
        summary = QMessageBox(self)
        summary.setIcon(QMessageBox.Icon.Information)
        summary.setWindowTitle("Scheme installed")
        summary.setText(f"{verb}: {title}")
        summary.setInformativeText(
            f"{result['locus_count']} {UNIT[page.kind]} are now in {where}."
            + ("" if found else "\n\nIt is not listed above; open the library folder to check it."))
        summary.setDetailedText(f"{path}\n\n" + ("\n".join(result["notes"]) or "Snapshot validated."))
        summary.setStandardButtons(QMessageBox.StandardButton.Ok)
        summary.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        # open(), never exec(): exec() blocks this call until someone clicks, which
        # stalls a scripted install and would stop a batch dead between schemes.
        # The box is still modal to this dialog; it just does not own the thread.
        # The reference keeps it alive until it closes and deletes itself.
        self._summary_box = summary
        summary.open()

    def _in_cgmlst_library(self, path) -> bool:
        resolved = Path(path).resolve()
        return any(base.resolve() in resolved.parents or base.resolve() == resolved
                   for base in cgmlst_library_roots(self.root) if base.exists())

    # -- lifecycle -------------------------------------------------------------

    def cancel(self):
        if self.worker and self.worker.isRunning():
            self.status.setText("Cancelling after the current network read…")
            self.worker.cancel()

    def confirm_stop(self) -> bool:
        """Ask before closing stops a download that may have run for many minutes.

        Closing used to cancel silently, so a scheme that was most of the way
        through was lost to a click meant to get the window out of the way. What
        has already been fetched is kept and the download resumes from there, and
        saying so is the difference between an alarming question and an easy one.
        """
        answer = QMessageBox.question(
            self, "Stop this download?",
            "A scheme is still downloading.\n\n"
            "Closing this window stops it. What has already been downloaded is kept, "
            "and starting the same scheme again continues from where it stopped "
            "rather than beginning again.",
            QMessageBox.StandardButton.Cancel | QMessageBox.StandardButton.Close,
            QMessageBox.StandardButton.Cancel)
        return answer == QMessageBox.StandardButton.Close

    def reject(self):
        if self.worker and self.worker.isRunning():
            if not self.confirm_stop():
                return
            self._close_pending = True
            self.cancel()
        else:
            super().reject()

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            if not self.confirm_stop():
                event.ignore()
                return
            self._close_pending = True
            self.cancel()
            event.ignore()
        else:
            event.accept()
