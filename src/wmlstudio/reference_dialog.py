"""Two scheme libraries — classical MLST and cgMLST — and the public catalogs behind them.

One tab per library, because a seven-locus MLST distance and a two-thousand-target
cgMLST distance are different quantities. Each tab lists what is installed, what is
catalogued and can be downloaded, and what a provider only permits under its own
terms; every row names the organism, the scheme, the size in that library's own
unit, the provider and the version. Nothing is downloaded until a row is selected
and its terms are shown.
"""

import threading
from pathlib import Path

from PySide6.QtCore import Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
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
from .reference_catalog import CGMLSTOrgCatalog, PasteurCatalog, PubMLSTCatalog
from .reference_index import scheme_library
from .sequence import AnalysisCancelled

KIND_TITLES = {"cgmlst": "cgMLST schemes", "mlst": "Classical MLST schemes"}
LIBRARY_NAMES = {"cgmlst": "cgMLST", "mlst": "classical MLST"}
UNIT = {"cgmlst": "targets", "mlst": "loci"}
SEPARATION = ("A classical seven-locus MLST distance and a 2,000-target cgMLST distance are different "
              "quantities: they never share a scale, an axis, a column or a threshold. That is why these "
              "are two libraries and two tabs, not one list.")
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
        lines.append(cgmlst_schemes.INTERPRETATION)
    return "\n\n".join(lines)


def _catalogued_detail(entry) -> str:
    """Everything the download button must show before a byte is transferred."""
    plan = cgmlst_schemes.download_plan(entry["key"])
    lines = [plan["title"],
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
            row = {"row_id": "path:" + str(resolved), "kind": kind, "state": "installed",
                   "status": "Installed · bundled" if packed else "Installed · your library",
                   "organism": entry["organism_label"], "scheme": entry["scheme_label"],
                   "count": entry["locus_count"], "target_set": source.get("target_set") or "",
                   "provider": entry["provider"], "version": entry["version"],
                   "title": entry["title"], "location": entry["path"], "path": entry["path"],
                   "key": source.get("catalog_key") or "", "entry": None,
                   "scheme_group": "", "requires_terms": False, "downloadable": False,
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
            "count": count, "target_set": "", "provider": str(entry.get("provider") or provider_key),
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
            # Requirement of the science, not of the layout: a core target set and an
            # accessory set are different quantities, so where a provider publishes
            # both the choice is explicit, and where it publishes one that is said in
            # words instead of shown as an empty menu.
            variants = QHBoxLayout()
            variants.addWidget(QLabel("Target set:"))
            self.variant = QComboBox()
            self.variant.setMinimumWidth(320)
            self.variant.activated.connect(self.choose_variant)
            variants.addWidget(self.variant)
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
                values.append({"core": "Core", "accessory": "Accessory"}.get(row["target_set"],
                                                                            "Not recorded"))
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

    # -- core versus accessory -------------------------------------------------

    def refresh_variants(self, row):
        """Offer the core/accessory choice where both exist; say so plainly where not.

        The choice is scoped to ONE provider's scheme family, because core and
        accessory are two target sets one provider publishes for one organism.
        Another provider's scheme for the same organism is a different scheme
        altogether, not the accessory half of this one, so it is named in words
        rather than dropped into the same menu.
        """
        if self.variant is None:
            return
        self.variant.blockSignals(True)
        self.variant.clear()
        if row is None or not (row["scheme_group"] or row["organism"]):
            self.variant.setEnabled(False)
            self.variant_note.setText(
                "Select a scheme to see whether a core and an accessory target set are both "
                "catalogued for it." if row is None else
                "This scheme records no organism, so no core or accessory pairing can be looked up "
                "for it. That is a gap in what the scheme records, not proof that none exists.")
            self.variant.blockSignals(False)
            return
        group = row["scheme_group"] or None
        variants = cgmlst_schemes.scheme_variants(row["organism"] or None, group=group)
        choices = [("Core target set", entry) for entry in variants["core"]]
        choices += [("Accessory / whole-genome set", entry) for entry in variants["accessory"]]
        for prefix, entry in choices:
            self.variant.addItem(f"{prefix} — {entry['title']}", entry["key"])
        self.variant.setEnabled(len(choices) > 1)
        index = self.variant.findData(row["key"])
        self.variant.setCurrentIndex(index if index >= 0 else 0)
        self.variant.blockSignals(False)
        note = variants["message"]
        if not choices:
            note = ("This scheme is not in the pinned catalogue, so no core or accessory pairing is "
                    "known for it. " + note)
        elif len(choices) == 1:
            note = f"{choices[0][0]} only. {note}"
        self.variant_note.setText(note + self._other_schemes(row, choices))

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
                            f"{entry['provider_name']}" for entry in others)
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
        self.setWindowTitle("Scheme libraries — classical MLST and cgMLST")
        self.resize(1080, 760)
        layout = QVBoxLayout(self)
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
        layout.addLayout(buttons)
        self.refresh_library()

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
        self.progress.setRange(0, 0)
        self.worker.start()

    def _progress(self, current, total, text):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(current)
        self.status.setText(text)

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

    def _in_cgmlst_library(self, path) -> bool:
        resolved = Path(path).resolve()
        return any(base.resolve() in resolved.parents or base.resolve() == resolved
                   for base in cgmlst_library_roots(self.root) if base.exists())

    # -- lifecycle -------------------------------------------------------------

    def cancel(self):
        if self.worker and self.worker.isRunning():
            self.status.setText("Cancelling after the current network read…")
            self.worker.cancel()

    def reject(self):
        if self.worker and self.worker.isRunning():
            self._close_pending = True
            self.cancel()
        else:
            super().reject()

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self._close_pending = True
            self.cancel()
            event.ignore()
        else:
            event.accept()
