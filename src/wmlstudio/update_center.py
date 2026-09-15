"""The Update tab: everything that can be installed or updated, with one button each.

Reference data reaches this application from several different places -- a species
panel that is downloaded, scheme libraries that are partly shipped and partly
fetched, AMR reference sets that belong to HYDRA, a characterization snapshot that
is staged into the build. Before this page a user had to know which menu owned
which one. Here they are one list: what is installed, which version or revision of
it, when that copy last changed on disk, how much room it takes or would cost to
fetch, whether an update is published, and the button that installs or updates it.

Three separations are the whole point of the page, because collapsing them is what
the user actually hit ("when I update instead of looking online it performs
analysis"):

* **Reading this computer** is not **asking a provider.** The rescan walks local
  folders and contacts nothing; the check makes one small network request for a
  version string. They are different buttons with different words, and a row says
  which of the two answered it.
* **A failed check is not an up-to-date result.** When the request does not get
  through, every row that depended on it says so and names the last check that did
  reach a provider, rather than falling back to silence that reads as "current".
* **A licence is not ours to accept.** A reference set whose provider's terms have
  to be read is reported as needing a decision, is never swept into "install and
  update everything", and its own button shows those terms and stops.

The state comes from `provisioning.report`, which probes rather than guesses. The
read-only probes run on this page's own background thread, so opening the tab
never raises the analysis progress dialog over the workspace; anything that
changes this computer runs through the window's own `launch_task`, where the
progress bar and Cancel behave exactly as they do for an analysis.
"""

from __future__ import annotations

import datetime
import importlib
import inspect
import json
import re
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import provisioning
from wmlstudio.background import FunctionWorker
from wmlstudio.ui_common import FlowLayout
from wmlstudio.widgets import button, label

INTRO = ("Everything WMLSTudio can install or update, in the order you meet it. Reading this "
         "computer and asking a provider are two different buttons: nothing is downloaded, and "
         "no server is contacted, until you press one that says so, and your own sequences are "
         "never uploaded.")

# How the four probe states read to somebody deciding whether to press the button.
INSTALLED_WORDS = {"ready": "Installed", "partial": "Partly installed",
                   "missing": "Not installed", "unusable": "Installed, not usable"}

ACTION_WORDS = {"ready": "Update…", "partial": "Complete…", "missing": "Install…",
                "unusable": "Repair…"}

# The window methods that already know how to install each thing, most specific
# first. A name that does not exist on the window is skipped, so a later agent can
# add `install_cgmlst_schemes` (or any other name listed here) and this page will
# use it without being edited.
ROUTES = {
    "species_panel": ("install_species_panel",),
    "mlst_schemes": ("install_mlst_schemes", "open_reference_manager"),
    "cgmlst_schemes": ("install_cgmlst_schemes", "open_cgmlst_catalog", "open_reference_manager"),
    "hydra_database": ("open_hydra_databases", "open_amr_databases"),
    "characterization": ("install_characterization_references",),
    "assembly_runtime": (),
    "blast_tools": (),
    "hydra_engine": (),
}

# Where the menu sends somebody for each entry, when there is a dedicated place.
MENU_ENTRIES = (
    ("species_panel", "Species reference panel…"),
    ("mlst_schemes", "MLST schemes…"),
    ("cgmlst_schemes", "cgMLST schemes…"),
    ("characterization", "Characterization references…"),
    ("hydra_database", "HYDRA databases…"),
)

#: What each column is really answering, for the header that has no room to say it.
COLUMN_NOTES = (
    "What this is, and whether an investigation can run without it.",
    "What this computer holds right now, read from its own folders.",
    "The version or revision this copy reports, and the day it last changed on disk.",
    "Room it takes here once installed, or what the download costs. Measured where it can be, "
    "and reported as unpublished rather than guessed where it cannot.",
    "Whether an update exists — and which question was answered to find out: this computer, "
    "what this build expects, or a provider over the internet.",
    "Install, update, or read the provider's terms first.",
)

#: What "an update" even means for something nobody downloads on its own. A part
#: that ships inside the application changes when the application does, and a
#: scheme catalogue is browsed one scheme at a time rather than version-compared;
#: neither may borrow the sentence that belongs to a set with a published release.
READY_WORDS = {
    "assembly_runtime": "Ships beside the application. It changes when WMLSTudio itself is "
                        "updated, not from here.",
    "blast_tools": "Ships beside the application. It changes when WMLSTudio itself is updated, "
                   "not from here.",
    "hydra_engine": "Ships inside this build. It changes when WMLSTudio itself is updated, not "
                    "from here.",
    "mlst_schemes": "Scheme catalogues are browsed and installed one scheme at a time; nothing "
                    "here compares your installed schemes against a published list.",
    "cgmlst_schemes": "Scheme catalogues are browsed and installed one scheme at a time; nothing "
                      "here compares your installed schemes against a published list.",
}

#: The reference sets whose provider publishes a release string this application
#: can compare against. Mirrors provisioning._database_steps deliberately: for
#: every other set, age is the only signal there is, and "we cannot tell whether
#: this is current" must never be written down as "it is".
PROVIDER_CHECKED = ("ncbi", "protein")

#: How many files a size measurement visits before it stops and says so. A cgMLST
#: library is hundreds of thousands of small files, and a page that has to stay
#: responsive is allowed to report "at least" rather than walk all of them.
MEASURE_LIMIT = 20000

#: A download size an item's own action already states, e.g. "about 120 MB". Read
#: from that sentence rather than repeated here, so the two can never disagree.
SIZE_PATTERN = re.compile(r"about\s+([\d.,]+)\s*(KB|MB|GB)", re.IGNORECASE)

#: Where "a provider answered us, on this day" is remembered. Beside the data
#: rather than in the project: it is a fact about this installation, not about one
#: investigation, and a failed check has to be able to name it.
CHECK_FILE = "Updates.ini"
CHECK_KEY = "last_successful_check"


def _stamp(path) -> str:
    """When this copy last changed on disk, as a plain date. Never a guess."""
    if not path:
        # Path("") is the working directory, and stamping *that* gave a thing
        # nobody has installed a plausible "last changed" date of today.
        return ""
    try:
        when = Path(path).stat().st_mtime
    except (OSError, TypeError, ValueError):
        return ""
    return datetime.datetime.fromtimestamp(when).strftime("%Y-%m-%d")


def _now() -> str:
    """This computer's clock, to the minute, in the words the page shows."""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def check_record(root, value=None):
    """Recall -- or, with a value, remember -- when a check last reached a provider.

    A failed check must be able to say when the answer it is still showing was
    last confirmed, and that has to outlive the window, so it is kept beside the
    data folder rather than in memory.
    """
    store = QSettings(str(Path(root) / CHECK_FILE), QSettings.Format.IniFormat)
    if value is not None:
        store.setValue(CHECK_KEY, json.dumps(value))
        store.sync()
        return value
    try:
        remembered = json.loads(store.value(CHECK_KEY, "") or "null")
    except (TypeError, ValueError):
        return None
    return remembered if isinstance(remembered, dict) else None


def measure_tree(path, *, limit=MEASURE_LIMIT):
    """Bytes under a folder, and whether every file in it was counted.

    Returns (bytes, complete). A measurement that stopped at the limit is reported
    as "at least" by the caller: an undercount announced as one is honest, an
    undercount presented as the size is not.
    """
    root = Path(path) if path else None
    if root is None or not root.exists():
        return None, True
    if root.is_file():
        try:
            return root.stat().st_size, True
        except OSError:
            return None, True
    total = counted = 0
    try:
        for item in root.rglob("*"):
            if counted >= limit:
                return total, False
            if item.is_file():
                counted += 1
                total += item.stat().st_size
    except OSError:
        # A folder that cannot be walked has an unknown size, not a size of zero.
        return (total or None), False
    return total, True


def human_size(count) -> str:
    """A byte count as somebody would say it out loud."""
    if not count:
        return ""
    for unit, size in (("GB", 1024 ** 3), ("MB", 1024 ** 2), ("KB", 1024)):
        if count >= size:
            value = count / size
            return f"{value:.1f} {unit}" if value < 10 else f"{value:.0f} {unit}"
    return f"{int(count)} bytes"


def installed_path(item) -> str:
    """The folder an item actually occupies, when the probe reported one."""
    detail = item.get("detail") or {}
    for key in ("path", "root", "roots"):
        value = detail.get(key)
        if isinstance(value, str) and value:
            return value
        if isinstance(value, list) and value:
            return str(value[0])
    return ""


def published_size(item) -> str:
    """The download size this item's own action states, when it states one."""
    action = item.get("action") or {}
    match = SIZE_PATTERN.search(str(action.get("detail") or ""))
    return f"{match.group(1)} {match.group(2).upper()}" if match else ""


def describe_size(item) -> str:
    """How much room this takes here, or what fetching it would cost. Measured, not guessed."""
    measured, complete = item.get("size_bytes"), item.get("size_complete", True)
    if measured and complete:
        return f"{human_size(measured)} here"
    if measured:
        return f"at least {human_size(measured)} here"
    published = published_size(item)
    if published:
        return f"about {published} to download"
    return "Not published"


def size_note(item) -> str:
    """Where that number came from, for the cell nobody has room to explain in."""
    measured, complete = item.get("size_bytes"), item.get("size_complete", True)
    if measured and complete:
        return "Measured on this computer, in the folder this item occupies."
    if measured:
        return ("Measured on this computer and stopped early: this folder holds more files than "
                f"one pass counts ({MEASURE_LIMIT:,}), so the real size is larger than the number "
                "shown.")
    if published_size(item):
        return "The download size this item's own installer states. Nothing has been measured."
    return ("Nobody publishes a size for this and nothing is installed here to measure, so no "
            "number is shown rather than a guessed one.")


def describe_version(item) -> str:
    """The version, revision or count this installation actually reports."""
    detail = item.get("detail") or {}
    key, parts = item.get("key", ""), []
    if item.get("state") == "missing":
        # The probe still reports the revision and the reference count this build
        # *expects*, and printing those beside "Not installed" read as though they
        # were here. What is expected is said as expected, or not at all.
        expected = detail.get("expected_revision")
        return ("Nothing installed"
                + (f" · this build expects revision {expected}" if expected else ""))
    if key == "species_panel":
        revision = detail.get("source_revision") or detail.get("expected_revision")
        if revision:
            parts.append(f"revision {revision}")
        if detail.get("reference_count"):
            parts.append(f"{detail['reference_count']} reference genomes")
    elif key in {"mlst_schemes", "cgmlst_schemes"}:
        if detail.get("count"):
            parts.append(f"{detail['count']} schemes")
            parts.append(f"{detail.get('installed_by_you', 0)} added by you")
    elif key == "hydra_database":
        catalogue = detail.get("catalogue") or {}
        if catalogue.get("total"):
            parts.append(f"{len(catalogue.get('installed') or [])} of {catalogue['total']} "
                         "reference sets installed")
        for name, version in sorted((detail.get("databases") or {}).items()):
            parts.append(f"{name} {version}")
    elif key == "characterization":
        if detail.get("manifest_format_version") is not None:
            parts.append(f"manifest format {detail['manifest_format_version']}")
        if detail.get("reference_digest"):
            parts.append(f"digest {str(detail['reference_digest'])[:8]}")
    elif key == "blast_tools":
        if detail.get("blastn_version"):
            parts.append(str(detail["blastn_version"]))
    elif key == "hydra_engine":
        for name in ("hydra_version", "version"):
            if detail.get(name):
                parts.append(str(detail[name]))
                break
    if not parts:
        for name in ("version", "revision"):
            if detail.get(name):
                parts.append(str(detail[name]))
    when = _stamp(installed_path(item))
    if when:
        parts.append(f"last changed {when}")
    return " · ".join(parts) or "Not recorded"


def failed_check_sentence(check, subject="a newer release") -> str:
    """A check that did not get through, said so that nobody reads it as "current"."""
    return (f"The check failed, so whether {subject} exists is unknown — "
            + str(check.get("error") or "the provider could not be reached") + " "
            + (f"A check last reached a provider on {check['last_success']}."
               if check.get("last_success") else
               "No check has ever reached a provider from this computer."))


def describe_update(item, check=None) -> str:
    """Whether an update exists for one requirement, and which question was answered.

    Being absent is read from this computer; being out of date can only come from
    whoever publishes the thing. The two are never written as one another, and a
    check that failed is never reported as nothing to do.
    """
    detail = item.get("detail") or {}
    if not item.get("ready"):
        action = item.get("action") or {}
        return action.get("label") or "Install the missing reference data"
    installed = str(detail.get("source_revision") or "")
    expected = str(detail.get("expected_revision") or "")
    if installed and expected and installed != expected:
        return (f"Update available: this build expects revision {expected}, and revision "
                f"{installed} is installed. Compared with this build, not with a provider.")
    if item.get("key") == "hydra_database":
        if not check:
            return "Installed. Not yet checked against what the provider publishes today."
        if not check.get("ok"):
            return failed_check_sentence(check, "a newer reference release")
        return str(check.get("message") or "Checked against the release the provider publishes.")
    if installed and expected:
        return (f"Matches revision {expected}, the one this build expects. No provider was asked; "
                "this is a comparison with what this build ships.")
    if item.get("key") in READY_WORDS:
        return READY_WORDS[item["key"]]
    return ("No version this application can compare against is published, so whether it is "
            "current cannot be confirmed here.")


def describe_set_update(row, step=None, check=None) -> str:
    """Whether one HYDRA reference set is current, and what its licence means for that."""
    name, provider = row.get("name", ""), row.get("provider") or "the provider"
    if not row.get("installed") and not row.get("open_licence"):
        # Reported as a decision, never quietly left out of the list: an absent set
        # reports nothing, and nothing reads exactly like a clean isolate.
        return (f"Needs your decision: {provider} publishes this under {row.get('licence')}, and "
                "this application will not accept those terms for you. It is never part of "
                "“install and update everything”; its own button shows the terms and asks.")
    if not row.get("installed") and row.get("download") != "automatic":
        return ("Not installed, and this engine cannot fetch it: install it by hand from "
                + (row.get("url") or provider) + ".")
    if not row.get("installed"):
        return f"Not installed, so nothing is reported for {row.get('purpose', 'this set')}."
    if check and not check.get("ok"):
        return failed_check_sentence(check, "a newer release of this set")
    if not check:
        return "Installed. Not yet checked against what the provider publishes today."
    if name in PROVIDER_CHECKED and step is not None and step.get("reason"):
        return str(step["reason"])
    if step is not None and step.get("action") in ("install", "update") and step.get("reason"):
        return str(step["reason"])
    return (f"Installed{', release ' + row['version'] if row.get('version') else ''}. {provider} "
            "publishes no version this application can compare against, so being current cannot "
            "be confirmed here.")


def decisions_needed(plan) -> list[dict]:
    """The reference sets left out of an automatic install because of their terms.

    "Install and update everything" fetches the core sets and whatever is already
    here. A set whose provider's terms have to be read is not among them, and this
    is the list that says so by name instead of letting it vanish from the plan.
    """
    rows = {row["name"]: row for row in (plan.get("catalogue") or {}).get("entries", [])}
    fetching = set(plan.get("databases") or ())
    waiting = []
    for step in plan.get("steps", []):
        row = rows.get(step.get("name", "")) if step.get("kind") == "database" else None
        if row is None or row.get("installed") or row.get("open_licence"):
            continue
        if row["name"] in fetching:
            # Someone asked for this one by name, so it is not being left out and
            # saying that it is would be its own kind of untruth.
            continue
        waiting.append({"name": row["name"], "title": row["title"], "licence": row["licence"],
                        "provider": row.get("provider", ""), "url": row.get("url", ""),
                        "size": row.get("size", ""), "purpose": row.get("purpose", "")})
    return waiting


def probe(data_root=None, scheme_paths=None, selected=None, *, cancelled=None, progress=None):
    """Read this computer: what is installed, how big it is. Contacts nothing.

    The size is measured rather than published wherever something is installed,
    because "how much room would this cost me" is the question a laptop user is
    actually asking, and a number copied from a release note is not an answer
    about this machine.
    """
    if progress:
        progress(0, 3, "Reading what is installed on this computer…")
    data = provisioning.report(data_root=data_root, scheme_paths=scheme_paths,
                               hydra_database_root=selected, verify=False, cancelled=cancelled)
    if progress:
        progress(1, 3, "Measuring how much room the installed reference data takes…")
    for item in data["items"]:
        if item.get("state") == "missing":
            # Nothing of this is installed, so whatever bytes sit under the folder
            # it would use belong to something else -- the cgMLST library and the
            # MLST library share one root, and one of them being absent must not
            # be reported as the other one's size.
            item["size_bytes"], item["size_complete"] = None, True
            continue
        measured, complete = measure_tree(installed_path(item))
        item["size_bytes"], item["size_complete"] = measured, complete
    if progress:
        progress(2, 3, "Listing every reference set this engine can search…")
    data["catalogue"] = provisioning.hydra_database_catalogue(data_root, selected, measure=True)
    if progress:
        progress(3, 3, data["summary"])
    return data


class UpdateCenter(QWidget):
    """The Update tab: installed state, published state, size, and one button each."""

    def __init__(self, window):
        super().__init__(window)
        self.window_ref = window
        self.report = None
        self.catalogue = {}
        self.plan = None
        self.steps = {}
        self.check = None
        self.worker = None
        # True while the window's own worker is installing everything, so the
        # button cannot be pressed twice by a probe finishing underneath it.
        self.installing = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(label(INTRO, "muted", True))
        self.check_line = label("", "small", True)
        layout.addWidget(self.check_line)
        self.table = QTableWidget(0, 6)
        # Short headers, because a clipped one cannot be read at 1000 px; the
        # sentence each column is actually answering lives in its tooltip.
        self.table.setHorizontalHeaderLabels(
            ["What", "Installed", "Version · date", "Size", "Update", ""])
        for column, explanation in enumerate(COLUMN_NOTES):
            header_item = QTableWidgetItem(self.table.horizontalHeaderItem(column).text())
            header_item.setToolTip(explanation)
            self.table.setHorizontalHeaderItem(column, header_item)
        self.table.verticalHeader().hide()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        for column in (1, 3, 5):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.table, 1)
        self.status = label("", "small", True)
        layout.addWidget(self.status)
        row = FlowLayout()
        # Two buttons, two words: one walks this computer's folders, the other
        # makes a network request. Conflating them is the bug this page exists for.
        self.rescan_button = button("Rescan what is installed here", self.refresh)
        self.rescan_button.setToolTip(
            "Reads this computer's own folders and nothing else. No server is contacted, nothing "
            "is downloaded and no analysis is run.")
        row.addWidget(self.rescan_button)
        self.check_button = button("Check online for updates", self.check_online, True)
        self.check_button.setToolTip(
            "Asks the providers, over the internet, which release they publish today. One small "
            "request for a version string; none of your sequences is sent, and nothing is "
            "downloaded or changed by asking.")
        row.addWidget(self.check_button)
        self.everything_button = button("Install and update everything…", self.plan_everything)
        self.everything_button.setToolTip(
            "Checks online, then shows you what is missing, what is stale and what is left alone "
            "and why, and only then installs and updates it. Reference data is published beside "
            "what you have, so an analysis already recorded keeps the snapshot it was run "
            "against. A set whose provider's terms have to be accepted is never included.")
        row.addWidget(self.everything_button)
        layout.addLayout(row)
        self.show_check(None)

    # --- asking, telling, and the confirmations a test can answer -------------
    def ask(self, title, text) -> bool:
        """Every yes/no on this page goes through here, so a test can answer it."""
        return QMessageBox.question(self, title, text) == QMessageBox.StandardButton.Yes

    def tell(self, title, text) -> None:
        """Every "here is how that is installed" goes through here, for the same reason."""
        QMessageBox.information(self, title, text)

    def confirm_terms(self, row) -> bool:
        """Show a provider's terms and stop. Nothing here accepts them for anybody."""
        return self.ask(
            f"{row.get('title', 'Reference set')} — the provider's own terms",
            "\n\n".join(filter(None, [
                f"{row.get('provider') or 'The provider'} publishes this reference set under: "
                f"{row.get('licence') or 'terms it does not record'}.",
                row.get("licence_note"),
                f"Terms and source: {row['url']}" if row.get("url") else "",
                "WMLSTudio cannot accept these terms for you, and never does so automatically. "
                "Download it only if your own use is allowed under them.",
                "Download this set now?"])))

    # --- reading this computer -----------------------------------------------
    def first_look(self):
        """Read the disk the first time this tab is opened, and only then.

        Opening the tab must not look like starting an analysis, so this runs on
        the page's own thread and never asks a provider anything. Re-reading is a
        button, not something the tab does behind a person's back.
        """
        return self.refresh() if self.report is None else False

    def refresh(self):
        """Re-probe this computer. Touches the disk, runs no tools, contacts nothing."""
        window = self.window_ref
        project = getattr(window, "project", None)
        selected = project.get_setting("hydra_database_root", "") if project is not None else ""
        data_root = str(getattr(window, "root", "") or "") or None
        scheme_paths = [str(path) for path in getattr(window, "scheme_paths", [])] or None

        def operation(cancelled, report_progress):
            # verify=False: the deep re-hash belongs to a deliberate check, not to
            # opening a tab. The probe sentence each row carries says which ran.
            return probe(data_root, scheme_paths, selected,
                         cancelled=cancelled, progress=report_progress)

        return self.run(operation, self.show_report,
                        "Reading what is installed on this computer. No server is contacted.")

    def run(self, operation, finished, message):
        """Run one read-only probe on this page's own thread.

        Deliberately not the window's `launch_task`: that raises the modal
        "analysis in progress" dialog over the workspace, which is how a tab that
        only reads a few folders came to look like it was analysing something.
        Anything that writes to this computer still goes through `launch_task`.
        """
        if self.worker is not None and self.worker.isRunning():
            self.status.setText("This page is already reading. One at a time.")
            return False
        self.worker = FunctionWorker(operation, self)
        self.worker.completed.connect(finished)
        self.worker.failed.connect(self.probe_failed)
        self.worker.progress.connect(lambda _percent, text: self.status.setText(str(text)))
        self.worker.finished.connect(self.probe_finished)
        for control in (self.rescan_button, self.check_button, self.everything_button):
            control.setEnabled(False)
        self.status.setText(message)
        self.worker.start()
        return True

    def probe_finished(self):
        for control in (self.rescan_button, self.check_button):
            control.setEnabled(True)
        # A probe finishing must not hand back a button an install is still using:
        # the plan arrives on this same signal, one step before the install starts.
        self.everything_button.setEnabled(not self.installing)

    def probe_failed(self, message):
        """A probe that could not finish says so here, and claims nothing about updates."""
        self.status.setText(f"This computer could not be read in full: {message}")

    def stop(self):
        """Let go of the background probe, so nothing outlives the window."""
        worker = self.worker
        if worker is not None and worker.isRunning():
            worker.cancel()
            worker.wait(5000)
        return True

    def show_report(self, report):
        if not isinstance(report, dict) or "items" not in report:
            return False
        self.report = report
        self.catalogue = report.get("catalogue") or {}
        self.render()
        self.status.setText(report.get("summary", ""))
        return True

    def clear(self):
        """The tab's own Clear: forget what is shown, read this computer again.

        The record of when a check last reached a provider is deliberately kept:
        clearing a view must not erase the only thing a failed check has to say.
        """
        self.report, self.plan, self.steps, self.check = None, None, {}, None
        self.catalogue = {}
        self.table.setRowCount(0)
        self.show_check(None)
        self.status.setText("Cleared. Reading this computer again; nothing online was forgotten.")
        return self.refresh()

    # --- asking the providers -------------------------------------------------
    def check_online(self):
        """The one network action on this page, and it is always a deliberate press."""
        window = self.window_ref
        project = getattr(window, "project", None)
        selected = project.get_setting("hydra_database_root", "") if project is not None else ""
        data_root = str(getattr(window, "root", "") or "") or None

        def operation(cancelled, report_progress):
            report_progress(0, 1, "Asking the providers which release they publish today…")
            return provisioning.installation_plan(data_root, selected, check_online=True,
                                                  cancelled=cancelled)

        return self.run(operation, self.checked,
                        "Asking the providers which release they publish today. None of your "
                        "sequences is sent.")

    def checked(self, plan):
        """Record what the provider said — or that it said nothing at all."""
        if not isinstance(plan, dict) or "release" not in plan:
            return False
        release = plan.get("release") or {}
        remembered = self.remembered_check()
        when = _now()
        if release.get("error"):
            # A failed check keeps the previous answer visible and names its date;
            # it must never leave the page reading as though nothing needs doing.
            self.show_check({"ok": False, "when": when, "error": release["error"],
                             "last_success": (remembered or {}).get("when", "")})
            self.status.setText("The check did not get through. Nothing was downloaded and "
                                "nothing on this computer was changed.")
            self.render()
            return False
        self.plan = plan
        self.steps = {step["key"]: step for step in plan.get("steps", [])}
        record = {"when": when, "latest": release.get("latest", ""),
                  "installed": release.get("installed", "")}
        root = getattr(self.window_ref, "root", None)
        if root:
            check_record(root, record)
        self.show_check({"ok": True, "when": when, "error": "",
                         "message": release.get("message", ""),
                         "latest": release.get("latest", ""),
                         "last_success": when})
        self.status.setText(plan.get("summary", ""))
        self.render()
        return True

    def remembered_check(self):
        """The last check that reached a provider, from a previous session if need be."""
        root = getattr(self.window_ref, "root", None)
        return check_record(root) if root else None

    def show_check(self, check):
        """One line at the top of the page: which question was last answered, and when."""
        self.check = check
        remembered = self.remembered_check() or {}
        if check is None:
            last = remembered.get("when", "")
            self.check_line.setText(
                "Nothing here has been compared with a provider in this session"
                + (f"; a check last reached one on {last}. " if last else
                   ". No check has ever reached a provider from this computer. ")
                + "Not checked is not the same as up to date — press “Check online for updates”.")
        elif not check.get("ok"):
            self.check_line.setText(
                f"The check on {check.get('when', 'today')} did not get through: "
                + str(check.get("error") or "the provider could not be reached") + " "
                + (f"A check last reached a provider on {check['last_success']}; nothing below "
                   "has been compared with one since."
                   if check.get("last_success") else
                   "No check has ever reached a provider from this computer.")
                + " A failed check is not an up-to-date result.")
        else:
            self.check_line.setText(
                f"Checked online {check.get('when', '')}. "
                + str(check.get("message") or "")
                + " Rows below that no provider publishes a version for say so instead.")
        return self.check_line.text()

    # --- the list itself ------------------------------------------------------
    def rows(self) -> list[dict]:
        """One row per installable thing: the requirements, then every reference set.

        The reference sets are listed under the database requirement rather than
        folded into it, because "a database is installed" says nothing about which
        of fifteen sets is there, and an absent set reports nothing at all.
        """
        listing = []
        for item in (self.report or {}).get("items", []):
            required = "required" if item.get("required") else "optional"
            listing.append({
                "key": item["key"], "title": f"{item['title']} ({required})",
                "installed": INSTALLED_WORDS.get(item.get("state", ""), item.get("state", "")),
                "version": describe_version(item), "size": describe_size(item),
                "size_note": size_note(item),
                "update": describe_update(item, self.check),
                "action": ACTION_WORDS.get(item.get("state", ""), "Install…"),
                "enabled": bool(self.route_for(item["key"]) or item.get("action")),
                "tooltip": "\n\n".join(filter(None, [item.get("reason"), item.get("probe"),
                                                     item.get("consequence")])),
            })
            if item["key"] == "hydra_database":
                listing.extend(self.database_rows())
        return listing

    def database_rows(self) -> list[dict]:
        """One row per HYDRA reference set, with its size, its licence and its state."""
        listing = []
        for row in self.catalogue.get("entries", []):
            name = row["name"]
            installed = ("Bundled with this release" if row.get("bundled") else
                         "Installed" if row.get("installed") else "Not installed")
            version = " · ".join(filter(None, [
                f"release {row['version']}" if row.get("version") else "",
                f"staged {str(row['staged'])[:10]}" if row.get("staged") else ""]))
            decision = not row.get("installed") and not row.get("open_licence")
            # The same distinction the requirement rows keep: room taken here is
            # measured, room a download would cost is what the provider publishes.
            size = (f"{row['size']} here" if row.get("installed") and row.get("bytes")
                    else f"{row['size']} to download" if row.get("size")
                    else "Not published")
            listing.append({
                "key": "database:" + name,
                "title": f"    ↳ {row['title']} (reference set)",
                "installed": installed,
                "version": version or ("Not installed" if not row.get("installed")
                                       else "Not recorded"),
                "size": size,
                "size_note": row.get("size_basis") or ("Nobody publishes a size for this set; it "
                                                       "is measured once it is installed."),
                "update": describe_set_update(row, self.steps.get("database:" + name), self.check),
                "action": ("Read the terms…" if decision else
                           "Update…" if row.get("installed") else
                           "Install…" if row.get("download") == "automatic" else
                           "How to install…"),
                "enabled": True,
                "tooltip": "\n\n".join(filter(None, [
                    row.get("purpose", ""), f"Published by {row.get('provider') or 'not recorded'}",
                    f"Licence: {row.get('licence', '')}", row.get("licence_note", ""),
                    row.get("size_basis", ""), row.get("url", "")])),
            })
        return listing

    def render(self):
        """Draw the list as it stands. Reads no disk and asks nobody anything."""
        listing = self.rows()
        self.table.setRowCount(len(listing))
        for index, row in enumerate(listing):
            for column, text in enumerate([row["title"], row["installed"], row["version"],
                                           row["size"], row["update"]]):
                cell = QTableWidgetItem(str(text))
                # A size has no room to explain itself in a cell, so its own cell
                # carries where the number came from rather than the row's reason.
                cell.setToolTip(row["size_note"] if column == 3
                                else row["tooltip"] or str(text))
                self.table.setItem(index, column, cell)
            action = button(row["action"],
                            lambda checked=False, key=row["key"]: self.start(key))
            action.setEnabled(bool(row["enabled"]))
            self.table.setCellWidget(index, 5, action)
        self.table.resizeRowsToContents()
        return listing

    # --- install and update everything ---------------------------------------
    def plan_everything(self):
        """Work out what one press would do, before anything is downloaded.

        Two steps on purpose. "Update everything" is otherwise a leap of faith:
        this reads the disk, asks NCBI for one version string and comes back with
        a list — what is missing, what is stale, what is left alone and why — and
        nothing is fetched until that list has been agreed to.
        """
        window = self.window_ref
        project = getattr(window, "project", None)
        selected = project.get_setting("hydra_database_root", "") if project is not None else ""
        data_root = str(getattr(window, "root", "") or "") or None

        def operation(cancelled, report_progress):
            report_progress(0, 1, "Checking what is installed, and which release is current…")
            return provisioning.installation_plan(data_root, selected, cancelled=cancelled)

        return self.run(operation, self.confirm_everything,
                        "Checking what is installed and what is missing…")

    @staticmethod
    def plan_text(plan):
        """The plan as a person reads it: what will happen to each thing, and why."""
        lines = [plan["summary"], ""]
        if plan["release"].get("error"):
            lines.insert(0, "The online check did not get through, so what follows is not a "
                            "statement that you are up to date.\n")
        for action, heading in (("install", "Will be installed"), ("update", "Will be updated"),
                                ("by hand", "Cannot be installed from here")):
            rows = [step for step in plan["steps"] if step["action"] == action]
            if not rows:
                continue
            lines.append(heading + ":")
            lines.extend(f"  • {step['title']}"
                         + (f" ({step['size']})" if step.get("size") else "")
                         + (f" — {step['licence']}" if step.get("licence") else "")
                         for step in rows)
            lines.append("")
        waiting = decisions_needed(plan)
        if waiting:
            lines.append("Needs your decision, and is not included here:")
            lines.extend(f"  • {row['title']} — {row['provider'] or 'the provider'} publishes it "
                         f"under {row['licence']}. This application will not accept those terms "
                         "for you; install it from its own row once you have read them."
                         for row in waiting)
            lines.append("")
        if plan["release"].get("message"):
            lines.append(plan["release"]["message"])
        lines.append("Your own sequences are never uploaded. Nothing already installed is "
                     "overwritten: new reference data is published beside it.")
        return "\n".join(lines)

    def confirm_everything(self, plan):
        """Show the plan, then run it only if the person says yes."""
        if not isinstance(plan, dict) or "steps" not in plan:
            return False
        self.checked(plan)
        self.status.setText(plan["summary"])
        if not plan["work"]:
            self.tell("Nothing to install", self.plan_text(plan))
            return False
        if not self.ask("Install and update everything", self.plan_text(plan)):
            self.status.setText("Nothing was downloaded. " + plan["summary"])
            return False
        window = self.window_ref
        project = getattr(window, "project", None)
        data_root = str(getattr(window, "root", "") or "") or None

        def operation(cancelled, progress):
            return provisioning.install_everything(data_root=data_root, project=project,
                                                   plan=plan, cancelled=cancelled,
                                                   progress=progress)

        launch = getattr(window, "launch_task", None)
        if callable(launch) and launch(operation, "update:everything", self.everything_finished):
            self.installing = True
            self.everything_button.setEnabled(False)
            self.status.setText("Installing and updating…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    def everything_finished(self, result):
        """Say what landed and what did not; a failure here is not a silent one."""
        self.installing = False
        self.everything_button.setEnabled(True)
        if not isinstance(result, dict):
            return
        waiting = decisions_needed(result.get("plan") or {})
        self.status.setText(result.get("summary", "")
                            + (f" {len(waiting)} reference set(s) were left out because their "
                               "provider's terms are yours to accept, not ours: "
                               + ", ".join(row["title"] for row in waiting) + "."
                               if waiting else "")
                            + " Press “Rescan what is installed here” to re-read this computer.")
        if result.get("failed"):
            self.tell("Some items were not installed", "\n\n".join(
                f"{failure['title']}: {failure['error']}\n{failure['state']}"
                for failure in result["failed"]))

    # --- installing one thing -------------------------------------------------
    def route_for(self, key):
        """The window method that already installs this thing, when there is one."""
        for name in ROUTES.get(key, ()):
            handler = getattr(self.window_ref, name, None)
            if callable(handler):
                return handler
        return None

    def item_for(self, key):
        for item in (self.report or {}).get("items", []):
            if item.get("key") == key:
                return item
        return None

    def set_for(self, name):
        for row in self.catalogue.get("entries", []):
            if row.get("name") == name:
                return row
        return None

    def start(self, key):
        """Install or update one thing, through the paths that already exist."""
        if str(key).startswith("database:"):
            return self.start_set(str(key).split(":", 1)[1])
        item = self.item_for(key) or {"key": key, "action": None}
        handler = self.route_for(key)
        if handler is not None:
            handler()
            self.status.setText("Press “Rescan what is installed here” once that finishes, to "
                                "see the new state.")
            return True
        if self.report is None:
            # Reached from the menu before the tab has read anything. "Nothing here
            # can install it" would be a claim about this computer that nobody has
            # looked at yet, so the page says what it actually knows.
            self.status.setText("This page has not read this computer yet, so it cannot say how "
                                "that one is installed. Press “Rescan what is installed here”.")
            return False
        action = item.get("action") or {}
        if action.get("automatic") and action.get("entry_point"):
            return self.run_entry_point(item, action)
        self.tell(item.get("title", "Not installed here"),
                  "\n\n".join(filter(None, [action.get("label"), action.get("detail"),
                                            item.get("reason"),
                                            f"Command: {action['command']}" if action.get("command")
                                            else ""])) or
                  "This part of the application is not installed, and nothing here can install it.")
        return False

    def start_set(self, name):
        """Install or update one HYDRA reference set, asking first where terms apply."""
        row = self.set_for(name)
        if row is None:
            self.status.setText(
                f"{name} is not in the catalogue this engine reads."
                if self.catalogue.get("entries") else
                "This page has not read this computer yet. Press “Rescan what is installed "
                "here” first.")
            return False
        if row.get("download") != "automatic":
            self.tell(row["title"],
                      "\n\n".join(filter(None, [
                          f"{row.get('provider') or 'The provider'} publishes no versioned "
                          "download this engine can fetch, so this set is installed by hand.",
                          f"Source: {row['url']}" if row.get("url") else "",
                          f"Licence: {row.get('licence', '')}", row.get("licence_note", ""),
                          "Until it is installed, nothing is reported for "
                          + row.get("purpose", "what it covers") + "."])))
            return False
        if not row.get("open_licence") and not self.confirm_terms(row):
            self.status.setText(
                f"Nothing was downloaded. {row['title']} needs its provider's terms accepted "
                "first, and this application will not accept them for you.")
            return False
        window = self.window_ref
        project = getattr(window, "project", None)
        root = Path(str(getattr(window, "root", "") or "."))
        download_root = root / "references" / "hydra"
        fetchable = set(self.catalogue.get("automatic") or ())
        # A new snapshot has to carry everything the store already holds, or the
        # next run silently loses whatever was left out of it.
        keep = sorted(({name} | set(self.catalogue.get("installed") or ())) & fetchable)

        def operation(cancelled, progress):
            return provisioning.update_hydra_databases(download_root, keep or [name],
                                                       project=project, cancelled=cancelled,
                                                       progress=progress)

        launch = getattr(window, "launch_task", None)
        if callable(launch) and launch(operation, "update:database:" + name, self.installed):
            self.status.setText(f"Downloading {row['title']} and everything the store already "
                                "holds, into a new snapshot beside it…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    def run_entry_point(self, item, action):
        """Run a provisioning installer on the existing worker, so Cancel works."""
        module_name, _, function_name = str(action["entry_point"]).partition(":")
        try:
            target = getattr(importlib.import_module(module_name), function_name)
            parameters = inspect.signature(target).parameters
        except (AttributeError, ImportError, TypeError, ValueError) as error:
            self.status.setText(f"That installer could not be started: {error}")
            return False
        if "cancelled" not in parameters or "progress" not in parameters:
            # Only the cancellable installers may run on the shared worker; the
            # others are configuration calls that belong to their own dialog.
            self.status.setText(action.get("detail") or "This one is installed from its own page.")
            return False
        argument = action.get("argument") or None

        def operation(cancelled, progress):
            if argument:
                return target(argument, cancelled=cancelled, progress=progress)
            return target(cancelled=cancelled, progress=progress)

        launch = getattr(self.window_ref, "launch_task", None)
        if callable(launch) and launch(operation, "update:" + item["key"], self.installed):
            self.status.setText(f"{item.get('title', 'Installing')}: working…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    def installed(self, result):
        """Say what landed, and leave re-reading the disk to a deliberate press."""
        detail = ""
        if isinstance(result, dict):
            for key in ("path", "root", "database_root", "species_count", "installed"):
                if result.get(key):
                    detail = f" {key.replace('_', ' ').capitalize()}: {result[key]}."
                    break
        self.status.setText(f"Finished.{detail} Press “Rescan what is installed here” to re-read "
                            "the state from disk.")


def open_update_center(window):
    """Show the window's single Update page, wherever the request came from.

    A menu entry, the command search and the tab itself all arrive here, and all
    three land on the same page: one place that says what is installed, not a tab
    and a dialog that can disagree with each other.
    """
    centre = getattr(window, "update_center", None)
    navigate = getattr(window, "navigate", None)
    if centre is not None and callable(navigate):
        navigate("update")
        centre.first_look()
    return centre
