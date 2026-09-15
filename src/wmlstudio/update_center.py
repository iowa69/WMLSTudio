"""One place for everything that can be installed or updated, with one button each.

Reference data reaches this application from several different places -- a species
panel that is downloaded, scheme libraries that are partly shipped and partly
fetched, AMR databases that belong to HYDRA, a characterization snapshot that is
staged into the build. Before this page a user had to know which menu owned which
one. Here they are one list: what is installed, which version or revision of it,
when that copy last changed on disk, whether an update is available, and the
button that installs or updates it.

The state comes from `provisioning.report`, which probes rather than guesses, and
every install runs through the window's own `launch_task`, so the progress bar and
Cancel behave exactly as they do for an analysis. Nothing here contacts a server
until a button is pressed.
"""

from __future__ import annotations

import datetime
import importlib
import inspect
from pathlib import Path

from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from wmlstudio import provisioning
from wmlstudio.widgets import button, label

INTRO = ("Everything WMLSTudio can install or update, in the order you meet it. Nothing is "
         "downloaded, and no server is contacted, until you press a button on this page. "
         "Your own sequences are never uploaded. “Install and update everything” shows you what "
         "is missing and what is stale before it fetches anything, and publishes new reference "
         "data beside what you already have.")

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


def _stamp(path) -> str:
    """When this copy last changed on disk, as a plain date. Never a guess."""
    try:
        when = Path(path).stat().st_mtime
    except (OSError, TypeError, ValueError):
        return ""
    return datetime.datetime.fromtimestamp(when).strftime("%Y-%m-%d")


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


def describe_version(item) -> str:
    """The version, revision or count this installation actually reports."""
    detail = item.get("detail") or {}
    key, parts = item.get("key", ""), []
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


def describe_update(item) -> str:
    """Whether an update is available, without inventing a check nobody ran."""
    detail = item.get("detail") or {}
    if not item.get("ready"):
        action = item.get("action") or {}
        return action.get("label") or "Install the missing reference data"
    installed = str(detail.get("source_revision") or "")
    expected = str(detail.get("expected_revision") or "")
    if installed and expected and installed != expected:
        return f"Update available: this build expects revision {expected}"
    catalogue = detail.get("catalogue") or {}
    if catalogue.get("available"):
        # Being installed is not being complete: the sets that are not here report
        # nothing, and that reads exactly like a clean isolate unless it is said.
        return (f"{len(catalogue['available'])} further reference sets are listed by name and can "
                "be chosen for download")
    return "Nothing is checked online until you press the button"


class UpdateCenter(QDialog):
    """The Update menu's page: installed state, versions, and one button each."""

    def __init__(self, window):
        super().__init__(window)
        self.window_ref = window
        self.report = None
        self.setWindowTitle("Updates · installed reference data")
        self.resize(1000, 560)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Updates and installed reference data", "title", True))
        layout.addWidget(label(INTRO, "muted", True))
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["What", "Installed", "Version · last changed", "Update", ""])
        self.table.verticalHeader().hide()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.table, 1)
        self.status = label("", "small", True)
        layout.addWidget(self.status)
        row = QHBoxLayout()
        self.check_button = button("Check what is installed", self.refresh, True)
        row.addWidget(self.check_button)
        self.everything_button = button("Install and update everything…", self.plan_everything)
        self.everything_button.setToolTip(
            "Checks what is installed, what is missing and what is stale, shows you the list, "
            "and only then installs and updates it. Reference data is published beside what you "
            "have, so an analysis already recorded keeps the snapshot it was run against.")
        row.addWidget(self.everything_button)
        row.addStretch()
        row.addWidget(button("Close", self.close))
        layout.addLayout(row)
        self.refresh()

    # --- reading the installed state ----------------------------------------
    def refresh(self):
        """Re-probe in the background: this touches the disk and runs no tools."""
        window = self.window_ref
        project = getattr(window, "project", None)
        selected = project.get_setting("hydra_database_root", "") if project is not None else ""
        data_root = str(getattr(window, "root", "") or "") or None
        scheme_paths = [str(path) for path in getattr(window, "scheme_paths", [])] or None

        def operation(cancelled, progress):
            # (done, total, message): the worker's own progress contract.
            progress(0, 1, "Checking what is installed on this computer…")
            # verify=False: the deep re-hash belongs to a deliberate check, not to
            # opening a menu. The probe sentence each row carries says which ran.
            return provisioning.report(data_root=data_root, scheme_paths=scheme_paths,
                                       hydra_database_root=selected, verify=False,
                                       cancelled=cancelled)

        launch = getattr(window, "launch_task", None)
        if callable(launch) and launch(operation, "provisioning", self.show_report):
            self.status.setText("Checking what is installed on this computer…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    def show_report(self, report):
        if not isinstance(report, dict) or "items" not in report:
            return
        self.report = report
        items = report["items"]
        self.table.setRowCount(len(items))
        for row, item in enumerate(items):
            required = "required" if item.get("required") else "optional"
            values = [f"{item['title']} ({required})",
                      INSTALLED_WORDS.get(item.get("state", ""), item.get("state", "")),
                      describe_version(item), describe_update(item)]
            for column, text in enumerate(values):
                cell = QTableWidgetItem(str(text))
                cell.setToolTip("\n\n".join(filter(None, [item.get("reason"), item.get("probe"),
                                                          item.get("consequence")])))
                self.table.setItem(row, column, cell)
            action = button(ACTION_WORDS.get(item.get("state", ""), "Install…"),
                            lambda checked=False, key=item["key"]: self.start(key))
            action.setEnabled(bool(self.route_for(item["key"]) or item.get("action")))
            self.table.setCellWidget(row, 4, action)
        self.table.resizeRowsToContents()
        self.status.setText(report.get("summary", ""))

    # --- install and update everything --------------------------------------
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

        def operation(cancelled, progress):
            progress(0, 1, "Checking what is installed, and which reference release is current…")
            return provisioning.installation_plan(data_root=data_root, selected=selected,
                                                  cancelled=cancelled)

        launch = getattr(window, "launch_task", None)
        if callable(launch) and launch(operation, "provisioning:plan", self.confirm_everything):
            self.status.setText("Checking what is installed and what is missing…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    @staticmethod
    def plan_text(plan):
        """The plan as a person reads it: what will happen to each thing, and why."""
        lines = [plan["summary"], ""]
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
        if plan["release"].get("message"):
            lines.append(plan["release"]["message"])
        lines.append("Your own sequences are never uploaded. Nothing already installed is "
                     "overwritten: new reference data is published beside it.")
        return "\n".join(lines)

    def confirm_everything(self, plan):
        """Show the plan, then run it only if the person says yes."""
        if not isinstance(plan, dict) or "steps" not in plan:
            return False
        self.status.setText(plan["summary"])
        if not plan["work"]:
            QMessageBox.information(self, "Nothing to install", self.plan_text(plan))
            return False
        if QMessageBox.question(self, "Install and update everything",
                                self.plan_text(plan)) != QMessageBox.StandardButton.Yes:
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
            self.everything_button.setEnabled(False)
            self.status.setText("Installing and updating…")
            return True
        self.status.setText("Another background task is running. Try again when it finishes.")
        return False

    def everything_finished(self, result):
        """Say what landed and what did not; a failure here is not a silent one."""
        self.everything_button.setEnabled(True)
        if not isinstance(result, dict):
            return
        self.status.setText(result.get("summary", "") + " Press “Check what is installed” to "
                            "re-read the state from disk.")
        if result.get("failed"):
            QMessageBox.warning(self, "Some items were not installed", "\n\n".join(
                f"{failure['title']}: {failure['error']}\n{failure['state']}"
                for failure in result["failed"]))

    # --- installing ---------------------------------------------------------
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

    def start(self, key):
        """Install or update one thing, through the paths that already exist."""
        item = self.item_for(key) or {"key": key, "action": None}
        handler = self.route_for(key)
        if handler is not None:
            handler()
            self.status.setText("Press “Check what is installed” once that finishes, to see "
                                "the new state here.")
            return True
        action = item.get("action") or {}
        if action.get("automatic") and action.get("entry_point"):
            return self.run_entry_point(item, action)
        QMessageBox.information(
            self, item.get("title", "Not installed here"),
            "\n\n".join(filter(None, [action.get("label"), action.get("detail"),
                                      item.get("reason"),
                                      f"Command: {action['command']}" if action.get("command")
                                      else ""])) or
            "This part of the application is not installed, and nothing here can install it.")
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
            for key in ("path", "root", "species_count", "installed"):
                if result.get(key):
                    detail = f" {key.replace('_', ' ').capitalize()}: {result[key]}."
                    break
        self.status.setText(f"Finished.{detail} Press “Check what is installed” to re-read "
                            "the state from disk.")


def open_update_center(window):
    """Show the window's single Update Centre, creating it the first time."""
    existing = getattr(window, "update_center", None)
    if existing is None:
        existing = UpdateCenter(window)
        window.update_center = existing
    else:
        existing.refresh()
    existing.show()
    existing.raise_()
    return existing
