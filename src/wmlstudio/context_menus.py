"""A reusable right-click system: what a click resolved to, and what may be done with it.

Two halves, kept apart on purpose:

* **Pure declaration.** `Selection` says what a right-click resolved to, `ActionSpec`
  says what an action is, and `plan_for` says which actions a view offers for that
  selection. None of it creates a widget, so it is testable without a menu.
* **Adapters.** One small class per kind of view turns a click point into a
  `Selection`, reading sample ids from the item data rather than the row number,
  because every table in this application is sortable and the visual row is not the
  data row.

Handler protocol
----------------
`ActionSpec.handler` is the *name* of a method on the main window, so this module
imports nothing from the window and can be unit tested on its own. The window
implements each one as::

    def context_add_to_report(self, selection: Selection) -> None: ...

Actions that carry a submenu receive the chosen entry as a second positional
argument::

    def context_export_selection(self, selection: Selection, fmt: str) -> None: ...
    def context_assign_organism_quick(self, selection: Selection,
                                      organism: tuple[str, str]) -> None: ...

A handler is responsible for its own confirmation and for refusing while a
background task runs; `action_state` marks such actions disabled with a reason,
it does not enforce them.

Two rules the wording must keep, because they are scientific claims:
archiving hides an isolate and destroys nothing, and an assigned organism is the
operator's record of identity, not a measurement this software made.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QModelIndex, QPoint, Qt
from PySide6.QtWidgets import (
    QGraphicsView,
    QTableView,
    QTableWidget,
    QTreeWidget,
    QWidget,
)

SEPARATOR = object()

USER_ROLE = Qt.ItemDataRole.UserRole

# Reasons shown as a tooltip on a disabled entry. An action that is meaningful but
# currently blocked is disabled with one of these; an action that is meaningless
# for the view is left out of the menu entirely.
REASONS = {
    "busy": "A background task is active. Cancel or finish it before changing your project.",
    "single": "Select exactly one isolate.",
    "no_file": "This isolate has no attached sequence file.",
    "bundled": "Bundled reference snapshots are read-only.",
    "running": "Wait for this sample's analysis to finish.",
    "none": "",
}

EXPORT_FORMATS = (
    ("CSV…", "csv"),
    ("TSV…", "tsv"),
    ("JSON…", "json"),
    ("HTML report…", "html"),
    ("PDF report…", "pdf"),
    ("Allelic profile table…", "profiles"),
    ("Evidence bundle (.json)…", "bundle"),
)


@dataclass(frozen=True)
class Selection:
    """What a right-click resolved to. Plain data; safe to hand to a worker."""

    view_id: str
    sample_ids: tuple[str, ...] = ()
    clicked_id: str | None = None
    folder: tuple[str, str] | None = None
    group_id: str | None = None
    paths: tuple[Path, ...] = ()
    on_blank: bool = False
    from_keyboard: bool = False

    @property
    def count(self) -> int:
        return len(self.sample_ids)

    @property
    def single(self) -> str | None:
        return self.sample_ids[0] if len(self.sample_ids) == 1 else None


@dataclass(frozen=True)
class ActionSpec:
    """One menu entry, declared without a widget and without a bound callable."""

    key: str
    title: str
    handler: str
    needs: str = "samples"          # samples | single | folder | paths | group | any
    group: int = 0                  # menu section; sections are separated automatically
    singular: str = ""              # wording used when the action targets one isolate
    destructive: bool = False
    mutates: bool = False           # writes to the project or the filesystem
    checkable: bool = False
    submenu: str = ""               # name of a submenu built at menu time
    enabled_when: object = None     # callable(window, selection) -> (bool, reason)
    views: frozenset[str] | None = None

    def format_title(self, selection) -> str:
        """Show the count whenever more than one isolate is affected."""
        count = getattr(selection, "count", 0)
        if "{n}" not in self.title:
            return self.title
        if count > 1:
            return self.title.replace("{n}", str(count))
        if self.singular:
            return self.singular
        return self.title.replace("{n} ", "").replace(" {n}", "")


def _scheme_removable(window, selection):
    """A scheme that ships with the application is read-only; an imported one is not."""
    root = getattr(window, "root", None)
    if root is None:
        return True, ""
    for path in selection.paths:
        if Path(root) not in Path(path).parents:
            return False, REASONS["bundled"]
    return True, ""


ACTIONS: dict[str, ActionSpec] = {spec.key: spec for spec in (
    # 0 — look at one thing
    ActionSpec("open_record", "Open isolate record…", "context_open_record", needs="single"),
    ActionSpec("open_folder", "Open containing folder", "context_open_folder", needs="single"),
    ActionSpec("open_scheme_folder", "Open scheme folder", "context_open_scheme_folder",
               needs="paths"),
    # 1 — send these isolates somewhere, explicitly
    ActionSpec("select_in_library", "Show {n} isolates in Isolate library",
               "context_select_in_library", group=1, singular="Show in Isolate library"),
    ActionSpec("add_to_comparison", "Add {n} isolates to comparison",
               "context_add_to_comparison", group=1, mutates=True,
               singular="Add to comparison"),
    ActionSpec("remove_from_comparison", "Remove {n} isolates from comparison",
               "context_remove_from_comparison", group=1, mutates=True,
               singular="Remove from comparison"),
    ActionSpec("add_to_characterization", "Add {n} isolates to characterization cohort",
               "context_add_to_characterization", group=1, mutates=True,
               singular="Add to characterization cohort"),
    ActionSpec("add_to_report", "Add {n} isolates to report", "context_add_to_report",
               group=1, mutates=True, singular="Add to report"),
    ActionSpec("remove_from_report", "Remove {n} isolates from report",
               "context_remove_from_report", group=1, mutates=True,
               singular="Remove from report"),
    # 2 — say what these isolates are
    ActionSpec("assign_organism_quick", "Assign genus / species", "context_assign_organism_quick",
               group=2, mutates=True, submenu="organisms"),
    ActionSpec("assign_organism", "Assign genus / species…", "context_assign_organism",
               group=2, mutates=True),
    ActionSpec("assign_scheme", "Assign typing scheme…", "context_assign_scheme",
               group=2, mutates=True),
    ActionSpec("rename_sample", "Rename isolate…", "context_rename_sample",
               needs="single", group=2, mutates=True),
    ActionSpec("rename_folder", "Rename this organism folder…", "context_rename_folder",
               needs="folder", group=2, mutates=True),
    ActionSpec("edit_annotations", "Edit annotations…", "context_edit_annotations",
               group=2, mutates=True),
    ActionSpec("add_collection", "Add {n} isolates to a collection…", "context_add_collection",
               group=2, mutates=True, singular="Add to a collection…"),
    # 3 — group them for review (a grouping is not a transmission inference)
    ActionSpec("highlight", "Highlight {n} isolates as a group…", "context_highlight",
               group=3, mutates=True, singular="Highlight as a group…"),
    ActionSpec("unhighlight", "Remove highlight from {n} isolates", "context_unhighlight",
               group=3, mutates=True, singular="Remove highlight"),
    ActionSpec("select_group", "Select this group's isolates", "context_select_group",
               needs="group", group=3),
    # 4 — bring more data in
    ActionSpec("add_files", "Add sequence files here…", "context_add_files",
               needs="any", group=4, mutates=True),
    ActionSpec("add_folder", "Add a sequence folder here…", "context_add_folder",
               needs="any", group=4, mutates=True),
    ActionSpec("import_scheme", "Import a scheme folder…", "context_import_scheme",
               needs="any", group=4, mutates=True),
    # 5 — take something away with you
    ActionSpec("export_selection", "Export {n} isolates…", "context_export_selection",
               group=5, singular="Export…", submenu="exports"),
    ActionSpec("copy_id", "Copy {n} sample IDs", "context_copy_id",
               group=5, singular="Copy sample ID"),
    ActionSpec("copy_path", "Copy file path", "context_copy_path", needs="single", group=5),
    ActionSpec("copy_scheme_path", "Copy scheme path", "context_copy_path",
               needs="paths", group=5),
    # 6 — destructive, last, never the default
    ActionSpec("archive", "Archive {n} isolates (keeps all evidence)…", "context_archive",
               group=6, destructive=True, mutates=True,
               singular="Archive this isolate (keeps all evidence)…"),
    ActionSpec("restore", "Restore {n} isolates from the archive", "context_restore",
               group=6, mutates=True, singular="Restore from the archive"),
    ActionSpec("remove", "Remove {n} isolates from this project…", "context_remove",
               group=6, destructive=True, mutates=True,
               singular="Remove this isolate from this project…"),
    ActionSpec("remove_scheme", "Remove this imported scheme…", "context_remove_scheme",
               needs="paths", group=6, destructive=True, mutates=True,
               enabled_when=_scheme_removable),
)}

_OPEN = ("open_record", "open_folder")
_ROUTING = ("select_in_library", "add_to_comparison", "add_to_characterization", "add_to_report")
_ASSIGN = ("assign_organism_quick", "assign_organism", "assign_scheme", "rename_sample",
           "edit_annotations", "add_collection")
_GROUPING = ("highlight", "unhighlight")
_ADDING = ("add_files", "add_folder")
_TAKING = ("export_selection", "copy_id", "copy_path")
_ARCHIVING = ("archive", "restore")

# The action set a plain sample table offers. Views deviate from it below.
SAMPLE_ACTIONS = _OPEN + _ROUTING + _ASSIGN + _GROUPING + _ADDING + _TAKING + _ARCHIVING

VIEW_ACTIONS: dict[str, tuple[str, ...]] = {
    "library": SAMPLE_ACTIONS + ("remove",),
    "overview.recent": SAMPLE_ACTIONS + ("remove",),
    "library.tree": (_OPEN + _ROUTING + ("assign_organism_quick", "assign_organism",
                     "rename_folder", "add_collection") + _GROUPING + _ADDING
                     + ("copy_id",) + _ARCHIVING),
    "schemes": ("open_scheme_folder", "copy_scheme_path", "import_scheme", "remove_scheme"),
    "compare.cohort": SAMPLE_ACTIONS + ("remove_from_comparison",),
    "compare.groups": ("select_group", "add_to_comparison", "add_to_report", "highlight",
                       "unhighlight", "copy_id"),
    "compare.profiles": SAMPLE_ACTIONS,
    "compare.graph": (("open_record",) + _ROUTING
                      + ("assign_organism_quick", "assign_organism")
                      + _GROUPING + ("select_group", "export_selection", "copy_id", "archive")),
    "characterization": SAMPLE_ACTIONS,
    "evidence.features": SAMPLE_ACTIONS,
    "evidence.amr": SAMPLE_ACTIONS,
    # Imported source rows are somebody else's evidence; this view never edits them.
    "evidence.hydra": ("open_record", "copy_id"),
    "reports": SAMPLE_ACTIONS + ("remove_from_report",),
    # A picker chooses; it never changes the project destructively.
    "picker.table": ("open_record", "assign_organism", "add_files", "copy_id"),
    "picker.folders": ("assign_organism", "add_files"),
}


def actions_for_view(view_id) -> tuple[ActionSpec, ...]:
    """Every action a view can ever offer, before the selection is considered."""
    keys = VIEW_ACTIONS.get(str(view_id))
    if keys is None:
        keys = tuple(key for key, spec in ACTIONS.items()
                     if spec.views is None or str(view_id) in spec.views)
    return tuple(ACTIONS[key] for key in dict.fromkeys(keys) if key in ACTIONS)


def is_applicable(spec, selection) -> bool:
    """Whether the entry belongs in the menu at all for this selection."""
    if selection.on_blank:
        return spec.needs == "any"
    if spec.needs == "any":
        return True
    if spec.needs == "folder":
        return selection.folder is not None
    if spec.needs == "paths":
        return bool(selection.paths)
    if spec.needs == "group":
        return bool(selection.group_id)
    # "samples" and "single" both need isolates. A single-isolate action stays
    # visible but disabled for a wider selection, so the scope is never narrowed
    # behind the user's back.
    return bool(selection.sample_ids)


def sample_state_reason(spec, samples) -> str:
    """Why the isolates themselves block this action, in the standard wording."""
    samples = list(samples or ())
    if spec.key in {"open_folder", "copy_path"} and any(not s.get("input_path") for s in samples):
        return REASONS["no_file"]
    if spec.mutates and any(s.get("status") == "running" for s in samples):
        return REASONS["running"]
    return ""


def action_state(spec, selection, *, window=None, busy=False, samples=None) -> tuple[bool, str]:
    """(enabled, reason). A blocked action keeps its place and explains itself."""
    if busy and spec.mutates:
        return False, REASONS["busy"]
    if spec.needs == "single" and selection.count != 1:
        return False, REASONS["single"]
    reason = sample_state_reason(spec, samples) if samples is not None else ""
    if reason:
        return False, reason
    if callable(spec.enabled_when):
        return spec.enabled_when(window, selection)
    return True, ""


def plan_for(view_id, selection) -> list:
    """The ordered menu for one view and one selection, with separators between groups.

    Pure: creates no widget and needs no QApplication.
    """
    applicable = [spec for spec in actions_for_view(view_id) if is_applicable(spec, selection)]
    applicable.sort(key=lambda spec: (spec.group, tuple(ACTIONS).index(spec.key)))
    plan: list = []
    previous = None
    for spec in applicable:
        if previous is not None and spec.group != previous:
            plan.append(SEPARATOR)
        plan.append(spec)
        previous = spec.group
    return plan


def handler_names() -> tuple[str, ...]:
    """Every method name the window must provide for the registry to be complete."""
    return tuple(sorted({spec.handler for spec in ACTIONS.values()}))


def missing_handlers(window) -> tuple[str, ...]:
    """Handler names the given window does not implement yet."""
    return tuple(name for name in handler_names() if not callable(getattr(window, name, None)))


def submenu_entries(name, window=None, selection=None) -> tuple[tuple[str, object], ...]:
    """(title, value) pairs for a submenu, built when the menu opens."""
    if name == "exports":
        return EXPORT_FORMATS
    if name == "organisms":
        return _organism_entries(window)
    return ()


def _organism_entries(window) -> tuple[tuple[str, object], ...]:
    """Genus/species pairs already used in this project, so one click covers the usual case."""
    samples = getattr(window, "current_samples", None) or ()
    from wmlstudio.ui_common import organism_for
    seen = {}
    for sample in samples:
        genus, species, _ = organism_for(sample)
        if not genus:
            continue
        seen.setdefault((genus, species or ""), f"{genus} {species}".strip())
    entries = [(title, pair) for pair, title in sorted(seen.items(), key=lambda item: item[1])]
    return tuple(entries[:12])


# --------------------------------------------------------------------------
# Adapters: a click point in one view becomes a Selection.
# --------------------------------------------------------------------------


def _dedup(values) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values if value))


class SelectionAdapter:
    """Base adapter. Subclasses resolve a point; everything else is shared."""

    def __init__(self, view_id, widget):
        self.view_id = str(view_id)
        self._widget = widget

    def widget(self) -> QWidget:
        return self._widget

    def viewport(self) -> QWidget:
        getter = getattr(self._widget, "viewport", None)
        return getter() if callable(getter) else self._widget

    def global_pos(self, point) -> QPoint:
        return self.viewport().mapToGlobal(point if point is not None else QPoint(0, 0))

    def keyboard_request(self, point) -> bool:
        """Shift+F10 and the Windows menu key arrive at the viewport origin, not on a row.

        Treat that as "act on what is already selected" instead of resolving an
        empty selection and offering a menu that cannot do anything.
        """
        if point is None:
            return True
        return point.isNull() or not self.viewport().rect().contains(point)

    def ensure_clicked_selected(self, point):
        """Right-click outside the selection replaces it, visibly, before the menu opens."""
        raise NotImplementedError

    def resolve(self, point) -> Selection:
        raise NotImplementedError


class TableWidgetAdapter(SelectionAdapter):
    """QTableWidget rows. Ids come from the id column's UserRole, never the row number."""

    def __init__(self, view_id, table, *, id_column=0):
        super().__init__(view_id, table)
        self.id_column = int(id_column)

    def id_at(self, row):
        item = self._widget.item(row, self.id_column)
        value = item.data(USER_ROLE) if item is not None else None
        return str(value) if value else None

    def ensure_clicked_selected(self, point):
        table = self._widget
        item = table.itemAt(point)
        if item is None:
            return None
        if item.row() not in {selected.row() for selected in table.selectedItems()}:
            table.selectRow(item.row())
        return item

    def resolve(self, point) -> Selection:
        table = self._widget
        keyboard = self.keyboard_request(point)
        item = None if keyboard else self.ensure_clicked_selected(point)
        rows = sorted({selected.row() for selected in table.selectedItems()})
        ids = _dedup(self.id_at(row) for row in rows)
        clicked = self.id_at(item.row()) if item is not None else None
        blank = item is None and not (keyboard and ids)
        return Selection(self.view_id, () if blank else ids, clicked,
                         on_blank=blank, from_keyboard=keyboard)


class PathTableAdapter(SelectionAdapter):
    """Rows that are files, not isolates: the scheme library reads its path column."""

    def __init__(self, view_id, table, *, path_column=2):
        super().__init__(view_id, table)
        self.path_column = int(path_column)

    def path_at(self, row):
        item = self._widget.item(row, self.path_column)
        text = item.text().strip() if item is not None else ""
        return Path(text) if text and text != "—" else None

    def ensure_clicked_selected(self, point):
        table = self._widget
        item = table.itemAt(point)
        if item is not None and item.row() not in {s.row() for s in table.selectedItems()}:
            table.selectRow(item.row())
        return item

    def resolve(self, point) -> Selection:
        table = self._widget
        keyboard = self.keyboard_request(point)
        item = None if keyboard else self.ensure_clicked_selected(point)
        rows = sorted({selected.row() for selected in table.selectedItems()})
        paths = tuple(dict.fromkeys(path for path in (self.path_at(row) for row in rows) if path))
        blank = item is None and not (keyboard and paths)
        return Selection(self.view_id, paths=() if blank else paths,
                         on_blank=blank, from_keyboard=keyboard)


class MatrixViewAdapter(SelectionAdapter):
    """QTableView over EvidenceMatrixModel: rows[i]['_sample_id'] follows the model's sort."""

    def __init__(self, view_id, view, *, id_field="_sample_id"):
        super().__init__(view_id, view)
        self.id_field = str(id_field)

    def _rows(self):
        model = self._widget.model()
        return getattr(model, "rows", None) or []

    def id_at(self, row):
        rows = self._rows()
        if not 0 <= row < len(rows):
            return None
        value = rows[row].get(self.id_field)
        return str(value) if value else None

    def _selected_rows(self):
        selection_model = self._widget.selectionModel()
        if selection_model is None:
            return []
        return sorted({index.row() for index in selection_model.selectedIndexes()})

    def ensure_clicked_selected(self, point):
        index = self._widget.indexAt(point)
        if not index.isValid():
            return QModelIndex()
        if index.row() not in self._selected_rows():
            self._widget.selectRow(index.row())
        return index

    def resolve(self, point) -> Selection:
        keyboard = self.keyboard_request(point)
        index = QModelIndex() if keyboard else self.ensure_clicked_selected(point)
        ids = _dedup(self.id_at(row) for row in self._selected_rows())
        clicked = self.id_at(index.row()) if index.isValid() else None
        blank = not index.isValid() and not (keyboard and ids)
        return Selection(self.view_id, () if blank else ids, clicked,
                         on_blank=blank, from_keyboard=keyboard)


class TreeWidgetAdapter(SelectionAdapter):
    """Organism folder trees. Decodes the tagged node data and expands a folder to its members."""

    def __init__(self, view_id, tree, *, dialect="library", folder_ids=None):
        super().__init__(view_id, tree)
        self.dialect = str(dialect)
        self.folder_ids = folder_ids

    def ensure_clicked_selected(self, point):
        tree = self._widget
        item = tree.itemAt(point)
        if item is not None and item not in tree.selectedItems():
            tree.setCurrentItem(item)
        return item

    @staticmethod
    def _descendant_ids(item) -> tuple[str, ...]:
        found: list[str] = []
        stack = [item]
        while stack:
            node = stack.pop()
            value = node.data(0, USER_ROLE)
            if isinstance(value, (tuple, list)) and len(value) == 2 and value[0] == "ids":
                found.extend(str(sample_id) for sample_id in value[1] or ())
            stack.extend(node.child(index) for index in range(node.childCount()))
        return _dedup(found)

    def _folder_members(self, item, genus, species) -> tuple[str, ...]:
        if callable(self.folder_ids):
            return _dedup(self.folder_ids(genus, species))
        return self._descendant_ids(item)

    def resolve(self, point) -> Selection:
        keyboard = self.keyboard_request(point)
        tree = self._widget
        item = (tree.currentItem() if keyboard else self.ensure_clicked_selected(point))
        if item is None:
            return Selection(self.view_id, on_blank=True, from_keyboard=keyboard)
        value = item.data(0, USER_ROLE)
        if self.dialect == "folders":
            return self._resolve_folder_node(item, value, keyboard)
        return self._resolve_library_node(item, value, keyboard)

    def _resolve_library_node(self, item, value, keyboard) -> Selection:
        if isinstance(value, (tuple, list)) and value:
            tag = value[0]
            if tag == "project":
                return Selection(self.view_id, paths=(Path(str(value[1])),),
                                 from_keyboard=keyboard)
            if tag == "ids":
                ids = _dedup(value[1] or ())
                return Selection(self.view_id, ids, ids[0] if len(ids) == 1 else None,
                                 from_keyboard=keyboard)
            if tag == "genus":
                genus = str(value[1])
                return Selection(self.view_id, self._folder_members(item, genus, ""),
                                 folder=(genus, ""), from_keyboard=keyboard)
            if tag == "organism":
                genus, species = str(value[1]), str(value[2]) if len(value) > 2 else ""
                return Selection(self.view_id, self._folder_members(item, genus, species),
                                 folder=(genus, species), from_keyboard=keyboard)
        # "All samples" and other untagged branches carry no scope of their own.
        return Selection(self.view_id, self._descendant_ids(item), on_blank=not item.childCount(),
                         from_keyboard=keyboard)

    def _resolve_folder_node(self, item, value, keyboard) -> Selection:
        if isinstance(value, (tuple, list)) and value:
            genus = str(value[0])
            species = str(value[1]) if len(value) > 1 else ""
            return Selection(self.view_id, self._folder_members(item, genus, species),
                             folder=(genus, species), from_keyboard=keyboard)
        return Selection(self.view_id, on_blank=True, from_keyboard=keyboard)


class GraphAdapter(SelectionAdapter):
    """The minimum spanning forest. A node that was not selected becomes the selection."""

    def node_key(self, item):
        key = getattr(item, "node_key", None)      # a GraphLabel points at its node
        if key is None:
            key = getattr(item, "key", None)       # a TreeNode carries the key itself
        return str(key) if key is not None and key in getattr(self._widget, "_members", {}) else None

    def ensure_clicked_selected(self, point):
        view = self._widget
        key = self.node_key(view.itemAt(point)) if point is not None else None
        if key is None:
            return None
        members = list(getattr(view, "_members", {}).get(key, ()))
        if members and not set(members) <= set(view.selected_ids()):
            view.select_ids(members)
        return key

    def resolve(self, point) -> Selection:
        view = self._widget
        keyboard = self.keyboard_request(point)
        key = None if keyboard else self.ensure_clicked_selected(point)
        ids = _dedup(view.selected_ids())
        members = _dedup(getattr(view, "_members", {}).get(key, ())) if key else ()
        blank = key is None and not ids
        return Selection(self.view_id, ids, members[0] if len(members) == 1 else None,
                         group_id=key, on_blank=blank, from_keyboard=keyboard)


def default_adapter(view_id, widget, **kwargs) -> SelectionAdapter:
    """Pick the adapter that matches the widget, so wiring a view is one line."""
    view_id = str(view_id)
    if isinstance(widget, QGraphicsView):
        return GraphAdapter(view_id, widget)
    if isinstance(widget, QTreeWidget):
        dialect = kwargs.pop("dialect", "folders" if view_id.startswith("picker.") else "library")
        return TreeWidgetAdapter(view_id, widget, dialect=dialect, **kwargs)
    if isinstance(widget, QTableWidget):
        if view_id == "schemes" or "path_column" in kwargs:
            return PathTableAdapter(view_id, widget, **kwargs)
        return TableWidgetAdapter(view_id, widget, **kwargs)
    if isinstance(widget, QTableView):
        return MatrixViewAdapter(view_id, widget, **kwargs)
    raise TypeError(f"No context-menu adapter for {type(widget).__name__} ({view_id}).")


def install_context_menu(widget, adapter, callback, **kwargs):
    """Route this view's right-click, menu key and Shift+F10 through one callback.

    `adapter` may be a ready adapter or a view id, in which case the matching one
    is built here. `callback(selection, global_position)` builds and shows the
    menu, which keeps the declarative half of this module free of QMenu.
    """
    if isinstance(adapter, str):
        adapter = default_adapter(adapter, widget, **kwargs)
    widget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    widget.customContextMenuRequested.connect(
        lambda point, a=adapter: callback(a.resolve(point), a.global_pos(point)))
    return adapter
