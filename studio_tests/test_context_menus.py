"""The reusable right-click system: what a click resolves to, and what it may do."""

from pathlib import Path

import pytest
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPoint, Qt
from PySide6.QtWidgets import (
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
)

from wmlstudio.context_menus import (
    ACTIONS,
    EXPORT_FORMATS,
    REASONS,
    SAMPLE_ACTIONS,
    SEPARATOR,
    VIEW_ACTIONS,
    GraphAdapter,
    MatrixViewAdapter,
    PathTableAdapter,
    Selection,
    TableWidgetAdapter,
    TreeWidgetAdapter,
    action_state,
    actions_for_view,
    default_adapter,
    handler_names,
    install_context_menu,
    missing_handlers,
    plan_for,
    sample_state_reason,
    submenu_entries,
)
from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import TreeView


def keys(plan):
    return [entry.key for entry in plan if entry is not SEPARATOR]


@pytest.fixture
def sample_table(qtbot):
    """A real sortable table built the way every table in the application is built."""
    table = make_table(["Sample", "Status"])
    qtbot.addWidget(table)
    table.resize(420, 320)
    table.setSortingEnabled(False)
    table.setRowCount(3)
    for row, (name, sample_id) in enumerate([("Alpha", "id-a"), ("Bravo", "id-b"),
                                             ("Charlie", "id-c")]):
        table.setItem(row, 0, cell(name, sample_id))
        table.setItem(row, 1, cell("Complete", sample_id))
    table.setSortingEnabled(True)
    table.show()
    qtbot.waitUntil(table.isVisible, timeout=2000)
    return table


def point_on_row(table, row):
    return table.visualItemRect(table.item(row, 0)).center()


def row_of(table, name):
    """make_table turns sorting on, so a name's visual row is decided by the view."""
    matches = table.findItems(name, Qt.MatchFlag.MatchExactly)
    assert matches, name
    return matches[0].row()


# --- the declarative half ------------------------------------------------


def test_the_library_offers_archive_and_remove_while_a_picker_offers_neither():
    selection = Selection("library", ("id-a", "id-b"))
    assert "archive" in keys(plan_for("library", selection))
    assert "remove" in keys(plan_for("library", selection))
    picker = keys(plan_for("picker.table", Selection("picker.table", ("id-a", "id-b"))))
    assert "archive" not in picker and "remove" not in picker
    folders = keys(plan_for("picker.folders", Selection("picker.folders", ("id-a",),
                                                        folder=("Klebsiella", "pneumoniae"))))
    assert "archive" not in folders and "remove" not in folders


def test_imported_source_rows_are_read_only():
    plan = keys(plan_for("evidence.hydra", Selection("evidence.hydra", ("id-a",))))
    assert plan == ["open_record", "copy_id"]


def test_a_right_click_on_empty_space_still_offers_something_to_do():
    for view_id in VIEW_ACTIONS:
        plan = plan_for(view_id, Selection(view_id, on_blank=True))
        assert all(entry is SEPARATOR or entry.needs == "any" for entry in plan)
    plan = plan_for("library", Selection("library", on_blank=True))
    assert keys(plan) == ["add_files", "add_folder"]


def test_sections_are_separated_and_destructive_actions_come_last():
    plan = plan_for("library", Selection("library", ("id-a",)))
    assert SEPARATOR in plan
    ordered = keys(plan)
    assert ordered[0] == "open_record"
    assert ordered[-1] == "remove"
    assert ordered.index("archive") > ordered.index("copy_id")
    assert plan[plan.index(ACTIONS["archive"]) - 1] is SEPARATOR


def test_titles_always_show_how_many_isolates_an_action_would_touch():
    many = Selection("library", ("a", "b", "c", "d", "e", "f", "g"))
    one = Selection("library", ("a",))
    assert ACTIONS["archive"].format_title(many) == "Archive 7 isolates (keeps all evidence)…"
    assert ACTIONS["archive"].format_title(one) == "Archive this isolate (keeps all evidence)…"
    assert ACTIONS["copy_id"].format_title(many) == "Copy 7 sample IDs"
    assert ACTIONS["copy_id"].format_title(one) == "Copy sample ID"
    assert ACTIONS["open_record"].format_title(many) == "Open isolate record…"


def test_a_single_isolate_action_is_disabled_with_a_reason_never_silently_narrowed():
    many = Selection("library", ("a", "b"))
    assert "rename_sample" in keys(plan_for("library", many)), "it must stay visible"
    enabled, reason = action_state(ACTIONS["rename_sample"], many)
    assert enabled is False
    assert reason == REASONS["single"]
    assert action_state(ACTIONS["rename_sample"], Selection("library", ("a",))) == (True, "")


def test_a_running_task_disables_the_actions_that_would_write_and_nothing_else():
    selection = Selection("library", ("a", "b"))
    blocked, reason = action_state(ACTIONS["archive"], selection, busy=True)
    assert blocked is False and reason == REASONS["busy"]
    assert action_state(ACTIONS["copy_id"], selection, busy=True) == (True, "")
    assert action_state(ACTIONS["export_selection"], selection, busy=True) == (True, "")
    assert action_state(ACTIONS["open_record"], Selection("library", ("a",)), busy=True) == (True, "")


def test_an_isolate_without_a_file_or_mid_analysis_blocks_the_action_with_a_reason():
    one = Selection("library", ("a",))
    running = [{"id": "a", "status": "running", "input_path": "/tmp/a.fasta"}]
    assert action_state(ACTIONS["archive"], one, samples=running) == (False, REASONS["running"])
    assert action_state(ACTIONS["copy_id"], one, samples=running) == (True, "")

    profile_only = [{"id": "a", "status": "completed", "input_path": ""}]
    blocked, reason = action_state(ACTIONS["open_folder"], one, samples=profile_only)
    assert blocked is False and reason == REASONS["no_file"]
    assert action_state(ACTIONS["open_record"], one, samples=profile_only) == (True, "")
    assert sample_state_reason(ACTIONS["archive"], []) == ""


def test_a_bundled_reference_snapshot_cannot_be_removed(tmp_path):
    class Window:
        root = tmp_path

    imported = Selection("schemes", paths=(tmp_path / "schemes" / "senterica",))
    bundled = Selection("schemes", paths=(Path("/opt/app/resources/schemes/senterica"),))
    assert action_state(ACTIONS["remove_scheme"], imported, window=Window()) == (True, "")
    enabled, reason = action_state(ACTIONS["remove_scheme"], bundled, window=Window())
    assert enabled is False
    assert reason == REASONS["bundled"]


def test_scheme_rows_are_files_and_offer_no_isolate_actions():
    plan = keys(plan_for("schemes", Selection("schemes", paths=(Path("/tmp/scheme"),))))
    assert plan == ["open_scheme_folder", "import_scheme", "copy_scheme_path", "remove_scheme"]
    assert "archive" not in plan


def test_group_rows_only_offer_actions_a_group_can_answer():
    selection = Selection("compare.groups", ("a", "b"), group_id="cluster-1")
    plan = keys(plan_for("compare.groups", selection))
    assert "select_group" in plan
    assert "remove" not in plan and "rename_sample" not in plan
    without_group = keys(plan_for("compare.groups", Selection("compare.groups", ("a",))))
    assert "select_group" not in without_group


def test_every_action_names_a_handler_the_window_must_implement():
    """Handlers live on the window; this module only declares their names."""
    assert len(handler_names()) >= 25
    for name in handler_names():
        assert name.startswith("context_") and name.isidentifier()

    class HalfDone:
        def context_copy_id(self, selection):
            return None

    missing = missing_handlers(HalfDone())
    assert "context_copy_id" not in missing
    assert "context_archive" in missing


def test_submenus_are_built_at_menu_time_from_the_project_itself():
    assert submenu_entries("exports") == EXPORT_FORMATS

    class Window:
        current_samples = [
            {"id": "a", "metadata": {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}}},
            {"id": "b", "metadata": {"organism": {"genus": "Klebsiella", "species": "pneumoniae"}}},
            {"id": "c", "metadata": {"organism": {"genus": "Escherichia", "species": "coli"}}},
            {"id": "d", "metadata": {}},
        ]

    entries = submenu_entries("organisms", Window())
    assert entries == (("Escherichia coli", ("Escherichia", "coli")),
                       ("Klebsiella pneumoniae", ("Klebsiella", "pneumoniae")))
    assert submenu_entries("nothing-like-this") == ()


def test_an_unknown_view_falls_back_to_the_standard_sample_action_set():
    fallback = [spec.key for spec in actions_for_view("something.new")]
    assert set(SAMPLE_ACTIONS) <= set(fallback)


# --- adapters ------------------------------------------------------------


def test_ids_come_from_the_item_data_even_after_the_user_sorts_the_table(sample_table):
    """Sorting is on in every table, so the visual row is not the data row."""
    adapter = TableWidgetAdapter("library", sample_table)
    sample_table.sortItems(0, Qt.SortOrder.AscendingOrder)
    assert sample_table.item(0, 0).text() == "Alpha"
    assert adapter.resolve(point_on_row(sample_table, 0)).sample_ids == ("id-a",)

    sample_table.sortItems(0, Qt.SortOrder.DescendingOrder)
    assert sample_table.item(0, 0).text() == "Charlie"
    selection = adapter.resolve(point_on_row(sample_table, 0))
    assert selection.sample_ids == ("id-c",), "row 0 is not data row 0 once sorted"
    assert selection.clicked_id == "id-c"
    assert selection.single == "id-c"
    assert selection.on_blank is False


def test_right_clicking_outside_the_selection_replaces_it_before_the_menu_opens(sample_table):
    adapter = TableWidgetAdapter("library", sample_table)
    alpha, charlie = row_of(sample_table, "Alpha"), row_of(sample_table, "Charlie")
    sample_table.selectRow(alpha)
    selection = adapter.resolve(point_on_row(sample_table, charlie))
    assert selection.sample_ids == ("id-c",)
    assert {item.row() for item in sample_table.selectedItems()} == {charlie}


def test_right_clicking_inside_the_selection_keeps_the_whole_selection(sample_table):
    adapter = TableWidgetAdapter("library", sample_table)
    alpha, charlie = row_of(sample_table, "Alpha"), row_of(sample_table, "Charlie")
    sample_table.selectRow(alpha)
    for column in range(sample_table.columnCount()):
        sample_table.item(charlie, column).setSelected(True)
    selection = adapter.resolve(point_on_row(sample_table, charlie))
    assert set(selection.sample_ids) == {"id-a", "id-c"}
    assert selection.count == 2
    assert selection.single is None
    assert {item.row() for item in sample_table.selectedItems()} == {alpha, charlie}


def test_a_click_on_empty_space_below_the_rows_acts_on_nothing(sample_table):
    adapter = TableWidgetAdapter("library", sample_table)
    sample_table.selectRow(1)
    blank = QPoint(5, sample_table.viewport().rect().bottom() - 2)
    assert sample_table.itemAt(blank) is None
    selection = adapter.resolve(blank)
    assert selection.on_blank is True
    assert selection.sample_ids == ()
    assert keys(plan_for("library", selection)) == ["add_files", "add_folder"]


def test_the_menu_key_falls_back_to_the_current_selection(sample_table):
    """Shift+F10 and the Windows menu key arrive at the viewport origin, not on a row."""
    adapter = TableWidgetAdapter("library", sample_table)
    sample_table.selectRow(1)
    selection = adapter.resolve(QPoint(0, 0))
    assert selection.from_keyboard is True
    assert selection.on_blank is False
    assert selection.sample_ids == ("id-b",)
    assert selection.clicked_id is None
    assert "archive" in keys(plan_for("library", selection))


def test_the_menu_key_with_nothing_selected_still_offers_the_blank_menu(sample_table):
    adapter = TableWidgetAdapter("library", sample_table)
    sample_table.clearSelection()
    selection = adapter.resolve(QPoint(0, 0))
    assert selection.from_keyboard is True
    assert selection.on_blank is True
    assert keys(plan_for("library", selection)) == ["add_files", "add_folder"]


def test_a_scheme_row_resolves_to_its_path_and_not_to_an_isolate(qtbot, tmp_path):
    table = QTableWidget(2, 3)
    qtbot.addWidget(table)
    table.resize(520, 220)
    for row, name in enumerate(["senterica", "kpneumoniae"]):
        table.setItem(row, 0, QTableWidgetItem(name))
        table.setItem(row, 1, QTableWidgetItem("Imported"))
        table.setItem(row, 2, QTableWidgetItem(str(tmp_path / "schemes" / name)))
    table.show()
    qtbot.waitUntil(table.isVisible, timeout=2000)
    adapter = PathTableAdapter("schemes", table)
    selection = adapter.resolve(table.visualItemRect(table.item(1, 0)).center())
    assert selection.paths == (tmp_path / "schemes" / "kpneumoniae",)
    assert selection.sample_ids == ()
    assert keys(plan_for("schemes", selection))[0] == "open_scheme_folder"


def test_a_matrix_view_resolves_through_the_models_own_row_order(qtbot):
    class Matrix(QAbstractTableModel):
        """Mirrors EvidenceMatrixModel: rows are dicts and sort() reorders them in place."""

        def __init__(self):
            super().__init__()
            self.headers = ["Isolate"]
            self.rows = [{"Isolate": "Alpha", "_sample_id": "id-a"},
                         {"Isolate": "Bravo", "_sample_id": "id-b"}]

        def rowCount(self, parent=QModelIndex()):
            return 0 if parent.isValid() else len(self.rows)

        def columnCount(self, parent=QModelIndex()):
            return 0 if parent.isValid() else len(self.headers)

        def data(self, index, role=Qt.ItemDataRole.DisplayRole):
            if index.isValid() and role == Qt.ItemDataRole.DisplayRole:
                return self.rows[index.row()][self.headers[index.column()]]
            return None

        def sort(self, column, order=Qt.SortOrder.AscendingOrder):
            self.layoutAboutToBeChanged.emit()
            self.rows.sort(key=lambda row: row["Isolate"],
                           reverse=order == Qt.SortOrder.DescendingOrder)
            self.layoutChanged.emit()

    view = QTableView()
    qtbot.addWidget(view)
    model = Matrix()
    view.setModel(model)
    view.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
    view.resize(360, 200)
    view.show()
    qtbot.waitUntil(view.isVisible, timeout=2000)
    model.sort(0, Qt.SortOrder.DescendingOrder)
    adapter = MatrixViewAdapter("evidence.features", view)
    point = view.visualRect(model.index(0, 0)).center()
    selection = adapter.resolve(point)
    assert selection.sample_ids == ("id-b",)
    assert selection.clicked_id == "id-b"


def library_tree(qtbot):
    tree = QTreeWidget()
    qtbot.addWidget(tree)
    tree.setHeaderHidden(True)
    tree.resize(320, 420)
    root = QTreeWidgetItem(["All samples  (3)"])
    root.setData(0, Qt.ItemDataRole.UserRole, None)
    tree.addTopLevelItem(root)
    genus = QTreeWidgetItem(["Klebsiella"])
    genus.setData(0, Qt.ItemDataRole.UserRole, ("genus", "Klebsiella"))
    root.addChild(genus)
    species = QTreeWidgetItem(["pneumoniae"])
    species.setData(0, Qt.ItemDataRole.UserRole, ("organism", "Klebsiella", "pneumoniae"))
    genus.addChild(species)
    st = QTreeWidgetItem(["ST 258  (2)"])
    st.setData(0, Qt.ItemDataRole.UserRole, ("ids", ["id-a", "id-b"]))
    species.addChild(st)
    other = QTreeWidgetItem(["Escherichia"])
    other.setData(0, Qt.ItemDataRole.UserRole, ("genus", "Escherichia"))
    root.addChild(other)
    coli = QTreeWidgetItem(["coli"])
    coli.setData(0, Qt.ItemDataRole.UserRole, ("organism", "Escherichia", "coli"))
    other.addChild(coli)
    lone = QTreeWidgetItem(["ST 131  (1)"])
    lone.setData(0, Qt.ItemDataRole.UserRole, ("ids", ["id-c"]))
    coli.addChild(lone)
    projects = QTreeWidgetItem(["Other saved projects"])
    tree.addTopLevelItem(projects)
    stored = QTreeWidgetItem(["Ward study  (4)"])
    stored.setData(0, Qt.ItemDataRole.UserRole, ("project", "/data/Ward study.wmlstudio"))
    projects.addChild(stored)
    tree.expandAll()
    tree.show()
    qtbot.waitUntil(tree.isVisible, timeout=2000)
    return tree


def tree_point(tree, text):
    matches = tree.findItems(text, Qt.MatchFlag.MatchExactly | Qt.MatchFlag.MatchRecursive)
    assert matches, text
    return tree.visualItemRect(matches[0]).center()


def test_an_organism_folder_expands_to_exactly_its_own_isolates(qtbot):
    tree = library_tree(qtbot)
    adapter = TreeWidgetAdapter("library.tree", tree)
    selection = adapter.resolve(tree_point(tree, "pneumoniae"))
    assert selection.folder == ("Klebsiella", "pneumoniae")
    assert selection.sample_ids == ("id-a", "id-b")
    assert "rename_folder" in keys(plan_for("library.tree", selection))


def test_a_genus_folder_covers_its_species_and_a_leaf_covers_one_isolate(qtbot):
    tree = library_tree(qtbot)
    adapter = TreeWidgetAdapter("library.tree", tree)
    genus = adapter.resolve(tree_point(tree, "Escherichia"))
    assert genus.folder == ("Escherichia", "")
    assert genus.sample_ids == ("id-c",)
    leaf = adapter.resolve(tree_point(tree, "ST 258  (2)"))
    assert leaf.folder is None
    assert leaf.sample_ids == ("id-a", "id-b")


def test_a_saved_project_node_resolves_to_a_path_not_to_isolates(qtbot):
    tree = library_tree(qtbot)
    adapter = TreeWidgetAdapter("library.tree", tree)
    selection = adapter.resolve(tree_point(tree, "Ward study  (4)"))
    assert selection.paths == (Path("/data/Ward study.wmlstudio"),)
    assert selection.sample_ids == ()


def test_the_cohort_pickers_folder_tree_uses_its_own_plain_tuples(qtbot):
    tree = QTreeWidget()
    qtbot.addWidget(tree)
    tree.resize(300, 300)
    everything = QTreeWidgetItem(["All isolates (3)"])
    everything.setData(0, Qt.ItemDataRole.UserRole, None)
    tree.addTopLevelItem(everything)
    genus = QTreeWidgetItem(["Klebsiella (2)"])
    genus.setData(0, Qt.ItemDataRole.UserRole, ("Klebsiella",))
    everything.addChild(genus)
    species = QTreeWidgetItem(["pneumoniae (2)"])
    species.setData(0, Qt.ItemDataRole.UserRole, ("Klebsiella", "pneumoniae"))
    genus.addChild(species)
    tree.expandAll()
    tree.show()
    qtbot.waitUntil(tree.isVisible, timeout=2000)
    members = {("Klebsiella", "pneumoniae"): ["id-a", "id-b"]}
    adapter = TreeWidgetAdapter("picker.folders", tree, dialect="folders",
                                folder_ids=lambda genus, species: members.get((genus, species), []))
    selection = adapter.resolve(tree_point(tree, "pneumoniae (2)"))
    assert selection.folder == ("Klebsiella", "pneumoniae")
    assert selection.sample_ids == ("id-a", "id-b")
    assert "archive" not in keys(plan_for("picker.folders", selection))


def test_a_graph_node_that_was_not_selected_becomes_the_selection(qtbot):
    # The graph is laid out but never shown: an offscreen repaint of a live scene
    # after this many torn-down widgets crashes Qt, and the adapter only needs the
    # view's geometry, which resize() already provides.
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(700, 500)
    results = [{"sample_id": key, "sample_name": f"Isolate {key}", "scheme_digest": "v1",
                "st": "258", "alleles": {"locus1": "1"},
                "calls": [{"locus": "locus1", "allele": "1", "status": "exact"}]}
               for key in ("a", "b")]
    view.draw_results(results, [{"source": "a", "target": "b", "distance": 1, "shared_loci": 1}])
    view.select_ids(["a"])
    adapter = GraphAdapter("compare.graph", view)

    node = view.nodes["b"]
    point = view.mapFromScene(node.sceneBoundingRect().center())
    selection = adapter.resolve(point)
    assert selection.sample_ids == ("b",)
    assert selection.group_id == "b"
    assert view.selected_ids() == ["b"]
    assert "select_group" in keys(plan_for("compare.graph", selection))
    assert "remove" not in keys(plan_for("compare.graph", selection))


def test_the_graph_menu_key_uses_whatever_is_already_selected(qtbot):
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(700, 500)
    results = [{"sample_id": key, "sample_name": f"Isolate {key}", "scheme_digest": "v1",
                "st": "258", "alleles": {"locus1": "1"},
                "calls": [{"locus": "locus1", "allele": "1", "status": "exact"}]}
               for key in ("a", "b")]
    view.draw_results(results, [])
    view.select_ids(["a", "b"])
    adapter = GraphAdapter("compare.graph", view)
    selection = adapter.resolve(QPoint(0, 0))
    assert selection.from_keyboard is True
    assert selection.sample_ids == ("a", "b")
    assert selection.group_id is None


def test_the_right_adapter_is_chosen_for_each_kind_of_view(qtbot, sample_table):
    assert isinstance(default_adapter("library", sample_table), TableWidgetAdapter)
    assert isinstance(default_adapter("schemes", sample_table), PathTableAdapter)
    tree = QTreeWidget()
    qtbot.addWidget(tree)
    assert default_adapter("library.tree", tree).dialect == "library"
    assert default_adapter("picker.folders", tree).dialect == "folders"
    view = QTableView()
    qtbot.addWidget(view)
    assert isinstance(default_adapter("evidence.amr", view), MatrixViewAdapter)
    graph = TreeView()
    qtbot.addWidget(graph)
    assert isinstance(default_adapter("compare.graph", graph), GraphAdapter)
    with pytest.raises(TypeError):
        default_adapter("mystery", QTreeWidgetItem(["x"]))


def test_installing_a_menu_routes_every_request_through_one_callback(sample_table):
    seen = []
    adapter = install_context_menu(sample_table, "library",
                                   lambda selection, position: seen.append((selection, position)))
    assert isinstance(adapter, TableWidgetAdapter)
    assert sample_table.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu
    sample_table.selectRow(1)
    sample_table.customContextMenuRequested.emit(point_on_row(sample_table, 1))
    assert len(seen) == 1
    selection, position = seen[0]
    assert selection.view_id == "library"
    assert selection.sample_ids == ("id-b",)
    assert not position.isNull()
    assert adapter.widget() is sample_table


# ---------------------------------------------------------------------------
# The window half: one menu builder, the generic handlers, and a soft report of
# whatever the per-view handlers have not landed yet.
# ---------------------------------------------------------------------------


@pytest.fixture
def studio(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "workspace")
    qtbot.addWidget(widget)
    yield widget
    widget.close()


def imported(studio, tmp_path, names=("alpha", "bravo")):
    for name in names:
        path = tmp_path / f"{name}.fasta"
        path.write_text(">c\nACGTACGTACGT\n", encoding="utf-8")
        studio.import_paths([path])
    return [sample["id"] for sample in studio.project.samples()]


def test_the_window_reports_the_handlers_it_still_lacks_instead_of_failing_to_build(studio):
    """A partially wired window must construct; the gap is reportable, never fatal."""
    report = studio.context_menu_report()
    assert isinstance(report, tuple)
    assert set(report) <= set(handler_names())
    for name in ("context_copy_id", "context_open_folder", "context_export_selection",
                 "context_copy_path"):
        assert name not in report
        assert callable(getattr(studio, name))


def test_a_menu_only_offers_entries_this_window_can_actually_run(studio, tmp_path):
    ids = imported(studio, tmp_path)
    selection = Selection("library", tuple(ids))
    plan = studio.context_menu_plan(selection)
    keys = [entry.key for entry in plan if entry is not SEPARATOR]
    assert "copy_id" in keys and "export_selection" in keys
    for entry in plan:
        if entry is not SEPARATOR:
            assert callable(getattr(studio, entry.handler)), entry.key
    assert plan and plan[0] is not SEPARATOR and plan[-1] is not SEPARATOR


def test_the_scheme_and_imported_source_views_are_wired_to_the_one_callback(studio):
    from PySide6.QtCore import Qt as QtCore_Qt
    for view_id, widget in (("schemes", studio.scheme_table),
                            ("evidence.hydra", studio.hydra_table)):
        assert studio._context_adapters[view_id].widget() is widget
        assert widget.contextMenuPolicy() == QtCore_Qt.ContextMenuPolicy.CustomContextMenu


def test_copying_identifiers_copies_exactly_the_selected_rows(studio, tmp_path):
    from PySide6.QtWidgets import QApplication
    ids = imported(studio, tmp_path)
    studio.context_copy_id(Selection("library", tuple(ids)))
    assert QApplication.clipboard().text() == "\n".join(ids)
    studio.context_copy_id(Selection("library", ()))
    assert QApplication.clipboard().text() == "\n".join(ids)


def test_opening_the_containing_folder_opens_the_users_own_input_directory(studio, tmp_path,
                                                                          monkeypatch):
    opened = []
    monkeypatch.setattr("wmlstudio.context_menus.QDesktopServices.openUrl", opened.append)
    ids = imported(studio, tmp_path, names=("alpha",))
    studio.context_open_folder(Selection("library", (ids[0],), clicked_id=ids[0]))
    # toLocalFile gives forward slashes on Windows, so compare locations, not spelling.
    assert [Path(url.toLocalFile()) for url in opened] == [tmp_path]


def test_opening_a_folder_for_an_isolate_without_a_file_says_so_rather_than_guessing(studio):
    messages = []
    studio.notify = messages.append
    studio.current_samples = [{"id": "x", "input_path": None}]
    studio.context_open_folder(Selection("library", ("x",)))
    assert messages and "no attached sequence file" in messages[-1]


def test_copying_a_path_reads_the_scheme_row_and_the_isolate_row(studio, tmp_path):
    from PySide6.QtWidgets import QApplication
    ids = imported(studio, tmp_path, names=("alpha",))
    studio.context_copy_path(Selection("library", (ids[0],)))
    assert QApplication.clipboard().text() == str(tmp_path / "alpha.fasta")
    studio.context_copy_path(Selection("schemes", paths=(Path("/schemes/practice_7"),)))
    assert QApplication.clipboard().text() == str(Path("/schemes/practice_7"))


def test_exporting_a_selection_writes_only_the_isolates_that_were_selected(studio, tmp_path,
                                                                           monkeypatch):
    import csv
    ids = imported(studio, tmp_path, names=("alpha", "bravo", "charlie"))
    destination = tmp_path / "selection.csv"
    monkeypatch.setattr("wmlstudio.context_menus.QFileDialog.getSaveFileName",
                        lambda *args, **kwargs: (str(destination), ""))
    messages = []
    studio.notify = messages.append
    studio.context_export_selection(Selection("library", (ids[0], ids[2])), "csv")
    with destination.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    exported = {row["sample_name"] for row in rows}
    assert exported == {"alpha", "charlie"}
    assert "2 explicitly selected isolates" in messages[-1]


def test_exporting_nothing_refuses_instead_of_exporting_the_whole_project(studio, monkeypatch):
    asked = []
    monkeypatch.setattr("wmlstudio.context_menus.QFileDialog.getSaveFileName",
                        lambda *args, **kwargs: asked.append(True) or ("", ""))
    messages = []
    studio.notify = messages.append
    studio.context_export_selection(Selection("library", ()), "csv")
    assert asked == []
    assert "Nothing is included automatically" in messages[-1]


def test_running_an_action_the_window_does_not_implement_does_nothing(studio):
    assert studio.run_context_action(None, Selection("library", ("a",))) is None
    assert studio.run_context_action(("context_not_a_real_handler", None),
                                     Selection("library", ("a",))) is None
