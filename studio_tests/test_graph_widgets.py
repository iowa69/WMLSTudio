"""Native graph interaction, stable presentation state, exports, and dark surfaces."""

import copy
import json
import math
from xml.etree import ElementTree

import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QImage, QPalette
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QGraphicsView,
    QLineEdit,
    QMenu,
    QStackedWidget,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from wmlstudio.theme import BACKGROUND, INK, STYLE, apply_dark_palette
from wmlstudio.widgets import TreeView, forest_layout


def records():
    rows = []
    for key, allele, location in [("a", "1", "Ward A"), ("b", "1", "Ward B"),
                                  ("c", "2", "Ward A"), ("d", None, None)]:
        rows.append({"sample_id": key, "sample_name": f"Isolate {key}", "scheme_digest": "same-version",
                     "st": "184" if allele == "1" else None,
                     "alleles": {"locus1": "1", "locus2": allele},
                     "calls": [{"locus": "locus1", "allele": "1", "status": "exact"},
                               {"locus": "locus2", "allele": allele,
                                "status": "exact" if allele else "missing"}],
                     "metadata": {"location": location}})
    return rows


def edges():
    return [{"source": "a", "target": "b", "distance": 0, "shared_loci": 2},
            {"source": "b", "target": "c", "distance": 1, "shared_loci": 2}]


@pytest.fixture
def graph(qtbot):
    widget = TreeView()
    qtbot.addWidget(widget)
    widget.resize(900, 600)
    widget.show()
    widget.draw_results(records(), edges())
    return widget


def test_force_layout_is_deterministic_and_disconnected_nodes_are_retained():
    first = forest_layout(["a", "b", "c", "d"], edges())
    second = forest_layout(["d", "c", "a", "b"], list(reversed(edges())))
    assert first == second
    assert set(first) == {"a", "b", "c", "d"}
    assert len(set(first.values())) == 4
    assert all(math.isfinite(coordinate) for point in first.values() for coordinate in point)


@pytest.mark.parametrize("count", [2, 7])
@pytest.mark.parametrize("size", [(590, 216), (790, 440)])
def test_fitted_small_cohort_keeps_every_label_inside_compact_view(qtbot, count, size):
    widget = TreeView()
    qtbot.addWidget(widget)
    widget.resize(*size)
    widget.show()
    rows = [{"sample_id": str(index), "sample_name": f"isolate{index}", "st": None}
            for index in range(count)]
    links = [{"source": str(index - 1), "target": str(index), "distance": 1, "shared_loci": 7}
             for index in range(1, count)]
    widget.draw_results(rows, links)
    widget.fit_tree()
    for title in widget.labels.values():
        projected = widget.mapFromScene(title.sceneBoundingRect()).boundingRect()
        assert widget.viewport().rect().contains(projected)
    if count == 2:
        assert not any(guide.isVisible() for guide in widget._label_guides.values())
        assert widget.transform().m11() * 48 >= 22


def test_selection_multiselect_signal_and_keyboard(graph, qtbot):
    observed = []
    graph.selectionChanged.connect(observed.append)
    graph.select_ids(["a", "c"])
    assert graph.selected_ids() == ["a", "c"]
    assert observed[-1] == ["a", "c"]
    graph.setFocus()
    qtbot.keyClick(graph, Qt.Key.Key_A, modifier=Qt.KeyboardModifier.ControlModifier)
    assert graph.selected_ids() == ["a", "b", "c", "d"]
    qtbot.keyClick(graph, Qt.Key.Key_Escape)
    assert graph.selected_ids() == []
    assert graph.dragMode() == QGraphicsView.DragMode.RubberBandDrag
    graph.set_interaction_mode("pan")
    assert graph.dragMode() == QGraphicsView.DragMode.ScrollHandDrag
    with pytest.raises(ValueError):
        graph.set_interaction_mode("invalid")


def test_native_click_ctrl_click_and_double_click(graph, qtbot):
    first = graph.mapFromScene(graph.nodes["a"].pos())
    second = graph.mapFromScene(graph.nodes["b"].pos())
    qtbot.mouseClick(graph.viewport(), Qt.MouseButton.LeftButton, pos=first)
    qtbot.mouseClick(graph.viewport(), Qt.MouseButton.LeftButton,
                     Qt.KeyboardModifier.ControlModifier, pos=second)
    assert graph.selected_ids() == ["a", "b"]
    with qtbot.waitSignal(graph.nodeActivated) as activated:
        qtbot.mouseDClick(graph.viewport(), Qt.MouseButton.LeftButton, pos=first)
    assert activated.args == ["a"]


def test_drag_updates_edges_and_refresh_preserves_layout_colors_labels(graph):
    a, _b, line, _text = graph.edges[0]
    old = line.line().p1()
    graph.nodes[a].setPos(333, 222)
    assert line.line().p1() != old
    graph.select_ids(["a", "b"])
    graph.set_selected_color("#FF4400")
    graph.set_node_label("a", "A custom display name")
    graph.draw_results(list(reversed(records())), list(reversed(edges())), cluster_threshold=0)
    assert graph.nodes[a].pos().x() == 333
    assert graph.nodes[a].pos().y() == 222
    assert graph.nodes["a"].slices == [("#ff4400", 1)]
    assert "A custom display name" in graph.labels["a"].text()
    assert graph.selected_ids() == ["a", "b"]
    assert any(line.pen().style() == Qt.PenStyle.DashLine for _, _, line, _ in graph.edges)


def test_metadata_coloring_missing_values_highlight_and_reset(graph):
    assert graph.available_color_fields() == ["cluster", "st", "metadata:location"]
    graph.set_color_by("metadata:location")
    assert set(graph.legend()) == {"Ward A", "Ward B", "Not recorded"}
    assert graph.nodes["a"].slices == graph.nodes["c"].slices
    assert graph.nodes["a"].slices != graph.nodes["b"].slices
    assert graph.highlight("ward b") == ["b"]
    assert graph.nodes["b"].highlighted
    assert not graph.nodes["a"].highlighted
    graph.highlight("")
    assert not graph.nodes["b"].highlighted
    graph.select_ids(["a"])
    graph.set_selected_color("magenta")
    assert graph.nodes["a"].slices == [("#ff00ff", 1)]
    graph.reset_colors(["a"])
    assert graph.nodes["a"].slices == graph.nodes["c"].slices


def test_nested_project_annotations_are_individual_color_fields(graph):
    rows = records()
    rows[0]["metadata"] = {"annotations": {"ward": "ICU", "collection_date": "2026-09-12"}}
    rows[1]["metadata"] = {"annotations": {"ward": "ICU"}}
    graph.draw_results(rows, edges())
    assert "metadata:annotations.ward" in graph.available_color_fields()
    graph.set_color_by("metadata:annotations.ward")
    assert graph.nodes["a"].slices == graph.nodes["b"].slices
    assert set(graph.legend()) == {"ICU", "Not recorded"}


def test_merge_requires_complete_exact_genotype_and_renders_metadata_pies(graph):
    original = copy.deepcopy(records())
    graph.set_color_by("metadata:location")
    graph.set_merge_identical(True)
    assert set(graph.nodes) == {"a", "c", "d"}
    assert graph.nodes["a"].count == 2
    assert len(graph.nodes["a"].slices) == 2
    assert len(graph.edges) == 1
    graph.select_ids(["b"])
    assert graph.selected_ids() == ["a", "b"]
    assert records() == original
    # Both have zero differences over one shared locus, but incomplete profiles stay separate.
    missing = copy.deepcopy(original[-1])
    missing.update(sample_id="e", sample_name="Other incomplete")
    graph.draw_results(original + [missing], edges() + [
        {"source": "d", "target": "e", "distance": 0, "shared_loci": 1}])
    assert "d" in graph.nodes and "e" in graph.nodes


@pytest.mark.parametrize("change", ["digest", "calls_missing", "calls_mixed", "allele_missing"])
def test_merge_does_not_guess_complete_evidence(graph, change):
    rows = records()[:2]
    if change == "digest":
        rows[1]["scheme_digest"] = "different-version"
    elif change == "calls_missing":
        rows[1].pop("calls")
    elif change == "calls_mixed":
        rows[1]["calls"][0]["status"] = "mixed"
    else:
        rows[1]["alleles"]["locus2"] = None
    graph.set_merge_identical(True)
    graph.draw_results(rows, edges()[:1])
    assert len(graph.nodes) == 2


def test_visibility_halos_and_state_roundtrip(graph):
    graph.set_halos_visible(True)
    assert graph._halos and all(halo.isVisible() for _, halo in graph._halos)
    graph.nodes["a"].setPos(123.5, 456.25)
    graph.set_node_colors({"a": "#abcdef"})
    graph.set_node_label("a", "Display alias")
    graph.set_color_by("st")
    graph.set_labels_visible(False)
    graph.set_edge_labels_visible(False)
    graph.set_show_st(False)
    state = json.loads(json.dumps(graph.export_state()))
    graph.reset_layout()
    graph.reset_colors()
    graph.restore_state(state)
    assert graph.nodes["a"].pos().x() == 123.5
    assert graph.nodes["a"].pos().y() == 456.25
    assert graph.nodes["a"].slices == [("#abcdef", 1)]
    assert graph.labels["a"].text() == "Display alias"
    assert all(not item.isVisible() for item in graph.labels.values())
    assert all(not text.isVisible() for _, _, _, text in graph.edges)
    assert graph.color_by == "st"
    assert graph.export_state() == state


def test_display_labels_do_not_overlap_when_nodes_are_close(graph):
    for index, node in enumerate(graph.nodes.values()):
        node.setPos(index * 38, 0)
    graph.update_edges()
    rectangles = [item.sceneBoundingRect() for item in graph.labels.values()]
    for index, rectangle in enumerate(rectangles):
        assert all(not rectangle.intersects(other) for other in rectangles[:index])
    assert any(item.isVisible() for item in graph._label_guides.values())


@pytest.mark.parametrize("value", [[float("nan"), 0], [0, float("inf")], [0], ["zero", 0], [0, 1e99]])
def test_invalid_layout_positions_are_rejected(graph, value):
    with pytest.raises(ValueError, match="position"):
        graph.restore_state({"version": 1, "positions": {"a": value}})


def test_invalid_colors_labels_fields_and_duplicate_ids_are_rejected(graph):
    for action in (lambda: graph.set_selected_color("not a color"),
                   lambda: graph.set_node_colors({"unknown": "red"}),
                   lambda: graph.set_color_by("metadata:unknown"),
                   lambda: graph.set_palette([]),
                   lambda: graph.set_node_label("a", "x" * 121),
                   lambda: graph.set_node_label("a", "x\ny"),
                   lambda: graph.restore_state({"version": 99}),
                   lambda: graph.draw_results(records() + records()[:1], edges())):
        graph.select_ids(["a"])
        with pytest.raises(ValueError):
            action()


def test_png_svg_graphml_and_forest_newick_exports(graph, tmp_path):
    graph.set_color_by("metadata:location")
    graph.set_merge_identical(True)
    graph.set_halos_visible(True)
    graph.save_png(tmp_path / "graph.png")
    graph.save_svg(tmp_path / "graph.svg")
    graph.save_graphml(tmp_path / "graph.graphml")
    graph.save_newick(tmp_path / "graph.nwk")
    image = QImage(str(tmp_path / "graph.png"))
    assert (image.width(), image.height()) == (1800, 1200)
    assert image.pixelColor(0, 0) == QColor(BACKGROUND)
    svg = ElementTree.parse(tmp_path / "graph.svg").getroot()
    assert svg.tag.endswith("svg")
    document = ElementTree.parse(tmp_path / "graph.graphml")
    namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
    assert len(document.findall(".//g:node", namespace)) == 3
    assert len(document.findall(".//g:edge", namespace)) == 1
    members = document.find(".//g:node[@id='a']/g:data[@key='members']", namespace)
    assert json.loads(members.text) == ["a", "b"]
    newick = (tmp_path / "graph.nwk").read_text()
    assert "not a phylogeny" in newick
    assert sum(line.endswith(";") for line in newick.splitlines()) == 2
    assert "Isolate a [a]" in newick and "Isolate b [b]" in newick


def test_dark_palette_covers_native_controls_dialogs_and_navigation(qapp, qtbot):
    original_palette, original_style = qapp.palette(), qapp.styleSheet()
    original_native = qapp.testAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs)
    try:
        apply_dark_palette(qapp)
        qapp.setStyleSheet(STYLE)
        assert qapp.palette().color(QPalette.ColorRole.Window) == QColor(BACKGROUND)
        assert qapp.palette().color(QPalette.ColorRole.Text) == QColor(INK)
        dialog = QDialog()
        qtbot.addWidget(dialog)
        layout = QVBoxLayout(dialog)
        stack = QStackedWidget()
        layout.addWidget(stack)
        graph = TreeView()
        stack.addWidget(graph)
        graph.draw_results(records(), edges())
        controls = QWidget()
        form = QVBoxLayout(controls)
        for control in (QLineEdit("Sample name"), QDoubleSpinBox(), QComboBox(), QTableWidget(2, 2)):
            form.addWidget(control)
        stack.addWidget(controls)
        menu = QMenu(dialog)
        menu.addAction("Export graph")
        dialog.resize(850, 600)
        dialog.show()
        for index in [0, 1, 0, 1, 0]:
            stack.setCurrentIndex(index)
            dialog.resize(850 + index * 70, 600 + index * 50)
            qapp.processEvents()
            dialog.repaint()
            image = dialog.grab().toImage()
            assert not image.isNull()
            assert image.pixelColor(2, 2).lightness() < 50
        image = graph.viewport().grab().toImage()
        assert image.pixelColor(0, 0) == QColor(BACKGROUND)
        assert graph.graphicsEffect() is None
        assert graph.viewportUpdateMode() == QGraphicsView.ViewportUpdateMode.FullViewportUpdate
    finally:
        qapp.setPalette(original_palette)
        qapp.setStyleSheet(original_style)
        qapp.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, original_native)
