"""Native graph interaction, stable presentation state, exports, and dark surfaces."""

import copy
import json
import math
from xml.etree import ElementTree

import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor, QContextMenuEvent, QImage, QPalette
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

import wmlstudio.widgets
from wmlstudio.theme import BACKGROUND, INK, STYLE, apply_dark_palette
from wmlstudio.widgets import (
    GRAPH_SUBTITLE,
    TreeView,
    forest_layout,
    graph_colors,
    render_side_by_side,
)


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
    # The ground of the picture is the ground of the theme in force, which is why
    # it is asked for here rather than written down: see test_graph_theme.py.
    assert image.pixelColor(0, 0) == QColor(graph_colors()["ground"])
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


def test_a_forest_can_be_handed_to_another_view_as_data_without_sharing_its_items(graph, qtbot):
    """A second view of one comparison must be a copy, not the same nodes twice.

    A detached window draws the same records, edges, threshold, groups and scale
    in a view of its own, so arranging one picture never moves a node in the
    other and neither can edit the evidence they were both built from.
    """
    graph.set_scale({"kind": "cgmlst", "title": "cgMLST", "target_word": "targets",
                     "targets": 2358, "caption": "cgMLST · kp · 2358 targets"})
    contents = graph.graph_contents()
    assert sorted(row["sample_id"] for row in contents["results"]) == ["a", "b", "c", "d"]
    assert contents["cluster_threshold"] == 1
    assert contents["scale"]["targets"] == 2358
    contents["results"][0]["sample_name"] = "Edited copy"
    contents["scale"]["targets"] = 7
    assert graph._results["a"]["sample_name"] == "Isolate a"
    assert graph.scale["targets"] == 2358
    other = TreeView()
    qtbot.addWidget(other)
    other.resize(900, 600)
    other.show_contents(graph.graph_contents())
    assert sorted(other.nodes) == sorted(graph.nodes)
    assert other.cluster_threshold == graph.cluster_threshold
    assert other.scale_caption() == graph.scale_caption()
    other.nodes["a"].setPos(500, 400)
    assert (graph.nodes["a"].pos().x(), graph.nodes["a"].pos().y()) != (500, 400)


def test_groups_labels_and_membership_are_readable_without_reaching_inside_the_view(graph):
    groups = graph.groups()
    assert [group["status"] for group in groups] == ["cluster", "singleton"]
    # Single linkage chains: a and c are one group through b, not by their own distance.
    assert graph.group_members("a") == graph.group_members("c") == ["a", "b", "c"]
    assert graph.group_members("d") == ["d"]
    # An isolate that is in no drawn group answers for itself, never for a group.
    assert graph.group_members("unknown") == ["unknown"]
    groups[0]["name"] = "Renamed elsewhere"
    assert graph.groups()[0]["name"] == "Group 1"
    assert graph.node_label("a") == "Isolate a"
    graph.set_node_label("a", "Ward B index case")
    assert graph.node_label("a") == "Ward B index case"
    assert graph.labels["a"].toolTip() == "Isolate a"
    with pytest.raises(ValueError, match="Unknown sample ID"):
        graph.node_label("unknown")


def test_scrolling_over_the_graph_zooms_it_and_says_the_reader_moved_the_view(graph):
    """The typing scale this view measures shadows the view's own zoom method.

    Calling the shadowed name turned every scroll over the graph into a type
    error instead of a zoom, so the magnification is asserted here rather than
    only the signal that follows it.
    """
    from PySide6.QtCore import QPoint, QPointF
    from PySide6.QtGui import QWheelEvent

    def wheel(delta):
        return QWheelEvent(QPointF(30, 30), QPointF(30, 30), QPoint(0, 0), QPoint(0, delta),
                           Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                           Qt.ScrollPhase.NoScrollPhase, False)

    adjustments = []
    graph.viewAdjusted.connect(lambda: adjustments.append(graph.transform().m11()))
    start = graph.transform().m11()
    graph.wheelEvent(wheel(120))
    assert graph.transform().m11() > start
    graph.wheelEvent(wheel(-120))
    assert graph.transform().m11() == pytest.approx(start)
    assert len(adjustments) == 2
    for _ in range(80):
        graph.wheelEvent(wheel(-120))
    assert 0.05 <= graph.transform().m11() <= 8
    # Fitting the graph is not the reader moving it, so it reports no adjustment.
    before = len(adjustments)
    graph.fit_tree()
    assert len(adjustments) == before


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
        assert image.pixelColor(0, 0) == QColor(graph_colors()["ground"])
        assert graph.graphicsEffect() is None
        assert graph.viewportUpdateMode() == QGraphicsView.ViewportUpdateMode.FullViewportUpdate
    finally:
        qapp.setPalette(original_palette)
        qapp.setStyleSheet(original_style)
        qapp.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, original_native)


@pytest.fixture
def quiet_graph(qtbot):
    """A laid-out but never shown forest: repainting one after many torn-down
    offscreen widgets crashes inside the graph label painter."""
    widget = TreeView()
    qtbot.addWidget(widget)
    widget.resize(900, 600)
    widget.draw_results(records(), edges())
    return widget


def test_an_export_can_name_which_tree_it_is_without_changing_the_default(quiet_graph, tmp_path):
    default_path, titled_path = tmp_path / "default.png", tmp_path / "titled.png"
    quiet_graph.save_image(default_path)
    quiet_graph.save_image(titled_path, title="Baseline forest",
                           subtitle="Frozen 2026-01-01 · not a phylogeny or transmission tree")
    assert default_path.stat().st_size > 0 and titled_path.stat().st_size > 0
    assert default_path.read_bytes() != titled_path.read_bytes()
    quiet_graph.save_svg(tmp_path / "titled.svg", title="Baseline forest")
    assert "Baseline forest" in (tmp_path / "titled.svg").read_text()
    quiet_graph.save_svg(tmp_path / "default.svg")
    assert "Allele-distance minimum spanning forest" in (tmp_path / "default.svg").read_text()
    assert GRAPH_SUBTITLE in (tmp_path / "default.svg").read_text()


def test_two_forests_can_be_shown_side_by_side_with_the_limit_in_the_picture(quiet_graph, qtbot, tmp_path):
    other = TreeView()
    qtbot.addWidget(other)
    other.resize(900, 600)
    other.draw_results(records()[:2], edges()[:1])
    image = render_side_by_side(quiet_graph, other, left_title="Baseline · 4 isolates",
                                right_title="Current · 2 isolates",
                                headline="Same reference, same link threshold", width=900, height=600)
    assert (image.width(), image.height()) == (1800, 676)
    for suffix in ("png", "jpg"):
        path = tmp_path / f"pair.{suffix}"
        assert image.save(str(path))
        assert path.stat().st_size > 0
    assert QImage(str(tmp_path / "pair.png")).width() == 1800


SNP_SCALE = {"kind": "snp", "title": "SKA2 split k-mer SNPs", "target_word": "shared split k-mers",
             "difference_word": "SNPs", "distance_phrase": "SNP-distance", "targets": 0,
             "caption": "SKA2 split k-mer SNPs · reference-free, k=31 · 2 isolates",
             "separation": "A SNP distance and an allele distance are different quantities."}
SNP_EDGE = {"source": "a", "target": "b", "distance": 3,
            "denominator_label": "19,880 split k-mers shared, 99.1% of the pair's combined set"}


def test_a_forest_that_is_not_allele_typing_takes_every_word_from_its_own_scale(quiet_graph, tmp_path):
    """A SNP edge carries no shared_loci at all: its denominator belongs to the pair.

    Reading one off the edge raised KeyError before this, so the tooltip is
    asserted here rather than only the wording, and every other place the view
    names a quantity is asserted beside it.
    """
    quiet_graph.show_contents({"results": records()[:2], "edges": [dict(SNP_EDGE)],
                               "cluster_threshold": -1, "scale": dict(SNP_SCALE)})
    assert "shared_loci" not in SNP_EDGE
    tooltip = quiet_graph._edge_tooltip(SNP_EDGE)
    assert tooltip.startswith("3 SNPs / 19,880 split k-mers shared, 99.1% of the pair's combined set")
    assert quiet_graph.edges[0][2].toolTip() == tooltip
    assert "SNP-distance minimum spanning forest; not a phylogeny." in quiet_graph._view_tooltip()
    assert quiet_graph.graph_subtitle().endswith("edge labels are SNPs")
    # Nothing is grouped below a negative threshold, so nothing is outlined and no
    # group is numbered: seven halos labelled "Group 1" read as seven clusters.
    assert quiet_graph._halos == [] and quiet_graph._halo_labels == {}
    assert [group["name"] for group in quiet_graph.groups()] == ["Not grouped", "Not grouped"]
    quiet_graph.save_svg(tmp_path / "snp.svg")
    written = (tmp_path / "snp.svg").read_text()
    assert "SNP-distance minimum spanning forest" in written
    assert "edge labels are SNPs" in written
    assert "single-link threshold: none set" in written
    assert "allele" not in written.casefold()
    quiet_graph.save_newick(tmp_path / "snp.nwk")
    assert "SNP-distance MST topology; not a phylogeny" in (tmp_path / "snp.nwk").read_text()
    quiet_graph.save_graphml(tmp_path / "snp.graphml")
    document = ElementTree.parse(tmp_path / "snp.graphml")
    namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
    assert [key.get("id") for key in document.findall(".//g:key[@for='edge']", namespace)] == [
        "distance", "unit", "shared_denominator"]
    assert {data.get("key"): data.text
            for data in document.findall(".//g:edge/g:data", namespace)} == {
        "distance": "3", "unit": "SNPs", "shared_denominator": SNP_EDGE["denominator_label"]}
    # An empty forest names what it is not showing rather than asking for profiles.
    quiet_graph.show_contents({"results": [], "edges": [], "scale": dict(SNP_SCALE)})
    assert any("No SKA2 split k-mer SNPs comparison is drawn yet." == item.text()
               for item in quiet_graph.canvas.items() if hasattr(item, "text"))


def test_an_allele_forest_keeps_every_word_it_had_before_a_scale_could_rename_them(graph, tmp_path):
    """The defaults are the wording this view has always used, so nothing moved."""
    assert graph._edge_tooltip(edges()[1]).startswith("1 differing alleles / 2 shared loci")
    assert "Allele-distance minimum spanning forest; not a phylogeny." in graph._view_tooltip()
    assert graph.graph_subtitle() == GRAPH_SUBTITLE
    assert [group["name"] for group in graph.groups()] == ["Group 1", "Group 2"]
    assert graph._halos and "Single-link group at ≤ 1 allele differences" in graph._halos[0][1].toolTip()
    graph.set_scale({"targets": 2358, "target_word": "targets", "caption": "cgMLST · kp · 2358 targets"})
    assert graph._edge_tooltip(edges()[1]).startswith("1 differing alleles / 2 of 2358 targets shared")
    graph.save_svg(tmp_path / "allele.svg")
    written = (tmp_path / "allele.svg").read_text()
    assert "cgMLST · kp · 2358 targets · allele-distance minimum spanning forest" in written
    assert GRAPH_SUBTITLE in written
    assert "single-link threshold: 1 of 2358 targets" in written


def test_a_pinned_colour_is_shared_and_pinning_the_same_key_twice_is_quiet(quiet_graph):
    quiet_graph.set_color_by("st")
    categories = quiet_graph.color_categories()
    assert set(categories) == {"ST 184", "ST unassigned"}
    emissions = []
    quiet_graph.legendChanged.connect(lambda values: emissions.append(dict(values)))
    assert quiet_graph.set_pinned_legend({"ST 184": "#123456"}) is True
    assert quiet_graph.legend()["ST 184"] == "#123456"
    assert len(emissions) == 1
    assert quiet_graph.set_pinned_legend({"ST 184": "#123456"}) is False
    assert len(emissions) == 1
    # A category this view does not show is ignored rather than invented.
    quiet_graph.set_pinned_legend({"ST 999": "#abcdef"})
    assert "ST 999" not in quiet_graph.legend()


def test_a_redraw_reuses_known_positions_and_only_a_new_isolate_is_laid_out_again(quiet_graph, monkeypatch):
    calls = []
    original = wmlstudio.widgets.forest_layout

    def counted(keys, links):
        calls.append(sorted(keys))
        return original(keys, links)

    monkeypatch.setattr(wmlstudio.widgets, "forest_layout", counted)
    quiet_graph.nodes["a"].setPos(40, 50)
    quiet_graph.update_edges()
    quiet_graph._redraw()
    assert calls == []
    assert quiet_graph.nodes["a"].pos().x() == 40 and quiet_graph.nodes["a"].pos().y() == 50
    quiet_graph.reset_layout()
    assert len(calls) == 1
    fresh = records() + [{"sample_id": "e", "sample_name": "Isolate e", "scheme_digest": "same-version",
                          "alleles": {"locus1": "1", "locus2": "9"}, "calls": [], "metadata": {}}]
    quiet_graph.draw_results(fresh, edges())
    assert len(calls) == 2 and "e" in calls[-1]


class RecordingMenu:
    """Stands in for the modal popup: records its entries and chooses one.

    The real QMenu.exec blocks on a native event loop, which an offscreen test
    can never dismiss, so the menu class the view constructs is replaced instead.
    """

    picked = ""

    def __init__(self, parent=None):
        self.entries = []

    def addAction(self, text):
        action = QAction(str(text))
        self.entries.append(action)
        return action

    def addSeparator(self):
        self.entries.append(None)

    def titles(self):
        return [action.text() for action in self.entries if action is not None]

    def exec(self, position):
        return next((action for action in self.entries
                     if action is not None and self.picked in action.text()), None)

    def deleteLater(self):
        pass


def test_the_graph_menu_can_be_extended_without_losing_its_own_entries(quiet_graph, monkeypatch):
    seen, chosen, menus = {}, [], []

    def extension(view, menu, node_key):
        seen["node_key"] = node_key
        menu.addSeparator()
        action = menu.addAction("Archive this isolate (keeps all evidence)…")
        return {action: lambda: chosen.append(node_key)}

    class Chooser(RecordingMenu):
        picked = "Archive"

        def __init__(self, parent=None):
            super().__init__(parent)
            menus.append(self)

    quiet_graph.context_extension = extension
    monkeypatch.setattr(wmlstudio.widgets, "QMenu", Chooser)
    quiet_graph.contextMenuEvent(_context_event(quiet_graph, quiet_graph.nodes["a"]))
    assert seen["node_key"] == "a"
    titles = menus[0].titles()
    assert "Select this threshold group" in titles and "Fit graph" in titles
    assert titles[-1] == "Archive this isolate (keeps all evidence)…"
    assert chosen == ["a"]


def test_an_extension_entry_never_shadows_the_views_own_actions(quiet_graph, monkeypatch):
    ran = []
    quiet_graph.context_extension = lambda view, menu, key: {menu.addAction("Archive…"): ran.append}

    class Chooser(RecordingMenu):
        picked = "Fit graph"

    monkeypatch.setattr(wmlstudio.widgets, "QMenu", Chooser)
    quiet_graph.contextMenuEvent(_context_event(quiet_graph, quiet_graph.nodes["a"]))
    assert ran == []


def _context_event(view, node):
    point = view.mapFromScene(node.sceneBoundingRect().center())
    return QContextMenuEvent(QContextMenuEvent.Reason.Mouse, point, view.mapToGlobal(point))


def test_graph_text_scale_resizes_lettering_without_moving_the_tree():
    """A bigger label must not become a different distance.

    Rescaling text for a high-resolution screen is a readability change. If it
    also moved nodes or changed an edge, the same data would draw as a different
    tree, so only the point size may follow the setting.
    """
    from wmlstudio.widgets import _graph_font, graph_text_scale, set_graph_text_scale
    previous = graph_text_scale()
    try:
        assert set_graph_text_scale(100) == 100
        base = _graph_font(14).pointSize()
        assert set_graph_text_scale(150) == 150
        assert _graph_font(14).pointSize() > base
        assert set_graph_text_scale(80) == 80
        assert _graph_font(14).pointSize() < base
        # Out-of-range requests are clamped, never applied and never raised at the user.
        assert set_graph_text_scale(500) == 150
        assert set_graph_text_scale(10) == 80
        # Even at the smallest scale a label stays legible rather than collapsing.
        assert _graph_font(6).pointSize() >= 6
    finally:
        set_graph_text_scale(previous)


def test_exported_image_keeps_the_standard_text_size_whatever_the_reader_chose(qtbot, tmp_path):
    """The Settings page promises a saved picture looks the same on every computer.

    Graph text size is a per-reader comfort setting. If it leaked into exports,
    two people reporting the same comparison would produce different figures.
    """
    from wmlstudio.widgets import TreeView, graph_text_scale, set_graph_text_scale
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(900, 600)
    view.draw_results([
        {"sample_id": "a", "sample_name": "A", "st": "1", "alleles": {"x": "1"}, "status": "completed"},
        {"sample_id": "b", "sample_name": "B", "st": "2", "alleles": {"x": "2"}, "status": "completed"},
    ], [], 10)
    previous = graph_text_scale()
    observed = []
    original = view.canvas.render
    try:
        # Record the scale in force at the moment the scene is painted into the
        # export, rather than comparing rendered bytes: node positions carry over
        # between redraws, so image equality would assert layout, not text size.
        view.canvas.render = lambda *args, **kwargs: (observed.append(graph_text_scale()),
                                                      original(*args, **kwargs))[1]
        set_graph_text_scale(150)
        view.save_image(tmp_path / "enlarged.png")
        assert observed == [100], f"export painted the scene at reader scale {observed}"
        # The reader's own setting survives the export untouched.
        assert graph_text_scale() == 150
    finally:
        del view.canvas.render
        set_graph_text_scale(previous)
