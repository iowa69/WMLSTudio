"""A forest in its own window: one identity, its own arrangement, the page left alone."""

import copy
import json
from xml.etree import ElementTree

import pytest
from PySide6.QtCore import QPoint, QPointF, QSettings, Qt
from PySide6.QtGui import QColor, QWheelEvent

import wmlstudio.graph_window as graph_window
from wmlstudio.graph_window import (
    EDIT_NOTE,
    HONESTY,
    MINIMUM_SIZE,
    GraphIdentity,
    GraphWindow,
    clamp_size,
    honesty_lines,
)
from wmlstudio.investigation import typing_scale
from wmlstudio.widgets import TreeView, graph_text_scale, set_graph_text_scale

CG_TARGETS = [f"target{index:04d}" for index in range(2358)]
MLST_LOCI = ["gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB"]


def cg_records():
    """Three isolates typed over a core-genome target set, one of them further away."""
    return [{"sample_id": key, "sample_name": f"Isolate {key}", "scheme": "kpneumoniae cgMLST",
             "scheme_digest": "cg-snapshot-1", "st": None, "primary_st": None,
             "alleles": {locus: "1" for locus in CG_TARGETS[:4]},
             "profile_loci": list(CG_TARGETS), "metadata": {"ward": ward}}
            for key, ward in (("a", "A"), ("b", "A"), ("c", "B"))]


def cg_edges():
    return [{"source": "a", "target": "b", "distance": 2, "shared_loci": 2301},
            {"source": "b", "target": "c", "distance": 30, "shared_loci": 2288}]


def mlst_records():
    return [{"sample_id": key, "sample_name": f"Isolate {key}", "scheme": "kpneumoniae 7-locus",
             "scheme_digest": "mlst-snapshot-1", "st": st, "primary_st": st,
             "alleles": {locus: "1" for locus in MLST_LOCI},
             "profile_loci": list(MLST_LOCI), "metadata": {"ward": ward}}
            for key, st, ward in (("a", "258", "A"), ("b", "258", "A"), ("c", "11", "B"))]


def mlst_edges():
    return [{"source": "a", "target": "b", "distance": 0, "shared_loci": 7},
            {"source": "b", "target": "c", "distance": 3, "shared_loci": 7}]


@pytest.fixture(autouse=True)
def forget_remembered_sizes(monkeypatch):
    """Window sizes are remembered between windows, never between tests."""
    monkeypatch.setattr(graph_window, "_SESSION_SIZES", {})


@pytest.fixture
def cg_view(qtbot):
    """A page's cgMLST forest: laid out with resize(), never shown."""
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(900, 600)
    view.set_scale(typing_scale(cg_records(), "cgmlst", scheme="kpneumoniae cgMLST"))
    view.draw_results(cg_records(), cg_edges(), 5)
    return view


@pytest.fixture
def mlst_view(qtbot):
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(900, 600)
    view.set_scale(typing_scale(mlst_records(), "mlst", scheme="kpneumoniae 7-locus"))
    view.draw_results(mlst_records(), mlst_edges(), 1)
    return view


@pytest.fixture
def window(qtbot, cg_view):
    detached = GraphWindow.from_view(cg_view, cohort="Ward B review",
                                     created="2026-09-15T08:30:00")
    qtbot.addWidget(detached)
    detached.resize(1100, 760)
    return detached


def test_a_window_draws_its_own_copy_and_never_takes_the_pages_graph_away(window, cg_view, qtbot):
    original = copy.deepcopy(cg_records())
    assert window.view is not cg_view
    assert sorted(window.view.nodes) == sorted(cg_view.nodes) == ["a", "b", "c"]
    assert window.view.cluster_threshold == cg_view.cluster_threshold == 5
    page_positions = {key: (node.pos().x(), node.pos().y()) for key, node in cg_view.nodes.items()}
    window.view.nodes["a"].setPos(321, 123)
    window.view.update_edges()
    window.view.set_node_colors({"a": "#ff0055"})
    window.view.set_node_label("a", "Index case")
    assert (cg_view.nodes["a"].pos().x(), cg_view.nodes["a"].pos().y()) == page_positions["a"]
    assert cg_view.nodes["a"].slices != [("#ff0055", 1)]
    assert cg_view.node_label("a") == "Isolate a"
    closes = []
    window.closed.connect(lambda: closes.append(True))
    window.close()
    assert closes == [True]
    # The page keeps its graph, its arrangement and its evidence after a window closes.
    assert {key: (node.pos().x(), node.pos().y())
            for key, node in cg_view.nodes.items()} == page_positions
    cg_view.fit_tree()
    assert cg_view.graph_contents()["results"] == original


def test_two_windows_on_one_comparison_keep_their_own_threshold_and_arrangement(cg_view, qtbot):
    first = GraphWindow.from_view(cg_view, cohort="Ward B review", created="2026-09-15T08:30:00")
    qtbot.addWidget(first)
    # The page moves on to a wider link threshold; the open window must not follow.
    cg_view.draw_results(cg_records(), cg_edges(), 40)
    second = GraphWindow.from_view(cg_view, cohort="Ward B review", created="2026-09-15T09:05:00")
    qtbot.addWidget(second)
    assert first.windowTitle() != second.windowTitle()
    assert "link ≤ 5 of 2358 targets" in first.windowTitle()
    assert "link ≤ 40 of 2358 targets" in second.windowTitle()
    assert (first.view.cluster_threshold, second.view.cluster_threshold) == (5, 40)
    assert len(first.view.groups()) == 2 and len(second.view.groups()) == 1
    first.view.nodes["a"].setPos(11, 22)
    first.view.update_edges()
    first.view.set_node_colors({"a": "#ff0055"})
    assert (second.view.nodes["a"].pos().x(), second.view.nodes["a"].pos().y()) != (11, 22)
    assert second.view.nodes["a"].slices != [("#ff0055", 1)]
    first.close()
    assert second.view.nodes["a"].slices != [("#ff0055", 1)]
    second.fit()
    assert sorted(second.view.nodes) == ["a", "b", "c"]


def test_a_classical_window_and_a_core_genome_window_never_borrow_each_others_numbers(
        qtbot, mlst_view, cg_view):
    classical = GraphWindow.from_view(mlst_view, cohort="Same three isolates")
    core = GraphWindow.from_view(cg_view, cohort="Same three isolates")
    for detached in (classical, core):
        qtbot.addWidget(detached)
    assert classical.kind_badge.text() == "Classical MLST"
    assert core.kind_badge.text() == "cgMLST"
    assert "7 loci" in classical.headline.text() and "2358" not in classical.headline.text()
    assert "2358 targets" in core.headline.text() and "7 loci" not in core.headline.text()
    assert classical.identity.threshold_words() == "link ≤ 1 of 7 loci"
    assert core.identity.threshold_words() == "link ≤ 5 of 2358 targets"
    assert classical.identity.export_subtitle() != core.identity.export_subtitle()
    assert classical.windowTitle() != core.windowTitle()
    for detached in (classical, core):
        assert not detached.separation.isHidden()
        assert "never share a scale" in detached.separation.text()


def test_the_window_is_a_resizable_window_that_remembers_the_size_it_was_left_at(
        qtbot, cg_view, mlst_view, tmp_path):
    preferences = QSettings(str(tmp_path / "Interface.ini"), QSettings.Format.IniFormat)
    first = GraphWindow.from_view(cg_view, cohort="Ward B review", preferences=preferences)
    qtbot.addWidget(first)
    assert first.windowFlags() & Qt.WindowType.Window
    assert (first.minimumWidth(), first.minimumHeight()) == MINIMUM_SIZE
    first.resize(980, 640)
    first.close()
    assert preferences.value("graph_window/cgmlst/size") == "980x640"
    reopened = GraphWindow.from_view(cg_view, cohort="Ward B review", preferences=preferences)
    qtbot.addWidget(reopened)
    assert (reopened.width(), reopened.height()) == clamp_size(980, 640,
                                                              reopened.available_screen_size())
    # A window for the other typing kind keeps its own size rather than inheriting one.
    other = GraphWindow.from_view(mlst_view, cohort="Ward B review", preferences=preferences)
    qtbot.addWidget(other)
    assert (other.width(), other.height()) == clamp_size(*graph_window.DEFAULT_SIZE,
                                                         other.available_screen_size())
    # A size larger than the screen is brought back inside it: an unreachable
    # window is worse than a small one.
    assert clamp_size(20000, 20000, (1024, 768)) == (1024, 768)
    assert clamp_size(10, 10, (1024, 768)) == MINIMUM_SIZE
    reopened.show()
    reopened.showMaximized()
    qtbot.waitUntil(reopened.isMaximized)
    reopened.showNormal()
    qtbot.waitUntil(lambda: not reopened.isMaximized())
    reopened.close()


def test_nodes_stay_where_they_are_put_and_their_edges_follow_them(window, qtbot):
    view = window.view
    source, _target, line, text = view.edges[0]
    before = (line.line().x1(), line.line().y1())
    view.nodes[source].setPos(400, 250)
    assert (line.line().x1(), line.line().y1()) == (400, 250) != before
    assert text.pos() != QPointF(0, 0)
    # Every control that redraws the forest leaves the node where it was put.
    window.toggles["st"].setChecked(False)
    window.color_field.setCurrentIndex(window.color_field.findData("st"))
    window.toggles["merge"].setChecked(True)
    assert (view.nodes[source].pos().x(), view.nodes[source].pos().y()) == (400, 250)
    assert (view.edges[0][2].line().x1(), view.edges[0][2].line().y1()) == (400, 250)
    state = json.loads(json.dumps(window.presentation_state()))
    assert state["view"]["positions"][source] == [400.0, 250.0]
    assert state["identity"] == {"kind": "cgmlst", "title": "cgMLST", "scheme": "kpneumoniae cgMLST",
                                 "targets": 2358, "target_word": "targets", "cohort": "Ward B review",
                                 "threshold": 5, "created": "2026-09-15 08:30", "note": "",
                                 "difference_word": "allele differences",
                                 "distance_phrase": "allele-distance"}
    restored = GraphWindow(window.identity)
    qtbot.addWidget(restored)
    restored.set_contents(view.graph_contents())
    restored.restore_presentation(state)
    assert (restored.view.nodes[source].pos().x(), restored.view.nodes[source].pos().y()) == (400, 250)


def test_colours_and_labels_can_be_set_per_node_and_per_group_without_touching_the_evidence(
        window, monkeypatch):
    monkeypatch.setattr("wmlstudio.graph_window.QColorDialog.getColor",
                        lambda *args, **kwargs: QColor("#ff0055"))
    monkeypatch.setattr("wmlstudio.graph_window.QInputDialog.getText",
                        lambda *args, **kwargs: ("Ward B index case", True))
    window.color_selected()
    assert "Select one or more isolates" in window.statusBar().currentMessage()
    window.view.select_ids(["a"])
    window.color_selected()
    assert window.view.nodes["a"].slices == [("#ff0055", 1)]
    assert window.view.nodes["b"].slices != [("#ff0055", 1)]
    # "a" and "b" are one single-link group at this threshold; "c" is not.
    window.color_group()
    assert window.view.nodes["b"].slices == [("#ff0055", 1)]
    assert window.view.nodes["c"].slices != [("#ff0055", 1)]
    assert "not an outbreak assignment" in window.statusBar().currentMessage()
    window.view.select_ids(["a"])
    window.rename_selected()
    assert window.view.node_label("a") == "Ward B index case"
    assert "Ward B index case" in window.view.labels["a"].text()
    # The original sample name is retained and still shown for the node.
    assert window.view.labels["a"].toolTip() == "Isolate a"
    assert window.view.graph_contents()["results"][0]["sample_name"] == "Isolate a"
    window.view.select_ids(["a", "b"])
    window.rename_selected()
    assert "exactly one isolate" in window.statusBar().currentMessage()
    window.reset_colors()
    assert window.view.nodes["a"].slices == window.view.nodes["b"].slices


def test_an_arrangement_is_handed_back_for_saving_without_a_redraw_storm(window, qtbot):
    with qtbot.waitSignal(window.stateChanged, timeout=2000) as handed_back:
        window.view.nodes["a"].setPos(150, 90)
    assert handed_back.args[0]["positions"]["a"] == [150.0, 90.0]
    saved = []
    window.stateChanged.connect(saved.append)
    for step in range(5):
        window.view.nodes["b"].setPos(step * 10, 40)
    qtbot.waitUntil(lambda: bool(saved), timeout=2000)
    assert len(saved) == 1
    assert saved[0]["positions"]["b"] == [40.0, 40.0]
    # An arrangement made in the last moment before closing is not lost.
    window.view.nodes["c"].setPos(77, 88)
    window.close()
    assert saved[-1]["positions"]["c"] == [77.0, 88.0]


def test_an_export_carries_this_windows_identity_at_the_standard_lettering_size(window, tmp_path):
    previous = graph_text_scale()
    observed = []
    original = window.view.canvas.render
    try:
        window.view.canvas.render = lambda *args, **kwargs: (observed.append(graph_text_scale()),
                                                             original(*args, **kwargs))[1]
        set_graph_text_scale(150)
        window.export(tmp_path / "forest.png")
        assert observed == [100], f"export painted the scene at reader scale {observed}"
        assert graph_text_scale() == 150
    finally:
        del window.view.canvas.render
        set_graph_text_scale(previous)
    window.export(tmp_path / "forest.svg")
    svg = (tmp_path / "forest.svg").read_text(encoding="utf-8")
    for expected in ("cgMLST", "kpneumoniae cgMLST", "2358 targets", "Ward B review",
                     "link ≤ 5 of 2358 targets, single linkage",
                     "not a phylogeny or transmission chain"):
        assert expected in svg, expected
    window.export(tmp_path / "forest.graphml")
    document = ElementTree.parse(tmp_path / "forest.graphml")
    namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
    interpretation = document.find(".//g:graph/g:data[@key='interpretation']", namespace).text
    assert "cgMLST" in interpretation and "Not a phylogeny" in interpretation
    window.export(tmp_path / "forest.nwk")
    assert "not a phylogeny" in (tmp_path / "forest.nwk").read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported graph export format"):
        window.export(tmp_path / "forest.txt")


def test_saving_from_the_window_writes_what_the_reader_chose_and_nothing_when_cancelled(
        window, tmp_path, monkeypatch):
    target = tmp_path / "chosen.svg"
    monkeypatch.setattr("wmlstudio.graph_window.QFileDialog.getSaveFileName",
                        lambda *args, **kwargs: (str(target), "SVG vector (*.svg)"))
    assert window.export_dialog() == target
    assert target.stat().st_size > 0
    assert "Saved chosen.svg · cgMLST" in window.statusBar().currentMessage()
    # A name typed without an extension is saved in the format that was chosen.
    monkeypatch.setattr("wmlstudio.graph_window.QFileDialog.getSaveFileName",
                        lambda *args, **kwargs: (str(tmp_path / "no extension"),
                                                 "GraphML nodes and edges (*.graphml)"))
    assert window.export_dialog() == tmp_path / "no extension.graphml"
    monkeypatch.setattr("wmlstudio.graph_window.QFileDialog.getSaveFileName",
                        lambda *args, **kwargs: ("", ""))
    assert window.export_dialog() is None
    assert sorted(path.name for path in tmp_path.iterdir()) == ["chosen.svg",
                                                               "no extension.graphml"]


def test_the_host_can_stop_an_export_from_overwriting_an_input_file(window, tmp_path):
    reads = tmp_path / "isolate-a.fasta"
    reads.write_text(">a\nACGT\n", encoding="utf-8")

    def refuse(path):
        raise ValueError("That is an input file; choose another destination.")

    window.check_output = refuse
    with pytest.raises(ValueError, match="input file"):
        window.export(reads)
    assert reads.read_text(encoding="utf-8") == ">a\nACGT\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["isolate-a.fasta"]


def test_the_graphs_honesty_stays_with_the_window(window):
    text = window.honesty.text()
    assert all(line in text for line in HONESTY)
    assert "single linkage" in text and "larger than the threshold" in text
    assert "loci the two isolates actually share" in text and "never a zero distance" in text
    assert "not a phylogeny" in text and "not a transmission chain" in text
    assert "never share a scale" in window.separation.text()
    assert EDIT_NOTE in window.statusBar().currentMessage()
    assert "picture only" in EDIT_NOTE and "Allele calls, distances" in EDIT_NOTE


def test_the_forest_is_fitted_to_the_window_until_the_reader_moves_the_view(window, qtbot):
    assert window.follow_window is True
    zoom = window.view.transform().m11()
    wheel = QWheelEvent(QPointF(20, 20), QPointF(20, 20), QPoint(0, 0), QPoint(0, 120),
                        Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                        Qt.ScrollPhase.NoScrollPhase, False)
    window.view.wheelEvent(wheel)
    assert window.view.transform().m11() > zoom
    assert window.follow_window is False
    chosen = window.view.transform().m11()
    window.resize(940, 700)
    qtbot.wait(150)
    assert window.view.transform().m11() == chosen
    window.fit()
    assert window.follow_window is True


def test_the_controls_show_what_the_view_shows_and_change_presentation_only(window, cg_view):
    fields = [window.color_field.itemData(index) for index in range(window.color_field.count())]
    assert fields == window.view.available_color_fields() == ["cluster", "st", "metadata:ward"]
    assert window.color_field.itemText(fields.index("metadata:ward")) == "Ward (recorded detail)"
    window.color_field.setCurrentIndex(fields.index("metadata:ward"))
    assert window.view.color_by == "metadata:ward"
    assert set(window.view.legend()) == {"A", "B"}
    assert cg_view.color_by == "cluster"
    assert all(box.isChecked() for name, box in window.toggles.items() if name != "merge")
    window.toggles["edges"].setChecked(False)
    assert all(not text.isVisible() for _a, _b, _line, text in window.view.edges)
    window.toggles["halos"].setChecked(False)
    assert window.view.show_halos is False
    window.toggles["labels"].setChecked(False)
    assert all(not item.isVisible() for item in window.view.labels.values())
    window.toggle_mode()
    assert window.view.interaction_mode == "pan"
    assert window.mode_button.text() == "Select isolates"
    assert window.view.graph_contents()["edges"] == cg_view.graph_contents()["edges"]


SNP_SCALE = {"kind": "snp", "title": "SKA2 split k-mer SNPs", "target_word": "shared split k-mers",
             "difference_word": "SNPs", "distance_phrase": "SNP-distance", "targets": 0,
             "scheme": "ska2:assembly-split-kmer-k31",
             "caption": "SKA2 split k-mer SNPs · reference-free, k=31 · 3 isolates",
             "denominator_note": "Every pair carries its own denominator: the split k-mers those "
                                 "two isolates share.",
             "separation": "A SNP distance, a classical MLST allele distance and a cgMLST target "
                           "distance are three different quantities."}


@pytest.fixture
def snp_view(qtbot):
    """A SNP forest: no cohort target count, a per-pair denominator, no link threshold."""
    view = TreeView()
    qtbot.addWidget(view)
    view.resize(900, 600)
    view.show_contents({
        "results": [{"sample_id": key, "sample_name": f"Isolate {key}", "st": None}
                    for key in ("a", "b", "c")],
        "edges": [{"source": "a", "target": "b", "distance": 3,
                   "denominator_label": "19,880 split k-mers shared, 99.1% of the pair's "
                                        "combined set"}],
        "cluster_threshold": -1, "scale": dict(SNP_SCALE)})
    return view


def test_a_snp_window_states_the_quantity_it_shows_and_never_calls_it_an_allele(snp_view, qtbot):
    window = GraphWindow.from_view(snp_view, cohort="Ward B review",
                                   created="2026-09-15T08:30:00")
    qtbot.addWidget(window)
    assert window.identity.kind == "snp"
    assert window.kind_badge.text() == "SKA2 split k-mer SNPs"
    # A per-pair denominator is not a missing target count, and is not reported as one.
    assert window.identity.caption() == SNP_SCALE["caption"]
    assert "target count not recorded" not in window.headline.text()
    assert window.identity.threshold_words() == "no link threshold set"
    assert window.windowTitle().startswith("SKA2 split k-mer SNPs forest · no link threshold set "
                                           "· Ward B review · ska2:assembly-split-kmer-k31")
    assert window.identity.export_title().endswith("· SNP-distance minimum spanning forest")
    assert "single linkage" not in window.identity.export_subtitle()
    assert "This is a layout of SNPs." in window.honesty.text()
    assert "own denominator" in window.honesty.text()
    assert "never a zero distance" in window.honesty.text()
    assert "allele" not in window.honesty.text().casefold()
    assert window.toggles["edges"].text() == "SNPs on edges"
    assert "SNPs" in window.toggles["edges"].toolTip()
    # A classical window opened beside it keeps every one of its own words.
    assert honesty_lines(GraphIdentity.from_scale(None)) == HONESTY
    window.close()


def test_an_identity_with_nothing_recorded_says_so_rather_than_guessing():
    identity = GraphIdentity.from_scale(None)
    assert identity.caption() == ("Unclassified typing · reference not recorded · "
                                  "target count not recorded")
    assert identity.threshold_words() == "link threshold not recorded"
    assert identity.cohort_words() == "cohort not named"
    assert "not a phylogeny or transmission chain" in identity.export_subtitle()
    # A threshold without a recorded target set never invents a denominator.
    bare = GraphIdentity.from_scale(None, threshold=3)
    assert bare.threshold_words() == "link ≤ 3 allele differences"
    named = GraphIdentity.from_scale(typing_scale(cg_records(), "cgmlst", scheme="kp cgMLST"),
                                     cohort="Ward B", threshold=10, created="2026-09-15T08:30:00")
    assert named.caption() == "cgMLST · kp cgMLST · 2358 targets"
    assert named.window_title().startswith("cgMLST forest · link ≤ 10 of 2358 targets · Ward B")
    assert named.created == "2026-09-15 08:30"
