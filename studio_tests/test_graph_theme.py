"""The forest follows the chosen theme, on screen and in every picture saved from it."""

import pytest
from PySide6.QtGui import QColor, QImage

from wmlstudio import theme
from wmlstudio.theme import BACKGROUND, INK, MUTED, PALETTE
from wmlstudio.widgets import (
    TreeView,
    contrast_ratio,
    graph_colors,
    halo_colors,
    readable_ink,
    render_side_by_side,
)


def records():
    return [{"sample_id": key, "sample_name": f"Isolate {key}", "scheme_digest": "same-version",
             "st": "184", "alleles": {"locus1": "1", "locus2": allele},
             "calls": [{"locus": "locus1", "allele": "1", "status": "exact"},
                       {"locus": "locus2", "allele": allele, "status": "exact"}],
             "metadata": {"location": "Ward A"}}
            for key, allele in [("a", "1"), ("b", "1"), ("c", "2")]]


def edges():
    return [{"source": "a", "target": "b", "distance": 0, "shared_loci": 2},
            {"source": "b", "target": "c", "distance": 1, "shared_loci": 2}]


@pytest.fixture
def theme_choice():
    """Choose a theme for one test, and put the reader's own choice back after it."""
    previous = theme.active_theme()
    yield lambda name: theme.set_active_theme(name)
    theme.set_active_theme(previous)


@pytest.fixture
def forest(qtbot):
    def build():
        view = TreeView()
        qtbot.addWidget(view)
        view.resize(900, 600)
        view.draw_results(records(), edges(), 1)
        return view
    return build


def _over(color, ground):
    """A translucent colour as it actually appears once painted on a ground."""
    color, ground = QColor(color), QColor(ground)
    weight = color.alphaF()
    return QColor(*[round(getattr(color, channel)() * weight + getattr(ground, channel)() * (1 - weight))
                    for channel in ("red", "green", "blue")])


def _strongest_contrast(image, ground, box):
    """How well the darkest or lightest lettering in a band reads on the ground."""
    left, top, right, bottom = box
    return max(contrast_ratio(image.pixelColor(x, y), ground)
               for y in range(top, bottom) for x in range(left, right, 2))


def test_the_same_forest_exported_under_two_themes_gets_two_different_grounds(
        forest, theme_choice, tmp_path):
    """A tree exported for a report used to come out dark whatever the reader chose.

    The ground was bound at import, so a light-theme workspace saved a black PNG
    onto a white report page. Both pictures are checked for lettering that can be
    read against their own ground, because a light picture with pale text would
    only be a different way of being unreadable.
    """
    theme_choice("dark")
    view = forest()
    view.save_image(tmp_path / "dark.png")
    theme_choice("light")
    # The picture is saved straight after a repaint, which is the moment a window
    # is in when the reader picks a theme and then exports: the canvas has already
    # been painted light while the pens on the items have not caught up yet.
    view.viewport().grab()
    view.save_image(tmp_path / "light.png")
    dark, light = QImage(str(tmp_path / "dark.png")), QImage(str(tmp_path / "light.png"))
    assert dark.pixelColor(0, 0) == QColor(theme.DARK["ground"])
    assert light.pixelColor(0, 0) == QColor(theme.LIGHT["ground"])
    assert dark.pixelColor(0, 0) != light.pixelColor(0, 0)
    # The title band of each picture, where the caption is drawn.
    assert _strongest_contrast(dark, theme.DARK["ground"], (35, 20, 700, 45)) >= 4.5
    assert _strongest_contrast(light, theme.LIGHT["ground"], (35, 20, 700, 45)) >= 4.5


def test_two_forests_saved_side_by_side_are_drawn_on_the_chosen_ground(
        forest, theme_choice, tmp_path):
    """A pair of trees is one picture, and it went dark on a light page too."""
    theme_choice("light")
    left, right = forest(), forest()
    image = render_side_by_side(left, right, left_title="Left", right_title="Right",
                                headline="Two layouts of the same cohort",
                                width=900, height=700)
    assert image.pixelColor(0, 0) == QColor(theme.LIGHT["ground"])
    assert _strongest_contrast(image, theme.LIGHT["ground"], (35, 16, 700, 40)) >= 4.5


def test_a_forest_already_on_screen_takes_a_theme_chosen_after_it_was_drawn(
        forest, theme_choice):
    """The canvas stayed black behind a pale window until the tree was redrawn.

    Nothing may be recalculated by the change: the same nodes, the same edges and
    the same numbers, on another ground.
    """
    theme_choice("dark")
    view = forest()
    assert view.backgroundBrush().color() == QColor(theme.DARK["ground"])
    before = {key: (node.pos().x(), node.pos().y()) for key, node in view.nodes.items()}
    theme_choice("light")
    # A repaint alone is enough for the canvas itself.
    assert view.viewport().grab().toImage().pixelColor(0, 0) == QColor(theme.LIGHT["ground"])
    view.refresh_theme()
    assert view.backgroundBrush().color() == QColor(theme.LIGHT["ground"])
    assert view.labels["a"].brush().color() == QColor(theme.LIGHT["ink"])
    assert view.edges[0][3].brush().color() == QColor(graph_colors()["edge_ink"])
    assert {key: (node.pos().x(), node.pos().y()) for key, node in view.nodes.items()} == before
    assert [edge["distance"] for edge in view.graph_contents()["edges"]] == [0, 1]


@pytest.mark.parametrize("name", ["dark", "slate", "light"])
def test_edge_numbers_group_outlines_and_the_selection_ring_read_on_every_theme(
        theme_choice, name):
    """Contrast on a tree is correctness: an unreadable edge number misreports a distance.

    Distance labels and the plain text around them are held to the readable-text
    ratio; rings and outlines, which are shapes rather than lettering, to the
    lower ratio a shape needs, and a group outline to being visible at all.
    """
    theme_choice(name)
    colors = graph_colors()
    ground = colors["ground"]
    for key in ("edge_ink", "ink", "muted"):
        assert contrast_ratio(colors[key], ground) >= 4.5, key
    for key in ("edge", "selected", "found", "singleton"):
        assert contrast_ratio(colors[key], ground) >= 3, key
    fill, border = halo_colors(PALETTE[0])
    assert contrast_ratio(_over(border, ground), ground) >= 1.5
    assert contrast_ratio(_over(fill, ground), ground) > 1


@pytest.mark.parametrize("name", ["dark", "slate", "light"])
def test_no_category_colour_disappears_into_the_ground_of_any_theme(theme_choice, name):
    """A node drawn in a pale category colour must still be a node on paper.

    Either the colour itself stands out, or the outline every node carries does.
    """
    theme_choice(name)
    colors = graph_colors()
    for color in PALETTE:
        assert (contrast_ratio(color, colors["ground"]) >= 3
                or contrast_ratio(colors["node_edge"], colors["ground"]) >= 3), color


@pytest.mark.parametrize("name", ["dark", "slate", "light"])
def test_the_count_on_a_merged_node_is_printed_in_ink_the_node_can_carry(theme_choice, name):
    """The count sits on the node, so it is read against the node, not the canvas.

    Printing it in the ground colour is right on a dark theme and invisible on a
    light one, where a pale mint node would carry a near-white numeral.
    """
    theme_choice(name)
    for color in PALETTE:
        assert contrast_ratio(readable_ink(color), color) >= 4.5, color
    theme_choice("dark")
    assert readable_ink(PALETTE[0]) == theme.DARK["ground"]
    theme_choice("light")
    assert readable_ink(PALETTE[0]) == theme.LIGHT["ink"]


def test_the_dark_theme_still_draws_the_forest_it_always_drew(theme_choice):
    """The colours the dark theme renders are the ones every saved picture has had.

    A theme selector must not quietly restyle the pictures already filed with an
    investigation, so the dark entries are pinned here, colour for colour.
    """
    theme_choice("dark")
    colors = graph_colors()
    assert (colors["ground"], colors["ink"], colors["muted"]) == (BACKGROUND, INK, MUTED)
    assert colors["edge"] == "#5F819B"
    assert colors["edge_ink"] == "#C7D8E8"
    assert colors["guide"] == "#36516C"
    assert colors["selected"] == "#E9F7FF"
    assert colors["found"] == "#FFD58A"
    assert colors["singleton"] == "#73869A"
    # The node outline was the ground itself, which is what separates two nodes
    # that overlap without drawing a ring on the canvas around each one.
    assert colors["node_edge"] == colors["ground"]
    fill, border = halo_colors("#48DCC0")
    assert (fill.alpha(), border.alpha()) == (15, 65)
    assert (fill.rgb() & 0xFFFFFF) == (QColor("#48DCC0").rgb() & 0xFFFFFF)
