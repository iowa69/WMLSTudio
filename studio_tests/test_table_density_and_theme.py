"""Dense tables and a chosen theme: the look of the workspace, and what it must never cost."""

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from wmlstudio import display, theme, ui_common
from wmlstudio.interface_settings import scaled_style_fragment

# Every colour the pre-token sheet wrote inline. If the token dict ever loses one
# of these, the dark theme quietly changes colour and nobody notices.
ORIGINAL_DARK_COLOURS = {
    "#E7EFF9", "#0B1220", "#0F1929", "#25354C", "#677C94", "#8FA9C2", "#A0B1C5",
    "#F2F7FF", "#183C3E", "#83E9D3", "#285755", "#121E30", "#293A53", "#132B35",
    "#28505A", "#122333", "#467786", "#16303D", "#48DCC0", "#302C24", "#6B593B",
    "#1B2A40", "#354A65", "#263B55", "#527392", "#30516A", "#244B59", "#48BFAE",
    "#C6FFF2", "#141F30", "#26364B", "#69E1CC", "#082B29", "#79EDD6", "#2AB69E",
    "#829B9F", "#254943", "#30554F", "#A3B4CA", "#1A2B42", "#183C43", "#28555A",
    "#85EDD7", "#101B2C", "#28526A", "#FFFFFF", "#55718F", "#131D2B", "#17253A",
    "#45617E", "#2B5067", "#203149", "#35526E", "#17243A", "#F5FBFF", "#2B3C54",
    "#25374E", "#24354A", "#21354D", "#6DD7C6", "#16253A", "#19293E", "#ABC0D6",
    "#293C55", "#101A2A", "#3A506D", "#577895", "#26384C", "#CCD9E8", "#28445D",
    "#F0F7FF", "#22324A", "#52708E", "#64809C", "#203C49", "#204C4C", "#36495F",
    "#152132", "#BDD0E2", "#152236", "#20374A", "#86EDD7", "#24394F", "#142739",
    "#C2D3E6",
}


@pytest.fixture(autouse=True)
def restored_appearance():
    """Put the shipped theme and density back: the sheet lives on the QApplication."""
    before = (theme.active_theme(), theme.active_density(), ui_common.table_density())
    yield
    theme.set_active_theme(before[0], before[1])
    ui_common.set_table_density(before[2])
    application = QApplication.instance()
    if application is not None:
        application.setStyleSheet(theme.active_style())


@pytest.fixture
def table(qtbot):
    widget = ui_common.make_table(["Isolate", "Organism evidence", "AMR evidence state"])
    qtbot.addWidget(widget)
    return widget


# ---------------------------------------------------------------------------
# The theme, now built from tokens rather than written out by hand.
# ---------------------------------------------------------------------------


def test_the_dark_theme_still_writes_every_colour_the_frozen_sheet_wrote():
    """The tokens were extracted from the old sheet; a dropped key would silently recolour it."""
    sheet = theme.build_style(theme.DARK, "spacious")
    for colour in ORIGINAL_DARK_COLOURS:
        assert colour in sheet, f"{colour} disappeared from the dark theme"
    # And the spacious density is the old spacing, declaration for declaration.
    assert "QTableView::item, QTableWidget::item { border-bottom: 1px solid #24354A; padding: 8px; }" in sheet
    assert "padding: 11px 8px; border: none; border-right: 1px solid #293C55;" in sheet
    assert "QPushButton, QToolButton { background: #1B2A40; border: 1px solid #354A65; " \
           "border-radius: 7px; padding: 9px 14px; font-weight: 600; }" in sheet


def test_every_theme_fills_in_every_token_the_sheet_and_the_palette_ask_for():
    """A palette missing one key would raise mid-render, or leave a literal $token on screen."""
    keys = set(theme.DARK)
    for name, tokens in theme.THEMES.items():
        assert set(tokens) == keys, f"{name} does not define the same colours as dark"
        for density in theme.TABLE_DENSITY:
            sheet = theme.build_style(tokens, density)
            assert "$" not in sheet, f"{name}/{density} left a token unsubstituted"
    assert set(theme.THEMES) == set(theme.THEME_LABELS) == set(theme.THEME_NOTES)
    assert set(theme.TABLE_DENSITY) == set(theme.DENSITY_LABELS)


def test_changing_the_theme_changes_the_sheet_so_the_window_actually_applies_it():
    """apply_application_style skips an unchanged string; an equal sheet would do nothing."""
    sheets = {name: theme.style_for(name, "compact") for name in theme.THEMES}
    assert len(set(sheets.values())) == len(sheets)
    assert len({theme.style_for("slate", density) for density in theme.TABLE_DENSITY}) == 3
    theme.set_active_theme("light")
    assert theme.active_style() == sheets["light"]
    assert theme.active_style() != sheets["dark"]


def test_an_unknown_theme_or_density_name_falls_back_instead_of_raising():
    """The name comes from an INI file a person can edit; it must never crash the window."""
    assert theme.theme_tokens("chartreuse") is theme.THEMES[theme.DEFAULT_THEME]
    assert theme.density_metrics("enormous") is theme.TABLE_DENSITY[theme.DEFAULT_DENSITY]
    theme.set_active_theme("chartreuse", "enormous")
    assert theme.active_theme() == theme.DEFAULT_THEME
    assert theme.active_density() == theme.DEFAULT_DENSITY


def test_the_light_theme_is_light_and_the_dark_themes_are_dark():
    """A theme named light that renders dark would be a lie in the one place it shows."""
    from PySide6.QtGui import QColor

    def lightness(hex_value):
        return QColor(hex_value).lightnessF()

    assert lightness(theme.LIGHT["ground"]) > 0.85
    assert lightness(theme.LIGHT["ink"]) < 0.25
    for name in ("dark", "slate"):
        assert lightness(theme.THEMES[name]["ground"]) < 0.25
        assert lightness(theme.THEMES[name]["ink"]) > 0.85
    # Slate is the lighter dark: its panels stand further off its ground than the
    # near-black theme's do, which is the whole reason it exists.
    assert (lightness(theme.SLATE["surface"]) - lightness(theme.SLATE["ground"])) > \
           (lightness(theme.DARK["surface"]) - lightness(theme.DARK["ground"]))
    assert lightness(theme.SLATE["ground"]) > lightness(theme.DARK["ground"])


def test_the_light_theme_points_at_an_arrow_that_can_be_seen_on_white():
    """The bundled chevron is pale: reused on a white field it leaves combo boxes blank."""
    from pathlib import Path

    from wmlstudio import theme as module
    icons = Path(module.__file__).resolve().parent / "resources/ui"
    assert theme.LIGHT["chevron_down"] != theme.DARK["chevron_down"]
    for tokens in theme.THEMES.values():
        for key in ("chevron_down", "chevron_up"):
            assert (icons / tokens[key]).exists(), tokens[key]


def test_the_text_size_rescaler_still_finds_every_font_size_in_every_theme():
    """scaled_style_fragment rewrites `font-size: Npx` by regex; a new spelling kills scaling."""
    for name in theme.THEMES:
        for density in theme.TABLE_DENSITY:
            sheet = theme.style_for(name, density)
            assert "font-size: 13px" in sheet
            bigger = scaled_style_fragment(sheet, 150)
            assert "font-size: 19.5px" in bigger
            assert "font-size: 13px" not in bigger


def test_the_qt_palette_follows_the_theme_it_is_given(qtbot):
    """Qt draws its own file and message dialogs from the palette, not from the sheet."""
    from PySide6.QtGui import QColor, QPalette

    application = QApplication.instance()
    before = application.palette()
    try:
        theme.apply_palette(application, theme.LIGHT)
        palette = application.palette()
        assert palette.color(QPalette.ColorRole.Window) == QColor(theme.LIGHT["ground"])
        assert palette.color(QPalette.ColorRole.Text) == QColor(theme.LIGHT["ink"])
        theme.apply_dark_palette(application)
        assert application.palette().color(QPalette.ColorRole.Window) == \
            QColor(theme.DARK["ground"])
    finally:
        application.setPalette(before)


# ---------------------------------------------------------------------------
# Tables: the reported complaint was that a screen holds four rows of evidence.
# ---------------------------------------------------------------------------


def test_a_new_table_is_compact_and_hides_the_row_number_strip(table):
    """39px rows and a row-number column were most of why a screenful held so little."""
    assert ui_common.table_density() == "compact"
    assert theme.DEFAULT_DENSITY == "compact"
    assert table.verticalHeader().defaultSectionSize() == 24
    assert table.verticalHeader().isHidden() is True
    assert table.showGrid() is True, "column rules replace the loose padding"
    assert table.property("tableDensity") == "compact"


def test_the_compact_sheet_pads_a_cell_far_less_than_the_sheet_that_was_reported():
    """8px cell and 11px header padding is what put ~330px of chrome above the first row."""
    compact = theme.style_for("slate", "compact")
    assert "padding: 3px 7px; }" in compact
    assert "padding: 5px 7px;" in compact
    # And the roomiest setting still offers the old spacing to anyone who wants it.
    assert "padding: 8px; }" in theme.style_for("slate", "spacious")


def test_a_column_is_never_opened_too_narrow_to_read_its_own_heading(table):
    """A heading elided to "R evidence st" leaves the reader unsure what the column holds."""
    header = table.horizontalHeader()
    for column in range(table.columnCount() - 1):
        assert header.sectionSize(column) >= ui_common.header_width(table, column)
    # The longest heading here is wider than the default column, so it had to grow.
    assert header.sectionSize(2) > theme.density_metrics("compact")["column_width"] or \
        ui_common.header_width(table, 2) <= theme.density_metrics("compact")["column_width"]


def test_a_very_long_heading_widens_its_column_only_up_to_a_limit(qtbot):
    """One enormous column would push every other one off the screen: the same defect again."""
    # The last section stretches to the viewport, so the long heading is measured
    # in the middle where the width this code chose is the width Qt keeps.
    widget = ui_common.make_table(["Isolate", "x" * 400, "Reviewed"])
    qtbot.addWidget(widget)
    assert ui_common.header_width(widget, 1) == ui_common.MAXIMUM_COLUMN
    assert widget.horizontalHeader().sectionSize(1) == ui_common.MAXIMUM_COLUMN


def test_the_columns_widen_to_what_the_rows_actually_hold_after_a_fill(table):
    """Fitting only the heading still hides a long value; fitting must follow the contents."""
    table.insertRow(0)
    for column, value in enumerate(("KP-014", "Klebsiella pneumoniae, provisional", "12")):
        table.setItem(0, column, ui_common.cell(value))
    narrow = table.horizontalHeader().sectionSize(1)
    ui_common.fit_table_columns(table)
    assert table.horizontalHeader().sectionSize(1) >= narrow
    assert table.horizontalHeader().sectionSize(1) >= ui_common.header_width(table, 1)


def test_a_column_a_caller_made_wider_keeps_its_width_through_a_fill(table):
    """A page that widened a column knows something about it that a refit does not."""
    table.horizontalHeader().resizeSection(1, 280)
    table.insertRow(0)
    table.setItem(0, 1, ui_common.cell("Klebsiella"))
    ui_common.fit_table_columns(table)
    assert table.horizontalHeader().sectionSize(1) == 280


def test_a_heading_wins_over_a_column_set_too_narrow_to_show_it(table):
    """Nothing is gained by a tidy narrow column whose heading cannot be read."""
    table.horizontalHeader().resizeSection(0, 20)
    ui_common.fit_table_columns(table)
    assert table.horizontalHeader().sectionSize(0) == ui_common.header_width(table, 0)


def test_a_matrix_too_wide_to_measure_is_left_alone(qtbot):
    """A cgMLST-sized table would cost more to measure than tidy columns are worth."""
    widget = ui_common.make_table([f"t{index}" for index in range(ui_common.AUTOFIT_COLUMN_LIMIT + 5)])
    qtbot.addWidget(widget)
    widths = [widget.horizontalHeader().sectionSize(c) for c in range(widget.columnCount())]
    ui_common.fit_table_columns(widget)
    assert [widget.horizontalHeader().sectionSize(c)
            for c in range(widget.columnCount())] == widths


def test_choosing_a_density_changes_tables_that_are_already_on_screen(qtbot):
    """A setting that only affects tables built later reads as a setting that did nothing."""
    from PySide6.QtWidgets import QVBoxLayout, QWidget

    host = QWidget()
    # The table is registered only through its parent: pytest-qt closes what it
    # was handed, and a child that Qt already deleted cannot be closed again.
    table = ui_common.make_table(["Isolate", "Organism evidence"])
    QVBoxLayout(host).addWidget(table)
    qtbot.addWidget(host)
    assert ui_common.restyle_tables(host, "spacious") == 1
    assert table.verticalHeader().defaultSectionSize() == 39
    assert table.showGrid() is False
    ui_common.restyle_tables(host, "compact")
    assert table.verticalHeader().defaultSectionSize() == 24


def test_set_table_density_refuses_a_name_it_does_not_know():
    """The name arrives from a settings file, so an unknown one must not build an empty table."""
    assert ui_common.set_table_density("roomy") == "roomy"
    assert ui_common.set_table_density("enormous") == theme.DEFAULT_DENSITY
    assert ui_common.set_table_density(None) == theme.DEFAULT_DENSITY


def test_a_quantity_is_right_aligned_so_a_column_of_numbers_can_be_compared():
    """Ragged digits are read one value at a time; a distance column is read down."""
    right = int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    assert int(ui_common.cell(12).textAlignment()) == right
    assert int(ui_common.cell("2358").textAlignment()) == right
    assert int(ui_common.cell("99.5%").textAlignment()) == right
    assert int(ui_common.cell("ST 258").textAlignment()) != right
    assert int(ui_common.cell("Klebsiella").textAlignment()) != right
    # Missing evidence stays where text sits: it is not a quantity and must not
    # line up in a column of numbers as though it were a zero.
    assert ui_common.cell(None).text() == "—"
    assert int(ui_common.cell(None).textAlignment()) != right
    assert ui_common.cell(None).toolTip() == "Not available"


def test_a_cell_can_still_be_aligned_by_the_caller(table):
    """Alignment is presentation: a caller that knows its column must be able to say so."""
    item = ui_common.cell(42, "KP-014", align=Qt.AlignmentFlag.AlignLeft)
    assert int(item.textAlignment()) == int(Qt.AlignmentFlag.AlignLeft)
    assert item.data(Qt.ItemDataRole.UserRole) == "KP-014"


# ---------------------------------------------------------------------------
# The saved choice. It lives beside the display keys, not with the project.
# ---------------------------------------------------------------------------


def test_the_chosen_theme_and_density_survive_a_restart(tmp_path):
    """A look the user had to hunt for twice is a look they will not choose again."""
    assert display.read_appearance(tmp_path) == {"theme": "", "density": ""}
    stored = display.write_appearance(theme="light", density="roomy", root=tmp_path)
    assert stored == {"theme": "light", "density": "roomy"}
    assert display.read_appearance(tmp_path) == stored
    # And they do not disturb the display keys they are filed beside.
    assert display.read_display_settings(tmp_path) == {
        "mode": "system", "scale_percent": 100, "rounding": "exact"}


def test_a_hand_edited_appearance_name_is_ignored_rather_than_applied(tmp_path):
    """Interface.ini is a plain text file on a shared machine; rubbish in it must not stick."""
    preferences = display.display_preferences(tmp_path)
    preferences.setValue("display/theme", "url(evil); background: red")
    preferences.setValue("display/table_density", "")
    preferences.sync()
    assert display.read_appearance(tmp_path) == {"theme": "", "density": ""}
    with pytest.raises(ValueError):
        display.write_appearance(theme="not a theme name", root=tmp_path)


def test_forgetting_a_choice_returns_to_the_shipped_look(tmp_path):
    display.write_appearance(theme="light", density="spacious", root=tmp_path)
    assert display.write_appearance(theme="", density="", root=tmp_path) == {
        "theme": "", "density": ""}
    assert display.read_appearance(tmp_path) == {"theme": "", "density": ""}


def test_the_appearance_page_says_what_a_theme_does_not_change():
    """A viewing preference must never read as something that touched the evidence."""
    assert "not what it found" in display.APPEARANCE_NOTICE
    assert "threshold" in display.APPEARANCE_NOTICE


# ---------------------------------------------------------------------------
# The Settings panel, on a real window.
# ---------------------------------------------------------------------------


@pytest.fixture
def window(qtbot, tmp_path):
    from wmlstudio.app import MainWindow
    widget = MainWindow(storage_root=tmp_path / "workspace")
    qtbot.addWidget(widget)
    yield widget
    widget.close()


def test_the_settings_page_offers_the_theme_and_the_table_density(window, tmp_path):
    """Both were asked for by name; a control nobody can find is a control that is absent."""
    from wmlstudio.interface_settings import InterfaceSettingsPanel
    panel = InterfaceSettingsPanel(window, display_root=tmp_path / "dataroot")
    offered = [panel.theme.itemData(index) for index in range(panel.theme.count())]
    assert offered == list(theme.THEMES)
    densities = [panel.density.itemData(index) for index in range(panel.density.count())]
    assert densities == list(theme.TABLE_DENSITY)
    assert panel.density.currentData() == ui_common.table_density()


def test_choosing_a_theme_saves_it_and_repaints_the_running_window(window, tmp_path):
    """The window skips a stylesheet that would not change; a theme that does nothing is a bug."""
    from wmlstudio.interface_settings import InterfaceSettingsPanel, restore_appearance
    panel = InterfaceSettingsPanel(window, display_root=tmp_path / "dataroot")
    before = QApplication.instance().styleSheet()
    panel.theme.setCurrentIndex(panel.theme.findData("light"))
    panel.density.setCurrentIndex(panel.density.findData("spacious"))
    assert display.read_appearance(tmp_path / "dataroot") == {
        "theme": "light", "density": "spacious"}
    assert theme.active_theme() == "light"
    assert ui_common.table_density() == "spacious"
    assert QApplication.instance().styleSheet() != before
    assert theme.LIGHT["ground"] in QApplication.instance().styleSheet()
    assert "spacious" in panel.preview_note.text() or "39 px rows" in panel.preview_note.text()
    # A restart reads the same two names back and puts them in force again.
    assert restore_appearance(tmp_path / "dataroot") == {
        "theme": "light", "density": "spacious"}
