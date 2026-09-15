"""Opaque, high-contrast native desktop surfaces and a matching Qt palette.

The sheet used to be one frozen string with every colour written inline, which
made a second theme impossible to add without rewriting it. It is now a template
of named tokens plus one dict of colours per theme, so a theme is a palette and
nothing else: no selector is duplicated, and a colour that is wrong is wrong in
exactly one place.

Density is the other axis. The reported complaint was that a table showed four
rows where a spreadsheet shows twenty, and the measured cause was padding, not
font size, so the paddings that decide how much chrome sits above the first row
live in `TABLE_DENSITY` and are chosen with the theme. Text size stays separate:
`interface_settings.scaled_style_fragment` rewrites every `font-size: Npx` by
regular expression, so that exact spelling must survive anything done here.
"""

from pathlib import Path
from string import Template

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette

_ICONS = (Path(__file__).resolve().parent / "resources/ui").as_posix()

# ---------------------------------------------------------------------------
# Palettes. Every theme defines every key: a missing one raises at build time
# rather than quietly falling back to a colour from a different theme.
# ---------------------------------------------------------------------------

DARK = {
    # Ground, panels and text. This is the original sheet's palette, unchanged,
    # so the dark theme still renders what it rendered before, colour for colour.
    "ink": "#E7EFF9", "ground": "#0B1220", "sidebar": "#0F1929", "rule": "#25354C",
    "ink_off": "#677C94", "eyebrow_ink": "#8FA9C2", "muted": "#A0B1C5",
    "metric_ink": "#F2F7FF", "surface": "#121E30", "card_border": "#293A53",
    "purpose_ink": "#C2D3E6", "group_title_ink": "#BDD0E2",
    # Accents.
    "accent": "#48DCC0", "accent_hover": "#79EDD6", "accent_press": "#2AB69E",
    "accent_dim": "#48BFAE", "on_accent": "#082B29", "focus_ring": "#69E1CC",
    # Badges, hero, drop target and the amber note.
    "badge_bg": "#183C3E", "badge_ink": "#83E9D3", "badge_border": "#285755",
    "hero_bg": "#132B35", "hero_border": "#28505A", "drop_bg": "#122333",
    "drop_border": "#467786", "drop_hover_bg": "#16303D",
    "note_bg": "#302C24", "note_border": "#6B593B",
    # Buttons.
    "button_bg": "#1B2A40", "button_hover_bg": "#263B55", "button_hover_border": "#527392",
    "button_press_bg": "#30516A", "button_on_bg": "#244B59", "button_on_ink": "#C6FFF2",
    "button_off_bg": "#141F30", "primary_off_ink": "#829B9F", "primary_off_bg": "#254943",
    "primary_off_border": "#30554F",
    # Sidebar navigation.
    "nav_ink": "#A3B4CA", "nav_hover_bg": "#1A2B42", "nav_on_bg": "#183C43",
    "nav_on_border": "#28555A", "nav_on_ink": "#85EDD7",
    # Fields, menus and steppers.
    "control_border": "#354A65", "control_border_off": "#26364B", "field_bg": "#101B2C",
    "field_hover_border": "#55718F", "field_off_bg": "#131D2B", "menu_bg": "#17253A",
    "menu_border": "#45617E", "menu_selection": "#2B5067", "stepper_bg": "#203149",
    "stepper_hover_bg": "#35526E", "selection": "#28526A", "on_selection": "#FFFFFF",
    # Tables and the rest of the item views.
    "row_alt": "#17243A", "on_row_selected": "#F5FBFF", "view_border": "#2B3C54",
    "gridline": "#25374E", "row_rule": "#24354A", "row_hover": "#21354D",
    "cell_focus": "#6DD7C6", "header_bg": "#16253A", "header_section_bg": "#19293E",
    "header_ink": "#ABC0D6", "header_rule": "#293C55",
    # Scrollbars, progress, menu bar, tooltips, checkboxes, tabs.
    "scrollbar_bg": "#101A2A", "scrollbar_handle": "#3A506D",
    "scrollbar_handle_hover": "#577895", "progress_track": "#26384C",
    "menubar_ink": "#CCD9E8", "menubar_selection": "#28445D", "tooltip_ink": "#F0F7FF",
    "tooltip_bg": "#22324A", "tooltip_border": "#52708E", "check_border": "#64809C",
    "check_hover_bg": "#203C49", "check_on_bg": "#204C4C", "check_off_border": "#36495F",
    "check_off_bg": "#152132", "tab_bg": "#152236", "tab_on_bg": "#20374A",
    "tab_on_ink": "#86EDD7", "tab_hover_bg": "#24394F", "workspace_tab_on_bg": "#142739",
    # Arrow artwork. The bundled chevrons are drawn pale, so a light theme has to
    # point at the ink-coloured copies or its combo boxes lose their arrows.
    "chevron_down": "chevron-down.svg", "chevron_up": "chevron-up.svg",
    # QPalette-only roles: Fusion draws Qt's own dialogs from these, and nothing
    # in the sheet can reach them.
    "bright_text": "#FF91A6", "link": "#80C8FF", "link_visited": "#C29AFF",
    "placeholder": "#7D91AA", "shade_light": "#3B4E68", "shade_midlight": "#2A3D57",
    "shade_mid": "#25364F", "shade_dark": "#070D17", "shade_shadow": "#03070D",
}

SLATE = {
    # The same shapes on a raised blue-grey. The near-black ground of DARK gives
    # a panel almost nothing to sit on; here the panel is several steps lighter
    # than the ground, which is what makes a dense table read as a sheet of rows.
    "ink": "#EDF2F8", "ground": "#1A2230", "sidebar": "#141B27", "rule": "#2C3849",
    "ink_off": "#7A8AA0", "eyebrow_ink": "#9AAEC4", "muted": "#AEBDCE",
    "metric_ink": "#F6FAFF", "surface": "#28344A", "card_border": "#3B4A63",
    "purpose_ink": "#C8D7E8", "group_title_ink": "#C6D6E6",
    "accent": "#45D6BC", "accent_hover": "#74E7D2", "accent_press": "#2BB29B",
    "accent_dim": "#43BFAB", "on_accent": "#05302B", "focus_ring": "#66DFC9",
    "badge_bg": "#1E4744", "badge_ink": "#8DEEDA", "badge_border": "#2F6460",
    "hero_bg": "#1D3A41", "hero_border": "#2F5D67", "drop_bg": "#1E3040",
    "drop_border": "#4E808E", "drop_hover_bg": "#22404C",
    "note_bg": "#3A3428", "note_border": "#776444",
    "button_bg": "#2B3749", "button_hover_bg": "#36455B", "button_hover_border": "#5E7C99",
    "button_press_bg": "#3E5C75", "button_on_bg": "#2C5866", "button_on_ink": "#CFFFF4",
    "button_off_bg": "#202A38", "primary_off_ink": "#86A09F", "primary_off_bg": "#2B5450",
    "primary_off_border": "#3A605B",
    "nav_ink": "#B2C1D3", "nav_hover_bg": "#27334A", "nav_on_bg": "#1E474C",
    "nav_on_border": "#2F6165", "nav_on_ink": "#8AEDD9",
    "control_border": "#435268", "control_border_off": "#323F52", "field_bg": "#161D29",
    "field_hover_border": "#607E9B", "field_off_bg": "#1C2431", "menu_bg": "#253044",
    "menu_border": "#4C6684", "menu_selection": "#325A73", "stepper_bg": "#2C3A50",
    "stepper_hover_bg": "#3E5A75", "selection": "#2F5A76", "on_selection": "#FFFFFF",
    "row_alt": "#2E3B52", "on_row_selected": "#F4FAFF", "view_border": "#3E4E68",
    "gridline": "#37455B", "row_rule": "#344157", "row_hover": "#323F57",
    "cell_focus": "#6BDCCB", "header_bg": "#222D3D", "header_section_bg": "#2A3648",
    "header_ink": "#C3D0DF", "header_rule": "#3B4A63",
    "scrollbar_bg": "#1B2331", "scrollbar_handle": "#445573",
    "scrollbar_handle_hover": "#5E7A99", "progress_track": "#303E53",
    "menubar_ink": "#D3DEEB", "menubar_selection": "#31506B", "tooltip_ink": "#F2F8FF",
    "tooltip_bg": "#2A3A52", "tooltip_border": "#5A7896", "check_border": "#6D89A5",
    "check_hover_bg": "#274450", "check_on_bg": "#26544F", "check_off_border": "#3D5066",
    "check_off_bg": "#1E2733", "tab_bg": "#202B3C", "tab_on_bg": "#2A3F52",
    "tab_on_ink": "#8AEDD9", "tab_hover_bg": "#2E3F54", "workspace_tab_on_bg": "#1F2F42",
    "chevron_down": "chevron-down.svg", "chevron_up": "chevron-up.svg",
    "bright_text": "#FF9AAE", "link": "#8ACEFF", "link_visited": "#C9A6FF",
    "placeholder": "#8496AD", "shade_light": "#4B5E78", "shade_midlight": "#3A4D67",
    "shade_mid": "#32445D", "shade_dark": "#101722", "shade_shadow": "#0A0F17",
}

LIGHT = {
    # Paper panels on a cool grey ground, for a bright room and for a screenshot
    # that has to survive being printed into a report. The teal is darkened until
    # it can carry text on white; the pale accents of the dark themes cannot.
    "ink": "#16202C", "ground": "#EEF2F7", "sidebar": "#FFFFFF", "rule": "#D5DEE8",
    "ink_off": "#97A5B5", "eyebrow_ink": "#5E7085", "muted": "#566878",
    "metric_ink": "#0E1822", "surface": "#FFFFFF", "card_border": "#DCE4EE",
    "purpose_ink": "#33465A", "group_title_ink": "#3D4F62",
    "accent": "#0E8C7A", "accent_hover": "#12A38E", "accent_press": "#0A6E60",
    "accent_dim": "#17A28D", "on_accent": "#FFFFFF", "focus_ring": "#0E8C7A",
    "badge_bg": "#DFF6F0", "badge_ink": "#0B6B5C", "badge_border": "#9FDCD0",
    "hero_bg": "#E4F4F1", "hero_border": "#A9D9D0", "drop_bg": "#F1F6FB",
    "drop_border": "#8FB4C4", "drop_hover_bg": "#E3F3F0",
    "note_bg": "#FDF4E0", "note_border": "#E0C387",
    "button_bg": "#FFFFFF", "button_hover_bg": "#EDF3F9", "button_hover_border": "#9FB0C4",
    "button_press_bg": "#DCE6F0", "button_on_bg": "#DCF2EC", "button_on_ink": "#0A5B4F",
    "button_off_bg": "#F2F5F9", "primary_off_ink": "#8FA8A3", "primary_off_bg": "#CFE3DE",
    "primary_off_border": "#BAD5CF",
    "nav_ink": "#46586B", "nav_hover_bg": "#E9F0F7", "nav_on_bg": "#DCF2EC",
    "nav_on_border": "#A6D9CD", "nav_on_ink": "#0A5B4F",
    "control_border": "#C6D0DD", "control_border_off": "#DFE5EC", "field_bg": "#FFFFFF",
    "field_hover_border": "#9FB0C4", "field_off_bg": "#F2F5F9", "menu_bg": "#FFFFFF",
    "menu_border": "#C6D0DD", "menu_selection": "#D8EAF7", "stepper_bg": "#F0F4F9",
    "stepper_hover_bg": "#E0E8F1", "selection": "#CBE4F5", "on_selection": "#0B1A26",
    "row_alt": "#F5F8FC", "on_row_selected": "#0B1A26", "view_border": "#D5DEE8",
    "gridline": "#E1E8F0", "row_rule": "#E7EDF4", "row_hover": "#EDF4FB",
    "cell_focus": "#0E8C7A", "header_bg": "#EAF0F6", "header_section_bg": "#EDF2F8",
    "header_ink": "#3D4F62", "header_rule": "#D5DEE8",
    "scrollbar_bg": "#EAEFF5", "scrollbar_handle": "#BCC8D6",
    "scrollbar_handle_hover": "#9BAABC", "progress_track": "#DDE5EE",
    "menubar_ink": "#33465A", "menubar_selection": "#DCE8F4", "tooltip_ink": "#0E1822",
    "tooltip_bg": "#FFFFFF", "tooltip_border": "#B6C4D4", "check_border": "#A6B4C4",
    # The bundled tick is drawn in pale mint, so a ticked box has to stay dark
    # enough to carry it: a pale ticked box would read as an empty one.
    "check_hover_bg": "#E0F2ED", "check_on_bg": "#0E8C7A", "check_off_border": "#DCE2EA",
    "check_off_bg": "#F2F5F9", "tab_bg": "#E7EDF4", "tab_on_bg": "#FFFFFF",
    "tab_on_ink": "#0A5B4F", "tab_hover_bg": "#EDF3F9", "workspace_tab_on_bg": "#FFFFFF",
    "chevron_down": "chevron-down-ink.svg", "chevron_up": "chevron-up-ink.svg",
    "bright_text": "#B3253F", "link": "#0B63B8", "link_visited": "#6B3FB0",
    "placeholder": "#97A5B5", "shade_light": "#FFFFFF", "shade_midlight": "#F0F4F9",
    "shade_mid": "#C6D0DD", "shade_dark": "#8794A5", "shade_shadow": "#6A7686",
}

THEMES = {"dark": DARK, "slate": SLATE, "light": LIGHT}
THEME_LABELS = {
    "dark": "Midnight — deepest dark",
    "slate": "Slate — lighter dark (recommended)",
    "light": "Daylight — light",
}
THEME_NOTES = {
    "dark": "Near-black ground. Best in a dark room.",
    "slate": "Raised blue-grey panels, stronger separation between page and table.",
    "light": "Paper panels for a bright room, or for a screenshot in a report.",
}
DEFAULT_THEME = "slate"

# The colours the graph, report and widget code imports directly. They stay bound
# to DARK: those modules read them once at import, so they cannot follow a theme
# chosen later, and a half-applied palette would be worse than a fixed one.
BACKGROUND = DARK["ground"]
SURFACE = DARK["surface"]
BORDER = DARK["view_border"]
INK = DARK["ink"]
TEAL = DARK["accent"]
MUTED = DARK["muted"]
PALETTE = ["#48DCC0", "#73ACFF", "#C29AFF", "#FFB976", "#FF86A4",
           "#B5D976", "#79D6EF", "#DBB0E8", "#E7D577", "#A3B5D2"]

# ---------------------------------------------------------------------------
# Density. Every number here was measured off a rendered screenshot: 39 px rows,
# 8 px cell padding and 11 px header padding put roughly 330 px of chrome above
# the first row of evidence, which is why the tables read as a brochure.
# ---------------------------------------------------------------------------

TABLE_DENSITY = {
    "compact": {
        "row_height": 24, "column_width": 130, "grid": True,
        "cell_padding": "3px 7px", "header_padding": "5px 7px",
        "button_padding": "6px 12px", "nav_padding": "8px 14px",
        "tab_padding": "6px 12px", "workspace_tab_padding": "6px 11px",
        "field_padding": "5px 8px", "strip_margins": (10, 2, 10, 2),
    },
    "roomy": {
        "row_height": 30, "column_width": 140, "grid": True,
        "cell_padding": "5px 8px", "header_padding": "7px 8px",
        "button_padding": "8px 13px", "nav_padding": "10px 15px",
        "tab_padding": "8px 14px", "workspace_tab_padding": "7px 12px",
        "field_padding": "7px", "strip_margins": (12, 3, 12, 3),
    },
    "spacious": {
        # What the application looked like before this work, kept so that the
        # change is reversible by anyone who preferred it.
        "row_height": 39, "column_width": 145, "grid": False,
        "cell_padding": "8px", "header_padding": "11px 8px",
        "button_padding": "9px 14px", "nav_padding": "12px 16px",
        "tab_padding": "9px 15px", "workspace_tab_padding": "8px 12px",
        "field_padding": "8px", "strip_margins": (12, 4, 12, 4),
    },
}

DENSITY_LABELS = {
    "compact": "Compact — most rows on screen (recommended)",
    "roomy": "Roomy — a little breathing room",
    "spacious": "Spacious — the widest spacing",
}
DEFAULT_DENSITY = "compact"


def theme_tokens(name=None):
    """The palette for a theme name, falling back to the default for anything else."""
    return THEMES.get(str(name or "").strip().lower(), THEMES[DEFAULT_THEME])


def density_metrics(name=None):
    """The measurements for a density name, falling back to the default."""
    return TABLE_DENSITY.get(str(name or "").strip().lower(), TABLE_DENSITY[DEFAULT_DENSITY])


def apply_palette(application, tokens=None):
    """Style native Fusion controls and Qt-owned dialogs, including disabled states."""
    tokens = tokens or DARK
    palette = QPalette()
    colors = {
        "Window": tokens["ground"], "WindowText": tokens["ink"], "Base": tokens["surface"],
        "AlternateBase": tokens["row_alt"], "ToolTipBase": tokens["tooltip_bg"],
        "ToolTipText": tokens["tooltip_ink"], "Text": tokens["ink"],
        "Button": tokens["button_bg"], "ButtonText": tokens["ink"],
        "BrightText": tokens["bright_text"], "Link": tokens["link"],
        "LinkVisited": tokens["link_visited"], "Highlight": tokens["selection"],
        "HighlightedText": tokens["on_row_selected"], "PlaceholderText": tokens["placeholder"],
        "Light": tokens["shade_light"], "Midlight": tokens["shade_midlight"],
        "Mid": tokens["shade_mid"], "Dark": tokens["shade_dark"],
        "Shadow": tokens["shade_shadow"],
    }
    for role, color in colors.items():
        palette.setColor(getattr(QPalette.ColorRole, role), QColor(color))
    for role in (QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText,
                 QPalette.ColorRole.WindowText, QPalette.ColorRole.PlaceholderText):
        palette.setColor(QPalette.ColorGroup.Disabled, role, QColor(tokens["ink_off"]))
    application.setPalette(palette)
    application.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, True)


def apply_dark_palette(application, tokens=None):
    """The original name, kept because the window and the graph tests call it.

    It stays on the dark palette when asked for nothing in particular. A caller
    that wants whatever theme the user chose should use `apply_theme`.
    """
    apply_palette(application, tokens if tokens is not None else DARK)


def apply_theme(application) -> str:
    """Put the active theme's QPalette on the application and return its sheet.

    The sheet is returned rather than set here so the caller keeps going through
    `ui_workbench.apply_application_style`, which skips an assignment that would
    change nothing: Qt re-polishes every widget of every window on each one.
    """
    apply_palette(application, theme_tokens(_active["theme"]))
    return active_style()


_TEMPLATE = """
QWidget { color: $ink; font-family: 'Segoe UI', 'Inter', 'DejaVu Sans'; font-size: 13px; }
QMainWindow, QWidget#main, QWidget#page, QStackedWidget, QDialog { background: $ground; }
QWidget#sidebar { background: $sidebar; border-right: 1px solid $rule; }
QLabel { background: transparent; }
QLabel:disabled { color: $ink_off; }
QLabel#brand { font-size: 22px; font-weight: 700; letter-spacing: -1px; }
QLabel#eyebrow { font-size: 10px; font-weight: 700; color: $eyebrow_ink; letter-spacing: 2px; }
QLabel#title { font-size: 28px; font-weight: 700; letter-spacing: -1px; }
QLabel#heroTitle { font-size: 33px; font-weight: 700; letter-spacing: -1px; }
QLabel#muted { color: $muted; }
QLabel#small { color: $muted; font-size: 11px; }
QLabel#metric { color: $metric_ink; font-size: 30px; font-weight: 600; }
QLabel#cardTitle { font-size: 15px; font-weight: 600; }
QLabel#badge { background: $badge_bg; color: $badge_ink; border: 1px solid $badge_border; border-radius: 12px; padding: 5px 10px; font-size: 11px; }
QFrame#card { background: $surface; border: 1px solid $card_border; border-radius: 12px; }
QFrame#hero { background: $hero_bg; border: 1px solid $hero_border; border-radius: 15px; }
QFrame#drop { background: $drop_bg; border: 1px dashed $drop_border; border-radius: 12px; }
QFrame#drop:hover { background: $drop_hover_bg; border-color: $accent; }
QFrame#note { background: $note_bg; border: 1px solid $note_border; border-radius: 9px; }
QPushButton, QToolButton { background: $button_bg; border: 1px solid $control_border; border-radius: 7px; padding: $button_padding; font-weight: 600; }
QPushButton:hover, QToolButton:hover { background: $button_hover_bg; border-color: $button_hover_border; }
QPushButton:pressed, QToolButton:pressed { background: $button_press_bg; }
QPushButton:checked, QToolButton:checked { background: $button_on_bg; border-color: $accent_dim; color: $button_on_ink; }
QPushButton:disabled, QToolButton:disabled { color: $ink_off; background: $button_off_bg; border-color: $control_border_off; }
QPushButton:focus, QToolButton:focus { border-color: $focus_ring; }
QPushButton#primary { color: $on_accent; background: $accent; border-color: $accent; }
QPushButton#primary:hover { background: $accent_hover; border-color: $accent_hover; }
QPushButton#primary:pressed { background: $accent_press; }
QPushButton#primary:disabled { color: $primary_off_ink; background: $primary_off_bg; border-color: $primary_off_border; }
QPushButton#nav { border: 1px solid transparent; background: transparent; text-align: left; padding: $nav_padding; color: $nav_ink; font-weight: 500; }
QPushButton#nav:hover { background: $nav_hover_bg; color: $ink; }
QPushButton#nav:checked { background: $nav_on_bg; border-color: $nav_on_border; color: $nav_on_ink; font-weight: 700; }
QPushButton#nav:focus { border-color: $accent; }
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QDateEdit, QDateTimeEdit, QTimeEdit { border: 1px solid $control_border; border-radius: 7px; background: $field_bg; padding: $field_padding; selection-background-color: $selection; selection-color: $on_selection; }
QLineEdit:hover, QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: $field_hover_border; }
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border-color: $accent; }
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled { background: $field_off_bg; color: $ink_off; border-color: $control_border_off; }
QComboBox { padding-right: 28px; }
QComboBox::drop-down { border-left: 1px solid $control_border; width: 25px; }
QComboBox::down-arrow { image: url('$icons/$chevron_down'); width: 12px; height: 12px; }
QComboBox QAbstractItemView { background: $menu_bg; border: 1px solid $menu_border; padding: 4px; color: $ink; selection-background-color: $menu_selection; selection-color: $on_selection; outline: 0; }
QSpinBox, QDoubleSpinBox { padding-right: 25px; }
QSpinBox::up-button, QDoubleSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::down-button { width: 22px; border-left: 1px solid $control_border; background: $stepper_bg; }
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover, QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover { background: $stepper_hover_bg; }
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow { image: url('$icons/$chevron_up'); width: 10px; height: 10px; }
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow { image: url('$icons/$chevron_down'); width: 10px; height: 10px; }
QAbstractItemView { background: $surface; alternate-background-color: $row_alt; color: $ink; selection-background-color: $selection; selection-color: $on_row_selected; border: 1px solid $view_border; outline: 0; }
QTableView, QTableWidget { gridline-color: $gridline; border: none; }
QTableView::item, QTableWidget::item { $cell_rule padding: $cell_padding; }
QAbstractItemView::item:selected { background: $selection; color: $on_row_selected; }
QAbstractItemView::item:hover:!selected { background: $row_hover; }
QAbstractItemView::item:focus { border: 1px solid $cell_focus; }
QTreeView::branch { background: transparent; }
QHeaderView { background: $header_bg; }
QHeaderView::section { background: $header_section_bg; color: $header_ink; font-size: 11px; font-weight: 600; padding: $header_padding; border: none; border-right: 1px solid $header_rule; border-bottom: 1px solid $control_border; }
QTableCornerButton::section { background: $header_section_bg; border: none; }
QScrollArea { border: none; background: $ground; }
QScrollBar:vertical { width: 11px; background: $scrollbar_bg; margin: 2px; }
QScrollBar:horizontal { height: 11px; background: $scrollbar_bg; margin: 2px; }
QScrollBar::handle:vertical { background: $scrollbar_handle; min-height: 30px; border-radius: 4px; }
QScrollBar::handle:horizontal { background: $scrollbar_handle; min-width: 30px; border-radius: 4px; }
QScrollBar::handle:hover { background: $scrollbar_handle_hover; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0px; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0px; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
QProgressBar { border: none; background: $progress_track; border-radius: 4px; max-height: 7px; color: transparent; }
QProgressBar::chunk { background: $accent; border-radius: 4px; }
QTextBrowser, QTextEdit, QPlainTextEdit { border: 1px solid $view_border; border-radius: 6px; background: $surface; color: $ink; padding: 10px; selection-background-color: $selection; selection-color: $on_selection; }
QMenuBar { background: $sidebar; color: $menubar_ink; border-bottom: 1px solid $rule; }
QMenuBar::item { background: transparent; padding: 7px 10px; }
QMenuBar::item:selected, QMenuBar::item:pressed { background: $menubar_selection; }
QMenu { background: $menu_bg; color: $ink; border: 1px solid $menu_border; padding: 5px; }
QMenu::item { padding: 8px 28px 8px 16px; border-radius: 4px; }
QMenu::item:selected { background: $menu_selection; color: $on_selection; }
QMenu::item:disabled { color: $ink_off; }
QMenu::separator { height: 1px; background: $control_border; margin: 5px 8px; }
QToolTip { color: $tooltip_ink; background: $tooltip_bg; border: 1px solid $tooltip_border; padding: 7px; }
QCheckBox, QRadioButton { spacing: 8px; background: transparent; }
QCheckBox::indicator, QAbstractItemView::indicator { width: 17px; height: 17px; border: 1px solid $check_border; border-radius: 4px; background: $field_bg; }
QCheckBox::indicator:hover, QAbstractItemView::indicator:hover { border-color: $badge_ink; background: $check_hover_bg; }
QCheckBox::indicator:checked, QAbstractItemView::indicator:checked { image: url('$icons/check.svg'); border-color: $accent; background: $check_on_bg; }
QCheckBox::indicator:indeterminate, QAbstractItemView::indicator:indeterminate { image: url('$icons/minus.svg'); border-color: $accent; background: $check_on_bg; }
QCheckBox::indicator:disabled { border-color: $check_off_border; background: $check_off_bg; }
QRadioButton::indicator { width: 18px; height: 18px; }
QGroupBox { border: 1px solid $control_border; border-radius: 7px; margin-top: 12px; padding-top: 10px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; color: $group_title_ink; }
QTabWidget::pane { border: 1px solid $control_border; background: $surface; }
QTabBar::tab { background: $tab_bg; color: $muted; padding: $tab_padding; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { background: $tab_on_bg; color: $tab_on_ink; border-bottom-color: $accent; }
QTabBar::tab:hover { background: $tab_hover_bg; }
/* Scoped to the workspace tab bar by id, and written without descendant
   selectors: Qt re-polishes every live widget on setStyleSheet, and a
   descendant rule makes that walk every ancestor chain. */
QTabWidget#workspaceTabs::pane { border: none; border-top: 1px solid $rule; background: $ground; }
QTabBar#workspaceTabBar { background: transparent; }
QTabBar#workspaceTabBar::tab { background: transparent; padding: $workspace_tab_padding; font-weight: 600; }
QTabBar#workspaceTabBar::tab:selected { background: $workspace_tab_on_bg; color: $tab_on_ink; }
QFrame#purposeStrip { background: $field_bg; border: 1px solid $rule; border-radius: 9px; }
QLabel#purpose { color: $purpose_ink; font-size: 12px; }
QLabel#cohortScope { color: $tab_on_ink; font-size: 11px; background: transparent; }
QPushButton#nextStep, QPushButton#pageGuide, QPushButton#cohortAdopt { padding: 4px 10px; font-size: 12px; }
QStatusBar { background: $sidebar; color: $muted; border-top: 1px solid $rule; font-size: 11px; }
QStatusBar::item { border: none; }
QSplitter::handle { background: $view_border; width: 2px; height: 2px; }
QSplitter::handle:hover { background: $accent_dim; }
QSizeGrip { background: $sidebar; }
"""


def build_style(tokens=None, density=DEFAULT_DENSITY):
    """One stylesheet from one palette and one set of measurements.

    Every `font-size: Npx` in the template is left in that exact spelling:
    `interface_settings.scaled_style_fragment` rewrites them by regular
    expression, and would silently stop scaling text if the form changed.
    """
    tokens = dict(theme_tokens(None) if tokens is None else tokens)
    metrics = density_metrics(density) if isinstance(density, str) else dict(density)
    values = dict(tokens)
    values.update({key: value for key, value in metrics.items() if isinstance(value, str)})
    # A grid line and a per-cell bottom border would draw two rules on the same
    # edge, which is what smudged the old sheet wherever both were turned on.
    values["cell_rule"] = ("" if metrics["grid"]
                           else f"border-bottom: 1px solid {tokens['row_rule']};")
    values["icons"] = _ICONS
    return Template(_TEMPLATE).substitute(values)


_active = {"theme": DEFAULT_THEME, "density": DEFAULT_DENSITY}
_cache = {}


def active_theme() -> str:
    return _active["theme"]


def active_density() -> str:
    return _active["density"]


def style_for(theme=None, density=None) -> str:
    """The sheet for one theme and density, built once and remembered."""
    theme = str(theme or DEFAULT_THEME).strip().lower()
    theme = theme if theme in THEMES else DEFAULT_THEME
    density = str(density or DEFAULT_DENSITY).strip().lower()
    density = density if density in TABLE_DENSITY else DEFAULT_DENSITY
    if (theme, density) not in _cache:
        _cache[(theme, density)] = build_style(THEMES[theme], density)
    return _cache[(theme, density)]


def active_style() -> str:
    return style_for(_active["theme"], _active["density"])


def set_active_theme(theme=None, density=None) -> str:
    """Choose the theme and/or density, and hand back the sheet to apply.

    The window applies stylesheets through `ui_workbench.apply_application_style`,
    which deliberately skips an assignment that would not change the string. A
    different theme or density produces a different string, so the change lands.
    """
    if theme is not None:
        name = str(theme).strip().lower()
        _active["theme"] = name if name in THEMES else DEFAULT_THEME
    if density is not None:
        name = str(density).strip().lower()
        _active["density"] = name if name in TABLE_DENSITY else DEFAULT_DENSITY
    return active_style()


STYLE = active_style()
