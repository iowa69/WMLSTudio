"""Opaque, high-contrast native desktop surfaces and a matching Qt palette."""

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette

BACKGROUND = "#0B1220"
SURFACE = "#121E30"
BORDER = "#2B3C54"
INK = "#E7EFF9"
TEAL = "#48DCC0"
MUTED = "#A0B1C5"
PALETTE = ["#48DCC0", "#73ACFF", "#C29AFF", "#FFB976", "#FF86A4",
           "#B5D976", "#79D6EF", "#DBB0E8", "#E7D577", "#A3B5D2"]


def apply_dark_palette(application):
    """Style native Fusion controls and Qt-owned dialogs, including disabled states."""
    palette = QPalette()
    colors = {
        "Window": BACKGROUND, "WindowText": INK, "Base": SURFACE,
        "AlternateBase": "#17243A", "ToolTipBase": "#22324A", "ToolTipText": INK,
        "Text": INK, "Button": "#1B2A40", "ButtonText": INK, "BrightText": "#FF91A6",
        "Link": "#80C8FF", "LinkVisited": "#C29AFF", "Highlight": "#28526A",
        "HighlightedText": "#F5FBFF", "PlaceholderText": "#7D91AA", "Light": "#3B4E68",
        "Midlight": "#2A3D57", "Mid": "#25364F", "Dark": "#070D17", "Shadow": "#03070D",
    }
    for role, color in colors.items():
        palette.setColor(getattr(QPalette.ColorRole, role), QColor(color))
    for role in (QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText,
                 QPalette.ColorRole.WindowText, QPalette.ColorRole.PlaceholderText):
        palette.setColor(QPalette.ColorGroup.Disabled, role, QColor("#677C94"))
    application.setPalette(palette)
    application.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, True)


STYLE = """
QWidget { color: #E7EFF9; font-family: 'Segoe UI', 'Inter', 'DejaVu Sans'; font-size: 13px; }
QMainWindow, QWidget#main, QWidget#page, QStackedWidget, QDialog { background: #0B1220; }
QWidget#sidebar { background: #0F1929; border-right: 1px solid #25354C; }
QLabel { background: transparent; }
QLabel:disabled { color: #677C94; }
QLabel#brand { font-size: 22px; font-weight: 700; letter-spacing: -1px; }
QLabel#eyebrow { font-size: 10px; font-weight: 700; color: #8FA9C2; letter-spacing: 2px; }
QLabel#title { font-size: 28px; font-weight: 700; letter-spacing: -1px; }
QLabel#heroTitle { font-size: 33px; font-weight: 700; letter-spacing: -1px; }
QLabel#muted { color: #A0B1C5; }
QLabel#small { color: #A0B1C5; font-size: 11px; }
QLabel#metric { color: #F2F7FF; font-size: 30px; font-weight: 600; }
QLabel#cardTitle { font-size: 15px; font-weight: 600; }
QLabel#badge { background: #183C3E; color: #83E9D3; border: 1px solid #285755; border-radius: 12px; padding: 5px 10px; font-size: 11px; }
QFrame#card { background: #121E30; border: 1px solid #293A53; border-radius: 12px; }
QFrame#hero { background: #132B35; border: 1px solid #28505A; border-radius: 15px; }
QFrame#drop { background: #122333; border: 1px dashed #467786; border-radius: 12px; }
QFrame#drop:hover { background: #16303D; border-color: #48DCC0; }
QFrame#note { background: #302C24; border: 1px solid #6B593B; border-radius: 9px; }
QPushButton, QToolButton { background: #1B2A40; border: 1px solid #354A65; border-radius: 7px; padding: 9px 14px; font-weight: 600; }
QPushButton:hover, QToolButton:hover { background: #263B55; border-color: #527392; }
QPushButton:pressed, QToolButton:pressed { background: #30516A; }
QPushButton:checked, QToolButton:checked { background: #244B59; border-color: #48BFAE; color: #C6FFF2; }
QPushButton:disabled, QToolButton:disabled { color: #677C94; background: #141F30; border-color: #26364B; }
QPushButton:focus, QToolButton:focus { border-color: #69E1CC; }
QPushButton#primary { color: #082B29; background: #48DCC0; border-color: #48DCC0; }
QPushButton#primary:hover { background: #79EDD6; border-color: #79EDD6; }
QPushButton#primary:pressed { background: #2AB69E; }
QPushButton#primary:disabled { color: #829B9F; background: #254943; border-color: #30554F; }
QPushButton#nav { border: 1px solid transparent; background: transparent; text-align: left; padding: 12px 16px; color: #A3B4CA; font-weight: 500; }
QPushButton#nav:hover { background: #1A2B42; color: #E7EFF9; }
QPushButton#nav:checked { background: #183C43; border-color: #28555A; color: #85EDD7; font-weight: 700; }
QPushButton#nav:focus { border-color: #48DCC0; }
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QDateEdit, QDateTimeEdit, QTimeEdit { border: 1px solid #354A65; border-radius: 7px; background: #101B2C; padding: 8px; selection-background-color: #28526A; selection-color: #FFFFFF; }
QLineEdit:hover, QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: #55718F; }
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border-color: #48DCC0; }
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled { background: #131D2B; color: #677C94; border-color: #26364B; }
QComboBox { padding-right: 28px; }
QComboBox::drop-down { border-left: 1px solid #354A65; width: 25px; }
QComboBox::down-arrow { image: url('@ICONS@/chevron-down.svg'); width: 12px; height: 12px; }
QComboBox QAbstractItemView { background: #17253A; border: 1px solid #45617E; padding: 4px; color: #E7EFF9; selection-background-color: #2B5067; selection-color: #FFFFFF; outline: 0; }
QSpinBox, QDoubleSpinBox { padding-right: 25px; }
QSpinBox::up-button, QDoubleSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::down-button { width: 22px; border-left: 1px solid #354A65; background: #203149; }
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover, QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover { background: #35526E; }
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow { image: url('@ICONS@/chevron-up.svg'); width: 10px; height: 10px; }
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow { image: url('@ICONS@/chevron-down.svg'); width: 10px; height: 10px; }
QAbstractItemView { background: #121E30; alternate-background-color: #17243A; color: #E7EFF9; selection-background-color: #28526A; selection-color: #F5FBFF; border: 1px solid #2B3C54; outline: 0; }
QTableView, QTableWidget { gridline-color: #25374E; border: none; }
QTableView::item, QTableWidget::item { border-bottom: 1px solid #24354A; padding: 8px; }
QAbstractItemView::item:selected { background: #28526A; color: #F5FBFF; }
QAbstractItemView::item:hover:!selected { background: #21354D; }
QAbstractItemView::item:focus { border: 1px solid #6DD7C6; }
QTreeView::branch { background: transparent; }
QHeaderView { background: #16253A; }
QHeaderView::section { background: #19293E; color: #ABC0D6; font-size: 11px; font-weight: 600; padding: 11px 8px; border: none; border-right: 1px solid #293C55; border-bottom: 1px solid #354A65; }
QTableCornerButton::section { background: #19293E; border: none; }
QScrollArea { border: none; background: #0B1220; }
QScrollBar:vertical { width: 11px; background: #101A2A; margin: 2px; }
QScrollBar:horizontal { height: 11px; background: #101A2A; margin: 2px; }
QScrollBar::handle:vertical { background: #3A506D; min-height: 30px; border-radius: 4px; }
QScrollBar::handle:horizontal { background: #3A506D; min-width: 30px; border-radius: 4px; }
QScrollBar::handle:hover { background: #577895; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0px; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0px; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
QProgressBar { border: none; background: #26384C; border-radius: 4px; max-height: 7px; color: transparent; }
QProgressBar::chunk { background: #48DCC0; border-radius: 4px; }
QTextBrowser, QTextEdit, QPlainTextEdit { border: 1px solid #2B3C54; border-radius: 6px; background: #121E30; color: #E7EFF9; padding: 10px; selection-background-color: #28526A; selection-color: #FFFFFF; }
QMenuBar { background: #0F1929; color: #CCD9E8; border-bottom: 1px solid #25354C; }
QMenuBar::item { background: transparent; padding: 7px 10px; }
QMenuBar::item:selected, QMenuBar::item:pressed { background: #28445D; }
QMenu { background: #17253A; color: #E7EFF9; border: 1px solid #45617E; padding: 5px; }
QMenu::item { padding: 8px 28px 8px 16px; border-radius: 4px; }
QMenu::item:selected { background: #2B5067; color: #FFFFFF; }
QMenu::item:disabled { color: #677C94; }
QMenu::separator { height: 1px; background: #354A65; margin: 5px 8px; }
QToolTip { color: #F0F7FF; background: #22324A; border: 1px solid #52708E; padding: 7px; }
QCheckBox, QRadioButton { spacing: 8px; background: transparent; }
QCheckBox::indicator, QAbstractItemView::indicator { width: 17px; height: 17px; border: 1px solid #64809C; border-radius: 4px; background: #101B2C; }
QCheckBox::indicator:hover, QAbstractItemView::indicator:hover { border-color: #83E9D3; background: #203C49; }
QCheckBox::indicator:checked, QAbstractItemView::indicator:checked { image: url('@ICONS@/check.svg'); border-color: #48DCC0; background: #204C4C; }
QCheckBox::indicator:indeterminate, QAbstractItemView::indicator:indeterminate { image: url('@ICONS@/minus.svg'); border-color: #48DCC0; background: #204C4C; }
QCheckBox::indicator:disabled { border-color: #36495F; background: #152132; }
QRadioButton::indicator { width: 18px; height: 18px; }
QGroupBox { border: 1px solid #354A65; border-radius: 7px; margin-top: 12px; padding-top: 10px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; color: #BDD0E2; }
QTabWidget::pane { border: 1px solid #354A65; background: #121E30; }
QTabBar::tab { background: #152236; color: #A0B1C5; padding: 9px 15px; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { background: #20374A; color: #86EDD7; border-bottom-color: #48DCC0; }
QTabBar::tab:hover { background: #24394F; }
/* Scoped to the workspace tab bar by id, and written without descendant
   selectors: Qt re-polishes every live widget on setStyleSheet, and a
   descendant rule makes that walk every ancestor chain. */
QTabWidget#workspaceTabs::pane { border: none; border-top: 1px solid #25354C; background: #0B1220; }
QTabBar#workspaceTabBar { background: transparent; }
QTabBar#workspaceTabBar::tab { background: transparent; padding: 8px 12px; font-weight: 600; }
QTabBar#workspaceTabBar::tab:selected { background: #142739; color: #86EDD7; }
QFrame#purposeStrip { background: #101B2C; border: 1px solid #25354C; border-radius: 9px; }
QLabel#purpose { color: #C2D3E6; font-size: 12px; }
QLabel#cohortScope { color: #86EDD7; font-size: 11px; background: transparent; }
QPushButton#nextStep, QPushButton#pageGuide, QPushButton#cohortAdopt { padding: 4px 10px; font-size: 12px; }
QStatusBar { background: #0F1929; color: #A0B1C5; border-top: 1px solid #25354C; font-size: 11px; }
QStatusBar::item { border: none; }
QSplitter::handle { background: #2B3C54; width: 2px; height: 2px; }
QSplitter::handle:hover { background: #48BFAE; }
QSizeGrip { background: #0F1929; }
""".replace("@ICONS@", (Path(__file__).resolve().parent / "resources/ui").as_posix())
