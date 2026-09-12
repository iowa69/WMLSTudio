"""Native widget theme; contrast and motion are controlled by the application."""

INK = "#203736"
TEAL = "#147D70"
MUTED = "#687A78"
PALETTE = ["#168879", "#6485C4", "#DC9368", "#9C7CBA", "#BD7485", "#71998B"]

STYLE = """
QWidget { color: #203736; font-family: 'Segoe UI', 'Inter', 'DejaVu Sans'; font-size: 13px; }
QMainWindow, QWidget#main, QWidget#page, QStackedWidget { background: #F6F8F7; }
QWidget#sidebar { background: #FFFFFF; border-right: 1px solid #E3EAE7; }
QLabel { background: transparent; }
QLabel#brand { font-size: 21px; font-weight: 700; letter-spacing: -1px; }
QLabel#eyebrow { font-size: 10px; font-weight: 700; color: #7B8D88; letter-spacing: 2px; }
QLabel#title { font-size: 29px; font-weight: 700; letter-spacing: -1px; }
QLabel#heroTitle { font-size: 34px; font-weight: 700; letter-spacing: -1px; }
QLabel#muted { color: #687A78; }
QLabel#small { color: #687A78; font-size: 11px; }
QLabel#metric { font-size: 30px; font-weight: 600; }
QLabel#cardTitle { font-size: 15px; font-weight: 600; }
QLabel#badge { background: #E5F4ED; color: #247963; border-radius: 12px; padding: 5px 10px; font-size: 11px; }
QFrame#card { background: white; border: 1px solid #E1E9E5; border-radius: 12px; }
QFrame#hero { background: #E9F3EE; border: 1px solid #DAE9DF; border-radius: 15px; }
QFrame#drop { background: #FFFFFF; border: 1px dashed #AABFB5; border-radius: 12px; }
QFrame#note { background: #FFF6E8; border: 1px solid #EDDEBB; border-radius: 9px; }
QPushButton { background: white; border: 1px solid #D5E0DA; border-radius: 7px; padding: 9px 14px; font-weight: 600; }
QPushButton:hover { background: #EEF5F1; border-color: #A9C4B7; }
QPushButton:pressed { background: #DDECE4; }
QPushButton:disabled { color: #9BA8A1; background: #EEF1EF; border-color: #E5EAE7; }
QPushButton:focus { border: 2px solid #228E80; }
QPushButton#primary { color: white; background: #187D6D; border-color: #187D6D; }
QPushButton#primary:hover { background: #126859; }
QPushButton#primary:disabled { background: #ADC9BC; border-color: #ADC9BC; }
QPushButton#nav { border: none; background: transparent; text-align: left; padding: 12px 17px; color: #677C73; font-weight: 500; }
QPushButton#nav:hover { background: #F0F5F2; }
QPushButton#nav:checked { background: #E7F2EC; color: #146D5F; font-weight: 700; }
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox { border: 1px solid #D8E3DC; border-radius: 7px; background: white; padding: 9px; selection-background-color: #218577; }
QComboBox::drop-down { border: none; width: 25px; }
QComboBox QAbstractItemView { background: white; selection-background-color: #E6F2EB; color: #203736; }
QTableWidget { background: white; border: none; gridline-color: #EDF1EE; alternate-background-color: #FBFCFB; selection-background-color: #E5F2ED; selection-color: #1B564A; }
QTableWidget::item { border-bottom: 1px solid #EDF1EE; padding: 8px; }
QTableWidget::item:selected { background: #E5F2ED; color: #1B564A; }
QHeaderView::section { background: #F8FAF8; color: #788A80; font-size: 10px; font-weight: 600; padding: 12px 8px; border: none; border-bottom: 1px solid #E5ECE7; }
QScrollArea { border: none; background: transparent; }
QScrollBar:vertical { width: 8px; background: transparent; margin: 2px; }
QScrollBar::handle:vertical { background: #CDDCD3; min-height: 30px; border-radius: 3px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0px; }
QProgressBar { border: none; background: #E1EBE5; border-radius: 4px; max-height: 7px; color: transparent; }
QProgressBar::chunk { background: #238776; border-radius: 4px; }
QTextBrowser { border: none; background: white; padding: 12px; }
QToolTip { color: #FFFFFF; background: #294D40; border: none; padding: 7px; }
QCheckBox { spacing: 8px; }
QCheckBox::indicator { width: 17px; height: 17px; }
QStatusBar { background: #FFFFFF; color: #7A8A82; border-top: 1px solid #E5ECE7; font-size: 11px; }
QSplitter::handle { background: #E3EBE5; width: 1px; }
"""
