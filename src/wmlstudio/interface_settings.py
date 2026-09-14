"""Simple application settings first; advanced controls stay behind disclosure."""

import re

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
    QScrollArea, QVBoxLayout, QWidget,
)

from wmlstudio import __version__
from wmlstudio.theme import STYLE
from wmlstudio.widgets import button, label


def interface_preferences(root):
    return QSettings(str(root / "Interface.ini"), QSettings.Format.IniFormat)


def scaled_style(percent):
    if isinstance(percent, bool) or not isinstance(percent, int) or not 80 <= percent <= 150:
        raise ValueError("Interface scale must be between 80% and 150%")
    factor = percent / 100
    return re.sub(r"font-size:\s*([0-9.]+)px", lambda match:
                  f"font-size: {round(float(match[1]) * factor, 2):g}px", STYLE)


class InterfaceSettingsDialog(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.preferences = interface_preferences(window.root)
        self.setWindowTitle("WMLSTudio settings")
        self.resize(640, 470)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Make the workspace comfortable", "title", True))
        form = QFormLayout()
        self.scale = QComboBox()
        for value in (80, 90, 100, 110, 125, 150):
            self.scale.addItem(f"{value}%" + (" · recommended" if value == 100 else ""), value)
        self.scale.setCurrentIndex(max(0, self.scale.findData(window.ui_scale)))
        self.scale.currentIndexChanged.connect(self.apply_scale)
        form.addRow("Interface text size", self.scale)
        self.motion = QCheckBox("Gentle logo animation")
        self.motion.setChecked(window.motion_enabled)
        self.motion.toggled.connect(window.set_motion)
        form.addRow("Motion", self.motion)
        self.policy = QComboBox()
        for title, value in (("Balanced", "balanced"), ("Faster throughput", "fast"), ("Low memory", "low_memory")):
            self.policy.addItem(title, value)
        self.policy.setCurrentIndex(max(0, self.policy.findData(self.preferences.value("resource_policy", "balanced"))))
        self.policy.currentIndexChanged.connect(lambda: self.preferences.setValue("resource_policy", self.policy.currentData()))
        form.addRow("Default analysis resources", self.policy)
        layout.addLayout(form)
        layout.addWidget(label("Changes apply immediately. Graph zoom remains independent: use the graph's zoom/fit controls. Each analysis and report reviews its own cohort.", "muted", True))
        self.advanced_button = button("Advanced settings ▸", self.toggle_advanced)
        self.advanced_button.setCheckable(True)
        layout.addWidget(self.advanced_button, alignment=Qt.AlignmentFlag.AlignLeft)
        self.advanced = QScrollArea()
        self.advanced.setWidgetResizable(True)
        panel = QWidget()
        detail = QVBoxLayout(panel)
        group = QGroupBox("References, files and reproducibility")
        fields = QVBoxLayout(group)
        fields.addWidget(label(f"WMLSTudio {__version__}\nPortable application data: {window.root}\nActive project: {window.project_path}", "small", True))
        fields.addWidget(button("Scheme references…", window.open_reference_manager))
        fields.addWidget(button("AMR / plasmid reference databases…", window.open_amr_databases))
        fields.addWidget(button("Characterization reference panel…", window.install_characterization_references))
        fields.addWidget(button("Open application data folder", window.open_data_folder))
        fields.addWidget(button("Save project backup copy…", window.save_project_copy))
        detail.addWidget(group)
        detail.addWidget(label("Advanced scientific thresholds belong to the reviewed analysis or investigation, not a hidden global default. Updates create immutable reference snapshots; previous results keep their original provenance.", "small", True))
        self.advanced.setWidget(panel)
        self.advanced.hide()
        layout.addWidget(self.advanced, 1)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)

    def apply_scale(self):
        self.window.set_ui_scale(self.scale.currentData())

    def toggle_advanced(self):
        showing = self.advanced_button.isChecked()
        self.advanced.setVisible(showing)
        self.advanced_button.setText("Advanced settings ▾" if showing else "Advanced settings ▸")
        self.resize(self.width(), 690 if showing else 470)
