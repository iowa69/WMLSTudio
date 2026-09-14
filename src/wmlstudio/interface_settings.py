"""Simple application settings first; advanced controls stay behind disclosure."""

import re

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from wmlstudio import __version__, display
from wmlstudio.theme import STYLE
from wmlstudio.widgets import button, label


def interface_preferences(root):
    return QSettings(str(root / "Interface.ini"), QSettings.Format.IniFormat)


def scaled_style_fragment(text, percent):
    """Rewrite every `font-size: Npx` in a stylesheet fragment at the chosen scale."""
    factor = percent / 100
    return re.sub(r"font-size:\s*([0-9.]+)px", lambda match:
                  f"font-size: {round(float(match[1]) * factor, 2):g}px", text)


def scaled_style(percent):
    if isinstance(percent, bool) or not isinstance(percent, int) or not 80 <= percent <= 150:
        raise ValueError("Interface scale must be between 80% and 150%")
    return scaled_style_fragment(STYLE, percent)


class InterfaceSettingsPanel(QWidget):
    """Display, text size and comfort controls, shared by the Settings tab and dialog.

    Two levers, kept apart because they do different things and the difference
    matters to someone deciding what to change: text size rewrites the stylesheet
    and applies at once, while the whole-interface scale is read by Qt when the
    application starts and therefore needs a restart. Neither changes a distance,
    a threshold or any exported image.
    """

    def __init__(self, window, *, references=True, display_root=None, parent=None):
        super().__init__(parent)
        self.window_ref = window
        self.preferences = interface_preferences(window.root)
        # Display keys live with the data root, not the project root: they are read
        # before any project is open, while the text scale stays per workspace.
        # `display_root` only exists so a test can redirect that file.
        self.display_root = display_root
        self.display_settings = display.read_display_settings(self.display_root)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        form = QFormLayout()
        self.scale = QComboBox()
        for value in display.TEXT_SCALE_CHOICES:
            self.scale.addItem(f"{value}%" + (" · recommended" if value == 100 else ""), value)
        self.scale.setCurrentIndex(max(0, self.scale.findData(window.ui_scale)))
        self.scale.currentIndexChanged.connect(self.apply_scale)
        form.addRow("Interface text size", self.scale)
        self.display_scale = QComboBox()
        self.display_scale.addItem("Follow my system setting (recommended)", 0)
        for value in display.available_scale_choices(self.available_size()):
            if value != 100:
                self.display_scale.addItem(f"{value}%", value)
        chosen = self.display_settings["scale_percent"] if self.display_settings["mode"] == "fixed" else 0
        self.display_scale.setCurrentIndex(max(0, self.display_scale.findData(chosen)))
        self.display_scale.currentIndexChanged.connect(self.apply_display_scale)
        form.addRow("Whole interface size", self.display_scale)
        self.graph_scale = QComboBox()
        for value in display.TEXT_SCALE_CHOICES:
            self.graph_scale.addItem(f"{value}%", value)
        self.graph_scale.setCurrentIndex(max(0, self.graph_scale.findData(self.saved_graph_scale())))
        self.graph_scale.currentIndexChanged.connect(self.apply_graph_scale)
        form.addRow("Graph text size", self.graph_scale)
        self.motion = QCheckBox("Gentle logo animation")
        self.motion.setChecked(bool(window.project.get_setting("motion", True)))
        self.motion.toggled.connect(window.set_motion)
        form.addRow("Motion", self.motion)
        self.policy = QComboBox()
        for title, value in (("Balanced", "balanced"), ("Faster throughput", "fast"), ("Low memory", "low_memory")):
            self.policy.addItem(title, value)
        self.policy.setCurrentIndex(max(0, self.policy.findData(self.preferences.value("resource_policy", "balanced"))))
        self.policy.currentIndexChanged.connect(lambda: self.preferences.setValue("resource_policy", self.policy.currentData()))
        form.addRow("Default analysis resources", self.policy)
        layout.addLayout(form)
        self.restart_notice = label(display.RESTART_NOTICE, "badge")
        self.restart_notice.setWordWrap(True)
        self.restart_notice.setVisible(False)
        layout.addWidget(self.restart_notice)
        self.hidden_notice = label(display.HIDDEN_CHOICES_NOTICE, "small", True)
        self.hidden_notice.setVisible(
            len(display.available_scale_choices(self.available_size())) < len(display.SCALE_CHOICES))
        layout.addWidget(self.hidden_notice)
        layout.addWidget(label(display.EXPORT_NOTICE, "small", True))
        layout.addWidget(self.build_preview())
        layout.addWidget(label("Text size applies immediately. Graph zoom stays independent: use the graph's own zoom and fit controls. Each analysis and report reviews its own cohort.", "muted", True))
        layout.addStretch(1)
        self.advanced_button = button("Advanced settings ▸", self.toggle_advanced)
        self.advanced_button.setCheckable(True)
        layout.addWidget(self.advanced_button, alignment=Qt.AlignmentFlag.AlignLeft)
        self.advanced = QScrollArea()
        self.advanced.setWidgetResizable(True)
        panel = QWidget()
        detail = QVBoxLayout(panel)
        self.rounding = QComboBox()
        for key, title in display.ROUNDING_LABELS.items():
            self.rounding.addItem(title, key)
        self.rounding.setCurrentIndex(max(0, self.rounding.findData(self.display_settings["rounding"])))
        self.rounding.currentIndexChanged.connect(self.apply_rounding)
        rounding_row = QFormLayout()
        rounding_row.addRow("Screen scale rounding", self.rounding)
        detail.addLayout(rounding_row)
        self.state_line = label(display.describe_display_state(self.resolved_state()), "small", True)
        detail.addWidget(self.state_line)
        detail.addWidget(label(display.RECOVERY_NOTICE, "small", True))
        if references:
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
        detail.addStretch()
        self.advanced.setWidget(panel)
        self.advanced.hide()
        layout.addWidget(self.advanced, 2)

    # --- helpers ------------------------------------------------------------
    def available_size(self):
        """The usable screen, so a scale the screen cannot show is never offered."""
        screen = self.window_ref.screen() if hasattr(self.window_ref, "screen") else None
        if screen is None:
            return (display.MINIMUM_WINDOW_WIDTH, display.MINIMUM_WINDOW_HEIGHT)
        return screen.availableGeometry().size()

    def saved_graph_scale(self):
        try:
            return display.clamp_text_scale(int(self.preferences.value("graph_scale", 100)))
        except (TypeError, ValueError):
            return 100

    def resolved_state(self):
        state = display.apply_display_settings(root=self.display_root)
        if state.get("reason") == "already_started":
            # Qt read its policy when this application started, so report what is
            # saved; the restart notice says when that is not yet what you see.
            state = {**state, **self.display_settings,
                     "reason": "preference" if self.display_settings["mode"] == "fixed" else "system"}
        return state

    def build_preview(self):
        """A sample row the user can judge before restarting anything."""
        frame = QFrame()
        frame.setObjectName("card")
        row = QHBoxLayout(frame)
        row.setContentsMargins(14, 10, 14, 10)
        self.preview_title = label("Isolate KP-014 · ST 258 · 7 of 7 loci called", "cardTitle")
        row.addWidget(self.preview_title, 1)
        self.preview_note = label("This is how a row of your evidence will read.", "small", True)
        row.addWidget(self.preview_note)
        self.preview_button = button("Example action")
        self.preview_button.setEnabled(False)
        row.addWidget(self.preview_button)
        return frame

    def refresh_preview(self):
        self.preview_note.setText(
            f"Text {self.window_ref.ui_scale}% · graph text {self.graph_scale.currentData()}% · "
            f"{self.display_scale.currentText()}")

    # --- the controls -------------------------------------------------------
    def apply_scale(self):
        self.window_ref.set_ui_scale(self.scale.currentData())
        self.refresh_preview()

    def apply_display_scale(self):
        percent = self.display_scale.currentData()
        self.display_settings = display.write_display_settings(
            mode="fixed" if percent else "system",
            scale_percent=percent if percent else None, root=self.display_root)
        self.restart_notice.setVisible(True)
        self.state_line.setText(display.describe_display_state(self.resolved_state()))
        self.refresh_preview()

    def apply_rounding(self):
        self.display_settings = display.write_display_settings(
            rounding=self.rounding.currentData(), root=self.display_root)
        self.restart_notice.setVisible(True)
        self.state_line.setText(display.describe_display_state(self.resolved_state()))

    def apply_graph_scale(self):
        percent = display.clamp_text_scale(self.graph_scale.currentData())
        self.preferences.setValue("graph_scale", percent)
        setter = getattr(self.window_ref, "set_graph_text_scale", None)
        if callable(setter):
            setter(percent)
        self.refresh_preview()

    def toggle_advanced(self):
        showing = self.advanced_button.isChecked()
        self.advanced.setVisible(showing)
        self.advanced_button.setText("Advanced settings ▾" if showing else "Advanced settings ▸")


class InterfaceSettingsDialog(QDialog):
    """The modal host for the same panel the Settings tab shows."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.setWindowTitle("WMLSTudio settings")
        self.resize(640, 560)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Make the workspace comfortable", "title", True))
        self.panel = InterfaceSettingsPanel(window, parent=self)
        layout.addWidget(self.panel, 1)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)

    # Kept so the existing entry points and tests keep addressing one object.
    @property
    def scale(self):
        return self.panel.scale

    @property
    def motion(self):
        return self.panel.motion

    @property
    def advanced(self):
        return self.panel.advanced

    @property
    def advanced_button(self):
        return self.panel.advanced_button

    def apply_scale(self):
        self.panel.apply_scale()

    def toggle_advanced(self):
        self.panel.toggle_advanced()
        self.resize(self.width(), 780 if self.panel.advanced_button.isChecked() else 560)
