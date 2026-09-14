"""User-adjustable display scaling: persisted, screen-aware, applied before the app exists."""

import json
import subprocess
import sys

import pytest
from PySide6.QtCore import QSize

from wmlstudio import display
from wmlstudio.interface_settings import scaled_style
from wmlstudio.theme import STYLE

# apply_display_settings has to run while no QApplication exists, which is exactly
# what a pytest-qt session cannot offer: the application object outlives every test.
# These probes therefore run in a fresh interpreter.
PROBE = """
import json, os, sys
from PySide6.QtCore import QSize
from PySide6.QtGui import QGuiApplication
from wmlstudio import display

root, override, screen = sys.argv[1], sys.argv[2], sys.argv[3]
available = QSize(*(int(part) for part in screen.split("x"))) if screen else None
state = display.apply_display_settings(int(override) if override else None, root=root,
                                       available_size=available)
state["QT_SCALE_FACTOR"] = os.environ.get("QT_SCALE_FACTOR")
state["policy"] = str(QGuiApplication.highDpiScaleFactorRoundingPolicy())
print(json.dumps(state))
"""


def probe(root, override="", environment=None, screen=""):
    env = {"QT_QPA_PLATFORM": "offscreen", "PATH": "/usr/bin:/bin"}
    env.update(environment or {})
    finished = subprocess.run(
        [sys.executable, "-c", PROBE, str(root), str(override), str(screen)],
        capture_output=True, text=True, env=env, timeout=120)
    assert finished.returncode == 0, finished.stderr
    return json.loads(finished.stdout.strip().splitlines()[-1])


def test_settings_fall_back_to_safe_values_when_the_file_is_missing(tmp_path):
    assert display.read_display_settings(tmp_path) == {
        "mode": "system", "scale_percent": 100, "rounding": "exact"}


def test_a_corrupt_or_unknown_preference_is_ignored_rather_than_applied(tmp_path):
    preferences = display.display_preferences(tmp_path)
    preferences.setValue("display/mode", "enormous")
    preferences.setValue("display/scale_percent", "not a number")
    preferences.setValue("display/rounding", "sideways")
    preferences.sync()
    assert display.read_display_settings(tmp_path) == {
        "mode": "system", "scale_percent": 100, "rounding": "exact"}


def test_chosen_display_settings_round_trip_through_the_interface_ini(tmp_path):
    stored = display.write_display_settings(mode="fixed", scale_percent=150,
                                            rounding="whole-steps", root=tmp_path)
    assert stored == {"mode": "fixed", "scale_percent": 150, "rounding": "whole-steps"}
    assert display.read_display_settings(tmp_path) == stored
    assert (tmp_path / "Interface.ini").exists()


def test_a_scale_outside_the_offered_choices_is_refused(tmp_path):
    with pytest.raises(ValueError):
        display.write_display_settings(scale_percent=137, root=tmp_path)
    with pytest.raises(ValueError):
        display.write_display_settings(mode="huge", root=tmp_path)
    with pytest.raises(ValueError):
        display.write_display_settings(rounding="sideways", root=tmp_path)
    assert display.read_display_settings(tmp_path)["scale_percent"] == 100


def test_a_scale_that_would_not_fit_the_screen_is_never_offered():
    """The window has a 1000x680 minimum in logical pixels, so scaling can exile it."""
    assert display.max_scale_percent(QSize(1920, 1080)) == 158
    assert display.available_scale_choices(QSize(1920, 1080)) == (100, 125, 150)
    assert display.max_scale_percent(QSize(3840, 2160)) == 317
    assert display.available_scale_choices(QSize(3840, 2160)) == display.SCALE_CHOICES
    assert display.max_scale_percent(QSize(0, 0)) == 100
    assert display.available_scale_choices(QSize(800, 600)) == (100,)
    assert display.max_scale_percent((1920, 1080)) == 158


def test_applying_settings_after_the_application_started_changes_nothing(qtbot, tmp_path,
                                                                        monkeypatch):
    """pytest-qt owns the QApplication; no policy and no environment may be touched."""
    display.write_display_settings(mode="fixed", scale_percent=200, root=tmp_path)
    monkeypatch.delenv("QT_SCALE_FACTOR", raising=False)
    state = display.apply_display_settings(root=tmp_path)
    assert state["applied"] is False
    assert state["reason"] == "already_started"
    assert state["scale_percent"] == 200
    assert "QT_SCALE_FACTOR" not in __import__("os").environ


def test_a_fixed_preference_sets_the_scale_factor_before_the_application_is_built(tmp_path):
    display.write_display_settings(mode="fixed", scale_percent=150, rounding="whole-steps",
                                   root=tmp_path)
    state = probe(tmp_path)
    assert state["applied"] is True
    assert state["reason"] == "preference"
    assert state["QT_SCALE_FACTOR"] == "1.5"
    assert "Round" in state["policy"] and "PreferFloor" not in state["policy"]


def test_following_the_system_setting_leaves_qt_scaling_alone(tmp_path):
    state = probe(tmp_path)
    assert state["applied"] is True
    assert state["reason"] == "system"
    assert state["QT_SCALE_FACTOR"] is None
    assert "PassThrough" in state["policy"]


def test_an_operator_environment_variable_always_wins_over_the_saved_preference(tmp_path):
    display.write_display_settings(mode="fixed", scale_percent=300, root=tmp_path)
    state = probe(tmp_path, environment={"QT_SCALE_FACTOR": "1.25"})
    assert state["reason"] == "environment"
    assert state["QT_SCALE_FACTOR"] == "1.25"


def test_the_command_line_override_is_the_documented_recovery_path(tmp_path):
    display.write_display_settings(mode="fixed", scale_percent=300, root=tmp_path)
    state = probe(tmp_path, override=100)
    assert state["reason"] == "override"
    assert state["QT_SCALE_FACTOR"] is None, "100% must not magnify anything"
    assert "--display-scale 100" in display.RECOVERY_NOTICE


def test_the_resolved_state_is_described_in_words_the_user_can_check():
    described = display.describe_display_state(
        {"reason": "preference", "scale_percent": 150, "rounding": "exact"})
    assert described == "Active display scale: 150% (from your saved preference) · rounding: exact"
    assert "system" in display.describe_display_state(
        {"reason": "system", "scale_percent": 100, "rounding": "exact"})


def test_a_scale_larger_than_the_screen_can_show_is_capped_not_obeyed(tmp_path):
    """A saved 300% on a 1080p laptop would push the window off the screen."""
    display.write_display_settings(mode="fixed", scale_percent=300, root=tmp_path)
    state = probe(tmp_path, screen="1920x1080")
    assert state["reason"] == "capped"
    assert state["scale_percent"] == 150
    assert state["QT_SCALE_FACTOR"] == "1.5"
    assert max(display.available_scale_choices(QSize(1920, 1080))) == 150


def test_a_running_application_short_circuits_before_any_capping(qtbot, tmp_path):
    display.write_display_settings(mode="fixed", scale_percent=300, root=tmp_path)
    environment = {}
    state = display.apply_display_settings(root=tmp_path, environment=environment,
                                           available_size=QSize(1920, 1080))
    assert state["reason"] == "already_started"
    assert environment == {}


def test_text_scaling_keeps_the_eighty_to_one_hundred_fifty_contract():
    assert display.validate_text_scale(80) == 80
    assert display.validate_text_scale(150) == 150
    for refused in (79, 151, 0, -100):
        with pytest.raises(ValueError):
            display.validate_text_scale(refused)
    with pytest.raises(ValueError):
        display.validate_text_scale(True)
    with pytest.raises(ValueError):
        display.validate_text_scale(100.0)
    assert display.clamp_text_scale(200) == 150
    assert display.clamp_text_scale(10) == 80
    assert display.clamp_text_scale("125") == 125
    assert display.clamp_text_scale(None) == 100


def test_the_existing_stylesheet_rescaler_holds_the_same_boundaries():
    """scaled_style had no coverage; its range is the contract display.py mirrors."""
    assert scaled_style(100) == STYLE
    assert "font-size: 19.5px" in scaled_style(150)
    assert "font-size: 10.4px" in scaled_style(80)
    for refused in (79, 151, True, 100.0, "100"):
        with pytest.raises(ValueError):
            scaled_style(refused)


def test_scaling_is_a_viewing_preference_and_says_so_where_the_user_reads_it():
    assert "not affected" in display.RESTART_NOTICE
    assert "same on every computer" in display.EXPORT_NOTICE
    assert "cannot fit the whole window" in display.HIDDEN_CHOICES_NOTICE
    assert set(display.ROUNDING) == set(display.ROUNDING_LABELS)
