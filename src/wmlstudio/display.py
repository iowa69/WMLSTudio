"""Screen scaling the user can adjust, applied before the application object exists.

Two different levers are kept apart because they do different things:

* **Text size** (80-150 %) rewrites the font sizes in the stylesheet. It changes
  text only, and it applies immediately.
* **Whole-interface scale** (100-300 %) sets Qt's own scale factor, so padding,
  row heights, icons and splitter handles grow with the text. Qt reads it when
  the QApplication is created, so it only takes effect after a restart.

Neither changes a distance, a threshold, an export or a rendered evidence image:
scaling is a viewing preference, never part of the evidence.

This module deliberately imports no QtWidgets, no theme and no widget code, so it
can be imported and called while no QApplication exists yet.
"""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QSettings, Qt
from PySide6.QtGui import QGuiApplication

from wmlstudio.paths import data_root

# Whole-interface scales offered in Settings. 100 % means "follow the system".
SCALE_CHOICES = (100, 125, 150, 175, 200, 250, 300)

# Text sizes offered in Settings, matching the existing interface scale combo.
TEXT_SCALE_CHOICES = (80, 90, 100, 110, 125, 150)
TEXT_SCALE_MINIMUM = 80
TEXT_SCALE_MAXIMUM = 150

ROUNDING = {
    "exact": Qt.HighDpiScaleFactorRoundingPolicy.PassThrough,
    "whole-steps": Qt.HighDpiScaleFactorRoundingPolicy.Round,
    "prefer-smaller": Qt.HighDpiScaleFactorRoundingPolicy.RoundPreferFloor,
}

ROUNDING_LABELS = {
    "exact": "Exact (use the size Windows reports)",
    "whole-steps": "Whole steps only (100 %, 200 %)",
    "prefer-smaller": "Prefer the smaller whole step",
}

# BaseWindow.setMinimumSize, in logical pixels. A fixed scale multiplies it, so a
# scale the screen cannot fit would leave the window partly off-screen.
MINIMUM_WINDOW_WIDTH = 1000
MINIMUM_WINDOW_HEIGHT = 680
# The same ceiling --window-size accepts, so the menu and the command line agree.
MAXIMUM_WINDOW_WIDTH = 7680
MAXIMUM_WINDOW_HEIGHT = 4320

# Window sizes offered by name, smallest first. These are window sizes, not screen
# resolutions: the screen keeps whatever resolution the operating system gave it.
WINDOW_SIZE_CHOICES = ((1000, 680), (1280, 800), (1366, 768), (1440, 900),
                       (1600, 1000), (1920, 1080), (2560, 1440))

DEFAULTS = {"mode": "system", "scale_percent": 100, "rounding": "exact"}

RESTART_NOTICE = ("Close and reopen WMLSTudio to use this size. "
                  "Your project and results are not affected.")
HIDDEN_CHOICES_NOTICE = ("Larger sizes are hidden because this screen cannot fit the whole "
                         "window at that scale.")
RECOVERY_NOTICE = ("If the window ever becomes too large to use, start WMLSTudio once with "
                   "--display-scale 100.")
EXPORT_NOTICE = ("Exported and reported images keep a fixed size so the picture looks the same "
                 "on every computer.")
# The words a person actually searches with when they cannot read the screen.
DISPLAY_SEARCH_TERMS = ("window size", "screen resolution", "text too small", "zoom",
                        "make everything bigger", "high DPI")
DISPLAY_TITLE = "Window size, text size and screen resolution"
DISPLAY_INTRO = ("Text too small, or the window too big for your screen? Everything that makes "
                 "WMLSTudio bigger, smaller or easier to read is on this page. None of it "
                 "changes a distance, a threshold or an exported result.")
WINDOW_SIZE_NOTICE = ("This resizes the WMLSTudio window only. Your screen keeps the resolution "
                      "your operating system gave it.")


def display_preferences(root=None) -> QSettings:
    """Display keys live beside the data root: they are read before a project opens."""
    base = Path(root) if root is not None else data_root()
    base.mkdir(parents=True, exist_ok=True)
    return QSettings(str(base / "Interface.ini"), QSettings.Format.IniFormat)


def _as_int(value, fallback):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return fallback


def read_display_settings(root=None) -> dict:
    """{'mode', 'scale_percent', 'rounding'}, falling back to safe values."""
    try:
        preferences = display_preferences(root)
        mode = str(preferences.value("display/mode", DEFAULTS["mode"]) or DEFAULTS["mode"])
        percent = _as_int(preferences.value("display/scale_percent", DEFAULTS["scale_percent"]),
                          DEFAULTS["scale_percent"])
        rounding = str(preferences.value("display/rounding", DEFAULTS["rounding"])
                       or DEFAULTS["rounding"])
    except (OSError, ValueError):
        return dict(DEFAULTS)
    if mode not in {"system", "fixed"}:
        mode = DEFAULTS["mode"]
    if percent not in SCALE_CHOICES:
        percent = DEFAULTS["scale_percent"]
    if rounding not in ROUNDING:
        rounding = DEFAULTS["rounding"]
    return {"mode": mode, "scale_percent": percent, "rounding": rounding}


def write_display_settings(*, mode=None, scale_percent=None, rounding=None, root=None) -> dict:
    """Persist the chosen display keys and return the stored settings."""
    current = read_display_settings(root)
    if mode is not None:
        if mode not in {"system", "fixed"}:
            raise ValueError("Display mode must be 'system' or 'fixed'.")
        current["mode"] = mode
    if scale_percent is not None:
        percent = _as_int(scale_percent, 0)
        if percent not in SCALE_CHOICES:
            raise ValueError(f"Display scale must be one of {', '.join(map(str, SCALE_CHOICES))}.")
        current["scale_percent"] = percent
    if rounding is not None:
        if rounding not in ROUNDING:
            raise ValueError("Rounding must be 'exact', 'whole-steps' or 'prefer-smaller'.")
        current["rounding"] = rounding
    preferences = display_preferences(root)
    preferences.setValue("display/mode", current["mode"])
    preferences.setValue("display/scale_percent", current["scale_percent"])
    preferences.setValue("display/rounding", current["rounding"])
    preferences.sync()
    return current


def validate_text_scale(percent) -> int:
    """Text sizes outside 80-150 % are refused, matching interface_settings.scaled_style."""
    if isinstance(percent, bool) or not isinstance(percent, int):
        raise ValueError("Interface text size must be a whole percentage.")
    if not TEXT_SCALE_MINIMUM <= percent <= TEXT_SCALE_MAXIMUM:
        raise ValueError(f"Interface text size must be between {TEXT_SCALE_MINIMUM}% "
                         f"and {TEXT_SCALE_MAXIMUM}%")
    return percent


def clamp_text_scale(percent) -> int:
    """The same range, applied to a value the user nudged past the end."""
    return max(TEXT_SCALE_MINIMUM, min(TEXT_SCALE_MAXIMUM, _as_int(percent, 100)))


def max_scale_percent(available_size) -> int:
    """The largest whole-interface scale whose minimum window still fits this screen."""
    width, height = _sides(available_size)
    if width <= 0 or height <= 0:
        return 100
    fit = min(width / MINIMUM_WINDOW_WIDTH, height / MINIMUM_WINDOW_HEIGHT)
    return max(100, int(fit * 100))


def available_scale_choices(available_size) -> tuple[int, ...]:
    """Offer only the scales this screen can actually show the whole window at."""
    ceiling = max_scale_percent(available_size)
    return tuple(percent for percent in SCALE_CHOICES if percent <= ceiling) or (100,)


def _sides(available_size):
    """(width, height) from a QSize, a QRect-like size or a plain pair."""
    width = getattr(available_size, "width", None)
    height = getattr(available_size, "height", None)
    width = width() if callable(width) else available_size[0]
    height = height() if callable(height) else available_size[1]
    return int(width), int(height)


def clamp_window_size(width, height, available_size=None):
    """A window size inside the supported range, and inside this screen when known."""
    width = max(MINIMUM_WINDOW_WIDTH, min(MAXIMUM_WINDOW_WIDTH, _as_int(width, 0)))
    height = max(MINIMUM_WINDOW_HEIGHT, min(MAXIMUM_WINDOW_HEIGHT, _as_int(height, 0)))
    if available_size is not None:
        screen_width, screen_height = _sides(available_size)
        if screen_width >= MINIMUM_WINDOW_WIDTH:
            width = min(width, screen_width)
        if screen_height >= MINIMUM_WINDOW_HEIGHT:
            height = min(height, screen_height)
    return width, height


def available_window_sizes(available_size=None):
    """Offer only the named sizes this screen can actually show in full."""
    if available_size is None:
        return WINDOW_SIZE_CHOICES
    screen_width, screen_height = _sides(available_size)
    fitting = tuple((width, height) for width, height in WINDOW_SIZE_CHOICES
                    if width <= screen_width and height <= screen_height)
    return fitting or ((MINIMUM_WINDOW_WIDTH, MINIMUM_WINDOW_HEIGHT),)


def describe_window_size(width, height) -> str:
    """How one size is written wherever it is offered or reported."""
    return f"{int(width)} × {int(height)}"


def window_size_token(width, height) -> str:
    """A size as one plain string, for a combo box entry or a stored preference.

    A tuple cannot be used there: Qt wraps arbitrary Python objects, and looking
    one up again with `findData` does not reliably match an equal tuple.
    """
    return f"{int(width)}x{int(height)}"


def parse_window_size(token):
    """'1280x800' back to (1280, 800), clamped. Anything else is None."""
    try:
        width, height = (int(part) for part in str(token).lower().split("x"))
    except (AttributeError, TypeError, ValueError):
        return None
    return clamp_window_size(width, height)


def read_window_size(root=None):
    """The window size the user last chose, or None when they never chose one."""
    try:
        preferences = display_preferences(root)
        width = _as_int(preferences.value("display/window_width", 0), 0)
        height = _as_int(preferences.value("display/window_height", 0), 0)
    except (OSError, ValueError):
        return None
    if width < MINIMUM_WINDOW_WIDTH or height < MINIMUM_WINDOW_HEIGHT:
        return None
    return clamp_window_size(width, height)


def write_window_size(width, height, *, root=None):
    """Persist a chosen window size beside the other display keys."""
    width, height = clamp_window_size(width, height)
    preferences = display_preferences(root)
    preferences.setValue("display/window_width", width)
    preferences.setValue("display/window_height", height)
    preferences.sync()
    return width, height


def clear_window_size(root=None) -> None:
    """Forget the chosen size, so the window sizes itself to the screen again."""
    preferences = display_preferences(root)
    preferences.remove("display/window_width")
    preferences.remove("display/window_height")
    preferences.sync()


def apply_display_settings(override_percent=None, *, root=None, environment=None,
                           available_size=None) -> dict:
    """Set Qt's scaling policy before QApplication exists. Reports what was applied.

    Returns {'applied', 'reason', 'mode', 'scale_percent', 'rounding', 'scale_factor'}.
    'reason' is one of already_started, environment, override, preference or system.
    """
    settings = read_display_settings(root)
    state = {"applied": False, "reason": "already_started", "scale_factor": None, **settings}
    if QCoreApplication.instance() is not None:
        # Qt has already read its scaling policy; changing it now would lie about
        # what the running interface is doing. Tests and embedders land here.
        return state
    environment = os.environ if environment is None else environment
    QGuiApplication.setHighDpiScaleFactorRoundingPolicy(ROUNDING[settings["rounding"]])
    state["applied"] = True
    if environment.get("QT_SCALE_FACTOR"):
        # An operator who set the variable outright always wins.
        state["reason"] = "environment"
        state["scale_factor"] = str(environment["QT_SCALE_FACTOR"])
        return state
    percent = None
    if override_percent is not None:
        percent = _as_int(override_percent, 0)
        if percent not in SCALE_CHOICES:
            raise ValueError(f"Display scale must be one of {', '.join(map(str, SCALE_CHOICES))}.")
        state["reason"] = "override"
        state["mode"] = "fixed" if percent != 100 else "system"
        state["scale_percent"] = percent
    elif settings["mode"] == "fixed":
        percent = settings["scale_percent"]
        state["reason"] = "preference"
    if percent is not None and available_size is not None:
        ceiling = max_scale_percent(available_size)
        if percent > ceiling:
            percent = max(c for c in SCALE_CHOICES if c <= ceiling)
            state["scale_percent"] = percent
            state["reason"] = "capped"
    if percent and percent != 100:
        state["scale_factor"] = f"{percent / 100:g}"
        environment["QT_SCALE_FACTOR"] = state["scale_factor"]
        return state
    # Nothing to magnify: Qt keeps following whatever the screen reports.
    if state["reason"] in {"already_started", "preference"}:
        state["reason"] = "system"
    state["mode"] = "system"
    state["scale_percent"] = 100
    return state


def describe_display_state(state) -> str:
    """One honest line for the Settings panel about what is actually in force."""
    sources = {
        "already_started": "already running, not changed",
        "environment": "from QT_SCALE_FACTOR in your environment",
        "override": "from the --display-scale option",
        "capped": "reduced to fit this screen",
        "preference": "from your saved preference",
        "system": "following your system setting",
    }
    reason = sources.get(state.get("reason", ""), "unknown source")
    percent = state.get("scale_percent", 100)
    rounding = state.get("rounding", DEFAULTS["rounding"])
    return f"Active display scale: {percent}% ({reason}) · rounding: {rounding}"
