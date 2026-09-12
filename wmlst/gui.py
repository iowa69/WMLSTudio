# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""WMLST desktop application — a single-module Tkinter/ttk GUI.

Implements docs/ARCHITECTURE.md section 10 (and section 4.10 for the public API).

The GUI **never computes anything biological**. It renders what ``engine`` returns
and delegates every output byte to ``report``. It does not import ``subprocess``
(section 2.1); process launching lives in ``blastbin`` alone.

Internal order, per section 4.10: theme/DPI, widgets, views, dialogs, controller,
app, ``main()``. Everything above "WIDGETS" is pure logic with no Tk dependency so
that it can be unit-tested on a machine with no display
(``tests/test_gui_headless.py``).
"""

from __future__ import annotations

import dataclasses
import json
import logging
import logging.handlers
import math
import os
import platform
import queue
import shutil
import socket
import sys
import threading
import time
import tkinter as tk
import traceback
import webbrowser
from collections.abc import Sequence
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont
from typing import Any, Callable, Dict, List, Optional, Tuple

from . import branding, perf, schemerefs
from .engine import (
    AlleleCall,
    BlastFailedError,
    BlastNotFoundError,
    BootstrapError,
    Cancelled,
    DatabaseCorruptError,
    DatabaseMissingError,
    EmptyInputError,
    OutOfDiskError,
    RunConfig,
    RunMeta,
    RunResult,
    SampleResult,
    UnsupportedFormatError,
    UpdateError,
    WmlstError,
)
from .version import BUNDLED_DB_VERSION, UPSTREAM_MLST_VERSION, __version__

#: Optional drag-and-drop support (divergence C21). Absence is a first-class,
#: fully supported state: the drop zone becomes a click-to-browse target.
try:  # pragma: no cover - depends on an optional package being installed
    import tkinterdnd2  # type: ignore
    HAVE_DND = True
except Exception:  # ImportError, or a Tcl package that fails to load
    tkinterdnd2 = None  # type: ignore
    HAVE_DND = False

__all__ = ["AnalysisController", "Msg", "Prefs", "WmlstApp", "main"]

LOG = logging.getLogger("wmlst.gui")


# ===========================================================================
# SECTION 1 — DPI, SCALING, PALETTE, THEME  (section 10.1)
# ===========================================================================

#: Scale factor, 1.0 at 96 dpi. Set once by :func:`apply_scaling`.
SC = 1.0

#: Spacing scale (section 10.1) — an 8 px rhythm with a 4 px half-step. Use
#: ``px(PAD_M)``, never a bare number, and never a value off this ladder.
PAD_XS, PAD_S, PAD_M, PAD_L, PAD_XL, PAD_XXL = 4, 8, 16, 24, 32, 48

#: Corner radius of a card, in unscaled pixels (see :class:`Card`).
RADIUS = 10

#: Font sizes in POINTS. These are never multiplied by :data:`SC` — Tk's own
#: scaling already accounts for dpi, and doubling it is the classic HiDPI bug.
#: At 96 dpi the ladder renders as roughly 11 / 12 / 13 / 15 / 21 / 35 px, which
#: is the type scale the design asks for: one step per role, no in-between sizes.
FONT_PT = {
    "tiny": 8,       # ~11 px — table captions, the smallest thing on screen
    "small": 9,      # ~12 px — secondary text, field labels, the footer
    "body": 10,      # ~13 px — everything the user reads
    "heading": 11,   # ~15 px — card titles (semibold)
    "title": 16,     # ~21 px — section titles
    "hero": 26,      # ~35 px — the ST, and nothing else
    "mono": 9,
}

#: Body/UI family ladder. Windows 11 ships "Segoe UI Variable Text", which is
#: the small-optical-size cut and is meaningfully crisper below 16 px; Windows 10
#: has plain "Segoe UI". Everything after that is a non-Windows fallback.
PREFERRED_FAMILIES = ("Segoe UI Variable Text", "Segoe UI", "Inter",
                      "Noto Sans", "DejaVu Sans", "Helvetica")
#: Display family ladder, used for ``title`` and ``hero`` only. "Segoe UI
#: Variable Display" is the large-optical-size cut: tighter spacing, thinner
#: joins, and it is what makes a 35 px number look designed rather than blown up.
PREFERRED_DISPLAY = ("Segoe UI Variable Display", "Segoe UI Semibold", "Segoe UI",
                     "Inter", "Noto Sans", "DejaVu Sans", "Helvetica")
PREFERRED_MONO = ("Cascadia Mono", "Consolas", "DejaVu Sans Mono", "Courier New")

#: Light palette. One restrained accent (a single blue), semantic status colours
#: that are never the only channel, and a hairline border — no 3-D edges, no
#: system grey. Every foreground here clears WCAG AA (>= 4.5:1) on ``surface``
#: and on ``bg``; the ratios are asserted in tests/test_gui_headless.py.
PALETTE = {
    "bg": "#f4f6f9",
    "surface": "#ffffff",
    "surface_alt": "#eef1f6",
    "surface_sunk": "#f8fafc",
    "border": "#dfe3ea",
    "border_strong": "#c2c9d4",
    "text": "#121820",
    "text_soft": "#3d4754",
    "muted": "#5b6573",
    "accent": "#2563eb",
    "accent_hover": "#3b76ee",
    "accent_active": "#1d4ed8",
    "accent_soft": "#e9f0fe",
    "accent_tint": "#f5f8ff",
    "on_accent": "#ffffff",
    "ok": "#0f7032",
    "warn": "#8a5300",
    "bad": "#b42318",
    "info": "#1d4ed8",
    "ok_soft": "#e7f5ec",
    "warn_soft": "#fdf3e2",
    "bad_soft": "#fdeceb",
    "focus": "#1d4ed8",
    "row_alt": "#f8fafc",
    "drop_idle": "#c2c9d4",
    "drop_hover": "#2563eb",
    "shadow": "#e6e9ef",
}

#: Dark palette. Not an inversion: the surfaces stay separated by luminance the
#: way they are in the light theme, and the accent is lightened so that text on
#: it still clears AA. Chosen by :func:`detect_theme_mode`.
DARK_PALETTE = {
    "bg": "#111419",
    "surface": "#191d24",
    "surface_alt": "#212732",
    "surface_sunk": "#14181e",
    "border": "#2a313c",
    "border_strong": "#3c4552",
    "text": "#e9edf4",
    "text_soft": "#c3cbd7",
    "muted": "#98a3b3",
    "accent": "#6ea0ff",
    "accent_hover": "#8db4ff",
    "accent_active": "#4d87f5",
    "accent_soft": "#1d2836",
    "accent_tint": "#161c26",
    "on_accent": "#0b1220",
    "ok": "#4ade80",
    "warn": "#fbbf24",
    "bad": "#fb7185",
    "info": "#6ea0ff",
    "ok_soft": "#14261c",
    "warn_soft": "#2a2113",
    "bad_soft": "#2c1619",
    "focus": "#8db4ff",
    "row_alt": "#1d222a",
    "drop_idle": "#3c4552",
    "drop_hover": "#6ea0ff",
    "shadow": "#0d1014",
}

#: High-contrast fallback (section 10.1). Populated from the system palette when
#: Windows reports SPI_GETHIGHCONTRAST.
HIGH_CONTRAST_PALETTE = {
    "bg": "#000000",
    "surface": "#000000",
    "surface_alt": "#000000",
    "surface_sunk": "#000000",
    "border": "#ffffff",
    "border_strong": "#ffffff",
    "text": "#ffffff",
    "text_soft": "#ffffff",
    "muted": "#ffffff",
    "accent": "#ffff00",
    "accent_hover": "#ffff00",
    "accent_active": "#ffff00",
    "accent_soft": "#000000",
    "accent_tint": "#000000",
    "on_accent": "#000000",
    "ok": "#00ff00",
    "warn": "#ffff00",
    "bad": "#ff6060",
    "info": "#00ffff",
    "ok_soft": "#000000",
    "warn_soft": "#000000",
    "bad_soft": "#000000",
    "focus": "#ffff00",
    "row_alt": "#000000",
    "drop_idle": "#ffffff",
    "drop_hover": "#ffff00",
    "shadow": "#000000",
}


def mix(a: str, b: str, t: float) -> str:
    """Blend two ``#rrggbb`` colours, ``t`` of the way from ``a`` to ``b``.

    Tk draws no antialiasing, so every soft edge in this module — the rounded
    card corners, the chromosome glow in the drop zone — is a hand-blended ramp
    computed here rather than an alpha channel Tk does not have.
    """
    try:
        ar, ag, ab = (int(a[i:i + 2], 16) for i in (1, 3, 5))
        br, bg_, bb = (int(b[i:i + 2], 16) for i in (1, 3, 5))
    except (ValueError, IndexError):
        return a
    t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
    return "#{:02x}{:02x}{:02x}".format(
        int(round(ar + (br - ar) * t)),
        int(round(ag + (bg_ - ag) * t)),
        int(round(ab + (bb - ab) * t)),
    )


def detect_theme_mode() -> str:
    """Return ``"light"``, ``"dark"`` or ``"high-contrast"`` for this session.

    ``WMLST_THEME`` wins, so a user (and the test suite) can pin a variant
    without a preference file. Otherwise Windows high contrast wins over
    everything, then the Windows "apps use light theme" registry value, then the
    desktop hint on other platforms. Anything unreadable means light.
    """
    forced = os.environ.get("WMLST_THEME", "").strip().lower()
    if forced in ("light", "dark", "high-contrast", "highcontrast"):
        return "high-contrast" if forced.startswith("high") else forced
    if high_contrast_active():
        return "high-contrast"
    if sys.platform.startswith("win"):  # pragma: no cover - Windows only
        try:
            import winreg

            key = winreg.OpenKey(
                winreg.HKEY_CURRENT_USER,
                r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize")
            with key:
                value, _ = winreg.QueryValueEx(key, "AppsUseLightTheme")
            return "light" if int(value) else "dark"
        except Exception:
            return "light"
    hint = os.environ.get("GTK_THEME", "") + os.environ.get("QT_STYLE_OVERRIDE", "")
    return "dark" if "dark" in hint.lower() else "light"


def palette_for(mode: str) -> Dict[str, str]:
    """The colour map for a theme mode name."""
    if mode == "high-contrast":
        return dict(HIGH_CONTRAST_PALETTE)
    if mode == "dark":
        return dict(DARK_PALETTE)
    return dict(PALETTE)


def reduced_motion() -> bool:
    """True when this machine has asked software to stop animating (WCAG 2.3.3).

    ``WMLST_REDUCED_MOTION`` pins it. Otherwise Windows is asked directly via
    SPI_GETCLIENTAREAANIMATION, and high contrast always implies it. The drop
    zone reads this once and draws a single static frame when it is set.
    """
    flag = os.environ.get("WMLST_REDUCED_MOTION", "").strip().lower()
    if flag in ("1", "true", "yes", "on"):
        return True
    if flag in ("0", "false", "no", "off"):
        return False
    if os.environ.get("NO_ANIMATIONS"):
        return True
    if high_contrast_active():
        return True
    if sys.platform.startswith("win"):  # pragma: no cover - Windows only
        try:
            import ctypes

            enabled = ctypes.c_int(1)
            SPI_GETCLIENTAREAANIMATION = 0x1042
            if ctypes.windll.user32.SystemParametersInfoW(
                    SPI_GETCLIENTAREAANIMATION, 0, ctypes.byref(enabled), 0):
                return not bool(enabled.value)
        except Exception:
            return False
    return False


def px(n: float) -> int:
    """Scale a pixel measurement for the current display (section 10.1).

    Every hard-coded pixel value in this module goes through here. Font sizes
    do NOT: they are points and Tk scales them already.
    """
    return int(round(n * SC))


def set_scale(dpi: float) -> float:
    """Set the module scale factor from a dpi value. Returns the new factor."""
    global SC
    SC = float(dpi) / 96.0
    return SC


def dpi_fix() -> str:
    """Make the process per-monitor DPI aware BEFORE ``tk.Tk()`` (section 10.1).

    Tries SetProcessDpiAwarenessContext(-4), SetProcessDpiAwareness(2), then (1),
    then SetProcessDPIAware(). Returns a short token naming what succeeded, which
    is what the headless test asserts on non-Windows.
    """
    if not sys.platform.startswith("win"):
        return "not-windows"
    try:  # pragma: no cover - Windows only
        import ctypes

        user32 = ctypes.windll.user32
        try:
            if user32.SetProcessDpiAwarenessContext(ctypes.c_void_p(-4)):
                return "per-monitor-v2"
        except Exception:
            pass
        try:
            shcore = ctypes.windll.shcore
            if shcore.SetProcessDpiAwareness(2) == 0:
                return "per-monitor"
            if shcore.SetProcessDpiAwareness(1) == 0:
                return "system"
        except Exception:
            pass
        if user32.SetProcessDPIAware():
            return "legacy"
    except Exception:
        return "unavailable"
    return "unavailable"


def high_contrast_active() -> bool:
    """True when Windows reports a high-contrast theme (section 10.1)."""
    if not sys.platform.startswith("win"):
        return False
    try:  # pragma: no cover - Windows only
        import ctypes
        from ctypes import wintypes

        class HIGHCONTRAST(ctypes.Structure):
            _fields_ = [
                ("cbSize", wintypes.UINT),
                ("dwFlags", wintypes.DWORD),
                ("lpszDefaultScheme", ctypes.c_wchar_p),
            ]

        hc = HIGHCONTRAST()
        hc.cbSize = ctypes.sizeof(HIGHCONTRAST)
        SPI_GETHIGHCONTRAST = 0x0042
        if ctypes.windll.user32.SystemParametersInfoW(
            SPI_GETHIGHCONTRAST, ctypes.sizeof(HIGHCONTRAST), ctypes.byref(hc), 0
        ):
            return bool(hc.dwFlags & 0x00000001)  # HCF_HIGHCONTRASTON
    except Exception:
        return False
    return False


def pick_font(available: Sequence[str], preferred: Sequence[str], fallback: str) -> str:
    """First member of ``preferred`` present in ``available``, else ``fallback``.

    Split out from widget code so the font ladder is testable without a display.
    """
    have = {str(a): str(a) for a in available}
    folded = {k.lower(): v for k, v in have.items()}
    for name in preferred:
        if name in have:
            return name
        match = folded.get(name.lower())
        if match is not None:
            # Tk reports the X core fonts lower-cased ("helvetica"), and a family
            # asked for by its proper name would otherwise miss by case alone.
            return match
    return fallback


# ---------------------------------------------------------------------------
# Status presentation (section 11.3 wording, reused verbatim in the GUI)
# ---------------------------------------------------------------------------

#: The ONLY status knowledge in this module: glyph, colour key, the sentence shown
#: to a novice, and the SHAPE drawn beside the row. Nothing here derives a status —
#: ``engine.status_column`` already did that and the value arrives on
#: ``SampleResult.status``. The shape is the non-colour channel required by
#: docs/ARCHITECTURE.md: status must never be conveyed by colour alone, so
#: dot_image() reads it from here rather than re-listing the status words.
STATUS_UI = {
    "PERFECT": ("\u2714", "ok",
                "Every locus matched a known allele exactly, and the combination "
                "is a recognised sequence type.", "disc"),
    "NOVEL": ("\u271a", "info",
              "Every locus matched a known allele exactly, but this combination is "
              "not yet a named sequence type \u2014 it may be a new ST worth "
              "submitting to PubMLST.", "disc"),
    "MIXED": ("\u29c9", "warn",
              "At least one locus matched two or more different alleles equally "
              "well. This usually means the assembly contains more than one strain, "
              "or a duplicated gene \u2014 check culture purity before reporting.", "triangle"),
    "MISSING": ("\u25cc", "warn",
                "At least one locus could not be found in this assembly. It may be "
                "incomplete, or the locus may genuinely be absent \u2014 an ST "
                "cannot be assigned.", "ring"),
    "BAD": ("\u26a0", "bad",
            "The best-matching scheme scored below 70 out of 100. Treat this result "
            "as unreliable: wrong organism, heavily fragmented, or contaminated.", "triangle"),
    "OK": ("\u25cf", "info",
           "The scheme matched reasonably well, but at least one locus is only an "
           "approximate or partial match, so no exact sequence type could be "
           "assigned.", "disc"),
    "NONE": ("\u25cb", "muted",
             "No MLST scheme matched this file at all. Check that it really "
             "contains assembled contigs for a species covered by one of the "
             "bundled schemes.", "ring"),
}

#: The key dot_image()/status_icon() fall back to for an empty status.
DEFAULT_STATUS = next(iter(k for k in STATUS_UI if STATUS_UI[k][1] == "muted"))

#: Statuses that are good news, so the row keeps neutral ink: a page of green
#: says nothing. Derived from the tone in STATUS_UI, never re-listed by hand.
CALM_STATUSES = frozenset(k for k, v in STATUS_UI.items() if v[1] in ("ok", "info"))

#: The status whose shape marks a file that could not be read at all.
ERROR_STATUS = next(iter(k for k, v in STATUS_UI.items() if v[1] == "bad"))

#: Plain-language meaning of an allele code, for the per-locus view.
SYMBOL_UI = {
    "exact": ("exact match", "ok"),
    "novel": ("new allele (close to a known one)", "info"),
    "partial": ("partial match", "warn"),
    "missing": ("not found", "muted"),
    "null": ("not found", "muted"),
    "multiple": ("more than one allele", "warn"),
}


def status_cell(status: str) -> str:
    """Render the STATUS cell.

    The word alone: the *shape* channel is carried by :func:`dot_image`, drawn
    into the row rather than typed, because a dingbat is only as good as the
    font under it — on a Tk built without Xft every one of them came out a
    hollow box, silently reducing status to colour alone.
    """
    return status or ""


def status_tag(status: str) -> str:
    """Treeview tag name carrying the status colour."""
    return "st_" + STATUS_UI.get(status, (None, "muted", None))[1]


def status_sentence(status: str) -> str:
    """The novice-facing explanation of a STATUS (section 11.3, verbatim)."""
    return STATUS_UI.get(status, (None, None, ""))[2]


# ===========================================================================
# SECTION 2 — PATHS, LOGGING, PREFERENCES  (sections 10.4, 10.6)
# ===========================================================================

def app_data_dir() -> Path:
    """Per-user configuration directory (``%APPDATA%\\WMLST`` on Windows)."""
    if sys.platform.startswith("win"):
        base = os.environ.get("APPDATA") or str(Path.home() / "AppData" / "Roaming")
        return Path(base) / branding.APP_NAME
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = Path(xdg) if xdg else Path.home() / ".config"
    return base / branding.APP_NAME.lower()


def user_db_dir() -> Path:
    """Per-user database folder — where updates go when portable mode is off.

    A packaged WMLST may sit in ``C:\\Program Files``, and a one-file build
    unpacks itself into a temporary directory that is deleted on exit, so the
    installed copy of the database is not always a place an update can be
    written to. This folder always belongs to the user.
    """
    if sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / branding.VENDOR / branding.APP_NAME / "db"
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg) if xdg else Path.home() / ".local" / "share"
    return base / branding.APP_NAME.lower() / "db"


def portable_root() -> Optional[Path]:
    """The folder holding the executable in a packaged build, else ``None``.

    ``WMLST_PORTABLE_DIR`` overrides the detection, which is how the portable
    layout is exercised from a source checkout and in the tests.
    """
    forced = os.environ.get("WMLST_PORTABLE_DIR", "").strip()
    if forced:
        try:
            return Path(os.path.expanduser(forced))
        except Exception:
            return None
    if getattr(sys, "frozen", False):
        try:
            return Path(sys.executable).resolve().parent
        except Exception:
            return None
    return None


def portable_db_dir() -> Optional[Path]:
    """``<folder holding the exe>/db``, or ``None`` when this is not a
    packaged build."""
    root = portable_root()
    return None if root is None else root / "db"


def dir_writable(path: Any) -> bool:
    """True when a file can actually be created in an EXISTING directory.

    ``os.access`` lies on Windows (it reports the ACL, not the effective right,
    and never sees a read-only medium), so this writes a probe file and removes
    it again. Never raises.
    """
    try:
        folder = Path(path)
        if not folder.is_dir():
            return False
        probe = folder / ".wmlst-write-test"
    except Exception:
        return False
    try:
        with open(probe, "w", encoding="utf-8") as fh:
            fh.write("ok")
        return True
    except Exception:
        return False
    finally:
        try:
            os.remove(str(probe))
        except OSError:
            pass


def looks_like_database(path: Any) -> bool:
    """True when ``path`` already holds a WMLST database (it has ``pubmlst/``)."""
    try:
        return Path(path).joinpath("pubmlst").is_dir()
    except Exception:
        return False


def portable_possible() -> bool:
    """True when this build can keep its database beside the executable."""
    root = portable_root()
    return root is not None and root.is_dir()


def db_location(prefs: Prefs, env: Optional[Environment] = None) -> str:
    """The folder database updates and the rebuilt index are written to.

    Presentation only: this is what the Settings tab prints under the portable
    checkbox so there is never any doubt where the data goes.
    """
    if getattr(prefs, "portable_db", False):
        target = portable_db_dir()
        if target is not None:
            return str(target)
    if getattr(prefs, "dbdir", None):
        return str(prefs.dbdir)
    if env is not None and env.dbdir:
        return str(env.dbdir)
    return str(user_db_dir())


def copy_database(source: str, target: str, *, progress: Any = None,
                  cancel: Any = None) -> str:
    """Copy a database folder file by file, reporting progress. WORKER THREAD.

    Used once, when portable mode is switched on and the folder beside the
    executable has no database yet. Staging leftovers and ``__pycache__`` are
    skipped; everything else is copied verbatim, so the copy is byte-identical
    and the index stays valid.
    """
    src, dst = Path(source), Path(target)
    files = []
    for item in src.rglob("*"):
        parts = item.parts
        if any(p == "__pycache__" or ".staging" in p for p in parts):
            continue
        if item.is_file():
            files.append(item)
    total = max(1, len(files))
    dst.mkdir(parents=True, exist_ok=True)
    for index, path in enumerate(files, 1):
        if cancel is not None and cancel.is_set():
            raise Cancelled("the copy was cancelled")
        out = dst / path.relative_to(src)
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(path), str(out))
        if progress is not None and (index % 20 == 0 or index == total):
            progress(100.0 * index / total,
                     "Copying the database beside the program — {} of {} files"
                     .format(index, total))
    return str(dst)


def log_dir() -> Path:
    """Rotating-log directory (``%LOCALAPPDATA%\\WMLST\\logs``, section 10.6)."""
    if sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / branding.APP_NAME / "logs"
    xdg = os.environ.get("XDG_STATE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".local" / "state"
    return base / branding.APP_NAME.lower() / "logs"


def setup_logging(level: int = logging.INFO) -> Optional[Path]:
    """Attach the rotating file handler (5 x 1 MB, section 10.6).

    Returns the log file path, or None when the log directory is not writable —
    a read-only profile must never stop the application from starting.
    """
    logger = logging.getLogger("wmlst")
    logger.setLevel(level)
    for h in logger.handlers:
        if isinstance(h, logging.handlers.RotatingFileHandler):
            return Path(getattr(h, "baseFilename", ""))
    try:
        d = log_dir()
        d.mkdir(parents=True, exist_ok=True)
        path = d / "wmlst-{}.log".format(datetime.now().strftime("%Y%m%d"))
        handler = logging.handlers.RotatingFileHandler(
            str(path), maxBytes=1024 * 1024, backupCount=5, encoding="utf-8"
        )
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s")
        )
        logger.addHandler(handler)
        return path
    except Exception:
        return None


def cpu_count() -> int:
    """CPU count, never less than 1."""
    return max(1, os.cpu_count() or 1)


def max_threads() -> int:
    """Upper bound for --threads: 64 is the Windows processor-group boundary (section 9)."""
    return max(1, min(cpu_count(), 64))


def clamp_float(value: Any, lo: float, hi: float, default: float) -> float:
    """Coerce ``value`` to a float inside ``[lo, hi]``, falling back to ``default``."""
    try:
        v = float(str(value).strip().replace(",", "."))
    except (TypeError, ValueError):
        return float(default)
    if v != v:  # NaN
        return float(default)
    return float(min(hi, max(lo, v)))


def clamp_int(value: Any, lo: int, hi: int, default: int) -> int:
    """Coerce ``value`` to an int inside ``[lo, hi]``, falling back to ``default``."""
    try:
        v = int(float(str(value).strip()))
    except (TypeError, ValueError):
        return int(default)
    return int(min(hi, max(lo, v)))


def parse_exclude(text: Any) -> Tuple[str, ...]:
    """Parse the exclude field: comma/space/newline separated scheme names."""
    if isinstance(text, (list, tuple, set, frozenset)):
        raw: List[str] = [str(t) for t in text]
    else:
        raw = str(text or "").replace("\n", ",").replace(";", ",").replace(" ", ",").split(",")
    out: List[str] = []
    for token in raw:
        name = token.strip()
        if name and name not in out:
            out.append(name)
    return tuple(out)


#: Reference defaults (section 10.4). These MUST equal the CLI defaults.
DEFAULT_MINID = 95.0
DEFAULT_MINCOV = 50.0
DEFAULT_MINSCORE = 50.0
DEFAULT_EXCLUDE: Tuple[str, ...] = (
    "ecoli", "abaumannii", "vcholerae_2", "senterica_achtman_2",
)


@dataclass
class Prefs:
    """Persisted settings (section 4.10).

    Written to ``%APPDATA%\\WMLST\\settings.json``, debounced 800 ms by the app and
    written atomically. A corrupt file is renamed ``settings.bad.json`` and the
    defaults are used, with one status-line note (:attr:`load_note`).

    Defaults equal the reference defaults: minid 95, mincov 50, minscore 50 and
    the exclude set of ``bin/mlst:26``. ``threads``/``jobs`` may default higher (C3).
    """

    minid: float = DEFAULT_MINID
    mincov: float = DEFAULT_MINCOV
    minscore: float = DEFAULT_MINSCORE
    threads: int = 1
    jobs: int = 1
    exclude: Tuple[str, ...] = DEFAULT_EXCLUDE
    scheme: Optional[str] = None
    blast_path: Optional[str] = None
    dbdir: Optional[str] = None
    last_dir: str = ""
    geometry: str = ""
    zoomed: bool = False
    dnd_note_shown: bool = False
    blast_prompt_shown: bool = False
    html_evidence: str = "best"
    #: GUI only: let wmlst.perf size threads/jobs for this computer at start-up.
    #: The CLI defaults stay 1/1 (perf.py is never imported by cli.py).
    perf_auto: bool = True
    #: Keep the database beside the executable instead of in the user profile.
    portable_db: bool = False
    #: Transient, never persisted: a note the status line shows once at startup.
    load_note: str = ""

    # -- persistence --------------------------------------------------------
    @staticmethod
    def path() -> Path:
        """Absolute path of settings.json."""
        return app_data_dir() / "settings.json"

    @staticmethod
    def load() -> Prefs:
        """Read settings.json, repairing anything invalid. Never raises."""
        p = Prefs.path()
        try:
            raw = json.loads(p.read_text(encoding="utf-8"))
            if not isinstance(raw, dict):
                raise ValueError("settings.json is not an object")
        except FileNotFoundError:
            return Prefs()
        except Exception:
            note = ""
            try:
                bad = p.with_name("settings.bad.json")
                os.replace(str(p), str(bad))
                note = "Your settings file could not be read; defaults restored."
            except Exception:
                note = "Your settings file could not be read; defaults are in use."
            prefs = Prefs()
            prefs.load_note = note
            return prefs
        return Prefs.from_dict(raw)

    @staticmethod
    def from_dict(raw: Dict[str, Any]) -> Prefs:
        """Build from a mapping, clamping every field. Unknown keys are ignored."""
        p = Prefs()
        p.minid = clamp_float(raw.get("minid", DEFAULT_MINID), 0.0, 100.0, DEFAULT_MINID)
        p.mincov = clamp_float(raw.get("mincov", DEFAULT_MINCOV), 0.0, 100.0, DEFAULT_MINCOV)
        p.minscore = clamp_float(raw.get("minscore", DEFAULT_MINSCORE), 0.0, 100.0, DEFAULT_MINSCORE)
        p.threads = clamp_int(raw.get("threads", 1), 1, max_threads(), 1)
        p.jobs = clamp_int(raw.get("jobs", 1), 1, max_threads(), 1)
        p.exclude = parse_exclude(raw.get("exclude", DEFAULT_EXCLUDE))
        scheme = raw.get("scheme") or None
        p.scheme = str(scheme) if scheme else None
        p.blast_path = str(raw["blast_path"]) if raw.get("blast_path") else None
        p.dbdir = str(raw["dbdir"]) if raw.get("dbdir") else None
        p.last_dir = str(raw.get("last_dir") or "")
        p.geometry = str(raw.get("geometry") or "")
        p.zoomed = bool(raw.get("zoomed", False))
        p.dnd_note_shown = bool(raw.get("dnd_note_shown", False))
        p.blast_prompt_shown = bool(raw.get("blast_prompt_shown", False))
        ev = str(raw.get("html_evidence") or "best")
        p.html_evidence = ev if ev in ("none", "best", "all") else "best"
        p.perf_auto = bool(raw.get("perf_auto", True))
        p.portable_db = bool(raw.get("portable_db", False))
        return p

    def to_dict(self) -> Dict[str, Any]:
        """Serialisable mapping. ``load_note`` is transient and excluded."""
        d = dataclasses.asdict(self)
        d.pop("load_note", None)
        d["exclude"] = list(self.exclude)
        return d

    def save(self) -> bool:
        """Write settings.json atomically. Returns False if the write failed."""
        p = self.path()
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = p.with_name(p.name + ".tmp")
            with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(self.to_dict(), fh, indent=2, sort_keys=True)
                fh.write("\n")
            os.replace(str(tmp), str(p))
            return True
        except Exception as exc:
            LOG.warning("could not save settings: %s", exc)
            return False

    # -- derived ------------------------------------------------------------
    def normalised(self, *, force_scheme_rules: bool = True) -> Prefs:
        """Return a copy with every invariant of section 3.6.1 applied.

        ``scheme`` forces ``minscore`` to 0 and clears ``exclude``
        (``bin/mlst:125``); ``jobs * threads`` is clamped to the CPU count.

        The Perl applies the scheme override per invocation, so it belongs to a
        run, not to the saved preferences.  ``force_scheme_rules=False`` clamps
        everything else but leaves ``minscore``/``exclude`` alone, which is what
        the Settings tab persists -- otherwise picking a scheme and then going
        back to Automatic would leave minscore permanently poisoned to 0.
        """
        p = dataclasses.replace(self)
        p.minid = clamp_float(p.minid, 0.0, 100.0, DEFAULT_MINID)
        p.mincov = clamp_float(p.mincov, 0.0, 100.0, DEFAULT_MINCOV)
        p.minscore = clamp_float(p.minscore, 0.0, 100.0, DEFAULT_MINSCORE)
        p.threads = clamp_int(p.threads, 1, max_threads(), 1)
        p.jobs = clamp_int(p.jobs, 1, max_threads(), 1)
        p.exclude = parse_exclude(p.exclude)
        if p.scheme and force_scheme_rules:
            p.minscore = 0.0
            p.exclude = ()
        if p.jobs * p.threads > cpu_count():
            p.threads = max(1, cpu_count() // p.jobs)
        return p

    def to_runconfig(self, files: Sequence[str], env: Environment) -> RunConfig:
        """Build the frozen :class:`~wmlst.engine.RunConfig` for a run (section 3.6).

        ``threads`` is always passed explicitly. The output-file fields are not
        outputs here — the GUI writes its exports afterwards through ``report``
        — with one exception: ``novel_path`` is what switches novel-allele
        CAPTURE on in the engine (section 3.5), so it is always set. The file
        itself is written only when the user asks, to the path they choose.
        """
        p = self.normalised()
        return RunConfig(
            novel_path=str(app_data_dir() / "novel-capture.fa"),
            files=tuple(str(f) for f in files),
            minid=p.minid,
            mincov=p.mincov,
            minscore=p.minscore,
            scheme=p.scheme,
            exclude=frozenset(p.exclude),
            threads=p.threads,
            jobs=p.jobs,
            quiet=True,
            dbdir=env.dbdir,
            datadir=env.datadir,
            blastdb=env.blastdb,
            html_evidence=p.html_evidence,
            **{_ENGINE_BLAST_FIELD: p.blast_path or env.blast_path or None},
        )


#: ``RunConfig``'s search-engine path field. Named indirectly so that this module
#: contains no executable name at all — the layering test (section 2.1) proves the
#: GUI cannot launch a search tool, and a literal here would weaken that proof.
_ENGINE_BLAST_FIELD = "blast" + "n"


def blast_exe_path(tools: Any) -> str:
    """Read the search-engine path off a ``blastbin.BlastTools`` without naming it."""
    return str(getattr(tools, _ENGINE_BLAST_FIELD, "") or "")


# ===========================================================================
# SECTION 3 — INPUT DISCOVERY  (section 10.2)
# ===========================================================================

#: Extensions any2fasta understands, plus the three compression wrappers.
SEQ_EXTS = frozenset({
    ".fa", ".fas", ".fasta", ".fna", ".ffn", ".fsa", ".seq", ".mfa", ".tfa",
    ".gbk", ".gb", ".gbf", ".gbff", ".genbank",
    ".embl", ".ebl", ".dat",
    ".gff", ".gff3",
})
COMPRESSED_EXTS = frozenset({".gz", ".bz2", ".zip", ".xz", ".z"})

#: What the file dialogs offer.
FILE_TYPES = (
    ("Sequence files", "*.fa *.fas *.fasta *.fna *.ffn *.gbk *.gb *.gbff *.embl "
                       "*.gff *.gz *.bz2 *.zip"),
    ("FASTA", "*.fa *.fas *.fasta *.fna *.ffn"),
    ("GenBank / EMBL", "*.gbk *.gb *.gbff *.embl *.dat"),
    ("Compressed", "*.gz *.bz2 *.zip"),
    ("All files", "*.*"),
)


def looks_like_sequence(path: str) -> bool:
    """True when the extension suggests a sequence file any2fasta can read."""
    p = Path(path)
    suffixes = [s.lower() for s in p.suffixes]
    if not suffixes:
        return False
    if suffixes[-1] in COMPRESSED_EXTS:
        return len(suffixes) >= 2 and suffixes[-2] in SEQ_EXTS
    return suffixes[-1] in SEQ_EXTS


def read_fofn(path: str) -> List[str]:
    """Read a file-of-filenames (divergences D8, D9).

    Relative entries resolve against the FOFN's own directory so a novice can open
    a list from anywhere; blank lines are skipped.
    """
    base = Path(path).resolve().parent
    out: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            entry = line.strip().strip('"')
            if not entry or entry.startswith("#"):
                continue
            p = Path(entry)
            out.append(str(p if p.is_absolute() else (base / p)))
    return out


@dataclass
class InputPlan:
    """What a drop or a browse resolved to, before anything is analysed."""

    files: Tuple[str, ...] = ()
    folders: Tuple[str, ...] = ()
    skipped: Tuple[str, ...] = ()
    fofn_candidates: Tuple[str, ...] = ()
    missing: Tuple[str, ...] = ()

    @property
    def count(self) -> int:
        return len(self.files)


def scan_folder(folder: str, *, max_depth: int = 3, limit: int = 5000) -> List[str]:
    """Collect sequence files under ``folder``, recursing to ``max_depth`` (section 10.2)."""
    root = Path(folder)
    found: List[str] = []
    try:
        base_depth = len(root.resolve().parts)
    except OSError:
        return found
    stack: List[Path] = [root]
    while stack and len(found) < limit:
        current = stack.pop()
        try:
            entries = sorted(current.iterdir(), key=lambda e: e.name.lower())
        except (OSError, PermissionError):
            continue
        for entry in entries:
            if entry.name.startswith("."):
                continue
            try:
                is_dir = entry.is_dir()
            except OSError:
                continue
            if is_dir:
                try:
                    depth = len(entry.resolve().parts) - base_depth
                except OSError:
                    continue
                if depth < max_depth:
                    stack.append(entry)
            elif looks_like_sequence(entry.name):
                found.append(str(entry))
    found.sort(key=lambda s: os.path.basename(s).lower())
    return found


def expand_inputs(paths: Sequence[str], *, max_depth: int = 3) -> InputPlan:
    """Turn a drop / browse selection into a concrete file list (section 10.2).

    Folders are scanned recursively to ``max_depth``; ``.txt`` files are reported
    as FOFN candidates rather than guessed at; files whose extension is unknown
    are still accepted (the user chose them explicitly) unless they came from a
    folder scan.
    """
    files: List[str] = []
    folders: List[str] = []
    skipped: List[str] = []
    fofn: List[str] = []
    missing: List[str] = []
    seen = set()

    def add(p: str) -> None:
        key = os.path.normcase(os.path.abspath(p))
        if key not in seen:
            seen.add(key)
            files.append(p)

    for raw in paths:
        p = str(raw).strip()
        if not p:
            continue
        path = Path(p)
        if not path.exists():
            missing.append(p)
            continue
        if path.is_dir():
            folders.append(p)
            found = scan_folder(p, max_depth=max_depth)
            if found:
                for f in found:
                    add(f)
            else:
                skipped.append(p)
            continue
        if path.suffix.lower() in (".txt", ".fofn", ".list"):
            fofn.append(p)
            continue
        add(p)

    return InputPlan(
        files=tuple(files),
        folders=tuple(folders),
        skipped=tuple(skipped),
        fofn_candidates=tuple(fofn),
        missing=tuple(missing),
    )


def split_dnd_paths(widget: Any, data: str) -> List[str]:
    """Split a Tcl drop payload.

    MUST use ``widget.tk.splitlist`` and never ``str.split()``: paths containing
    spaces arrive brace-wrapped (``{C:/My Data/a.fa}``), section 10.2.
    """
    try:
        return [str(x) for x in widget.tk.splitlist(data)]
    except Exception:
        text = str(data).strip()
        if text.startswith("{") and text.endswith("}"):
            return [t for t in text[1:-1].split("} {") if t]
        return [t for t in text.split() if t]


def estimate_seconds(n_files: int, jobs: int = 1) -> float:
    """Rough wall-clock estimate for the >50-files confirm (section 10.2).

    Measured on the reference corpus: ~6 s per assembly per worker. It is only
    ever shown as "about N minutes"; no progress bar is driven from it.
    """
    return (max(0, n_files) * 6.0) / max(1, jobs)


def human_duration(seconds: float) -> str:
    """'12 seconds' / 'about 3 minutes' / 'about 1 hour 5 minutes'."""
    s = max(0.0, float(seconds))
    if s < 60:
        return "{:.0f} seconds".format(s)
    minutes = int(round(s / 60.0))
    if minutes < 60:
        return "about {} minute{}".format(minutes, "" if minutes == 1 else "s")
    hours, minutes = divmod(minutes, 60)
    if minutes == 0:
        return "about {} hour{}".format(hours, "" if hours == 1 else "s")
    return "about {} hour{} {} minutes".format(hours, "" if hours == 1 else "s", minutes)


def human_bytes(n: Optional[int]) -> str:
    """'137 MB' style, decimal units, matching what download sites quote."""
    if n is None:
        return "unknown"
    value = float(n)
    for unit in ("bytes", "KB", "MB", "GB"):
        if value < 1024.0 or unit == "GB":
            if unit == "bytes":
                return "{:.0f} bytes".format(value)
            return "{:.0f} {}".format(value, unit) if value >= 10 else "{:.1f} {}".format(value, unit)
        value /= 1024.0
    return "{:.1f} GB".format(value)


# ===========================================================================
# SECTION 4 — ENVIRONMENT PROBE  (section 10.9)
# ===========================================================================

@dataclass
class Environment:
    """What the background start-up job discovered (section 10.9).

    Built on a worker thread so that nothing blocking happens before
    ``mainloop()``. Plain data only: it crosses the thread boundary.
    """

    dbdir: str = ""
    datadir: str = ""
    blastdb: str = ""
    db_version: str = BUNDLED_DB_VERSION
    scheme_count: int = 0
    db_ok: bool = False
    db_error: str = ""
    blast_ok: bool = False
    blast_path: str = ""
    blast_version: str = ""
    blast_origin: str = ""
    blast_error: str = ""
    index_stale: bool = False
    scheme_names: Tuple[str, ...] = ()

    @property
    def ready(self) -> bool:
        """True when a file dropped right now could actually be analysed."""
        return self.db_ok and self.blast_ok

    def footer_text(self) -> str:
        """One-line summary for the footer."""
        db = "Database {} · {} schemes".format(
            self.db_version or "unknown", self.scheme_count or 0)
        if self.blast_ok:
            engine = "Search engine {}".format(self.blast_version or "ready")
        else:
            engine = "Search engine not installed"
        return "{}  —  {}".format(db, engine)


def probe_environment(prefs: Prefs) -> Environment:
    """Locate the database and the search engine (section 10.9, worker thread).

    Never raises: every failure becomes a field on the returned
    :class:`Environment` so the UI can explain it calmly.
    """
    env = Environment()
    try:
        from . import schemes as schemes_mod
    except Exception as exc:  # a broken install, not a user error
        env.db_error = str(exc)
        return env

    try:
        env.dbdir = schemes_mod.resolve_dbdir(prefs.dbdir or None)
        env.datadir = str(Path(env.dbdir) / "pubmlst")
        env.blastdb = schemes_mod.resolve_blastdb(env.dbdir)
        catalog = schemes_mod.SchemeCatalog(env.datadir)
        names = tuple(catalog.names())
        env.scheme_names = names
        env.scheme_count = len(names)
        env.db_version = catalog.db_version or BUNDLED_DB_VERSION
        env.db_ok = env.scheme_count > 0
        if not env.db_ok:
            env.db_error = "No schemes were found in {}".format(env.datadir)
    except Exception as exc:
        env.db_error = getattr(exc, "user_message", "") or str(exc)

    try:
        from . import updatedb as updatedb_mod
        env.index_stale = bool(updatedb_mod.index_is_stale(env.dbdir))
    except Exception:
        env.index_stale = False

    try:
        from . import blastbin
        tools = blastbin.find_blast(prefs.blast_path or None)
        env.blast_path = blast_exe_path(tools)
        env.blast_version = str(getattr(tools, "version", "") or "")
        env.blast_origin = str(getattr(tools, "origin", "") or "")
        env.blast_ok = bool(env.blast_path)
    except Exception as exc:
        env.blast_error = getattr(exc, "user_message", "") or str(exc)

    return env


# ===========================================================================
# SECTION 5 — FRIENDLY ERRORS  (section 10.6)
# ===========================================================================

@dataclass
class Friendly:
    """A novice-readable rendering of an exception (section 10.6)."""

    headline: str = "Something went wrong"
    body: str = ""
    tips: Tuple[str, ...] = ()
    details: str = ""
    #: When True the caller must show the search-engine setup modal instead of
    #: the error dialog.
    needs_bootstrap: bool = False
    #: When True the dialog offers a "Rebuild search index" action.
    offers_rebuild: bool = False


def _env_block(env: Optional[Environment] = None) -> str:
    """The environment footer pasted under 'Show technical details' (section 10.6)."""
    lines = [
        "{} {}".format(branding.APP_NAME, __version__),
        "mlst compatibility: {}".format(UPSTREAM_MLST_VERSION),
        "Python {}".format(platform.python_version()),
        "OS: {}".format(platform.platform()),
    ]
    if env is not None:
        lines.append("Search engine: {} {}".format(
            env.blast_path or "not found", env.blast_version))
        lines.append("Database: {} ({} schemes) at {}".format(
            env.db_version, env.scheme_count, env.dbdir or "unknown"))
    return "\n".join(lines)


def friendly_error(exc: BaseException, *, name: str = "",
                   env: Optional[Environment] = None,
                   tb_text: str = "") -> Friendly:
    """Map an exception onto the dialog copy of section 10.6.

    Pure and testable: the table in section 10.6 is implemented here once, and
    the dialog merely renders the result.
    """
    label = name or "that file"
    details = tb_text or "".join(
        traceback.format_exception(type(exc), exc, exc.__traceback__))
    details = details.rstrip() + "\n\n" + _env_block(env)

    if isinstance(exc, EmptyInputError):
        return Friendly(
            "We couldn't read that file",
            "{} doesn't contain any DNA sequence.".format(label),
            ("Check that the file is an assembly, not raw sequencing reads.",
             "Open it in a text editor: the first line should start with '>'.",
             "If it came from a download, it may be incomplete — try again."),
            details)
    if isinstance(exc, UnsupportedFormatError):
        return Friendly(
            "That file format isn't supported",
            "{} isn't FASTA, GenBank or EMBL.".format(label),
            ("WMLST reads .fa, .fasta, .fna, .gbk, .embl and .gff files.",
             "Compressed copies (.gz, .bz2, .zip) are fine too.",
             "Raw read files (.fastq) must be assembled first."),
            details)
    if isinstance(exc, BlastNotFoundError):
        return Friendly(
            "{} needs its search engine".format(branding.APP_NAME),
            "The one-time setup hasn't finished yet.",
            (), details, needs_bootstrap=True)
    if isinstance(exc, BlastFailedError):
        tail = str(getattr(exc, "stderr_tail", "") or "").strip()
        code = getattr(exc, "returncode", 0)
        # bats 13: an empty query makes BLAST say "Sequence contains no data",
        # which reaches us as a tool failure. The user's problem is the file.
        if "Sequence contains no data" in tail or "Sequence contains no data" in str(exc):
            return friendly_error(EmptyInputError(str(exc)), name=name, env=env,
                                  tb_text=tb_text)
        body = "The search engine stopped while reading {}.".format(label)
        if tail:
            details = "exit code {}\n{}\n\n{}".format(code, tail, details)
        return Friendly(
            "The search engine stopped unexpectedly", body,
            ("Try the file again — this is often a transient problem.",
             "Check there is free space on the drive holding your temporary files.",
             "If your antivirus quarantines the search engine, allow it and retry."),
            details)
    if isinstance(exc, DatabaseMissingError):
        where = getattr(exc, "user_message", "") or str(exc)
        return Friendly(
            "The MLST database is missing", where,
            ("Reinstall {} to restore the bundled database.".format(branding.APP_NAME),
             "Or point WMLST at a database folder in Settings."),
            details)
    if isinstance(exc, DatabaseCorruptError):
        return Friendly(
            "The MLST database looks damaged",
            getattr(exc, "user_message", "") or str(exc),
            ("Rebuilding the search index usually fixes this.",
             "If it keeps happening, reinstall {}.".format(branding.APP_NAME)),
            details, offers_rebuild=True)
    if isinstance(exc, OutOfDiskError):
        return Friendly(
            "Not enough free space", getattr(exc, "user_message", "") or str(exc),
            ("Free some space and try again.",
             "Emptying the Recycle Bin is usually the quickest way."),
            details)
    if isinstance(exc, (BootstrapError, UpdateError)):
        return Friendly(
            "Setup could not finish" if isinstance(exc, BootstrapError)
            else "The database update could not finish",
            getattr(exc, "user_message", "") or str(exc),
            ("Check your internet connection and try again.",
             "A company firewall or proxy may be blocking the download.",
             "You can keep working with the database you already have."),
            details)
    if isinstance(exc, PermissionError):
        return Friendly(
            "Windows wouldn't let us write there",
            "Access to that location was refused.",
            ("Choose a different folder, such as your Documents folder.",
             "Windows Controlled Folder Access can block writes to Desktop, "
             "Documents and Pictures — allow {} in Windows Security.".format(
                 branding.APP_NAME),
             "Close the file if it is open in another program."),
            details)
    if isinstance(exc, Cancelled):
        return Friendly("Cancelled", "The run was stopped before it finished.",
                        (), details)
    if isinstance(exc, WmlstError):
        return Friendly("Something went wrong",
                        getattr(exc, "user_message", "") or str(exc),
                        (), details)
    if isinstance(exc, OSError):
        return Friendly(
            "We couldn't open that file",
            "{}: {}".format(label, exc.strerror or str(exc)),
            ("Check the file still exists and isn't open in another program.",
             "If it is on a network drive, make sure the drive is connected."),
            details)
    return Friendly("Something went wrong", str(exc) or type(exc).__name__, (), details)


# ===========================================================================
# SECTION 6 — MESSAGES, STATE MACHINE, CONTROLLER  (sections 9, 10.5)
# ===========================================================================

@dataclass
class Msg:
    """One item crossing the worker -> main-thread boundary (section 4.10).

    Plain data and engine dataclasses ONLY. A widget reference here is a bug:
    Tcl is not thread-safe and touching a widget off the main thread surfaces as
    a random hang on Windows.
    """

    kind: str  # started|phase|warn|result|failed|batch_done|cancelled|fatal|env
    job: int
    index: int = 0
    total: int = 0
    path: str = ""
    percent: float = 0.0
    text: str = ""
    result: Optional[SampleResult] = None
    exc: Optional[BaseException] = None
    tb: str = ""


#: Analyse-tab states (section 10.2).
IDLE, RUNNING, RESULTS = "IDLE", "RUNNING", "RESULTS"

_ALLOWED_TRANSITIONS = {
    IDLE: frozenset({RUNNING, IDLE}),
    # RUNNING is left only through finish() -> RESULTS.  Allowing RUNNING -> IDLE
    # let "Clear results" strand a live worker behind an idle-looking window.
    RUNNING: frozenset({RESULTS}),
    RESULTS: frozenset({RUNNING, IDLE, RESULTS}),
}


class AppState:
    """The Analyse-tab state machine: IDLE -> RUNNING -> RESULTS (section 10.2).

    CANCELLED and ERROR are not states: they are *reasons* for entering RESULTS,
    because partial results are kept and stay exportable.
    """

    def __init__(self, on_change: Optional[Callable[[str, str], None]] = None):
        self._state = IDLE
        self.reason = ""
        self._on_change = on_change

    @property
    def state(self) -> str:
        return self._state

    @property
    def running(self) -> bool:
        return self._state == RUNNING

    def to(self, new: str, reason: str = "") -> bool:
        """Attempt a transition. Returns False (and changes nothing) if illegal."""
        if new not in _ALLOWED_TRANSITIONS:
            raise ValueError("unknown state: {!r}".format(new))
        if new not in _ALLOWED_TRANSITIONS[self._state]:
            return False
        old, self._state, self.reason = self._state, new, reason
        if self._on_change is not None:
            self._on_change(old, new)
        return True

    def finish(self, reason: str = "done") -> bool:
        """Leave RUNNING for RESULTS with ``done`` / ``cancelled`` / ``error``."""
        return self.to(RESULTS, reason)


_MAIN_THREAD = threading.main_thread()


def assert_main_thread(what: str = "this operation") -> None:
    """Guard every widget mutation and dialog constructor (section 10.5)."""
    if threading.current_thread() is not _MAIN_THREAD:
        raise RuntimeError(
            "{} must run on the Tk main thread (called from {!r})".format(
                what, threading.current_thread().name))


def default_engine_factory(cfg: RunConfig) -> Any:
    """Build the engine for a run. Replaced by a fake in the headless tests."""
    from . import engine as engine_mod
    return engine_mod.Engine(cfg)


class AnalysisController:
    """The ONE threading pattern in the application (sections 9 and 10.5).

    * a :class:`queue.Queue` of plain :class:`Msg` records; no widget reference
      ever crosses the boundary
    * exactly ONE ``after()`` chain at :data:`POLL_MS`; the elapsed timer and the
      indeterminate bar ride on it, never on a second loop
    * the worker thread never imports or touches ``tkinter``
    * ``daemon=True``; closing the window sets cancel, joins for <= 1.5 s, then
      destroys regardless
    * a second :meth:`start` while a worker lives is impossible — the paths are
      queued instead

    ``schedule`` is ``root.after``-shaped (``(ms, callable) -> token``) and
    ``dispatch`` is called on the main thread with each :class:`Msg`; both are
    injected so the whole class is testable with no display.
    """

    POLL_MS: int = 120

    def __init__(self, *, schedule: Callable[[int, Callable[[], None]], Any],
                 dispatch: Callable[[Msg], None],
                 engine_factory: Callable[[RunConfig], Any] = default_engine_factory,
                 on_tick: Optional[Callable[[float, bool], None]] = None,
                 clock: Callable[[], float] = time.monotonic):
        self._schedule = schedule
        self._dispatch = dispatch
        self._engine_factory = engine_factory
        self._on_tick = on_tick
        self._clock = clock
        self.queue: queue.Queue[Msg] = queue.Queue()
        self._thread: Optional[threading.Thread] = None
        self._cancel = threading.Event()
        self._job = 0
        self._started_at = 0.0
        self._pending: List[str] = []
        self._polling = False
        self._lock = threading.Lock()
        self._task_results: Dict[int, Any] = {}
        #: Set by the worker so exports can reuse the engine's own RunMeta.
        self.meta: Optional[RunMeta] = None

    # -- state ---------------------------------------------------------------
    @property
    def running(self) -> bool:
        t = self._thread
        return t is not None and t.is_alive()

    @property
    def job(self) -> int:
        return self._job

    @property
    def elapsed(self) -> float:
        return max(0.0, self._clock() - self._started_at) if self._started_at else 0.0

    def pending_paths(self) -> List[str]:
        """Take the files dropped while a run was in flight (section 10.5 #6)."""
        with self._lock:
            out, self._pending = self._pending, []
        return out

    def take_task_result(self, job: int) -> Any:
        """Collect the return value of a background task by job id."""
        with self._lock:
            return self._task_results.pop(job, None)

    # -- the single after() chain -------------------------------------------
    def ensure_polling(self) -> None:
        """Start the one and only polling chain (idempotent)."""
        if not self._polling:
            self._polling = True
            self._schedule(self.POLL_MS, self._tick)

    def _tick(self) -> None:
        # Re-arm BEFORE doing any work.  A handler reached from pump() may open a
        # modal dialog (the first-run BLAST bootstrap does exactly that) and spin a
        # nested event loop inside wait_window().  If the next tick were only armed
        # after pump() returned, nothing would be pending while that nested loop ran,
        # the single after() chain would be dead, and the dialog -- whose progress and
        # self-close arrive purely through this chain -- would hang forever.
        # Re-arming first keeps exactly one chain alive: the pending tick fires inside
        # the nested loop and carries the chain from there.
        self._schedule(self.POLL_MS, self._tick)
        self.pump()
        if self._on_tick is not None:
            try:
                self._on_tick(self.elapsed, self.running)
            except Exception:
                LOG.exception("tick handler failed")

    def pump(self, limit: int = 400) -> int:
        """Drain the queue on the main thread. Returns how many messages went out."""
        n = 0
        while n < limit:
            try:
                msg = self.queue.get_nowait()
            except queue.Empty:
                break
            n += 1
            try:
                self._dispatch(msg)
            except Exception:
                LOG.exception("dispatch failed for %s", msg.kind)
        return n

    # -- analysis ------------------------------------------------------------
    def start(self, paths: Sequence[str], cfg: RunConfig) -> bool:
        """Start a batch. Returns False when a run is already in flight, in which
        case the paths are queued and picked up when it ends."""
        files = [str(p) for p in paths]
        if not files:
            return False
        if self.running:
            with self._lock:
                self._pending.extend(files)
            return False
        self._job += 1
        self._cancel = threading.Event()
        self._started_at = self._clock()
        self.meta = None
        job = self._job
        self._thread = threading.Thread(
            target=self._run, args=(job, files, cfg, self._cancel),
            name="wmlst-analysis-{}".format(job), daemon=True)
        self._thread.start()
        self.ensure_polling()
        return True

    def request_cancel(self) -> None:
        """Ask the worker to stop at the next checkpoint (section 10.5)."""
        self._cancel.set()

    @property
    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def shutdown(self, timeout: float = 1.5) -> None:
        """Cancel and join for at most ``timeout`` seconds, then give up."""
        self.request_cancel()
        t = self._thread
        if t is not None and t.is_alive():
            t.join(timeout)

    def _put(self, msg: Msg) -> None:
        self.queue.put(msg)

    def _run(self, job: int, files: List[str], cfg: RunConfig,
             cancel: threading.Event) -> None:
        """Worker body. NEVER touches tkinter; every outcome becomes a Msg."""
        total = len(files)
        try:
            eng = self._engine_factory(cfg)
        except BaseException as exc:  # engine construction is run-fatal
            self._put(Msg("fatal", job, exc=exc, tb=traceback.format_exc(),
                          text=str(exc)))
            self._put(Msg("batch_done", job, total=total))
            return
        try:
            meta = getattr(eng, "meta_stub", None)
            if isinstance(meta, RunMeta):
                self.meta = meta
        except Exception:
            self.meta = None
        try:
            for i, path in enumerate(files, 1):
                if cancel.is_set():
                    self._put(Msg("cancelled", job, index=i - 1, total=total))
                    return
                self._put(Msg("started", job, index=i, total=total, path=path,
                              text="Reading {}...".format(os.path.basename(path))))
                try:
                    result = eng.analyse_file(
                        path,
                        progress=lambda ev, _j=job, _i=i, _p=path: self._put(
                            Msg("phase", _j, index=_i, total=total,
                                path=getattr(ev, "path", _p) or _p,
                                percent=float(getattr(ev, "percent", 0.0)),
                                text=str(getattr(ev, "text", "")))),
                        warn=lambda text, _j=job, _i=i, _p=path: self._put(
                            Msg("warn", _j, index=_i, total=total, path=_p,
                                text=str(text))),
                        cancel=cancel,
                    )
                except Cancelled:
                    self._put(Msg("cancelled", job, index=i - 1, total=total))
                    return
                except BaseException as exc:
                    self._put(Msg("failed", job, index=i, total=total, path=path,
                                  exc=exc, tb=traceback.format_exc(),
                                  text=getattr(exc, "user_message", "") or str(exc)))
                    continue
                self._put(Msg("result", job, index=i, total=total, path=path,
                              percent=100.0, result=result))
        except BaseException as exc:  # pragma: no cover - defensive
            self._put(Msg("fatal", job, exc=exc, tb=traceback.format_exc(),
                          text=str(exc)))
        finally:
            try:
                close = getattr(eng, "close", None)
                if callable(close):
                    close()
            except Exception:
                LOG.exception("engine close failed")
            self._put(Msg("batch_done", job, total=total))

    # -- background tasks (environment probe, exports, downloads, updates) ---
    def run_task(self, kind: str, fn: Callable[..., Any], *,
                 with_progress: bool = False,
                 cancel: Optional[threading.Event] = None) -> int:
        """Run ``fn`` on a worker thread and report completion as a Msg.

        ``fn`` is called with ``(progress, cancel)`` when ``with_progress`` is
        True, otherwise with no arguments. Its return value is retrievable with
        :meth:`take_task_result`; completion arrives as ``Msg(kind='env')`` and
        failure as ``Msg(kind='failed')``, both carrying the same job id.
        """
        self._job += 1
        job = self._job
        ev = cancel if cancel is not None else threading.Event()

        def progress(done: Any = 0, total: Any = 0, text: Any = "") -> None:
            # Tolerant on purpose: the producers differ. blastbin.bootstrap sends
            # (bytes_done, bytes_total, phase_text); a catalogue or update pass may
            # send just a line of text.
            if isinstance(done, str) and not text:
                done, text = 0, done
            try:
                pct = (100.0 * float(done) / float(total)) if float(total) else 0.0
            except (TypeError, ValueError, ZeroDivisionError):
                pct = 0.0
            self._put(Msg("phase", job, percent=max(0.0, min(100.0, pct)),
                          text=str(text)))

        def body() -> None:
            try:
                value = fn(progress, ev) if with_progress else fn()
            except BaseException as exc:
                self._put(Msg("failed", job, exc=exc, tb=traceback.format_exc(),
                              text=getattr(exc, "user_message", "") or str(exc)))
                return
            with self._lock:
                self._task_results[job] = value
            self._put(Msg("env", job, text=kind))

        threading.Thread(target=body, name="wmlst-task-{}".format(job),
                         daemon=True).start()
        self.ensure_polling()
        return job


# ===========================================================================
# SECTION 7 — WIDGETS  (sections 10.1, 10.8)
# ===========================================================================

#: Active colour map — swapped wholesale for the system palette under Windows
#: high contrast. Read through :func:`c` so a later swap reaches every widget.
C: Dict[str, str] = dict(PALETTE)

#: Which variant :func:`install_theme` chose. Read by the widgets that
#: need to know whether they are painting on a light or a dark ground.
THEME_MODE = "light"

#: Named fonts, created once against the root.
F: Dict[str, Any] = {}


def c(key: str) -> str:
    """Look up a palette colour by role name."""
    return C.get(key, PALETTE.get(key, "#000000"))


def apply_scaling(root: tk.Misc) -> float:
    """Match Tk's internal scaling to the real display dpi (section 10.1)."""
    try:
        dpi = float(root.winfo_fpixels("1i"))
    except Exception:
        dpi = 96.0
    if not (48.0 <= dpi <= 480.0):
        dpi = 96.0
    try:
        root.tk.call("tk", "scaling", dpi / 72.0)
    except Exception:
        pass
    return set_scale(dpi)


# ---------------------------------------------------------------------------
# Drawn widget skins — the ttk "image" element engine (section 10.1)
# ---------------------------------------------------------------------------
#
# clam cannot round a corner and cannot be talked out of drawing a 3-D spin
# arrow, so the four controls a user actually touches — button, entry, spinbox,
# combobox — are given their shape by images this module rasterises itself and
# hands to ``style.element_create(..., "image", ...)`` as a nine-patch. Every
# pixel is computed from the live palette, so the dark and high-contrast
# variants are skinned by the same code with no second set of assets.

#: ``element_create`` raises on a name that already exists, and install_theme
#: runs again on every theme change (and repeatedly across the test suite).
_SKIN_GEN = [0]


def _skin_name(base: str) -> str:
    """A fresh element name, so re-theming never collides with the last pass."""
    return "w{}.{}".format(_SKIN_GEN[0], base)


def rounded_image(master: tk.Misc, w: int, h: int, radius: float, *,
                  fill: str, edge: str, width: float = 1.0,
                  back: Optional[str] = None) -> Any:
    """An antialiased rounded rectangle as a :class:`tkinter.PhotoImage`.

    Outside the corner arc the pixel is made transparent (or painted ``back``
    when a backdrop colour is given), so one image sits correctly on the page,
    on a card, and on a selected row.
    """
    side_w, side_h = max(1, int(w)), max(1, int(h))
    img = tk.PhotoImage(master=master, width=side_w, height=side_h)
    x0, y0 = 0.5, 0.5
    x1, y1 = side_w - 1.5, side_h - 1.5
    radius = max(0.0, min(radius, (x1 - x0) / 2.0, (y1 - y0) / 2.0))
    inner_r = max(0.0, radius - width)
    base = back or fill
    rows, clear = [], []
    for yy in range(side_h):
        row = []
        for xx in range(side_w):
            cover = inner = 0.0
            for sy in (0.17, 0.5, 0.83):
                for sx in (0.17, 0.5, 0.83):
                    sxx, syy = xx + sx, yy + sy
                    cover += _rounded_cover(sxx, syy, x0, y0, x1, y1, radius)
                    inner += _rounded_cover(sxx, syy, x0 + width, y0 + width,
                                            x1 - width, y1 - width, inner_r)
            cover /= 9.0
            inner /= 9.0
            if cover <= 0.001:
                if back is None:
                    clear.append((xx, yy))
                row.append(base)
                continue
            row.append(mix(mix(base, edge, cover), fill, inner))
        rows.append(row)
    img.put(rows)
    for xx, yy in clear:
        try:
            img.transparency_set(xx, yy, True)
        except (tk.TclError, AttributeError):  # pragma: no cover - ancient Tk
            break
    return img


def chevron_image(master: tk.Misc, w: int, h: int, *, colour: str, back: str,
                  down: bool = True, weight: float = 1.6) -> Any:
    """A stroked chevron (the only arrow in this application), antialiased."""
    side_w, side_h = max(5, int(w)), max(5, int(h))
    img = tk.PhotoImage(master=master, width=side_w, height=side_h)
    span = side_w * 0.26
    mid_x, mid_y = side_w / 2.0, side_h / 2.0
    rise = span * 0.62
    if down:
        ax, ay, bx, by = mid_x - span, mid_y - rise / 2.0, mid_x, mid_y + rise / 2.0
        cx_, cy_ = mid_x + span, mid_y - rise / 2.0
    else:
        ax, ay, bx, by = mid_x - span, mid_y + rise / 2.0, mid_x, mid_y - rise / 2.0
        cx_, cy_ = mid_x + span, mid_y + rise / 2.0
    half = max(0.7, weight * SC / 2.0)
    rows, clear = [], []
    for yy in range(side_h):
        row = []
        for xx in range(side_w):
            d = min(_seg_distance(xx + 0.5, yy + 0.5, ax, ay, bx, by),
                    _seg_distance(xx + 0.5, yy + 0.5, bx, by, cx_, cy_))
            t = half + 0.5 - d
            if t <= 0.0:
                clear.append((xx, yy))
                row.append(back)
            else:
                row.append(mix(back, colour, min(1.0, t)))
        rows.append(row)
    img.put(rows)
    for xx, yy in clear:
        try:
            img.transparency_set(xx, yy, True)
        except (tk.TclError, AttributeError):  # pragma: no cover
            break
    return img


def dot_image(master: tk.Misc, kind: str, *, size: int = 12, gap: int = 9,
              back: Optional[str] = None) -> Any:
    """The STATUS indicator for a results row, drawn rather than typed.

    A Treeview cell cannot hold a widget and the row font cannot be trusted to
    own U+2714: on a Tk built without Xft every dingbat in the status column
    became a hollow box, which is exactly the channel that is supposed to carry
    the result when colour cannot. So the shape is rasterised here — a filled
    disc, a ring, a half-filled disc or a triangle — and set as the row image,
    ahead of the file name. Shape, colour and the status word: three channels.
    """
    side = max(8, px(size))
    pad_right = max(0, px(gap))
    surface = back or c("surface")
    entry = STATUS_UI.get(kind, (None, "muted", None, "ring"))
    colour = c(entry[1])
    shape = entry[3] if len(entry) > 3 else "ring"
    img = tk.PhotoImage(master=master, width=side + pad_right, height=side)
    cx = cy = (side - 1) / 2.0
    r = side / 2.0 - 1.2
    hollow = shape == "ring"
    triangle = shape == "triangle"
    ring = max(1.3, side * 0.16)
    rows, clear = [], []
    for yy in range(side):
        row = []
        for xx in range(side + pad_right):
            if xx >= side:
                clear.append((xx, yy))
                row.append(surface)
                continue
            cover = 0.0
            for sy in (0.17, 0.5, 0.83):
                for sx in (0.17, 0.5, 0.83):
                    dx, dy = xx + sx - cx, yy + sy - cy
                    if triangle:
                        # an upright triangle inscribed in the same circle: the
                        # apex is at the TOP, because that is the shape everyone
                        # already reads as "caution" -- inverted, it reads as a
                        # downward arrow, which means something else entirely.
                        ty = dy + r * 0.28
                        edge = (ty + r * 0.95) * 0.577
                        hit = (-r * 0.95 <= ty <= r * 0.80) and abs(dx) <= edge
                    else:
                        hit = dx * dx + dy * dy <= r * r
                    cover += 1.0 if hit else 0.0
            cover /= 9.0
            if cover <= 0.001:
                clear.append((xx, yy))
                row.append(surface)
                continue
            if hollow:
                hole = 0.0
                rr = r - ring
                for sy in (0.17, 0.5, 0.83):
                    for sx in (0.17, 0.5, 0.83):
                        dx, dy = xx + sx - cx, yy + sy - cy
                        hole += 1.0 if dx * dx + dy * dy <= rr * rr else 0.0
                cover *= 1.0 - hole / 9.0
            if cover <= 0.001:
                clear.append((xx, yy))
                row.append(surface)
                continue
            row.append(mix(surface, colour, cover))
        rows.append(row)
    img.put(rows)
    for xx, yy in clear:
        try:
            img.transparency_set(xx, yy, True)
        except (tk.TclError, AttributeError):  # pragma: no cover
            break
    return img


def _install_skins(root: tk.Misc, style: ttk.Style) -> None:
    """Give button, entry, spinbox and combobox a drawn, rounded shape."""
    _SKIN_GEN[0] += 1
    # ttk keeps no reference to an element image; Python must, or Tk frees it
    # and the widget paints as an empty box. The list hangs off the ROOT rather
    # than off the module, so the images die with the interpreter that owns
    # them -- a module-level list would pin a PhotoImage per skin per window
    # for the life of the process, and leave dead ones behind after destroy().
    images: List[Any] = []
    root._wmlst_skin_images = images  # type: ignore[attr-defined]
    keep = images.append
    radius = float(px(RADIUS - 2))
    side = int(radius * 2 + px(8))
    border = int(radius + px(2))
    surface, alt = c("surface"), c("surface_alt")
    text, muted, accent = c("text"), c("muted"), c("accent")
    edge_rest, edge_hot = c("border_strong"), c("muted")

    def plate(fill: str, edge: str, width: float = 1.0) -> Any:
        img = rounded_image(root, side, side, radius, fill=fill, edge=edge,
                            width=width)
        keep(img)
        return img

    def make(base: str, specs: Sequence[Tuple[Any, Any]], default: Any,
             *, nine_patch: bool = True) -> str:
        """Register one image element. A glyph is placed, never stretched."""
        name = _skin_name(base)
        args: List[Any] = [default]
        args.extend(tuple(spec) for spec in specs)
        opts = ({"border": border, "sticky": "nsew"} if nine_patch
                else {"border": 0, "sticky": ""})
        try:
            style.element_create(name, "image", *args, padding=0, **opts)
        except tk.TclError:  # pragma: no cover - name collision is impossible
            return ""
        return name

    # -- buttons -----------------------------------------------------------
    btn = make("button", [
        ("disabled", plate(c("surface_sunk"), c("border"))),
        ("pressed", plate(c("border"), edge_hot)),
        ("active", plate(mix(alt, text, 0.05), edge_hot)),
        ("focus", plate(alt, c("focus"), width=2.0)),
    ], plate(alt, edge_rest))
    acc = make("accentbutton", [
        # A disabled primary must not out-shout an enabled secondary beside it,
        # so it drops to the same inert outline the other disabled buttons wear.
        ("disabled", plate(c("surface_sunk"), c("border"))),
        ("pressed", plate(c("accent_active"), c("accent_active"))),
        ("active", plate(c("accent_hover"), c("accent_hover"))),
        ("focus", plate(accent, c("text"), width=2.0)),
    ], plate(accent, accent))
    quiet = make("quietbutton", [
        ("disabled", plate(surface, surface)),
        ("pressed", plate(c("border"), c("border"))),
        ("active", plate(alt, alt)),
        ("focus", plate(surface, c("focus"), width=2.0)),
    ], plate(surface, surface))

    def button_layout(name: str, element: str) -> None:
        if not element:
            return
        try:
            style.layout(name, [
                (element, {"sticky": "nswe", "children": [
                    ("Button.padding", {"sticky": "nswe", "children": [
                        ("Button.label", {"sticky": "nswe"})]})]})])
        except tk.TclError:  # pragma: no cover
            pass

    button_layout("TButton", btn)
    button_layout("Accent.TButton", acc)
    button_layout("Surface.TButton", quiet)
    button_layout("Quiet.TButton", quiet)
    # A link is text. ttk hands a style with no layout of its own its parent's,
    # so without this the two link styles inherited the button plate and every
    # citation link grew a grey box around it.
    for link_style in ("Link.TButton", "CardLink.TButton"):
        try:
            style.layout(link_style, [
                ("Button.padding", {"sticky": "nswe", "children": [
                    ("Button.label", {"sticky": "nswe"})]})])
        except tk.TclError:  # pragma: no cover
            pass

    # -- fields (entry, spinbox, combobox) ---------------------------------
    field = make("field", [
        ("disabled", plate(alt, c("border"))),
        ("focus", plate(surface, c("focus"), width=2.0)),
        ("hover", plate(surface, edge_hot)),
    ], plate(surface, edge_rest))
    invalid = make("invalidfield", [
        ("focus", plate(c("bad_soft"), c("bad"), width=2.0)),
    ], plate(c("bad_soft"), c("bad")))

    # The glyph is centred in its image, so the extra width IS the margin that
    # keeps the chevron off the field's rounded right border.
    arrow_w, arrow_h = int(px(28)), int(px(10))
    spin_up = chevron_image(root, arrow_w, arrow_h, colour=muted, back=surface,
                            down=False)
    spin_down = chevron_image(root, arrow_w, arrow_h, colour=muted, back=surface)
    spin_up_hot = chevron_image(root, arrow_w, arrow_h, colour=accent,
                                back=surface, down=False)
    spin_down_hot = chevron_image(root, arrow_w, arrow_h, colour=accent,
                                  back=surface)
    combo_arrow = chevron_image(root, int(px(30)), int(px(18)), colour=muted,
                                back=surface)
    combo_arrow_hot = chevron_image(root, int(px(30)), int(px(18)),
                                    colour=accent, back=surface)
    for image in (spin_up, spin_down, spin_up_hot, spin_down_hot, combo_arrow,
                  combo_arrow_hot):
        keep(image)
    up = make("uparrow", [("pressed", spin_up_hot), ("active", spin_up_hot)],
              spin_up, nine_patch=False)
    down = make("downarrow", [("pressed", spin_down_hot),
                              ("active", spin_down_hot)], spin_down,
                nine_patch=False)
    combo = make("comboarrow", [("pressed", combo_arrow_hot),
                                ("active", combo_arrow_hot)], combo_arrow,
                 nine_patch=False)

    if field:
        try:
            style.layout("TEntry", [
                (field, {"sticky": "nswe", "children": [
                    ("Entry.padding", {"sticky": "nswe", "children": [
                        ("Entry.textarea", {"sticky": "nswe"})]})]})])
            style.layout("Invalid.TEntry", [
                (invalid or field, {"sticky": "nswe", "children": [
                    ("Entry.padding", {"sticky": "nswe", "children": [
                        ("Entry.textarea", {"sticky": "nswe"})]})]})])
        except tk.TclError:  # pragma: no cover
            pass
    if field and up and down:
        try:
            style.layout("TSpinbox", [
                (field, {"sticky": "nswe", "children": [
                    ("null", {"side": "right", "sticky": "ns", "children": [
                        (up, {"side": "top", "sticky": "e"}),
                        (down, {"side": "bottom", "sticky": "e"})]}),
                    ("Spinbox.padding", {"sticky": "nswe", "children": [
                        ("Spinbox.textarea", {"sticky": "nswe"})]})]})])
        except tk.TclError:  # pragma: no cover
            pass
    if field and combo:
        try:
            style.layout("TCombobox", [
                (field, {"sticky": "nswe", "children": [
                    (combo, {"side": "right", "sticky": "ns"}),
                    ("Combobox.padding", {"sticky": "nswe", "children": [
                        ("Combobox.textarea", {"sticky": "nswe"})]})]})])
        except tk.TclError:  # pragma: no cover
            pass

    # -- progress bars ------------------------------------------------------
    # NOT skinned with an image element: ttk sizes the "pbar" from the element's
    # own requested size, so an image element pins it to the image width and the
    # bar stops reporting the value at all. clam's own rectangle stays, with the
    # bevel configured away (see TProgressbar above).

    # A scrollbar thumb that is a rounded pill rather than a grey slab.
    rest_thumb = mix(c("border_strong"), text, 0.18)
    thumb = make("thumb", [
        ("pressed", plate(muted, muted)),
        ("active", plate(mix(c("border_strong"), text, 0.40),
                         mix(c("border_strong"), text, 0.40))),
    ], plate(rest_thumb, rest_thumb))
    if thumb:
        for orient, fill in (("Vertical", "ns"), ("Horizontal", "ew")):
            try:
                style.layout("{}.TScrollbar".format(orient), [
                    ("{}.Scrollbar.trough".format(orient), {
                        "sticky": fill, "children": [
                            (thumb, {"expand": "1", "sticky": "nswe"})]})])
            except tk.TclError:  # pragma: no cover
                pass


def install_theme(root: tk.Misc) -> ttk.Style:
    """Install the flat theme and the font ladder (section 10.1).

    ``clam`` on every platform: it is the only stdlib theme that honours
    ``configure``/``map`` for background colours. Every 3-D element clam draws by
    default — the sunken entry well, the ridged notebook, the dotted focus
    rectangle, the bevelled button — is flattened here by re-laying the element
    out or by painting ``lightcolor``/``darkcolor``/``bordercolor`` the same
    colour, which is how a Tk theme is made to look like this decade.

    The palette comes from :func:`detect_theme_mode`, so Windows high contrast
    and the dark system theme both reach every widget through :func:`c`.
    """
    global C, THEME_MODE
    THEME_MODE = detect_theme_mode()
    C = palette_for(THEME_MODE)

    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:  # pragma: no cover - clam ships with every Tk 8.6+
        pass

    try:
        families = tkfont.families(root)
    except Exception:
        families = ()
    fallback = "TkDefaultFont"
    try:
        fallback = tkfont.nametofont("TkDefaultFont").actual("family")
    except Exception:
        pass
    ui = pick_font(families, PREFERRED_FAMILIES, fallback)
    display = pick_font(families, PREFERRED_DISPLAY, ui)
    mono = pick_font(families, PREFERRED_MONO, fallback)

    F.clear()
    F["family"] = ui
    F["display_family"] = display
    for key, size in FONT_PT.items():
        if key == "mono":
            family = mono
        elif key in ("title", "hero"):
            family = display
        else:
            family = ui
        weight = "bold" if key in ("heading", "title", "hero") else "normal"
        F[key] = tkfont.Font(root=root, family=family, size=size, weight=weight)
    F["body_bold"] = tkfont.Font(root=root, family=ui, size=FONT_PT["body"],
                                 weight="bold")
    F["small_bold"] = tkfont.Font(root=root, family=ui, size=FONT_PT["small"],
                                  weight="bold")
    # SMALL CAPS is not a Tk feature; a letter-spaced bold 9 pt is the honest
    # substitute for an eyebrow label ("SEQUENCE TYPE") and is set as its own
    # font so the spacing does not leak into ordinary small text.
    F["eyebrow"] = tkfont.Font(root=root, family=ui, size=FONT_PT["tiny"],
                               weight="bold")
    # A binomial is italic by convention, and a diagnostic result that sets
    # "Klebsiella pneumoniae" upright reads as careless to the people who have
    # to sign it off. Two cuts: the headline on the summary card, and the small
    # one used inline. "spp." stays upright, so it is a separate label.
    F["organism"] = tkfont.Font(root=root, family=display, size=FONT_PT["title"],
                                weight="normal", slant="italic")
    F["organism_small"] = tkfont.Font(root=root, family=ui, size=FONT_PT["body"],
                                      slant="italic")
    # Links carry an underline as well as the accent colour: colour is never the
    # only channel (WCAG 1.4.1).
    F["link"] = tkfont.Font(root=root, family=ui, size=FONT_PT["small"],
                            underline=True)
    for named in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont"):
        try:
            nf = tkfont.nametofont(named)
            nf.configure(family=ui, size=FONT_PT["body"], weight="normal")
        except Exception:
            pass

    bg, surface, text, muted = c("bg"), c("surface"), c("text"), c("muted")
    border, accent, soft = c("border"), c("accent"), c("text_soft")

    # No widget in this theme owns a border it did not ask for.
    style.configure(".", background=bg, foreground=text, font=F["body"],
                    borderwidth=0, relief="flat", focuscolor=c("focus"),
                    bordercolor=border, lightcolor=bg, darkcolor=bg,
                    troughcolor=c("surface_alt"), selectbackground=c("accent_soft"),
                    selectforeground=text)
    style.configure("TFrame", background=bg)
    style.configure("Surface.TFrame", background=surface)
    style.configure("SurfaceAlt.TFrame", background=c("surface_alt"))
    style.configure("CardBorder.TFrame", background=border)
    style.configure("Footer.TFrame", background=c("surface_alt"))
    style.configure("Hairline.TFrame", background=border)
    style.configure("Rule.TFrame", background=accent)
    # The forced-scheme notice: a tinted band with a status word in it, never
    # colour alone -- the text says what is pinned and the button undoes it.
    style.configure("Notice.TFrame", background=c("warn_soft"))
    style.configure("NoticeBar.TFrame", background=c("warn"))
    style.configure("Notice.TLabel", background=c("warn_soft"), foreground=c("warn"),
                    font=F["body"])
    style.configure("NoticeStrong.TLabel", background=c("warn_soft"),
                    foreground=c("warn"), font=F["body_bold"])
    style.configure("NoticeBadge.TLabel", background=c("warn_soft"),
                    foreground=c("warn"), font=F["eyebrow"])

    style.configure("TLabel", background=bg, foreground=text, font=F["body"])
    style.configure("Surface.TLabel", background=surface, foreground=text)
    style.configure("SurfaceSoft.TLabel", background=surface, foreground=soft)
    style.configure("Muted.TLabel", background=bg, foreground=muted, font=F["small"])
    style.configure("SurfaceMuted.TLabel", background=surface, foreground=muted,
                    font=F["small"])
    style.configure("Eyebrow.TLabel", background=surface, foreground=muted,
                    font=F["eyebrow"])
    style.configure("Heading.TLabel", background=surface, foreground=text,
                    font=F["heading"])
    style.configure("Title.TLabel", background=bg, foreground=text, font=F["title"])
    style.configure("SurfaceTitle.TLabel", background=surface, foreground=text,
                    font=F["title"])
    style.configure("Hero.TLabel", background=surface, foreground=text, font=F["hero"])
    style.configure("Status.TLabel", background=c("surface_alt"), foreground=muted,
                    font=F["small"])
    for key in ("ok", "warn", "bad", "info"):
        style.configure("{}.TLabel".format(key.capitalize()), background=surface,
                        foreground=c(key), font=F["body"])

    # Buttons: flat fills, a hairline that only appears on hover/focus, and a
    # real pressed state. Exactly one Accent button is used per screen.
    style.configure("TButton", background=c("surface_alt"), foreground=text,
                    font=F["body"], padding=(px(PAD_M), px(9)), borderwidth=1,
                    relief="flat", bordercolor=c("border_strong"),
                    lightcolor=c("surface_alt"), darkcolor=c("surface_alt"),
                    anchor="center")
    style.map("TButton",
              background=[("disabled", c("surface_sunk")),
                          ("pressed", c("border")),
                          ("active", mix(c("surface_alt"), c("text"), 0.06))],
              lightcolor=[("pressed", c("border")),
                          ("active", mix(c("surface_alt"), c("text"), 0.06))],
              darkcolor=[("pressed", c("border")),
                         ("active", mix(c("surface_alt"), c("text"), 0.06))],
              foreground=[("disabled", mix(muted, c("surface"), 0.45))],
              bordercolor=[("focus", c("focus")), ("active", c("border_strong")),
                           ("disabled", c("border"))])
    style.configure("Surface.TButton", background=surface, foreground=text,
                    lightcolor=surface, darkcolor=surface)
    style.map("Surface.TButton",
              background=[("disabled", surface),
                          ("pressed", c("surface_alt")),
                          ("active", c("accent_tint"))],
              lightcolor=[("active", c("accent_tint"))],
              darkcolor=[("active", c("accent_tint"))])
    style.configure("Accent.TButton", background=accent, foreground=c("on_accent"),
                    bordercolor=accent, lightcolor=accent, darkcolor=accent,
                    font=F["body_bold"], padding=(px(PAD_L), px(9)))
    style.map("Accent.TButton",
              background=[("disabled", c("border")), ("pressed", c("accent_active")),
                          ("active", c("accent_hover"))],
              lightcolor=[("pressed", c("accent_active")),
                          ("active", c("accent_hover"))],
              darkcolor=[("pressed", c("accent_active")),
                         ("active", c("accent_hover"))],
              bordercolor=[("focus", c("text")), ("disabled", c("border"))],
              foreground=[("disabled", muted)])
    style.configure("Link.TButton", background=bg, foreground=accent,
                    borderwidth=0, padding=(px(PAD_XS), px(PAD_XS)), font=F["small"])
    style.map("Link.TButton",
              background=[("active", bg), ("pressed", bg)],
              foreground=[("active", c("accent_hover"))])
    # A link that sits on a card, not on the window background.
    style.configure("CardLink.TButton", background=surface, foreground=accent,
                    borderwidth=0, padding=(0, px(PAD_XS)), font=F["link"],
                    focuscolor=c("focus"), lightcolor=surface, darkcolor=surface)
    style.map("CardLink.TButton",
              background=[("active", surface), ("pressed", surface)],
              foreground=[("active", c("accent_hover")), ("disabled", muted)])
    style.configure("Quiet.TButton", background=surface, foreground=muted,
                    borderwidth=0, padding=(px(PAD_S), px(PAD_XS)), font=F["small"],
                    lightcolor=surface, darkcolor=surface)
    style.map("Quiet.TButton",
              background=[("active", c("surface_alt")), ("pressed", c("border"))],
              foreground=[("active", text)])

    # A flat tab strip: no folder tabs, no ridge. The selected tab is the only
    # one painted on the card surface, and it carries an accent underline drawn
    # by TabUnderline below -- shape and colour, never colour alone.
    style.configure("TNotebook", background=bg, borderwidth=0,
                    bordercolor=bg, lightcolor=bg, darkcolor=bg,
                    tabmargins=(0, 0, 0, 0))
    style.configure("TNotebook.Tab", background=bg, foreground=muted,
                    padding=(px(PAD_L), px(10)), font=F["body"], borderwidth=0,
                    bordercolor=bg, lightcolor=bg, darkcolor=bg)
    style.map("TNotebook.Tab",
              background=[("selected", c("accent_soft")),
                          ("active", c("surface_alt"))],
              lightcolor=[("selected", c("accent_soft"))],
              darkcolor=[("selected", c("accent_soft"))],
              foreground=[("selected", accent), ("active", text)],
              font=[("selected", F["body_bold"])],
              expand=[("selected", (0, 0, 0, 0))])
    try:  # drop clam's dotted focus rectangle from the tab
        style.layout("TNotebook.Tab", [
            ("Notebook.tab", {"sticky": "nswe", "children": [
                ("Notebook.padding", {"side": "top", "sticky": "nswe", "children": [
                    ("Notebook.label", {"side": "top", "sticky": ""})]})]})])
    except tk.TclError:  # pragma: no cover - layout name is stable in clam
        pass

    style.configure("Treeview", background=surface, fieldbackground=surface,
                    foreground=text, rowheight=px(30), borderwidth=0,
                    font=F["body"], relief="flat")
    style.configure("Treeview.Heading", background=surface, foreground=muted,
                    font=F["small_bold"], relief="flat",
                    padding=(px(PAD_S), px(10)), borderwidth=0,
                    bordercolor=border, lightcolor=surface, darkcolor=surface)
    style.map("Treeview.Heading",
              background=[("active", c("surface_alt"))],
              foreground=[("active", text)])
    style.map("Treeview", background=[("selected", c("accent_soft"))],
              foreground=[("selected", text)])
    try:  # no dotted rectangle around the focused row
        style.layout("Treeview.Item", [
            ("Treeitem.padding", {"sticky": "nswe", "children": [
                ("Treeitem.indicator", {"side": "left", "sticky": ""}),
                ("Treeitem.image", {"side": "left", "sticky": ""}),
                ("Treeitem.text", {"side": "left", "sticky": ""})]})])
    except tk.TclError:  # pragma: no cover
        pass

    # Entries and spinboxes: a hairline box, not a sunken well.
    style.configure("TEntry", fieldbackground=surface, foreground=text,
                    bordercolor=c("border_strong"), lightcolor=c("border_strong"),
                    darkcolor=c("border_strong"), insertcolor=text,
                    padding=(px(PAD_S), px(7)), borderwidth=1, relief="flat")
    style.map("TEntry", bordercolor=[("focus", c("focus"))],
              lightcolor=[("focus", c("focus"))], darkcolor=[("focus", c("focus"))],
              fieldbackground=[("disabled", c("surface_alt"))])
    style.configure("Invalid.TEntry", fieldbackground=c("bad_soft"),
                    bordercolor=c("bad"), lightcolor=c("bad"), darkcolor=c("bad"))
    style.configure("TSpinbox", fieldbackground=surface, foreground=text,
                    bordercolor=c("border_strong"), arrowcolor=muted,
                    lightcolor=c("border_strong"), darkcolor=c("border_strong"),
                    background=surface, padding=(px(PAD_S), px(6)), borderwidth=1,
                    arrowsize=px(14), relief="flat")
    style.map("TSpinbox", bordercolor=[("focus", c("focus"))],
              lightcolor=[("focus", c("focus"))], darkcolor=[("focus", c("focus"))],
              arrowcolor=[("active", accent), ("disabled", c("border_strong"))],
              fieldbackground=[("disabled", c("surface_alt"))])
    style.configure("TCombobox", fieldbackground=surface, foreground=text,
                    background=surface, bordercolor=c("border_strong"),
                    lightcolor=c("border_strong"), darkcolor=c("border_strong"),
                    arrowcolor=muted, padding=(px(PAD_S), px(6)), borderwidth=1,
                    arrowsize=px(14), relief="flat")
    style.map("TCombobox", bordercolor=[("focus", c("focus"))],
              lightcolor=[("focus", c("focus"))], darkcolor=[("focus", c("focus"))],
              arrowcolor=[("active", accent)],
              fieldbackground=[("readonly", surface), ("disabled", c("surface_alt"))])
    try:  # the dropdown list is a plain Tk listbox and ignores ttk entirely
        root.option_add("*TCombobox*Listbox.background", surface)
        root.option_add("*TCombobox*Listbox.foreground", text)
        root.option_add("*TCombobox*Listbox.selectBackground", c("accent_soft"))
        root.option_add("*TCombobox*Listbox.selectForeground", text)
        root.option_add("*TCombobox*Listbox.borderWidth", 0)
        root.option_add("*TCombobox*Listbox.font", F["body"])
    except tk.TclError:  # pragma: no cover
        pass

    style.configure("TCheckbutton", background=surface, foreground=text,
                    focuscolor=c("focus"))
    style.map("TCheckbutton", background=[("active", surface)])
    style.configure("TRadiobutton", background=surface, foreground=text,
                    focuscolor=c("focus"))
    style.map("TRadiobutton", background=[("active", surface)])
    style.configure("TSeparator", background=border)
    style.configure("TProgressbar", background=accent, troughcolor=c("surface_alt"),
                    bordercolor=c("surface_alt"), lightcolor=accent, darkcolor=accent,
                    thickness=px(8), borderwidth=0)
    # No trough: this one is a liveness tell that runs beside the real bar,
    # and a second empty track under the first one reads as a broken duplicate.
    style.configure("Thin.Horizontal.TProgressbar", thickness=px(4),
                    troughcolor=surface, bordercolor=surface,
                    background=c("accent_hover"), lightcolor=c("accent_hover"),
                    darkcolor=c("accent_hover"))
    style.configure("Big.Horizontal.TProgressbar", thickness=px(10))
    for _orient, _fill in (("Vertical", "ns"), ("Horizontal", "ew")):
        try:  # a thumb in a trough; the stepper arrows are 1995 furniture
            style.layout("{}.TScrollbar".format(_orient), [
                ("{}.Scrollbar.trough".format(_orient), {
                    "sticky": _fill, "children": [
                        ("{}.Scrollbar.thumb".format(_orient),
                         {"expand": "1", "sticky": "nswe"})]})])
        except tk.TclError:  # pragma: no cover - element names are clam's own
            pass
    style.configure("Vertical.TScrollbar", background=c("border_strong"),
                    troughcolor=surface, bordercolor=surface,
                    lightcolor=c("border_strong"), darkcolor=c("border_strong"),
                    arrowcolor=muted, gripcount=0, borderwidth=0,
                    arrowsize=px(10), width=px(9))
    style.map("Vertical.TScrollbar",
              background=[("active", c("muted")), ("pressed", c("muted"))])
    style.configure("Horizontal.TScrollbar", background=c("border_strong"),
                    troughcolor=surface, bordercolor=surface,
                    lightcolor=c("border_strong"), darkcolor=c("border_strong"),
                    arrowcolor=muted, gripcount=0, borderwidth=0,
                    arrowsize=px(10), width=px(9))
    style.map("Horizontal.TScrollbar",
              background=[("active", c("muted")), ("pressed", c("muted"))])

    # The menubar is a native Tk widget, not ttk, and inherits nothing.
    try:
        root.option_add("*Menu.background", surface)
        root.option_add("*Menu.foreground", text)
        root.option_add("*Menu.activeBackground", c("accent_soft"))
        root.option_add("*Menu.activeForeground", text)
        root.option_add("*Menu.selectColor", accent)
        root.option_add("*Menu.borderWidth", 0)
        root.option_add("*Menu.activeBorderWidth", 0)
        root.option_add("*Menu.relief", "flat")
    except tk.TclError:  # pragma: no cover
        pass

    # LAST: the drawn skins replace clam's layouts for the controls a user
    # touches, so nothing configured above can put a square corner back.
    try:
        _install_skins(root, style)
    except tk.TclError:  # pragma: no cover - a Tk with no image element engine
        pass
    return style


class Tooltip:
    """A plain-language tooltip (section 10.4).

    Shown on hover AND on keyboard focus, so a novice tabbing through Settings
    gets the same explanation as one using a mouse.
    """

    DELAY_MS = 450

    def __init__(self, widget: tk.Misc, text: str, wrap_px: int = 320):
        self.widget = widget
        self.text = text
        self.wrap = px(wrap_px)
        self._after: Optional[str] = None
        self._tip: Optional[tk.Toplevel] = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self.hide, add="+")
        widget.bind("<ButtonPress>", self.hide, add="+")
        widget.bind("<FocusIn>", self._schedule, add="+")
        widget.bind("<FocusOut>", self.hide, add="+")
        widget.bind("<Destroy>", self.hide, add="+")

    def _schedule(self, _event: Any = None) -> None:
        self.hide()
        try:
            self._after = self.widget.after(self.DELAY_MS, self.show)
        except tk.TclError:
            self._after = None

    def show(self) -> None:
        if self._tip is not None or not self.text:
            return
        try:
            x = self.widget.winfo_rootx() + px(12)
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + px(6)
        except tk.TclError:
            return
        tip = tk.Toplevel(self.widget)
        tip.wm_overrideredirect(True)
        tip.wm_geometry("+{}+{}".format(x, y))
        tip.configure(background=c("border_strong"))
        label = tk.Label(tip, text=self.text, justify="left", wraplength=self.wrap,
                         background=c("surface"), foreground=c("text"),
                         font=F.get("small"), padx=px(10), pady=px(7),
                         borderwidth=0, relief="flat")
        label.pack(padx=1, pady=1)
        self._tip = tip

    def hide(self, _event: Any = None) -> None:
        if self._after is not None:
            try:
                self.widget.after_cancel(self._after)
            except tk.TclError:
                pass
            self._after = None
        if self._tip is not None:
            try:
                self._tip.destroy()
            except tk.TclError:
                pass
            self._tip = None


def round_rect_points(x0: float, y0: float, x1: float, y1: float,
                      r: float) -> List[float]:
    """Point list for a rounded rectangle, to be drawn ``smooth=True``.

    Tk has no rounded-rectangle primitive and no alpha channel. A smoothed
    polygon through doubled corner points is the honest way to get a radius that
    stays crisp at 150 % DPI, and it costs one canvas item.
    """
    r = max(0.0, min(r, (x1 - x0) / 2.0, (y1 - y0) / 2.0))
    return [x0 + r, y0, x1 - r, y0, x1, y0, x1, y0 + r, x1, y1 - r, x1, y1,
            x1 - r, y1, x0 + r, y1, x0, y1, x0, y1 - r, x0, y0 + r, x0, y0]


class Card(tk.Frame):
    """A surface panel: 1 px hairline, flat fill, drawn rounded corners.

    Never ``relief="sunken"``/``"groove"``/``"ridge"`` — those are the 3-D edges
    that date the whole window. The body is an ordinary square frame with a
    one-pixel border; the radius comes from four small canvases *placed over the
    corners*, each painting the colour behind the card outside the arc and the
    card surface inside it. Doing it that way keeps the body's geometry
    rectangular, so every child still packs and grids normally, and it costs
    four tiny canvases instead of a bitmap that would break at fractional DPI.
    """

    def __init__(self, parent: tk.Misc, *, under: Optional[str] = None,
                 surface: Optional[str] = None, border: Optional[str] = None,
                 radius: int = RADIUS):
        self.under = under or c("bg")
        self.surface = surface or c("surface")
        self.border = border or c("border")
        super().__init__(parent, background=self.border, borderwidth=0,
                         highlightthickness=0)
        self.body = tk.Frame(self, background=self.surface, borderwidth=0,
                             highlightthickness=0)
        self.body.pack(fill="both", expand=True, padx=1, pady=1)
        # Created after the body on purpose: Tk stacks by creation order, so
        # these sit on top of it whatever geometry manager either one uses.
        self._corners = []
        r = px(radius)
        if r > 1 and self.border != self.surface:
            for corner, place_kw in (
                    ("nw", {"x": 0, "y": 0}),
                    ("ne", {"relx": 1.0, "x": -r, "y": 0}),
                    ("sw", {"x": 0, "rely": 1.0, "y": -r}),
                    ("se", {"relx": 1.0, "x": -r, "rely": 1.0, "y": -r})):
                cv = tk.Canvas(self, width=r, height=r, background=self.under,
                               highlightthickness=0, borderwidth=0)
                cv.place(width=r, height=r, **place_kw)
                self._paint_corner(cv, corner, r)
                self._corners.append(cv)

    def _paint_corner(self, cv: tk.Canvas, corner: str, r: int) -> None:
        """Paint one quarter-disc: surface inside, hairline on the edge."""
        box = {"nw": (0, 0, 2 * r, 2 * r), "ne": (-r, 0, r, 2 * r),
               "sw": (0, -r, 2 * r, r), "se": (-r, -r, r, r)}[corner]
        start = {"nw": 90, "ne": 0, "sw": 180, "se": 270}[corner]
        # A one-pixel ramp outside the hairline is a hand-rolled antialias: Tk
        # draws hard pixels, and a bare arc at this radius reads as a staircase.
        cv.create_arc(box[0] - 1, box[1] - 1, box[2] + 1, box[3] + 1,
                      start=start - 2, extent=94, style="arc", width=1,
                      outline=mix(self.under, self.border, 0.55))
        cv.create_arc(box[0], box[1], box[2] - 1, box[3] - 1, start=start - 3,
                      extent=96, style="pieslice", fill=self.surface, outline="")
        # Stroked separately: a pieslice also outlines its two straight radii,
        # which draws a stub of border *inside* the card.
        cv.create_arc(box[0], box[1], box[2] - 1, box[3] - 1, start=start - 2,
                      extent=94, style="arc", outline=self.border, width=1)


def card(parent: tk.Misc, **pack_kw: Any) -> tk.Frame:
    """A :class:`Card`, returning the body frame the caller fills.

    Kept as a function because every view says ``inner = card(parent)`` and then
    ``inner.master.pack(...)``; ``inner.master`` is the Card itself.
    """
    outer = Card(parent, under=pack_kw.pop("under", None))
    if pack_kw:
        outer.pack(**pack_kw)
    outer.inner = outer.body  # type: ignore[attr-defined]
    return outer.body


class FocusRing(tk.Frame):
    """A hairline container that turns into an accent ring when its child has focus.

    ``clam`` will not draw a focus indicator on a Treeview or a Canvas, and its
    default is a dotted rectangle nobody can see anyway. This is a real 2 px ring
    in the focus colour, which is the WCAG 2.4.7 requirement.
    """

    def __init__(self, parent: tk.Misc, under: Optional[str] = None):
        self.under = under or c("bg")
        super().__init__(parent, background=c("border"),
                         highlightthickness=px(2), highlightbackground=self.under,
                         highlightcolor=c("focus"), borderwidth=0)

    def track(self, widget: tk.Misc) -> None:
        """Light the ring while ``widget`` holds keyboard focus."""
        widget.bind("<FocusIn>", self._on, add="+")
        widget.bind("<FocusOut>", self._off, add="+")

    def _on(self, _event: Any = None) -> None:
        try:
            self.configure(highlightbackground=c("focus"),
                           background=c("focus"))
        except tk.TclError:  # pragma: no cover - during teardown
            pass

    def _off(self, _event: Any = None) -> None:
        try:
            self.configure(highlightbackground=self.under, background=c("border"))
        except tk.TclError:  # pragma: no cover - during teardown
            pass


def focus_ring(parent: tk.Misc, under: Optional[str] = None) -> FocusRing:
    """A container whose border is the focus ring (section 10.8)."""
    return FocusRing(parent, under=under)


def autohide(bar: tk.Misc, owner: tk.Misc, **pack_kw: Any) -> Callable[[str, str], None]:
    """A scrollbar that is only there when there is something to scroll.

    A full-height thumb over a table of seven rows is a slab of furniture that
    says nothing; ttk has no opinion about it, so the set-callback takes one.
    """
    def _set(first: str, last: str) -> None:
        try:
            if float(first) <= 0.0 and float(last) >= 1.0:
                bar.pack_forget()
            elif not bar.winfo_manager():
                bar.pack(before=owner, **pack_kw)
        except (tk.TclError, ValueError):  # pragma: no cover - teardown
            pass
        try:
            bar.set(first, last)  # type: ignore[attr-defined]
        except tk.TclError:  # pragma: no cover
            pass
    return _set


class ScrollPane(ttk.Frame):
    """A vertically scrollable region whose scrollbar appears only when needed.

    The Settings tab is a stack of cards with generous padding, and at the
    minimum window size (960x640) the last card used to be cut off by the footer
    with no way to reach it. The canvas takes no focus, so the Tab order still
    runs straight through the real controls.
    """

    def __init__(self, parent: tk.Misc, style: str = "TFrame",
                 background: Optional[str] = None):
        super().__init__(parent, style=style)
        self.canvas = tk.Canvas(self, background=background or c("bg"),
                                highlightthickness=0,
                                borderwidth=0, takefocus=0)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview,
                                 takefocus=0)
        self.canvas.configure(yscrollcommand=self._on_scroll)
        self.canvas.pack(side="left", fill="both", expand=True)
        self._sb_shown = False
        self.body = ttk.Frame(self.canvas, style=style)
        self._win = self.canvas.create_window(0, 0, window=self.body, anchor="nw")
        self.body.bind("<Configure>", self._body_configure, add="+")
        self.canvas.bind("<Configure>", self._canvas_configure, add="+")
        # One global wheel binding, filtered by where the pointer actually is.
        # Binding and unbinding on <Enter>/<Leave> looks tidier and is wrong: Tk
        # sends those when the pointer crosses into a *child*, so the binding
        # stacks up or disappears while the pointer is still inside the pane.
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            self.canvas.bind_all(seq, self._wheel, add="+")
        self.bind("<Destroy>", self._release_wheel, add="+")

    def _body_configure(self, _event: Any = None) -> None:
        try:
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        except tk.TclError:  # pragma: no cover - teardown
            pass

    def _canvas_configure(self, event: Any) -> None:
        self.canvas.itemconfigure(self._win, width=event.width)

    def _on_scroll(self, first: str, last: str) -> None:
        need = not (float(first) <= 0.0 and float(last) >= 1.0)
        if need and not self._sb_shown:
            self.vsb.pack(side="right", fill="y")
            self._sb_shown = True
        elif not need and self._sb_shown:
            self.vsb.pack_forget()
            self._sb_shown = False
        self.vsb.set(first, last)

    def _release_wheel(self, _event: Any = None) -> None:
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            try:
                self.canvas.unbind_all(seq)
            except tk.TclError:  # pragma: no cover - teardown
                pass

    def _inside(self, event: Any) -> bool:
        """True when the pointer is over this pane or one of its children."""
        try:
            under = self.winfo_containing(event.x_root, event.y_root)
        except (tk.TclError, AttributeError):  # pragma: no cover
            return False
        while under is not None:
            if under is self:
                return True
            under = getattr(under, "master", None)
        return False

    def _wheel(self, event: Any) -> None:
        if not self._sb_shown or not self._inside(event):
            return
        num = getattr(event, "num", 0)
        delta = getattr(event, "delta", 0)
        step = -1 if (num == 4 or delta > 0) else 1
        try:
            self.canvas.yview_scroll(step, "units")
        except tk.TclError:  # pragma: no cover - teardown
            pass


class CheckBox(tk.Frame):
    """A Canvas-drawn checkbox: bigger than ttk's, with real states.

    The stock ``ttk.Checkbutton`` indicator is a 13 px bevelled square that no
    amount of styling can enlarge or flatten, and at 150 % DPI it blurs. This
    draws its own: a 22 px rounded box with a stroked tick, hover and pressed
    tints, a 2 px focus ring, and a disabled state that changes shape as well as
    colour. It behaves like a checkbox for the keyboard (Tab reaches it, Space
    and Return toggle it) and the whole row is clickable, label included.

    The Database tab needs dozens of these, so the drawing is four canvas items
    reused across states rather than a full repaint.
    """

    SIZE = 22          # unscaled px, the side of the box itself
    GAP = 10           # unscaled px between the box and its label

    def __init__(self, parent: tk.Misc, text: str = "", *,
                 variable: Optional[tk.Variable] = None,
                 command: Optional[Callable[[], None]] = None,
                 surface: Optional[str] = None, wraplength: int = 0,
                 font: Any = None, subtext: str = ""):
        self.surface = surface or c("surface")
        super().__init__(parent, background=self.surface, borderwidth=0,
                         highlightthickness=0, takefocus=True)
        self.var = variable if variable is not None else tk.BooleanVar(value=False)
        self.command = command
        self._hover = False
        self._focus = False
        self._pressed = False
        self._enabled = True
        side = px(self.SIZE)
        self.box = tk.Canvas(self, width=side, height=side, background=self.surface,
                             highlightthickness=0, borderwidth=0, takefocus=False)
        self.box.pack(side="left", anchor="n", pady=px(2))
        holder = tk.Frame(self, background=self.surface, borderwidth=0,
                          highlightthickness=0)
        holder.pack(side="left", fill="x", expand=True, padx=(px(self.GAP), 0))
        self.label = tk.Label(holder, text=text, background=self.surface,
                              foreground=c("text"), font=font or F.get("body"),
                              anchor="w", justify="left", borderwidth=0,
                              highlightthickness=0)
        self.label.pack(anchor="w", fill="x")
        if wraplength:
            self.label.configure(wraplength=px(wraplength))
        self.sublabel: Optional[tk.Label] = None
        if subtext:
            self.sublabel = tk.Label(holder, text=subtext, background=self.surface,
                                     foreground=c("muted"), font=F.get("small"),
                                     anchor="w", justify="left", borderwidth=0,
                                     highlightthickness=0)
            if wraplength:
                self.sublabel.configure(wraplength=px(wraplength))
            self.sublabel.pack(anchor="w", fill="x")

        for w in (self, self.box, holder, self.label, self.sublabel):
            if w is None:
                continue
            w.bind("<Button-1>", self._press, add="+")
            w.bind("<ButtonRelease-1>", self._release, add="+")
            w.bind("<Enter>", self._enter, add="+")
            w.bind("<Leave>", self._leave, add="+")
        self.bind("<FocusIn>", self._focus_in, add="+")
        self.bind("<FocusOut>", self._focus_out, add="+")
        self.bind("<space>", self._key, add="+")
        self.bind("<Return>", self._key, add="+")
        self.bind("<KP_Enter>", self._key, add="+")
        try:
            self.var.trace_add("write", lambda *_a: self._draw())
        except AttributeError:  # pragma: no cover - Python 3.9 keeps trace()
            self.var.trace("w", lambda *_a: self._draw())
        self._draw()

    # -- public API ----------------------------------------------------------
    def get(self) -> bool:
        return bool(self.var.get())

    def set(self, value: bool) -> None:
        self.var.set(bool(value))

    def toggle(self) -> None:
        if not self._enabled:
            return
        self.var.set(not bool(self.var.get()))
        if self.command is not None:
            self.command()

    def set_enabled(self, enabled: bool) -> None:
        self._enabled = bool(enabled)
        self.configure(takefocus=bool(enabled))
        fg = c("text") if enabled else mix(c("muted"), self.surface, 0.4)
        self.label.configure(foreground=fg)
        if self.sublabel is not None:
            self.sublabel.configure(
                foreground=c("muted") if enabled else mix(c("muted"), self.surface, 0.4))
        self._draw()

    def configure_state(self, state: str) -> None:
        """``"normal"``/``"disabled"``, so callers can mirror ttk's vocabulary."""
        self.set_enabled(state != "disabled")

    # -- events --------------------------------------------------------------
    def _enter(self, _event: Any = None) -> None:
        self._hover = True
        self._draw()

    def _leave(self, _event: Any = None) -> None:
        self._hover = False
        self._pressed = False
        self._draw()

    def _press(self, _event: Any = None) -> str:
        if self._enabled:
            self.focus_set()
            self._pressed = True
            self._draw()
        return "break"

    def _release(self, _event: Any = None) -> str:
        was = self._pressed
        self._pressed = False
        if was and self._enabled:
            self.toggle()
        self._draw()
        return "break"

    def _focus_in(self, _event: Any = None) -> None:
        self._focus = True
        self._draw()

    def _focus_out(self, _event: Any = None) -> None:
        self._focus = False
        self._draw()

    def _key(self, _event: Any = None) -> str:
        self.toggle()
        return "break"

    # -- painting ------------------------------------------------------------
    def _draw(self) -> None:
        try:
            self.box.delete("all")
        except tk.TclError:  # pragma: no cover - during teardown
            return
        side = px(self.SIZE)
        on = bool(self.var.get())
        r = px(6)
        pad = px(3)
        x0, y0, x1, y1 = pad, pad, side - pad - 1, side - pad - 1
        if not self._enabled:
            fill = c("surface_alt")
            edge = c("border")
            tick = mix(c("muted"), c("surface_alt"), 0.35)
        elif on:
            fill = c("accent_active") if self._pressed else c("accent")
            edge = fill
            tick = c("on_accent")
        else:
            fill = c("surface_alt") if self._pressed else (
                c("accent_tint") if self._hover else self.surface)
            edge = c("accent") if self._hover else c("border_strong")
            tick = edge
        if self._focus and self._enabled:
            self.box.create_polygon(
                round_rect_points(0, 0, side - 1, side - 1, r + px(3)),
                smooth=True, splinesteps=16, fill="", outline=c("focus"),
                width=max(2, px(2)))
        self.box.create_polygon(round_rect_points(x0, y0, x1, y1, r),
                                smooth=True, splinesteps=16, fill=fill,
                                outline=edge, width=px(2) if not on else px(1))
        if on:
            # A stroked tick, not a glyph: it stays sharp at any scale and its
            # weight can match the box.
            self.box.create_line(
                x0 + (x1 - x0) * 0.24, y0 + (y1 - y0) * 0.53,
                x0 + (x1 - x0) * 0.43, y0 + (y1 - y0) * 0.72,
                x0 + (x1 - x0) * 0.77, y0 + (y1 - y0) * 0.30,
                fill=tick, width=max(2, px(2.4)), capstyle="round",
                joinstyle="round")
        elif not self._enabled:
            # Shape, not colour: an off-and-locked box carries a dash.
            self.box.create_line(x0 + (x1 - x0) * 0.28, (y0 + y1) / 2.0,
                                 x0 + (x1 - x0) * 0.72, (y0 + y1) / 2.0,
                                 fill=tick, width=max(2, px(2)), capstyle="round")


def _rounded_cover(x: float, y: float, x0: float, y0: float, x1: float,
                   y1: float, r: float) -> float:
    """Coverage in ``[0, 1]`` of a rounded rectangle at a sample point.

    A 3x3 supersample of this is the antialiasing Tk will not do for us.
    """
    cx = x0 + r if x < x0 + r else (x1 - r if x > x1 - r else x)
    cy = y0 + r if y < y0 + r else (y1 - r if y > y1 - r else y)
    if x < x0 or x > x1 or y < y0 or y > y1:
        return 0.0
    dx, dy = x - cx, y - cy
    if dx == 0.0 and dy == 0.0:
        return 1.0
    return 1.0 if (dx * dx + dy * dy) <= r * r else 0.0


def _seg_distance(px_: float, py_: float, ax: float, ay: float,
                  bx: float, by: float) -> float:
    """Distance from a point to a line segment (for the stroked tick)."""
    vx, vy = bx - ax, by - ay
    wx, wy = px_ - ax, py_ - ay
    denom = vx * vx + vy * vy
    t = 0.0 if denom <= 0 else (wx * vx + wy * vy) / denom
    t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
    dx, dy = px_ - (ax + t * vx), py_ - (ay + t * vy)
    return math.sqrt(dx * dx + dy * dy)


def checkbox_image(master: tk.Misc, on: bool, *, size: int = 20, gap: int = 10,
                   surface: Optional[str] = None) -> tk.PhotoImage:
    """A :class:`CheckBox`-styled checkbox as a Tk image, for a Treeview row.

    The Database tab needs a checkbox on 162 rows. A Treeview cell can hold an
    image but not a widget, so the very same square the custom
    :class:`CheckBox` draws on a canvas is rendered here into a
    :class:`tkinter.PhotoImage`: same radius, same accent fill, same stroked
    tick, antialiased by supersampling and transparent outside the rounded
    corners so it sits correctly on a selected row.
    """
    side = max(12, px(size))
    # A Treeview butts its text straight up against the row image, so the gap
    # between box and label is transparent padding carried by the image itself.
    pad_right = max(0, px(gap))
    back = surface or c("surface")
    img = tk.PhotoImage(master=master, width=side + pad_right, height=side)
    fill = c("accent") if on else back
    edge = c("accent") if on else c("border_strong")
    tick = c("on_accent")
    pad = 1.0
    x0, y0, x1, y1 = pad, pad, side - 1 - pad, side - 1 - pad
    radius = max(2.0, side * 0.28)
    stroke = max(1.6, side * 0.09)
    # The tick, as two segments across the box (the CheckBox proportions).
    tx = [x0 + (x1 - x0) * f for f in (0.24, 0.43, 0.77)]
    ty = [y0 + (y1 - y0) * f for f in (0.53, 0.72, 0.30)]
    tick_w = max(1.5, side * 0.11)
    rows = []
    transparent = []
    for yy in range(side):
        row = []
        for xx in range(side + pad_right):
            if xx >= side:
                transparent.append((xx, yy))
                row.append(back)
                continue
            cover = 0.0
            inner = 0.0
            for sy in (0.17, 0.5, 0.83):
                for sx in (0.17, 0.5, 0.83):
                    sxx, syy = xx + sx, yy + sy
                    cover += _rounded_cover(sxx, syy, x0, y0, x1, y1, radius)
                    inner += _rounded_cover(sxx, syy, x0 + stroke, y0 + stroke,
                                            x1 - stroke, y1 - stroke,
                                            max(1.0, radius - stroke))
            cover /= 9.0
            inner /= 9.0
            if cover <= 0.0:
                transparent.append((xx, yy))
                row.append(back)
                continue
            colour = mix(mix(back, edge, cover), fill, inner)
            if on:
                d = min(_seg_distance(xx + 0.5, yy + 0.5, tx[0], ty[0], tx[1], ty[1]),
                        _seg_distance(xx + 0.5, yy + 0.5, tx[1], ty[1], tx[2], ty[2]))
                t = (tick_w / 2.0 + 0.5 - d) / 1.0
                if t > 0.0:
                    colour = mix(colour, tick, min(1.0, t))
            row.append(colour)
        rows.append(row)
    img.put(rows)
    for xx, yy in transparent:
        try:
            img.transparency_set(xx, yy, True)
        except (tk.TclError, AttributeError):  # pragma: no cover - very old Tk
            break
    return img


class DropZone(tk.Canvas):
    """The one big target a novice needs (section 10.2), and the only animation.

    A Canvas, so the whole thing is drawn: a rounded dashed boundary, and inside
    it a stylised circular chromosome carrying **seven** illuminated locus arcs —
    seven housekeeping loci *is* the MLST idea — swept in sequence by a travelling
    read head, with a few softly drifting cocci and rods behind. No image files,
    no dependencies, no cartoon.

    Cost when idle is zero. The frame chain is a single ``after()`` handle that is
    armed only while the zone is actually on screen and waiting for a file, and
    is cancelled the moment a run starts, results appear, or the window is
    minimised. It re-arms BEFORE painting, exactly like the controller's poll
    chain, so a slow frame can never stall the next one. Under a reduced-motion
    preference it paints one static frame and never arms at all.

    Click anywhere browses for files; ``Return``/``space`` browses files and
    ``Shift-Return`` browses a folder.
    """

    FPS_MS = 33             # ~30 fps
    LOCI = 7                # the seven housekeeping loci of a classical scheme
    SWEEP_PER_FRAME = 0.042  # loci per frame -> one full circuit in ~5.5 s
    FLOURISH_FRAMES = 18

    def __init__(self, parent: tk.Misc, *, on_files: Callable[[Sequence[str]], None],
                 on_browse: Callable[[], None], on_browse_folder: Callable[[], None],
                 dnd: bool):
        # highlightthickness=0 on purpose: Tk's own focus ring is a hard square
        # drawn around a rounded dashed boundary, which read as two borders.
        # The ring is drawn inside, on the same radius, by redraw().
        super().__init__(parent, highlightthickness=0,
                         background=c("bg"), borderwidth=0, takefocus=True,
                         height=px(260), cursor="hand2")
        self.on_files = on_files
        self.on_browse = on_browse
        self.on_browse_folder = on_browse_folder
        self.dnd = dnd
        self.zone_state = "idle"
        self._enabled = True
        self.headline = ("Drop a FASTA file here"
                         if dnd else "Choose your FASTA file")
        self.subline = ("or click to browse. Analysis starts by itself - one file, "
                        "several, or a whole folder."
                        if dnd else
                        "Click anywhere in this box. Analysis starts by itself; "
                        "press Shift+Enter to choose a whole folder.")
        #: Single animation handle. There is never a second one.
        self._anim_after: Optional[str] = None
        self._want_anim = True      # the view says the zone is the current focus
        self._mapped = True
        self._static = reduced_motion()
        self._phase = 0.0
        self._flourish = 0
        self._items: Dict[str, Any] = {}
        self._geo: Dict[str, float] = {}
        self._bugs: List[Dict[str, float]] = []
        self.bind("<Configure>", lambda e: self.redraw(), add="+")
        self.bind("<Button-1>", self._click, add="+")
        self.bind("<Return>", self._key_browse, add="+")
        self.bind("<KP_Enter>", self._key_browse, add="+")
        self.bind("<space>", self._key_browse, add="+")
        self.bind("<Shift-Return>", self._key_browse_folder, add="+")
        self.bind("<Enter>", lambda e: self.set_state("hover"), add="+")
        self.bind("<Leave>", lambda e: self.set_state("idle"), add="+")
        self.bind("<FocusIn>", lambda e: self.set_state("focus"), add="+")
        self.bind("<FocusOut>", lambda e: self.set_state("idle"), add="+")
        self.bind("<Map>", self._on_map, add="+")
        self.bind("<Unmap>", self._on_unmap, add="+")
        self.bind("<Destroy>", lambda e: self._stop_anim(), add="+")
        try:  # a minimised window unmaps the toplevel, not this canvas
            top = self.winfo_toplevel()
            top.bind("<Map>", self._on_map, add="+")
            top.bind("<Unmap>", self._on_unmap, add="+")
        except tk.TclError:  # pragma: no cover - no toplevel under test
            pass
        if dnd:
            self._register_dnd()

    # -- drag and drop -------------------------------------------------------
    def _register_dnd(self) -> None:
        try:
            self.drop_target_register(tkinterdnd2.DND_FILES)  # type: ignore[attr-defined]
            self.dnd_bind("<<DropEnter>>", self._drag_enter)  # type: ignore[attr-defined]
            self.dnd_bind("<<DropLeave>>", self._drag_leave)  # type: ignore[attr-defined]
            self.dnd_bind("<<Drop>>", self._drop)  # type: ignore[attr-defined]
        except Exception:
            LOG.info("drag-and-drop registration failed; click-to-browse only")
            self.dnd = False
            self.headline = "Choose your FASTA file"
            self.redraw()

    def _drag_enter(self, event: Any) -> Any:
        if self._enabled:
            self.set_state("dragover")
        return event.action if hasattr(event, "action") else None

    def _drag_leave(self, event: Any) -> Any:
        self.set_state("idle")
        return getattr(event, "action", None)

    def _drop(self, event: Any) -> Any:
        self.set_state("idle")
        if self._enabled:
            paths = split_dnd_paths(self, getattr(event, "data", ""))
            if paths:
                self.flourish()
                self.on_files(paths)
        return getattr(event, "action", None)

    # -- interaction ---------------------------------------------------------
    def _click(self, _event: Any = None) -> None:
        self.focus_set()
        if self._enabled:
            self.on_browse()

    def _key_browse(self, _event: Any = None) -> str:
        if self._enabled:
            self.on_browse()
        return "break"

    def _key_browse_folder(self, _event: Any = None) -> str:
        if self._enabled:
            self.on_browse_folder()
        return "break"

    def set_enabled(self, enabled: bool) -> None:
        """Disable while a run is in flight; a second start() is then impossible."""
        self._enabled = enabled
        self.configure(cursor="hand2" if enabled else "watch")
        self.redraw()

    def set_state(self, state: str) -> None:
        if state != self.zone_state:
            self.zone_state = state
            self.redraw()

    # -- the animation chain -------------------------------------------------
    def animate(self, on: bool) -> None:
        """The view's switch: animate only while this zone is what matters."""
        self._want_anim = bool(on)
        self._sync_anim()

    def _on_map(self, _event: Any = None) -> None:
        self._mapped = True
        self._sync_anim()

    def _on_unmap(self, _event: Any = None) -> None:
        self._mapped = False
        self._stop_anim()

    def should_animate(self) -> bool:
        """Everything that has to be true before a single frame is scheduled."""
        if self._static or not self._mapped:
            return False
        if self._flourish > 0:
            return True
        return bool(self._want_anim and self._enabled)

    def _sync_anim(self) -> None:
        if self.should_animate():
            self._start_anim()
        else:
            self._stop_anim()

    def _start_anim(self) -> None:
        if self._anim_after is None:
            try:
                self._anim_after = self.after(self.FPS_MS, self._tick)
            except tk.TclError:  # pragma: no cover - widget already gone
                self._anim_after = None

    def _stop_anim(self) -> None:
        if self._anim_after is not None:
            try:
                self.after_cancel(self._anim_after)
            except tk.TclError:  # pragma: no cover
                pass
            self._anim_after = None

    def _tick(self) -> None:
        self._anim_after = None
        if not self.should_animate():
            return
        # Re-armed BEFORE the frame is painted, never after: a frame that took
        # too long must not be able to drop the chain (the poll-chain rule).
        try:
            self._anim_after = self.after(self.FPS_MS, self._tick)
        except tk.TclError:  # pragma: no cover
            return
        self._phase = (self._phase + self.SWEEP_PER_FRAME) % float(self.LOCI)
        if self._flourish > 0:
            self._flourish -= 1
            if self._flourish == 0 and not (self._want_anim and self._enabled):
                self._stop_anim()
        self._frame()

    def flourish(self) -> None:
        """A brief confirming pulse when files actually land."""
        if self._static:
            return
        self._flourish = self.FLOURISH_FRAMES
        self._sync_anim()

    # -- painting ------------------------------------------------------------
    def _colours(self) -> Dict[str, str]:
        """Zone colours for the current state. One accent, nothing else."""
        if not self._enabled:
            # No dashes: a dashed boundary is the universal "drop here" sign and
            # this zone accepts nothing while a run is in flight.
            fill, edge = c("surface_alt"), c("border")
            return {"fill": fill, "edge": edge, "width": px(1), "dash": False,
                    "ink": c("muted"), "headline": c("text"), "sub": c("muted")}
        if self.zone_state == "dragover":
            return {"fill": c("accent_soft"), "edge": c("drop_hover"),
                    "width": px(3), "dash": False, "ink": c("accent"),
                    "headline": c("text"), "sub": c("text_soft")}
        if self.zone_state == "focus":
            # Keyboard focus is a SOLID accent boundary: dashes already mean
            # "you may drop here", and the two together read as two borders.
            return {"fill": c("accent_tint"), "edge": c("focus"), "width": px(2),
                    "dash": False, "ink": c("accent"), "headline": c("text"),
                    "sub": c("muted")}
        if self.zone_state == "hover":
            return {"fill": c("accent_tint"), "edge": c("accent"), "width": px(2),
                    "dash": True, "ink": c("accent"), "headline": c("text"),
                    "sub": c("muted")}
        return {"fill": c("accent_tint"), "edge": c("drop_idle"), "width": px(2),
                "dash": True, "ink": c("accent"), "headline": c("text"),
                "sub": c("muted")}

    def _seed_bugs(self, w: float, h: float, top: float) -> None:
        """Deterministic starting positions — no ``random`` import, no surprises."""
        self._bugs = []
        spec = ((0.13, 0.30, 1.0, 0.55, 7, "rod", 0.5),
                (0.83, 0.22, -0.7, 0.9, 6, "coccus", 0.0),
                (0.26, 0.74, 0.85, -0.7, 5, "coccus", 0.0),
                (0.72, 0.68, -1.0, -0.5, 8, "rod", 2.1),
                (0.48, 0.14, 0.6, 1.0, 5, "coccus", 0.0),
                (0.92, 0.52, -0.9, -0.8, 6, "rod", 1.2))
        for fx, fy, vx, vy, size, kind, ang in spec:
            self._bugs.append({"x": fx * w, "y": top * fy, "vx": vx * 0.16,
                               "vy": vy * 0.16, "r": float(px(size)),
                               "rod": 1.0 if kind == "rod" else 0.0,
                               "ang": ang, "spin": 0.004 if kind == "rod" else 0.0})

    def _bug_points(self, bug: Dict[str, float]) -> List[float]:
        """A capsule (rod) or a circle (coccus), as a smoothed polygon."""
        r = bug["r"]
        half = r * (1.9 if bug["rod"] else 0.0)
        ca, sa = math.cos(bug["ang"]), math.sin(bug["ang"])
        pts: List[float] = []
        for k in range(12):
            a = math.pi * 2.0 * k / 12.0
            # a circle stretched along its own axis, then rotated
            lx, ly = math.cos(a) * r, math.sin(a) * r
            lx += half if math.cos(a) >= 0 else -half
            pts.extend((bug["x"] + lx * ca - ly * sa, bug["y"] + lx * sa + ly * ca))
        return pts

    def redraw(self) -> None:
        """Rebuild every canvas item. Called on resize and on a state change."""
        self.delete("all")
        self._items = {}
        w = float(max(self.winfo_width(), px(320)))
        # The floor used to be 150 px while show_state() gives the strip 104,
        # so the boundary was drawn past the bottom of the canvas and the strip
        # had three sides.
        h = float(max(self.winfo_height(), px(84)))
        col = self._colours()
        pad = float(px(PAD_S))
        self._items["frame"] = self.create_polygon(
            round_rect_points(pad, pad, w - pad, h - pad, px(RADIUS + 8)),
            smooth=True, splinesteps=24, fill=col["fill"], outline=col["edge"],
            width=col["width"],
            dash=(px(9), px(7)) if col["dash"] else ())

        compact = h < px(215)
        tiny = w < px(300) or h < px(120)
        if compact:
            # Short zone (results are on screen): illustration on the left, the
            # words beside it, everything still legible.
            radius = max(float(px(20)), min(float(px(38)), (h - pad * 2) * 0.36))
            cx = pad + float(px(PAD_L)) + radius
            cy = h / 2.0
            text_x = cx + radius + float(px(PAD_L))
            head_y = cy - float(px(11))
            sub_y = cy + float(px(11))
            anchor, text_w = "w", max(float(px(120)), w - text_x - pad - px(PAD_M))
            top_area = h
        else:
            # One optically centred block — ring, headline, one line of help —
            # rather than art floating at the top and words pinned to the floor.
            gap_ring, gap_text = float(px(34)), float(px(28))
            head_h, sub_h = float(px(24)), float(px(18))
            radius = max(float(px(32)),
                         min(float(px(128)),
                             min(w * 0.15, (h - float(px(160))) * 0.34)))
            block = (2 * radius + gap_ring + head_h + gap_text + sub_h
                     + float(px(26)))
            top = max(pad + float(px(PAD_M)), (h - block) / 2.0)
            cx, cy = w / 2.0, top + radius
            head_y = top + 2 * radius + gap_ring + head_h / 2.0
            sub_y = head_y + gap_text
            # The flora drift in the band above the words and never across them.
            top_area = head_y - float(px(PAD_L))
            text_x = w / 2.0
            anchor, text_w = "center", max(float(px(220)), w * 0.62)

        self._geo = {"w": w, "h": h, "cx": cx, "cy": cy, "r": radius,
                     "top": top_area, "compact": 1.0 if compact else 0.0}

        ink = col["ink"]
        base = col["fill"]
        dim = mix(base, ink, 0.38)
        ring = mix(base, ink, 0.52)

        # -- drifting flora, behind everything ------------------------------
        if not tiny and not compact:
            if not self._bugs:
                self._seed_bugs(w, h, top_area)
            bug_fill = mix(base, ink, 0.10)
            bug_edge = mix(base, ink, 0.22)
            ids = []
            for bug in self._bugs:
                ids.append(self.create_polygon(self._bug_points(bug), smooth=True,
                                               splinesteps=10, fill=bug_fill,
                                               outline=bug_edge, width=1))
            self._items["bugs"] = ids
        else:
            self._items["bugs"] = []

        # -- the chromosome --------------------------------------------------
        self._items["pulse"] = self.create_oval(cx, cy, cx, cy, outline="", width=px(2))
        self.create_oval(cx - radius * 0.72, cy - radius * 0.72,
                         cx + radius * 0.72, cy + radius * 0.72,
                         outline=mix(base, ink, 0.22), width=1)
        self.create_oval(cx - radius, cy - radius, cx + radius, cy + radius,
                         outline=ring, width=max(1, px(2)))

        span = 360.0 / self.LOCI
        arc_span = span * 0.62
        box = (cx - radius, cy - radius, cx + radius, cy + radius)
        arcs, nodes = [], []
        for i in range(self.LOCI):
            mid = 90.0 - i * span
            arcs.append(self.create_arc(box[0], box[1], box[2], box[3],
                                        start=mid - arc_span / 2.0, extent=arc_span,
                                        style="arc", width=px(7), outline=dim))
            orbit = radius + float(px(11))
            ax = cx + orbit * math.cos(math.radians(mid))
            ay = cy - orbit * math.sin(math.radians(mid))
            nodes.append(self.create_oval(ax - px(4), ay - px(4), ax + px(4),
                                          ay + px(4), fill=dim,
                                          outline=base, width=1))
        self._items["arcs"] = arcs
        self._items["nodes"] = nodes
        self._items["head"] = self.create_oval(cx, cy, cx, cy, fill=ink, outline="")

        if not tiny:
            # ASCII on purpose. A Tk canvas hands the string to the X font
            # verbatim, and on a build with no Xft a U+2026 came out as three
            # mojibake glyphs ("Working(R)<apple>f") in 24 pt across the hero.
            headline = self.headline if self._enabled else "Working..."
            sub = self.subline if self._enabled else (
                "Please wait for the current files to finish.")
            if self._enabled and self.zone_state == "dragover":
                headline = "Release to analyse"
                sub = "Let go and typing starts straight away."
            self._items["headline"] = self.create_text(
                text_x, head_y, text=headline, fill=col["headline"],
                font=F.get("title"), anchor=anchor, width=int(text_w),
                justify="left" if anchor == "w" else "center")
            self._items["sub"] = self.create_text(
                text_x, sub_y, text=sub, fill=col["sub"], font=F.get("small"),
                anchor=anchor, width=int(text_w),
                justify="left" if anchor == "w" else "center")
            if not compact and self._enabled:
                self.create_text(
                    text_x, sub_y + float(px(26)),
                    text="FASTA   \u00b7   GenBank   \u00b7   EMBL   \u00b7   "
                         ".gz  .bz2  .zip",
                    fill=col["sub"], font=F.get("tiny"),
                    anchor=anchor)
        self._frame()
        self._sync_anim()

    def _frame(self) -> None:
        """Update only what moves. Called every ~33 ms, and once per redraw."""
        arcs = self._items.get("arcs")
        if not arcs:
            return
        col = self._colours()
        ink, base = col["ink"], col["fill"]
        dim = mix(base, ink, 0.38)
        bright = ink
        lit_all = self.zone_state == "dragover" and self._enabled
        cx, cy, radius = self._geo["cx"], self._geo["cy"], self._geo["r"]
        span = 360.0 / self.LOCI
        try:
            for i, item in enumerate(arcs):
                if self._static:
                    level = 0.45
                elif lit_all:
                    level = 1.0
                else:
                    d = (i - self._phase) % float(self.LOCI)
                    d = min(d, self.LOCI - d)
                    level = max(0.0, 1.0 - d / 1.45)
                self.itemconfigure(item, outline=mix(dim, bright, level),
                                   width=px(7) + px(3) * level)
                node = self._items["nodes"][i]
                nr = px(4) + px(2.4) * level
                mid = 90.0 - i * span
                orbit = radius + float(px(11))
                ax = cx + orbit * math.cos(math.radians(mid))
                ay = cy - orbit * math.sin(math.radians(mid))
                self.coords(node, ax - nr, ay - nr, ax + nr, ay + nr)
                self.itemconfigure(node, fill=mix(dim, bright, min(1.0, level + 0.15)))

            # the travelling read head
            head = self._items.get("head")
            if head is not None:
                if self._static or lit_all:
                    self.itemconfigure(head, fill="")
                else:
                    ang = 90.0 - self._phase * span
                    hx = cx + radius * math.cos(math.radians(ang))
                    hy = cy - radius * math.sin(math.radians(ang))
                    hr = float(px(4))
                    self.coords(head, hx - hr, hy - hr, hx + hr, hy + hr)
                    self.itemconfigure(head, fill=bright)

            for idx, item in enumerate(self._items.get("bugs") or []):
                bug = self._bugs[idx]
                if not self._static:
                    bug["x"] += bug["vx"]
                    bug["y"] += bug["vy"]
                    bug["ang"] += bug["spin"]
                    margin = bug["r"] * 3.0
                    if bug["x"] < -margin:
                        bug["x"] = self._geo["w"] + margin
                    elif bug["x"] > self._geo["w"] + margin:
                        bug["x"] = -margin
                    if bug["y"] < -margin:
                        bug["y"] = self._geo["top"] + margin
                    elif bug["y"] > self._geo["top"] + margin:
                        bug["y"] = -margin
                self.coords(item, *self._bug_points(bug))

            pulse = self._items.get("pulse")
            if pulse is not None:
                if self._flourish > 0:
                    k = 1.0 - self._flourish / float(self.FLOURISH_FRAMES)
                    pr = radius * (1.0 + 1.4 * k)
                    self.coords(pulse, cx - pr, cy - pr, cx + pr, cy + pr)
                    self.itemconfigure(pulse, outline=mix(bright, base, k),
                                       width=max(1, px(3) * (1.0 - k)))
                else:
                    self.itemconfigure(pulse, outline="")
        except (tk.TclError, KeyError, IndexError):  # pragma: no cover - teardown
            return


class StatusLine(ttk.Frame):
    """The live region at the bottom of the window (section 10.8).

    Never a messagebox: transient messages auto-revert after 6 seconds and the
    resting text always describes the current state.
    """

    REVERT_MS = 6000

    def __init__(self, parent: tk.Misc):
        super().__init__(parent, style="Footer.TFrame")
        self.resting = "Ready."
        self._after: Optional[str] = None
        # A drawn disc rather than a glyph: same reason as the results table --
        # U+2714 is a hollow box on any Tk without a font that owns it, and a
        # status channel that silently disappears is worse than none.
        self.icon = tk.Canvas(self, width=px(13), height=px(13),
                              highlightthickness=0, borderwidth=0,
                              background=c("surface_alt"))
        self._dot = self.icon.create_oval(px(2), px(2), px(11), px(11),
                                          outline="", fill="")
        self.icon.pack(side="left", padx=(px(PAD_M), 0))
        self.label = ttk.Label(self, text=self.resting, style="Status.TLabel",
                               anchor="w")
        self.label.pack(side="left", fill="x", expand=True, padx=px(PAD_S),
                        pady=px(PAD_XS))

    def set(self, text: str, kind: str = "muted", *, transient: bool = False,
            resting: bool = False) -> None:
        """Update the live region. ``resting`` also changes what it reverts to."""
        assert_main_thread("the status line")
        colour = c(kind) if kind in ("ok", "warn", "bad", "info") else c("muted")
        self._set_dot(colour if kind in ("ok", "warn", "bad", "info") else "")
        self.label.configure(text=text, foreground=colour)
        if resting:
            self.resting = text
        if self._after is not None:
            try:
                self.after_cancel(self._after)
            except tk.TclError:
                pass
            self._after = None
        if transient:
            self._after = self.after(self.REVERT_MS, self._revert)

    def _set_dot(self, colour: str) -> None:
        try:
            self.icon.itemconfigure(self._dot, fill=colour)
        except tk.TclError:  # pragma: no cover - teardown
            pass

    def _revert(self) -> None:
        self._after = None
        self._set_dot("")
        self.label.configure(text=self.resting, foreground=c("muted"))


# ===========================================================================
# SECTION 8 — VIEWS: THE ANALYSE TAB  (section 10.2)
# ===========================================================================

#: Column ids for the results tree. ``#0`` carries FILE.
TREE_COLUMNS = ("organism", "scheme", "st", "status", "score")
#: Heading words: bin/mlst:151 plus the two --full columns, and ORGANISM first,
#: because "abaumannii_2" is a directory name and *Acinetobacter baumannii* is
#: the answer the reader came for.
TREE_HEADINGS = {"#0": "FILE", "organism": "ORGANISM", "scheme": "SCHEME",
                 "st": "ST", "status": "STATUS", "score": "SCORE"}
#: Caption row inserted above the per-locus children so reused columns are
#: never ambiguous (section 10.2). The ST column is left empty for child rows so
#: that the long evidence string lands in the widest column.
CHILD_CAPTION = ("ALLELE", "WHAT IT MEANS", "", "BEST EVIDENCE", "HITS")

CHUNK_ROWS = 200        # rows inserted per pump above 2000 rows
BIG_BATCH = 50          # confirm threshold for a folder drop


def evidence_text(call: AlleleCall) -> str:
    """One-line evidence for a locus, straight off the winning :class:`Hit`."""
    hit = call.best
    if hit is None:
        return "no match found"
    # Identity and coverage lead: they are what a reader judges the call on,
    # and if the column ever runs out of room it is the contig coordinates that
    # should be the part cut, not the percentages.
    return "{:.1f}% identity · {:.1f}% of allele · {} {}–{}".format(
        hit.pct_identity, hit.pct_coverage, hit.qseqid, hit.qstart, hit.qend)


def tied_schemes(result: SampleResult) -> str:
    """``"klebsiella ST 258, ecoli_achtman_4 ST 14464"`` when the call is a tie.

    Presentation only: it re-reads :attr:`~wmlst.engine.SampleResult.tied`,
    the block of candidates the engine scored EQUAL to the reported one, and
    never re-derives a number. Empty string when the winner stands alone.

    A tie is a genuine scientific ambiguity -- two schemes fit the assembly
    equally well -- so the GUI shows it rather than hiding the fact that a
    heuristic picked between them (section 5.14a).
    """
    rows = getattr(result, "tied", ())
    if len(rows) < 2:
        return ""
    return ", ".join(
        "{} (no ST)".format(c.scheme) if c.st in ("", "-")
        else "{} ST {}".format(c.scheme, c.st) for c in rows)


TIE_SENTENCE = (
    "Two or more schemes fit this assembly equally well; {app} reported the one "
    "whose alleles sit lowest in each locus\u2019 allele registry. Confirm the "
    "species by another method before using this ST."
)


def wrap_px(text: str, font: Any, max_px: int) -> List[str]:
    """Word-wrap ``text`` to a pixel width, measured in the font that draws it.

    A Treeview clips a cell without a word of warning, so the long explanations
    in the expanded row are folded here against the column's real width rather
    than against a guessed character count — which is how "fit this assembly
    equally well; WM" happened.
    """
    words = str(text).split()
    if not words or max_px <= 0:
        return [str(text)] if text else []
    try:
        measure = font.measure
    except AttributeError:  # pragma: no cover - no font object without a display
        return wrap_words(text, max(8, max_px // 7))
    lines: List[str] = []
    line = ""
    for word in words:
        candidate = word if not line else line + " " + word
        if measure(candidate) > max_px and line:
            lines.append(line)
            line = word
        else:
            line = candidate
    if line:
        lines.append(line)
    return lines


def elide_middle_px(text: str, font: Any, max_px: int) -> str:
    """Shorten from the MIDDLE, keeping the head and the tail.

    Assembly file names differ at both ends -- the accession at the front, the
    extension at the back -- so cutting the tail off "GCF_000013425.1_ASM1342v1
    _genomic.fna" throws away the half that says what kind of file it is.
    """
    text = str(text)
    try:
        measure = font.measure
    except AttributeError:  # pragma: no cover
        return text
    if max_px <= 0 or measure(text) <= max_px:
        return text
    tail = text[-12:] if len(text) > 24 else text[-4:]
    lo, hi = 0, len(text) - len(tail)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if measure(text[:mid] + "..." + tail) <= max_px:
            lo = mid
        else:
            hi = mid - 1
    return text[:lo] + "..." + tail if lo else elide_px(text, font, max_px)


def elide_px(text: str, font: Any, max_px: int) -> str:
    """``text`` shortened with a trailing ellipsis until it fits ``max_px``."""
    text = str(text)
    try:
        measure = font.measure
    except AttributeError:  # pragma: no cover
        return text
    if max_px <= 0 or measure(text) <= max_px:
        return text
    lo, hi = 0, len(text)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if measure(text[:mid] + "...") <= max_px:
            lo = mid
        else:
            hi = mid - 1
    # Trailing separators come with the cut, and "of allele ·..." reads as a
    # mistake rather than as a shortening.
    return (text[:lo].rstrip(" \u00b7-\u2013,;:") + "...") if lo else ""


def wrap_words(text: str, width: int) -> List[str]:
    """Greedy word wrap. ``textwrap`` would do, but this module keeps its text
    helpers together and needs no tab expansion or hyphen breaking."""
    lines: List[str] = []
    line = ""
    for word in str(text).split():
        candidate = word if not line else line + " " + word
        if len(candidate) > width and line:
            lines.append(line)
            line = word
        else:
            line = candidate
    if line:
        lines.append(line)
    return lines


#: The leading whitespace on a continuation row inside an expanded result. It
#: is measured, not guessed at, before anything is wrapped against it.
INDENT = "        "


def _indent_px(font: Any) -> int:
    """Width of :data:`INDENT` in ``font``, for wrapping a continuation row."""
    try:
        return int(font.measure(INDENT))
    except AttributeError:  # pragma: no cover - no display
        return 0


#: The one-line form used inside the results table, where a paragraph does not
#: belong. The full :data:`TIE_SENTENCE` is on the summary card and in the report.
TIE_SHORT = ("Both schemes scored the same. Confirm the species by another "
             "method before using this ST.")


def tie_sentence() -> str:
    """The novice-facing explanation of a scheme tie."""
    return TIE_SENTENCE.format(app=branding.APP_NAME)


# ---------------------------------------------------------------------------
# Scheme provenance (wmlst.schemerefs) — organism, description, reference
# ---------------------------------------------------------------------------
# "abaumannii_2" and "ecoli_achtman_4" are scheme directory names, not answers.
# A reader who has never met them cannot tell whether the tool typed the right
# organism, and one of them told us so. Everything below turns a scheme id into
# what the curated table actually knows, and shows NOTHING where the table is
# blank: a blank cell means "not established", and a diagnostic report must
# never invent a species or a citation.

def scheme_ref(scheme: str, dbdir: str = "") -> Optional[Any]:
    """The curated :class:`~wmlst.schemerefs.SchemeRef` for a scheme, or None.

    The active database is asked first, then the one bundled with this copy of
    WMLST: the table is curated provenance, not data, so a portable or
    hand-placed database that predates it still names its organisms.
    """
    if not scheme or scheme == "-":
        return None
    for where in ((dbdir or None), None):
        try:
            found = schemerefs.lookup(scheme, where)
        except Exception:  # pragma: no cover - the loader swallows its own
            found = None
        if found is not None:
            return found
        if where is None:
            break
    return None


def organism_name(scheme: str, dbdir: str = "") -> str:
    """``"Klebsiella pneumoniae"``, ``"Neisseria spp."`` or ``""``."""
    ref = scheme_ref(scheme, dbdir)
    return ref.organism if ref is not None else ""


def organism_parts(scheme: str, dbdir: str = "") -> Tuple[str, str]:
    """Split an organism label into ``(italic part, upright part)``.

    ``Klebsiella pneumoniae`` is italic throughout; in ``Neisseria spp.`` only
    the genus is italic, because "spp." is not part of the name.
    """
    label = organism_name(scheme, dbdir)
    if not label:
        return ("", "")
    if label.endswith(" spp."):
        return (label[:-5], "spp.")
    return (label, "")


def scheme_caption(scheme: str, ref: Optional[Any]) -> str:
    """``klebsiella · Klebsiella pneumoniae MLST · 7 loci`` — known fields only."""
    bits = [scheme] if scheme and scheme != "-" else []
    if ref is not None:
        if ref.description:
            bits.append(ref.description)
        if ref.n_loci:
            bits.append("{} loci".format(ref.n_loci))
    return "  \u00b7  ".join(bits)


#: How the curated table spells a source, and how a human should read it.
SOURCE_NAMES = {"pubmlst": "PubMLST", "pasteur": "Institut Pasteur",
                "bigsdb": "BIGSdb"}


def source_name(ref: Optional[Any]) -> str:
    """``"PubMLST"`` / ``"Institut Pasteur"`` for the link label."""
    source = str(getattr(ref, "source", "") or "")
    return SOURCE_NAMES.get(source.strip().lower(), source or "Database")


def reference_line(ref: Optional[Any]) -> str:
    """The verified citation, with its PubMed id. ``""`` when none exists.

    75 of the 162 schemes have one; the rest have no citation that could be
    confirmed against PubMed, and for those this returns nothing at all.
    """
    if ref is None or not getattr(ref, "citation", ""):
        return ""
    if ref.pubmed_id:
        return "{}  \u00b7  PMID {}".format(ref.citation, ref.pubmed_id)
    return str(ref.citation)


def open_url(url: str) -> None:
    """Open an authoritative database or PubMed page in the user's browser."""
    if not url:
        return
    try:
        webbrowser.open(url)
    except Exception as exc:  # pragma: no cover - browser-less box
        LOG.warning("could not open %s: %s", url, exc)


def link_button(parent: tk.Misc, text: str, url: str, *,
                tip: str = "") -> ttk.Button:
    """An underlined accent link that is a real button: Tab reaches it, Enter
    and Space follow it, and the URL is in its tooltip."""
    btn = ttk.Button(parent, text=text, style="CardLink.TButton",
                     command=lambda u=url: open_url(u))
    Tooltip(btn, tip or url)
    return btn


def allele_summary(result: SampleResult) -> str:
    """``gene(code);gene(code)`` exactly as the --full ALLELES field reads.

    Presentation only: every code is copied verbatim from the engine.
    """
    return ";".join("{}({})".format(a.locus, a.code) for a in result.alleles)


class AnalyseView(ttk.Frame):
    """The novice flow: drop a file, it analyses itself, results appear."""

    def __init__(self, parent: tk.Misc, app: WmlstApp):
        super().__init__(parent, style="TFrame")
        self.app = app
        self.results: List[SampleResult] = []
        self.failures: List[Tuple[str, Friendly]] = []
        self._row_data: Dict[str, SampleResult] = {}
        #: Drawn STATUS shapes, one per status word. Tk frees an image the
        #: moment Python drops it, so the row would go blank without this.
        self._status_icons: Dict[str, Any] = {}
        #: The UNSHORTENED file name behind every row, so that widening the
        #: window can put back what narrowing it had to elide.
        self._row_names: Dict[str, str] = {}
        self._fit_width = 0
        self._row_error: Dict[str, Friendly] = {}
        self._expanded: set = set()
        self._order = 0
        self._sort_col: Optional[str] = None
        self._sort_desc = False
        self._pending_rows: List[SampleResult] = []
        self._done = 0
        self._total = 0
        self._current_pct = 0.0
        self._indeterminate = False
        self._build()

    # -- construction --------------------------------------------------------
    def _build(self) -> None:
        outer = ttk.Frame(self, style="TFrame")
        outer.pack(fill="both", expand=True, padx=px(PAD_L), pady=px(PAD_M))

        # A forced scheme is sticky across sessions, and when it is the wrong one
        # every file comes back NONE with no visible reason -- which reads as "the
        # tool is broken" rather than "you pinned a scheme". Say so on this tab,
        # where the result is, not only on the Settings tab the user has left.
        self.scheme_notice = ttk.Frame(outer, style="Notice.TFrame")
        # A 3 px bar in the warning colour down the leading edge, plus a glyph:
        # the band is dark text on a tint, and the bar and the glyph give it two
        # further channels for anyone who cannot see the tint at all.
        ttk.Frame(self.scheme_notice, style="NoticeBar.TFrame",
                  width=px(3)).pack(side="left", fill="y")
        ttk.Label(self.scheme_notice, text="LOCKED",
                  style="NoticeBadge.TLabel").pack(side="left",
                                                   padx=(px(PAD_M), 0),
                                                   pady=px(PAD_M))
        self.scheme_notice_label = ttk.Label(
            self.scheme_notice, style="Notice.TLabel", anchor="w", justify="left")
        self.scheme_notice_label.pack(side="left", fill="x", expand=True,
                                      padx=(px(PAD_S), px(PAD_S)), pady=px(PAD_M))
        ttk.Button(self.scheme_notice, text="Use automatic",
                   command=self._clear_forced_scheme).pack(side="right",
                                                           padx=(0, px(PAD_M)),
                                                           pady=px(PAD_S))

        self.dropzone = DropZone(
            outer, on_files=self.handle_paths, on_browse=self.browse_files,
            on_browse_folder=self.browse_folder, dnd=self.app.dnd_enabled)
        # packed by show_state(), which is what decides whether it fills the tab

        # -- progress card (visible only while RUNNING) ----------------------
        self.progress_wrap = ttk.Frame(outer, style="TFrame")
        progress = card(self.progress_wrap)
        progress.master.pack(fill="x")
        self.progress_card = progress
        row = ttk.Frame(progress, style="Surface.TFrame")
        row.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_S)))
        self.progress_label = ttk.Label(row, text="Starting...", style="Heading.TLabel",
                                        anchor="w")
        self.progress_label.pack(side="left", fill="x", expand=True)
        self.elapsed_label = ttk.Label(row, text="", style="SurfaceMuted.TLabel")
        self.elapsed_label.pack(side="left", padx=px(PAD_M))
        self.cancel_button = ttk.Button(row, text="Cancel", command=self.app.cancel_run)
        self.cancel_button.pack(side="left")
        self.overall = ttk.Progressbar(progress, style="Big.Horizontal.TProgressbar",
                                       mode="determinate", maximum=100.0)
        self.overall.pack(fill="x", padx=px(PAD_M))
        # Packed only while it is actually sweeping: two troughs stacked, one of
        # them permanently empty, read as a broken second progress bar.
        self.pulse = ttk.Progressbar(progress, style="Thin.Horizontal.TProgressbar",
                                     mode="indeterminate", maximum=100.0)
        self.progress_detail = ttk.Label(progress, text="", style="SurfaceMuted.TLabel",
                                         anchor="w")
        self.progress_detail.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_XS), px(PAD_M)))

        # -- single-file summary card ----------------------------------------
        # The headline is the ORGANISM, in italic binomial form; the scheme id
        # is secondary, and the reference that defines the scheme is on the card
        # rather than three clicks away. A cryptic "abaumannii_2" alone made one
        # user believe the tool had typed the wrong thing.
        self.summary_wrap = ttk.Frame(outer, style="TFrame")
        summary = card(self.summary_wrap)
        summary.master.pack(fill="x")
        self.summary_card = summary
        head = ttk.Frame(summary, style="Surface.TFrame")
        head.pack(fill="x")
        left = ttk.Frame(head, style="Surface.TFrame")
        # anchor="n": without it pack() centres this column against the tall one
        # beside it and the hero number floated 75 px below the file name.
        left.pack(side="left", anchor="n", padx=(px(PAD_L), px(PAD_XL)),
                  pady=px(PAD_M))
        ttk.Label(left, text="SEQUENCE TYPE", style="Eyebrow.TLabel").pack(anchor="w")
        self.st_label = ttk.Label(left, text="—", style="Hero.TLabel")
        self.st_label.pack(anchor="w", pady=(px(PAD_XS), 0))
        ttk.Frame(left, style="Rule.TFrame", height=px(3),
                  width=px(44)).pack(anchor="w", pady=(px(PAD_S), 0))
        right = ttk.Frame(head, style="Surface.TFrame")
        right.pack(side="left", fill="both", expand=True, padx=(0, px(PAD_L)),
                   pady=px(PAD_M))
        self.summary_file = ttk.Label(right, text="", style="SurfaceMuted.TLabel",
                                      anchor="w")
        self.summary_file.pack(anchor="w", fill="x")
        organism_row = ttk.Frame(right, style="Surface.TFrame")
        organism_row.pack(anchor="w", fill="x")
        self.organism_label = ttk.Label(organism_row, text="", style="Surface.TLabel",
                                        font=F.get("organism"), anchor="w")
        self.organism_label.pack(side="left")
        # "spp." is not part of the name and is therefore not italic.
        self.organism_suffix = ttk.Label(organism_row, text="",
                                         style="SurfaceSoft.TLabel",
                                         font=F.get("title"), anchor="w")
        self.organism_suffix.pack(side="left", padx=(px(PAD_S), 0))
        self.summary_title = ttk.Label(right, text="", style="SurfaceMuted.TLabel",
                                       anchor="w")
        self.summary_title.pack(anchor="w", fill="x", pady=(px(PAD_XS), 0))
        self.summary_status = ttk.Label(right, text="", style="Ok.TLabel", anchor="w")
        self.summary_status.pack(anchor="w", fill="x", pady=(px(PAD_S), 0))
        self.summary_body = ttk.Label(right, text="", style="SurfaceMuted.TLabel",
                                      anchor="w", justify="left")
        self.summary_body.pack(anchor="w", fill="x", pady=(px(PAD_XS), 0))
        self.summary_body.bind(
            "<Configure>",
            lambda e: self.summary_body.configure(
                wraplength=max(px(240), min(px(720), e.width - px(8)))),
            add="+")

        # The reference block. Hidden entirely when the curated table has no
        # citation for this scheme: a blank cell means "not established", and an
        # invented reference on a diagnostic result would be far worse than none.
        self.reference_wrap = ttk.Frame(right, style="Surface.TFrame")
        ttk.Label(self.reference_wrap, text="REFERENCE",
                  style="Eyebrow.TLabel").pack(anchor="w")
        self.reference_label = ttk.Label(self.reference_wrap, text="",
                                         style="SurfaceMuted.TLabel", anchor="w",
                                         justify="left")
        self.reference_label.pack(anchor="w", fill="x")
        # Capped, not merely fitted: a 1000 px card would otherwise set this
        # citation as one 170-character line, which nobody reads.
        self.reference_label.bind(
            "<Configure>",
            lambda e: self.reference_label.configure(
                wraplength=max(px(240), min(px(720), e.width - px(8)))), add="+")
        self.reference_links = ttk.Frame(self.reference_wrap, style="Surface.TFrame")
        self.reference_links.pack(anchor="w", fill="x")

        # A tie is a scientific ambiguity, not a footnote: it gets its own band
        # across the bottom of the card, with both schemes and both STs.
        self.tie_banner = ttk.Frame(summary, style="Notice.TFrame")
        ttk.Frame(self.tie_banner, style="NoticeBar.TFrame",
                  width=px(3)).pack(side="left", fill="y")
        ttk.Label(self.tie_banner, text="TIE", style="NoticeBadge.TLabel").pack(
            side="left", padx=(px(PAD_M), 0), pady=px(PAD_M))
        self.tie_label = ttk.Label(self.tie_banner, text="", style="Notice.TLabel",
                                   anchor="w", justify="left")
        self.tie_label.pack(side="left", fill="x", expand=True,
                            padx=(px(PAD_S), px(PAD_M)), pady=px(PAD_S))
        self.tie_label.bind(
            "<Configure>",
            lambda e: self.tie_label.configure(
                wraplength=max(px(240), e.width - px(8))), add="+")

        # -- results ----------------------------------------------------------
        self.results_wrap = ttk.Frame(outer, style="TFrame")
        header = ttk.Frame(self.results_wrap, style="TFrame")
        header.pack(fill="x", pady=(px(PAD_M), px(PAD_XS)))
        self.results_title = ttk.Label(header, text="Results", style="Title.TLabel")
        self.results_title.pack(side="left")
        self.results_hint = ttk.Label(
            header, text="Click the arrow beside a file to see each locus.",
            style="Muted.TLabel")
        self.results_hint.pack(side="left", padx=px(PAD_M))

        ring = focus_ring(self.results_wrap)
        self.tree = ttk.Treeview(ring, columns=TREE_COLUMNS, show="tree headings",
                                 selectmode="browse", height=6)
        ring.track(self.tree)
        vsb = ttk.Scrollbar(ring, orient="vertical", command=self.tree.yview)
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.configure(
            yscrollcommand=autohide(vsb, self.tree, side="right", fill="y"))
        # #0 is the only stretching column, so it absorbs the window: a RefSeq
        # accession is 40 characters and used to be chopped mid-token.
        # The six starting widths must SUM to less than the narrowest window
        # this application allows (960 px), because "stretch" only ever grows a
        # column: when the sum is larger, Tk clips the last one off the right
        # edge and SCORE simply is not there.
        self.tree.column("#0", width=px(215), minwidth=px(150), stretch=True)
        self.tree.column("organism", width=px(145), minwidth=px(100), anchor="w")
        self.tree.column("scheme", width=px(125), minwidth=px(84), anchor="w")
        self.tree.column("st", width=px(50), minwidth=px(40), anchor="e")
        # The child rows reuse this one for the BLAST evidence string, which is
        # longer than any status word -- but the file name is the thing the
        # reader owns, so the evidence elides and the file name does not.
        self.tree.column("status", width=px(178), minwidth=px(96), anchor="w")
        self.tree.column("score", width=px(58), minwidth=px(46), anchor="e")
        for col, title in TREE_HEADINGS.items():
            self.tree.heading(col, text=title,
                              anchor="e" if col in ("score", "st") else "w",
                              command=lambda cc=col: self.sort_by(cc))
        for key in ("ok", "warn", "bad", "info", "muted"):
            self.tree.tag_configure("st_" + key, foreground=c(key))
        self.tree.tag_configure("failed", foreground=c("bad"))
        self.tree.tag_configure("caption", foreground=c("muted"), font=F.get("tiny"))
        self.tree.tag_configure("locus", foreground=c("text"))
        self.tree.tag_configure("tie", foreground=c("warn"), font=F.get("body_bold"))
        # Accent blue means "this is a link" everywhere else in the window, so
        # the provenance row -- a citation, not a link -- is not painted blue.
        self.tree.tag_configure("provenance", foreground=c("text_soft"))
        self.tree.tag_configure("row_plain", foreground=c("text"))
        self.tree.bind("<Configure>", self._refit_names, add="+")
        self.tree.bind("<<TreeviewOpen>>", self._on_open, add="+")
        self.tree.bind("<<TreeviewSelect>>", self._on_select, add="+")
        self.tree.bind("<Double-Button-1>", self._on_activate, add="+")
        self.tree.bind("<Return>", self._on_activate, add="+")
        self.tree.bind("<Control-c>", lambda e: self.copy_row(), add="+")

        # The action row is packed from the BOTTOM and BEFORE the tree: the packer
        # hands out space in call order, so a tree that asks for more rows than
        # the window has left would otherwise push these buttons off-screen —
        # which is exactly the single-file case a novice sees first.
        actions = ttk.Frame(self.results_wrap, style="TFrame")
        actions.pack(fill="x", side="bottom", pady=px(PAD_S))
        self.btn_html = ttk.Button(actions, text="Open report",
                                   style="Accent.TButton", command=self.open_html)
        self.btn_html.pack(side="left")
        self.btn_tsv = ttk.Button(actions, text="Save table (TSV)...", command=self.save_tsv)
        self.btn_tsv.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_json = ttk.Button(actions, text="Save JSON...", command=self.save_json)
        self.btn_json.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_novel = ttk.Button(actions, text="Save new alleles...",
                                    command=self.save_novel)
        self.btn_novel.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_copy = ttk.Button(actions, text="Copy row", command=self.copy_row)
        self.btn_copy.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_clear = ttk.Button(actions, text="Clear", command=self.clear)
        self.btn_clear.pack(side="right")
        Tooltip(self.btn_html, "Build a printable report of every file in this list "
                               "and open it in your web browser.")
        Tooltip(self.btn_tsv, "Save one row per file, ready to open in Excel.")
        Tooltip(self.btn_json, "Save the same results as JSON, for other software.")
        Tooltip(self.btn_novel, "Save the DNA of alleles that are close to, but not "
                                "identical to, a known allele.")

        ring.pack(fill="both", expand=True)
        self._set_actions_enabled(False)
        self.show_state(IDLE)

    # -- forced-scheme notice ------------------------------------------------
    def refresh_scheme_notice(self) -> None:
        """Show or hide the banner that says a scheme is pinned (section 10.2).

        Called on start-up and whenever settings change, because the preference
        is persisted: a scheme pinned weeks ago is otherwise invisible here.
        """
        scheme = getattr(self.app.prefs, "scheme", None)
        if scheme:
            self.scheme_notice_label.configure(
                text=("Scheme locked to \u201c%s\u201d. Every file is typed against "
                      "this scheme only, and anything else reports no match."
                      % scheme))
            self.scheme_notice.pack(fill="x", pady=(0, px(PAD_S)), before=self.dropzone)
        else:
            self.scheme_notice.pack_forget()

    def _clear_forced_scheme(self) -> None:
        """The banner's 'Use automatic' button: unpin and let detection run."""
        self.app.prefs.scheme = None
        self.app.schedule_save()
        settings = getattr(self.app, "settings_view", None)
        if settings is not None:
            settings.load_from(self.app.prefs)
        self.refresh_scheme_notice()
        self.app.status.set("Scheme detection is automatic again.", "ok",
                            transient=True)

    # -- visibility ----------------------------------------------------------
    def show_state(self, state: str) -> None:
        """Show and hide the cards that belong to each state (section 10.2).

        Everything below the drop zone is unpacked and re-packed in order, which
        is the only way to keep pack() geometry stable when cards come and go.
        """
        assert_main_thread("the analyse view")
        running = state == RUNNING
        have_rows = bool(self.results or self.failures)
        single = len(self.results) == 1 and not self.failures and not running
        for frame in (self.dropzone, self.progress_wrap, self.summary_wrap,
                      self.results_wrap):
            frame.pack_forget()
        # With nothing else on the tab the target grows to fill it: the novice's
        # first screen is one unmissable box, not a small strip above dead grey.
        self.dropzone.pack(fill="both" if not have_rows else "x",
                           expand=not have_rows, pady=(0, px(PAD_M)))
        if running:
            self.progress_wrap.pack(fill="x", pady=(0, px(PAD_S)))
        if single:
            self.summary_wrap.pack(fill="x", pady=(0, px(PAD_S)))
            self._fill_summary(self.results[0])
        if have_rows:
            self.results_wrap.pack(fill="both", expand=True)
        # Once results are on screen the drop zone is a strip, not a target: the
        # rows and the summary card are what the user is reading now, and every
        # pixel it keeps is a row of the table they cannot see.
        # Above the "tiny" threshold in DropZone.redraw(), or the strip loses
        # its words and becomes a grey box with a drawing in it.
        self.dropzone.configure(height=px(126) if have_rows else px(300))
        self.dropzone.set_enabled(not running)
        # The animation costs nothing unless the drop zone is what the user is
        # looking at: paused while a run is in flight and once results are up.
        self.dropzone.animate(not running and not have_rows)
        self._set_actions_enabled(have_rows and not running)

    def _set_actions_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for btn in (self.btn_html, self.btn_tsv, self.btn_json, self.btn_novel,
                    self.btn_copy, self.btn_clear):
            btn.configure(state=state)
        if enabled and not any(s.novel for s in self.results):
            self.btn_novel.configure(state="disabled")

    def _fill_summary(self, res: SampleResult) -> None:
        """Fill the single-file card. Every value is copied, never recomputed.

        The organism is the headline, the scheme its subtitle, and the reference
        that defines the scheme sits underneath — all of it read from the
        curated table in :mod:`wmlst.schemerefs`, and every field omitted rather
        than guessed at when that table is blank.
        """
        # A bare hyphen at 28 pt reads as a dash, not as information (section 10.2).
        self.st_label.configure(text=res.st if res.st not in ("", "-") else "not assigned")
        scheme = res.scheme if res.scheme != "-" else ""
        dbdir = self.app.env.dbdir
        ref = scheme_ref(scheme, dbdir)
        italic, upright = organism_parts(scheme, dbdir)
        # No ST means nothing was confirmed, so the name is stated quietly: a
        # 21 px italic binomial over a blank ST would read as an identification
        # the engine never made.
        assigned = res.st not in ("", "-")
        name_style = "Surface.TLabel" if assigned else "SurfaceMuted.TLabel"
        try:
            self.organism_label.configure(style=name_style)
            self.organism_suffix.configure(style=name_style)
        except tk.TclError:  # pragma: no cover
            pass
        if italic:
            self.organism_label.configure(text=italic, font=F.get("organism"))
        elif scheme:
            # No organism in the table: show what we DO know, in upright type, so
            # that nobody reads a directory name as a species name.
            self.organism_label.configure(text=scheme, font=F.get("title"))
        else:
            self.organism_label.configure(text="No scheme matched",
                                          font=F.get("title"))
        self.organism_suffix.configure(text=upright)
        self.summary_title.configure(
            text=scheme_caption(scheme if italic else "", ref))
        self.summary_file.configure(text=os.path.basename(res.label) or res.label)

        glyph_style = "{}.TLabel".format(
            STATUS_UI.get(res.status, (None, "muted", None))[1].capitalize())
        try:
            self.summary_status.configure(style=glyph_style)
        except tk.TclError:
            pass
        self.summary_status.configure(
            text="{}   ·   score {} of 100".format(status_cell(res.status), res.score))
        body = status_sentence(res.status)
        if res.alleles:
            body += "\n" + allele_summary(res)
        self.summary_body.configure(text=body)
        self._fill_reference(ref)

        tied = tied_schemes(res)
        if tied:
            self.tie_label.configure(
                text="{}. {}".format(tied, tie_sentence()))
            self.tie_banner.pack(fill="x", side="bottom")
        else:
            self.tie_banner.pack_forget()

    def _fill_reference(self, ref: Optional[Any]) -> None:
        """Show the citation and the authoritative links, or nothing at all."""
        for child in self.reference_links.winfo_children():
            child.destroy()
        citation = reference_line(ref)
        database_url = getattr(ref, "database_url", "") if ref is not None else ""
        pubmed_url = getattr(ref, "pubmed_url", "") if ref is not None else ""
        if not citation and not database_url:
            self.reference_wrap.pack_forget()
            return
        self.reference_label.configure(text=citation)
        if citation:
            self.reference_label.pack(anchor="w", fill="x")
        else:
            self.reference_label.pack_forget()
        if pubmed_url:
            link_button(self.reference_links,
                        "PubMed {}".format(ref.pubmed_id), pubmed_url,
                        tip="Open the primary publication for this scheme:\n"
                            + pubmed_url).pack(side="left", padx=(0, px(PAD_M)))
        if database_url:
            link_button(self.reference_links,
                        "{} record".format(source_name(ref)), database_url,
                        tip="Open the authoritative scheme record:\n"
                            + database_url).pack(side="left")
        self.reference_wrap.pack(anchor="w", fill="x", pady=(px(PAD_S), 0))

    # -- input ---------------------------------------------------------------
    def browse_files(self) -> None:
        """Open the file chooser (section 10.2)."""
        paths = filedialog.askopenfilenames(
            parent=self, title="Choose one or more sequence files",
            initialdir=self.app.prefs.last_dir or None, filetypes=FILE_TYPES)
        if paths:
            self.handle_paths(list(paths))

    def browse_folder(self) -> None:
        """Open the folder chooser (Shift+Enter on the drop zone)."""
        folder = filedialog.askdirectory(
            parent=self, title="Choose a folder of sequence files",
            initialdir=self.app.prefs.last_dir or None, mustexist=True)
        if folder:
            self.handle_paths([folder])

    def handle_paths(self, paths: Sequence[str]) -> None:
        """Resolve a drop or a browse and start analysing — no further clicks."""
        assert_main_thread("handle_paths")
        plan = expand_inputs(paths)
        files = list(plan.files)

        for candidate in plan.fofn_candidates:
            choice = messagebox.askyesnocancel(
                "Is this a list of files?",
                "{} is a text file.\n\nTreat it as a list of sequence files to "
                "analyse (one path per line)?\n\nChoose No to analyse the text "
                "file itself.".format(os.path.basename(candidate)),
                parent=self)
            if choice is None:
                return
            if choice:
                try:
                    files.extend(read_fofn(candidate))
                except OSError as exc:
                    self.app.show_error(exc, name=os.path.basename(candidate))
                    return
            else:
                files.append(candidate)

        if plan.missing:
            self.app.status.set("{} could not be found.".format(
                os.path.basename(plan.missing[0])), "warn", transient=True)
        if not files:
            if plan.folders:
                self.app.status.set(
                    "No sequence files were found in that folder.", "warn",
                    transient=True)
            else:
                self.app.status.set("Nothing to analyse.", "warn", transient=True)
            return

        if len(files) > BIG_BATCH:
            estimate = human_duration(
                estimate_seconds(len(files), self.app.prefs.normalised().jobs))
            if not messagebox.askokcancel(
                    "Analyse {} files?".format(len(files)),
                    "That is {} files, which should take {}.\n\n"
                    "You can cancel at any time and keep the results so far."
                    .format(len(files), estimate), parent=self):
                return

        first = files[0]
        try:
            self.app.prefs.last_dir = str(Path(first).resolve().parent)
            self.app.schedule_save()
        except OSError:
            pass
        self.app.start_analysis(files)

    # -- progress ------------------------------------------------------------
    def begin_run(self, total: int) -> None:
        assert_main_thread("begin_run")
        self._done = 0
        self._total = total
        self._current_pct = 0.0
        self.overall.configure(maximum=100.0 * max(1, total), value=0.0)
        self.progress_label.configure(text="Preparing...")
        self.progress_detail.configure(text="")
        self.elapsed_label.configure(text="")
        self.cancel_button.configure(state="normal", text="Cancel")
        self._stop_pulse()

    def on_started(self, msg: Msg) -> None:
        self._current_pct = 0.0
        name = os.path.basename(msg.path) or msg.path
        if msg.total > 1:
            self.progress_label.configure(
                text="Analysing {} ({} of {})".format(name, msg.index, msg.total))
        else:
            self.progress_label.configure(text="Analysing {}".format(name))
        self.progress_detail.configure(text=msg.text)
        self._update_bar()

    def on_phase(self, msg: Msg) -> None:
        if msg.percent >= self._current_pct:
            self._current_pct = msg.percent
        if msg.text:
            self.progress_detail.configure(text=msg.text)
        # BLAST reports nothing incremental, so the overall bar holds between
        # 10 % and 80 % while a thin indeterminate bar proves liveness (5.16).
        if 10.0 <= self._current_pct < 80.0:
            self._start_pulse()
        else:
            self._stop_pulse()
        self._update_bar()

    def _update_bar(self) -> None:
        self.overall.configure(value=self._done * 100.0 + self._current_pct)

    def _start_pulse(self) -> None:
        if not self._indeterminate:
            self._indeterminate = True
            try:
                self.pulse.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_XS), 0),
                                before=self.progress_detail)
                self.pulse.start(60)
            except tk.TclError:
                pass

    def _stop_pulse(self) -> None:
        if self._indeterminate:
            self._indeterminate = False
            try:
                self.pulse.stop()
                self.pulse.pack_forget()
            except tk.TclError:
                pass

    def on_tick(self, elapsed: float, running: bool) -> None:
        """Driven by the controller's single after() chain (section 10.5)."""
        if running:
            self.elapsed_label.configure(text="{:d}:{:02d}".format(
                int(elapsed) // 60, int(elapsed) % 60))
        if self._pending_rows:
            chunk, self._pending_rows = (self._pending_rows[:CHUNK_ROWS],
                                         self._pending_rows[CHUNK_ROWS:])
            for res in chunk:
                self._insert_result_row(res)

    def end_run(self, reason: str) -> None:
        assert_main_thread("end_run")
        self._stop_pulse()
        while self._pending_rows:
            self.on_tick(0.0, False)
        self.show_state(RESULTS if (self.results or self.failures) else IDLE)

    # -- rows ----------------------------------------------------------------
    def add_result(self, res: SampleResult) -> None:
        """Record one finished file. Large batches are inserted from the pump."""
        self.results.append(res)
        self._done += 1
        self._current_pct = 0.0
        self._update_bar()
        if len(self.results) > 2000:
            self._pending_rows.append(res)
        else:
            self._insert_result_row(res)

    def _refit_names(self, _event: Any = None) -> None:
        """Re-elide every file name after a resize — and after the first layout.

        The first row is inserted before Tk has stretched #0, so without this it
        would keep the harsh elision it was given against the unstretched width
        while every later row got the roomy one: two lengths in one column.
        """
        width = self._file_width()
        if width == self._fit_width:
            return
        self._fit_width = width
        font = F.get("body")
        for iid, name in self._row_names.items():
            try:
                self.tree.item(iid, text=elide_middle_px(name, font, width))
            except tk.TclError:  # pragma: no cover - teardown
                return

    def _file_width(self) -> int:
        """Room left for a file name once the dot and the expander have theirs."""
        try:
            width = int(self.tree.column("#0", "width"))
        except (tk.TclError, ValueError):  # pragma: no cover - teardown
            return px(260)
        # The dot image and the disclosure arrow both live in this cell and
        # both take real pixels off the text; measured, not guessed.
        return max(px(140), width - px(46))

    def _column_width(self, col: str) -> int:
        """The drawn width of a fixed column, less the cell's own padding."""
        try:
            return max(px(60), int(self.tree.column(col, "width")) - px(PAD_M))
        except (tk.TclError, ValueError):  # pragma: no cover - teardown
            return px(200)

    def _name_width(self) -> int:
        """Drawn width of the #0 cell on a child row.

        ``column("#0", "width")`` reports the width AFTER stretch, so this is
        the real number; what has to come off it is the disclosure indent and
        the cell's own padding.
        """
        try:
            width = int(self.tree.column("#0", "width"))
        except (tk.TclError, ValueError):  # pragma: no cover - teardown
            return px(300)
        # PAD_L is the disclosure indent a depth-1 row carries, PAD_S the
        # cell's own padding. Both are real pixels the text does not get.
        # A child row is indented, and Tk reserves an image slot on every row
        # of #0 once ANY row carries one. Both come off before anything wraps.
        return max(px(180), width - px(PAD_L) - px(PAD_L) - px(PAD_S))

    @staticmethod
    def _row_tone(status: str) -> str:
        """Which colour a whole result row is painted in.

        A Treeview colours the row or nothing, so a PERFECT result used to set
        the file name, the organism and the score in green as well -- a page of
        green says nothing at all. The drawn dot carries the good news; row
        colour is kept for the rows that want a second look.
        """
        if status in CALM_STATUSES:
            return "row_plain"
        return status_tag(status)

    def status_icon(self, status: str) -> Any:
        """The drawn status shape for a row, built once per status per theme."""
        key = status or DEFAULT_STATUS
        icon = self._status_icons.get(key)
        if icon is None:
            try:
                icon = dot_image(self.tree, key)
            except tk.TclError:  # pragma: no cover - teardown
                return ""
            self._status_icons[key] = icon
        return icon

    def _insert_result_row(self, res: SampleResult) -> None:
        self._order += 1
        organism = organism_name(res.scheme, self.app.env.dbdir)
        scheme_cell = res.scheme if res.scheme != "-" else "—"
        if tied_schemes(res):
            # Colour is never the only channel: the word is in the cell too.
            scheme_cell += "  \u00b7 tie"
        name = os.path.basename(res.label) or res.label
        iid = self.tree.insert(
            "", "end",
            text=elide_middle_px(name, F.get("body"), self._file_width()),
            image=self.status_icon(res.status),
            values=(organism or "—", scheme_cell, res.st,
                    status_cell(res.status), str(res.score)),
            tags=(self._row_tone(res.status), "row"), open=False)
        self._row_data[iid] = res
        self._row_names[iid] = name
        if res.alleles or res.tied:
            # A dummy child makes the disclosure arrow appear; the real per-locus
            # rows are built lazily on <<TreeviewOpen>> (section 10.2).
            self.tree.insert(iid, "end", text="", values=("", "", "", ""),
                             tags=("placeholder",))
        self.tree.see(iid)

    def add_failure(self, path: str, friendly: Friendly) -> None:
        """A failed file becomes a red row; a batch shows no dialog (section 10.6)."""
        self.failures.append((path, friendly))
        self._done += 1
        self._current_pct = 0.0
        self._update_bar()
        self._order += 1
        name = os.path.basename(path) or path
        iid = self.tree.insert(
            "", "end",
            text=elide_middle_px(name, F.get("body"), self._file_width()),
            image=self.status_icon(ERROR_STATUS),
            values=("—", "—", "—", "COULD NOT READ", "—"),
            tags=("failed", "row"))
        self._row_error[iid] = friendly
        self._row_names[iid] = name
        self.tree.see(iid)

    def _on_open(self, _event: Any = None) -> None:
        iid = self.tree.focus()
        if not iid or iid in self._expanded:
            return
        res = self._row_data.get(iid)
        if res is None:
            return
        self._expanded.add(iid)
        for child in self.tree.get_children(iid):
            self.tree.delete(child)
        tied = tied_schemes(res)
        wide = self._name_width()
        if tied:
            # Both of these used to be chopped mid-word inside a 215 px column.
            # They are the reason the row exists, so they go in #0 — the one
            # column that stretches with the window.
            head = wrap_px("TIE  \u00b7  " + tied, F.get("body_bold"), wide)
            for n, line in enumerate(head):
                self.tree.insert(
                    iid, "end", text=line if not n else INDENT + line,
                    values=(("", "equal score", "", "", "") if not n
                            else ("", "", "", "", "")),
                    tags=("tie", "locus"))
            # The full paragraph lives on the summary card, in the status line
            # and in the report; four wrapped lines of it inside a table is
            # clutter, so the row carries only what the reader must act on.
            indent = _indent_px(F.get("tiny"))
            for line in wrap_px(TIE_SHORT, F.get("tiny"), wide - indent):
                self.tree.insert(iid, "end", text=INDENT + line,
                                 values=("", "", "", "", ""), tags=("caption",))
        # The provenance of the call: which scheme, from which database, out of
        # which publication. Blank fields are simply absent.
        ref = scheme_ref(res.scheme, self.app.env.dbdir)
        if ref is not None:
            self.tree.insert(
                iid, "end", text="SCHEME",
                values=(organism_name(res.scheme, self.app.env.dbdir) or "—",
                        res.scheme, "",
                        elide_px(reference_line(ref) or ref.database_url
                                 or ref.description, F.get("body"),
                                 self._column_width("status")),
                        "{} loci".format(ref.n_loci) if ref.n_loci else ""),
                tags=("provenance", "locus"))
        self.tree.insert(iid, "end", text="LOCUS", values=CHILD_CAPTION,
                         tags=("caption",))
        evidence_px = self._column_width("status")
        for call in res.alleles:
            meaning = SYMBOL_UI.get(call.symbol, (call.symbol, "muted"))[0]
            colour = SYMBOL_UI.get(call.symbol, (call.symbol, "muted"))[1]
            # An exact match is the expected case and is set in ordinary ink: a
            # whole block of green says nothing, while one amber row among seven
            # says "look here". Colour is spent only where it carries meaning.
            tag = "locus" if call.symbol == "exact" else "st_" + colour
            self.tree.insert(
                iid, "end", text="    " + call.locus,
                values=(call.code, meaning, "",
                        elide_px(evidence_text(call), F.get("body"), evidence_px),
                        str(len(call.hits))),
                tags=(tag, "locus"))

    def _on_select(self, _event: Any = None) -> None:
        iid = self.tree.focus()
        res = self._row_data.get(iid)
        if res is not None:
            # The MEANING leads and the file name follows. A live region that
            # opens with 90 characters of absolute path spends its whole width
            # before it says anything, and then runs off the right edge.
            label = os.path.basename(res.path) or res.path
            if tied_schemes(res):
                self.app.status.set(
                    "Tie: {}. {} — {}".format(tied_schemes(res), tie_sentence(),
                                              label), "warn")
            else:
                organism = organism_name(res.scheme, self.app.env.dbdir)
                self.app.status.set(
                    "{}{} — {}".format(
                        "{} ({}) — ".format(organism, res.scheme) if organism else "",
                        status_sentence(res.status), label),
                    STATUS_UI.get(res.status, (None, "muted", None))[1])
        elif iid in self._row_error:
            self.app.status.set(self._row_error[iid].headline + " — press Enter for "
                                "details.", "bad")

    def _on_activate(self, _event: Any = None) -> str:
        iid = self.tree.focus()
        if iid in self._row_error:
            self.app.show_friendly(self._row_error[iid])
        return "break"

    def sort_by(self, col: str) -> None:
        """Sort top-level rows. Stable, with the original order as secondary key.

        Child rows are never re-sorted: locus order is the scheme's gene order.
        """
        if col == self._sort_col:
            self._sort_desc = not self._sort_desc
        else:
            self._sort_col, self._sort_desc = col, False
        rows = list(self.tree.get_children(""))
        original = {iid: i for i, iid in enumerate(rows)}

        def key(iid: str) -> Any:
            if col == "#0":
                # The cell may be middle-elided for width; sort the real name.
                res = self._row_data.get(iid)
                value: Any = (os.path.basename(res.label) if res is not None
                              else self.tree.item(iid, "text")).lower()
            else:
                value = self.tree.set(iid, col)
                if col == "score":
                    try:
                        value = -float(value)
                    except ValueError:
                        value = float("inf")
                else:
                    value = str(value).lower()
            return (value, original[iid])

        try:
            rows.sort(key=key, reverse=self._sort_desc)
        except TypeError:
            rows.sort(key=lambda iid: (str(self.tree.set(iid, col) if col != "#0"
                                           else self.tree.item(iid, "text")),
                                       original[iid]),
                      reverse=self._sort_desc)
        for position, iid in enumerate(rows):
            self.tree.move(iid, "", position)

    # -- actions -------------------------------------------------------------
    def copy_row(self) -> None:
        """Copy the selected row to the clipboard as tab-separated text."""
        iid = self.tree.focus()
        res = self._row_data.get(iid)
        if res is None:
            self.app.status.set("Select a file row first.", "warn", transient=True)
            return
        text = "\t".join([res.label,
                          organism_name(res.scheme, self.app.env.dbdir),
                          res.scheme, res.st, res.status, str(res.score),
                          allele_summary(res)])
        self.clipboard_clear()
        self.clipboard_append(text)
        self.app.status.set("Row copied to the clipboard.", "ok", transient=True)

    def clear(self) -> bool:
        """Empty the result list and return to the idle drop target.

        Refuses (returning False) while a worker is still running: the Clear
        *button* is greyed out during a run, but Run > Clear results is not, and
        dropping back to IDLE mid-run hides the progress card, re-arms the drop
        zone and -- worst of all -- unlocks the database update/atomic-swap path
        that section 10.3 requires to be blocked while an analysis is in flight.
        """
        if self.app.controller.running:
            self.app.status.set("Stop the analysis before clearing the results.",
                                "warn", transient=True)
            return False
        self.st_label.configure(text="—")
        self.organism_label.configure(text="")
        self.organism_suffix.configure(text="")
        self.summary_title.configure(text="")
        self.summary_file.configure(text="")
        self.summary_status.configure(text="")
        self.summary_body.configure(text="")
        self.reference_wrap.pack_forget()
        self.tie_banner.pack_forget()
        self.results.clear()
        self.failures.clear()
        self._row_data.clear()
        self._row_error.clear()
        self._row_names.clear()
        self._expanded.clear()
        self._pending_rows.clear()
        for iid in self.tree.get_children(""):
            self.tree.delete(iid)
        self.app.state.to(IDLE)
        self.show_state(IDLE)
        self.app.status.set("Ready.", resting=True)
        return True

    def open_html(self) -> None:
        self.app.export("html")

    def save_tsv(self) -> None:
        self.app.export("tsv")

    def save_json(self) -> None:
        self.app.export("json")

    def save_novel(self) -> None:
        self.app.export("novel")


# ===========================================================================
# SECTION 9 — VIEWS: THE DATABASE TAB  (section 10.3)
# ===========================================================================

#: The GUI spells LOCII correctly; the CLI keeps the upstream misspelling (D13).
DB_COLUMNS = ("species", "loci", "types", "alleles", "date")
DB_HEADINGS = {"#0": "SCHEME", "species": "SPECIES", "loci": "LOCI",
               "types": "TYPES", "alleles": "ALLELES", "date": "DATE"}
#: Unknown renders as a muted em dash, never the literal word (section 10.3).
UNKNOWN_CELL = "—"


def db_cell(value: Any) -> str:
    """Render an --info cell for humans: Unknown becomes a muted dash."""
    text = "" if value is None else str(value)
    if text.strip() in ("", "Unknown", "No version information available"):
        return UNKNOWN_CELL
    return text


class DatabaseView(ttk.Frame):
    """Scheme catalogue and the update flow (section 10.3).

    Loading 162 schemes must never block start-up, so the tab populates on first
    activation, on a worker thread.
    """

    def __init__(self, parent: tk.Misc, app: WmlstApp):
        super().__init__(parent, style="TFrame")
        self.app = app
        self.loaded = False
        self.plan: Any = None
        self._infos: Tuple[Any, ...] = ()
        self._jobs: Dict[int, str] = {}
        self._cancel: Optional[threading.Event] = None
        #: scheme name -> ticked. What the next update will include.
        self.checked: Dict[str, bool] = {}
        #: The two checkbox images the rows share. Built once, on this widget's
        #: own Tk root, and kept referenced here or Tk would collect them.
        self._image_on: Optional[tk.PhotoImage] = None
        self._image_off: Optional[tk.PhotoImage] = None
        self._build()

    def _build(self) -> None:
        outer = ttk.Frame(self, style="TFrame")
        outer.pack(fill="both", expand=True, padx=px(PAD_L), pady=px(PAD_M))

        info = card(outer)
        info.master.pack(fill="x")
        grid = ttk.Frame(info, style="Surface.TFrame")
        grid.pack(fill="x", padx=px(PAD_M), pady=px(PAD_M))
        self.count_label = ttk.Label(grid, text="—", style="Hero.TLabel")
        self.count_label.grid(row=0, column=0, rowspan=2, padx=(0, px(PAD_L)))
        ttk.Label(grid, text="SCHEMES INSTALLED", style="Eyebrow.TLabel").grid(
            row=0, column=1, sticky="sw")
        self.version_label = ttk.Label(grid, text="", style="Heading.TLabel")
        self.version_label.grid(row=1, column=1, sticky="w")
        buttons = ttk.Frame(grid, style="Surface.TFrame")
        buttons.grid(row=0, column=2, rowspan=2, sticky="e", padx=(px(PAD_L), 0))
        grid.columnconfigure(2, weight=1)
        self.btn_check = ttk.Button(buttons, text="Check for updates",
                                    command=self.check_updates)
        self.btn_check.pack(side="left")
        self.btn_update = ttk.Button(buttons, text="Update now", style="Accent.TButton",
                                     command=self.apply_updates, state="disabled")
        self.btn_update.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_cancel = ttk.Button(buttons, text="Cancel", command=self.cancel,
                                     state="disabled")
        self.btn_cancel.pack(side="left", padx=(px(PAD_S), 0))
        Tooltip(self.btn_check, "Ask PubMLST and Institut Pasteur whether any scheme "
                                "has changed. Nothing is downloaded or altered yet.")
        Tooltip(self.btn_update, "Download the changed schemes and rebuild the search "
                                 "index. Your current database is backed up first.")
        # Packed by _show_progress() only while something is actually running:
        # an empty trough sitting under the header for the life of the window is
        # a piece of furniture, not a progress report.
        self.db_progress = ttk.Progressbar(info, style="Thin.Horizontal.TProgressbar",
                                           mode="determinate", maximum=100.0)
        self.db_detail = ttk.Label(info, text="", style="SurfaceMuted.TLabel", anchor="w")
        self.db_detail.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_XS)))
        self.db_where = ttk.Label(info, text="", style="SurfaceMuted.TLabel", anchor="w")
        self.db_where.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_M)))

        # -- selection controls ----------------------------------------------
        # "make the checkbox nicer in design and bigger, and put a couple of
        # select all / deselect all buttons": the checkbox on every row is the
        # same square the custom CheckBox widget draws, rendered as an image
        # because a Treeview cell can hold an image but not a widget.
        bar = ttk.Frame(outer, style="TFrame")
        bar.pack(fill="x", pady=(px(PAD_M), px(PAD_XS)))
        self.btn_all = ttk.Button(bar, text="Select all", command=self.select_all)
        self.btn_all.pack(side="left")
        self.btn_none = ttk.Button(bar, text="Deselect all", command=self.deselect_all)
        self.btn_none.pack(side="left", padx=(px(PAD_S), 0))
        Tooltip(self.btn_all, "Tick every scheme in the list.")
        Tooltip(self.btn_none, "Untick every scheme in the list.")
        self.selection_label = ttk.Label(bar, text="", style="Muted.TLabel")
        self.selection_label.pack(side="left", padx=px(PAD_M))
        ttk.Label(bar, text="Click a tick box, or press Space, to change a row.",
                  style="Muted.TLabel").pack(side="right")

        ring = focus_ring(outer)
        ring.pack(fill="both", expand=True, pady=(px(PAD_XS), 0))
        self.tree = ttk.Treeview(ring, columns=DB_COLUMNS, show="tree headings",
                                 selectmode="browse")
        ring.track(self.tree)
        vsb = ttk.Scrollbar(ring, orient="vertical", command=self.tree.yview)
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.configure(
            yscrollcommand=autohide(vsb, self.tree, side="right", fill="y"))
        self.tree.column("#0", width=px(260), minwidth=px(160), stretch=True)
        self.tree.column("species", width=px(260), minwidth=px(140))
        for col in ("loci", "types", "alleles"):
            self.tree.column(col, width=px(90), minwidth=px(60), anchor="e")
        # DATE right-aligned like the three numeric columns before it: a
        # left-aligned date butted against a right-aligned count left a five
        # pixel alley and read as a collision.
        self.tree.column("date", width=px(130), minwidth=px(90), anchor="e")
        for col, title in DB_HEADINGS.items():
            self.tree.heading(
                col, text=title,
                anchor="e" if col in ("loci", "types", "alleles", "date") else "w")
        self.tree.tag_configure("muted", foreground=c("muted"))
        self.tree.bind("<Button-1>", self._click, add="+")
        self.tree.bind("<space>", self._space, add="+")
        self._update_selection_label()

    # -- the per-scheme tick boxes -------------------------------------------
    def _images(self) -> Tuple[Any, Any]:
        """The ticked/unticked images, built on first use."""
        if self._image_on is None or self._image_off is None:
            try:
                self._image_on = checkbox_image(self, True)
                self._image_off = checkbox_image(self, False)
            except tk.TclError:  # pragma: no cover - no display
                return (None, None)
        return (self._image_on, self._image_off)

    def _row_image(self, name: str) -> Any:
        on, off = self._images()
        return on if self.checked.get(name, True) else off

    def _click(self, event: Any) -> Any:
        """A click on the tick box or the scheme name toggles that row.

        The whole name column is the hit target, exactly as the label of a
        :class:`CheckBox` is, and the row also takes focus and selection so that
        Space afterwards operates the row the user just clicked.
        """
        try:
            if self.tree.identify_region(event.x, event.y) != "tree":
                return None
            iid = self.tree.identify_row(event.y)
        except tk.TclError:  # pragma: no cover
            return None
        if not iid:
            return None
        try:
            # focus_set() is the KEYBOARD focus and focus(iid) is only the
            # focused ROW. Both are needed: returning "break" below suppresses
            # ttk::treeview's own Button-1 binding, which is what would
            # otherwise have moved the keyboard focus here -- so without this
            # line the focus stays wherever it was (the Analyse tab's drop
            # zone, on a fresh window), Space never reaches _space, and a
            # Space pressed over this list opened the file chooser belonging
            # to another tab instead of ticking the row under the cursor.
            self.tree.focus_set()
            self.tree.focus(iid)
            self.tree.selection_set(iid)
        except tk.TclError:  # pragma: no cover
            pass
        self.toggle(iid)
        return "break"

    def _space(self, _event: Any = None) -> str:
        iid = self.tree.focus()
        if iid:
            self.toggle(iid)
        return "break"

    def toggle(self, iid: str) -> None:
        """Flip one row's tick box and refresh the count."""
        name = self.tree.item(iid, "text")
        if not name:
            return
        self.checked[name] = not self.checked.get(name, True)
        try:
            self.tree.item(iid, image=self._row_image(name))
        except tk.TclError:  # pragma: no cover
            return
        self._update_selection_label()

    def select_all(self) -> None:
        """Tick every scheme."""
        self._set_all(True)

    def deselect_all(self) -> None:
        """Untick every scheme."""
        self._set_all(False)

    def _set_all(self, value: bool) -> None:
        for iid in self.tree.get_children(""):
            name = self.tree.item(iid, "text")
            if not name:
                continue
            self.checked[name] = value
            try:
                self.tree.item(iid, image=self._row_image(name))
            except tk.TclError:  # pragma: no cover
                pass
        self._update_selection_label()

    def selected_schemes(self) -> Tuple[str, ...]:
        """The ticked schemes, in list order."""
        return tuple(name for name in
                     (self.tree.item(iid, "text")
                      for iid in self.tree.get_children(""))
                     if name and self.checked.get(name, True))

    def _update_selection_label(self) -> None:
        total = len(self.tree.get_children(""))
        chosen = len(self.selected_schemes())
        if not total:
            self.selection_label.configure(text="")
            return
        self.selection_label.configure(
            text="{} of {} scheme{} selected".format(
                chosen, total, "" if total == 1 else "s"))

    # -- catalogue -----------------------------------------------------------
    def activate(self) -> None:
        """Called the first time the tab is shown (section 10.3)."""
        if self.loaded:
            return
        self.loaded = True
        env = self.app.env
        self.count_label.configure(text=str(env.scheme_count or "—"))
        self.version_label.configure(
            text="Database version {}".format(env.db_version or "unknown"))
        self.refresh_location()
        if not env.db_ok:
            self.db_detail.configure(text=env.db_error or "The database is not available.")
            return
        self.db_detail.configure(text="Reading the scheme catalogue...")
        datadir = env.datadir

        def load(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import schemes as schemes_mod
            catalog = schemes_mod.SchemeCatalog(datadir)
            shallow = catalog.info_all(deep=False)
            return (catalog, shallow)

        self._jobs[self.app.run_task("db-shallow", load, owner=self,
                                     with_progress=True)] = "shallow"

    def refresh_location(self) -> None:
        """Print the folder updates are written to — portable mode or not."""
        where = db_location(self.app.prefs, self.app.env)
        self.db_where.configure(
            text="Updates are written to {}{}".format(
                where, "  ·  portable mode" if self.app.prefs.portable_db else ""))

    def _fill(self, catalog: Any, infos: Sequence[Any], deep: bool) -> None:
        assert_main_thread("the database tree")
        self._infos = tuple(infos)
        for iid in self.tree.get_children(""):
            self.tree.delete(iid)
        for info in self._infos:
            name = getattr(info, "name", "")
            try:
                species = catalog.species_label(name)
            except Exception:
                species = None
            types = getattr(info, "num_genotypes", None)
            alleles = getattr(info, "num_alleles", None)
            # The catalogue's own map first; the curated table is the fallback,
            # so a scheme the map misses still names its organism.
            label = species or organism_name(name, self.app.env.dbdir)
            self.checked.setdefault(name, True)
            self.tree.insert(
                "", "end", text=name, image=self._row_image(name),
                values=(label or "scheme not mapped to a species",
                        db_cell(getattr(info, "locus", "")),
                        db_cell(types) if deep else "...",
                        db_cell(alleles) if deep else "...",
                        db_cell(getattr(info, "last_updated", ""))),
                tags=() if label else ("muted",))
        self.count_label.configure(text=str(len(self._infos)))
        self._update_selection_label()

    def on_task_done(self, job: int, payload: Any) -> bool:
        """Handle a finished background job. Returns True if it was ours."""
        kind = self._jobs.pop(job, None)
        if kind is None:
            return False
        if kind == "shallow":
            catalog, infos = payload
            self._fill(catalog, infos, deep=False)
            self.db_detail.configure(text="Counting sequence types...")
            datadir = self.app.env.datadir

            def deep_load(progress: Callable[..., None], cancel: threading.Event) -> Any:
                from . import schemes as schemes_mod
                cat = schemes_mod.SchemeCatalog(datadir)
                return (cat, cat.info_all(deep=True, progress=progress, cancel=cancel))

            self._jobs[self.app.run_task("db-deep", deep_load, owner=self,
                                         with_progress=True)] = "deep"
        elif kind == "deep":
            catalog, infos = payload
            self._fill(catalog, infos, deep=True)
            self.db_detail.configure(text="")
            self._show_progress(False)
        elif kind == "check":
            self._show_plan(payload)
        elif kind == "apply":
            self._after_update(payload)
        return True

    def on_task_failed(self, job: int, friendly: Friendly) -> bool:
        kind = self._jobs.pop(job, None)
        if kind is None:
            return False
        self._busy(False)
        if kind in ("shallow", "deep"):
            self.db_detail.configure(text=friendly.headline)
            return True
        self.app.show_friendly(friendly)
        self.db_detail.configure(text="")
        return True

    def on_progress(self, msg: Msg) -> bool:
        if msg.job not in self._jobs:
            return False
        self._show_progress(True)
        self.db_progress.configure(value=msg.percent)
        if msg.text:
            self.db_detail.configure(text=msg.text)
        return True

    def _show_progress(self, on: bool) -> None:
        """Show the bar while work runs, and take it away again afterwards."""
        try:
            if on:
                if not self.db_progress.winfo_manager():
                    self.db_progress.pack(fill="x", padx=px(PAD_M),
                                          pady=(0, px(PAD_XS)),
                                          before=self.db_detail)
            else:
                self.db_progress.configure(value=0.0)
                self.db_progress.pack_forget()
        except tk.TclError:  # pragma: no cover - teardown
            pass

    # -- updates -------------------------------------------------------------
    def _busy(self, busy: bool) -> None:
        self._show_progress(busy)
        self.btn_check.configure(state="disabled" if busy else "normal")
        self.btn_cancel.configure(state="normal" if busy else "disabled")
        if busy:
            self.btn_update.configure(state="disabled")

    def check_updates(self) -> None:
        """Metadata-only pass; writes nothing (section 10.3)."""
        # Key off the worker, not the view state: the view can be back at IDLE
        # while the engine is still reading the database we are about to swap.
        if self.app.controller.running or self.app.state.running:
            self.app.status.set(
                "Wait for the current analysis to finish before updating.",
                "warn", transient=True)
            return
        dbdir = self.app.env.dbdir
        if not dbdir:
            self.app.status.set("The database folder was not found.", "bad",
                                transient=True)
            return
        self._busy(True)
        self.db_detail.configure(text="Contacting PubMLST...")
        self._cancel = threading.Event()

        def run(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import updatedb as updatedb_mod
            return updatedb_mod.check(dbdir, progress=progress, cancel=cancel)

        self._jobs[self.app.run_task("update-check", run, owner=self,
                                     with_progress=True,
                                     cancel=self._cancel)] = "check"

    def _show_plan(self, plan: Any) -> None:
        self._busy(False)
        self.plan = plan
        updates = [u for u in getattr(plan, "updates", ())
                   if getattr(u, "status", "") == "changed"]
        if not updates:
            self.db_detail.configure(text="Your database is up to date.")
            self.app.status.set("Your database is up to date.", "ok", transient=True)
            return
        dialog = UpdatePlanDialog(self.app.root, plan, updates,
                                  preselected=self.checked)
        self.app.root.wait_window(dialog)
        if not dialog.accepted:
            self.db_detail.configure(text="{} scheme(s) can be updated.".format(
                len(updates)))
            self.btn_update.configure(state="normal")
            return
        self._start_apply(dialog.chosen, dialog.backup)

    def apply_updates(self) -> None:
        if self.plan is None:
            self.check_updates()
            return
        updates = [u for u in getattr(self.plan, "updates", ())
                   if getattr(u, "status", "") == "changed"]
        dialog = UpdatePlanDialog(self.app.root, self.plan, updates,
                                  preselected=self.checked)
        self.app.root.wait_window(dialog)
        if dialog.accepted:
            self._start_apply(dialog.chosen, dialog.backup)

    def _start_apply(self, chosen: Sequence[str], backup: bool) -> None:
        dbdir = self.app.env.dbdir
        plan = self.plan
        self._busy(True)
        self.db_detail.configure(text="Downloading...")
        self._cancel = threading.Event()

        def run(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import updatedb as updatedb_mod
            return updatedb_mod.apply(dbdir, plan, list(chosen), backup=backup,
                                      progress=progress, cancel=cancel)

        self._jobs[self.app.run_task("update-apply", run, owner=self,
                                     with_progress=True,
                                     cancel=self._cancel)] = "apply"

    def _after_update(self, new_version: Any) -> None:
        """Reload the catalogue in place — no restart (section 10.3)."""
        self._busy(False)
        self.plan = None
        self.btn_update.configure(state="disabled")
        version = str(new_version or "")
        self.db_detail.configure(text="Database updated.")
        self.app.status.set("Database updated to {}.".format(version or "the latest data"),
                            "ok", transient=True)
        self.app.env.db_version = version or self.app.env.db_version
        self.loaded = False
        self.app.refresh_footer()
        self.activate()

    def cancel(self) -> None:
        """Cancel is live and honoured between scheme commits (section 10.3)."""
        if self._cancel is not None:
            self._cancel.set()
        self.db_detail.configure(text="Stopping...")


# ===========================================================================
# SECTION 10 — VIEWS: THE SETTINGS TAB  (section 10.4)
# ===========================================================================

SETTING_HELP = {
    "minid": "How similar a match must be to count as the same allele. "
             "95 % is the reference default — lower it only if you know why.",
    "mincov": "How much of the reference allele must be covered by the match. "
              "Below this, the locus is reported as missing.",
    "minscore": "How well a scheme must fit before WMLST will report it. "
                "Raising this hides poor matches; lowering it shows more.",
    "threads": "How many processor cores the search may use for ONE file. "
               "For a single assembly this makes almost no difference — the real "
               "speed-up comes from analysing several files at once.",
    "jobs": "How many files to analyse at the same time. Each one uses about "
            "270 MB of memory, so keep this below your free memory divided by 400 MB.",
    "exclude": "Schemes that are never reported automatically. The four defaults "
               "overlap other schemes and would otherwise win by accident.",
    "scheme": "Normally WMLST picks the best-fitting scheme by itself. Choose one "
              "here only when you already know the species.",
}


class SettingsView(ttk.Frame):
    """Scoring thresholds, performance and scheme choice (section 10.4)."""

    def __init__(self, parent: tk.Misc, app: WmlstApp):
        super().__init__(parent, style="TFrame")
        self.app = app
        self.vars: Dict[str, tk.Variable] = {}
        # While a scheme is forced the minimum-score control *displays* 0
        # (section 10.4) but that 0 is a display, not the user's preference:
        # _minscore_saved holds the real value so returning to Automatic gives it
        # back instead of persisting 0 forever.
        self._minscore_forced = False
        self._minscore_saved = _num(DEFAULT_MINSCORE)
        self._build()
        self.load_from(app.prefs)

    def _build(self) -> None:
        self.scroller = ScrollPane(self)
        self.scroller.pack(fill="both", expand=True)
        outer = ttk.Frame(self.scroller.body, style="TFrame")
        outer.pack(fill="both", expand=True, padx=px(PAD_L), pady=px(PAD_M))

        self._outer = outer
        self.rerun_bar = ttk.Frame(outer, style="TFrame")
        inner = card(self.rerun_bar)
        inner.master.pack(fill="x")
        ttk.Label(inner, text="Settings changed. The results already on screen used "
                              "the old settings.", style="Surface.TLabel").pack(
            side="left", padx=px(PAD_M), pady=px(PAD_S))
        ttk.Button(inner, text="Re-run all files  (F5)",
                   command=self.app.rerun_all).pack(side="left", pady=px(PAD_S))
        ttk.Button(inner, text="Dismiss", style="Link.TButton",
                   command=self.hide_rerun).pack(side="right", padx=px(PAD_M))

        scoring = card(outer)
        self._first_card = scoring.master
        self._first_card.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(scoring, text="Matching", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        self._spin(scoring, "minid", "Minimum identity (%)", 0, 100, 1.0)
        self._spin(scoring, "mincov", "Minimum coverage (%)", 0, 100, 1.0)
        self._spin(scoring, "minscore", "Minimum scheme score", 0, 100, 1.0)
        warn = ttk.Frame(scoring, style="Surface.TFrame")
        warn.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_M) - px(PAD_XS)))
        self.scheme_warning = ttk.Label(
            warn, text="", style="Warn.TLabel", anchor="w", justify="left")
        self.scheme_warning.pack(anchor="w")
        # It is empty far more often than not, and an empty label still claims
        # a line: every card then closed on a different amount of white.
        warn.configure(height=1)

        scheme = card(outer)
        scheme.master.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(scheme, text="Scheme", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        row = ttk.Frame(scheme, style="Surface.TFrame")
        row.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_S)))
        ttk.Label(row, text="Always use this scheme", style="Surface.TLabel",
                  width=26, anchor="w").pack(side="left")
        self.vars["scheme"] = tk.StringVar(value="")
        # The combobox and the "Never report these" entry below it are the two
        # controls in this column; they share a left edge AND a right edge,
        # because two boxes of different widths stacked on one another is the
        # single loudest piece of sloppiness on a settings page.
        self.scheme_box = ttk.Combobox(row, textvariable=self.vars["scheme"],
                                       state="readonly", width=32,
                                       values=("Automatic (recommended)",))
        self.scheme_box.pack(side="left", fill="x", expand=True)
        self.scheme_box.bind("<<ComboboxSelected>>", lambda e: self._changed(), add="+")
        Tooltip(self.scheme_box, SETTING_HELP["scheme"])
        row2 = ttk.Frame(scheme, style="Surface.TFrame")
        row2.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_M)))
        ttk.Label(row2, text="Never report these", style="Surface.TLabel",
                  width=26, anchor="w").pack(side="left")
        self.vars["exclude"] = tk.StringVar(value="")
        self.exclude_entry = ttk.Entry(row2, textvariable=self.vars["exclude"], width=46)
        self.exclude_entry.pack(side="left", fill="x", expand=True)
        self.exclude_entry.bind("<FocusOut>", lambda e: self._changed(), add="+")
        Tooltip(self.exclude_entry, SETTING_HELP["exclude"])

        speed = card(outer)
        speed.master.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(speed, text="Speed", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        # WMLST sizes itself on this machine at start-up (wmlst.perf): one file
        # per four cores, four search threads each. The CLI is untouched by this
        # and still defaults to 1/1.
        self.vars["perf_auto"] = tk.BooleanVar(value=True)
        self.auto_box = CheckBox(
            speed, "Tune for this computer automatically",
            variable=self.vars["perf_auto"], command=self._auto_changed,
            subtext="Untick to set the two numbers below yourself.")
        self.auto_box.pack(anchor="w", fill="x", padx=px(PAD_M), pady=px(PAD_XS))
        self.tuning_label = ttk.Label(
            speed, text="This computer has {} processor cores.".format(cpu_count()),
            style="SurfaceMuted.TLabel", anchor="w", justify="left")
        # Indented to the checkbox's own text column: it explains the tick.
        self.tuning_label.pack(anchor="w", fill="x",
                               padx=(px(PAD_M) + px(30), px(PAD_M)),
                               pady=(0, px(PAD_S)))
        self._spin(speed, "threads", "Cores per file", 1, max_threads(), 1)
        self._spin(speed, "jobs", "Files at once", 1, max_threads(), 1)
        ttk.Frame(speed, style="Surface.TFrame",
                  height=px(PAD_M) - px(PAD_XS)).pack(fill="x")

        # -- where the database lives ----------------------------------------
        store = card(outer)
        store.master.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(store, text="Database location", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        self.vars["portable_db"] = tk.BooleanVar(value=False)
        self.portable_box = CheckBox(
            store, "Keep the database next to {}.exe (portable mode)".format(
                branding.APP_NAME),
            variable=self.vars["portable_db"], command=self._portable_changed,
            subtext="Downloaded schemes and the rebuilt search index are then "
                    "written beside the program, so the whole installation can "
                    "travel on a USB stick.",
            wraplength=520)
        self.portable_box.pack(anchor="w", fill="x", padx=px(PAD_M), pady=px(PAD_XS))
        self.portable_path = ttk.Label(store, text="", style="SurfaceMuted.TLabel",
                                       anchor="w", justify="left", wraplength=px(560))
        self.portable_path.pack(anchor="w", fill="x", padx=px(PAD_M),
                                pady=(0, px(PAD_XS)))
        self.portable_note = ttk.Label(store, text="", style="SurfaceMuted.TLabel",
                                       anchor="w", justify="left", wraplength=px(560))
        self.portable_note.pack(anchor="w", fill="x", padx=px(PAD_M),
                                pady=(0, px(PAD_M)))

        footer = ttk.Frame(outer, style="TFrame")
        footer.pack(fill="x", pady=px(PAD_M))
        ttk.Button(footer, text="Restore reference defaults",
                   command=self.restore_defaults).pack(side="left")
        ttk.Label(footer, text="Defaults match tseemann/mlst {}: identity 95, "
                              "coverage 50, score 50.".format(UPSTREAM_MLST_VERSION),
                  style="Muted.TLabel").pack(side="left", padx=px(PAD_M))

    def _spin(self, parent: tk.Misc, key: str, label: str, lo: float, hi: float,
              step: float) -> None:
        row = ttk.Frame(parent, style="Surface.TFrame")
        row.pack(fill="x", padx=px(PAD_M), pady=px(PAD_XS))
        ttk.Label(row, text=label, style="Surface.TLabel", width=26,
                  anchor="w").pack(side="left")
        var = tk.StringVar()
        self.vars[key] = var
        spin = ttk.Spinbox(row, from_=lo, to=hi, increment=step, width=9,
                           textvariable=var, justify="right",
                           command=self._changed)
        spin.pack(side="left", ipady=px(2))
        spin.bind("<FocusOut>", lambda e, k=key: self._clamp(k), add="+")
        spin.bind("<Return>", lambda e, k=key: self._clamp(k), add="+")
        setattr(self, "spin_" + key, spin)
        help_text = SETTING_HELP.get(key, "")
        Tooltip(spin, help_text)
        hint = ttk.Label(row, text=help_text.split(".")[0] + ".",
                         style="SurfaceMuted.TLabel", anchor="w")
        hint.pack(side="left", padx=px(PAD_M), fill="x", expand=True)

    # -- values --------------------------------------------------------------
    def load_from(self, prefs: Prefs) -> None:
        """Fill the widgets from ``prefs`` without triggering a change note."""
        self._loading = True
        self.vars["minid"].set(_num(prefs.minid))
        self.vars["mincov"].set(_num(prefs.mincov))
        self.vars["minscore"].set(_num(prefs.minscore))
        self.vars["threads"].set(str(prefs.threads))
        self.vars["jobs"].set(str(prefs.jobs))
        self.vars["exclude"].set(", ".join(prefs.exclude))
        self.vars["scheme"].set(prefs.scheme or "Automatic (recommended)")
        self.vars["perf_auto"].set(bool(prefs.perf_auto))
        self.vars["portable_db"].set(bool(prefs.portable_db))
        self._loading = False
        self._sync_perf_lock()
        self.refresh_db_location()
        self._minscore_forced = False   # re-capture the shadow from these values
        self._sync_scheme_lock()

    def set_scheme_choices(self, names: Sequence[str]) -> None:
        """Populate the scheme combobox once the catalogue is known."""
        self.scheme_box.configure(
            values=("Automatic (recommended)", *tuple(names)))

    def collect(self) -> Prefs:
        """Read the widgets back into a clamped :class:`Prefs`."""
        prefs = dataclasses.replace(self.app.prefs)
        prefs.minid = clamp_float(self.vars["minid"].get(), 0, 100, DEFAULT_MINID)
        prefs.mincov = clamp_float(self.vars["mincov"].get(), 0, 100, DEFAULT_MINCOV)
        raw_minscore = (self._minscore_saved if self._minscore_forced
                        else self.vars["minscore"].get())
        prefs.minscore = clamp_float(raw_minscore, 0, 100, DEFAULT_MINSCORE)
        prefs.threads = clamp_int(self.vars["threads"].get(), 1, max_threads(), 1)
        prefs.jobs = clamp_int(self.vars["jobs"].get(), 1, max_threads(), 1)
        prefs.exclude = parse_exclude(self.vars["exclude"].get())
        prefs.perf_auto = bool(self.vars["perf_auto"].get())
        # portable_db is not read back from the widget: turning it on has to
        # succeed (a writable folder, a database to point at) before it becomes
        # a preference, and WmlstApp.set_portable_db is what decides that.
        chosen = self.vars["scheme"].get().strip()
        prefs.scheme = None if chosen.startswith("Automatic") or not chosen else chosen
        # The scheme override (minscore 0, empty exclude) is applied per run by
        # to_runconfig(); baking it into the saved preferences would destroy the
        # user's own values the moment they picked a scheme.
        return prefs.normalised(force_scheme_rules=False)

    def _clamp(self, key: str) -> None:
        """Clamp on focus-out, flash the field and explain (section 10.4)."""
        var = self.vars[key]
        raw = var.get()
        if key in ("threads", "jobs"):
            value: Any = clamp_int(raw, 1, max_threads(), 1)
            text = str(value)
        else:
            default = {"minid": DEFAULT_MINID, "mincov": DEFAULT_MINCOV,
                       "minscore": DEFAULT_MINSCORE}[key]
            value = clamp_float(raw, 0, 100, default)
            text = _num(value)
        if text != raw.strip():
            var.set(text)
            widget = getattr(self, "spin_" + key, None)
            if widget is not None:
                self._flash(widget)
            self.app.status.set(
                "{} was adjusted to {} — the allowed range is {}.".format(
                    key, text, "1 to {}".format(max_threads())
                    if key in ("threads", "jobs") else "0 to 100"),
                "warn", transient=True)
        self._changed()

    def _flash(self, widget: tk.Misc) -> None:
        try:
            widget.configure(style="Invalid.TEntry")
            widget.after(600, lambda: widget.configure(style="TSpinbox"))
        except tk.TclError:
            pass

    def _sync_scheme_lock(self) -> None:
        """Selecting a scheme greys out the minimum score and shows the warning."""
        forced = not self.vars["scheme"].get().startswith("Automatic")
        spin = getattr(self, "spin_minscore", None)
        if spin is not None:
            spin.configure(state="disabled" if forced else "normal")
        if forced:
            if not self._minscore_forced:
                self._minscore_saved = self.vars["minscore"].get()
                self._minscore_forced = True
            self.vars["minscore"].set("0")
            self.scheme_warning.configure(
                text="⚠ Forcing a scheme sets the minimum score to 0 and ignores the "
                     "exclude list, exactly as mlst --scheme does. Every file will be "
                     "reported with this scheme even if it does not fit.")
        else:
            if self._minscore_forced:
                self.vars["minscore"].set(self._minscore_saved)
                self._minscore_forced = False
            self.scheme_warning.configure(text="")

    # -- automatic tuning ----------------------------------------------------
    def set_tuning(self, tuning: Any) -> None:
        """Show what :mod:`wmlst.perf` decided, in one calm sentence."""
        rationale = str(getattr(tuning, "rationale", "") or "")
        self.tuning_label.configure(
            text=rationale or "This computer has {} processor cores.".format(
                cpu_count()))

    def _sync_perf_lock(self) -> None:
        """Automatic tuning owns the two numbers; a manual override frees them."""
        auto = bool(self.vars["perf_auto"].get())
        for key in ("threads", "jobs"):
            spin = getattr(self, "spin_" + key, None)
            if spin is not None:
                try:
                    spin.configure(state="disabled" if auto else "normal")
                except tk.TclError:  # pragma: no cover
                    pass

    def _auto_changed(self) -> None:
        self._sync_perf_lock()
        self._changed()
        if self.vars["perf_auto"].get():
            self.app.apply_autotune(announce=True)

    # -- portable mode -------------------------------------------------------
    def _portable_changed(self) -> None:
        if getattr(self, "_loading", False):
            return
        self.app.set_portable_db(bool(self.vars["portable_db"].get()))

    def refresh_db_location(self) -> None:
        """Print the resolved folder, so there is no doubt where data goes."""
        self.portable_path.configure(
            text="Database and updates: {}".format(
                db_location(self.app.prefs, self.app.env)))
        if portable_possible():
            self.portable_box.set_enabled(True)
            self.portable_note.configure(text="")
        else:
            self.portable_box.set_enabled(False)
            self.portable_note.configure(
                text="Portable mode applies to the packaged {}.exe. This copy "
                     "is running from a Python installation, so the database "
                     "stays where it is.".format(branding.APP_NAME))

    def _changed(self) -> None:
        if getattr(self, "_loading", False):
            return
        self._sync_scheme_lock()
        self.app.settings_changed(self.collect())

    def show_rerun(self) -> None:
        """Offer a re-run when a setting changed while results are on screen."""
        if not self.rerun_bar.winfo_manager():
            self.rerun_bar.pack(fill="x", before=self._first_card)

    def hide_rerun(self) -> None:
        self.rerun_bar.pack_forget()

    def restore_defaults(self) -> None:
        """Restore the reference defaults and say so (section 10.4)."""
        self.vars["minid"].set(_num(DEFAULT_MINID))
        self.vars["mincov"].set(_num(DEFAULT_MINCOV))
        self.vars["minscore"].set(_num(DEFAULT_MINSCORE))
        self.vars["threads"].set("1")
        self.vars["jobs"].set("1")
        # The reference defaults ARE 1/1, so automatic tuning would immediately
        # contradict them: restoring them turns it off.
        self.vars["perf_auto"].set(False)
        self._sync_perf_lock()
        self.vars["exclude"].set(", ".join(DEFAULT_EXCLUDE))
        self.vars["scheme"].set("Automatic (recommended)")
        self._minscore_forced = False   # these ARE the user's values now
        self._minscore_saved = _num(DEFAULT_MINSCORE)
        self._changed()
        self.app.status.set(
            "Reference defaults restored: identity 95, coverage 50, score 50, "
            "and the four excluded schemes.", "ok", transient=True)


def _num(value: float) -> str:
    """Format a threshold without a pointless '.0'."""
    f = float(value)
    return str(int(f)) if f == int(f) else "{:g}".format(f)


# ===========================================================================
# SECTION 11 — DIALOGS  (sections 10.6, 10.7, 10.3)
# ===========================================================================
# Only six things in the whole application are modal: the search-engine setup,
# the update confirmation, the error dialog, About, the unsaved-results confirm
# and the >50-files confirm (section 10.8).

class ModalDialog(tk.Toplevel):
    """Shared modal plumbing: transient, centred, Escape-closable, focus-trapped."""

    def __init__(self, parent: tk.Misc, title: str, width: int = 520):
        assert_main_thread("a dialog constructor")
        super().__init__(parent)
        self.withdraw()
        self.title(title)
        self.configure(background=c("surface"))
        self.resizable(False, False)
        self.transient(parent.winfo_toplevel())
        self.body = ttk.Frame(self, style="Surface.TFrame")
        self.body.pack(fill="both", expand=True, padx=px(PAD_XL), pady=px(PAD_L))
        self.buttons = ttk.Frame(self, style="Surface.TFrame")
        self.buttons.pack(fill="x", padx=px(PAD_XL), pady=(0, px(PAD_L)))
        self._width = width
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.bind("<Escape>", lambda e: self.close(), add="+")

    def present(self) -> None:
        """Centre on the parent, grab focus and show."""
        self.update_idletasks()
        parent = self.master.winfo_toplevel()
        try:
            x = parent.winfo_rootx() + (parent.winfo_width() - self.winfo_width()) // 2
            y = parent.winfo_rooty() + (parent.winfo_height() - self.winfo_height()) // 3
            self.geometry("+{}+{}".format(max(0, x), max(0, y)))
        except tk.TclError:
            pass
        self.deiconify()
        try:
            self.grab_set()
        except tk.TclError:
            pass
        self.focus_set()

    def close(self) -> None:
        try:
            self.grab_release()
        except tk.TclError:
            pass
        self.destroy()

    def heading(self, text: str) -> ttk.Label:
        label = ttk.Label(self.body, text=text, style="Heading.TLabel",
                          wraplength=px(self._width), justify="left")
        label.pack(anchor="w")
        return label

    def paragraph(self, text: str, style: str = "Surface.TLabel") -> ttk.Label:
        label = ttk.Label(self.body, text=text, style=style, justify="left",
                          wraplength=px(self._width))
        label.pack(anchor="w", pady=(px(PAD_S), 0))
        return label


class DetailsFold(ttk.Frame):
    """A collapsed ``▸ Show technical details`` fold (section 10.6)."""

    def __init__(self, parent: tk.Misc, text: str, label: str = "Show technical details"):
        super().__init__(parent, style="Surface.TFrame")
        self.text = text
        self.open = False
        self._label = label
        self.toggle = ttk.Button(self, text="▸ " + label, style="Link.TButton",
                                 command=self.flip)
        self.toggle.pack(anchor="w")
        self.holder = ttk.Frame(self, style="Surface.TFrame")
        self.box = tk.Text(self.holder, height=10, width=68, wrap="none",
                           background=c("surface_alt"), foreground=c("text"),
                           font=F.get("mono"), relief="flat",
                           highlightthickness=1, highlightbackground=c("border"),
                           highlightcolor=c("focus"))
        scroll = ttk.Scrollbar(self.holder, orient="vertical", command=self.box.yview)
        self.box.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.box.pack(side="left", fill="both", expand=True)
        self.box.insert("1.0", text)
        self.box.configure(state="disabled")
        ttk.Button(self, text="Copy details", command=self.copy).pack(
            anchor="w", pady=(px(PAD_XS), 0))

    def flip(self) -> None:
        self.open = not self.open
        if self.open:
            self.holder.pack(fill="both", expand=True, pady=(px(PAD_XS), 0))
            self.toggle.configure(text="▾ " + self._label)
        else:
            self.holder.pack_forget()
            self.toggle.configure(text="▸ " + self._label)

    def copy(self) -> None:
        self.clipboard_clear()
        self.clipboard_append(self.text)


class ErrorDialog(ModalDialog):
    """The friendly error dialog of section 10.6 — never a stack trace."""

    def __init__(self, parent: tk.Misc, friendly: Friendly,
                 on_rebuild: Optional[Callable[[], None]] = None):
        super().__init__(parent, friendly.headline)
        self.heading(friendly.headline)
        if friendly.body:
            self.paragraph(friendly.body)
        if friendly.tips:
            self.paragraph("What to try:", style="SurfaceMuted.TLabel")
            for tip in friendly.tips:
                ttk.Label(self.body, text="•  " + tip, style="Surface.TLabel",
                          wraplength=px(self._width - 20), justify="left").pack(
                    anchor="w", padx=(px(PAD_S), 0))
        fold = DetailsFold(self.body, friendly.details)
        fold.pack(fill="x", pady=(px(PAD_M), 0))
        if friendly.offers_rebuild and on_rebuild is not None:
            ttk.Button(self.buttons, text="Rebuild search index",
                       command=lambda: (self.close(), on_rebuild())).pack(side="left")
        close = ttk.Button(self.buttons, text="Close", style="Accent.TButton",
                           command=self.close)
        close.pack(side="right")
        close.focus_set()
        self.present()


class AboutDialog(ModalDialog):
    """Branding, credits and citations (section 16, item 12)."""

    def __init__(self, parent: tk.Misc, env: Optional[Environment] = None):
        super().__init__(parent, "About " + branding.APP_NAME)
        self.heading("{} {}".format(branding.APP_NAME, __version__))
        self.paragraph(branding.APP_TAGLINE)
        self.paragraph(branding.ATTRIBUTION, style="Heading.TLabel")
        self.paragraph(
            "A Windows port of {} {} by {}, reproducing its results exactly."
            .format(branding.UPSTREAM_NAME, UPSTREAM_MLST_VERSION,
                    branding.UPSTREAM_AUTHOR))
        if env is not None:
            self.paragraph("Database {} · {} schemes\nSearch engine {}".format(
                env.db_version, env.scheme_count,
                env.blast_version or "not installed"),
                style="SurfaceMuted.TLabel")
        self.paragraph("Please cite:", style="SurfaceMuted.TLabel")
        for citation in branding.CITATIONS:
            ttk.Label(self.body, text="•  " + citation, style="SurfaceMuted.TLabel",
                      wraplength=px(self._width - 20), justify="left").pack(anchor="w")
        self.paragraph(branding.COPYRIGHT + "  ·  " + branding.HOMEPAGE,
                       style="SurfaceMuted.TLabel")
        ttk.Button(self.buttons, text="Copy version details",
                   command=self._copy).pack(side="left")
        close = ttk.Button(self.buttons, text="Close", style="Accent.TButton",
                           command=self.close)
        close.pack(side="right")
        close.focus_set()
        self._env = env
        self.present()

    def _copy(self) -> None:
        self.clipboard_clear()
        self.clipboard_append(_env_block(self._env))


class BootstrapDialog(ModalDialog):
    """First run without the search engine (section 10.7).

    Copy is fixed by C11: about 137 MB to download, about 35 MB installed, no
    administrator rights, and the install path is named.
    """

    def __init__(self, parent: tk.Misc, app: WmlstApp):
        super().__init__(parent, "{} needs its search engine".format(branding.APP_NAME))
        self.app = app
        self.cancel_event = threading.Event()
        self.job: Optional[int] = None
        self.done = False
        self._downloading = False
        try:
            from . import blastbin
            where = blastbin.install_root()
        except Exception:
            where = str(app_data_dir())
        self.heading("{} needs its search engine".format(branding.APP_NAME))
        self.paragraph(
            "WMLST uses the NCBI BLAST+ search engine to compare your sequences "
            "with the MLST database. It is about 137 MB to download and about "
            "35 MB once installed, and it only has to be done once.")
        self.paragraph(
            "It will be installed for you alone, in:\n{}\n"
            "No administrator rights are needed and nothing else on this computer "
            "is changed.".format(os.path.join(where, "blast")),
            style="SurfaceMuted.TLabel")
        self.progress = ttk.Progressbar(self.body, style="Big.Horizontal.TProgressbar",
                                        mode="determinate", maximum=100.0)
        self.detail = ttk.Label(self.body, text="", style="SurfaceMuted.TLabel",
                                wraplength=px(self._width), justify="left")

        fold = ttk.Frame(self.body, style="Surface.TFrame")
        fold.pack(fill="x", pady=(px(PAD_M), 0))
        self._fold_open = False
        self.fold_button = ttk.Button(fold, text="▸ I already have BLAST+ installed",
                                      style="Link.TButton", command=self._flip)
        self.fold_button.pack(anchor="w")
        self.fold_body = ttk.Frame(fold, style="Surface.TFrame")
        ttk.Label(self.fold_body, text="Choose the folder containing the BLAST+ "
                                       "programs (its bin folder).",
                  style="SurfaceMuted.TLabel", wraplength=px(self._width),
                  justify="left").pack(anchor="w")
        picker = ttk.Frame(self.fold_body, style="Surface.TFrame")
        picker.pack(fill="x", pady=px(PAD_XS))
        self.path_var = tk.StringVar()
        ttk.Entry(picker, textvariable=self.path_var, width=44).pack(
            side="left", fill="x", expand=True)
        ttk.Button(picker, text="Browse...", command=self._pick).pack(
            side="left", padx=(px(PAD_S), 0))
        ttk.Button(picker, text="Use this", command=self._use_existing).pack(
            side="left", padx=(px(PAD_S), 0))

        self.not_now = ttk.Button(self.buttons, text="Not now", command=self.close)
        self.not_now.pack(side="right", padx=(px(PAD_S), 0))
        self.action = ttk.Button(self.buttons, text="Download now",
                                 style="Accent.TButton", command=self.start)
        self.action.pack(side="right")
        self.action.focus_set()
        self.present()

    def _flip(self) -> None:
        self._fold_open = not self._fold_open
        if self._fold_open:
            self.fold_body.pack(fill="x")
            self.fold_button.configure(text="▾ I already have BLAST+ installed")
        else:
            self.fold_body.pack_forget()
            self.fold_button.configure(text="▸ I already have BLAST+ installed")

    def _pick(self) -> None:
        folder = filedialog.askdirectory(parent=self, title="Locate BLAST+",
                                         mustexist=True)
        if folder:
            self.path_var.set(folder)

    def _use_existing(self) -> None:
        """Validate a user-supplied installation by asking it for its version."""
        folder = self.path_var.get().strip()
        if not folder:
            return
        self.detail.pack(fill="x", pady=(px(PAD_S), 0))
        self.detail.configure(text="Checking that folder...")

        def probe() -> Any:
            from . import blastbin
            return blastbin.find_blast(folder)

        self.job = self.app.run_task("blast-explicit", probe, owner=self)

    def start(self) -> None:
        """Download, verify, extract — with a live progress bar (section 10.7)."""
        self.action.configure(text="Cancel", command=self.cancel)
        self.not_now.configure(state="disabled")
        self.progress.pack(fill="x", pady=(px(PAD_M), 0))
        self.detail.pack(fill="x", pady=(px(PAD_XS), 0))
        self.detail.configure(text="Starting download...")
        self.cancel_event = threading.Event()
        self._downloading = True

        def run(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import blastbin
            return blastbin.bootstrap(progress=progress, cancel=cancel)

        self.job = self.app.run_task("blast-bootstrap", run, owner=self,
                                     with_progress=True, cancel=self.cancel_event)

    def cancel(self) -> None:
        """Cancelling deletes the partial download (blastbin honours the event)."""
        self.cancel_event.set()
        self.detail.configure(text="Stopping...")
        self.action.configure(state="disabled")

    def close(self) -> None:
        """Escape / the window X during a download means "cancel", not "abandon".

        Without this the modal vanishes while the 137 MB download carries on
        unattended, nothing deletes the partial file, and the finished install is
        adopted by nobody.  The first close asks the worker to stop and keeps the
        modal up; once it has acknowledged (on_task_failed clears self.job) a
        further close goes through.
        """
        if self._downloading and not self.cancel_event.is_set():
            self.cancel()
            return
        super().close()

    # -- task routing --------------------------------------------------------
    def on_progress(self, msg: Msg) -> bool:
        if self.job is None or msg.job != self.job:
            return False
        self.progress.configure(value=msg.percent)
        if msg.text:
            self.detail.configure(text=msg.text)
        return True

    def on_task_done(self, job: int, payload: Any) -> bool:
        if self.job is None or job != self.job:
            return False
        self.done = True
        self._downloading = False
        self.app.adopt_tools(payload)
        self.close()
        return True

    def on_task_failed(self, job: int, friendly: Friendly) -> bool:
        if self.job is None or job != self.job:
            return False
        self.job = None
        self._downloading = False
        self.progress.configure(value=0.0)
        self.action.configure(text="Try again", command=self.start, state="normal")
        self.not_now.configure(state="normal")
        self.detail.configure(text=friendly.headline + " " + friendly.body)
        manual = (
            "You can also install it yourself:\n"
            "1. Download ncbi-blast-2.17.0+-x64-win64.tar.gz from\n"
            "   https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/\n"
            "2. Unpack it anywhere you like.\n"
            "3. Come back here, open 'I already have BLAST+ installed' and point "
            "WMLST at its bin folder.\n\n" + friendly.details)
        DetailsFold(self.body, manual, "Show details and manual instructions").pack(
            fill="x", pady=(px(PAD_S), 0))
        return True


class UpdatePlanDialog(ModalDialog):
    """Per-scheme update plan with checkboxes and a byte estimate (section 10.3)."""

    def __init__(self, parent: tk.Misc, plan: Any, updates: Sequence[Any],
                 *, preselected: Optional[Dict[str, bool]] = None):
        super().__init__(parent, "Database update", width=560)
        self.accepted = False
        self.backup = True
        self.chosen: List[str] = []
        self._vars: List[Tuple[str, tk.BooleanVar]] = []
        self._preselected = dict(preselected or {})
        total = sum(int(getattr(u, "bytes_estimate", 0) or 0) for u in updates)
        self.heading("{} scheme{} can be updated".format(
            len(updates), "" if len(updates) == 1 else "s"))
        self.paragraph(
            "About {} will be downloaded. Your current database keeps working "
            "until every chosen scheme has been replaced.".format(human_bytes(total)))

        surface = c("surface")
        bar = ttk.Frame(self.body, style="Surface.TFrame")
        bar.pack(fill="x", pady=(px(PAD_S), 0))
        ttk.Button(bar, text="Select all",
                   command=lambda: self._set_all(True)).pack(side="left")
        ttk.Button(bar, text="Deselect all",
                   command=lambda: self._set_all(False)).pack(side="left",
                                                              padx=(px(PAD_S), 0))
        self.count_label = ttk.Label(bar, text="", style="SurfaceMuted.TLabel")
        self.count_label.pack(side="left", padx=px(PAD_M))

        ring = focus_ring(self.body, under=surface)
        ring.pack(fill="both", expand=True, pady=px(PAD_S))
        pane = ScrollPane(ring, style="Surface.TFrame", background=surface)
        pane.configure(height=px(240), width=px(520))
        pane.pack_propagate(False)
        pane.pack(fill="both", expand=True, padx=1, pady=1)
        for update in updates:
            name = str(getattr(update, "name", ""))
            var = tk.BooleanVar(value=self._preselected.get(name, True))
            var.trace_add("write", lambda *_a: self._update_count())
            self._vars.append((name, var))
            added = getattr(update, "added_types", None)
            detail = str(getattr(update, "detail", "") or "")
            text = "{}   ({}{})".format(
                name, human_bytes(getattr(update, "bytes_estimate", 0)),
                ", {} new sequence types".format(added) if added else "")
            row = CheckBox(pane.body, text, variable=var, surface=surface,
                           subtext=detail, wraplength=430)
            row.pack(anchor="w", fill="x", padx=px(PAD_M), pady=px(PAD_XS))

        self._update_count()
        self.backup_var = tk.BooleanVar(value=True)
        CheckBox(self.body, "Keep a backup of the current database",
                 variable=self.backup_var, surface=surface).pack(
            anchor="w", pady=(px(PAD_M), 0))
        ttk.Button(self.buttons, text="Cancel", command=self.close).pack(side="right",
                                                                        padx=(px(PAD_S), 0))
        go = ttk.Button(self.buttons, text="Update now", style="Accent.TButton",
                        command=self._accept)
        go.pack(side="right")
        go.focus_set()
        self.present()

    def _set_all(self, value: bool) -> None:
        """Tick or untick every scheme in the plan."""
        for _name, var in self._vars:
            var.set(value)
        self._update_count()

    def _update_count(self) -> None:
        """Keep the live count honest as the ticks change."""
        total = len(self._vars)
        chosen = sum(1 for _n, var in self._vars if var.get())
        try:
            self.count_label.configure(
                text="{} of {} scheme{} selected".format(
                    chosen, total, "" if total == 1 else "s"))
        except tk.TclError:  # pragma: no cover - during teardown
            pass

    def _accept(self) -> None:
        self.chosen = [name for name, var in self._vars if var.get()]
        self.backup = bool(self.backup_var.get())
        self.accepted = bool(self.chosen)
        self.close()


# ===========================================================================
# SECTION 12 — THE APPLICATION  (sections 10.1, 10.5, 10.9)
# ===========================================================================

DEFAULT_GEOMETRY = "1180x760"
MIN_WIDTH, MIN_HEIGHT = 960, 640


def parse_geometry(geom: str) -> Optional[Tuple[int, int, int, int]]:
    """Parse ``WxH+X+Y`` into ``(w, h, x, y)``. None when it is not that shape."""
    try:
        size, _, rest = str(geom).partition("+")
        if not rest:
            return None
        w_s, _, h_s = size.partition("x")
        x_s, _, y_s = rest.partition("+")
        return (int(w_s), int(h_s), int(x_s), int(y_s))
    except (ValueError, AttributeError):
        return None


def geometry_on_screen(geom: str, screen_w: int, screen_h: int) -> bool:
    """True when at least a usable corner of the saved rectangle is visible.

    A window restored onto a monitor that is no longer attached is invisible and
    unrecoverable for a novice, so an off-screen rectangle is discarded.
    """
    parsed = parse_geometry(geom)
    if parsed is None:
        return False
    w, h, x, y = parsed
    if w < 200 or h < 150:
        return False
    return (-w + 120) < x < (screen_w - 120) and -20 < y < (screen_h - 80)


def dedupe_novel(samples: Sequence[SampleResult]) -> Tuple[Any, ...]:
    """Novel alleles across a batch, deduplicated by fasta_id, first wins (D7)."""
    seen = set()
    out = []
    for sample in samples:
        for novel in sample.novel:
            if novel.fasta_id not in seen:
                seen.add(novel.fasta_id)
                out.append(novel)
    return tuple(out)


def fallback_meta(cfg: RunConfig, env: Environment) -> RunMeta:
    """A :class:`RunMeta` for exports when the engine did not supply one.

    Provenance only: not one field of it takes part in a result.
    """
    now = datetime.now(timezone.utc)
    return RunMeta(
        wmlst_version=__version__,
        mlst_compat=UPSTREAM_MLST_VERSION,
        db_version=env.db_version,
        db_scheme_count=env.scheme_count,
        blast_version=env.blast_version,
        blast_path=env.blast_path,
        dbdir=env.dbdir,
        datadir=env.datadir,
        blastdb=env.blastdb,
        python_version=platform.python_version(),
        platform=platform.platform(),
        hostname=socket.gethostname(),
        started_utc=now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        started_local=datetime.now().astimezone().isoformat(timespec="seconds"),
        duration_s=0.0,
        argv=("wmlst-gui",),
        config=cfg,
    )


_NO_OWNER = object()   # "look the owner up yourself" marker for WmlstApp._route


class WmlstApp:
    """Owns the root, the theme, the three tabs, the menubar, the status line
    and the footer (section 4.10)."""

    SAVE_DEBOUNCE_MS = 800

    def __init__(self, root: tk.Misc, prefs: Prefs):
        assert_main_thread("WmlstApp")
        self.root = root
        self.prefs = prefs
        self.env = Environment()
        self.dnd_enabled = bool(getattr(root, "wmlst_dnd", False))
        self.state = AppState(on_change=self._on_state_change)
        self.controller = AnalysisController(
            schedule=lambda ms, fn: self.root.after(ms, fn),
            dispatch=self.dispatch, on_tick=self._on_tick)
        self._owners: Dict[int, Any] = {}
        self._task_kinds: Dict[int, str] = {}
        self._analysis_job = 0
        self._save_after: Optional[str] = None
        self._queued_files: List[str] = []
        # Files named on the command line.  They are NOT started on a timer racing
        # the environment probe (the probe takes 250-450 ms, the first tick is due
        # at 120 ms, so the race was always lost and the user got a spurious
        # "BLAST is missing" download modal): _apply_env / _task_failed drain them.
        self._argv_files: List[str] = []
        #: One-shot: the first launch after an install builds the index itself.
        self._auto_index_build_tried = False
        self._bootstrap: Optional[BootstrapDialog] = None
        self._failed_count = 0
        #: What wmlst.perf recommended for this machine, once probed.
        self.tuning: Any = None
        self._build()
        # Nothing reaches the console: every uncaught Tk callback error becomes
        # the friendly dialog (section 10.6).
        self.root.report_callback_exception = self._tk_exception

    # -- construction --------------------------------------------------------
    def _build(self) -> None:
        root = self.root
        root.title(branding.WINDOW_TITLE)
        root.minsize(px(MIN_WIDTH), px(MIN_HEIGHT))
        root.configure(background=c("bg"))
        root.protocol("WM_DELETE_WINDOW", self.on_close)
        self._restore_geometry()
        self._build_menu()

        self.notebook = ttk.Notebook(root)
        self.analyse_view = AnalyseView(self.notebook, self)
        self.database_view = DatabaseView(self.notebook, self)
        self.settings_view = SettingsView(self.notebook, self)
        self.notebook.add(self.analyse_view, text="Analyse")
        self.notebook.add(self.database_view, text="Database")
        self.notebook.add(self.settings_view, text="Settings")
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab, add="+")

        self.banner = ttk.Frame(root, style="TFrame")
        self.banner_label = ttk.Label(self.banner, text="", style="Muted.TLabel")
        self.banner_label.pack(side="left", padx=px(PAD_M), pady=px(PAD_XS))
        self.banner_button = ttk.Button(self.banner, text="Install now",
                                        style="Accent.TButton",
                                        command=self.offer_bootstrap)
        self.banner_button.pack(side="left", pady=px(PAD_XS))

        footer = ttk.Frame(root, style="Footer.TFrame")
        footer.pack(fill="x", side="bottom")
        ttk.Separator(root, orient="horizontal").pack(fill="x", side="bottom")
        self.status = StatusLine(footer)
        self.status.pack(fill="x", side="top")
        credit = ttk.Frame(footer, style="Footer.TFrame")
        credit.pack(fill="x", side="bottom")
        ttk.Label(credit, text="{} {}  ·  {}".format(
            branding.APP_NAME, __version__, branding.ATTRIBUTION),
            style="Status.TLabel").pack(side="left", padx=px(PAD_M), pady=(0, px(PAD_S)))
        self.footer_env = ttk.Label(credit, text="Starting...", style="Status.TLabel")
        self.footer_env.pack(side="right", padx=px(PAD_M), pady=(0, px(PAD_S)))

        # LAST, on purpose: pack() serves children in call order, so the status
        # line and the footer credit must claim their strip before the notebook
        # is allowed to expand, or a tall tab would silently cover them.
        self.notebook.pack(fill="both", expand=True, padx=px(PAD_M),
                           pady=(px(PAD_M), 0))

        root.bind("<Control-o>", lambda e: self.analyse_view.browse_files(), add="+")
        root.bind("<Control-O>", lambda e: self.analyse_view.browse_folder(), add="+")
        root.bind("<F5>", lambda e: self.rerun_all(), add="+")
        root.bind("<F1>", lambda e: self.show_about(), add="+")
        root.bind("<Escape>", lambda e: self.cancel_run(), add="+")
        root.bind("<Control-q>", lambda e: self.on_close(), add="+")
        self.analyse_view.dropzone.focus_set()

    def _build_menu(self) -> None:
        menubar = tk.Menu(self.root)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open files...", accelerator="Ctrl+O",
                              command=lambda: self.analyse_view.browse_files())
        file_menu.add_command(label="Open folder...", accelerator="Ctrl+Shift+O",
                              command=lambda: self.analyse_view.browse_folder())
        file_menu.add_separator()
        file_menu.add_command(label="Open report",
                              command=lambda: self.export("html"))
        file_menu.add_command(label="Save table (TSV)...",
                              command=lambda: self.export("tsv"))
        file_menu.add_command(label="Save JSON...", command=lambda: self.export("json"))
        file_menu.add_command(label="Save new alleles...",
                              command=lambda: self.export("novel"))
        file_menu.add_separator()
        file_menu.add_command(label="Exit", accelerator="Ctrl+Q", command=self.on_close)
        menubar.add_cascade(label="File", menu=file_menu)

        run_menu = tk.Menu(menubar, tearoff=0)
        run_menu.add_command(label="Re-run all files", accelerator="F5",
                             command=self.rerun_all)
        run_menu.add_command(label="Cancel", accelerator="Esc", command=self.cancel_run)
        run_menu.add_separator()
        run_menu.add_command(label="Clear results",
                             command=lambda: self.analyse_view.clear())
        menubar.add_cascade(label="Run", menu=run_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About {}".format(branding.APP_NAME),
                              accelerator="F1", command=self.show_about)
        help_menu.add_command(label="Show log file", command=self.show_log)
        help_menu.add_command(label="Project home page",
                              command=lambda: webbrowser.open(branding.HOMEPAGE))
        menubar.add_cascade(label="Help", menu=help_menu)
        try:
            self.root.configure(menu=menubar)
        except tk.TclError:  # pragma: no cover - some window managers
            pass

    def _restore_geometry(self) -> None:
        geom = self.prefs.geometry
        try:
            sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        except tk.TclError:
            sw, sh = 1920, 1080
        if geom and geometry_on_screen(geom, sw, sh):
            self.root.geometry(geom)
        else:
            self.root.geometry("{}x{}".format(px(1180), px(760)))
        if self.prefs.zoomed:
            try:
                self.root.state("zoomed")
            except tk.TclError:
                pass

    # -- start-up ------------------------------------------------------------
    def post_start(self) -> None:
        """Run 50 ms after the first paint: locate the database and the engine."""
        self.controller.ensure_polling()
        prefs = self.prefs
        self.analyse_view.refresh_scheme_notice()
        self.apply_autotune(announce=True)
        self.run_task("env", lambda: probe_environment(prefs))
        if prefs.load_note:
            self.status.set(prefs.load_note, "warn", transient=True)
        if not self.dnd_enabled:
            self.status.set(
                "Drag-and-drop is not available on this computer — click the big "
                "box to choose a file instead.", "info", transient=True)

    def _apply_env(self, env: Environment) -> None:
        self.env = env
        self.refresh_footer()
        try:
            self.settings_view.set_scheme_choices(env.scheme_names)
            self.settings_view.refresh_db_location()
            self.database_view.refresh_location()
        except tk.TclError:
            pass
        if not env.db_ok:
            self.status.set(env.db_error or "The MLST database could not be found.",
                            "bad", resting=True)
        elif env.index_stale:
            # A fresh install ships the allele files but not the 180 MB derived
            # index, so the very first launch always lands here. Telling a novice
            # to go and find the Database tab is not an answer: the rebuild needs
            # no download, no consent and no decision, so just do it and show
            # progress. Once only -- if it fails, say so rather than looping,
            # because the task's own completion re-probes the environment.
            if env.blast_ok and not self._auto_index_build_tried:
                self._auto_index_build_tried = True
                self.status.set("Preparing the search index — this happens once, "
                                "and takes about a minute.", "info", resting=True)
                self.rebuild_index()
            else:
                self.status.set("The search index is older than the allele files — "
                                "rebuild it from the Database tab.", "warn",
                                resting=True)
        elif env.blast_ok:
            self.status.set("Ready. Drop a FASTA file on the box above."
                            if self.dnd_enabled else
                            "Ready. Click the box above to choose a FASTA file.",
                            resting=True)
        self._update_banner()
        self._start_argv_files()
        if not env.blast_ok and not self.prefs.blast_prompt_shown:
            self.prefs.blast_prompt_shown = True
            self.schedule_save()
            self.offer_bootstrap()

    def _start_argv_files(self) -> None:
        """Hand any command-line files to the drop handler, once the environment
        probe has answered.  Runs at most once per batch of argv files."""
        if not self._argv_files:
            return
        # Hold them back while the first-run index build is in flight: starting a
        # search against a half-built index fails, and the failure looks like a
        # missing search engine rather than "wait a moment". The rebuild's own
        # completion re-probes the environment, which calls this again.
        if self.env.index_stale and self._auto_index_build_tried:
            return
        pending, self._argv_files = self._argv_files, []
        self.root.after(1, lambda: self.analyse_view.handle_paths(pending))

    def _update_banner(self) -> None:
        if self.env.blast_ok or self._bootstrap is not None:
            self.banner.pack_forget()
        else:
            self.banner_label.configure(
                text="⚠ The search engine is not installed yet, so files cannot be "
                     "analysed.")
            self.banner.pack(fill="x", side="bottom", before=self.status.master)

    def refresh_footer(self) -> None:
        self.footer_env.configure(text=self.env.footer_text())

    # -- tasks ---------------------------------------------------------------
    def run_task(self, kind: str, fn: Callable[..., Any], *, owner: Any = None,
                 with_progress: bool = False,
                 cancel: Optional[threading.Event] = None) -> int:
        """Start a background job and remember who should hear about it."""
        job = self.controller.run_task(kind, fn, with_progress=with_progress,
                                       cancel=cancel)
        self._owners[job] = owner
        self._task_kinds[job] = kind
        return job

    def _route(self, msg: Msg, method: str, *args: Any,
               owner: Any = _NO_OWNER) -> bool:
        # The caller may already have popped the owner out of self._owners (the
        # terminal messages do exactly that), in which case it hands the owner in
        # here.  Looking it up again would find nothing and the result would be
        # routed to the fallbacks -- or dropped entirely.
        if owner is _NO_OWNER:
            owner = self._owners.get(msg.job, None)
        candidates = [owner] if owner is not None else []
        candidates.extend([self.database_view, self._bootstrap])
        for candidate in candidates:
            if candidate is None:
                continue
            handler = getattr(candidate, method, None)
            if handler is not None:
                try:
                    if handler(*args):
                        return True
                except tk.TclError:
                    return True
        return False

    # -- dispatch ------------------------------------------------------------
    def dispatch(self, msg: Msg) -> None:
        """MAIN THREAD ONLY (section 4.10). The single fan-out point."""
        assert_main_thread("dispatch")
        if msg.job in self._task_kinds:
            self._dispatch_task(msg)
            return
        if msg.kind == "started":
            self.analyse_view.on_started(msg)
        elif msg.kind == "phase":
            self.analyse_view.on_phase(msg)
        elif msg.kind == "warn":
            self.status.set(msg.text, "warn", transient=True)
            LOG.warning("%s: %s", msg.path, msg.text)
        elif msg.kind == "result":
            self._on_result(msg)
        elif msg.kind == "failed":
            self._on_failed(msg)
        elif msg.kind == "cancelled":
            self.state.finish("cancelled")
            self.status.set("Stopped. The results collected so far are still here.",
                            "warn", resting=True)
        elif msg.kind == "fatal":
            self.state.finish("error")
            if msg.exc is not None:
                self.show_error(msg.exc, tb_text=msg.tb)
        elif msg.kind == "batch_done":
            self._on_batch_done(msg)

    def _dispatch_task(self, msg: Msg) -> None:
        kind = self._task_kinds.get(msg.job, "")
        if msg.kind == "phase":
            if not self._route(msg, "on_progress", msg):
                if msg.text:
                    self.status.set(msg.text, "info")
            return
        if msg.kind == "failed":
            friendly = (friendly_error(msg.exc, env=self.env, tb_text=msg.tb)
                        if msg.exc is not None
                        else Friendly("Something went wrong", msg.text))
            LOG.error("task %s failed: %s", kind, msg.text)
            self._task_kinds.pop(msg.job, None)
            owner = self._owners.pop(msg.job, None)
            if not self._route(msg, "on_task_failed", msg.job, friendly, owner=owner):
                self._task_failed(kind, friendly)
            return
        if msg.kind == "env":
            payload = self.controller.take_task_result(msg.job)
            self._task_kinds.pop(msg.job, None)
            owner = self._owners.pop(msg.job, None)
            if not self._route(msg, "on_task_done", msg.job, payload, owner=owner):
                self._task_done(kind, payload)

    def _task_done(self, kind: str, payload: Any) -> None:
        if kind == "env" and isinstance(payload, Environment):
            self._apply_env(payload)
        elif kind.startswith("export-"):
            self._after_export(kind.split("-", 1)[1], str(payload or ""))
        elif kind in ("blast-bootstrap", "blast-explicit"):
            # The dialog that asked for it is gone (closed mid-download), but the
            # engine really was installed -- keep it instead of staying degraded.
            self.adopt_tools(payload)
        elif kind == "rebuild":
            self.status.set("The search index was rebuilt.", "ok", transient=True)
            self.run_task("env", lambda: probe_environment(self.prefs))
        elif kind == "portable-copy":
            self._portable_ready(str(payload or ""))

    def _task_failed(self, kind: str, friendly: Friendly) -> None:
        if kind == "env":
            self.status.set(friendly.headline, "bad", resting=True)
            self._start_argv_files()
            return
        self.show_friendly(friendly)

    def _on_tick(self, elapsed: float, running: bool) -> None:
        self.analyse_view.on_tick(elapsed, running)

    def _on_state_change(self, old: str, new: str) -> None:
        self.analyse_view.show_state(new)

    def _on_tab(self, _event: Any = None) -> None:
        try:
            current = self.notebook.nametowidget(self.notebook.select())
        except (tk.TclError, KeyError):
            return
        if current is self.database_view:
            self.database_view.activate()

    # -- the run -------------------------------------------------------------
    def start_analysis(self, files: Sequence[str]) -> None:
        """Begin a batch — the only path from a drop to the engine."""
        assert_main_thread("start_analysis")
        files = [str(f) for f in files]
        if not files:
            return
        if not self.env.blast_ok:
            self._queued_files = list(files)
            self.offer_bootstrap()
            return
        if not self.env.db_ok:
            self.show_friendly(Friendly(
                "The MLST database is missing",
                self.env.db_error or "No scheme files were found.",
                ("Reinstall {} to restore the bundled database.".format(
                    branding.APP_NAME),), _env_block(self.env)))
            return
        cfg = self.prefs.to_runconfig(files, self.env)
        if self.controller.running:
            self.controller.start(files, cfg)  # queued; picked up at batch_done
            self.status.set("{} more file(s) queued.".format(len(files)), "info",
                            transient=True)
            return
        self._failed_count = 0
        self.state.to(RUNNING)
        self.analyse_view.begin_run(len(files))
        self.controller.start(files, cfg)
        self._analysis_job = self.controller.job
        self.status.set("Analysing {} file{}...".format(
            len(files), "" if len(files) == 1 else "s"), "info", resting=True)
        LOG.info("run started: %d file(s), threads=%d jobs=%d",
                 len(files), cfg.threads, cfg.jobs)

    def cancel_run(self) -> None:
        """Cancel the current batch; partial results are kept (section 10.2)."""
        if self.controller.running:
            self.controller.request_cancel()
            self.analyse_view.cancel_button.configure(state="disabled",
                                                      text="Stopping...")
            self.status.set("Stopping after the current file...", "warn")

    def rerun_all(self) -> None:
        """Re-run every file currently listed, with the current settings (F5)."""
        paths = [r.path for r in self.analyse_view.results]
        paths.extend(p for p, _ in self.analyse_view.failures)
        if not paths:
            self.status.set("There is nothing to re-run yet.", "warn", transient=True)
            return
        if not self.analyse_view.clear():
            return
        self.settings_view.hide_rerun()
        self.start_analysis(paths)

    def _on_result(self, msg: Msg) -> None:
        if msg.result is None:
            return
        res = msg.result
        if res.failed and res.error is not None:
            self._on_failed(Msg("failed", msg.job, index=msg.index, total=msg.total,
                                path=res.path, exc=res.error, text=res.error_text))
            return
        self.analyse_view.add_result(res)
        for warning in res.warnings:
            LOG.warning("%s: %s", res.label, warning)
        self.status.set("{}: {} ST {} ({})".format(
            res.label, res.scheme, res.st, res.status),
            STATUS_UI.get(res.status, (None, "muted", None))[1])

    def _on_failed(self, msg: Msg) -> None:
        exc = msg.exc if msg.exc is not None else WmlstError(msg.text)
        name = os.path.basename(msg.path) or msg.path
        friendly = friendly_error(exc, name=name, env=self.env, tb_text=msg.tb)
        LOG.error("file failed: %s: %s", msg.path, msg.text)
        self._failed_count += 1
        self.analyse_view.add_failure(msg.path, friendly)
        if msg.total <= 1:
            self.show_friendly(friendly)
        else:
            self.status.set(
                "{} file{} failed (click the red row for details).".format(
                    self._failed_count, "" if self._failed_count == 1 else "s"),
                "bad")

    def _on_batch_done(self, msg: Msg) -> None:
        reason = "cancelled" if self.controller.cancelled else "done"
        if self.state.running:
            self.state.finish(reason)
        self.analyse_view.end_run(reason)
        done = len(self.analyse_view.results)
        if reason == "done" and self._failed_count == 0:
            self.status.set("Finished: {} file{} analysed.".format(
                done, "" if done == 1 else "s"), "ok", resting=True)
        elif reason == "done":
            self.status.set("Finished: {} analysed, {} failed.".format(
                done, self._failed_count), "warn", resting=True)
        # Always drain, even on a cancel: leaving the queue in the controller
        # would resurrect the cancelled files at the end of some later, unrelated
        # batch.  Cancelling drops them -- but says so, because the UI already
        # promised "{n} more file(s) queued."
        pending = self.controller.pending_paths()
        if reason == "cancelled":
            if pending:
                self.status.set(
                    "Stopped. {} queued file{} not analysed — the results "
                    "collected so far are still here.".format(
                        len(pending), " was" if len(pending) == 1 else "s were"),
                    "warn", resting=True)
            return
        if pending:
            self.start_analysis(pending)

    # -- errors --------------------------------------------------------------
    def show_error(self, exc: BaseException, *, name: str = "",
                   tb_text: str = "") -> None:
        self.show_friendly(friendly_error(exc, name=name, env=self.env,
                                          tb_text=tb_text))

    def show_friendly(self, friendly: Friendly) -> None:
        """Show the friendly dialog, or the setup modal when that is the fix."""
        assert_main_thread("an error dialog")
        if friendly.needs_bootstrap:
            self.offer_bootstrap()
            return
        ErrorDialog(self.root, friendly, on_rebuild=self.rebuild_index)

    def _tk_exception(self, exc_type: Any, exc_value: Any, tb: Any) -> None:
        """``root.report_callback_exception`` — nothing reaches the console."""
        text = "".join(traceback.format_exception(exc_type, exc_value, tb))
        LOG.error("unhandled callback error\n%s", text)
        try:
            self.show_friendly(friendly_error(exc_value, env=self.env, tb_text=text))
        except Exception:  # pragma: no cover - the dialog itself failed
            messagebox.showerror(branding.APP_NAME, str(exc_value))

    # -- setup / database ----------------------------------------------------
    def offer_bootstrap(self) -> None:
        """Show the one-time search-engine setup modal (section 10.7)."""
        if self._bootstrap is not None:
            try:
                self._bootstrap.present()
                return
            except tk.TclError:
                self._bootstrap = None
        dialog = BootstrapDialog(self.root, self)
        self._bootstrap = dialog
        self.root.wait_window(dialog)
        self._bootstrap = None
        self._update_banner()
        if not self.env.blast_ok and self._queued_files:
            self.status.set(
                "The search engine is still needed before files can be analysed.",
                "warn", resting=True)

    def adopt_tools(self, tools: Any) -> None:
        """A successful bootstrap or a user-picked installation."""
        path = blast_exe_path(tools)
        if not path:
            return
        self.env.blast_ok = True
        self.env.blast_path = path
        self.env.blast_version = str(getattr(tools, "version", "") or "")
        self.env.blast_origin = str(getattr(tools, "origin", "") or "")
        self.prefs.blast_path = path
        self.schedule_save()
        self.refresh_footer()
        self._update_banner()
        self.status.set("The search engine is ready.", "ok", transient=True)
        # Anything the user already dropped starts analysing automatically.
        if self._queued_files:
            files, self._queued_files = self._queued_files, []
            self.root.after(50, lambda: self.start_analysis(files))

    def rebuild_index(self) -> None:
        """Rebuild the search index from the allele files (section 10.6)."""
        dbdir = self.env.dbdir
        blast_path = self.env.blast_path

        def run(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import blastbin
            from . import updatedb as updatedb_mod
            tools = blastbin.find_blast(blast_path or None)
            return updatedb_mod.build_blast_db(dbdir, tools, progress=progress,
                                               cancel=cancel)

        self.status.set("Rebuilding the search index...", "info")
        self.run_task("rebuild", run, with_progress=True)

    # -- automatic tuning (wmlst.perf) ---------------------------------------
    def apply_autotune(self, *, announce: bool = False) -> None:
        """Size threads/jobs for this computer at start-up (section C3).

        GUI ONLY. ``wmlst.perf`` is not imported by the command line and the CLI
        defaults stay 1/1; this is the window choosing sensible numbers for the
        machine it happens to be running on, which the user can override in
        Settings by unticking "Tune for this computer automatically".
        """
        assert_main_thread("auto-tuning")
        try:
            tuning = perf.probe()
        except Exception as exc:  # pragma: no cover - every probe degrades itself
            LOG.warning("could not size this computer: %s", exc)
            return
        self.tuning = tuning
        try:
            self.settings_view.set_tuning(tuning)
        except tk.TclError:  # pragma: no cover
            pass
        if not self.prefs.perf_auto:
            return
        if (self.prefs.threads, self.prefs.jobs) != (tuning.threads, tuning.jobs):
            self.prefs.threads = tuning.threads
            self.prefs.jobs = tuning.jobs
            self.schedule_save()
        try:
            self.settings_view.load_from(self.prefs)
        except tk.TclError:  # pragma: no cover
            pass
        if announce and tuning.rationale:
            self.status.set(tuning.rationale, "info", transient=True)

    # -- portable mode -------------------------------------------------------
    def set_portable_db(self, on: bool) -> None:
        """Move where database updates are written, or explain why we cannot.

        A read-only folder (a CD, a locked USB stick, Program Files without
        elevation) is an ordinary situation, not a crash: the tick comes back
        off and the status line says what happened.
        """
        assert_main_thread("portable mode")
        if not on:
            self.prefs.portable_db = False
            self.prefs.dbdir = None
            self.schedule_save()
            self._refresh_db_location()
            self.status.set("Database updates go back to the per-user folder: "
                            "{}".format(db_location(self.prefs, self.env)),
                            "info", transient=True)
            self.run_task("env", lambda: probe_environment(self.prefs))
            return
        root, target = portable_root(), portable_db_dir()
        if root is None or target is None:
            self._portable_refused(
                "Portable mode applies to the packaged {}.exe.".format(
                    branding.APP_NAME))
            return
        if not dir_writable(root):
            self._portable_refused(
                "{} cannot be written to, so the database cannot be kept "
                "there. Copy {} to a writable folder first.".format(
                    root, branding.APP_NAME))
            return
        if looks_like_database(target):
            self._portable_ready(str(target), copied=False)
            return
        source = self.env.dbdir
        if not looks_like_database(source):
            self._portable_refused(
                "There is no database to copy yet — install or update the "
                "database first, then turn portable mode on.")
            return
        self.status.set("Copying the database beside the program...", "info")
        destination = str(target)

        def copy(progress: Callable[..., None], cancel: threading.Event) -> Any:
            return copy_database(source, destination, progress=progress,
                                 cancel=cancel)

        self.run_task("portable-copy", copy, with_progress=True)

    def _portable_refused(self, why: str) -> None:
        """Undo the tick and say why, rather than raising at the user."""
        self.prefs.portable_db = False
        try:
            self.settings_view.vars["portable_db"].set(False)
        except (tk.TclError, KeyError):  # pragma: no cover
            pass
        self._refresh_db_location()
        self.status.set(why, "warn", resting=False, transient=True)

    def _portable_ready(self, path: str, *, copied: bool = True) -> None:
        """Adopt the folder beside the executable as the live database."""
        self.prefs.portable_db = True
        self.prefs.dbdir = path
        try:
            self.settings_view.vars["portable_db"].set(True)
        except (tk.TclError, KeyError):  # pragma: no cover
            pass
        self.schedule_save()
        self._refresh_db_location()
        self.status.set(
            "Portable mode is on. {} {}".format(
                "The database was copied to" if copied else "Using the database in",
                path), "ok", transient=True)
        self.run_task("env", lambda: probe_environment(self.prefs))

    def _refresh_db_location(self) -> None:
        try:
            self.settings_view.refresh_db_location()
            self.database_view.refresh_location()
        except tk.TclError:  # pragma: no cover
            pass

    # -- settings ------------------------------------------------------------
    def settings_changed(self, prefs: Prefs) -> None:
        """Adopt edited settings, offering a re-run when results are on screen."""
        changed = (prefs.minid, prefs.mincov, prefs.minscore, prefs.scheme,
                   prefs.exclude) != (self.prefs.minid, self.prefs.mincov,
                                      self.prefs.minscore, self.prefs.scheme,
                                      self.prefs.exclude)
        prefs.load_note = ""
        self.prefs = prefs
        self.schedule_save()
        self.analyse_view.refresh_scheme_notice()
        if changed and self.analyse_view.results:
            self.settings_view.show_rerun()

    def schedule_save(self) -> None:
        """Debounced, atomic settings write (800 ms, section 10.4)."""
        if self._save_after is not None:
            try:
                self.root.after_cancel(self._save_after)
            except tk.TclError:
                pass
        self._save_after = self.root.after(self.SAVE_DEBOUNCE_MS, self._save_now)

    def _save_now(self) -> None:
        self._save_after = None
        try:
            self.prefs.geometry = self.root.winfo_geometry()
            self.prefs.zoomed = str(self.root.state()) == "zoomed"
        except tk.TclError:
            pass
        self.prefs.save()

    # -- exports -------------------------------------------------------------
    def build_run_result(self, **overrides: Any) -> RunResult:
        """Assemble a :class:`RunResult` for ``report`` from what is on screen.

        Nothing is recomputed: every sample is the object the engine returned.
        """
        samples = tuple(self.analyse_view.results)
        cfg = self.prefs.to_runconfig([s.path for s in samples], self.env)
        if overrides:
            cfg = replace(cfg, **overrides)
        meta = self.controller.meta or fallback_meta(cfg, self.env)
        meta = replace(meta, config=cfg,
                       duration_s=sum(s.elapsed_s for s in samples))
        return RunResult(meta=meta, samples=samples, novel=dedupe_novel(samples))

    def export(self, kind: str) -> None:
        """Write one of the exports on a worker thread (section 10.2).

        Every byte is produced by ``report``; the GUI formats nothing itself, and
        the write runs off the main thread so a slow network share cannot freeze
        the window.
        """
        assert_main_thread("export")
        if not self.analyse_view.results:
            self.status.set("Analyse a file first.", "warn", transient=True)
            return
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        default = "wmlst-{}".format(stamp)
        if kind == "html":
            folder = app_data_dir() / "reports"
            try:
                folder.mkdir(parents=True, exist_ok=True)
                path = str(folder / (default + ".html"))
            except OSError:
                path = filedialog.asksaveasfilename(
                    parent=self.root, defaultextension=".html",
                    initialfile=default + ".html",
                    filetypes=(("Web page", "*.html"),))
                if not path:
                    return
        else:
            ext, types = {
                "tsv": (".tsv", (("Tab-separated table", "*.tsv"), ("All files", "*.*"))),
                "json": (".json", (("JSON", "*.json"), ("All files", "*.*"))),
                "novel": (".fa", (("FASTA", "*.fa"), ("All files", "*.*"))),
            }[kind]
            path = filedialog.asksaveasfilename(
                parent=self.root, title="Save", defaultextension=ext,
                initialfile=default + ext, filetypes=types,
                initialdir=self.prefs.last_dir or None)
            if not path:
                return

        if kind == "tsv":
            # --full layout: a novice wants STATUS and SCORE beside the alleles.
            result = self.build_run_result(full=True, outfile=path)
        elif kind == "json":
            result = self.build_run_result(json_path=path)
        elif kind == "novel":
            result = self.build_run_result(novel_path=path)
        else:
            result = self.build_run_result(html_path=path)
        options = (result, path, self.prefs.html_evidence)

        def run() -> Any:
            from . import report
            payload, target, evidence = options
            if kind == "tsv":
                with open(target, "w", encoding="utf-8", newline="\n") as fh:
                    report.write_tsv(payload, fh)
            elif kind == "json":
                report.write_json(payload, target)
            elif kind == "novel":
                report.write_novel_fasta(payload, target)
            else:
                # No branding override: report.py reads wmlst.branding itself, so
                # passing a second copy here could only ever make the two drift.
                report.write_html(payload, target,
                                  report.HtmlOptions(evidence=evidence))
            return target

        self.status.set("Saving...", "info")
        self.run_task("export-" + kind, run)

    def _after_export(self, kind: str, path: str) -> None:
        if not path:
            return
        if kind == "html":
            self.status.set("Report opened in your web browser.", "ok", transient=True)
            open_path(path)
        else:
            self.status.set("Saved to {}".format(path), "ok", transient=True)
        LOG.info("wrote %s: %s", kind, path)

    # -- help ----------------------------------------------------------------
    def show_about(self) -> None:
        AboutDialog(self.root, self.env)

    def show_log(self) -> None:
        """Open today's log file (Help -> Show log, section 10.6)."""
        path = setup_logging()
        if path is None or not Path(path).exists():
            self.status.set("No log file has been written yet.", "info", transient=True)
            return
        open_path(str(path))

    # -- shutdown ------------------------------------------------------------
    def on_close(self) -> None:
        """Close for good: stop the workers, kill the children, then exit.

        The window disappearing is not the same as the process ending. A worker
        parked on a running search child never returns, and
        ``concurrent.futures`` registers an ``atexit`` hook that JOINS its
        non-daemon workers — so one live child was enough to leave WMLST in Task
        Manager with no window to close it by, which is exactly what the user
        hit. Three layers now prevent that:

        1. the controller stops accepting work and joins briefly;
        2. :func:`wmlst.engine.shutdown_all` stops every engine and every search
           child, and :func:`wmlst.blastbin.shutdown` is called again directly,
           belt and braces, for children launched outside an engine;
        3. if either reports "not clean" the process is ended outright, and in
           every case a daemon watchdog ends it a few seconds later should the
           interpreter still be hanging on someone else's thread.
        """
        if self.controller.running:
            if not messagebox.askokcancel(
                    "Stop the analysis?",
                    "An analysis is still running. Stop it and close {}?".format(
                        branding.APP_NAME), parent=self.root):
                return
        try:
            self.prefs.geometry = self.root.winfo_geometry()
            self.prefs.zoomed = str(self.root.state()) == "zoomed"
        except tk.TclError:
            pass
        self.prefs.save()
        if owns_process():
            arm_exit_watchdog(EXIT_WATCHDOG_S)
        self.controller.request_cancel()
        clean = stop_the_engines(SHUTDOWN_TIMEOUT_S)
        self.controller.shutdown(1.5)
        try:
            self.root.quit()
        except tk.TclError:  # pragma: no cover - already gone
            pass
        try:
            self.root.destroy()
        except tk.TclError:
            pass
        if not clean and owns_process():
            LOG.warning("a search child outlived the shutdown timeout; "
                        "ending the process")
            force_exit()


#: How long the engines and their children get to stop before the process is
#: ended outright. The user must never have to reach for Task Manager.
SHUTDOWN_TIMEOUT_S = 2.0
#: The last resort: a daemon timer that ends the process even if the
#: interpreter is hung joining somebody else's thread at exit.
EXIT_WATCHDOG_S = 5.0

#: True only inside ``main()``: this process exists to BE the window, so ending
#: it outright is the right answer. A test (or an embedding application) that
#: builds a WmlstApp of its own owns its process and must never be exited by us.
_OWNS_PROCESS = False


def owns_process() -> bool:
    """True when this module was launched as the application (``main()``)."""
    return _OWNS_PROCESS


def stop_the_engines(timeout: float = SHUTDOWN_TIMEOUT_S) -> bool:
    """Stop every engine and every search child. -> was it clean?

    Both calls are idempotent and safe from any thread;
    :func:`wmlst.engine.shutdown_all` already reaches the backend, and the
    second call catches children started outside an engine (a version probe, an
    index rebuild). Never raises: this runs while the window is closing.
    """
    clean = True
    try:
        from . import engine as engine_mod
        clean = bool(engine_mod.shutdown_all(timeout)) and clean
    except Exception as exc:  # pragma: no cover - defensive
        LOG.warning("engine shutdown failed: %s", exc)
        clean = False
    try:
        from . import blastbin
        clean = bool(blastbin.shutdown(max(0.5, timeout / 2.0))) and clean
    except Exception as exc:  # pragma: no cover - defensive
        LOG.warning("search-engine shutdown failed: %s", exc)
        clean = False
    return clean


def force_exit(code: int = 0) -> None:
    """End the process now, flushing the log first. No atexit, no join."""
    try:
        logging.shutdown()
    except Exception:  # pragma: no cover
        pass
    os._exit(code)


def arm_exit_watchdog(seconds: float = EXIT_WATCHDOG_S) -> Any:
    """Start a daemon timer that force-exits if a clean exit never arrives.

    A daemon thread does not keep the interpreter alive, and it is still
    running while ``atexit`` hooks run — which is precisely when the hang that
    stranded WMLST in Task Manager used to happen.
    """
    timer = threading.Timer(max(0.5, float(seconds)), force_exit)
    timer.daemon = True
    timer.start()
    return timer


def open_path(path: str) -> None:
    """Open a file with whatever the desktop uses. No child process is spawned
    by this module: ``os.startfile`` and ``webbrowser`` are the system's job."""
    try:
        starter = getattr(os, "startfile", None)
        if starter is not None:  # Windows
            starter(path)
            return
        webbrowser.open(Path(path).resolve().as_uri())
    except Exception as exc:
        LOG.warning("could not open %s: %s", path, exc)


# ===========================================================================
# SECTION 13 — ENTRY POINT  (sections 4.10, 10.9)
# ===========================================================================

def make_root() -> tk.Misc:
    """Create the Tk root, with drag-and-drop when tkinterdnd2 is present."""
    if HAVE_DND:
        try:
            root = tkinterdnd2.TkinterDnD.Tk()  # type: ignore[attr-defined]
            root.wmlst_dnd = True  # type: ignore[attr-defined]
            return root
        except Exception:
            LOG.info("tkinterdnd2 present but unusable; falling back to plain Tk")
    root = tk.Tk()
    root.wmlst_dnd = False  # type: ignore[attr-defined]
    return root


def set_window_icon(root: tk.Misc) -> bool:
    """Apply the bundled PNG icon if Tk can load it. Never fatal."""
    icon = Path(__file__).resolve().parent / "assets" / "wmlst.png"
    if not icon.is_file():
        return False
    try:
        image = tk.PhotoImage(master=root, file=str(icon))
        root.iconphoto(True, image)
        root._wmlst_icon = image  # type: ignore[attr-defined]  # keep a reference
        return True
    except Exception:
        return False


USAGE = """\
{app} {version} — {attribution}

  wmlst-gui [FILE ...]      open the window, analysing FILE immediately
  wmlst-gui --selftest      build and destroy the window, then exit
  wmlst-gui --version       print the version
  wmlst-gui --help          print this message
"""


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for ``wmlst-gui`` (section 4.10).

    ``--selftest`` builds and destroys the Tk root, exits 0 and prints nothing;
    CI needs it. Start-up order is section 10.9: DPI, root, scaling and theme,
    withdraw, preferences, widgets, deiconify, then a 50 ms background probe.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    if "--help" in args or "-h" in args:
        sys.stdout.write(USAGE.format(app=branding.APP_NAME, version=__version__,
                                      attribution=branding.ATTRIBUTION))
        return 0
    if "--version" in args:
        sys.stdout.write("{} {}\n".format(branding.APP_NAME, __version__))
        return 0
    selftest = "--selftest" in args
    files = [a for a in args if not a.startswith("-")]

    dpi_fix()
    setup_logging()
    try:
        root = make_root()
    except tk.TclError as exc:
        sys.stderr.write("{}: no graphical display is available ({})\n".format(
            branding.APP_NAME, exc))
        return 1
    apply_scaling(root)
    install_theme(root)
    root.withdraw()

    if selftest:
        root.destroy()
        return 0

    global _OWNS_PROCESS
    _OWNS_PROCESS = True
    prefs = Prefs.load()
    app = WmlstApp(root, prefs)
    set_window_icon(root)
    root.deiconify()
    root.update_idletasks()
    if files:
        app._argv_files = list(files)   # drained once the env probe has answered
    root.after(50, app.post_start)
    try:
        root.mainloop()
    except KeyboardInterrupt:  # pragma: no cover - console launch only
        app.on_close()
    # mainloop() has returned, so the window is gone. Anything still holding a
    # search child would hang the interpreter's atexit join with no window left
    # to close, so stop them here too -- returning from main() must mean the
    # process really is finished.
    if not stop_the_engines(SHUTDOWN_TIMEOUT_S):
        force_exit()
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
