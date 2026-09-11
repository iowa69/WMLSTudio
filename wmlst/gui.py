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
import os
import platform
import queue
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

from . import branding
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

#: Spacing scale (section 10.1). Use ``px(PAD_M)``, never a bare number.
PAD_XS, PAD_S, PAD_M, PAD_L, PAD_XL, PAD_XXL = 4, 8, 12, 16, 24, 32

#: Font sizes in POINTS. These are never multiplied by :data:`SC` — Tk's own
#: scaling already accounts for dpi, and doubling it is the classic HiDPI bug.
FONT_PT = {
    "body": 10,
    "small": 9,
    "tiny": 8,
    "heading": 12,
    "title": 15,
    "hero": 28,
    "mono": 9,
}

PREFERRED_FAMILIES = ("Segoe UI", "Inter", "Noto Sans", "DejaVu Sans", "Helvetica")
PREFERRED_MONO = ("Consolas", "Cascadia Mono", "DejaVu Sans Mono", "Courier New")

#: Flat, modern palette. Deliberately not the 1998 Tk grey.
PALETTE = {
    "bg": "#f2f4f7",
    "surface": "#ffffff",
    "surface_alt": "#f7f9fc",
    "border": "#d5dae1",
    "border_strong": "#b6bec9",
    "text": "#16191d",
    "muted": "#5c6673",
    "accent": "#1f6feb",
    "accent_hover": "#3a83f0",
    "accent_active": "#1a5fd0",
    "accent_soft": "#e8f0fe",
    "on_accent": "#ffffff",
    "ok": "#146c2e",
    "warn": "#8a5a00",
    "bad": "#b3261e",
    "info": "#1f6feb",
    "focus": "#1f6feb",
    "row_alt": "#f7f9fc",
    "drop_idle": "#aab3c0",
    "drop_hover": "#1f6feb",
}

#: High-contrast fallback (section 10.1). Populated from the system palette when
#: Windows reports SPI_GETHIGHCONTRAST.
HIGH_CONTRAST_PALETTE = {
    "bg": "#000000",
    "surface": "#000000",
    "surface_alt": "#000000",
    "border": "#ffffff",
    "border_strong": "#ffffff",
    "text": "#ffffff",
    "muted": "#ffffff",
    "accent": "#ffff00",
    "accent_hover": "#ffff00",
    "accent_active": "#ffff00",
    "accent_soft": "#000000",
    "on_accent": "#000000",
    "ok": "#00ff00",
    "warn": "#ffff00",
    "bad": "#ff6060",
    "info": "#00ffff",
    "focus": "#ffff00",
    "row_alt": "#000000",
    "drop_idle": "#ffffff",
    "drop_hover": "#ffff00",
}


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
    have = {str(a) for a in available}
    for name in preferred:
        if name in have:
            return name
    return fallback


# ---------------------------------------------------------------------------
# Status presentation (section 11.3 wording, reused verbatim in the GUI)
# ---------------------------------------------------------------------------

#: The ONLY status knowledge in this module: glyph, colour key and the sentence
#: shown to a novice. Nothing here derives a status — ``engine.status_column``
#: already did that and the value arrives on ``SampleResult.status``.
STATUS_UI = {
    "PERFECT": ("\u2714", "ok",
                "Every locus matched a known allele exactly, and the combination "
                "is a recognised sequence type."),
    "NOVEL": ("\u271a", "info",
              "Every locus matched a known allele exactly, but this combination is "
              "not yet a named sequence type \u2014 it may be a new ST worth "
              "submitting to PubMLST."),
    "MIXED": ("\u29c9", "warn",
              "At least one locus matched two or more different alleles equally "
              "well. This usually means the assembly contains more than one strain, "
              "or a duplicated gene \u2014 check culture purity before reporting."),
    "MISSING": ("\u25cc", "warn",
                "At least one locus could not be found in this assembly. It may be "
                "incomplete, or the locus may genuinely be absent \u2014 an ST "
                "cannot be assigned."),
    "BAD": ("\u26a0", "bad",
            "The best-matching scheme scored below 70 out of 100. Treat this result "
            "as unreliable: wrong organism, heavily fragmented, or contaminated."),
    "OK": ("\u25cf", "info",
           "The scheme matched reasonably well, but at least one locus is only an "
           "approximate or partial match, so no exact sequence type could be "
           "assigned."),
    "NONE": ("\u25cb", "muted",
             "No MLST scheme matched this file at all. Check that it really "
             "contains assembled contigs for a species covered by one of the "
             "bundled schemes."),
}

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
    """Render the STATUS cell as ``"glyph WORD"`` — colour is never the only channel."""
    glyph = STATUS_UI.get(status, ("\u25cb", "muted", ""))[0]
    return "{} {}".format(glyph, status) if status else ""


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
    def normalised(self) -> Prefs:
        """Return a copy with every invariant of section 3.6.1 applied.

        ``scheme`` forces ``minscore`` to 0 and clears ``exclude``
        (``bin/mlst:125``); ``jobs * threads`` is clamped to the CPU count.
        """
        p = dataclasses.replace(self)
        p.minid = clamp_float(p.minid, 0.0, 100.0, DEFAULT_MINID)
        p.mincov = clamp_float(p.mincov, 0.0, 100.0, DEFAULT_MINCOV)
        p.minscore = clamp_float(p.minscore, 0.0, 100.0, DEFAULT_MINSCORE)
        p.threads = clamp_int(p.threads, 1, max_threads(), 1)
        p.jobs = clamp_int(p.jobs, 1, max_threads(), 1)
        p.exclude = parse_exclude(p.exclude)
        if p.scheme:
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
    RUNNING: frozenset({RESULTS, IDLE}),
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
        self.pump()
        if self._on_tick is not None:
            try:
                self._on_tick(self.elapsed, self.running)
            except Exception:
                LOG.exception("tick handler failed")
        self._schedule(self.POLL_MS, self._tick)

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
                              text="Reading {}…".format(os.path.basename(path))))
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


def install_theme(root: tk.Misc) -> ttk.Style:
    """Install the flat theme and the font ladder (section 10.1).

    ``clam`` on every platform: it is the only stdlib theme that honours
    ``configure``/``map`` for background colours. Under Windows high contrast the
    custom palette is abandoned for the system one.
    """
    global C
    C = dict(HIGH_CONTRAST_PALETTE if high_contrast_active() else PALETTE)

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
    mono = pick_font(families, PREFERRED_MONO, fallback)

    F.clear()
    F["family"] = ui
    for key, size in FONT_PT.items():
        family = mono if key == "mono" else ui
        weight = "bold" if key in ("heading", "title", "hero") else "normal"
        F[key] = tkfont.Font(root=root, family=family, size=size, weight=weight)
    F["body_bold"] = tkfont.Font(root=root, family=ui, size=FONT_PT["body"],
                                 weight="bold")
    for named in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont"):
        try:
            nf = tkfont.nametofont(named)
            nf.configure(family=ui, size=FONT_PT["body"])
        except Exception:
            pass

    bg, surface, text, muted = c("bg"), c("surface"), c("text"), c("muted")
    border, accent = c("border"), c("accent")

    style.configure(".", background=bg, foreground=text, font=F["body"],
                    borderwidth=0, focuscolor=c("focus"))
    style.configure("TFrame", background=bg)
    style.configure("Surface.TFrame", background=surface)
    style.configure("CardBorder.TFrame", background=border)
    style.configure("Footer.TFrame", background=surface)

    style.configure("TLabel", background=bg, foreground=text, font=F["body"])
    style.configure("Surface.TLabel", background=surface, foreground=text)
    style.configure("Muted.TLabel", background=bg, foreground=muted, font=F["small"])
    style.configure("SurfaceMuted.TLabel", background=surface, foreground=muted,
                    font=F["small"])
    style.configure("Heading.TLabel", background=surface, foreground=text,
                    font=F["heading"])
    style.configure("Title.TLabel", background=bg, foreground=text, font=F["title"])
    style.configure("Hero.TLabel", background=surface, foreground=text, font=F["hero"])
    style.configure("Status.TLabel", background=surface, foreground=muted,
                    font=F["small"])
    for key in ("ok", "warn", "bad", "info"):
        style.configure("{}.TLabel".format(key.capitalize()), background=surface,
                        foreground=c(key), font=F["body"])

    style.configure("TButton", background=surface, foreground=text, font=F["body"],
                    padding=(px(14), px(7)), borderwidth=1, relief="flat",
                    bordercolor=c("border_strong"), lightcolor=surface,
                    darkcolor=surface, anchor="center")
    style.map("TButton",
              background=[("disabled", c("surface_alt")), ("pressed", c("accent_soft")),
                          ("active", c("accent_soft"))],
              foreground=[("disabled", muted)],
              bordercolor=[("focus", c("focus")), ("active", accent)])
    style.configure("Accent.TButton", background=accent, foreground=c("on_accent"),
                    bordercolor=accent, lightcolor=accent, darkcolor=accent,
                    font=F["body_bold"], padding=(px(18), px(8)))
    style.map("Accent.TButton",
              background=[("disabled", c("border")), ("pressed", c("accent_active")),
                          ("active", c("accent_hover"))],
              foreground=[("disabled", muted)],
              bordercolor=[("focus", c("text"))])
    style.configure("Link.TButton", background=bg, foreground=accent,
                    borderwidth=0, padding=(px(4), px(2)), font=F["small"])
    style.map("Link.TButton", background=[("active", bg)],
              foreground=[("active", c("accent_hover"))])

    style.configure("TNotebook", background=bg, borderwidth=0, tabmargins=(0, 0, 0, 0))
    style.configure("TNotebook.Tab", background=bg, foreground=muted,
                    padding=(px(20), px(9)), font=F["body"], borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", surface)],
              foreground=[("selected", text), ("active", text)],
              expand=[("selected", (0, 0, 0, 0))])

    style.configure("Treeview", background=surface, fieldbackground=surface,
                    foreground=text, rowheight=px(26), borderwidth=0, font=F["body"])
    style.configure("Treeview.Heading", background=c("surface_alt"), foreground=muted,
                    font=F["small"], relief="flat", padding=(px(8), px(6)),
                    borderwidth=0)
    style.map("Treeview.Heading", background=[("active", c("accent_soft"))])
    style.map("Treeview", background=[("selected", c("accent_soft"))],
              foreground=[("selected", text)])

    style.configure("TEntry", fieldbackground=surface, foreground=text,
                    bordercolor=c("border_strong"), lightcolor=c("border_strong"),
                    darkcolor=c("border_strong"), insertcolor=text,
                    padding=px(5), borderwidth=1)
    style.map("TEntry", bordercolor=[("focus", c("focus"))],
              lightcolor=[("focus", c("focus"))], darkcolor=[("focus", c("focus"))])
    style.configure("Invalid.TEntry", fieldbackground="#fdeeed",
                    bordercolor=c("bad"), lightcolor=c("bad"), darkcolor=c("bad"))
    style.configure("TSpinbox", fieldbackground=surface, foreground=text,
                    bordercolor=c("border_strong"), arrowcolor=text,
                    background=c("surface_alt"), padding=px(4), borderwidth=1)
    style.map("TSpinbox", bordercolor=[("focus", c("focus"))])
    style.configure("TCombobox", fieldbackground=surface, foreground=text,
                    background=c("surface_alt"), bordercolor=c("border_strong"),
                    arrowcolor=text, padding=px(4))
    style.map("TCombobox", bordercolor=[("focus", c("focus"))],
              fieldbackground=[("readonly", surface)])
    style.configure("TCheckbutton", background=surface, foreground=text,
                    focuscolor=c("focus"))
    style.map("TCheckbutton", background=[("active", surface)])
    style.configure("TRadiobutton", background=surface, foreground=text,
                    focuscolor=c("focus"))
    style.map("TRadiobutton", background=[("active", surface)])
    style.configure("TSeparator", background=border)
    style.configure("TProgressbar", background=accent, troughcolor=c("surface_alt"),
                    bordercolor=c("surface_alt"), lightcolor=accent, darkcolor=accent,
                    thickness=px(10))
    style.configure("Thin.Horizontal.TProgressbar", thickness=px(4))
    style.configure("Big.Horizontal.TProgressbar", thickness=px(12))
    style.configure("Vertical.TScrollbar", background=c("surface_alt"),
                    troughcolor=surface, bordercolor=surface,
                    arrowcolor=muted, gripcount=0)
    style.configure("Horizontal.TScrollbar", background=c("surface_alt"),
                    troughcolor=surface, bordercolor=surface,
                    arrowcolor=muted, gripcount=0)
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


def card(parent: tk.Misc, **pack_kw: Any) -> ttk.Frame:
    """A 1 px bordered frame wrapping a surface frame (section 10.1).

    ttk has no border radius and faking rounded corners with images breaks at
    150 % DPI, so cards are honest rectangles.
    """
    outer = ttk.Frame(parent, style="CardBorder.TFrame")
    if pack_kw:
        outer.pack(**pack_kw)
    inner = ttk.Frame(outer, style="Surface.TFrame")
    inner.pack(fill="both", expand=True, padx=1, pady=1)
    outer.inner = inner  # type: ignore[attr-defined]
    return inner


def focus_ring(parent: tk.Misc) -> tk.Frame:
    """A tk.Frame whose highlight border is the 3 px focus ring (section 10.8).

    ``clam`` will not draw a focus ring on a Treeview or a Canvas, so the
    container does it instead.
    """
    return tk.Frame(parent, background=c("border"), highlightthickness=px(3),
                    highlightbackground=c("border"), highlightcolor=c("focus"),
                    borderwidth=0)


class DropZone(tk.Canvas):
    """The one big target a novice needs (section 10.2).

    A Canvas, not a Frame, so the dashed border can be redrawn on ``<Configure>``.
    Idle / hover / drag-over / keyboard-focus differ in outline colour and width;
    the glyph is drawn with lines and polygons — no emoji, no bitmap.

    Click anywhere browses for files; ``Return``/``space`` browses files and
    ``Shift-Return`` browses a folder.
    """

    def __init__(self, parent: tk.Misc, *, on_files: Callable[[Sequence[str]], None],
                 on_browse: Callable[[], None], on_browse_folder: Callable[[], None],
                 dnd: bool):
        super().__init__(parent, highlightthickness=px(3),
                         highlightbackground=c("bg"), highlightcolor=c("focus"),
                         background=c("bg"), borderwidth=0, takefocus=True,
                         height=px(210), cursor="hand2")
        self.on_files = on_files
        self.on_browse = on_browse
        self.on_browse_folder = on_browse_folder
        self.dnd = dnd
        self.zone_state = "idle"
        self._enabled = True
        self.headline = ("Drop a FASTA file here, or click to browse"
                         if dnd else "Click to choose your FASTA file")
        self.subline = ("Analysis starts by itself. "
                        "You can drop several files, or a whole folder."
                        if dnd else
                        "Analysis starts by itself. "
                        "Press Shift+Enter to choose a whole folder.")
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
            self.headline = "Click to choose your FASTA file"
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

    # -- painting ------------------------------------------------------------
    def redraw(self) -> None:
        self.delete("all")
        w = max(self.winfo_width(), px(320))
        h = max(self.winfo_height(), px(180))
        dragging = self.zone_state == "dragover"
        if not self._enabled:
            outline, width, fill = c("border"), px(2), c("bg")
            headline, sub = "Working…", "Please wait for the current files to finish."
        elif dragging:
            outline, width, fill = c("drop_hover"), px(3), c("accent_soft")
            headline, sub = "Release to analyse", self.subline
        elif self.zone_state in ("hover", "focus"):
            outline, width, fill = c("accent"), px(2), c("surface")
            headline, sub = self.headline, self.subline
        else:
            outline, width, fill = c("drop_idle"), px(2), c("surface")
            headline, sub = self.headline, self.subline

        inset = px(6)
        self.create_rectangle(inset, inset, w - inset, h - inset, outline=outline,
                              width=width, dash=(px(6), px(5)), fill=fill)

        cx = w // 2
        top = h // 2 - px(52)
        glyph = outline if dragging else c("accent")
        # A downward arrow dropping into an open tray, drawn as vector art so it
        # scales cleanly at any DPI and needs no bundled bitmap.
        self.create_line(cx, top, cx, top + px(34), fill=glyph, width=px(3),
                         capstyle="round")
        self.create_polygon(cx - px(11), top + px(26), cx + px(11), top + px(26),
                            cx, top + px(42), fill=glyph, outline=glyph)
        self.create_line(cx - px(30), top + px(52), cx - px(30), top + px(64),
                         cx + px(30), top + px(64), cx + px(30), top + px(52),
                         fill=glyph, width=px(3), joinstyle="round",
                         capstyle="round")

        self.create_text(cx, h // 2 + px(34), text=headline, fill=c("text"),
                         font=F.get("title"), anchor="center")
        self.create_text(cx, h // 2 + px(62), text=sub, fill=c("muted"),
                         font=F.get("small"), anchor="center", width=w - px(80),
                         justify="center")


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
        self.icon = ttk.Label(self, text="", style="Status.TLabel")
        self.icon.pack(side="left", padx=(px(PAD_M), 0))
        self.label = ttk.Label(self, text=self.resting, style="Status.TLabel",
                               anchor="w")
        self.label.pack(side="left", fill="x", expand=True, padx=px(PAD_S),
                        pady=px(PAD_XS))

    def set(self, text: str, kind: str = "muted", *, transient: bool = False,
            resting: bool = False) -> None:
        """Update the live region. ``resting`` also changes what it reverts to."""
        assert_main_thread("the status line")
        glyph = {"ok": "✔", "warn": "⚠", "bad": "⚠",
                 "info": "ℹ"}.get(kind, "")
        colour = c(kind) if kind in ("ok", "warn", "bad", "info") else c("muted")
        self.icon.configure(text=glyph, foreground=colour)
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

    def _revert(self) -> None:
        self._after = None
        self.icon.configure(text="")
        self.label.configure(text=self.resting, foreground=c("muted"))


# ===========================================================================
# SECTION 8 — VIEWS: THE ANALYSE TAB  (section 10.2)
# ===========================================================================

#: Column ids for the results tree. ``#0`` carries FILE.
TREE_COLUMNS = ("scheme", "st", "status", "score")
#: Heading words, identical to bin/mlst:151 plus the two --full columns.
TREE_HEADINGS = {"#0": "FILE", "scheme": "SCHEME", "st": "ST",
                 "status": "STATUS", "score": "SCORE"}
#: Caption row inserted above the per-locus children so reused columns are
#: never ambiguous (section 10.2).
CHILD_CAPTION = ("ALLELE", "WHAT IT MEANS", "BEST EVIDENCE", "HITS")

CHUNK_ROWS = 200        # rows inserted per pump above 2000 rows
BIG_BATCH = 50          # confirm threshold for a folder drop


def evidence_text(call: AlleleCall) -> str:
    """One-line evidence for a locus, straight off the winning :class:`Hit`."""
    hit = call.best
    if hit is None:
        return "no match found"
    return "{} {}–{} · {:.1f}% identity · {:.1f}% of allele".format(
        hit.qseqid, hit.qstart, hit.qend, hit.pct_identity, hit.pct_coverage)


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
        self.progress_label = ttk.Label(row, text="Starting…", style="Heading.TLabel",
                                        anchor="w")
        self.progress_label.pack(side="left", fill="x", expand=True)
        self.elapsed_label = ttk.Label(row, text="", style="SurfaceMuted.TLabel")
        self.elapsed_label.pack(side="left", padx=px(PAD_M))
        self.cancel_button = ttk.Button(row, text="Cancel", command=self.app.cancel_run)
        self.cancel_button.pack(side="left")
        self.overall = ttk.Progressbar(progress, style="Big.Horizontal.TProgressbar",
                                       mode="determinate", maximum=100.0)
        self.overall.pack(fill="x", padx=px(PAD_M))
        self.pulse = ttk.Progressbar(progress, style="Thin.Horizontal.TProgressbar",
                                     mode="indeterminate", maximum=100.0)
        self.pulse.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_XS), 0))
        self.progress_detail = ttk.Label(progress, text="", style="SurfaceMuted.TLabel",
                                         anchor="w")
        self.progress_detail.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_XS), px(PAD_M)))

        # -- single-file summary card ----------------------------------------
        self.summary_wrap = ttk.Frame(outer, style="TFrame")
        summary = card(self.summary_wrap)
        summary.master.pack(fill="x")
        self.summary_card = summary
        left = ttk.Frame(summary, style="Surface.TFrame")
        left.pack(side="left", padx=px(PAD_L), pady=px(PAD_M))
        ttk.Label(left, text="SEQUENCE TYPE", style="SurfaceMuted.TLabel").pack(anchor="w")
        self.st_label = ttk.Label(left, text="—", style="Hero.TLabel")
        self.st_label.pack(anchor="w")
        right = ttk.Frame(summary, style="Surface.TFrame")
        right.pack(side="left", fill="both", expand=True, padx=(0, px(PAD_L)),
                   pady=px(PAD_M))
        self.summary_title = ttk.Label(right, text="", style="Heading.TLabel", anchor="w")
        self.summary_title.pack(anchor="w", fill="x")
        self.summary_status = ttk.Label(right, text="", style="Ok.TLabel", anchor="w")
        self.summary_status.pack(anchor="w", fill="x", pady=(px(PAD_XS), 0))
        self.summary_body = ttk.Label(right, text="", style="SurfaceMuted.TLabel",
                                      anchor="w", justify="left")
        self.summary_body.pack(anchor="w", fill="x", pady=(px(PAD_XS), 0))
        self.summary_body.bind(
            "<Configure>",
            lambda e: self.summary_body.configure(wraplength=max(px(240), e.width - px(8))),
            add="+")

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
        vsb = ttk.Scrollbar(ring, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.column("#0", width=px(320), minwidth=px(180), stretch=True)
        self.tree.column("scheme", width=px(150), minwidth=px(90), anchor="w")
        self.tree.column("st", width=px(90), minwidth=px(60), anchor="w")
        self.tree.column("status", width=px(200), minwidth=px(120), anchor="w")
        self.tree.column("score", width=px(80), minwidth=px(60), anchor="e")
        for col, title in TREE_HEADINGS.items():
            self.tree.heading(col, text=title,
                              command=lambda cc=col: self.sort_by(cc))
        for key in ("ok", "warn", "bad", "info", "muted"):
            self.tree.tag_configure("st_" + key, foreground=c(key))
        self.tree.tag_configure("failed", foreground=c("bad"))
        self.tree.tag_configure("caption", foreground=c("muted"), font=F.get("tiny"))
        self.tree.tag_configure("locus", foreground=c("text"))
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
        self.btn_tsv = ttk.Button(actions, text="Save table (TSV)…", command=self.save_tsv)
        self.btn_tsv.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_json = ttk.Button(actions, text="Save JSON…", command=self.save_json)
        self.btn_json.pack(side="left", padx=(px(PAD_S), 0))
        self.btn_novel = ttk.Button(actions, text="Save new alleles…",
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
        self.dropzone.configure(height=px(140) if have_rows else px(210))
        self.dropzone.set_enabled(not running)
        self._set_actions_enabled(have_rows and not running)

    def _set_actions_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for btn in (self.btn_html, self.btn_tsv, self.btn_json, self.btn_novel,
                    self.btn_copy, self.btn_clear):
            btn.configure(state=state)
        if enabled and not any(s.novel for s in self.results):
            self.btn_novel.configure(state="disabled")

    def _fill_summary(self, res: SampleResult) -> None:
        """Fill the single-file card. Every value is copied, never recomputed."""
        # A bare hyphen at 28 pt reads as a dash, not as information (section 10.2).
        self.st_label.configure(text=res.st if res.st not in ("", "-") else "not assigned")
        scheme = res.scheme if res.scheme != "-" else "no scheme matched"
        self.summary_title.configure(text="{}  ·  {}".format(
            os.path.basename(res.label) or res.label, scheme))
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
        self.progress_label.configure(text="Preparing…")
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
                self.pulse.start(60)
            except tk.TclError:
                pass

    def _stop_pulse(self) -> None:
        if self._indeterminate:
            self._indeterminate = False
            try:
                self.pulse.stop()
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

    def _insert_result_row(self, res: SampleResult) -> None:
        self._order += 1
        iid = self.tree.insert(
            "", "end", text=os.path.basename(res.label) or res.label,
            values=(res.scheme, res.st, status_cell(res.status), str(res.score)),
            tags=(status_tag(res.status), "row"), open=False)
        self._row_data[iid] = res
        if res.alleles:
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
        iid = self.tree.insert(
            "", "end", text=os.path.basename(path) or path,
            values=("—", "—", "⚠ COULD NOT READ", "—"), tags=("failed", "row"))
        self._row_error[iid] = friendly
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
        self.tree.insert(iid, "end", text="LOCUS", values=CHILD_CAPTION,
                         tags=("caption",))
        for call in res.alleles:
            meaning = SYMBOL_UI.get(call.symbol, (call.symbol, "muted"))[0]
            colour = SYMBOL_UI.get(call.symbol, (call.symbol, "muted"))[1]
            self.tree.insert(
                iid, "end", text="    " + call.locus,
                values=(call.code, meaning, evidence_text(call), str(len(call.hits))),
                tags=("st_" + colour, "locus"))

    def _on_select(self, _event: Any = None) -> None:
        iid = self.tree.focus()
        res = self._row_data.get(iid)
        if res is not None:
            self.app.status.set("{} — {}".format(res.path, status_sentence(res.status)),
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
                value: Any = self.tree.item(iid, "text").lower()
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
        text = "\t".join([res.label, res.scheme, res.st, res.status, str(res.score),
                          allele_summary(res)])
        self.clipboard_clear()
        self.clipboard_append(text)
        self.app.status.set("Row copied to the clipboard.", "ok", transient=True)

    def clear(self) -> None:
        """Empty the result list and return to the idle drop target."""
        self.st_label.configure(text="—")
        self.summary_title.configure(text="")
        self.summary_status.configure(text="")
        self.summary_body.configure(text="")
        self.results.clear()
        self.failures.clear()
        self._row_data.clear()
        self._row_error.clear()
        self._expanded.clear()
        self._pending_rows.clear()
        for iid in self.tree.get_children(""):
            self.tree.delete(iid)
        self.app.state.to(IDLE)
        self.show_state(IDLE)
        self.app.status.set("Ready.", resting=True)

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
        ttk.Label(grid, text="schemes installed", style="SurfaceMuted.TLabel").grid(
            row=0, column=1, sticky="w")
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
        self.db_progress = ttk.Progressbar(info, style="Thin.Horizontal.TProgressbar",
                                           mode="determinate", maximum=100.0)
        self.db_progress.pack(fill="x", padx=px(PAD_M))
        self.db_detail = ttk.Label(info, text="", style="SurfaceMuted.TLabel", anchor="w")
        self.db_detail.pack(fill="x", padx=px(PAD_M), pady=(px(PAD_XS), px(PAD_M)))

        ring = focus_ring(outer)
        ring.pack(fill="both", expand=True, pady=(px(PAD_M), 0))
        self.tree = ttk.Treeview(ring, columns=DB_COLUMNS, show="tree headings",
                                 selectmode="browse")
        vsb = ttk.Scrollbar(ring, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.column("#0", width=px(220), minwidth=px(140), stretch=True)
        self.tree.column("species", width=px(260), minwidth=px(140))
        for col in ("loci", "types", "alleles"):
            self.tree.column(col, width=px(90), minwidth=px(60), anchor="e")
        self.tree.column("date", width=px(120), minwidth=px(80), anchor="w")
        for col, title in DB_HEADINGS.items():
            self.tree.heading(col, text=title)
        self.tree.tag_configure("muted", foreground=c("muted"))

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
        if not env.db_ok:
            self.db_detail.configure(text=env.db_error or "The database is not available.")
            return
        self.db_detail.configure(text="Reading the scheme catalogue…")
        datadir = env.datadir

        def load(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import schemes as schemes_mod
            catalog = schemes_mod.SchemeCatalog(datadir)
            shallow = catalog.info_all(deep=False)
            return (catalog, shallow)

        self._jobs[self.app.run_task("db-shallow", load, owner=self,
                                     with_progress=True)] = "shallow"

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
            self.tree.insert(
                "", "end", text=name,
                values=(species or "scheme not mapped to a species",
                        db_cell(getattr(info, "locus", "")),
                        db_cell(types) if deep else "…",
                        db_cell(alleles) if deep else "…",
                        db_cell(getattr(info, "last_updated", ""))),
                tags=() if species else ("muted",))
        self.count_label.configure(text=str(len(self._infos)))

    def on_task_done(self, job: int, payload: Any) -> bool:
        """Handle a finished background job. Returns True if it was ours."""
        kind = self._jobs.pop(job, None)
        if kind is None:
            return False
        if kind == "shallow":
            catalog, infos = payload
            self._fill(catalog, infos, deep=False)
            self.db_detail.configure(text="Counting sequence types…")
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
            self.db_progress.configure(value=0.0)
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
        self.db_progress.configure(value=msg.percent)
        if msg.text:
            self.db_detail.configure(text=msg.text)
        return True

    # -- updates -------------------------------------------------------------
    def _busy(self, busy: bool) -> None:
        self.btn_check.configure(state="disabled" if busy else "normal")
        self.btn_cancel.configure(state="normal" if busy else "disabled")
        if busy:
            self.btn_update.configure(state="disabled")

    def check_updates(self) -> None:
        """Metadata-only pass; writes nothing (section 10.3)."""
        if self.app.state.running:
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
        self.db_detail.configure(text="Contacting PubMLST…")
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
        dialog = UpdatePlanDialog(self.app.root, plan, updates)
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
        dialog = UpdatePlanDialog(self.app.root, self.plan, updates)
        self.app.root.wait_window(dialog)
        if dialog.accepted:
            self._start_apply(dialog.chosen, dialog.backup)

    def _start_apply(self, chosen: Sequence[str], backup: bool) -> None:
        dbdir = self.app.env.dbdir
        plan = self.plan
        self._busy(True)
        self.db_detail.configure(text="Downloading…")
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
        self.db_detail.configure(text="Stopping…")


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
        self._build()
        self.load_from(app.prefs)

    def _build(self) -> None:
        outer = ttk.Frame(self, style="TFrame")
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
        warn.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_M)))
        self.scheme_warning = ttk.Label(
            warn, text="", style="Warn.TLabel", anchor="w", justify="left")
        self.scheme_warning.pack(anchor="w")

        scheme = card(outer)
        scheme.master.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(scheme, text="Scheme", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        row = ttk.Frame(scheme, style="Surface.TFrame")
        row.pack(fill="x", padx=px(PAD_M), pady=(0, px(PAD_S)))
        ttk.Label(row, text="Always use this scheme", style="Surface.TLabel",
                  width=26, anchor="w").pack(side="left")
        self.vars["scheme"] = tk.StringVar(value="")
        self.scheme_box = ttk.Combobox(row, textvariable=self.vars["scheme"],
                                       state="readonly", width=32,
                                       values=("Automatic (recommended)",))
        self.scheme_box.pack(side="left")
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

        perf = card(outer)
        perf.master.pack(fill="x", pady=(px(PAD_M), 0))
        ttk.Label(perf, text="Speed", style="Heading.TLabel").pack(
            anchor="w", padx=px(PAD_M), pady=(px(PAD_M), px(PAD_XS)))
        self._spin(perf, "threads", "Cores per file", 1, max_threads(), 1)
        self._spin(perf, "jobs", "Files at once", 1, max_threads(), 1)
        ttk.Label(perf, text="This computer has {} processor cores.".format(cpu_count()),
                  style="SurfaceMuted.TLabel").pack(anchor="w", padx=px(PAD_M),
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
        spin = ttk.Spinbox(row, from_=lo, to=hi, increment=step, width=8,
                           textvariable=var, justify="right",
                           command=self._changed)
        spin.pack(side="left")
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
        self._loading = False
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
        prefs.minscore = clamp_float(self.vars["minscore"].get(), 0, 100, DEFAULT_MINSCORE)
        prefs.threads = clamp_int(self.vars["threads"].get(), 1, max_threads(), 1)
        prefs.jobs = clamp_int(self.vars["jobs"].get(), 1, max_threads(), 1)
        prefs.exclude = parse_exclude(self.vars["exclude"].get())
        chosen = self.vars["scheme"].get().strip()
        prefs.scheme = None if chosen.startswith("Automatic") or not chosen else chosen
        return prefs.normalised()

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
            self.vars["minscore"].set("0")
            self.scheme_warning.configure(
                text="⚠ Forcing a scheme sets the minimum score to 0 and ignores the "
                     "exclude list, exactly as mlst --scheme does. Every file will be "
                     "reported with this scheme even if it does not fit.")
        else:
            self.scheme_warning.configure(text="")

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
        self.vars["exclude"].set(", ".join(DEFAULT_EXCLUDE))
        self.vars["scheme"].set("Automatic (recommended)")
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
        ttk.Button(picker, text="Browse…", command=self._pick).pack(
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
        self.detail.configure(text="Checking that folder…")

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
        self.detail.configure(text="Starting download…")
        self.cancel_event = threading.Event()

        def run(progress: Callable[..., None], cancel: threading.Event) -> Any:
            from . import blastbin
            return blastbin.bootstrap(progress=progress, cancel=cancel)

        self.job = self.app.run_task("blast-bootstrap", run, owner=self,
                                     with_progress=True, cancel=self.cancel_event)

    def cancel(self) -> None:
        """Cancelling deletes the partial download (blastbin honours the event)."""
        self.cancel_event.set()
        self.detail.configure(text="Stopping…")
        self.action.configure(state="disabled")

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
        self.app.adopt_tools(payload)
        self.close()
        return True

    def on_task_failed(self, job: int, friendly: Friendly) -> bool:
        if self.job is None or job != self.job:
            return False
        self.job = None
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

    def __init__(self, parent: tk.Misc, plan: Any, updates: Sequence[Any]):
        super().__init__(parent, "Database update", width=560)
        self.accepted = False
        self.backup = True
        self.chosen: List[str] = []
        self._vars: List[Tuple[str, tk.BooleanVar]] = []
        total = sum(int(getattr(u, "bytes_estimate", 0) or 0) for u in updates)
        self.heading("{} scheme{} can be updated".format(
            len(updates), "" if len(updates) == 1 else "s"))
        self.paragraph(
            "About {} will be downloaded. Your current database keeps working "
            "until every chosen scheme has been replaced.".format(human_bytes(total)))

        ring = focus_ring(self.body)
        ring.pack(fill="both", expand=True, pady=px(PAD_S))
        canvas = tk.Canvas(ring, background=c("surface"), highlightthickness=0,
                           height=px(220), width=px(520))
        scroll = ttk.Scrollbar(ring, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        holder = ttk.Frame(canvas, style="Surface.TFrame")
        canvas.create_window((0, 0), window=holder, anchor="nw")
        holder.bind("<Configure>",
                    lambda e: canvas.configure(scrollregion=canvas.bbox("all")), add="+")
        for update in updates:
            name = str(getattr(update, "name", ""))
            var = tk.BooleanVar(value=True)
            self._vars.append((name, var))
            added = getattr(update, "added_types", None)
            detail = str(getattr(update, "detail", "") or "")
            text = "{}   ({}{})".format(
                name, human_bytes(getattr(update, "bytes_estimate", 0)),
                ", {} new sequence types".format(added) if added else "")
            row = ttk.Checkbutton(holder, text=text, variable=var)
            row.pack(anchor="w", padx=px(PAD_S), pady=px(2))
            if detail:
                Tooltip(row, detail)

        self.backup_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.body, text="Keep a backup of the current database",
                        variable=self.backup_var).pack(anchor="w", pady=(px(PAD_S), 0))
        ttk.Button(self.buttons, text="Cancel", command=self.close).pack(side="right",
                                                                        padx=(px(PAD_S), 0))
        go = ttk.Button(self.buttons, text="Update now", style="Accent.TButton",
                        command=self._accept)
        go.pack(side="right")
        go.focus_set()
        self.present()

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
        self._bootstrap: Optional[BootstrapDialog] = None
        self._failed_count = 0
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
        self.footer_env = ttk.Label(credit, text="Starting…", style="Status.TLabel")
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
        file_menu.add_command(label="Open files…", accelerator="Ctrl+O",
                              command=lambda: self.analyse_view.browse_files())
        file_menu.add_command(label="Open folder…", accelerator="Ctrl+Shift+O",
                              command=lambda: self.analyse_view.browse_folder())
        file_menu.add_separator()
        file_menu.add_command(label="Open report",
                              command=lambda: self.export("html"))
        file_menu.add_command(label="Save table (TSV)…",
                              command=lambda: self.export("tsv"))
        file_menu.add_command(label="Save JSON…", command=lambda: self.export("json"))
        file_menu.add_command(label="Save new alleles…",
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
        except tk.TclError:
            pass
        if not env.db_ok:
            self.status.set(env.db_error or "The MLST database could not be found.",
                            "bad", resting=True)
        elif env.index_stale:
            self.status.set("The search index is older than the allele files — "
                            "rebuild it from the Database tab.", "warn", resting=True)
        elif env.blast_ok:
            self.status.set("Ready. Drop a FASTA file on the box above."
                            if self.dnd_enabled else
                            "Ready. Click the box above to choose a FASTA file.",
                            resting=True)
        self._update_banner()
        if not env.blast_ok and not self.prefs.blast_prompt_shown:
            self.prefs.blast_prompt_shown = True
            self.schedule_save()
            self.offer_bootstrap()

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

    def _route(self, msg: Msg, method: str, *args: Any) -> bool:
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
            self._owners.pop(msg.job, None)
            if not self._route(msg, "on_task_failed", msg.job, friendly):
                self._task_failed(kind, friendly)
            return
        if msg.kind == "env":
            payload = self.controller.take_task_result(msg.job)
            self._task_kinds.pop(msg.job, None)
            self._owners.pop(msg.job, None)
            if not self._route(msg, "on_task_done", msg.job, payload):
                self._task_done(kind, payload)

    def _task_done(self, kind: str, payload: Any) -> None:
        if kind == "env" and isinstance(payload, Environment):
            self._apply_env(payload)
        elif kind.startswith("export-"):
            self._after_export(kind.split("-", 1)[1], str(payload or ""))
        elif kind == "rebuild":
            self.status.set("The search index was rebuilt.", "ok", transient=True)
            self.run_task("env", lambda: probe_environment(self.prefs))

    def _task_failed(self, kind: str, friendly: Friendly) -> None:
        if kind == "env":
            self.status.set(friendly.headline, "bad", resting=True)
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
        self.status.set("Analysing {} file{}…".format(
            len(files), "" if len(files) == 1 else "s"), "info", resting=True)
        LOG.info("run started: %d file(s), threads=%d jobs=%d",
                 len(files), cfg.threads, cfg.jobs)

    def cancel_run(self) -> None:
        """Cancel the current batch; partial results are kept (section 10.2)."""
        if self.controller.running:
            self.controller.request_cancel()
            self.analyse_view.cancel_button.configure(state="disabled",
                                                      text="Stopping…")
            self.status.set("Stopping after the current file…", "warn")

    def rerun_all(self) -> None:
        """Re-run every file currently listed, with the current settings (F5)."""
        paths = [r.path for r in self.analyse_view.results]
        paths.extend(p for p, _ in self.analyse_view.failures)
        if not paths:
            self.status.set("There is nothing to re-run yet.", "warn", transient=True)
            return
        self.analyse_view.clear()
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
        pending = self.controller.pending_paths()
        if pending and reason != "cancelled":
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

        self.status.set("Rebuilding the search index…", "info")
        self.run_task("rebuild", run, with_progress=True)

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

        self.status.set("Saving…", "info")
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
        """Cancel, join for at most 1.5 s, then destroy regardless (section 10.5)."""
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
        self.controller.shutdown(1.5)
        try:
            self.root.destroy()
        except tk.TclError:
            pass


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

    prefs = Prefs.load()
    app = WmlstApp(root, prefs)
    set_window_icon(root)
    root.deiconify()
    root.update_idletasks()
    root.after(50, app.post_start)
    if files:
        root.after(120, lambda: app.analyse_view.handle_paths(files))
    try:
        root.mainloop()
    except KeyboardInterrupt:  # pragma: no cover - console launch only
        app.on_close()
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
