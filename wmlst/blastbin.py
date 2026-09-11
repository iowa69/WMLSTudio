# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""BLAST+ discovery, bootstrap, process-launch policy and Windows path safety.

This is the ONLY module in WMLST permitted to import :mod:`subprocess`
(docs/ARCHITECTURE.md sections 2.1, 4.4, 8, 9).

Everything here is written so that it imports and behaves sanely on Linux — every
Windows-only mechanism (``CREATE_NO_WINDOW``, ``GetShortPathNameW``,
``SetErrorMode``, the ``.exe`` suffix, the bootstrap downloader) is guarded on
``os.name == 'nt'`` and degrades to the obvious POSIX behaviour.
"""

from __future__ import annotations

import contextlib
import hashlib
import itertools
import logging
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import urllib.error
import urllib.request
import uuid
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

__all__ = [
    "BLAST_BYTES",
    "BLAST_MD5",
    "BLAST_MEMBERS",
    "BLAST_MIN_VERSION",
    "BLAST_OUTFMT",
    "BLAST_PINNED_VERSION",
    "BLAST_URL",
    "CREATE_NO_WINDOW",
    "IS_WINDOWS",
    "BlastTools",
    "ansi_safe",
    "blastn_argv",
    "bootstrap",
    "child_env",
    "find_blast",
    "install_root",
    "job_dir",
    "makeblastdb_argv",
    "probe_version",
    "run_blastn",
    "run_makeblastdb",
    "run_tool",
    "temp_root",
    "terminate_all",
    "verify_sha256sums",
    "write_sha256sums",
]

_LOG = logging.getLogger("wmlst.blastbin")

IS_WINDOWS = os.name == "nt"

#: Windows process-creation flag that keeps a console from flashing (section 8.2).
CREATE_NO_WINDOW = 0x08000000
_STARTF_USESHOWWINDOW = 1
_SW_HIDE = 0

#: NTSTATUS codes that mean "the loader could not start the child" (section 8.3).
_STATUS_DLL_NOT_FOUND = 0xC0000135
_STATUS_DLL_INIT_FAILED = 0xC0000142
_VC_REDIST_URL = "https://aka.ms/vs/17/release/vc_redist.x64.exe"

BLAST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/"
    "ncbi-blast-2.17.0+-x64-win64.tar.gz"
)
BLAST_MD5 = "dcd973097407a2910061ff4fb51b09fb"
BLAST_BYTES = 143_400_333
BLAST_MEMBERS: Tuple[str, ...] = (
    "blastn.exe",
    "makeblastdb.exe",
    "nghttp2.dll",
    "blastdbcmd.exe",
    "blastn.exe.manifest",
    "makeblastdb.exe.manifest",
    "blastdbcmd.exe.manifest",
    "LICENSE",
    "BLAST_PRIVACY",
)

#: The build WMLST is pinned to; anything else is usable but warned about.
BLAST_PINNED_VERSION = "2.17.0+"
#: Below this, seeding and HSP boundaries are not the ones the goldens were made with.
BLAST_MIN_VERSION: Tuple[int, int, int] = (2, 9, 0)

#: The nine columns the engine parses (section 5.5). One string, one argv element.
BLAST_OUTFMT = "6 sseqid slen length nident qseqid qstart qend qseq sstrand"

#: Sentinel written LAST by bootstrap(), so a half-extracted tree is never trusted.
INSTALL_OK = ".wmlst-install-ok"
SHA256SUMS = "SHA256SUMS"

_EXC_CACHE: Dict[str, type] = {}
_PROCS: set = set()
_PROCS_LOCK = threading.Lock()
_JOB_COUNTER = itertools.count(1)


def _exc(name: str) -> type:
    """Resolve an exception class from ``wmlst.engine`` lazily.

    ``engine`` imports this module (section 2.1), so a module-level
    ``from wmlst.engine import ...`` would be a circular import.
    """
    cls = _EXC_CACHE.get(name)
    if cls is None:
        from wmlst import engine  # local import: see docstring

        cls = getattr(engine, name)
        _EXC_CACHE[name] = cls
    return cls


def __getattr__(name: str):  # PEP 562 — re-export the engine exceptions lazily
    if name in (
        "WmlstError",
        "Cancelled",
        "BlastNotFoundError",
        "BlastFailedError",
        "BootstrapError",
        "OutOfDiskError",
    ):
        return _exc(name)
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


@dataclass(frozen=True)
class BlastTools:
    """A resolved BLAST+ installation (section 4.4)."""

    blastn: str
    makeblastdb: str
    blastdbcmd: Optional[str]
    version: str
    version_tuple: Tuple[int, ...]
    origin: str


# ---------------------------------------------------------------------------
# Parent-process hardening (section 8.2)
# ---------------------------------------------------------------------------

def _harden_parent() -> None:
    """``SetErrorMode`` so a child loader failure cannot pop a modal dialog."""
    if not IS_WINDOWS:
        return
    try:
        import ctypes

        SEM_FAILCRITICALERRORS = 0x0001
        SEM_NOGPFAULTERRORBOX = 0x0002
        ctypes.windll.kernel32.SetErrorMode(  # type: ignore[attr-defined]
            SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX
        )
    except Exception:  # pragma: no cover - never fatal
        _LOG.debug("SetErrorMode failed", exc_info=True)


_harden_parent()


# ---------------------------------------------------------------------------
# 4.4  run_tool — the one process launcher
# ---------------------------------------------------------------------------

def _startupinfo():
    if not IS_WINDOWS:
        return None
    si = subprocess.STARTUPINFO()  # type: ignore[attr-defined]
    si.dwFlags |= _STARTF_USESHOWWINDOW
    si.wShowWindow = _SW_HIDE
    return si


def run_tool(
    argv: Sequence[str],
    *,
    cwd: Optional[str] = None,
    env: Optional[Dict[str, str]] = None,
    timeout: Optional[float] = None,
    cancel: Optional[threading.Event] = None,
) -> subprocess.CompletedProcess:
    """Launch a child process under the mandatory policy of section 4.4/8.2.

    ``args`` is always a LIST and ``shell=False``; stdin is ``DEVNULL`` (a windowed
    parent has no valid std handles); stdout/stderr are pipes decoded as UTF-8 with
    ``errors='replace'`` and universal newlines.  On Windows the child gets
    ``CREATE_NO_WINDOW`` plus ``STARTUPINFO``/``SW_HIDE`` and nothing else — never
    ``DETACHED_PROCESS`` or ``CREATE_NEW_CONSOLE``.

    Raises ``subprocess.TimeoutExpired`` after ``kill()`` **then** ``wait()`` (the
    wait is required on Windows or the handle leaks and the job dir cannot be
    removed), and ``engine.Cancelled`` when `cancel` is set.
    """
    args = [str(a) for a in argv]
    kwargs: Dict[str, object] = {
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "close_fds": True,
        "cwd": cwd,
        "env": env,
        "shell": False,
        "encoding": "utf-8",
        "errors": "replace",
    }
    if IS_WINDOWS:
        kwargs["creationflags"] = CREATE_NO_WINDOW
        kwargs["startupinfo"] = _startupinfo()

    _LOG.debug("run_tool: %r (cwd=%r)", args, cwd)
    try:
        proc = subprocess.Popen(args, **kwargs)  # type: ignore[arg-type]
    except OSError as exc:
        raise _exc("BlastNotFoundError")(
            "Could not start %s: %s" % (args[0], exc),
            user_message="WMLST could not start '%s'." % (os.path.basename(args[0]),),
        ) from exc

    with _PROCS_LOCK:
        _PROCS.add(proc)
    deadline = None if timeout is None else time.monotonic() + timeout
    try:
        while True:
            if cancel is not None and cancel.is_set():
                _terminate(proc)
                raise _exc("Cancelled")("Cancelled before %s finished" % (args[0],))
            if deadline is None:
                slice_s = None if cancel is None else 0.2
            else:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    proc.kill()
                    proc.wait()
                    raise subprocess.TimeoutExpired(args, timeout)
                slice_s = remaining if cancel is None else min(0.2, remaining)
            try:
                out, err = proc.communicate(timeout=slice_s)
            except subprocess.TimeoutExpired:
                continue
            return subprocess.CompletedProcess(args, proc.returncode, out, err)
    finally:
        with _PROCS_LOCK:
            _PROCS.discard(proc)


def _terminate(proc) -> None:
    """terminate() -> wait(5) -> kill() -> wait() (section 9)."""
    try:
        proc.terminate()
    except OSError:  # pragma: no cover - already gone
        return
    try:
        proc.wait(timeout=5)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        proc.kill()
        proc.wait()
    except OSError:  # pragma: no cover
        pass


def terminate_all() -> int:
    """Stop every child this process still owns. -> how many were signalled."""
    with _PROCS_LOCK:
        procs = list(_PROCS)
    for proc in procs:
        if proc.poll() is None:
            _terminate(proc)
    return len(procs)


# ---------------------------------------------------------------------------
# 4.4  Paths: ansi_safe, temp_root, job_dir
# ---------------------------------------------------------------------------

_ANSI_MAX = 200


def _short_path(path: str) -> Optional[str]:
    """``GetShortPathNameW``; None when 8.3 generation is off for the volume."""
    if not IS_WINDOWS:  # pragma: no cover - Windows-only API
        return None
    try:
        import ctypes
        from ctypes import wintypes

        GetShortPathNameW = ctypes.windll.kernel32.GetShortPathNameW  # type: ignore[attr-defined]
        GetShortPathNameW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
        GetShortPathNameW.restype = wintypes.DWORD
        need = GetShortPathNameW(path, None, 0)
        if not need:
            return None
        buf = ctypes.create_unicode_buffer(need)
        got = GetShortPathNameW(path, buf, need)
        if not got:
            return None
        short = buf.value
        # 8.3 name generation can be disabled per volume; then short == long.
        return None if short == path else short
    except Exception:  # pragma: no cover
        _LOG.debug("GetShortPathNameW failed for %r", path, exc_info=True)
        return None


def ansi_safe(path: str) -> Optional[str]:
    """Return a form of `path` that ``blastn.exe`` can actually open (section 8.5).

    ``blastn.exe`` has no ``<activeCodePage>UTF-8</>`` in its manifest, so it opens
    files through ``CreateFileA`` and cannot open a path outside the system ANSI
    code page nor one at/beyond MAX_PATH; the ``\\\\?\\`` prefix does not help,
    being a wide-API feature.  Ladder: ACP test -> ``GetShortPathNameW`` -> None,
    which tells the caller to stage the file under :func:`temp_root`.
    On POSIX every path is fine, so the resolved path is returned unchanged.
    """
    resolved = os.path.abspath(path)
    if not IS_WINDOWS:
        return resolved
    try:
        resolved.encode("mbcs")
        ok_acp = True
    except (UnicodeEncodeError, LookupError):
        ok_acp = False
    if ok_acp and len(resolved) < _ANSI_MAX:
        return resolved
    short = _short_path(resolved)
    if short is not None and len(short) < _ANSI_MAX:
        try:
            short.encode("ascii")
        except UnicodeEncodeError:
            return None
        return short
    return None


def _usable_dir(path: str) -> bool:
    """True when `path` exists (or can be made), is ANSI-safe, and is writable."""
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        return False
    if ansi_safe(path) is None:
        return False
    # The probe name must be unique per *call*, not per process: --jobs N runs
    # this concurrently in several threads and a shared name makes them delete
    # each other's file, so the loser's os.remove() raises ENOENT.
    probe = os.path.join(
        path, "wmlst-probe-%d-%d-%s.tmp"
        % (os.getpid(), threading.get_ident(), uuid.uuid4().hex[:8]),
    )
    try:
        with open(probe, "wb") as fh:
            fh.write(b"wmlst")
        os.remove(probe)
    except OSError:
        return False
    return True


def temp_root() -> str:
    """A short, ANSI-safe, writable scratch directory (section 4.4).

    Ladder: ``%WMLST_TMPDIR%`` -> the system temp dir -> its 8.3 short name ->
    ``%LOCALAPPDATA%\\IOWA-Tech\\WMLST\\tmp`` -> ``C:\\ProgramData\\IOWA-Tech\\WMLST\\tmp``
    -> ``%SystemDrive%\\WMLST-tmp``.  Each candidate is validated by writing and
    deleting a probe file.
    """
    candidates: List[str] = []
    env_tmp = os.environ.get("WMLST_TMPDIR")
    if env_tmp:
        candidates.append(env_tmp)
    sys_tmp = tempfile.gettempdir()
    candidates.append(sys_tmp)
    if IS_WINDOWS:
        short = _short_path(os.path.abspath(sys_tmp))
        if short:
            candidates.append(short)
        candidates.append(os.path.join(install_root(), "tmp"))
        programdata = os.environ.get("ProgramData", r"C:\ProgramData")
        candidates.append(os.path.join(programdata, "IOWA-Tech", "WMLST", "tmp"))
        sysdrive = os.environ.get("SystemDrive", "C:")
        candidates.append(os.path.join(sysdrive + os.sep, "WMLST-tmp"))

    for cand in candidates:
        if _usable_dir(cand):
            return os.path.abspath(cand)
    raise _exc("WmlstError")(
        "No writable temporary directory found (tried: %s)" % (", ".join(candidates),),
        user_message="WMLST could not find a usable temporary folder. "
        "Set WMLST_TMPDIR to a folder you can write to.",
    )


def _on_rm_error(func, path, exc_info) -> None:
    """rmtree hardening: clear the read-only bit and retry; then log and continue."""
    for _ in range(5):
        try:
            os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
        except OSError:
            pass
        try:
            func(path)
            return
        except OSError:
            time.sleep(0.1)
    _LOG.warning("could not remove %s during temp cleanup", path)


def _rmtree(path: str) -> None:
    if sys.version_info >= (3, 12):
        shutil.rmtree(path, onexc=lambda f, p, e: _on_rm_error(f, p, e))
    else:  # pragma: no cover - exercised on 3.9-3.11
        shutil.rmtree(path, onerror=_on_rm_error)


@contextlib.contextmanager
def job_dir(root: Optional[str] = None):
    """Yield a private scratch directory and remove it afterwards (sections 4.4, 9).

    ONE DIRECTORY PER CONCURRENTLY-EXECUTING JOB: the Perl reuses a single dir for
    the whole run, which is safe only because it is strictly serial.  Cleanup never
    aborts a scan — a directory that refuses to die is logged and left behind.
    """
    base = root or temp_root()
    while True:
        path = os.path.join(base, "wmlst-%d-%d" % (os.getpid(), next(_JOB_COUNTER)))
        try:
            os.mkdir(path)
            break
        except FileExistsError:  # pragma: no cover - counter collision
            continue
    try:
        yield path
    finally:
        try:
            _rmtree(path)
        except Exception:  # pragma: no cover - cleanup must never raise
            _LOG.warning("temp cleanup failed for %s", path, exc_info=True)


# ---------------------------------------------------------------------------
# 4.4  Child environment
# ---------------------------------------------------------------------------

def child_env(blastdb_dir: str) -> Dict[str, str]:
    """The environment every BLAST+ child gets (section 4.4).

    Disabling the usage report is what stops the Windows Defender Firewall prompt
    appearing over the Tkinter window on first run, and is the right default for
    patient-derived isolates.
    """
    env = os.environ.copy()
    safe = ansi_safe(blastdb_dir) if blastdb_dir else None
    env["BLASTDB"] = safe if safe is not None else os.path.abspath(blastdb_dir or ".")
    env["BLAST_USAGE_REPORT"] = "false"
    env["NCBI_USAGE_REPORT_ENABLED"] = "0"
    env["NCBI_DONT_USE_NCBIRC"] = "1"
    env["NCBI_DONT_USE_LOCAL_CONFIG"] = "1"
    env.pop("BLASTDB_LMDB_MAP_SIZE", None)
    env.pop("NCBI_CONFIG_OVERRIDES", None)
    return env


# ---------------------------------------------------------------------------
# 4.4 / 5.5  argv builders
# ---------------------------------------------------------------------------

def _fmt_num(value: float) -> str:
    """Perl's default numeric stringification: 95.0 -> '95', 95.5 -> '95.5'."""
    return "%.15g" % (float(value),)


def _argv_path(path: str, what: str) -> str:
    safe = ansi_safe(path)
    if safe is None:
        raise _exc("WmlstError")(
            "%s path cannot be passed to blastn: %r" % (what, path),
            user_message="The folder used for %s has characters or a length that "
            "BLAST+ cannot open. Set WMLST_TMPDIR to a short, plain-ASCII folder."
            % (what,),
        )
    return safe


def blastn_argv(
    tools: BlastTools,
    *,
    query: str,
    out: str,
    db_basename: str,
    threads: int,
    minid: float,
) -> List[str]:
    """Exact transliteration of ``bin/mlst:309-315`` (section 5.5).

    No quote characters anywhere: the outfmt spec is ONE argv element, not a
    shell-quoted string.  `db_basename` MUST be the bare index basename, with the
    caller setting ``cwd`` and ``BLASTDB`` to the directory that holds it — BLAST
    splits the ``-db`` value on whitespace internally, so an argv list does not
    rescue a database path containing a space.

    Every flag is load-bearing; do not add, remove or reorder any of them, and
    never pass ``-task``, ``-mt_mode`` or ``-lcase_masking``.
    """
    return [
        # argv[0] is handed to CreateProcessW (a wide API) and is never opened
        # as a data file, so the CreateFileA-motivated ansi_safe() gate must NOT
        # be applied to it (section 5.5 spells the argv out with a raw exe path).
        # probe_version() already launches this very path raw; gating it here
        # made discovery accept an install that every later search refused.
        tools.blastn,
        "-query",
        _argv_path(query, "query"),
        "-out",
        _argv_path(out, "output"),
        "-db",
        db_basename,
        "-num_threads",
        str(int(threads)),
        "-ungapped",
        "-dust",
        "no",
        "-word_size",
        "32",
        "-max_target_seqs",
        "100000",
        "-perc_identity",
        _fmt_num(minid),
        "-evalue",
        "1E-20",
        "-outfmt",
        BLAST_OUTFMT,
    ]


def makeblastdb_argv(tools: BlastTools, fasta: str) -> List[str]:
    """``scripts/mlst-make_blast_db:21``. The title is exactly ``PubMLST``.

    No vendor string: the title is embedded in the index bytes (section 4.2).

    Only the BASENAME of `fasta` ever reaches the command line: makeblastdb
    splits the ``-in`` value on whitespace exactly like blastn splits ``-db``
    (an absolute ``C:\\Program Files\\...`` path exits 1 with "Please provide a
    database name using -out"), so the caller MUST set ``cwd`` to the directory
    that holds the FASTA.  That is necessary but NOT sufficient — see
    :func:`run_makeblastdb`, which is what callers should use.
    """
    return [
        # argv[0] raw: see blastn_argv(). CreateProcessW, not CreateFileA.
        tools.makeblastdb,
        "-hash_index",
        "-in",
        os.path.basename(fasta),
        "-dbtype",
        "nucl",
        "-title",
        "PubMLST",
        "-parse_seqids",
    ]


_WHITESPACE_RE = re.compile(r"\s")


def _has_whitespace(path: str) -> bool:
    return bool(_WHITESPACE_RE.search(path))


def _scratch_root_without_whitespace() -> str:
    """A writable scratch directory whose absolute path has no whitespace."""
    candidates: List[str] = []
    try:
        root = temp_root()
    except Exception:  # pragma: no cover - temp_root already explains itself
        root = None
    if root:
        candidates.append(root)
        if IS_WINDOWS:  # pragma: no cover - Windows-only API
            short = _short_path(root)
            if short:
                candidates.append(short)
    if IS_WINDOWS:  # pragma: no cover - Windows-only fallbacks
        sysdrive = os.environ.get("SystemDrive", "C:")
        candidates.append(os.path.join(sysdrive + os.sep, "WMLST-tmp"))
    else:
        candidates.append(tempfile.gettempdir())
    for cand in candidates:
        resolved = os.path.abspath(cand)
        if _has_whitespace(resolved):
            continue
        if _usable_dir(resolved):
            return resolved
    raise _exc("WmlstError")(
        "No whitespace-free temporary directory found (tried: %s)"
        % (", ".join(candidates) if candidates else "nowhere",),
        user_message="WMLST could not find a temporary folder without spaces in "
        "its path, which makeblastdb needs. Set WMLST_TMPDIR to a folder such "
        "as C:\\WMLST-tmp.",
    )


def _makeblastdb_build_dir(directory: str) -> Tuple[str, bool]:
    """Where makeblastdb may run for an index destined for `directory`.

    Returns ``(build_dir, relocate)``; when `relocate` is True the caller must
    move the finished index files into `directory` itself.
    """
    resolved = os.path.abspath(directory)
    if not _has_whitespace(resolved) and ansi_safe(resolved) is not None:
        return resolved, False
    if IS_WINDOWS:  # pragma: no cover - Windows-only API
        short = _short_path(resolved)
        if short and not _has_whitespace(short) and ansi_safe(short) is not None:
            # The very same directory under its 8.3 alias: nothing to move.
            return short, False
    root = _scratch_root_without_whitespace()
    return tempfile.mkdtemp(prefix="wmlst-mkdb-", dir=root), True


def run_makeblastdb(
    tools: BlastTools,
    fasta: str,
    *,
    timeout: Optional[float] = None,
    cancel: Optional[threading.Event] = None,
) -> subprocess.CompletedProcess:
    """Index `fasta` in place, next to itself, whatever its path (sections 4.9, 8.8).

    makeblastdb splits on whitespace TWICE: once in the ``-in`` value, and again
    in the absolute path under which it re-opens the database it has just built
    for the metadata pass.  So under ``C:\\Program Files\\...`` an absolute ``-in``
    exits 1 ("Please provide a database name using -out") and a bare ``-in``
    with ``cwd`` set exits 2 ("No alias or index file found ... [C:\\Program]")
    after writing 11 of the 12 files — the ``.njs`` is never produced.  Neither
    ``-out``, nor dropping ``-parse_seqids``/``-hash_index``, rescues it: no argv
    arrangement can write a database into a directory whose path has a space.

    The index files are fully relocatable (section 8.8), so when the destination
    cannot be indexed directly this builds under a whitespace-free scratch root
    and moves the artifacts into place.  Returns the CompletedProcess unchanged;
    exit-code interpretation stays the caller's job.
    """
    fasta = os.path.abspath(fasta)
    directory = os.path.dirname(fasta) or os.curdir
    name = os.path.basename(fasta)
    build_dir, relocate = _makeblastdb_build_dir(directory)
    argv = makeblastdb_argv(tools, name)
    if not relocate:
        return run_tool(argv, cwd=build_dir, env=child_env(build_dir),
                        timeout=timeout, cancel=cancel)
    _LOG.debug("makeblastdb staged into %r for %r", build_dir, directory)
    try:
        shutil.copyfile(fasta, os.path.join(build_dir, name))
        completed = run_tool(argv, cwd=build_dir, env=child_env(build_dir),
                             timeout=timeout, cancel=cancel)
        if completed.returncode == 0:
            for entry in sorted(os.listdir(build_dir)):
                if entry == name or not entry.startswith(name + "."):
                    continue
                dest = os.path.join(directory, entry)
                if os.path.lexists(dest):
                    os.remove(dest)
                shutil.move(os.path.join(build_dir, entry), dest)
        return completed
    finally:
        try:
            _rmtree(build_dir)
        except Exception:  # pragma: no cover - cleanup must never raise
            _LOG.warning("could not remove %s", build_dir, exc_info=True)


def run_blastn(
    tools: BlastTools,
    *,
    query: str,
    out: str,
    blastdb: str,
    threads: int = 1,
    minid: float = 95.0,
    timeout: Optional[float] = 900.0,
    cancel: Optional[threading.Event] = None,
) -> subprocess.CompletedProcess:
    """Convenience wrapper that gets ``-db``/``cwd``/``BLASTDB`` right (section 5.5).

    `blastdb` is the full index prefix, e.g. ``.../db/blast/mlst.fa``; only its
    basename reaches the command line.  Returns the CompletedProcess unchanged —
    exit-code interpretation is the engine's job (section 8.6).
    """
    db_dir = os.path.dirname(os.path.abspath(blastdb)) or os.curdir
    # Diagnose a missing index HERE, while we still know what is missing. cwd is
    # the database directory, so if it does not exist Popen fails with a bare
    # "No such file or directory: <db dir>" attributed to the blastn path -- which
    # reads as "BLAST is not installed" and sends the user off to download 137 MB
    # that will not help. The index is derived and rebuildable; say so.
    if not os.path.isdir(db_dir):
        raise _exc("DatabaseMissingError")(
            "the BLAST index directory does not exist: %s" % db_dir,
            user_message=("The search index has not been built yet. WMLST can "
                          "rebuild it from the allele files."),
        )
    argv = blastn_argv(
        tools,
        query=query,
        out=out,
        db_basename=os.path.basename(blastdb),
        threads=threads,
        minid=minid,
    )
    return run_tool(
        argv,
        cwd=_argv_path(db_dir, "database"),
        env=child_env(db_dir),
        timeout=timeout,
        cancel=cancel,
    )


# ---------------------------------------------------------------------------
# 4.4  Version probing
# ---------------------------------------------------------------------------

_VERSION_RE = re.compile(r"(\d+)\.(\d+)\.(\d+)(\+?)", re.ASCII)


def probe_version(blastn: str) -> Tuple[str, Tuple[int, ...]]:
    """Run ``blastn -version`` under :func:`run_tool`. -> ('2.17.0+', (2, 17, 0)).

    A loader failure (``0xC0000135`` / ``0xC0000142``) with empty output raises
    ``BootstrapError`` naming the VC++ redistributable (sections 4.4, 8.3).
    """
    proc = run_tool([blastn, "-version"], timeout=60)
    text = (proc.stdout or "") + (proc.stderr or "")
    rc = proc.returncode
    unsigned = rc & 0xFFFFFFFF if isinstance(rc, int) and rc < 0 else rc
    if unsigned in (_STATUS_DLL_NOT_FOUND, _STATUS_DLL_INIT_FAILED) and not text.strip():
        raise _exc("BootstrapError")(
            "blastn could not start (0x%08X): a required DLL is missing" % (unsigned,),
            user_message=(
                "BLAST+ could not start because the Microsoft Visual C++ runtime is "
                "missing. Install it from %s and try again." % (_VC_REDIST_URL,)
            ),
        )
    match = _VERSION_RE.search(text)
    if proc.returncode != 0 or not match:
        raise _exc("BlastNotFoundError")(
            "Could not determine blastn version from %r (exit %s): %r"
            % (blastn, proc.returncode, text.strip()[:200]),
            user_message="WMLST could not run BLAST+ at '%s'." % (blastn,),
        )
    major, minor, patch = (int(match.group(i)) for i in (1, 2, 3))
    version = "%d.%d.%d%s" % (major, minor, patch, match.group(4) or "+")
    return version, (major, minor, patch)


# ---------------------------------------------------------------------------
# 4.4  Discovery ladder
# ---------------------------------------------------------------------------

def _exe(name: str) -> str:
    return name + ".exe" if IS_WINDOWS else name


def install_root() -> str:
    """``%LOCALAPPDATA%\\IOWA-Tech\\WMLST`` on Windows, ``~/.local/share/wmlst`` elsewhere."""
    if IS_WINDOWS:
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser(
            r"~\AppData\Local"
        )
        return os.path.join(base, "IOWA-Tech", "WMLST")
    base = os.environ.get("XDG_DATA_HOME") or os.path.join(
        os.path.expanduser("~"), ".local", "share"
    )
    return os.path.join(base, "wmlst")


def write_sha256sums(directory: str) -> str:
    """Write ``SHA256SUMS`` for every regular file in `directory`. -> its path."""
    lines = []
    for name in sorted(os.listdir(directory)):
        if name in (SHA256SUMS, INSTALL_OK):
            continue
        full = os.path.join(directory, name)
        if not os.path.isfile(full):
            continue
        lines.append("%s  %s\n" % (_sha256_file(full), name))
    path = os.path.join(directory, SHA256SUMS)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(lines)
    return path


def verify_sha256sums(directory: str) -> bool:
    """True when every file listed in ``SHA256SUMS`` is present and matches."""
    path = os.path.join(directory, SHA256SUMS)
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            rows = fh.read().splitlines()
    except OSError:
        return False
    if not rows:
        return False
    for row in rows:
        row = row.strip()
        if not row:
            continue
        parts = row.split(None, 1)
        if len(parts) != 2:
            return False
        digest, name = parts[0], parts[1].lstrip("*").strip()
        full = os.path.join(directory, name)
        if not os.path.isfile(full) or _sha256_file(full) != digest:
            _LOG.warning("checksum mismatch for %s", full)
            return False
    return True


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _tools_from_dir(directory: str, origin: str) -> Optional[BlastTools]:
    """Build a :class:`BlastTools` from `directory` if a runnable blastn lives there."""
    blastn = os.path.join(directory, _exe("blastn"))
    if not os.path.isfile(blastn):
        return None
    return _tools_from_blastn(blastn, origin)


def _tools_from_blastn(blastn: str, origin: str) -> Optional[BlastTools]:
    blastn = os.path.abspath(blastn)
    directory = os.path.dirname(blastn)
    try:
        version, vtuple = probe_version(blastn)
    except Exception as exc:
        _LOG.debug("candidate %s rejected: %s", blastn, exc)
        return None
    if vtuple < BLAST_MIN_VERSION:
        _LOG.warning(
            "ignoring blastn %s at %s: WMLST needs %s or newer",
            version,
            blastn,
            ".".join(str(n) for n in BLAST_MIN_VERSION),
        )
        return None
    if version != BLAST_PINNED_VERSION:
        # Informational, not a warning: any build at or above BLAST_MIN_VERSION
        # is supported, and the one below rejects anything older. As a warning it
        # survived --quiet and broke stderr byte-parity with upstream on every
        # distro-packaged BLAST -- Debian and Ubuntu ship 2.12, not the pinned
        # 2.17.0+. The neighbouring "Found blastn" line is already an info line.
        _LOG.info(
            "blastn %s found at %s; WMLST is validated against %s",
            version,
            blastn,
            BLAST_PINNED_VERSION,
        )
    makeblastdb = os.path.join(directory, _exe("makeblastdb"))
    blastdbcmd = os.path.join(directory, _exe("blastdbcmd"))
    if not os.path.isfile(makeblastdb):
        found = shutil.which("makeblastdb", path=directory) or shutil.which("makeblastdb")
        makeblastdb = found or makeblastdb
    return BlastTools(
        blastn=blastn,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd if os.path.isfile(blastdbcmd) else None,
        version=version,
        version_tuple=vtuple,
        origin=origin,
    )


def _frozen_bundle_dir() -> Optional[str]:
    if not getattr(sys, "frozen", False):
        return None
    base = getattr(sys, "_MEIPASS", None) or os.path.dirname(os.path.abspath(sys.executable))
    return os.path.join(base, "blast")


def _is_world_writable(directory: str) -> bool:
    """POSIX-only binary-planting guard for the PATH rung (section 4.4).

    Windows has no POSIX mode bits: CPython synthesises them from the file
    attributes (``attributes_to_mode()`` in Python/fileutils.c), so EVERY
    directory without FILE_ATTRIBUTE_READONLY reports ``0o40777`` and every
    read-only one reports ``0o40555``.  Testing S_IWOTH there does not measure
    permissiveness at all — it inverts it, rejecting every ordinary BLAST+
    install directory on PATH and accepting read-only ones.  Report False on
    Windows and leave the (platform-neutral) CWD check to do the work; a real
    Windows guard needs a DACL inspection, which belongs in its own function.
    """
    if IS_WINDOWS:
        return False
    try:
        mode = os.stat(directory).st_mode
    except OSError:  # pragma: no cover
        return True
    return bool(mode & stat.S_IWOTH)


def _conda_dirs() -> List[str]:
    prefix = os.environ.get("CONDA_PREFIX")
    dirs: List[str] = []
    if not prefix:
        return dirs
    if IS_WINDOWS:
        dirs.append(os.path.join(prefix, "Library", "bin"))
        dirs.append(os.path.join(prefix, "Scripts"))
    else:
        dirs.append(os.path.join(prefix, "bin"))
    envs = os.path.join(os.path.dirname(os.path.dirname(prefix)), "envs")
    if os.path.isdir(envs):
        for name in sorted(os.listdir(envs)):
            sib = os.path.join(envs, name)
            dirs.append(
                os.path.join(sib, "Library", "bin") if IS_WINDOWS
                else os.path.join(sib, "bin")
            )
    return dirs


def find_blast(explicit: Optional[str] = None) -> BlastTools:
    """Find a usable BLAST+ installation; first hit wins (section 4.4).

    Ladder: `explicit` -> frozen bundle -> ``%WMLST_BLAST_DIR%`` -> the per-user
    install root (only with a valid ``.wmlst-install-ok`` and ``SHA256SUMS``) ->
    ``PATH`` (rejecting a directory that is the CWD or world-writable, which is the
    binary-planting guard) -> conda prefixes.  Raises ``BlastNotFoundError``.
    """
    tried: List[str] = []

    if explicit:
        tried.append(explicit)
        path = explicit
        if os.path.isdir(path):
            tools = _tools_from_dir(path, "explicit")
        else:
            tools = _tools_from_blastn(path, "explicit") if os.path.isfile(path) else None
            if tools is None:
                which = shutil.which(explicit)
                tools = _tools_from_blastn(which, "explicit") if which else None
        if tools is not None:
            return tools

    bundle = _frozen_bundle_dir()
    if bundle:
        tried.append(bundle)
        tools = _tools_from_dir(bundle, "bundled")
        if tools is not None:
            return tools

    env_dir = os.environ.get("WMLST_BLAST_DIR")
    if env_dir:
        tried.append(env_dir)
        tools = _tools_from_dir(env_dir, "env")
        if tools is not None:
            return tools

    appdata = os.path.join(install_root(), "blast")
    tried.append(appdata)
    if os.path.isfile(os.path.join(appdata, INSTALL_OK)) and verify_sha256sums(appdata):
        tools = _tools_from_dir(appdata, "appdata")
        if tools is not None:
            return tools
    elif os.path.isdir(appdata):
        _LOG.warning(
            "ignoring %s: missing %s or failed checksum verification",
            appdata,
            INSTALL_OK,
        )

    which = shutil.which("blastn")
    if which:
        which = os.path.realpath(which)
        parent = os.path.dirname(which)
        tried.append(which)
        if os.path.abspath(parent) == os.path.abspath(os.getcwd()):
            _LOG.warning("ignoring blastn in the current directory: %s", which)
        elif _is_world_writable(parent):
            _LOG.warning("ignoring blastn in a world-writable directory: %s", which)
        else:
            tools = _tools_from_blastn(which, "path")
            if tools is not None:
                return tools

    for directory in _conda_dirs():
        tried.append(directory)
        tools = _tools_from_dir(directory, "conda")
        if tools is not None:
            return tools

    raise _exc("BlastNotFoundError")(
        "blastn not found (looked in: %s)" % (", ".join(tried) if tried else "nowhere",),
        user_message=(
            "WMLST could not find the BLAST+ search engine. "
            + (
                "Use Database -> Install BLAST+ to download it automatically."
                if IS_WINDOWS
                else "Install ncbi-blast+ with your package manager, or set WMLST_BLAST_DIR."
            )
        ),
    )


# ---------------------------------------------------------------------------
# 4.4  Bootstrap downloader (Windows payload)
# ---------------------------------------------------------------------------

_MD5LINE_RE = re.compile(r"^([0-9a-f]{32})\s+(\S+)$", re.ASCII)
_BAD_CHARS_RE = re.compile(r'[:*?"<>|]', re.ASCII)


def _http_get(url: str, timeout: float = 60.0):
    req = urllib.request.Request(url, headers={"User-Agent": "WMLST"})
    return urllib.request.urlopen(req, timeout=timeout)


def _expected_md5(url: str) -> str:
    """Step 1: fetch ``<url>.md5`` and check that it names the archive we want."""
    try:
        with _http_get(url + ".md5") as resp:
            text = resp.read(4096).decode("ascii", "replace").strip()
    except urllib.error.URLError as exc:
        raise _exc("BootstrapError")(
            "Could not fetch %s.md5: %s" % (url, exc),
            user_message="WMLST could not reach the NCBI download server.",
        ) from exc
    match = _MD5LINE_RE.match(text.splitlines()[0] if text else "")
    if not match:
        raise _exc("BootstrapError")(
            "Unparseable md5 file: %r" % (text[:120],),
            user_message="The BLAST+ checksum file from NCBI was not readable.",
        )
    if os.path.basename(match.group(2)) != os.path.basename(url):
        raise _exc("BootstrapError")(
            "md5 file names %r, expected %r"
            % (match.group(2), os.path.basename(url)),
            user_message="The BLAST+ checksum file from NCBI did not match the download.",
        )
    return match.group(1)


def _content_length(url: str) -> Optional[int]:
    """Step 2: HEAD the archive; a mismatch is a warning, never an abort."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "WMLST"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.headers.get("Content-Length")
            return int(raw) if raw else None
    except (urllib.error.URLError, ValueError, OSError):
        return None


def _safe_member(member: tarfile.TarInfo) -> bool:
    """Step 7's guard: reject anything that is not a plain file in a plain place."""
    name = member.name.replace("\\", "/")
    if not member.isfile():
        return False
    if member.issym() or member.islnk() or member.ischr() or member.isblk() or member.isfifo():
        return False
    if name.startswith("/") or name.startswith("../") or "/../" in name or name == "..":
        return False
    if re.match(r"^[A-Za-z]:", name, re.ASCII):
        return False
    if _BAD_CHARS_RE.search(name):
        return False
    return os.path.basename(name) in BLAST_MEMBERS


def bootstrap(
    dest: Optional[str] = None,
    *,
    progress: Optional[Callable[[int, int, str], None]] = None,
    cancel: Optional[threading.Event] = None,
) -> BlastTools:
    """Download, verify, selectively extract and install BLAST+ (section 4.4).

    Download -> verify MD5 -> selective extract -> probe -> atomic swap, with the
    ``.wmlst-install-ok`` sentinel written LAST so a half-extracted tree from an AV
    quarantine is never trusted.  `progress` is
    ``Callable[[bytes_done, bytes_total, phase_text], None]``; `cancel` is honoured
    between chunks and between members.

    urllib is mandatory here (not PowerShell): urllib attaches no Mark-of-the-Web,
    so the extracted ``blastn.exe`` never needs ``Unblock-File``.
    """
    root = os.path.abspath(dest or install_root())
    blast_dir = os.path.join(root, "blast")
    dl_dir = os.path.join(root, ".dl")
    stage_dir = os.path.join(root, ".stage")
    name = os.path.basename(BLAST_URL)

    def report(done: int, total: int, text: str) -> None:
        if progress is not None:
            progress(done, total, text)

    def check_cancel() -> None:
        if cancel is not None and cancel.is_set():
            raise _exc("Cancelled")("BLAST+ installation cancelled")

    os.makedirs(dl_dir, exist_ok=True)
    if os.path.isdir(stage_dir):
        _rmtree(stage_dir)
    os.makedirs(stage_dir, exist_ok=True)

    report(0, BLAST_BYTES, "Checking the download")
    expected_md5 = _expected_md5(BLAST_URL)
    if expected_md5 != BLAST_MD5:
        raise _exc("BootstrapError")(
            "NCBI now publishes md5 %s for %s; WMLST is pinned to %s"
            % (expected_md5, name, BLAST_MD5),
            user_message="The BLAST+ download on the NCBI server has changed. "
            "WMLST will not install an unverified build.",
        )
    length = _content_length(BLAST_URL)
    if length is not None and length != BLAST_BYTES:
        _LOG.warning("Content-Length %s != expected %s", length, BLAST_BYTES)

    total = length or BLAST_BYTES
    free = shutil.disk_usage(root).free
    needed = BLAST_BYTES + 40 * 1024 * 1024
    if free < needed:
        raise _exc("OutOfDiskError")(
            "%d bytes free at %s, need %d" % (free, root, needed),
            user_message="There is not enough free disk space to install BLAST+ "
            "(about %d MB is needed)." % (needed // (1024 * 1024),),
        )

    archive = os.path.join(dl_dir, name)
    part = archive + ".part"
    digest = hashlib.md5()
    done = 0
    report(0, total, "Downloading BLAST+ (about 137 MB)")
    try:
        with _http_get(BLAST_URL, timeout=120) as resp, open(part, "wb") as fh:
            while True:
                check_cancel()
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                digest.update(chunk)
                done += len(chunk)
                report(done, total, "Downloading BLAST+ (about 137 MB)")
    except BaseException:
        _quiet_remove(part)
        raise
    if digest.hexdigest() != BLAST_MD5:
        got = digest.hexdigest()
        _quiet_remove(part)
        raise _exc("BootstrapError")(
            "md5 mismatch: got %s, expected %s" % (got, BLAST_MD5),
            user_message="The BLAST+ download was damaged in transit "
            "(checksum %s, expected %s). Please try again." % (got, BLAST_MD5),
        )
    os.replace(part, archive)

    report(done, total, "Installing BLAST+")
    wanted = set(BLAST_MEMBERS)
    extracted: List[str] = []
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            check_cancel()
            if not _safe_member(member):
                continue
            base = os.path.basename(member.name.replace("\\", "/"))
            if base not in wanted:
                continue
            src = tar.extractfile(member)
            if src is None:  # pragma: no cover - isfile() already checked
                continue
            target = os.path.join(stage_dir, base)
            with src, open(target, "wb") as out:
                shutil.copyfileobj(src, out, 1 << 20)
            if (member.mode & 0o111) or base.lower().endswith((".exe", ".dll")):
                os.chmod(target, 0o755)
            extracted.append(base)
            wanted.discard(base)
            report(done, total, "Installing BLAST+ (%s)" % (base,))

    missing = [m for m in ("blastn.exe", "makeblastdb.exe", "nghttp2.dll") if m not in extracted]
    if IS_WINDOWS and missing:
        raise _exc("BootstrapError")(
            "archive did not contain %s" % (", ".join(missing),),
            user_message="The BLAST+ download was incomplete (%s missing)."
            % (", ".join(missing),),
        )

    write_sha256sums(stage_dir)

    staged_blastn = os.path.join(stage_dir, _exe("blastn"))
    if os.path.isfile(staged_blastn):
        probe_version(staged_blastn)

    if os.path.isdir(blast_dir):
        old = blast_dir + ".old-%d" % (os.getpid(),)
        os.replace(blast_dir, old)
        try:
            _rmtree(old)
        except Exception:  # pragma: no cover
            _LOG.warning("could not remove %s", old)
    os.replace(stage_dir, blast_dir)
    with open(os.path.join(blast_dir, INSTALL_OK), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(BLAST_PINNED_VERSION + "\n")

    if not os.environ.get("WMLST_KEEP_DOWNLOAD"):
        try:
            _rmtree(dl_dir)
        except Exception:  # pragma: no cover
            _LOG.warning("could not remove %s", dl_dir)

    report(total, total, "BLAST+ installed")
    tools = _tools_from_dir(blast_dir, "appdata")
    if tools is None:
        raise _exc("BootstrapError")(
            "installed tree at %s is not runnable" % (blast_dir,),
            user_message="BLAST+ was installed but would not run.",
        )
    return tools


def _quiet_remove(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass
