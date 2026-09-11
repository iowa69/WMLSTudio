# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Shared pytest fixtures for the whole WMLST suite (docs/ARCHITECTURE.md 13.4).

Three jobs:

1. Make ``import wmlst`` work from a bare checkout, so a contributor can run
   ``python -m pytest`` before ``pip install -e .``.
2. Find a real ``blastn``/``makeblastdb`` -- on PATH, via ``WMLST_BLAST_DIR``,
   or in an unpacked ``ncbi-blast-*`` tree lying about the machine -- and put it
   on PATH for the whole session, because most of the suite locates BLAST with
   a plain ``shutil.which("blastn")``.
3. Turn the marker vocabulary declared in pyproject.toml into automatic skips,
   so a machine without BLAST, without a display, without Windows or without
   network reports SKIPPED rather than a wall of red.

Set ``WMLST_TEST_NO_PATH_INJECT=1`` to disable job 2 and test exactly what is
on the ambient PATH.
"""

from __future__ import annotations

import glob
import os
import shutil
import sys

import pytest

#: Repository root, resolved from this file rather than from the cwd so the
#: suite behaves identically under `pytest`, `pytest tests/`, and an IDE runner.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

#: Executable base names, with the Windows suffix applied where it belongs.
_EXE = ".exe" if os.name == "nt" else ""
BLASTN = "blastn" + _EXE
MAKEBLASTDB = "makeblastdb" + _EXE


def _candidate_blast_dirs():
    """Yield directories that might contain blastn, best guess first.

    The globs matter for local development: the reference BLAST+ 2.17.0 tree is
    routinely unpacked into a scratch directory rather than installed, and the
    exact path is session-specific, so it is searched for rather than pinned.
    """
    explicit = [
        os.environ.get("WMLST_TEST_BLAST_DIR"),
        os.environ.get("WMLST_BLAST_DIR"),
    ]
    for path in explicit:
        if path:
            yield path

    patterns = [
        os.path.join(REPO_ROOT, "vendor", "ncbi-blast-*", "bin"),
        os.path.join(REPO_ROOT, os.pardir, "ncbi-blast-*", "bin"),
        os.path.join(REPO_ROOT, "packaging", "payload", "blast"),
        os.path.join(
            os.path.expanduser("~"), "AppData", "Local", "IOWA-Tech", "WMLST", "blast"
        ),
        os.path.join(os.path.expanduser("~"), ".local", "share", "wmlst", "blast"),
        os.path.join(os.sep + "tmp", "claude-*", "*", "*", "scratchpad", "ncbi-blast-*", "bin"),
        os.path.join(os.sep + "tmp", "ncbi-blast-*", "bin"),
        os.path.join(os.sep + "opt", "ncbi-blast-*", "bin"),
    ]
    for pattern in patterns:
        for path in sorted(glob.glob(pattern)):
            yield path


def find_blast_dir():
    """Return a directory holding both blastn and makeblastdb, or None.

    Checks PATH first so a system install always wins over a scratch copy.
    """
    found = shutil.which("blastn")
    if found and shutil.which("makeblastdb"):
        return os.path.dirname(os.path.abspath(found))
    for directory in _candidate_blast_dirs():
        if not directory or not os.path.isdir(directory):
            continue
        if os.path.isfile(os.path.join(directory, BLASTN)) and os.path.isfile(
            os.path.join(directory, MAKEBLASTDB)
        ):
            return os.path.abspath(directory)
    return None


# The PATH injection happens at import time, not inside a fixture, because
# modules such as test_cli_bats_parity call shutil.which() at collection time to
# decide whether to skip. A session fixture would run too late.
_BLAST_DIR = None
if not os.environ.get("WMLST_TEST_NO_PATH_INJECT"):
    _BLAST_DIR = find_blast_dir()
    if _BLAST_DIR:
        os.environ["PATH"] = _BLAST_DIR + os.pathsep + os.environ.get("PATH", "")
        os.environ.setdefault("WMLST_BLAST_DIR", _BLAST_DIR)
else:  # pragma: no cover - opt-out path
    _BLAST_DIR = find_blast_dir()

# NCBI's BLAST+ phones home with a usage beacon unless this is set. A test run
# must never make an unexpected outbound request.
os.environ.setdefault("BLAST_USAGE_REPORT", "false")
os.environ.setdefault("PYTHONUTF8", "1")


# ---------------------------------------------------------------------------
# Path fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def repo_root():
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def db_dir(repo_root):
    """Absolute path to ``db/`` -- the committed PubMLST snapshot (section 6.1)."""
    return os.path.join(repo_root, "db")


@pytest.fixture(scope="session")
def pubmlst_dir(db_dir):
    """Absolute path to ``db/pubmlst/``: 162 scheme directories."""
    return os.path.join(db_dir, "pubmlst")


@pytest.fixture(scope="session")
def blast_index(db_dir):
    """Absolute path to the derived ``db/blast/mlst.fa``.

    Skips when the index has not been built; it is gitignored and rebuildable
    with ``wmlst --make-blast-db``.
    """
    path = os.path.join(db_dir, "blast", "mlst.fa")
    if not os.path.isfile(path + ".nsq") and not os.path.isfile(path + ".nin"):
        pytest.skip("db/blast is not built - run: python -m wmlst --make-blast-db")
    return path


@pytest.fixture(scope="session")
def data_dir(repo_root):
    """Absolute path to ``tests/data/`` -- the 13 input fixtures."""
    return os.path.join(repo_root, "tests", "data")


@pytest.fixture(scope="session")
def golden_dir(repo_root):
    """Absolute path to ``tests/golden/`` -- real Perl mlst 2.35.0 output."""
    return os.path.join(repo_root, "tests", "golden")


@pytest.fixture(scope="session")
def fixture_path(data_dir):
    """Return a callable mapping a fixture base name to its absolute path."""

    def _resolve(name):
        path = os.path.join(data_dir, name)
        if not os.path.exists(path):
            raise AssertionError("no such test fixture: %s" % path)
        return path

    return _resolve


@pytest.fixture(scope="session")
def golden_bytes(golden_dir):
    """Return a callable mapping a golden file name to its exact bytes.

    Bytes, never text: the golden corpus is a byte-identity surface and
    decoding it would hide exactly the CRLF bug section 8.4 is about.
    """

    def _read(name):
        path = os.path.join(golden_dir, name)
        if not os.path.isfile(path):
            raise AssertionError("no such golden file: %s" % path)
        with open(path, "rb") as handle:
            return handle.read()

    return _read


# ---------------------------------------------------------------------------
# BLAST fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def blast_dir():
    """Directory holding blastn and makeblastdb, or skip the test."""
    if _BLAST_DIR is None:
        pytest.skip(
            "no blastn/makeblastdb found. Put NCBI BLAST+ on PATH or set "
            "WMLST_TEST_BLAST_DIR to an ncbi-blast-*/bin directory."
        )
    return _BLAST_DIR


@pytest.fixture(scope="session")
def blastn_path(blast_dir):
    """Absolute path to a real ``blastn`` executable."""
    return os.path.join(blast_dir, BLASTN)


@pytest.fixture(scope="session")
def makeblastdb_path(blast_dir):
    """Absolute path to a real ``makeblastdb`` executable."""
    return os.path.join(blast_dir, MAKEBLASTDB)


@pytest.fixture
def clean_env(monkeypatch):
    """A child environment with every WMLST/BLAST variable removed.

    Use it whenever a test asserts default behaviour: a developer with
    MLST_DBDIR exported would otherwise get a silent pass (divergence D15).
    """
    for name in (
        "WMLST_DBDIR", "MLST_DBDIR", "WMLST_BLAST_DIR", "WMLST_TEST_BLAST_DIR",
        "WMLST_JOBS", "WMLST_QUIET", "BLASTDB",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("BLAST_USAGE_REPORT", "false")
    return os.environ


# ---------------------------------------------------------------------------
# Marker plumbing
# ---------------------------------------------------------------------------
def pytest_configure(config):
    """Register the markers when the suite runs without pyproject.toml.

    ``--strict-markers`` is on, and a stand-alone ``pytest tests/test_x.py``
    from a source tarball may not pick the ini up.
    """
    for marker, help_text in (
        ("golden", "byte-for-byte comparison against real Perl mlst 2.35.0 output"),
        ("needs_blast", "requires a working blastn/makeblastdb"),
        ("needs_db", "requires the committed PubMLST database in db/pubmlst"),
        ("windows", "only meaningful on Windows"),
        ("gui", "drives tkinter; skipped without a display"),
        ("network", "talks to PubMLST/Pasteur/NCBI"),
        ("slow", "takes more than ~10 s"),
    ):
        config.addinivalue_line("markers", "%s: %s" % (marker, help_text))


def _has_display():
    if os.name == "nt" or sys.platform == "darwin":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def pytest_collection_modifyitems(config, items):
    """Apply the automatic skips the marker vocabulary implies."""
    no_blast = pytest.mark.skip(
        reason="no blastn found (set WMLST_TEST_BLAST_DIR or install NCBI BLAST+)"
    )
    no_db = pytest.mark.skip(reason="db/pubmlst is missing from this checkout")
    not_windows = pytest.mark.skip(reason="Windows-only behaviour")
    no_display = pytest.mark.skip(reason="no display; set DISPLAY to run the GUI tests")
    no_network = pytest.mark.skip(reason="network tests need WMLST_TEST_NETWORK=1")

    db_present = os.path.isdir(os.path.join(REPO_ROOT, "db", "pubmlst"))
    display = _has_display()
    network = os.environ.get("WMLST_TEST_NETWORK") == "1"

    for item in items:
        keywords = item.keywords
        if "needs_blast" in keywords and _BLAST_DIR is None:
            item.add_marker(no_blast)
        if "needs_db" in keywords and not db_present:
            item.add_marker(no_db)
        if "windows" in keywords and os.name != "nt":
            item.add_marker(not_windows)
        if "gui" in keywords and not display:
            item.add_marker(no_display)
        if "network" in keywords and not network:
            item.add_marker(no_network)


def pytest_report_header(config):
    """Print what the suite actually found, so a skip is never a mystery."""
    blast = _BLAST_DIR or "NOT FOUND"
    db = os.path.join(REPO_ROOT, "db", "pubmlst")
    index = os.path.join(REPO_ROOT, "db", "blast", "mlst.fa.nsq")
    return [
        "WMLST: repo=%s" % REPO_ROOT,
        "WMLST: blast=%s" % blast,
        "WMLST: db/pubmlst=%s  blast index=%s"
        % (
            "present" if os.path.isdir(db) else "MISSING",
            "built" if os.path.isfile(index) else "not built",
        ),
    ]
