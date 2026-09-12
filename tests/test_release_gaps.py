# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Release-gap regressions found by the integration check.

Two independent defects, both invisible to the rest of the suite because both
only bite outside a source checkout or outside the happy path:

* **Shipping.** ``db/scheme_refs.tsv`` is loaded by :mod:`wmlst.schemerefs`,
  which degrades *silently* to an empty table when the file is absent. It was
  declared in neither ``pyproject.toml`` nor ``MANIFEST.in``, so every wheel
  and sdist shipped without it and every organism label came back blank with
  no error anywhere. The frozen bundle was unaffected (it takes the whole
  ``db`` tree), which is exactly why nothing caught it.

* **Shutdown error surfacing.** ``find_blast`` walks a ladder of candidates and
  rejects any candidate whose ``probe_version`` raises. Once
  ``blastbin.shutdown`` latches, ``run_tool`` refuses to launch and raises
  ``Cancelled`` -- which the ladder swallowed as "candidate is bad", tried every
  remaining entry, and finally raised ``BlastNotFoundError`` whose
  ``user_message`` tells the user to go and install BLAST+. Closing the window
  mid-run is not a broken installation and must not be reported as one.
"""

from __future__ import annotations

import os
import re

import pytest

from wmlst import blastbin, engine, schemerefs

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

#: Every table shipped loose in ``db/``. They are declared in two places and
#: have to stay in lock-step, which is the only reason this list exists.
DB_TABLES = ("VERSION.txt", "scheme_species_map.tab", "scheme_refs.tsv",
             "schemes.manifest.tsv")


def _package_data() -> str:
    """Raw body of ``[tool.setuptools.package-data]``'s ``wmlst_db = [...]``.

    Sliced out of the text rather than parsed: ``tomllib`` landed in 3.11 and
    this suite runs on the 3.9 floor, and a TOML parser is the one dependency
    that cannot be back-filled under the stdlib-only rule.
    """
    with open(os.path.join(REPO, "pyproject.toml"), encoding="utf-8") as fh:
        text = fh.read()
    head = text.index("[tool.setuptools.package-data]")
    match = re.search(r"^wmlst_db\s*=\s*\[(.*?)\]", text[head:],
                      re.DOTALL | re.MULTILINE)
    assert match, "pyproject.toml has no package-data wmlst_db list"
    return match.group(1)


def _manifest() -> str:
    with open(os.path.join(REPO, "MANIFEST.in"), encoding="utf-8") as fh:
        return fh.read()


def _spec() -> str:
    with open(os.path.join(REPO, "packaging", "wmlst.spec"), encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# Shipping
# ---------------------------------------------------------------------------
def test_scheme_refs_is_declared_as_wheel_package_data():
    """Without this line the wheel has no scheme_refs.tsv and nothing warns."""
    assert '"%s"' % (schemerefs.REFS_FILENAME,) in _package_data(), (
        "%s missing from [tool.setuptools.package-data].wmlst_db; the wheel "
        "would ship without it and schemerefs would silently degrade to blank "
        "genus/species/citation for all 162 schemes" % (schemerefs.REFS_FILENAME,))


def test_scheme_refs_is_declared_in_the_sdist_manifest():
    assert "include db/%s" % (schemerefs.REFS_FILENAME,) in _manifest()


def test_every_loose_db_table_is_declared_in_both_places():
    """The sibling tables in db/ must not drift apart again."""
    data, manifest = _package_data(), _manifest()
    for name in DB_TABLES:
        assert '"%s"' % (name,) in data, "%s not in package-data" % (name,)
        assert "include db/%s" % (name,) in manifest, "%s not in MANIFEST" % (name,)


def test_the_shipped_table_actually_exists_where_the_loader_looks():
    """A declaration is worthless if the file moved."""
    assert os.path.isfile(os.path.join(REPO, "db", schemerefs.REFS_FILENAME))
    assert len(schemerefs.load_refs(os.path.join(REPO, "db"))) == 162


def test_the_new_leaf_modules_are_pinned_as_hidden_imports():
    """perf/schemerefs are reached by name from the GUI; PyInstaller needs them."""
    spec = _spec()
    for name in ("wmlst.perf", "wmlst.schemerefs"):
        assert '"%s"' % (name,) in spec, "%s missing from the spec HIDDEN list" % (name,)


# ---------------------------------------------------------------------------
# Shutdown error surfacing
# ---------------------------------------------------------------------------
@pytest.fixture
def latched():
    """blastbin in its shut-down state, always reset afterwards."""
    blastbin.shutdown(0.5)
    assert blastbin.is_shutting_down()
    try:
        yield
    finally:
        blastbin.reset_shutdown()
    assert not blastbin.is_shutting_down()


def test_find_blast_reports_cancelled_not_missing_blast(latched):
    """The regression: shutdown must never look like a broken installation."""
    with pytest.raises(engine.Cancelled):
        blastbin.find_blast()


def test_find_blast_with_an_explicit_path_also_reports_cancelled(latched):
    """The explicit branch runs before the ladder and must latch too."""
    with pytest.raises(engine.Cancelled):
        blastbin.find_blast(os.path.join(REPO, "no", "such", "blastn"))


def test_a_cancelled_probe_does_not_silently_reject_the_candidate(latched):
    """_tools_from_blastn must propagate Cancelled, not return None."""
    with pytest.raises(engine.Cancelled):
        blastbin._tools_from_blastn(os.path.join(REPO, "blastn"), "explicit")


def test_cancelled_is_not_a_blast_not_found_error():
    """The two must stay distinguishable: the GUI branches on the difference."""
    assert not issubclass(engine.Cancelled, engine.BlastNotFoundError)
    assert not issubclass(engine.BlastNotFoundError, engine.Cancelled)


def test_a_genuinely_bad_candidate_is_still_rejected_quietly():
    """The broad reject path is intact when we are NOT shutting down."""
    assert not blastbin.is_shutting_down()
    assert blastbin._tools_from_blastn(
        os.path.join(REPO, "definitely", "not", "blastn"), "explicit") is None


def test_find_blast_works_again_after_reset_shutdown():
    """reset_shutdown un-latches the ladder (the GUI's aborted-close path)."""
    blastbin.shutdown(0.5)
    blastbin.reset_shutdown()
    assert not blastbin.is_shutting_down()
    try:
        blastbin.find_blast()
    except Exception as exc:  # no BLAST+ on this box is fine; Cancelled is not
        assert not isinstance(exc, engine.Cancelled), \
            "still latched after reset_shutdown()"
