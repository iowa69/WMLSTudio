# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Database location and on-disk layout — docs/ARCHITECTURE.md sections 4.5, 6.

Two things are checked here:

* ``resolve_dbdir`` / ``resolve_blastdb`` honour the documented ladder,
  including the ``MLST_DBDIR`` drop-in compatibility variable (divergence D15).
* the shipped ``db/`` tree still satisfies every invariant the updater and the
  typing engine rely on (section 6.1): 162 scheme directories, 1,108 ``.tfa``
  files, verbatim upstream bytes with no CR and no blank lines, and a manifest
  that resolves every one of them.

Runnable as ``python3 -m pytest tests/test_dbpath.py`` or ``python3 tests/test_dbpath.py``.
"""

from __future__ import annotations

import glob
import json
import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst import updatedb as U
from wmlst.engine import DatabaseMissingError
from wmlst.version import BUNDLED_DB_VERSION

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(REPO, "db")
SCHEMES = os.path.join(DB, "pubmlst")

schemes_module = pytest.importorskip("wmlst.schemes",
                                     reason="wmlst/schemes.py is not present yet")

SCHEME_DIRS = sorted(n for n in os.listdir(SCHEMES)
                     if os.path.isdir(os.path.join(SCHEMES, n))
                     and not n.startswith("."))


# ---------------------------------------------------------------------------
# 4.5 resolve_dbdir
# ---------------------------------------------------------------------------

def _clear_env(monkeypatch):
    for name in ("WMLST_DBDIR", "MLST_DBDIR"):
        monkeypatch.delenv(name, raising=False)


def test_an_explicit_directory_wins(monkeypatch):
    _clear_env(monkeypatch)
    with tempfile.TemporaryDirectory() as tmp:
        monkeypatch.setenv("WMLST_DBDIR", DB)
        assert os.path.realpath(schemes_module.resolve_dbdir(tmp)) == \
            os.path.realpath(tmp)


def test_wmlst_dbdir_is_honoured(monkeypatch):
    _clear_env(monkeypatch)
    with tempfile.TemporaryDirectory() as tmp:
        monkeypatch.setenv("WMLST_DBDIR", tmp)
        assert os.path.realpath(schemes_module.resolve_dbdir()) == os.path.realpath(tmp)


def test_mlst_dbdir_is_honoured_for_drop_in_compatibility(monkeypatch):
    """D15: upstream's test.sh:41 sets MLST_DBDIR."""
    _clear_env(monkeypatch)
    with tempfile.TemporaryDirectory() as tmp:
        monkeypatch.setenv("MLST_DBDIR", tmp)
        assert os.path.realpath(schemes_module.resolve_dbdir()) == os.path.realpath(tmp)


def test_wmlst_dbdir_beats_mlst_dbdir(monkeypatch):
    _clear_env(monkeypatch)
    with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
        monkeypatch.setenv("WMLST_DBDIR", first)
        monkeypatch.setenv("MLST_DBDIR", second)
        assert os.path.realpath(schemes_module.resolve_dbdir()) == \
            os.path.realpath(first)


def test_a_missing_directory_raises_with_upstream_wording(monkeypatch):
    """4.5: 'Database directory does not exist: <path>' mirrors bin/mlst:78."""
    _clear_env(monkeypatch)
    missing = os.path.join(REPO, "no-such-database-directory")
    with pytest.raises(DatabaseMissingError) as excinfo:
        schemes_module.resolve_dbdir(missing)
    assert "Database directory does not exist" in str(excinfo.value)
    assert missing in str(excinfo.value)


def test_the_source_checkout_is_found_without_any_environment(monkeypatch):
    _clear_env(monkeypatch)
    resolved = schemes_module.resolve_dbdir()
    assert os.path.isdir(os.path.join(resolved, "pubmlst"))


def test_resolve_blastdb_points_at_the_shipped_index():
    stem = schemes_module.resolve_blastdb(DB)
    assert os.path.basename(stem) == "mlst.fa"
    assert stem.endswith(os.path.join("blast", "mlst.fa"))


def test_resolve_blastdb_never_returns_a_relative_path():
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "pubmlst"))
        stem = schemes_module.resolve_blastdb(tmp)
        assert os.path.isabs(stem) or stem.startswith(tmp)


# ---------------------------------------------------------------------------
# 6.1 the shipped tree
# ---------------------------------------------------------------------------

def test_the_tree_has_162_schemes_and_1108_allele_files():
    assert len(SCHEME_DIRS) == 162
    assert len(glob.glob(os.path.join(SCHEMES, "*", "*.tfa"))) == 1108


def test_every_scheme_has_a_profile_table_and_an_info_file():
    for name in SCHEME_DIRS:
        assert os.path.isfile(os.path.join(SCHEMES, name, name + ".txt")), name
        assert os.path.isfile(os.path.join(SCHEMES, name, name + "_info.json")), name


def test_no_scheme_ships_a_database_version_file():
    """D12/4.5: --info's DATE column must stay 'Unknown' for all 162."""
    assert glob.glob(os.path.join(SCHEMES, "*", "database_version.txt")) == []


def test_version_txt_matches_the_pinned_bundle_version():
    assert U.db_version(DB) == BUNDLED_DB_VERSION


def test_the_manifest_resolves_every_directory_on_disk():
    refs = U.load_manifest(DB)
    assert sorted(r.name for r in refs) == SCHEME_DIRS
    assert all(r.resolved for r in refs)
    assert all(r.api.startswith(("https://rest.pubmlst.org/",
                                 "https://bigsdb.pasteur.fr/api/")) for r in refs)


def test_info_locus_count_matches_the_allele_file_count():
    for name in SCHEME_DIRS:
        info = json.loads(open(os.path.join(SCHEMES, name, name + "_info.json"),
                               "rb").read())
        tfas = glob.glob(os.path.join(SCHEMES, name, "*.tfa"))
        assert info["locus"] == len(tfas), name


def test_every_payload_file_is_free_of_cr_and_blank_lines():
    """6.1: verbatim upstream bytes; the corpus is hygienic and must stay so."""
    for path in sorted(glob.glob(os.path.join(SCHEMES, "*", "*.tfa")) +
                       glob.glob(os.path.join(SCHEMES, "*", "*.txt"))):
        with open(path, "rb") as fh:
            data = fh.read()
        assert b"\r" not in data, path
        assert b"\n\n" not in data, path
        assert data.endswith(b"\n"), path


def test_every_locus_name_is_a_safe_windows_filename():
    """7.6: a ':' in a locus name would create an NTFS alternate data stream."""
    for path in glob.glob(os.path.join(SCHEMES, "*", "*.tfa")):
        U._validate_locus_name(os.path.basename(path)[:-4])


def test_the_only_unparseable_loci_are_the_three_documented_ones():
    """5.19(b): the loci behind halobacteria.EF-2_N, rpoB'_N and mhominis_3.p120'_N.

    The spec names them by seqid (``<locus>_<allele>``); on disk the files are
    ``EF-2.tfa``, ``rpoB'.tfa`` and ``p120'.tfa``.
    """
    odd = []
    for path in sorted(glob.glob(os.path.join(SCHEMES, "*", "*.tfa"))):
        locus = os.path.basename(path)[:-4]
        if U._unparseable_loci([locus]):
            odd.append("%s.%s" % (os.path.basename(os.path.dirname(path)), locus))
    assert sorted(odd) == ["halobacteria.EF-2", "halobacteria.rpoB'",
                           "mhominis_3.p120'"]


def test_no_scheme_has_a_duplicate_seqid():
    for name in SCHEME_DIRS:
        seen = set()
        for path in glob.glob(os.path.join(SCHEMES, name, "*.tfa")):
            with open(path, "rb") as fh:
                for line in fh:
                    if line.startswith(b">"):
                        seqid = line[1:].strip()
                        assert seqid not in seen, (name, seqid)
                        seen.add(seqid)


def test_no_defline_contains_a_pipe():
    """A '|' in a defline would change how BLAST parses the seqid."""
    for path in glob.glob(os.path.join(SCHEMES, "*", "*.tfa")):
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    assert b"|" not in line, path


def test_profile_headers_start_with_ST():
    for name in SCHEME_DIRS:
        with open(os.path.join(SCHEMES, name, name + ".txt"), "rb") as fh:
            header = fh.readline().rstrip(b"\n").split(b"\t")
        assert header[0] == b"ST", name
        assert len(header) >= 2, name


def test_mgenitalium_has_a_header_only_profile_table():
    """6.2: zero profiles is legal and must not be repaired."""
    path = os.path.join(SCHEMES, "mgenitalium", "mgenitalium.txt")
    with open(path, "rb") as fh:
        data = fh.read()
    assert data.count(b"\n") == 1
    assert len(data) == 59
    assert U.scheme_info(DB, "mgenitalium")["last_updated"] == U.NO_VERSION


def test_aphagocytophilum_has_eight_header_genes_and_seven_files():
    """5.19(a): the one scheme where len(genes) != count(*.tfa). Never 'fix' it."""
    with open(os.path.join(SCHEMES, "aphagocytophilum",
                           "aphagocytophilum.txt"), "rb") as fh:
        header = fh.readline().decode().rstrip("\n").split("\t")
    loci = [c for c in header if c not in U.PROFILE_NON_LOCUS]
    assert len(loci) == 8 and loci[-1] == "MLST_cluster"
    assert len(glob.glob(os.path.join(SCHEMES, "aphagocytophilum", "*.tfa"))) == 7


# ---------------------------------------------------------------------------
# 6.4 scheme_species_map.tab
# ---------------------------------------------------------------------------

def test_species_map_shape():
    path = os.path.join(DB, "scheme_species_map.tab")
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert len(lines) == 133
    assert lines[0].split("\t") == ["#SCHEME", "GENUS", "SPECIES"]
    rows = [line.split("\t") for line in lines[1:]]
    assert len(rows) == 132
    assert len({r[0] for r in rows}) == 129
    assert any(r[1] != r[1].strip() for r in rows), "some GENUS values are padded"


def test_species_map_is_never_the_source_of_scheme_names():
    """6.4: 80 directories are unmapped and 47 rows name vanished directories."""
    path = os.path.join(DB, "scheme_species_map.tab")
    with open(path, encoding="utf-8") as fh:
        mapped = {line.split("\t")[0].strip() for line in fh.read().splitlines()[1:]}
    unmapped = set(SCHEME_DIRS) - mapped
    stale = mapped - set(SCHEME_DIRS)
    assert len(unmapped) == 80
    assert len(stale) == 47


# ---------------------------------------------------------------------------
# 4.9 index staleness
# ---------------------------------------------------------------------------

def test_index_is_stale_is_false_for_a_freshly_built_tree():
    if not os.path.isfile(os.path.join(DB, "blast", "mlst.fa.nsq")):
        pytest.skip("db/blast has not been built in this checkout")
    assert U.index_is_stale(DB) is False


def test_index_is_stale_when_the_index_is_absent():
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "pubmlst", "tiny"))
        assert U.index_is_stale(tmp) is True


def test_index_is_stale_on_an_empty_tree_rather_than_raising():
    with tempfile.TemporaryDirectory() as tmp:
        assert U.index_is_stale(tmp) is True


if __name__ == "__main__":
    raise SystemExit(pytest.main([os.path.abspath(__file__), "-q"]))
