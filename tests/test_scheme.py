# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Scheme / SchemeCatalog / DB-path tests. docs/ARCHITECTURE.md 4.5, 6, 5.19.

Every assertion here is measured against the shipped `db/pubmlst` snapshot
2025-12-29 and against the Perl semantics of `MLST::Scheme`.

Runnable as ``python3 -m pytest tests/test_scheme.py`` or ``python3
tests/test_scheme.py``.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst.engine import DatabaseMissingError
from wmlst.schemes import (
    HEADER_EXCLUDE_RE,
    Scheme,
    SchemeCatalog,
    perl_split,
    perl_truthy,
    resolve_blastdb,
    resolve_datadir,
    resolve_dbdir,
)

DBDIR = resolve_dbdir()
DATADIR = resolve_datadir(DBDIR)
CATALOG = SchemeCatalog(DATADIR)


# ---------------------------------------------------------------------------
# Perl emulation
# ---------------------------------------------------------------------------
def test_perl_truthy_treats_the_string_zero_as_false():
    """Section 0.1 item 5 -- eleven shipped schemes depend on this."""
    assert not perl_truthy(None)
    assert not perl_truthy("")
    assert not perl_truthy("0")
    assert perl_truthy("00")
    assert perl_truthy("0.0")
    assert perl_truthy("-")
    assert perl_truthy("1")


def test_perl_split_drops_trailing_empty_fields():
    assert perl_split("\t", "a\tb\t\t") == ["a", "b"]
    assert perl_split("\t", "a\t\tb") == ["a", "", "b"]
    assert perl_split("\t", "") == []


def test_header_exclude_regex_is_exactly_six_anchored_names():
    for name in ("ST", "mlst_clade", "clonal_complex", "species", "CC",
                 "Lineage"):
        assert HEADER_EXCLUDE_RE.match(name)
    # Case-sensitive, anchored, and deliberately NOT including MLST_cluster.
    for name in ("st", "MLST_cluster", "STx", "cc", "clonal_complexes"):
        assert not HEADER_EXCLUDE_RE.match(name), name


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
def test_resolve_dbdir_honours_mlst_dbdir():
    """Divergence D15 / test.sh:41 -- the legacy variable still works."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        old = os.environ.get("MLST_DBDIR")
        os.environ["MLST_DBDIR"] = td
        try:
            assert resolve_dbdir() == os.path.abspath(td)
        finally:
            if old is None:
                os.environ.pop("MLST_DBDIR", None)
            else:
                os.environ["MLST_DBDIR"] = old


def test_resolve_dbdir_explicit_wins_over_the_environment():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        old = os.environ.get("WMLST_DBDIR")
        os.environ["WMLST_DBDIR"] = os.path.join(td, "nope")
        try:
            assert resolve_dbdir(td) == os.path.abspath(td)
        finally:
            if old is None:
                os.environ.pop("WMLST_DBDIR", None)
            else:
                os.environ["WMLST_DBDIR"] = old


def test_resolve_dbdir_raises_with_upstream_wording():
    try:
        resolve_dbdir(os.path.join(DBDIR, "definitely-not-here"))
    except DatabaseMissingError as exc:
        assert "Database directory does not exist:" in str(exc)
    else:
        raise AssertionError("expected DatabaseMissingError")


def test_resolve_blastdb_points_at_the_shipped_index():
    stem = resolve_blastdb(DBDIR)
    assert os.path.basename(stem) == "mlst.fa"
    # The shipped tree really does hold the built index.
    assert os.path.isfile(stem + ".nin")


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------
def test_catalogue_finds_162_schemes_sorted():
    names = CATALOG.names()
    assert len(names) == 162
    assert list(names) == sorted(names)
    assert "sepidermidis" in names and "saureus" in names


def test_membership_is_case_sensitive_and_not_a_path_test():
    """Section 5.2 step 12 -- NTFS would let 'SAUREUS' through otherwise."""
    assert "saureus" in CATALOG
    assert "SAUREUS" not in CATALOG
    assert "no_such_scheme" not in CATALOG


def test_getitem_raises_keyerror_for_an_unknown_scheme():
    try:
        CATALOG["no_such_scheme"]
    except KeyError:
        pass
    else:
        raise AssertionError("expected KeyError")


def test_db_version():
    assert CATALOG.db_version == "2025-12-29"


# ---------------------------------------------------------------------------
# Scheme
# ---------------------------------------------------------------------------
def test_sepidermidis_shape():
    s = CATALOG["sepidermidis"]
    assert s.genes == ("arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL")
    assert s.num_genes == 7
    assert s.num_genotypes == 1253
    assert s.num_alleles == 601


def test_the_known_st_184_signature():
    s = CATALOG["sepidermidis"]
    calls = {"arcC": "16", "aroE": "1", "gtr": "2", "mutS": "1",
             "pyrR": "2", "tpiA": "1", "yqiL": "1"}
    sig = s.signature_of(calls)
    assert sig == "16/1/2/1/2/1/1"
    assert s.sequence_type(sig) == "184"


def test_signature_of_uses_perl_falsiness_and_gene_order():
    s = CATALOG["sepidermidis"]
    calls = {"arcC": "16", "aroE": "", "gtr": "0", "yqiL": "1"}
    assert s.signature_of(calls) == "16/-/-/-/-/-/1"


def test_sequence_type_returns_hyphen_for_an_unknown_signature():
    s = CATALOG["sepidermidis"]
    assert s.sequence_type("9999/9999/9999/9999/9999/9999/9999") == "-"


def test_aphagocytophilum_has_eight_genes_and_seven_tfa_files():
    """Section 5.19a / B9 -- the single scheme where these disagree."""
    s = CATALOG["aphagocytophilum"]
    assert s.genes == ("pheS", "glyA", "fumC", "mdh", "sucA", "dnaN", "atpA",
                       "MLST_cluster")
    assert s.num_genes == 8
    tfa = [f for f in os.listdir(s.path) if f.endswith(".tfa")]
    assert len(tfa) == 7


def test_only_aphagocytophilum_has_that_mismatch():
    odd = []
    for name in CATALOG.names():
        s = CATALOG[name]
        tfa = len([f for f in os.listdir(s.path) if f.endswith(".tfa")])
        if tfa != s.num_genes:
            odd.append(name)
    assert odd == ["aphagocytophilum"], odd


def test_null_alleles_become_hyphens_in_the_genotype_table():
    """Section 5.19c -- literal '0' values in eleven shipped profile tables."""
    s = CATALOG["sepidermidis"]
    with_null = [sig for sig in s.genotypes if sig.split("/").count("-")]
    assert with_null, "sepidermidis is one of the eleven null-allele schemes"


def test_genotype_columns_are_taken_by_index_not_by_name():
    """B12 / 5.19a -- aphagocytophilum's 8th field is clonal_complex."""
    s = CATALOG["aphagocytophilum"]
    with open(s.tab_file, encoding="utf-8", newline=None) as fh:
        fh.readline()
        first = fh.readline().rstrip("\r\n").split("\t")
    sig = "/".join(first[i] if first[i] not in ("", "0") else "-"
                   for i in range(1, 9))
    assert sig in s.genotypes
    assert s.genotypes[sig] == first[0]


def test_last_updated_is_unknown_for_every_shipped_scheme():
    """C12 / D12 -- --info's DATE column is a byte-identity surface."""
    assert all(CATALOG[n].last_updated == "Unknown" for n in CATALOG.names())


def test_lazy_loading_does_not_read_every_profile():
    """162 eager loads would cost seconds and ~90 MB for nothing."""
    s = Scheme(DATADIR, "sepidermidis")
    assert s._genotypes is None
    assert s.genes                  # touching .genes must not load the profile
    assert s._genotypes is None
    assert s.genotypes              # touching .genotypes must
    assert s._genotypes is not None


def test_every_scheme_parses():
    for name in CATALOG.names():
        s = CATALOG[name]
        assert s.num_genes >= 1, name
        assert isinstance(s.genes, tuple)


# ---------------------------------------------------------------------------
# Species map (C13 / 6.4)
# ---------------------------------------------------------------------------
def test_species_map_is_many_to_one_aware():
    assert len(CATALOG.species_rows("campylobacter")) > 1
    assert CATALOG.species_label("campylobacter") == "Campylobacter coli/jejuni"


def test_unmapped_schemes_return_none_and_are_not_invented():
    unmapped = [n for n in CATALOG.names() if not CATALOG.species_rows(n)]
    assert len(unmapped) == 80, len(unmapped)
    assert CATALOG.species_label(unmapped[0]) is None


def test_genus_only_rows_render_as_sp():
    assert CATALOG.species_label("achromobacter") == "Achromobacter sp."


def test_fields_are_stripped():
    """Several GENUS values carry a trailing space (6.4)."""
    for name in CATALOG.names():
        for genus, species in CATALOG.species_rows(name):
            assert genus == genus.strip()
            assert species == species.strip()


# ---------------------------------------------------------------------------
# SchemeInfo (6.2)
# ---------------------------------------------------------------------------
def test_info_reads_the_canonical_eight_keys():
    info = CATALOG.info("sepidermidis")
    assert info.name == "sepidermidis"
    assert info.locus == 7
    assert info.source == "pubmlst"
    assert info.api.startswith("https://rest.pubmlst.org/")
    assert info.authenticated is False
    assert info.species == "epidermidis"


def test_info_all_covers_every_scheme():
    infos = CATALOG.info_all()
    assert len(infos) == 162
    assert [i.name for i in infos] != []


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS %s" % name)
            except AssertionError as exc:
                failed += 1
                print("FAIL %s: %s" % (name, exc))
    print("%s" % ("ALL PASS" if not failed else "%d FAILED" % failed))
    sys.exit(1 if failed else 0)
