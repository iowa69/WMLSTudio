# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""``db/scheme_refs.tsv`` and :mod:`wmlst.schemerefs`.

Two halves:

* the **data** -- the shipped table must cover all 162 bundled schemes, agree
  with each scheme's ``_info.json`` on locus count and source, carry a
  clickable database URL for every row, and never carry a malformed PMID or a
  citation without one. A blank genus/species/citation is legal and expected.
* the **loader** -- caching, a missing file degrading to empty, short rows,
  CRLF, duplicate rows, and the ``organism_label`` contract.

Runnable as ``python3 -m pytest tests/test_schemerefs.py`` or as
``python3 tests/test_schemerefs.py``.
"""

from __future__ import annotations

import json
import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst import schemerefs as R

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(REPO, "db")
SCHEMES_DIR = os.path.join(DB, "pubmlst")
TABLE = os.path.join(DB, R.REFS_FILENAME)

#: Schemes the brief singles out: the ESKAPE set plus the high-volume
#: public-health schemes. Each of these must carry a verified citation.
PRIORITY = (
    "saureus", "abaumannii", "abaumannii_2", "ecoli_achtman_4", "efaecium",
    "klebsiella", "paeruginosa", "ecloacae", "senterica_achtman_2",
    "listeria_2", "campylobacter", "spneumoniae", "cdifficile",
)


def bundled_schemes():
    """Every scheme directory under ``db/pubmlst``."""
    return sorted(
        name for name in os.listdir(SCHEMES_DIR)
        if os.path.isdir(os.path.join(SCHEMES_DIR, name))
        and name != "__pycache__"
    )


@pytest.fixture(scope="module")
def refs():
    R.clear_cache()
    table = R.load_refs(DB)
    assert table, "db/%s parsed to nothing" % R.REFS_FILENAME
    return table


# ---------------------------------------------------------------------------
# The shipped table
# ---------------------------------------------------------------------------
def test_table_exists_and_is_lf_utf8():
    assert os.path.isfile(TABLE), "db/%s is missing" % R.REFS_FILENAME
    with open(TABLE, "rb") as fh:
        raw = fh.read()
    assert b"\r" not in raw, "scheme_refs.tsv contains CR; LF only (section 0)"
    raw.decode("utf-8")  # raises on a bad byte
    assert raw.endswith(b"\n")


def test_header_names_every_column_in_order():
    with open(TABLE, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
    assert header.startswith("#")
    assert tuple(header[1:].split("\t")) == R.COLUMNS


def test_covers_every_bundled_scheme_exactly_once(refs):
    names = bundled_schemes()
    assert len(names) == 162, "expected 162 bundled schemes, found %d" % len(names)
    assert sorted(refs) == names


def test_rows_are_in_scheme_order():
    rows = []
    with open(TABLE, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.split("\t")[0])
    assert rows == sorted(rows), "rows must be sorted by scheme name"


def test_every_row_has_a_description_and_field_count(refs):
    with open(TABLE, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            assert len(cols) == len(R.COLUMNS), (
                "line %d has %d columns, expected %d"
                % (lineno, len(cols), len(R.COLUMNS)))
    for ref in refs.values():
        assert ref.description, "%s has no description" % ref.scheme
        assert "MLST" in ref.description


def test_n_loci_and_source_agree_with_info_json(refs):
    for name, ref in sorted(refs.items()):
        path = os.path.join(SCHEMES_DIR, name, name + "_info.json")
        with open(path, encoding="utf-8") as fh:
            info = json.load(fh)
        assert ref.n_loci == info["locus"], name
        assert ref.source == info["source"], name
        tfa = [f for f in os.listdir(os.path.join(SCHEMES_DIR, name))
               if f.endswith(".tfa")]
        assert ref.n_loci == len(tfa), "%s: n_loci != .tfa count" % name


def test_every_row_has_a_clickable_database_url(refs):
    for name, ref in sorted(refs.items()):
        assert ref.database_url.startswith("https://"), name
        if ref.source == "pasteur":
            assert "bigsdb.pasteur.fr" in ref.database_url, name
        else:
            assert "pubmlst.org" in ref.database_url, name
        assert "page=schemeInfo" in ref.database_url, name


def test_genus_and_species_are_clean(refs):
    """No trailing space, no lower-case genus, no species duplicating genus.

    ``scheme_species_map.tab`` ships ``"Enterococcus "`` and
    ``"Pseudomonas "`` with a trailing space and ``afumigatus`` as a species;
    none of that may survive into this table.
    """
    for name, ref in sorted(refs.items()):
        for field in (ref.genus, ref.species):
            assert field == field.strip(), name
            assert "  " not in field, name
        if ref.genus:
            assert ref.genus[0].isupper(), name
        if ref.species:
            assert ref.genus, "%s has a species but no genus" % name
            # A bare lower-case epithet: no "afumigatus", no scheme name.
            assert re.fullmatch(r"[a-z][a-z-]+", ref.species), (
                "%s: %r is not a species epithet" % (name, ref.species))
            assert ref.species != name, name


def test_species_is_never_guessed_where_the_scheme_spans_a_genus(refs):
    """Genus-wide and above-genus schemes carry a blank species on purpose."""
    for name in ("neisseria", "borrelia", "vibrio", "leptospira", "brucella",
                 "campylobacter", "aeromonas", "achromobacter"):
        assert refs[name].species == "", name
    for name in ("chlamydiales", "halobacteria", "mycobacteria_2",
                 "llactis_phage"):
        assert refs[name].genus == "", name
        assert refs[name].organism == "", name


def test_pubmed_ids_are_well_formed_and_paired_with_a_citation(refs):
    for name, ref in sorted(refs.items()):
        if ref.pubmed_id:
            assert re.fullmatch(r"[1-9][0-9]{6,8}", ref.pubmed_id), (
                "%s: %r is not a PMID" % (name, ref.pubmed_id))
            assert ref.citation, "%s has a PMID but no citation" % name
        if ref.citation:
            assert ref.pubmed_id, "%s has a citation but no PMID" % name
            # "Author A, ... Title. Journal Year;vol:pages."
            assert re.search(r"\b(19|20)\d{2}\b", ref.citation), name
            assert ref.citation.endswith("."), name


def test_pubmed_ids_are_unique_per_publication(refs):
    """A PMID may be shared only by schemes that really share a paper."""
    shared = {}
    for name, ref in sorted(refs.items()):
        if ref.pubmed_id:
            shared.setdefault(ref.pubmed_id, []).append(name)
    for pmid, names in sorted(shared.items()):
        if len(names) > 1:
            assert set(names) in (
                # the same Achtman 7-locus scheme, mirrored in two databases
                {"salmonella", "senterica_achtman_2"},
                # cdiphtheriae is an alias of diphtheria_3 (schemes.manifest)
                {"cdiphtheriae", "diphtheria_3"},
                # one paper describes both the MLST and the expanded eMLST
                {"mhominis", "mhominis_3"},
            ), "PMID %s reused by %s" % (pmid, names)


def test_priority_schemes_all_carry_a_verified_citation(refs):
    missing = [s for s in PRIORITY if not refs[s].has_citation]
    assert not missing, "priority schemes without a citation: %s" % missing


def test_known_rows_are_exactly_right(refs):
    sa = refs["saureus"]
    assert (sa.genus, sa.species) == ("Staphylococcus", "aureus")
    assert sa.organism == "Staphylococcus aureus"
    assert sa.n_loci == 7
    assert sa.pubmed_id == "10698988"
    assert "Enright" in sa.citation

    # The two A. baumannii schemes are different papers, not the same one.
    assert refs["abaumannii"].pubmed_id != refs["abaumannii_2"].pubmed_id
    assert "Oxford" in refs["abaumannii"].description
    assert "Pasteur" in refs["abaumannii_2"].description

    # ecoli_achtman_4 is absent from scheme_species_map.tab; it must be here.
    ec = refs["ecoli_achtman_4"]
    assert (ec.genus, ec.species) == ("Escherichia", "coli")

    # Trailing-space quirks of the legacy map must not reappear.
    assert refs["efaecium"].genus == "Enterococcus"
    assert refs["paeruginosa"].genus == "Pseudomonas"

    # hsuis is Helicobacter suis; the legacy map says "Haematopinus".
    assert refs["hsuis"].genus == "Helicobacter"

    # A genus-wide scheme renders with "spp.".
    assert refs["neisseria"].organism == "Neisseria spp."


def test_citation_counts_are_reported(refs):
    """Not an assertion about a number; a guard that the split is sane."""
    cited = R.cited_schemes(DB)
    assert 25 <= len(cited) <= len(refs)
    assert all(r.has_citation and r.pubmed_id for r in cited)
    sys.stderr.write(
        "\nscheme_refs.tsv: %d/%d rows carry a verified citation, %d blank\n"
        % (len(cited), len(refs), len(refs) - len(cited)))


# ---------------------------------------------------------------------------
# The loader
# ---------------------------------------------------------------------------
def test_lookup_and_organism_label_against_the_real_db(monkeypatch):
    monkeypatch.setenv("WMLST_DBDIR", DB)
    R.clear_cache()
    assert R.lookup("saureus").scheme == "saureus"
    assert R.organism_label("saureus") == "Staphylococcus aureus"
    assert R.organism_label("chlamydiales") == ""
    assert R.lookup("no_such_scheme") is None
    assert R.organism_label("no_such_scheme") == ""
    assert R.organism_label("") == ""
    assert R.describe("saureus").startswith("Staphylococcus aureus")
    assert R.describe("no_such_scheme") == "no_such_scheme"


def test_lookup_strips_whitespace(refs):
    assert R.lookup("  saureus  ", DB) == refs["saureus"]
    assert R.lookup("saureus", DB) is R.lookup("  saureus  ", DB)


def test_missing_file_degrades_to_empty(tmp_path):
    R.clear_cache()
    empty = str(tmp_path)
    assert R.load_refs(empty) == {}
    assert R.lookup("saureus", empty) is None
    assert R.organism_label("saureus", empty) == ""
    assert R.cited_schemes(empty) == []


def test_unreadable_file_degrades_to_empty(tmp_path):
    """A directory where the table should be is not a crash."""
    os.mkdir(os.path.join(str(tmp_path), R.REFS_FILENAME))
    R.clear_cache()
    assert R.load_refs(str(tmp_path)) == {}


def test_result_is_cached_per_directory(tmp_path):
    path = os.path.join(str(tmp_path), R.REFS_FILENAME)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("#" + "\t".join(R.COLUMNS) + "\n")
        fh.write("x\tGenus\tspecies\tGenus species MLST\t7\tpubmlst\t"
                 "https://pubmlst.org/\t\t\n")
    R.clear_cache()
    first = R.load_refs(str(tmp_path))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("")
    assert R.load_refs(str(tmp_path)) is first, "second load re-read the file"
    R.clear_cache()
    assert R.load_refs(str(tmp_path)) == {}


def test_parse_refs_tolerates_junk():
    text = (
        "#" + "\t".join(R.COLUMNS) + "\r\n"
        "\r\n"
        "   \n"
        "# a comment\n"
        "short\tGenus\tspecies\tGenus species MLST\r\n"   # 4 of 9 columns
        "wide\tG\ts\tG s MLST\t7\tpubmlst\thttps://x/\t1\tA. 2000.\textra\n"
        "\tOrphan\t\t\t\t\t\t\t\n"                        # no scheme name
        "dup\tA\tb\tA b MLST\t3\tpubmlst\thttps://x/\t\t\n"
        "dup\tZ\tz\tZ z MLST\t9\tpubmlst\thttps://x/\t\t\n"
        "bad\tG\ts\tG s MLST\tnot-a-number\tpubmlst\thttps://x/\t\t\n"
    )
    table = R.parse_refs(text)
    assert sorted(table) == ["bad", "dup", "short", "wide"]
    assert table["short"].n_loci == 0
    assert table["short"].source == ""
    assert table["short"].organism == "Genus species"
    assert table["wide"].citation == "A. 2000."
    assert table["dup"].genus == "A", "the first row for a scheme must win"
    assert table["bad"].n_loci == 0


def test_parse_refs_of_nothing_is_empty():
    assert R.parse_refs("") == {}
    assert R.parse_refs("#only\ta\theader\n") == {}


def test_schemeref_helpers():
    ref = R.SchemeRef(scheme="x", genus="Genus", species="species",
                      pubmed_id="12345678")
    assert ref.organism == "Genus species"
    assert ref.pubmed_url == "https://pubmed.ncbi.nlm.nih.gov/12345678/"
    assert ref.has_citation is False
    assert R.SchemeRef(scheme="x", genus="Genus").organism == "Genus spp."
    assert R.SchemeRef(scheme="x").organism == ""
    assert R.SchemeRef(scheme="x").pubmed_url == ""


def test_module_is_a_stdlib_only_leaf():
    """No gui, no engine internals, no third-party import (section 13.1)."""
    import ast

    src = os.path.join(REPO, "wmlst", "schemerefs.py")
    tree = ast.parse(open(src, encoding="utf-8").read(), filename=src)
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert "wmlst.gui" not in imported
    assert "wmlst.engine" not in imported
    assert not any(m.startswith("wmlst.") and m != "wmlst.schemes"
                   for m in imported)
    assert "subprocess" not in imported


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
