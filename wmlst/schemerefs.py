# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Per-scheme provenance: organism, scheme description and literature reference.

WMLST calls a scheme; a reader of the result needs to know what that scheme
*is*. ``db/scheme_refs.tsv`` answers four questions for every one of the 162
bundled schemes:

* which **genus** and **species** the scheme types,
* the scheme's **full name / description** (``Acinetobacter baumannii MLST
  (Oxford)`` rather than the bare directory name ``abaumannii``),
* a **database_url** the reader can click through to the authoritative
  PubMLST / Institut Pasteur record,
* and, where the primary MLST publication could be **verified**, its PubMed
  identifier and a short citation.

The table supersedes ``db/scheme_species_map.tab`` as a *presentation* source.
That file stays exactly as upstream ships it because ``report.py`` renders it
and ``tests/test_dbpath.py`` pins its bytes; this module never touches it.

Two rules govern the contents, and the loader below assumes both:

1. **A blank cell means "not established", never "none".** Most schemes have no
   ``pubmed_id``: PubMLST's REST API exposes no per-scheme citation, so every
   reference in the table was verified by hand against PubMed. A wrong citation
   on a diagnostic report is far worse than an absent one.
2. **A missing or unreadable table degrades to empty.** :func:`lookup` returns
   ``None`` and :func:`organism_label` returns ``""``; nothing here raises, so a
   results view can call it unconditionally.

This module is deliberately a leaf: stdlib only, no ``engine`` internals, no
``gui``. It reads one small text file and caches the parse.

Public API::

    SchemeRef          frozen dataclass, one row of the table
    load_refs(dbdir)   -> Mapping[str, SchemeRef], cached per directory
    lookup(scheme)     -> SchemeRef | None
    organism_label(s)  -> "Staphylococcus aureus" | "Neisseria spp." | ""
    clear_cache()      forget every parse (tests, and after an updatedb run)
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional

#: Basename of the table inside the database directory.
REFS_FILENAME = "scheme_refs.tsv"

#: Column order of ``db/scheme_refs.tsv``. The header line repeats it, prefixed
#: with ``#``; the loader trusts this constant rather than the header so that a
#: hand-edited file with a mangled comment still parses.
COLUMNS = (
    "scheme",
    "genus",
    "species",
    "description",
    "n_loci",
    "source",
    "database_url",
    "pubmed_id",
    "citation",
)

__all__ = [
    "COLUMNS",
    "REFS_FILENAME",
    "SchemeRef",
    "clear_cache",
    "load_refs",
    "lookup",
    "organism_label",
]


@dataclass(frozen=True)
class SchemeRef:
    """One row of ``db/scheme_refs.tsv``.

    Every string field is already stripped. ``genus``, ``species``,
    ``pubmed_id`` and ``citation`` are ``""`` when the value could not be
    established from a reliable source -- that is a deliberate, meaningful
    value and must be rendered as "unknown", never guessed at.
    """

    scheme: str
    genus: str = ""
    species: str = ""
    description: str = ""
    n_loci: int = 0
    source: str = ""
    database_url: str = ""
    pubmed_id: str = ""
    citation: str = ""

    # -- convenience, for the results view ----------------------------------
    @property
    def organism(self) -> str:
        """``"Genus species"``, ``"Genus spp."`` or ``""``.

        A genus-wide scheme (``neisseria``, ``borrelia``, ``vibrio``, ...)
        carries no species on purpose and renders as ``Genus spp.``.
        """
        if self.genus and self.species:
            return "%s %s" % (self.genus, self.species)
        if self.genus:
            return "%s spp." % self.genus
        return ""

    @property
    def pubmed_url(self) -> str:
        """Canonical PubMed URL, or ``""`` when no PMID was verified."""
        if not self.pubmed_id:
            return ""
        return "https://pubmed.ncbi.nlm.nih.gov/%s/" % self.pubmed_id

    @property
    def has_citation(self) -> bool:
        """True when a primary publication was verified for this scheme."""
        return bool(self.citation)


# ---------------------------------------------------------------------------
# Database directory
# ---------------------------------------------------------------------------
def _fallback_dbdir() -> str:
    """``<repo>/db`` next to the installed package -- the last-resort guess."""
    return os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "db")


def default_dbdir() -> str:
    """Locate the database directory without ever raising.

    Delegates to :func:`wmlst.schemes.resolve_dbdir` so there is one ladder in
    the product, and falls back to the source-checkout / frozen locations when
    that import or that call fails (no database installed, for instance).
    """
    try:
        from wmlst.schemes import resolve_dbdir
    except Exception:
        resolve_dbdir = None  # type: ignore[assignment]
    if resolve_dbdir is not None:
        try:
            return resolve_dbdir()
        except Exception:
            pass
    for value in (os.environ.get("WMLST_DBDIR"), os.environ.get("MLST_DBDIR")):
        if value:
            return os.path.abspath(os.path.expanduser(value))
    base = getattr(sys, "_MEIPASS", None)
    if base:
        frozen = os.path.join(base, "wmlst_db")
        if os.path.isdir(frozen):
            return frozen
    return _fallback_dbdir()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def _to_int(text: str) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def parse_refs(text: str) -> Dict[str, SchemeRef]:
    """Parse the whole table from a string. Never raises.

    Blank lines and ``#`` comment lines (including the header) are skipped, as
    are rows with an empty first column. A short row is padded with ``""``; a
    long row keeps only the known columns. The first row for a scheme wins, so
    a duplicate cannot silently overwrite a curated entry.
    """
    table = {}  # type: Dict[str, SchemeRef]
    for raw in text.splitlines():
        line = raw.rstrip("\r\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cols = [c.strip() for c in line.split("\t")]
        if len(cols) < len(COLUMNS):
            cols += [""] * (len(COLUMNS) - len(cols))
        name = cols[0]
        if not name or name in table:
            continue
        table[name] = SchemeRef(
            scheme=name,
            genus=cols[1],
            species=cols[2],
            description=cols[3],
            n_loci=_to_int(cols[4]),
            source=cols[5],
            database_url=cols[6],
            pubmed_id=cols[7],
            citation=cols[8],
        )
    return table


#: dbdir (absolute, normalised) -> parsed table. Populated on first use.
_CACHE = {}  # type: Dict[str, Dict[str, SchemeRef]]


def clear_cache() -> None:
    """Drop every cached parse. Call after the updater rewrites the table."""
    _CACHE.clear()


def _key(dbdir: Optional[str]) -> str:
    return os.path.normcase(os.path.abspath(dbdir or default_dbdir()))


def load_refs(dbdir: Optional[str] = None) -> Mapping[str, SchemeRef]:
    """Return ``{scheme: SchemeRef}`` for ``<dbdir>/scheme_refs.tsv``.

    Cached per directory. A missing, unreadable or empty file yields an empty
    mapping -- this function does not raise, whatever is on disk.
    """
    key = _key(dbdir)
    cached = _CACHE.get(key)
    if cached is not None:
        return cached
    text = ""
    try:
        with open(os.path.join(key, REFS_FILENAME), encoding="utf-8",
                  errors="replace", newline=None) as fh:
            text = fh.read()
    except OSError:
        text = ""
    table = parse_refs(text)
    _CACHE[key] = table
    return table


def lookup(scheme: str, dbdir: Optional[str] = None) -> Optional[SchemeRef]:
    """The :class:`SchemeRef` for ``scheme``, or ``None`` when unlisted."""
    if not scheme:
        return None
    return load_refs(dbdir).get(scheme.strip())


def organism_label(scheme: str, dbdir: Optional[str] = None) -> str:
    """``"Staphylococcus aureus"``, ``"Neisseria spp."`` or ``""``.

    ``""`` means "this scheme's organism is not established" -- either the
    scheme is unlisted, or its genus was deliberately left blank because the
    scheme spans a rank above genus (``chlamydiales``, ``halobacteria``).
    """
    ref = lookup(scheme, dbdir)
    return ref.organism if ref is not None else ""


def describe(scheme: str, dbdir: Optional[str] = None) -> str:
    """Best available one-line name for ``scheme``.

    The curated description when there is one, else the organism label, else
    the scheme name itself -- so a caller can print this unconditionally.
    """
    ref = lookup(scheme, dbdir)
    if ref is None:
        return scheme or ""
    return ref.description or ref.organism or ref.scheme


def cited_schemes(dbdir: Optional[str] = None) -> List[SchemeRef]:
    """Every scheme with a verified primary publication, in scheme order."""
    refs = load_refs(dbdir)
    return [refs[name] for name in sorted(refs) if refs[name].has_citation]
