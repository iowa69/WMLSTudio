# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Scheme / SchemeCatalog / database path resolution.

Ported from perl5/MLST/Scheme.pm:1-120, perl5/MLST/Downloader.pm (layout only)
and bin/mlst:96-130 (catalogue listing, --list/--longlist/--info).

Python port of ``MLST::Scheme`` and ``MLST::PubMLST`` plus the database
location ladder. Implements docs/ARCHITECTURE.md sections 4.5 and 6.

Nothing in here does biology: it reads the profile tables and hands the
engine gene lists, signatures and sequence types. Every value it returns is
a ``str`` -- allele numbers and STs are never coerced to ``int`` (section 3,
typing rules; ``"2.002"`` and ``"007"`` are real shipped values).
"""

from __future__ import annotations

import bisect
import json
import os
import re
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .version import BUNDLED_DB_VERSION

__all__ = [
    "DB_VERSION_FALLBACK",
    "HEADER_EXCLUDE_RE",
    "UNKNOWN_DATE",
    "Scheme",
    "SchemeCatalog",
    "SchemeInfo",
    "perl_split",
    "perl_truthy",
    "resolve_blastdb",
    "resolve_datadir",
    "resolve_dbdir",
]

def _errors():
    """The exception classes, imported lazily.

    Section 3.8 puts the hierarchy in ``engine.py`` and section 2.1 makes
    ``engine`` import ``schemes`` -- so a module-level import here would be a
    cycle. Deferring it to the raise sites keeps the module graph acyclic
    without duplicating the classes.
    """
    from . import engine
    return engine


def _db_missing(message: str, user_message: str):
    return _errors().DatabaseMissingError(message, user_message=user_message)


#: Returned by :attr:`Scheme.last_updated` when no usable date is on disk.
UNKNOWN_DATE = "Unknown"

#: Used when ``db/VERSION.txt`` is absent or unreadable.
DB_VERSION_FALLBACK = BUNDLED_DB_VERSION

#: Scheme.pm:46. Six names, anchored, case-SENSITIVE. Do NOT add MLST_cluster
#: -- `aphagocytophilum` really does get 8 "genes" out of this (section 5.19a).
HEADER_EXCLUDE_RE = re.compile(
    r"^(ST|mlst_clade|clonal_complex|species|CC|Lineage)$", re.ASCII)

#: ``>gapA_417`` / ``>penA-2.002``: the allele number at the end of a ``.tfa``
#: header. Same separator class and same decimal allowance as the sseqid regex
#: in ``engine.HIT_RE``; greedy ``.*`` so ``glmU_1_12`` yields ``12``.
_TFA_HEADER_RE = re.compile(r"^>.*[_-](\d+(?:\.\d+)?)$", re.ASCII)

#: Allowed characters in a `database_version.txt` date (Scheme.pm:30).
_DATE_BAD_RE = re.compile(r"[^\d-]", re.ASCII)

#: Profile rows whose first column starts with "ST" are header rows
#: (Scheme.pm:54 -- an unanchored-to-line-end ``m/^ST/``).
_PROFILE_HEADER_RE = re.compile(r"^ST", re.ASCII)


# ---------------------------------------------------------------------------
# Perl emulation helpers
# ---------------------------------------------------------------------------
def perl_truthy(value) -> bool:
    """True exactly when Perl would consider ``value`` true.

    ``undef``, ``""`` and the string ``"0"`` are all false (section 0.1 item 5).
    Eleven shipped schemes carry literal ``0`` allele values, so a plain
    ``if value is not None`` gets this wrong (section 5.19c).
    """
    if value is None:
        return False
    if value == "" or value == "0":
        return False
    return True


def perl_split(sep: str, text: str) -> List[str]:
    """``split /sep/, $text`` with Perl's default LIMIT: trailing empty fields
    are discarded (section 6; ``sep`` is treated as a literal string).
    """
    parts = text.split(sep)
    while parts and parts[-1] == "":
        parts.pop()
    return parts


def _first_line(path: str) -> Optional[str]:
    """First line of ``path``, chomped of ``\\r`` and ``\\n`` (divergence D5).

    Returns ``None`` when the path is not a regular file or cannot be read.
    """
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8", errors="replace",
                  newline=None) as fh:
            line = fh.readline()
    except OSError:
        return None
    return line.rstrip("\r\n")


# ---------------------------------------------------------------------------
# 4.5  Database path resolution
# ---------------------------------------------------------------------------
def _meipass_db() -> Optional[str]:
    base = getattr(sys, "_MEIPASS", None)
    if not base:
        return None
    return os.path.join(base, "wmlst_db")


def _wheel_db() -> Optional[str]:
    try:
        from importlib import resources
    except Exception:  # pragma: no cover - importlib always present on 3.9+
        return None
    try:
        files = getattr(resources, "files", None)
        if files is None:  # pragma: no cover - 3.9 without the new API
            return None
        return str(files("wmlst_db"))
    except Exception:
        return None


def _checkout_db() -> str:
    return os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "db")


def resolve_dbdir(explicit: Optional[str] = None) -> str:
    """Locate the database root. Section 4.5.

    First hit wins:

    1. ``explicit``
    2. ``%WMLST_DBDIR%``
    3. ``%MLST_DBDIR%`` -- drop-in compatibility, divergence D15
    4. frozen: ``sys._MEIPASS/wmlst_db``
    5. installed wheel: ``importlib.resources.files('wmlst_db')``
    6. source checkout: ``<repo>/db``

    Raises :class:`~wmlst.engine.DatabaseMissingError` when the winning
    candidate is not a directory, with upstream's wording (bin/mlst:78).
    """
    candidate = None
    for value in (explicit,
                  os.environ.get("WMLST_DBDIR"),
                  os.environ.get("MLST_DBDIR")):
        if value:
            candidate = value
            break
    if candidate is None:
        for getter in (_meipass_db, _wheel_db):
            value = getter()
            if value and os.path.isdir(value):
                candidate = value
                break
    if candidate is None:
        candidate = _checkout_db()

    candidate = os.path.abspath(os.path.expanduser(candidate))
    if not os.path.isdir(candidate):
        raise _db_missing(
            "Database directory does not exist: %s" % candidate,
            "The MLST allele database was not found at %s." % candidate)
    return candidate


def resolve_datadir(dbdir: str) -> str:
    """``<dbdir>/pubmlst``. Section 4.5."""
    return os.path.join(dbdir, "pubmlst")


def _local_appdata_blastdir() -> str:
    base = os.environ.get("LOCALAPPDATA")
    if not base:
        base = os.path.join(os.path.expanduser("~"), ".local", "share")
    return os.path.join(base, "IOWA-Tech", "WMLST", "blast")


#: The 12 BLAST v5 index extensions makeblastdb -hash_index produces.
_INDEX_EXT = (".ndb", ".nhd", ".nhi", ".nhr", ".nin", ".njs",
              ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto")


def _index_complete(stem: str) -> bool:
    return all(os.path.isfile(stem + ext) for ext in _INDEX_EXT)


def resolve_blastdb(dbdir: str, log=None) -> str:
    """Locate (or choose a home for) the BLAST index stem. Section 4.5.

    Normally ``<dbdir>/blast/mlst.fa``. When that directory holds no complete
    index AND cannot be written to, fall back to the per-user location under
    ``%LOCALAPPDATA%``. A complete index in the writable fallback is preferred
    over an incomplete one in the shipped tree.

    ``log`` is an optional one-argument callable used for the
    ``BLAST database directory is read-only; using <path>`` notice.
    """
    primary_dir = os.path.join(dbdir, "blast")
    primary = os.path.join(primary_dir, "mlst.fa")
    fallback = os.path.join(_local_appdata_blastdir(), "mlst.fa")

    if _index_complete(primary):
        return primary
    if _index_complete(fallback):
        return fallback

    writable = os.path.isdir(primary_dir) and os.access(primary_dir, os.W_OK)
    if not writable and not os.path.isdir(primary_dir):
        writable = os.path.isdir(dbdir) and os.access(dbdir, os.W_OK)
    if writable:
        return primary
    if log is not None:
        log("BLAST database directory is read-only; using %s" % fallback)
    return fallback


# ---------------------------------------------------------------------------
# 4.5  SchemeInfo
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SchemeInfo:
    """The cheap, ``*_info.json``-derived view of a scheme. Section 4.5.

    Used by the GUI catalogue and the HTML report. NEVER consulted by the
    typing path -- no number in here can reach a call.
    """

    name: str
    description: str = ""
    locus: int = 0
    download_date: str = ""
    last_updated: str = ""
    source: str = ""
    api: str = ""
    authenticated: bool = False
    genus: Optional[str] = None
    species: Optional[str] = None
    species_rows: Tuple[Tuple[str, str], ...] = ()
    num_genes: int = 0
    genes: Tuple[str, ...] = ()
    extra: Tuple[Tuple[str, str], ...] = ()


# ---------------------------------------------------------------------------
# 4.5  Scheme
# ---------------------------------------------------------------------------
class Scheme:
    """One ``<datadir>/<name>/`` directory. Section 4.5.

    ``<name>.txt`` is parsed lazily on first property access and then cached:
    eagerly loading all 162 profile tables costs several seconds and ~90 MB
    for no benefit, since a run touches at most a handful.
    """

    __slots__ = (
        "_allele_index",
        "_genes",
        "_genotypes",
        "_last_updated",
        "_num_alleles",
        "dir",
        "name",
    )

    def __init__(self, dir: str, name: str):
        self.dir = dir
        self.name = name
        self._genes = None          # type: Optional[Tuple[str, ...]]
        self._genotypes = None      # type: Optional[Dict[str, str]]
        self._num_alleles = None    # type: Optional[int]
        self._last_updated = None   # type: Optional[str]
        self._allele_index = None   # type: Optional[Dict[str, Tuple[float, ...]]]

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return "Scheme(%r)" % (self.name,)

    # -- files --------------------------------------------------------------
    @property
    def path(self) -> str:
        """The scheme's directory."""
        return os.path.join(self.dir, self.name)

    @property
    def tab_file(self) -> str:
        """``<dir>/<name>/<name>.txt`` (Scheme.pm:17)."""
        return os.path.join(self.dir, self.name, self.name + ".txt")

    def locus_file(self, gene: str) -> str:
        """``<dir>/<name>/<gene>.tfa`` -- one locus' allele FASTA."""
        return os.path.join(self.dir, self.name, gene + ".tfa")

    @property
    def info_file(self) -> str:
        """``<dir>/<name>/<name>_info.json``. Section 6.2."""
        return os.path.join(self.dir, self.name, self.name + "_info.json")

    def _open_tab(self):
        try:
            return open(self.tab_file, encoding="utf-8",
                        errors="replace", newline=None)
        except OSError as exc:
            raise _db_missing(
                "Could not open scheme file: %s" % self.tab_file,
                "The profile table for scheme '%s' is missing."
                % self.name) from exc

    # -- genes --------------------------------------------------------------
    @property
    def genes(self) -> Tuple[str, ...]:
        """Locus names in profile-header order (Scheme.pm:38-47).

        The first line of ``<name>.txt``, chomped, split on a single tab,
        with exactly the six :data:`HEADER_EXCLUDE_RE` names removed. Empty
        fields from consecutive tabs are kept as loci named ``''``.
        """
        if self._genes is None:
            with self._open_tab() as fh:
                header = fh.readline()
            header = header.rstrip("\r\n")
            row = perl_split("\t", header)
            self._genes = tuple(
                c for c in row if not HEADER_EXCLUDE_RE.match(c))
        return self._genes

    @property
    def num_genes(self) -> int:
        """``len(self.genes)``. This is the ``n`` of the scoring formula."""
        return len(self.genes)

    # -- genotypes ----------------------------------------------------------
    @property
    def genotypes(self) -> Dict[str, str]:
        """``signature -> ST`` (Scheme.pm:49-65).

        Columns ``1..num_genes`` are taken **by index**, never by name
        (section 5.19b/B12), so `aphagocytophilum`'s 8th "gene" reads the
        `clonal_complex` column. Perl falsiness turns ``''`` and ``'0'`` into
        ``'-'``. The LAST row wins on a duplicate signature.
        """
        if self._genotypes is None:
            out = {}  # type: Dict[str, str]
            n = self.num_genes
            with self._open_tab() as fh:
                for line in fh:
                    if _PROFILE_HEADER_RE.match(line):
                        continue
                    col = perl_split("\t", line.rstrip("\r\n"))
                    sig = "/".join(
                        col[i] if (i < len(col) and perl_truthy(col[i])) else "-"
                        for i in range(1, n + 1))
                    out[sig] = col[0] if col else ""
            self._genotypes = out
        return self._genotypes

    @property
    def num_genotypes(self) -> int:
        """Number of DISTINCT signatures, not of rows (Scheme.pm:67-71)."""
        return len(self.genotypes)

    @property
    def num_alleles(self) -> int:
        """Distinct ``<gene>_<value>`` strings derived FROM THE SIGNATURES
        (Scheme.pm:73-86), not from the ``.tfa`` files.

        It therefore counts the ``-`` pseudo-allele as a real one. Replicated
        literally because ``--info`` is a byte-identity surface (section 12.5).
        """
        if self._num_alleles is None:
            genes = self.genes
            seen = set()
            for sig in self.genotypes:
                parts = perl_split("/", sig)
                for i, value in enumerate(parts):
                    gene = genes[i] if i < len(genes) else ""
                    seen.add(gene + "_" + value)
            self._num_alleles = len(seen)
        return self._num_alleles

    @property
    def last_updated(self) -> str:
        """First line of ``database_version.txt``, else ``'Unknown'``.

        ``'Unknown'`` also when the date holds any character outside
        ``[0-9-]`` (Scheme.pm:23-32). No shipped scheme has the file, so all
        162 return ``'Unknown'`` -- which is what the ``--info`` golden
        expects (C12/D12).
        """
        if self._last_updated is None:
            date = _first_line(
                os.path.join(self.dir, self.name, "database_version.txt"))
            if date is None or _DATE_BAD_RE.search(date):
                self._last_updated = UNKNOWN_DATE
            else:
                self._last_updated = date
        return self._last_updated

    # -- the two hot functions ---------------------------------------------
    def signature_of(self, calls: Mapping[str, str]) -> str:
        """``'/'.join(calls[g] or '-' for g in self.genes)`` (Scheme.pm:88-95).

        Perl falsiness: a missing key, ``''`` and ``'0'`` all render ``'-'``.
        Gene order is profile-header order. Section 5.9 step 39.
        """
        genes = self.genes
        out = []
        for gene in genes:
            value = calls.get(gene)
            out.append(value if perl_truthy(value) else "-")
        return "/".join(out)

    def sequence_type(self, signature: str) -> str:
        """``genotypes[sig] || '-'`` (Scheme.pm:84-87). Section 5.9 step 40.

        Perl falsiness again: a stored ST of ``''`` or ``'0'`` reads as
        ``'-'``.
        """
        st = self.genotypes.get(signature)
        return st if perl_truthy(st) else "-"

    # -- metadata -----------------------------------------------------------
    def read_info(self) -> Dict[str, object]:
        """The raw ``<name>_info.json`` mapping, or ``{}``. Section 6.2."""
        try:
            with open(self.info_file, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            return {}
        return data if isinstance(data, dict) else {}

    # -- allele registry depth (section 5.14a, tie-break only) --------------
    @property
    def allele_index(self) -> Dict[str, Tuple[float, ...]]:
        """``gene -> sorted tuple of the allele numbers in <gene>.tfa``.

        Read from the ``.tfa`` files, NOT from the profile table: the profile
        only lists alleles that appear in some ST, while the ``.tfa`` is the
        registry as PubMLST issued it. Numbering is NOT dense -- withdrawn
        alleles leave holes (``ecoli_achtman_4.adk`` ships 1914 alleles whose
        highest number is 2060) -- so the empirical rank over these numbers is
        used rather than ``number / count``.

        Loci with no ``.tfa`` are simply absent from the mapping; the only one
        in the shipped DB is ``aphagocytophilum.MLST_cluster``, the pseudo-gene
        section 5.19a deliberately keeps in the header. Unreadable files are
        absent too -- this feeds a tie-break, never a call, so it degrades to
        "no opinion" rather than raising.

        Lazy and cached. Nothing reads it unless two schemes tie on score.
        """
        if self._allele_index is None:
            index = {}  # type: Dict[str, Tuple[float, ...]]
            for gene in self.genes:
                nums = _read_allele_numbers(self.locus_file(gene))
                if nums:
                    index[gene] = nums
            self._allele_index = index
        return self._allele_index

    def allele_percentile(self, gene: str, number: str) -> Optional[float]:
        """Where ``number`` sits in ``gene``'s allele registry, in ``(0, 1]``.

        ``bisect_right(ids, n) / len(ids)`` -- the fraction of the locus'
        registered alleles numbered at or below this one. Allele numbers are
        issued in order of first observation, so a LOW percentile means "an
        allele this scheme's own species has been reporting since the scheme
        opened" and a HIGH one means "a rare, recently registered variant".

        ``None`` when the locus has no registry, or the code is not a number.
        """
        ids = self.allele_index.get(gene)
        if not ids:
            return None
        try:
            value = float(number)
        except (TypeError, ValueError):
            return None
        return bisect.bisect_right(ids, value) / float(len(ids))


def _read_allele_numbers(path: str) -> Tuple[float, ...]:
    """Every allele number in one ``.tfa``, sorted ascending.

    The header regex mirrors :data:`~wmlst.engine.HIT_RE`'s allele group
    (``[_-]`` separator, one optional decimal point, ``re.ASCII``) so the
    ``.tfa`` and the BLAST index agree on what an allele number is. Headers
    that do not match are skipped, exactly as their BLAST hits would be
    (section 5.6, the three unparseable loci).
    """
    out = []  # type: List[float]
    try:
        with open(path, encoding="utf-8", errors="replace", newline=None) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                match = _TFA_HEADER_RE.match(line.rstrip("\r\n"))
                if match is None:
                    continue
                try:
                    out.append(float(match.group(1)))
                except ValueError:      # pragma: no cover - regex forbids it
                    continue
    except OSError:
        return ()
    out.sort()
    return tuple(out)


# ---------------------------------------------------------------------------
# 4.5  SchemeCatalog
# ---------------------------------------------------------------------------
_CANONICAL_INFO_KEYS = ("name", "description", "locus", "download_date",
                        "last_updated", "source", "API", "authenticated")


class SchemeCatalog:
    """Every scheme directory under ``datadir``. Section 4.5.

    ``Scheme`` objects are created once and cached, but each one stays lazy,
    so building the catalogue is a single ``listdir``.
    """

    def __init__(self, datadir: str):
        self.datadir = datadir
        self._names = None          # type: Optional[Tuple[str, ...]]
        self._schemes = {}          # type: Dict[str, Scheme]
        self._species = None        # type: Optional[Dict[str, Tuple]]
        self._db_version = None     # type: Optional[str]

    # -- discovery ----------------------------------------------------------
    def names(self) -> Tuple[str, ...]:
        """Sorted scheme names (PubMLST.pm:20-31, divergence D1).

        Entries of ``datadir`` that are directories and whose name does not
        start with ``'.'``. Upstream uses raw ``readdir`` order; WMLST sorts.
        """
        if self._names is None:
            try:
                entries = os.listdir(self.datadir)
            except OSError as exc:
                raise _db_missing(
                    "Database directory does not exist: %s" % self.datadir,
                    "The MLST allele database was not found at %s."
                    % self.datadir) from exc
            self._names = tuple(sorted(
                e for e in entries
                if not e.startswith(".") and not e.startswith("_")
                and os.path.isdir(os.path.join(self.datadir, e))))
        return self._names

    def __contains__(self, name) -> bool:
        """Set membership, case-SENSITIVE.

        Never a path-existence test: NTFS is case-insensitive, so
        ``Path('db/pubmlst/SAUREUS').exists()`` is True on Windows and would
        let ``--scheme SAUREUS`` through (section 5.2 step 12).
        """
        return name in set(self.names())

    def __getitem__(self, name: str) -> Scheme:
        scheme = self._schemes.get(name)
        if scheme is None:
            if name not in self:
                raise KeyError(name)
            scheme = Scheme(self.datadir, name)
            self._schemes[name] = scheme
        return scheme

    def get(self, name: str) -> Optional[Scheme]:
        """``self[name]`` or ``None``."""
        try:
            return self[name]
        except KeyError:
            return None

    def __iter__(self):
        return iter(self.names())

    def __len__(self) -> int:
        return len(self.names())

    # -- versions -----------------------------------------------------------
    @property
    def dbdir(self) -> str:
        """The parent of ``datadir``."""
        return os.path.dirname(os.path.abspath(self.datadir))

    @property
    def db_version(self) -> str:
        """First line of ``<dbdir>/VERSION.txt``, else the bundled fallback."""
        if self._db_version is None:
            line = _first_line(os.path.join(self.dbdir, "VERSION.txt"))
            self._db_version = line.strip() if line else DB_VERSION_FALLBACK
        return self._db_version

    # -- species map --------------------------------------------------------
    def _species_map(self) -> Dict[str, Tuple[Tuple[str, str], ...]]:
        """``db/scheme_species_map.tab`` as scheme -> all matching rows.

        Many-to-one aware: three schemes have more than one row (C13).
        """
        if self._species is None:
            table = {}  # type: Dict[str, List[Tuple[str, str]]]
            path = os.path.join(self.dbdir, "scheme_species_map.tab")
            try:
                with open(path, encoding="utf-8", errors="replace",
                          newline=None) as fh:
                    for line in fh:
                        line = line.rstrip("\r\n")
                        if not line or line.startswith("#"):
                            continue
                        cols = line.split("\t")
                        name = cols[0].strip()
                        genus = cols[1].strip() if len(cols) > 1 else ""
                        species = cols[2].strip() if len(cols) > 2 else ""
                        if not name:
                            continue
                        table.setdefault(name, []).append((genus, species))
            except OSError:
                table = {}
            self._species = {k: tuple(v) for k, v in table.items()}
        return self._species

    def species_rows(self, name: str) -> Tuple[Tuple[str, str], ...]:
        """All ``(genus, species)`` rows for ``name``; may be empty. Sec 4.5.

        80 of the 162 shipped directories are unmapped. Never fabricate a row.
        """
        return self._species_map().get(name, ())

    def species_label(self, name: str) -> Optional[str]:
        """Human species line, or ``None`` when the scheme is unmapped.

        0 rows -> None; 1 row -> ``'Genus species'`` or ``'Genus sp.'``;
        n rows same genus -> ``'Genus a/b'``; n rows different genera ->
        ``'Genus a / Other b'``.
        """
        rows = self.species_rows(name)
        if not rows:
            return None
        if len(rows) == 1:
            genus, species = rows[0]
            return ("%s %s" % (genus, species)) if species else ("%s sp." % genus)
        genera = {g for g, _ in rows}
        if len(genera) == 1:
            genus = rows[0][0]
            parts = [s for _, s in rows if s]
            if not parts:
                return "%s sp." % genus
            return "%s %s" % (genus, "/".join(parts))
        return " / ".join(
            ("%s %s" % (g, s)) if s else ("%s sp." % g) for g, s in rows)

    # -- info ---------------------------------------------------------------
    def info(self, name: str, *, deep: bool = False) -> SchemeInfo:
        """Build a :class:`SchemeInfo` for ``name``. Section 4.5/6.2.

        Cheap by default: only ``<name>_info.json`` and the species map are
        read. ``deep=True`` additionally loads the profile header so
        ``genes``/``num_genes`` are populated.
        """
        scheme = self[name]
        raw = scheme.read_info()
        rows = self.species_rows(name)
        genus = rows[0][0] if len(rows) == 1 else None
        species = (rows[0][1] or None) if len(rows) == 1 else None
        extra = tuple(
            (k, str(v)) for k, v in raw.items()
            if k not in _CANONICAL_INFO_KEYS)
        try:
            locus = int(raw.get("locus") or 0)
        except (TypeError, ValueError):
            locus = 0
        genes = scheme.genes if deep else ()
        return SchemeInfo(
            name=str(raw.get("name") or name),
            description=str(raw.get("description") or ""),
            locus=locus,
            download_date=str(raw.get("download_date") or ""),
            last_updated=str(raw.get("last_updated") or ""),
            source=str(raw.get("source") or ""),
            api=str(raw.get("API") or ""),
            authenticated=bool(raw.get("authenticated")),
            genus=genus,
            species=species,
            species_rows=rows,
            num_genes=len(genes),
            genes=tuple(genes),
            extra=extra,
        )

    def info_all(self, *, deep: bool = False, progress=None,
                 cancel=None) -> Tuple[SchemeInfo, ...]:
        """:meth:`info` for every scheme, in sorted name order. Section 4.5.

        ``progress(done, total, name)`` is called after each scheme and
        ``cancel`` is any object with ``is_set()``; a set event raises
        :class:`~wmlst.engine.Cancelled`.
        """
        Cancelled = _errors().Cancelled

        names = self.names()
        total = len(names)
        out = []
        for i, name in enumerate(names, 1):
            if cancel is not None and cancel.is_set():
                raise Cancelled("Catalogue scan cancelled")
            out.append(self.info(name, deep=deep))
            if progress is not None:
                progress(i, total, name)
        return tuple(out)
