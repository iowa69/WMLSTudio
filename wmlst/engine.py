# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""WMLST engine — BLAST hits to allele calls to scheme scores to a SampleResult.

This module is the only place biology happens. Everything downstream
(`report.py`, `cli.py`, `gui.py`) consumes the frozen dataclasses defined below
and MUST NOT re-derive any number from them.

Ported from bin/mlst:300-412 (``find_mlst``: hit parsing, allele calling,
signature scoring, the tie warning) and bin/mlst:160-205 (``main``'s per-file
loop and ``status_column``), with ``perl5/MLST/Scheme.pm`` reached through
``wmlst.schemes``.

See docs/ARCHITECTURE.md sections 3 and 5 for the normative contract.

>>> ------------------------------------------------------------------------
>>> SECTION A: SHARED TYPE CONTRACT  --  FROZEN.
>>> Every other module imports these. Changing a field name or type here is a
>>> cross-module break: update docs/ARCHITECTURE.md section 3 in the same change
>>> or do not make the change at all.
>>> ------------------------------------------------------------------------
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Optional, Tuple

__all__ = [
    "DEFAULT_EXCLUDE",
    "SEEN",
    "SEP",
    "AlleleCall",
    "BlastFailedError",
    "BlastNotFoundError",
    "BootstrapError",
    "Cancelled",
    "DatabaseCorruptError",
    "DatabaseMissingError",
    "EmptyInputError",
    "Hit",
    "NovelAllele",
    "OutOfDiskError",
    "RunConfig",
    "RunMeta",
    "RunResult",
    "SampleResult",
    "SchemeScore",
    "UnsupportedFormatError",
    "UpdateError",
    "WmlstError",
    "allele_count",
    "classify",
    "out_sep",
]

# Python 3.10+ gets slots for free; 3.9 is the floor and must still import.
if sys.version_info >= (3, 10):
    def _frozen(cls):
        return dataclass(frozen=True, slots=True)(cls)
else:  # pragma: no cover - exercised only on 3.9
    def _frozen(cls):
        return dataclass(frozen=True)(cls)


#: Signature field separator used by upstream mlst.
SEP = "/"

#: Sentinel stored in the novel-allele map once an exact hit is seen for a locus.
SEEN = "N"

#: Upstream's default --exclude list (bin/mlst:26).
DEFAULT_EXCLUDE = frozenset(
    {"ecoli", "abaumannii", "vcholerae_2", "senterica_achtman_2"}
)


# ---------------------------------------------------------------------------
# 3.1  Hit
# ---------------------------------------------------------------------------
@_frozen
class Hit:
    """One line of ``blastn -outfmt '6 sseqid slen length nident qseqid qstart
    qend qseq sstrand'`` that matched the sseqid regex AND passed the mincov gate.

    Field order mirrors the BLAST column order EXACTLY. The single highest-risk
    transcription error in this port is swapping ``slen`` (column 2) with
    ``length`` (column 3): ``slen`` is the ALLELE length, ``length`` is the
    ALIGNMENT length.
    """

    sseqid: str          # col 1, verbatim, e.g. "saureus.arcC_1"
    scheme: str          # regex group 1 - everything before the FIRST '.'
    locus: str           # regex group 2 - may contain '_' (e.g. "glmU_1")
    allele: str          # regex group 3 - STRING; may be "12" or "2.002". NEVER int().
    slen: int            # col 2 - Perl $hlen: length of the reference allele
    length: int          # col 3 - Perl $alen: length of the alignment
    nident: int          # col 4
    qseqid: str          # col 5 - contig id, UNTRUSTED user data
    qstart: int          # col 6 - 1-based inclusive, as BLAST reports
    qend: int            # col 7 - 1-based inclusive
    qseq: str            # col 8 - aligned query sequence, UPPERCASE from BLAST
    sstrand: str         # col 9 - exactly "plus" or "minus"
    pct_identity: float  # 100.0 * nident / length
    pct_coverage: float  # 100.0 * nident / slen  <- what --mincov is tested against
    call_kind: str       # "exact" | "novel" | "partial"
    index: int           # 1-based order of appearance in mlst.bls


# ---------------------------------------------------------------------------
# 3.2  AlleleCall
# ---------------------------------------------------------------------------
@_frozen
class AlleleCall:
    """One locus of the reported scheme, in PROFILE-HEADER gene order."""

    locus: str
    code: str            # "16" | "~16" | "16?" | "-" | "0" | "1,2" | "2.002"
    symbol: str          # exact|novel|partial|missing|null|multiple
    best: Optional[Hit]  # the hit that SET the code; None for missing/null
    hits: Tuple[Hit, ...] = ()   # every surviving hit for this locus, by Hit.index


def classify(code: str) -> str:
    """Derive :attr:`AlleleCall.symbol` from ``code`` alone, in this exact order.

    ``"0"`` is tested before the numeric case: a ``0`` code is the ``/-`` -> ``/0``
    rewrite, not an allele named zero.
    """
    if code == "-":
        return "missing"
    if code == "0":
        return "null"
    if "," in code:
        return "multiple"
    if code.startswith("~"):
        return "novel"
    if code.endswith("?"):
        return "partial"
    return "exact"


# ---------------------------------------------------------------------------
# 3.3  SchemeScore
# ---------------------------------------------------------------------------
@_frozen
class SchemeScore:
    """One row of the candidate table, including the sentinel."""

    scheme: str
    st: str
    signature: str       # "/"-joined codes AFTER the /- -> /0 rewrite
    score: int           # final integer score, post +10-if-ST
    num_loci: int        # 0 for the sentinel
    n_exact: int = 0
    n_novel: int = 0     # count of "~" in the PRE-rewrite signature
    n_partial: int = 0   # count of "?"
    n_missing: int = 0   # count of "-" in the POST-rewrite signature
    n_null: int = 0      # count of "0" produced by the rewrite
    excluded: bool = False
    below_minscore: bool = False


# ---------------------------------------------------------------------------
# 3.4  NovelAllele
# ---------------------------------------------------------------------------
@_frozen
class NovelAllele:
    scheme: str
    locus: str
    md5: str             # 32 lowercase hex, computed AFTER revcom
    seq: str             # captured qseq, reverse-complemented iff sstrand == "minus"
    source_label: str    # the sample's LABEL (not its path)
    fasta_id: str        # f"{scheme}.{locus}-{md5}"
    nearest_allele: str
    length_bp: int


# ---------------------------------------------------------------------------
# 3.5  SampleResult
# ---------------------------------------------------------------------------
@_frozen
class SampleResult:
    """Everything about one input file."""

    path: str            # the argv string, verbatim, UNTRUSTED. JSON "filename".
    label: str           # the FILE column and the JSON "id". UNTRUSTED.
    scheme: str          # winning scheme name, or "-"
    st: str              # sequence type STRING, or "-"
    signature: str
    score: int
    status: str          # PERFECT|NONE|NOVEL|MIXED|MISSING|BAD|OK
    alleles: Tuple[AlleleCall, ...] = ()
    candidates: Tuple[SchemeScore, ...] = ()
    novel: Tuple[NovelAllele, ...] = ()
    warnings: Tuple[str, ...] = ()
    n_contigs: int = 0
    total_bp: int = 0
    hits_seen: int = 0
    hits_kept: int = 0
    elapsed_s: float = 0.0
    failed: bool = False
    error: Optional[WmlstError] = None
    error_text: str = ""
    #: Loci of the winning scheme that had an EXACT hit, i.e. upstream's
    #: ``$nov{$LABEL}{$sch}{$gene} = $SEEN`` sentinel (bin/mlst:352). Only
    #: populated when ``--novel`` is in force; it exists so the run-level
    #: merge can reproduce Perl's per-LABEL blocking across input files that
    #: share a label. Never reported.
    novel_seen: Tuple[str, ...] = ()


# ---------------------------------------------------------------------------
# 3.6  RunConfig
# ---------------------------------------------------------------------------
@_frozen
class RunConfig:
    """Every knob, resolved. Built once by cli.py or gui.py; the engine never
    reads argv, os.environ or a settings file."""

    files: Tuple[str, ...] = ()
    minid: float = 95.0
    mincov: float = 50.0
    minscore: float = 50.0
    scheme: Optional[str] = None
    exclude: frozenset = DEFAULT_EXCLUDE
    full: bool = False
    legacy: bool = False
    csv: bool = False
    label: Optional[str] = None
    nopath: bool = False
    outfile: Optional[str] = None
    json_path: Optional[str] = None
    novel_path: Optional[str] = None
    html_path: Optional[str] = None
    evidence_tsv_path: Optional[str] = None
    threads: int = 1
    jobs: int = 1
    quiet: bool = False
    debug: bool = False
    blast_timeout_s: float = 900.0
    dbdir: str = ""
    datadir: str = ""
    blastdb: str = ""
    blastn: Optional[str] = None
    repair_locus_ids: bool = False
    html_evidence: str = "best"
    html_max_novel_bp: int = 100_000


def out_sep(cfg: RunConfig) -> str:
    """Field separator for the tabular outputs."""
    return "," if cfg.csv else "\t"


def allele_count(cfg: RunConfig, catalog=None) -> int:
    """Width of the sentinel signature: 7, or the forced scheme's locus count."""
    if cfg.scheme and catalog is not None:
        return catalog[cfg.scheme].num_genes
    return 7


# ---------------------------------------------------------------------------
# 3.7  RunMeta / RunResult
# ---------------------------------------------------------------------------
@_frozen
class RunMeta:
    wmlst_version: str
    mlst_compat: str
    db_version: str
    db_scheme_count: int
    blast_version: str
    blast_path: str
    dbdir: str
    datadir: str
    blastdb: str
    python_version: str
    platform: str
    hostname: str
    started_utc: str
    started_local: str
    duration_s: float
    argv: Tuple[str, ...]
    config: RunConfig


@_frozen
class RunResult:
    meta: RunMeta
    samples: Tuple[SampleResult, ...] = ()
    novel: Tuple[NovelAllele, ...] = ()


# ---------------------------------------------------------------------------
# 3.8  Exception hierarchy
# ---------------------------------------------------------------------------
class WmlstError(Exception):
    """Base. ``.user_message`` is a one-line, novice-readable sentence."""

    user_message: str = "Something went wrong."

    def __init__(self, message: str = "", user_message: str = ""):
        super().__init__(message or user_message or self.user_message)
        if user_message:
            self.user_message = user_message
        elif message:
            self.user_message = message


class Cancelled(WmlstError):
    user_message = "The run was cancelled."


class EmptyInputError(WmlstError):
    user_message = "That file contains no sequence data."


class UnsupportedFormatError(WmlstError):
    user_message = "That file is not a FASTA, GenBank, EMBL or GFF file."


class BlastNotFoundError(WmlstError):
    user_message = "The BLAST+ search engine could not be found."


class BlastFailedError(WmlstError):
    user_message = "The BLAST+ search engine failed to run."

    def __init__(self, message: str = "", returncode: int = 0,
                 stderr_tail: str = "", user_message: str = ""):
        super().__init__(message, user_message)
        self.returncode = returncode
        self.stderr_tail = stderr_tail


class DatabaseMissingError(WmlstError):
    user_message = "The MLST allele database is missing."


class DatabaseCorruptError(WmlstError):
    user_message = "The MLST allele database looks damaged."


class OutOfDiskError(WmlstError):
    user_message = "There is not enough free disk space."


class BootstrapError(WmlstError):
    user_message = "Setting up the BLAST+ search engine failed."


class UpdateError(WmlstError):
    user_message = "Updating the MLST database failed."


# >>> ------------------------------------------------------------------------
# >>> END OF FROZEN TYPE CONTRACT. The algorithm implementation follows below.
# >>> ------------------------------------------------------------------------


# ===========================================================================
# SECTION B: THE ALGORITHM  --  docs/ARCHITECTURE.md section 5, steps 20-52.
# ===========================================================================
import hashlib
import os
import platform as _platform
import re
import socket
import threading
import time
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from typing import Callable, Dict, List

from .schemes import (
    Scheme,
    SchemeCatalog,
    perl_truthy,
    resolve_blastdb,
    resolve_datadir,
    resolve_dbdir,
)
from .version import UPSTREAM_MLST_VERSION, __version__

__all__ += [
    "HIT_RE",
    "HIT_RE_REPAIR",
    "PHASE_BLAST",
    "PHASE_CONVERT",
    "PHASE_PARSE",
    "PHASE_SCORE",
    "PHASE_WEIGHT",
    "Engine",
    "ProgressEvent",
    "ProgressFn",
    "WarnFn",
    "analyse_file",
    "analyse_files",
    "build_signature",
    "collect_calls",
    "parse_blast",
    "parse_blast_full",
    "revcom",
    "score_signature",
    "sort_duplicate_codes",
    "status_column",
]


# ---------------------------------------------------------------------------
# 4.6 / 5.16  Progress
# ---------------------------------------------------------------------------
PHASE_CONVERT = "convert"
PHASE_BLAST = "blast"
PHASE_PARSE = "parse"
PHASE_SCORE = "score"

#: Share of a file's 0..100 bar each phase owns (section 5.16).
PHASE_WEIGHT = {"convert": 10.0, "blast": 70.0, "parse": 12.0, "score": 8.0}


@_frozen
class ProgressEvent:
    """One progress tick for one input file. Section 4.6 / 5.16.

    ``text`` is novice-facing wording owned by the engine (C17); the GUI and
    CLI render it verbatim and never invent their own.
    """

    phase: str
    percent: float
    text: str
    path: str


ProgressFn = Callable[[ProgressEvent], None]
WarnFn = Callable[[str], None]


# ---------------------------------------------------------------------------
# 4.6 / 5.6  The sseqid regex
# ---------------------------------------------------------------------------
#: bin/mlst:320-325. ``re.ASCII`` IS MANDATORY (section 0.1 item 1) and there
#: is deliberately no ``$`` anchor. Greedy backtracking in group 2 is
#: load-bearing: ``leptospira.glmU_1_12`` must yield gene ``glmU_1``, num ``12``.
HIT_RE = re.compile(
    r'^(\w+)\.(\w+)[_-](\d+(?:\.\d+)?)'
    r'\t(\d+)\t(\d+)\t(\d+)'
    r'\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)',
    re.ASCII)

#: ``--repair-locus-ids`` (C19/D11). Relaxes ONLY the locus group so the 108
#: alleles behind ``halobacteria.rpoB'_N`` and friends parse. Default OFF:
#: results with this on do NOT match tseemann/mlst.
HIT_RE_REPAIR = re.compile(
    r'^([^\t.]+)\.([^\t]+?)[_-](\d+(?:\.\d+)?)'
    r'\t(\d+)\t(\d+)\t(\d+)'
    r'\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)',
    re.ASCII)

#: bin/mlst:344 -- a prior code containing either of these is NOT comma-appended.
_TILDE_Q_RE = re.compile(r"[~?]", re.ASCII)

#: bin/mlst:279 / 5.15 step 52.
_ALL_DIGITS_RE = re.compile(r"\d+", re.ASCII)
_NON_DIGIT_RE = re.compile(r"\D", re.ASCII)

#: bin/mlst:178 / 5.15 step 51.
_NUMERIC_DUP_RE = re.compile(r"[\d.,]+", re.ASCII)

#: Perl tr/ATGCatgc/TACGtacg/. IUPAC codes are reversed but NOT complemented
#: (B7); "fixing" this changes every novel-allele MD5.
_REVCOM_TABLE = str.maketrans("ATGCatgc", "TACGtacg")


def revcom(dna: str) -> str:
    """Reverse-complement, Perl-faithfully (bin/mlst:287-292).

    Reverse FIRST, then translate. Only ``ATGCatgc`` are complemented; every
    other character (``N``, ``R``, ``Y``, ...) is carried through unchanged.
    """
    return dna[::-1].translate(_REVCOM_TABLE)


# ---------------------------------------------------------------------------
# 5.6  parse_blast (steps 28-33)
# ---------------------------------------------------------------------------
def parse_blast_full(text: str, cfg: RunConfig, dbg: WarnFn = None):
    """Parse ``mlst.bls``. Steps 28-33. Returns
    ``(hits, excluded_hits, seen, kept)``.

    ``hits`` survived the mincov gate, the ``--scheme`` filter and
    ``--exclude``; ``excluded_hits`` are the ones ``--exclude`` removed, kept
    only so the HTML report can explain *why* a scheme is absent (B10). They
    can never reach a call.

    ``seen`` is Perl's ``$res_count``: lines that matched the regex, counted
    BEFORE the mincov gate. Lines that do not match are silently skipped and
    are not counted (three shipped loci never parse -- section 5.19b).
    """
    pattern = HIT_RE_REPAIR if cfg.repair_locus_ids else HIT_RE
    forced = cfg.scheme
    exclude = cfg.exclude or frozenset()
    mincov = cfg.mincov

    hits = []           # type: List[Hit]
    excluded = []       # type: List[Hit]
    seen = 0
    kept = 0

    for raw in text.split("\n"):
        line = raw.rstrip("\r\n")
        if not line:
            continue
        m = pattern.match(line)
        if not m:
            if dbg is not None and line:
                dbg("Unparseable BLAST sseqid, dropping line: %s"
                    % line.split("\t", 1)[0])
            continue
        (sch, gene, num, hlen_s, alen_s, nident_s,
         qid, qstart_s, qend_s, qseq, sstrand) = m.groups()
        seen += 1
        hlen = int(hlen_s)
        alen = int(alen_s)
        nident = int(nident_s)

        if dbg is not None:
            dbg("[%d] %s:%s-%s(%s) | %s %s %s | id=%d/%d | cov=%d/%d | seq=%s"
                % (seen, qid, qstart_s, qend_s, sstrand, sch, gene, num,
                   nident, alen, alen, hlen, qseq))

        # Step 32. nident/hlen, NOT alen/hlen (B1). Never cross-multiply:
        # both sides are IEEE-754 doubles and boundary cases flip if you do.
        if not (nident / hlen >= mincov / 100):
            continue

        if forced and sch != forced:
            if dbg is not None:
                dbg("Skipping %s.%s.%s allele as user specified --scheme %s"
                    % (sch, gene, num, forced))
            continue

        dropped = sch in exclude
        if dropped and dbg is not None:
            dbg("Excluding %s.%s.%s due to --exclude option" % (sch, gene, num))

        if alen == hlen:
            call_kind = "exact" if nident == hlen else "novel"
        else:
            call_kind = "partial"

        hit = Hit(
            sseqid=line.split("\t", 1)[0],
            scheme=sch, locus=gene, allele=num,
            slen=hlen, length=alen, nident=nident,
            qseqid=qid, qstart=int(qstart_s), qend=int(qend_s),
            qseq=qseq, sstrand=sstrand,
            pct_identity=100.0 * nident / alen,
            pct_coverage=100.0 * nident / hlen,
            call_kind=call_kind,
            index=seen,
        )
        if dropped:
            excluded.append(hit)
            continue
        kept += 1
        hits.append(hit)

    return hits, excluded, seen, kept


def parse_blast(text: str, cfg: RunConfig, dbg: WarnFn = None):
    """``-> (hits, seen, kept)``. Steps 27-33. See :func:`parse_blast_full`."""
    hits, _excluded, seen, kept = parse_blast_full(text, cfg, dbg)
    return hits, seen, kept


# ---------------------------------------------------------------------------
# 5.7  Building calls (steps 34-37)
# ---------------------------------------------------------------------------
class _Walk:
    """Mutable accumulator for one pass over the surviving hits."""

    __slots__ = ("best", "by_gene", "messages", "nov", "res", "warnings")

    def __init__(self):
        self.res = {}       # type: Dict[str, Dict[str, str]]
        self.nov = {}       # type: Dict[str, Dict[str, str]]
        self.best = {}      # type: Dict[str, Dict[str, Hit]]
        self.by_gene = {}   # type: Dict[str, Dict[str, List[Hit]]]
        self.messages = []  # type: List[str]
        self.warnings = []  # type: List[str]


def _walk_hits(hits: Iterable[Hit], capture_novel: bool = False) -> _Walk:
    """The exact/inexact loop of bin/mlst:342-364. Steps 34-37.

    Hits MUST arrive in BLAST output order; every first-wins rule below
    depends on it.
    """
    w = _Walk()
    for hit in hits:
        sch, gene, num = hit.scheme, hit.locus, hit.allele
        res = w.res.setdefault(sch, {})
        nov = w.nov.setdefault(sch, {})
        best = w.best.setdefault(sch, {})
        w.by_gene.setdefault(sch, {}).setdefault(gene, []).append(hit)

        if hit.slen == hit.length and hit.nident == hit.slen:
            # Step 34, exact branch.
            prior = res.get(gene)
            if gene in res and not _TILDE_Q_RE.search(prior):
                line = ("WARNING: found additional exact allele match %s.%s-%s"
                        % (sch, gene, num))
                w.messages.append(line)
                w.warnings.append(line)
                res[gene] = prior + "," + num
                # `best` stays the hit that created the FIRST component.
            else:
                w.messages.append("Found exact allele match %s.%s-%s"
                                  % (sch, gene, num))
                res[gene] = num
                best[gene] = hit
            # Unconditional: an exact match wipes any captured novel sequence.
            nov[gene] = SEEN
        else:
            # Step 35, inexact branch. Perl `||=`, so a stored "0" would be
            # overwritten -- keep the falsiness, do not use `in`.
            label = ("~" + num) if hit.length == hit.slen else (num + "?")
            if not perl_truthy(res.get(gene)):
                res[gene] = label
                best[gene] = hit
            if capture_novel and hit.length == hit.slen:
                seq = revcom(hit.qseq) if hit.sstrand == "minus" else hit.qseq
                if not perl_truthy(nov.get(gene)):
                    nov[gene] = seq
    return w


def collect_calls(hits: Iterable[Hit]) -> Dict[str, Dict[str, str]]:
    """``-> dict[scheme][gene] -> code``. Steps 34-37 (pure; no novel capture)."""
    return _walk_hits(hits, capture_novel=False).res


# ---------------------------------------------------------------------------
# 5.9  Signature, ST, rewrite and score (steps 39-45)
# ---------------------------------------------------------------------------
def build_signature(scheme: Scheme, calls: Mapping[str, str]):
    """``-> (pre_rewrite_signature, st, post_rewrite_signature)``. Steps 39-41.

    The ST lookup uses the PRE-rewrite signature because that is the form the
    genotype table stores. The ``/-`` -> ``/0`` rewrite is slash-anchored, so
    a null FIRST locus keeps its ``-`` and its full ``-1.0`` penalty (B2), and
    it runs BEFORE scoring, so a rewritten null costs nothing (B3).
    """
    pre = scheme.signature_of(calls)
    st = scheme.sequence_type(pre)
    if st != "-" and "/-" in pre:
        return pre, st, pre.replace("/-", "/0")
    return pre, st, pre


def score_signature(signature: str, num_loci: int, st: str) -> int:
    """EXACT port of bin/mlst:382-392. Step 42.

    Plain ``float``, this exact operation order, ``int()`` truncation toward
    zero. ``Decimal``, ``Fraction``, ``round()``, ``math.floor()``, ``//`` and
    any algebraic simplification are forbidden: ``int(s*90/n)`` differs from
    ``int(s*(90/n))`` for 31 input shapes, and exact-rational arithmetic
    differs from IEEE-754 for 79 of the 986 tuples in tests/golden/scores.tsv.
    """
    n = num_loci
    s = float(n)
    s -= 0.3 * signature.count("~")
    s -= 0.5 * signature.count("?")
    s -= 1.0 * signature.count("-")
    score = int(s * 90 / n)
    if st != "-":
        score += 10
    return score


def sort_duplicate_codes(code: str) -> str:
    """bin/mlst:177-184. Step 51.

    Applies ONLY when the field holds a comma AND matches ``^[\\d.,]+$``, so
    ``~5,6`` and ``5?,6`` keep their discovery order (B5).
    """
    if "," not in code:
        return code
    if not _NUMERIC_DUP_RE.fullmatch(code):
        return code
    try:
        return ",".join(sorted(code.split(","), key=float))
    except ValueError:
        return code


def status_column(st: str, score: int, codes: Sequence[str]) -> str:
    """EXACT port of bin/mlst:275-284. Step 52.

    THE ORDER OF THE TESTS IS THE SPECIFICATION. ``PERFECT`` is decided from
    the ST alone, so a scheme can be PERFECT while carrying rewritten ``0``
    codes; ``MISSING`` tests whether a code *contains* ``-``, not equals it.
    """
    if _ALL_DIGITS_RE.fullmatch(st):
        return "PERFECT"
    if score <= 0:
        return "NONE"
    if not any(_NON_DIGIT_RE.search(c) for c in codes):
        return "NOVEL"
    if any("," in c for c in codes):
        return "MIXED"
    if any("-" in c for c in codes):
        return "MISSING"
    if score < 70:
        return "BAD"
    return "OK"


# ---------------------------------------------------------------------------
# 5.8 / 5.9 / 5.14  Candidate table
# ---------------------------------------------------------------------------
def _score_one(scheme: Scheme, name: str, calls: Mapping[str, str],
               minscore: float, excluded: bool, dbg: WarnFn = None) -> SchemeScore:
    """One row of the candidate table. Steps 39-45."""
    pre, st, sig = build_signature(scheme, calls)
    n = scheme.num_genes
    score = score_signature(sig, n, st)
    if dbg is not None:
        dbg("SCORE=%d\t%s\t%s\t%s\t(%d genes)" % (score, name, st, sig, n))
    fields = sig.split("/")
    return SchemeScore(
        scheme=name, st=st, signature=sig, score=score, num_loci=n,
        n_exact=sum(1 for c in fields if classify(c) == "exact"),
        n_novel=pre.count("~"),
        n_partial=sig.count("?"),
        n_missing=sig.count("-"),
        n_null=sum(1 for c in fields if c == "0"),
        excluded=excluded,
        below_minscore=(not excluded) and not (score >= minscore),
    )


def _build_candidates(catalog: SchemeCatalog, res: Mapping[str, Mapping[str, str]],
                      excluded_res: Mapping[str, Mapping[str, str]],
                      cfg: RunConfig, allele_cnt: int, dbg: WarnFn = None,
                      warn: WarnFn = None):
    """``-> (kept, all_rows)``. Steps 38-45.

    ``kept`` is upstream's ``@sig`` after the ``--minscore`` filter, sentinel
    first, then ``sorted(res)`` (D1/C8) -- so a stable sort on ``-score``
    makes the lexicographically smallest scheme name win a tie and lets the
    sentinel win a 0-0 tie.

    Scheme names come from the BLAST index (the sseqid prefix) but the
    :class:`~wmlst.schemes.Scheme` objects come from the datadir, so a
    mismatched ``--blastdb``/``--datadir`` pair -- or a datadir caught
    mid-update -- can name a scheme that has no profile. Such a name is
    dropped with a warning (the same guard the excluded rows already use)
    instead of raising ``KeyError`` out of the whole batch: D10 says a
    per-file problem must not abort the run.
    """
    sentinel = SchemeScore(
        scheme=(cfg.scheme or "-"), st="-",
        signature="/".join(["-"] * allele_cnt),
        score=0, num_loci=0, n_missing=allele_cnt)
    kept = [sentinel]           # type: List[SchemeScore]
    dropped = []                # type: List[SchemeScore]
    for name in sorted(res):
        if name not in catalog:
            line = ("WARNING: BLAST index names scheme '%s' which is not in %s"
                    " - ignoring it" % (name, cfg.datadir or "the datadir"))
            if warn is not None:
                warn(line)
            if dbg is not None:
                dbg(line)
            continue
        row = _score_one(catalog[name], name, res[name], cfg.minscore,
                         False, dbg)
        (dropped if row.below_minscore else kept).append(row)
    excluded_rows = [
        _score_one(catalog[name], name, excluded_res[name], cfg.minscore,
                   True, dbg)
        for name in sorted(excluded_res) if name in catalog
    ]
    return kept, dropped, excluded_rows


def _order_candidates(kept: List[SchemeScore], dropped: List[SchemeScore],
                      excluded_rows: List[SchemeScore]):
    """Step 46 plus the evidence tail.

    Python's sort is stable, so the insertion order set up by
    :func:`_build_candidates` breaks ties. Filtered and excluded rows are
    appended AFTER every kept row so index 0 is always the winner -- a
    below-minscore row can out-score the sentinel and must never displace it.
    """
    kept = sorted(kept, key=lambda c: -c.score)
    dropped = sorted(dropped, key=lambda c: -c.score)
    excluded_rows = sorted(excluded_rows, key=lambda c: -c.score)
    return kept, tuple(kept) + tuple(dropped) + tuple(excluded_rows)


# ---------------------------------------------------------------------------
# Deferred, defensive imports of the owner-B modules
# ---------------------------------------------------------------------------
def _import_blastbin():
    """The BLAST backend (`wmlst.blastbin`, owner B).

    Imported lazily so `engine` stays importable -- and unit-testable -- in a
    tree where `blastbin.py` has not landed yet, and so `report.py` can import
    the dataclasses without dragging BLAST discovery in.
    """
    from . import blastbin
    return blastbin


def _import_any2fasta():
    """The converter (`wmlst.any2fasta`, owner B). Lazy for the same reason."""
    from . import any2fasta
    return any2fasta


def _basename(path: str) -> str:
    """Basename that splits on BOTH separators, whatever the host OS.

    A run on Linux may be handed a Windows path out of a saved GUI session,
    and `--nopath` must still produce the file name (step 21).
    """
    return re.split(r"[\\/]", path)[-1]


def _utc_now():
    return datetime.now(timezone.utc)


# ---------------------------------------------------------------------------
# 4.6  Engine
# ---------------------------------------------------------------------------
class Engine:
    """Runs section 5 over the files in ``cfg``. Section 4.6.

    Stateless apart from a cached :class:`~wmlst.schemes.SchemeCatalog`, so a
    single instance is safe to share across the ``cfg.jobs`` worker threads.

    ``backend`` and ``converter`` exist purely as injection seams for tests
    and for integrations that supply their own BLAST launcher; leave them
    ``None`` and the real :mod:`wmlst.blastbin` / :mod:`wmlst.any2fasta` are
    used. ``engine`` itself never launches a process (section 2.1):
    every child is started inside ``blastbin``.
    """

    def __init__(self, cfg: RunConfig, tools=None, *, catalog=None,
                 backend=None, converter=None, argv: Sequence[str] = ()):
        self.cfg = cfg
        self._tools = tools
        self._catalog = catalog
        self._backend = backend
        self._converter = converter
        self.argv = tuple(argv)
        # Re-entrant: `tools` resolves `backend` while holding it.
        self._lock = threading.RLock()
        self._started = _utc_now()
        self._started_monotonic = time.time()

        # 3.6.1 -- programmer errors, not user errors.
        assert not (cfg.legacy and cfg.scheme is None), \
            "RunConfig.legacy requires a scheme (3.6.1 invariant 1)"
        assert cfg.scheme is None or cfg.minscore == 0.0, \
            "RunConfig.scheme forces minscore=0 (3.6.1 invariant 2)"
        assert cfg.scheme is None or not cfg.exclude, \
            "RunConfig.scheme clears exclude (3.6.1 invariant 2)"
        assert cfg.label is None or len(cfg.files) <= 1, \
            "RunConfig.label forbidden with >1 file (3.6.1 invariant 3)"
        assert cfg.threads >= 1 and cfg.jobs >= 1

    # -- lazily resolved collaborators -------------------------------------
    @property
    def catalog(self) -> SchemeCatalog:
        """The scheme catalogue, built once from ``cfg.datadir``."""
        with self._lock:
            if self._catalog is None:
                datadir = self.cfg.datadir or resolve_datadir(
                    self.cfg.dbdir or resolve_dbdir())
                self._catalog = SchemeCatalog(datadir)
            return self._catalog

    @property
    def backend(self):
        """The BLAST launcher module (:mod:`wmlst.blastbin` by default)."""
        with self._lock:
            if self._backend is None:
                try:
                    self._backend = _import_blastbin()
                except ImportError as exc:
                    raise BlastNotFoundError(
                        "wmlst.blastbin is not available: %s" % exc,
                        user_message=("The BLAST+ search engine could not be "
                                      "started because part of WMLST is "
                                      "missing.")) from exc
            return self._backend

    @property
    def converter(self):
        """The FASTA converter module (:mod:`wmlst.any2fasta` by default)."""
        with self._lock:
            if self._converter is None:
                try:
                    self._converter = _import_any2fasta()
                except ImportError as exc:
                    raise WmlstError(
                        "wmlst.any2fasta is not available: %s" % exc,
                        user_message=("WMLST cannot read sequence files "
                                      "because part of it is missing.")) from exc
            return self._converter

    @property
    def tools(self):
        """The resolved :class:`BlastTools`, discovered on first use."""
        with self._lock:
            if self._tools is None:
                self._tools = self.backend.find_blast(self.cfg.blastn)
            return self._tools

    def close(self) -> None:
        """Release cached collaborators. Section 4.6."""
        with self._lock:
            self._catalog = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    # -- metadata -----------------------------------------------------------
    def _blast_identity(self):
        try:
            tools = self.tools
        except WmlstError:
            return "", ""
        return getattr(tools, "version", ""), getattr(tools, "blastn", "")

    def _make_meta(self, duration_s: float) -> RunMeta:
        cfg = self.cfg
        version, path = self._blast_identity()
        try:
            n_schemes = len(self.catalog.names())
            db_version = self.catalog.db_version
        except WmlstError:
            n_schemes, db_version = 0, ""
        started = self._started
        return RunMeta(
            wmlst_version=__version__,
            mlst_compat=UPSTREAM_MLST_VERSION,
            db_version=db_version,
            db_scheme_count=n_schemes,
            blast_version=version,
            blast_path=path,
            dbdir=cfg.dbdir,
            datadir=cfg.datadir,
            blastdb=cfg.blastdb,
            python_version=_platform.python_version(),
            platform=_platform.platform(),
            hostname=socket.gethostname(),
            started_utc=started.strftime("%Y-%m-%dT%H:%M:%SZ"),
            started_local=started.astimezone().isoformat(timespec="seconds"),
            duration_s=duration_s,
            argv=self.argv,
            config=cfg,
        )

    @property
    def meta_stub(self) -> RunMeta:
        """:class:`RunMeta` with ``duration_s == 0.0``, for pre-run display."""
        return self._make_meta(0.0)

    # -- helpers ------------------------------------------------------------
    @staticmethod
    def _check_cancel(cancel):
        if cancel is not None and cancel.is_set():
            raise Cancelled("The run was cancelled")

    @staticmethod
    def _emit(progress, phase, percent, text, path):
        if progress is not None:
            progress(ProgressEvent(phase=phase, percent=percent, text=text,
                                   path=path))

    def label_for(self, path: str, label=None) -> str:
        """``--label``, else the basename under ``--nopath``, else the path.

        Step 21, with Perl falsiness: ``--label ''`` and ``--label 0`` both
        fall through to the filename.
        """
        if perl_truthy(label):
            return label
        if perl_truthy(self.cfg.label):
            return self.cfg.label
        return _basename(path) if self.cfg.nopath else path

    def _allele_cnt(self) -> int:
        """Step 18: 7, or the forced scheme's locus count."""
        if self.cfg.scheme:
            try:
                return self.catalog[self.cfg.scheme].num_genes
            except (KeyError, WmlstError):
                return 7
        return 7

    def _blastdb(self) -> str:
        if self.cfg.blastdb:
            return self.cfg.blastdb
        dbdir = self.cfg.dbdir or resolve_dbdir()
        return resolve_blastdb(dbdir)

    # -- the per-file pipeline ---------------------------------------------
    def _failed(self, path: str, label: str, exc: WmlstError,
                elapsed: float, allele_cnt: int) -> SampleResult:
        """A ``failed=True`` row (divergence D10).

        It never reaches the compat TSV/CSV/JSON -- upstream produces no row
        for a file that kills it -- but the GUI, the HTML report and the
        evidence TSV need to show what happened.
        """
        return SampleResult(
            path=path, label=label, scheme="-", st="-",
            signature="/".join(["-"] * allele_cnt), score=0, status="NONE",
            failed=True, error=exc, error_text=str(exc), elapsed_s=elapsed)

    def analyse_file(self, path: str, *, label=None, progress: ProgressFn = None,
                     warn: WarnFn = None, msg: WarnFn = None,
                     dbg: WarnFn = None, cancel=None) -> SampleResult:
        """The whole per-file pipeline, steps 20-52. Section 4.6.

        NEVER raises for a per-file problem: an unreadable input, a converter
        failure or a non-zero ``blastn`` returns ``SampleResult(failed=True)``
        so a batch can continue. It DOES raise for run-fatal conditions
        (:class:`BlastNotFoundError`, :class:`DatabaseMissingError`) and for
        :class:`Cancelled`.

        ``warn`` receives Perl ``wrn()`` lines (never suppressed by
        ``--quiet``); ``msg`` receives ``msg()`` lines (suppressed); ``dbg``
        receives ``--debug`` lines. All three are optional.
        """
        cfg = self.cfg
        t0 = time.time()
        lbl = self.label_for(path, label)
        allele_cnt = self._allele_cnt()
        self._check_cancel(cancel)

        if msg is not None and lbl != path:
            msg("Using label '%s' for file %s" % (lbl, path))

        # Step 20: readability BEFORE the directory test, exactly as upstream.
        if not os.access(path, os.R_OK):
            return self._failed(
                path, lbl, WmlstError("Unable to read from '%s'" % path,
                                      user_message="That file could not be read."),
                time.time() - t0, allele_cnt)
        if os.path.isdir(path):
            return self._failed(
                path, lbl,
                WmlstError("'%s' seems to be a directory, not a file" % path,
                           user_message="That is a folder, not a file."),
                time.time() - t0, allele_cnt)

        backend = self.backend
        blastdb = self._blastdb()
        blastdb_dir = os.path.dirname(blastdb)
        db_basename = os.path.basename(blastdb)
        tools = self.tools

        self._emit(progress, PHASE_CONVERT, 0.0,
                   "Reading %s…" % _basename(path), path)

        with backend.job_dir() as jd:
            fna = os.path.join(jd, "mlst.fna")
            bls = os.path.join(jd, "mlst.bls")

            # Step 24: an input path blastn's ANSI file API cannot open is
            # staged into the (guaranteed-safe) job directory first. The
            # DISPLAYED filename always stays the user's original.
            src = path
            if hasattr(backend, "ansi_safe") and backend.ansi_safe(path) is None:
                staged = os.path.join(jd, "input.fa")
                with open(path, "rb") as fin, open(staged, "wb") as fout:
                    while True:
                        chunk = fin.read(1 << 20)
                        if not chunk:
                            break
                        fout.write(chunk)
                src = staged

            # Steps 25-26.
            try:
                n_contigs, total_bp, _fmt = self.converter.convert_to_file(src, fna)
            except WmlstError as exc:
                return self._failed(path, lbl, exc, time.time() - t0, allele_cnt)
            except OSError as exc:
                return self._failed(
                    path, lbl,
                    WmlstError(str(exc), user_message="That file could not be read."),
                    time.time() - t0, allele_cnt)

            self._emit(progress, PHASE_CONVERT, 10.0,
                       "Read {0} contigs ({1:,} bp).".format(n_contigs, total_bp),
                       path)
            self._check_cancel(cancel)

            # Step 23. Upstream is safe only because a failed command kills
            # the process; WMLST continues past a bad file, so a stale .bls
            # from the previous file would be parsed as this file's hits.
            try:
                os.remove(bls)
            except OSError:
                pass

            try:
                n_schemes = len(self.catalog.names())
            except WmlstError:
                n_schemes = 0
            self._emit(progress, PHASE_BLAST, 10.0,
                       "Searching %d schemes…" % n_schemes, path)

            # Step 27.
            query = fna
            out = bls
            if hasattr(backend, "ansi_safe"):
                query = backend.ansi_safe(fna) or fna
                out = backend.ansi_safe(bls) or bls
            argv = backend.blastn_argv(
                tools, query=query, out=out, db_basename=db_basename,
                threads=cfg.threads, minid=cfg.minid)
            try:
                proc = backend.run_tool(
                    argv, cwd=blastdb_dir, env=backend.child_env(blastdb_dir),
                    timeout=cfg.blast_timeout_s, cancel=cancel)
            except Cancelled:
                raise
            except WmlstError as exc:
                return self._failed(path, lbl, exc, time.time() - t0, allele_cnt)

            stderr = (getattr(proc, "stderr", "") or "")
            if isinstance(stderr, bytes):
                stderr = stderr.decode("utf-8", "replace")
            if proc.returncode != 0:
                tail = "\n".join(stderr.strip().splitlines()[-5:])
                return self._failed(
                    path, lbl,
                    BlastFailedError(
                        tail or ("blastn exited with status %s" % proc.returncode),
                        returncode=proc.returncode, stderr_tail=tail,
                        user_message=("BLAST+ could not search this file."
                                      if not tail else tail)),
                    time.time() - t0, allele_cnt)

            self._emit(progress, PHASE_BLAST, 80.0, "Search finished.", path)
            self._check_cancel(cancel)

            # Step 28: universal newlines AND an explicit rstrip. A stray
            # "\r" on the last column silently turns sstrand into "plus\r".
            try:
                with open(bls, encoding="utf-8", errors="replace",
                          newline=None) as fh:
                    text = fh.read()
            except FileNotFoundError:
                text = ""

        self._emit(progress, PHASE_PARSE, 80.0, "Reading matches…", path)
        return self._finish(path, lbl, text, allele_cnt, n_contigs, total_bp,
                            t0, progress, warn, msg, dbg)

    def analyse_blast_text(self, path: str, text: str, *, label=None,
                           n_contigs: int = 0, total_bp: int = 0,
                           progress: ProgressFn = None, warn: WarnFn = None,
                           msg: WarnFn = None, dbg: WarnFn = None) -> SampleResult:
        """Steps 29-52 over already-captured ``mlst.bls`` text. Section 5.6-5.15.

        The half of the pipeline that launches no child process: useful for
        re-deriving a call from stored BLAST output and for testing the
        algorithm without BLAST+ installed.
        """
        return self._finish(path, self.label_for(path, label), text,
                            self._allele_cnt(), n_contigs, total_bp,
                            time.time(), progress, warn, msg, dbg)

    def _finish(self, path, lbl, text, allele_cnt, n_contigs, total_bp, t0,
                progress=None, warn=None, msg=None, dbg=None) -> SampleResult:
        """Steps 29-52: parse, call, score, choose, post-process."""
        cfg = self.cfg
        catalog = self.catalog

        hits, excluded_hits, seen, kept = parse_blast_full(text, cfg, dbg)
        self._emit(progress, PHASE_PARSE, 92.0, "%d alleles matched." % kept, path)
        self._emit(progress, PHASE_SCORE, 92.0, "Choosing the best scheme…",
                   path)

        walk = _walk_hits(hits, capture_novel=bool(cfg.novel_path))
        if msg is not None:
            for line in walk.messages:
                msg(line)
        excluded_walk = _walk_hits(excluded_hits, capture_novel=False)

        kept_rows, dropped_rows, excluded_rows = _build_candidates(
            catalog, walk.res, excluded_walk.res, cfg, allele_cnt, dbg, warn)
        kept_rows, all_rows = _order_candidates(kept_rows, dropped_rows,
                                                excluded_rows)

        # Step 47: every equal-first, via wrn() -- NOT suppressed by --quiet.
        warnings = list(walk.warnings)
        top = kept_rows[0]
        for other in kept_rows[1:]:
            if other.score == top.score:
                line = ("WARNING: %s(%s)==%s(%s) score=%s %s"
                        % (top.scheme, top.st, other.scheme, other.st,
                           str(int(top.score)), path))
                warnings.append(line)
                if warn is not None:
                    warn(line)

        winner = top
        if cfg.scheme:
            assert winner.scheme == cfg.scheme, (
                "BUG: got back %s despite --scheme %s" % (winner.scheme, cfg.scheme))

        # Steps 50-51.
        if winner.scheme == "-":
            codes = []      # deliberately empty even though the signature is -/-/...
            genes = ()
        else:
            codes = [sort_duplicate_codes(c) for c in winner.signature.split(SEP)]
            genes = catalog[winner.scheme].genes

        status = status_column(winner.st, winner.score, codes)

        # Step 36: attach evidence, in profile-header gene order.
        alleles = []
        if winner.scheme != "-":
            best = walk.best.get(winner.scheme, {})
            by_gene = walk.by_gene.get(winner.scheme, {})
            for i, gene in enumerate(genes):
                code = codes[i] if i < len(codes) else "-"
                symbol = classify(code)
                alleles.append(AlleleCall(
                    locus=gene, code=code, symbol=symbol,
                    best=(None if symbol in ("missing", "null")
                          else best.get(gene)),
                    hits=tuple(sorted(by_gene.get(gene, ()),
                                      key=lambda h: h.index))))

        # Steps 49 + 5.17: prune novel alleles to the winning scheme.
        # bin/mlst:233-241 walks `keys %{$nov{$fname}{$sch}}`, i.e. EVERY locus
        # the BLAST index produced a full-length inexact hit for -- not only
        # the loci the profile header lists. A locus present in the index but
        # missing from the header (a custom --blastdb, or an upstream update
        # that ships a .tfa before the profile catches up) still gets its
        # novel sequence written. Header order first, then the extras sorted,
        # so D7's deterministic record order survives.
        novel = []
        novel_seen = ()
        if cfg.novel_path and winner.scheme != "-":
            nov_map = walk.nov.get(winner.scheme, {})
            novel_seen = tuple(sorted(g for g in nov_map if nov_map[g] == SEEN))
            extra = sorted(g for g in nov_map if g not in genes)
            for gene in list(genes) + extra:
                seq = nov_map.get(gene)
                if seq is None or seq == SEEN:
                    continue
                digest = hashlib.md5(seq.encode("ascii", "replace")).hexdigest()
                hit = walk.best.get(winner.scheme, {}).get(gene)
                novel.append(NovelAllele(
                    scheme=winner.scheme, locus=gene, md5=digest, seq=seq,
                    source_label=lbl,
                    fasta_id="%s.%s-%s" % (winner.scheme, gene, digest),
                    nearest_allele=(hit.allele if hit is not None else ""),
                    length_bp=len(seq)))

        signature = (SEP.join(codes) if winner.scheme != "-"
                     else winner.signature)
        self._emit(progress, PHASE_SCORE, 100.0,
                   ("%s ST %s" % (winner.scheme, winner.st))
                   if winner.scheme != "-" else "No scheme matched.", path)

        return SampleResult(
            path=path, label=lbl, scheme=winner.scheme, st=winner.st,
            signature=signature, score=winner.score, status=status,
            alleles=tuple(alleles), candidates=all_rows, novel=tuple(novel),
            novel_seen=novel_seen,
            warnings=tuple(warnings), n_contigs=n_contigs, total_bp=total_bp,
            hits_seen=seen, hits_kept=kept, elapsed_s=time.time() - t0,
            failed=False, error=None, error_text="")

    @staticmethod
    def _merge_novel(results) -> List[NovelAllele]:
        """bin/mlst's per-LABEL %nov, collapsed across the run. See analyse()."""
        live = [(i, r) for i, r in enumerate(results)
                if r is not None and not r.failed]

        # The last file to name a label decides which scheme survives its map.
        final_scheme = {}
        for _i, res in live:
            final_scheme[res.label] = res.scheme

        # Walk backwards: a file contributes only while it -- and every later
        # file sharing its label -- won that surviving scheme.
        contributes = {}
        still_ok = {}
        for i, res in reversed(live):
            ok = (still_ok.get(res.label, True)
                  and res.scheme != "-"
                  and res.scheme == final_scheme[res.label])
            contributes[i] = ok
            still_ok[res.label] = ok

        blocked = set()          # (label, locus) an exact hit closed off
        for i, res in live:
            if contributes[i]:
                for locus in res.novel_seen:
                    blocked.add((res.label, locus))

        taken = set()            # (label, locus) already filled -- `||=`
        seen_ids = set()         # $seen{$id}++
        novel = []
        for i, res in live:
            if not contributes[i]:
                continue
            for allele in res.novel:
                key = (res.label, allele.locus)
                if key in blocked or key in taken:
                    continue
                taken.add(key)
                if allele.fasta_id in seen_ids:
                    continue
                seen_ids.add(allele.fasta_id)
                novel.append(allele)
        return novel

    # -- the run ------------------------------------------------------------
    def analyse(self, *, progress: ProgressFn = None, warn: WarnFn = None,
                msg: WarnFn = None, dbg: WarnFn = None, cancel=None) -> RunResult:
        """Every file in ``cfg.files``, then :class:`RunResult`. Section 4.6.

        Files run on a ``cfg.jobs``-wide thread pool and the results are
        re-sorted into argv order before the tuple is built (section 9);
        ``RunResult.samples`` is NEVER in completion order.
        """
        cfg = self.cfg
        self._started = _utc_now()
        t0 = time.time()
        files = list(cfg.files)
        results = [None] * len(files)   # type: List[Any]

        def one(i_path):
            i, p = i_path
            return i, self.analyse_file(p, progress=progress, warn=warn,
                                        msg=msg, dbg=dbg, cancel=cancel)

        if cfg.jobs > 1 and len(files) > 1:
            with ThreadPoolExecutor(max_workers=cfg.jobs) as pool:
                for i, res in pool.map(one, list(enumerate(files))):
                    results[i] = res
        else:
            for i, p in enumerate(files):
                results[i] = self.analyse_file(p, progress=progress, warn=warn,
                                               msg=msg, dbg=dbg, cancel=cancel)

        # 5.17: collapse the per-sample novel maps the way bin/mlst does.
        #
        # Upstream's %nov is keyed by LABEL, not by file (bin/mlst:85, 362),
        # so two inputs that share a label share ONE map, and:
        #   * bin/mlst:170-172 deletes, after EVERY file, every scheme in that
        #     label's map that is not that file's winning scheme -- so only
        #     the trailing run of same-label files that all won the same
        #     scheme can contribute anything;
        #   * bin/mlst:352 stores the $SEEN sentinel unconditionally, so an
        #     exact hit in ANY contributing file blocks that locus outright;
        #   * bin/mlst:362 is `||=`, so among the contributors the FIRST
        #     sequence for a locus wins;
        #   * bin/mlst:240's `$seen{$id}++` then dedups identical sequences
        #     across labels.
        # Distinct labels (the normal case) are unaffected by all of this.
        # Records stay in (sample argv index, scheme, gene) order -- D7.
        novel = self._merge_novel(results)
        if cfg.novel_path and msg is not None:
            msg("Found %d novel alleles" % len(novel))

        return RunResult(meta=self._make_meta(time.time() - t0),
                         samples=tuple(results), novel=tuple(novel))


# ---------------------------------------------------------------------------
# Module-level entry points
# ---------------------------------------------------------------------------
def analyse_file(path: str, cfg: RunConfig, catalog: SchemeCatalog = None,
                 meta=None, *, label=None, tools=None, backend=None,
                 converter=None, progress: ProgressFn = None,
                 warn: WarnFn = None, msg: WarnFn = None, dbg: WarnFn = None,
                 cancel=None) -> SampleResult:
    """Type one file. Convenience wrapper over :meth:`Engine.analyse_file`.

    ``meta`` is accepted and ignored: a :class:`RunMeta` describes a run, and
    a single-file call has no run to describe. It is in the signature so
    callers that already hold one need no special case.
    """
    engine = Engine(cfg, tools, catalog=catalog, backend=backend,
                    converter=converter)
    return engine.analyse_file(path, label=label, progress=progress, warn=warn,
                               msg=msg, dbg=dbg, cancel=cancel)


def analyse_files(cfg: RunConfig, catalog: SchemeCatalog = None, *, tools=None,
                  backend=None, converter=None, argv: Sequence[str] = (),
                  progress: ProgressFn = None, warn: WarnFn = None,
                  msg: WarnFn = None, dbg: WarnFn = None,
                  cancel=None) -> RunResult:
    """Type every file in ``cfg.files``. Wrapper over :meth:`Engine.analyse`."""
    engine = Engine(cfg, tools, catalog=catalog, backend=backend,
                    converter=converter, argv=argv)
    return engine.analyse(progress=progress, warn=warn, msg=msg, dbg=dbg,
                          cancel=cancel)
