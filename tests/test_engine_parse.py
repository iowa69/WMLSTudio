# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""BLAST parsing, call building and the winner choice.
docs/ARCHITECTURE.md 5.6, 5.7, 5.8, 5.9, 5.14, 5.15, 5.17, 5.19.

The heavyweight case parses ``tests/golden/example.blast.tsv`` -- 1,423 real
``blastn`` lines captured from the shipped database -- and must derive
sepidermidis ST 184, score 100, PERFECT, with no BLAST+ installed.

Runnable as ``python3 -m pytest tests/test_engine_parse.py`` or ``python3
tests/test_engine_parse.py``.
"""

from __future__ import annotations

import ast
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst import engine as engine_mod
from wmlst.engine import (
    HIT_RE,
    HIT_RE_REPAIR,
    SEEN,
    Engine,
    RunConfig,
    build_signature,
    collect_calls,
    parse_blast,
    parse_blast_full,
    revcom,
    sort_duplicate_codes,
    status_column,
)
from wmlst.schemes import SchemeCatalog, resolve_datadir, resolve_dbdir

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "golden")
PKG = os.path.join(os.path.dirname(HERE), "wmlst")

DATADIR = resolve_datadir(resolve_dbdir())
CATALOG = SchemeCatalog(DATADIR)

BASE_CFG = RunConfig(minid=95.0, mincov=50.0, minscore=50.0)


def _line(sseqid, slen, length, nident, qseqid="c1", qstart=1, qend=10,
          qseq="ACGT", sstrand="plus"):
    return "\t".join(str(x) for x in (sseqid, slen, length, nident, qseqid,
                                      qstart, qend, qseq, sstrand))


# ---------------------------------------------------------------------------
# 5.6 -- the sseqid regex
# ---------------------------------------------------------------------------
def test_regex_splits_on_the_first_dot():
    m = HIT_RE.match(_line("saureus.arcC_1", 456, 456, 456))
    assert m.group(1) == "saureus"
    assert m.group(2) == "arcC"
    assert m.group(3) == "1"


def test_greedy_backtracking_is_load_bearing():
    """`leptospira.glmU_1_12` must give gene 'glmU_1', num '12' (5.6)."""
    m = HIT_RE.match(_line("leptospira.glmU_1_12", 1, 1, 1))
    assert (m.group(2), m.group(3)) == ("glmU_1", "12")


def test_allele_may_be_decimal_and_stays_a_string():
    """B6 -- ngstar.penA_2.002 (test.sh:120)."""
    m = HIT_RE.match(_line("ngstar.penA_2.002", 1, 1, 1))
    assert (m.group(1), m.group(2), m.group(3)) == ("ngstar", "penA", "2.002")


def test_locus_separator_may_be_a_hyphen():
    m = HIT_RE.match(_line("abc.def-7", 1, 1, 1))
    assert (m.group(2), m.group(3)) == ("def", "7")


def test_two_dots_never_parse():
    assert HIT_RE.match(_line("a.b.c_1", 1, 1, 1)) is None


def test_there_is_no_dollar_anchor():
    """Trailing extra columns are tolerated (5.6)."""
    assert HIT_RE.match(_line("s.g_1", 1, 1, 1) + "\textra\tcolumns")


def test_the_three_shipped_unparseable_loci_still_do_not_parse():
    """B8 / 5.19b -- 108 alleles are lost here; reproduce, do not fix."""
    for sseqid in ("halobacteria.EF-2_1", "halobacteria.rpoB'_1",
                   "mhominis_3.p120'_1"):
        assert HIT_RE.match(_line(sseqid, 1, 1, 1)) is None, sseqid


def test_repair_locus_ids_recovers_them():
    """C19 / D11 -- opt-in, and results then do NOT match tseemann/mlst."""
    for sseqid, gene in (("halobacteria.rpoB'_1", "rpoB'"),
                         ("mhominis_3.p120'_1", "p120'")):
        m = HIT_RE_REPAIR.match(_line(sseqid, 1, 1, 1))
        assert m is not None and m.group(2) == gene, sseqid


def test_hit_re_is_ascii():
    assert HIT_RE.flags & re.ASCII
    assert HIT_RE_REPAIR.flags & re.ASCII


def test_ascii_regex_changes_behaviour():
    """Section 0.1 item 1 -- Unicode \\w would accept a non-ASCII scheme name."""
    unicode_re = re.compile(HIT_RE.pattern)
    line = _line("s\u00e4ureus.arcC_1", 1, 1, 1)
    assert unicode_re.match(line) is not None
    assert HIT_RE.match(line) is None


#: The modules that transliterate a Perl regex and where Unicode-vs-byte
#: semantics can therefore change a call or a byte-identity field:
#: `engine` (the sseqid regex, status_column, the duplicate sort), `schemes`
#: (the profile-header filter, the database_version.txt date guard),
#: `any2fasta` (the format sniffers) and `cli` (option parsing).
#: `report` and `updatedb` contain only WMLST-original patterns -- URLs and
#: scheme-description keywords -- which section 0.1's "every PORTED regex"
#: rule does not reach.
PARITY_MODULES = ("engine.py", "schemes.py", "any2fasta.py", "cli.py")


def _compiled_regexes(path):
    """Yield ``(lineno, pattern, flag_source)`` for each ``re.compile`` call."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), path)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "compile"
                and isinstance(func.value, ast.Name) and func.value.id == "re"):
            continue
        parts = []
        for arg in node.args[:1]:
            for sub in ast.walk(arg):
                if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                    parts.append(sub.value)
        flags = [ast.dump(a) for a in node.args[1:]]
        flags += [ast.dump(k.value) for k in node.keywords]
        yield node.lineno, "".join(parts), " ".join(flags)


def test_all_ported_regexes_are_ascii():
    """Acceptance-checklist item 5 / section 0.1 item 1.

    Python's ``\\w``, ``\\d``, ``\\s`` and ``\\b`` are Unicode-aware by
    default; Perl gives them byte semantics here. Every ported pattern that
    uses one MUST pass ``re.ASCII``, or a non-ASCII contig id or scheme name
    silently changes what matches.
    """
    sensitive = re.compile(r"\\[wWdDsSbB]")
    offenders = []
    for fname in PARITY_MODULES:
        path = os.path.join(PKG, fname)
        if not os.path.isfile(path):
            continue                      # module not written yet
        for lineno, pattern, flags in _compiled_regexes(path):
            if not pattern or not sensitive.search(pattern):
                continue
            if "ASCII" not in flags:
                offenders.append("%s:%d %r -- add re.ASCII"
                                 % (fname, lineno, pattern[:60]))
    assert not offenders, "ported regexes missing re.ASCII: %s" % offenders


def test_every_regex_in_the_algorithm_modules_is_ascii():
    """Stricter rule for the two modules that decide a call.

    In `engine` and `schemes` there is no such thing as a non-ported regex,
    so the flag is required unconditionally -- no pattern-content heuristic
    to argue with.
    """
    offenders = []
    for fname in ("engine.py", "schemes.py"):
        for lineno, pattern, flags in _compiled_regexes(os.path.join(PKG, fname)):
            if "ASCII" not in flags:
                offenders.append("%s:%d %r" % (fname, lineno, pattern[:60]))
    assert not offenders, offenders


# ---------------------------------------------------------------------------
# 5.6 -- the mincov gate and hit derivation
# ---------------------------------------------------------------------------
def test_mincov_gates_on_nident_over_slen_not_length_over_slen():
    """B1 -- the highest-value single assertion in this file.

    slen 100, a 60 bp alignment with 40 identities: length/slen is 60 % and
    would pass a 50 % gate, but nident/slen is 40 % and must not.
    """
    cfg = RunConfig(mincov=50.0)
    hits, seen, kept = parse_blast(_line("s.g_1", 100, 60, 40), cfg)
    assert seen == 1 and kept == 0 and hits == []
    hits, seen, kept = parse_blast(_line("s.g_1", 100, 60, 55), cfg)
    assert kept == 1 and hits[0].pct_coverage == 55.0


def test_mincov_boundary_is_inclusive():
    cfg = RunConfig(mincov=50.0)
    assert parse_blast(_line("s.g_1", 100, 50, 50), cfg)[2] == 1
    assert parse_blast(_line("s.g_1", 100, 49, 49), cfg)[2] == 0


def test_seen_counts_parsed_lines_before_the_gate():
    """Perl increments $res_count before the mincov test (step 30)."""
    text = "\n".join([_line("s.g_1", 100, 10, 10),
                      "garbage line that never parses",
                      _line("s.g_2", 100, 100, 100)])
    _hits, seen, kept = parse_blast(text, RunConfig(mincov=50.0))
    assert (seen, kept) == (2, 1)


def test_call_kind_derivation():
    """Section 3.1.1."""
    cfg = RunConfig(mincov=0.0)
    exact = parse_blast(_line("s.g_1", 100, 100, 100), cfg)[0][0]
    novel = parse_blast(_line("s.g_1", 100, 100, 97), cfg)[0][0]
    partial = parse_blast(_line("s.g_1", 100, 60, 60), cfg)[0][0]
    assert (exact.call_kind, novel.call_kind, partial.call_kind) == (
        "exact", "novel", "partial")


def test_slen_and_length_are_not_transposed():
    """The single highest-risk transcription error in this port (3.1)."""
    hit = parse_blast(_line("s.g_1", 500, 300, 290), RunConfig(mincov=0.0))[0][0]
    assert hit.slen == 500 and hit.length == 300 and hit.nident == 290


def test_scheme_filter_and_exclude():
    text = "\n".join([_line("aaa.g_1", 10, 10, 10), _line("bbb.g_1", 10, 10, 10)])
    cfg = RunConfig(mincov=0.0, scheme="aaa", minscore=0.0,
                    exclude=frozenset())
    hits, seen, kept = parse_blast(text, cfg)
    assert seen == 2 and kept == 1 and hits[0].scheme == "aaa"

    cfg = RunConfig(mincov=0.0, exclude=frozenset({"bbb"}))
    hits, excluded, seen, kept = parse_blast_full(text, cfg)
    assert kept == 1 and hits[0].scheme == "aaa"
    assert [h.scheme for h in excluded] == ["bbb"]


def test_crlf_never_reaches_sstrand():
    """Section 0.1 item 4 -- a trailing \\r breaks the minus-strand test."""
    hit = parse_blast(_line("s.g_1", 10, 10, 10, sstrand="minus") + "\r",
                      RunConfig(mincov=0.0))[0][0]
    assert hit.sstrand == "minus"


# ---------------------------------------------------------------------------
# 5.7 -- building calls
# ---------------------------------------------------------------------------
def _walk(lines, capture_novel=False):
    hits = parse_blast("\n".join(lines), RunConfig(mincov=0.0))[0]
    return engine_mod._walk_hits(hits, capture_novel=capture_novel)


def test_exact_then_exact_appends_with_a_comma():
    w = _walk([_line("s.g_1", 10, 10, 10), _line("s.g_2", 10, 10, 10)])
    assert w.res["s"]["g"] == "1,2"
    assert any("WARNING: found additional exact allele match s.g-2" == m
               for m in w.warnings)


def test_a_prior_approximate_code_is_overwritten_not_appended():
    """Step 34 -- the '[~?]' guard."""
    w = _walk([_line("s.g_5", 10, 10, 9), _line("s.g_1", 10, 10, 10)])
    assert w.res["s"]["g"] == "1"
    w = _walk([_line("s.g_5", 10, 6, 6), _line("s.g_1", 10, 10, 10)])
    assert w.res["s"]["g"] == "1"


def test_a_comma_list_keeps_appending():
    """'3,3' contains no [~?] so a third exact hit appends again."""
    w = _walk([_line("s.g_3", 10, 10, 10)] * 3)
    assert w.res["s"]["g"] == "3,3,3"


def test_inexact_is_first_wins():
    """Step 35 -- Perl `||=`."""
    w = _walk([_line("s.g_5", 10, 10, 9), _line("s.g_6", 10, 10, 9)])
    assert w.res["s"]["g"] == "~5"


def test_tilde_is_a_prefix_and_question_a_suffix():
    w = _walk([_line("s.a_5", 10, 10, 9), _line("s.b_6", 10, 6, 6)])
    assert w.res["s"]["a"] == "~5"
    assert w.res["s"]["b"] == "6?"


def test_collect_calls_matches_the_walk():
    lines = [_line("s.a_5", 10, 10, 9), _line("s.b_1", 10, 10, 10)]
    hits = parse_blast("\n".join(lines), RunConfig(mincov=0.0))[0]
    assert collect_calls(hits) == {"s": {"a": "~5", "b": "1"}}


def test_best_hit_tracks_the_code_that_won():
    w = _walk([_line("s.g_5", 10, 10, 9), _line("s.g_1", 10, 10, 10),
               _line("s.g_2", 10, 10, 10)])
    assert w.res["s"]["g"] == "1,2"
    assert w.best["s"]["g"].allele == "1"      # the FIRST comma component
    assert [h.allele for h in w.by_gene["s"]["g"]] == ["5", "1", "2"]


# ---------------------------------------------------------------------------
# 5.7 / 5.17 -- novel capture
# ---------------------------------------------------------------------------
def test_novel_capture_only_for_full_length_hits():
    w = _walk([_line("s.a_5", 10, 10, 9, qseq="ACGTACGTAC"),
               _line("s.b_6", 10, 6, 6, qseq="ACGTAC")], capture_novel=True)
    assert w.nov["s"]["a"] == "ACGTACGTAC"
    assert "b" not in w.nov["s"]


def test_novel_capture_reverse_complements_a_minus_hit():
    w = _walk([_line("s.a_5", 4, 4, 3, qseq="ATGC", sstrand="minus")],
              capture_novel=True)
    assert w.nov["s"]["a"] == "GCAT"


def test_an_exact_hit_writes_the_seen_sentinel_unconditionally():
    """Step 34 -- SEEN overwrites an already-captured sequence."""
    w = _walk([_line("s.a_5", 4, 4, 3, qseq="ATGC"),
               _line("s.a_1", 4, 4, 4, qseq="ATGC")], capture_novel=True)
    assert w.nov["s"]["a"] == SEEN


def test_seen_is_truthy_so_a_later_inexact_hit_cannot_overwrite_it():
    w = _walk([_line("s.a_1", 4, 4, 4, qseq="ATGC"),
               _line("s.a_9", 4, 4, 3, qseq="TTTT")], capture_novel=True)
    assert w.nov["s"]["a"] == SEEN


def test_revcom_does_not_complement_iupac_codes():
    """B7 -- 'fixing' this changes every novel-allele MD5."""
    assert revcom("ATGC") == "GCAT"
    assert revcom("atgc") == "gcat"
    assert revcom("ATGCN") == "NGCAT"
    assert revcom("ATRYGC") == "GCYRAT"


# ---------------------------------------------------------------------------
# 5.9 -- signature, ST, rewrite
# ---------------------------------------------------------------------------
def test_the_null_rewrite_is_slash_anchored():
    """B2 -- a null FIRST locus keeps its '-' and its full penalty."""
    scheme = CATALOG["sepidermidis"]
    pre, st, post = build_signature(scheme, {"arcC": "16", "aroE": "1",
                                             "gtr": "2", "mutS": "1",
                                             "pyrR": "2", "tpiA": "1",
                                             "yqiL": "1"})
    assert (pre, st, post) == ("16/1/2/1/2/1/1", "184", "16/1/2/1/2/1/1")

    # Synthetic: the rewrite only fires when an ST was found.
    assert "-/-".replace("/-", "/0") == "-/0"


def test_rewrite_requires_both_conditions():
    class _Fake:
        genes = ("a", "b", "c")

        def signature_of(self, calls):
            return "/".join(calls.get(g) or "-" for g in self.genes)

        def sequence_type(self, sig):
            return "7" if sig == "-/1/-" else "-"

    fake = _Fake()
    pre, st, post = build_signature(fake, {"b": "1"})
    assert (pre, st, post) == ("-/1/-", "7", "-/1/0")
    pre, st, post = build_signature(fake, {"b": "2"})
    assert (pre, st, post) == ("-/2/-", "-", "-/2/-")


def test_a_real_null_allele_scheme_produces_a_zero_code():
    """Section 5.19c -- eleven shipped schemes carry literal 0 values.

    `efaecium` ST 1421 has a null in column 6, so the rewrite fires and the
    row types PERFECT while rendering `gdh(0)`.
    """
    scheme = CATALOG["efaecium"]
    null_sig = next(s for s in scheme.genotypes
                    if "-" in s and not s.startswith("-"))
    calls = dict(zip(scheme.genes, null_sig.split("/")))
    pre, st, post = build_signature(scheme, calls)
    assert st != "-"
    assert "-" in pre.split("/")
    assert "0" in post.split("/")
    assert status_column(st, 100, post.split("/")) == "PERFECT"


def test_a_leading_null_is_never_rewritten_on_real_data():
    """B2 -- sepidermidis ST 781 has its null in column 1, so it stays '-'.

    Both of sepidermidis's null-bearing profiles are of this shape, which is
    exactly why the rewrite must be slash-anchored and not a bare '-' -> '0'.
    """
    scheme = CATALOG["sepidermidis"]
    leading = next(s for s in scheme.genotypes if s.startswith("-"))
    calls = dict(zip(scheme.genes, leading.split("/")))
    _pre, st, post = build_signature(scheme, calls)
    assert st == scheme.genotypes[leading]
    assert post.startswith("-/")
    assert "0" not in post.split("/")


# ---------------------------------------------------------------------------
# 5.6-5.15 -- the real 1,423-line BLAST capture
# ---------------------------------------------------------------------------
def _engine(**cfg_kw):
    kw = {"minid": 95.0, "mincov": 50.0, "minscore": 50.0, "datadir": DATADIR}
    kw.update(cfg_kw)
    return Engine(RunConfig(**kw), catalog=CATALOG)


def test_golden_blast_capture_yields_sepidermidis_184():
    """The contract case: parse real blastn output, derive the call."""
    with open(os.path.join(GOLDEN, "example.blast.tsv"),
              encoding="utf-8", newline=None) as fh:
        text = fh.read()
    res = _engine().analyse_blast_text("example.fna", text)
    assert res.scheme == "sepidermidis"
    assert res.st == "184"
    assert res.score == 100
    assert res.status == "PERFECT"
    assert res.signature == "16/1/2/1/2/1/1"
    assert [(a.locus, a.code) for a in res.alleles] == [
        ("arcC", "16"), ("aroE", "1"), ("gtr", "2"), ("mutS", "1"),
        ("pyrR", "2"), ("tpiA", "1"), ("yqiL", "1")]
    assert all(a.symbol == "exact" for a in res.alleles)
    assert res.hits_seen == 1423
    assert res.hits_kept == 640
    assert not res.failed


def test_golden_capture_full_row_renders_as_the_golden_tsv():
    with open(os.path.join(GOLDEN, "example.blast.tsv"),
              encoding="utf-8", newline=None) as fh:
        text = fh.read()
    res = _engine(nopath=True).analyse_blast_text("example.fna", text)
    row = "\t".join([res.label, res.scheme, res.st, res.status,
                     str(int(res.score)),
                     ";".join("%s(%s)" % (a.locus, a.code) for a in res.alleles)])
    with open(os.path.join(GOLDEN, "full.out"), encoding="utf-8",
              newline=None) as fh:
        expected = fh.read().splitlines()[1]
    assert row == expected


def test_the_winner_has_a_best_hit_for_every_locus():
    with open(os.path.join(GOLDEN, "example.blast.tsv"),
              encoding="utf-8", newline=None) as fh:
        text = fh.read()
    res = _engine().analyse_blast_text("example.fna", text)
    for call in res.alleles:
        assert call.best is not None, call.locus
        assert call.best.scheme == "sepidermidis"
        assert call.best.allele == call.code
        assert call.best in call.hits


def test_candidates_are_sorted_and_index_zero_is_the_winner():
    with open(os.path.join(GOLDEN, "example.blast.tsv"),
              encoding="utf-8", newline=None) as fh:
        text = fh.read()
    res = _engine().analyse_blast_text("example.fna", text)
    kept = [c for c in res.candidates if not c.below_minscore and not c.excluded]
    assert kept[0].scheme == res.scheme
    assert [c.score for c in kept] == sorted((c.score for c in kept),
                                             reverse=True)
    assert any(c.scheme == "-" for c in kept), "the sentinel must be present"


# ---------------------------------------------------------------------------
# 5.8 / 5.14 -- sentinel and tie-breaking
# ---------------------------------------------------------------------------
def test_no_hits_gives_the_sentinel_and_none():
    res = _engine().analyse_blast_text("none.fa", "")
    assert (res.scheme, res.st, res.score, res.status) == ("-", "-", 0, "NONE")
    assert res.alleles == ()
    assert res.signature == "-/-/-/-/-/-/-"


def test_forced_scheme_widens_the_sentinel_signature():
    """Step 18 -- allele_cnt becomes the forced scheme's locus count."""
    res = _engine(scheme="aphagocytophilum", minscore=0.0,
                  exclude=frozenset()).analyse_blast_text("x.fa", "")
    assert res.scheme == "aphagocytophilum"
    assert res.signature == "/".join(["-"] * 8)
    assert res.status == "NONE"
    assert [a.code for a in res.alleles] == ["-"] * 8


def test_tie_is_broken_by_the_lexicographically_smallest_scheme_name():
    """C8 / D1 -- deterministic where upstream flips a coin.

    Two schemes, identical perfect calls, identical scores: the smaller name
    wins and BOTH tie warnings are emitted.
    """
    class _Fake:
        def __init__(self, name):
            self.name = name
            self.genes = ("a",)
            self.num_genes = 1

        def signature_of(self, calls):
            return calls.get("a") or "-"

        def sequence_type(self, sig):
            return "-"

    fake_catalog = {"zzz": _Fake("zzz"), "aaa": _Fake("aaa")}
    res = {"zzz": {"a": "1"}, "aaa": {"a": "1"}}
    cfg = RunConfig(minscore=50.0)
    kept, dropped, excl = engine_mod._build_candidates(
        fake_catalog, res, {}, cfg, 7)
    kept, _all = engine_mod._order_candidates(kept, dropped, excl)
    assert kept[0].scheme == "aaa"
    assert [c.scheme for c in kept] == ["aaa", "zzz", "-"]


def test_the_sentinel_wins_a_zero_zero_tie():
    """C8 -- the sentinel is inserted first and is never minscore-filtered."""
    res = _engine(minscore=0.0).analyse_blast_text("x.fa", "")
    assert res.scheme == "-"


def test_below_minscore_rows_are_kept_as_evidence_but_never_win():
    """A dropped row can out-score the sentinel; it must not displace it."""
    with open(os.path.join(GOLDEN, "example.blast.tsv"),
              encoding="utf-8", newline=None) as fh:
        text = fh.read()
    res = _engine().analyse_blast_text("example.fna", text)
    dropped = [c for c in res.candidates if c.below_minscore]
    assert dropped, "example.fna produces sub-threshold candidates"
    assert max(c.score for c in dropped) > 0
    assert res.candidates[0].scheme == "sepidermidis"


# ---------------------------------------------------------------------------
# 5.15 -- post-processing
# ---------------------------------------------------------------------------
def test_duplicate_codes_are_numerically_resorted_in_the_result():
    """Step 51 flows into AlleleCall.code and SampleResult.signature."""
    assert sort_duplicate_codes("10,2") == "2,10"
    assert status_column("-", 90, ["2,10"]) == "MIXED"


def test_failed_rows_are_flagged_and_carry_no_call():
    res = _engine().analyse_file(os.path.join(HERE, "no-such-file.fa"))
    assert res.failed is True
    assert res.scheme == "-" and res.st == "-" and res.score == 0
    assert "Unable to read from" in res.error_text


def test_label_falls_back_through_perl_falsiness():
    """Step 21 -- --label '' and --label 0 are ignored."""
    eng = _engine(nopath=True)
    assert eng.label_for("/tmp/dir/x.fa") == "x.fa"
    assert eng.label_for("/tmp/dir/x.fa", label="") == "x.fa"
    assert eng.label_for("/tmp/dir/x.fa", label="0") == "x.fa"
    assert eng.label_for("/tmp/dir/x.fa", label="S1") == "S1"
    assert _engine().label_for("/tmp/dir/x.fa") == "/tmp/dir/x.fa"
    assert _engine(nopath=True).label_for(r"C:\data\x.fa") == "x.fa"


def test_debug_lines_match_the_upstream_format():
    """Step 31."""
    lines = []
    parse_blast(_line("s.g_1", 10, 8, 8, qseqid="ctg1", qstart=5, qend=12,
                      qseq="ACGTACGT", sstrand="minus"),
                RunConfig(mincov=0.0), lines.append)
    assert lines[0] == ("[1] ctg1:5-12(minus) | s g 1 | id=8/8 | cov=8/10 "
                        "| seq=ACGTACGT")


# ---------------------------------------------------------------------------
# End to end, through the real BLAST+ and the real converter
# ---------------------------------------------------------------------------
def _blast_available():
    import shutil
    if shutil.which("blastn") is None:
        return False
    try:
        from wmlst import any2fasta, blastbin  # noqa: F401
    except ImportError:
        return False
    return os.path.isfile(os.path.join(
        os.path.dirname(DATADIR), "blast", "mlst.fa.nin"))


def test_end_to_end_example_fna():
    """The whole of steps 20-52 against the shipped index (section 16 item 9).

    Skipped, not failed, when BLAST+ or the index is absent -- this file's
    other 40-odd tests already cover the algorithm without them.
    """
    if not _blast_available():
        print("SKIP test_end_to_end_example_fna: no blastn / no built index")
        return
    from wmlst.schemes import resolve_blastdb
    from wmlst.schemes import resolve_dbdir as _rdb
    dbdir = _rdb()
    cfg = RunConfig(minid=95.0, mincov=50.0, minscore=50.0, nopath=True,
                    dbdir=dbdir, datadir=DATADIR,
                    blastdb=resolve_blastdb(dbdir), threads=1)
    res = Engine(cfg, catalog=CATALOG).analyse_file(
        os.path.join(HERE, "data", "example.fna"))
    assert not res.failed, res.error_text
    assert (res.label, res.scheme, res.st, res.score, res.status) == (
        "example.fna", "sepidermidis", "184", 100, "PERFECT")
    assert ";".join("%s(%s)" % (a.locus, a.code) for a in res.alleles) == (
        "arcC(16);aroE(1);gtr(2);mutS(1);pyrR(2);tpiA(1);yqiL(1)")
    assert res.n_contigs == 42
    assert res.hits_seen == 1423 and res.hits_kept == 640


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
