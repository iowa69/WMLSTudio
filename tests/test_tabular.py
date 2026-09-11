# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Byte-identity tests for the compat emitters in ``wmlst.report``.

Covers docs/ARCHITECTURE.md section 12. Fixtures are built BY HAND from the
frozen dataclasses in ``wmlst.engine`` so this suite does not depend on the
engine's algorithm being finished; the expectations are the real-Perl goldens
in ``tests/golden/``.

Runs under ``python3 -m pytest tests/test_tabular.py`` and as a plain script
(``python3 tests/test_tabular.py``).
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wmlst import report
from wmlst.engine import (
    AlleleCall,
    Hit,
    NovelAllele,
    RunConfig,
    RunMeta,
    RunResult,
    SampleResult,
    classify,
)

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = os.path.join(HERE, "golden")


# ---------------------------------------------------------------------------
# fixture builders
# ---------------------------------------------------------------------------

def read_golden(name: str) -> str:
    with open(os.path.join(GOLDEN, name), encoding="utf-8", newline="") as fh:
        return fh.read()


def parse_alleles(text: str):
    """'arcC(16);aroE(1)' -> tuple[AlleleCall, ...] in the given order."""
    calls = []
    if not text:
        return tuple(calls)
    for item in text.split(";"):
        locus, _, rest = item.partition("(")
        code = rest[:-1]
        calls.append(AlleleCall(locus=locus, code=code, symbol=classify(code),
                                best=None, hits=()))
    return tuple(calls)


def mk_sample(path, scheme, st, status, score, alleles="", label=None, failed=False,
              error_text="", candidates=(), novel=(), warnings=()):
    calls = parse_alleles(alleles)
    signature = "/".join(c.code for c in calls)
    return SampleResult(
        path=path, label=label if label is not None else path, scheme=scheme, st=st,
        signature=signature, score=score, status=status, alleles=calls,
        candidates=candidates, novel=novel, warnings=warnings,
        n_contigs=1, total_bp=1000, hits_seen=0, hits_kept=0, elapsed_s=0.5,
        failed=failed, error=None, error_text=error_text,
    )


def mk_meta(cfg: RunConfig) -> RunMeta:
    return RunMeta(
        wmlst_version="1.0.0", mlst_compat="2.35.0", db_version="2025-12-29",
        db_scheme_count=162, blast_version="2.17.0+", blast_path="/usr/bin/blastn",
        dbdir=os.path.join(os.path.dirname(HERE), "db"),
        datadir=os.path.join(os.path.dirname(HERE), "db", "pubmlst"),
        blastdb=os.path.join(os.path.dirname(HERE), "db", "blast", "mlst.fa"),
        python_version="3.12.0", platform="Linux-x86_64", hostname="testhost",
        started_utc="2026-01-01T00:00:00Z", started_local="2026-01-01T01:00:00+01:00",
        duration_s=1.25, argv=("mlst", "--full", "example.fna"), config=cfg,
    )


def mk_result(samples, cfg=None, novel=()):
    cfg = cfg or RunConfig()
    return RunResult(meta=mk_meta(cfg), samples=tuple(samples), novel=tuple(novel))


def tsv(result) -> str:
    buf = io.StringIO()
    report.write_tsv(result, buf)
    return buf.getvalue()


# The known-good calls, transcribed from the golden corpus' fixtures.
EXAMPLE = ("sepidermidis", "184", "PERFECT", 100,
           "arcC(16);aroE(1);gtr(2);mutS(1);pyrR(2);tpiA(1);yqiL(1)")
MESSY = ("hparasuis", "-", "MISSING", 51,
         "atpD(~66);infB(~80);mdh(~83);rpoB(~34);6pgd(114?);g3pd(-);frdB(~117)")
MIXED = ("mgenitalium", "-", "MIXED", 90,
         "MLST_adk(7);MLST_atpA(1,1);MLST_gmk(1);MLST_gyrB(1,1);MLST_pgm(3,3);MLST_ppa(1,1)")
NOVEL = ("mgenitalium", "-", "NOVEL", 90,
         "MLST_adk(7);MLST_atpA(1);MLST_gmk(1);MLST_gyrB(1);MLST_pgm(3);MLST_ppa(1)")
NONE = ("-", "-", "NONE", 0, "")
ISSUE146 = ("salmonella", "64", "PERFECT", 100,
            "aroC(10);dnaN(14);hemD(15);hisD(31);purE(25);sucA(20);thrA(33)")
EQUALITY = ("salmonella", "3529", "PERFECT", 100,
            "aroC(684);dnaN(603);hemD(554);hisD(836);purE(690);sucA(648);thrA(700)")
LEPTO = ("leptospira_2", "-", "OK", 74,
         "adk_2(85);glmU_2(~24);icdA_2(~46);lipL32_2(20);lipL41_2(39);"
         "mreA_2(~56);pntA_2(~26)")
FORCED = ("saureus", "-", "NONE", 0,
          "arcC(-);aroE(-);glpF(-);gmk(-);pta(-);tpi(-);yqiL(-)")


# ---------------------------------------------------------------------------
# 12.1  print_row
# ---------------------------------------------------------------------------

def test_print_row_plain():
    buf = io.StringIO()
    report.print_row(["a", "b", "c"], "\t", buf)
    assert buf.getvalue() == "a\tb\tc\n"


def test_print_row_quotes_on_separator_only():
    # A comma is NOT quoted in TSV mode and a tab is NOT quoted in CSV mode:
    # the trigger is the CURRENT separator (section 12.1).
    buf = io.StringIO()
    report.print_row(["x,y"], "\t", buf)
    assert buf.getvalue() == "x,y\n"
    buf = io.StringIO()
    report.print_row(["x\ty"], ",", buf)
    assert buf.getvalue() == "x\ty\n"
    buf = io.StringIO()
    report.print_row(["x,y"], ",", buf)
    assert buf.getvalue() == '"x,y"\n'


def test_print_row_backslash_doubling_order():
    # quote doubling FIRST, then backslash doubling: a"b\c -> "a""b\\c"
    buf = io.StringIO()
    report.print_row(['a"b\\c'], "\t", buf)
    assert buf.getvalue() == '"a""b\\\\c"\n'


def test_print_row_backslash_alone_is_not_quoted():
    buf = io.StringIO()
    report.print_row(["a\\b"], "\t", buf)
    assert buf.getvalue() == "a\\b\n"


def test_print_row_newline_does_not_trigger_quoting():
    buf = io.StringIO()
    report.print_row(["a\nb"], "\t", buf)
    assert buf.getvalue() == "a\nb\n"


def test_print_row_pinned_csv_case():
    # test.sh:127-130 -- `mlst --csv mixed.fa.zip` line 0 contains ,"MLST_gyrB(1,1)",
    result = mk_result([mk_sample("mixed.fa.zip", *MIXED)], RunConfig(csv=True))
    line = tsv(result)
    assert ',"MLST_gyrB(1,1)",' in line


# ---------------------------------------------------------------------------
# 12.2 / 12.3  the golden rows
# ---------------------------------------------------------------------------

def test_golden_default():
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], RunConfig())
    assert tsv(result) == read_golden("default.out")


def test_golden_csv():
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], RunConfig(csv=True))
    assert tsv(result) == read_golden("csv.out")


def test_golden_full_family():
    cases = {
        "full.out": ("example.fna", EXAMPLE, {}),
        "nopath.out": ("example.fna", EXAMPLE, {"nopath": True}),
        "gz.out": ("example.fna.gz", EXAMPLE, {}),
        "gbk.out": ("example.gbk.gz", EXAMPLE, {}),
        "messy.out": ("messy.fa", MESSY, {}),
        "mixedzip.out": ("mixed.fa.zip", MIXED, {}),
        "novelfa.out": ("novel.fa", NOVEL, {}),
        "nonefa.out": ("none.fa", NONE, {}),
        "issue146.out": ("issue146.fa", ISSUE146, {}),
        "equality.out": ("equality.fa.gz", EQUALITY, {}),
        "novelbz2.out": ("novel.fasta.bz2", LEPTO, {}),
    }
    for golden, (path, spec, extra) in sorted(cases.items()):
        cfg = RunConfig(full=True, **extra)
        result = mk_result([mk_sample(path, *spec)], cfg)
        assert tsv(result) == read_golden(golden), golden


def test_golden_label():
    cfg = RunConfig(full=True, label="MYSAMPLE")
    result = mk_result([mk_sample("example.fna", *EXAMPLE, label="MYSAMPLE")], cfg)
    assert tsv(result) == read_golden("label.out")


def test_golden_scheme_forced():
    # --scheme forces minscore to 0 and clears --exclude (3.6.1).
    cfg = RunConfig(full=True, scheme="saureus", minscore=0.0, exclude=frozenset())
    result = mk_result([mk_sample("example.fna", *FORCED)], cfg)
    assert tsv(result) == read_golden("scheme_forced.out")


def test_golden_legacy():
    cfg = RunConfig(legacy=True, scheme="sepidermidis", minscore=0.0,
                    exclude=frozenset())
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], cfg)
    assert tsv(result) == read_golden("legacy.out")


def test_legacy_header_without_any_sample():
    # No usable sample to take the gene order from: fall back to the profile
    # header on disk so the --legacy header is still correct.
    cfg = RunConfig(legacy=True, scheme="sepidermidis", minscore=0.0,
                    exclude=frozenset(),
                    datadir=os.path.join(os.path.dirname(HERE), "db", "pubmlst"))
    result = mk_result([mk_sample("null.fa", "-", "-", "NONE", 0, failed=True,
                                  error_text="empty")], cfg)
    assert tsv(result) == read_golden("legacy.out").splitlines(True)[0]


def test_failed_rows_are_excluded_from_compat_output():
    # emptyfa.out and nullfa.out are header-only: upstream died before the row.
    cfg = RunConfig(full=True)
    for path, golden in (("empty.fa", "emptyfa.out"), ("null.fa", "nullfa.out")):
        result = mk_result(
            [mk_sample(path, "-", "-", "NONE", 0, failed=True, error_text="empty")], cfg)
        assert tsv(result) == read_golden(golden), golden


def test_golden_multi():
    # upstream aborts on null.fa, so its golden holds only the two rows that
    # were printed before the failure (D10).
    cfg = RunConfig(full=True)
    result = mk_result([
        mk_sample("example.fna", *EXAMPLE),
        mk_sample("novel.fa", *NOVEL),
        mk_sample("null.fa", "-", "-", "NONE", 0, failed=True, error_text="empty"),
    ], cfg)
    assert tsv(result) == read_golden("multi.out")


def test_no_scheme_full_row_has_trailing_separator():
    cfg = RunConfig(full=True)
    result = mk_result([mk_sample("none.fa", *NONE)], cfg)
    assert tsv(result).splitlines(True)[1] == "none.fa\t-\t-\tNONE\t0\t\n"


def test_no_scheme_default_row_has_three_fields():
    result = mk_result([mk_sample("none.fa", *NONE)], RunConfig())
    assert tsv(result) == "none.fa\t-\t-\n"


def test_no_carriage_returns_anywhere():
    cfg = RunConfig(full=True)
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], cfg)
    assert "\r" not in tsv(result)


def test_score_is_written_as_a_bare_integer():
    cfg = RunConfig(full=True)
    result = mk_result([mk_sample("messy.fa", *MESSY)], cfg)
    assert "\t51\t" in tsv(result)
    assert "51.0" not in tsv(result)


# ---------------------------------------------------------------------------
# 12.4  JSON
# ---------------------------------------------------------------------------

def test_json_matches_golden_structurally():
    cfg = RunConfig(json_path="x.json")
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], cfg)
    ours = json.loads(report.json_text(result))
    theirs = json.loads(read_golden("example.json"))
    assert ours == theirs


def test_json_formatting_is_perl_pretty():
    # indent=3, " : " separator, ascii-only -- proven by re-encoding the golden
    # document in ITS OWN key order and getting the golden bytes back.
    raw = read_golden("example.json")
    import collections
    decoded = json.loads(raw, object_pairs_hook=collections.OrderedDict)
    reencoded = json.dumps(decoded, indent=3, separators=(",", " : "),
                           ensure_ascii=True)
    # Perl's JSON::PP appends a newline that Python's json.dumps does not;
    # C5 mandates no trailing newline, so compare against the stripped golden.
    assert reencoded == raw.rstrip("\n")
    assert raw.count("\n") == reencoded.count("\n") + 1


def test_json_has_one_trailing_newline_and_fixed_key_order():
    cfg = RunConfig(json_path="x.json")
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], cfg)
    text = report.json_text(result)
    # Perl's to_json(pretty=>1) ends with a single LF; so do we.
    assert text.endswith("]\n")
    assert not text.endswith("\n\n")
    keys = re.findall(r'^      "(\w+)" :', text, re.M)
    assert keys == ["id", "filename", "scheme", "sequence_type", "alleles"]


def test_json_alleles_null_when_no_scheme():
    result = mk_result([mk_sample("none.fa", *NONE)], RunConfig())
    obj = json.loads(report.json_text(result))
    assert obj[0]["alleles"] is None


def test_json_every_value_is_a_string():
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], RunConfig())
    obj = json.loads(report.json_text(result))[0]
    assert obj["sequence_type"] == "184" and isinstance(obj["sequence_type"], str)
    assert all(isinstance(v, str) for v in obj["alleles"].values())


def test_json_alleles_are_in_gene_order():
    result = mk_result([mk_sample("example.fna", *EXAMPLE)], RunConfig())
    text = report.json_text(result)
    order = re.findall(r'^         "(\w+)" :', text, re.M)
    assert order == ["arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]


def test_json_honours_label_and_path(tmpdir=None):
    cfg = RunConfig(label="MYSAMPLE")
    result = mk_result(
        [mk_sample("tests/data/example.fna", *EXAMPLE, label="MYSAMPLE")], cfg)
    obj = json.loads(report.json_text(result))[0]
    assert obj["id"] == "MYSAMPLE"
    assert obj["filename"] == "tests/data/example.fna"


# ---------------------------------------------------------------------------
# 12.6  novel FASTA
# ---------------------------------------------------------------------------

def novel_from_golden(name: str):
    """Rebuild NovelAllele fixtures from a golden FASTA's records."""
    alleles = []
    text = read_golden(name)
    for block in text.split(">"):
        if not block.strip():
            continue
        header, seq = block.split("\n", 1)
        seq = seq.strip()
        fasta_id, _, label = header.partition(" ")
        scheme, _, rest = fasta_id.partition(".")
        locus, _, md5 = rest.rpartition("-")
        alleles.append(NovelAllele(
            scheme=scheme, locus=locus, md5=md5, seq=seq, source_label=label,
            fasta_id=fasta_id, nearest_allele="1", length_bp=len(seq)))
    return alleles


def test_novel_fasta_is_byte_identical(tmp_path=None):
    import tempfile
    for golden in ("messy_novel.fa", "lepto_novel.fa"):
        alleles = novel_from_golden(golden)
        assert alleles, golden
        result = mk_result([mk_sample("x.fa", *NONE)], RunConfig(novel_path="n.fa"),
                           novel=alleles)
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "novel.fa")
            report.write_novel_fasta(result, out)
            with open(out, encoding="utf-8", newline="") as fh:
                assert fh.read() == read_golden(golden), golden


def test_novel_fasta_ids_are_md5_of_the_sequence():
    for golden in ("messy_novel.fa", "lepto_novel.fa"):
        for allele in novel_from_golden(golden):
            assert hashlib.md5(allele.seq.encode("ascii")).hexdigest() == allele.md5


def test_novel_fasta_dedups_by_id_first_wins():
    import tempfile
    base = novel_from_golden("messy_novel.fa")[0]
    dup = NovelAllele(scheme=base.scheme, locus=base.locus, md5=base.md5,
                      seq=base.seq, source_label="OTHER.fa", fasta_id=base.fasta_id,
                      nearest_allele=base.nearest_allele, length_bp=base.length_bp)
    result = mk_result([mk_sample("x.fa", *NONE)], RunConfig(), novel=[base, dup])
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "n.fa")
        report.write_novel_fasta(result, out)
        text = open(out, encoding="utf-8").read()
    assert text.count(">") == 1
    assert "OTHER.fa" not in text


def test_novel_fasta_empty_file_is_still_created():
    import tempfile
    result = mk_result([mk_sample("x.fa", *NONE)], RunConfig())
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "n.fa")
        report.write_novel_fasta(result, out)
        assert os.path.isfile(out)
        assert open(out, encoding="utf-8").read() == ""


# ---------------------------------------------------------------------------
# 12.5  --list / --longlist / --info
# ---------------------------------------------------------------------------

class FakeScheme:
    def __init__(self, name, genes, genotypes, alleles):
        self.name = name
        self.genes = genes
        self.num_genes = len(genes)
        self.num_genotypes = genotypes
        self.num_alleles = alleles
        self.last_updated = "Unknown"


class FakeCatalog:
    def __init__(self, schemes):
        self._s = {s.name: s for s in schemes}

    def names(self):
        return tuple(sorted(self._s))

    def __getitem__(self, name):
        return self._s[name]


CATALOG = FakeCatalog([
    FakeScheme("saureus", ("arcC", "aroE", "glpF"), 5, 30),
    FakeScheme("abaumannii", ("gltA", "gyrB"), 2, 9),
])


def test_write_list():
    buf = io.StringIO()
    report.write_list(CATALOG, buf)
    assert buf.getvalue() == "abaumannii saureus\n"


def test_write_longlist():
    buf = io.StringIO()
    report.write_longlist(CATALOG, buf, "\t")
    assert buf.getvalue() == "abaumannii\tgltA\tgyrB\nsaureus\tarcC\taroE\tglpF\n"


def test_write_info_has_six_columns_and_keeps_LOCII():
    buf = io.StringIO()
    report.write_info(CATALOG, buf, "\t")
    lines = buf.getvalue().splitlines()
    assert lines[0] == "SCHEME\tLOCII\tTYPES\tALLELES\tDATE\tLOCII_NAMES"
    assert lines[1] == "abaumannii\t2\t2\t9\tUnknown\tgltA gyrB"
    assert lines[2] == "saureus\t3\t5\t30\tUnknown\tarcC aroE glpF"
    assert all(len(line.split("\t")) == 6 for line in lines)


# ---------------------------------------------------------------------------
# 12.7  the WMLST-only evidence dump
# ---------------------------------------------------------------------------

def mk_hit(**kw):
    base = {"sseqid": "sepidermidis.arcC_16", "scheme": "sepidermidis", "locus": "arcC",
                "allele": "16", "slen": 465, "length": 465, "nident": 465, "qseqid": "contig_1",
                "qstart": 100, "qend": 564, "qseq": "ACGT", "sstrand": "plus",
                "pct_identity": 100.0, "pct_coverage": 100.0, "call_kind": "exact", "index": 1}
    base.update(kw)
    return Hit(**base)


def test_evidence_tsv_columns_and_rows():
    import tempfile
    hit = mk_hit()
    call = AlleleCall(locus="arcC", code="16", symbol="exact", best=hit, hits=(hit,))
    sample = SampleResult(path="example.fna", label="example.fna",
                          scheme="sepidermidis", st="184", signature="16",
                          score=100, status="PERFECT", alleles=(call,))
    result = mk_result([sample], RunConfig())
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "e.tsv")
        report.write_evidence_tsv(result, out)
        lines = open(out, encoding="utf-8").read().splitlines()
    assert lines[0].split("\t") == list(report.EVIDENCE_COLUMNS)
    row = lines[1].split("\t")
    assert row[0] == "example.fna" and row[2] == "arcC" and row[5] == "yes"
    assert row[6] == "contig_1" and row[9] == "plus"
    assert row[11] == "465" and row[12] == "465"
    assert row[13] == "100.0" and row[14] == "100.0"


def test_evidence_tsv_includes_failed_samples():
    import tempfile
    result = mk_result([mk_sample("bad.fa", "-", "-", "NONE", 0, failed=True,
                                  error_text="The input appears to be empty")],
                       RunConfig())
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "e.tsv")
        report.write_evidence_tsv(result, out)
        text = open(out, encoding="utf-8").read()
    assert "bad.fa" in text and "failed" in text


def test_inexact_percentages_are_floored_not_rounded():
    # C14: an inexact call must never print 100.0 %.
    assert report._pct_text(99.99, exact=False) == "99.9"
    assert report._pct_text(99.99, exact=True) == "100.0"
    assert report._pct_text(86.666, exact=False) == "86.6"


# ---------------------------------------------------------------------------

def _main():
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("ok   " + name)
            except AssertionError as exc:
                failures += 1
                print("FAIL " + name + ": " + str(exc))
            except Exception as exc:
                failures += 1
                print("ERROR " + name + ": " + repr(exc))
    print("{0} failure(s)".format(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
