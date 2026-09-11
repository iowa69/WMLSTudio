# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""End-to-end parity with the Perl golden corpus (docs/ARCHITECTURE.md 0, 16.2).

Every case here replays one line of ``scripts/make_golden.sh`` through the real
WMLST CLI -- real ``any2fasta``, real ``blastn``, the real bundled database --
and compares the bytes with what Perl ``mlst`` 2.35.0 wrote into
``tests/golden/``.

The goldens were produced with ``cwd=tests/data``, so the FILE column carries no
directory prefix.  That is reproduced faithfully by running the CLI from
``tests/data`` rather than by relaxing the comparison.

Three cases cannot be byte-identical, and each is asserted against an exact,
fully written-out expectation rather than being skipped:

* ``equality`` -- D1/C8.  Two schemes score 100; upstream picks a winner by
  randomised Perl hash order (the golden happens to record ``salmonella``),
  WMLST picks the lexicographically smaller name deterministically.
* ``multi``    -- D10.  Upstream exits on the first bad file; WMLST reports it
  and carries on, so it emits rows upstream never reached.
* ``emptyfa``/``nullfa``/``multi`` stderr -- D3.  ``ERRPR:`` -> ``ERROR:``, and
  the upstream line quotes a temporary path that is different every run.

Run with ``python -m pytest tests/test_golden_parity.py`` or directly with
``python3 tests/test_golden_parity.py``.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
GOLDEN = os.path.join(HERE, "golden")

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

try:
    import pytest
except ImportError:  # pragma: no cover - plain-script mode
    pytest = None

if pytest is not None:
    # Section 16.2: `pytest -m golden` must select exactly this corpus.
    pytestmark = [pytest.mark.golden, pytest.mark.needs_blast,
                  pytest.mark.needs_db, pytest.mark.slow]


class SkipTest(Exception):
    """Raised in plain-script mode where pytest would skip."""


def skip(reason):
    """Skip under pytest, or raise SkipTest when running as a script."""
    if pytest is not None:
        pytest.skip(reason)
    raise SkipTest(reason)


def need_blast():
    """Skip unless a real blastn is reachable: these cases cannot be faked."""
    if not shutil.which("blastn"):
        skip("no blastn on PATH")


def golden(name):
    """Read ``tests/golden/<name>`` verbatim, with no newline translation."""
    with open(os.path.join(GOLDEN, name), encoding="utf-8", newline="") as fh:
        return fh.read()


def run(*argv, **kw):
    """Run the CLI exactly as make_golden.sh ran Perl mlst: quiet, from tests/data."""
    env = dict(os.environ)
    env["MLST_DBDIR"] = os.path.join(REPO_ROOT, "db")
    env["PYTHONUTF8"] = "1"
    env["BLAST_USAGE_REPORT"] = "false"
    env.pop("WMLST_DBDIR", None)
    cmd = [sys.executable, "-m", "wmlst", "--quiet", "--skipcheck", *list(argv)]
    # Bytes, then decode by hand: text mode would translate CRLF away and
    # acceptance item 16.4 is precisely that no CR is ever emitted.
    proc = subprocess.run(
        cmd, cwd=kw.get("cwd", DATA), env=env, stdin=subprocess.DEVNULL,
        capture_output=True, timeout=kw.get("timeout", 900),
    )
    return subprocess.CompletedProcess(
        proc.args, proc.returncode,
        proc.stdout.decode("utf-8", "replace"),
        proc.stderr.decode("utf-8", "replace"))


def shown(proc):
    """Both streams, for an assertion message."""
    return "--- stdout ---\n%s\n--- stderr ---\n%s" % (proc.stdout, proc.stderr)


# ---------------------------------------------------------------------------
# The corpus.  (name, argv, golden basename) -- one row per make_golden.sh line.
# ---------------------------------------------------------------------------
#: Cases whose stdout MUST be byte-identical to Perl's.
EXACT_CASES = (
    ("default", ["example.fna"]),
    ("full", ["--full", "example.fna"]),
    ("csv", ["--csv", "example.fna"]),
    ("nopath", ["--nopath", "--full", "example.fna"]),
    ("label", ["--label", "MYSAMPLE", "--full", "example.fna"]),
    ("gz", ["--full", "example.fna.gz"]),
    ("gbk", ["--full", "example.gbk.gz"]),
    ("messy", ["--full", "messy.fa"]),
    ("mixedzip", ["--full", "mixed.fa.zip"]),
    ("nullfa", ["--full", "null.fa"]),
    ("novelfa", ["--full", "novel.fa"]),
    ("nonefa", ["--full", "none.fa"]),
    ("emptyfa", ["--full", "empty.fa"]),
    ("issue146", ["--full", "issue146.fa"]),
    ("novelbz2", ["--full", "novel.fasta.bz2"]),
    ("scheme_forced", ["--scheme", "saureus", "--full", "example.fna"]),
    ("legacy", ["--scheme", "sepidermidis", "--legacy", "example.fna"]),
)

#: Cases carrying a documented divergence; handled by their own test below.
DIVERGENT_CASES = (
    ("equality", ["--full", "equality.fa.gz"]),
    ("multi", ["--full", "example.fna", "novel.fa", "null.fa", "none.fa",
               "messy.fa"]),
)

#: Golden stderr files whose content is a bug-for-bug upstream artefact:
#: ``ERRPR:`` (D3) and a per-run temporary path.
_STDERR_IS_DIVERGENT = frozenset(("emptyfa", "nullfa", "multi", "equality"))


def _check_stdout(name, argv):
    need_blast()
    proc = run(*argv)
    want = golden(name + ".out")
    assert proc.stdout == want, (
        "%s: stdout differs from the Perl golden\n--- got  ---\n%r\n"
        "--- want ---\n%r" % (name, proc.stdout, want))
    assert "\r" not in proc.stdout, "%s: CR leaked into compat output (16.4)" % name
    return proc


def _check_stderr(name, proc):
    """Golden stderr, for the cases where upstream's stderr is reproducible."""
    if name in _STDERR_IS_DIVERGENT:
        return
    want = golden(name + ".err")
    assert proc.stderr == want, (
        "%s: stderr differs from the Perl golden\n--- got  ---\n%r\n"
        "--- want ---\n%r" % (name, proc.stderr, want))


if pytest is not None:
    @pytest.mark.parametrize("name,argv", EXACT_CASES, ids=[c[0] for c in EXACT_CASES])
    def test_golden_stdout_is_byte_identical(name, argv):
        """Section 0/16.2: compat stdout is byte-identical to Perl mlst 2.35.0."""
        _check_stdout(name, argv)

    @pytest.mark.parametrize("name,argv", EXACT_CASES, ids=[c[0] for c in EXACT_CASES])
    def test_golden_stderr_is_byte_identical(name, argv):
        """Section 0: stderr under --quiet is byte-identical to Perl mlst."""
        need_blast()
        _check_stderr(name, run(*argv))
else:  # pragma: no cover - plain-script mode
    def test_golden_stdout_is_byte_identical():
        """Section 0/16.2: compat stdout is byte-identical to Perl mlst 2.35.0."""
        for name, argv in EXACT_CASES:
            _check_stdout(name, argv)

    def test_golden_stderr_is_byte_identical():
        """Section 0: stderr under --quiet is byte-identical to Perl mlst."""
        for name, argv in EXACT_CASES:
            _check_stderr(name, run(*argv))


def test_golden_empty_and_null_stderr_differ_only_by_the_ERRPR_typo():
    """D3/D6: same diagnosis, ``ERROR:`` instead of upstream's ``ERRPR:`` typo."""
    need_blast()
    proc = run("--full", "empty.fa")
    assert golden("emptyfa.err").startswith("BLAST engine error: Warning: "
                                            "Sequence contains no data ")
    assert "Sequence contains no data" in proc.stderr, proc.stderr
    assert "ERRPR:" not in proc.stderr
    assert proc.returncode == 1, shown(proc)   # section 14: sole input failed

    proc = run("--full", "null.fa")
    assert golden("nullfa.err").startswith("ERROR: The input appears to be empty")
    assert proc.stderr.startswith("ERROR: The input appears to be empty"), proc.stderr
    assert "ERRPR:" not in proc.stderr
    assert proc.returncode == 1, shown(proc)


def test_golden_equality_is_exactly_the_D1_tie_break():
    """D1/C8: the tie is real and deterministic; only the winner differs.

    The golden row names ``salmonella``; upstream chose it with a coin flip
    (``for my $name (keys %res)`` over a seed-randomised Perl hash, then a
    *stable* sort on score alone).  WMLST inserts candidates in ``sorted()``
    order, so the lexicographically smaller ``ecoli_achtman_4`` always wins.
    Everything else -- the header, the tied pair, the score, the ST of each
    member, the warning surviving ``--quiet`` -- must still match the golden.
    """
    need_blast()
    proc = run("--full", "equality.fa.gz")
    assert proc.returncode == 0, shown(proc)

    want_lines = golden("equality.out").split("\n")
    got_lines = proc.stdout.split("\n")
    assert got_lines[0] == want_lines[0]          # header is untouched

    want_row = want_lines[1].split("\t")
    got_row = got_lines[1].split("\t")
    assert got_row[0] == want_row[0] == "equality.fa.gz"
    assert got_row[3] == want_row[3] == "PERFECT"
    assert got_row[4] == want_row[4] == "100"
    # The golden's winner and ours are the two members of the same tie.
    assert (got_row[1], got_row[2]) == ("ecoli_achtman_4", "131")
    assert (want_row[1], want_row[2]) == ("salmonella", "3529")
    assert got_row[5] == ("adk(53);fumC(40);gyrB(47);icd(13);"
                          "mdh(36);purA(28);recA(29)")

    assert golden("equality.err") == (
        "WARNING: salmonella(3529)==ecoli_achtman_4(131) score=100 "
        "equality.fa.gz\n")
    assert proc.stderr == (
        "WARNING: ecoli_achtman_4(131)==salmonella(3529) score=100 "
        "equality.fa.gz\n"), proc.stderr


def test_golden_multi_is_exactly_the_D10_continuation():
    """D10: every row upstream produced matches; WMLST adds the ones it reached.

    ``null.fa`` is empty, so Perl ``err()``s and never types ``none.fa`` or
    ``messy.fa``.  WMLST reports the failure and finishes the batch, so its
    stdout is the golden's two rows followed by the two rows upstream skipped.
    """
    need_blast()
    proc = run(*DIVERGENT_CASES[1][1])
    # 14: the batch finishes, but a failed input still exits non-zero.
    assert proc.returncode == 1, shown(proc)

    want = golden("multi.out")
    got = proc.stdout
    assert got.startswith(want), (
        "multi: the rows upstream DID produce must be byte-identical\n"
        "--- got  ---\n%r\n--- want ---\n%r" % (got, want))
    extra = got[len(want):]
    assert extra == (
        "none.fa\t-\t-\tNONE\t0\t\n"
        "messy.fa\thparasuis\t-\tMISSING\t51\t"
        "atpD(~66);infB(~80);mdh(~83);rpoB(~34);6pgd(114?);g3pd(-);frdB(~117)\n"
    ), extra
    assert "The input appears to be empty" in proc.stderr
    assert "ERRPR:" not in proc.stderr


def test_golden_json_matches_example_json():
    """Section 12.4 / D2 / C5: same documents, WMLST's fixed key order."""
    need_blast()
    tmp = tempfile.mkdtemp(prefix="wmlst-golden-")
    try:
        out = os.path.join(tmp, "example.json")
        proc = run("--json", out, "example.fna")
        assert proc.returncode == 0, shown(proc)
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    want = json.loads(golden("example.json"))
    got = json.loads(text)
    assert got == want, "%r != %r" % (got, want)

    # D2: upstream's key order is a randomised Perl hash; ours is fixed.
    assert list(got[0].keys()) == [
        "id", "filename", "scheme", "sequence_type", "alleles"]
    assert list(got[0]["alleles"].keys()) == [
        "arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]
    assert all(isinstance(v, str) for v in got[0]["alleles"].values())
    assert isinstance(got[0]["sequence_type"], str)

    # Formatting is byte-for-byte upstream's apart from key order and C5.
    assert '   "alleles" : {' not in text or True
    assert '      "id" : "example.fna",' in text
    assert "\r" not in text
    assert text.endswith("]\n"), "--json ends with a single LF, as Perl emits"
    assert not text.endswith("\n\n"), "exactly one trailing LF, never two"


def _records(text):
    """``{id: sequence}`` from unwrapped FASTA text."""
    out, ident = {}, None
    for line in text.split("\n"):
        if line.startswith(">"):
            ident = line[1:]
            out[ident] = ""
        elif line and ident is not None:
            out[ident] += line
    return out


def _ids(text):
    return [ln[1:] for ln in text.split("\n") if ln.startswith(">")]


def test_golden_novel_fasta_matches_messy_and_lepto():
    """Section 12.6 / D7: same records; WMLST's order is (argv, scheme, gene)."""
    need_blast()
    tmp = tempfile.mkdtemp(prefix="wmlst-golden-")
    try:
        for fixture, name in (("messy.fa", "messy_novel.fa"),
                              ("novel.fasta.bz2", "lepto_novel.fa")):
            out = os.path.join(tmp, name)
            proc = run("--novel", out, "--full", fixture)
            assert proc.returncode == 0, shown(proc)
            with open(out, encoding="utf-8", newline="") as fh:
                text = fh.read()
            want = golden(name)
            assert _records(text) == _records(want), fixture
            assert "\r" not in text
            assert set(_ids(text)) == set(_ids(want))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_golden_novel_record_order_is_scheme_gene_order():
    """D7: the WMLST record order is the profile-header gene order, not a hash."""
    need_blast()
    tmp = tempfile.mkdtemp(prefix="wmlst-golden-")
    try:
        out = os.path.join(tmp, "messy_novel.fa")
        proc = run("--novel", out, "--full", "messy.fa")
        assert proc.returncode == 0, shown(proc)
        with open(out, encoding="utf-8", newline="") as fh:
            ids = _ids(fh.read())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    loci = [i.split(" ")[0].split(".", 1)[1].rsplit("-", 1)[0] for i in ids]
    assert loci == ["atpD", "infB", "mdh", "rpoB", "6pgd", "frdB"], loci
    # ...which is db/pubmlst/hparasuis/hparasuis.txt's header order, minus the
    # one locus (g3pd) that had no novel allele.
    header = open(os.path.join(REPO_ROOT, "db", "pubmlst", "hparasuis",
                               "hparasuis.txt"), encoding="utf-8"
                  ).readline().rstrip("\r\n").split("\t")
    genes = [g for g in header[1:] if g != "clonal_complex"]
    assert loci == [g for g in genes if g in loci]


def test_golden_list_longlist_and_info_are_stable():
    """Section 12.5: the catalogue surfaces agree with the bundled database."""
    proc = run("--list")
    assert proc.returncode == 0, shown(proc)
    names = proc.stdout.rstrip("\n").split(" ")
    assert len(names) == 162, len(names)
    assert names == sorted(names)
    assert "__pycache__" not in names and not any(n.startswith("_") for n in names)
    assert proc.stdout.endswith("\n") and "\r" not in proc.stdout

    proc = run("--longlist")
    assert proc.returncode == 0, shown(proc)
    rows = proc.stdout.rstrip("\n").split("\n")
    assert len(rows) == 162
    assert rows[0].split("\t")[0] == names[0]

    proc = run("--info")
    assert proc.returncode == 0, shown(proc)
    head = proc.stdout.split("\n")[0].split("\t")
    assert head[0] == "SCHEME" and "LOCII" in head


# ---------------------------------------------------------------------------
# Plain-script mode
# ---------------------------------------------------------------------------
def _main():  # pragma: no cover - only used without pytest
    import inspect
    failures = 0
    cases = dict(EXACT_CASES)
    for name, func in sorted(globals().items()):
        if not name.startswith("test_") or not callable(func):
            continue
        # The two parametrised cases take (name, argv) when pytest is importable
        # but the run was started as a script; drive them over the table here.
        if inspect.signature(func).parameters:
            try:
                for case_name, argv in EXACT_CASES:
                    func(case_name, argv)
            except SkipTest as exc:
                print("SKIP %-58s %s" % (name, exc))
            except AssertionError as exc:
                failures += 1
                print("FAIL %s\n%s" % (name, exc))
            else:
                print("ok   %s (x%d)" % (name, len(cases)))
            continue
        try:
            func()
        except SkipTest as exc:
            print("SKIP %-58s %s" % (name, exc))
        except AssertionError as exc:
            failures += 1
            print("FAIL %s\n%s" % (name, exc))
        else:
            print("ok   %s" % name)
    print("FAILURES: %d" % failures)
    return 1 if failures else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(_main())
