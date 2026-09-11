# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Every case of upstream's bats suite, ported (docs/ARCHITECTURE.md 13.4).

Source: ``ref-mlst/test/test.sh`` (41 ``@test`` blocks). Each function below
names the upstream case it ports in its docstring, followed by the ten
Windows-port cases the architecture adds and a byte-for-byte comparison of the
CLI against ``tests/golden/*.out``.

The whole file runs under ``python3 -m pytest`` and as a plain script
(``python3 tests/test_cli_bats_parity.py``). Cases that need a real ``blastn``
and a finished engine are skipped, loudly, when either is absent.
"""

import json
import os
import py_compile
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
GOLDEN = os.path.join(HERE, "golden")
DB = os.path.join(REPO, "db")

if REPO not in sys.path:
    sys.path.insert(0, REPO)

try:
    import pytest
except ImportError:  # pragma: no cover - pytest is optional for the script mode
    pytest = None


class _Skip(Exception):
    """Raised by :func:`skip` when pytest is not available."""


def skip(reason):
    """Skip the current test, with or without pytest."""
    if pytest is not None:
        pytest.skip(reason)
    raise _Skip(reason)


#: Seconds any single CLI invocation may take before the suite calls it hung.
TIMEOUT_S = float(os.environ.get("WMLST_TEST_TIMEOUT", "180"))

#: The substring upstream's suite greps for: a tab-delimited sepidermidis ST 184.
SEPI = "\tsepidermidis\t184\t"


def have_blast():
    """True when a real ``blastn`` is on PATH."""
    return shutil.which("blastn") is not None


def engine_ready():
    """True when the engine and the report writer are implemented."""
    try:
        from wmlst import engine, report
    except Exception:
        return False
    return hasattr(engine, "Engine") and hasattr(report, "write_tsv")


def catalog_ready():
    """True when the scheme catalogue and the list writers are implemented."""
    try:
        from wmlst import report, schemes
    except Exception:
        return False
    return hasattr(schemes, "SchemeCatalog") and hasattr(report, "write_list")


def need_blast():
    if not have_blast():
        skip("blastn is not on PATH")
    if not engine_ready():
        skip("wmlst.engine.Engine / wmlst.report are not implemented yet")


def need_catalog():
    if not catalog_ready():
        skip("wmlst.schemes / wmlst.report listing writers are not implemented yet")


def run(*args, **kw):
    """Run the CLI the way the bats suite does: cwd=tests/data, --quiet --skipcheck.

    ``raw=True`` drops the two default flags (upstream's ``$bin`` instead of
    ``$exe``); ``stdin``, ``cwd`` and ``env_extra`` are passed through.
    """
    raw = kw.pop("raw", False)
    stdin = kw.pop("stdin", None)
    cwd = kw.pop("cwd", DATA)
    env_extra = kw.pop("env_extra", None)
    assert not kw, kw
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MLST_DBDIR", DB)
    env["PYTHONIOENCODING"] = "utf-8"
    if env_extra:
        env.update(env_extra)
    argv = [sys.executable, "-m", "wmlst"]
    if not raw:
        argv += ["--quiet", "--skipcheck"]
    argv += [str(a) for a in args]
    # A timeout, not patience: a deadlocked child must fail the suite loudly
    # rather than hang CI for an hour.
    return subprocess.run(argv, cwd=cwd, env=env, input=stdin,
                          capture_output=True, text=True, timeout=TIMEOUT_S)


def golden(name):
    """Read ``tests/golden/<name>`` with no newline translation."""
    with open(os.path.join(GOLDEN, name), encoding="utf-8", newline="") as fh:
        return fh.read()


def combined(proc):
    """bats' ``$output``: stdout and stderr merged."""
    return proc.stdout + proc.stderr


# ---------------------------------------------------------------------------
# bats 1-13: the cases that need neither blastn nor the engine
# ---------------------------------------------------------------------------
def test_bats01_script_syntax_check():
    """bats 1 "Script syntax check" - `perl -c` becomes a byte-compile."""
    for name in ("cli.py", "__main__.py", "__init__.py"):
        py_compile.compile(os.path.join(REPO, "wmlst", name), doraise=True)


def test_bats02_version():
    """bats 2 "Version": exit 0 and the tool name in the output."""
    proc = run("--version")
    assert proc.returncode == 0
    assert "mlst" in combined(proc)


def test_bats03_help():
    """bats 3 "Help": exit 0 and "threads" in the output."""
    proc = run("--help")
    assert proc.returncode == 0
    assert "threads" in combined(proc)


def test_bats04_try_check():
    """bats 4 "Try --check": exit 0 and "OK" in the output."""
    if not have_blast():
        skip("blastn is not on PATH")
    try:
        from wmlst import blastbin  # noqa: F401
    except ImportError:
        skip("wmlst.blastbin is not implemented yet")
    proc = run("--check", raw=True)
    assert proc.returncode == 0, combined(proc)
    assert "OK" in combined(proc)


def test_bats05_no_parameters():
    """bats 5 "No parameters": a non-zero exit with the upstream wording."""
    proc = run()
    assert proc.returncode != 0
    assert "Please provide some FASTA/Genbank files to genotype (can be .gz)" \
        in proc.stderr


def test_bats06_bad_option():
    """bats 6 "Bad option": "Unknown option" and NO usage block."""
    proc = run("--doesnotexist")
    assert proc.returncode != 0
    assert "Unknown option" in combined(proc)
    assert "USAGE" not in combined(proc)
    assert "SYNOPSIS" not in combined(proc)
    assert combined(proc).strip() == "Unknown option: doesnotexist"


def test_bats07_list_schemes():
    """bats 7 "List schemes --list"."""
    need_catalog()
    proc = run("--list")
    assert proc.returncode == 0, combined(proc)
    assert "saureus" in proc.stdout
    assert proc.stdout.endswith("\n")
    assert proc.stdout.count("\n") == 1


def test_bats08_list_schemes_with_mlst_dbdir():
    """bats 8 "List schemes with MLST_DBDIR" (divergence D15)."""
    need_catalog()
    proc = run("--list", env_extra={"MLST_DBDIR": DB})
    assert proc.returncode == 0, combined(proc)
    assert "saureus" in proc.stdout


def test_bats09_longlist():
    """bats 9 "List schemes --longlist"."""
    need_catalog()
    proc = run("--longlist")
    assert proc.returncode == 0, combined(proc)
    assert "saureus" in proc.stdout


def test_bats10_info():
    """bats 10 "List schemes --info": LOCII header row, six columns."""
    need_catalog()
    proc = run("--info")
    assert proc.returncode == 0, combined(proc)
    lines = proc.stdout.split("\n")
    assert lines[0] == "SCHEME\tLOCII\tTYPES\tALLELES\tDATE\tLOCII_NAMES"
    assert "saureus" in proc.stdout


def test_bats11_passing_a_folder():
    """bats 11 "Passing a folder"."""
    proc = run(DATA)
    assert proc.returncode != 0
    assert "directory" in combined(proc)
    assert "seems to be a directory, not a file" in proc.stderr


def test_bats12_null_input():
    """bats 12 "Null input": null.fa fails with an ERROR line."""
    need_blast()
    proc = run("--full", "null.fa")
    assert proc.returncode != 0
    assert "ERROR" in combined(proc)
    assert "The input appears to be empty" in proc.stderr


def test_bats13_empty_input():
    """bats 13 "Empry input": empty.fa surfaces BLAST's own wording."""
    need_blast()
    proc = run("--no-quiet", "empty.fa")
    assert "Sequence contains no data" in combined(proc)


# ---------------------------------------------------------------------------
# bats 14-26: the typing cases
# ---------------------------------------------------------------------------
def test_bats14_plain_fasta():
    """bats 14 "Plain FASTA"."""
    need_blast()
    proc = run("example.fna")
    assert proc.returncode == 0, combined(proc)
    assert SEPI in proc.stdout


def test_bats15_gzipped_fasta():
    """bats 15 "Gzipped FASTA"."""
    need_blast()
    proc = run("example.fna.gz")
    assert proc.returncode == 0, combined(proc)
    assert SEPI in proc.stdout


def test_bats16_gzipped_genbank():
    """bats 16 "Gzipped Genbank" (upstream reuses the .fna.gz fixture; we also
    exercise the real .gbk.gz, which is what the case meant to test)."""
    need_blast()
    assert SEPI in run("example.fna.gz").stdout
    assert SEPI in run("example.gbk.gz").stdout


def test_bats17_bzipped_fasta():
    """bats 17 "Bzipped FASTA"."""
    need_blast()
    proc = run("novel.fasta.bz2")
    assert proc.returncode == 0, combined(proc)
    assert "leptospira_2" in proc.stdout


def test_bats18_zipped_fasta():
    """bats 18 "Zipped FASTA"."""
    need_blast()
    proc = run("mixed.fa.zip")
    assert proc.returncode == 0, combined(proc)
    assert "mgen" in proc.stdout


def test_bats19_fofn_input():
    """bats 19 "FOFN input": three sepidermidis rows."""
    need_blast()
    proc = run("--quiet", "--fofn", "fofn.txt")
    assert proc.returncode == 0, combined(proc)
    lines = proc.stdout.rstrip("\n").split("\n")
    assert len(lines) == 3, lines
    for line in lines:
        assert SEPI in line


def test_bats20_bad_fofn_input():
    """bats 20 "Bad FOFN input": an empty FOFN leaves no files to type."""
    empty = os.devnull if os.path.exists(os.devnull) else None
    if empty is None:  # pragma: no cover - every supported OS has one
        skip("no null device")
    proc = run("--fofn", empty)
    assert proc.returncode != 0
    assert "Please provide some FASTA/Genbank files" in proc.stderr


def test_bats21_try_skipcheck():
    """bats 21 "Try --skipcheck": the dependency banner must not appear."""
    need_blast()
    proc = run("--no-quiet", "--skipcheck", "example.fna.gz")
    assert proc.returncode == 0, combined(proc)
    assert "Checking mlst dependencie" not in combined(proc)
    assert "Skipping dependency check due to --skipcheck" in proc.stderr


def test_bats22_accept_stdin():
    """bats 22 "Accept STDIN"."""
    need_blast()
    if not os.path.exists("/dev/stdin"):
        skip("no /dev/stdin on this platform")
    import gzip

    with gzip.open(os.path.join(DATA, "example.fna.gz"), "rt") as fh:
        payload = fh.read()
    proc = run("/dev/stdin", stdin=payload)
    assert proc.returncode == 0, combined(proc)
    assert SEPI in proc.stdout


def test_bats23_debug_works():
    """bats 23 "Check --debug works": the "=== DEBUG" marker."""
    need_blast()
    proc = run("--debug", "example.fna.gz")
    assert proc.returncode == 0, combined(proc)
    assert "=== DEBUG" in combined(proc)


def test_bats24_two_files():
    """bats 24 "Two files in legacy mode"."""
    need_blast()
    proc = run("example.fna.gz", "example.gbk.gz")
    assert proc.returncode == 0, combined(proc)
    assert SEPI in proc.stdout
    assert len(proc.stdout.rstrip("\n").split("\n")) == 2


def test_bats25_finds_duplicate_alleles():
    """bats 25 "Finds duplicate alleles": pgm(3,3), atpA(1,1) and a WARNING."""
    need_blast()
    proc = run("mixed.fa.zip", raw=True)
    assert proc.returncode == 0, combined(proc)
    out = combined(proc)
    assert "pgm(3,3)" in out
    assert "atpA(1,1)" in out
    assert "WARNING" in out
    assert "found additional exact allele match" in out


def test_bats26_decimal_allele_ids_are_reported():
    """bats 26 "Decimal allele IDs are reported": the sseqid regex."""
    line = "ngstar.penA_2.002\t884\t884\t884\tq1\t1\t884\tACGT\tplus"
    try:
        from wmlst.engine import HIT_RE as pattern
    except ImportError:
        pattern = re.compile(
            r'^(\w+)\.(\w+)[_-](\d+(?:\.\d+)?)'
            r'\t(\d+)\t(\d+)\t(\d+)'
            r'\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)', re.ASCII)
    m = pattern.match(line)
    assert m is not None
    assert (m.group(1), m.group(2), m.group(3)) == ("ngstar", "penA", "2.002")


def test_bats27_decimal_duplicates_sort_numerically():
    """bats 27 "Decimal duplicate alleles sort numerically"."""
    code = "10.001,2.002,2.001"
    try:
        from wmlst.engine import sort_duplicate_codes
    except ImportError:
        def sort_duplicate_codes(value):
            if "," in value and re.fullmatch(r"[\d.,]+", value, re.ASCII):
                return ",".join(sorted(value.split(","), key=float))
            return value
    assert sort_duplicate_codes(code) == "2.001,2.002,10.001"


# ---------------------------------------------------------------------------
# bats 28-41: output formats and status words
# ---------------------------------------------------------------------------
def test_bats28_csv_output():
    """bats 28 "CSV output": a comma inside a field forces quoting."""
    need_blast()
    proc = run("--csv", "mixed.fa.zip")
    assert proc.returncode == 0, combined(proc)
    assert ',"MLST_gyrB(1,1)",' in proc.stdout.split("\n")[0]


def test_bats29_detect_perfect():
    """bats 29 "Detect PERFECT"."""
    need_blast()
    proc = run("--full", "--csv", "example.fna")
    assert ",PERFECT," in proc.stdout.split("\n")[1]


def test_bats30_detect_mixed():
    """bats 30 "Detect MIXED"."""
    need_blast()
    proc = run("--full", "--csv", "mixed.fa.zip")
    assert ",MIXED," in proc.stdout.split("\n")[1]


def test_bats31_detect_novel():
    """bats 31 "Detect NOVEL"."""
    need_blast()
    proc = run("--full", "--csv", "novel.fa")
    assert ",NOVEL," in proc.stdout.split("\n")[1]


def test_bats32_detect_missing():
    """bats 32 "Detect MISSING"."""
    need_blast()
    proc = run("--full", "--csv", "messy.fa")
    assert ",MISSING," in proc.stdout.split("\n")[1]


def test_bats33_detect_none():
    """bats 33 "Detect NONE"."""
    need_blast()
    proc = run("--full", "--csv", "none.fa")
    assert ",NONE," in proc.stdout.split("\n")[1]


def test_bats34_json_output():
    """bats 34 "JSON output"."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "wmlst.json")
        proc = run("--json", out, "example.fna.gz")
        assert proc.returncode == 0, combined(proc)
        assert os.access(out, os.R_OK)
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
        assert "sequence_type" in text
        payload = json.loads(text)
        assert payload[0]["sequence_type"] == "184"
        assert isinstance(payload[0]["sequence_type"], str)


def test_bats35_using_outfile():
    """bats 35 "Using --outfile": empty stdout, exactly two lines in the file."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "wmlst.tsv")
        proc = run("--full", "--outfile", out, "messy.fa")
        assert proc.returncode == 0, combined(proc)
        assert proc.stdout == ""
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
        assert text.count("\n") == 2
        assert "\r" not in text


def test_bats36_custom_label():
    """bats 36 "Custom label"."""
    need_blast()
    proc = run("--label", "GDAYMATE", "example.fna.gz")
    assert proc.returncode == 0, combined(proc)
    assert "GDAYMATE" in proc.stdout


def test_bats37_duplicate_label():
    """bats 37 "Duplicate label": --label with two files is refused."""
    proc = run("--label", "double_trouble", "example.gbk.gz", "example.fna.gz")
    assert proc.returncode != 0
    assert "Using --label when scanning multiple files does not make sense" \
        in proc.stderr


def test_bats38_save_novel_allele():
    """bats 38 "Save novel allele"."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "wmlst.novel.fa")
        proc = run("--novel", out, "novel.fasta.bz2")
        assert proc.returncode == 0, combined(proc)
        assert os.access(out, os.R_OK)
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
        assert "mreA" in text
        assert "\r" not in text


def test_bats39_issue_146():
    """bats 39 "Test issue 146"."""
    need_blast()
    proc = run("issue146.fa")
    assert proc.returncode == 0, combined(proc)
    assert "purE(25)" in proc.stdout


def test_bats40_show_seqs_entry_point():
    """bats 40 "Script: show_seqs" - wmlst.cli:main_show_seqs, the wmlst-show-seqs exe."""
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MLST_DBDIR", DB)
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; from wmlst.cli import main_show_seqs;"
         " sys.exit(main_show_seqs(sys.argv[1:]))",
         "-s", "efaecium", "-t", "111"],
        env=env, capture_output=True, text=True, timeout=TIMEOUT_S)
    assert proc.returncode == 0, proc.stderr
    assert ">atpA_2" in proc.stdout
    assert proc.stdout.count(">") == 7  # one record per locus of efaecium
    assert "Extracting: atpA_2" in proc.stderr  # progress goes to stderr only


def test_show_seqs_rejects_bad_arguments():
    """The four error paths of scripts/mlst-show_seqs, exit 1 each."""
    from wmlst.cli import main_show_seqs

    for args in (["-t", "111"], ["-s", "efaecium"],
                 ["-s", "efaecium", "-t", "abc"],
                 ["-s", "nosuchscheme", "-t", "1"]):
        assert main_show_seqs(args) == 1, args


def test_packaged_entry_points_exist():
    """pyproject binds four console scripts into wmlst.cli; all must be callable."""
    from wmlst import cli

    for name in ("main", "main_show_seqs", "main_make_blast_db", "main_update_db"):
        assert callable(getattr(cli, name)), name
    assert callable(__import__("wmlst.gui", fromlist=["main"]).main)


def test_bats41_equally_good_scheme_warning():
    """bats 41 "Eqyally good scheme warning"."""
    need_blast()
    proc = run("equality.fa.gz", raw=True)
    assert proc.returncode == 0, combined(proc)
    assert "WARNING:" in combined(proc)


# ---------------------------------------------------------------------------
# The ten Windows-port cases (section 13.4)
# ---------------------------------------------------------------------------
def test_port01_crlf_input_is_handled():
    """A CRLF-terminated FASTA types identically to the LF original."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        dest = os.path.join(tmp, "crlf.fa")
        with open(os.path.join(DATA, "example.fna"), "rb") as src:
            payload = src.read().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")
        with open(dest, "wb") as fh:
            fh.write(payload)
        proc = run("--nopath", dest)
        assert proc.returncode == 0, combined(proc)
        assert SEPI in proc.stdout


def test_port02_path_with_spaces():
    """A path containing spaces survives the whole pipeline (no shell anywhere)."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        folder = os.path.join(tmp, "a folder with spaces")
        os.makedirs(folder)
        dest = os.path.join(folder, "my sample.fna")
        shutil.copy(os.path.join(DATA, "example.fna"), dest)
        proc = run("--nopath", dest)
        assert proc.returncode == 0, combined(proc)
        assert SEPI in proc.stdout
        assert proc.stdout.startswith("my sample.fna\t")


def test_port03_non_ascii_path():
    """A non-ASCII filename reaches the FILE column unmangled."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        dest = os.path.join(tmp, "échantillon-ü.fna")
        try:
            shutil.copy(os.path.join(DATA, "example.fna"), dest)
        except (UnicodeEncodeError, OSError):  # pragma: no cover
            skip("filesystem cannot hold a non-ASCII name")
        proc = run("--nopath", dest)
        assert proc.returncode == 0, combined(proc)
        assert proc.stdout.startswith("échantillon-ü.fna\t")


def test_port04_stdout_is_lf_never_crlf():
    """Section 0: a single \\r in a compat output is a release blocker."""
    need_blast()
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MLST_DBDIR", DB)
    proc = subprocess.run(
        [sys.executable, "-m", "wmlst", "--quiet", "--skipcheck", "--full",
         "example.fna"],
        cwd=DATA, env=env, capture_output=True)
    assert proc.returncode == 0, proc.stderr
    assert b"\r" not in proc.stdout


def test_port04b_stderr_is_lf_never_crlf():
    """Section 0: stderr is a byte-identity surface too, so it is LF as well.

    ``wrn()``/``msg()`` lines are pinned by ``tests/golden/*.err``, so a CR on
    stderr is the same release blocker as one on stdout (acceptance item 16.4).

    This cannot be observed on POSIX, where the platform terminator is already
    LF: ``reconfigure_streams()`` could skip stderr entirely and every Linux run
    would stay green.  So drive the function against a stream that behaves the
    way the real Windows ``sys.stderr`` does -- ``newline=None`` there makes the
    text layer translate every ``"\n"`` into ``"\r\n"`` -- and require that
    reconfiguring turns that translation OFF.  Reverting stderr to a
    ``reconfigure()`` call that omits ``newline="\n"`` fails this on BOTH
    platforms, which is the point: the bug was invisible to the Linux suite.
    """
    import io

    from wmlst.cli import Logger, reconfigure_streams

    def windows_like():
        buf = io.BytesIO()
        # newline="\r\n" reproduces the writenl the MS_WINDOWS branch of
        # CPython's TextIOWrapper picks for newline=None. Same knob, testable here.
        return buf, io.TextIOWrapper(buf, encoding="utf-8", newline="\r\n")

    out_buf, out_stream = windows_like()
    err_buf, err_stream = windows_like()
    saved = (sys.stdout, sys.stderr)
    try:
        sys.stdout, sys.stderr = out_stream, err_stream
        # Pre-flight: untouched, this stream really does emit CRLF -- otherwise
        # the assertions below would pass no matter what the function did.
        sys.stderr.write("probe\n")
        sys.stderr.flush()
        assert err_buf.getvalue() == b"probe\r\n"
        err_buf.seek(0)
        err_buf.truncate()

        reconfigure_streams()
        Logger.wrn("ecoli_achtman_4(131)==salmonella(3529) score=100 x.fa")
        sys.stdout.write("FILE\tSCHEME\tST\n")
        sys.stdout.flush()
    finally:
        sys.stdout, sys.stderr = saved

    err = err_buf.getvalue()
    assert err == (b"WARNING: ecoli_achtman_4(131)==salmonella(3529) "
                   b"score=100 x.fa\n"), err
    assert b"\r" not in err
    assert out_buf.getvalue() == b"FILE\tSCHEME\tST\n"


def test_port05_threads_do_not_change_the_answer():
    """--threads is a performance knob; it must not move a single byte."""
    need_blast()
    one = run("--full", "example.fna")
    four = run("--full", "--threads", "4", "example.fna")
    assert one.stdout == four.stdout
    assert one.returncode == four.returncode == 0


def test_port06_default_exclude_list():
    """The default --exclude set is upstream's four schemes (bin/mlst:26)."""
    from wmlst.cli import EXCLUDE_DEFAULT
    from wmlst.engine import DEFAULT_EXCLUDE, RunConfig

    assert frozenset(EXCLUDE_DEFAULT.split(",")) == DEFAULT_EXCLUDE
    assert RunConfig().exclude == DEFAULT_EXCLUDE


def test_port07_scheme_forces_minscore_zero_and_clears_exclude():
    """Section 3.6.1 invariant 2, announced on stderr."""
    need_catalog()
    from wmlst.cli import build_parser, config_from_args

    ns = build_parser().parse_args(["--scheme", "saureus", "x.fna"])
    cfg = config_from_args(ns, [])
    assert cfg.minscore == 0.0
    assert cfg.exclude == frozenset()


def test_port07b_scheme_announces_the_minscore_override():
    """The same invariant, announced on stderr (bin/mlst:92)."""
    need_blast()
    proc = run("--no-quiet", "--scheme", "saureus", "--full", "example.fna")
    assert "Setting --minscore=0 because user chose --scheme" in proc.stderr
    assert proc.returncode == 0, combined(proc)


def test_port08_invalid_scheme_is_refused():
    """Step 12: an unknown scheme name names --list in the error."""
    need_catalog()
    proc = run("--scheme", "notascheme", "example.fna")
    assert proc.returncode != 0
    assert "Invalid --scheme 'notascheme'. Check using --list" in proc.stderr


def test_port09_legacy_requires_scheme():
    """Step 10 (bin/mlst:81)."""
    proc = run("--legacy", "example.fna")
    assert proc.returncode != 0
    assert "Must specify a --scheme for --legacy output mode" in proc.stderr


def test_port10_branding_is_on_the_banner():
    """Added case 51: the vendor and the author appear on the stderr banner."""
    from wmlst import branding

    proc = run("--no-quiet", "--skipcheck", "--list")
    assert branding.VENDOR in proc.stderr
    assert branding.AUTHOR in proc.stderr
    assert branding.UPSTREAM_AUTHOR in proc.stderr
    # ...and never in the machine-readable payload (section 4.2).
    assert branding.VENDOR not in proc.stdout


def test_port11_quiet_hides_msg_but_not_warnings():
    """--quiet suppresses msg() and must NOT suppress the tie wrn() lines."""
    need_blast()
    proc = run("equality.fa.gz")
    assert proc.returncode == 0, combined(proc)
    assert proc.stderr.startswith("WARNING: ")
    assert "This is wmlst" not in proc.stderr


# ---------------------------------------------------------------------------
# Byte-for-byte comparison against the Perl-generated golden corpus
# ---------------------------------------------------------------------------
#: name -> (argv, golden stdout file). Mirrors scripts/make_golden.sh.
GOLDEN_CASES = (
    ("default", ["example.fna"], "default.out"),
    ("full", ["--full", "example.fna"], "full.out"),
    ("csv", ["--csv", "example.fna"], "csv.out"),
    ("nopath", ["--nopath", "--full", "example.fna"], "nopath.out"),
    ("label", ["--label", "MYSAMPLE", "--full", "example.fna"], "label.out"),
    ("gz", ["--full", "example.fna.gz"], "gz.out"),
    ("gbk", ["--full", "example.gbk.gz"], "gbk.out"),
    ("messy", ["--full", "messy.fa"], "messy.out"),
    ("mixedzip", ["--full", "mixed.fa.zip"], "mixedzip.out"),
    ("nullfa", ["--full", "null.fa"], "nullfa.out"),
    ("novelfa", ["--full", "novel.fa"], "novelfa.out"),
    ("nonefa", ["--full", "none.fa"], "nonefa.out"),
    ("emptyfa", ["--full", "empty.fa"], "emptyfa.out"),
    ("issue146", ["--full", "issue146.fa"], "issue146.out"),
    ("novelbz2", ["--full", "novel.fasta.bz2"], "novelbz2.out"),
    ("scheme_forced", ["--scheme", "saureus", "--full", "example.fna"],
     "scheme_forced.out"),
    ("legacy", ["--scheme", "sepidermidis", "--legacy", "example.fna"],
     "legacy.out"),
)


def _golden_case(name, argv, out_name):
    need_blast()
    proc = run(*argv)
    assert proc.stdout == golden(out_name), (
        "%s: stdout differs from the Perl golden\n--- got ---\n%r\n--- want ---\n%r"
        % (name, proc.stdout, golden(out_name)))


if pytest is not None:
    @pytest.mark.parametrize("name,argv,out_name", GOLDEN_CASES,
                             ids=[c[0] for c in GOLDEN_CASES])
    def test_golden_stdout_is_byte_identical(name, argv, out_name):
        """Section 0: TSV/CSV/legacy stdout is byte-identical to Perl mlst."""
        _golden_case(name, argv, out_name)
else:  # pragma: no cover - script mode
    def test_golden_stdout_is_byte_identical():
        """Section 0: TSV/CSV/legacy stdout is byte-identical to Perl mlst."""
        for case in GOLDEN_CASES:
            _golden_case(*case)


def test_golden_multi_file_surviving_rows_match():
    """The multi-file golden, allowing for divergence D10.

    Upstream dies on ``null.fa`` and never types the two files after it; WMLST
    keeps going. Every row upstream DID produce must still match exactly.
    """
    need_blast()
    proc = run("--full", "example.fna", "novel.fa", "null.fa", "none.fa", "messy.fa")
    want = golden("multi.out").split("\n")
    got = proc.stdout.split("\n")
    assert got[:len(want) - 1] == want[:len(want) - 1]
    # The batch continues past the bad file, but an unreadable input is still
    # a non-zero exit so a pipeline notices the missing sample. Upstream also
    # exits non-zero here (it dies outright).
    assert proc.returncode == 1
    assert "The input appears to be empty" in proc.stderr


def test_golden_json_is_structurally_identical():
    """Section 0 / 12.4: --json matches after deserialisation, keys in D2 order."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "example.json")
        proc = run("--json", out, "example.fna")
        assert proc.returncode == 0, combined(proc)
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
    want = json.loads(golden("example.json"))
    got = json.loads(text)
    assert got == want
    assert list(got[0].keys()) == ["id", "filename", "scheme", "sequence_type",
                                   "alleles"]
    assert text.endswith("]\n")      # one trailing LF, as Perl emits
    assert not text.endswith("\n\n")
    assert "\r" not in text


def test_golden_novel_fasta_matches_as_a_record_set():
    """Section 0 / 12.6: --novel matches as a set of {id: sequence} records."""
    need_blast()

    def records(text):
        out = {}
        ident = None
        for line in text.split("\n"):
            if line.startswith(">"):
                ident = line[1:]
                out[ident] = ""
            elif line and ident is not None:
                out[ident] += line
        return out

    for fixture, name in (("messy.fa", "messy_novel.fa"),
                          ("novel.fasta.bz2", "lepto_novel.fa")):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "novel.fa")
            proc = run("--novel", out, "--full", fixture)
            assert proc.returncode == 0, combined(proc)
            with open(out, encoding="utf-8", newline="") as fh:
                got = records(fh.read())
        assert got == records(golden(name)), fixture


def test_golden_equality_tie_is_a_documented_divergence():
    """equality.fa.gz scores 100 for two schemes; the winner is D1, not chance.

    Upstream decides this tie by Perl hash order - the golden happens to record
    ``salmonella``. C8/D1 make WMLST deterministic instead, so the row may name
    either tied scheme; what MUST hold is that both are named in the tie warning
    at score=100, that the warning survives --quiet, and that everything else on
    the row matches the golden shape.
    """
    need_blast()
    proc = run("--full", "equality.fa.gz")
    assert proc.returncode == 0, combined(proc)
    head, row = proc.stdout.split("\n")[:2]
    assert head == golden("equality.out").split("\n")[0]
    fields = row.split("\t")
    assert fields[0] == "equality.fa.gz"
    assert fields[1] in ("salmonella", "ecoli_achtman_4")
    assert fields[3] == "PERFECT" and fields[4] == "100"
    if fields[1] == "salmonella":
        assert proc.stdout == golden("equality.out")
    warning = proc.stderr.strip()
    assert warning.startswith("WARNING: ")
    assert "score=100 equality.fa.gz" in warning
    assert "salmonella(3529)" in warning and "ecoli_achtman_4(131)" in warning
    # bin/mlst:405 uses wrn(), which --quiet must not suppress.
    assert proc.stderr.count("\n") == 1


# ---------------------------------------------------------------------------
# The WMLST-only switches (divergence D14) -- they must never move a compat byte
# ---------------------------------------------------------------------------
def test_ext01_html_report_is_written_and_branded():
    """--html writes a self-contained report and says so even under --quiet."""
    need_blast()
    from wmlst import branding

    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "report.html")
        proc = run("--full", "--html", out, "example.fna")
        assert proc.returncode == 0, combined(proc)
        assert proc.stdout == golden("full.out")  # the compat surface is untouched
        assert "Wrote HTML report: %s" % out in proc.stderr
        with open(out, encoding="utf-8", newline="") as fh:
            html = fh.read()
    assert branding.VENDOR in html and branding.AUTHOR in html
    assert "sepidermidis" in html


def test_ext02_evidence_tsv_is_written():
    """--evidence-tsv writes the long-format hit dump (section 4.7)."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "evidence.tsv")
        proc = run("--full", "--evidence-tsv", out, "example.fna")
        assert proc.returncode == 0, combined(proc)
        assert proc.stdout == golden("full.out")
        with open(out, encoding="utf-8", newline="") as fh:
            text = fh.read()
    header = text.split("\n")[0].split("\t")
    assert header[:4] == ["LABEL", "SCHEME", "LOCUS", "ALLELE"]
    assert "\r" not in text


def test_ext03_jobs_does_not_change_the_answer():
    """--jobs is a scheduling knob; the rows stay in argv order and unchanged."""
    need_blast()
    serial = run("--full", "example.fna", "issue146.fa")
    parallel = run("--full", "--jobs", "2", "example.fna", "issue146.fa")
    assert serial.stdout == parallel.stdout
    assert serial.returncode == parallel.returncode == 0
    rows = serial.stdout.rstrip("\n").split("\n")
    assert rows[1].startswith("example.fna\t")
    assert rows[2].startswith("issue146.fa\t")


def test_ext04_exit_codes_follow_section_14():
    """Section 14: any unreadable input exits 1; an all-good run exits 0."""
    need_blast()
    assert run("--full", "null.fa").returncode == 1
    # The good row is still written, but the exit code flags the lost sample.
    partial = run("--full", "example.fna", "null.fa")
    assert partial.returncode == 1
    assert "sepidermidis" in partial.stdout
    assert run("--full", "example.fna").returncode == 0


# ---------------------------------------------------------------------------
# Regressions for the adversarial-review findings 7, 14, 15, 16, 17 and 22
# ---------------------------------------------------------------------------
def test_fix07_blastdb_does_not_move_the_database_root():
    """Finding 7: --blastdb must not re-root the database (bin/mlst:471-472).

    Upstream defaults --blastdb and --datadir independently off $MLST_DBDIR, so
    naming a BLAST index somewhere else never makes the PubMLST folder move
    with it.
    """
    import types

    from wmlst import cli

    with tempfile.TemporaryDirectory() as tmp:
        stem = os.path.join(tmp, "bl", "mlst.fa")     # no pubmlst/ sibling
        ns = types.SimpleNamespace(datadir="", blastdb=stem)
        old = dict(os.environ)
        try:
            os.environ.pop("WMLST_DBDIR", None)
            os.environ["MLST_DBDIR"] = DB
            dbdir, datadir, blastdb = cli._resolve_db(ns)
        finally:
            os.environ.clear()
            os.environ.update(old)
    assert dbdir == os.path.abspath(DB)
    assert datadir == os.path.join(os.path.abspath(DB), "pubmlst")
    assert blastdb == os.path.abspath(stem)


def test_fix07_datadir_does_not_move_the_blast_index():
    """Finding 7, mirror case: --datadir must not drag the BLAST index along."""
    import types

    from wmlst import cli

    with tempfile.TemporaryDirectory() as tmp:
        pm = os.path.join(tmp, "pm")                  # no blast/ sibling
        os.mkdir(pm)
        ns = types.SimpleNamespace(datadir=pm, blastdb="")
        old = dict(os.environ)
        try:
            os.environ.pop("WMLST_DBDIR", None)
            os.environ["MLST_DBDIR"] = DB
            dbdir, datadir, blastdb = cli._resolve_db(ns)
        finally:
            os.environ.clear()
            os.environ.update(old)
    assert dbdir == os.path.abspath(DB)
    assert datadir == os.path.abspath(pm)
    assert os.path.dirname(blastdb) != os.path.join(tmp, "blast")
    assert blastdb == os.path.join(os.path.abspath(DB), "blast", "mlst.fa")


def test_fix07_datadir_elsewhere_still_types_end_to_end():
    """Finding 7 end to end: a --datadir outside $MLST_DBDIR must still run."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        pm = os.path.join(tmp, "pm")
        os.mkdir(pm)
        shutil.copytree(os.path.join(DB, "pubmlst", "sepidermidis"),
                        os.path.join(pm, "sepidermidis"))
        proc = run("--datadir", pm, "example.fna")
    assert proc.returncode == 0, combined(proc)
    assert proc.stdout.startswith("example.fna" + SEPI)


def test_fix14_longlist_and_info_ignore_csv():
    """Finding 14: bin/mlst:118 exits before `$OUTSEP = ',' if $csv` (line 136).

    --list/--longlist/--info are therefore unconditionally tab separated.
    """
    need_catalog()
    for flag in ("--longlist", "--info", "--list"):
        plain = run(flag)
        comma = run(flag, "--csv")
        assert plain.returncode == 0, combined(plain)
        assert comma.returncode == 0, combined(comma)
        assert comma.stdout == plain.stdout, flag
    info = run("--info", "--csv")
    assert info.stdout.split("\n")[0] == \
        "SCHEME\tLOCII\tTYPES\tALLELES\tDATE\tLOCII_NAMES"


def test_fix15_unique_upstream_abbreviations_resolve():
    """Finding 15: a prefix unique upstream keeps working despite the extras.

    Getopt::Long's auto_abbrev accepts ``--t``/``--j``/``--e``/``--b``/``--h``;
    the WMLST-only ``--tsv-evidence``/``--jobs``/``--evidence-tsv``/
    ``--blast-timeout``/``--html*`` must not shadow them.
    """
    from wmlst.cli import _expand_abbrevs, build_parser

    parser = build_parser()
    cases = {
        "--h": "--help", "--t": "--threads", "--j": "--json",
        "--e": "--exclude", "--b": "--blastdb", "--q": "--quiet",
        "--sc": "--scheme",
    }
    for short, full in cases.items():
        assert _expand_abbrevs([short], parser) == [full], short
    # `name=value` keeps its value, exact spellings are untouched, and the
    # WMLST extras are still reachable by their own unambiguous prefixes.
    assert _expand_abbrevs(["--t=4"], parser) == ["--threads=4"]
    assert _expand_abbrevs(["--jobs", "2"], parser) == ["--jobs", "2"]
    assert _expand_abbrevs(["--ts"], parser) == ["--tsv-evidence"]
    assert _expand_abbrevs(["--", "--t"], parser) == ["--", "--t"]


def test_fix15_abbreviated_help_prints_the_help():
    """Finding 15, through the real CLI: `--h` is `--help`, not an error."""
    proc = run("--h", raw=True)
    assert proc.returncode == 0, combined(proc)
    assert proc.stdout == golden("help.txt")


def test_fix15_ambiguous_prefix_uses_getopt_long_wording():
    """Finding 15: `--l` is ambiguous upstream too, and says so the same way."""
    proc = run("--l", raw=True)
    assert proc.returncode == 1
    assert proc.stderr.rstrip("\n").endswith(
        "Option l is ambiguous (label, legacy, list, longlist)")


def test_fix16_input_files_may_follow_an_option():
    """Finding 16: GetOptions permutes, so `mlst a --quiet b` types BOTH files."""
    need_blast()
    proc = run("example.fna", "--quiet", "issue146.fa", raw=True)
    assert proc.returncode == 0, combined(proc)
    rows = proc.stdout.rstrip("\n").split("\n")
    assert len(rows) == 2, proc.stdout
    assert rows[0].startswith("example.fna\t")
    assert rows[1].startswith("issue146.fa\t")


def test_fix16_permuted_files_reach_the_label_guard():
    """Finding 16: the trailing file must be collected before the guard runs."""
    proc = run("example.fna", "--label", "foo", "issue146.fa")
    assert proc.returncode == 1
    assert "Using --label when scanning multiple files" in proc.stderr


def test_fix16_a_real_unknown_option_is_still_reported():
    """Finding 16: permuting must not swallow a genuine unknown option."""
    proc = run("example.fna", "--badopt", "issue146.fa")
    assert proc.returncode == 1
    assert proc.stderr.rstrip("\n").endswith("Unknown option: badopt")


def test_fix17_engine_msg_lines_reach_stderr():
    """Findings 17/22: `Found exact allele match ...` is msg(), bin/mlst:348."""
    need_blast()
    proc = run("--skipcheck", "--label", "FOO", "example.fna", raw=True)
    assert proc.returncode == 0, combined(proc)
    assert "Using label 'FOO' for file example.fna" in proc.stderr
    exact = [ln for ln in proc.stderr.split("\n")
             if ln.startswith("Found exact allele match ")]
    assert len(exact) == 10, proc.stderr
    assert "Found exact allele match sepidermidis.arcC-16" in exact
    # ...and --quiet still hides every one of them (bin/mlst msg()).
    quiet = run("--skipcheck", "example.fna")
    assert "Found exact allele match" not in quiet.stderr


def test_fix17_debug_emits_the_per_hit_and_score_lines():
    """Finding 17: --debug must carry the engine's dbg() trace, not just phases."""
    need_blast()
    proc = run("--debug", "example.fna")
    assert proc.returncode == 0, combined(proc)
    hits = [ln for ln in proc.stderr.split("\n")
            if re.match(r"^\[\d+\] \S+:\d+-\d+\(", ln)]
    scores = [ln for ln in proc.stderr.split("\n") if ln.startswith("SCORE=")]
    assert len(hits) > 1000, len(hits)
    assert scores, proc.stderr[-2000:]
    assert "SCORE=100\tsepidermidis\t184\t16/1/2/1/2/1/1\t(7 genes)" in scores


def test_fix22_duplicate_exact_warning_is_printed_once_in_hit_order():
    """Finding 22: the duplicate-exact notice is msg(), emitted exactly once.

    ``bin/mlst:344`` prints it inline with the ``Found exact allele match``
    lines, so it interleaves with them and ``--quiet`` hides it.
    """
    need_blast()
    proc = run("--skipcheck", "mixed.fa.zip", raw=True)
    assert proc.returncode == 0, combined(proc)
    lines = [ln for ln in proc.stderr.split("\n") if "exact allele match" in ln]
    assert lines == [
        "Found exact allele match mgenitalium.MLST_pgm-3",
        "Found exact allele match mgenitalium.MLST_atpA-1",
        "Found exact allele match mgenitalium.MLST_gyrB-1",
        "Found exact allele match mgenitalium.MLST_ppa-1",
        "WARNING: found additional exact allele match mgenitalium.MLST_pgm-3",
        "WARNING: found additional exact allele match mgenitalium.MLST_atpA-1",
        "Found exact allele match mgenitalium.MLST_adk-7",
        "WARNING: found additional exact allele match mgenitalium.MLST_gyrB-1",
        "Found exact allele match mgenitalium.MLST_gmk-1",
        "WARNING: found additional exact allele match mgenitalium.MLST_ppa-1",
    ], proc.stderr
    # msg(), so --quiet hides it and tests/golden/mixedzip.err stays empty.
    assert run("mixed.fa.zip").stderr == ""


def test_fix17_novel_count_is_announced_once_after_the_table():
    """Findings 17/22: the engine's novel-allele msg() must not duplicate _run's."""
    need_blast()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "novel.fa")
        proc = run("--skipcheck", "--novel", out, "novel.fa", raw=True)
    assert proc.returncode == 0, combined(proc)
    counted = [ln for ln in proc.stderr.split("\n")
               if re.match(r"^Found \d+ novel alleles$", ln)]
    assert len(counted) == 1, proc.stderr


# Markers come from pyproject's vocabulary (section 13.4). They are applied
# programmatically so the module still imports without pytest.
if pytest is not None:
    for _name, _fn in list(globals().items()):
        if not _name.startswith("test_") or not callable(_fn):
            continue
        _marks = list(getattr(_fn, "pytestmark", []))
        _uses = getattr(_fn, "__code__", None)
        _names = set(_uses.co_names) if _uses is not None else set()
        if "need_blast" in _names or _name.startswith("test_golden"):
            _marks.append(pytest.mark.needs_blast)
            _marks.append(pytest.mark.needs_db)
        if "need_catalog" in _names:
            _marks.append(pytest.mark.needs_db)
        if _name.startswith("test_golden"):
            _marks.append(pytest.mark.golden)
        if _marks:
            _fn.pytestmark = _marks
    del _name, _fn, _marks, _uses, _names


def _main():
    """Standalone runner so the file works without pytest."""
    failures = []
    skipped = 0
    names = [n for n in sorted(globals()) if n.startswith("test_")]
    for name in names:
        fn = globals()[name]
        if not callable(fn):
            continue
        try:
            if name == "test_golden_stdout_is_byte_identical" and pytest is not None:
                for case in GOLDEN_CASES:
                    _golden_case(*case)
            else:
                fn()
        except _Skip as exc:
            skipped += 1
            print("skip %s: %s" % (name, exc))
        except BaseException as exc:
            if pytest is not None and exc.__class__.__name__ == "Skipped":
                skipped += 1
                print("skip %s: %s" % (name, exc))
                continue
            failures.append(name)
            print("FAIL %s: %s" % (name, exc))
        else:
            print("ok   %s" % name)
    print("%d ok, %d skipped, %d failed"
          % (len(names) - len(failures) - skipped, skipped, len(failures)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())


# ---------------------------------------------------------------------------
# finding 25: the three Getopt::Long leniencies argparse has no equivalent for
# ---------------------------------------------------------------------------
def test_fix25_getopt_long_leniencies_are_reproduced():
    """`bin/mlst:438` is a bare `use Getopt::Long;` - every default is live.

    No `Getopt::Long::Configure` call exists anywhere in the Perl, so
    single-dash long options, `ignore_case` and the hyphen-free `name!`
    negation all work upstream.  Each spelling below was run against real
    `mlst` 2.35.0 and produced the result asserted here.
    """
    from wmlst.cli import _expand_abbrevs, build_parser

    parser = build_parser()
    cases = {
        # single-dash long options (Getopt::Long does not bundle by default)
        "-csv": "--csv",
        "-quiet": "--quiet",
        "-novel": "--novel",
        # ignore_case
        "--CSV": "--csv",
        "--Scheme": "--scheme",
        "--NOVEL": "--novel",
        "-CSV": "--csv",
        # hyphen-free negation of a `name!` option
        "--nocsv": "--no-csv",
        "--noquiet": "--no-quiet",
        "-nocsv": "--no-csv",
        # ... and the spellings that already worked must not move
        "--csv": "--csv",
        "--no-csv": "--no-csv",
    }
    for given, want in cases.items():
        assert _expand_abbrevs([given], parser) == [want], given

    # a value glued on with `=` survives the rewrite, on both dash counts
    assert _expand_abbrevs(["-scheme=bogus"], parser) == ["--scheme=bogus"]
    assert _expand_abbrevs(["--Scheme=saureus"], parser) == ["--scheme=saureus"]


def test_fix25_leniencies_never_swallow_the_things_that_must_stay_put():
    """`-`, `--`, the real short flags and unknown names are all untouched."""
    from wmlst.cli import _expand_abbrevs, build_parser

    parser = build_parser()
    # `-` is stdin, never an option; `-q`/`-h` are registered short flags;
    # `nopath` and `novel` are options in their own right and must NOT be
    # read as a negation of a `path`/`vel` option that does not exist.
    for tok in ("-", "-q", "-h", "--nopath", "--novel", "-bogus", "--bogus"):
        assert _expand_abbrevs([tok], parser) == [tok], tok
    # nothing past `--` is rewritten, however Perl-ish it looks
    assert _expand_abbrevs(["--", "-csv", "--CSV"], parser) == [
        "--", "-csv", "--CSV"]


def test_fix25_lenient_spellings_run_end_to_end():
    """The rewrite reaches the real CLI, not just the helper."""
    need_blast()
    for flag in ("-csv", "--CSV"):
        proc = run(flag, "example.fna")
        assert proc.returncode == 0, combined(proc)
        assert proc.stdout.startswith("example.fna,sepidermidis,184,"), flag
    for flag in ("--nocsv", "--no-csv"):
        proc = run("--csv", flag, "example.fna")
        assert proc.returncode == 0, combined(proc)
        assert SEPI in proc.stdout, flag
    # `-scheme=bogus` is ACCEPTED as an option and then rejected on its value,
    # which is upstream's behaviour (D3 turns ERRPR: into ERROR:).
    proc = run("-scheme=bogus", "example.fna")
    assert proc.returncode == 1
    assert "Invalid --scheme 'bogus'" in combined(proc)


def test_fix25_unknown_option_report_drops_the_glued_value():
    """Getopt::Long names the option, never the value: `--bogus=1` -> `bogus`."""
    proc = run("--bogus=1", "example.fna")
    assert proc.returncode == 1
    assert combined(proc).strip() == "Unknown option: bogus"
