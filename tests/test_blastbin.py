# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Tests for wmlst/blastbin.py (docs/ARCHITECTURE.md 4.4, 5.5, 8, 9).

Everything that does not need a real BLAST+ runs everywhere; the tests that do
skip themselves when no ``blastn`` can be found (set ``WMLST_BLAST_DIR`` or put
one on ``$PATH``).

Runnable two ways:  ``python3 -m pytest tests/test_blastbin.py``
                    ``python3 tests/test_blastbin.py``
"""

from __future__ import annotations

import hashlib
import io
import os
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from wmlst import any2fasta, blastbin
from wmlst.blastbin import BlastTools
from wmlst.engine import (
    BlastNotFoundError,
    BootstrapError,
    Cancelled,
)

DATA = os.path.join(REPO, "tests", "data")
GOLDEN = os.path.join(REPO, "tests", "golden")
BLASTDB = os.path.join(REPO, "db", "blast", "mlst.fa")

FAKE = BlastTools(
    blastn="/opt/blast/bin/blastn",
    makeblastdb="/opt/blast/bin/makeblastdb",
    blastdbcmd=None,
    version="2.17.0+",
    version_tuple=(2, 17, 0),
    origin="explicit",
)


def _tools():
    """A real BlastTools, or skip."""
    try:
        return blastbin.find_blast(os.environ.get("WMLST_BLASTN"))
    except Exception as exc:  # BlastNotFoundError and friends
        raise unittest.SkipTest("no usable blastn: %s" % exc) from exc


def _need_db():
    if not os.path.isfile(BLASTDB + ".nin"):
        raise unittest.SkipTest("db/blast/mlst.fa index not built")


# ---------------------------------------------------------------------------
# 5.5 — the blastn command line
# ---------------------------------------------------------------------------

def test_blastn_argv_is_an_exact_transliteration():
    argv = blastbin.blastn_argv(
        FAKE, query="/tmp/job/mlst.fna", out="/tmp/job/mlst.bls",
        db_basename="mlst.fa", threads=1, minid=95.0)
    assert argv == [
        "/opt/blast/bin/blastn",
        "-query", "/tmp/job/mlst.fna",
        "-out", "/tmp/job/mlst.bls",
        "-db", "mlst.fa",
        "-num_threads", "1",
        "-ungapped",
        "-dust", "no",
        "-word_size", "32",
        "-max_target_seqs", "100000",
        "-perc_identity", "95",
        "-evalue", "1E-20",
        "-outfmt", "6 sseqid slen length nident qseqid qstart qend qseq sstrand",
    ], argv


def test_blastn_argv_has_no_quotes_and_forbidden_flags():
    argv = blastbin.blastn_argv(
        FAKE, query="/tmp/q", out="/tmp/o", db_basename="mlst.fa",
        threads=8, minid=90.0)
    for item in argv:
        assert "'" not in item and '"' not in item, item
    assert argv.count("-outfmt") == 1
    assert argv[-1] == blastbin.BLAST_OUTFMT
    for forbidden in ("-task", "-mt_mode", "-lcase_masking"):
        assert forbidden not in argv


def test_perc_identity_is_stringified_like_perl():
    def perc(value):
        return blastbin.blastn_argv(
            FAKE, query="/tmp/q", out="/tmp/o", db_basename="d",
            threads=1, minid=value)[-5]

    assert perc(95.0) == "95"
    assert perc(95.5) == "95.5"
    assert perc(100) == "100"
    assert perc(50.0) == "50"
    assert perc(99.99) == "99.99"


def test_makeblastdb_argv():
    argv = blastbin.makeblastdb_argv(FAKE, "/db/blast/mlst.fa")
    assert argv == [
        "/opt/blast/bin/makeblastdb", "-hash_index", "-in", "/db/blast/mlst.fa",
        "-dbtype", "nucl", "-title", "PubMLST", "-parse_seqids",
    ], argv
    # The title is embedded in the index bytes: no vendor string (section 4.2).
    assert "IOWA" not in " ".join(argv)


def test_child_env():
    env = blastbin.child_env("/db/blast")
    assert env["BLASTDB"] == os.path.abspath("/db/blast")
    assert env["BLAST_USAGE_REPORT"] == "false"
    assert env["NCBI_USAGE_REPORT_ENABLED"] == "0"
    assert env["NCBI_DONT_USE_NCBIRC"] == "1"
    assert env["NCBI_DONT_USE_LOCAL_CONFIG"] == "1"
    assert "BLASTDB_LMDB_MAP_SIZE" not in env
    assert "NCBI_CONFIG_OVERRIDES" not in env


def test_child_env_drops_inherited_lmdb_settings():
    os.environ["BLASTDB_LMDB_MAP_SIZE"] = "1000000"
    os.environ["NCBI_CONFIG_OVERRIDES"] = "x"
    try:
        env = blastbin.child_env("/db/blast")
        assert "BLASTDB_LMDB_MAP_SIZE" not in env
        assert "NCBI_CONFIG_OVERRIDES" not in env
    finally:
        os.environ.pop("BLASTDB_LMDB_MAP_SIZE", None)
        os.environ.pop("NCBI_CONFIG_OVERRIDES", None)


# ---------------------------------------------------------------------------
# 4.4 — paths
# ---------------------------------------------------------------------------

def test_ansi_safe_returns_an_absolute_path():
    got = blastbin.ansi_safe("tests/data/example.fna")
    assert got is not None and os.path.isabs(got)
    if not blastbin.IS_WINDOWS:
        assert got == os.path.abspath("tests/data/example.fna")


def test_temp_root_is_writable_and_absolute():
    root = blastbin.temp_root()
    assert os.path.isabs(root) and os.path.isdir(root)
    probe = os.path.join(root, "wmlst-selftest.tmp")
    with open(probe, "w") as fh:
        fh.write("x")
    os.remove(probe)


def test_wmlst_tmpdir_is_honoured():
    with tempfile.TemporaryDirectory() as tmp:
        target = os.path.join(tmp, "scratch")
        os.environ["WMLST_TMPDIR"] = target
        try:
            assert blastbin.temp_root() == os.path.abspath(target)
        finally:
            os.environ.pop("WMLST_TMPDIR", None)


def test_job_dir_is_unique_per_job_and_is_removed():
    seen = []
    with blastbin.job_dir() as a, blastbin.job_dir() as b:
        assert a != b
        assert os.path.isdir(a) and os.path.isdir(b)
        with open(os.path.join(a, "mlst.bls"), "w") as fh:
            fh.write("x")
        seen += [a, b]
    for path in seen:
        assert not os.path.exists(path), path


def test_job_dir_cleanup_never_raises():
    with blastbin.job_dir() as path:
        sub = os.path.join(path, "ro")
        os.mkdir(sub)
        with open(os.path.join(sub, "f"), "w") as fh:
            fh.write("x")
        os.chmod(os.path.join(sub, "f"), stat.S_IREAD)
    assert not os.path.exists(path)


# ---------------------------------------------------------------------------
# 4.4 / 8.2 — run_tool policy
# ---------------------------------------------------------------------------

def test_run_tool_pipes_and_devnull_stdin():
    code = (
        "import sys;"
        "sys.stdout.write('OUT');"
        "sys.stderr.write('ERR');"
        "sys.stdout.write(repr(sys.stdin.read()))"
    )
    proc = blastbin.run_tool([sys.executable, "-c", code], timeout=60)
    assert proc.returncode == 0
    assert proc.stdout == "OUT''"
    assert proc.stderr == "ERR"


def test_run_tool_never_uses_a_shell():
    # An argv list means metacharacters are inert (section 8.2).
    code = "import sys; sys.stdout.write(sys.argv[1])"
    nasty = "a b & c | d ; e $(f) 'g' \"h\""
    proc = blastbin.run_tool([sys.executable, "-c", code, nasty], timeout=60)
    assert proc.stdout == nasty


def test_run_tool_decodes_crlf_and_bad_bytes():
    code = (
        "import sys;"
        "sys.stdout.buffer.write(b'a\\r\\nb\\r\\n\\xff\\xfe')"
    )
    proc = blastbin.run_tool([sys.executable, "-c", code], timeout=60)
    assert proc.stdout.startswith("a\nb\n")
    assert "\ufffd" in proc.stdout


def test_run_tool_passes_env_and_cwd():
    code = "import os,sys; sys.stdout.write(os.environ['BLAST_USAGE_REPORT']+'|'+os.getcwd())"
    with tempfile.TemporaryDirectory() as tmp:
        real = os.path.realpath(tmp)
        proc = blastbin.run_tool(
            [sys.executable, "-c", code], cwd=real, env=blastbin.child_env(real),
            timeout=60)
    assert proc.stdout == "false|" + real


def test_run_tool_timeout_kills_the_child():
    start = time.monotonic()
    try:
        blastbin.run_tool([sys.executable, "-c", "import time; time.sleep(30)"], timeout=1.0)
    except subprocess.TimeoutExpired:
        pass
    else:
        raise AssertionError("no TimeoutExpired")
    assert time.monotonic() - start < 20


def test_run_tool_cancel_raises_cancelled_and_leaves_no_orphan():
    cancel = threading.Event()
    threading.Timer(0.4, cancel.set).start()
    try:
        blastbin.run_tool(
            [sys.executable, "-c", "import time; time.sleep(30)"], cancel=cancel)
    except Cancelled:
        pass
    else:
        raise AssertionError("no Cancelled")
    assert blastbin.terminate_all() == 0


def test_run_tool_missing_executable_is_a_wmlst_error():
    try:
        blastbin.run_tool([os.path.join(REPO, "no-such-binary-xyz")])
    except BlastNotFoundError as exc:
        assert "could not start" in exc.user_message.lower()
    else:
        raise AssertionError("no BlastNotFoundError")


# ---------------------------------------------------------------------------
# 4.4 — checksums and the extraction guard
# ---------------------------------------------------------------------------

def test_sha256sums_roundtrip_and_tamper_detection():
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "a.bin"), "wb") as fh:
            fh.write(b"hello")
        with open(os.path.join(tmp, "b.bin"), "wb") as fh:
            fh.write(b"world")
        blastbin.write_sha256sums(tmp)
        assert blastbin.verify_sha256sums(tmp)
        with open(os.path.join(tmp, "b.bin"), "wb") as fh:
            fh.write(b"tampered")
        assert not blastbin.verify_sha256sums(tmp)
        os.remove(os.path.join(tmp, "b.bin"))
        assert not blastbin.verify_sha256sums(tmp)


def test_verify_sha256sums_is_false_without_the_file():
    with tempfile.TemporaryDirectory() as tmp:
        assert not blastbin.verify_sha256sums(tmp)


def test_safe_member_rejects_everything_dangerous():
    def info(name, **kw):
        ti = tarfile.TarInfo(name)
        ti.type = kw.get("type", tarfile.REGTYPE)
        ti.size = 1
        return ti

    assert blastbin._safe_member(info("ncbi-blast-2.17.0+/bin/blastn.exe"))
    assert blastbin._safe_member(info("ncbi-blast-2.17.0+/LICENSE"))
    assert not blastbin._safe_member(info("ncbi-blast-2.17.0+/bin/blastp.exe"))
    assert not blastbin._safe_member(info("../blastn.exe"))
    assert not blastbin._safe_member(info("a/../../blastn.exe"))
    assert not blastbin._safe_member(info("/etc/blastn.exe"))
    assert not blastbin._safe_member(info("C:/windows/blastn.exe"))
    assert not blastbin._safe_member(info("bin/bl:astn.exe"))
    link = info("ncbi-blast-2.17.0+/bin/blastn.exe", type=tarfile.SYMTYPE)
    assert not blastbin._safe_member(link)
    dev = info("ncbi-blast-2.17.0+/bin/blastn.exe", type=tarfile.CHRTYPE)
    assert not blastbin._safe_member(dev)


def test_bootstrap_constants_are_pinned():
    assert blastbin.BLAST_URL.endswith("ncbi-blast-2.17.0+-x64-win64.tar.gz")
    assert "/2.17.0/" in blastbin.BLAST_URL and "LATEST" not in blastbin.BLAST_URL
    assert blastbin.BLAST_URL.startswith("https://")
    assert blastbin.BLAST_MD5 == "dcd973097407a2910061ff4fb51b09fb"
    assert blastbin.BLAST_BYTES == 143_400_333
    assert "nghttp2.dll" in blastbin.BLAST_MEMBERS
    assert "blastn.exe" in blastbin.BLAST_MEMBERS


def test_install_root_shape():
    root = blastbin.install_root()
    assert os.path.isabs(root)
    if blastbin.IS_WINDOWS:
        assert root.endswith(os.path.join("IOWA-Tech", "WMLST"))
    else:
        assert root.endswith("wmlst")


# ---------------------------------------------------------------------------
# 4.4 — bootstrap, exercised offline against a synthetic mirror
# ---------------------------------------------------------------------------

class _FakeResponse(io.BytesIO):
    def __init__(self, data):
        io.BytesIO.__init__(self, data)
        self.headers = {"Content-Length": str(len(data))}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _patch_mirror(payload, md5_body):
    """Serve `payload` for BLAST_URL and `md5_body` for BLAST_URL.md5."""
    def fake_get(url, timeout=60.0):
        if url.endswith(".md5"):
            return _FakeResponse(md5_body.encode("ascii"))
        return _FakeResponse(payload)
    return fake_get


def test_bootstrap_rejects_a_corrupt_download_and_deletes_the_part_file():
    original_get, original_len = blastbin._http_get, blastbin._content_length
    blastbin._http_get = _patch_mirror(
        b"not the real archive",
        "%s  %s\n" % (blastbin.BLAST_MD5, os.path.basename(blastbin.BLAST_URL)))
    blastbin._content_length = lambda url: None
    try:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                blastbin.bootstrap(tmp)
            except BootstrapError as exc:
                assert blastbin.BLAST_MD5 in str(exc)
            else:
                raise AssertionError("no BootstrapError")
            leftovers = os.listdir(os.path.join(tmp, ".dl"))
            assert leftovers == [], leftovers
    finally:
        blastbin._http_get, blastbin._content_length = original_get, original_len


def test_bootstrap_refuses_an_md5_file_that_names_another_archive():
    original_get = blastbin._http_get
    blastbin._http_get = _patch_mirror(b"x", "%s  other.tar.gz\n" % blastbin.BLAST_MD5)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                blastbin.bootstrap(tmp)
            except BootstrapError as exc:
                assert "other.tar.gz" in str(exc)
            else:
                raise AssertionError("no BootstrapError")
    finally:
        blastbin._http_get = original_get


def test_bootstrap_happy_path_against_a_synthetic_archive():
    if blastbin.IS_WINDOWS:
        raise unittest.SkipTest("the synthetic payload is a POSIX shell blastn")
    with tempfile.TemporaryDirectory() as tmp:
        src_dir = os.path.join(tmp, "src")
        os.mkdir(src_dir)
        stub = os.path.join(src_dir, "blastn")
        with open(stub, "w", newline="\n") as fh:
            fh.write("#!/bin/sh\necho 'blastn: 2.17.0+'\n")
        os.chmod(stub, 0o755)
        licence = os.path.join(src_dir, "LICENSE")
        with open(licence, "w") as fh:
            fh.write("public domain\n")
        evil = os.path.join(src_dir, "evil")
        with open(evil, "w") as fh:
            fh.write("pwned\n")

        archive = os.path.join(src_dir, "payload.tar.gz")
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(stub, arcname="ncbi-blast-2.17.0+/bin/blastn")
            tar.add(stub, arcname="ncbi-blast-2.17.0+/bin/makeblastdb")
            tar.add(licence, arcname="ncbi-blast-2.17.0+/LICENSE")
            tar.add(evil, arcname="../blastn")
            tar.add(evil, arcname="ncbi-blast-2.17.0+/bin/blastp")
        with open(archive, "rb") as fh:
            payload = fh.read()
        digest = hashlib.md5(payload).hexdigest()

        saved = (blastbin._http_get, blastbin._content_length,
                 blastbin.BLAST_MD5, blastbin.BLAST_MEMBERS)
        blastbin._http_get = _patch_mirror(
            payload, "%s  %s\n" % (digest, os.path.basename(blastbin.BLAST_URL)))
        blastbin._content_length = lambda url: len(payload)
        blastbin.BLAST_MD5 = digest
        blastbin.BLAST_MEMBERS = ("blastn", "makeblastdb", "LICENSE")
        events = []
        try:
            root = os.path.join(tmp, "root")
            tools = blastbin.bootstrap(
                root, progress=lambda d, t, text: events.append((d, t, text)))
        finally:
            (blastbin._http_get, blastbin._content_length,
             blastbin.BLAST_MD5, blastbin.BLAST_MEMBERS) = saved

        blast_dir = os.path.join(root, "blast")
        assert tools.version == "2.17.0+" and tools.origin == "appdata"
        assert sorted(os.listdir(blast_dir)) == [
            blastbin.INSTALL_OK, "LICENSE", blastbin.SHA256SUMS, "blastn", "makeblastdb"]
        assert blastbin.verify_sha256sums(blast_dir)
        assert not os.path.exists(os.path.join(root, ".stage"))
        assert not os.path.exists(os.path.join(root, ".dl"))
        assert not os.path.exists(os.path.join(tmp, "blastn"))  # never escaped
        assert events and events[-1][2] == "BLAST+ installed"


# ---------------------------------------------------------------------------
# 4.4 — discovery, with a real BLAST+
# ---------------------------------------------------------------------------

def test_find_blast_and_probe_version():
    tools = _tools()
    assert os.path.isabs(tools.blastn) and os.path.isfile(tools.blastn)
    assert tools.version_tuple >= blastbin.BLAST_MIN_VERSION
    version, vtuple = blastbin.probe_version(tools.blastn)
    assert version == tools.version and vtuple == tools.version_tuple
    assert version.endswith("+")


def test_find_blast_explicit_file_and_directory():
    tools = _tools()
    by_file = blastbin.find_blast(tools.blastn)
    by_dir = blastbin.find_blast(os.path.dirname(tools.blastn))
    assert by_file.blastn == by_dir.blastn == tools.blastn
    assert by_file.origin == "explicit"


def test_find_blast_raises_when_nothing_is_installed():
    saved_path = os.environ.get("PATH")
    saved_dir = os.environ.get("WMLST_BLAST_DIR")
    saved_conda = os.environ.get("CONDA_PREFIX")
    os.environ["PATH"] = os.path.join(REPO, "no-such-dir")
    os.environ.pop("WMLST_BLAST_DIR", None)
    os.environ.pop("CONDA_PREFIX", None)
    try:
        try:
            blastbin.find_blast("/definitely/not/here/blastn")
        except BlastNotFoundError as exc:
            assert "BLAST+" in exc.user_message
        else:
            raise AssertionError("no BlastNotFoundError")
    finally:
        if saved_path is not None:
            os.environ["PATH"] = saved_path
        if saved_dir is not None:
            os.environ["WMLST_BLAST_DIR"] = saved_dir
        if saved_conda is not None:
            os.environ["CONDA_PREFIX"] = saved_conda


def test_world_writable_guard():
    with tempfile.TemporaryDirectory() as tmp:
        os.chmod(tmp, 0o777)
        assert blastbin._is_world_writable(tmp)
        os.chmod(tmp, 0o755)
        assert not blastbin._is_world_writable(tmp)


# ---------------------------------------------------------------------------
# 5.5 — the real thing, against the golden BLAST output
# ---------------------------------------------------------------------------

def test_real_blastn_reproduces_the_golden_output():
    tools = _tools()
    _need_db()
    golden = os.path.join(GOLDEN, "example.blast.tsv")
    if not os.path.isfile(golden):
        raise unittest.SkipTest("no golden BLAST output")
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, "mlst.fna")
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(os.path.join(DATA, "example.fna"), fna)
        proc = blastbin.run_blastn(
            tools, query=fna, out=bls, blastdb=BLASTDB, threads=1, minid=95.0)
        assert proc.returncode == 0, proc.stderr
        with open(bls, "rb") as fh:
            got = fh.read()
    with open(golden, "rb") as fh:
        expected = fh.read()
    assert got.count(b"\n") == 1423
    assert got == expected


def test_blast_exit_3_on_a_sequence_with_no_data():
    tools = _tools()
    _need_db()
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, "mlst.fna")
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(os.path.join(DATA, "empty.fa"), fna)
        proc = blastbin.run_blastn(tools, query=fna, out=bls, blastdb=BLASTDB)
    # Section 8.6: surface the string verbatim; test.sh:62-65 greps for it.
    assert proc.returncode == 3, proc.returncode
    assert "Sequence contains no data" in proc.stderr


def test_bad_database_name_is_exit_2():
    tools = _tools()
    _need_db()
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, "mlst.fna")
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(os.path.join(DATA, "none.fa"), fna)
        proc = blastbin.run_blastn(
            tools, query=fna, out=bls,
            blastdb=os.path.join(os.path.dirname(BLASTDB), "no_such_db.fa"))
    assert proc.returncode == 2, proc.returncode
    assert "BLAST Database error" in proc.stderr


def test_makeblastdb_builds_a_usable_index():
    tools = _tools()
    with blastbin.job_dir() as jd:
        fa = os.path.join(jd, "tiny.fa")
        with open(fa, "w", newline="\n") as fh:
            fh.write(">saureus.arcC_1\n" + "ACGTACGTAC" * 12 + "\n")
        proc = blastbin.run_tool(
            blastbin.makeblastdb_argv(tools, fa),
            cwd=jd, env=blastbin.child_env(jd), timeout=300)
        assert proc.returncode == 0, proc.stderr
        assert os.path.isfile(fa + ".nin")
        query = os.path.join(jd, "q.fna")
        out = os.path.join(jd, "q.bls")
        shutil.copy(fa, query)
        run = blastbin.run_blastn(tools, query=query, out=out, blastdb=fa)
        assert run.returncode == 0, run.stderr
        with open(out, encoding="utf-8") as fh:
            row = fh.readline().rstrip("\r\n").split("\t")
        assert len(row) == 9 and row[0] == "saureus.arcC_1" and row[8] == "plus"


# ---------------------------------------------------------------------------
# plain-python runner
# ---------------------------------------------------------------------------

def _main():
    failures = skipped = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
        except unittest.SkipTest as exc:
            skipped += 1
            print("SKIP %s (%s)" % (name, exc))
        except Exception as exc:
            failures += 1
            import traceback
            print("FAIL %s: %s" % (name, exc))
            traceback.print_exc()
        else:
            print("ok   %s" % name)
    print("\n%d failed, %d skipped" % (failures, skipped))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
