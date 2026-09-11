# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""The Windows runtime matrix (docs/ARCHITECTURE.md 8, 9, 13.4).

Windows 11 is the target, Linux is where this runs in CI, so every case here is
written to execute on both: the genuinely Windows-only APIs are exercised with
``blastbin.IS_WINDOWS`` forced on and the OS calls stubbed, and the rest
(process policy, path shapes, CRLF, cancellation, thread invariance) is real on
either platform.

Runnable two ways:  ``python3 -m pytest tests/test_windows.py``
                    ``python3 tests/test_windows.py``
"""

from __future__ import annotations

import io
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import time
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from wmlst import any2fasta, blastbin
from wmlst.engine import Cancelled

DATA = os.path.join(REPO, "tests", "data")
GOLDEN = os.path.join(REPO, "tests", "golden")
BLASTDB = os.path.join(REPO, "db", "blast", "mlst.fa")
SOURCES = ("any2fasta.py", "blastbin.py")

try:  # pytest marker when running under pytest; harmless otherwise
    import pytest

    pytestmark = pytest.mark.windows
except ImportError:  # pragma: no cover
    pytest = None


def _tools():
    try:
        return blastbin.find_blast(os.environ.get("WMLST_BLASTN"))
    except Exception as exc:
        raise unittest.SkipTest("no usable blastn: %s" % exc) from exc


def _need_db():
    if not os.path.isfile(BLASTDB + ".nin"):
        raise unittest.SkipTest("db/blast/mlst.fa index not built")


def _link_index(dest_dir):
    """Make the 12 .n* index files visible in `dest_dir` as cheaply as possible."""
    os.makedirs(dest_dir, exist_ok=True)
    for name in sorted(os.listdir(os.path.dirname(BLASTDB))):
        if not name.startswith(os.path.basename(BLASTDB) + "."):
            continue  # mlst.fa itself is deliberately not shipped (C10)
        src = os.path.join(os.path.dirname(BLASTDB), name)
        dst = os.path.join(dest_dir, name)
        for attempt in (os.symlink, os.link, shutil.copy):
            try:
                attempt(src, dst)
                break
            except OSError:
                continue
        else:  # pragma: no cover
            raise unittest.SkipTest("cannot stage the BLAST index")
    return os.path.join(dest_dir, os.path.basename(BLASTDB))


def _blast_to_string(tools, query_src, db, threads=1, query_name="mlst.fna"):
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, query_name)
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(query_src, fna)
        proc = blastbin.run_blastn(tools, query=fna, out=bls, blastdb=db,
                                   threads=threads)
        assert proc.returncode == 0, proc.stderr
        with open(bls, "rb") as fh:
            return fh.read()


# ---------------------------------------------------------------------------
# 8.2 — process policy
# ---------------------------------------------------------------------------

class _FakePopen:
    """Captures the kwargs run_tool() would hand to subprocess.Popen."""

    captured: dict = {}   # noqa: RUF012 - a deliberate shared capture slot

    def __init__(self, args, **kwargs):
        _FakePopen.captured = dict(kwargs)
        _FakePopen.captured["args"] = args
        self.returncode = 0

    def communicate(self, timeout=None):
        return "", ""

    def poll(self):
        return 0


def test_windows_launch_flags_hide_every_console():
    """CREATE_NO_WINDOW + STARTUPINFO/SW_HIDE, and nothing else (section 8.2)."""
    saved = (blastbin.IS_WINDOWS, blastbin._startupinfo, subprocess.Popen)
    sentinel = object()
    blastbin.IS_WINDOWS = True
    blastbin._startupinfo = lambda: sentinel
    subprocess.Popen = _FakePopen
    try:
        blastbin.run_tool(["blastn.exe", "-version"])
    finally:
        blastbin.IS_WINDOWS, blastbin._startupinfo, subprocess.Popen = saved
    cap = _FakePopen.captured
    assert cap["creationflags"] == 0x08000000
    assert cap["startupinfo"] is sentinel
    # Never combined with DETACHED_PROCESS (0x8) or CREATE_NEW_CONSOLE (0x10).
    assert cap["creationflags"] & 0x00000008 == 0
    assert cap["creationflags"] & 0x00000010 == 0
    assert cap["stdin"] == subprocess.DEVNULL
    assert cap["stdout"] == subprocess.PIPE and cap["stderr"] == subprocess.PIPE
    assert cap["shell"] is False and isinstance(cap["args"], list)
    assert cap["close_fds"] is True
    assert cap["encoding"] == "utf-8" and cap["errors"] == "replace"


def test_no_creationflags_on_posix():
    saved = subprocess.Popen
    subprocess.Popen = _FakePopen
    try:
        blastbin.run_tool([sys.executable, "-c", "pass"])
    finally:
        subprocess.Popen = saved
    if not blastbin.IS_WINDOWS:
        assert "creationflags" not in _FakePopen.captured


def test_no_shell_strings_anywhere_in_the_runtime():
    for name in SOURCES:
        with open(os.path.join(REPO, "wmlst", name), encoding="utf-8") as fh:
            text = fh.read()
        assert "shell=True" not in text, name
        assert "os.system" not in text, name
        assert "NamedTemporaryFile(" not in text, name


def test_only_blastbin_imports_subprocess():
    with open(os.path.join(REPO, "wmlst", "any2fasta.py"), encoding="utf-8") as fh:
        assert not re.search(r"^import subprocess", fh.read(), re.M)


def test_downloader_is_urllib_not_powershell():
    # urllib attaches no Mark-of-the-Web, so blastn.exe needs no Unblock-File.
    with open(os.path.join(REPO, "wmlst", "blastbin.py"), encoding="utf-8") as fh:
        text = fh.read()
    assert "urllib.request" in text
    for banned in ("powershell", "Invoke-WebRequest", "curl", "certutil"):
        assert banned not in text


def test_run_tool_under_a_windowed_parent_has_no_std_handles():
    """A --windowed build has no valid stdio; DEVNULL is what keeps the child alive."""
    code = "import sys; sys.stdout.write(repr(sys.stdin.read()))"
    proc = blastbin.run_tool([sys.executable, "-c", code], timeout=60)
    assert proc.stdout == "''"


# ---------------------------------------------------------------------------
# 8.3 — bootstrap payload and the VC++ failure mode
# ---------------------------------------------------------------------------

def test_dll_load_failure_names_the_vc_redistributable():
    saved = blastbin.run_tool
    for code in (0xC0000135, 0xC0000142, -1073741515, -1073741502):
        blastbin.run_tool = (
            lambda argv, code=code, **kw: subprocess.CompletedProcess(argv, code, "", "")
        )
        try:
            blastbin.probe_version("blastn.exe")
        except Exception as exc:
            assert "vc_redist" in exc.user_message, exc.user_message
            assert "WinError" not in exc.user_message
        else:
            blastbin.run_tool = saved
            raise AssertionError("no BootstrapError for 0x%08X" % (code & 0xFFFFFFFF,))
        finally:
            blastbin.run_tool = saved


def test_nghttp2_is_in_the_extracted_payload():
    # blastn.exe statically imports it; forgetting it is a hard STATUS_DLL_NOT_FOUND.
    assert "nghttp2.dll" in blastbin.BLAST_MEMBERS
    assert "ncbi-vdb-md.dll" not in blastbin.BLAST_MEMBERS
    assert "blastp.exe" not in blastbin.BLAST_MEMBERS


def test_install_root_is_per_user_and_not_roaming():
    root = blastbin.install_root()
    assert "Roaming" not in root
    assert "OneDrive" not in root


# ---------------------------------------------------------------------------
# 8.4 — text, encoding, newlines
# ---------------------------------------------------------------------------

def test_crlf_blast_output_parses_with_a_clean_sstrand():
    """The single highest-value line of defence in the port (section 8.4)."""
    row = ("saureus.arcC_1\t456\t456\t456\tcontig1\t1\t456\t"
           + "A" * 456 + "\tminus")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "mlst.bls")
        with open(path, "wb") as fh:
            fh.write((row + "\r\n").encode("ascii") * 2)
        with open(path, encoding="utf-8", errors="replace", newline=None) as fh:
            fields = [line.rstrip("\r\n").split("\t") for line in fh]
    assert len(fields) == 2
    for parts in fields:
        assert len(parts) == 9
        assert parts[8] == "minus", repr(parts[8])
        assert "\r" not in parts[8]


def test_converted_query_is_lf_only_even_from_crlf_input():
    raw, _ = None, None
    with tempfile.TemporaryDirectory() as tmp:
        dest = os.path.join(tmp, "mlst.fna")
        any2fasta.convert_to_file(os.path.join(DATA, "issue146.fa"), dest)
        with open(dest, "rb") as fh:
            raw = fh.read()
    assert b"\r" not in raw


def test_crlf_input_and_lf_input_give_the_same_blast_hits():
    """issue146.fa is CRLF; converting it to LF must not move a single hit."""
    tools = _tools()
    _need_db()
    src = os.path.join(DATA, "issue146.fa")
    with tempfile.TemporaryDirectory() as tmp:
        lf_src = os.path.join(tmp, "issue146.lf.fa")
        with open(src, "rb") as fh, open(lf_src, "wb") as out:
            out.write(fh.read().replace(b"\r\n", b"\n"))
        a = _blast_to_string(tools, src, BLASTDB, threads=4)
        b = _blast_to_string(tools, lf_src, BLASTDB, threads=4)
    assert a == b and a.count(b"\n") > 1000


def test_a_crlf_file_of_filenames_reads_clean():
    with tempfile.TemporaryDirectory() as tmp:
        fofn = os.path.join(tmp, "fofn.txt")
        with open(fofn, "wb") as fh:
            fh.write(b"example.fna\r\nexample.fna.gz\r\n\r\n")
        with open(fofn, encoding="utf-8", errors="replace", newline=None) as fh:
            names = [line.rstrip("\r\n") for line in fh if line.strip()]
    assert names == ["example.fna", "example.fna.gz"]


# ---------------------------------------------------------------------------
# 8.5 / 8.7 — awkward paths
# ---------------------------------------------------------------------------

def test_database_directory_containing_a_space():
    """-db is a bare basename with cwd and BLASTDB set (section 5.5)."""
    tools = _tools()
    _need_db()
    golden = os.path.join(GOLDEN, "example.blast.tsv")
    if not os.path.isfile(golden):
        raise unittest.SkipTest("no golden BLAST output")
    with tempfile.TemporaryDirectory() as tmp:
        db = _link_index(os.path.join(tmp, "Program Files", "MLST db"))
        got = _blast_to_string(tools, os.path.join(DATA, "example.fna"), db)
    with open(golden, "rb") as fh:
        assert got == fh.read()


def test_building_a_blast_db_under_a_path_containing_a_space():
    """C:\\Program Files is the DEFAULT Windows install location (sections 4.9, 8.8).

    The sibling test above only LINKS a prebuilt index into a directory with a
    space and queries it; it never runs makeblastdb, which is exactly why this
    was invisible.  makeblastdb splits on whitespace twice -- in the -in value
    (exit 1, "Please provide a database name using -out") and again in the
    absolute path it re-opens the finished database under (exit 2, "No alias or
    index file found ... [C:\\Program]", .njs never written) -- so updatedb goes
    through blastbin.run_makeblastdb(), which builds in a whitespace-free
    scratch directory and moves the relocatable index into place.
    """
    from wmlst import updatedb

    tools = _tools()
    with tempfile.TemporaryDirectory() as tmp:
        dbdir = os.path.join(tmp, "Program Files", "MLST db")
        scheme = os.path.join(dbdir, "pubmlst", "saureus")
        os.makedirs(scheme)
        assert " " in dbdir
        with open(os.path.join(scheme, "arcC.tfa"), "w", newline="\n") as fh:
            fh.write(">arcC_1\n" + "ACGTACGTAC" * 12 + "\n")
        with open(os.path.join(scheme, "aroE.tfa"), "w", newline="\n") as fh:
            fh.write(">aroE_1\n" + "TTGACCAGTC" * 12 + "\n")
        with open(os.path.join(scheme, "saureus.txt"), "w", newline="\n") as fh:
            fh.write("ST\tarcC\taroE\n1\t1\t1\n")

        fasta = updatedb.build_blast_db(dbdir, tools)

        assert os.path.isfile(fasta)
        missing = [ext for ext in updatedb.BLAST_INDEX_EXTENSIONS
                   if not os.path.isfile(fasta + "." + ext)]
        assert not missing, missing
        # The staging directory must be gone, not left in the destination.
        assert not os.path.isdir(
            os.path.join(dbdir, updatedb.BLAST_DIRNAME, updatedb.STAGING_DIRNAME))
        # And the relocated index must answer a real query.
        query = os.path.join(tmp, "q.fna")
        shutil.copy(fasta, query)
        with blastbin.job_dir() as jd:
            out = os.path.join(jd, "q.bls")
            run = blastbin.run_blastn(tools, query=query, out=out, blastdb=fasta)
            assert run.returncode == 0, run.stderr
            with open(out, encoding="utf-8") as fh:
                hits = [line.split("\t")[0] for line in fh]
    assert "saureus.arcC_1" in hits and "saureus.aroE_1" in hits


def test_path_rung_survives_the_st_mode_windows_synthesises():
    """Windows reports 0o40777 for every writable directory (section 4.4).

    CPython's attributes_to_mode() makes st_mode a function of the file
    ATTRIBUTES, not of any ACL: a plain directory is 0o40777 and a read-only one
    0o40555.  Testing S_IWOTH there inverted the binary-planting guard and made
    the PATH rung -- the only rung a `pip install wmlst` user with NCBI's own
    BLAST+ can reach -- reject every ordinary install directory.
    """
    which = shutil.which("blastn")
    if not which:
        raise unittest.SkipTest("no blastn on PATH")
    bindir = os.path.dirname(os.path.realpath(which))
    target = os.path.abspath(bindir)

    class _VersionPopen:
        def __init__(self, args, **kwargs):
            self.returncode = 0

        def communicate(self, timeout=None):
            return "blastn: 2.17.0+\n", ""

        def poll(self):
            return 0

    real_stat = blastbin.os.stat

    def windows_stat(path, *a, **kw):
        st = real_stat(path, *a, **kw)
        if os.path.abspath(path) == target:
            return os.stat_result([stat.S_IFDIR | 0o111 | 0o666, *tuple(st)[1:]])
        return st

    saved = (blastbin.IS_WINDOWS, blastbin._startupinfo, subprocess.Popen,
             blastbin.os.stat)
    saved_env = {k: os.environ.get(k)
                 for k in ("PATH", "WMLST_BLAST_DIR", "CONDA_PREFIX")}
    try:
        blastbin.IS_WINDOWS = True
        blastbin._startupinfo = lambda: None
        subprocess.Popen = _VersionPopen
        blastbin.os.stat = windows_stat
        os.environ.pop("WMLST_BLAST_DIR", None)
        os.environ.pop("CONDA_PREFIX", None)
        os.environ["PATH"] = bindir
        assert windows_stat(bindir).st_mode & stat.S_IWOTH  # the trap
        assert not blastbin._is_world_writable(bindir)      # ... now disarmed
        tools = blastbin.find_blast()
    finally:
        (blastbin.IS_WINDOWS, blastbin._startupinfo, subprocess.Popen,
         blastbin.os.stat) = saved
        for key, value in saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    assert tools.origin == "path", tools.origin


def test_query_filename_with_shell_metacharacters():
    tools = _tools()
    _need_db()
    with tempfile.TemporaryDirectory() as tmp:
        nasty = os.path.join(tmp, "we ird '&^%!()+ name.fna")
        shutil.copy(os.path.join(DATA, "none.fa"), nasty)
        with blastbin.job_dir() as jd:
            fna = os.path.join(jd, "mlst.fna")
            bls = os.path.join(jd, "mlst.bls")
            any2fasta.convert_to_file(nasty, fna)
            proc = blastbin.run_blastn(tools, query=fna, out=bls, blastdb=BLASTDB)
        assert proc.returncode == 0, proc.stderr


def test_non_ascii_input_path():
    tools = _tools()
    _need_db()
    with tempfile.TemporaryDirectory() as tmp:
        try:
            folder = os.path.join(tmp, "Ren\u00e9 \u00dcmlaut \u4e2d\u6587")
            os.mkdir(folder)
        except (UnicodeEncodeError, OSError) as exc:  # pragma: no cover
            raise unittest.SkipTest(
                "filesystem cannot hold non-ASCII names") from exc
        src = os.path.join(folder, "pr\u00f8ve.fna")
        shutil.copy(os.path.join(DATA, "none.fa"), src)
        # The job dir lives under the ASCII temp root, so blastn never sees this path.
        with blastbin.job_dir() as jd:
            fna = os.path.join(jd, "mlst.fna")
            bls = os.path.join(jd, "mlst.bls")
            any2fasta.convert_to_file(src, fna)
            proc = blastbin.run_blastn(tools, query=fna, out=bls, blastdb=BLASTDB)
        assert proc.returncode == 0, proc.stderr


def test_input_from_a_path_longer_than_260_characters():
    with tempfile.TemporaryDirectory() as tmp:
        deep = tmp
        while len(deep) < 300:
            deep = os.path.join(deep, "a_very_long_directory_component")
        try:
            os.makedirs(deep, exist_ok=True)
        except OSError as exc:  # pragma: no cover
            raise unittest.SkipTest(
                "cannot create a 300-char path: %s" % exc) from exc
        src = os.path.join(deep, "sample.fna")
        shutil.copy(os.path.join(DATA, "empty.fa"), src)
        assert len(src) > 260
        out = io.StringIO()
        assert any2fasta.convert(src, out) == (1, 0, "FASTA")


def test_ansi_safe_ladder():
    if blastbin.IS_WINDOWS:  # pragma: no cover - exercised on the target platform
        safe = blastbin.ansi_safe(os.environ.get("LOCALAPPDATA", "C:\\"))
        assert safe is None or len(safe) < 200
    else:
        assert blastbin.ansi_safe("/tmp/x") == "/tmp/x"
    # temp_root() must stay short enough for <root>\wmlst-<pid>-<n>\mlst.bls.
    assert len(blastbin.temp_root()) <= 60


def test_temp_root_is_ansi_safe():
    assert blastbin.ansi_safe(blastbin.temp_root()) is not None


# ---------------------------------------------------------------------------
# 9 — concurrency, cancellation, stale files
# ---------------------------------------------------------------------------

def test_thread_count_does_not_change_the_output_bytes():
    tools = _tools()
    _need_db()
    one = _blast_to_string(tools, os.path.join(DATA, "example.fna"), BLASTDB, threads=1)
    four = _blast_to_string(tools, os.path.join(DATA, "example.fna"), BLASTDB, threads=4)
    assert one == four and one.count(b"\n") == 1423


def test_every_job_gets_a_private_empty_directory():
    """Sharing mlst.fna/mlst.bls is the classic port bug (section 9)."""
    seen = set()
    for _ in range(4):
        with blastbin.job_dir() as jd:
            assert os.listdir(jd) == []
            assert jd not in seen
            seen.add(jd)
            with open(os.path.join(jd, "mlst.bls"), "w") as fh:
                fh.write("stale")
    for jd in seen:
        assert not os.path.exists(jd)


def test_cancelling_a_real_blastn_leaves_no_orphan():
    tools = _tools()
    _need_db()
    cancel = threading.Event()
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, "mlst.fna")
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(os.path.join(DATA, "issue146.fa"), fna)
        threading.Timer(0.6, cancel.set).start()
        start = time.monotonic()
        try:
            blastbin.run_blastn(tools, query=fna, out=bls, blastdb=BLASTDB,
                                threads=1, cancel=cancel)
        except Cancelled:
            pass
        else:
            raise unittest.SkipTest("blastn finished before the cancel fired")
        assert time.monotonic() - start < 20
    # The registry is empty and the partial .bls went with the job dir.
    assert blastbin.terminate_all() == 0
    assert not os.path.exists(jd)


def test_child_environment_disables_the_usage_report():
    """No Defender Firewall prompt over the Tkinter window (section 4.4)."""
    code = (
        "import os,sys;"
        "sys.stdout.write('|'.join("
        "os.environ.get(k, 'MISSING') for k in "
        "('BLAST_USAGE_REPORT','NCBI_USAGE_REPORT_ENABLED',"
        "'NCBI_DONT_USE_NCBIRC','NCBI_DONT_USE_LOCAL_CONFIG','BLASTDB')))"
    )
    db_dir = os.path.dirname(BLASTDB)
    proc = blastbin.run_tool([sys.executable, "-c", code],
                            env=blastbin.child_env(db_dir), timeout=60)
    assert proc.stdout == "false|0|1|1|" + os.path.abspath(db_dir)


# ---------------------------------------------------------------------------
# 8.6 — exit codes
# ---------------------------------------------------------------------------

def test_empty_sequence_is_exit_3_with_the_verbatim_message():
    tools = _tools()
    _need_db()
    with blastbin.job_dir() as jd:
        fna = os.path.join(jd, "mlst.fna")
        bls = os.path.join(jd, "mlst.bls")
        any2fasta.convert_to_file(os.path.join(DATA, "empty.fa"), fna)
        proc = blastbin.run_blastn(tools, query=fna, out=bls, blastdb=BLASTDB)
    assert proc.returncode == 3
    assert "BLAST engine error: Warning: Sequence contains no data" in proc.stderr
    with open(os.path.join(GOLDEN, "emptyfa.err"), encoding="utf-8") as fh:
        # the golden line keeps its trailing space; so must ours
        assert fh.readline().rstrip("\r\n") == proc.stderr.splitlines()[0]


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
