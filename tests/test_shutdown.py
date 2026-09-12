# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Shutdown, cancellation and "no orphaned blastn" (docs/ARCHITECTURE.md 9).

The bug these tests exist for: closing the window left the process alive, so
only Task Manager could end it.  The mechanism, reproduced verbatim by
:func:`test_a_parked_pool_worker_hangs_the_interpreter_at_exit`:

* ``Engine.analyse`` runs ``--jobs`` files on a ``ThreadPoolExecutor``;
* since Python 3.9 its worker threads are NOT daemon threads, and
  ``concurrent.futures.thread`` registers an ``atexit`` hook that JOINS every
  one of them at interpreter shutdown;
* a worker blocked on a running ``blastn`` is therefore never released -- the
  window closes, the join blocks forever and the process survives.

The cure is that a shutdown reaches the CHILDREN: ``blastbin.shutdown()`` kills
every tracked child and refuses new launches, which releases the workers, and
``Engine.shutdown()`` / ``engine.shutdown_all()`` cancel the run and drain the
pool under a deadline, reporting whether everything stopped so a GUI can decide
about a last-resort force exit.

Run with ``python -m pytest tests/test_shutdown.py``.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import CancelledError

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
DBDIR = os.path.join(REPO_ROOT, "db")

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from wmlst import blastbin
from wmlst import engine as engine_mod
from wmlst.engine import Cancelled, Engine, RunConfig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
@pytest.fixture(autouse=True)
def _clean_shutdown_state():
    """Never leave the module latched: the rest of the suite launches children."""
    blastbin.reset_shutdown()
    yield
    blastbin.reset_shutdown()
    blastbin.terminate_all(timeout=2.0)


def _sleeper(seconds: float = 60.0):
    """A child that outlives any test unless something stops it."""
    return [sys.executable, "-c", "import time; time.sleep(%r)" % (seconds,)]


def _wait_for(predicate, timeout: float = 10.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    return predicate()


def _tracked() -> int:
    with blastbin._PROCS_LOCK:
        return len(blastbin._PROCS)


def _tracked_pids():
    with blastbin._PROCS_LOCK:
        return sorted(proc.pid for proc in blastbin._PROCS)


def _search_pids():
    """PIDs of the tracked children that are real ``blastn`` SEARCHES.

    Not every child is one: resolving ``tools`` runs ``blastn -version``
    first, and stopping during that probe would prove nothing about a worker
    parked on a search.
    """
    with blastbin._PROCS_LOCK:
        procs = list(blastbin._PROCS)
    return sorted(proc.pid for proc in procs
                  if "-query" in [str(a) for a in (proc.args or ())])


def _child_names(pid: int):
    """Names of the live children of `pid`. POSIX only; () elsewhere."""
    if not os.path.isdir("/proc"):
        return ()
    try:
        out = subprocess.run(["ps", "-o", "comm=", "--ppid", str(pid)],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             timeout=10)
    except (OSError, subprocess.SubprocessError):  # pragma: no cover
        return ()
    return tuple(line.strip() for line in out.stdout.decode().splitlines()
                 if line.strip())


def _pid_is_gone(pid: int) -> bool:
    """True when `pid` has exited.

    NOT ``os.kill(pid, 0)`` on Windows. That is not a liveness probe there:
    CPython maps os.kill to OpenProcess+TerminateProcess for every signal
    except CTRL_C_EVENT/CTRL_BREAK_EVENT, so the "probe" would kill whatever
    it touched. It is also wrong in the other direction -- a process that has
    already exited keeps an openable PID for as long as anyone holds a handle
    to it, and Popen holds one until it is reaped, so a dead child reported as
    alive. Ask the kernel whether the process object is signalled instead.
    """
    if os.path.isdir("/proc"):
        return not os.path.isdir("/proc/%d" % pid)
    if os.name == "nt":  # pragma: no cover - Windows only
        import ctypes

        SYNCHRONIZE = 0x00100000
        PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
        WAIT_OBJECT_0 = 0
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(
            SYNCHRONIZE | PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not handle:
            return True  # no such process, or it is already reaped
        try:
            # A process object becomes signalled the moment the process exits,
            # regardless of who still holds a handle to it.
            return kernel32.WaitForSingleObject(handle, 0) == WAIT_OBJECT_0
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except OSError:
        return True
    return False


def _run_config(files, **kw):
    return RunConfig(files=tuple(files), dbdir=DBDIR, **kw)


def _blastn_copies(tmp_path, n):
    """`n` distinct copies of a genome big enough to keep blastn busy."""
    import shutil

    src = os.path.join(DATA, "issue146.fa")
    out = []
    for i in range(n):
        dst = str(tmp_path / ("sample%d.fa" % i))
        shutil.copyfile(src, dst)
        out.append(dst)
    return out


# ---------------------------------------------------------------------------
# 1. The mechanism, reproduced -- this is what the user hit
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_a_parked_pool_worker_hangs_the_interpreter_at_exit():
    """Control: a non-daemon pool worker blocks exit even from a daemon thread.

    No WMLST code is involved -- this is the CPython behaviour the fix works
    around, asserted so nobody 'simplifies' the fix away later.
    """
    script = (
        "import concurrent.futures as cf, threading, time\n"
        "def worker():\n"
        "    pool = cf.ThreadPoolExecutor(max_workers=2)\n"
        "    pool.submit(time.sleep, 60).result()\n"
        "threading.Thread(target=worker, daemon=True).start()\n"
        "time.sleep(0.5)\n"
        "print('main returning', flush=True)\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", script],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        with pytest.raises(subprocess.TimeoutExpired):
            proc.communicate(timeout=5)
    finally:
        proc.kill()
        proc.wait()


# ---------------------------------------------------------------------------
# 2. blastbin: shutdown reaches the children
# ---------------------------------------------------------------------------
def test_shutdown_refuses_new_launches():
    blastbin.shutdown(1.0)
    assert blastbin.is_shutting_down() is True
    with pytest.raises(Cancelled):
        blastbin.run_tool(_sleeper(5))
    assert _tracked() == 0


def test_reset_shutdown_allows_launching_again():
    blastbin.shutdown(1.0)
    blastbin.reset_shutdown()
    assert blastbin.is_shutting_down() is False
    proc = blastbin.run_tool([sys.executable, "-c", "print('alive')"], timeout=30)
    assert proc.returncode == 0
    assert proc.stdout.strip() == "alive"


def test_shutdown_kills_a_running_child_and_releases_its_worker():
    """The heart of the fix: the parked worker returns, and fast."""
    raised = []

    def body():
        try:
            blastbin.run_tool(_sleeper(120))
        except BaseException as exc:   # recorded, then asserted
            raised.append(exc)

    worker = threading.Thread(target=body, name="parked-worker", daemon=True)
    worker.start()
    assert _wait_for(lambda: _tracked() == 1, 10.0), "the child was never tracked"
    pids = _tracked_pids()

    start = time.monotonic()
    clean = blastbin.shutdown(5.0)
    elapsed = time.monotonic() - start

    assert clean is True
    assert elapsed < 5.0
    worker.join(5.0)
    assert not worker.is_alive(), "the worker was not released"
    assert raised and isinstance(raised[0], Cancelled)
    assert _tracked() == 0
    assert blastbin.terminate_all() == 0
    for pid in pids:
        assert _wait_for(lambda p=pid: _pid_is_gone(p), 5.0), (
            "child %d survived shutdown" % pid)


def test_shutdown_is_clean_and_quick_with_nothing_running():
    start = time.monotonic()
    assert blastbin.shutdown(5.0) is True
    assert time.monotonic() - start < 1.0


def test_shutdown_never_waits_longer_than_its_timeout():
    """Even with several children, the grace is shared, not multiplied."""
    raised = []

    def body():
        try:
            blastbin.run_tool(_sleeper(120))
        except BaseException as exc:   # asserted after the join
            raised.append(type(exc))

    threads = []
    for i in range(4):
        t = threading.Thread(target=body, name="parked-%d" % i, daemon=True)
        t.start()
        threads.append(t)
    assert _wait_for(lambda: _tracked() == 4, 10.0)
    start = time.monotonic()
    clean = blastbin.shutdown(6.0)
    elapsed = time.monotonic() - start
    assert clean is True
    assert elapsed < 6.0, "shutdown took %.1f s" % elapsed
    for t in threads:
        t.join(5.0)
        assert not t.is_alive()
    assert raised == [Cancelled] * 4


# ---------------------------------------------------------------------------
# 3. run_tool's wait loop: bounded, cancel-observant, unchanged otherwise
# ---------------------------------------------------------------------------
def test_cancel_is_honoured_promptly_with_no_timeout_set():
    """timeout=None used to mean an UNBOUNDED wait when cancel was None."""
    cancel = threading.Event()
    threading.Timer(0.3, cancel.set).start()
    start = time.monotonic()
    with pytest.raises(Cancelled):
        blastbin.run_tool(_sleeper(120), timeout=None, cancel=cancel)
    assert time.monotonic() - start < 5.0
    assert _tracked() == 0


def test_cancel_set_before_the_call_stops_at_once():
    cancel = threading.Event()
    cancel.set()
    start = time.monotonic()
    with pytest.raises(Cancelled):
        blastbin.run_tool(_sleeper(120), cancel=cancel)
    assert time.monotonic() - start < 5.0
    assert _tracked() == 0


def test_the_bounded_poll_still_collects_output_correctly():
    """The polling loop must not truncate or mangle a normal run."""
    code = ("import sys, time; time.sleep(0.5); "
            "sys.stdout.write('out\\n'); sys.stderr.write('err\\n')")
    proc = blastbin.run_tool([sys.executable, "-c", code], timeout=None)
    assert proc.returncode == 0
    assert proc.stdout == "out\n"
    assert proc.stderr == "err\n"
    assert _tracked() == 0


def test_timeout_still_raises_timeout_expired_and_reaps():
    start = time.monotonic()
    with pytest.raises(subprocess.TimeoutExpired):
        blastbin.run_tool(_sleeper(120), timeout=0.5)
    assert time.monotonic() - start < 5.0
    assert _tracked() == 0


def test_terminate_all_counts_what_it_tracked():
    assert blastbin.terminate_all() == 0
    started = threading.Event()

    def body():
        try:
            started.set()
            blastbin.run_tool(_sleeper(120))
        except BaseException:          # the terminate is the point
            pass

    t = threading.Thread(target=body, daemon=True)
    t.start()
    assert _wait_for(lambda: _tracked() == 1, 10.0)
    assert blastbin.terminate_all(timeout=5.0) == 1
    t.join(5.0)
    assert not t.is_alive()


# ---------------------------------------------------------------------------
# 4. engine: cancellation view, pool naming, shutdown contract
# ---------------------------------------------------------------------------
def test_engine_shutdown_before_a_run_makes_analyse_raise_cancelled():
    eng = Engine(_run_config([os.path.join(DATA, "example.fna")]))
    assert eng.shutdown(2.0) is True
    with pytest.raises(Cancelled):
        eng.analyse()


def test_engine_shutdown_is_clean_when_idle():
    eng = Engine(_run_config([os.path.join(DATA, "example.fna")]))
    start = time.monotonic()
    assert eng.shutdown(5.0) is True
    assert time.monotonic() - start < 2.0


def test_shutdown_all_is_clean_when_nothing_runs():
    Engine(_run_config([os.path.join(DATA, "example.fna")]))
    start = time.monotonic()
    assert engine_mod.shutdown_all(5.0) is True
    assert time.monotonic() - start < 2.0


def test_an_engine_side_abort_never_sets_the_callers_cancel_event():
    """The GUI reads its own Event to mean 'the user pressed Cancel'."""
    eng = Engine(_run_config([os.path.join(DATA, "example.fna")]))
    caller = threading.Event()
    view = eng._cancel_view(caller)
    view.set()
    assert view.is_set() is True
    assert caller.is_set() is False
    caller.set()
    assert eng._cancel_view(caller).is_set() is True


def test_the_cancel_view_is_not_wrapped_twice():
    eng = Engine(_run_config([os.path.join(DATA, "example.fna")]))
    view = eng._cancel_view(None)
    assert eng._cancel_view(view) is view


@pytest.mark.needs_db
def test_job_threads_carry_the_wmlst_name_prefix():
    """A named pool is what makes a stuck worker identifiable in a dump."""
    names = []

    class Spy(Engine):
        def analyse_file(self, path, **kw):
            names.append(threading.current_thread().name)
            return super().analyse_file(path, **kw)

    missing = [os.path.join(DATA, "no-such-file-%d.fa" % i) for i in range(4)]
    eng = Spy(_run_config(missing, jobs=4))
    result = eng.analyse()
    assert len(result.samples) == 4
    assert all(sample.failed for sample in result.samples)
    assert names and all(n.startswith("wmlst-job") for n in names), names


@pytest.mark.needs_db
def test_results_stay_in_argv_order_through_the_pool():
    missing = [os.path.join(DATA, "no-such-file-%d.fa" % i) for i in range(6)]
    eng = Engine(_run_config(missing, jobs=4))
    result = eng.analyse()
    assert [s.path for s in result.samples] == missing


@pytest.mark.needs_db
def test_shutdown_mid_pool_raises_wmlst_cancelled_not_futures_cancellederror():
    """Callers catch ``engine.Cancelled``; a dropped future must not leak out.

    ``shutdown`` drains the pool with ``cancel_futures=True``, so the queued
    jobs die as ``concurrent.futures.CancelledError`` -- a different class,
    which every caller (CLI and GUI alike) would report as a crash.
    """

    class Slow(Engine):
        def analyse_file(self, path, **kw):
            time.sleep(0.4)
            return super().analyse_file(path, **kw)

    missing = [os.path.join(DATA, "no-such-file-%d.fa" % i) for i in range(12)]
    eng = Slow(_run_config(missing, jobs=2))
    outcome = {}

    def body():
        try:
            outcome["result"] = eng.analyse()
        except BaseException as exc:   # asserted below
            outcome["exc"] = exc

    runner = threading.Thread(target=body, daemon=True)
    runner.start()
    time.sleep(0.5)
    assert eng.shutdown(5.0) is True
    runner.join(10.0)
    assert not runner.is_alive()
    assert isinstance(outcome.get("exc"), Cancelled), outcome


@pytest.mark.needs_db
def test_a_dropped_future_surfaces_as_wmlst_cancelled():
    """The same conversion, forced rather than raced.

    ``_drain_pool`` uses ``cancel_futures=True``; a job that never ran comes
    back as ``concurrent.futures.CancelledError``, and ``_map_files`` must
    translate it, because the CLI and the GUI both catch ``engine.Cancelled``
    and would otherwise report a crash.
    """
    eng = Engine(_run_config([os.path.join(DATA, "example.fna")] * 2, jobs=2))

    def one(i_path):
        raise CancelledError("the pool dropped this job")

    with pytest.raises(Cancelled):
        eng._map_files(["a", "b"], one, eng._cancel_view(None))


# ---------------------------------------------------------------------------
# 5. The real thing: a --jobs 4 run against the bundled DB, stopped mid-flight
# ---------------------------------------------------------------------------
@pytest.mark.needs_blast
@pytest.mark.needs_db
@pytest.mark.slow
def test_shutdown_stops_a_jobs_4_run_and_leaves_no_orphan_blastn(tmp_path):
    files = _blastn_copies(tmp_path, 4)
    eng = Engine(_run_config(files, jobs=4))
    outcome = {}

    def body():
        try:
            outcome["result"] = eng.analyse()
        except BaseException as exc:   # asserted below
            outcome["exc"] = exc

    runner = threading.Thread(target=body, name="wmlst-analysis-test",
                              daemon=True)
    runner.start()

    # Wait for TWO concurrent searches: that is the pool actually parked on
    # children, which is the state the bug needs.
    if not _wait_for(lambda: len(_search_pids()) >= 2, 180.0):
        pytest.skip("blastn searches never overlapped; nothing to shut down")
    pids = _tracked_pids()
    assert len(_search_pids()) >= 2

    start = time.monotonic()
    clean = eng.shutdown(10.0)
    elapsed = time.monotonic() - start

    assert clean is True, "shutdown reported %r after %.1f s" % (clean, elapsed)
    assert elapsed < 10.0, "shutdown took %.1f s" % elapsed

    runner.join(10.0)
    assert not runner.is_alive(), "analyse() never returned"
    assert isinstance(outcome.get("exc"), Cancelled), outcome

    assert _tracked() == 0
    assert blastbin.terminate_all() == 0
    for pid in pids:
        assert _wait_for(lambda p=pid: _pid_is_gone(p), 10.0), (
            "orphaned blastn %d" % pid)
    assert "blastn" not in "".join(_child_names(os.getpid()))

    # Every pool worker really finished: nothing is left for atexit to join.
    assert _wait_for(
        lambda: not [t for t in threading.enumerate()
                     if t.name.startswith("wmlst-job") and t.is_alive()],
        10.0), [t.name for t in threading.enumerate()]


@pytest.mark.needs_blast
@pytest.mark.needs_db
@pytest.mark.slow
def test_the_interpreter_exits_after_shutdown_all(tmp_path):
    """End to end, in a real child process: the GUI shape, then exit.

    A daemon thread runs a ``--jobs 4`` analysis; the main thread asks for a
    shutdown and returns, exactly like ``WM_DELETE_WINDOW`` -> ``destroy()`` ->
    ``main()`` returning.  Before the fix this child never exited.
    """
    files = _blastn_copies(tmp_path, 4)
    pidfile = str(tmp_path / "child-pids.txt")
    script = """
import json, os, sys, threading, time
sys.path.insert(0, %(repo)r)
from wmlst import blastbin
from wmlst import engine as engine_mod

cfg = engine_mod.RunConfig(files=tuple(%(files)r), dbdir=%(dbdir)r, jobs=4)
eng = engine_mod.Engine(cfg)

def body():
    try:
        eng.analyse()
    except BaseException:
        pass

threading.Thread(target=body, daemon=True).start()
deadline = time.time() + 180
pids = []
while time.time() < deadline:
    with blastbin._PROCS_LOCK:
        procs = list(blastbin._PROCS)
    pids = sorted(p.pid for p in procs
                  if "-query" in [str(a) for a in (p.args or ())])
    if len(pids) >= 2:
        break
    time.sleep(0.05)
open(%(pidfile)r, "w").write(json.dumps(pids))
clean = engine_mod.shutdown_all(10.0)
print("clean=%%s" %% clean, flush=True)
""" % {"repo": REPO_ROOT, "files": files, "dbdir": DBDIR, "pidfile": pidfile}

    env = dict(os.environ)
    env.setdefault("BLAST_USAGE_REPORT", "false")
    start = time.monotonic()
    proc = subprocess.run([sys.executable, "-c", script], timeout=120,
                          capture_output=True, env=env)
    elapsed = time.monotonic() - start

    assert proc.returncode == 0, proc.stderr.decode()
    assert b"clean=True" in proc.stdout, proc.stdout + proc.stderr
    assert elapsed < 120

    import json

    with open(pidfile) as fh:
        pids = json.load(fh)
    assert len(pids) >= 2, "the child never ran two blastn searches"
    for pid in pids:
        assert _wait_for(lambda p=pid: _pid_is_gone(p), 10.0), (
            "orphaned blastn %d survived its parent" % pid)


@pytest.mark.needs_blast
@pytest.mark.needs_db
@pytest.mark.slow
@pytest.mark.skipif(os.name == "nt", reason="POSIX signal delivery")
def test_ctrl_c_stops_a_jobs_4_cli_run(tmp_path):
    """The same hang bit the CLI: Ctrl-C with --jobs > 1 could not exit.

    ``KeyboardInterrupt`` lands in the main thread, but the pool's ``__exit__``
    joined the workers, which were parked on ``blastn``. Now the abort path
    releases them, so the process really does end.
    """
    files = _blastn_copies(tmp_path, 4)
    env = dict(os.environ)
    env.setdefault("BLAST_USAGE_REPORT", "false")
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.Popen(
        [sys.executable, "-m", "wmlst", "--jobs", "4", *files],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    try:
        if not _wait_for(
                lambda: "blastn" in " ".join(_child_names(proc.pid)), 180.0):
            pytest.skip("no blastn child to interrupt")
        proc.send_signal(signal.SIGINT)
        start = time.monotonic()
        proc.communicate(timeout=120)
        # 5 s is deliberately below the ~5.6 s a single blastn needs for
        # issue146.fa: waiting the children out instead of stopping them
        # cannot pass this.
        assert time.monotonic() - start < 5.0
        assert "blastn" not in " ".join(_child_names(proc.pid))
    finally:
        if proc.poll() is None:  # pragma: no cover - only if the fix regressed
            proc.kill()
            proc.communicate()


@pytest.mark.needs_blast
@pytest.mark.needs_db
@pytest.mark.slow
def test_jobs_4_is_byte_identical_to_jobs_1(tmp_path):
    """The fix must not disturb the parallel path's output (section 9)."""
    files = ["example.fna", "issue146.fa", "messy.fa", "novel.fa"]
    base = [sys.executable, "-m", "wmlst", "--quiet", *files]
    env = dict(os.environ)
    env.setdefault("BLAST_USAGE_REPORT", "false")
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")

    def run(jobs):
        proc = subprocess.run([*base, "--jobs", str(jobs)], cwd=DATA,
                              capture_output=True, timeout=900, env=env)
        assert proc.returncode == 0, proc.stderr.decode()
        return proc.stdout

    assert run(4) == run(1)
