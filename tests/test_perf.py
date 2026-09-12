# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Tests for :mod:`wmlst.perf`, the GUI's start-up auto-tuner.

Three things are pinned here:

1. The sizing rule the user asked for -- one file per four cores, four BLAST
   threads each -- at every boundary that matters (1, 2, 3, 4, 7, 8, 12, 16,
   32, 64 cores), plus the invariants that hold at *every* core count.
2. The memory clamp, which may only ever reduce ``jobs``.
3. That this module changes no CLI default. ``RunConfig`` stays 1/1
   (docs/ARCHITECTURE.md C3) and nothing outside the GUI may import perf.

No test here touches BLAST, the database or a display.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path

import pytest

from wmlst import perf

REPO_ROOT = Path(__file__).resolve().parents[1]

MB = 1024 * 1024
GB = 1024 * MB


# ---------------------------------------------------------------------------
# The rule: one file per four cores, four threads each
# ---------------------------------------------------------------------------
#: (cores, expected threads, expected jobs). The three the user wrote out by
#: hand -- 4 -> 1 file, 8 -> 2 files, 16 -> 4 files -- plus every boundary
#: around them.
SIZING = [
    (1, 1, 1),
    (2, 2, 1),
    (3, 3, 1),
    (4, 4, 1),
    (7, 4, 1),
    (8, 4, 2),
    (12, 4, 3),
    (16, 4, 4),
    (32, 4, 8),
    (64, 4, 16),
]


@pytest.mark.parametrize(("cores", "threads", "jobs"), SIZING)
def test_sizing_rule(cores, threads, jobs):
    """4 cores -> 1 file, 8 -> 2 files, 16 -> 4 files, 4 threads each."""
    tuning = perf.recommend(cores, memory_bytes=64 * GB)
    assert (tuning.threads, tuning.jobs) == (threads, jobs), tuning.rationale


@pytest.mark.parametrize(("cores", "threads", "jobs"), SIZING)
def test_sizing_rule_is_identical_with_memory_unknown(cores, threads, jobs):
    """An unreadable memory figure must not change the sizing."""
    tuning = perf.recommend(cores, memory_bytes=None)
    # memory_bytes=None probes the host, which on a normal machine has plenty;
    # the clamp may only ever lower jobs, never raise it or touch threads.
    assert tuning.threads == threads
    assert 1 <= tuning.jobs <= jobs


def test_the_three_examples_the_user_gave():
    """The verbatim requirement: 8 cores -> 2 files, 16 cores -> 4 files."""
    assert perf.recommend(4, 64 * GB).jobs == 1
    assert perf.recommend(8, 64 * GB).jobs == 2
    assert perf.recommend(16, 64 * GB).jobs == 4
    for cores in (4, 8, 16):
        assert perf.recommend(cores, 64 * GB).threads == perf.THREADS_PER_FILE


@pytest.mark.parametrize("cores", range(1, 129))
def test_invariants_hold_at_every_core_count(cores):
    """threads >= 1, jobs >= 1, and threads*jobs never oversubscribes."""
    tuning = perf.recommend(cores, memory_bytes=256 * GB)
    assert tuning.threads >= 1
    assert tuning.jobs >= 1
    assert tuning.threads * tuning.jobs <= cores
    assert tuning.total_threads == tuning.threads * tuning.jobs


@pytest.mark.parametrize("cores", [4, 8, 12, 16, 32, 64, 128, 256])
def test_a_multiple_of_four_uses_every_core(cores):
    """On an exact multiple of 4 the tuning spends the whole machine."""
    tuning = perf.recommend(cores, memory_bytes=1024 * GB)
    assert tuning.threads * tuning.jobs == cores
    assert tuning.jobs == cores // 4


@pytest.mark.parametrize("cores", [1, 2, 3])
def test_below_four_cores_is_one_file_with_every_core(cores):
    """Under 4 cores: one file, threads = core count."""
    tuning = perf.recommend(cores, memory_bytes=64 * GB)
    assert tuning.jobs == 1
    assert tuning.threads == cores


@pytest.mark.parametrize("cores", [0, -1, -64])
def test_a_nonsense_core_count_floors_at_one(cores):
    """Zero or negative cores can never yield 0 threads or 0 jobs."""
    tuning = perf.recommend(cores, memory_bytes=64 * GB)
    assert tuning == perf.Tuning(1, 1, tuning.rationale)
    assert "1 core detected" in tuning.rationale


def test_a_non_integer_core_count_does_not_raise():
    """A GUI spinbox hands over strings; sizing must survive one."""
    assert perf.recommend("16", memory_bytes=64 * GB).jobs == 4
    assert perf.recommend(None, memory_bytes=64 * GB).jobs >= 1


# ---------------------------------------------------------------------------
# The memory clamp
# ---------------------------------------------------------------------------
def test_memory_clamps_jobs_down():
    """16 cores but 600 MB free -> 2 files, still 4 threads each."""
    tuning = perf.recommend(16, memory_bytes=600 * MB)
    assert tuning.jobs == 2
    assert tuning.threads == 4
    assert "memory" in tuning.rationale


def test_memory_never_raises_jobs():
    """A terabyte of RAM on 4 cores is still one file."""
    assert perf.recommend(4, memory_bytes=1024 * GB).jobs == 1
    assert perf.recommend(2, memory_bytes=1024 * GB).jobs == 1
    assert perf.recommend(7, memory_bytes=1024 * GB).jobs == 1


@pytest.mark.parametrize(
    ("free", "expected_jobs"),
    [
        (0, 1),                       # no headroom at all: one file still runs
        (1, 1),
        (perf.MEMORY_PER_FILE_BYTES, 1),
        (2 * perf.MEMORY_PER_FILE_BYTES, 2),
        (3 * perf.MEMORY_PER_FILE_BYTES, 3),
        (4 * perf.MEMORY_PER_FILE_BYTES, 4),
        (16 * perf.MEMORY_PER_FILE_BYTES, 4),   # clamp above the CPU answer
    ],
)
def test_memory_clamp_boundaries_on_sixteen_cores(free, expected_jobs):
    """The clamp is min(cpu answer, free // 300 MB), floored at 1."""
    assert perf.recommend(16, memory_bytes=free).jobs == expected_jobs


def test_memory_clamp_never_touches_threads():
    """Only jobs is clamped; each file still gets its four threads."""
    for free in (0, 1, 300 * MB, 900 * MB, 64 * GB):
        assert perf.recommend(64, memory_bytes=free).threads == 4


def test_a_starved_container_still_runs_one_file():
    """2 cores and 50 MB free is still a runnable 1x1 tuning."""
    tuning = perf.recommend(2, memory_bytes=50 * MB)
    assert (tuning.threads, tuning.jobs) == (2, 1)
    assert tuning.threads * tuning.jobs <= 2


def test_jobs_allowed_by_memory():
    """The helper is None-in/None-out and never returns 0."""
    assert perf.jobs_allowed_by_memory(None) is None
    assert perf.jobs_allowed_by_memory(0) == 1
    assert perf.jobs_allowed_by_memory(-1) == 1
    assert perf.jobs_allowed_by_memory("not a number") is None
    assert perf.jobs_allowed_by_memory(4 * perf.MEMORY_PER_FILE_BYTES) == 4
    assert perf.jobs_allowed_by_memory(4 * perf.MEMORY_PER_FILE_BYTES - 1) == 3


def test_memory_per_file_is_about_three_hundred_megabytes():
    """The documented figure; the GUI copy quotes it."""
    assert perf.MEMORY_PER_FILE_BYTES == 300 * MB


# ---------------------------------------------------------------------------
# The rationale sentence the GUI displays
# ---------------------------------------------------------------------------
def test_rationale_is_the_sentence_from_the_specification():
    """Exactly the example wording, so the GUI can show it verbatim."""
    tuning = perf.recommend(16, memory_bytes=64 * GB)
    assert tuning.rationale == "16 cores detected: 4 files at a time, 4 threads each"


def test_rationale_is_singular_on_one_core():
    """No '1 cores', no '1 files', no '1 threads'."""
    assert perf.recommend(1, memory_bytes=64 * GB).rationale == (
        "1 core detected: 1 file at a time, 1 thread each")


def test_rationale_names_the_memory_when_it_clamped():
    """A user must be able to see why they got fewer files than cores/4."""
    tuning = perf.recommend(32, memory_bytes=1200 * MB)
    assert tuning.jobs == 4
    assert tuning.rationale == (
        "32 cores detected but only 1.2 GB of memory free: "
        "4 files at a time, 4 threads each")


@pytest.mark.parametrize("cores", [1, 2, 3, 4, 7, 8, 12, 16, 32, 64])
def test_rationale_is_one_short_sentence(cores):
    """One line, no newline, short enough for a status bar."""
    for free in (None, 0, 512 * MB, 64 * GB):
        text = perf.recommend(cores, memory_bytes=free).rationale
        assert text and "\n" not in text and "\r" not in text
        assert len(text) <= 100
        assert text[0].isdigit()


# ---------------------------------------------------------------------------
# detect_cpu / detect_memory
# ---------------------------------------------------------------------------
def test_detect_cpu_never_returns_zero_or_none():
    """The contract: three integers, all >= 1."""
    info = perf.detect_cpu()
    assert isinstance(info, perf.CpuInfo)
    for value in (info.logical, info.physical_guess, info.usable):
        assert isinstance(value, int)
        assert value >= 1


def test_detect_cpu_respects_the_affinity_mask():
    """usable <= logical, and physical_guess is never invented above logical."""
    info = perf.detect_cpu()
    assert info.usable <= info.logical
    assert info.physical_guess <= info.logical
    affinity = getattr(os, "sched_getaffinity", None)
    if affinity is not None:
        assert info.usable <= len(affinity(0))


def test_detect_cpu_prefers_affinity_over_cpu_count(monkeypatch):
    """A taskset-restricted process must be sized for its mask, not the host."""
    monkeypatch.setattr(perf.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(perf, "_affinity_count", lambda: 8)
    monkeypatch.setattr(perf, "_cgroup_cpu_quota", lambda: None)
    info = perf.detect_cpu()
    assert info.logical == 64
    assert info.usable == 8
    assert perf.recommend(info.usable, memory_bytes=64 * GB).jobs == 2


def test_detect_cpu_respects_a_cgroup_quota(monkeypatch):
    """docker --cpus=4 on a 64-core host must size as four cores."""
    monkeypatch.setattr(perf.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(perf, "_affinity_count", lambda: 64)
    monkeypatch.setattr(perf, "_cgroup_cpu_quota", lambda: 4)
    info = perf.detect_cpu()
    assert info.usable == 4
    assert perf.recommend(info.usable, memory_bytes=64 * GB).jobs == 1


def test_detect_cpu_survives_a_broken_environment(monkeypatch):
    """cpu_count() None, no affinity, no cgroup: one core, not a crash."""
    monkeypatch.setattr(perf.os, "cpu_count", lambda: None)
    monkeypatch.setattr(perf, "_affinity_count", lambda: None)
    monkeypatch.setattr(perf, "_cgroup_cpu_quota", lambda: None)
    monkeypatch.setattr(perf, "_physical_cores_linux", lambda: None)
    monkeypatch.setattr(perf, "_physical_cores_windows", lambda: None)
    assert perf.detect_cpu() == perf.CpuInfo(1, 1, 1)


def test_a_nonsense_topology_falls_back_to_logical(monkeypatch):
    """A shared /proc reporting more physical cores than logical is ignored."""
    monkeypatch.setattr(perf.os, "cpu_count", lambda: 4)
    monkeypatch.setattr(perf, "_affinity_count", lambda: 4)
    monkeypatch.setattr(perf, "_cgroup_cpu_quota", lambda: None)
    monkeypatch.setattr(perf, "_physical_cores_linux", lambda: 96)
    monkeypatch.setattr(perf, "_physical_cores_windows", lambda: None)
    assert perf.detect_cpu().physical_guess == 4


def test_detect_memory_returns_plausible_numbers():
    """Either field may be None, but a number must be a positive int."""
    info = perf.detect_memory()
    assert isinstance(info, perf.MemoryInfo)
    for value in (info.total, info.available):
        assert value is None or (isinstance(value, int) and value >= 0)
    if info.total is not None and info.available is not None:
        assert info.available <= info.total


def test_unreadable_probes_degrade_to_unknown(monkeypatch):
    """Every memory source failing means "unknown", not an exception."""
    monkeypatch.setattr(perf, "_read_text", lambda path: None)
    monkeypatch.setattr(perf, "_read_int", lambda path: None)
    monkeypatch.setattr(perf, "_sysconf_memory", lambda: perf.MemoryInfo(None, None))
    if os.name != "nt":
        assert perf.detect_memory() == perf.MemoryInfo(None, None)
    # And a sizing with no memory information is pure CPU sizing.
    assert perf.recommend(16, memory_bytes=None).jobs == 4


def test_smallest_limit_wins(monkeypatch):
    """A cgroup limit beats the host figure /proc/meminfo reports."""
    monkeypatch.setattr(perf, "_meminfo_memory",
                        lambda: perf.MemoryInfo(64 * GB, 60 * GB))
    monkeypatch.setattr(perf, "_cgroup_memory",
                        lambda: perf.MemoryInfo(2 * GB, 1 * GB))
    if os.name != "nt":
        assert perf.detect_memory() == perf.MemoryInfo(2 * GB, 1 * GB)


def test_human_bytes():
    """The status line formatting used inside the rationale."""
    assert perf._human_bytes(1200 * MB) == "1.2 GB"
    assert perf._human_bytes(512 * MB) == "512 MB"
    assert perf._human_bytes(0) == "0 bytes"
    assert perf._human_bytes(-5) == "0 bytes"


def test_cgroup_v2_cpu_quota_is_parsed(tmp_path, monkeypatch):
    """`cpu.max` is "<quota> <period>"; a fraction rounds up, "max" is None."""
    path = tmp_path / "cpu.max"
    monkeypatch.setattr(perf, "_CGROUP_V2_CPU", str(path))
    monkeypatch.setattr(perf, "_CGROUP_V1_QUOTA", str(tmp_path / "absent"))
    monkeypatch.setattr(perf, "_CGROUP_V1_PERIOD", str(tmp_path / "absent"))

    path.write_text("400000 100000\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() == 4
    path.write_text("150000 100000\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() == 2
    path.write_text("max 100000\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() is None
    path.write_text("garbage\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() is None


def test_cgroup_v1_cpu_quota_is_parsed(tmp_path, monkeypatch):
    """cfs_quota_us / cfs_period_us, with -1 meaning unlimited."""
    quota = tmp_path / "cpu.cfs_quota_us"
    period = tmp_path / "cpu.cfs_period_us"
    monkeypatch.setattr(perf, "_CGROUP_V2_CPU", str(tmp_path / "absent"))
    monkeypatch.setattr(perf, "_CGROUP_V1_QUOTA", str(quota))
    monkeypatch.setattr(perf, "_CGROUP_V1_PERIOD", str(period))

    period.write_text("100000\n", encoding="utf-8")
    quota.write_text("800000\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() == 8
    quota.write_text("-1\n", encoding="utf-8")
    assert perf._cgroup_cpu_quota() is None


def test_meminfo_is_parsed(tmp_path, monkeypatch):
    """MemAvailable is preferred; kB is converted to bytes."""
    path = tmp_path / "meminfo"
    path.write_text(
        "MemTotal:       16305416 kB\n"
        "MemFree:          500000 kB\n"
        "MemAvailable:    8000000 kB\n"
        "Buffers:          100000 kB\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(perf, "_read_text",
                        lambda p: path.read_text(encoding="utf-8")
                        if p == "/proc/meminfo" else None)
    info = perf._meminfo_memory()
    assert info.total == 16305416 * 1024
    assert info.available == 8000000 * 1024


def test_read_int_handles_max_and_junk(tmp_path, monkeypatch):
    """cgroup files say "max" for "no limit"; that is not an integer."""
    path = tmp_path / "value"
    path.write_text("max\n", encoding="utf-8")
    assert perf._read_int(str(path)) is None
    path.write_text("12345\n", encoding="utf-8")
    assert perf._read_int(str(path)) == 12345
    assert perf._read_int(str(tmp_path / "absent")) is None


# ---------------------------------------------------------------------------
# probe(): the cheap start-up call
# ---------------------------------------------------------------------------
def test_probe_returns_a_usable_tuning():
    """What the GUI shows the moment its window opens."""
    tuning = perf.probe(refresh=True)
    assert isinstance(tuning, perf.Tuning)
    assert tuning.threads >= 1
    assert tuning.jobs >= 1
    assert tuning.rationale
    assert tuning.threads * tuning.jobs <= perf.detect_cpu().usable


def test_probe_is_cached_and_refreshable():
    """Repeated calls are free; refresh=True re-reads."""
    first = perf.probe(refresh=True)
    assert perf.probe() is first
    assert perf.probe(refresh=True) == first


def test_probe_launches_no_subprocess_and_runs_no_benchmark():
    """The start-up probe must not delay the window: no loop, no process."""
    source = (REPO_ROOT / "wmlst" / "perf.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert alias.name.split(".")[0] not in {
                    "subprocess", "multiprocessing", "timeit", "threading"}
        elif isinstance(node, ast.ImportFrom) and node.module:
            assert node.module.split(".")[0] not in {
                "subprocess", "multiprocessing", "timeit", "threading"}
    assert "time.sleep" not in source
    # No benchmark: nothing in this module measures elapsed time.
    assert "perf_counter" not in source


def test_probe_is_fast():
    """A hard ceiling, so a slow probe can never creep onto the GUI's path."""
    import time

    start = time.perf_counter()
    for _ in range(20):
        perf.probe(refresh=True)
    elapsed = time.perf_counter() - start
    assert elapsed < 1.0, "20 probes took %.3f s" % elapsed


# ---------------------------------------------------------------------------
# perf.py changes nothing about the CLI (docs/ARCHITECTURE.md C3)
# ---------------------------------------------------------------------------
def test_the_cli_defaults_are_untouched():
    """RunConfig stays threads=1, jobs=1: the golden files depend on it."""
    from wmlst import engine

    config = engine.RunConfig()
    assert config.threads == 1
    assert getattr(config, "jobs", 1) == 1


def test_no_engine_or_cli_module_imports_perf():
    """perf is GUI-only; the engine must stay a pure function of RunConfig."""
    for name in ("engine", "cli", "report", "schemes", "blastbin",
                 "any2fasta", "updatedb", "version", "branding"):
        path = REPO_ROOT / "wmlst" / (name + ".py")
        if not path.is_file():
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and (node.module or ""):
                assert not node.module.split(".")[-1] == "perf", (
                    "%s.py imports perf" % name)
                for alias in node.names:
                    assert alias.name != "perf", "%s.py imports perf" % name
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    assert not alias.name.endswith("perf"), (
                        "%s.py imports perf" % name)


def test_perf_imports_nothing_from_wmlst():
    """perf.py sits outside the section 2.1 rank graph: it depends on nothing."""
    tree = ast.parse((REPO_ROOT / "wmlst" / "perf.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert not node.level, "perf.py has a relative import"
            assert not (node.module or "").startswith("wmlst")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                assert not alias.name.startswith("wmlst")
