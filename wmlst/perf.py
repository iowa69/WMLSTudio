# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Cheap start-up auto-tuner: how many files at once, how many BLAST threads each.

**This module is for the GUI. It does NOT change any CLI default.**
``RunConfig.threads`` and ``RunConfig.jobs`` remain ``1``/``1`` everywhere
(docs/ARCHITECTURE.md C3), because the byte-compatibility contract of section 0
is verified against upstream's single-threaded invocation. Nothing here is
imported by :mod:`wmlst.cli` and nothing here is consulted by
:mod:`wmlst.engine`; the engine stays a pure function of the ``RunConfig`` it is
handed. The GUI calls :func:`probe` once when its window opens, shows the
numbers, and passes whatever the user finally chose explicitly.

The sizing rule, in one line: **one file per four cores, four BLAST threads
each.**

=========  =========  ======  ==========================================
cores      threads    jobs    meaning
=========  =========  ======  ==========================================
1          1          1       one file, one thread
2          2          1       one file, two threads
3          3          1       one file, three threads
4          4          1       one file, four threads
8          4          2       two files at a time, four threads each
16         4          4       four files at a time, four threads each
64         4          16      sixteen files at a time, four threads each
=========  =========  ======  ==========================================

Below four cores the machine cannot fill a whole slot, so it runs a single file
with as many threads as it has. ``threads * jobs`` never exceeds the usable core
count, and both are always at least ``1``.

``jobs`` is then clamped -- **only ever downwards** -- by free memory, at
:data:`MEMORY_PER_FILE_BYTES` (300 MB) per concurrent file. A 32-core container
holding 1 GB of headroom gets three files, not eight.

Cost: this module reads ``os.sched_getaffinity``/``os.cpu_count``, one or two
small pseudo-files under ``/proc`` and ``/sys/fs/cgroup``, or one ``ctypes``
call on Windows. There is no benchmark loop, no timing run and no subprocess,
so calling :func:`probe` on the GUI's start-up path costs well under a
millisecond and cannot delay the window appearing.

Every probe degrades to a safe answer rather than raising: an unreadable
``/proc``, a locked-down container, a missing ``sched_getaffinity`` and a
``cpu_count()`` of ``None`` all end at "one core, one file, one thread".

Byte sizes printed in :attr:`Tuning.rationale` use 1024-based units with the
familiar ``KB``/``MB``/``GB`` labels, matching how desktop tools label them.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

__all__ = [
    "CORES_PER_FILE",
    "MEMORY_PER_FILE_BYTES",
    "THREADS_PER_FILE",
    "CpuInfo",
    "MemoryInfo",
    "Tuning",
    "detect_cpu",
    "detect_memory",
    "jobs_allowed_by_memory",
    "probe",
    "recommend",
]

#: Cores consumed by one concurrent input file. The user's rule: multiples of 4.
CORES_PER_FILE = 4

#: ``blastn -num_threads`` for each concurrent file. Same number, by design:
#: a slot is four cores and it spends them on four threads.
THREADS_PER_FILE = 4

#: Rough resident cost of one concurrent file (BLAST database pages, the hit
#: table and the contig sequences). Deliberately generous: over-subscribing
#: memory makes a run swap, which is far worse than running one file fewer.
MEMORY_PER_FILE_BYTES = 300 * 1024 * 1024

#: cgroup v2 and v1 control files, read in that order.
_CGROUP_V2_CPU = "/sys/fs/cgroup/cpu.max"
_CGROUP_V1_QUOTA = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
_CGROUP_V1_PERIOD = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
_CGROUP_V2_MEM_MAX = "/sys/fs/cgroup/memory.max"
_CGROUP_V2_MEM_CUR = "/sys/fs/cgroup/memory.current"
_CGROUP_V1_MEM_MAX = "/sys/fs/cgroup/memory/memory.limit_in_bytes"
_CGROUP_V1_MEM_CUR = "/sys/fs/cgroup/memory/memory.usage_in_bytes"

#: A cgroup "no limit" sentinel is sometimes a huge number rather than "max".
_CGROUP_UNLIMITED = 1 << 62


# ---------------------------------------------------------------------------
# Value types
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CpuInfo:
    """What the machine has, and what this process is actually allowed to use.

    :param logical: logical processors the OS reports (hyperthreads included).
    :param physical_guess: best guess at physical cores. A *guess*: it equals
        ``logical`` whenever the topology cannot be read, which is the honest
        answer rather than a fabricated ``logical // 2``.
    :param usable: cores this process may really run on -- CPU affinity and any
        cgroup quota applied. This is the number :func:`recommend` sizes from.

    All three are integers ``>= 1``; none is ever ``0`` or ``None``.
    """

    logical: int
    physical_guess: int
    usable: int


@dataclass(frozen=True)
class MemoryInfo:
    """Physical memory, as far as it can be read. Either field may be ``None``.

    :param total: total physical (or cgroup-limited) bytes, or ``None``.
    :param available: bytes that could be handed out without swapping, or
        ``None`` when no source could be read.
    """

    total: Optional[int]
    available: Optional[int]


@dataclass(frozen=True)
class Tuning:
    """A concurrency recommendation, with a sentence explaining itself.

    :param threads: ``blastn -num_threads`` for each file. Always ``>= 1``.
    :param jobs: files to process concurrently. Always ``>= 1``.
    :param rationale: one short sentence, ready to drop into a GUI status line,
        e.g. ``"16 cores detected: 4 files at a time, 4 threads each"``.

    ``threads * jobs`` never exceeds the usable core count.
    """

    threads: int
    jobs: int
    rationale: str

    @property
    def total_threads(self) -> int:
        """Peak worker threads this tuning asks for: ``threads * jobs``."""
        return self.threads * self.jobs


# ---------------------------------------------------------------------------
# Small readers. Every one of them swallows its own errors.
# ---------------------------------------------------------------------------
def _read_text(path: str) -> Optional[str]:
    """Return the contents of a small pseudo-file, or ``None`` if unreadable."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            return handle.read()
    except (OSError, ValueError):
        return None


def _read_int(path: str) -> Optional[int]:
    """Return a whole-file integer, or ``None`` (also for the literal ``max``)."""
    text = _read_text(path)
    if text is None:
        return None
    text = text.strip()
    if not text or text == "max":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _cgroup_cpu_quota() -> Optional[int]:
    """Whole cores permitted by a cgroup quota, or ``None`` when unlimited.

    Rounds a fractional quota **up**, so a 1.5-core limit reports 2 rather than
    collapsing to 1: the quota throttles, it does not fail, and reporting 1 on a
    1.9-core allowance would waste nearly half the allocation.
    """
    text = _read_text(_CGROUP_V2_CPU)
    if text is not None:
        parts = text.split()
        if len(parts) == 2 and parts[0] != "max":
            try:
                quota, period = int(parts[0]), int(parts[1])
            except ValueError:
                quota = period = 0
            if quota > 0 and period > 0:
                return max(1, -(-quota // period))

    quota = _read_int(_CGROUP_V1_QUOTA)
    period = _read_int(_CGROUP_V1_PERIOD)
    if quota is not None and period is not None and quota > 0 and period > 0:
        return max(1, -(-quota // period))
    return None


def _affinity_count() -> Optional[int]:
    """Cores in this process's CPU affinity mask, or ``None`` where unsupported.

    ``os.sched_getaffinity`` exists on Linux and is the only call that sees a
    ``taskset``/``cpuset`` restriction; ``os.cpu_count()`` reports the whole
    machine regardless.
    """
    getter = getattr(os, "sched_getaffinity", None)
    if getter is None:
        return None
    try:
        return len(getter(0))
    except OSError:
        return None


def _physical_cores_linux() -> Optional[int]:
    """Physical cores from ``/proc/cpuinfo`` topology, or ``None``.

    Counts distinct ``(physical id, core id)`` pairs, which collapses the
    hyperthreads of one core into one entry.
    """
    text = _read_text("/proc/cpuinfo")
    if not text:
        return None
    pairs = set()
    package = core = None
    for line in text.splitlines():
        if ":" not in line:
            if package is not None and core is not None:
                pairs.add((package, core))
            package = core = None
            continue
        key, _, value = line.partition(":")
        key = key.strip()
        value = value.strip()
        if key == "physical id":
            package = value
        elif key == "core id":
            core = value
    if package is not None and core is not None:
        pairs.add((package, core))
    return len(pairs) or None


def _physical_cores_windows() -> Optional[int]:
    """Physical cores from ``GetLogicalProcessorInformation``, or ``None``.

    Counts ``SYSTEM_LOGICAL_PROCESSOR_INFORMATION`` records whose relationship
    is ``RelationProcessorCore`` (0). Any failure -- and any answer that fails
    the sanity check in :func:`detect_cpu` -- falls back to the logical count.
    """
    if os.name != "nt":
        return None
    try:
        import ctypes
        from ctypes import wintypes

        class _Union(ctypes.Union):
            _fields_ = [("Reserved", ctypes.c_ulonglong * 2)]

        class _Info(ctypes.Structure):
            _fields_ = [
                ("ProcessorMask", ctypes.POINTER(ctypes.c_ulong)),
                ("Relationship", wintypes.DWORD),
                ("Union", _Union),
            ]

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        length = wintypes.DWORD(0)
        kernel32.GetLogicalProcessorInformation(None, ctypes.byref(length))
        if length.value <= 0:
            return None
        count = length.value // ctypes.sizeof(_Info)
        if count <= 0:
            return None
        buffer = (_Info * count)()
        if not kernel32.GetLogicalProcessorInformation(
            ctypes.byref(buffer), ctypes.byref(length)
        ):
            return None
        cores = sum(1 for entry in buffer if entry.Relationship == 0)
        return cores or None
    except Exception:  # a cosmetic figure may never break start-up
        return None


def _windows_memory() -> MemoryInfo:
    """``(total, available)`` physical bytes via ``GlobalMemoryStatusEx``."""
    try:
        import ctypes
        from ctypes import wintypes

        class _MemStatus(ctypes.Structure):
            _fields_ = [
                ("dwLength", wintypes.DWORD),
                ("dwMemoryLoad", wintypes.DWORD),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        status = _MemStatus()
        status.dwLength = ctypes.sizeof(_MemStatus)
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        if not kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return MemoryInfo(None, None)
        return MemoryInfo(int(status.ullTotalPhys), int(status.ullAvailPhys))
    except Exception:  # an unreadable probe means "unknown"
        return MemoryInfo(None, None)


def _meminfo_memory() -> MemoryInfo:
    """``(total, available)`` from ``/proc/meminfo``, in bytes."""
    text = _read_text("/proc/meminfo")
    if not text:
        return MemoryInfo(None, None)
    values = {}
    for line in text.splitlines():
        key, _, rest = line.partition(":")
        fields = rest.split()
        if not fields:
            continue
        try:
            amount = int(fields[0])
        except ValueError:
            continue
        unit = fields[1].lower() if len(fields) > 1 else "kb"
        values[key.strip()] = amount * 1024 if unit == "kb" else amount
    total = values.get("MemTotal")
    available = values.get("MemAvailable")
    if available is None:
        free = values.get("MemFree")
        if free is not None:
            # Pre-3.14 kernels have no MemAvailable; reclaimable cache is the
            # documented approximation.
            available = free + values.get("Cached", 0) + values.get("Buffers", 0)
    return MemoryInfo(total, available)


def _sysconf_memory() -> MemoryInfo:
    """``(total, available)`` from ``os.sysconf`` page counts, in bytes."""
    sysconf = getattr(os, "sysconf", None)
    names = getattr(os, "sysconf_names", {})
    if sysconf is None:
        return MemoryInfo(None, None)

    def _get(name: str) -> Optional[int]:
        if name not in names:
            return None
        try:
            value = sysconf(name)
        except (OSError, ValueError):
            return None
        return value if isinstance(value, int) and value > 0 else None

    page = _get("SC_PAGE_SIZE") or _get("SC_PAGESIZE")
    if page is None:
        return MemoryInfo(None, None)
    total_pages = _get("SC_PHYS_PAGES")
    free_pages = _get("SC_AVPHYS_PAGES")
    return MemoryInfo(
        total_pages * page if total_pages else None,
        free_pages * page if free_pages else None,
    )


def _cgroup_memory() -> MemoryInfo:
    """``(limit, limit - current)`` from a cgroup memory controller."""
    limit = _read_int(_CGROUP_V2_MEM_MAX)
    current = _read_int(_CGROUP_V2_MEM_CUR)
    if limit is None:
        limit = _read_int(_CGROUP_V1_MEM_MAX)
        current = _read_int(_CGROUP_V1_MEM_CUR)
    if limit is None or limit <= 0 or limit >= _CGROUP_UNLIMITED:
        return MemoryInfo(None, None)
    if current is None or current < 0:
        return MemoryInfo(limit, None)
    return MemoryInfo(limit, max(0, limit - current))


def _human_bytes(count: int) -> str:
    """Format a byte count for a status line: ``"1.2 GB"``, ``"512 MB"``."""
    count = max(0, int(count))
    for unit, size in (("GB", 1024 ** 3), ("MB", 1024 ** 2), ("KB", 1024)):
        if count >= size:
            scaled = count / float(size)
            if scaled >= 100:
                return "%d %s" % (int(scaled), unit)
            return "%.1f %s" % (scaled, unit)
    return "%d bytes" % count


# ---------------------------------------------------------------------------
# The public probes
# ---------------------------------------------------------------------------
def detect_cpu() -> CpuInfo:
    """Return :class:`CpuInfo` for this machine and this process.

    ``usable`` respects **both** the CPU affinity mask (``taskset``, a Docker
    ``--cpuset-cpus``, a scheduler pin) and a cgroup CPU quota (a Kubernetes
    ``limits.cpu``, ``docker --cpus``); ``os.cpu_count()`` sees neither and is
    only the fallback. The result is never ``0`` and never ``None``.
    """
    logical = os.cpu_count() or 0
    if logical < 1:
        logical = _affinity_count() or 1
    logical = max(1, logical)

    usable = _affinity_count()
    if usable is None or usable < 1:
        usable = logical
    quota = _cgroup_cpu_quota()
    if quota is not None and quota >= 1:
        usable = min(usable, quota)
    usable = max(1, min(usable, logical))

    physical = _physical_cores_linux() if os.name == "posix" else None
    if physical is None:
        physical = _physical_cores_windows()
    # A topology answer outside [1, logical] is nonsense from a container's
    # shared /proc; prefer the honest logical count to a fabricated one.
    if physical is None or physical < 1 or physical > logical:
        physical = logical

    return CpuInfo(logical=logical, physical_guess=physical, usable=usable)


def detect_memory() -> MemoryInfo:
    """Return :class:`MemoryInfo`, using the tightest limit any probe reports.

    Sources, all cheap: ``GlobalMemoryStatusEx`` on Windows, ``/proc/meminfo``
    and the cgroup memory controller on Linux, ``os.sysconf`` page counts
    elsewhere. When two sources disagree the smaller wins, because a cgroup
    limit is real even though ``/proc/meminfo`` shows the whole host. Both
    fields are ``None`` when nothing could be read.
    """
    if os.name == "nt":
        return _windows_memory()

    best = _meminfo_memory()
    if best.total is None and best.available is None:
        best = _sysconf_memory()

    cgroup = _cgroup_memory()
    total = _smallest(best.total, cgroup.total)
    available = _smallest(best.available, cgroup.available)
    return MemoryInfo(total, available)


def _smallest(left: Optional[int], right: Optional[int]) -> Optional[int]:
    """The smaller of two optional byte counts, ignoring ``None``."""
    if left is None:
        return right
    if right is None:
        return left
    return min(left, right)


def jobs_allowed_by_memory(memory_bytes: Optional[int]) -> Optional[int]:
    """How many concurrent files ``memory_bytes`` affords, or ``None``.

    ``None`` in, ``None`` out: unknown memory must never influence the sizing.
    The answer is never below ``1`` -- one file always runs, even on a machine
    reporting no headroom at all, because refusing to work is not an option a
    typing tool has.
    """
    if memory_bytes is None:
        return None
    try:
        usable = int(memory_bytes)
    except (TypeError, ValueError):
        return None
    if usable < 0:
        usable = 0
    return max(1, usable // MEMORY_PER_FILE_BYTES)


def recommend(
    cpu_count: Optional[int] = None,
    memory_bytes: Optional[int] = None,
) -> Tuning:
    """Size a run: one file per four cores, four BLAST threads each.

    :param cpu_count: cores to size for. ``None`` probes the machine with
        :func:`detect_cpu` and uses ``usable``. Values below ``1`` are raised
        to ``1``.
    :param memory_bytes: free memory to size against. ``None`` probes with
        :func:`detect_memory`; when that also comes back unknown, no memory
        clamp is applied at all.
    :returns: a :class:`Tuning`.

    Four cores give one file with four threads; every further complete group of
    four cores adds another concurrent file. Under four cores, one file takes
    every core as threads. Memory can only ever *reduce* ``jobs``, never raise
    it, and ``threads * jobs`` never exceeds ``cpu_count``.
    """
    if cpu_count is None:
        cores = detect_cpu().usable
    else:
        try:
            cores = int(cpu_count)
        except (TypeError, ValueError):
            cores = 1
    cores = max(1, cores)

    if cores < CORES_PER_FILE:
        threads = cores
        jobs = 1
    else:
        threads = THREADS_PER_FILE
        jobs = cores // CORES_PER_FILE

    if memory_bytes is None:
        memory_bytes = detect_memory().available

    allowed = jobs_allowed_by_memory(memory_bytes)
    clamped_to = None
    if allowed is not None and allowed < jobs:
        jobs = allowed
        clamped_to = memory_bytes

    # Invariants, restated as code rather than trusted to the branches above.
    threads = max(1, min(threads, cores))
    jobs = max(1, jobs)
    if threads * jobs > cores:
        jobs = max(1, cores // threads)

    return Tuning(threads=threads, jobs=jobs,
                  rationale=_rationale(cores, threads, jobs, clamped_to))


def _rationale(cores: int, threads: int, jobs: int,
               clamped_to: Optional[int]) -> str:
    """One sentence a status bar can show verbatim."""
    core_word = "core" if cores == 1 else "cores"
    file_word = "file" if jobs == 1 else "files"
    thread_word = "thread" if threads == 1 else "threads"
    tail = "%d %s at a time, %d %s each" % (jobs, file_word, threads, thread_word)
    if clamped_to is not None:
        return "%d %s detected but only %s of memory free: %s" % (
            cores, core_word, _human_bytes(clamped_to), tail)
    return "%d %s detected: %s" % (cores, core_word, tail)


#: Memoised result of the start-up probe. The core count and the memory limit
#: do not change under a running GUI in any way worth re-reading on every
#: repaint; ``probe(refresh=True)`` re-reads on demand.
_CACHED_PROBE: Optional[Tuning] = None


def probe(refresh: bool = False) -> Tuning:
    """The start-up call: detect, size, and return a :class:`Tuning`.

    This is what the GUI runs when its window opens. It reads a core count and
    a memory figure and does nothing else -- no benchmark, no timing loop, no
    subprocess -- so it cannot delay the window. The result is cached; pass
    ``refresh=True`` to probe again.
    """
    global _CACHED_PROBE
    if refresh or _CACHED_PROBE is None:
        cpu = detect_cpu()
        _CACHED_PROBE = recommend(cpu.usable, detect_memory().available)
    return _CACHED_PROBE
