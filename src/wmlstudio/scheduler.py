"""Portable CPU/RAM admission control for cancellable per-isolate operations.

Memory is a conservative reservation, not an operating-system process limit.
Callbacks run serially on the coordinator; operations must not access Qt or a
project database. Every operation receives its explicit native-tool allocation.
"""

from __future__ import annotations

import math
import os
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from wmlstudio.sequence import AnalysisCancelled

GIB = 1024 ** 3
#: The least memory a single sample is ever scheduled with. Below this the work
#: genuinely cannot proceed, so the user is told rather than left with a run that
#: stops without a result.
MINIMUM_SAMPLE_MEMORY_GB = 2
#: How far a per-sample request may be scaled down to fit the machine. Within this
#: the request is read as a generous default; beyond it, as a genuine requirement.
REDUCIBLE_MEMORY_FACTOR = 3


@dataclass(frozen=True)
class HardwareSnapshot:
    cpus: int
    available_memory: int | None
    total_memory: int | None = None
    source: str = "detected"


def detect_hardware():
    cpus = max(1, os.cpu_count() or 1)
    available = total = None
    try:
        cpus = min(cpus, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    if os.name == "nt":
        try:
            import ctypes
            from ctypes import wintypes

            class MemoryStatus(ctypes.Structure):
                _fields_ = [("length", wintypes.DWORD), ("load", wintypes.DWORD),
                            *[(name, ctypes.c_ulonglong) for name in
                              ("total", "available", "page_total", "page_available",
                               "virtual_total", "virtual_available", "extended")]]
            status = MemoryStatus()
            status.length = ctypes.sizeof(status)
            kernel = ctypes.windll.kernel32
            kernel.GlobalMemoryStatusEx.argtypes = [ctypes.POINTER(MemoryStatus)]
            kernel.GlobalMemoryStatusEx.restype = wintypes.BOOL
            if kernel.GlobalMemoryStatusEx(ctypes.byref(status)):
                total, available = int(status.total), int(status.available)
            process_mask, system_mask = ctypes.c_size_t(), ctypes.c_size_t()
            kernel.GetCurrentProcess.restype = wintypes.HANDLE
            kernel.GetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
            kernel.GetProcessAffinityMask.restype = wintypes.BOOL
            if kernel.GetProcessAffinityMask(kernel.GetCurrentProcess(), ctypes.byref(process_mask), ctypes.byref(system_mask)):
                cpus = min(cpus, process_mask.value.bit_count() or 1)
        except (AttributeError, OSError, ValueError):
            pass
    else:
        try:
            values = {line.split(":")[0]: int(line.split()[1]) * 1024
                      for line in Path("/proc/meminfo").read_text().splitlines() if ":" in line}
            total, available = values.get("MemTotal"), values.get("MemAvailable")
        except (OSError, ValueError, IndexError):
            pass
        # Respect common cgroup-v2 limits as well as process CPU affinity.
        try:
            quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
            if quota != "max":
                cpus = min(cpus, max(1, int(quota) // int(period)))
        except (OSError, ValueError, ZeroDivisionError):
            pass
        try:
            limit = Path("/sys/fs/cgroup/memory.max").read_text().strip()
            if limit != "max":
                current = int(Path("/sys/fs/cgroup/memory.current").read_text())
                available = min(available if available is not None else int(limit), max(0, int(limit) - current))
                total = min(total if total is not None else int(limit), int(limit))
        except (OSError, ValueError):
            pass
    return HardwareSnapshot(max(1, cpus), available, total)


@dataclass(frozen=True)
class ResourcePlan:
    threads_per_sample: int
    memory_gb: int
    max_parallel: int
    cpu_budget: int
    memory_budget_gb: float
    reserve_gb: float
    policy: str = "balanced"
    memory_detected: bool = True
    #: True when the per-sample memory request was lowered to fit this computer.
    reduced_for_memory: bool = False

    def to_dict(self):
        return asdict(self)

    def validate(self):
        for name in ("threads_per_sample", "memory_gb", "max_parallel", "cpu_budget"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("memory_budget_gb", "reserve_gb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.max_parallel * self.threads_per_sample > self.cpu_budget:
            raise ValueError("Concurrent samples exceed the CPU budget")
        if self.max_parallel * self.memory_gb > self.memory_budget_gb:
            raise ValueError("Concurrent samples exceed the memory reservation budget")
        if self.reserve_gb < 0:
            raise ValueError("Reserved memory must not be negative")
        return self


def plan_resources(*, threads_per_sample=4, memory_gb=8, cpu_budget=None,
                   memory_budget_gb=None, max_parallel=None, policy="balanced", hardware=None):
    hardware = hardware or detect_hardware()
    if policy not in {"balanced", "fast", "low_memory"}:
        raise ValueError("Unknown resource policy")
    for name, value in (("threads_per_sample", threads_per_sample), ("memory_gb", memory_gb),
                        ("cpu_budget", cpu_budget), ("max_parallel", max_parallel)):
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 1):
            raise ValueError(f"{name} must be a positive integer")
    cpus = min(hardware.cpus, cpu_budget or hardware.cpus)
    threads = min(threads_per_sample, cpus)
    reserve = 0.0
    if hardware.available_memory is not None:
        available = hardware.available_memory / GIB
        reserve = max(1.0, available * {"balanced": .20, "fast": .10, "low_memory": .25}[policy])
        budget = max(0.0, available - reserve)
    else:
        budget = float(memory_gb)  # Unknown RAM: one job only, never guessed parallelism.
    if memory_budget_gb is not None:
        if isinstance(memory_budget_gb, bool) or not isinstance(memory_budget_gb, (float, int)) or not 0 < memory_budget_gb < 1e6:
            raise ValueError("Memory budget must be a finite positive GiB value")
        budget = min(budget, float(memory_budget_gb))
    # A modest laptop must still be able to work. Asking for 8 GiB per sample on a
    # machine with 6 GiB free used to raise, the assembly task failed, and because
    # typing only starts after a successful assembly the run stopped with no result
    # and no explanation the user could connect to memory. Run one smaller job
    # instead of refusing, and say what was reduced.
    if budget < memory_gb:
        # A modest shortfall means the default per-sample figure was a little
        # generous for this computer, so run one smaller job. A request many times
        # larger than the machine is a real requirement that cannot be met, and
        # quietly handing it a fraction would only fail later and less clearly.
        if budget < MINIMUM_SAMPLE_MEMORY_GB or memory_gb > budget * REDUCIBLE_MEMORY_FACTOR:
            raise ValueError(
                f"Only {budget:.1f} GiB is safely available for jobs and each sample requests "
                f"{memory_gb} GiB. Close other programs to free memory, or lower the memory "
                f"each sample may use.")
    per_sample = min(memory_gb, max(MINIMUM_SAMPLE_MEMORY_GB, int(budget)))
    parallel = min(cpus // threads, int(budget // per_sample), max_parallel or cpus)
    if policy == "low_memory" or hardware.available_memory is None:
        parallel = min(parallel, 1)
    parallel = max(1, parallel)
    return ResourcePlan(threads, per_sample, parallel, cpus, round(budget, 3),
                        round(reserve, 3), policy, hardware.available_memory is not None,
                        per_sample < memory_gb).validate()


def resource_plan(value=None, **defaults):
    if isinstance(value, ResourcePlan):
        return value.validate()
    if value:
        return ResourcePlan(**value).validate()
    return plan_resources(**defaults)


# --- how much of this computer to use, said in two numbers -----------------
# Two explicit numbers instead of a policy name, which is the model WMLST uses
# and people read without help: JOBS is how many samples run at the same time,
# THREADS is how many CPU threads each of those samples may hand to its native
# tools. Both are clamped to [1, min(CPU count, 64)] -- 64 is the Windows
# processor-group boundary -- and jobs x threads is clamped to the CPU count,
# because asking for more threads than the machine has makes every sample
# slower rather than faster. Nothing here changes a result: it is throughput.

#: Neither number may exceed this, whatever the machine reports.
THREAD_LIMIT = 64
#: Automatic sizing: threads one concurrent sample occupies, and the threads it
#: is then told to use. The same number by design -- a slot is four threads and
#: it spends them on four threads.
CPUS_PER_JOB = 4
THREADS_PER_JOB = 4


@dataclass(frozen=True)
class WorkSize:
    """Samples at a time, threads each, and the machine the numbers were sized for."""

    jobs: int
    threads: int
    cpus: int
    automatic: bool = False
    # "", "cpu" or "memory": what actually held these numbers down, so the
    # interface can say why rather than showing an unexplained smaller number.
    limited_by: str = ""

    @property
    def total_threads(self) -> int:
        return self.jobs * self.threads

    def to_dict(self):
        return asdict(self)


def cpu_count(hardware=None) -> int:
    """CPU threads this process may actually use. Never 0, never None."""
    return max(1, (hardware or detect_hardware()).cpus)


def thread_ceiling(cpus=None) -> int:
    """The largest value either number may take on this computer."""
    return max(1, min(_whole(cpus, 0) or cpu_count(), THREAD_LIMIT))


def _whole(value, fallback=1) -> int:
    """Coerce a spin box or a stored preference to a whole number."""
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return int(fallback)


def clamp_worksize(jobs, threads, cpus=None, *, automatic=False, limited_by="") -> WorkSize:
    """Both numbers into range, then jobs x threads into the CPU count."""
    cpus = max(1, _whole(cpus, 0) or cpu_count())
    ceiling = thread_ceiling(cpus)
    jobs = min(ceiling, max(1, _whole(jobs)))
    threads = min(ceiling, max(1, _whole(threads)))
    if jobs * threads > cpus:
        # The per-sample threads come down, never the number of samples: the
        # user asked for that many samples at once and still gets them.
        threads = max(1, cpus // jobs)
        limited_by = limited_by or "cpu"
    return WorkSize(jobs, threads, cpus, bool(automatic), str(limited_by))


def auto_worksize(hardware=None, *, memory_gb=1) -> WorkSize:
    """Size this computer: one sample per four CPU threads, four threads each.

    Below four threads the machine cannot fill a whole slot, so one sample takes
    every thread. Free memory only ever reduces the number of samples, and an
    unreadable memory figure means one sample at a time rather than a guess.
    """
    hardware = hardware or detect_hardware()
    cpus = max(1, hardware.cpus)
    if cpus < CPUS_PER_JOB:
        size = clamp_worksize(1, cpus, cpus, automatic=True)
    else:
        size = clamp_worksize(cpus // CPUS_PER_JOB, THREADS_PER_JOB, cpus, automatic=True)
    if hardware.available_memory is None:
        return replace(size, jobs=1, limited_by="memory")
    available = hardware.available_memory / GIB
    budget = max(0.0, available - max(1.0, available * .20))
    affordable = int(budget // max(1, memory_gb))
    if affordable < size.jobs:
        size = replace(size, jobs=max(1, affordable), limited_by="memory")
    return size


def describe_worksize(size) -> str:
    """The plain sentence the Settings page shows: what this choice actually means."""
    samples = "sample" if size.jobs == 1 else "samples"
    threads = "thread" if size.threads == 1 else "threads"
    sentence = (f"{size.jobs} {samples} at a time, {size.threads} {threads} each, using "
                f"{size.total_threads} of your {size.cpus} CPU threads.")
    if size.limited_by == "memory":
        sentence += " Free memory, not the processor, is what holds this down."
    elif size.limited_by == "cpu":
        sentence += (" The threads for each sample came down so that the two numbers "
                     "together fit this computer.")
    return sentence


def worksize_constraint(cpus=None) -> str:
    """The rule in words, to stand beside the two numbers it governs."""
    cpus = max(1, _whole(cpus, 0) or cpu_count())
    return (f"Samples at a time × threads each can never be more than the {cpus} CPU threads "
            f"this computer has. Ask for more and the threads for each sample are reduced, "
            f"never the samples. Each number may be 1 to {thread_ceiling(cpus)}.")


def plan_for(jobs, threads=None, *, memory_gb=1, policy="balanced", hardware=None):
    """Turn the two numbers into the admission plan run_bounded already enforces."""
    hardware = hardware or detect_hardware()
    size = jobs if isinstance(jobs, WorkSize) else clamp_worksize(jobs, threads, hardware.cpus)
    return plan_resources(threads_per_sample=size.threads, memory_gb=memory_gb,
                          max_parallel=size.jobs, policy=policy, hardware=hardware)


#: The Settings page's two numbers, for runs that state none of their own. None
#: means "nothing has been chosen", which is why importing this module changes
#: no existing behaviour: only a caller that sets it is affected.
_DEFAULT_WORKSIZE = None


def set_default_worksize(size):
    """Record the chosen size application-wide, or clear it with None."""
    global _DEFAULT_WORKSIZE
    if size is None:
        _DEFAULT_WORKSIZE = None
    elif isinstance(size, WorkSize):
        _DEFAULT_WORKSIZE = clamp_worksize(size.jobs, size.threads, size.cpus,
                                           automatic=size.automatic)
    else:
        _DEFAULT_WORKSIZE = clamp_worksize(size["jobs"], size["threads"], size.get("cpus"),
                                           automatic=size.get("automatic", False))
    return _DEFAULT_WORKSIZE


def default_worksize():
    """The chosen size, or None when the user has never chosen one."""
    return _DEFAULT_WORKSIZE


def resources_for_run(plan=None, *, memory_gb=1):
    plan = plan or {}
    if _DEFAULT_WORKSIZE is not None and not plan.get("resource_plan") and "threads" not in plan:
        # This run named no resources of its own, so the Settings page answers.
        return plan_for(_DEFAULT_WORKSIZE, memory_gb=plan.get("memory_gb", memory_gb))
    return resource_plan(plan.get("resource_plan"), threads_per_sample=plan.get("threads", 4),
                         memory_gb=plan.get("memory_gb", memory_gb), policy=plan.get("resource_policy", "balanced"))


def run_bounded(items, operation, plan, *, cancelled=None, on_started=None,
                on_result=None, on_error=None, progress=None, memory_probe=None):
    """Run bounded work; deliver callbacks serially and retain completed results.

    operation(item, allocation, cancellation_callback, progress_callback) must
    pass allocation.threads_per_sample to native tools and observe cancellation.
    progress_callback accepts (done, total, message). No later sample is admitted
    after cancellation. Existing jobs are drained before this function returns.
    """
    plan = resource_plan(plan)
    items = list(items)
    if not items:
        return []
    memory_probe = memory_probe or detect_hardware
    initial = memory_probe()
    cpu_budget = min(plan.cpu_budget, initial.cpus)
    threads = min(plan.threads_per_sample, cpu_budget)
    memory_budget = min(plan.memory_budget_gb, max(0.0, initial.available_memory / GIB - plan.reserve_gb)) if initial.available_memory is not None else plan.memory_budget_gb
    parallel = min(plan.max_parallel, cpu_budget // threads, int(memory_budget // plan.memory_gb))
    if parallel < 1:
        raise RuntimeError("Available RAM fell below the reviewed per-sample allocation. No new job was started; free memory and retry.")
    plan = replace(plan, cpu_budget=cpu_budget, threads_per_sample=threads, max_parallel=parallel,
                   memory_budget_gb=memory_budget).validate()
    stop = threading.Event()
    def is_cancelled():
        return stop.is_set() or bool(cancelled and cancelled())
    completed, running, results = 0, {}, []
    next_index = 0
    message_lock = threading.Lock()
    notices = []
    def report(done, total, message):
        with message_lock:
            notices[:] = [str(message)]
    with ThreadPoolExecutor(max_workers=min(plan.max_parallel, len(items)), thread_name_prefix="WMLSTudio-isolate") as pool:
        try:
            while next_index < len(items) or running:
                while not is_cancelled() and next_index < len(items) and len(running) < plan.max_parallel:
                    observed = memory_probe()
                    if observed.available_memory is not None and observed.available_memory / GIB < plan.memory_gb + plan.reserve_gb:
                        if not running:
                            raise RuntimeError("Available RAM fell below the reviewed per-sample allocation. No new job was started; free memory and retry.")
                        break
                    item = items[next_index]
                    if on_started:
                        on_started(item)
                    future = pool.submit(operation, item, plan, is_cancelled, report)
                    running[future] = item
                    next_index += 1
                if not running:
                    break
                ready, _ = wait(running, timeout=.1, return_when=FIRST_COMPLETED)
                with message_lock:
                    notice = notices.pop() if notices else None
                if progress and notice:
                    progress(completed, len(items), notice)
                for future in sorted(ready, key=lambda f: items.index(running[f])):
                    item = running.pop(future)
                    try:
                        result = future.result()
                    except AnalysisCancelled:
                        stop.set()
                        if on_error:
                            on_error(item, AnalysisCancelled())
                    except Exception as error:
                        if on_error:
                            on_error(item, error)
                        else:
                            raise
                    else:
                        results.append((item, result))
                        if on_result:
                            on_result(item, result)
                    completed += 1
                    if progress:
                        progress(completed, len(items), f"{completed}/{len(items)} samples finished · up to {plan.max_parallel} concurrent")
            if is_cancelled():
                raise AnalysisCancelled("Queue cancelled; completed results are retained and pending samples were not started.")
        finally:
            stop.set()
    return results
