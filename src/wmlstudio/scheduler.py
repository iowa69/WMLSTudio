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
    parallel = min(cpus // threads, int(budget // memory_gb), max_parallel or cpus)
    if policy == "low_memory" or hardware.available_memory is None:
        parallel = min(parallel, 1)
    if parallel < 1:
        raise ValueError(f"Only {budget:.1f} GiB is safely available for jobs; each sample requests {memory_gb} GiB. Reduce per-sample memory or free RAM.")
    return ResourcePlan(threads, memory_gb, parallel, cpus, round(budget, 3),
                        round(reserve, 3), policy, hardware.available_memory is not None).validate()


def resource_plan(value=None, **defaults):
    if isinstance(value, ResourcePlan):
        return value.validate()
    if value:
        return ResourcePlan(**value).validate()
    return plan_resources(**defaults)


def resources_for_run(plan=None, *, memory_gb=1):
    plan = plan or {}
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
