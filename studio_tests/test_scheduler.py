"""Resource arithmetic plus real overlapping/cancellable bounded operations."""

import threading

import pytest

from wmlstudio.scheduler import GIB, HardwareSnapshot, ResourcePlan, plan_resources, run_bounded
from wmlstudio.sequence import AnalysisCancelled


def hardware(cpus=16, available_gb=64):
    return HardwareSnapshot(cpus, int(available_gb * GIB), 64 * GIB, "test fixture")


def test_sixteen_cores_allocate_four_real_four_thread_jobs():
    plan = plan_resources(threads_per_sample=4, memory_gb=8, hardware=hardware())
    assert plan.max_parallel == 4
    assert plan.max_parallel * plan.threads_per_sample == 16
    assert plan.max_parallel * plan.memory_gb <= plan.memory_budget_gb


def test_memory_affinity_and_low_memory_policy_cap_admission():
    plan = plan_resources(threads_per_sample=4, memory_gb=8, hardware=hardware(16, 24))
    assert plan.max_parallel == 2
    assert plan_resources(threads_per_sample=4, memory_gb=3, hardware=hardware(2)).threads_per_sample == 2
    assert plan_resources(policy="low_memory", hardware=hardware()).max_parallel == 1
    assert plan_resources(hardware=HardwareSnapshot(32, None)).max_parallel == 1
    with pytest.raises(ValueError, match="safely available"):
        plan_resources(memory_gb=8, hardware=hardware(16, 8))


@pytest.mark.parametrize("field,value", [("cpu_budget", 0), ("threads_per_sample", True),
                                        ("memory_budget_gb", float("nan")), ("reserve_gb", float("inf"))])
def test_saved_invalid_resource_allocations_are_rejected(field, value):
    values = ResourcePlan(4, 8, 2, 16, 32, 2).to_dict()
    values[field] = value
    with pytest.raises(ValueError):
        ResourcePlan(**values).validate()


def test_real_tasks_overlap_within_cpu_and_ram_budget_and_callbacks_are_serial():
    plan = plan_resources(threads_per_sample=4, memory_gb=8, hardware=hardware())
    lock, barrier = threading.Lock(), threading.Barrier(4)
    active = maximum = 0
    coordinator = threading.get_ident()
    callbacks = []
    def operation(item, allocation, cancelled, progress):
        nonlocal active, maximum
        with lock:
            active += 1
            maximum = max(maximum, active)
            assert active * allocation.threads_per_sample <= allocation.cpu_budget
            assert active * allocation.memory_gb <= allocation.memory_budget_gb
        barrier.wait(timeout=5)
        progress(1, 1, f"complete {item}")
        with lock:
            active -= 1
        return item * 2
    def completed(item, value):
        assert threading.get_ident() == coordinator
        callbacks.append((item, value))
    output = run_bounded(range(8), operation, plan, on_result=completed, memory_probe=hardware)
    assert maximum == 4
    assert sorted(output) == [(index, index * 2) for index in range(8)]
    assert sorted(callbacks) == sorted(output)


def test_cancel_stops_pending_admission_and_drains_all_running_jobs():
    plan = plan_resources(threads_per_sample=4, memory_gb=8, max_parallel=2, hardware=hardware())
    barrier, cancel = threading.Barrier(2), threading.Event()
    started, failed = [], []
    def operation(item, allocation, cancelled, progress):
        barrier.wait(timeout=5)
        cancel.set()
        if cancelled():
            raise AnalysisCancelled()
    with pytest.raises(AnalysisCancelled):
        run_bounded(range(10), operation, plan, cancelled=cancel.is_set,
                    on_started=started.append, on_error=lambda item, error: failed.append(item), memory_probe=hardware)
    assert started == [0, 1]
    assert sorted(failed) == [0, 1]


def test_one_bad_sample_does_not_drop_other_results():
    plan = plan_resources(max_parallel=1, hardware=hardware())
    failed = []
    def operation(item, allocation, cancelled, progress):
        if item == "bad":
            raise ValueError("bad input")
        return item
    output = run_bounded(["bad", "good"], operation, plan,
                         on_error=lambda item, error: failed.append((item, str(error))), memory_probe=hardware)
    assert output == [("good", "good")]
    assert failed == [("bad", "bad input")]


def test_ram_pressure_before_start_fails_without_launching_any_input():
    plan = plan_resources(hardware=hardware())
    started = []
    with pytest.raises(RuntimeError, match="No new job"):
        run_bounded([1], lambda *args: None, plan, on_started=started.append,
                    memory_probe=lambda: hardware(16, 2))
    assert not started


def test_saved_plan_is_clamped_to_actual_cpu_affinity_at_launch():
    plan = plan_resources(threads_per_sample=8, hardware=hardware(32))
    observed = []
    def operation(item, allocation, cancelled, progress):
        observed.append(allocation)
    run_bounded([1], operation, plan, memory_probe=lambda: hardware(2))
    assert observed[0].threads_per_sample == 2
    assert observed[0].cpu_budget == 2 and observed[0].max_parallel == 1
