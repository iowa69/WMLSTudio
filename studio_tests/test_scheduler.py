"""Resource arithmetic plus real overlapping/cancellable bounded operations."""

import threading

import pytest

from wmlstudio import scheduler
from wmlstudio.scheduler import GIB, HardwareSnapshot, ResourcePlan, plan_resources, run_bounded
from wmlstudio.sequence import AnalysisCancelled


def hardware(cpus=16, available_gb=64):
    return HardwareSnapshot(cpus, int(available_gb * GIB), 64 * GIB, "test fixture")


@pytest.fixture
def no_default_worksize():
    """The chosen size is application-wide; no test may leak one into the next."""
    previous = scheduler.default_worksize()
    scheduler.set_default_worksize(None)
    yield
    scheduler.set_default_worksize(previous)


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
    # This used to refuse. A capable machine with 8 GiB free is exactly the laptop
    # this software is meant to run on, and refusing it stopped the pipeline with
    # no result, so the default request is now trimmed to fit instead.
    fitted = plan_resources(memory_gb=8, hardware=hardware(16, 8))
    assert fitted.max_parallel >= 1 and fitted.reduced_for_memory is True
    assert fitted.max_parallel * fitted.memory_gb <= fitted.memory_budget_gb


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


# ---------------------------------------------------------------------------
# The two numbers a person actually sets: samples at a time, threads each.
# This is the legacy WMLST model, which the user asked us to copy, so the rule
# it encodes is asserted here rather than left to the widget that shows it.
# ---------------------------------------------------------------------------


def test_both_numbers_stay_inside_one_and_the_processor_group_boundary():
    assert scheduler.thread_ceiling(16) == 16
    assert scheduler.thread_ceiling(128) == scheduler.THREAD_LIMIT == 64
    size = scheduler.clamp_worksize(0, 0, 8)
    assert (size.jobs, size.threads) == (1, 1)
    assert scheduler.clamp_worksize(999, 1, 8).jobs == 8
    assert scheduler.clamp_worksize("4", "2.0", 16).to_dict() == {
        "jobs": 4, "threads": 2, "cpus": 16, "automatic": False, "limited_by": ""}
    assert scheduler.clamp_worksize(None, None, 16).jobs == 1


def test_asking_for_more_than_the_machine_has_reduces_threads_not_samples():
    """WMLST's own rule: jobs x threads is clamped to the CPU count by threads."""
    size = scheduler.clamp_worksize(8, 8, 16)
    assert (size.jobs, size.threads, size.total_threads) == (8, 2, 16)
    assert size.limited_by == "cpu"
    assert scheduler.clamp_worksize(4, 4, 16).limited_by == ""
    for jobs, threads, cpus in [(3, 7, 9), (5, 5, 16), (2, 64, 4), (7, 3, 7)]:
        size = scheduler.clamp_worksize(jobs, threads, cpus)
        assert size.jobs * size.threads <= cpus, (jobs, threads, cpus)
        assert size.jobs >= 1 and size.threads >= 1


@pytest.mark.parametrize("cpus,expected", [(1, (1, 1)), (2, (1, 2)), (3, (1, 3)), (4, (1, 4)),
                                           (8, (2, 4)), (16, (4, 4)), (64, (16, 4))])
def test_automatic_sizing_is_one_sample_per_four_threads_four_threads_each(cpus, expected):
    size = scheduler.auto_worksize(hardware(cpus))
    assert (size.jobs, size.threads) == expected
    assert size.automatic is True
    assert size.total_threads <= cpus


def test_memory_only_ever_reduces_the_number_of_samples():
    tight = scheduler.auto_worksize(hardware(16, available_gb=3), memory_gb=1)
    assert tight.jobs == 2 and tight.threads == 4, "four slots fit the CPU, two fit the RAM"
    assert tight.limited_by == "memory"
    scarce = scheduler.auto_worksize(hardware(16, available_gb=1.5), memory_gb=1)
    assert scarce.jobs == 1 and scarce.limited_by == "memory"
    # Unknown memory is never guessed at: one sample at a time.
    unknown = scheduler.auto_worksize(HardwareSnapshot(32, None))
    assert unknown.jobs == 1 and unknown.limited_by == "memory"
    assert scheduler.auto_worksize(hardware(16, 64)).jobs == 4


def test_the_choice_is_described_in_one_plain_sentence():
    sentence = scheduler.describe_worksize(scheduler.clamp_worksize(4, 2, 16))
    assert sentence == "4 samples at a time, 2 threads each, using 8 of your 16 CPU threads."
    assert scheduler.describe_worksize(scheduler.clamp_worksize(1, 1, 8)).startswith(
        "1 sample at a time, 1 thread each")
    assert "Free memory" in scheduler.describe_worksize(
        scheduler.auto_worksize(HardwareSnapshot(32, None)))
    assert "came down" in scheduler.describe_worksize(scheduler.clamp_worksize(8, 8, 16))


def test_the_constraint_between_the_two_numbers_is_stated_in_words():
    words = scheduler.worksize_constraint(16)
    assert "16 CPU threads" in words
    assert "never the samples" in words
    assert "1 to 16" in words


def test_the_two_numbers_become_the_admission_plan_run_bounded_already_enforces():
    plan = scheduler.plan_for(4, 2, memory_gb=1, hardware=hardware())
    assert (plan.max_parallel, plan.threads_per_sample) == (4, 2)
    assert plan.max_parallel * plan.threads_per_sample <= plan.cpu_budget
    plan.validate()
    size = scheduler.auto_worksize(hardware())
    assert scheduler.plan_for(size, memory_gb=1, hardware=hardware()).max_parallel == size.jobs
    with pytest.raises(ValueError, match="safely available"):
        scheduler.plan_for(4, 4, memory_gb=64, hardware=hardware(16, 8))


def test_a_run_that_states_no_resources_of_its_own_uses_the_chosen_numbers(no_default_worksize):
    machine = scheduler.detect_hardware()
    before = scheduler.resources_for_run({})
    scheduler.set_default_worksize(scheduler.clamp_worksize(2, 1, machine.cpus))
    chosen = scheduler.resources_for_run({})
    assert (chosen.max_parallel, chosen.threads_per_sample) == (min(2, machine.cpus), 1)
    # A run that states its own resources still wins over the saved default.
    stated = scheduler.resources_for_run({"threads": 1, "memory_gb": 1})
    assert stated.threads_per_sample == 1
    explicit = scheduler.resources_for_run({"resource_plan": before.to_dict()})
    assert explicit.to_dict() == before.to_dict()
    scheduler.set_default_worksize(None)
    assert scheduler.default_worksize() is None
    assert scheduler.resources_for_run({}).to_dict() == before.to_dict()


@pytest.mark.parametrize("cpus,available_gb", [(2, 7), (4, 6), (8, 5)])
def test_a_modest_computer_gets_a_smaller_job_rather_than_a_refusal(cpus, available_gb):
    """Refusing here stopped the pipeline with no result and no usable explanation.

    plan_resources used to raise whenever free memory was below the 8 GiB a sample
    asks for. The assembly task then failed, and because typing only begins after a
    successful assembly the run ended with result None. An ordinary 8 GiB laptop hit
    this, and so did CI, which is how it was found.
    """
    plan = plan_resources(hardware=hardware(cpus, available_gb))
    assert plan.max_parallel >= 1
    assert plan.reduced_for_memory is True
    assert scheduler.MINIMUM_SAMPLE_MEMORY_GB <= plan.memory_gb < 8
    # The reduced plan must still be internally consistent, not merely non-raising.
    assert plan.max_parallel * plan.memory_gb <= plan.memory_budget_gb
    assert plan.max_parallel * plan.threads_per_sample <= plan.cpu_budget


def test_a_roomy_computer_is_unchanged_by_the_small_machine_path():
    plan = plan_resources(hardware=hardware(16, 30))
    assert plan.memory_gb == 8 and plan.reduced_for_memory is False and plan.max_parallel > 1


def test_too_little_memory_says_so_in_words_the_user_can_act_on():
    with pytest.raises(ValueError) as failure:
        plan_resources(hardware=hardware(2, 1))
    message = str(failure.value)
    # It has to name memory as the cause and give the user something to do, because
    # the old failure surfaced as a run that simply produced no result.
    assert "memory" in message.lower() and "available" in message.lower()
    assert "free memory" in message.lower() or "lower the memory" in message.lower()


def test_a_request_far_beyond_the_machine_is_refused_rather_than_quietly_shrunk():
    """Shrinking a real requirement by tenfold would only fail later, less clearly."""
    with pytest.raises(ValueError, match="safely available"):
        plan_resources(memory_gb=64, hardware=hardware(16, 8))
    # Two GiB free is genuinely too little to assemble, so it is refused too.
    with pytest.raises(ValueError, match="safely available"):
        plan_resources(hardware=hardware(2, 3))
