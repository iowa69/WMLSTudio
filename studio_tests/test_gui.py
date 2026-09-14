"""Native workflow integration tests using real Qt queues and project files."""

import gzip
import hashlib
import json
import threading
import time

import pytest
from PySide6.QtGui import QImage

from wmlstudio import jobs
from wmlstudio.app import MainWindow
from wmlstudio.project import Project
from wmlstudio.scheduler import GIB, HardwareSnapshot, plan_resources
from wmlstudio.sequence import AnalysisCancelled


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, message: errors.append(str(message)))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    widget.test_errors = errors
    qtbot.addWidget(widget)
    widget.show()
    yield widget
    if widget.worker and widget.worker.isRunning():
        widget.worker.cancel()
        qtbot.waitUntil(lambda: not widget.worker.isRunning(), timeout=10000)
        qtbot.waitUntil(lambda: widget.run_button.isEnabled(), timeout=10000)
    widget.close()


def finish_queue(qtbot, window):
    qtbot.waitUntil(
        lambda: window.worker is not None and not window.worker.isRunning()
        and window.run_button.isEnabled(), timeout=15000,
    )


def fasta(tmp_path, name="isolate.fasta"):
    path = tmp_path / name
    path.write_text(">contig\nACGTACGTACGT\n", encoding="ascii")
    return path


def test_import_and_actual_read_qc_preserve_original_files(window, qtbot, tmp_path):
    assembly = fasta(tmp_path)
    reads = tmp_path / "reads.fastq.gz"
    reads.write_bytes(gzip.compress(b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\n!!!!\n"))
    unsupported = tmp_path / "notes.txt"
    unsupported.write_text("not a sequence")
    original = {path: path.read_bytes() for path in (assembly, reads, unsupported)}
    window.import_paths([assembly, reads, unsupported, tmp_path / "absent.fasta"])
    assert len(window.project.samples()) == 2
    assert window.sample_table.rowCount() == 2
    assert {sample["status"] for sample in window.project.samples()} == {"queued"}
    window.scheme_combo.setCurrentIndex(0)
    window.start_analysis()
    finish_queue(qtbot, window)
    assert window.test_errors == []
    samples = window.project.samples()
    assert {sample["status"] for sample in samples} == {"completed"}
    by_kind = {sample["result"]["kind"]: sample["result"] for sample in samples}
    assert by_kind["fastq"]["status"] == "qc_only"
    assert by_kind["fastq"]["qc"]["records"] == 2
    assert by_kind["fastq"]["qc"]["q30_percent"] == 50.0
    assert by_kind["fastq"]["alleles"] == {}
    assert by_kind["fasta"]["qc"]["total_bases"] == 12
    assert window.comparison_results() == []
    for path, data in original.items():
        assert path.read_bytes() == data
    assert by_kind["fastq"]["input_sha256"] == hashlib.sha256(original[reads]).hexdigest()


def test_demo_runs_real_typing_and_builds_conservative_interactive_forest(window, qtbot, tmp_path):
    window.load_demo()
    finish_queue(qtbot, window)
    assert window.test_errors == []
    samples = window.project.samples()
    assert len(samples) == 7
    assert all(sample["status"] == "completed" for sample in samples)
    results = {sample["name"]: sample["result"] for sample in samples}
    assert results["Practice_A01"]["st"] == "1"
    assert results["Practice_A03"]["st"] == "2"
    assert results["Practice_partial"]["status"] == "incomplete"
    assert len({result["scheme_digest"] for result in results.values()}) == 1
    window.cohort_ids = {sample['id'] for sample in samples}
    window.navigate(2)
    window.refresh_comparison()
    assert len(window.distance_rows) == 21
    assert sum(not row["comparable"] for row in window.distance_rows) == 6
    assert len(window.tree.nodes) == 7
    assert len(window.tree.edges) == 5
    pair = next(row for row in window.distance_rows if
                {row["source_name"], row["target_name"]} ==
                {"Practice_A01", "Practice_A02"})
    assert pair["distance"] == 0 and pair["shared_loci"] == 7
    source, _target, line, _label = window.tree.edges[0]
    old_position = line.line().p1()
    node = window.tree.nodes[source]
    node.setPos(node.pos().x() + 80, node.pos().y() + 35)
    assert line.line().p1() != old_position
    destination = tmp_path / "forest.png"
    window.tree.save_png(destination)
    image = QImage(str(destination))
    assert not image.isNull()
    assert image.width() == 1800 and image.height() == 1200


def test_practice_project_can_be_revisited_without_duplicate_samples(window, qtbot):
    window.load_demo()
    finish_queue(qtbot, window)
    original_ids = {sample["id"] for sample in window.project.samples()}
    window.load_demo()
    finish_queue(qtbot, window)
    assert {sample["id"] for sample in window.project.samples()} == original_ids
    assert len(window.project.samples()) == 7


def test_project_switch_preserves_results_and_failed_open_keeps_current_project(window, qtbot, tmp_path):
    assembly = fasta(tmp_path)
    window.import_paths([assembly])
    original_path = window.project_path
    sample_id = window.project.samples()[0]["id"]
    window.project.set_metadata(sample_id, {"location": "Dublin"})
    window.start_analysis()
    finish_queue(qtbot, window)
    assert window.switch_project(tmp_path / "second.wmlstudio") is True
    assert window.project.samples() == []
    assert window.sample_table.rowCount() == 0
    assert window.switch_project(original_path) is True
    restored = window.project.get_sample(sample_id)
    assert restored["metadata"] == {"location": "Dublin"}
    assert restored["result"]["qc"]["total_bases"] == 12
    invalid = tmp_path / "not-a-project.wmlstudio"
    invalid.write_bytes(b"not sqlite")
    assert window.switch_project(invalid) is False
    assert window.project_path == original_path
    assert window.project.get_sample(sample_id)["status"] == "completed"
    assert len(window.test_errors) == 1


def test_filtered_view_exports_all_sample_records(window, tmp_path, monkeypatch):
    first = fasta(tmp_path, "first.fasta")
    second = fasta(tmp_path, "second.fasta")
    window.import_paths([first, second])
    first_id = window.project.samples()[0]["id"]
    window.project.set_result(first_id, {
        "sample_name": "first.fasta", "status": "complete", "scheme": "test",
        "scheme_digest": "abc", "st": "7", "alleles": {"locus": "1"},
    })
    window.refresh()
    window.search.setText("first")
    assert window.sample_table.rowCount() == 1
    destination = tmp_path / "all-results.json"
    monkeypatch.setattr(
        "wmlstudio.app.QFileDialog.getSaveFileName", lambda *args: (str(destination), "JSON"),
    )
    window.export_project("json")
    exported = json.loads(destination.read_text())["samples"]
    assert len(exported) == 2
    assert {row["job_status"] for row in exported} == {"queued", "completed"}
    assert "all 2 project samples" in window.progress_text.text()
    assert window.test_errors == []


def test_failed_reanalysis_does_not_enter_comparison(window, tmp_path):
    window.import_paths([fasta(tmp_path)])
    sample_id = window.project.samples()[0]["id"]
    window.project.set_result(sample_id, {
        "sample_name": "isolate", "status": "complete", "scheme": "test",
        "scheme_digest": "abc", "alleles": {"locus": "1"},
    })
    window.cohort_ids = {sample_id}
    assert len(window.comparison_results()) == 1
    window.sample_failed(sample_id, "Input changed")
    assert window.comparison_results() == []
    window.refresh_comparison()
    assert window.distance_rows == []
    assert window.tree.nodes == {}


def test_invalid_fastq_fails_one_sample_without_losing_next_result(window, qtbot, tmp_path):
    invalid = tmp_path / "truncated.fastq"
    invalid.write_text("@read\nACGT\n+\nII\n")
    assembly = fasta(tmp_path)
    window.import_paths([invalid, assembly])
    window.start_analysis()
    finish_queue(qtbot, window)
    rows = {sample["name"]: sample for sample in window.project.samples()}
    assert rows["truncated"]["status"] == "failed"
    assert "truncated FASTQ" in rows["truncated"]["error"]
    assert rows["isolate"]["status"] == "completed"
    assert window.run_button.isEnabled()
    assert window.test_errors == []


@pytest.mark.parametrize("close_window", [False, True])
def test_cancellation_keeps_unstarted_samples_and_defers_close(
    window, qtbot, tmp_path, monkeypatch, close_window,
):
    entered = threading.Event()

    def cancellable_inspection(path, max_reads, cancelled):
        entered.set()
        deadline = time.monotonic() + 10
        while not cancelled():
            if time.monotonic() > deadline:
                raise RuntimeError("Test cancellation timed out")
            threading.Event().wait(0.005)
        raise AnalysisCancelled("Cancelled safely")

    monkeypatch.setattr(jobs, "inspect_sequence", cancellable_inspection)
    allocation = plan_resources(threads_per_sample=1, memory_gb=1, max_parallel=1,
                                hardware=HardwareSnapshot(4, 8 * GIB))
    monkeypatch.setattr("wmlstudio.scheduler.resources_for_run", lambda *args, **kwargs: allocation)
    window.import_paths([fasta(tmp_path, "one.fasta"), fasta(tmp_path, "two.fasta")])
    path = window.project_path
    window.start_analysis()
    qtbot.waitUntil(entered.is_set, timeout=5000)
    qtbot.waitUntil(lambda: window.project.samples()[0]["status"] == "running", timeout=5000)
    assert window.switch_project(tmp_path / "not-opened.wmlstudio") is False
    assert not (tmp_path / "not-opened.wmlstudio").exists()
    if close_window:
        window.close()
        assert window.closing_after_cancel is True
        finish_queue(qtbot, window)
        qtbot.waitUntil(lambda: not window.isVisible(), timeout=5000)
        with Project(path) as project:
            samples = project.samples()
    else:
        window.cancel_analysis()
        finish_queue(qtbot, window)
        samples = window.project.samples()
        assert window.isVisible()
    assert [sample["status"] for sample in samples] == ["interrupted", "queued"]
    assert samples[0]["result"] is None
    assert window.test_errors == []


def test_parallel_cancellation_interrupts_only_admitted_samples(window, qtbot, tmp_path, monkeypatch):
    entered = set()
    lock = threading.Lock()

    def inspection(path, max_reads, cancelled):
        with lock:
            entered.add(str(path))
        deadline = time.monotonic() + 10
        while not cancelled():
            if time.monotonic() > deadline:
                raise RuntimeError("Test cancellation timed out")
            threading.Event().wait(0.005)
        raise AnalysisCancelled("Cancelled safely")

    allocation = plan_resources(threads_per_sample=1, memory_gb=1, max_parallel=2,
                                hardware=HardwareSnapshot(4, 8 * GIB))
    monkeypatch.setattr("wmlstudio.scheduler.resources_for_run", lambda *args, **kwargs: allocation)
    monkeypatch.setattr(jobs, "inspect_sequence", inspection)
    window.import_paths([fasta(tmp_path, name) for name in ("one.fasta", "two.fasta", "three.fasta")])
    window.start_analysis()
    qtbot.waitUntil(lambda: len(entered) == 2, timeout=5000)
    qtbot.waitUntil(lambda: sum(s["status"] == "running" for s in window.project.samples()) == 2,
                    timeout=5000)
    window.cancel_analysis()
    finish_queue(qtbot, window)
    assert len(entered) == 2
    assert [sample["status"] for sample in window.project.samples()] == ["interrupted", "interrupted", "queued"]
    assert all(sample["result"] is None for sample in window.project.samples())
    assert window.test_errors == []


def test_explicit_analysis_cohort_does_not_inherit_global_selection(window, qtbot, tmp_path):
    window.import_paths([fasta(tmp_path, "first.fasta"), fasta(tmp_path, "second.fasta")])
    first, second = window.project.samples()
    window.selection_ids = {first["id"]}
    window.start_analysis(sample_ids={second["id"]})
    finish_queue(qtbot, window)
    assert window.project.get_sample(first["id"])["status"] == "queued"
    assert window.project.get_sample(second["id"])["status"] == "completed"
    assert window.selection_ids == {first["id"]}


def test_hydra_import_is_persisted_and_does_not_leak_across_projects(window, tmp_path, monkeypatch):
    upstream = {
        "hydra_version": "1.4.0", "command": "hydra run --format json", "databases": [],
        "parameters": {}, "samples": [{"sample": "HYDRA_ONLY_<script>", "hits": []}],
    }
    source = tmp_path / "hydra.json"
    source.write_text(json.dumps(upstream))
    monkeypatch.setattr(
        "wmlstudio.app.QFileDialog.getOpenFileName", lambda *args: (str(source), "JSON"),
    )
    original_project = window.project_path
    window.import_hydra()
    saved = window.project.get_setting("hydra_report")
    assert saved["samples"][0]["sample"] == "HYDRA_ONLY_<script>"
    assert saved["import_provenance"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert "HYDRA_ONLY_<script>" in window.hydra_view.toPlainText()
    assert "<script>" not in window.hydra_view.toHtml()
    assert window.project.samples() == []
    assert window.switch_project(tmp_path / "clean.wmlstudio") is True
    assert window.project.get_setting("hydra_report") is None
    assert "HYDRA_ONLY_" not in window.hydra_view.toPlainText()
    source.unlink()
    assert window.switch_project(original_project) is True
    assert window.project.get_setting("hydra_report") == saved
    assert "HYDRA_ONLY_<script>" in window.hydra_view.toPlainText()
    assert window.test_errors == []


def test_motion_preference_is_restored_when_switching_projects(window, tmp_path):
    original = window.project_path
    window.motion.setChecked(False)
    assert not window.helix.timer.isActive()
    assert window.switch_project(tmp_path / "default-motion.wmlstudio") is True
    assert window.motion.isChecked() is True
    assert window.motion_enabled is True
    assert window.helix.timer.isActive()
    assert window.switch_project(original) is True
    assert window.motion.isChecked() is False
    assert window.motion_enabled is False
    assert not window.helix.timer.isActive()
