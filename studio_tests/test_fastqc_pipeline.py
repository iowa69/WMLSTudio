import threading

import pytest
from PySide6.QtWidgets import QMessageBox

from wmlstudio.app import MainWindow
from wmlstudio.sequence import AnalysisCancelled


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    errors = []
    monkeypatch.setattr(MainWindow, "error", lambda self, text: errors.append(text))
    widget = MainWindow(storage_root=tmp_path / "workspace")
    qtbot.addWidget(widget)
    path = tmp_path / "reads.fastq"
    path.write_text("@one\nACGT\n+\nIIII\n")
    widget.import_paths([path])
    widget.show()
    monkeypatch.setattr(widget, "review_run_plan", lambda *args, **kwargs: {"fastqc": True})
    yield widget
    if widget.worker and widget.worker.isRunning():
        widget.cancel_analysis()
        qtbot.waitUntil(lambda: not widget.worker.isRunning(), timeout=10000)
        qtbot.wait(50)
    assert not errors
    widget.close()


def wait_idle(qtbot, window):
    qtbot.waitUntil(lambda: window.worker is not None and not window.worker.isRunning()
                    and window.run_button.isEnabled(), timeout=15000)


@pytest.mark.parametrize("flag,answer,typed", [
    ("WARN", None, True), ("FAIL", QMessageBox.StandardButton.Yes, True),
    ("FAIL", QMessageBox.StandardButton.No, False),
])
def test_fastqc_precedes_typing_and_failed_module_requires_review(window, qtbot, monkeypatch, flag, answer, typed):
    calls = []
    def fastqc(project, samples, output, allocation, **options):
        calls.append("fastqc")
        result = {"engine": "FastQC", "reports": [{"qc_status": flag}]}
        project.update_metadata(samples[0]["id"], {"fastqc": result})
        return [{"sample_id": samples[0]["id"], "result": result}]
    monkeypatch.setattr("wmlstudio.fastqc_dialog.run_project_fastqc", fastqc)
    if answer is not None:
        monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: answer)
    else:
        monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: pytest.fail("WARN must not be relabeled FAIL"))
    original = window.begin_typing
    def typing(samples, scheme):
        calls.append("typing")
        return original(samples, scheme)
    monkeypatch.setattr(window, "begin_typing", typing)
    window.start_analysis(confirm=True)
    wait_idle(qtbot, window)
    assert calls == (["fastqc", "typing"] if typed else ["fastqc"])
    sample = window.project.samples()[0]
    assert sample["metadata"]["fastqc"]["reports"][0]["qc_status"] == flag
    assert sample["status"] == ("completed" if typed else "queued")


def test_fastqc_cancellation_never_starts_assembly_or_typing(window, qtbot, monkeypatch):
    entered = threading.Event()
    def fastqc(project, samples, output, allocation, **options):
        entered.set()
        while not options["cancelled"]():
            threading.Event().wait(.005)
        raise AnalysisCancelled()
    monkeypatch.setattr("wmlstudio.fastqc_dialog.run_project_fastqc", fastqc)
    monkeypatch.setattr(window, "begin_typing", lambda *args: pytest.fail("Cancelled FastQC must not start typing"))
    monkeypatch.setattr(window, "assemble_pairs", lambda *args: pytest.fail("Cancelled FastQC must not start assembly"))
    window.start_analysis(confirm=True)
    qtbot.waitUntil(entered.is_set, timeout=5000)
    window.cancel_analysis()
    wait_idle(qtbot, window)
    assert window.project.samples()[0]["status"] == "queued"
    assert window._run_plan == {}
