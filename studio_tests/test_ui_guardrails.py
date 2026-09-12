"""Guard against data loss and misleading state in native application workflows."""

import json
from types import SimpleNamespace

import pytest

from wmlstudio.app import MainWindow
from wmlstudio.export import write_csv, write_html, write_json, write_tsv


@pytest.mark.parametrize("arguments, expected", [([], []), (["--window-size", "1080x720"], [(1080, 720)])])
def test_entrypoint_preserves_screen_aware_size_without_explicit_override(monkeypatch, arguments, expected):
    import wmlstudio.app as desktop

    sizes = []
    application = SimpleNamespace(
        setApplicationName=lambda value: None, setOrganizationName=lambda value: None,
        setStyle=lambda value: None, setStyleSheet=lambda value: None, exec=lambda: 0,
    )
    window = SimpleNamespace(resize=lambda width, height: sizes.append((width, height)), show=lambda: None)
    monkeypatch.setattr(desktop, "QApplication", SimpleNamespace(instance=lambda: application))
    monkeypatch.setattr(desktop, "MainWindow", lambda *args: window)
    assert desktop.main(arguments) == 0
    assert sizes == expected


@pytest.mark.parametrize("writer", [write_csv, write_html, write_json, write_tsv])
def test_export_cannot_overwrite_input(writer, tmp_path):
    sequence = tmp_path / "isolate.fasta"
    sequence.write_text(">c\nACGT\n")
    original = sequence.read_bytes()
    with pytest.raises(ValueError, match="different output"):
        writer([{"sample_name": "isolate", "input_path": str(sequence)}], sequence)
    assert sequence.read_bytes() == original


def test_project_lock_and_safe_outputs(qtbot, tmp_path):
    window = MainWindow(storage_root=tmp_path)
    qtbot.addWidget(window)
    path = window.project_path
    with pytest.raises(ValueError, match="already open"):
        MainWindow.lock_project(path)
    with pytest.raises(ValueError, match="different output"):
        window.check_output(path)
    assert window.switch_project(path)
    window.close()
    lock = MainWindow.lock_project(path)
    lock.unlock()


def test_native_pdf_and_atomic_project_copy(qtbot, tmp_path, monkeypatch):
    window = MainWindow(storage_root=tmp_path / "data")
    qtbot.addWidget(window)
    source = tmp_path / "sample.fasta"
    source.write_text(">c\nACGT\n")
    window.import_paths([source])
    pdf = tmp_path / "report.pdf"
    window.write_pdf_report(pdf)
    assert pdf.read_bytes().startswith(b"%PDF-")
    copied = tmp_path / "copy.wmlstudio"
    monkeypatch.setattr("wmlstudio.app.QFileDialog.getSaveFileName", lambda *a: (str(copied), ""))
    window.save_project_copy()
    assert window.switch_project(copied)
    assert len(window.project.samples()) == 1
    assert window.project.samples()[0]["input_path"] == str(source)
    window.close()


def test_hydra_native_summary_does_not_merge_local_calls(qtbot, tmp_path):
    window = MainWindow(storage_root=tmp_path)
    qtbot.addWidget(window)
    report = {
        "hydra_version": "1.4.0", "samples": [{"sample": "<isolate>", "species": {"name": "Test organism"},
        "mlst": {"sequence_type": "42"}, "summary": {"amr_genes": 1}, "hits": [{"gene": "<script>",
        "element_type": "AMR", "database": "synthetic", "identity_pct": 100, "coverage_pct": 100}]}],
        "parameters": {}, "databases": [], "import_provenance": {"sha256": "a" * 64},
    }
    window.project.set_setting("hydra_report", report)
    window.restore_hydra()
    assert window.hydra_table.rowCount() == 1
    assert window.hydra_table.item(0, 3).text() == "1"
    assert "<script>" in window.hydra_view.toPlainText()
    assert window.project.samples() == []
    assert json.loads(json.dumps(window.project.get_setting("hydra_report")))["samples"][0]["mlst"]["sequence_type"] == "42"
    window.close()
