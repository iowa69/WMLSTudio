"""Reproduce native screenshots and PDF using the synthetic practice project.

Usage: QT_QPA_PLATFORM=offscreen uv run python studio_scripts/capture_desktop.py
Artifacts are generated in results/2026-09-12_validation; no real input is changed.
"""

import argparse
import tempfile
import time
from pathlib import Path

from PySide6.QtWidgets import QApplication

from wmlstudio.app import MainWindow
from wmlstudio.theme import STYLE


def capture(output: Path):
    app = QApplication.instance() or QApplication([])
    app.setApplicationName("WMLSTudio")
    app.setStyle("Fusion")
    app.setStyleSheet(STYLE)
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="wmlstudio-capture-") as folder:
        window = MainWindow(storage_root=Path(folder))
        window.resize(1380, 940)
        window.show()
        window.load_demo()
        deadline = time.monotonic() + 30
        while window.worker.isRunning() or not window.run_button.isEnabled():
            app.processEvents()
            if time.monotonic() > deadline:
                window.cancel_analysis()
                window.worker.wait(10000)
                raise RuntimeError("Practice analysis timed out")
        window.set_motion(False)
        for index, name in [(0, "overview"), (1, "samples"), (2, "comparison"), (3, "schemes"), (5, "reports")]:
            window.navigate(index)
            if index == 1:
                window.sample_table.selectRow(0)
            app.processEvents()
            window.grab().save(str(output / f"{name}.png"))
        window.write_pdf_report(output / "practice-report.pdf")
        window.resize(1080, 720)
        window.navigate(0)
        app.processEvents()
        window.grab().save(str(output / "compact-overview.png"))
        window.close()
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/2026-09-12_validation"))
    print(capture(parser.parse_args().output))
