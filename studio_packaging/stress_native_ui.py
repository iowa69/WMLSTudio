"""Exercise the real seven-page native workspace and capture repaint/size evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import time
from pathlib import Path

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication, QScrollArea

from wmlstudio.app import MainWindow
from wmlstudio.theme import STYLE, apply_dark_palette


def settle(app, milliseconds=120):
    loop = QEventLoop()
    QTimer.singleShot(milliseconds, loop.quit)
    loop.exec()
    app.processEvents()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    apply_dark_palette(app)
    app.setStyleSheet(STYLE)
    results = {"platform": app.platformName(), "captures": [], "navigation_count": 0}
    with tempfile.TemporaryDirectory(prefix="wmlstudio-native-audit-") as work:
        window = MainWindow(storage_root=Path(work))
        errors = []
        window.error = lambda message: errors.append(str(message))
        window.show()
        window.load_demo()
        deadline = time.monotonic() + 90
        while window.worker and window.worker.isRunning():
            if time.monotonic() > deadline:
                window.cancel_analysis()
                raise RuntimeError("Practice analysis did not finish during the UI audit")
            settle(app, 50)
        settle(app, 200)
        window.set_motion(False)
        for width, height in [(1380, 940), (1080, 720)]:
            window.resize(width, height)
            settle(app)
            for cycle in range(3):
                for index in [0, 2, 1, 4, 3, 6, 5]:
                    window.nav_buttons[index].click()
                    settle(app)
                    results["navigation_count"] += 1
                    if window.pages.currentIndex() != index:
                        raise AssertionError("Navigation did not select the requested page")
                    if cycle != 2:
                        continue
                    name = window.nav_names[index].lower().replace(" ", "-")
                    destination = args.output / f"{width}x{height}-{index}-{name}.png"
                    screenshot = None
                    method = "widget-render"
                    if app.platformName() == "windows":
                        screenshot = window.screen().grabWindow(int(window.winId()))
                        method = "native-window-surface"
                    if screenshot is None or screenshot.isNull():
                        screenshot = window.grab()
                        method = "widget-render"
                    if not screenshot.save(str(destination)):
                        raise OSError(f"Could not capture {destination}")
                    page = window.pages.currentWidget()
                    details = {"page": index, "name": name, "requested_size": [width, height],
                               "actual_size": [window.width(), window.height()], "capture": str(destination),
                               "capture_method": method, "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
                               "graphics_effect": str(window.pages.graphicsEffect())}
                    if isinstance(page, QScrollArea):
                        details.update(viewport_width=page.viewport().width(), content_width=page.widget().width(),
                                       minimum_content_width=page.widget().minimumSizeHint().width(),
                                       horizontal_overflow=page.horizontalScrollBar().maximum(),
                                       vertical_scroll=page.verticalScrollBar().maximum())
                    results["captures"].append(details)
                    print(json.dumps(details), flush=True)
        results["errors"] = errors
        results["sample_count"] = len(window.project.samples())
        results["graph_nodes"] = len(window.tree.nodes)
        window.close()
        settle(app, 50)
    (args.output / "audit.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return 1 if results["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
