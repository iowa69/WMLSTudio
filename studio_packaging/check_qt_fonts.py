"""CI-only font visibility gate; never redistribute Windows system fonts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PySide6.QtGui import QFont, QFontDatabase, QRawFont
from PySide6.QtWidgets import QApplication


def check_fonts(family):
    app = QApplication.instance() or QApplication([])
    families = QFontDatabase.families()
    if family not in families:
        raise RuntimeError(f"Headless Qt did not discover required system font: {family}")
    font = QRawFont.fromFont(QFont(family, 14))
    if not font.isValid() or not all(font.supportsCharacter(ord(char)) for char in "WMLSTudio genome ABCabc123"):
        raise RuntimeError(f"Headless Qt cannot render ordinary UI characters in {family}")
    return {"platform": app.platformName(), "requested_family": family,
            "resolved_family": font.familyName(), "family_count": len(families), "ascii_glyphs": "passed"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", default="Segoe UI")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = check_fonts(args.family)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
