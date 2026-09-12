"""Resolve resources independently of working directory and installed Python."""

import sys
from pathlib import Path

from PySide6.QtCore import QStandardPaths


def resource_root() -> Path:
    return Path(__file__).resolve().parent / "resources"


def data_root() -> Path:
    if getattr(sys, "frozen", False):
        root = Path(sys.executable).resolve().parent / "Data"
    else:
        root = Path(QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation))
    root.mkdir(parents=True, exist_ok=True)
    return root


def scheme_locations(root: Path) -> list[Path]:
    locations = [resource_root() / "schemes", root / "schemes"]
    return sorted({p for base in locations if base.is_dir() for p in base.iterdir()
                   if p.is_dir()}, key=lambda p: p.name.casefold())
