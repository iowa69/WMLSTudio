"""Explicit build step: preserve version-matched Python/Qt license texts."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "WMLSTudio-build-notices"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read(2 * 1024 * 1024)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=ROOT / "studio_packaging/generated_notices")
    parser.add_argument("--python-version", default=platform.python_version(),
                        help="Target interpreter version when staging for a separate build runtime")
    args = parser.parse_args()
    qt_version = version("PySide6-Essentials")
    python_version = args.python_version
    if not all(re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value)
               for value in (qt_version, python_version)):
        raise ValueError("Stable Python and PySide6 releases are required for versioned notices")
    requests = [("Python-LICENSE.txt",
                 f"https://raw.githubusercontent.com/python/cpython/v{python_version}/LICENSE")]
    for repository in ("pyside/pyside-setup", "qt/qtbase"):
        listing = json.loads(fetch(
            f"https://api.github.com/repos/{repository}/contents/LICENSES?ref=v{qt_version}"))
        for entry in listing:
            if entry["type"] == "file":
                filename = Path(entry["name"])
                if filename.name != entry["name"]:
                    raise ValueError("Invalid license filename")
                requests.append((f"{repository.split('/')[1]}/{filename.name}",
                                 entry["download_url"]))
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest = {"python_version": python_version, "qt_version": qt_version, "files": []}

    def retrieve(item):
        relative, url = item
        return relative, url, fetch(url)

    with ThreadPoolExecutor(max_workers=8) as pool:
        for relative, url, content in pool.map(retrieve, requests):
            path = destination / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            manifest["files"].append({"path": relative, "source": url,
                                      "sha256": hashlib.sha256(content).hexdigest()})
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Preserved {len(manifest['files'])} version-matched license texts in {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
