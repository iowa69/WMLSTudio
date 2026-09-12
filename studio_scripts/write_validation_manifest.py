"""Capture the exact local source and environment behind a development handoff."""

import hashlib
import json
import platform
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    files = [ROOT / "pyproject.toml", ROOT / "uv.lock"]
    for folder in ("src/wmlstudio", "studio_tests", "studio_scripts", "studio_packaging"):
        files.extend(sorted((ROOT / folder).glob("*.py")))
    files.extend(sorted((ROOT / "studio_packaging").glob("*.spec")))
    entries = [{"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
               for path in files]
    manifest = {
        "created_utc": datetime.now(UTC).isoformat(), "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": {name: version(name) for name in ("wmlstudio", "PySide6-Essentials", "pyahocorasick", "pytest", "pyinstaller")},
        "source_sha256": hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest(),
        "source_files": entries,
        "interpretation": "Records local source and environment. Windows/Wine execution is recorded separately; no clean Windows 11 acceptance is implied.",
    }
    output = ROOT / "results/2026-09-12_validation/source-manifest.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(output)


if __name__ == "__main__":
    main()
