"""Explicit build step: preserve version-matched Python/Qt license texts."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import platform
import re
import sys
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import distribution, distributions, version
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "WMLSTudio-build-notices"})
    parsed = urlsplit(url)
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if parsed.scheme == "https" and parsed.hostname == "api.github.com" and parsed.port in (None, 443) and token:
        # Unredirected headers are deliberately not forwarded by urllib to a
        # redirected URL, including a raw-content/provider host.
        request.add_unredirected_header("Authorization", "Bearer " + token)
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read(2 * 1024 * 1024)


def preserve_installed_notices(destination, manifest):
    """Retain transitive-helper attributions without guessing the freeze graph."""
    for package in distributions():
        name = str(package.metadata.get("Name", "unknown")).lower().replace("_", "-")
        manifest["packages"][name] = package.version
        for relative in package.files or []:
            if not any(word in str(relative).lower() for word in ("license", "copying", "notice")):
                continue
            source = Path(package.locate_file(relative))
            if not source.is_file() or ".." in relative.parts:
                continue
            target = destination / "packages" / name / str(relative)
            target.parent.mkdir(parents=True, exist_ok=True)
            content = source.read_bytes()
            target.write_bytes(content)
            manifest["files"].append({"path": target.relative_to(destination).as_posix(),
                                      "source": f"installed wheel: {name}=={package.version}",
                                      "sha256": hashlib.sha256(content).hexdigest()})


def preserve_windows_python_notices(destination, manifest):
    """The Windows interpreter's combined notices exceed CPython's source LICENSE."""
    candidates = (Path(sys.base_prefix) / "LICENSE.txt", Path(sys.executable).resolve().parent / "LICENSE.txt")
    source = next((path for path in candidates if path.is_file()), None)
    if source is None:
        raise ValueError("The target Windows Python interpreter is missing its bundled LICENSE.txt")
    content = source.read_bytes()
    target = destination / "Python-Windows-runtime-LICENSE.txt"
    target.write_bytes(content)
    manifest["files"].append({"path": target.name,
        "source": f"Target Windows Python {platform.python_version()} interpreter LICENSE.txt (bundled runtime notices)",
        "sha256": hashlib.sha256(content).hexdigest()})


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
    manifest = {"python_version": python_version, "qt_version": qt_version, "files": [], "packages": {}}

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
    # Include every isolated build distribution: pandas/setuptools may pull
    # incidental helper modules into PYZ through imports not listed by the app.
    preserve_installed_notices(destination, manifest)
    if platform.system() == "Windows":
        preserve_windows_python_notices(destination, manifest)
    # Pyrodigal requires corresponding source; Pyskani's source archive also
    # preserves embedded skani and vendored Rust dependency license texts.
    for source_package in ("pyrodigal", "pyskani"):
        if source_package not in manifest["packages"]:
            continue
        package_version = manifest["packages"][source_package]
        listing = json.loads(fetch(f"https://pypi.org/pypi/{source_package}/{package_version}/json"))
        source_entry = next(item for item in listing["urls"] if item["packagetype"] == "sdist")
        target = destination / "sources" / source_entry["filename"]
        target.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(source_entry["url"], headers={"User-Agent": "WMLSTudio-build-notices"})
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read(64 * 1024 * 1024)
        actual = hashlib.sha256(content).hexdigest()
        if actual != source_entry["digests"]["sha256"]:
            raise ValueError(f"{source_package} source archive hash does not match the provider record")
        target.write_bytes(content)
        manifest["files"].append({"path": target.relative_to(destination).as_posix(),
                                  "source": source_entry["url"], "sha256": actual})
    hydra_commit = "6d36c109491c16544e8919fe6962b4b62e97d3d7"
    hydra_url = f"https://github.com/iowa69/hydra/archive/{hydra_commit}.zip"
    with urllib.request.urlopen(hydra_url, timeout=60) as response:
        hydra_source = response.read(16 * 1024 * 1024)
    with zipfile.ZipFile(io.BytesIO(hydra_source)) as archive:
        if archive.testzip() is not None:
            raise ValueError("HYDRA corresponding-source archive failed integrity verification")
    hydra_target = destination / "sources" / f"hydra-{hydra_commit}.zip"
    hydra_target.parent.mkdir(parents=True, exist_ok=True)
    hydra_target.write_bytes(hydra_source)
    manifest["files"].append({"path": hydra_target.relative_to(destination).as_posix(),
                              "source": hydra_url, "source_commit": hydra_commit,
                              "sha256": hashlib.sha256(hydra_source).hexdigest()})
    # Include the exact application/build sources accompanying this combined
    # binary, including local uncommitted changes; never include sample data.
    source_files = {ROOT / name for name in ("README.md", "LICENSE", "pyproject.toml", "uv.lock")}
    for folder in ("src/wmlstudio", "studio_tests", "studio_scripts", "studio_packaging", "docs", ".github/workflows"):
        for path in (ROOT / folder).rglob("*"):
            relative = path.relative_to(ROOT)
            if (path.is_file() and not any(part in {"generated_notices", "__pycache__"} for part in relative.parts)
                    and path.suffix in {".py", ".md", ".svg", ".spec", ".ps1", ".sh", ".patch", ".yml", ".txt"}
                    and not ("resources" in relative.parts and "ui" not in relative.parts)):
                source_files.add(path)
    source_target = destination / "sources" / "wmlstudio-application-source.zip"
    with zipfile.ZipFile(source_target, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(source_files):
            if path.is_file():
                archive.write(path, arcname="WMLSTudio/" + path.relative_to(ROOT).as_posix())
    manifest["files"].append({"path": source_target.relative_to(destination).as_posix(),
                              "source": "Exact local application/build source at staging time; sample data excluded",
                              "sha256": hashlib.sha256(source_target.read_bytes()).hexdigest()})
    if platform.system() == "Windows":
        qt_runtime = Path(distribution("PySide6-Essentials").locate_file("PySide6"))
        manifest["native_runtime"] = []
        for name in ("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"):
            path = qt_runtime / name
            content = path.read_bytes()
            if content[:2] != b"MZ":
                raise ValueError(f"The official Qt wheel has no valid native runtime: {name}")
            manifest["native_runtime"].append({"name": name, "sha256": hashlib.sha256(content).hexdigest(),
                "provider": f"Official PySide6-Essentials {qt_version} Windows wheel",
                "bundle_path": f"Tools/blast/bin/{name}"})
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Preserved {len(manifest['files'])} version-matched license texts in {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
