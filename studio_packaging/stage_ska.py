"""Explicit, verified SKA2 staging. Never download at analysis/runtime startup."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.5.1"
REVISION = "fcf9413d2768dc6538d31f664a4bf310651449e7"
LINUX_URL = "https://github.com/bacpop/ska.rust/releases/download/v0.5.1/ska-v0.5.1-ubuntu-latest-stable.tar.gz"
LINUX_SHA256 = "a48ff20fe5723669ef09dc2aa59ca19ef6a230d68b64c34d9e62b816224b5e1e"
SOURCE_URL = f"https://codeload.github.com/bacpop/ska.rust/tar.gz/{REVISION}"
SOURCE_SHA256 = "ab4de86d31ce44055a64cfccc088d577f252e9b8fff9ca86472e0f8a4a049ce1"


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def safe_relative(value):
    relative = PurePosixPath(value)
    if not relative.parts or relative.is_absolute() or ".." in relative.parts or "\\" in value or ":" in value:
        raise ValueError(f"Unsafe SKA2 package path: {value}")
    return Path(*relative.parts)


def verify(directory, platform=None):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("tool") != "SKA2" or manifest.get("version") != VERSION
            or manifest.get("source_commit") != REVISION or platform and manifest.get("platform") != platform):
        raise ValueError("SKA2 source, version or platform mismatch")
    windows = manifest.get("platform") == "windows-x64"
    required = {"ska.exe" if windows else "ska", "LICENSE", "NOTICE"}
    if windows:
        required |= {"Cargo.lock", "ska-0.5.1-corresponding-source.tar.gz", "smoke-test.json"}
    else:
        required.add("ska-pinned-upstream-source.tar.gz")
    seen = set()
    for entry in manifest.get("files", []):
        relative = safe_relative(entry["path"])
        folded = relative.as_posix().casefold()
        if folded in seen:
            raise ValueError("Duplicated SKA2 manifest file")
        seen.add(folded)
        path = directory / relative
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(directory):
            raise ValueError("SKA2 manifest file is missing or unsafe")
        if digest(path) != entry["sha256"]:
            raise ValueError(f"SKA2 payload checksum mismatch: {relative}")
    if not {name.casefold() for name in required}.issubset(seen):
        raise ValueError("SKA2 required binary/source/license files are missing")
    if windows and json.loads((directory / "smoke-test.json").read_text()).get("passed") is not True:
        raise ValueError("Native Windows SKA2 has no successful synthetic smoke evidence")
    return manifest


def fetch(url, checksum, target):
    target = Path(target)
    if not target.is_file():
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".ska-download-", delete=False) as handle:
            temporary = Path(handle.name)
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "WMLSTudio-build"}), timeout=120) as response:
                    shutil.copyfileobj(response, handle)
            except BaseException:
                handle.close()
                temporary.unlink(missing_ok=True)
                raise
        if digest(temporary) != checksum:
            temporary.unlink()
            raise ValueError("Official SKA2 download checksum mismatch")
        os.replace(temporary, target)
    if digest(target) != checksum:
        raise ValueError("Cached SKA2 archive checksum mismatch")
    return target


def stage(destination, platform, source=None, cache=None):
    destination = Path(destination).resolve()
    if destination.exists():
        return verify(destination, platform)
    if platform == "windows-x64":
        if source is None:
            raise ValueError("Windows SKA2 requires the verified native workflow artifact --source")
        source = Path(source).resolve()
        manifest = verify(source, platform)
    elif platform != "linux-x64":
        raise ValueError("Unsupported SKA2 platform")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".ska-stage-") as temporary:
        target = Path(temporary) / "ska2"
        target.mkdir()
        if platform == "windows-x64":
            # Copy only manifest-listed files. Do not distribute unrelated files
            # present beside an artifact or silently accept symlinked payloads.
            for entry in manifest["files"]:
                relative = safe_relative(entry["path"])
                (target / relative).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / relative, target / relative)
        else:
            cache = Path(cache or ROOT / "artifacts/ska2-upstream")
            archive = fetch(LINUX_URL, LINUX_SHA256, cache / "ska-v0.5.1-ubuntu-latest-stable.tar.gz")
            with tarfile.open(archive) as bundle:
                seen = set()
                for member in bundle.getmembers():
                    relative = safe_relative(member.name)
                    if not member.isfile() or relative.as_posix().casefold() in seen:
                        raise ValueError("Unexpected SKA2 Linux archive entry")
                    seen.add(relative.as_posix().casefold())
                    output = target / relative
                    output.parent.mkdir(parents=True, exist_ok=True)
                    with bundle.extractfile(member) as original, output.open("wb") as handle:
                        shutil.copyfileobj(original, handle)
                    output.chmod(member.mode & 0o777)
            upstream = fetch(SOURCE_URL, SOURCE_SHA256, cache / "ska-fcf9413-source.tar.gz")
            shutil.copy2(upstream, target / "ska-pinned-upstream-source.tar.gz")
            manifest = {"tool": "SKA2", "version": VERSION, "platform": platform,
                        "source_commit": REVISION, "source_url": SOURCE_URL,
                        "source_archive_sha256": SOURCE_SHA256,
                        "binary_archive_url": LINUX_URL, "binary_archive_sha256": LINUX_SHA256,
                        "license": "Apache-2.0; original upstream LICENSE/NOTICE retained",
                        "build_provenance": "Official upstream Linux release binary, not the separately compiled Windows toolchain",
                        "files": [{"path": path.relative_to(target).as_posix(), "bytes": path.stat().st_size,
                                   "sha256": digest(path)} for path in sorted(target.rglob("*")) if path.is_file()]}
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        verify(target, platform)
        os.replace(target, destination)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("windows-x64", "linux-x64"), required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--destination", type=Path)
    args = parser.parse_args()
    destination = args.destination or ROOT / "src/wmlstudio/resources/tools/ska2" / args.platform
    result = stage(destination, args.platform, args.source)
    print(f"Verified SKA2 {result['version']} {result['platform']} ({len(result['files'])} files)")
