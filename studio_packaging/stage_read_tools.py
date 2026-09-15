"""Explicit, verified fastp staging. Never download at analysis or startup time.

Upstream fastp publishes an official hash-pinned Linux x86-64 binary and no
Windows binary at all. Linux therefore stages the official artifact; Windows
stages only a separately built, reviewed native artifact supplied with
``--source``, exactly as the Windows SKA2 and SKESA payloads are handled. When
no verified Windows artifact exists the package simply ships without fastp and
the application says so, rather than substituting a different program.
"""

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
VERSION = "1.3.7"
SOURCE_COMMIT = "8a2397b6628ae14127efdb7566f67fc05f9aea56"
LINUX_URL = f"https://opengene.org/fastp/fastp.{VERSION}"
LINUX_SHA256 = "d54016cb118628f4d6645798bdf4e66b4e5263ce1c9838668ea3d595522fb495"
SOURCE_URL = f"https://codeload.github.com/OpenGene/fastp/tar.gz/{SOURCE_COMMIT}"
SOURCE_SHA256 = "7c758d46aee6044549175e3e1fb37dd1ae32b5d6aca7a49a9f521f6e72e108c6"
SOURCE_ARCHIVE = f"fastp-{SOURCE_COMMIT[:7]}-source.tar.gz"
#: Texts copied out of the pinned source archive so the licence travels with the binary.
SOURCE_TEXTS = {"LICENSE": "LICENSE", "README.md": "fastp-README.md"}


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def safe_relative(value):
    relative = PurePosixPath(value)
    if not relative.parts or relative.is_absolute() or ".." in relative.parts or "\\" in value or ":" in value:
        raise ValueError(f"Unsafe fastp package path: {value}")
    return Path(*relative.parts)


def verify(directory, platform=None):
    """Confirm a staged directory is the reviewed payload, file by file."""
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("tool") != "fastp" or manifest.get("version") != VERSION
            or manifest.get("source_commit") != SOURCE_COMMIT
            or platform and manifest.get("platform") != platform):
        raise ValueError("fastp source, version or platform mismatch")
    windows = manifest.get("platform") == "windows-x64"
    required = {"fastp.exe" if windows else "fastp", "LICENSE", SOURCE_ARCHIVE}
    if windows:
        required.add("smoke-test.json")
    seen = set()
    for entry in manifest.get("files", []):
        relative = safe_relative(entry["path"])
        folded = relative.as_posix().casefold()
        if folded in seen:
            raise ValueError("Duplicated fastp manifest file")
        seen.add(folded)
        path = directory / relative
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(directory):
            raise ValueError("fastp manifest file is missing or unsafe")
        if digest(path) != entry["sha256"]:
            raise ValueError(f"fastp payload checksum mismatch: {relative}")
    if not {name.casefold() for name in required}.issubset(seen):
        raise ValueError("fastp required binary, source or license files are missing")
    if windows and json.loads((directory / "smoke-test.json").read_text()).get("passed") is not True:
        raise ValueError("Native Windows fastp has no successful trimming smoke evidence")
    return manifest


def fetch(url, checksum, target):
    target = Path(target)
    if not target.is_file():
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".fastp-download-", delete=False) as handle:
            temporary = Path(handle.name)
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "WMLSTudio-build"})
                with urllib.request.urlopen(request, timeout=180) as response:
                    shutil.copyfileobj(response, handle, 1024 * 1024)
            except BaseException:
                handle.close()
                temporary.unlink(missing_ok=True)
                raise
        if digest(temporary) != checksum:
            temporary.unlink()
            raise ValueError(f"Official fastp download checksum mismatch: {url}")
        os.replace(temporary, target)
    if digest(target) != checksum:
        raise ValueError(f"Cached fastp archive checksum mismatch: {target}")
    return target


def _extract_texts(archive, target):
    """Take only the licence and readme the manifest promises, never the whole tree."""
    root = f"fastp-{SOURCE_COMMIT}"
    found = {}
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            relative = safe_relative(member.name)
            if not member.isfile() or relative.parts[0] != root or len(relative.parts) != 2:
                continue
            name = SOURCE_TEXTS.get(relative.parts[1])
            if name is None or name in found:
                continue
            with bundle.extractfile(member) as source, (target / name).open("wb") as handle:
                shutil.copyfileobj(source, handle)
            found[name] = True
    missing = sorted(set(SOURCE_TEXTS.values()) - set(found))
    if missing:
        raise ValueError(f"Pinned fastp source archive is missing required texts: {missing}")


def stage(destination, platform, source=None, cache=None):
    """Stage one platform's payload atomically; a valid existing snapshot is reused."""
    destination = Path(destination).resolve()
    if destination.exists():
        return verify(destination, platform)
    if platform == "windows-x64":
        if source is None:
            raise ValueError(
                "Upstream fastp publishes no official Windows binary. Supply a reviewed native "
                "Windows artifact with --source, or build without fastp: the application reports "
                "read trimming as unavailable rather than substituting another program.")
        source = Path(source).resolve()
        manifest = verify(source, platform)
    elif platform != "linux-x64":
        raise ValueError("Unsupported fastp platform")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".fastp-stage-") as temporary:
        target = Path(temporary) / "fastp"
        target.mkdir()
        if platform == "windows-x64":
            # Copy only manifest-listed files: never redistribute whatever else
            # happens to sit beside an artifact, and never follow a symlink.
            for entry in manifest["files"]:
                relative = safe_relative(entry["path"])
                (target / relative).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source / relative, target / relative)
        else:
            cache = Path(cache or ROOT / "artifacts/fastp-upstream")
            binary = fetch(LINUX_URL, LINUX_SHA256, cache / f"fastp.{VERSION}")
            archive = fetch(SOURCE_URL, SOURCE_SHA256, cache / SOURCE_ARCHIVE)
            shutil.copy2(binary, target / "fastp")
            (target / "fastp").chmod(0o755)
            shutil.copy2(archive, target / SOURCE_ARCHIVE)
            _extract_texts(archive, target)
            manifest = {
                "tool": "fastp", "version": VERSION, "platform": platform,
                "source_commit": SOURCE_COMMIT, "source_url": SOURCE_URL,
                "source_archive_sha256": SOURCE_SHA256,
                "binary_url": LINUX_URL, "binary_sha256": LINUX_SHA256,
                "license": "MIT; the original upstream LICENSE text is retained beside the binary",
                "build_provenance": ("Official upstream Linux x86-64 release binary from opengene.org, "
                                     "hash-pinned and unmodified; not a WMLSTudio build."),
                "scope": ("Adapter trimming and quality filtering of paired FASTQ reads. Trimming is "
                          "not an isolate validation, a purity check, or a species assignment."),
                "files": [{"path": path.relative_to(target).as_posix(), "bytes": path.stat().st_size,
                           "sha256": digest(path)} for path in sorted(target.rglob("*")) if path.is_file()],
            }
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        verify(target, platform)
        os.replace(target, destination)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("windows-x64", "linux-x64"), required=True)
    parser.add_argument("--source", type=Path, help="Reviewed native Windows fastp artifact directory")
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    destination = args.destination or ROOT / "src/wmlstudio/resources/tools/fastp" / args.platform
    result = stage(destination, args.platform, args.source, args.cache)
    total = sum(entry["bytes"] for entry in result["files"])
    print(f"Verified fastp {result['version']} {result['platform']}: "
          f"{len(result['files'])} files, {total:,} bytes in {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
