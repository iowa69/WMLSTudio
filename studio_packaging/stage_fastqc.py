"""Explicit build-time staging of original FastQC and an app-local Temurin JRE.

No runtime download, system Java installation, or modifications to either
upstream program. Exact corresponding source archives and legal trees travel
with the portable package. Only FASTQ operation is exposed by WMLSTudio.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import posixpath
import shutil
import stat
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
FASTQC_VERSION = "0.12.1"
FASTQC_COMMIT = "e7ef390bf10382f60786bdd0cf28abd4f8683ffd"
JRE_VERSION = "17.0.20.1+1"
RELEASE = "https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.20.1%2B1/"
PACKAGES = {
    "fastqc": ("fastqc_v0.12.1.zip", "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip",
               "5f4dba8780231a25a6b8e11ab2c238601920c9704caa5458d9de559575d58aa7"),
    "fastqc_source": ("FastQC-e7ef390.tar.gz", f"https://codeload.github.com/s-andrews/FastQC/tar.gz/{FASTQC_COMMIT}",
                      "69c51ed13083eb5fa85c2d79245e3b16b70eeba7987e33f52fe6f522ce3293e2"),
    "java_source": ("OpenJDK17U-jdk-sources_17.0.20.1_1.tar.gz", RELEASE + "OpenJDK17U-jdk-sources_17.0.20.1_1.tar.gz",
                    "21e2a065d244ab048e737f21af5d1fc74daaeb6707de36477ead8db1dca71214"),
    "windows-x64": ("OpenJDK17U-jre_x64_windows_hotspot_17.0.20.1_1.zip", RELEASE + "OpenJDK17U-jre_x64_windows_hotspot_17.0.20.1_1.zip",
                    "bc21a93923103cdaac93ee337b0ae4365e739fde36df823dd456bc67c8a9d352"),
    "linux-x64": ("OpenJDK17U-jre_x64_linux_hotspot_17.0.20.1_1.tar.gz", RELEASE + "OpenJDK17U-jre_x64_linux_hotspot_17.0.20.1_1.tar.gz",
                  "0b2b640e3046b64c8ec504de0ab9d91bb5610182bda21fad454681ce54d45a62"),
}


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def acquire(key, cache, *, offline=False):
    name, url, checksum = PACKAGES[key]
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    destination = cache / name
    if not destination.is_file():
        if offline:
            raise ValueError(f"Required pinned archive is not cached: {destination}")
        with tempfile.NamedTemporaryFile(dir=cache, prefix=".download-", delete=False) as handle:
            temporary = Path(handle.name)
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "WMLSTudio-build/1"})
                with urllib.request.urlopen(request, timeout=120) as response:
                    shutil.copyfileobj(response, handle)
            except BaseException:
                handle.close()
                temporary.unlink(missing_ok=True)
                raise
        if digest(temporary) != checksum:
            temporary.unlink()
            raise ValueError(f"Pinned archive checksum mismatch: {name}")
        os.replace(temporary, destination)
    if digest(destination) != checksum:
        raise ValueError(f"Pinned archive checksum mismatch: {name}")
    return destination


def _relative(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError(f"Unsafe archive path: {name}")
    return Path(*path.parts[1:])


def extract(archive, destination):
    """Strip one root; materialize internal license symlinks as regular files."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    seen = set()

    def output(name):
        relative = _relative(name)
        if not relative.parts:
            return None
        key = relative.as_posix().casefold()
        if key in seen:
            raise ValueError(f"Duplicate archive path: {name}")
        seen.add(key)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as bundle:
            for entry in bundle.infolist():
                _relative(entry.orig_filename)
                _relative(entry.filename)
                if entry.is_dir():
                    continue
                if stat.S_ISLNK(entry.external_attr >> 16):
                    raise ValueError("ZIP symlinks are not supported")
                target = output(entry.filename)
                if target is not None:
                    with bundle.open(entry) as source, target.open("wb") as handle:
                        shutil.copyfileobj(source, handle)
    else:
        with tarfile.open(archive) as bundle:
            for entry in bundle.getmembers():
                _relative(entry.name)
                if entry.isdir():
                    continue
                if entry.issym() or entry.islnk():
                    target_name = (posixpath.join(posixpath.dirname(entry.name), entry.linkname)
                                   if entry.issym() else entry.linkname)
                    normalized = posixpath.normpath(target_name)
                    _relative(normalized)
                    if PurePosixPath(normalized).parts[0] != PurePosixPath(entry.name).parts[0]:
                        raise ValueError("Archive link escapes its root")
                elif not entry.isfile():
                    raise ValueError("Special archive entries are not supported")
                target = output(entry.name)
                if target is not None:
                    with bundle.extractfile(entry) as source, target.open("wb") as handle:
                        shutil.copyfileobj(source, handle)
                    target.chmod(entry.mode & 0o777)


def verify(destination, platform=None):
    destination = Path(destination)
    manifest = json.loads((destination / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("fastqc_version") != FASTQC_VERSION or manifest.get("java_version") != JRE_VERSION
            or platform and manifest.get("platform") != platform):
        raise ValueError("FastQC/JRE staging version or platform mismatch")
    if not manifest.get("files"):
        raise ValueError("Empty FastQC/JRE manifest")
    for entry in manifest["files"]:
        relative = PurePosixPath(entry["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe staged manifest path")
        path = destination / Path(*relative.parts)
        if not path.is_file() or digest(path) != entry["sha256"]:
            raise ValueError(f"FastQC/JRE staged checksum mismatch: {entry['path']}")
    return manifest


def stage(destination, platform, cache, *, offline=False):
    if platform not in ("windows-x64", "linux-x64"):
        raise ValueError("Unsupported FastQC target platform")
    destination = Path(destination).resolve()
    if destination.exists():
        return verify(destination, platform)
    archives = {key: acquire(key, cache, offline=offline)
                for key in ("fastqc", "fastqc_source", "java_source", platform)}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".fastqc-stage-") as temporary:
        target = Path(temporary) / "fastqc"
        extract(archives["fastqc"], target)
        extract(archives[platform], target / "jre")
        sources = target / "sources"
        sources.mkdir()
        for key in ("fastqc_source", "java_source"):
            shutil.copy2(archives[key], sources / archives[key].name)
        required = ["uk/ac/babraham/FastQC/FastQCApplication.class", "LICENSE", "jre/release",
                    "jre/NOTICE", "jre/legal/java.base/LICENSE", "htsjdk.jar", "jbzip2-0.9.jar",
                    "jre/bin/java.exe" if platform == "windows-x64" else "jre/bin/java"]
        if any(not (target / item).is_file() for item in required):
            raise ValueError("Pinned FastQC/JRE archive is incomplete")
        files = [{"path": path.relative_to(target).as_posix(), "size": path.stat().st_size,
                  "sha256": digest(path)} for path in sorted(target.rglob("*")) if path.is_file()]
        manifest = {"tool": "FastQC", "fastqc_version": FASTQC_VERSION, "fastqc_commit": FASTQC_COMMIT,
                    "java_provider": "Eclipse Temurin OpenJDK", "java_version": JRE_VERSION,
                    "platform": platform, "scope": "Original FastQC FASTQ checks; no automatic trimming",
                    "licenses": "FastQC GPL-3.0-or-later; OpenJDK GPL-2.0 with Classpath Exception; bundled component notices retained",
                    "archives": [{"name": PACKAGES[key][0], "url": PACKAGES[key][1],
                                  "sha256": PACKAGES[key][2]} for key in archives], "files": files}
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        os.replace(target, destination)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("windows-x64", "linux-x64"), required=True)
    parser.add_argument("--destination", type=Path, default=ROOT / "src/wmlstudio/resources/tools/fastqc")
    parser.add_argument("--cache", type=Path, default=ROOT / "artifacts/fastqc-downloads")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    result = stage(args.destination, args.platform, args.cache, offline=args.offline)
    print(f"Verified FastQC {result['fastqc_version']} / Temurin {result['java_version']}: {len(result['files'])} files")
