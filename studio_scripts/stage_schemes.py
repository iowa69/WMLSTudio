"""Stage an attributed, checksummed scheme snapshot; never invoked at app startup."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
PINNED_REVISION = "6cad46ffd9dfddfaa55f7993cf391f80a70556d3"
SOURCE_REPOSITORY = "https://github.com/iowa69/WMLST"
ARCHIVE_URL = f"https://codeload.github.com/iowa69/WMLST/tar.gz/{PINNED_REVISION}"
ALLOWED_SUFFIXES = {".tfa", ".fasta", ".fa", ".fna", ".txt", ".json"}
MAX_ARCHIVE_BYTES = 512 * 1024 * 1024


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def unpack_schemes(archive: Path, destination: Path) -> Path:
    """Extract only regular scheme data, with path and decompression bounds."""
    destination.mkdir(parents=True, exist_ok=True)
    total = 0
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            parts = PurePosixPath(member.name).parts
            if len(parts) < 5 or parts[1:3] != ("db", "pubmlst"):
                continue
            relative = PurePosixPath(*parts[3:])
            if ".." in relative.parts or "\\" in str(relative) or ":" in str(relative):
                raise ValueError(f"Unsafe archive path: {member.name}")
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"Non-regular archive entry: {member.name}")
            if len(relative.parts) != 2 or relative.suffix.lower() not in ALLOWED_SUFFIXES:
                continue
            total += member.size
            if total > MAX_ARCHIVE_BYTES:
                raise ValueError("Scheme archive exceeds the unpacked size limit")
            output = destination.joinpath(*relative.parts)
            output.parent.mkdir(parents=True, exist_ok=True)
            source = bundle.extractfile(member)
            if source is None:
                raise ValueError(f"Missing archive content: {member.name}")
            with source, output.open("xb") as target:
                shutil.copyfileobj(source, target)
    return destination


def local_provenance(source: Path) -> dict:
    result = {"method": "local", "path": str(source.resolve())}
    try:
        result["repository_revision"] = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        result["working_tree_data_modified"] = bool(subprocess.check_output(
            ["git", "-C", str(source), "status", "--porcelain", "--", "."],
            text=True, stderr=subprocess.DEVNULL,
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        result["repository_revision"] = None
    return result


def stage_schemes(source: Path, destination: Path, provenance: dict | None = None) -> dict:
    source, destination = source.resolve(), destination.resolve()
    if not source.is_dir():
        raise ValueError(f"Scheme source directory does not exist: {source}")
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError("Source and destination must be separate, non-nested directories")
    files = []
    schemes = []
    for directory in sorted(source.iterdir()):
        if directory.is_symlink():
            raise ValueError(f"Scheme source must not contain symlinks: {directory}")
        if not directory.is_dir():
            continue
        data_files = [item for item in sorted(directory.iterdir())
                      if item.suffix.lower() in ALLOWED_SUFFIXES]
        if not any(item.suffix.lower() in {".tfa", ".fasta", ".fa", ".fna"}
                   for item in data_files):
            continue
        metadata = {}
        for item in data_files:
            if item.is_symlink() or not item.is_file():
                raise ValueError(f"Scheme data must be regular files: {item}")
            if item.name.endswith("_info.json"):
                metadata = json.loads(item.read_text(encoding="utf-8"))
            files.append({"path": item.relative_to(source).as_posix(),
                          "bytes": item.stat().st_size, "sha256": sha256(item)})
        schemes.append({"id": directory.name, "source_metadata": metadata})
    if not schemes:
        raise ValueError("No scheme allele FASTA files were found")
    allowed = {entry["path"] for entry in files} | {"manifest.json"}
    if destination.exists():
        for existing in destination.rglob("*"):
            if existing.is_symlink():
                raise ValueError(f"Staging destination contains a symlink: {existing}")
            if existing.is_file() and existing.relative_to(destination).as_posix() not in allowed:
                raise ValueError(f"Staging would leave stale data: {existing}; use a fresh directory")
    destination.mkdir(parents=True, exist_ok=True)
    for entry in files:
        output = destination / entry["path"]
        output.parent.mkdir(parents=True, exist_ok=True)
        if not output.exists() or sha256(output) != entry["sha256"]:
            shutil.copy2(source / entry["path"], output)
        if sha256(output) != entry["sha256"]:
            raise ValueError(f"Copy verification failed: {output}")
    digest = hashlib.sha256(json.dumps(files, sort_keys=True,
                                      separators=(",", ":")).encode()).hexdigest()
    manifest = {
        "format_version": 1,
        "created_utc": datetime.now(UTC).isoformat(),
        "source_repository": SOURCE_REPOSITORY,
        "provenance": provenance if provenance is not None else local_provenance(source),
        "snapshot_sha256": digest,
        "scheme_count": len(schemes),
        "file_count": len(files),
        "total_bytes": sum(entry["bytes"] for entry in files),
        "data_provider": "PubMLST, University of Oxford; individual scheme curators",
        "terms_url": "https://pubmlst.org/terms-conditions",
        "redistribution_review": "Submission dates not independently verified; see bundled notices",
        "schemes": schemes,
        "files": files,
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--source", type=Path, help="Local directory of scheme directories")
    inputs.add_argument("--download", action="store_true", help="Fetch the pinned GitHub snapshot")
    parser.add_argument("--destination", type=Path,
                        default=ROOT / "src/wmlstudio/resources/schemes")
    args = parser.parse_args(argv)
    if args.source:
        manifest = stage_schemes(args.source, args.destination)
    else:
        with tempfile.TemporaryDirectory(prefix="wmlstudio-schemes-") as temporary:
            archive = Path(temporary) / "snapshot.tar.gz"
            request = urllib.request.Request(ARCHIVE_URL, headers={"User-Agent": "WMLSTudio-build"})
            total = 0
            with urllib.request.urlopen(request, timeout=120) as response, archive.open("wb") as out:
                while chunk := response.read(1024 * 1024):
                    total += len(chunk)
                    if total > MAX_ARCHIVE_BYTES:
                        raise ValueError("Scheme archive exceeds the download size limit")
                    out.write(chunk)
            source = unpack_schemes(archive, Path(temporary) / "schemes")
            manifest = stage_schemes(source, args.destination, {
                "method": "pinned_github_archive", "repository_revision": PINNED_REVISION,
                "archive_url": ARCHIVE_URL, "archive_sha256": sha256(archive),
            })
    print(json.dumps({key: manifest[key] for key in (
        "scheme_count", "file_count", "total_bytes", "snapshot_sha256")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
