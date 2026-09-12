"""Explicit build step: stage SHA-256-pinned official native NCBI BLAST+ tools."""

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
VERSION = "2.17.0"
PACKAGES = {
    "windows-x64": {"archive": "ncbi-blast-2.17.0+-x64-win64.tar.gz",
                    "sha256": "ccde8788641e8f4137536aaadedfeac2f3599dbbc6166e701b5d89d19fa79038",
                    "suffix": ".exe"},
    "linux-x64": {"archive": "ncbi-blast-2.17.0+-x64-linux.tar.gz",
                  "sha256": "3888112d8207831aa47371d93583c601f058f88b5db22dc782438b039a3a411b",
                  "suffix": ""},
}
PROGRAMS = ("blastn", "blastp", "blastx", "tblastn", "makeblastdb", "blastdbcmd")


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_skesa_bundle(source):
    """Require the reviewed native payload, exact source and redistribution texts."""
    source = Path(source).resolve()
    manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("tool") != "SKESA" or manifest.get("platform") != "windows-x64-ucrt"
            or manifest.get("source_commit") != "c1413581e4f37211892d3c4310d01f3d9a9b3490"):
        raise ValueError("SKESA is not the reviewed native Windows payload.")
    recorded = manifest.get("files", [])
    required = {"skesa.exe", "SKESA-LICENSE.txt", "AGPL-3.0.txt", "skesa-2.4.0-wmlstudio-source.tar.gz"}
    paths = set()
    for item in recorded:
        relative = PurePosixPath(item["path"])
        path = (source / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(source):
            raise ValueError("Unsafe SKESA manifest path.")
        if not path.is_file() or digest(path) != item["sha256"]:
            raise ValueError(f"SKESA payload integrity check failed: {relative}")
        paths.add(relative.as_posix())
    if not required.issubset(paths) or not any(name.startswith("dependency-licenses/") for name in paths):
        raise ValueError("SKESA payload is missing required binaries, sources or license texts.")
    evidence = json.loads((source / "smoke-test.json").read_text(encoding="utf-8"))
    if not evidence.get("all_contigs_match_reference") or not evidence.get("native_adapter", {}).get("passed"):
        raise ValueError("SKESA payload has no successful native assembly smoke evidence.")
    return manifest


def stage_skesa_bundle(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    manifest = verify_skesa_bundle(source)
    if destination.exists():
        if verify_skesa_bundle(destination) != manifest:
            raise ValueError("Existing SKESA payload differs; choose a fresh staging directory.")
        return manifest
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="wmlstudio-skesa-stage-", dir=destination.parent) as work:
        target = Path(work) / "payload"
        shutil.copytree(source, target)
        verify_skesa_bundle(target)
        os.replace(target, destination)
    return manifest


def package_url(platform):
    return f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{VERSION}/{PACKAGES[platform]['archive']}"


def stage(platform, destination, archive=None):
    package = PACKAGES[platform]
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        manifest_path = destination / "manifest.json"
        if manifest_path.is_file():
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
            recorded = previous.get("files", [])
            required = {f"bin/{name}{package['suffix']}" for name in PROGRAMS} | {"LICENSE", "README"}
            safe = (isinstance(recorded, list) and recorded
                    and all(isinstance(item, dict) and isinstance(item.get("path"), str)
                            and not PurePosixPath(item["path"]).is_absolute()
                            and ".." not in PurePosixPath(item["path"]).parts for item in recorded))
            if (safe and required.issubset({item["path"] for item in recorded})
                    and previous.get("archive_sha256") == package["sha256"] and all(
                (destination / item["path"]).is_file()
                and digest(destination / item["path"]) == item["sha256"]
                for item in recorded)):
                return previous
        raise ValueError(f"Destination already exists or has a different snapshot: {destination}. Choose a new staging directory.")
    with tempfile.TemporaryDirectory(prefix="wmlstudio-blast-stage-", dir=destination.parent) as temporary:
        temporary = Path(temporary)
        source = Path(archive) if archive else temporary / package["archive"]
        if archive is None:
            request = urllib.request.Request(package_url(platform), headers={"User-Agent": "WMLSTudio-build"})
            with urllib.request.urlopen(request, timeout=60) as response, source.open("wb") as output:
                shutil.copyfileobj(response, output, 1024 * 1024)
        actual = digest(source)
        if actual != package["sha256"]:
            raise ValueError(f"BLAST archive SHA-256 mismatch: expected {package['sha256']}, got {actual}")
        target = temporary / "blast"
        target.mkdir()
        programs = {name + package["suffix"] for name in PROGRAMS}
        files = []
        with tarfile.open(source, "r:gz") as bundle:
            for member in bundle:
                parts = PurePosixPath(member.name).parts
                if not parts or parts[0] != f"ncbi-blast-{VERSION}+" or not member.isfile():
                    continue
                relative = PurePosixPath(*parts[1:])
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("Unsafe BLAST archive member.")
                filename = relative.name
                chosen = (relative.as_posix() in {"LICENSE", "BLAST_PRIVACY", "README", "ncbi_package_info"}
                          or len(relative.parts) == 2 and relative.parts[0] == "bin"
                          and (filename in programs or filename.removesuffix(".manifest") in programs
                               or filename.lower().endswith(".dll"))
                          or relative.parts[0] == "lib")
                if not chosen:
                    continue
                output = target / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                handle = bundle.extractfile(member)
                if handle is None:
                    raise ValueError(f"Cannot extract {member.name}")
                with handle, output.open("wb") as destination_file:
                    shutil.copyfileobj(handle, destination_file)
                if relative.parts[0] == "bin" and not package["suffix"]:
                    output.chmod(0o755)
                files.append({"path": relative.as_posix(), "size": output.stat().st_size,
                              "sha256": digest(output)})
        missing = [name for name in programs if not (target / "bin" / name).is_file()]
        if missing:
            raise ValueError(f"Official archive is missing required programs: {missing}")
        manifest = {"tool": "NCBI BLAST+", "version": VERSION, "platform": platform,
                    "source_url": package_url(platform), "archive_sha256": actual,
                    "license": "NCBI public-domain notice and bundled component terms; see LICENSE",
                    "files": sorted(files, key=lambda item: item["path"])}
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        os.replace(target, destination)
    return manifest


def stage_hydra_database(source, destination):
    """Stage only the reviewed NCBI starter data; never redistribute other providers implicitly."""
    from wmlstudio.hydra_runtime import _database_provenance, installed_databases

    source, destination = Path(source).resolve(), Path(destination).resolve()
    entries = installed_databases(source)
    if set(entries) != {"ncbi", "protein"}:
        raise ValueError("The portable starter must contain exactly the NCBI nucleotide and protein databases.")
    provenance = _database_provenance(source, ["ncbi", "protein"])
    if destination.exists():
        previous = destination / "snapshot_provenance.json"
        recorded = json.loads(previous.read_text()) if previous.is_file() else {}
        if recorded.get("manifest_sha256") == provenance["manifest_sha256"]:
            current = _database_provenance(destination, ["ncbi", "protein"])
            if current["databases"] == provenance["databases"] == recorded.get("databases"):
                return current
        raise ValueError(f"Starter staging already contains a different snapshot: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".hydra-starter-stage-", dir=destination.parent) as temporary:
        target = Path(temporary) / "starter"
        target.mkdir()
        for relative in ("nucl/ncbi", "prot/protein", "mutation"):
            if (source / relative).is_dir():
                shutil.copytree(source / relative, target / relative)
        shutil.copy2(source / "manifest.json", target / "manifest.json")
        provenance["redistribution_scope"] = "NCBI AMRFinderPlus reference data only; provider terms and attribution retained separately."
        (target / "snapshot_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        os.replace(target, destination)
    return provenance


def download_hydra_starter(destination):
    """An explicit build-time download; a valid already-staged snapshot is reused."""
    from wmlstudio.hydra_runtime import _database_provenance, update_databases

    destination = Path(destination).resolve()
    if destination.exists():
        provenance_path = destination / "snapshot_provenance.json"
        if provenance_path.is_file():
            recorded = json.loads(provenance_path.read_text(encoding="utf-8"))
            current = _database_provenance(destination, ["ncbi", "protein"])
            if current["manifest_sha256"] == recorded["manifest_sha256"] and current["databases"] == recorded["databases"]:
                return current
        raise ValueError("An existing HYDRA starter staging directory failed snapshot verification.")
    with tempfile.TemporaryDirectory(prefix="wmlstudio-build-reference-") as temporary:
        result = update_databases(Path(temporary) / "ncbi", ["protein", "ncbi"],
                                  progress=lambda _i, _n, message: print(message, flush=True))
        return stage_hydra_database(result["database_root"], destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=PACKAGES, required=True)
    parser.add_argument("--archive", type=Path, help="Already downloaded official archive, still hash checked")
    parser.add_argument("--destination", type=Path, default=ROOT / "src/wmlstudio/resources/tools/blast")
    reference = parser.add_mutually_exclusive_group()
    reference.add_argument("--hydra-source", type=Path, help="Stage an already validated NCBI-only HYDRA starter snapshot")
    reference.add_argument("--download-hydra-starter", action="store_true", help="Explicitly download/build the NCBI-only starter (or verify/reuse one already staged)")
    parser.add_argument("--hydra-destination", type=Path, default=ROOT / "src/wmlstudio/resources/hydra/starter")
    parser.add_argument("--skesa-source", type=Path, help="Verified native Windows SKESA artifact directory")
    parser.add_argument("--require-skesa", action="store_true", help="Fail unless a verified native Windows assembler is staged")
    args = parser.parse_args()
    manifest = stage(args.platform, args.destination, args.archive)
    print(f"Staged BLAST+ {VERSION} ({args.platform}): {len(manifest['files'])} verified files in {args.destination}")
    if args.hydra_source:
        provenance = stage_hydra_database(args.hydra_source, args.hydra_destination)
        print(f"Staged NCBI-only HYDRA starter with manifest SHA-256 {provenance['manifest_sha256']}")
    elif args.download_hydra_starter:
        provenance = download_hydra_starter(args.hydra_destination)
        print(f"Staged NCBI-only HYDRA starter with manifest SHA-256 {provenance['manifest_sha256']}")
    skesa = ROOT / "src/wmlstudio/resources/tools/skesa"
    if args.skesa_source:
        stage_skesa_bundle(args.skesa_source, skesa)
    if args.require_skesa:
        verify_skesa_bundle(skesa)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
