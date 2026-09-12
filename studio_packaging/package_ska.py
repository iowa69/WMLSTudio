"""Package a native SKA2 build with its pinned source, lock and crate sources.

Run only in the explicit native-tool build workflow, never at application startup.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

REVISION = "fcf9413d2768dc6538d31f664a4bf310651449e7"
VERSION = "0.5.1"


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, capture_output=True,
                          encoding="utf-8", errors="replace", timeout=180, **kwargs)


def smoke(binary):
    """Known synthetic isolated SNPs, plus an identical pair; no patient data."""
    randomizer = random.Random(18091)
    sequence = "".join(randomizer.choices("ACGT", k=24000))
    mutated = list(sequence)
    for position in (500, 3500, 7500, 14500, 21000):
        mutated[position] = next(base for base in "ACGT" if base != sequence[position])
    with tempfile.TemporaryDirectory(prefix="SKA synthetic é ") as directory:
        work = Path(directory)
        for name, bases in (("a", sequence), ("b", sequence), ("c", "".join(mutated))):
            (work / f"{name}.fa").write_text(f">contig\n{bases}\n", encoding="ascii")
        files = work / "inputs.tsv"
        files.write_text("".join(f"{name}\t{work / (name + '.fa')}\n" for name in ("a", "b", "c")), encoding="utf-8")
        build = run([str(binary), "build", "-f", str(files), "-o", str(work / "cohort"),
                     "-k", "31", "--threads", "2"])
        distance = run([str(binary), "distance", str(work / "cohort.skf"), "--threads", "2"])
        rows = distance.stdout.strip().splitlines()
        if len(rows) != 4:
            raise ValueError(f"Unexpected SKA2 pairwise output: {distance.stdout}")
        evidence = {"seed": 18091, "bases": len(sequence), "isolated_snps": 5,
                    "distance_tsv": distance.stdout, "build_log": build.stderr,
                    "passed": False}
        # Column order and exact result checked explicitly against the CLI format.
        header = rows[0].split("\t")
        snp_column = next((index for index, value in enumerate(header) if "snp" in value.lower()), None)
        if snp_column is None:
            raise ValueError(f"SKA2 distance output has no SNP column: {header}")
        observed = {}
        for row in rows[1:]:
            values = row.split("\t")
            observed[tuple(sorted(values[:2]))] = float(values[snp_column])
        if observed != {("a", "b"): 0.0, ("a", "c"): 5.0, ("b", "c"): 5.0}:
            raise ValueError(f"SKA2 did not reproduce the known isolated SNPs: {observed}")
        evidence["passed"] = True
        return evidence


def package(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    revision = run(["git", "rev-parse", "HEAD"], cwd=source).stdout.strip()
    if revision != REVISION:
        raise ValueError("Wrong upstream SKA2 source revision")
    if destination.exists():
        raise ValueError("Choose a fresh SKA2 output directory; existing bundles are not overwritten")
    binary = source / "target/release/ska.exe"
    if not binary.is_file() or not (source / "Cargo.lock").is_file() or not (source / "vendor").is_dir():
        raise ValueError("Native binary, dependency lock and vendored corresponding sources are required")
    evidence = smoke(binary)
    destination.mkdir(parents=True)
    shutil.copy2(binary, destination / "ska.exe")
    for name in ("LICENSE", "NOTICE", "README.md", "Cargo.lock"):
        shutil.copy2(source / name, destination / name)
    (destination / "smoke-test.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    archive = destination / "ska-0.5.1-corresponding-source.tar.gz"
    with tarfile.open(archive, "w:gz") as output:
        for name in ("src", "vendor", "Cargo.toml", "Cargo.lock", "LICENSE", "NOTICE", "README.md"):
            output.add(source / name, arcname="ska-0.5.1/" + name)
    manifest = {"tool": "SKA2", "version": VERSION, "platform": "windows-x64",
                "source_commit": REVISION, "source_url": "https://github.com/bacpop/ska.rust",
                "license": "Apache-2.0; vendored crates retain their own source licenses",
                "rustc": run(["rustc", "--version"]).stdout.strip(),
                "cargo": run(["cargo", "--version"]).stdout.strip(),
                "build": "cargo build --release --locked; RUSTFLAGS=-C target-feature=+crt-static",
                "lock_sha256": sha256(source / "Cargo.lock"),
                "files": [{"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
                          for path in sorted(destination.iterdir()) if path.is_file()],
                "validation": "Known synthetic 0/5-SNP pairwise distances and Unicode/space-containing paths"}
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(json.dumps(package(args.source, args.destination), indent=2))


if __name__ == "__main__":
    main()
