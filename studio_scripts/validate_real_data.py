"""Opt-in, restartable validation of local sequences; never modifies raw inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.sequence import file_sha256, inspect_sequence
from wmlstudio.typing import call_assembly, load_scheme

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLY_SUFFIXES = (".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz",
                     ".fasta.bz2", ".fa.bz2", ".fna.bz2")


def code_digest() -> str:
    digest = hashlib.sha256()
    paths = [Path(__file__), *sorted((ROOT / "src/wmlstudio").glob("*.py"))]
    for path in paths:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def save_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".checkpoint")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, action="append", default=[])
    parser.add_argument("--assembly-root", type=Path, help="Recursively discover a bounded cohort")
    parser.add_argument("--max-assemblies", type=int, default=10)
    parser.add_argument("--scheme", type=Path, help="Local scheme for all selected assemblies")
    parser.add_argument("--expected-st", help="Known answer, allowed only with one assembly")
    parser.add_argument("--fastq", type=Path, action="append", default=[])
    parser.add_argument("--max-reads", type=int, default=10000)
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts/real-data-validation.json")
    parser.add_argument("--resume", action="store_true", help="Reuse exact successful input/code hashes")
    args = parser.parse_args(argv)
    if args.max_assemblies < 1 or args.max_reads < 1:
        parser.error("max-assemblies and max-reads must be positive")
    assemblies = list(args.assembly)
    if args.assembly_root:
        if not args.assembly_root.is_dir():
            parser.error("assembly-root does not exist")
        assemblies += sorted(path for path in args.assembly_root.rglob("*")
                             if path.is_file() and path.name.lower().endswith(ASSEMBLY_SUFFIXES))
    assemblies = list(dict.fromkeys(path.resolve() for path in assemblies))[:args.max_assemblies]
    if assemblies and not args.scheme:
        parser.error("--scheme is required for assembly typing")
    if not assemblies and not args.fastq:
        parser.error("select --assembly, --assembly-root, or --fastq")
    if args.expected_st is not None and len(assemblies) != 1:
        parser.error("--expected-st requires exactly one selected assembly")
    fingerprint = code_digest()
    previous = {}
    if args.resume and args.output.exists():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
    reusable = {}
    if previous.get("code_sha256") == fingerprint:
        reusable = {record["key"]: record for record in previous.get("records", [])
                    if record.get("validation_status") == "passed"}
    start = time.perf_counter()
    scheme = load_scheme(args.scheme) if args.scheme else None
    scheme_seconds = time.perf_counter() - start
    report = {
        "format_version": 1, "created_utc": datetime.now(UTC).isoformat(),
        "application_version": __version__, "code_sha256": fingerprint,
        "platform": platform.platform(), "python": sys.version,
        "scheme": {"path": str(args.scheme.resolve()), "name": scheme.name,
                   "sha256": scheme.digest, "loci": scheme.locus_count,
                   "alleles": scheme.allele_count} if scheme else None,
        "scheme_load_seconds": round(scheme_seconds, 6),
        "notes": ["Local inputs are read only; no sequence data are uploaded.",
                  "FASTQ QC scans a bounded prefix; input SHA-256 covers all on-disk bytes.",
                  "A passing run without an expected ST proves execution, not typing accuracy."],
        "records": [],
    }
    inputs = [(path, "assembly") for path in assemblies]
    inputs += [(path.resolve(), "fastq") for path in args.fastq]
    for path, kind in inputs:
        expected = args.expected_st if kind == "assembly" else None
        key_data = {"path": str(path), "kind": kind, "expected_st": expected,
                    "max_reads": args.max_reads if kind == "fastq" else None,
                    "scheme_digest": scheme.digest if scheme and kind == "assembly" else None}
        key = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
        record = {"key": key, **key_data, "resumed": False}
        started = time.perf_counter()
        try:
            digest = file_sha256(path)
            if key in reusable and reusable[key].get("input_sha256") == digest:
                record = {**reusable[key], "resumed": True}
            else:
                result = call_assembly(path, scheme) if kind == "assembly" else inspect_sequence(
                    path, max_reads=args.max_reads)
                record.update(input_sha256=digest, result=result)
                if result.get("input_sha256") != digest:
                    raise ValueError("Input changed between fingerprinting and analysis")
                if kind == "fastq" and result.get("kind") != "fastq":
                    raise ValueError("The selected FASTQ input is not FASTQ")
                if expected is not None and str(result.get("st")) != str(expected):
                    raise ValueError(f"Expected ST {expected}, observed {result.get('st')}")
                record["validation_status"] = "passed"
                record["elapsed_seconds"] = round(time.perf_counter() - started, 6)
        except Exception as error:
            record.update(validation_status="failed", error=f"{type(error).__name__}: {error}",
                          elapsed_seconds=round(time.perf_counter() - started, 6))
        report["records"].append(record)
        report["passed"] = sum(item["validation_status"] == "passed" for item in report["records"])
        report["failed"] = len(report["records"]) - report["passed"]
        save_report(args.output, report)
        print(f"{record['validation_status']}: {path.name} ({record['elapsed_seconds']:.3f} s)",
              flush=True)
    print(f"Report: {args.output.resolve()}")
    return 1 if report["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
