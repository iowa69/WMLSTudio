"""Batch access to the same scientific engine used by the native application."""

import argparse
import json
import sys
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.comparison import minimum_spanning_forest, pairwise_distances
from wmlstudio.export import ensure_separate_destination, export_results, write_distances
from wmlstudio.sequence import inspect_sequence, validate_read_pair
from wmlstudio.typing import call_assembly, load_scheme


def main(argv=None):
    parser = argparse.ArgumentParser(prog="wmlstudio-cli", description="Reproducible WMLSTudio sequence analysis")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    qc = commands.add_parser("qc", help="Validate FASTA or inspect a bounded prefix of Phred+33 FASTQ")
    qc.add_argument("inputs", nargs="+", type=Path)
    qc.add_argument("--max-reads", type=int, default=100000)
    qc.add_argument("--output", "-o", type=Path)
    qc.add_argument("--format", choices=["json", "csv", "tsv", "html"], default="json")
    typing = commands.add_parser("type", help="Exact known-allele matching on assemblies")
    typing.add_argument("inputs", nargs="+", type=Path)
    typing.add_argument("--scheme", required=True, type=Path)
    typing.add_argument("--output", "-o", type=Path)
    typing.add_argument("--format", choices=["json", "csv", "tsv", "html"], default="json")
    pair = commands.add_parser("check-pair", help="Validate paired FASTQ record identities")
    pair.add_argument("r1", type=Path)
    pair.add_argument("r2", type=Path)
    pair.add_argument("--max-reads", type=int, default=100000)
    compare = commands.add_parser("compare", help="Compare an exported result JSON file")
    compare.add_argument("input", type=Path)
    compare.add_argument("--min-overlap", type=float, default=0.95)
    compare.add_argument("--output", "-o", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "check-pair":
            print(json.dumps(validate_read_pair(args.r1, args.r2, args.max_reads), indent=2))
            return 0
        if args.command == "compare":
            document = json.loads(args.input.read_text(encoding="utf-8"))
            if not isinstance(document, dict) or not isinstance(document.get("samples"), list):
                raise ValueError("Expected a WMLSTudio JSON export containing a samples list.")
            if not all(isinstance(sample, dict) for sample in document["samples"]):
                raise ValueError("Every sample in the comparison input must be an object.")
            results = [r for r in document["samples"] if r.get("alleles") and r.get("job_status", "completed") == "completed"]
            pairs = pairwise_distances(results, args.min_overlap)
            if args.output:
                ensure_separate_destination(args.output, [args.input])
                write_distances(pairs, args.output, args.min_overlap)
            else:
                print(json.dumps({"pairs": pairs, "forest": minimum_spanning_forest(results, args.min_overlap)}, indent=2))
            return 0
        scheme = load_scheme(args.scheme) if args.command == "type" else None
        results = []
        failed = False
        for path in args.inputs:
            try:
                result = call_assembly(path, scheme) if scheme else inspect_sequence(path, args.max_reads)
                result["software_version"] = __version__
                results.append(result)
            except (ValueError, OSError) as exc:
                failed = True
                results.append({"sample_name": path.name, "input_path": str(path), "status": "failed", "error": str(exc)})
                print(f"{path.name}: {exc}", file=sys.stderr)
        if args.output:
            export_results(results, args.output, args.format)
        else:
            print(json.dumps({"application": f"WMLSTudio {__version__}", "samples": results}, indent=2))
        return 1 if failed else 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"WMLSTudio: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
