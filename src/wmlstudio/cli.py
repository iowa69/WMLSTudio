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


def parse_organism(value):
    """'Klebsiella pneumoniae' -> ('Klebsiella', 'pneumoniae'); a bare genus keeps an empty species."""
    parts = (value or "").split(maxsplit=1)
    if not parts:
        return None
    return parts[0], parts[1] if len(parts) > 1 else ""


def describe_modules(modules):
    """What each organism-specific module is, which taxa it covers, and what it does not establish."""
    from wmlstudio.organism_modules import BOUNDARY
    return {"boundary": BOUNDARY,
            "modules": [{"key": key, "title": module.title, "column_title": module.column_title,
                         "purpose": module.purpose, "genera": sorted(module.match.genera),
                         "species": sorted(module.match.species),
                         "excluded_species": sorted(module.match.exclude_species),
                         "reference_sections": list(module.manifest_sections),
                         "limitations": list(module.limitations)}
                        for key, module in sorted(modules.items())]}


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
    fastqc = commands.add_parser("fastqc", help="Original FastQC complete-file reports with bundled Java; no trimming")
    fastqc.add_argument("inputs", nargs="+", type=Path)
    fastqc.add_argument("--output", "-o", type=Path, required=True)
    fastqc.add_argument("--threads", type=int, default=2)
    fastqc.add_argument("--memory-gb", type=int, default=1)
    characterization = commands.add_parser("characterize", help="Independent ANI, virulence and explicitly linked accessory evidence")
    characterization.add_argument("input", nargs="?", type=Path)
    characterization.add_argument("--references", type=Path)
    characterization.add_argument("--no-species", action="store_true")
    characterization.add_argument("--no-virulence", action="store_true")
    characterization.add_argument("--module", action="append", metavar="KEY",
        help="Run an organism-specific typing module; repeat for several. See --list-modules")
    characterization.add_argument("--list-modules", action="store_true",
        help="Print the registered organism-specific modules, their taxa and their stated limits")
    characterization.add_argument("--organism", metavar="\"Genus species\"",
        help="Record which taxa this isolate belongs to; modules are labelled applicable or off-panel, never skipped")
    characterization.add_argument("--hydra-report", type=Path)
    characterization.add_argument("--hydra-sample")
    characterization.add_argument("--threads", type=int, default=2)
    characterization.add_argument("--output", "-o", type=Path)
    ska = commands.add_parser("ska", help="Research-only native split-kmer SNP evidence; not a transmission decision")
    ska.add_argument("inputs", nargs="+", type=Path)
    ska.add_argument("--output", "-o", type=Path, required=True)
    ska.add_argument("--threads", type=int, default=4)
    ska.add_argument("--k", type=int, default=31)
    ska.add_argument("--min-shared-fraction", type=float, default=.95)
    ska.add_argument("--reference", type=Path)
    ska.add_argument("--annotation", type=Path)
    ska.add_argument("--annotation-reference-sha256", help="Explicit user-declared reference binding for GFF3 lacking embedded FASTA")
    ska.add_argument("--mask-bed", type=Path)
    ska.add_argument("--coding-only", action="store_true")
    ska.add_argument("--no-repeat-mask", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "ska":
            from wmlstudio.sequence import file_sha256, sample_name
            from wmlstudio.ska_runtime import run_ska
            sources = [path.resolve() for path in args.inputs]
            if len(set(sources)) != len(sources):
                raise ValueError("Select distinct assembly files for SKA2")
            samples = [{"id": f"cli-{index+1}", "name": sample_name(path), "input_path": str(path),
                        "input_sha256": file_sha256(path)} for index, path in enumerate(sources)]
            result = run_ska(samples, args.output, threads=args.threads, k=args.k,
                min_shared_fraction=args.min_shared_fraction, reference_path=args.reference,
                annotation_path=args.annotation, annotation_reference_sha256=args.annotation_reference_sha256,
                mask_bed=args.mask_bed, coding_only=args.coding_only, repeat_mask=not args.no_repeat_mask)
            print(json.dumps(result, indent=2))
            return 0
        if args.command == "fastqc":
            from wmlstudio.fastqc import run_fastqc
            result = run_fastqc(args.inputs, args.output, threads=args.threads, memory_gb=args.memory_gb)
            print(json.dumps(result, indent=2))
            return 0
        if args.command == "characterize":
            from wmlstudio.characterization import characterize_assembly
            from wmlstudio.characterization_refs import bundled_reference_root
            from wmlstudio.organism_modules import registered_modules
            modules = registered_modules()
            if args.list_modules:
                print(json.dumps(describe_modules(modules), indent=2))
                return 0
            if args.input is None:
                raise ValueError("Provide an assembly to characterize, or use --list-modules.")
            if args.threads < 1:
                raise ValueError("Threads must be positive")
            unknown = [key for key in (args.module or []) if key not in modules]
            if unknown:
                raise ValueError(f"Unknown organism module {', '.join(sorted(unknown))}. "
                                 f"Available: {', '.join(sorted(modules))}.")
            if args.output:
                ensure_separate_destination(args.output, [args.input, *([args.hydra_report] if args.hydra_report else [])])
            result = characterize_assembly(args.input, args.references or bundled_reference_root(),
                args.hydra_report, threads=args.threads, hydra_sample_name=args.hydra_sample,
                species=not args.no_species, virulence=not args.no_virulence,
                modules={key: True for key in (args.module or [])},
                organism=parse_organism(args.organism))
            encoded = json.dumps(result, indent=2)
            if args.output:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(encoded + "\n", encoding="utf-8")
            else:
                print(encoded)
            return 1 if result["status"] == "failed" else 0
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
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as exc:
        print(f"WMLSTudio: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
