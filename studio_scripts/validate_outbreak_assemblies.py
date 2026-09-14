"""Plan/execute a bounded native SKESA technical test on explicitly reviewed reads.

This deliberately retains existing FastQC FAIL flags. It is an engine/input-
integrity experiment, not a QC-pass, taxonomy certification, or clinical assay.
Requires an explicit --execute switch after inspecting the saved plan.
"""

import argparse
import json
import os
import shutil
import time
from pathlib import Path

from wmlstudio.assembly import run_skesa
from wmlstudio.scheduler import plan_resources, run_bounded
from wmlstudio.sequence import file_sha256, inspect_sequence
from wmlstudio.typing import call_assembly, load_scheme


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--fastqc-validation", type=Path, required=True)
    parser.add_argument("--skesa", type=Path, required=True)
    parser.add_argument("--schemes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parallel", type=int, default=3)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    allocation = plan_resources(threads_per_sample=4, memory_gb=8, max_parallel=args.parallel)
    qc = json.loads(args.fastqc_validation.read_text(encoding="utf-8"))
    by_accession = {row["accession"]: row for row in qc["records"]}
    cases = [("kpneumoniae", "SRR26465495", "klebsiella"),
             ("abaumannii", "SRR26465518", "abaumannii"),
             ("efaecium", "DRR428002", "efaecium")]
    plan = {"status": "planned", "technical_validation_only": True,
            "read_quality_status": "FastQC FAIL flags retained; explicit technical-test review permits untrimmed assembly",
            "resource_allocation": allocation.to_dict(), "skesa": str(args.skesa.resolve()),
            "skesa_sha256": file_sha256(args.skesa), "free_disk_bytes": shutil.disk_usage(args.output).free,
            "cases": [{"directory_label": folder, "accession": accession,
                       "scheme_comparator": scheme, "fastqc_flags": [report["modules"] for report in by_accession[accession]["fastqc"]["reports"]]}
                      for folder, accession, scheme in cases],
            "domain_triage": {"status": "not_run", "reason": "QuickClade container runtime/reference unavailable; no certified domain/species assignment is inferred from a directory name or ST."},
            "quast": "Run separately on final contigs using installed QUAST5.3.0; no completeness claim from N50 alone."}
    if plan["free_disk_bytes"] < 10 * 1024**3:
        raise ValueError("At least 10 GiB free scratch space is required")
    plan_path = args.output / "validation-plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in plan.items() if key != "cases"}, indent=2), flush=True)
    if not args.execute:
        return
    records = []
    def assemble(case, resources, cancelled, progress):
        folder, accession, scheme_name = case
        originals = [args.input_root / folder / f"{accession}_{mate}.fastq.gz" for mate in (1, 2)]
        before = [file_sha256(path, cancelled) for path in originals]
        if before != by_accession[accession]["input_sha256_before"]:
            raise ValueError("Input FASTQs differ from the reviewed FastQC experiment")
        work = args.output / accession
        aliases = work / "original read copies é"
        aliases.mkdir(parents=True, exist_ok=True)
        copied = []
        for mate, (source, checksum) in enumerate(zip(originals, before, strict=True), 1):
            target = aliases / f"raw R{mate} é.fastq.gz"
            if not target.exists():
                shutil.copy2(source, target)
            if file_sha256(target, cancelled) != checksum:
                raise ValueError("Unicode-named raw copy differs from its original")
            copied.append(target)
        destination = work / "assembly output é"
        result_path = destination / "assembly.json"
        resumed = result_path.is_file()
        if resumed:
            result = json.loads(result_path.read_text(encoding="utf-8"))
            contigs = Path(result["assembly_path"])
            if not contigs.is_file() or file_sha256(contigs) != result["provenance"]["assembly_sha256"]:
                raise ValueError("Existing assembler result failed resume checksum validation")
        else:
            result = run_skesa(*copied, destination, executable=args.skesa,
                               threads=resources.threads_per_sample, memory_gb=resources.memory_gb,
                               cancelled=cancelled, progress=progress)
        contigs = Path(result["assembly_path"])
        metrics = inspect_sequence(contigs, cancelled=cancelled)
        selected_scheme = args.schemes / scheme_name
        typing = call_assembly(contigs, load_scheme(selected_scheme), cancelled=cancelled) if selected_scheme.is_dir() else {
            "status": "not_run", "reason": f"Scheme directory unavailable: {scheme_name}"}
        after = [file_sha256(path, cancelled) for path in originals]
        if after != before or [file_sha256(path, cancelled) for path in copied] != before:
            raise ValueError("Original inputs or Unicode read copies changed")
        record = {"accession": accession, "resumed": resumed, "directory_label": folder,
                  "input_sha256_before": before, "input_sha256_after": after,
                  "fastqc_modules": [report["modules"] for report in by_accession[accession]["fastqc"]["reports"]],
                  "assembly": result, "native_contiguity_qc": metrics, "exact_scheme_call": typing,
                  "interpretation": "Technical engine/identity check on flagged reads; MLST agreement is not independent species confirmation, purity or clinical validation."}
        (work / "validation.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        return record
    started = time.monotonic()
    run_bounded(cases, assemble, allocation, on_result=lambda case, result: records.append(result),
                progress=lambda done, total, message: print(f"[{done}/{total}] {message}", flush=True))
    report = {**plan, "status": "completed", "elapsed_seconds": round(time.monotonic()-started, 3),
              "platform": os.name, "records": records}
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "seconds": report["elapsed_seconds"],
                      "results": [{"accession": row["accession"], "st": row["exact_scheme_call"].get("st"),
                                   "qc": row["native_contiguity_qc"]["qc"]} for row in records]}, indent=2))


if __name__ == "__main__":
    main()
