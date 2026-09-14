"""Explicit, restartable full-read FastQC evidence on selected local isolates.

Folder organism names are inventory labels, not independently verified species
or truth labels. Source FASTQs are never altered. This is a measured bounded
three-isolate test, not a claim about throughput for hundreds of isolates.
"""

import argparse
import json
import os
import time
from pathlib import Path

from wmlstudio.fastqc import run_fastqc
from wmlstudio.scheduler import plan_resources, run_bounded
from wmlstudio.sequence import file_sha256, validate_read_pair

CASES = (
    ("Klebsiella pneumoniae (directory label)", "kpneumoniae", "SRR26465495"),
    ("Acinetobacter baumannii (directory label)", "abaumannii", "SRR26465518"),
    ("Enterococcus faecium (directory label)", "efaecium", "DRR428002"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--fastqc-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parallel", type=int, default=3)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    allocation = plan_resources(threads_per_sample=2, memory_gb=1, max_parallel=args.parallel)
    records = []
    def analyse(case, resources, cancelled, progress):
        label, directory, accession = case
        reads = [args.input_root / directory / f"{accession}_{mate}.fastq.gz" for mate in (1, 2)]
        before = [file_sha256(path, cancelled) for path in reads]
        # This limit deliberately exceeds any plausible isolate run; both EOFs
        # and complete-file status are asserted, never inferred from a prefix.
        pair = validate_read_pair(*reads, max_reads=2**63-1, cancelled=cancelled)
        if not pair["complete_file"] or pair["sampled"]:
            raise ValueError("The intended complete-file pair check was not completed")
        destination = args.output / accession
        saved = destination / "fastqc-result.json"
        resumed = saved.is_file()
        if resumed:
            result = json.loads(saved.read_text(encoding="utf-8"))
            if [entry["sha256"] for entry in result["inputs"]] != before:
                raise ValueError("Existing validation belongs to different input bytes")
        else:
            result = run_fastqc(reads, destination, threads=resources.threads_per_sample,
                                memory_gb=resources.memory_gb, root=args.fastqc_root,
                                cancelled=cancelled, progress=progress)
        after = [file_sha256(path, cancelled) for path in reads]
        if after != before:
            raise ValueError("An original read input changed during validation")
        return {"inventory_label": label, "accession": accession, "pair_validation": pair,
                "input_sha256_before": before, "input_sha256_after": after,
                "fastqc": result, "resumed": resumed}
    started = time.monotonic()
    usage_before = None
    if os.name != "nt":
        import resource
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    run_bounded(CASES, analyse, allocation, on_result=lambda case, result: records.append(result),
                progress=lambda done, total, message: print(f"[{done}/{total}] {message}", flush=True))
    report = {"status": "completed", "elapsed_seconds": round(time.monotonic()-started, 3),
              "resource_plan": allocation.to_dict(), "records": records,
              "claim_scope": "Three real paired-read samples, complete FastQC and pair-identity validation; no trimming; no hundreds-sample throughput claim."}
    if usage_before is not None:
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        report["child_cpu_seconds"] = round(usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime, 3)
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "records"}, indent=2))


if __name__ == "__main__":
    main()
