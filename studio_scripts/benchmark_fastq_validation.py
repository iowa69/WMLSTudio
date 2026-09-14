"""Measure a bounded real FASTQ parser prefix, before/after an equivalent check.

This is not an assembler/cohort throughput claim. All original file bytes are
hashed before and after; observed records must have identical content hashes.
"""

import argparse
import hashlib
import json
import statistics
import time

import wmlstudio.sequence as sequence


class LegacyQualityCheck:
    def search(self, text):
        return any(ord(char) < 33 or ord(char) > 126 for char in text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--records", type=int, default=100000)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if not 1 <= args.records <= 1000000 or not 1 <= args.runs <= 10:
        raise ValueError("Keep benchmark bounded:1..1,000,000 records and1..10 runs")
    original_hash = sequence.file_sha256(args.input)
    optimized = sequence._INVALID_QUALITY
    results = []
    try:
        for repeat in range(args.runs):
            for name in (("legacy", "optimized") if repeat % 2 == 0 else ("optimized", "legacy")):
                sequence._INVALID_QUALITY = LegacyQualityCheck() if name == "legacy" else optimized
                started = time.perf_counter()
                checksum, count = hashlib.sha256(), 0
                with sequence.SequenceReader(args.input) as reader:
                    if reader.kind != "fastq":
                        raise ValueError("Benchmark input must be FASTQ")
                    for record in reader:
                        checksum.update((record.name + "\0" + record.sequence + "\0" + record.quality + "\n").encode("ascii"))
                        count += 1
                        if count >= args.records:
                            break
                    complete = reader.complete
                result = {"method": name, "repeat": repeat, "seconds": time.perf_counter()-started,
                          "records": count, "record_sha256": checksum.hexdigest(), "complete_file": complete}
                results.append(result)
                print(json.dumps(result), flush=True)
    finally:
        sequence._INVALID_QUALITY = optimized
    if len({(row["records"], row["record_sha256"]) for row in results}) != 1:
        raise ValueError("Optimized and original parsing results differ")
    if sequence.file_sha256(args.input) != original_hash:
        raise ValueError("Original FASTQ bytes changed")
    medians = {name: statistics.median(row["seconds"] for row in results if row["method"] == name)
               for name in ("legacy", "optimized")}
    report = {"input_sha256": original_hash, "unchanged": True, "records_limit": args.records,
              "results": results, "median_seconds": medians, "observed_ratio": medians["legacy"] / medians["optimized"],
              "scope": "Bounded parser-only benchmark, identical records; not end-to-end assembly or hundreds-isolate throughput."}
    from pathlib import Path
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
