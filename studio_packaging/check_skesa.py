"""Real SKESA smoke assembly against a deterministic synthetic ground truth.

This is an executable/platform test, not validation on clinical isolates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def check(executable, output=None, adapter_source=None, expected_version="2.4.0"):
    executable = Path(executable).resolve()
    version = subprocess.check_output([str(executable), "--version"], stderr=subprocess.STDOUT,
                                      text=True, timeout=20).strip()
    if expected_version not in version:
        raise RuntimeError(f"Unexpected SKESA version: {version}")
    rng = random.Random(41027)
    genome = "".join(rng.choices("ACGT", k=12000))
    with tempfile.TemporaryDirectory(prefix="wmlstudio-skesa-smoke-") as temporary:
        directory = Path(temporary)
        with (directory / "mate1.fastq").open("w", newline="\n") as first, \
                (directory / "mate2.fastq").open("w", newline="\n") as second:
            for number, start in enumerate(range(0, len(genome) - 400 + 1, 5)):
                for mate, handle, sequence in (
                    (1, first, genome[start:start + 150]),
                    (2, second, reverse_complement(genome[start + 250:start + 400])),
                ):
                    handle.write(f"@pair{number}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        command = [str(executable), "--reads", "mate1.fastq,mate2.fastq", "--cores", "2",
                   "--memory", "4", "--min_contig", "200", "--contigs_out", "contigs.fasta"]
        environment = os.environ.copy()
        if os.name == "nt":
            system = environment.get("SystemRoot", r"C:\Windows")
            environment["PATH"] = str(executable.parent) + os.pathsep + system + r"\System32"
        completed = subprocess.run(command, cwd=directory, capture_output=True, text=True,
                                   env=environment, timeout=180, check=False)
        if completed.returncode:
            raise RuntimeError(f"Real SKESA smoke exited {completed.returncode}:\n"
                               f"{completed.stdout}\n{completed.stderr}")
        contigs = []
        for line in (directory / "contigs.fasta").read_text().splitlines():
            if line.startswith(">"):
                contigs.append("")
            elif line.strip() and contigs:
                contigs[-1] += line.strip().upper()
            elif line.strip():
                raise RuntimeError("Assembler output is not FASTA.")
        if not contigs or max(map(len, contigs)) < 9000:
            raise RuntimeError("Synthetic assembly did not recover a >=9 kbp contig.")
        if any(not seq or (seq not in genome and reverse_complement(seq) not in genome) for seq in contigs):
            raise RuntimeError("A synthetic contig disagrees with the known reference sequence.")
        report = {"tool_version": version, "test_kind": "synthetic platform smoke; not clinical validation",
                  "reference_sha256": hashlib.sha256(genome.encode()).hexdigest(),
                  "reference_bases": len(genome), "pairs": number + 1, "contigs": len(contigs),
                  "assembled_bases": sum(map(len, contigs)), "largest_contig": max(map(len, contigs)),
                  "all_contigs_match_reference": True, "command": command,
                  "stderr": completed.stderr[-12000:]}
        if adapter_source:
            sys.path.insert(0, str(Path(adapter_source).resolve() / "src"))
            from wmlstudio.assembly import run_skesa
            from wmlstudio.sequence import iter_sequences

            # Exercise a real native process from a Unicode working directory.
            reads_directory = directory / "Unicode reads λ, paired"
            reads_directory.mkdir()
            original_reads = []
            for mate in (1, 2):
                path = reads_directory / f"isolate, λ_R{mate}.fastq"
                shutil.copy2(directory / f"mate{mate}.fastq", path)
                original_reads.append(path)
            before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in original_reads]
            result = run_skesa(*original_reads, reads_directory / "assembled λ",
                               executable=executable, threads=2, memory_gb=4)
            assembled = [record.sequence for record in iter_sequences(result["assembly_path"])]
            if not assembled or max(map(len, assembled)) < 9000 or any(
                    seq not in genome and reverse_complement(seq) not in genome for seq in assembled):
                raise RuntimeError("Native adapter assembly does not match synthetic ground truth.")
            after = [hashlib.sha256(path.read_bytes()).hexdigest() for path in original_reads]
            if before != after or result["pairing"]["records_checked"] != number + 1:
                raise RuntimeError("Native adapter changed inputs or did not validate every pair.")
            report["native_adapter"] = {
                "passed": True, "unicode_comma_paths": True, "original_inputs_unchanged": True,
                "assembly_sha256": result["provenance"]["assembly_sha256"],
                "read_pairs": result["pairing"]["records_checked"], "qc": result["qc"],
            }
    if output:
        Path(output).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable")
    parser.add_argument("--output")
    parser.add_argument("--adapter-source", help="WMLSTudio checkout to test its real native process adapter")
    parser.add_argument("--expected-version", default="2.4.0")
    arguments = parser.parse_args()
    check(arguments.executable, arguments.output, arguments.adapter_source, arguments.expected_version)
