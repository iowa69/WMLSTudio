"""Re-run the bounded local regression cohort against the prior WMLST evidence.

Usage: uv run python studio_scripts/validate_local_cohort.py
Requires the original local input paths and artifacts/wmlst-paeruginosa-reference.json.
Results/checkpoints are written outside the input directories; no raw data changes.
"""

import json
from pathlib import Path

from validate_real_data import main as validate

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/2026-09-12_validation"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    checks = [validate([
        "--assembly", str(ROOT.parent / "wmlst/tests/data/example.fna"),
        "--scheme", str(ROOT / "src/wmlstudio/resources/schemes/sepidermidis"),
        "--expected-st", "184", "--output", str(OUTPUT / "sepidermidis.json"),
    ])]
    reference = json.loads((ROOT / "artifacts/wmlst-paeruginosa-reference.json").read_text())
    comparison = []
    for expected in reference:
        accession = Path(expected["filename"]).parents[1].name
        output = OUTPUT / f"{accession}.json"
        checks.append(validate([
            "--assembly", expected["filename"], "--scheme", str(ROOT / "src/wmlstudio/resources/schemes/paeruginosa"),
            "--expected-st", expected["sequence_type"], "--output", str(output),
        ]))
        report = json.loads(output.read_text())
        observed = report["records"][0]["result"]
        alleles_agree = observed["alleles"] == expected["alleles"]
        if not alleles_agree:
            checks.append(1)
        comparison.append({"accession": accession, "expected_st": expected["sequence_type"],
                           "observed_st": observed["st"], "all_alleles_agree": alleles_agree,
                           "input_sha256": observed["input_sha256"], "scheme_digest": observed["scheme_digest"]})
    previous = json.loads((ROOT / "artifacts/paeruginosa-and-fastq-validation.json").read_text())
    fastq = next(record["path"] for record in previous["records"] if record["kind"] == "fastq")
    checks.append(validate(["--fastq", fastq, "--max-reads", "10000", "--output", str(OUTPUT / "real-fastq.json")]))
    destination = OUTPUT / "reference-concordance.json"
    temporary = destination.with_suffix(".tmp")
    temporary.write_text(json.dumps({"comparator": "WMLST 1.2.1 / BLAST 2.12.0, same scheme snapshot",
                                     "scope": "Three P. aeruginosa assemblies; regression concordance, not population accuracy",
                                     "comparisons": comparison, "passed": not any(checks)}, indent=2) + "\n")
    temporary.replace(destination)
    return int(any(checks))


if __name__ == "__main__":
    raise SystemExit(main())
