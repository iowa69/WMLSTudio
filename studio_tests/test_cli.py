"""Command-line contracts tested through the actual process exit and JSON output."""

import csv
import gzip
import hashlib
import json
import subprocess
import sys

import pytest

from wmlstudio import __version__

ARC1 = "AACCGTACGTTAG"
ARC2 = "AACCGTTCGTTAG"
GYR1 = "TTGGCATACCTGA"


@pytest.fixture
def typing_inputs(tmp_path):
    scheme = tmp_path / "known_scheme"
    scheme.mkdir()
    (scheme / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_2\n{ARC2}\n")
    (scheme / "gyrB.tfa").write_text(f">gyrB_1\n{GYR1}\n")
    (scheme / "profiles.tsv").write_text("ST\tarcA\tgyrB\n7\t1\t1\n8\t2\t1\n")
    first = tmp_path / "isolate_A.fasta"
    first.write_text(f">one\n{ARC1}\n>two\n{GYR1}\n")
    second = tmp_path / "isolate_B.fasta"
    second.write_text(f">one\n{ARC2}\n>two\n{GYR1}\n")
    partial = tmp_path / "partial.fasta"
    partial.write_text(f">one\n{ARC1}\n")
    return scheme, first, second, partial


def cli(*arguments):
    return subprocess.run(
        [sys.executable, "-m", "wmlstudio.cli", *map(str, arguments)],
        capture_output=True, text=True, timeout=20, check=False,
    )


def test_version_and_help_are_runnable():
    version = cli("--version")
    assert version.returncode == 0
    assert version.stdout.strip() == __version__
    help_result = cli("--help")
    assert help_result.returncode == 0
    assert all(command in help_result.stdout for command in ("qc", "type", "check-pair", "compare"))
    assert "fastqc" in help_result.stdout and "characterize" in help_result.stdout


def test_characterization_cli_not_run_is_not_negative_and_protects_input(tmp_path):
    source = tmp_path / "input.fasta"
    source.write_text(">contig\nACGTACGT\n")
    original = source.read_bytes()
    result = cli("characterize", source, "--no-species", "--no-virulence")
    assert result.returncode == 0, result.stderr
    output = json.loads(result.stdout)
    assert output["status"] == "not_run"
    assert output["species_evidence"]["status"] == "not_run"
    assert output["drug_associations"]["status"] == "not_run"
    invalid = cli("characterize", source, "--output", source)
    assert invalid.returncode == 1 and "Traceback" not in invalid.stderr
    assert source.read_bytes() == original


def test_fastqc_cli_rejects_invalid_resources_without_traceback(tmp_path):
    result = cli("fastqc", tmp_path / "absent.fastq", "--output", tmp_path / "reports", "--threads", "0")
    assert result.returncode == 1
    assert "thread allocation" in result.stderr and "Traceback" not in result.stderr
    assert not (tmp_path / "reports").exists()


def test_type_produces_expected_sts_and_provenance(typing_inputs):
    scheme, first, second, _ = typing_inputs
    completed = cli("type", first, second, "--scheme", scheme)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    document = json.loads(completed.stdout)
    assert document["application"] == f"WMLSTudio {__version__}"
    results = document["samples"]
    assert [result["st"] for result in results] == ["7", "8"]
    assert [result["status"] for result in results] == ["complete", "complete"]
    assert results[0]["alleles"] == {"arcA": "1", "gyrB": "1"}
    assert results[1]["alleles"] == {"arcA": "2", "gyrB": "1"}
    assert results[0]["input_sha256"] == hashlib.sha256(first.read_bytes()).hexdigest()
    assert len(results[0]["scheme_digest"]) == 64
    assert results[0]["scheme_digest"] == results[1]["scheme_digest"]


def test_type_batch_failure_returns_nonzero_and_keeps_valid_results(tmp_path, typing_inputs):
    scheme, first, _, _ = typing_inputs
    malformed = tmp_path / "bad.fasta"
    malformed.write_text(">empty_contig\n")
    missing = tmp_path / "missing.fasta"
    completed = cli("type", malformed, first, missing, "--scheme", scheme)
    assert completed.returncode == 1
    results = json.loads(completed.stdout)["samples"]
    assert [result["status"] for result in results] == ["failed", "complete", "failed"]
    assert results[1]["st"] == "7"
    assert "empty sequence" in results[0]["error"]
    assert "bad.fasta" in completed.stderr
    assert "missing.fasta" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_type_does_not_assign_st_to_raw_reads(tmp_path, typing_inputs):
    scheme, first, _, _ = typing_inputs
    reads = tmp_path / "disguised.fasta"
    reads.write_text(f"@r\n{ARC1}\n+\n{'I' * len(ARC1)}\n")
    completed = cli("type", reads, first, "--scheme", scheme)
    assert completed.returncode == 1
    results = json.loads(completed.stdout)["samples"]
    assert results[0]["status"] == "failed"
    assert "FASTA assembly" in results[0]["error"]
    assert results[1]["st"] == "7"


def test_qc_fastq_sampling_keeps_full_compressed_input_hash(tmp_path):
    reads = tmp_path / "reads.fastq.gz"
    raw = b"@first\nACGT\n+\nIIII\n@second\nTGCA\n+\n!!!!\n"
    reads.write_bytes(gzip.compress(raw))
    completed = cli("qc", reads, "--max-reads", "1")
    assert completed.returncode == 0, completed.stderr
    result = json.loads(completed.stdout)["samples"][0]
    assert result["kind"] == "fastq"
    assert result["qc"]["records"] == 1
    assert result["qc"]["sampled"] is True
    assert result["qc"]["q30_percent"] == 100
    assert result["input_sha256"] == hashlib.sha256(reads.read_bytes()).hexdigest()
    assert any("not validated" in note for note in result["notes"])
    assert any("Phred+33" in note for note in result["notes"])


def test_qc_batch_invalid_limit_and_malformed_file_report_failures(tmp_path):
    good = tmp_path / "good.fastq"
    good.write_text("@r\nACGT\n+\nIIII\n")
    bad = tmp_path / "bad.fastq"
    bad.write_text("@r\nACGT\n+\nII\n")
    batch = cli("qc", bad, good)
    assert batch.returncode == 1
    rows = json.loads(batch.stdout)["samples"]
    assert rows[0]["status"] == "failed"
    assert rows[1]["qc"]["records"] == 1
    invalid_limit = cli("qc", good, "--max-reads", "0")
    assert invalid_limit.returncode == 1
    assert "positive integer" in invalid_limit.stderr


def test_check_pair_full_prefix_and_mismatch(tmp_path):
    first, second = tmp_path / "one.fq", tmp_path / "two.fq"
    first.write_text("@a/1\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n")
    second.write_text("@a/2\nTGCA\n+\nIIII\n@b/2\nTGCA\n+\nIIII\n")
    full = cli("check-pair", first, second)
    assert full.returncode == 0
    assert json.loads(full.stdout)["status"] == "verified"
    prefix = cli("check-pair", first, second, "--max-reads", "1")
    assert prefix.returncode == 0
    assert json.loads(prefix.stdout)["status"] == "prefix_verified"
    second.write_text("@a/2\nTGCA\n+\nIIII\n@different/2\nTGCA\n+\nIIII\n")
    mismatch = cli("check-pair", first, second)
    assert mismatch.returncode == 1
    assert mismatch.stdout == ""
    assert "identifiers differ at record 2" in mismatch.stderr


def test_compare_real_typing_export_keeps_missing_denominators(tmp_path, typing_inputs):
    scheme, first, second, partial = typing_inputs
    exported = tmp_path / "typing.json"
    typed = cli("type", first, second, partial, "--scheme", scheme, "--output", exported)
    assert typed.returncode == 0, typed.stderr
    compared = cli("compare", exported)
    assert compared.returncode == 0, compared.stderr
    document = json.loads(compared.stdout)
    assert len(document["pairs"]) == 3
    known = [pair for pair in document["pairs"] if pair["comparable"]]
    assert len(known) == 1
    assert known[0]["distance"] == 1
    assert known[0]["shared_loci"] == 2
    assert known[0]["total_loci"] == 2
    assert len(document["forest"]) == 1
    incomplete = [pair for pair in document["pairs"] if not pair["comparable"]]
    assert all(pair["distance"] is None and pair["shared_loci"] == 1 for pair in incomplete)
    assert all(pair["overlap"] == 0.5 for pair in incomplete)
    relaxed = cli("compare", exported, "--min-overlap", "0.5")
    relaxed_document = json.loads(relaxed.stdout)
    assert len(relaxed_document["forest"]) == 2
    assert all(pair["comparable"] for pair in relaxed_document["pairs"])


def test_compare_output_records_requested_overlap_and_rejects_invalid_threshold(
    tmp_path, typing_inputs,
):
    scheme, first, second, _ = typing_inputs
    exported = tmp_path / "typing.json"
    assert cli("type", first, second, "--scheme", scheme, "--output", exported).returncode == 0
    destination = tmp_path / "distances.json"
    completed = cli("compare", exported, "--min-overlap", "0.75", "--output", destination)
    assert completed.returncode == 0
    assert completed.stdout == ""
    report = json.loads(destination.read_text())
    assert report["min_overlap"] == 0.75
    assert report["pairs"][0]["distance"] == 1
    assert "shared" in report["missing_policy"]
    invalid = cli("compare", exported, "--min-overlap", "nan")
    assert invalid.returncode == 1
    assert "finite number" in invalid.stderr


def test_compare_excludes_failed_stored_results_and_mixed_samples(tmp_path, typing_inputs):
    scheme, first, second, _ = typing_inputs
    typed = cli("type", first, second, "--scheme", scheme)
    document = json.loads(typed.stdout)
    stale = dict(document["samples"][0], sample_name="stale", job_status="failed")
    mixed = dict(document["samples"][0], sample_name="mixed", status="mixed")
    document["samples"].extend([stale, mixed])
    source = tmp_path / "stored.json"
    source.write_text(json.dumps(document))
    completed = cli("compare", source)
    assert completed.returncode == 0
    pairs = json.loads(completed.stdout)["pairs"]
    assert len(pairs) == 3
    assert all("stale" not in {pair["source"], pair["target"]} for pair in pairs)
    mixed_pairs = [pair for pair in pairs if "mixed" in {pair["source"], pair["target"]}]
    assert len(mixed_pairs) == 2
    assert all(not pair["comparable"] and pair["distance"] is None for pair in mixed_pairs)


@pytest.mark.parametrize("format", ["json", "csv", "tsv", "html"])
def test_type_writes_requested_export_format(tmp_path, typing_inputs, format):
    scheme, first, _, _ = typing_inputs
    destination = tmp_path / f"report.{format}"
    completed = cli("type", first, "--scheme", scheme, "--format", format, "--output", destination)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == ""
    if format == "json":
        assert json.loads(destination.read_text())["samples"][0]["st"] == "7"
    elif format in {"csv", "tsv"}:
        with destination.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="," if format == "csv" else "\t"))
        assert rows[0]["st"] == "7"
        assert rows[0]["allele:arcA"] == "1"
    else:
        text = destination.read_text()
        assert text.startswith("<!doctype html>")
        assert "isolate_A" in text and "Sequence type</dt><dd>7" in text


def test_cli_top_level_failures_have_no_traceback(tmp_path, typing_inputs):
    _, first, _, _ = typing_inputs
    missing_scheme = cli("type", first, "--scheme", tmp_path / "absent")
    assert missing_scheme.returncode == 1
    assert "Scheme directory does not exist" in missing_scheme.stderr
    invalid_json = tmp_path / "invalid.json"
    invalid_json.write_text("this is not JSON")
    invalid_comparison = cli("compare", invalid_json)
    assert invalid_comparison.returncode == 1
    assert "WMLSTudio:" in invalid_comparison.stderr
    assert "Traceback" not in missing_scheme.stderr + invalid_comparison.stderr


@pytest.mark.parametrize("document", [
    {"samples": [None]},
    {"samples": ["not a result"]},
    {"samples": {"sample": {"alleles": {}}}},
    {"samples": None},
])
def test_compare_rejects_malformed_record_shapes_without_traceback(tmp_path, document):
    source = tmp_path / "bad_shape.json"
    source.write_text(json.dumps(document))
    completed = cli("compare", source)
    assert completed.returncode == 1
    assert completed.stdout == ""
    assert "WMLSTudio:" in completed.stderr
    assert "Traceback" not in completed.stderr


@pytest.mark.parametrize("format", ["json", "csv", "tsv", "html"])
def test_cli_export_cannot_overwrite_original_sequence(typing_inputs, format):
    scheme, first, _, _ = typing_inputs
    original = first.read_bytes()
    completed = cli("type", first, "--scheme", scheme, "--format", format, "--output", first)
    assert completed.returncode == 1
    assert first.read_bytes() == original
    assert "Traceback" not in completed.stderr


def test_cli_lists_organism_modules_with_their_taxa_and_stated_limits():
    listed = cli("characterize", "--list-modules")
    assert listed.returncode == 0, listed.stderr
    document = json.loads(listed.stdout)
    modules = {entry["key"]: entry for entry in document["modules"]}
    assert {"sccmec", "klebsiella_locus_st", "klebsiella_capsule"} <= set(modules)
    assert modules["sccmec"]["genera"] == ["Staphylococcus"]
    assert modules["sccmec"]["species"] == ["aureus"]
    assert modules["klebsiella_capsule"]["genera"] == ["Klebsiella"]
    for entry in modules.values():
        assert "does not establish" in entry["purpose"]
        assert entry["limitations"]
    # The claim boundary is printed where a command-line user actually sees it.
    for project in ("Kleborate", "Kaptive", "AMRFinderPlus", "SCCmecFinder"):
        assert project in document["boundary"]
    assert "none is equivalent to those tools" in document["boundary"]


def test_cli_module_selection_records_applicability_and_never_fakes_a_negative(tmp_path):
    source = tmp_path / "assembly.fasta"
    source.write_text(">contig\nACGTACGTACGT\n")
    completed = cli("characterize", source, "--no-species", "--no-virulence",
                    "--module", "sccmec", "--organism", "Staphylococcus aureus")
    assert completed.returncode == 0, completed.stderr
    output = json.loads(completed.stdout)
    assert output["sccmec"]["applicability"] == "recommended"
    assert output["sccmec"]["status"] == "not_run"
    # A snapshot without the panel is diagnosable, not a negative SCCmec result.
    assert "sccmec.targets" in output["sccmec"]["reason"]
    assert "not_detected" not in json.dumps(output["sccmec"])
    assert output["klebsiella_capsule"]["applicability"] == "off_panel"
    assert output["klebsiella_capsule"]["reason"] == "Assay not selected."


def test_cli_rejects_an_unknown_module_key_without_a_traceback(tmp_path):
    source = tmp_path / "assembly.fasta"
    source.write_text(">contig\nACGT\n")
    rejected = cli("characterize", source, "--module", "kleborate")
    assert rejected.returncode == 1
    assert "Unknown organism module kleborate" in rejected.stderr
    assert "klebsiella_locus_st" in rejected.stderr
    assert "Traceback" not in rejected.stderr
    assert rejected.stdout == ""
    without_input = cli("characterize")
    assert without_input.returncode == 1
    assert "--list-modules" in without_input.stderr and "Traceback" not in without_input.stderr
