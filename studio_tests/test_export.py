import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from wmlstudio.export import (
    export_results,
    write_csv,
    write_distances,
    write_html,
    write_json,
    write_tsv,
)


@pytest.fixture
def records():
    return [{
        "id": "stable-id", "name": '=HYPERLINK("malicious")', "status": "completed",
        "input_path": '/local/<img src=x onerror="alert(1)">.fasta',
        "metadata": {"site": "Dublin & Cork", "note": "<script>alert(1)</script>"},
        "result": {
            "sample_name": "old name", "status": "complete", "scheme": "Demo",
            "scheme_digest": "abc123", "st": "1", "qc": {"contigs": 2},
            "notes": ["<script>alert(2)</script>"], "alleles": {"xyz": None, "abc": "7"},
            "calls": [{"locus": "abc", "status": "exact"}],
        },
    }, {
        "id": "queued-id", "name": "Pending isolate", "status": "queued",
        "input_path": "/local/pending.fasta", "metadata": {}, "result": None,
    }]


@pytest.mark.parametrize("writer,delimiter,suffix", [(write_csv, ",", ".csv"), (write_tsv, "\t", ".tsv")])
def test_tabular_export_retains_all_records_and_guards_formula(records, tmp_path, writer, delimiter, suffix):
    destination = tmp_path / f"result{suffix}"
    writer(records, destination)
    with destination.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    assert len(rows) == 2
    assert rows[0]["sample_name"].startswith("'=HYPERLINK")
    assert rows[0]["sample_id"] == "stable-id"
    assert rows[0]["allele:abc"] == "7"
    assert rows[0]["allele:xyz"] == ""
    assert rows[0]["status"] == "complete"
    assert rows[1]["job_status"] == "queued"
    assert json.loads(rows[0]["metadata"])["site"] == "Dublin & Cork"


@pytest.mark.parametrize("name", ["+SUM(1,2)", "-1+2", "@SUM(1,2)", "  =1+1", "\t=1+1", "\n=1+1"])
def test_formula_protection_handles_whitespace(name, tmp_path):
    destination = tmp_path / "protected.csv"
    write_csv([{"sample_name": name, "st": 12}], destination)
    with destination.open(encoding="utf-8-sig", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["sample_name"] == "'" + name
    assert row["st"] == "12"


def test_json_preserves_exact_text_and_provenance(records, tmp_path):
    destination = tmp_path / "results.json"
    export_results(records, destination)
    document = json.loads(destination.read_text())
    assert document["format_version"] == 1
    assert document["samples"][0]["sample_name"] == records[0]["name"]
    assert document["samples"][0]["scheme_digest"] == "abc123"
    assert document["samples"][1]["job_status"] == "queued"


def test_html_escapes_user_content_and_is_self_contained(records, tmp_path):
    destination = tmp_path / "report.html"
    write_html(records, destination)
    report = destination.read_text()
    assert "<script>" not in report
    assert "<img src=" not in report
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in report
    assert "Dublin &amp; Cork" in report
    assert "Pending isolate" in report
    assert "abc123" in report
    assert "@media print" in report
    assert "does not establish an outbreak" in report
    assert "http://" not in report and "https://" not in report


def test_failed_atomic_replace_preserves_existing_export_and_cleans_temp(tmp_path):
    destination = tmp_path / "existing.json"
    destination.write_text("previous data")
    with patch("wmlstudio.export.os.replace", side_effect=OSError("simulated full drive")):
        with pytest.raises(OSError, match="full drive"):
            write_json([{"sample_name": "new"}], destination)
    assert destination.read_text() == "previous data"
    assert list(tmp_path.iterdir()) == [destination]


def test_invalid_result_never_overwrites_existing_export(tmp_path):
    destination = tmp_path / "existing.json"
    destination.write_text("previous data")
    with pytest.raises(ValueError):
        write_json([{"bad": float("nan")}], destination)
    assert destination.read_text() == "previous data"


def test_unknown_format_is_rejected_without_creating_file(tmp_path):
    destination = tmp_path / "output.fasta"
    with pytest.raises(ValueError, match="export format"):
        export_results([], destination)
    assert not Path(destination).exists()


def test_distance_export_identifies_pairs_and_preserves_denominators(tmp_path):
    destination = tmp_path / "distances.json"
    pairs = [{"source": "a", "target": "b", "distance": None, "shared_loci": 0,
              "total_loci": 7, "comparable": False, "reason": "No shared calls"}]
    write_distances(pairs, destination, min_overlap=0.8)
    document = json.loads(destination.read_text())
    assert document["pairs"] == pairs
    assert "samples" not in document
    assert document["min_overlap"] == 0.8
    assert document["metric"] == "allele differences"
    assert "union of profile loci" in document["missing_policy"]
    assert "scheme fingerprints" in document["scheme_requirement"]
