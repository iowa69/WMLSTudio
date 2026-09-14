import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from wmlstudio.export import (
    REPORT_PRESETS,
    export_results,
    review_report_html,
    write_csv,
    write_distances,
    write_html,
    write_json,
    write_tsv,
)
from wmlstudio.simple_report import SIMPLE_REPORT_PRESET

# Every key the "Customize sections…" dialog renders with options[key]; a preset
# missing one of them raises KeyError the moment a user opens that dialog.
DIALOG_KEYS = ("title", "investigation", "qc", "amr", "virulence", "plasmid_hypotheses",
               "drug_associations", "graph", "graph_jpeg", "provenance")


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


def test_every_report_preset_carries_every_option_the_section_dialog_renders():
    for key, preset in REPORT_PRESETS.items():
        assert [name for name in DIALOG_KEYS if name in preset] == list(DIALOG_KEYS), key
    # The simple preset is defined once, and appended last so no saved template
    # and no stored combo position can be repointed by adding a preset.
    assert REPORT_PRESETS["simple"] == SIMPLE_REPORT_PRESET
    assert list(REPORT_PRESETS)[-1] == "simple"
    assert not any(preset.get("layout") for key, preset in REPORT_PRESETS.items() if key != "simple")


def test_the_one_page_layout_is_reached_through_the_ordinary_report_writer(records):
    simple = review_report_html(records, selected_ids={"stable-id"}, options=REPORT_PRESETS["simple"])
    assert "<title>Simple outbreak summary</title>" in simple
    assert "What this report does not tell you" in simple
    assert "No comparison has been built for these isolates" in simple
    # The detailed shape is untouched: it still loops one section per isolate.
    detailed = review_report_html(records, selected_ids={"stable-id"}, options=REPORT_PRESETS["cohort"])
    assert "At a glance" in detailed and "What this report does not tell you" not in detailed


def test_the_report_writer_threads_the_image_type_and_the_scope_sentence(records):
    picture = b"\xff\xd8\xff\xe0 pretend image bytes"
    note = "No isolates were chosen for this report, so it covers all 2 isolates in the project."
    threaded = review_report_html(records, selected_ids={"stable-id"}, graph_png=picture,
                                  graph_mime="image/jpeg", scope_note=note, scope_implicit=True)
    assert "data:image/jpeg;base64," in threaded
    assert note in threaded
    assert "This scope was not chosen for this report." in threaded
    # Every existing caller keeps PNG, the counted scope line, and no extra notice.
    default = review_report_html(records, selected_ids={"stable-id"}, graph_png=picture)
    assert "data:image/png;base64," in default
    assert "1 explicitly selected isolate(s)" in default
    assert "This scope was not chosen" not in default


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
