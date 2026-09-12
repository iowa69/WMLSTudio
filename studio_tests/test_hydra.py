import hashlib
import json

import pytest

from wmlstudio.hydra import HydraImportError, load_hydra_report


def payload():
    # Field names mirror SampleResult.as_dict()/Hit.as_row() at HYDRA 1.4.0.
    hit = {"sample": "KP001", "database": "ncbi", "gene": "blaKPC-2",
           "element_type": "AMR", "element_subtype": "AMR", "class": "BETA-LACTAM",
           "method": "EXACTP", "resolution": "COMPLETE", "identity_pct": 100.0,
           "coverage_pct": 100.0, "depth": None, "allele_fraction": None,
           "primary": True, "note": ""}
    return {"hydra_version": "1.4.0", "command": "hydra run -a KP001.fasta --format json",
            "databases": ["ncbi", "card"], "parameters": {"min_identity": 90},
            "samples": [{"sample": "KP001", "input_type": "assembly", "inputs": ["KP001.fasta"],
                         "species": {"name": "Klebsiella pneumoniae", "confidence": "strong"},
                         "mlst": {"scheme": "klebsiella", "sequence_type": "258",
                                  "alleles": {"gapA": "3"}, "source": "assembly"},
                         "typing": [], "scores": {}, "qc": {}, "warnings": [],
                         "runtime_seconds": 16.3,
                         "hits": [hit, {**hit, "database": "card", "primary": False}]}]}


def write_report(tmp_path, report):
    path = tmp_path / "hydra.json"
    path.write_text(json.dumps(report), encoding="utf-8")
    return path


def test_import_preserves_upstream_calls_and_provenance(tmp_path):
    original = payload()
    path = write_report(tmp_path, original)
    result = load_hydra_report(path)
    assert result["samples"][0]["summary"]["amr_genes"] == 1
    assert result["samples"][0]["summary"]["total_hits"] == 2
    assert result["samples"][0]["hits"] == original["samples"][0]["hits"]
    assert result["samples"][0]["mlst"] == original["samples"][0]["mlst"]
    assert result["samples"][0]["runtime_seconds"] == 16.3
    assert result["import_provenance"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert result["import_provenance"]["analysis_origin"] == "external_hydra_report"
    assert result["import_warnings"] == []


def test_unknown_version_has_warning_and_preserves_extra_fields(tmp_path):
    report = payload()
    report["hydra_version"] = "2.0.0"
    report["future_evidence"] = {"note": "retained"}
    result = load_hydra_report(write_report(tmp_path, report))
    assert "not been validated" in result["import_warnings"][0]
    assert result["future_evidence"] == report["future_evidence"]


def test_missing_primary_is_not_counted_as_confirmed(tmp_path):
    report = payload()
    del report["samples"][0]["hits"][0]["primary"]
    result = load_hydra_report(write_report(tmp_path, report))
    assert result["samples"][0]["summary"]["amr_genes"] == 0
    assert result["import_warnings"]


def test_catalogued_heteroresistance_is_distinguished_from_other_variants(tmp_path):
    report = payload()
    report["samples"][0]["hits"][0].update(
        method="POINTR", resolution="POINT", allele_fraction=0.21,
        note="G2576T; HETERORESISTANT", depth=100)
    report["samples"][0]["hits"][1].update(
        method="VARIANTR", resolution="POINT", primary=True, allele_fraction=0.3)
    summary = load_hydra_report(write_report(tmp_path, report))["samples"][0]["summary"]
    assert summary["point_mutations"] == 1
    assert summary["heteroresistant_sites"] == 1


@pytest.mark.parametrize("field,value", [
    ("primary", "false"), ("identity_pct", 101), ("coverage_pct", -1),
    ("allele_fraction", 2), ("depth", -2), ("sample", "OTHER"), ("gene", []),
    ("note", {}), ("class", []),
])
def test_invalid_evidence_fails_whole_import(tmp_path, field, value):
    report = payload()
    report["samples"][0]["hits"][0][field] = value
    with pytest.raises(HydraImportError):
        load_hydra_report(write_report(tmp_path, report))


@pytest.mark.parametrize("report", [[], {"samples": []}, {"hydra_version": "1.4.0", "samples": {}},
                                   {"hydra_version": "1.4.0", "samples": [{"sample": "a"}]}])
def test_unsupported_shape(tmp_path, report):
    with pytest.raises(HydraImportError):
        load_hydra_report(write_report(tmp_path, report))


def test_duplicate_samples_rejected(tmp_path):
    report = payload()
    report["samples"].append(report["samples"][0])
    with pytest.raises(HydraImportError, match="Duplicate HYDRA sample"):
        load_hydra_report(write_report(tmp_path, report))


@pytest.mark.parametrize("text", ['{"hydra_version":null,"hydra_version":"1.4.0","samples":[]}',
                                  '{"hydra_version":"1.4.0","samples":[],"bad":NaN}',
                                  '{"hydra_version":"1.4.0","samples":[],"bad":1e999}', "<html>"])
def test_ambiguous_or_non_json_rejected(tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text)
    with pytest.raises(HydraImportError):
        load_hydra_report(path)


def test_source_paths_are_not_opened(tmp_path):
    report = payload()
    report["samples"][0]["inputs"] = ["https://not-a-real-location/input.fasta", "../absent.fasta"]
    assert load_hydra_report(write_report(tmp_path, report))["samples"][0]["inputs"]
