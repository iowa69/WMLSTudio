import csv
from pathlib import Path

import pytest

from wmlstudio import metadata_ingest
from wmlstudio.organism_assignments import (
    ASSIGNMENT_BASIS,
    ASSIGNMENT_CONFIDENCE,
    TEMPLATE_COLUMNS,
    apply_assignments,
    read_assignments,
    write_template,
)
from wmlstudio.project import Project
from wmlstudio.storage import QUARANTINE_ROOT, import_samples

GENERA = [("Klebsiella", "pneumoniae"), ("Escherichia", "coli"), ("Staphylococcus", "aureus"),
          ("Enterococcus", "faecium"), ("Acinetobacter", "baumannii"), ("Pseudomonas", "aeruginosa"),
          ("Salmonella", "enterica"), ("Streptococcus", "pneumoniae"), ("Listeria", "monocytogenes"),
          ("Neisseria", "meningitidis"), ("Haemophilus", "influenzae"), ("Campylobacter", "jejuni"),
          ("Clostridioides", "difficile"), ("Enterobacter", "cloacae"), ("Serratia", "marcescens"),
          ("Citrobacter", "freundii"), ("Proteus", "mirabilis"), ("Morganella", "morganii"),
          ("Providencia", "stuartii"), ("Aeromonas", "hydrophila")]


def sequence(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">contig\n{text}\n")
    return path


@pytest.fixture
def project(tmp_path):
    with Project(tmp_path / "study.wmlstudio") as project:
        yield project


@pytest.fixture
def cohort(project, tmp_path):
    """Twenty unfiled isolates, as a flat drop of mixed files would arrive."""
    assignments = []
    for index in range(len(GENERA)):
        path = sequence(tmp_path / "inbox" / f"isolate{index:02d}.fasta", "ACGT" * (index + 1))
        assignments.append({"path": path, "typing_mode": "auto", "quarantine": "awaiting_identification",
                            "organism_evidence": {"status": "quarantined",
                                                  "quarantine_reason": "awaiting_identification",
                                                  "basis": "none", "confidence": "unresolved"}})
    identifiers = import_samples(project, assignments, tmp_path / "managed")
    return identifiers, tmp_path / "managed"


def write_table(path, rows, columns=("sample_id", "genus", "species"), delimiter=","):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(columns)
        writer.writerows(rows)
    return path


def test_a_spreadsheet_round_trip_files_twenty_isolates_into_twenty_folders(project, cohort, tmp_path):
    identifiers, root = cohort
    template = write_template(tmp_path / "organisms.csv", [
        {"sample_id": sid, "sample_name": project.get_sample(sid)["name"]} for sid in identifiers])
    assert template.read_bytes().startswith(b"\xef\xbb\xbf")
    with template.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == list(TEMPLATE_COLUMNS)
    for row, (genus, species) in zip(rows, GENERA, strict=True):
        row["genus"], row["species"] = genus, species
    edited = write_table(tmp_path / "edited.csv",
                         [[row["sample_id"], row["genus"], row["species"]] for row in rows])
    parsed, problems = read_assignments(edited, samples=project.samples())
    assert problems == [] and len(parsed) == 20
    report = apply_assignments(project, parsed, storage_root=root)
    assert report["assigned"] == 20 and report["quarantined"] == 0
    assert len(report["filing"]["moved"]) == 20
    folders = {Path(sample["input_path"]).relative_to(root).parts[:2] for sample in project.samples()}
    assert folders == {(genus, species) for genus, species in GENERA}
    assert not (root / QUARANTINE_ROOT / "Awaiting_identification").exists()
    for sample in project.samples():
        assert sample["id"] in Path(sample["input_path"]).parts
        assert Path(sample["metadata"]["workflow"]["source_path"]).exists()


def test_a_spreadsheet_assignment_is_never_shown_as_genomic_evidence(project, cohort, tmp_path):
    identifiers, root = cohort
    table = write_table(tmp_path / "one.csv", [[identifiers[0], "Klebsiella", "pneumoniae"]])
    parsed, problems = read_assignments(table, samples=project.samples())
    apply_assignments(project, parsed, storage_root=root)
    evidence = project.get_sample(identifiers[0])["metadata"]["organism_evidence"]
    assert evidence["basis"] == ASSIGNMENT_BASIS == "csv_import"
    assert evidence["confidence"] == ASSIGNMENT_CONFIDENCE == "unresolved"
    assert evidence["status"] == "confirmed" and evidence["confirmed_by"] == "user"
    assert evidence["accepted"] == {"genus": "Klebsiella", "species": "pneumoniae"}


def test_an_ambiguous_sample_name_is_reported_and_nothing_is_written(project, tmp_path):
    for index in range(2):
        import_samples(project, [{"path": sequence(tmp_path / f"dir{index}" / "isolate.fasta",
                                                   "ACGT" * (index + 1)), "typing_mode": "auto"}],
                       tmp_path / "managed")
    before = [sample["input_path"] for sample in project.samples()]
    table = write_table(tmp_path / "names.csv", [["isolate", "Klebsiella", "pneumoniae"]],
                        columns=("sample_name", "genus", "species"))
    parsed, problems = read_assignments(table, samples=project.samples())
    assert parsed == [] and len(problems) == 1
    assert "ambiguous" in problems[0] and "use sample_id" in problems[0]
    assert [sample["input_path"] for sample in project.samples()] == before


def test_an_unknown_identifier_is_a_problem_rather_than_a_guessed_match(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "mixed.csv", [
        [identifiers[0], "Klebsiella", "pneumoniae"],
        ["not-a-real-sample-id", "Escherichia", "coli"]])
    parsed, problems = read_assignments(table, samples=project.samples())
    assert [row["sample_id"] for row in parsed] == [identifiers[0]]
    assert any("never replaced by a guessed name" in problem for problem in problems)


def test_a_row_without_a_genus_or_a_quarantine_reason_is_refused(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "blank.csv", [[identifiers[0], "", ""]])
    parsed, problems = read_assignments(table, samples=project.samples())
    assert parsed == []
    assert "give a genus, or a quarantine reason" in problems[0]


def test_a_row_cannot_claim_an_organism_and_a_quarantine_at_the_same_time(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "both.csv", [[identifiers[0], "Klebsiella", "", "low_confidence"]],
                        columns=("sample_id", "genus", "species", "quarantine"))
    parsed, problems = read_assignments(table, samples=project.samples())
    assert parsed == [] and "cannot be both assigned" in problems[0]


def test_an_unknown_quarantine_reason_is_refused_by_name(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "bad.csv", [[identifiers[0], "", "", "probably_novel_species"]],
                        columns=("sample_id", "genus", "species", "quarantine"))
    parsed, problems = read_assignments(table, samples=project.samples())
    assert parsed == [] and "quarantine must be one of" in problems[0]


def test_a_quarantine_row_moves_a_filed_isolate_back_to_needs_review(project, cohort, tmp_path):
    identifiers, root = cohort
    first = write_table(tmp_path / "a.csv", [[identifiers[0], "Klebsiella", "pneumoniae"]])
    apply_assignments(project, read_assignments(first, samples=project.samples())[0], storage_root=root)
    assert (root / "Klebsiella" / "pneumoniae").is_dir()
    second = write_table(tmp_path / "b.csv", [[identifiers[0], "", "", "user_deferred"]],
                         columns=("sample_id", "genus", "species", "quarantine"))
    report = apply_assignments(project, read_assignments(second, samples=project.samples())[0],
                               storage_root=root)
    assert report["quarantined"] == 1 and report["assigned"] == 0
    filed = Path(project.get_sample(identifiers[0])["input_path"]).relative_to(root)
    assert filed.parts[:2] == (QUARANTINE_ROOT, "User_deferred")
    assert not (root / "Klebsiella").exists()


def test_the_same_isolate_cannot_be_assigned_twice_in_one_table(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "dupe.csv", [[identifiers[0], "Klebsiella", "pneumoniae"],
                                                [identifiers[0], "Escherichia", "coli"]])
    parsed, problems = read_assignments(table, samples=project.samples())
    assert len(parsed) == 1 and "already assigned on row 2" in problems[0]


def test_a_tab_separated_table_and_a_file_column_are_both_accepted(project, tmp_path):
    first = sequence(tmp_path / "inbox" / "a.fasta", "ACGT")
    second = sequence(tmp_path / "inbox" / "b.fasta", "AACC")
    table = write_table(tmp_path / "drop.tsv", [[first.name, "Klebsiella", "pneumoniae"],
                                                [str(second), "Escherichia", "coli"]],
                        columns=("file", "genus", "species"), delimiter="\t")
    parsed, problems = read_assignments(table, paths=[first, second])
    assert problems == []
    assert [Path(row["path"]).name for row in parsed] == ["a.fasta", "b.fasta"]
    assert all(row["sample_id"] is None for row in parsed)
    with pytest.raises(ValueError, match="import review"):
        apply_assignments(project, parsed)


def test_two_identity_columns_are_refused_rather_than_silently_ranked(project, cohort, tmp_path):
    identifiers, _ = cohort
    table = write_table(tmp_path / "two.csv", [[identifiers[0], "isolate00", "Klebsiella"]],
                        columns=("sample_id", "sample_name", "genus"))
    with pytest.raises(ValueError, match="exactly one identity column"):
        read_assignments(table, samples=project.samples())


def test_a_table_with_no_identity_column_is_refused(tmp_path):
    table = write_table(tmp_path / "none.csv", [["Klebsiella", "pneumoniae"]],
                        columns=("genus", "species"))
    with pytest.raises(ValueError, match="sample_id, sample_name or file"):
        read_assignments(table)


def test_organism_columns_are_still_refused_by_the_epidemiology_import():
    for field in ("genus", "species", "organism", "workflow"):
        with pytest.raises(ValueError, match="Protected field"):
            metadata_ingest.annotation_field(field)
    assert {"organism", "genus", "species"} <= metadata_ingest.PROTECTED_FIELDS
