import gzip
import hashlib
import json
import random

import pytest

from wmlstudio.sequence import AnalysisCancelled, SequenceError
from wmlstudio.typing import SchemeError, call_assembly, load_scheme, reverse_complement

ARC1 = "AACCGTACGTTAG"
ARC2 = "AACCGTTCGTTAG"
GYR1 = "TTGGCATACCTGA"
GYR2 = "TTGGCATTCCTGA"


@pytest.fixture
def schema_path(tmp_path):
    directory = tmp_path / "local_scheme"
    directory.mkdir()
    (directory / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_2\n{ARC2}\n")
    (directory / "gyrB.fasta").write_text(f">gyrB_1\n{GYR1}\n>gyrB_2\n{GYR2}\n")
    (directory / "profiles.tsv").write_text("ST\tarcA\tgyrB\tclonal_complex\n1\t1\t1\tA\n2\t2\t1\tB\n")
    (directory / "scheme.json").write_text(json.dumps({"name": "Known truth", "source": "synthetic"}))
    return directory


def assembly(tmp_path, sequence, name="assembly.fa"):
    path = tmp_path / name
    path.write_text(f">contig1\n{sequence}\n")
    return path


def test_schema_contract(schema_path):
    scheme = load_scheme(schema_path)
    assert scheme.name == "Known truth"
    assert scheme.loci == ("arcA", "gyrB")
    assert scheme.allele_count == 4
    assert scheme.profile_count == 2
    assert scheme.profiles[("1", "1")] == ("1",)
    assert len(scheme.digest) == 64
    assert scheme.scheme_digest == scheme.digest


def test_known_st_forward_and_reverse_with_exact_coordinates(tmp_path, schema_path):
    path = assembly(tmp_path, "NNN" + ARC1 + "NNNNN" + reverse_complement(GYR1) + "NN")
    result = call_assembly(path, load_scheme(schema_path))
    assert result["status"] == "complete"
    assert result["st"] == "1"
    assert result["alleles"] == {"arcA": "1", "gyrB": "1"}
    assert result["input_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    arc, gyr = result["calls"]
    assert arc["hits"][0]["start"] == 4
    assert arc["hits"][0]["end"] == 16
    assert arc["hits"][0]["strand"] == "+"
    assert gyr["hits"][0]["start"] == 22
    assert gyr["hits"][0]["end"] == 34
    assert gyr["hits"][0]["strand"] == "-"


def test_known_st_in_compressed_fasta(tmp_path, schema_path):
    path = tmp_path / "assembly.fa.gz"
    path.write_bytes(gzip.compress(f">a\n{ARC2}\n>b\n{GYR1}\n".encode()))
    result = call_assembly(path, load_scheme(schema_path))
    assert result["st"] == "2"
    assert result["qc"]["records"] == 2


def test_complete_new_combination_is_novel_profile_only(tmp_path, schema_path):
    result = call_assembly(assembly(tmp_path, ARC1 + "NNN" + GYR2), load_scheme(schema_path))
    assert result["status"] == "novel_profile"
    assert result["st"] is None
    assert result["alleles"] == {"arcA": "1", "gyrB": "2"}


def test_missing_exact_match_is_not_called_novel_allele(tmp_path, schema_path):
    result = call_assembly(assembly(tmp_path, ARC1 + "NNN" + GYR1[:-1] + "C"), load_scheme(schema_path))
    assert result["status"] == "incomplete"
    assert result["st"] is None
    assert result["alleles"]["gyrB"] is None
    assert result["calls"][1]["status"] == "missing"
    assert any("do not establish novel alleles" in note for note in result["notes"])


def test_different_alleles_at_separate_positions_are_mixed(tmp_path, schema_path):
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + ARC2 + "NN" + GYR1), load_scheme(schema_path))
    assert result["status"] == "mixed"
    assert result["st"] is None
    assert result["alleles"]["arcA"] is None
    assert result["calls"][0]["candidates"] == ["1", "2"]


def test_repeated_same_allele_is_ambiguous_copy_number(tmp_path, schema_path):
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["status"] == "ambiguous"
    assert result["alleles"]["arcA"] is None
    assert result["calls"][0]["hit_count"] == 2


def test_identical_sequence_aliases_are_ambiguous(tmp_path, schema_path):
    (schema_path / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_2\n{ARC1}\n")
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["status"] == "ambiguous"
    assert result["st"] is None
    assert result["calls"][0]["status"] == "ambiguous"


def test_contained_reference_alleles_are_ambiguous(tmp_path, schema_path):
    (schema_path / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_2\n{ARC1[2:-2]}\n")
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["calls"][0]["status"] == "ambiguous"


def test_hits_beyond_display_limit_still_prevent_false_st(tmp_path, schema_path):
    result = call_assembly(assembly(tmp_path, (ARC1 + "NN") * 100 + ARC2 + "NN" + GYR1), load_scheme(schema_path))
    arc = result["calls"][0]
    assert arc["hit_count"] == 101
    assert len(arc["hits"]) == 50
    assert arc["evidence_truncated"] is True
    assert arc["candidates"] == ["1", "2"]
    assert arc["status"] == "mixed"
    assert result["st"] is None


def test_does_not_match_across_contigs(tmp_path, schema_path):
    path = tmp_path / "fragmented.fa"
    path.write_text(f">one\n{ARC1[:5]}\n>two\n{ARC1[5:]}\n>three\n{GYR1}\n")
    result = call_assembly(path, load_scheme(schema_path))
    assert result["alleles"]["arcA"] is None
    assert result["status"] == "incomplete"


def test_duplicate_contig_ids_cannot_collapse_into_single_exact_call(tmp_path, schema_path):
    path = tmp_path / "duplicate_contigs.fa"
    path.write_text(f">same first\n{ARC1}\n>same second\n{ARC1}\n>other\n{GYR1}\n")
    with pytest.raises(SequenceError, match="duplicate FASTA identifier 'same'"):
        call_assembly(path, load_scheme(schema_path))


def test_match_across_internal_processing_chunk(tmp_path, schema_path):
    prefix = "N" * (65536 - 5)
    result = call_assembly(assembly(tmp_path, prefix + ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["st"] == "1"
    assert result["calls"][0]["hits"][0]["start"] == len(prefix) + 1
    assert result["calls"][0]["hits"][0]["end"] == len(prefix) + len(ARC1)


def test_large_schema_seed_matches_are_verified_as_full_alleles(tmp_path):
    rng = random.Random(501)
    root = tmp_path / 'large'
    root.mkdir()
    sequences = [''.join(rng.choice('ACGT') for _ in range(101)) for _ in range(31)]
    # Equal central seed, different complete sequence: not shared-locus evidence.
    sequences[1] = 'A' * 35 + sequences[0][35:66] + 'C' * 35
    for index, sequence in enumerate(sequences):
        (root / f'locus{index:02}.fa').write_text(f'>locus{index:02}_1\n{sequence}\n')
    prefix = 'N' * (65536 - 60)
    path = assembly(tmp_path, prefix + reverse_complement(sequences[0]))
    result = call_assembly(path, load_scheme(root))
    assert result['parameters']['index'] == 'seed-verified'
    assert result['alleles']['locus00'] == '1'
    assert result['alleles']['locus01'] is None
    assert sum(value is not None for value in result['alleles'].values()) == 1
    hit = result['calls'][0]['hits'][0]
    assert (hit['start'], hit['end'], hit['strand']) == (len(prefix) + 1, len(prefix) + 101, '-')


def test_empty_profile_directory_supports_cgmlst(tmp_path, schema_path):
    (schema_path / "profiles.tsv").unlink()
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["status"] == "profile_unavailable"
    assert result["st"] is None
    assert result["alleles"] == {"arcA": "1", "gyrB": "1"}


def test_multiple_sts_for_same_profile_not_silently_assigned(tmp_path, schema_path):
    (schema_path / "profiles.tsv").write_text("ST\tarcA\tgyrB\n1\t1\t1\n99\t1\t1\n")
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + GYR1), load_scheme(schema_path))
    assert result["status"] == "ambiguous"
    assert result["st"] is None
    assert any("multiple STs" in note for note in result["notes"])


def test_identical_references_between_loci_are_ambiguous(tmp_path, schema_path):
    (schema_path / "gyrB.fasta").write_text(f">gyrB_1\n{ARC1}\n>gyrB_2\n{GYR2}\n")
    result = call_assembly(assembly(tmp_path, ARC1), load_scheme(schema_path))
    assert result["status"] == "ambiguous"
    assert result["alleles"] == {"arcA": None, "gyrB": None}


def test_schema_digest_tracks_reference_profiles_and_metadata(schema_path):
    initial = load_scheme(schema_path).digest
    assert load_scheme(schema_path).digest == initial
    (schema_path / "unrelated.log").write_text("does not affect biological sources")
    assert load_scheme(schema_path).digest == initial
    (schema_path / "scheme.json").write_text('{"name":"changed"}')
    changed_metadata = load_scheme(schema_path).digest
    assert changed_metadata != initial
    (schema_path / "profiles.tsv").write_text("ST\tarcA\tgyrB\n3\t1\t1\n")
    changed_profile = load_scheme(schema_path).digest
    assert changed_profile != changed_metadata
    (schema_path / "arcA.tfa").write_text(f">arcA_1\n{ARC2}\n>arcA_2\n{ARC1}\n")
    assert load_scheme(schema_path).digest != changed_profile


@pytest.mark.parametrize("profile,match", [
    ("ST\tarcA\n1\t1\n", "lacks loci"),
    ("ST\tarcA\tgyrB\n1\t1\n", "row width"),
    ("ST\tarcA\tgyrB\n1\t1\t1\n1\t2\t1\n", "conflicting"),
    ("ST\tarcA\tarcA\tgyrB\n1\t1\t1\t1\n", "duplicate"),
])
def test_invalid_profiles_rejected(schema_path, profile, match):
    (schema_path / "profiles.tsv").write_text(profile)
    with pytest.raises(SchemeError, match=match):
        load_scheme(schema_path)


def test_incomplete_profiles_excluded(schema_path):
    (schema_path / "profiles.tsv").write_text("ST\tarcA\tgyrB\n1\t1\t1\n2\t0\t1\n")
    scheme = load_scheme(schema_path)
    assert scheme.profile_count == 1
    assert any("incomplete" in note for note in scheme.notes)


def test_ambiguous_reference_dna_excluded_from_matching(tmp_path, schema_path):
    (schema_path / "arcA.tfa").write_text(">arcA_1\nACGTN\n")
    scheme = load_scheme(schema_path)
    assert any("ambiguous DNA" in note for note in scheme.notes)
    result = call_assembly(assembly(tmp_path, "ACGTANNNNACGTNNNN" + GYR1), scheme)
    assert result["status"] == "incomplete"
    assert result["alleles"]["arcA"] is None


def test_unavailable_reference_profiles_are_visible_and_never_assigned(tmp_path, schema_path):
    (schema_path / "profiles.tsv").write_text("ST\tarcA\tgyrB\n1\t9\t1\n2\t1\t1\n")
    scheme = load_scheme(schema_path)
    assert any("arcA_9" in note and "cannot be assigned" in note for note in scheme.notes)
    result = call_assembly(assembly(tmp_path, ARC1 + "NN" + GYR1), scheme)
    assert result["st"] == "2"


def test_duplicate_reference_identifiers_rejected(schema_path):
    (schema_path / "arcA.tfa").write_text(f">arcA_1\n{ARC1}\n>arcA_1\n{ARC2}\n")
    with pytest.raises(SchemeError, match="duplicate"):
        load_scheme(schema_path)


def test_fastq_is_not_typed_as_assembly(tmp_path, schema_path):
    path = tmp_path / "reads.fastq"
    path.write_text(f"@r\n{ARC1}\n+\n{'I' * len(ARC1)}\n")
    with pytest.raises(SequenceError, match="FASTA assembly"):
        call_assembly(path, load_scheme(schema_path))


def test_cancellation_and_progress(tmp_path, schema_path):
    with pytest.raises(AnalysisCancelled):
        load_scheme(schema_path, cancelled=lambda: True)
    scheme = load_scheme(schema_path)
    path = assembly(tmp_path, ARC1 + "NN" + GYR1)
    with pytest.raises(AnalysisCancelled):
        call_assembly(path, scheme, cancelled=lambda: True)
    updates = []
    result = call_assembly(path, scheme, progress=lambda *args: updates.append(args))
    assert result["st"] == "1"
    assert any("Indexing" in message for _, _, message in updates)
    assert any("Scanned" in message for _, _, message in updates)


def test_palindromic_allele_is_one_hit(tmp_path):
    directory = tmp_path / "palindrome"
    directory.mkdir()
    (directory / "pal.fa").write_text(">pal_1\nACGTACGT\n")
    result = call_assembly(assembly(tmp_path, "NNACGTACGTNN"), load_scheme(directory))
    assert result["alleles"] == {"pal": "1"}
    assert result["calls"][0]["hit_count"] == 1
    assert result["calls"][0]["hits"][0]["strand"] == "both"
