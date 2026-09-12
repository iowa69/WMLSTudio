import bz2
import gzip
import hashlib
import re

import pytest

from wmlstudio.sequence import (
    AnalysisCancelled,
    SequenceError,
    inspect_sequence,
    iter_sequences,
    read_pair_identity,
    sample_name,
    validate_read_pair,
)


@pytest.mark.parametrize("compression", [None, "gzip", "bzip2"])
def test_fasta_qc_and_hash(tmp_path, compression):
    raw = b">long description\nACGTACGTNN\n>short\ngcry\n"
    encoded = gzip.compress(raw) if compression == "gzip" else (
        bz2.compress(raw) if compression == "bzip2" else raw
    )
    path = tmp_path / "sample.data"
    path.write_bytes(encoded)
    result = inspect_sequence(path)
    assert result["kind"] == "fasta"
    assert result["input_sha256"] == hashlib.sha256(encoded).hexdigest()
    qc = result["qc"]
    assert qc["records"] == 2
    assert qc["total_bases"] == 14
    assert qc["n50"] == 10
    assert qc["gc_percent"] == 60
    assert qc["ambiguous_bases"] == 4
    assert qc["n_bases"] == 2
    assert qc["n_percent"] == pytest.approx(100 * 2 / 14)
    assert qc["min_length"] == 4
    assert qc["max_length"] == 10
    assert qc["complete_file"] is True
    assert qc["sampled"] is False
    assert qc["mean_quality"] is None


def test_wrapped_fastq_and_phred_qc(tmp_path):
    path = tmp_path / "reads.fq"
    path.write_text("@r1/1\nACG\nTN\n+r1/1\n!5\n?II\n@r2/1\nGC\n+\n55\n")
    result = inspect_sequence(path)
    assert result["kind"] == "fastq"
    qc = result["qc"]
    assert qc["records"] == 2
    assert qc["total_bases"] == 7
    assert qc["mean_quality"] == pytest.approx(170 / 7)
    assert qc["q20_percent"] == pytest.approx(600 / 7)
    assert qc["q30_percent"] == pytest.approx(300 / 7)
    assert qc["n50"] is None
    assert qc["quality_encoding"] == "Phred+33 (assumed)"


def test_fastq_cap_does_not_claim_tail_validation(tmp_path):
    path = tmp_path / "reads.fastq"
    raw = b"@r1\nACGT\n+\nIIII\n@malformed_tail\nINVALID\n"
    path.write_bytes(raw)
    result = inspect_sequence(path, max_reads=1)
    assert result["qc"]["records"] == 1
    assert result["qc"]["sampled"] is True
    assert result["qc"]["complete_file"] is False
    assert result["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert any("not validated" in note for note in result["notes"])
    with pytest.raises(SequenceError):
        inspect_sequence(path, max_reads=2)


def test_fastq_exactly_at_cap_is_complete(tmp_path):
    path = tmp_path / "reads.fq"
    path.write_text("@r1\nACGT\n+\nIIII\n")
    result = inspect_sequence(path, max_reads=1)
    assert result["qc"]["complete_file"] is True
    assert result["qc"]["sampled"] is False


@pytest.mark.parametrize("raw,match", [
    ("", "empty"),
    ("ACGT\n", "expected"),
    (">\nACGT\n", "nonempty"),
    (">a\n>b\nACGT\n", "empty sequence"),
    (">a\nAC GTX\n", "invalid DNA"),
    (">a\nACGU\n", "invalid DNA"),
    (">a\nACGT\n>a same_id\nACGT\n", "duplicate"),
    ("@r\nACGT\n", "missing '+'"),
    ("@r\nACGT\n+\nIII", "truncated"),
    ("@r\nACGT\n+\nIIIII\n", "lengths differ"),
    ("@r\nACGT\n+different\nIIII\n", "does not match"),
    ("@r\nACGT\n+\nI I!\n", "ASCII"),
    ("@r\n\n+\n\n", "blank"),
    ("@r\nACGT\n+\nIIII\nBADHEADER\n", "nonempty"),
])
def test_malformed_inputs_are_rejected(tmp_path, raw, match):
    path = tmp_path / "bad.sequence"
    path.write_text(raw)
    with pytest.raises(SequenceError, match=re.escape(match)):
        inspect_sequence(path)


@pytest.mark.parametrize("compress", [gzip.compress, bz2.compress])
def test_truncated_compressed_input_rejected(tmp_path, compress):
    path = tmp_path / "truncated.data"
    path.write_bytes(compress(b">record\nACGT\n")[:-5])
    with pytest.raises(SequenceError):
        inspect_sequence(path)


def test_mislabelled_compression_rejected(tmp_path):
    path = tmp_path / "plain.fasta.gz"
    path.write_text(">a\nACGT\n")
    with pytest.raises(SequenceError, match="compression header"):
        inspect_sequence(path)


def test_all_ambiguous_sequence_has_no_gc_estimate(tmp_path):
    path = tmp_path / "ambiguous.fa"
    path.write_text(">a\nNNRY\n")
    assert inspect_sequence(path)["qc"]["gc_percent"] is None


@pytest.mark.parametrize("max_reads", [0, -1, 1.5, True, None])
def test_read_cap_must_be_positive_integer(tmp_path, max_reads):
    with pytest.raises(ValueError, match="positive integer"):
        inspect_sequence(tmp_path / "irrelevant", max_reads=max_reads)


def test_cancellation(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text(">r\nACGT\n")
    with pytest.raises(AnalysisCancelled):
        inspect_sequence(path, cancelled=lambda: True)


@pytest.mark.parametrize("header,expected", [
    ("read/1", ("read", 1)),
    ("@read/2 description", ("read", 2)),
    ("instrument:lane:x:y 1:N:0:AGTC", ("instrument:lane:x:y", 1)),
    ("instrument:lane:x:y 2:Y:0:AGTC", ("instrument:lane:x:y", 2)),
    ("sample_R1", ("sample_R1", None)),
    ("read/12", ("read/12", None)),
])
def test_mate_identification_never_guesses_from_filenames(header, expected):
    assert read_pair_identity(header) == expected


def test_conflicting_mate_header_rejected():
    with pytest.raises(SequenceError, match="Conflicting"):
        read_pair_identity("read/1 2:N:0:ACGT")


def test_pair_validation_full_and_sampled(tmp_path):
    first, second = tmp_path / "a.fq", tmp_path / "b.fq"
    first.write_text("@a/1\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n")
    second.write_text("@a/2\nTGCA\n+\nIIII\n@b/2\nTGCA\n+\nIIII\n")
    assert validate_read_pair(first, second)["status"] == "verified"
    prefix = validate_read_pair(first, second, max_reads=1)
    assert prefix["sampled"] is True
    assert prefix["status"] == "prefix_verified"
    second.write_text("@a/2\nTGCA\n+\nIIII\n@c/2\nTGCA\n+\nIIII\n")
    with pytest.raises(SequenceError, match="identifiers differ at record 2"):
        validate_read_pair(first, second)


def test_pair_reversed_and_count_mismatch(tmp_path):
    first, second = tmp_path / "a.fq", tmp_path / "b.fq"
    first.write_text("@a/1\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n")
    second.write_text("@a/2\nTGCA\n+\nIIII\n")
    with pytest.raises(SequenceError, match="different numbers"):
        validate_read_pair(first, second, max_reads=1)
    with pytest.raises(SequenceError, match="reversed"):
        validate_read_pair(second, first)


def test_unmarked_pairs_not_asserted_verified(tmp_path):
    first, second = tmp_path / "a.fq", tmp_path / "b.fq"
    first.write_text("@a\nACGT\n+\nIIII\n")
    second.write_text("@a\nTGCA\n+\nIIII\n")
    assert validate_read_pair(first, second)["status"] == "unmarked"


def test_iter_records_and_sample_name(tmp_path):
    path = tmp_path / "name.FASTA.GZ"
    path.write_bytes(gzip.compress(b">contig description\naCgT\n"))
    record = list(iter_sequences(path))[0]
    assert record.name == "contig description"
    assert record.identifier == "contig"
    assert record.sequence == "ACGT"
    assert sample_name(path) == "name"
