"""Streaming nucleotide input, explicit QC sampling, and read-pair validation.

Sequence files are never modified. FASTQ qualities are interpreted as Phred+33;
the encoding cannot reliably be inferred from the observed character range.
"""

from __future__ import annotations

import bz2
import gzip
import hashlib
import io
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, TextIO

CancelCallback = Callable[[], bool] | None
DNA = frozenset("ACGTRYSWKMBDHVN")


class SequenceError(ValueError):
    """Input is empty, malformed, unsupported, or changed during analysis."""


class AnalysisCancelled(Exception):
    """A user cancellation was observed at a safe boundary."""


def check_cancelled(cancelled: CancelCallback) -> None:
    if cancelled is not None and cancelled():
        raise AnalysisCancelled("Analysis cancelled.")


def sample_name(path: str | Path) -> str:
    name = Path(path).name
    if name.lower().endswith((".gz", ".bz2")):
        name = name.rsplit(".", 1)[0]
    if name.lower().endswith((".fasta", ".fna", ".fa", ".fas", ".tfa", ".fastq", ".fq")):
        name = name.rsplit(".", 1)[0]
    return name


def file_sha256(path: str | Path, cancelled: CancelCallback = None) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            check_cancelled(cancelled)
            digest.update(chunk)
    check_cancelled(cancelled)
    return digest.hexdigest()


def file_signature(path: str | Path) -> tuple[int, int, int]:
    stat = Path(path).stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    name: str
    sequence: str
    quality: str | None = None

    @property
    def identifier(self) -> str:
        return self.name.split()[0]


class SequenceReader:
    """One-record-at-a-time reader, with one header of look-ahead.

    ``complete`` is true when EOF has been observed after the last yielded record.
    A bounded FASTQ scan does not validate records beyond the requested prefix.
    Both four-line and wrapped FASTQ records are supported.
    """

    def __init__(self, path: str | Path, cancelled: CancelCallback = None):
        self.path = Path(path)
        self.cancelled = cancelled
        self.kind = ""
        self.complete = False
        self._line_number = 0
        self._handle: TextIO | None = None
        self._raw: io.BufferedReader | None = None
        self._pending = ""

    def __enter__(self) -> SequenceReader:
        check_cancelled(self.cancelled)
        try:
            self._raw = self.path.open("rb")
            magic = self._raw.peek(3)[:3]
            if magic.startswith(b"\x1f\x8b"):
                binary = gzip.GzipFile(fileobj=self._raw)
            elif magic.startswith(b"BZh"):
                binary = bz2.BZ2File(self._raw)
            else:
                if self.path.suffix.lower() in {".gz", ".bz2"}:
                    raise SequenceError(f"{self.path.name}: compression header does not match suffix.")
                binary = self._raw
            self._handle = io.TextIOWrapper(binary, encoding="ascii", newline=None)
            self._pending = self._readline()
            while self._pending and not self._pending.strip():
                self._pending = self._readline()
            if not self._pending:
                raise SequenceError(f"{self.path.name}: sequence file is empty.")
            if self._pending.startswith(">"):
                self.kind = "fasta"
            elif self._pending.startswith("@"):
                self.kind = "fastq"
            else:
                raise self._error("expected a FASTA '>' or FASTQ '@' header")
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *exc: object) -> None:
        if self._handle is not None:
            self._handle.close()
        if self._raw is not None:
            self._raw.close()

    def _error(self, message: str) -> SequenceError:
        return SequenceError(f"{self.path.name}, line {self._line_number}: {message}.")

    def _readline(self) -> str:
        check_cancelled(self.cancelled)
        assert self._handle is not None
        try:
            line = self._handle.readline()
        except (UnicodeDecodeError, OSError, EOFError) as error:
            raise self._error(f"cannot decode complete sequence data ({error})") from error
        if line:
            self._line_number += 1
        return line

    def _sequence_line(self, line: str) -> str:
        sequence = line.rstrip("\r\n").upper()
        invalid = set(sequence) - DNA
        if invalid:
            shown = ", ".join(repr(x) for x in sorted(invalid)[:5])
            raise self._error(f"invalid DNA characters: {shown}")
        return sequence

    def __iter__(self) -> Iterator[SequenceRecord]:
        identifiers: set[str] = set()
        while self._pending:
            check_cancelled(self.cancelled)
            header = self._pending.rstrip("\r\n")
            marker = ">" if self.kind == "fasta" else "@"
            if not header.startswith(marker) or not header[1:].strip():
                raise self._error(f"expected a nonempty {marker!r} header")
            name = header[1:].strip()
            identifier = name.split()[0]
            if self.kind == "fasta":
                if identifier in identifiers:
                    raise self._error(f"duplicate FASTA identifier {identifier!r}")
                identifiers.add(identifier)
                parts: list[str] = []
                line = self._readline()
                while line and not line.startswith(">"):
                    if line.rstrip("\r\n"):
                        parts.append(self._sequence_line(line))
                    line = self._readline()
                sequence = "".join(parts)
                if not sequence:
                    raise self._error(f"empty sequence for {identifier!r}")
                self._pending = line
                self.complete = not bool(line)
                yield SequenceRecord(name, sequence)
                continue
            parts = []
            line = self._readline()
            while line and not line.startswith("+"):
                if not line.rstrip("\r\n"):
                    raise self._error("blank FASTQ sequence line")
                parts.append(self._sequence_line(line))
                line = self._readline()
            if not line:
                raise self._error(f"truncated FASTQ record {identifier!r}: missing '+' separator")
            repeated_header = line[1:].rstrip("\r\n")
            if repeated_header and repeated_header != name:
                raise self._error("FASTQ '+' header does not match the sequence header")
            sequence = "".join(parts)
            if not sequence:
                raise self._error(f"empty FASTQ sequence for {identifier!r}")
            qualities: list[str] = []
            quality_length = 0
            while quality_length < len(sequence):
                line = self._readline()
                if not line:
                    raise self._error(f"truncated FASTQ quality for {identifier!r}")
                quality = line.rstrip("\r\n")
                if not quality or any(ord(char) < 33 or ord(char) > 126 for char in quality):
                    raise self._error("FASTQ quality must use ASCII characters 33 through 126")
                quality_length += len(quality)
                if quality_length > len(sequence):
                    raise self._error(f"FASTQ sequence/quality lengths differ for {identifier!r}")
                qualities.append(quality)
            self._pending = self._readline()
            self.complete = not bool(self._pending)
            yield SequenceRecord(name, sequence, "".join(qualities))


def iter_sequences(
    path: str | Path, cancelled: CancelCallback = None
) -> Iterator[SequenceRecord]:
    with SequenceReader(path, cancelled) as reader:
        yield from reader


class QCAccumulator:
    """Streaming base/quality counts plus a length histogram for exact N50."""

    def __init__(self) -> None:
        self.lengths: Counter[int] = Counter()
        self.bases: Counter[str] = Counter()
        self.qualities: Counter[str] = Counter()

    def add(self, record: SequenceRecord) -> None:
        self.lengths[len(record.sequence)] += 1
        self.bases.update(record.sequence)
        if record.quality is not None:
            self.qualities.update(record.quality)

    def result(self, kind: str, complete: bool, limit: int | None) -> dict:
        records = sum(self.lengths.values())
        total = sum(length * count for length, count in self.lengths.items())
        canonical = sum(self.bases[x] for x in "ACGT")
        n50 = None
        if kind == "fasta":
            cumulative = 0
            for length, count in sorted(self.lengths.items(), reverse=True):
                cumulative += length * count
                if cumulative * 2 >= total:
                    n50 = length
                    break
        quality_bases = sum(self.qualities.values())
        return {
            "records": records,
            "total_bases": total,
            "min_length": min(self.lengths, default=0),
            "max_length": max(self.lengths, default=0),
            "mean_length": total / records if records else 0,
            "n50": n50,
            "gc_percent": 100 * (self.bases["G"] + self.bases["C"]) / canonical
            if canonical else None,
            "gc_denominator": "ACGT bases",
            "acgt_bases": canonical,
            "n_bases": self.bases["N"],
            "n_percent": 100 * self.bases["N"] / total if total else 0,
            "ambiguous_bases": total - canonical,
            "ambiguous_percent": 100 * (total - canonical) / total if total else 0,
            "mean_quality": sum((ord(q) - 33) * n for q, n in self.qualities.items())
            / quality_bases if quality_bases else None,
            "q20_percent": 100 * sum(n for q, n in self.qualities.items() if ord(q) - 33 >= 20)
            / quality_bases if quality_bases else None,
            "q30_percent": 100 * sum(n for q, n in self.qualities.items() if ord(q) - 33 >= 30)
            / quality_bases if quality_bases else None,
            "quality_encoding": "Phred+33 (assumed)" if kind == "fastq" else None,
            "sampled": not complete,
            "complete_file": complete,
            "records_limit": limit if kind == "fastq" else None,
        }


def inspect_sequence(
    path: str | Path, max_reads: int = 100000, cancelled: CancelCallback = None
) -> dict:
    """Inspect all FASTA records or at most the first ``max_reads`` FASTQ records.

    Hashing always covers the complete on-disk input, including compression bytes.
    A sampled QC result makes no claim about validity of the unparsed remainder.
    """
    if not isinstance(max_reads, int) or isinstance(max_reads, bool) or max_reads < 1:
        raise ValueError("max_reads must be a positive integer.")
    path = Path(path).resolve()
    signature = file_signature(path)
    qc = QCAccumulator()
    with SequenceReader(path, cancelled) as reader:
        for count, record in enumerate(reader, 1):
            qc.add(record)
            if reader.kind == "fastq" and count >= max_reads:
                break
        kind = reader.kind
        complete = reader.complete
    digest = file_sha256(path, cancelled)
    if file_signature(path) != signature:
        raise SequenceError(f"{path.name}: input changed during analysis; run it again.")
    notes = []
    if kind == "fastq":
        notes.append("Read quality assumes Phred+33 encoding.")
        notes.append("Reads require assembly before this exact assembly typing engine can be used.")
    if not complete:
        notes.append(
            f"QC examined the first {max_reads:,} reads only; subsequent records were not validated. "
            "SHA-256 covers the complete input file."
        )
    return {
        "sample_name": sample_name(path),
        "input_path": str(path),
        "kind": kind,
        "qc": qc.result(kind, complete, max_reads),
        "input_sha256": digest,
        "notes": notes,
    }


def read_pair_identity(header: str) -> tuple[str, int | None]:
    """Extract a read identity and explicit /1,/2 or CASAVA mate; never guess."""
    fields = header.removeprefix("@").split()
    if not fields:
        raise SequenceError("Empty read identifier.")
    identifier = fields[0]
    mate = None
    suffix = re.search(r"/([12])$", identifier)
    if suffix:
        mate = int(suffix.group(1))
        identifier = identifier[:-2]
    if len(fields) > 1 and re.match(r"^[12]:[YN]:\d+:", fields[1]):
        casava_mate = int(fields[1][0])
        if mate is not None and mate != casava_mate:
            raise SequenceError("Conflicting read mate indicators in FASTQ header.")
        mate = casava_mate
    return identifier, mate


def validate_read_pair(
    first: str | Path,
    second: str | Path,
    max_reads: int = 100000,
    cancelled: CancelCallback = None,
) -> dict:
    """Validate corresponding read IDs, explicit mate direction and prefix counts.

    Unmarked matching identifiers are reported as unmarked, never proved paired.
    Matching filenames alone are not evidence of a biological read pair.
    """
    if not isinstance(max_reads, int) or isinstance(max_reads, bool) or max_reads < 1:
        raise ValueError("max_reads must be a positive integer.")
    if Path(first).resolve() == Path(second).resolve():
        raise SequenceError("Read-pair inputs must be different files.")
    count = 0
    explicit = True
    with SequenceReader(first, cancelled) as left, SequenceReader(second, cancelled) as right:
        if left.kind != "fastq" or right.kind != "fastq":
            raise SequenceError("Read-pair validation requires two FASTQ files.")
        left_iter, right_iter = iter(left), iter(right)
        for _ in range(max_reads):
            a, b = next(left_iter, None), next(right_iter, None)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise SequenceError("Read-pair files contain different numbers of records.")
            left_id, left_mate = read_pair_identity(a.name)
            right_id, right_mate = read_pair_identity(b.name)
            if left_id != right_id:
                raise SequenceError(f"Read-pair identifiers differ at record {count + 1}.")
            if left_mate not in {None, 1} or right_mate not in {None, 2}:
                raise SequenceError(f"Read-pair mate indicators are reversed or invalid at record {count + 1}.")
            if (left_mate is None) != (right_mate is None):
                raise SequenceError(f"Only one read has a mate indicator at record {count + 1}.")
            explicit = explicit and left_mate == 1 and right_mate == 2
            count += 1
            if left.complete != right.complete:
                raise SequenceError("Read-pair files contain different numbers of records.")
        complete = left.complete and right.complete
    return {
        "records_checked": count,
        "sampled": not complete,
        "complete_file": complete,
        "explicit_mates": explicit,
        "status": "verified" if complete and explicit else "prefix_verified" if explicit else "unmarked",
    }
