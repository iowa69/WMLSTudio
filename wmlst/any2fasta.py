# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Pure-Python replacement for ``any2fasta`` 0.9.0 as invoked by ``mlst``.

Upstream ``bin/mlst:307`` runs ``any2fasta -q FILE > mlst.fna``.  ``-q`` only, so
``-n`` / ``-l`` / ``-u`` / ``-s`` / ``-g`` / ``-p`` / ``-k`` are all OFF and the
converter is a strict **pass-through**: no case folding, no character filtering,
no re-wrapping (docs/ARCHITECTURE.md C9, D6, section 4.3).

Implements docs/ARCHITECTURE.md section 4.3; text/encoding policy from section 8.4.

Divergences from the Perl, both mandated by section 8.4:

* input is decoded with universal newlines, so a CRLF source file is converted to
  an LF FASTA (the Perl copies CRLF through verbatim).  BLAST is insensitive to
  the difference; verified on ``tests/data/issue146.fa``.
* a UTF-8 BOM is stripped rather than being fed to the format sniffer.
"""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import os
import re
import sys
import warnings
import zipfile
from collections.abc import Iterable, Iterator
from typing import Callable, Dict, List, Optional, Tuple

__all__ = [
    "FORMATS",
    "VERSION",
    "convert",
    "convert_to_file",
    "detect_format",
    "open_text",
    "revcom",
    "sniff_compression",
]

#: The any2fasta release this module reproduces.
VERSION = "0.9.0"

#: Formats WMLST accepts, in detection order (any2fasta:88-99, section 4.3).
FORMATS: Tuple[str, ...] = ("GENBANK", "UNIPROT", "EMBL", "FASTA", "FASTQ", "GFF")

_EXC_CACHE: Dict[str, type] = {}


def _exc(name: str) -> type:
    """Resolve an exception class from ``wmlst.engine`` lazily.

    ``engine`` imports this module (section 2.1), so a module-level
    ``from wmlst.engine import ...`` would be a circular import.
    """
    cls = _EXC_CACHE.get(name)
    if cls is None:
        from wmlst import engine  # local import: see docstring

        cls = getattr(engine, name)
        _EXC_CACHE[name] = cls
    return cls


def __getattr__(name: str):  # PEP 562 — re-export the engine exceptions lazily
    if name in ("WmlstError", "EmptyInputError", "UnsupportedFormatError"):
        return _exc(name)
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


# ---------------------------------------------------------------------------
# Compression
# ---------------------------------------------------------------------------

_MAGIC: Tuple[Tuple[bytes, str], ...] = (
    (b"\x1f\x8b\x08", "gzip"),
    (b"BZh", "bzip2"),
    (b"\xfd7zXZ\x00", "xz"),
    (b"PK\x03\x04", "zip"),
    (b"\x28\xb5\x2f\xfd", "zstd"),
)

_STDIN_READER: Optional[io.BufferedReader] = None


def _stdin_reader() -> io.BufferedReader:
    """A peekable buffered view of ``sys.stdin``, created at most once."""
    global _STDIN_READER
    buf = getattr(sys.stdin, "buffer", sys.stdin)
    if _STDIN_READER is None or _STDIN_READER.raw is not buf:  # type: ignore[attr-defined]
        _STDIN_READER = io.BufferedReader(buf)  # type: ignore[arg-type]
    return _STDIN_READER


def sniff_compression(path: str) -> str:
    """Magic-byte sniff of `path`; NEVER extension-based (section 4.3).

    -> one of 'plain' | 'gzip' | 'bzip2' | 'xz' | 'zip' | 'zstd'.
    ``path == '-'`` peeks at stdin without consuming it.
    """
    if path == "-":
        head = bytes(_stdin_reader().peek(8)[:8])
    else:
        # Opened through a BufferedReader and *peeked*, never read: `path` may
        # name a pipe or FIFO (`/dev/stdin`, a process substitution), where a
        # consuming read would swallow bytes that open_text() must still see.
        fh = open(path, "rb", buffering=0)
        try:
            head = bytes(io.BufferedReader(fh).peek(8)[:8])
        finally:
            fh.close()
    for magic, name in _MAGIC:
        if head.startswith(magic):
            return name
    return "plain"


def _sniff_stream(reader: io.BufferedReader) -> str:
    """Magic-byte sniff of an already-open buffered stream, without consuming it."""
    head = bytes(reader.peek(8)[:8])
    for magic, name in _MAGIC:
        if head.startswith(magic):
            return name
    return "plain"


class _TextSource:
    """An iterator of decoded lines that owns and closes its whole stream stack."""

    def __init__(self, text: io.TextIOBase, closers: Iterable[object]) -> None:
        self._text = text
        self._closers = list(closers)

    def __iter__(self) -> Iterator[str]:
        return iter(self._text)

    def __next__(self) -> str:
        return next(self._text)  # type: ignore[arg-type]

    def readline(self) -> str:
        return self._text.readline()

    def read(self, size: int = -1) -> str:
        return self._text.read(size)

    def close(self) -> None:
        try:
            self._text.detach()
        except Exception:
            pass
        for obj in self._closers:
            close = getattr(obj, "close", None)
            if close is not None:
                try:
                    close()
                except Exception:
                    pass
        self._closers = []

    def __enter__(self) -> _TextSource:
        return self

    def __exit__(self, *exc) -> bool:
        self.close()
        return False


def _zstd_stream(binary):
    """Open a zstd stream, or explain why we cannot (section 4.3)."""
    try:
        from compression.zstd import ZstdFile  # Python 3.14+
    except ImportError as exc:
        raise _exc("UnsupportedFormatError")(
            "zstd-compressed input needs Python 3.14 or newer; "
            "decompress the file first (zstd -d) and try again"
        ) from exc
    return ZstdFile(binary, "rb")


def open_text(path: str):
    """Open `path` as a line iterator of decoded text (sections 4.3, 8.4).

    Universal newlines, ``utf-8`` with ``errors='replace'`` (``utf-8-sig`` when a
    BOM is present — the ``utf-8-sig`` codec is a no-op on BOM-less input, so it
    is used unconditionally).  Decompresses gzip/bzip2/xz/zip/zstd transparently.
    ``zip`` reads ``namelist()[0]`` ONLY, matching ``IO::Uncompress::Unzip``, and
    warns when the archive holds more than one member.
    ``path == '-'`` reads ``sys.stdin.buffer``.
    """
    closers: List[object] = []

    if path == "-":
        raw = _stdin_reader()
    else:
        # One open, then peek: re-opening to sniff would lose the first bytes of
        # a non-seekable path such as /dev/stdin.
        handle = open(path, "rb", buffering=0)
        closers.append(handle)
        raw = io.BufferedReader(handle)
        closers.insert(0, raw)
    comp = _sniff_stream(raw)

    if comp == "plain":
        stream = raw
    elif comp == "gzip":
        stream = gzip.GzipFile(fileobj=raw, mode="rb")
        closers.insert(0, stream)
    elif comp == "bzip2":
        stream = bz2.BZ2File(raw, "rb")
        closers.insert(0, stream)
    elif comp == "xz":
        stream = lzma.LZMAFile(raw, "rb")
        closers.insert(0, stream)
    elif comp == "zstd":
        stream = _zstd_stream(raw)
        closers.insert(0, stream)
    elif comp == "zip":
        if not raw.seekable():  # zipfile needs random access; a pipe has none
            raw = io.BytesIO(raw.read())
        zf = zipfile.ZipFile(raw)
        names = zf.namelist()
        if not names:
            zf.close()
            for obj in closers:
                getattr(obj, "close", lambda: None)()
            raise _exc("EmptyInputError")("The input appears to be empty")
        if len(names) > 1:
            warnings.warn(
                "zip archive %r holds %d members; only %r is read"
                % (path, len(names), names[0]),
                RuntimeWarning,
                stacklevel=2,
            )
        stream = zf.open(names[0], "r")
        closers.insert(0, zf)
        closers.insert(0, stream)
    else:  # pragma: no cover - _MAGIC and this branch list are the same set
        stream = raw

    text = io.TextIOWrapper(
        stream, encoding="utf-8-sig", errors="replace", newline=None
    )
    return _TextSource(text, closers)


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

_DETECTORS: Tuple[Tuple[str, re.Pattern[str]], ...] = (
    ("GENBANK", re.compile(r"^LOCUS[ \t]", re.ASCII)),
    ("UNIPROT", re.compile(r"^ID[ \t]+Reviewed", re.ASCII)),
    ("EMBL", re.compile(r"^ID[ \t]", re.ASCII)),
    ("FASTA", re.compile(r"^>\S", re.ASCII)),
    ("FASTQ", re.compile(r"^@\S", re.ASCII)),
    ("GFF", re.compile(r"^##gff", re.ASCII)),
)


def detect_format(first_line: str) -> str:
    """Classify a file from its FIRST LINE ONLY, first match wins (section 4.3).

    Raises ``UnsupportedFormatError`` otherwise.  The upstream misspelling
    'Unfamilar' is preserved (any2fasta:159); the trailing newline is trimmed
    because the caller renders the message as ``ERROR: <msg>``.
    """
    for name, pattern in _DETECTORS:
        if pattern.match(first_line):
            return name
    raise _exc("UnsupportedFormatError")(
        "Unfamilar format with first line: " + first_line.rstrip("\r\n")
    )


# ---------------------------------------------------------------------------
# Parsers.  Each takes an iterator of lines (terminators included) and a text
# stream, and returns (n_sequences, total_bp).
# ---------------------------------------------------------------------------

_BLANK_RE = re.compile(r"^\s*$", re.ASCII)
_WS_RE = re.compile(r"\s", re.ASCII)
_HSPACE_RE = re.compile(r"[ \t]", re.ASCII)
_TRAIL_NUM_RE = re.compile(r"\d+$", re.ASCII)
_LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)", re.ASCII)
_ID_RE = re.compile(r"^ID[ \t]+([^;\t ]+)", re.ASCII)
_SQ_RE = re.compile(r"^SQ\s", re.ASCII)


def _parse_fasta(lines: Iterator[str], out) -> Tuple[int, int]:
    """any2fasta:202-225 — drop blank lines, emit every other line verbatim."""
    count = 0
    total = 0
    for line in lines:
        if _BLANK_RE.match(line):
            continue
        if line.startswith(">"):
            count += 1
        else:
            total += len(line.rstrip("\n"))
        out.write(line)
    return count, total


def _parse_fastq(lines: Iterator[str], out) -> Tuple[int, int]:
    """any2fasta:229-240 — strict 4-line stride; incomplete trailing record dropped."""
    count = 0
    total = 0
    group: List[str] = []
    for line in lines:
        group.append(line)
        if len(group) == 4:
            count, total = _emit_fastq(group, out, count, total)
            group = []
    # Perl's loop guard is `$i < $#lines`, i.e. a record needs the id line plus
    # at least one more line to exist; quality lines are never read.
    if len(group) >= 2:
        count, total = _emit_fastq(group, out, count, total)
    return count, total


def _emit_fastq(group: List[str], out, count: int, total: int) -> Tuple[int, int]:
    head = group[0]
    seq = group[1]
    out.write(">" + head[1:])
    out.write(seq)
    return count + 1, total + len(seq.rstrip("\n"))


def _parse_gff(lines: Iterator[str], out) -> Tuple[int, int]:
    """any2fasta:244-255 — skip to the ##FASTA payload, then run the FASTA parser."""
    for line in lines:
        if line.startswith(">"):
            return _parse_fasta(_prepend(line, lines), out)
    return 0, 0


def _parse_genbank(lines: Iterator[str], out) -> Tuple[int, int]:
    """any2fasta:280-323 — ORIGIN opens the block, ``//`` flushes the record."""
    acc = ""
    dna: List[str] = []
    in_seq = False
    count = 0
    total = 0
    for raw in lines:
        line = raw[:-1] if raw.endswith("\n") else raw  # chomp
        if line.startswith("//"):
            out.write(">" + acc + "\n")
            out.write("".join(dna))
            count += 1
            in_seq = False
            dna = []
            acc = ""
            continue
        if line.startswith("ORIGIN"):
            in_seq = True
            continue
        if in_seq:
            # A coordinate of 10 or more digits eats a base here. That is an
            # upstream bug (any2fasta:314) and is reproduced deliberately.
            chunk = _WS_RE.sub("", line[10:])
            dna.append(chunk + "\n")
            total += len(chunk)
        else:
            m = _LOCUS_RE.match(line)
            if m:
                acc = m.group(1)
    return count, total


def _parse_embl(lines: Iterator[str], out) -> Tuple[int, int]:
    """any2fasta:327-371 — SQ opens the block, ``//`` flushes the record."""
    acc = ""
    dna: List[str] = []
    in_seq = False
    count = 0
    total = 0
    for raw in lines:
        line = raw[:-1] if raw.endswith("\n") else raw  # chomp
        if line.startswith("//"):
            out.write(">" + acc + "\n")
            out.write("".join(dna))
            count += 1
            in_seq = False
            dna = []
            acc = ""
            continue
        if _SQ_RE.match(line):
            in_seq = True
            continue
        if in_seq:
            chunk = _TRAIL_NUM_RE.sub("", _HSPACE_RE.sub("", line))
            dna.append(chunk + "\n")
            total += len(chunk)
        else:
            m = _ID_RE.match(line)
            if m:
                acc = m.group(1)
    return count, total


_PARSERS: Dict[str, Callable[[Iterator[str], object], Tuple[int, int]]] = {
    "GENBANK": _parse_genbank,
    "UNIPROT": _parse_embl,
    "EMBL": _parse_embl,
    "FASTA": _parse_fasta,
    "FASTQ": _parse_fastq,
    "GFF": _parse_gff,
}


def _prepend(first: str, rest: Iterator[str]) -> Iterator[str]:
    yield first
    yield from rest


# ---------------------------------------------------------------------------
# Public conversion API
# ---------------------------------------------------------------------------


def convert(src: str, out) -> Tuple[int, int, str]:
    """Convert `src` to FASTA on the text stream `out` (section 4.3).

    `out` MUST have been opened ``newline='\\n'``, ``encoding='utf-8'``.
    -> ``(n_sequences, total_bp, format_name)``.
    Raises ``EmptyInputError('The input appears to be empty')`` when there is no
    first line, and ``EmptyInputError("No sequences found in '<src>'")`` when the
    parser produced no records.  ``UnsupportedFormatError`` on an unknown first line.
    """
    with open_text(src) as fh:
        lines = iter(fh)
        try:
            first = next(lines)
        except StopIteration:
            first = ""
        # Perl tests `if (not $header)`, and the string "0" is false in Perl, so a
        # file whose entire content is a bare `0` is reported as empty (any2fasta:135).
        if first == "" or first == "0":
            raise _exc("EmptyInputError")("The input appears to be empty")
        fmt = detect_format(first)
        count, total = _PARSERS[fmt](_prepend(first, lines), out)
    if count == 0:
        raise _exc("EmptyInputError")("No sequences found in '%s'" % (src,))
    return count, total, fmt


def convert_to_file(src: str, dest: str) -> Tuple[int, int, str]:
    """Convert `src` into the new file `dest`, then flush, fsync and CLOSE it.

    Implements section 5.4 step 25.  The handle must be closed before ``blastn``
    is launched or Windows answers with ``ERROR_SHARING_VIOLATION (32)``;
    ``tempfile.NamedTemporaryFile`` must not be used for the same reason.
    """
    fh = open(dest, "w", encoding="utf-8", errors="replace", newline="\n")
    try:
        result = convert(src, fh)
        fh.flush()
        try:
            os.fsync(fh.fileno())
        except (OSError, ValueError):  # pragma: no cover - e.g. odd filesystems
            pass
    finally:
        fh.close()
    return result


_COMPLEMENT = str.maketrans("ATGCatgc", "TACGtacg")


def revcom(dna: str) -> str:
    """Perl-faithful reverse complement (``bin/mlst:287-292``, section 4.3).

    Reverse FIRST, then ``tr/ATGCatgc/TACGtacg/``.  IUPAC ambiguity codes are
    reversed but NOT complemented — upstream behaviour, and it MUST NOT be fixed.
    """
    return dna[::-1].translate(_COMPLEMENT)
