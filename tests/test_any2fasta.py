# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Acceptance tests for wmlst/any2fasta.py (docs/ARCHITECTURE.md 4.3, 8.4, 13.4).

The contract is byte-identity with the real Perl ``any2fasta 0.9.0`` run as
``any2fasta -q FILE``.  The expected bytes are pinned as SHA-256 digests taken
from the real converter on all 13 fixtures; when a real ``any2fasta`` is
reachable (``$ANY2FASTA``, or on ``$PATH``, with ``perl`` available) the live
byte comparison is run as well.

Runnable two ways:  ``python3 -m pytest tests/test_any2fasta.py``
                    ``python3 tests/test_any2fasta.py``
"""

from __future__ import annotations

import hashlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from wmlst import any2fasta
from wmlst.engine import EmptyInputError, UnsupportedFormatError, WmlstError

DATA = os.path.join(REPO, "tests", "data")

# name -> (sha256 of the expected FASTA, n_sequences, total_bp, format)
# Digests were taken from the real any2fasta 0.9.0; issue146.fa is CRLF input and
# its digest is therefore of the LF-normalised bytes (section 8.4 mandates
# universal-newline decoding, verified not to change any BLAST hit).
EXPECTED = {
    "empty.fa": ("c1c7fae320b011b56d205607173ea36a9642f5c7577eb8d4177015fd380df7ec", 1, 0, "FASTA"),
    "equality.fa.gz": ("335cd52c20fe7f6e30fcfcefd1c1d81da7bd0efec92aa0e7eef2bac5a4181eaa", 158, 5052677, "FASTA"),
    "example.fna": ("a18f883d81d4a09350526706698a9e2a0637e3240f2acd638a60fa854c8ded39", 42, 2476164, "FASTA"),
    "example.fna.gz": ("a18f883d81d4a09350526706698a9e2a0637e3240f2acd638a60fa854c8ded39", 42, 2476164, "FASTA"),
    "example.gbk.gz": ("0e8584e6fe3cf62199be32fddeab6793bd065c2e04c1bf8ca5f522ccbbb5eef8", 42, 2476164, "GENBANK"),
    "issue146.fa": ("8a046c79315712fef10ff2051c86e1b34b24e0f87365f61bc72d4131cabe45c7", 1, 4954218, "FASTA"),
    "messy.fa": ("2c95b88e35e65b093a129aa48a2524d3e0e74bbc69555cb15537f703a3371a85", 71, 2304338, "FASTA"),
    "mixed.fa.zip": ("6d4da1b76913aa43ef8961cb2770716e786d0f039aac46dc11e9c6746c42b06a", 2, 832003, "FASTA"),
    "none.fa": ("2bc224029a4fcd58108434641fed15beb00438f379d9fb580cd1a29e6251035e", 1, 59993, "FASTA"),
    "novel.fa": ("fa9a890d79267aa0f3dbe55cc86587cee5fc676c1c389cf67bd9838f1b25f9a2", 1, 582003, "FASTA"),
    "novel.fasta.bz2": ("98ed325751307695e362731efdc804ae563f4723adc1e48a147b129f04c39726", 4, 5051077, "FASTA"),
}

#: The two fixtures that are errors, not conversions.
EXPECTED_ERRORS = {
    "null.fa": (EmptyInputError, "The input appears to be empty"),
    "fofn.txt": (UnsupportedFormatError, "Unfamilar format with first line: example.fna"),
}

#: CRLF input is normalised to LF, so it cannot be compared byte-for-byte.
CRLF_FIXTURES = ("issue146.fa",)


def _convert_bytes(src):
    with tempfile.TemporaryDirectory() as tmp:
        dest = os.path.join(tmp, "mlst.fna")
        result = any2fasta.convert_to_file(src, dest)
        with open(dest, "rb") as fh:
            return fh.read(), result


def _write(tmp, name, data, mode="w"):
    path = os.path.join(tmp, name)
    with open(path, mode, newline="" if mode == "w" else None) as fh:
        fh.write(data)
    return path


# ---------------------------------------------------------------------------
# The acceptance test: all 13 fixtures
# ---------------------------------------------------------------------------

def test_all_fixtures_match_the_real_any2fasta():
    for name, (digest, nseq, bp, fmt) in sorted(EXPECTED.items()):
        raw, (got_n, got_bp, got_fmt) = _convert_bytes(os.path.join(DATA, name))
        assert hashlib.sha256(raw).hexdigest() == digest, "%s: bytes differ" % name
        assert (got_n, got_bp, got_fmt) == (nseq, bp, fmt), "%s: %r" % (
            name, (got_n, got_bp, got_fmt))


def test_the_two_error_fixtures():
    for name, (exc_type, message) in sorted(EXPECTED_ERRORS.items()):
        try:
            _convert_bytes(os.path.join(DATA, name))
        except exc_type as exc:
            assert str(exc) == message, "%s: %r" % (name, str(exc))
        else:
            raise AssertionError("%s did not raise %s" % (name, exc_type.__name__))


def test_live_byte_identity_against_real_any2fasta():
    """Byte-for-byte against the real Perl converter, when one is available."""
    exe = os.environ.get("ANY2FASTA") or shutil.which("any2fasta")
    if not exe or not shutil.which("perl"):
        raise unittest.SkipTest("no real any2fasta/perl available")
    for name in sorted(EXPECTED):
        src = os.path.join(DATA, name)
        proc = subprocess.run(
            ["perl", exe, "-q", src] if exe.endswith((".pl", "any2fasta")) else [exe, "-q", src],
            capture_output=True, check=True)
        expected = proc.stdout
        if name in CRLF_FIXTURES:
            expected = expected.replace(b"\r\n", b"\n")
        got, _ = _convert_bytes(src)
        assert got == expected, "%s differs from the real any2fasta" % name


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

def test_detect_format_table():
    cases = [
        ("LOCUS       NZ_AHMY02000075  683 bp\n", "GENBANK"),
        ("LOCUS\tX\n", "GENBANK"),
        ("ID   Reviewed; 259 AA.\n", "UNIPROT"),
        ("ID   K02675; SV 1; linear\n", "EMBL"),
        (">contig1 some description\n", "FASTA"),
        ("@read1\n", "FASTQ"),
        ("##gff-version 3\n", "GFF"),
    ]
    for line, expected in cases:
        assert any2fasta.detect_format(line) == expected, line


def test_defline_with_a_space_after_the_angle_bracket_is_not_fasta():
    # ^>\S -- '> contig1' must NOT match (section 4.3).
    for bad in ("> contig1\n", "IDx foo\n", "LOCUSX y\n", "example.fna\n", "\n"):
        try:
            any2fasta.detect_format(bad)
        except UnsupportedFormatError as exc:
            assert str(exc).startswith("Unfamilar format with first line: ")
        else:
            raise AssertionError("accepted %r" % bad)


def test_formats_constant():
    assert any2fasta.FORMATS == ("GENBANK", "UNIPROT", "EMBL", "FASTA", "FASTQ", "GFF")


# ---------------------------------------------------------------------------
# Compression sniffing
# ---------------------------------------------------------------------------

def test_sniff_compression_is_magic_based_not_extension_based():
    expect = {
        "example.fna": "plain", "example.fna.gz": "gzip", "example.gbk.gz": "gzip",
        "equality.fa.gz": "gzip", "novel.fasta.bz2": "bzip2", "mixed.fa.zip": "zip",
        "null.fa": "plain", "empty.fa": "plain", "fofn.txt": "plain",
    }
    for name, comp in sorted(expect.items()):
        assert any2fasta.sniff_compression(os.path.join(DATA, name)) == comp, name


def test_xz_roundtrip():
    import lzma
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.fa.xz")
        with lzma.open(path, "wb") as fh:
            fh.write(b">a\nACGT\n")
        assert any2fasta.sniff_compression(path) == "xz"
        out = io.StringIO()
        assert any2fasta.convert(path, out) == (1, 4, "FASTA")
        assert out.getvalue() == ">a\nACGT\n"


def test_zip_reads_only_the_first_member_and_warns():
    import zipfile
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "two.zip")
        with zipfile.ZipFile(path, "w") as zf:
            zf.writestr("a.fa", ">a\nACGT\n")
            zf.writestr("b.fa", ">b\nTTTT\n")
        out = io.StringIO()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            any2fasta.convert(path, out)
        assert out.getvalue() == ">a\nACGT\n"
        assert any(issubclass(w.category, RuntimeWarning) for w in caught)


# ---------------------------------------------------------------------------
# Truncated / corrupt compressed input (finding 05, finding 18)
# ---------------------------------------------------------------------------

def _corrupt_inputs(tmp):
    """name -> path, one per way a compressed stream can fail to decompress."""
    import bz2 as _bz2
    import gzip as _gzip
    import lzma as _lzma
    import zipfile as _zipfile
    paths = {}

    body = b">a\n" + b"ACGT" * 4096 + b"\n"

    def _trunc(name, data, keep):
        path = os.path.join(tmp, name)
        with open(path, "wb") as fh:
            fh.write(data[:keep])
        paths[name] = path

    _trunc("trunc.fa.gz", _gzip.compress(body), 40)
    _trunc("trunc.fa.bz2", _bz2.compress(body), 60)
    _trunc("trunc.fa.xz", _lzma.compress(body), 80)

    zpath = os.path.join(tmp, "whole.zip")
    with _zipfile.ZipFile(zpath, "w", _zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("a.fa", body)
    with open(zpath, "rb") as fh:
        _trunc("trunc.fa.zip", fh.read(), 60)

    # A damaged deflate body: gzip's own header and CRC are intact enough that
    # the failure comes out of zlib, not as gzip.BadGzipFile.
    blob = bytearray(_gzip.compress(body))
    blob[30] ^= 0xFF
    path = os.path.join(tmp, "bad-deflate.fa.gz")
    with open(path, "wb") as fh:
        fh.write(bytes(blob))
    paths["bad-deflate.fa.gz"] = path
    return paths


def test_truncated_or_corrupt_compressed_input_raises_wmlst_error():
    # Finding 05: EOFError / zlib.error / lzma.LZMAError / zipfile.BadZipFile
    # are neither WmlstError nor OSError, so if they escape convert() they kill
    # the whole batch with a traceback instead of failing one file
    # (docs/ARCHITECTURE.md D10, and the GUI shows the stack trace).
    with tempfile.TemporaryDirectory() as tmp:
        for name, path in sorted(_corrupt_inputs(tmp).items()):
            try:
                any2fasta.convert(path, io.StringIO())
            except WmlstError as exc:
                assert str(exc).startswith("Could not read '"), "%s: %r" % (name, str(exc))
                assert exc.user_message == "That file is truncated or corrupt."
            except Exception as exc:  # pragma: no cover - the bug being fixed
                raise AssertionError(
                    "%s escaped as %s: %s"
                    % (name, type(exc).__name__, exc)) from exc
            else:
                raise AssertionError("%s did not raise" % name)


def test_corrupt_input_still_fails_per_file_through_convert_to_file():
    with tempfile.TemporaryDirectory() as tmp:
        path = _corrupt_inputs(tmp)["trunc.fa.gz"]
        dest = os.path.join(tmp, "mlst.fna")
        try:
            any2fasta.convert_to_file(path, dest)
        except WmlstError:
            pass
        else:
            raise AssertionError("no WmlstError")
        # engine.py only catches WmlstError/OSError; the output handle must
        # still have been closed, hence removable (section 5.4).
        os.remove(dest)


def _is_closed(obj):
    closed = getattr(obj, "closed", None)
    if closed is not None:
        return bool(closed)
    return getattr(obj, "fp", None) is None  # zipfile.ZipFile


def _open_text_closers(exc):
    """The `closers` list of the open_text frame in `exc`'s traceback."""
    tb = exc.__traceback__
    while tb is not None:
        if tb.tb_frame.f_code.co_name == "open_text":
            return tb.tb_frame.f_locals.get("closers")
        tb = tb.tb_next
    return None


def test_open_text_closes_its_handles_when_it_raises():
    # Finding 18: the raised exception is stored on the failed SampleResult
    # (engine.py), and its __traceback__ pins open_text's frame -- so an
    # unclosed handle in `closers` stays open for the whole run.
    with tempfile.TemporaryDirectory() as tmp:
        path = _corrupt_inputs(tmp)["trunc.fa.zip"]  # ZipFile() raises in open_text
        try:
            any2fasta.open_text(path)
        except Exception as exc:
            closers = _open_text_closers(exc)
            assert closers, "open_text frame not found in the traceback"
            assert all(_is_closed(o) for o in closers), (
                "leaked: %r" % [o for o in closers if not _is_closed(o)])
        else:
            raise AssertionError("open_text did not raise")


def test_open_text_closes_its_handles_when_zstd_is_unavailable():
    # The branch that actually bites: on Python < 3.14 every .zst input raises
    # UnsupportedFormatError, which the engine keeps -- 200 of them used to mean
    # 200 held file descriptors.
    def _boom(binary):
        raise UnsupportedFormatError("zstd-compressed input needs Python 3.14 or newer")

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa.zst")
        with open(path, "wb") as fh:
            fh.write(b"\x28\xb5\x2f\xfd" + b"\x00" * 32)
        assert any2fasta.sniff_compression(path) == "zstd"
        saved = any2fasta._zstd_stream
        any2fasta._zstd_stream = _boom
        try:
            kept = []
            for _ in range(20):
                try:
                    any2fasta.convert(path, io.StringIO())
                except UnsupportedFormatError as exc:
                    kept.append(exc)  # engine.py stores the exception; so do we
            assert len(kept) == 20
            for exc in kept:
                closers = _open_text_closers(exc)
                assert closers, "open_text frame not found in the traceback"
                assert all(_is_closed(o) for o in closers)
            if os.path.isdir("/proc/self/fd"):  # Linux: prove no accumulation
                assert len(os.listdir("/proc/self/fd")) < 40
        finally:
            any2fasta._zstd_stream = saved


# ---------------------------------------------------------------------------
# Parser behaviour
# ---------------------------------------------------------------------------

def test_pass_through_does_not_uppercase_or_filter():
    # C9: strict pass-through -- no case change, no character substitution,
    # no re-wrapping.
    payload = ">m ixed desc\nacgtRYKMxxx-*\nACGTN\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "p.fa", payload)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
    assert out.getvalue() == payload
    assert (n, bp, fmt) == (1, 18, "FASTA")


def test_blank_lines_are_dropped_but_nothing_else_is():
    payload = ">a\n\nACGT\n   \n\t\nTTTT\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "b.fa", payload)
        out = io.StringIO()
        any2fasta.convert(path, out)
    assert out.getvalue() == ">a\nACGT\nTTTT\n"


def test_missing_final_newline_is_preserved():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "n.fa", ">a\nACGT")
        out = io.StringIO()
        any2fasta.convert(path, out)
    assert out.getvalue() == ">a\nACGT"


def test_crlf_input_is_normalised_to_lf():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "c.fa", ">a\r\nACGT\r\nTTTT\r\n")
        out = io.StringIO()
        n, bp, _ = any2fasta.convert(path, out)
    assert out.getvalue() == ">a\nACGT\nTTTT\n"
    assert (n, bp) == (1, 8)


def test_utf8_bom_is_stripped():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "bom.fa")
        with open(path, "wb") as fh:
            fh.write(b"\xef\xbb\xbf>a\nACGT\n")
        out = io.StringIO()
        assert any2fasta.convert(path, out) == (1, 4, "FASTA")
        assert out.getvalue() == ">a\nACGT\n"


def test_zero_byte_file_is_empty_input():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "z.fa", "")
        try:
            any2fasta.convert(path, io.StringIO())
        except EmptyInputError as exc:
            assert str(exc) == "The input appears to be empty"
        else:
            raise AssertionError("no EmptyInputError")


def test_bare_zero_is_empty_because_perl_says_so():
    # `if (not $header)` — the string "0" is false in Perl (any2fasta:135).
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "zero.fa", "0")
        try:
            any2fasta.convert(path, io.StringIO())
        except EmptyInputError as exc:
            assert str(exc) == "The input appears to be empty"
        else:
            raise AssertionError("no EmptyInputError")


def test_genbank_parser():
    gbk = (
        "LOCUS       ACC1  12 bp    DNA     linear   CON 23-NOV-2017\n"
        "VERSION     ACC1.3\n"
        "ORIGIN\n"
        "        1 acgtac gtacgt\n"
        "//\n"
        "LOCUS       ACC2  4 bp    DNA\n"
        "ORIGIN\n"
        "        1 tttt\n"
        "//\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "g.gbk", gbk)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
    # -g is OFF, so the accession is LOCUS's, never VERSION's.
    assert out.getvalue() == ">ACC1\nacgtacgtacgt\n>ACC2\ntttt\n"
    assert (n, bp, fmt) == (2, 16, "GENBANK")


def test_genbank_record_without_a_terminator_is_dropped():
    gbk = "LOCUS       ACC1  4 bp\nORIGIN\n        1 acgt\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "g2.gbk", gbk)
        try:
            any2fasta.convert(path, io.StringIO())
        except EmptyInputError as exc:
            assert str(exc).startswith("No sequences found in ")
        else:
            raise AssertionError("no EmptyInputError")


def test_embl_parser_strips_spaces_then_a_trailing_coordinate():
    embl = (
        "ID   K02675; SV 1; linear; genomic DNA; STD; UNC; 12 BP.\n"
        "SQ   Sequence 12 BP; 3 A; 3 C; 3 G; 3 T; 0 other;\n"
        "     agtcgcttt taa        12\n"
        "//\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "e.embl", embl)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
    assert out.getvalue() == ">K02675\nagtcgcttttaa\n"
    assert (n, bp, fmt) == (1, 12, "EMBL")


def test_uniprot_is_routed_to_the_embl_parser():
    rec = "ID   Reviewed; 5 AA.\nSQ   SEQUENCE   5 AA;\n     MGIFD\n//\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "u.dat", rec)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
    assert fmt == "UNIPROT" and out.getvalue() == ">Reviewed\nMGIFD\n"
    assert (n, bp) == (1, 5)


def test_gff_uses_the_embedded_fasta_only():
    gff = (
        "##gff-version 3\n"
        "ctg1\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=x\n"
        "##FASTA\n"
        ">ctg1\nACGTACGTA\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "a.gff", gff)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
    assert out.getvalue() == ">ctg1\nACGTACGTA\n"
    assert (n, bp, fmt) == (1, 9, "GFF")


def test_gff_without_an_embedded_fasta_is_empty():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "b.gff", "##gff-version 3\nctg1\tx\tCDS\t1\t9\t.\t+\t0\tID=x\n")
        try:
            any2fasta.convert(path, io.StringIO())
        except EmptyInputError as exc:
            assert str(exc).startswith("No sequences found in ")
        else:
            raise AssertionError("no EmptyInputError")


def test_fastq_stride_and_incomplete_trailing_record():
    fq = "@r1 d\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "r.fq", fq)
        out = io.StringIO()
        n, bp, fmt = any2fasta.convert(path, out)
        assert out.getvalue() == ">r1 d\nACGT\n>r2\nTTTT\n"
        assert (n, bp, fmt) == (2, 8, "FASTQ")

        path = _write(tmp, "r2.fq", "@r1\nACGT\n+\n!!!!\n@r2\n")
        out = io.StringIO()
        assert any2fasta.convert(path, out)[0] == 1

        path = _write(tmp, "r3.fq", "@r1\n")
        try:
            any2fasta.convert(path, io.StringIO())
        except EmptyInputError:
            pass
        else:
            raise AssertionError("a lone id line must produce no records")


def test_messy_single_long_line_survives_unchanged():
    raw, (n, _bp, _) = _convert_bytes(os.path.join(DATA, "messy.fa"))
    longest = max(len(line) for line in raw.split(b"\n"))
    assert longest == 282005, longest
    assert n == 71


def test_stdin_is_supported():
    original = sys.stdin
    any2fasta._STDIN_READER = None
    try:
        class _Fake:
            buffer = io.BytesIO(b">a\nACGT\n")
        sys.stdin = _Fake()
        out = io.StringIO()
        assert any2fasta.convert("-", out) == (1, 4, "FASTA")
        assert out.getvalue() == ">a\nACGT\n"
    finally:
        sys.stdin = original
        any2fasta._STDIN_READER = None


# ---------------------------------------------------------------------------
# revcom
# ---------------------------------------------------------------------------

def test_revcom_is_perl_faithful():
    assert any2fasta.revcom("ATGC") == "GCAT"
    assert any2fasta.revcom("atgc") == "gcat"
    assert any2fasta.revcom("AaTtGgCc") == "gGcCaAtT"
    # IUPAC codes are reversed but NOT complemented (bin/mlst:287-292).
    assert any2fasta.revcom("ACGTRYN") == "NYRACGT"
    assert any2fasta.revcom("") == ""


# ---------------------------------------------------------------------------
# convert_to_file
# ---------------------------------------------------------------------------

def test_convert_to_file_closes_the_handle_and_writes_lf():
    with tempfile.TemporaryDirectory() as tmp:
        src = _write(tmp, "s.fa", ">a\r\nACGT\r\n")
        dest = os.path.join(tmp, "mlst.fna")
        any2fasta.convert_to_file(src, dest)
        with open(dest, "rb") as fh:
            assert fh.read() == b">a\nACGT\n"
        # The file must be closed, hence removable, before blastn is launched.
        os.remove(dest)


def test_no_named_temporary_file_anywhere():
    # Section 5.4: NamedTemporaryFile's O_TEMPORARY share mode gives blastn
    # ERROR_SHARING_VIOLATION (32) on Windows.
    with open(os.path.join(REPO, "wmlst", "any2fasta.py"), encoding="utf-8") as fh:
        assert "NamedTemporaryFile(" not in fh.read()


# ---------------------------------------------------------------------------
# plain-python runner
# ---------------------------------------------------------------------------

def _main():
    failures = 0
    skipped = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
        except unittest.SkipTest as exc:
            skipped += 1
            print("SKIP %s (%s)" % (name, exc))
        except Exception as exc:
            failures += 1
            import traceback
            print("FAIL %s: %s" % (name, exc))
            traceback.print_exc()
        else:
            print("ok   %s" % name)
    print("\n%d failed, %d skipped" % (failures, skipped))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
