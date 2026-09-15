"""Auditable exact nucleotide allele calls from local MLST/cgMLST schemas.

This engine does not infer novel alleles, align approximate matches, assemble reads,
or infer a species. A missing exact match is consequently labelled ``missing``.
Coordinates are one-based, inclusive, on the input contig; strands are explicit.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import ahocorasick

from . import __version__
from .sequence import (
    CancelCallback,
    QCAccumulator,
    SequenceError,
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
    sample_name,
)

ProgressCallback = Callable[[int, int, str], None] | None
ALLELE_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".tfa"}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class SchemeError(ValueError):
    """A local scheme cannot be interpreted without making an unsafe assumption."""


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def _natural_key(value: str) -> tuple:
    return tuple((0, int(part)) if part.isdigit() else (1, part.lower())
                 for part in re.split(r"(\d+)", value))


@dataclass
class Scheme:
    name: str
    path: Path
    loci: tuple[str, ...]
    alleles: dict[str, dict[str, str]]
    profiles: dict[tuple[str, ...], tuple[str, ...]]
    metadata: dict
    digest: str
    notes: list[str] = field(default_factory=list)
    #: False when this scheme was loaded to be checked rather than to be typed
    #: against, so its allele sequences were never kept. Nothing may match
    #: against such a scheme: every locus would silently report no hit.
    sequences_loaded: bool = True
    _automaton: object | None = field(default=None, init=False, repr=False, compare=False)

    @property
    def locus_count(self) -> int:
        return len(self.loci)

    @property
    def allele_count(self) -> int:
        return sum(len(values) for values in self.alleles.values())

    @property
    def profile_count(self) -> int:
        return sum(len(values) for values in self.profiles.values())

    @property
    def scheme_digest(self) -> str:
        return self.digest


def _allele_identifier(header: str, locus: str) -> str:
    identifier = header.split()[0]
    for prefix in (f"{locus}_", f"{locus}-"):
        if identifier.startswith(prefix):
            allele = identifier[len(prefix):]
            break
    else:
        allele = identifier if identifier.isdecimal() else ""
    if not allele or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", allele):
        raise SchemeError(
            f"{locus}: allele header {identifier!r} must be a number or '{locus}_<allele>'."
        )
    if allele in {"0", "N", "NA"}:
        raise SchemeError(f"{locus}: {allele!r} is a reserved missing-allele identifier.")
    return allele


def _source_digest(files: list[Path], root: Path, cancelled: CancelCallback,
                   progress: ProgressCallback = None) -> str:
    ordered = sorted(files, key=lambda file: file.relative_to(root).as_posix())
    digest = hashlib.sha256(b"WMLSTudio-schema-v1\0")
    for number, path in enumerate(ordered, 1):
        check_cancelled(cancelled)
        name = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(name).to_bytes(8, "big"))
        digest.update(name)
        digest.update(bytes.fromhex(file_sha256(path, cancelled)))
        if progress:
            progress(number, len(ordered), f"Fingerprinting {path.name}")
    return digest.hexdigest()


def load_scheme(path: str | Path, cancelled: CancelCallback = None, *,
                progress: ProgressCallback = None, sequences: bool = True) -> Scheme:
    """Load a flat local allele FASTA directory and optional tab-separated profiles.

    Locus names are allele-file stems. Profiles use an ST column and one column
    per locus; extra metadata columns are allowed. Allele-only cgMLST directories
    are supported, and never acquire an inferred sequence type.

    ``sequences=False`` checks a scheme without keeping it. Every file is still
    read, parsed and fingerprinted and every error is still raised, but the allele
    sequences are discarded as they are read: only the identifiers are kept, which
    is all the remaining checks and the digest need. A 2,358-target cgMLST scheme
    is several gigabytes of FASTA, and holding it to confirm a freshly downloaded
    copy cost about 1.5 times its own size in memory — on a laptop with 8 GB that
    is the difference between a check that takes a minute and one that swaps for
    an hour and reads as a frozen download. A scheme loaded this way can never be
    typed against: ``_automaton`` refuses it by name rather than quietly matching
    nothing.

    ``progress`` is called as ``progress(done, total, message)`` per file. This
    pass reads the whole scheme twice, so a caller that shows nothing here leaves
    a full progress bar standing still for minutes.
    """
    root = Path(path).resolve()
    if not root.is_dir():
        raise SchemeError(f"Scheme directory does not exist: {root}")
    files = sorted((file for file in root.iterdir() if file.is_file()), key=lambda p: p.name)
    allele_files = []
    for file in files:
        name = file.with_suffix("") if file.suffix.lower() in {".gz", ".bz2"} else file
        if name.suffix.lower() in ALLELE_SUFFIXES:
            allele_files.append((file, name.stem))
    if not allele_files:
        raise SchemeError("No allele FASTA files (.tfa, .fa, .fasta, .fna) in scheme directory.")
    signatures = {file: file_signature(file) for file in files}
    alleles: dict[str, dict[str, str]] = {}
    included: list[Path] = []
    ambiguous_references: list[str] = []
    for number, (file, locus) in enumerate(allele_files, 1):
        check_cancelled(cancelled)
        if progress:
            progress(number, len(allele_files), f"Reading locus {locus}")
        if locus in alleles:
            raise SchemeError(f"More than one allele FASTA file defines locus {locus!r}.")
        values: dict[str, str] = {}
        try:
            with SequenceReader(file, cancelled) as reader:
                if reader.kind != "fasta":
                    raise SchemeError(f"{file.name}: scheme alleles must be FASTA sequences.")
                for record in reader:
                    identifier = _allele_identifier(record.name, locus)
                    if identifier in values:
                        raise SchemeError(f"{file.name}: duplicate allele ID {identifier!r}.")
                    if set(record.sequence) - set("ACGT"):
                        ambiguous_references.append(f"{locus}_{identifier}")
                    # Checking a download needs the identifiers, not the bases; the
                    # sequence is dropped here rather than after gigabytes of it
                    # have been accumulated.
                    values[identifier] = record.sequence if sequences else ""
        except SequenceError as error:
            raise SchemeError(str(error)) from error
        alleles[locus] = values
        included.append(file)

    metadata: dict = {}
    metadata_files = [file for file in files if file.name == "scheme.json"
                      or file.name.endswith("_info.json")]
    if len(metadata_files) > 1:
        preferred = [file for file in metadata_files if file.name == "scheme.json"]
        if not preferred:
            raise SchemeError("Multiple scheme metadata files found; retain a single *_info.json.")
        metadata_files = preferred
    if metadata_files:
        try:
            metadata = json.loads(metadata_files[0].read_text(encoding="utf-8"))
            if not isinstance(metadata, dict):
                raise ValueError("expected a JSON object")
        except (UnicodeDecodeError, ValueError) as error:
            raise SchemeError(f"Cannot read {metadata_files[0].name}: {error}") from error
        included.extend(metadata_files)

    profile_candidates: list[tuple[Path, list[str], int]] = []
    for file in files:
        if file.suffix.lower() not in {".txt", ".tsv"}:
            continue
        try:
            with file.open(encoding="utf-8-sig", newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"), [])
        except (UnicodeDecodeError, csv.Error) as error:
            raise SchemeError(f"Cannot read potential profile file {file.name}: {error}") from error
        header = [value.strip() for value in header]
        profile_fields = {"st", "sequence_type", "sequence type"}
        if isinstance(metadata.get("profile_field"), str):
            profile_fields.add(metadata["profile_field"].casefold())
        st_columns = [index for index, value in enumerate(header)
                      if value.lower() in profile_fields]
        if st_columns:
            if len(st_columns) != 1 or len(set(header)) != len(header):
                raise SchemeError(f"{file.name}: duplicate or ambiguous profile columns.")
            missing = set(alleles) - set(header)
            if missing:
                raise SchemeError(f"{file.name}: profile table lacks loci: {', '.join(sorted(missing))}")
            profile_candidates.append((file, header, st_columns[0]))
    if len(profile_candidates) > 1:
        raise SchemeError("Multiple ST profile tables found; use one profile table per scheme.")

    loci = tuple(sorted(alleles, key=_natural_key))
    profiles: dict[tuple[str, ...], tuple[str, ...]] = {}
    notes: list[str] = []
    if metadata.get("access_notice"):
        notes.append(str(metadata["access_notice"]))
    if ambiguous_references:
        notes.append(
            f"Excluded {len(ambiguous_references)} reference alleles with ambiguous DNA from "
            f"exact matching (including {', '.join(ambiguous_references[:5])}). "
            "Their ambiguous bases are never treated as wildcards."
        )
    if profile_candidates:
        file, header, st_column = profile_candidates[0]
        loci = tuple(column for column in header if column in alleles)
        columns = [header.index(locus) for locus in loci]
        accumulated: dict[tuple[str, ...], set[str]] = {}
        st_profiles: dict[str, tuple[str, ...]] = {}
        incomplete = 0
        unavailable_references: set[str] = set()
        with file.open(encoding="utf-8-sig", newline="") as handle:
            rows = csv.reader(handle, delimiter="\t")
            next(rows)
            for line, row in enumerate(rows, 2):
                check_cancelled(cancelled)
                if not row or not any(value.strip() for value in row):
                    continue
                if len(row) != len(header):
                    raise SchemeError(f"{file.name}, line {line}: profile row width differs from header.")
                st = row[st_column].strip()
                if not st:
                    raise SchemeError(f"{file.name}, line {line}: missing sequence type.")
                profile = tuple(row[column].strip() for column in columns)
                if any(value in {"", "0", "-", "N", "NA", "?"} for value in profile):
                    incomplete += 1
                    continue
                for locus, allele in zip(loci, profile, strict=True):
                    if allele not in alleles[locus]:
                        unavailable_references.add(f"{locus}_{allele}")
                if st in st_profiles and st_profiles[st] != profile:
                    raise SchemeError(f"{file.name}: ST {st} is assigned conflicting allele profiles.")
                st_profiles[st] = profile
                accumulated.setdefault(profile, set()).add(st)
        profiles = {profile: tuple(sorted(sts, key=_natural_key))
                    for profile, sts in accumulated.items()}
        if incomplete:
            notes.append(f"Ignored {incomplete} incomplete profile-table rows for ST assignment.")
        if unavailable_references:
            notes.append(
                f"Profile table references {len(unavailable_references)} alleles with no reference "
                f"sequence (including {', '.join(sorted(unavailable_references)[:5])}); "
                "STs requiring those alleles cannot be assigned."
            )
        included.append(file)
    else:
        notes.append("No ST profile table supplied; exact allele calls are available without ST assignment.")
    digest = _source_digest(included, root, cancelled, progress)
    if any(file_signature(file) != signatures[file] for file in included):
        raise SchemeError("Scheme files changed while loading; load the scheme again.")
    return Scheme(
        name=str(metadata.get("name") or root.name), path=root, loci=loci, alleles=alleles,
        profiles=profiles, metadata=metadata, digest=digest, notes=notes,
        sequences_loaded=sequences,
    )


def _automaton(scheme: Scheme, cancelled: CancelCallback, progress: ProgressCallback):
    if not scheme.sequences_loaded:
        # A scheme loaded to be checked holds allele identifiers and no bases.
        # Matching against it would find nothing and report every locus missing,
        # which is a wrong answer rather than an error, so it is refused here.
        raise SchemeError(
            f"{scheme.name} was loaded to be checked, not to be typed against, so its reference "
            "sequences were never read. Load the scheme again before matching.")
    if scheme._automaton is not None:
        return scheme._automaton
    automaton = ahocorasick.Automaton()
    for index, locus in enumerate(scheme.loci, 1):
        check_cancelled(cancelled)
        for number, (allele, sequence) in enumerate(scheme.alleles[locus].items()):
            if number % 1024 == 0:
                check_cancelled(cancelled)
            if set(sequence) - set("ACGT"):
                continue
            reverse = reverse_complement(sequence)
            strands = [(sequence, "both")] if reverse == sequence else [(sequence, "+"), (reverse, "-")]
            for word, strand in strands:
                if len(scheme.loci) > 30:
                    offset = max(0, (len(word) - 31) // 2)
                    seed = word[offset:offset + 31]
                    values = automaton.get(seed, [])
                    values.append((locus, allele, strand, len(sequence), word, offset))
                    automaton.add_word(seed, values)
                else:
                    values = automaton.get(word, [])
                    values.append((locus, allele, strand, len(sequence)))
                    automaton.add_word(word, values)
        if progress:
            progress(index, len(scheme.loci), f"Indexing locus {locus}")
    if not len(automaton):
        raise SchemeError("This scheme contains no unambiguous reference sequences for exact matching.")
    automaton.make_automaton()
    check_cancelled(cancelled)
    scheme._automaton = automaton
    return automaton


def _summarize_locus(locus: str, matches: list[dict]) -> dict:
    candidates = sorted({match["allele"] for match in matches}, key=_natural_key)
    locations = {(match["contig"], match["start"], match["end"]) for match in matches}
    status, allele, reason = "missing", None, "No complete exact reference allele match."
    if len(candidates) == 1 and len(locations) == 1:
        status, allele, reason = "exact", candidates[0], "One unambiguous exact allele match."
    elif matches:
        status, reason = "ambiguous", "Multiple exact candidates overlap, or the locus has multiple copies."
        # Different alleles at disjoint positions are a mixed-call signal, not a
        # diagnosis of contamination. Every hit participates in this decision.
        if len(candidates) > 1:
            ordered = sorted(matches, key=lambda hit: (hit["contig"], hit["start"], hit["end"]))
            earliest: dict[str, tuple[str, int]] = {}
            for hit in ordered:
                for other_allele, (other_contig, other_end) in earliest.items():
                    if other_allele != hit["allele"] and (
                        other_contig != hit["contig"] or other_end < hit["start"]
                    ):
                        status = "mixed"
                        reason = "Different exact allele candidates occur at distinct positions."
                        break
                if status == "mixed":
                    break
                previous = earliest.get(hit["allele"])
                if previous is None or (hit["contig"], hit["end"]) < previous:
                    earliest[hit["allele"]] = (hit["contig"], hit["end"])
    if status == "exact" and any(hit["shared_reference"] for hit in matches):
        status, allele = "ambiguous", None
        reason = "The matched reference sequence is identical between different loci."
    evidence_limit = 50
    return {
        "locus": locus, "allele": allele, "status": status, "candidates": candidates,
        "hit_count": len(matches), "hits": matches[:evidence_limit],
        "evidence_truncated": len(matches) > evidence_limit, "reason": reason,
    }


def call_assembly(
    path: str | Path,
    scheme: Scheme,
    cancelled: CancelCallback = None,
    progress: ProgressCallback = None,
) -> dict:
    """Call complete exact alleles on both strands of each FASTA contig.

    Known STs require every locus to be unambiguous and the full allele vector to
    map to exactly one ST. No matches are allowed across separate contigs.
    """
    path = Path(path).resolve()
    signature = file_signature(path)
    automaton = _automaton(scheme, cancelled, progress)
    matches: dict[str, list[dict]] = {locus: [] for locus in scheme.loci}
    qc = QCAccumulator()
    with SequenceReader(path, cancelled) as reader:
        if reader.kind != "fasta":
            raise SequenceError("Exact allele calling requires a FASTA assembly; FASTQ reads need assembly first.")
        for contig_number, record in enumerate(reader, 1):
            check_cancelled(cancelled)
            qc.add(record)
            # Chunks bound cancellation latency even on a very large contig with
            # no matches. The iterator carries automaton state across chunks.
            iterator = automaton.iter("")
            for start in range(0, len(record.sequence), 65536):
                check_cancelled(cancelled)
                chunk = record.sequence[start:start + 65536]
                iterator.set(chunk, False)
                for hit_number, (end, values) in enumerate(iterator):
                    if hit_number % 4096 == 0:
                        check_cancelled(cancelled)
                    shared = len({value[0] for value in values}) > 1
                    for value in values:
                        locus, allele, strand, length = value[:4]
                        hit_start, hit_end = end - length + 2, end + 1
                        reference_shared = shared
                        if len(value) == 6:
                            full_sequence, offset = value[4:]
                            seed_length = min(31, length)
                            start_zero = end - seed_length + 1 - offset
                            if start_zero < 0 or not record.sequence.startswith(full_sequence, start_zero):
                                continue
                            hit_start, hit_end = start_zero + 1, start_zero + length
                            reference_shared = shared and any(
                                other[0] != locus and other[4] == full_sequence for other in values)
                        matches[locus].append({
                            "allele": allele, "contig": record.identifier,
                            "start": hit_start, "end": hit_end,
                            "strand": strand, "shared_reference": reference_shared,
                        })
            if progress:
                progress(contig_number, 0, f"Scanned contig {record.identifier}")
    digest = file_sha256(path, cancelled)
    if file_signature(path) != signature:
        raise SequenceError(f"{path.name}: input changed during analysis; run it again.")
    calls = [_summarize_locus(locus, matches[locus]) for locus in scheme.loci]
    alleles = {call["locus"]: call["allele"] for call in calls}
    call_statuses = {call["status"] for call in calls}
    st = None
    notes = list(scheme.notes)
    notes.append("Exact nucleotide matching only; missing matches do not establish novel alleles.")
    if "mixed" in call_statuses:
        status = "mixed"
        notes.append("Distinct allele candidates were found; review assembly duplication or mixed input.")
    elif "ambiguous" in call_statuses:
        status = "ambiguous"
    elif "missing" in call_statuses:
        status = "incomplete"
    elif not scheme.profiles:
        status = "profile_unavailable"
    else:
        vector = tuple(alleles[locus] for locus in scheme.loci)
        sts = scheme.profiles.get(vector, ())
        if len(sts) == 1:
            st, status = sts[0], "complete"
        elif len(sts) > 1:
            status = "ambiguous"
            notes.append(f"This exact profile maps to multiple STs: {', '.join(sts)}.")
        else:
            status = "novel_profile"
            notes.append("All alleles are known, but their complete combination is absent from this local profile table.")
    return {
        "sample_name": sample_name(path), "input_path": str(path), "kind": "fasta",
        "qc": qc.result("fasta", True, None), "scheme": scheme.name,
        "scheme_digest": scheme.digest, "scheme_metadata": dict(scheme.metadata),
        "st": st, "status": status, "alleles": alleles, "calls": calls,
        "input_sha256": digest, "notes": notes,
        "parameters": {"method": "exact-nucleotide", "strands": "both",
                       "index": "seed-verified" if len(scheme.loci) > 30 else "full-allele",
                       "coordinates": "1-based inclusive", "evidence_limit_per_locus": 50},
        "engine_version": __version__,
    }
