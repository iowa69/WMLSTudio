"""Explicit assembly/read associations without replacing any sample or evidence.

Filename suggestions are only review aids. Full pair structure/QC and SHA-256
are validated in a background task; biological identity remains user-confirmed.
"""

from __future__ import annotations

import copy
import hashlib
import re
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio.pairing_dialog import mate_hint
from wmlstudio.sequence import (
    QCAccumulator,
    SequenceError,
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
    read_pair_identity,
)


def _assembly_key(path):
    name = Path(path).name.casefold().removesuffix(".gz").removesuffix(".bz2")
    name = re.sub(r"\.(?:fasta|fna|fa|fas)$", "", name)
    return re.sub(r"[_.-](?:contigs|assembly|assembled|scaffolds)$", "", name)


def suggest_read_attachments(assemblies, reads):
    """Return unique exact-name suggestions, withholding every ambiguous mapping."""
    suggestions, ambiguous, missing = [], {}, []
    for assembly in assemblies:
        workflow = assembly.get("metadata", {}).get("workflow", {})
        keys = {_assembly_key(path) for path in (assembly.get("input_path"), workflow.get("source_path")) if path}
        matching = {1: [], 2: []}
        for read in reads:
            paths = [read.get("input_path"), read.get("metadata", {}).get("workflow", {}).get("source_path")]
            hints = {mate_hint(path) for path in paths if path}
            for direction in (1, 2):
                if any((key, direction) in hints for key in keys):
                    matching[direction].append(read)
        if len(matching[1]) == len(matching[2]) == 1:
            suggestions.append({"sample_id": assembly["id"], "read1_id": matching[1][0]["id"],
                                "read2_id": matching[2][0]["id"], "basis": "Exact normalized filename; user confirmation required"})
        elif matching[1] or matching[2]:
            ambiguous[assembly["id"]] = {str(key): [item["id"] for item in value] for key, value in matching.items()}
        else:
            missing.append(assembly["id"])
    for read_id in {item[key] for item in suggestions for key in ("read1_id", "read2_id")}:
        conflicts = [item for item in suggestions if read_id in (item["read1_id"], item["read2_id"])]
        if len(conflicts) > 1:
            for item in conflicts:
                ambiguous[item["sample_id"]] = {"conflicting_read": read_id}
    return {"suggestions": [item for item in suggestions if item["sample_id"] not in ambiguous],
            "ambiguous": ambiguous, "missing": missing}


def validate_read_attachment(sample, read1, read2, *, read_sample_ids=None, cancelled=None, progress=None):
    """Validate complete FASTQ pairs/QC without trimming, staging or changing reads."""
    paths = [Path(path).expanduser().resolve() for path in (sample["input_path"], read1, read2)]
    if any(not path.is_file() for path in paths):
        raise SequenceError("The assembly and both original FASTQ files must exist.")
    if any(paths[left].samefile(paths[right]) for left, right in ((0, 1), (0, 2), (1, 2))):
        raise SequenceError("Choose three different files: one assembly and two FASTQ mates.")
    signatures = [file_signature(path) for path in paths]
    with SequenceReader(paths[0], cancelled) as assembly:
        if assembly.kind != "fasta":
            raise SequenceError("Attach reads to a FASTA assembly, not another read file.")
        if not sum(1 for _ in assembly):
            raise SequenceError("The assembly contains no valid FASTA records.")
    accumulators, count, explicit = [QCAccumulator(), QCAccumulator()], 0, True
    with SequenceReader(paths[1], cancelled) as first, SequenceReader(paths[2], cancelled) as second:
        if first.kind != "fastq" or second.kind != "fastq":
            raise SequenceError("Read attachment requires two FASTQ inputs.")
        left, right = iter(first), iter(second)
        while True:
            check_cancelled(cancelled)
            a, b = next(left, None), next(right, None)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise SequenceError("FASTQ mates have different numbers of records.")
            aid, am = read_pair_identity(a.name)
            bid, bm = read_pair_identity(b.name)
            if aid != bid or am not in (None, 1) or bm not in (None, 2) or (am is None) != (bm is None):
                raise SequenceError(f"FASTQ identifiers or mate directions disagree at pair {count + 1}.")
            for accumulator, record in zip(accumulators, (a, b), strict=True):
                accumulator.add(record)
            count += 1
            explicit = explicit and am == 1 and bm == 2
            if progress and count % 10000 == 0:
                progress(0, 1, f"Validating complete pairing · {count:,} read pairs")
    if not count:
        raise SequenceError("FASTQ attachment contains no read pairs.")
    hashes = [file_sha256(path, cancelled) for path in paths]
    if signatures != [file_signature(path) for path in paths]:
        raise SequenceError("An input changed during read-attachment validation; nothing was linked.")
    source_ids = read_sample_ids or [None, None]
    if len(source_ids) != 2:
        raise ValueError("Supply exactly two optional source-read sample IDs.")
    identity = hashlib.sha256((sample["id"] + "\0" + "\0".join(hashes[1:])).encode()).hexdigest()
    return {"attachment_id": identity, "sample_id": sample["id"], "assembly_path": str(paths[0]),
            "assembly_sha256_at_link": hashes[0], "assembly_signature": list(signatures[0]),
            "reads": [{"read_id": hashlib.sha256((identity + str(index)).encode()).hexdigest(),
                       "path": str(path), "sha256": digest, "signature": list(signature),
                       "mate": index, "source_sample_id": source_ids[index - 1],
                       "qc": accumulators[index - 1].result("fastq", True, None)}
                      for index, (path, digest, signature) in enumerate(zip(paths[1:], hashes[1:], signatures[1:], strict=True), 1)],
            "pairing": {"records_checked": count, "complete_file": True, "sampled": False,
                        "explicit_mates": explicit, "status": "verified" if explicit else "unmarked"},
            "linked_at": datetime.now(UTC).isoformat(), "identity_basis": "Explicit user association; filenames/read IDs do not prove assembly provenance",
            "preprocessing": "None; originals unchanged. Native structural/Phred+33 QC, not FastQC or fastp."}


def attach_read_pair(project, sample_id, evidence):
    """Serial metadata-only commit; preserve assembly input/result and all stable IDs."""
    if evidence.get("sample_id") != sample_id:
        raise ValueError("Read evidence belongs to a different stable sample ID.")
    reads = evidence.get("reads", [])
    pairing = evidence.get("pairing", {})
    if (len(reads) != 2 or [read.get("mate") for read in reads] != [1, 2]
            or not pairing.get("complete_file") or pairing.get("sampled")
            or pairing.get("records_checked", 0) < 1
            or any(not re.fullmatch("[0-9a-f]{64}", read.get("sha256", "")) for read in reads)):
        raise ValueError("Read attachment requires complete validated pairing and both SHA-256 identities.")
    identity = hashlib.sha256((sample_id + "\0" + "\0".join(read["sha256"] for read in reads)).encode()).hexdigest()
    if evidence.get("attachment_id") != identity:
        raise ValueError("Read-attachment identity does not match the validated files.")
    sample = project.get_sample(sample_id)
    if Path(sample["input_path"]).resolve() != Path(evidence["assembly_path"]).resolve():
        raise ValueError("The assembly input changed while reads were being validated.")
    records = [{"path": evidence["assembly_path"], "signature": evidence["assembly_signature"]}, *evidence["reads"]]
    for record in records:
        if list(file_signature(record["path"])) != record["signature"]:
            raise ValueError("A validated input changed before attachment; validate it again.")
    for read in evidence["reads"]:
        if read.get("source_sample_id"):
            source = project.get_sample(read["source_sample_id"])
            if Path(source["input_path"]).resolve() != Path(read["path"]).resolve():
                raise ValueError("A source read record changed during validation.")
            linked = source.get("metadata", {}).get("workflow", {}).get("paired_with")
            if linked and linked != sample_id:
                raise ValueError("A selected read record is already attached to another isolate.")
    with project.transaction():
        current = project.get_sample(sample_id)
        if current["input_path"] != sample["input_path"]:
            raise ValueError("Assembly changed before the read association was saved.")
        metadata = current["metadata"]
        previous = metadata.get("reads")
        if previous and previous.get("attachment_id") == evidence["attachment_id"]:
            return previous["attachment_id"]
        if previous:
            project.record_history(sample_id, "read_attachment_superseded", {"evidence": previous})
        metadata["reads"] = copy.deepcopy(evidence)
        project.set_metadata(sample_id, metadata)
        project.record_history(sample_id, "reads_attached", {"attachment_id": evidence["attachment_id"],
                               "assembly_sha256": evidence["assembly_sha256_at_link"], "reads": evidence["reads"]})
        for read in reads:
            source_id = read.get("source_sample_id")
            if source_id:
                project.update_metadata(source_id, {"workflow": {"source_kind": "read_mate", "paired_with": sample_id,
                    "read_attachment_id": identity, "mate": read["mate"]}})
                project.record_history(source_id, "read_attached_to_assembly", {
                    "sample_id": sample_id, "attachment_id": identity, "mate": read["mate"], "sha256": read["sha256"]})
    return evidence["attachment_id"]
