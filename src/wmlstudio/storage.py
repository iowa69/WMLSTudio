"""Managed sequence copies, organism assignments, and reversible organisation.

Original user files are never renamed, moved, or overwritten. Only a verified
managed copy may be removed after its replacement and database update succeed.
"""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
import unicodedata
import uuid
from pathlib import Path

from wmlstudio.project import Project
from wmlstudio.sequence import check_cancelled, file_sha256, file_signature, sample_name

_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(10)),
             *(f"LPT{i}" for i in range(10))}


def safe_component(value: object, fallback: str = "Unknown") -> str:
    """Produce a bounded Windows-safe path component, never a traversal segment."""
    text = unicodedata.normalize("NFKC", str(value or ""))
    text = re.sub(r'[<>:"/\\|?*\x00-\x1f\x7f]', "_", text).strip(" .")
    text = re.sub(r"\s+", " ", text)
    if not text or text in {".", ".."}:
        text = fallback
    if text.split(".", 1)[0].upper() in _RESERVED:
        text = "_" + text
    return text[:72].rstrip(" .") or fallback


def _assignment(value: dict) -> tuple[dict, dict]:
    mode = value.get("typing_mode", "auto")
    if mode not in {"auto", "manual", "unknown"}:
        raise ValueError(f"Unsupported organism assignment mode: {mode}")
    organism = {"genus": str(value.get("genus", "")).strip(),
                "species": str(value.get("species", "")).strip()}
    if mode == "manual" and not organism["genus"]:
        raise ValueError("Manual organism assignment requires a genus.")
    workflow = {"typing_mode": mode, "scheme_path": value.get("scheme_path") or None}
    return organism, workflow


def assign_organism(
    project: Project, sample_ids, genus: str, species: str = "", scheme_path=None,
    typing_mode: str = "manual",
) -> None:
    organism, workflow = _assignment({"genus": genus, "species": species,
                                      "typing_mode": typing_mode, "scheme_path": scheme_path})
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            sample = project.get_sample(sample_id)
            if sample["status"] == "running":
                raise ValueError("Wait for this sample's analysis before changing its assignment.")
            old = sample["metadata"]
            changed = old.get("organism") != organism or any(
                old.get("workflow", {}).get(key) != value for key, value in workflow.items())
            project.update_metadata(sample_id, {"organism": organism, "workflow": workflow})
            if changed and sample.get("result") is not None and not sample.get("profile_only"):
                project.invalidate_result(sample_id, "Organism or typing scheme assignment changed")


def _relink_digest(sample):
    """Prefer the current scientific input fingerprint over older raw-source hashes."""
    metadata = sample.get("metadata", {})
    workflow = metadata.get("workflow", {})
    result_hash = (sample.get("result") or {}).get("input_sha256")
    assembly_hash = metadata.get("assembly", {}).get("provenance", {}).get("assembly_sha256")
    choices = [(result_hash, "result.input_sha256"), (assembly_hash, "assembly.provenance.assembly_sha256"),
               (workflow.get("managed_sha256"), "workflow.managed_sha256"),
               (workflow.get("relinked_input_sha256"), "workflow.relinked_input_sha256")]
    if workflow.get("source_kind") != "assembly":
        choices.append((workflow.get("source_sha256"), "workflow.source_sha256"))
    for digest, source in choices:
        if digest:
            if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
                raise ValueError(f"The stored {source} is not a complete SHA-256; an evidence-preserving relink is unavailable.")
            return digest.lower(), source
    raise ValueError("No complete input SHA-256 is recorded. Import the file as a new sample instead of relinking evidence.")


def relink_input(project: Project, sample_id: str, new_path, *, cancelled=None) -> Path:
    """Relink an explicitly chosen byte-identical file, preserving all profiles.

    Intended for a worker: hashes the whole file and checks it did not change.
    No copy, deletion, directory scan, retyping or result replacement occurs.
    A file outside the known sample-owned storage location becomes unmanaged;
    explicit relinking never grants permission to delete the chosen file later.
    """
    sample = project.get_sample(sample_id)
    if sample.get("profile_only") or not sample.get("input_path"):
        raise ValueError("A profile-only sample has no sequence input to relink.")
    if sample["status"] == "running":
        raise ValueError("Wait for this sample's analysis before relinking its input.")
    expected, fingerprint_source = _relink_digest(sample)
    target = Path(new_path).expanduser().resolve()
    if not target.is_file():
        raise FileNotFoundError(f"Replacement sequence input does not exist: {target}")
    check_cancelled(cancelled)
    signature = file_signature(target)
    actual = file_sha256(target, cancelled)
    check_cancelled(cancelled)
    if file_signature(target) != signature:
        raise ValueError("The replacement input changed during hashing; no relink was performed.")
    if actual != expected:
        raise ValueError("The replacement file has a different SHA-256. Existing sample evidence was not changed.")
    with project.transaction():
        check_cancelled(cancelled)
        if file_signature(target) != signature:
            raise ValueError("The replacement input changed after hashing; no relink was performed.")
        current = project.get_sample(sample_id)
        if (current["status"] == "running" or current["input_path"] != sample["input_path"]
                or _relink_digest(current)[0] != expected):
            raise ValueError("Sample input or evidence changed while hashing; retry the relink.")
        workflow = current["metadata"].get("workflow", {})
        managed = False
        storage_root = workflow.get("storage_root")
        if workflow.get("managed") and storage_root:
            root = Path(storage_root).expanduser().resolve()
            managed = (target.is_relative_to(root) and sample_id in target.relative_to(root).parts
                       and target != Path(workflow.get("source_path") or current["input_path"]).resolve())
        project.set_input_path(sample_id, target)
        project.update_metadata(sample_id, {"workflow": {
            "managed": managed, "managed_sha256": expected if managed else None,
            "storage_root": storage_root if managed else None,
            "relinked_input_sha256": expected,
        }})
        if workflow.get("source_kind") == "assembly":
            project.update_metadata(sample_id, {"assembly": {"assembly_path": str(target)}})
        project.record_history(sample_id, "input_relinked", {
            "from": sample["input_path"], "to": str(target), "sha256": expected,
            "fingerprint_source": fingerprint_source, "profiles_preserved": True,
        })
    return target


def _target(root: Path, sample_id: str, original: Path, organism: dict, st, append_st: bool) -> Path:
    genus = safe_component(organism.get("genus"), "Unknown_genus")
    species = safe_component(organism.get("species"), "Unknown_species")
    st_name = "ST_" + safe_component(st) if st else "ST_unassigned"
    identifier = safe_component(sample_id)
    # Preserve compression and format suffixes while appending ST to the stem.
    suffixes = "".join(original.suffixes[-2:] if original.suffix.lower() in {".gz", ".bz2"}
                       else original.suffixes[-1:])
    stem = original.name[:-len(suffixes)] if suffixes else original.name
    filename = safe_component(stem, "sample") + ("_ST_" + safe_component(st) if st and append_st else "")
    filename += re.sub(r'[^A-Za-z0-9.]', "_", suffixes)
    target = root / genus / species / st_name / identifier / filename
    if not target.resolve().is_relative_to(root):
        raise ValueError("Managed storage path escapes its selected root.")
    return target


def _copy_atomic(source: Path, destination: Path, cancelled=None) -> str:
    check_cancelled(cancelled)
    if destination.exists():
        raise FileExistsError(f"A managed file already exists at {destination}")
    before = file_signature(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        digest = hashlib.sha256()
        with source.open("rb") as original, tempfile.NamedTemporaryFile(
            mode="wb", prefix=".copy-", suffix=".tmp", dir=destination.parent, delete=False,
        ) as target:
            temporary = Path(target.name)
            while chunk := original.read(1024 * 1024):
                check_cancelled(cancelled)
                digest.update(chunk)
                target.write(chunk)
            target.flush()
            os.fsync(target.fileno())
        check_cancelled(cancelled)
        if file_signature(source) != before:
            raise ValueError(f"Source file changed while it was being copied: {source}")
        # Windows rename refuses an existing destination and works on NTFS and
        # removable drives. POSIX link supplies the same no-overwrite guarantee.
        if os.name == "nt":
            os.rename(temporary, destination)
        else:
            os.link(temporary, destination)
            temporary.unlink()
        temporary = None
        return digest.hexdigest()
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def import_samples(
    project: Project, assignments, storage_root=None, managed: bool = True, append_st: bool = False,
    cancelled=None, progress=None,
) -> list[str]:
    """Stage a complete import before one database transaction; cancel rolls it back.

    Assignments contain path, typing_mode, genus, species and optional scheme_path.
    The function is synchronous and intended to run in an application's worker.
    """
    assignments = list(assignments)
    root = Path(storage_root or project.path.parent / "Sequences").expanduser().resolve()
    staged = []
    created = []
    committed = False
    try:
        for index, assignment in enumerate(assignments):
            check_cancelled(cancelled)
            original = Path(assignment["path"]).expanduser().resolve()
            if not original.is_file():
                raise FileNotFoundError(f"Sequence input does not exist: {original}")
            organism, workflow = _assignment(assignment)
            sample_id = uuid.uuid4().hex
            destination = _target(root, sample_id, original, organism, None, append_st) if managed else original
            if progress:
                progress(index, len(assignments), f"Importing {original.name}")
            if managed:
                digest = _copy_atomic(original, destination, cancelled)
                created.append(destination)
            else:
                digest = file_sha256(original, cancelled)
            workflow.update({"managed": managed, "append_st": append_st,
                             "storage_root": str(root) if managed else None,
                             "source_path": str(original), "source_sha256": digest,
                             "managed_sha256": digest if managed else None})
            staged.append((sample_id, destination, assignment.get("name") or sample_name(original),
                           {"organism": organism, "workflow": workflow}))
        check_cancelled(cancelled)
        with project.transaction():
            for sample_id, destination, name, metadata in staged:
                check_cancelled(cancelled)
                project.add_sample(destination, name, sample_id=sample_id)
                project.set_metadata(sample_id, metadata)
        committed = True
        if progress:
            progress(len(assignments), len(assignments), "Import complete")
        return [entry[0] for entry in staged]
    except BaseException:
        if not committed:
            for destination in created:
                destination.unlink(missing_ok=True)
        raise


def import_managed(
    project: Project, paths, storage_root, genus: str = "", species: str = "",
    append_st: bool = False, cancelled=None,
) -> list[str]:
    mode = "manual" if genus else "auto"
    return import_samples(project, [{"path": path, "genus": genus, "species": species,
                                    "typing_mode": mode} for path in paths],
                          storage_root, append_st=append_st, cancelled=cancelled)


def organize_sample(
    project: Project, sample_id: str, storage_root=None, append_st: bool | None = None,
    cancelled=None,
) -> Path:
    """Copy to the assigned genus/species/ST location, then update the project.

    Only an existing managed copy is eligible for cleanup. Original paths remain
    in workflow metadata. If cleanup fails, its path is recorded for recovery.
    """
    sample = project.get_sample(sample_id)
    metadata = sample["metadata"]
    workflow = metadata.get("workflow", {})
    source = Path(sample["input_path"]).resolve()
    root = Path(storage_root or workflow.get("storage_root") or
                project.path.parent / "Sequences").expanduser().resolve()
    append = bool(workflow.get("append_st", False)) if append_st is None else append_st
    result = sample.get("result") or {}
    st = result.get("st") if sample["status"] == "completed" and result.get("status") == "complete" else None
    original = Path(workflow.get("source_path") or source)
    organism = dict(metadata.get("organism", {}))
    identified = result.get("identification", {}) or {}
    if not organism.get("genus"):
        detected = identified.get("organism") or {}
        detected = detected if isinstance(detected, dict) else {}
        organism = {"genus": detected.get("genus") or identified.get("genus", ""),
                    "species": detected.get("species") or identified.get("species", "")}
    # An assembled isolate must never inherit its original FASTQ extension.
    naming_source = (Path(safe_component(sample["name"]) + ".fasta")
                     if workflow.get("source_kind") == "assembly" else original)
    destination = _target(root, sample_id, naming_source, organism, st, append)
    check_cancelled(cancelled)
    if source == destination:
        return destination
    created = False
    if destination.exists():
        digest = file_sha256(source, cancelled)
        if file_sha256(destination, cancelled) != digest:
            raise FileExistsError("A different file already occupies the requested managed destination.")
    else:
        digest = _copy_atomic(source, destination, cancelled)
        created = True
    try:
        check_cancelled(cancelled)
        if result.get("input_sha256") and result["input_sha256"] != digest:
            raise ValueError("Input contents changed since typing; rerun before organising by ST.")
        with project.transaction():
            # An input cannot be changed behind a running analysis.
            current = project.get_sample(sample_id)
            if current["status"] == "running" or current["input_path"] != str(source):
                raise ValueError("Sample input changed or analysis is running; retry organisation later.")
            project.set_input_path(sample_id, destination)
            project.update_metadata(sample_id, {"workflow": {
                "managed": True, "storage_root": str(root), "append_st": append,
                "source_path": str(original), "source_sha256": workflow.get("source_sha256") or digest,
                "managed_sha256": digest,
            }})
            if workflow.get("source_kind") == "assembly":
                project.update_metadata(sample_id, {"assembly": {"assembly_path": str(destination)}})
    except BaseException:
        if created:
            destination.unlink(missing_ok=True)
        raise
    previous_root = workflow.get("storage_root")
    if (workflow.get("managed") and previous_root and source != original.resolve()
            and source.is_relative_to(Path(previous_root).resolve())):
        try:
            if file_sha256(source) == workflow.get("managed_sha256"):
                source.unlink()
        except OSError:
            project.record_history(sample_id, "managed_copy_retained", {"path": str(source)})
    return destination
