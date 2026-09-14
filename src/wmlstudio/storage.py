"""Managed sequence copies, organism assignments, and reversible organisation.

Original user files are never renamed, moved, or overwritten. Only a verified
managed copy may be removed after its replacement and database update succeed.

A folder path is a filing decision, never a taxonomic claim. When identification
declines to decide, the file is filed under the quarantine tree rather than being
pushed into a genus folder the evidence does not support.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
import unicodedata
import uuid
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio.project import Project
from wmlstudio.sequence import (
    AnalysisCancelled,
    check_cancelled,
    file_sha256,
    file_signature,
    sample_name,
)

_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(10)),
             *(f"LPT{i}" for i in range(10))}

# "The app declined to decide" is a different statement from the user's explicit
# "I do not know", which keeps using Unknown_genus/Unknown_species. A leading
# underscore survives safe_component and no PubMLST scheme name may start with
# one, so this bucket can never collide with a real genus folder.
QUARANTINE_ROOT = "_Unresolved"
QUARANTINE_BUCKETS = {
    "awaiting_identification": "Awaiting_identification",
    "low_confidence": "Low_confidence",
    "conflicting_evidence": "Conflicting_evidence",
    "reads_not_assembled": "Reads_not_assembled",
    "not_in_reference_panel": "Not_in_reference_panel",
    "user_deferred": "User_deferred",
}
QUARANTINE_README = """Needs review — these files are not filed by organism.

WMLSTudio put a file here because it declined to decide what the organism is,
not because the organism is new, unusual, or wrong. A folder name in this program
is where a copy is stored; it is never a laboratory identification.

  Awaiting_identification  No organism evidence has been reviewed yet.
  Low_confidence           Evidence exists but is below the configured filing floor.
  Conflicting_evidence     Competing references were too close, or a second organism
                           covered a substantial part of the assembly.
  Reads_not_assembled      Raw reads cannot be identified until they are assembled.
  Not_in_reference_panel   No installed reference was close enough. This is a limit of
                           the installed panel, not evidence of a novel species.
  User_deferred            You chose to review this file later.

Assign an organism in WMLSTudio to move a file out of this folder. Your original
files are untouched; everything here is a managed copy.
"""
ORGANISM_EVIDENCE_VERSION = 1
_EVIDENCE_STATUS = {"proposed", "confirmed", "quarantined"}
_EVIDENCE_BYTES = 64 * 1024

# What a project files without asking. Only a whole-genome comparison against the
# installed panel is strong enough to move a file into a genus folder on its own;
# everything weaker is imported all the same and waits in the needs-review tree.
INTAKE_DEFAULT_POLICY = {"auto_confirm": True, "min_confidence": "genomic_reference_supported"}


def default_storage_root(project: Project) -> Path:
    """The managed folder that belongs to this project file, beside it on disk.

    Derived from the project's own current path, so a project that moved with its
    folder computes the folder it moved to rather than the one it was created in.
    """
    return project.path.with_suffix(".files")


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


def quarantine_token(value) -> str | None:
    """Validate a caller-supplied quarantine reason; unknown tokens are refused."""
    if value in (None, "", False):
        return None
    token = str(value)
    if token not in QUARANTINE_BUCKETS:
        raise ValueError(f"Unknown quarantine bucket: {token}")
    return token


def organism_evidence(value) -> dict | None:
    """Validate a decision record before it is written into sample metadata.

    The record states which engine produced a label and how strong that label is.
    It is bounded on purpose: the full ANI and MLST payloads stay where the
    engines wrote them, and this is only the filing decision.
    """
    if value in (None, {}):
        return None
    if not isinstance(value, dict):
        raise TypeError("Organism evidence must be a mapping.")
    evidence = dict(value)
    evidence.setdefault("format_version", ORGANISM_EVIDENCE_VERSION)
    if evidence["format_version"] != ORGANISM_EVIDENCE_VERSION:
        raise ValueError("Unsupported organism evidence format version.")
    status = evidence.setdefault("status", "proposed")
    if status not in _EVIDENCE_STATUS:
        raise ValueError(f"Unsupported organism evidence status: {status}")
    quarantine_token(evidence.get("quarantine_reason"))
    encoded = json.dumps(evidence, sort_keys=True, ensure_ascii=False, allow_nan=False)
    if len(encoded.encode()) > _EVIDENCE_BYTES:
        raise ValueError("Organism evidence exceeds the 64 KiB decision-record bound.")
    return evidence


def review_bucket(metadata: dict) -> str | None:
    """The needs-review folder a stored decision record implies, or None when filed.

    Only an accepted organism leaves the review tree. A proposal nobody has
    accepted yet is not a decision, so it keeps waiting under its own bucket
    rather than drifting into a genus folder the next time files are re-filed.
    An unknown reason still needs review; it never becomes "no reason".
    """
    metadata = metadata if isinstance(metadata, dict) else {}
    evidence = metadata.get("organism_evidence")
    if not isinstance(evidence, dict):
        return None
    if evidence.get("status") == "confirmed":
        return None
    if evidence.get("status") not in {"quarantined", "proposed"}:
        return None
    token = str(evidence.get("quarantine_reason") or "awaiting_identification")
    return token if token in QUARANTINE_BUCKETS else "awaiting_identification"


def assign_organism(
    project: Project, sample_ids, genus: str, species: str = "", scheme_path=None,
    typing_mode: str = "manual",
) -> None:
    """Record an organism assignment, invalidating only a result it can invalidate.

    Under typing_mode 'auto' the scheme is derived from the organism, so changing
    the organism can change which panel applies and the result must be re-earned.
    When the user pinned an explicit scheme_path in manual mode, the ST belongs to
    that scheme's digest rather than to the display label, so correcting a
    mis-called genus must not destroy a valid, still-applicable result.
    """
    organism, workflow = _assignment({"genus": genus, "species": species,
                                      "typing_mode": typing_mode, "scheme_path": scheme_path})
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            sample = project.get_sample(sample_id)
            if sample["status"] == "running":
                raise ValueError("Wait for this sample's analysis before changing its assignment.")
            old = sample["metadata"]
            old_workflow = old.get("workflow", {})
            scheme_bound = (old_workflow.get("typing_mode") == "manual"
                            and bool(old_workflow.get("scheme_path")))
            workflow_changed = any(old_workflow.get(key) != value for key, value in workflow.items())
            organism_changed = old.get("organism") != organism
            changed = workflow_changed or (organism_changed and not scheme_bound)
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


def _target(root: Path, sample_id: str, original: Path, organism: dict, st, append_st: bool,
            *, quarantine: str | None = None) -> Path:
    identifier = safe_component(sample_id)
    if quarantine:
        st = None  # A quarantined file has no organism, so it has no ST folder.
    # Preserve compression and format suffixes while appending ST to the stem.
    suffixes = "".join(original.suffixes[-2:] if original.suffix.lower() in {".gz", ".bz2"}
                       else original.suffixes[-1:])
    stem = original.name[:-len(suffixes)] if suffixes else original.name
    filename = safe_component(stem, "sample") + ("_ST_" + safe_component(st) if st and append_st else "")
    filename += re.sub(r'[^A-Za-z0-9.]', "_", suffixes)
    if quarantine:
        bucket = QUARANTINE_BUCKETS.get(quarantine)
        if bucket is None:
            raise ValueError(f"Unknown quarantine bucket: {quarantine}")
        target = root / QUARANTINE_ROOT / bucket / identifier / filename
    else:
        genus = safe_component(organism.get("genus"), "Unknown_genus")
        if genus.casefold() == QUARANTINE_ROOT.casefold():
            # No organism is named after the needs-review folder; keep them apart
            # so a browsable quarantine tree can never be mistaken for a taxon.
            genus += "_genus"
        species = safe_component(organism.get("species"), "Unknown_species")
        st_name = "ST_" + safe_component(st) if st else "ST_unassigned"
        target = root / genus / species / st_name / identifier / filename
    if not target.resolve().is_relative_to(root):
        raise ValueError("Managed storage path escapes its selected root.")
    return target


def preview_target(root, sample_id: str, original, genus: str = "", species: str = "", *,
                   st=None, append_st: bool = False, quarantine=None) -> Path:
    """The exact path an import would create, computed without touching the disk.

    A review dialog renders its "goes to" column from here, so what the user is
    shown and what import_samples creates are produced by the same code.
    """
    return _target(Path(root).expanduser().resolve(), sample_id, Path(original),
                   {"genus": genus, "species": species}, st, append_st,
                   quarantine=quarantine_token(quarantine))


def _quarantine_notice(root: Path) -> None:
    """Explain the quarantine tree once, in the folder the user will open."""
    notice = root / QUARANTINE_ROOT / "README.txt"
    if notice.exists():
        return
    try:
        notice.parent.mkdir(parents=True, exist_ok=True)
        with notice.open("x", encoding="utf-8", newline="\r\n") as handle:
            handle.write(QUARANTINE_README)
    except OSError:
        pass  # An explanatory file is never worth failing an import for.


def _prune_empty(root: Path, directory: Path) -> None:
    """Remove directories vacated by a move, upward to but never including root."""
    root = Path(root).resolve()
    try:
        current = Path(directory).resolve()
    except OSError:
        return
    while current != root and current.is_relative_to(root):
        if current.is_symlink() or not current.is_dir():
            return
        try:
            if any(current.iterdir()):
                return
            current.rmdir()
        except OSError:
            return
        current = current.parent


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


def known_digests(project: Project) -> dict[str, str]:
    """Every input fingerprint this project already holds, mapped to its isolate.

    Both the user's original and the managed copy are recorded, so a file
    re-offered from either location is recognised as one this project already has.
    """
    known: dict[str, str] = {}
    for sample in project.samples():
        workflow = (sample.get("metadata") or {}).get("workflow", {})
        for digest in (workflow.get("source_sha256"), workflow.get("managed_sha256")):
            if isinstance(digest, str) and digest:
                known.setdefault(digest.lower(), sample["id"])
    return known


def sample_for_digest(project: Project, sha256: str) -> str | None:
    """The isolate this project already holds for these exact bytes, or None.

    Identity is the file's SHA-256, never its name: the same assembly offered
    twice under two names is one isolate, and two different files that happen to
    share a name are not.
    """
    if not isinstance(sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", sha256):
        raise ValueError("A duplicate check needs a complete SHA-256.")
    return known_digests(project).get(sha256.lower())


def import_samples(
    project: Project, assignments, storage_root=None, managed: bool = True, append_st: bool = False,
    cancelled=None, progress=None, *, duplicates: str = "allow", notes=None,
) -> list[str]:
    """Stage a complete import before one database transaction; cancel rolls it back.

    Assignments contain path, typing_mode, genus, species and optional scheme_path,
    quarantine and organism_evidence. Identification is expected to have already
    run against the user's original file, so a reviewed verdict is filed once and
    a misidentified file is never copied into a genus folder and moved back.
    The function is synchronous and intended to run in an application's worker.

    duplicates='skip' drops a file whose bytes already exist in this project and
    appends a record of it to notes; 'allow' keeps every file and is the default.
    """
    if duplicates not in {"allow", "skip"}:
        raise ValueError(f"Unsupported duplicate import policy: {duplicates}")
    assignments = list(assignments)
    root = Path(storage_root or default_storage_root(project)).expanduser().resolve()
    known = known_digests(project) if duplicates == "skip" else {}
    staged = []
    created = []
    skipped = []
    committed = False
    try:
        for index, assignment in enumerate(assignments):
            check_cancelled(cancelled)
            original = Path(assignment["path"]).expanduser().resolve()
            if not original.is_file():
                raise FileNotFoundError(f"Sequence input does not exist: {original}")
            organism, workflow = _assignment(assignment)
            quarantine = quarantine_token(assignment.get("quarantine"))
            evidence = organism_evidence(assignment.get("organism_evidence"))
            sample_id = uuid.uuid4().hex
            destination = (_target(root, sample_id, original, organism, None, append_st,
                                   quarantine=quarantine) if managed else original)
            if progress:
                progress(index, len(assignments), f"Importing {original.name}")
            if managed:
                if quarantine:
                    _quarantine_notice(root)
                digest = _copy_atomic(original, destination, cancelled)
            else:
                digest = file_sha256(original, cancelled)
            if duplicates == "skip" and digest.lower() in known:
                if managed:
                    destination.unlink(missing_ok=True)
                    _prune_empty(root, destination.parent)
                skipped.append({"path": str(original), "sha256": digest,
                                "duplicate_of": known[digest.lower()]})
                continue
            if managed:
                created.append(destination)
            workflow.update({"managed": managed, "append_st": append_st,
                             "storage_root": str(root) if managed else None,
                             "source_path": str(original), "source_sha256": digest,
                             "managed_sha256": digest if managed else None})
            metadata = {"organism": organism, "workflow": workflow}
            if evidence is not None:
                metadata["organism_evidence"] = evidence
            staged.append((sample_id, destination, assignment.get("name") or sample_name(original),
                           metadata, quarantine, evidence))
            known.setdefault(digest.lower(), sample_id)
        check_cancelled(cancelled)
        with project.transaction():
            for sample_id, destination, name, metadata, quarantine, evidence in staged:
                check_cancelled(cancelled)
                project.add_sample(destination, name, sample_id=sample_id)
                project.set_metadata(sample_id, metadata)
                if evidence is not None:
                    project.record_history(sample_id, "organism_identified", {
                        "basis": evidence.get("basis"), "confidence": evidence.get("confidence"),
                        "status": evidence.get("status"), "proposed": evidence.get("proposed")})
                if quarantine:
                    project.record_history(sample_id, "organism_quarantined", {
                        "reason": quarantine, "bucket": QUARANTINE_BUCKETS[quarantine],
                        "path": str(destination)})
            for entry in skipped:
                project.record_history(None, "duplicate_import_skipped", entry)
        committed = True
        if notes is not None:
            notes.extend(skipped)
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


def filing_policy(project: Project | None = None, override=None) -> dict:
    """What this project files by itself and what it holds back for a person.

    A project that has never been given a policy files only whole-genome
    comparison supported calls without asking. Everything weaker — a genus-only
    reference, a complex this panel cannot split, an MLST panel match, nothing at
    all — is still imported, and waits in the needs-review tree under its reason.
    """
    from wmlstudio.organism_id import (
        SETTING_AUTO_CONFIRM,
        SETTING_MIN_CONFIDENCE,
        resolve_policy,
    )
    if override is not None:
        return resolve_policy(override)
    stored = {}
    if project is not None:
        for key, field in ((SETTING_AUTO_CONFIRM, "auto_confirm"),
                           (SETTING_MIN_CONFIDENCE, "min_confidence")):
            value = project.get_setting(key, None)
            if value is not None:
                stored[field] = value
    return resolve_policy({**INTAKE_DEFAULT_POLICY, **stored})


def set_filing_policy(project: Project, *, auto_confirm=None, min_confidence=None) -> dict:
    """Record the filing policy on the project and return what now applies.

    Raising the floor never retroactively re-files anything: it decides what the
    next intake may file without a person looking at it.
    """
    from wmlstudio.organism_id import SETTING_AUTO_CONFIRM, SETTING_MIN_CONFIDENCE
    requested = dict(filing_policy(project))
    if auto_confirm is not None:
        requested["auto_confirm"] = bool(auto_confirm)
    if min_confidence is not None:
        requested["min_confidence"] = str(min_confidence)
    effective = filing_policy(None, requested)
    with project.transaction():
        project.set_setting(SETTING_AUTO_CONFIRM, effective["auto_confirm"])
        project.set_setting(SETTING_MIN_CONFIDENCE, effective["min_confidence"])
        project.record_history(None, "filing_policy_changed", effective)
    return effective


def _unidentified(path: Path, name, reason: str, bucket: str = "awaiting_identification") -> dict:
    """An assignment for a file nothing has decided about yet; it still imports."""
    assignment = {"path": str(path), "genus": "", "species": "", "typing_mode": "auto",
                  "quarantine": bucket,
                  "organism_evidence": {"status": "quarantined", "quarantine_reason": bucket,
                                        "basis": "none", "confidence": "unresolved",
                                        "proposed": {"genus": "", "species": ""},
                                        "accepted": {"genus": "", "species": ""},
                                        "reason": reason}}
    if name:
        assignment["name"] = name
    return assignment


def _phase_progress(progress, base: int, span: int, total: int):
    """Report one stage of intake against one overall bar."""
    if progress is None:
        return None

    def report(done, count, message):
        scaled = base + (span * int(done)) // max(int(count), 1)
        progress(min(scaled, total), total, message)
    return report


def _identification_report(identify, species_panel_root, kpsc_panel_root, verdicts,
                           identifiers, needs_review) -> dict:
    """State plainly whether identification could run, and what it could not do.

    A missing reference panel is a provisioning fact the user can act on. It is
    reported as one, rather than leaving every sample sitting in needs-review with
    no explanation of why nothing was proposed.
    """
    panels = [str(root) for root in (species_panel_root, kpsc_panel_root) if root]
    errors, seen = [], set()
    for verdict in verdicts:
        for error in verdict.get("errors") or []:
            message = f"{error.get('tier')}: {error.get('message')}"
            if message not in seen:
                seen.add(message)
                errors.append(message)
    available = bool(identify and panels)
    if not identify:
        notice = ("Identification was not run, so every imported file is waiting in Needs "
                  "review until an organism is decided.")
    elif not panels:
        notice = ("No species reference panel is installed, so no genomic comparison was "
                  "attempted. Every file was imported and is waiting in Needs review. Install "
                  "the species reference panel, then identify these samples again to file them "
                  "by organism.")
    elif needs_review:
        notice = (f"{len(needs_review)} of {len(identifiers)} imported files are waiting in "
                  "Needs review: the installed references did not support an organism firmly "
                  "enough to file them, so none was filed into a genus folder on a guess.")
    else:
        notice = ""
    return {"available": available, "attempted": bool(identify), "panels": panels,
            "errors": errors, "notice": notice}


def intake_samples(
    project: Project, paths, *, storage_root=None, species_panel_root=None, kpsc_panel_root=None,
    scheme_paths=(), policy=None, names=None, managed: bool = True, append_st: bool = False,
    identify: bool = True, cancelled=None, progress=None,
) -> dict:
    """Load files once: identify them, then file and record each of them exactly once.

    This is the whole of "add samples to the project". Identification runs on the
    user's own files before a byte is copied, the verdict decides the folder, and
    the copy is made straight into that folder, so nothing is imported to one place
    and moved to another. A file this cannot identify is imported all the same and
    waits under ``_Unresolved`` with the reason it is there; it is never blocked,
    and it is never pushed into a genus folder the evidence does not support.

    Bytes already in the project are recognised and skipped, so offering the same
    file again — from the same folder or a copy of it — never produces a second
    isolate. Returns a report naming what was imported, what was skipped and which
    isolate it duplicates, what is waiting for review, and whether identification
    could run at all, so a caller can say why rather than appearing to do nothing.
    """
    from wmlstudio.organism_id import assignment_for, identify_batch
    inputs = []
    for value in paths:
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Sequence input does not exist: {path}")
        inputs.append(path)
    labels = {str(Path(key).expanduser().resolve()): value
              for key, value in dict(names or {}).items()}
    effective = filing_policy(project, policy)
    root = Path(storage_root or default_storage_root(project)).expanduser().resolve()
    total = max(2 * len(inputs), 1)
    if not inputs:
        return {"imported": [], "skipped": [], "verdicts": [], "filed": {}, "needs_review": [],
                "policy": effective, "storage_root": str(root),
                "identification": {"available": False, "attempted": False, "panels": [],
                                   "errors": [], "notice": "No files were offered."}}

    known = known_digests(project)
    skipped: list[dict] = []
    representatives: list[Path] = []
    digests: dict[str, str] = {}
    twins: dict[str, list[Path]] = {}
    for index, path in enumerate(inputs):
        check_cancelled(cancelled)
        if progress:
            progress(index, total, f"Checking {path.name}")
        digest = file_sha256(path, cancelled).lower()
        if digest in known:
            skipped.append({"path": str(path), "sha256": digest, "duplicate_of": known[digest]})
        elif digest in digests:
            twins.setdefault(digest, []).append(path)
        else:
            digests[digest] = str(path)
            representatives.append(path)

    verdicts: list[dict] = []
    assignments: list[dict] = []
    if representatives and identify:
        verdicts = identify_batch(
            representatives, species_panel_root=species_panel_root,
            kpsc_panel_root=kpsc_panel_root, scheme_paths=list(scheme_paths), cancelled=cancelled,
            progress=_phase_progress(progress, 0, len(inputs), total))
    by_path = {str(Path(verdict["input_path"]).expanduser().resolve()): verdict
               for verdict in verdicts}
    for path in representatives:
        verdict = by_path.get(str(path))
        if verdict is None:
            assignments.append(_unidentified(
                path, labels.get(str(path)),
                "Identification was not run for this file, so no organism evidence was reviewed."
                if not identify else
                "Identification did not complete for this file, so no organism was proposed."))
        else:
            assignments.append(assignment_for(verdict, policy=effective,
                                              name=labels.get(str(path))))

    notes: list[dict] = []
    identifiers = import_samples(
        project, assignments, root, managed=managed, append_st=append_st, cancelled=cancelled,
        progress=_phase_progress(progress, len(inputs), len(inputs), total),
        duplicates="skip", notes=notes)
    skipped.extend(notes)

    filed, needs_review = {}, []
    imported_by_digest = {}
    for sample_id in identifiers:
        sample = project.get_sample(sample_id)
        filed[sample_id] = sample["input_path"]
        workflow = sample["metadata"].get("workflow", {})
        digest = str(workflow.get("source_sha256") or "").lower()
        if digest:
            imported_by_digest[digest] = sample_id
        if review_bucket(sample["metadata"]):
            needs_review.append(sample_id)
    for digest, paths_for_digest in twins.items():
        for path in paths_for_digest:
            skipped.append({"path": str(path), "sha256": digest,
                            "duplicate_of": imported_by_digest.get(digest, "")})

    identification = _identification_report(
        identify, species_panel_root, kpsc_panel_root, verdicts, identifiers, needs_review)
    project.record_history(None, "samples_intake", {
        "policy": effective, "storage_root": str(root), "imported": len(identifiers),
        "skipped": len(skipped), "needs_review": len(needs_review),
        "identification_available": identification["available"]})
    if progress:
        progress(total, total, "Import complete")
    return {"imported": identifiers, "skipped": skipped, "verdicts": verdicts,
            "filed": filed, "needs_review": needs_review, "policy": effective,
            "storage_root": str(root), "identification": identification}


def _reidentification_targets(project: Project, sample_ids):
    """Which samples a fresh identification may act on, and why the rest may not."""
    samples = {sample["id"]: sample for sample in project.samples()}
    requested = list(dict.fromkeys(sample_ids)) if sample_ids is not None else None
    for sample_id in requested or ():
        if sample_id not in samples:
            raise KeyError(f"Unknown sample: {sample_id}")
    targets, skipped = [], []
    for sample in samples.values():
        if requested is not None and sample["id"] not in requested:
            continue
        waiting = review_bucket(sample["metadata"])
        if sample.get("profile_only") or not sample.get("input_path"):
            reason = "A profile-only sample has no sequence file to identify."
        elif sample["status"] == "running":
            reason = "Wait for this sample's analysis before identifying it again."
        elif not Path(sample["input_path"]).is_file():
            reason = "This sample's sequence file is not where the project recorded it."
        elif waiting is None:
            reason = ("This sample's organism was already accepted. Change it with a manual "
                      "reassignment rather than a fresh identification.")
        else:
            targets.append(sample)
            continue
        if requested is not None:
            skipped.append((sample["id"], reason))
    return targets, skipped


def reidentify_samples(
    project: Project, sample_ids=None, *, species_panel_root=None, kpsc_panel_root=None,
    scheme_paths=(), policy=None, storage_root=None, cancelled=None, progress=None,
) -> dict:
    """Identify samples that are waiting for review again, then file what is decided.

    This is what "install the reference panel and try again" means: the files are
    already in the project, so nothing is re-imported and nothing is copied twice.
    A sample whose organism a person already accepted is left alone — a fresh
    comparison never overwrites a decision somebody made.
    """
    from wmlstudio.organism_id import apply_policy, identify_batch, proposed_destination
    effective = filing_policy(project, policy)
    targets, skipped = _reidentification_targets(project, sample_ids)
    paths = [Path(sample["input_path"]) for sample in targets]
    verdicts = identify_batch(paths, species_panel_root=species_panel_root,
                              kpsc_panel_root=kpsc_panel_root, scheme_paths=list(scheme_paths),
                              cancelled=cancelled, progress=progress) if paths else []
    by_path = {str(Path(verdict["input_path"]).expanduser().resolve()): verdict
               for verdict in verdicts}
    accepted, waiting = [], []
    for sample in targets:
        check_cancelled(cancelled)
        verdict = by_path.get(str(Path(sample["input_path"]).expanduser().resolve()))
        if verdict is None:
            skipped.append((sample["id"], "Identification did not complete for this sample."))
            continue
        decided = apply_policy(verdict, effective)
        bucket, genus, species = proposed_destination(decided, policy=effective)
        if bucket is None and genus:
            confirm_organism(project, [sample["id"]], genus, species, evidence=decided,
                             typing_mode="manual", basis=decided.get("basis") or "genomic_ani",
                             confidence=decided.get("confidence") or "unresolved",
                             confirmed_by=decided.get("confirmed_by") or "auto_policy")
            accepted.append((sample["id"], genus, species))
        else:
            project.update_metadata(sample["id"], {
                "organism_evidence": organism_evidence(decided)})
            project.record_history(sample["id"], "organism_identified", {
                "basis": decided.get("basis"), "confidence": decided.get("confidence"),
                "status": decided.get("status"), "proposed": decided.get("proposed"),
                "bucket": bucket})
            waiting.append((sample["id"], bucket))
    moved = refile_samples(project, [sample["id"] for sample in targets],
                           storage_root=storage_root, cancelled=cancelled)
    identification = _identification_report(
        True, species_panel_root, kpsc_panel_root, verdicts,
        [sample["id"] for sample in targets], [sample_id for sample_id, _ in waiting])
    return {"accepted": accepted, "waiting": waiting, "skipped": [*skipped, *moved["skipped"]],
            "moved": moved["moved"], "verdicts": verdicts, "policy": effective,
            "identification": identification}


def missing_managed_copies(project: Project) -> list[dict]:
    """Managed copies this project expects and cannot find. Stats files, hashes none."""
    missing = []
    for sample in project.samples():
        workflow = (sample.get("metadata") or {}).get("workflow", {})
        if not workflow.get("managed") or not sample.get("input_path"):
            continue
        if Path(sample["input_path"]).is_file():
            continue
        missing.append({"sample_id": sample["id"], "name": sample["name"],
                        "path": sample["input_path"],
                        "storage_root": str(workflow.get("storage_root") or ""),
                        "managed_sha256": str(workflow.get("managed_sha256") or "")})
    return missing


def _relocation_candidates(project: Project, storage_root, recorded: str) -> list[Path]:
    """Folders a managed tree could have moved to, most explicit first."""
    roots, seen = [], set()
    options = [storage_root, default_storage_root(project)]
    if recorded:
        options.append(project.path.parent / Path(recorded).name)
    for option in options:
        if not option:
            continue
        path = Path(option).expanduser().resolve()
        if str(path) not in seen:
            seen.add(str(path))
            roots.append(path)
    return roots


def relocate_managed_storage(project: Project, *, storage_root=None, cancelled=None,
                             progress=None) -> dict:
    """Re-point managed copies that travelled with the project, verified byte for byte.

    A portable project is reopened from wherever its folder now is — another
    computer, another drive letter — and the absolute paths recorded at import no
    longer exist. A candidate in the same place under the new root is adopted only
    when its SHA-256 equals the fingerprint recorded for that sample's managed copy,
    so a file that merely has the right name is never adopted as this isolate's
    evidence. Nothing is copied, moved or deleted; only the recorded location changes.
    """
    outstanding = missing_managed_copies(project)
    report: dict[str, list] = {"relocated": [], "unresolved": []}
    for index, entry in enumerate(outstanding):
        check_cancelled(cancelled)
        if progress:
            progress(index, len(outstanding), f"Looking for {Path(entry['path']).name}")
        sample_id, expected = entry["sample_id"], entry["managed_sha256"].lower()
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            report["unresolved"].append((sample_id, "No managed-copy SHA-256 is recorded for this "
                                                    "sample, so no replacement can be verified."))
            continue
        recorded_root = entry["storage_root"]
        if not recorded_root:
            report["unresolved"].append((sample_id, "No managed storage location is recorded for "
                                                    "this sample."))
            continue
        try:
            relative = Path(entry["path"]).relative_to(Path(recorded_root))
        except ValueError:
            report["unresolved"].append((sample_id, "The recorded file is not inside the recorded "
                                                    "managed storage location."))
            continue
        adopted = None
        for root in _relocation_candidates(project, storage_root, recorded_root):
            check_cancelled(cancelled)
            candidate = root / relative
            if not candidate.is_file() or candidate.is_symlink():
                continue
            if file_sha256(candidate, cancelled) != expected:
                report["unresolved"].append((sample_id, f"A file exists at {candidate} but its "
                                                        "contents differ from the managed copy "
                                                        "recorded for this sample; it was not adopted."))
                adopted = False
                break
            with project.transaction():
                current = project.get_sample(sample_id)
                if current["status"] == "running" or current["input_path"] != entry["path"]:
                    raise ValueError("Sample input changed or analysis is running; retry later.")
                project.set_input_path(sample_id, candidate)
                project.update_metadata(sample_id, {"workflow": {
                    "managed": True, "storage_root": str(root), "managed_sha256": expected}})
                if current["metadata"].get("workflow", {}).get("source_kind") == "assembly":
                    project.update_metadata(sample_id, {"assembly": {"assembly_path": str(candidate)}})
                project.record_history(sample_id, "managed_storage_relocated", {
                    "from": entry["path"], "to": str(candidate), "sha256": expected,
                    "verified": "sha256 matches the recorded managed copy"})
            report["relocated"].append((sample_id, str(candidate)))
            adopted = True
            break
        if adopted is None:
            report["unresolved"].append((sample_id, "The managed copy was not found under any "
                                                    "folder this project knows about."))
    if progress:
        progress(len(outstanding), len(outstanding), "Storage check complete")
    return report


def _filing(project: Project, sample: dict, storage_root, append_st, quarantine) -> dict:
    """Shared destination computation, so a preview can never drift from a move."""
    metadata = sample["metadata"]
    workflow = metadata.get("workflow", {})
    source = Path(sample["input_path"]).resolve()
    root = Path(storage_root or workflow.get("storage_root") or
                default_storage_root(project)).expanduser().resolve()
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
    bucket = review_bucket(metadata) if quarantine is None else quarantine_token(quarantine)
    # An assembled isolate must never inherit its original FASTQ extension.
    naming_source = (Path(safe_component(sample["name"]) + ".fasta")
                     if workflow.get("source_kind") == "assembly" else original)
    destination = _target(root, sample["id"], naming_source, organism, st, append, quarantine=bucket)
    return {"root": root, "source": source, "original": original, "organism": organism,
            "st": st, "append_st": append, "quarantine": bucket, "destination": destination,
            "workflow": workflow, "result": result}


def plan_filing(project: Project, sample_id: str, *, storage_root=None, append_st: bool | None = None,
                quarantine=None) -> dict:
    """Where a sample would be filed, computed without touching a single file.

    This is what a review dialog shows before anything is copied: the same
    function organise and re-file use, so the preview and the move cannot differ.
    """
    sample = project.get_sample(sample_id)
    evidence = sample["metadata"].get("organism_evidence") or {}
    plan = {"sample_id": sample_id, "name": sample["name"], "eligible": False, "reason": "",
            "current": None, "destination": None, "relative": "", "changed": False,
            "quarantine": None, "genus": "", "species": "", "st": None, "root": None,
            "basis": evidence.get("basis") or "", "confidence": evidence.get("confidence") or "",
            "status": evidence.get("status") or ""}
    if sample.get("profile_only") or not sample.get("input_path"):
        plan["reason"] = "A profile-only sample has no sequence file to file."
        return plan
    if not sample["metadata"].get("workflow", {}).get("managed"):
        # Organising means moving WMLSTudio's own copy. A file the user chose to
        # keep where it is stays there; correcting its label never adopts it.
        plan["reason"] = ("This file is linked where you keep it, so only its label changed. "
                          "WMLSTudio organises its own managed copies.")
        return plan
    filing = _filing(project, sample, storage_root, append_st, quarantine)
    plan.update(eligible=True, current=filing["source"], destination=filing["destination"],
                relative=filing["destination"].relative_to(filing["root"]).as_posix(),
                changed=filing["source"] != filing["destination"], quarantine=filing["quarantine"],
                genus=filing["organism"].get("genus", ""), species=filing["organism"].get("species", ""),
                st=filing["st"], root=filing["root"])
    if sample["status"] == "running":
        plan.update(eligible=False, reason="Wait for this sample's analysis before re-filing it.")
    return plan


def organize_sample(
    project: Project, sample_id: str, storage_root=None, append_st: bool | None = None,
    cancelled=None, *, quarantine=None,
) -> Path:
    """Copy to the assigned genus/species/ST location, then update the project.

    Only an existing managed copy is eligible for cleanup. Original paths remain
    in workflow metadata. If cleanup fails, its path is recorded for recovery.
    An explicit quarantine token forces the needs-review tree; when it is omitted
    the bucket recorded in metadata['organism_evidence'] decides, so a sample the
    app declined to identify is never pushed into an unsupported genus folder.
    """
    sample = project.get_sample(sample_id)
    filing = _filing(project, sample, storage_root, append_st, quarantine)
    workflow, result = filing["workflow"], filing["result"]
    source, original, root = filing["source"], filing["original"], filing["root"]
    append, destination = filing["append_st"], filing["destination"]
    if filing["quarantine"]:
        _quarantine_notice(root)
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
                _prune_empty(Path(previous_root).resolve(), source.parent)
        except OSError:
            project.record_history(sample_id, "managed_copy_retained", {"path": str(source)})
    return destination


def refile_samples(project: Project, sample_ids, *, storage_root=None, append_st: bool | None = None,
                   cancelled=None, progress=None) -> dict:
    """Move managed copies whose planned destination no longer matches where they are.

    A sample that is running, profile-only or unmanaged is skipped with a reason
    rather than failed, so correcting twenty labels never aborts on the first one.
    """
    identifiers = list(dict.fromkeys(sample_ids))
    report: dict[str, list] = {"moved": [], "unchanged": [], "skipped": []}
    for index, sample_id in enumerate(identifiers):
        check_cancelled(cancelled)
        if progress:
            progress(index, len(identifiers), f"Re-filing {index + 1} of {len(identifiers)}")
        plan = plan_filing(project, sample_id, storage_root=storage_root, append_st=append_st)
        if not plan["eligible"]:
            report["skipped"].append((sample_id, plan["reason"]))
            continue
        if not plan["changed"]:
            report["unchanged"].append(sample_id)
            continue
        try:
            destination = organize_sample(project, sample_id, storage_root, append_st, cancelled)
        except AnalysisCancelled:
            raise
        except (OSError, ValueError) as error:
            report["skipped"].append((sample_id, str(error)))
            continue
        report["moved"].append((sample_id, str(destination)))
    if progress:
        progress(len(identifiers), len(identifiers), "Re-filing complete")
    return report


def _evidence_patch(evidence, **changes) -> dict:
    record = organism_evidence(dict(evidence or {})) or {}
    record.setdefault("format_version", ORGANISM_EVIDENCE_VERSION)
    record.update(changes)
    return organism_evidence(record)


def confirm_organism(project: Project, sample_ids, genus: str, species: str = "", *,
                     evidence=None, scheme_path=None, typing_mode: str = "manual",
                     basis: str = "user_assigned", confidence: str = "unresolved",
                     confirmed_by: str = "user") -> None:
    """Accept an organism for filing, recording who accepted it and on what basis.

    Files are not moved here; the caller runs refile_samples in a worker so a slow
    copy never blocks the interface. A spreadsheet or a click is an operator
    statement, so its confidence stays 'unresolved' and never becomes genomic evidence.
    """
    assign_organism(project, sample_ids, genus, species, scheme_path, typing_mode)
    accepted = {"genus": str(genus).strip(), "species": str(species).strip()}
    stamp = datetime.now(UTC).isoformat()
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            current = project.get_sample(sample_id)["metadata"].get("organism_evidence")
            record = _evidence_patch(evidence if evidence is not None else current,
                                     status="confirmed", accepted=accepted,
                                     quarantine_reason=None, confirmed_by=confirmed_by,
                                     confirmed_utc=stamp)
            proposed = record.get("proposed") or {}
            record.setdefault("basis", basis)
            record.setdefault("confidence", confidence)
            # A label the operator typed over the proposal is not supported by the
            # evidence that produced the proposal, so it loses that evidence's basis.
            if ((proposed.get("genus", ""), proposed.get("species", "")) != (accepted["genus"], accepted["species"])
                    or record["basis"] in {"", "none", None}):
                record.update(basis=basis, confidence=confidence)
            project.update_metadata(sample_id, {"organism_evidence": record})
            project.record_history(sample_id, "organism_confirmed", {
                "accepted": accepted, "basis": record.get("basis"),
                "confidence": record.get("confidence"), "confirmed_by": confirmed_by})


def quarantine_samples(project: Project, sample_ids, reason: str = "user_deferred", *,
                       evidence=None) -> None:
    """Send samples to the needs-review tree without claiming anything about them."""
    token = quarantine_token(reason)
    if token is None:
        raise ValueError("Quarantine requires an explicit reason.")
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            sample = project.get_sample(sample_id)
            if sample["status"] == "running":
                raise ValueError("Wait for this sample's analysis before changing its assignment.")
            current = sample["metadata"].get("organism_evidence")
            record = _evidence_patch(evidence if evidence is not None else current,
                                     status="quarantined", quarantine_reason=token,
                                     confirmed_by=None, confirmed_utc=None)
            record.setdefault("basis", "none")
            record.setdefault("confidence", "unresolved")
            project.update_metadata(sample_id, {"organism_evidence": record})
            project.record_history(sample_id, "organism_quarantined", {
                "reason": token, "bucket": QUARANTINE_BUCKETS[token]})


def reassign_organism(project: Project, sample_ids, genus: str, species: str = "", *,
                      scheme_path=None, typing_mode: str = "manual", storage_root=None,
                      basis: str = "user_assigned", cancelled=None, progress=None) -> dict:
    """Confirm a corrected organism and then move the managed copies to match it.

    This is the whole of a manual override: without the move, a corrected label
    would leave the file sitting in a folder that contradicts it.
    """
    confirm_organism(project, sample_ids, genus, species, scheme_path=scheme_path,
                     typing_mode=typing_mode, basis=basis)
    return refile_samples(project, sample_ids, storage_root=storage_root,
                          cancelled=cancelled, progress=progress)


def _managed_copy(sample: dict, *, cancelled=None) -> tuple[Path | None, str]:
    """Decide whether this app created the file at input_path and may delete it."""
    workflow = (sample.get("metadata") or {}).get("workflow", {})
    if not sample.get("input_path"):
        return None, "This sample has no sequence file."
    if not workflow.get("managed"):
        return None, "This file was linked, not copied here; WMLSTudio will not delete it."
    storage_root = workflow.get("storage_root")
    if not storage_root:
        return None, "No managed storage location is recorded for this sample."
    path = Path(sample["input_path"]).resolve()
    root = Path(storage_root).expanduser().resolve()
    if not path.is_relative_to(root) or sample["id"] not in path.relative_to(root).parts:
        return None, "This file is outside the managed storage location for this sample."
    if path == Path(workflow.get("source_path") or path).expanduser().resolve():
        return None, "This is your original file, not a managed copy."
    if not path.is_file() or path.is_symlink():
        return None, "The managed copy is missing or is not a regular file."
    if file_sha256(path, cancelled) != workflow.get("managed_sha256"):
        return None, "The file no longer matches the recorded managed copy; it was left in place."
    return path, ""


def managed_copy_for(project: Project, sample_id: str, *, cancelled=None) -> Path | None:
    """The verified managed copy this app may delete, or None when there is none.

    None never means "nothing is there": it means nothing WMLSTudio is entitled to
    remove is there. A user's original file is never returned from here.
    """
    return _managed_copy(project.get_sample(sample_id), cancelled=cancelled)[0]


def delete_managed_copy(project: Project, sample_id: str, *, cancelled=None) -> dict:
    """Reclaim the space a managed copy uses, relinking the sample to the original.

    Refused unless the user's original still exists with its recorded hash, because
    deleting the only remaining copy of the evidence is never a cleanup operation.
    """
    sample = project.get_sample(sample_id)
    if sample["status"] == "running":
        raise ValueError("Wait for this sample's analysis before removing its managed copy.")
    path, reason = _managed_copy(sample, cancelled=cancelled)
    if path is None:
        return {"deleted": False, "path": None, "reason": reason}
    workflow = sample["metadata"].get("workflow", {})
    original = Path(workflow.get("source_path") or "").expanduser()
    if not original.is_file() or file_sha256(original, cancelled) != workflow.get("source_sha256"):
        project.record_history(sample_id, "managed_copy_retained", {
            "path": str(path), "reason": "original_unavailable"})
        return {"deleted": False, "path": str(path),
                "reason": "Your original file is missing or has changed, so this copy is the only "
                          "remaining evidence and was kept."}
    root = Path(workflow["storage_root"]).expanduser().resolve()
    with project.transaction():
        current = project.get_sample(sample_id)
        if current["status"] == "running" or current["input_path"] != str(path):
            raise ValueError("Sample input changed or analysis is running; retry later.")
        project.set_input_path(sample_id, original.resolve())
        project.update_metadata(sample_id, {"workflow": {
            "managed": False, "managed_sha256": None, "storage_root": None}})
    path.unlink(missing_ok=True)
    _prune_empty(root, path.parent)
    project.record_history(sample_id, "managed_copy_removed", {
        "path": str(path), "relinked_to": str(original.resolve())})
    return {"deleted": True, "path": str(path), "reason": ""}


def remove_samples(project: Project, sample_ids, *, delete_managed_copy: bool = False,
                   cancelled=None) -> dict:
    """Remove samples, optionally deleting only copies this app made and verified.

    The user's original file is never touched, whatever the caller asks for.
    """
    report: dict[str, list] = {"removed": [], "deleted": [], "retained": []}
    for sample_id in dict.fromkeys(sample_ids):
        check_cancelled(cancelled)
        sample = project.get_sample(sample_id)
        if sample["status"] == "running":
            raise ValueError("Wait for this sample's analysis before removing it.")
        path, reason = (_managed_copy(sample, cancelled=cancelled) if delete_managed_copy
                        else (None, "The managed copy was kept."))
        root = Path(sample["metadata"].get("workflow", {}).get("storage_root") or "").expanduser()
        with project.transaction():
            project.remove_sample(sample_id)
        report["removed"].append(sample_id)
        if path is None:
            if delete_managed_copy:
                report["retained"].append((sample_id, reason))
                project.record_history(sample_id, "managed_copy_retained", {"reason": reason})
            continue
        try:
            path.unlink()
        except OSError as error:
            report["retained"].append((sample_id, str(error)))
            project.record_history(sample_id, "managed_copy_retained",
                                   {"path": str(path), "reason": str(error)})
            continue
        if root.is_dir():
            _prune_empty(root.resolve(), path.parent)
        report["deleted"].append(str(path))
        project.record_history(sample_id, "managed_copy_removed", {"path": str(path)})
    return report
