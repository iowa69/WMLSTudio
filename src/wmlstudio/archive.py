"""Reversible archiving: hide an isolate from the working views, lose nothing.

Archiving never deletes a sample row, a stored result, an allele call, a history
entry, or a file on disk. It writes a dated, optionally justified note into the
sample's metadata; restoring removes that note and leaves every other record
exactly as it was. ``Project.samples()``, the audit trail and full exports keep
returning archived isolates, so the evidence chain stays complete. Only the
working views drop them, through ``active_samples``.

Discarding evidence is a separate, explicit operation: ``Project.remove_sample``,
which itself keeps the user's original input file and records what it deleted.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from typing import Any

from wmlstudio.project import Project

FORMAT_VERSION = 1


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def archive_record(sample: Mapping[str, Any]) -> dict[str, Any] | None:
    """The sample's archive note (``at``, ``reason``), or None when it is active.

    Only an object written by this module counts. A stray ``archived`` annotation
    of some other shape leaves the isolate visible, because the safe failure is a
    sample the user can still see.
    """
    if not isinstance(sample, Mapping):
        return None
    record = (sample.get("metadata") or {}).get("archived")
    return dict(record) if isinstance(record, Mapping) else None


def is_archived(sample: Mapping[str, Any]) -> bool:
    """True when the sample is marked out of the working set."""
    return archive_record(sample) is not None


def active_samples(samples: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """The isolates the user is currently working with; the stored project is unchanged."""
    return [sample for sample in samples if not is_archived(sample)]


def archived_samples(samples: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """The archived isolates, still fully present in the project with their evidence."""
    return [sample for sample in samples if is_archived(sample)]


def archive_samples(project: Project, sample_ids: Iterable[str], reason: str = "") -> list[str]:
    """Mark isolates as out of the working set, keeping every record they carry.

    Results, allele calls, history and files are untouched; the isolate stays in
    ``Project.samples()`` and in any export. A running analysis is refused, as an
    organism assignment is. Returns the identifiers actually changed.
    """
    text = str(reason or "").strip()
    changed: list[str] = []
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            sample = project.get_sample(sample_id)
            if sample["status"] == "running":
                raise ValueError("Wait for this sample's analysis before archiving it.")
            if is_archived(sample):
                continue
            record = {"at": _now(), "reason": text, "format_version": FORMAT_VERSION}
            project.update_metadata(sample_id, {"archived": record})
            project.record_history(sample_id, "sample_archived", record)
            changed.append(sample_id)
    return changed


def restore_samples(project: Project, sample_ids: Iterable[str]) -> list[str]:
    """Return archived isolates to the working views with their evidence intact.

    The archive note is removed rather than set to a falsey value, so a restored
    isolate is indistinguishable from one that was never archived. Returns the
    identifiers actually changed.
    """
    changed: list[str] = []
    with project.transaction():
        for sample_id in dict.fromkeys(sample_ids):
            sample = project.get_sample(sample_id)
            if sample["status"] == "running":
                raise ValueError("Wait for this sample's analysis before restoring it.")
            if not is_archived(sample):
                continue
            metadata = dict(sample["metadata"])
            record = metadata.pop("archived", None)
            project.set_metadata(sample_id, metadata)  # update_metadata merges, it cannot delete
            project.record_history(sample_id, "sample_restored", {"archived": record})
            changed.append(sample_id)
    return changed
