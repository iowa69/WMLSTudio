"""Transactional, single-file project storage for the desktop application.

Sequence inputs are referenced by absolute path, not copied into the project.
Result snapshots remain available if an input is moved or disconnected.
"""

from __future__ import annotations

import json
import sqlite3
import threading
import uuid
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 4
APPLICATION_ID = 0x574D4C53  # WMLS
STATUSES = frozenset({"queued", "running", "failed", "completed", "interrupted"})

# Classical MLST and core-genome MLST answer different questions and are never
# interchangeable, so every stored profile is filed under one of these and a
# profile that carries no scheme at all is left unclassified rather than sorted
# into a tree it does not belong in.
TYPING_KINDS = ("mlst", "cgmlst", "unclassified")
# The same locus count the analysis planner already uses to choose between the
# exact-match caller and the full-CDS cgMLST caller, so a stored profile is
# classified the way it was produced.
CGMLST_LOCUS_FLOOR = 30


class ProjectError(ValueError):
    """The selected file is not a supported WMLSTudio project."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True)


def classify_typing(result: Mapping[str, Any]) -> dict[str, str]:
    """Say whether a stored profile is classical MLST or cgMLST, and on what basis.

    The basis travels with the answer because the evidence differs in strength:
    a caller or a scheme that named itself is a record, while a locus count is an
    inference from the size of the panel. A result with no allelic profile stays
    'unclassified'; it is never pushed into one of the two typing views on a guess.
    A whole-genome (wgMLST) scheme is grouped with cgMLST, and the basis says so.
    """
    if not isinstance(result, Mapping):
        raise TypeError("A stored analysis must be a mapping.")
    declared = str(result.get("analysis_kind") or "").strip().casefold()
    if declared in {"mlst", "cgmlst"}:
        return {"kind": declared, "basis": "recorded in the stored result"}
    method = str((result.get("parameters") or {}).get("method") or "").casefold()
    if "cgmlst" in method:
        return {"kind": "cgmlst", "basis": f"allele caller {method!r} recorded in the result"}
    metadata = result.get("scheme_metadata")
    scheme_type = (str(metadata.get("type") or "").strip().casefold()
                   if isinstance(metadata, Mapping) else "")
    if scheme_type in {"mlst", "cgmlst", "wgmlst"}:
        return {"kind": "mlst" if scheme_type == "mlst" else "cgmlst",
                "basis": f"scheme metadata declares type {scheme_type!r}"}
    alleles = result.get("alleles")
    count = len(alleles) if isinstance(alleles, Mapping) else 0
    if not count or not result.get("scheme_digest"):
        return {"kind": "unclassified", "basis": "no allelic profile is stored"}
    if count > CGMLST_LOCUS_FLOOR:
        return {"kind": "cgmlst",
                "basis": f"{count} loci, above the {CGMLST_LOCUS_FLOOR}-locus core-genome floor"}
    return {"kind": "mlst", "basis": f"{count} loci, at classical scheme size"}


def typing_kind(result: Mapping[str, Any]) -> str:
    """'mlst', 'cgmlst' or 'unclassified' for one stored analysis result."""
    return classify_typing(result)["kind"]


class Project:
    """A portable SQLite project with stable sample identifiers.

    Opening a project converts persisted running jobs to interrupted jobs. They
    must be explicitly rerun; opening a file never starts an analysis.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path).expanduser().resolve()
        self._lock = threading.RLock()
        self._closed = False
        self._depth = 0
        self._connection = sqlite3.connect(
            str(self.path), timeout=10, isolation_level=None, check_same_thread=False
        )
        self._connection.row_factory = sqlite3.Row
        try:
            self._initialize()
            with self.transaction():
                self._connection.execute(
                    "UPDATE samples SET status = 'interrupted', error = ?, updated_at = ? "
                    "WHERE status = 'running'",
                    ("Analysis was interrupted before completion. Rerun this sample.", _now()),
                )
        except BaseException:
            self._connection.close()
            self._closed = True
            raise

    def _initialize(self) -> None:
        try:
            tables = {
                row[0]
                for row in self._connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            version = self._connection.execute("PRAGMA user_version").fetchone()[0]
            application = self._connection.execute("PRAGMA application_id").fetchone()[0]
        except sqlite3.DatabaseError as exc:
            raise ProjectError("This file is not a readable SQLite project.") from exc
        if not tables and version == 0 and application == 0:
            with self.transaction():
                self._connection.execute(
                    "CREATE TABLE samples ("
                    "id TEXT PRIMARY KEY, name TEXT NOT NULL, input_path TEXT NOT NULL, "
                    "status TEXT NOT NULL CHECK(status IN "
                    "('queued','running','failed','completed','interrupted')), "
                    "error TEXT NOT NULL DEFAULT '', metadata TEXT NOT NULL DEFAULT '{}', "
                    "result TEXT, created_at TEXT NOT NULL, updated_at TEXT NOT NULL)"
                )
                self._connection.execute(
                    "CREATE TABLE settings (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
                )
                self._create_history()
                self._create_analyses()
                self._connection.execute(f"PRAGMA application_id = {APPLICATION_ID}")
                self._connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            return
        if application != APPLICATION_ID:
            raise ProjectError("This SQLite file is not a WMLSTudio project.")
        if version not in {1, 2, 3, SCHEMA_VERSION}:
            raise ProjectError(
                f"Project schema version {version} is unsupported; "
                f"this application supports version {SCHEMA_VERSION}."
            )
        expected = {
            "samples": {
                "id", "name", "input_path", "status", "error", "metadata", "result",
                "created_at", "updated_at",
            },
            "settings": {"key", "value"},
        }
        for table, columns in expected.items():
            actual = {
                row[1] for row in self._connection.execute(f"PRAGMA table_info({table})")
            }
            if not columns <= actual:
                raise ProjectError(f"Project table {table!r} is missing required columns.")
        if version in {1, 2, 3}:
            with self.transaction():
                if version == 1:
                    self._create_history()
                if version in {1, 2}:
                    self._create_analyses()
                    for row in self._connection.execute("SELECT id, result FROM samples WHERE result IS NOT NULL"):
                        result = json.loads(row["result"])
                        if result.get("scheme_digest"):
                            self.set_analysis(row["id"], result)
                self._classify_stored_analyses()
                self._connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
                self.record_history(None, "schema_migrated", {"from": version, "to": SCHEMA_VERSION})
        else:
            columns = {row[1] for row in self._connection.execute("PRAGMA table_info(history)")}
            if not {"id", "sample_id", "action", "details", "created_at"} <= columns:
                raise ProjectError("Project history table is missing required columns.")
            columns = {row[1] for row in self._connection.execute("PRAGMA table_info(analyses)")}
            if not {"sample_id", "scheme_digest", "result", "typing_kind", "updated_at"} <= columns:
                raise ProjectError("Project analysis table is missing required columns.")

    def _create_history(self) -> None:
        self._connection.execute(
            "CREATE TABLE history (id INTEGER PRIMARY KEY AUTOINCREMENT, sample_id TEXT, "
            "action TEXT NOT NULL, details TEXT NOT NULL, created_at TEXT NOT NULL)"
        )
        self._connection.execute("CREATE INDEX history_sample ON history(sample_id, id)")

    def _create_analyses(self) -> None:
        self._connection.execute(
            "CREATE TABLE analyses (sample_id TEXT NOT NULL, scheme_digest TEXT NOT NULL, "
            "result TEXT NOT NULL, typing_kind TEXT NOT NULL DEFAULT '', updated_at TEXT NOT NULL, "
            "PRIMARY KEY(sample_id, scheme_digest))"
        )

    def _classify_stored_analyses(self) -> None:
        """Label every stored profile MLST or cgMLST once, from its own payload.

        Nothing is rewritten: the scientific result is untouched and only the
        column that says which of the two typing views it belongs to is filled in.
        """
        columns = {row[1] for row in self._connection.execute("PRAGMA table_info(analyses)")}
        if "typing_kind" not in columns:
            self._connection.execute(
                "ALTER TABLE analyses ADD COLUMN typing_kind TEXT NOT NULL DEFAULT ''"
            )
        rows = self._connection.execute(
            "SELECT sample_id, scheme_digest, result FROM analyses WHERE typing_kind = ''"
        ).fetchall()
        for row in rows:
            self._connection.execute(
                "UPDATE analyses SET typing_kind = ? WHERE sample_id = ? AND scheme_digest = ?",
                (typing_kind(json.loads(row["result"])), row["sample_id"], row["scheme_digest"]),
            )

    def _check_open(self) -> None:
        if self._closed:
            raise RuntimeError("This project is closed.")

    @contextmanager
    def transaction(self) -> Iterator[Project]:
        """Group writes atomically, with savepoints for nested operations."""
        with self._lock:
            self._check_open()
            depth = self._depth
            savepoint = f"wmlstudio_{depth}"
            self._connection.execute("BEGIN IMMEDIATE" if depth == 0 else f"SAVEPOINT {savepoint}")
            self._depth += 1
            try:
                yield self
                self._connection.execute("COMMIT" if depth == 0 else f"RELEASE {savepoint}")
            except BaseException:
                if depth == 0:
                    self._connection.execute("ROLLBACK")
                else:
                    self._connection.execute(f"ROLLBACK TO {savepoint}")
                    self._connection.execute(f"RELEASE {savepoint}")
                raise
            finally:
                self._depth -= 1

    def add_sample(
        self, path: str | Path, name: str | None = None, *, sample_id: str | None = None,
    ) -> str:
        source = Path(path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"Sequence input does not exist: {source}")
        sample_name = source.name if name is None else name.strip()
        if not sample_name:
            raise ValueError("Sample name cannot be empty.")
        sample_id = sample_id or uuid.uuid4().hex
        if not isinstance(sample_id, str) or not sample_id.strip():
            raise ValueError("Sample identifier cannot be empty.")
        timestamp = _now()
        with self.transaction():
            self._connection.execute(
                "INSERT INTO samples "
                "(id, name, input_path, status, created_at, updated_at) VALUES (?, ?, ?, ?, ?, ?)",
                (sample_id, sample_name, str(source), "queued", timestamp, timestamp),
            )
            self.record_history(sample_id, "sample_imported", {"input_path": str(source)})
        return sample_id

    @staticmethod
    def _decode(row: sqlite3.Row) -> dict[str, Any]:
        sample = dict(row)
        sample["metadata"] = json.loads(sample["metadata"])
        sample["result"] = json.loads(sample["result"]) if sample["result"] is not None else None
        sample["missing_input"] = bool(sample["input_path"]) and not Path(sample["input_path"]).is_file()
        sample["profile_only"] = not bool(sample["input_path"])
        return sample

    def add_profile(
        self, name: str, result: Mapping[str, Any], metadata: Mapping[str, Any] | None = None,
        *, sample_id: str | None = None,
    ) -> str:
        """Import a profile without requiring or fabricating a sequence input file."""
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Profile sample name cannot be empty.")
        identifier = sample_id or uuid.uuid4().hex
        details = dict(metadata or {})
        details["workflow"] = {**details.get("workflow", {}), "source_kind": "profile", "managed": False}
        timestamp = _now()
        with self.transaction():
            self._connection.execute(
                "INSERT INTO samples (id,name,input_path,status,metadata,created_at,updated_at) "
                "VALUES (?,?,'','queued',?,?,?)",
                (identifier, name.strip(), _json(details), timestamp, timestamp),
            )
            self.set_result(identifier, result)
            self.record_history(identifier, "profile_imported", {"source": details.get("provenance", {})})
        return identifier

    def create_collection(self, name: str) -> str:
        if not str(name).strip():
            raise ValueError("Collection name cannot be empty.")
        with self.transaction():
            entries = self.get_setting("collections", [])
            for entry in entries:
                if entry["name"].casefold() == name.strip().casefold():
                    return entry["id"]
            identifier = uuid.uuid4().hex
            entries.append({"id": identifier, "name": name.strip(), "sample_ids": []})
            self.set_setting("collections", entries)
            return identifier

    def collections(self) -> list[dict[str, Any]]:
        entries = self.get_setting("collections", [])
        existing = {sample["id"] for sample in self.samples()}
        return [{**entry, "sample_ids": [sid for sid in entry["sample_ids"] if sid in existing]}
                for entry in entries]

    def set_collection_members(self, collection_id: str, sample_ids, *, add: bool = True) -> None:
        identifiers = list(dict.fromkeys(sample_ids))
        with self.transaction():
            entries = self.collections()
            collection = next((entry for entry in entries if entry["id"] == collection_id), None)
            if collection is None:
                raise KeyError(f"Unknown collection: {collection_id}")
            for sample_id in identifiers:
                self.get_sample(sample_id)
            collection["sample_ids"] = (list(dict.fromkeys([*collection["sample_ids"], *identifiers]))
                                        if add else [sid for sid in collection["sample_ids"] if sid not in identifiers])
            self.set_setting("collections", entries)

    def rename_collection(self, collection_id: str, name: str) -> None:
        if not str(name).strip():
            raise ValueError("Collection name cannot be empty.")
        with self.transaction():
            entries = self.collections()
            found = False
            for entry in entries:
                if entry["id"] == collection_id:
                    entry["name"] = name.strip()
                    found = True
                elif entry["name"].casefold() == name.strip().casefold():
                    raise ValueError("A collection with that name already exists.")
            if not found:
                raise KeyError(f"Unknown collection: {collection_id}")
            self.set_setting("collections", entries)

    def delete_collection(self, collection_id: str) -> None:
        with self.transaction():
            entries = self.collections()
            filtered = [entry for entry in entries if entry["id"] != collection_id]
            if len(filtered) == len(entries):
                raise KeyError(f"Unknown collection: {collection_id}")
            self.set_setting("collections", filtered)

    def samples(self) -> list[dict[str, Any]]:
        with self._lock:
            self._check_open()
            rows = self._connection.execute(
                "SELECT * FROM samples ORDER BY created_at, id"
            ).fetchall()
            return [self._decode(row) for row in rows]

    def comparison_revision(self) -> tuple:
        """Small invalidation key; never decode whole allelic profiles on the GUI thread."""
        with self._lock:
            self._check_open()
            rows = self._connection.execute(
                'SELECT s.id, s.status, s.updated_at, MAX(a.updated_at) '
                'FROM samples s LEFT JOIN analyses a ON a.sample_id = s.id '
                'GROUP BY s.id, s.status, s.updated_at ORDER BY s.id'
            ).fetchall()
            return tuple(tuple(row) for row in rows)

    def get_sample(self, sample_id: str) -> dict[str, Any]:
        with self._lock:
            self._check_open()
            row = self._connection.execute(
                "SELECT * FROM samples WHERE id = ?", (sample_id,)
            ).fetchone()
            if row is None:
                raise KeyError(f"Unknown sample: {sample_id}")
            return self._decode(row)

    def _update(self, sample_id: str, assignments: str, values: tuple[Any, ...]) -> None:
        with self.transaction():
            cursor = self._connection.execute(
                f"UPDATE samples SET {assignments}, updated_at = ? WHERE id = ?",
                (*values, _now(), sample_id),
            )
            if cursor.rowcount != 1:
                raise KeyError(f"Unknown sample: {sample_id}")

    def set_result(self, sample_id: str, result: Mapping[str, Any], *, kind: str | None = None) -> None:
        """Record the sample's current headline result and file it under its typing kind.

        The headline result is whichever analysis ran most recently; the kind-addressed
        copy in the analyses table is what keeps a classical ST and a cgMLST profile
        separately retrievable afterwards.
        """
        if not isinstance(result, Mapping):
            raise TypeError("An analysis result must be a mapping.")
        snapshot = dict(result)
        snapshot["sample_id"] = sample_id
        encoded = _json(snapshot)
        with self.transaction():
            previous = self.get_sample(sample_id)["result"]
            if previous is not None:
                self.record_history(sample_id, "result_superseded", {"result": previous})
            self._update(sample_id, "result = ?, status = 'completed', error = ''", (encoded,))
            if snapshot.get("scheme_digest"):
                self.set_analysis(sample_id, snapshot, kind=kind)
            self.record_history(sample_id, "analysis_completed", {
                "scheme": snapshot.get("scheme"), "scheme_digest": snapshot.get("scheme_digest"),
                "input_sha256": snapshot.get("input_sha256"), "st": snapshot.get("st"),
                "typing_kind": kind or typing_kind(snapshot),
            })

    def set_analysis(self, sample_id: str, result: Mapping[str, Any], *, kind: str | None = None) -> None:
        """Store a scheme result under its own typing kind, replacing nothing else.

        The kind is derived from the result itself unless the caller states it,
        and it is kept beside the profile rather than inside it, so the stored
        scientific payload is exactly what the caller produced.
        """
        snapshot = dict(result)
        digest = snapshot.get("scheme_digest")
        if not isinstance(digest, str) or not digest:
            raise ValueError("A stored scheme analysis requires a nonempty scheme fingerprint.")
        if kind is not None and kind not in TYPING_KINDS:
            raise ValueError(f"Unknown typing kind {kind!r}; expected one of {', '.join(TYPING_KINDS)}.")
        snapshot["sample_id"] = sample_id
        encoded = _json(snapshot)
        stored_kind = kind or typing_kind(snapshot)
        with self.transaction():
            self.get_sample(sample_id)
            previous = self._connection.execute(
                "SELECT result FROM analyses WHERE sample_id = ? AND scheme_digest = ?",
                (sample_id, digest),
            ).fetchone()
            if previous is not None and previous["result"] != encoded:
                self.record_history(sample_id, "scheme_analysis_superseded", {
                    "scheme_digest": digest, "result": json.loads(previous["result"]),
                })
            self._connection.execute(
                "INSERT INTO analyses (sample_id, scheme_digest, result, typing_kind, updated_at) "
                "VALUES (?, ?, ?, ?, ?) ON CONFLICT(sample_id, scheme_digest) DO UPDATE SET "
                "result = excluded.result, typing_kind = excluded.typing_kind, "
                "updated_at = excluded.updated_at",
                (sample_id, digest, encoded, stored_kind, _now()),
            )

    def analysis_results(self, sample_id: str) -> list[dict[str, Any]]:
        with self._lock:
            self._check_open()
            self.get_sample(sample_id)
            return [json.loads(row["result"]) for row in self._connection.execute(
                "SELECT result FROM analyses WHERE sample_id = ? ORDER BY updated_at, rowid",
                (sample_id,),
            )]

    def analyses_by_kind(self, sample_id: str) -> dict[str, list[dict[str, Any]]]:
        """Every stored profile for one sample, split into its MLST and cgMLST sets.

        Both keys are always present, and both are always lists: a sample with an
        ST and no core-genome profile reads as an empty cgMLST list, never as an
        absent one. Full allelic profiles are decoded here; call it off the GUI thread.
        """
        grouped: dict[str, list[dict[str, Any]]] = {kind: [] for kind in TYPING_KINDS}
        with self._lock:
            self._check_open()
            self.get_sample(sample_id)
            for row in self._connection.execute(
                "SELECT result, typing_kind FROM analyses WHERE sample_id = ? "
                "ORDER BY updated_at, rowid", (sample_id,),
            ):
                grouped.setdefault(row["typing_kind"] or "unclassified", []).append(
                    json.loads(row["result"]))
        return grouped

    def latest_analysis(self, sample_id: str, kind: str) -> dict[str, Any] | None:
        """The most recently stored MLST or cgMLST profile for one sample, or None.

        None means this sample has no profile of that kind: it never means the
        other kind's profile can stand in for it.
        """
        if kind not in TYPING_KINDS:
            raise ValueError(f"Unknown typing kind {kind!r}; expected one of {', '.join(TYPING_KINDS)}.")
        with self._lock:
            self._check_open()
            self.get_sample(sample_id)
            row = self._connection.execute(
                "SELECT result FROM analyses WHERE sample_id = ? AND typing_kind = ? "
                "ORDER BY updated_at DESC, scheme_digest DESC LIMIT 1", (sample_id, kind),
            ).fetchone()
            return None if row is None else json.loads(row["result"])

    def results_of_kind(self, kind: str) -> list[dict[str, Any]]:
        """Every stored profile of one typing kind across the project, newest last.

        Built for the separate MLST and cgMLST comparisons: each row is the stored
        result with its sample identity attached, and a sample contributes one row
        per scheme it was typed against. Decodes whole profiles, so run it in a worker.
        """
        if kind not in TYPING_KINDS:
            raise ValueError(f"Unknown typing kind {kind!r}; expected one of {', '.join(TYPING_KINDS)}.")
        with self._lock:
            self._check_open()
            rows = self._connection.execute(
                "SELECT a.result, s.id, s.name, s.status FROM analyses a "
                "JOIN samples s ON s.id = a.sample_id WHERE a.typing_kind = ? "
                "ORDER BY a.updated_at, a.scheme_digest", (kind,),
            ).fetchall()
        return [dict(json.loads(row["result"]), sample_id=row["id"], sample_name=row["name"],
                     job_status=row["status"]) for row in rows]

    def typing_kind_index(self) -> dict[str, dict[str, int]]:
        """How many MLST and cgMLST profiles each sample has, without decoding any.

        Cheap enough for a table column: it reads only the kind column, never a
        locus array, so a sample list can show what a sample already has.
        """
        index: dict[str, dict[str, int]] = {}
        with self._lock:
            self._check_open()
            rows = self._connection.execute(
                "SELECT a.sample_id, a.typing_kind, COUNT(*) AS total FROM analyses a "
                "JOIN samples s ON s.id = a.sample_id GROUP BY a.sample_id, a.typing_kind"
            ).fetchall()
        for row in rows:
            counts = index.setdefault(row["sample_id"], {kind: 0 for kind in TYPING_KINDS})
            counts[row["typing_kind"] or "unclassified"] = row["total"]
        return index

    def analysis_summaries(self, sample_id: str) -> list[dict[str, Any]]:
        """Lightweight labels/fingerprints without decoding locus/call arrays in Python.

        SQLite JSON support is present in packaged CPython. Older external
        SQLite builds use a correctness-preserving fallback, not guessed labels.
        """
        with self._lock:
            self._check_open()
            if self._connection.execute('SELECT 1 FROM samples WHERE id = ?', (sample_id,)).fetchone() is None:
                raise KeyError(f'Unknown sample: {sample_id}')
            fields = ('scheme', 'input_sha256', 'status', 'analysis_kind', 'kind', 'scheme_path', 'st')
            expressions = ', '.join(f"json_extract(a.result, '$.{field}') AS {field}" for field in fields)
            try:
                rows = self._connection.execute(
                    'SELECT a.scheme_digest, a.updated_at, a.typing_kind, ' + expressions + ', '
                    "CASE WHEN json_type(a.result, '$.alleles') = 'object' "
                    "THEN (SELECT COUNT(*) FROM json_each(a.result, '$.alleles')) ELSE 0 END AS locus_count "
                    'FROM analyses a WHERE a.sample_id = ? ORDER BY a.updated_at, a.scheme_digest',
                    (sample_id,),
                ).fetchall()
                return [dict(row) for row in rows]
            except sqlite3.OperationalError as error:
                if not any(message in str(error).lower() for message in ('no such function: json_', 'no such table: json_each')):
                    raise
                summaries = []
                for row in self._connection.execute(
                    'SELECT scheme_digest, updated_at, typing_kind, result FROM analyses WHERE sample_id = ? ORDER BY updated_at, rowid',
                    (sample_id,),
                ):
                    result = json.loads(row['result'])
                    summaries.append({**{field: result.get(field) for field in fields},
                                      'scheme_digest': row['scheme_digest'], 'updated_at': row['updated_at'],
                                      'typing_kind': row['typing_kind'],
                                      'locus_count': len(result['alleles']) if isinstance(result.get('alleles'), dict) else 0})
                return summaries

    def rename_sample(self, sample_id: str, name: str) -> None:
        """Change the display name only; the input file and every stored result are untouched.

        Allele calls, analyses and history stay attached to the sample identifier,
        which never changes.
        """
        text = str(name).strip()
        if not text:
            raise ValueError("Sample name cannot be empty.")
        with self.transaction():
            previous = self.get_sample(sample_id)["name"]
            if previous == text:
                return
            self._update(sample_id, "name = ?", (text,))
            self.record_history(sample_id, "sample_renamed", {"from": previous, "to": text})

    def set_status(self, sample_id: str, status: str, error: str = "") -> None:
        if status not in STATUSES:
            raise ValueError(f"Unknown job status {status!r}.")
        self._update(sample_id, "status = ?, error = ?", (status, str(error)))

    def set_metadata(self, sample_id: str, metadata: Mapping[str, Any]) -> None:
        if not isinstance(metadata, Mapping):
            raise TypeError("Sample metadata must be a mapping.")
        encoded = _json(dict(metadata))
        with self.transaction():
            previous = self.get_sample(sample_id)["metadata"]
            self._update(sample_id, "metadata = ?", (encoded,))
            if previous != metadata:
                self.record_history(sample_id, "metadata_changed", {"before": previous})

    def update_metadata(self, sample_id: str, patch: Mapping[str, Any]) -> None:
        """Merge top-level fields and nested objects without dropping other metadata."""
        def merge(target: dict, changes: Mapping) -> dict:
            for key, value in changes.items():
                if isinstance(value, Mapping) and isinstance(target.get(key), dict):
                    target[key] = merge(target[key], value)
                else:
                    target[key] = value
            return target

        with self.transaction():
            self.set_metadata(sample_id, merge(self.get_sample(sample_id)["metadata"], patch))

    def set_input_path(self, sample_id: str, path: str | Path) -> None:
        source = Path(path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"Sequence input does not exist: {source}")
        with self.transaction():
            previous = self.get_sample(sample_id)["input_path"]
            self._update(sample_id, "input_path = ?", (str(source),))
            if previous != str(source):
                self.record_history(sample_id, "input_relocated", {"from": previous, "to": str(source)})

    def invalidate_result(self, sample_id: str, reason: str) -> None:
        """Archive a result when its assignment changes, then queue a fresh analysis."""
        with self.transaction():
            previous = self.get_sample(sample_id)["result"]
            self.record_history(sample_id, "result_invalidated", {"reason": reason, "result": previous})
            self._update(sample_id, "result = NULL, status = 'queued', error = ''", ())

    def record_history(self, sample_id: str | None, action: str, details: Mapping[str, Any]) -> None:
        with self.transaction():
            self._connection.execute(
                "INSERT INTO history (sample_id, action, details, created_at) VALUES (?, ?, ?, ?)",
                (sample_id, str(action), _json(dict(details)), _now()),
            )

    def history(self, sample_id: str | None = None) -> list[dict[str, Any]]:
        with self._lock:
            self._check_open()
            if sample_id is None:
                rows = self._connection.execute("SELECT * FROM history ORDER BY id").fetchall()
            else:
                rows = self._connection.execute(
                    "SELECT * FROM history WHERE sample_id = ? ORDER BY id", (sample_id,),
                ).fetchall()
            return [{**dict(row), "details": json.loads(row["details"])} for row in rows]

    def remove_sample(self, sample_id: str) -> None:
        """Remove the stored record and result; never delete the input file.

        The row and every stored scheme analysis are captured in history before
        the delete, so restore_removed_sample can rebuild the discarded evidence.
        """
        with self.transaction():
            sample = self.get_sample(sample_id)
            analyses = self.analysis_results(sample_id)
            kinds = {row["scheme_digest"]: row["typing_kind"] for row in self._connection.execute(
                "SELECT scheme_digest, typing_kind FROM analyses WHERE sample_id = ?", (sample_id,))}
            cursor = self._connection.execute("DELETE FROM samples WHERE id = ?", (sample_id,))
            if cursor.rowcount != 1:
                raise KeyError(f"Unknown sample: {sample_id}")
            self.record_history(sample_id, "sample_removed",
                                {"sample": sample, "analyses": analyses, "format_version": 2,
                                 "analysis_kinds": kinds})
            self._connection.execute("DELETE FROM analyses WHERE sample_id = ?", (sample_id,))

    def restore_removed_sample(self, history_id: int) -> str:
        """Rebuild a removed sample from its own removal record; nothing is invented.

        Only a format_version 2 ``sample_removed`` entry carries the analyses it
        deleted. An older entry is refused rather than restored without them, so a
        restored isolate never looks like it had fewer scheme results than it had.
        """
        with self.transaction():
            row = self._connection.execute(
                "SELECT * FROM history WHERE id = ?", (history_id,)
            ).fetchone()
            if row is None:
                raise KeyError(f"Unknown history entry: {history_id}")
            details = json.loads(row["details"])
            sample = details.get("sample")
            if (row["action"] != "sample_removed" or details.get("format_version") != 2
                    or not isinstance(sample, dict)):
                raise ValueError("This history entry does not hold a restorable removed sample.")
            sample_id = str(sample.get("id") or "")
            if not sample_id:
                raise ValueError("The stored removal record carries no sample identifier.")
            if self._connection.execute(
                "SELECT 1 FROM samples WHERE id = ?", (sample_id,)
            ).fetchone() is not None:
                raise ValueError(f"Sample {sample_id} is already in this project; nothing was restored.")
            status, error = sample.get("status"), str(sample.get("error") or "")
            if status not in STATUSES:
                raise ValueError(f"The stored removal record has an unknown job status {status!r}.")
            if status == "running":
                status = "interrupted"
                error = "Analysis was interrupted before completion. Rerun this sample."
            result, timestamp = sample.get("result"), _now()
            self._connection.execute(
                "INSERT INTO samples (id, name, input_path, status, error, metadata, result, "
                "created_at, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (sample_id, str(sample.get("name") or sample_id), str(sample.get("input_path") or ""),
                 status, error, _json(sample.get("metadata") or {}),
                 None if result is None else _json(result),
                 str(sample.get("created_at") or timestamp), timestamp),
            )
            analyses = details.get("analyses") or []
            # A kind a caller stated explicitly is part of the discarded evidence.
            # An older record carries none, so those analyses classify themselves again.
            kinds = details.get("analysis_kinds") or {}
            for analysis in analyses:
                recorded = kinds.get(str(analysis.get("scheme_digest")))
                self.set_analysis(sample_id, analysis,
                                  kind=recorded if recorded in TYPING_KINDS else None)
            self.record_history(sample_id, "sample_restored_from_history",
                                {"history_id": row["id"], "analyses": len(analyses)})
        return sample_id

    def get_setting(self, key: str, default: Any = None) -> Any:
        with self._lock:
            self._check_open()
            row = self._connection.execute(
                "SELECT value FROM settings WHERE key = ?", (key,),
            ).fetchone()
            return default if row is None else json.loads(row["value"])

    def set_setting(self, key: str, value: Any) -> None:
        encoded = _json(value)
        with self.transaction():
            self._connection.execute(
                "INSERT INTO settings (key, value) VALUES (?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value = excluded.value", (key, encoded),
            )

    def close(self) -> None:
        with self._lock:
            if not self._closed:
                if self._depth:
                    raise RuntimeError("Cannot close a project during a transaction.")
                self._connection.close()
                self._closed = True

    def __enter__(self) -> Project:
        self._check_open()
        return self

    def __exit__(self, *_exc: Any) -> None:
        self.close()
