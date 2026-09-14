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

SCHEMA_VERSION = 3
APPLICATION_ID = 0x574D4C53  # WMLS
STATUSES = frozenset({"queued", "running", "failed", "completed", "interrupted"})


class ProjectError(ValueError):
    """The selected file is not a supported WMLSTudio project."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True)


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
        if version not in {1, 2, SCHEMA_VERSION}:
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
        if version in {1, 2}:
            with self.transaction():
                if version == 1:
                    self._create_history()
                self._create_analyses()
                for row in self._connection.execute("SELECT id, result FROM samples WHERE result IS NOT NULL"):
                    result = json.loads(row["result"])
                    if result.get("scheme_digest"):
                        self.set_analysis(row["id"], result)
                self._connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
                self.record_history(None, "schema_migrated", {"from": version, "to": SCHEMA_VERSION})
        else:
            columns = {row[1] for row in self._connection.execute("PRAGMA table_info(history)")}
            if not {"id", "sample_id", "action", "details", "created_at"} <= columns:
                raise ProjectError("Project history table is missing required columns.")
            columns = {row[1] for row in self._connection.execute("PRAGMA table_info(analyses)")}
            if not {"sample_id", "scheme_digest", "result", "updated_at"} <= columns:
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
            "result TEXT NOT NULL, updated_at TEXT NOT NULL, PRIMARY KEY(sample_id, scheme_digest))"
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

    def set_result(self, sample_id: str, result: Mapping[str, Any]) -> None:
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
                self.set_analysis(sample_id, snapshot)
            self.record_history(sample_id, "analysis_completed", {
                "scheme": snapshot.get("scheme"), "scheme_digest": snapshot.get("scheme_digest"),
                "input_sha256": snapshot.get("input_sha256"), "st": snapshot.get("st"),
            })

    def set_analysis(self, sample_id: str, result: Mapping[str, Any]) -> None:
        """Store a secondary scheme result without replacing primary MLST/ST or job state."""
        snapshot = dict(result)
        digest = snapshot.get("scheme_digest")
        if not isinstance(digest, str) or not digest:
            raise ValueError("A stored scheme analysis requires a nonempty scheme fingerprint.")
        snapshot["sample_id"] = sample_id
        encoded = _json(snapshot)
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
                "INSERT INTO analyses (sample_id, scheme_digest, result, updated_at) VALUES (?, ?, ?, ?) "
                "ON CONFLICT(sample_id, scheme_digest) DO UPDATE SET "
                "result = excluded.result, updated_at = excluded.updated_at",
                (sample_id, digest, encoded, _now()),
            )

    def analysis_results(self, sample_id: str) -> list[dict[str, Any]]:
        with self._lock:
            self._check_open()
            self.get_sample(sample_id)
            return [json.loads(row["result"]) for row in self._connection.execute(
                "SELECT result FROM analyses WHERE sample_id = ? ORDER BY updated_at, scheme_digest",
                (sample_id,),
            )]

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
                    'SELECT a.scheme_digest, a.updated_at, ' + expressions + ', '
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
                    'SELECT scheme_digest, updated_at, result FROM analyses WHERE sample_id = ? ORDER BY updated_at, scheme_digest',
                    (sample_id,),
                ):
                    result = json.loads(row['result'])
                    summaries.append({**{field: result.get(field) for field in fields},
                                      'scheme_digest': row['scheme_digest'], 'updated_at': row['updated_at'],
                                      'locus_count': len(result['alleles']) if isinstance(result.get('alleles'), dict) else 0})
                return summaries

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
        """Remove the stored record and result; never delete the input file."""
        with self.transaction():
            sample = self.get_sample(sample_id)
            cursor = self._connection.execute("DELETE FROM samples WHERE id = ?", (sample_id,))
            if cursor.rowcount != 1:
                raise KeyError(f"Unknown sample: {sample_id}")
            self.record_history(sample_id, "sample_removed", {"sample": sample})
            self._connection.execute("DELETE FROM analyses WHERE sample_id = ?", (sample_id,))

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
