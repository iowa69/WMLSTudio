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

SCHEMA_VERSION = 1
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
                self._connection.execute(f"PRAGMA application_id = {APPLICATION_ID}")
                self._connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            return
        if application != APPLICATION_ID:
            raise ProjectError("This SQLite file is not a WMLSTudio project.")
        if version != SCHEMA_VERSION:
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

    def add_sample(self, path: str | Path, name: str | None = None) -> str:
        source = Path(path).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"Sequence input does not exist: {source}")
        sample_name = source.name if name is None else name.strip()
        if not sample_name:
            raise ValueError("Sample name cannot be empty.")
        sample_id = uuid.uuid4().hex
        timestamp = _now()
        with self.transaction():
            self._connection.execute(
                "INSERT INTO samples "
                "(id, name, input_path, status, created_at, updated_at) VALUES (?, ?, ?, ?, ?, ?)",
                (sample_id, sample_name, str(source), "queued", timestamp, timestamp),
            )
        return sample_id

    @staticmethod
    def _decode(row: sqlite3.Row) -> dict[str, Any]:
        sample = dict(row)
        sample["metadata"] = json.loads(sample["metadata"])
        sample["result"] = json.loads(sample["result"]) if sample["result"] is not None else None
        sample["missing_input"] = not Path(sample["input_path"]).is_file()
        return sample

    def samples(self) -> list[dict[str, Any]]:
        with self._lock:
            self._check_open()
            rows = self._connection.execute(
                "SELECT * FROM samples ORDER BY created_at, id"
            ).fetchall()
            return [self._decode(row) for row in rows]

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
        self._update(sample_id, "result = ?, status = 'completed', error = ''", (_json(snapshot),))

    def set_status(self, sample_id: str, status: str, error: str = "") -> None:
        if status not in STATUSES:
            raise ValueError(f"Unknown job status {status!r}.")
        self._update(sample_id, "status = ?, error = ?", (status, str(error)))

    def set_metadata(self, sample_id: str, metadata: Mapping[str, Any]) -> None:
        if not isinstance(metadata, Mapping):
            raise TypeError("Sample metadata must be a mapping.")
        self._update(sample_id, "metadata = ?", (_json(dict(metadata)),))

    def remove_sample(self, sample_id: str) -> None:
        """Remove the stored record and result; never delete the input file."""
        with self.transaction():
            cursor = self._connection.execute("DELETE FROM samples WHERE id = ?", (sample_id,))
            if cursor.rowcount != 1:
                raise KeyError(f"Unknown sample: {sample_id}")

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
