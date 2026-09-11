# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
# Ported from scripts/mlst-make_blast_db:1-25 and bin/mlst:76-95
"""WMLST database updater — PubMLST / Pasteur BIGSdb refresh and BLAST index build.

Implements docs/ARCHITECTURE.md section 4.9 (public API), section 6 (the database
on disk) and section 7 (the update algorithm).

Design points that are contractual rather than incidental:

* Change detection keys on a CONTENT HASH (:func:`content_sha256`), never on the
  upstream ``last_updated`` date — that field advances daily while the anonymous
  payload is frozen (section 7.2).
* Resolution keys on the TRIPLE ``(source, db, scheme_id)``; ``pubmlst_yersinia_seqdef``
  exists on both hosts with different content (section 6.3).
* ``.tfa`` and ``<scheme>.txt`` are verbatim upstream bytes, written in binary mode
  so the text layer can never insert a ``\\r`` (section 6.1).
* Every commit is per-scheme and atomic, so an interrupted run always leaves a
  self-consistent tree (section 7.3).

Runtime dependencies are stdlib only; ``subprocess`` is never imported here —
``makeblastdb`` is launched through :mod:`wmlst.blastbin` (section 4.4).
"""

from __future__ import annotations

import concurrent.futures
import functools
import hashlib
import importlib
import inspect
import json
import logging
import os
import random
import re
import shutil
import ssl
import stat
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone

from .branding import AUTHOR, HOMEPAGE, VENDOR
from .engine import Cancelled, DatabaseMissingError, UpdateError
from .version import __version__

__all__ = [
    "REST_ROOTS",
    "USER_AGENT",
    "SchemeRef",
    "SchemeUpdate",
    "UpdatePlan",
    "apply",
    "build_blast_db",
    "build_manifest",
    "check",
    "content_sha256",
    "db_version",
    "discover_new_schemes",
    "export_bundle",
    "import_bundle",
    "index_is_stale",
    "journal_path",
    "load_manifest",
    "local_content_sha256",
    "manifest_path",
    "read_journal",
    "rollback",
    "scheme_info",
    "stamp_db_version",
    "sweep_staging",
    "write_manifest",
    "write_scheme_info",
]

log = logging.getLogger("wmlst.updatedb")

# --------------------------------------------------------------------------
# 7.1 Endpoints and 7.4 etiquette
# --------------------------------------------------------------------------

REST_ROOTS = {"pubmlst": "https://rest.pubmlst.org",
              "pasteur": "https://bigsdb.pasteur.fr/api"}

USER_AGENT = "WMLST/{version} ({vendor}; {author}; {url})".format(
    version=__version__, vendor=VENDOR, author=AUTHOR, url=HOMEPAGE)

#: In-flight requests per host (PubMLST and Pasteur are counted separately).
MAX_WORKERS_PER_HOST = 4
#: Per-worker politeness delay, seconds.
REQUEST_DELAY = 0.1
#: Attempts per request, and the exponential backoff schedule in seconds.
MAX_ATTEMPTS = 5
BACKOFF_BASE = 1.0
BACKOFF_JITTER = 0.25
#: Only these statuses are retried. 400/401/403/404 are terminal.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
CONNECT_TIMEOUT = 15.0
READ_TIMEOUT = 120.0
#: Consecutive failures against one host before it is paused, and for how long.
BREAKER_FAILURES = 10
BREAKER_PAUSE = 60.0
BREAKER_MAX_PAUSES = 3
#: Every hop of every fetch must stay on one of these hosts, over https. A
#: hostile or compromised API could otherwise redirect the run onto plaintext
#: http, or hand back a ``loci`` list pointing at an internal service (D-note).
ALLOWED_NETLOCS = frozenset(
    urllib.parse.urlsplit(root).netloc for root in REST_ROOTS.values())
#: Hard ceiling on one response body. The largest bundled payload is ~2.4 MB
#: (``helicobacter/atpA.tfa``); this leaves ample headroom while keeping a
#: server that streams forever from exhausting memory.
MAX_BODY = 64 * 1024 * 1024

# --------------------------------------------------------------------------
# 6.x on-disk layout
# --------------------------------------------------------------------------

PUBMLST_DIRNAME = "pubmlst"
BLAST_DIRNAME = "blast"
BLAST_FASTA = "mlst.fa"
VERSION_FILE = "VERSION.txt"
MANIFEST_FILE = "schemes.manifest.tsv"
SPECIES_MAP_FILE = "scheme_species_map.tab"
JOURNAL_FILE = ".wmlst-update.jsonl"
STAGING_DIRNAME = ".staging"
ROLLBACK_DIRNAME = ".rollback"
#: ``blast`` is installed by renaming a SIBLING directory over it, so a crash or
#: a failed rename can never leave a half-replaced index behind (section 7.3).
BLAST_STAGING_PREFIX = BLAST_DIRNAME + ".staging."
BLAST_OLD_PREFIX = BLAST_DIRNAME + ".old."
INFO_SUFFIX = "_info.json"
PROFILE_SUFFIX = ".txt"
ALLELE_SUFFIX = ".tfa"
VERSION_STAMP_FILE = "database_version.txt"

MANIFEST_HEADER = "#SCHEME\tSOURCE\tDB\tSCHEME_ID\tNLOCI\tALIAS_OF"

#: Written into a manifest row whose upstream identity could not be established.
UNRESOLVED = "unresolved"

#: The 12 files a v5 ``makeblastdb -hash_index -parse_seqids`` run must produce.
BLAST_INDEX_EXTENSIONS = ("ndb", "nhd", "nhi", "nhr", "nin", "njs",
                          "nog", "nos", "not", "nsq", "ntf", "nto")

#: The six non-locus column names of a profile header (MLST/Scheme.pm:46).
PROFILE_NON_LOCUS = frozenset(
    {"ST", "mlst_clade", "clonal_complex", "species", "CC", "Lineage"})

#: ``last_updated`` for a scheme upstream never versioned (e.g. mgenitalium).
NO_VERSION = "No version information available"

#: NTFS-hostile characters. A ``:`` would silently create an alternate data stream.
# Windows reserves these in a filename; NUL and the other control codes are
# rejected too, because open() raises on an embedded NUL rather than returning a
# clean error, and a control byte in a path is never a legitimate scheme name.
_FORBIDDEN_IN_LOCUS = set(':*?"<>|\\/') | {chr(c) for c in range(0x20)} | {"\x7f"}

# Windows refuses to create a file with one of these stems, with or without an
# extension, so a scheme so named would fail only once it reached a real install.
_RESERVED_WINDOWS_STEMS = frozenset(
    ["CON", "PRN", "AUX", "NUL"]
    + ["COM%d" % i for i in range(1, 10)]
    + ["LPT%d" % i for i in range(1, 10)]
)

_API_RE = re.compile(
    r"^https://(?:rest\.pubmlst\.org|bigsdb\.pasteur\.fr/api)"
    r"/db/([^/]+)/schemes/(\d+)$", re.ASCII)

_LOCUS_WORD_RE = re.compile(r"^\w+$", re.ASCII)

#: A zip member such as ``C:/evil.txt`` is drive-relative on Windows and escapes
#: any directory it is joined onto; ``C:\\Windows\\...`` is outright absolute.
_DRIVE_RE = re.compile(r"^[A-Za-z]:", re.ASCII)

#: The three documented dirname violations of the ``_N`` suffix rule (6.3).
#: Hardcoded because they cannot be inferred; used only to build a manifest when
#: a bundled ``_info.json`` is missing, never to override one that is present.
DIRNAME_OVERRIDES = {
    "salmonella": ("pubmlst", "pubmlst_salmonella_seqdef", "2"),
    "mgenitalium": ("pubmlst", "pubmlst_mgenitalium_seqdef", "2"),
    "cdiphtheriae": ("pasteur", "pubmlst_diphtheria_seqdef", "3"),
}

#: Schemes whose database is absent from the host's ``/db`` group listing yet
#: answers HTTP 200 when addressed directly (6.3). Never prune these.
UNLISTED_DATABASES = frozenset({
    "pubmlst_ecoli_achtman_seqdef", "pubmlst_halobacteria_seqdef",
    "pubmlst_mamphoriforme_seqdef", "pubmlst_mcatarrhalis_achtman_seqdef",
    "pubmlst_senterica_achtman_seqdef", "pubmlst_streptothermophilus_seqdef",
})

#: 7.5 discovery filters.
_MLST_DESC_RE = re.compile(r"\bMLST\b", re.IGNORECASE)
_NOT_MLST_RE = re.compile(r"(cgMLST|rMLST|wgMLST|Core\s+\d+|TCS-\d+)", re.IGNORECASE)
DISCOVERY_LOCUS_RANGE = (2, 15)

if sys.version_info >= (3, 10):
    def _frozen(cls):
        return dataclass(frozen=True, slots=True)(cls)
else:  # pragma: no cover - exercised only on 3.9
    def _frozen(cls):
        return dataclass(frozen=True)(cls)


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def _utc_now():
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _utc_now().replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _today() -> str:
    return _utc_now().strftime("%Y-%m-%d")


@functools.lru_cache(maxsize=64)
def _callback_arity(callback) -> int:
    """How many positional arguments a progress hook wants (2 or 3)."""
    try:
        parameters = inspect.signature(callback).parameters.values()
    except (TypeError, ValueError):
        return 2
    count = 0
    for parameter in parameters:
        if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
            return 3
        if parameter.kind in (inspect.Parameter.POSITIONAL_ONLY,
                              inspect.Parameter.POSITIONAL_OR_KEYWORD):
            count += 1
    return 3 if count >= 3 else 2


def _emit(progress, fraction: float, message: str) -> None:
    """Report progress as a fraction in [0, 1] plus a human message.

    Two shapes are supported and told apart by the callback's own signature:
    ``(fraction, message)``, and the ``(done, total, text)`` shape that
    ``blastbin.bootstrap`` and the GUI's task runner already use — the latter
    receives ``(fraction, 1.0, message)``, which renders as the same percentage.
    A callback that raises is logged and ignored: a broken GUI hook must never
    abort a database update.
    """
    if progress is None:
        return
    fraction = max(0.0, min(1.0, float(fraction)))
    try:
        if _callback_arity(progress) >= 3:
            progress(fraction, 1.0, message)
        else:
            progress(fraction, message)
    except Exception:  # pragma: no cover - defensive
        log.debug("progress callback raised", exc_info=True)


def _cancelled(cancel) -> bool:
    if cancel is None:
        return False
    is_set = getattr(cancel, "is_set", None)
    return bool(is_set()) if callable(is_set) else bool(cancel)


def _check_cancel(cancel) -> None:
    if _cancelled(cancel):
        raise Cancelled("The database update was cancelled.")


def _pubmlst_dir(dbdir: str) -> str:
    return os.path.join(dbdir, PUBMLST_DIRNAME)


def _scheme_dir(dbdir: str, name: str) -> str:
    return os.path.join(_pubmlst_dir(dbdir), name)


def journal_path(dbdir: str) -> str:
    """-> the path of the newline-delimited JSON update journal (section 7.3)."""
    return os.path.join(dbdir, JOURNAL_FILE)


def manifest_path(dbdir: str) -> str:
    """-> the path of ``db/schemes.manifest.tsv`` (section 6.3)."""
    return os.path.join(dbdir, MANIFEST_FILE)


def db_version(dbdir: str) -> str:
    """-> the first line of ``db/VERSION.txt``, or '' when it is absent (6.1)."""
    try:
        with open(os.path.join(dbdir, VERSION_FILE), encoding="utf-8") as fh:
            return fh.readline().strip()
    except OSError:
        return ""


def stamp_db_version(dbdir: str, version: str) -> str:
    """Write ``db/VERSION.txt`` as ``<version>\\n`` in binary mode. -> version (6.1)."""
    _write_bytes(os.path.join(dbdir, VERSION_FILE),
                 (version + "\n").encode("utf-8"))
    return version


def _on_rm_error(func, path, _exc) -> None:
    """Clear the read-only bit and retry: Windows refuses to unlink +R files.

    A locked file (antivirus, backup agent, Explorer preview) usually frees up
    within a few hundred milliseconds, so the retry loop is worth its cost.
    """
    for _ in range(5):
        try:
            os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
        except OSError:
            pass
        try:
            func(path)
            return
        except OSError:
            time.sleep(0.1)
    log.warning("could not remove %s", path)


def _rmtree(path: str, *, required: bool = False) -> None:
    """``shutil.rmtree`` hardened for Windows (mirrors ``blastbin._rmtree``).

    ``shutil.rmtree(..., ignore_errors=True)`` gives up on the first ``EACCES``
    and says nothing, so a read-only or briefly locked survivor is silently
    reused as if it were a fresh directory. With ``required=True`` a directory
    that still exists afterwards is an actionable error instead.
    """
    if os.path.isdir(path):
        if sys.version_info >= (3, 12):
            shutil.rmtree(path, onexc=lambda f, p, e: _on_rm_error(f, p, e))
        else:  # pragma: no cover - exercised on 3.9-3.11
            shutil.rmtree(path, onerror=_on_rm_error)
    if required and os.path.isdir(path):
        raise UpdateError(
            "Could not clear the staging folder %s. Close any program holding "
            "files there (antivirus, backup, Explorer), clear the read-only "
            "attribute, delete the folder and try again." % path)


def _write_bytes(path: str, data: bytes) -> None:
    """Atomic-ish write: temp file in the same directory, then :func:`os.replace`."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    tmp = os.path.join(directory, ".wmlst-tmp-%d-%s" % (os.getpid(), os.path.basename(path)))
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _dir_size(path: str) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return total


def _scheme_names_on_disk(dbdir: str):
    root = _pubmlst_dir(dbdir)
    if not os.path.isdir(root):
        raise DatabaseMissingError(
            "Database directory does not exist: %s" % root)
    return tuple(sorted(
        (n for n in os.listdir(root)
         if not n.startswith(".") and not n.startswith("_")
         and os.path.isdir(os.path.join(root, n))),
        key=lambda s: s.encode("utf-8")))


def _tfa_names(scheme_dir: str):
    return tuple(sorted(
        (n for n in os.listdir(scheme_dir) if n.endswith(ALLELE_SUFFIX)),
        key=lambda s: s.encode("utf-8")))


def _is_within(base: str, path: str) -> bool:
    """True when ``path`` resolves inside ``base`` (never follows to a parent)."""
    base = os.path.abspath(base)
    path = os.path.abspath(path)
    return path == base or path.startswith(base + os.sep) or (
        os.altsep is not None and path.startswith(base + os.altsep))


def _locus_of_url(url: str) -> str:
    """-> the locus name carried by a BIGSdb ``/loci/<name>`` URL, percent-decoded."""
    return urllib.parse.unquote(url.rstrip("/").rsplit("/", 1)[-1])


def _field_names(meta) -> tuple:
    return tuple(_locus_of_url(u) for u in meta.get("fields", ()) or ())


def _derive_dirname(db: str, scheme_id: str) -> str:
    """The ``_N`` suffix rule of section 6.3 — a VALIDATION aid, never the resolver."""
    base = db[len("pubmlst_"):] if db.startswith("pubmlst_") else db
    if base.endswith("_seqdef"):
        base = base[:-len("_seqdef")]
    return base if scheme_id == "1" else "%s_%s" % (base, scheme_id)


# --------------------------------------------------------------------------
# 4.9 dataclasses
# --------------------------------------------------------------------------

@_frozen
class SchemeRef:
    """One row of ``db/schemes.manifest.tsv`` (sections 4.9, 6.3).

    ``alias_of`` names the canonical directory when two directories hold the
    same upstream payload (``cdiphtheriae`` -> ``diphtheria_3``).
    """

    name: str
    source: str
    db: str
    scheme_id: str
    nloci: int
    alias_of: str | None

    @property
    def resolved(self) -> bool:
        """False for a manifest row whose upstream identity is unknown."""
        return self.source in REST_ROOTS and bool(self.db) and bool(self.scheme_id)

    @property
    def api(self) -> str:
        """-> the scheme document URL (section 7.1). '' when unresolved."""
        if not self.resolved:
            return ""
        return "%s/db/%s/schemes/%s" % (REST_ROOTS[self.source], self.db, self.scheme_id)

    @property
    def triple(self) -> tuple:
        """-> ``(source, db, scheme_id)`` — the only legal resolution key (6.3)."""
        return (self.source, self.db, self.scheme_id)


@_frozen
class SchemeUpdate:
    """The verdict for one scheme after a metadata pass (section 4.9)."""

    name: str
    local_date: str
    remote_date: str
    status: str
    bytes_estimate: int
    added_types: int | None
    detail: str


@_frozen
class UpdatePlan:
    """Everything :func:`check` learned; the input to :func:`apply` (section 4.9)."""

    updates: tuple
    new_schemes: tuple
    total_bytes: int
    checked_at: str
    run_id: str

    def by_name(self, name: str):
        """-> the :class:`SchemeUpdate` for ``name``, or None."""
        for upd in self.updates:
            if upd.name == name:
                return upd
        return None

    @property
    def changed(self) -> tuple:
        """-> the updates whose status is ``changed``."""
        return tuple(u for u in self.updates if u.status == "changed")


# --------------------------------------------------------------------------
# HTTP layer — 7.4 etiquette
# --------------------------------------------------------------------------

class _Response:
    __slots__ = ("body", "headers", "status", "url")

    def __init__(self, status, headers, body, url):
        self.status = status
        self.headers = headers
        self.body = body
        self.url = url

    @property
    def content_type(self) -> str:
        return (self.headers.get("Content-Type") or "").lower()


class _HostState:
    __slots__ = ("failures", "lock", "paused_until", "pauses")

    def __init__(self):
        self.lock = threading.Lock()
        self.failures = 0
        self.paused_until = 0.0
        self.pauses = 0


def check_fetch_url(url: str) -> str:
    """Refuse any URL that is not https on a known BIGSdb host. -> ``url``.

    Applied to the FIRST request and to every redirect hop, and to the ``loci``
    URLs the scheme document hands back, which are attacker-controlled data:
    without it an https fetch can finish over cleartext http, or on an internal
    host of the server's choosing.
    """
    parts = urllib.parse.urlsplit(url)
    if parts.scheme != "https" or parts.netloc not in ALLOWED_NETLOCS:
        raise UpdateError(
            "Refusing to fetch %s: only https on %s is allowed."
            % (url, ", ".join(sorted(ALLOWED_NETLOCS))))
    return url


class _PinnedRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Re-apply :func:`check_fetch_url` to every 30x target (section 7.4)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_fetch_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


class _Fetcher:
    """Polite, retrying, cancellable HTTP GET pool (section 7.4).

    4 in-flight requests per host, 100 ms per-worker delay, 5 attempts with
    1/2/4/8/16 s exponential backoff and +-25 % jitter, ``Retry-After`` honoured,
    a 10-consecutive-failure circuit breaker that pauses a host for 60 s and
    aborts the run on the third pause. TLS verification is always on.
    """

    def __init__(self, *, cancel=None, user_agent: str = USER_AGENT,
                 workers: int = MAX_WORKERS_PER_HOST, delay: float = REQUEST_DELAY,
                 opener=None, sleep=time.sleep):
        self.cancel = cancel
        self.user_agent = user_agent
        self.workers = max(1, int(workers))
        self.delay = float(delay)
        self._sleep = sleep
        self._hosts = {}
        self._hosts_lock = threading.Lock()
        self.bytes_read = 0
        self.requests = 0
        self._counter_lock = threading.Lock()
        if opener is not None:
            self._opener = opener
        else:
            context = ssl.create_default_context()
            self._opener = urllib.request.build_opener(
                urllib.request.HTTPSHandler(context=context),
                _PinnedRedirectHandler())
            self._opener.addheaders = []

    # -- internals ---------------------------------------------------------
    def _state(self, host: str) -> _HostState:
        with self._hosts_lock:
            state = self._hosts.get(host)
            if state is None:
                state = self._hosts[host] = _HostState()
            return state

    def _wait_for_host(self, host: str) -> None:
        state = self._state(host)
        while True:
            _check_cancel(self.cancel)
            with state.lock:
                remaining = state.paused_until - time.monotonic()
            if remaining <= 0:
                return
            self._sleep(min(remaining, 1.0))

    def _note_success(self, host: str) -> None:
        state = self._state(host)
        with state.lock:
            state.failures = 0

    def _note_failure(self, host: str) -> None:
        state = self._state(host)
        with state.lock:
            state.failures += 1
            if state.failures < BREAKER_FAILURES:
                return
            state.failures = 0
            state.pauses += 1
            pauses = state.pauses
            state.paused_until = time.monotonic() + BREAKER_PAUSE
        if pauses >= BREAKER_MAX_PAUSES:
            raise UpdateError(
                "%s is not responding reliably; the update was stopped and can be "
                "resumed later." % host)
        log.warning("pausing %s for %.0f s after repeated failures", host, BREAKER_PAUSE)

    def _open(self, url: str) -> _Response:
        check_fetch_url(url)
        request = urllib.request.Request(url, headers={
            "User-Agent": self.user_agent,
            "Accept": "*/*",
            "Accept-Encoding": "identity",
        })
        # urllib applies one timeout to connect and to each read; the read
        # budget is the binding one (largest payload ~2.4 MB).
        with self._opener.open(request, timeout=READ_TIMEOUT) as handle:
            body = handle.read(MAX_BODY + 1)
            if len(body) > MAX_BODY:
                raise UpdateError(
                    "%s returned more than %d bytes; refusing to buffer it."
                    % (url, MAX_BODY))
            status = getattr(handle, "status", None) or handle.getcode()
            return _Response(int(status), dict(handle.headers), body, url)

    # -- public ------------------------------------------------------------
    def get(self, url: str) -> _Response:
        """GET ``url`` with the full retry / breaker policy. -> :class:`_Response`."""
        host = urllib.parse.urlsplit(url).netloc
        last = ""
        for attempt in range(1, MAX_ATTEMPTS + 1):
            _check_cancel(self.cancel)
            self._wait_for_host(host)
            if self.delay:
                self._sleep(self.delay)
            retry_after = None
            try:
                response = self._open(url)
            except urllib.error.HTTPError as exc:
                body = b""
                try:
                    body = exc.read()
                except Exception:  # pragma: no cover - defensive
                    pass
                if exc.code not in RETRY_STATUS:
                    raise HttpStatusError(exc.code, url, body) from exc
                last = "HTTP %d" % exc.code
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                self._note_failure(host)
            except (urllib.error.URLError, OSError, TimeoutError) as exc:
                last = str(getattr(exc, "reason", None) or exc)
                self._note_failure(host)
            else:
                self._note_success(host)
                with self._counter_lock:
                    self.requests += 1
                    self.bytes_read += len(response.body)
                return response
            if attempt == MAX_ATTEMPTS:
                break
            delay = BACKOFF_BASE * (2 ** (attempt - 1))
            if retry_after:
                try:
                    delay = max(delay, float(retry_after))
                except ValueError:
                    pass
            delay *= 1.0 + random.uniform(-BACKOFF_JITTER, BACKOFF_JITTER)
            log.debug("retry %d/%d for %s in %.1fs (%s)",
                      attempt, MAX_ATTEMPTS, url, delay, last)
            self._sleep(delay)
        raise UpdateError(
            "Could not download %s after %d attempts (%s). Check your internet "
            "connection and try again." % (url, MAX_ATTEMPTS, last or "unknown error"))

    def get_json(self, url: str):
        """GET ``url`` and parse it as JSON. -> the decoded object."""
        response = self.get(url)
        try:
            return json.loads(response.body.decode("utf-8"))
        except (UnicodeDecodeError, ValueError) as exc:
            raise UpdateError(
                "%s did not return valid JSON (%s)." % (url, exc)) from exc

    def gather(self, urls, *, progress=None, base: float = 0.0, span: float = 0.0,
               message: str = "Downloading"):
        """GET every URL, 4 workers per host. -> ``dict[url] = _Response``."""
        urls = list(dict.fromkeys(urls))
        if not urls:
            return {}
        hosts = {urllib.parse.urlsplit(u).netloc for u in urls}
        results = {}
        done = 0
        limits = {h: threading.Semaphore(self.workers) for h in hosts}

        def task(url):
            host = urllib.parse.urlsplit(url).netloc
            with limits[host]:
                return url, self.get(url)

        max_workers = self.workers * max(1, len(hosts))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(task, u) for u in urls]
            try:
                for future in concurrent.futures.as_completed(futures):
                    url, response = future.result()
                    results[url] = response
                    done += 1
                    if span:
                        _emit(progress, base + span * done / len(urls),
                              "%s (%d/%d)" % (message, done, len(urls)))
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        return results


class HttpStatusError(UpdateError):
    """A terminal (never retried) HTTP status: 400/401/403/404 and friends."""

    def __init__(self, code: int, url: str, body: bytes = b""):
        self.code = int(code)
        self.url = url
        self.body = body
        text = ""
        try:
            payload = json.loads(body.decode("utf-8"))
            if isinstance(payload, dict):
                text = str(payload.get("message", ""))
        except Exception:
            text = body[:200].decode("utf-8", "replace")
        self.message_text = text
        super().__init__("HTTP %d for %s%s" % (code, url, (": " + text) if text else ""))
        self.user_message = (
            "The database server refused a request (HTTP %d). Nothing was changed."
            % self.code)

    @property
    def is_undefined_scheme(self) -> bool:
        """True for the 404 that means 'this scheme no longer exists upstream'."""
        return self.code == 404 and "has not been defined" in self.message_text


# --------------------------------------------------------------------------
# 4.9 content_sha256
# --------------------------------------------------------------------------

def content_sha256(profiles_csv: bytes, alleles) -> str:
    """The change oracle of section 7.2 / 4.9.

    ``h = sha256(); h.update(b'WMLST-DBv1\\n'); h.update(profiles_csv)`` then, for
    each locus in byte-sorted order, ``h.update(locus + b'\\x00'); h.update(bytes)``.
    Upstream's ``last_updated`` advances daily while the anonymous payload is
    frozen, so only this hash can decide whether anything really changed.
    """
    digest = hashlib.sha256()
    digest.update(b"WMLST-DBv1\n")
    digest.update(profiles_csv)
    for locus in sorted(alleles, key=lambda s: s.encode("utf-8")):
        digest.update(locus.encode("utf-8") + b"\x00")
        digest.update(alleles[locus])
    return digest.hexdigest()


def local_content_sha256(dbdir: str, name: str) -> str | None:
    """Recompute :func:`content_sha256` from the files on disk (sections 4.9, 7.2).

    Authoritative even for a pristine bundle that carries no
    ``wmlst_content_sha256`` ledger key yet. -> None when the scheme is absent.
    """
    scheme_dir = _scheme_dir(dbdir, name)
    profiles = os.path.join(scheme_dir, name + PROFILE_SUFFIX)
    if not os.path.isdir(scheme_dir) or not os.path.isfile(profiles):
        return None
    with open(profiles, "rb") as fh:
        profiles_csv = fh.read()
    alleles = {}
    for filename in _tfa_names(scheme_dir):
        with open(os.path.join(scheme_dir, filename), "rb") as fh:
            alleles[filename[:-len(ALLELE_SUFFIX)]] = fh.read()
    return content_sha256(profiles_csv, alleles)


# --------------------------------------------------------------------------
# 6.2 _info.json
# --------------------------------------------------------------------------

def scheme_info(dbdir: str, name: str) -> dict:
    """-> the parsed ``<scheme>_info.json``, or ``{}`` when missing (section 6.2)."""
    path = os.path.join(_scheme_dir(dbdir, name), name + INFO_SUFFIX)
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _info_bytes(info: dict) -> bytes:
    """Serialise ``_info.json``: 8 canonical keys first, 2-space indent, LF (6.2).

    ``ensure_ascii=True`` because the bundled files escape their non-ASCII
    descriptions (``MLST#2 (Bek\u0151)``); writing raw UTF-8 there would make
    every re-download differ from the shipped bytes for no benefit.
    """
    ordered = {}
    for key in ("name", "description", "locus", "download_date", "last_updated",
                "source", "API", "authenticated"):
        if key in info:
            ordered[key] = info[key]
    for key, value in info.items():
        if key not in ordered:
            ordered[key] = value
    text = json.dumps(ordered, indent=2, ensure_ascii=True) + "\n"
    return text.encode("utf-8")


def write_scheme_info(path: str, info: dict) -> None:
    """Write a ``<scheme>_info.json`` in the bundled byte format (section 6.2)."""
    _write_bytes(path, _info_bytes(info))


def _build_info(ref: SchemeRef, meta: dict, *, digest: str, remote_date: str,
                previous: dict, unparseable) -> dict:
    """Assemble the 8 canonical keys plus the four WMLST ledger keys (6.2)."""
    info = {
        "name": ref.name,
        "description": meta.get("description") or previous.get("description", "MLST"),
        "locus": int(meta.get("locus_count", ref.nloci)),
        "download_date": _today(),
        "last_updated": remote_date,
        "source": ref.source,
        "API": ref.api,
        "authenticated": bool(previous.get("authenticated", False)),
        "wmlst_content_sha256": digest,
        "wmlst_verified_unchanged_at": _now_iso(),
        "wmlst_upstream_last_updated_at_verification": remote_date,
    }
    if unparseable:
        info["wmlst_unparseable_loci"] = list(unparseable)
    return info


# --------------------------------------------------------------------------
# 6.3 the manifest
# --------------------------------------------------------------------------

def _ref_from_info(name: str, info: dict) -> SchemeRef:
    match = _API_RE.match(str(info.get("API", "")))
    if match:
        db, scheme_id = match.group(1), match.group(2)
        source = "pasteur" if "pasteur" in info["API"] else "pubmlst"
        declared = str(info.get("source") or source)
        if declared != source:
            log.warning("%s: _info.json source %r disagrees with its API host; using %r",
                        name, declared, source)
        nloci = int(info.get("locus") or 0)
        return SchemeRef(name, source, db, scheme_id, nloci, None)
    if name in DIRNAME_OVERRIDES:
        source, db, scheme_id = DIRNAME_OVERRIDES[name]
        return SchemeRef(name, source, db, scheme_id, int(info.get("locus") or 0), None)
    return SchemeRef(name, UNRESOLVED, "", "", int(info.get("locus") or 0), None)


def _apply_aliases(refs):
    """Point every duplicate ``(source, db, scheme_id)`` at one canonical name.

    The canonical member is the directory the ``_N`` rule would generate; the
    others become aliases and are fetched once, materialised many times (7.3).
    """
    groups = {}
    for ref in refs:
        if ref.resolved:
            groups.setdefault(ref.triple, []).append(ref.name)
    canonical = {}
    for triple, names in groups.items():
        if len(names) < 2:
            continue
        expected = _derive_dirname(triple[1], triple[2])
        names = sorted(names, key=lambda s: s.encode("utf-8"))
        primary = expected if expected in names else names[0]
        for other in names:
            if other != primary:
                canonical[other] = primary
    out = []
    for ref in refs:
        alias = canonical.get(ref.name)
        out.append(ref if alias is None else SchemeRef(
            ref.name, ref.source, ref.db, ref.scheme_id, ref.nloci, alias))
    return tuple(out)


def load_manifest(dbdir: str) -> tuple:
    """-> ``tuple[SchemeRef, ...]`` from ``db/schemes.manifest.tsv`` (section 4.9).

    Falls back to parsing every ``<scheme>_info.json`` ``API`` field with
    ``^https://(?:rest\\.pubmlst\\.org|bigsdb\\.pasteur\\.fr/api)/db/([^/]+)/schemes/(\\d+)$``.
    A scheme directory present on disk but absent from the manifest is derived
    the same way, so a hand-added directory still resolves. Resolution keys on
    the TRIPLE ``(source, db, scheme_id)``, never on ``db`` alone: the database
    ``pubmlst_yersinia_seqdef`` exists on both hosts with different content.
    """
    on_disk = set(_scheme_names_on_disk(dbdir))
    refs = {}
    path = manifest_path(dbdir)
    if os.path.isfile(path):
        with open(path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) != 6:
                    raise UpdateError(
                        "%s line %d is not a 6-column row." % (path, lineno))
                name, source, db, scheme_id, nloci, alias = (p.strip() for p in parts)
                _validate_scheme_name(name, "%s line %d" % (path, lineno))
                try:
                    nloci_int = int(nloci) if nloci not in ("", "-") else 0
                except ValueError as exc:
                    raise UpdateError(
                        "%s line %d has a non-numeric NLOCI."
                        % (path, lineno)) from exc
                refs[name] = SchemeRef(
                    name,
                    source if source != UNRESOLVED else UNRESOLVED,
                    "" if db == "-" else db,
                    "" if scheme_id == "-" else scheme_id,
                    nloci_int,
                    None if alias in ("", "-") else alias)
    for name in sorted(on_disk, key=lambda s: s.encode("utf-8")):
        if name not in refs:
            refs[name] = _ref_from_info(name, scheme_info(dbdir, name))
    ordered = [refs[n] for n in sorted(refs, key=lambda s: s.encode("utf-8"))]
    if any(r.alias_of for r in ordered):
        return tuple(ordered)
    return _apply_aliases(ordered)


def build_manifest(dbdir: str, *, verify: bool = False, progress=None,
                   cancel=None, fetcher=None) -> tuple:
    """Derive the manifest from the bundled ``_info.json`` files (section 6.3).

    With ``verify=True`` every triple is confirmed against the live API and
    ``NLOCI`` is taken from the upstream ``locus_count``; a scheme that cannot be
    resolved is recorded as ``unresolved`` rather than guessed. Absence from the
    host's ``/db`` listing is never a reason to prune — only a 404 is, and even
    then the local directory is left untouched.
    """
    names = _scheme_names_on_disk(dbdir)
    refs = _apply_aliases(tuple(
        _ref_from_info(name, scheme_info(dbdir, name)) for name in names))
    if not verify:
        return refs
    fetcher = fetcher or _Fetcher(cancel=cancel)
    urls = {}
    for ref in refs:
        if ref.resolved and ref.alias_of is None:
            urls.setdefault(ref.api, []).append(ref)
    _emit(progress, 0.0, "Verifying %d scheme identities" % len(urls))
    responses = fetcher.gather(list(urls), progress=progress, base=0.0, span=1.0,
                               message="Verifying schemes")
    meta_by_api = {}
    for url, response in responses.items():
        try:
            meta_by_api[url] = json.loads(response.body.decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            meta_by_api[url] = None
    out = []
    for ref in refs:
        source_ref = ref
        if ref.alias_of:
            source_ref = next(r for r in refs if r.name == ref.alias_of)
        meta = meta_by_api.get(source_ref.api)
        if meta is None:
            out.append(ref)
            continue
        nloci = int(meta.get("locus_count") or ref.nloci)
        out.append(SchemeRef(ref.name, ref.source, ref.db, ref.scheme_id,
                             nloci, ref.alias_of))
    return tuple(out)


def write_manifest(dbdir: str, refs, path: str | None = None) -> str:
    """Write ``db/schemes.manifest.tsv``: header + one row per scheme (6.3). -> path."""
    target = path or manifest_path(dbdir)
    lines = [MANIFEST_HEADER]
    for ref in sorted(refs, key=lambda r: r.name.encode("utf-8")):
        lines.append("\t".join((
            ref.name,
            ref.source if ref.source else UNRESOLVED,
            ref.db or "-",
            str(ref.scheme_id) if ref.scheme_id else "-",
            str(int(ref.nloci or 0)),
            ref.alias_of or "-")))
    _write_bytes(target, ("\n".join(lines) + "\n").encode("utf-8"))
    return target


# --------------------------------------------------------------------------
# 7.3 journal
# --------------------------------------------------------------------------

class _Journal:
    """Newline-delimited JSON, fsynced after every record (section 7.3)."""

    def __init__(self, path: str, run_id: str):
        self.path = path
        self.run_id = run_id
        self._lock = threading.Lock()

    def record(self, event: str, **fields) -> None:
        payload = {"ts": _now_iso(), "run_id": self.run_id, "event": event}
        payload.update(fields)
        line = json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n"
        with self._lock:
            try:
                os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
                with open(self.path, "ab") as fh:
                    fh.write(line.encode("utf-8"))
                    fh.flush()
                    os.fsync(fh.fileno())
            except OSError:  # pragma: no cover - a journal must never abort a run
                log.debug("could not append to the update journal", exc_info=True)


def read_journal(dbdir: str) -> tuple:
    """-> ``tuple[dict, ...]`` of journal records, oldest first (section 7.3)."""
    path = journal_path(dbdir)
    out = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    out.append(json.loads(line))
                except ValueError:
                    continue
    except OSError:
        return ()
    return tuple(out)


def committed_schemes(dbdir: str, run_id: str | None = None) -> frozenset:
    """-> the names already committed, for ``--resume`` (section 7.3)."""
    names = set()
    for record in read_journal(dbdir):
        if run_id is not None and record.get("run_id") != run_id:
            continue
        if record.get("event") == "scheme_commit":
            names.add(record.get("scheme"))
    names.discard(None)
    return frozenset(names)


def sweep_staging(dbdir: str) -> tuple:
    """Repair a crashed run: finish half-done swaps, delete orphans (section 7.3).

    A crash between ``live -> <scheme>.old.<pid>`` and ``staging -> <scheme>``
    leaves the live directory missing; that case is repaired by renaming the
    ``.old`` copy back. -> the tuple of paths removed or restored.
    """
    touched = []
    for entry in sorted(os.listdir(dbdir) if os.path.isdir(dbdir) else ()):
        path = os.path.join(dbdir, entry)
        if not os.path.isdir(path):
            continue
        if entry.startswith(BLAST_OLD_PREFIX):
            live = os.path.join(dbdir, BLAST_DIRNAME)
            if not os.path.isdir(live):
                os.replace(path, live)
                touched.append(live)
            else:
                _rmtree(path)
                touched.append(path)
        elif entry.startswith(BLAST_STAGING_PREFIX):
            _rmtree(path)
            touched.append(path)
    root = _pubmlst_dir(dbdir)
    if not os.path.isdir(root):
        return tuple(touched)
    for entry in sorted(os.listdir(root)):
        path = os.path.join(root, entry)
        if not os.path.isdir(path):
            continue
        if ".old." in entry:
            live = os.path.join(root, entry.split(".old.", 1)[0])
            if not os.path.isdir(live):
                os.replace(path, live)
                touched.append(live)
            else:
                _rmtree(path)
                touched.append(path)
    staging = os.path.join(root, STAGING_DIRNAME)
    if os.path.isdir(staging):
        for entry in sorted(os.listdir(staging)):
            target = os.path.join(staging, entry)
            _rmtree(target)
            touched.append(target)
    return tuple(touched)


# --------------------------------------------------------------------------
# 7.6 validation gates
# --------------------------------------------------------------------------

def _require_text_plain(response: _Response, what: str) -> None:
    if response.status != 200:
        raise UpdateError("%s returned HTTP %d." % (what, response.status))
    if not response.content_type.startswith("text/plain"):
        raise UpdateError("%s returned %r, not text/plain."
                          % (what, response.content_type or "no content type"))
    body = response.body
    if body[:1] == b"{" and b'"message"' in body[:200]:
        raise UpdateError("%s returned a JSON message instead of data." % what)


def _validate_locus_name(locus: str) -> None:
    """Reject a locus name that cannot be a safe filename (section 7.6)."""
    if not locus:
        raise UpdateError("A locus with an empty name was returned upstream.")
    bad = set(locus) & _FORBIDDEN_IN_LOCUS
    if bad:
        raise UpdateError(
            "Locus %r contains a character that is illegal in a Windows filename "
            "(%s)." % (locus, "".join(sorted(bad))))
    if locus in (".", "..") or locus.startswith("."):
        raise UpdateError("Locus %r is not a usable file name." % locus)


def _validate_scheme_name(name: str, where: str = "") -> None:
    """Reject a scheme name that cannot be one safe path component (7.6).

    The name becomes a directory under ``db/pubmlst`` and the stem of
    ``<scheme>.txt`` / ``<scheme>_info.json``, so it is a path component in
    everything :func:`_materialise` writes. ``db/schemes.manifest.tsv`` is an
    ordinary user-writable file, and ``--dbdir`` / ``%WMLST_DBDIR%`` can point
    at a tree the user did not author, so the NAME column is untrusted input.
    """
    prefix = ("%s: " % where) if where else ""
    if not name:
        raise UpdateError("%sa scheme with an empty name was supplied." % prefix)
    bad = set(name) & _FORBIDDEN_IN_LOCUS
    if bad:
        raise UpdateError(
            "%sscheme %r contains a character that is illegal in a Windows "
            "filename (%s)." % (prefix, name, "".join(sorted(bad))))
    if name.startswith(".") or name.startswith("_"):
        raise UpdateError(
            "%sscheme %r is not a usable directory name." % (prefix, name))
    if os.path.isabs(name) or os.path.basename(name) != name:
        raise UpdateError(
            "%sscheme %r is not a single path component." % (prefix, name))
    if name.split(".")[0].upper() in _RESERVED_WINDOWS_STEMS:
        raise UpdateError(
            "%sscheme %r is a reserved Windows device name." % (prefix, name))
    if name.endswith((" ", ".")):
        # Windows silently strips these, so two distinct names would collide.
        raise UpdateError(
            "%sscheme %r ends with a space or a dot." % (prefix, name))


def validate_tfa(locus: str, response: _Response) -> None:
    """Every gate of section 7.6 for one ``alleles_fasta`` payload."""
    what = "The allele file for %s" % locus
    _require_text_plain(response, what)
    body = response.body
    if not body:
        raise UpdateError("%s is empty." % what)
    if body[:1] != b">":
        raise UpdateError("%s does not start with '>'." % what)
    if b"\r" in body:
        raise UpdateError("%s contains a carriage return." % what)
    header_re = re.compile(r"^>" + re.escape(locus) + r"_\S+$", re.ASCII)
    records = 0
    for raw in body.split(b"\n"):
        if raw == b"" :
            continue
        if raw.startswith(b"<"):
            raise UpdateError("%s looks like an HTML page, not FASTA." % what)
        if raw.startswith(b">"):
            records += 1
            line = raw.decode("utf-8", "replace")
            if not header_re.match(line):
                raise UpdateError("%s has a malformed header %r." % (what, line))
    if body.count(b"\n\n"):
        raise UpdateError("%s contains a blank line." % what)
    if records == 0:
        raise UpdateError("%s contains no sequences." % what)


def validate_profiles(name: str, response: _Response, loci, fields) -> int:
    """Gates of section 7.6 for ``<scheme>.txt``. -> the number of data rows.

    Zero data rows is legal (``mgenitalium`` has a 59-byte header-only table).
    """
    what = "The profile table for %s" % name
    _require_text_plain(response, what)
    body = response.body
    if b"\r" in body:
        raise UpdateError("%s contains a carriage return." % what)
    lines = body.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    if not lines:
        raise UpdateError("%s is empty." % what)
    header = lines[0].decode("utf-8", "replace").split("\t")
    if len(header) < 2 or header[0] != "ST":
        raise UpdateError("%s does not begin with an 'ST' column." % what)
    declared = set(header) - PROFILE_NON_LOCUS - set(fields)
    if declared != set(loci):
        missing = sorted(set(loci) - declared)
        extra = sorted(declared - set(loci))
        raise UpdateError(
            "%s does not match the scheme definition (missing %s, unexpected %s)."
            % (what, missing or "nothing", extra or "nothing"))
    for number, raw in enumerate(lines[1:], 2):
        if raw == b"":
            raise UpdateError("%s has a blank line at line %d." % (what, number))
        if len(raw.decode("utf-8", "replace").split("\t")) != len(header):
            raise UpdateError(
                "%s line %d has the wrong number of columns." % (what, number))
    return len(lines) - 1


def _validate_scheme(name: str, staged_dir: str, meta: dict, loci, live_dir: str,
                     *, allow_shrink: bool) -> None:
    """Per-scheme gates of section 7.6, run just before the atomic swap."""
    tfas = _tfa_names(staged_dir)
    locus_count = int(meta.get("locus_count") or len(loci))
    if len(tfas) != locus_count:
        raise UpdateError("%s: staged %d allele files but upstream declares %d."
                          % (name, len(tfas), locus_count))
    for locus in loci:
        if not os.path.isfile(os.path.join(staged_dir, locus + ALLELE_SUFFIX)):
            raise UpdateError("%s: no allele file was produced for %s." % (name, locus))
    staged_bytes = _dir_size(staged_dir)
    if staged_bytes <= 0:
        raise UpdateError("%s: nothing was downloaded." % name)
    if not allow_shrink and os.path.isdir(live_dir):
        live_bytes = _dir_size(live_dir)
        if live_bytes and staged_bytes < live_bytes // 2:
            raise UpdateError(
                "%s: the download is less than half the size of the installed "
                "scheme (%d vs %d bytes); refusing to replace it."
                % (name, staged_bytes, live_bytes))


def _unparseable_loci(loci):
    """Loci whose names the sseqid regex cannot parse (section 5.19(b))."""
    return tuple(line for line in loci if not _LOCUS_WORD_RE.match(line))


# --------------------------------------------------------------------------
# 4.9 check
# --------------------------------------------------------------------------

def _remote_date(meta: dict) -> str:
    value = meta.get("last_updated")
    if value in (None, ""):
        return NO_VERSION
    return str(value)


def _local_rows(dbdir: str, name: str) -> int | None:
    path = os.path.join(_scheme_dir(dbdir, name), name + PROFILE_SUFFIX)
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return None
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    return max(0, len(lines) - 1)


def check(dbdir: str, refs=None, *, progress=None, cancel=None, fetcher=None) -> UpdatePlan:
    """Metadata-only pass over every scheme (sections 4.9, 7.2). Writes nothing.

    162 GETs, ~230 KB, ~15 s. Each scheme is classified as ``up-to-date``,
    ``changed`` (a payload fetch is needed to know whether the content really
    moved), ``retired-upstream``, ``quarantined`` or ``error``. A scheme is never
    pruned for being absent from the host's ``/db`` group listing — only a 404
    ``"...has not been defined"`` justifies ``retired-upstream``, and even then
    the local directory is left untouched.
    """
    refs = tuple(refs) if refs is not None else load_manifest(dbdir)
    fetcher = fetcher or _Fetcher(cancel=cancel)
    run_id = "%s-%s" % (_utc_now().strftime("%Y%m%dT%H%M%SZ"),
                        hashlib.sha256(os.urandom(16)).hexdigest()[:8])
    checked_at = _now_iso()

    fetch_refs = [r for r in refs if r.resolved and r.alias_of is None]
    unresolved = [r for r in refs if not r.resolved]
    total = max(1, len(fetch_refs))
    _emit(progress, 0.0, "Checking %d schemes for updates" % len(refs))

    metas = {}
    errors = {}
    lock = threading.Lock()
    done = [0]

    def one(ref):
        try:
            meta = fetcher.get_json(ref.api)
        except HttpStatusError as exc:
            with lock:
                errors[ref.api] = exc
        except Cancelled:
            raise
        except UpdateError as exc:
            with lock:
                errors[ref.api] = exc
        else:
            with lock:
                metas[ref.api] = meta
        with lock:
            done[0] += 1
            position = done[0]
        _emit(progress, position / total, "Checked %d/%d schemes" % (position, total))

    hosts = {urllib.parse.urlsplit(r.api).netloc for r in fetch_refs}
    limits = {h: threading.Semaphore(MAX_WORKERS_PER_HOST) for h in hosts}

    def guarded(ref):
        with limits[urllib.parse.urlsplit(ref.api).netloc]:
            one(ref)

    if fetch_refs:
        workers = MAX_WORKERS_PER_HOST * max(1, len(hosts))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(guarded, r) for r in fetch_refs]
            try:
                for future in concurrent.futures.as_completed(futures):
                    future.result()
            except BaseException:
                for future in futures:
                    future.cancel()
                raise

    updates = []
    for ref in sorted(refs, key=lambda r: r.name.encode("utf-8")):
        info = scheme_info(dbdir, ref.name)
        local_date = str(info.get("last_updated") or "Unknown")
        if not ref.resolved:
            updates.append(SchemeUpdate(
                ref.name, local_date, "", "error", 0, None,
                "No upstream mapping is recorded for this scheme."))
            continue
        api = ref.api
        if api in errors:
            exc = errors[api]
            if isinstance(exc, HttpStatusError) and exc.is_undefined_scheme:
                updates.append(SchemeUpdate(
                    ref.name, local_date, "", "retired-upstream", 0, None,
                    "Upstream no longer defines this scheme; the installed copy "
                    "was kept."))
            else:
                updates.append(SchemeUpdate(
                    ref.name, local_date, "", "error", 0, None, str(exc)))
            continue
        meta = metas.get(api)
        if meta is None:
            updates.append(SchemeUpdate(
                ref.name, local_date, "", "error", 0, None,
                "No answer from %s." % api))
            continue
        remote = _remote_date(meta)
        locus_count = int(meta.get("locus_count") or 0)
        if ref.nloci and locus_count and locus_count != ref.nloci:
            updates.append(SchemeUpdate(
                ref.name, local_date, remote, "quarantined", 0, None,
                "Upstream now has %d loci, the installed scheme has %d; confirm "
                "before updating." % (locus_count, ref.nloci)))
            continue
        verified = info.get("wmlst_upstream_last_updated_at_verification")
        local_dir = _scheme_dir(dbdir, ref.name)
        if not os.path.isdir(local_dir):
            updates.append(SchemeUpdate(
                ref.name, local_date, remote, "changed", 0, None,
                "Not installed locally."))
            continue
        if verified and str(verified) == remote:
            updates.append(SchemeUpdate(
                ref.name, local_date, remote, "up-to-date", 0, 0,
                "Upstream has not changed since the last verification."))
            continue
        rows = _local_rows(dbdir, ref.name)
        records = meta.get("records")
        added = None
        if rows is not None and isinstance(records, int):
            added = records - rows
        updates.append(SchemeUpdate(
            ref.name, local_date, remote, "changed", _dir_size(local_dir), added,
            "Upstream reports %s; the content hash will be checked on download."
            % remote))

    total_bytes = sum(u.bytes_estimate for u in updates if u.status == "changed")
    _emit(progress, 1.0, "Checked %d schemes" % len(updates))
    if unresolved:
        log.warning("%d schemes have no upstream mapping in the manifest",
                    len(unresolved))
    return UpdatePlan(tuple(updates), (), total_bytes, checked_at, run_id)


# --------------------------------------------------------------------------
# 7.3 the atomic swap
# --------------------------------------------------------------------------

def _swap_dir(live: str, staged: str, backup: str | None) -> None:
    """Replace ``live`` with ``staged`` (section 7.3).

    ``os.replace`` cannot atomically replace a non-empty directory on Windows,
    so the live copy is renamed aside first and deleted only after the staging
    directory is in place. A crash between the two renames is repaired by
    :func:`sweep_staging`.
    """
    os.makedirs(os.path.dirname(os.path.abspath(live)), exist_ok=True)
    old = None
    if os.path.isdir(live):
        if backup:
            if os.path.isdir(backup):
                _rmtree(backup)
            os.makedirs(os.path.dirname(os.path.abspath(backup)), exist_ok=True)
            shutil.copytree(live, backup)
        old = "%s.old.%d" % (live, os.getpid())
        if os.path.isdir(old):
            _rmtree(old)
        os.replace(live, old)
    try:
        os.replace(staged, live)
    except OSError:
        if old and not os.path.isdir(live):
            os.replace(old, live)
        raise
    if old:
        _rmtree(old)


def _materialise(dbdir: str, ref: SchemeRef, meta: dict, profiles: bytes,
                 alleles, digest: str, *, backup: bool, write_version_files: bool,
                 allow_shrink: bool) -> None:
    """Stage, validate and atomically commit one scheme directory (7.3, 7.6)."""
    _validate_scheme_name(ref.name)
    root = _pubmlst_dir(dbdir)
    staging_root = os.path.join(root, STAGING_DIRNAME)
    staged = os.path.join(staging_root, "%s.%d" % (ref.name, os.getpid()))
    _rmtree(staged, required=True)
    os.makedirs(staged, exist_ok=True)
    live = os.path.join(root, ref.name)
    if not (_is_within(root, live) and _is_within(staging_root, staged)):
        raise UpdateError("Scheme %r would be written outside %s." % (ref.name, root))
    try:
        for locus in sorted(alleles, key=lambda s: s.encode("utf-8")):
            _validate_locus_name(locus)
            with open(os.path.join(staged, locus + ALLELE_SUFFIX), "wb") as fh:
                fh.write(alleles[locus])
        with open(os.path.join(staged, ref.name + PROFILE_SUFFIX), "wb") as fh:
            fh.write(profiles)
        remote = _remote_date(meta)
        info = _build_info(ref, meta, digest=digest, remote_date=remote,
                           previous=scheme_info(dbdir, ref.name),
                           unparseable=_unparseable_loci(sorted(alleles)))
        write_scheme_info(os.path.join(staged, ref.name + INFO_SUFFIX), info)
        if write_version_files:
            _write_bytes(os.path.join(staged, VERSION_STAMP_FILE),
                         (remote + "\n").encode("utf-8"))
        _validate_scheme(ref.name, staged, meta, sorted(alleles), live,
                         allow_shrink=allow_shrink)
        backup_dir = (os.path.join(dbdir, ROLLBACK_DIRNAME, ref.name)
                      if backup else None)
        _swap_dir(live, staged, backup_dir)
    finally:
        _rmtree(staged)


def _touch_info(dbdir: str, ref: SchemeRef, meta: dict, digest: str) -> None:
    """``date-moved-content-identical``: refresh the ledger, rewrite nothing else."""
    _validate_scheme_name(ref.name)
    info = dict(scheme_info(dbdir, ref.name))
    remote = _remote_date(meta)
    info.setdefault("name", ref.name)
    info.setdefault("description", meta.get("description", "MLST"))
    info.setdefault("locus", int(meta.get("locus_count") or ref.nloci))
    info.setdefault("source", ref.source)
    info.setdefault("API", ref.api)
    info.setdefault("authenticated", False)
    info["download_date"] = _today()
    info["last_updated"] = remote
    info["wmlst_content_sha256"] = digest
    info["wmlst_verified_unchanged_at"] = _now_iso()
    info["wmlst_upstream_last_updated_at_verification"] = remote
    write_scheme_info(
        os.path.join(_scheme_dir(dbdir, ref.name), ref.name + INFO_SUFFIX), info)


# --------------------------------------------------------------------------
# 4.9 apply
# --------------------------------------------------------------------------

def _chosen_names(plan: UpdatePlan, chosen) -> tuple:
    if chosen is None:
        picked = [u.name for u in plan.changed]
    else:
        picked = []
        for item in chosen:
            picked.append(item.name if isinstance(item, SchemeUpdate) else str(item))
    return tuple(sorted(set(picked), key=lambda s: s.encode("utf-8")))


def apply(dbdir: str, plan: UpdatePlan, chosen, *, backup: bool = True,
          write_version_files: bool = False, progress=None, cancel=None,
          refs=None, allow_shrink: bool = False, force_blast: bool = False,
          resume: bool = False, tools=None, fetcher=None) -> str:
    """Download, validate, commit per scheme, rebuild the index (sections 4.9, 7).

    ``chosen`` is an iterable of scheme names or :class:`SchemeUpdate` objects
    (``None`` means every ``changed`` scheme in ``plan``). Alias pairs
    (``cdiphtheriae`` / ``diphtheria_3``) are fetched once and materialised twice.
    A scheme whose freshly downloaded content hashes to the installed value is
    recorded as ``date-moved-content-identical``: only its ``_info.json`` ledger
    is touched, the payload files are not rewritten and the BLAST index is not
    invalidated. Cancellation is honoured BETWEEN scheme commits, so the tree is
    never left inconsistent. -> the new ``db/VERSION.txt`` value.
    """
    refs = tuple(refs) if refs is not None else load_manifest(dbdir)
    by_name = {r.name: r for r in refs}
    names = _chosen_names(plan, chosen)
    journal = _Journal(journal_path(dbdir), plan.run_id)
    sweep_staging(dbdir)
    already = committed_schemes(dbdir, plan.run_id) if resume else frozenset()

    unknown = [n for n in names if n not in by_name]
    if unknown:
        raise UpdateError("Unknown scheme(s): %s" % ", ".join(sorted(unknown)))

    # One fetch per upstream identity; aliases ride along (7.3).
    groups = {}
    for name in names:
        if name in already:
            journal.record("scheme_skip", scheme=name, reason="already committed")
            continue
        ref = by_name[name]
        if not ref.resolved:
            journal.record("scheme_error", scheme=name, reason="unresolved")
            raise UpdateError(
                "%s has no upstream mapping; fix db/%s first." % (name, MANIFEST_FILE))
        source = by_name[ref.alias_of] if ref.alias_of and ref.alias_of in by_name else ref
        groups.setdefault(source.triple, {"source": source, "members": []})
        groups[source.triple]["members"].append(ref)
    journal.record("run_start", schemes=list(names), backup=bool(backup),
                   version=__version__)
    fetcher = fetcher or _Fetcher(cancel=cancel)
    dirty = []
    unchanged = []
    total_groups = max(1, len(groups))
    index = 0
    try:
        for triple in sorted(groups, key=lambda t: groups[t]["source"].name.encode("utf-8")):
            _check_cancel(cancel)
            group = groups[triple]
            source = group["source"]
            members = sorted(group["members"], key=lambda r: r.name.encode("utf-8"))
            base = index / total_groups
            span = 1.0 / total_groups
            index += 1
            label = "/".join(m.name for m in members)
            _emit(progress, base, "Downloading %s" % label)
            journal.record("scheme_fetch_begin", scheme=source.name,
                           members=[m.name for m in members], api=source.api)
            try:
                meta = fetcher.get_json(source.api)
                locus_urls = list(meta.get("loci") or ())
                if not locus_urls:
                    raise UpdateError("%s: upstream returned no loci." % source.name)
                loci = {}
                for url in locus_urls:
                    check_fetch_url(str(url))
                    locus = _locus_of_url(url)
                    _validate_locus_name(locus)
                    loci[url + "/alleles_fasta"] = locus
                profiles_url = source.api + "/profiles_csv"
                responses = fetcher.gather(
                    [profiles_url, *list(loci)],
                    progress=progress, base=base, span=span * 0.9,
                    message="Downloading %s" % label)
                fields = _field_names(meta)
                locus_names = [loci[u] for u in loci]
                rows = validate_profiles(source.name, responses[profiles_url],
                                         locus_names, fields)
                alleles = {}
                for url, locus in loci.items():
                    validate_tfa(locus, responses[url])
                    alleles[locus] = responses[url].body
                profiles = responses[profiles_url].body
                digest = content_sha256(profiles, alleles)
                for member in members:
                    _check_cancel(cancel)
                    if local_content_sha256(dbdir, member.name) == digest:
                        _touch_info(dbdir, member, meta, digest)
                        unchanged.append(member.name)
                        journal.record("scheme_commit", scheme=member.name,
                                       status="date-moved-content-identical",
                                       sha256=digest, rows=rows)
                        continue
                    _materialise(dbdir, member, meta, profiles, alleles, digest,
                                 backup=backup,
                                 write_version_files=write_version_files,
                                 allow_shrink=allow_shrink)
                    dirty.append(member.name)
                    journal.record("scheme_commit", scheme=member.name,
                                   status="changed", sha256=digest, rows=rows)
            except Cancelled:
                raise
            except UpdateError as exc:
                journal.record("scheme_error", scheme=source.name, reason=str(exc))
                raise
            _emit(progress, base + span, "Updated %s" % label)
    except Cancelled:
        journal.record("run_end", status="cancelled", changed=dirty,
                       unchanged=unchanged)
        raise
    except BaseException as exc:
        journal.record("run_end", status="error", changed=dirty,
                       unchanged=unchanged, reason=str(exc))
        raise

    if backup:
        _write_rollback_meta(dbdir, plan.run_id, dirty)
    # Rebuild the index BEFORE stamping the version. Stamping first would leave a
    # failed rebuild advertising a fresh date over a stale or absent index, and
    # the next run would then see nothing to do.
    if dirty or force_blast or index_is_stale(dbdir):
        _emit(progress, 0.95, "Rebuilding the BLAST database")
        build_blast_db(dbdir, tools, progress=progress, cancel=cancel)
    version = stamp_db_version(dbdir, _today())
    journal.record("run_end", status="ok", changed=dirty, unchanged=unchanged,
                   db_version=version)
    _emit(progress, 1.0, "Database updated (%d changed, %d already current)"
          % (len(dirty), len(unchanged)))
    return version


# --------------------------------------------------------------------------
# 7.3 rollback
# --------------------------------------------------------------------------

def _rollback_dir(dbdir: str) -> str:
    return os.path.join(dbdir, ROLLBACK_DIRNAME)


def _write_rollback_meta(dbdir: str, run_id: str, changed) -> None:
    root = _rollback_dir(dbdir)
    if not os.path.isdir(root):
        return
    meta = {"run_id": run_id, "created_at": _now_iso(),
            "previous_db_version": _rollback_previous_version(dbdir),
            "schemes": sorted(changed)}
    _write_bytes(os.path.join(root, "rollback.json"),
                 (json.dumps(meta, indent=2, sort_keys=True) + "\n").encode("utf-8"))


def _rollback_previous_version(dbdir: str) -> str:
    root = _rollback_dir(dbdir)
    stash = os.path.join(root, VERSION_FILE)
    if os.path.isfile(stash):
        with open(stash, encoding="utf-8") as fh:
            return fh.readline().strip()
    current = db_version(dbdir)
    if current:
        _write_bytes(stash, (current + "\n").encode("utf-8"))
    return current


def rollback(dbdir: str) -> str:
    """Restore the copies kept by the last :func:`apply` (section 7.3).

    -> the ``db/VERSION.txt`` value after the restore. The BLAST index is left
    stale on purpose so :func:`index_is_stale` forces a rebuild.
    """
    root = _rollback_dir(dbdir)
    if not os.path.isdir(root):
        raise UpdateError("There is no rollback copy to restore.")
    restored = []
    for entry in sorted(os.listdir(root)):
        path = os.path.join(root, entry)
        if not os.path.isdir(path):
            continue
        live = _scheme_dir(dbdir, entry)
        staged = os.path.join(_pubmlst_dir(dbdir), STAGING_DIRNAME,
                              "%s.restore.%d" % (entry, os.getpid()))
        os.makedirs(os.path.dirname(staged), exist_ok=True)
        _rmtree(staged, required=True)
        shutil.copytree(path, staged)
        _swap_dir(live, staged, None)
        restored.append(entry)
    version = _rollback_previous_version(dbdir)
    if version:
        stamp_db_version(dbdir, version)
    meta_path = os.path.join(root, "rollback.json")
    meta = {}
    if os.path.isfile(meta_path):
        try:
            with open(meta_path, encoding="utf-8") as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = {}
    meta["restored_at"] = _now_iso()
    meta["restored"] = restored
    _rmtree(root)
    os.makedirs(root, exist_ok=True)
    _write_bytes(meta_path, (json.dumps(meta, indent=2, sort_keys=True) + "\n")
                 .encode("utf-8"))
    _Journal(journal_path(dbdir), meta.get("run_id", "rollback")).record(
        "rollback", schemes=restored, db_version=version)
    return version or db_version(dbdir)


# --------------------------------------------------------------------------
# 4.9 build_blast_db / index_is_stale
# --------------------------------------------------------------------------

def _blastbin():
    """Lazily import :mod:`wmlst.blastbin` — the only module allowed subprocess."""
    try:
        return importlib.import_module("wmlst.blastbin")
    except ImportError as exc:  # pragma: no cover - packaging failure
        raise UpdateError(
            "The BLAST+ helper module is unavailable (%s); the search index "
            "cannot be rebuilt." % exc) from exc


def concat_fasta(dbdir: str, out_path: str, *, progress=None, cancel=None) -> int:
    """Write the concatenated ``mlst.fa`` (port of ``scripts/mlst-make_blast_db``).

    Deterministic: schemes in byte order, ``*.tfa`` within a scheme in byte
    order (the shell version uses readdir + ``LC_COLLATE``, which is neither
    stable nor portable). Each header becomes ``b'>' + scheme + b'.' + line[1:]``
    on BYTES. Whole records whose header contains ``not a locus`` are dropped —
    upstream's ``grep -v`` orphans their sequence lines onto the previous record
    (divergence D4). -> the number of sequence records written.
    """
    schemes = _scheme_names_on_disk(dbdir)
    records = 0
    total = max(1, len(schemes))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "wb") as out:
        for position, scheme in enumerate(schemes, 1):
            _check_cancel(cancel)
            prefix = b">" + scheme.encode("utf-8") + b"."
            scheme_dir = os.path.join(_pubmlst_dir(dbdir), scheme)
            for filename in _tfa_names(scheme_dir):
                with open(os.path.join(scheme_dir, filename), "rb") as fh:
                    data = fh.read()
                keep = True
                for line in data.split(b"\n"):
                    if line == b"":
                        continue
                    if line.startswith(b">"):
                        keep = b"not a locus" not in line
                        if keep:
                            records += 1
                            out.write(prefix + line[1:] + b"\n")
                        continue
                    if keep:
                        out.write(line + b"\n")
            _emit(progress, position / total,
                  "Collecting alleles (%d/%d schemes)" % (position, total))
    return records


def build_blast_db(dbdir: str, tools=None, *, progress=None, cancel=None) -> str:
    """Rebuild ``db/blast/mlst.fa`` and its v5 index (section 4.9). -> the FASTA path.

    Builds into the sibling ``blast.staging.<pid>/``, runs ``makeblastdb
    -hash_index -in <fa> -dbtype nucl -title PubMLST -parse_seqids`` there,
    verifies all 12 ``.n*`` files, then installs the whole directory with ONE
    rename: 13 individual renames could be interrupted half-way and leave a new
    ``.nin`` indexing the old ``.nsq``, which types every sample as NONE without
    any error. The title is exactly ``PubMLST``; it is embedded in the index
    bytes and a vendor string there would break byte-comparison with the
    reference distribution.
    """
    blastbin = _blastbin()
    if tools is None:
        tools = blastbin.find_blast()
    blast_dir = os.path.join(dbdir, BLAST_DIRNAME)
    # A SIBLING of blast/, never a child: the install is then one directory
    # rename instead of 13 file renames, so an interruption can never leave the
    # new .nin indexing the old .nsq. The pid keeps two runs out of each
    # other's way, and sweep_staging clears anything a crash leaves behind.
    staging = os.path.join(dbdir, BLAST_STAGING_PREFIX + str(os.getpid()))
    _rmtree(staging, required=True)
    os.makedirs(staging, exist_ok=True)
    staged_fa = os.path.join(staging, BLAST_FASTA)
    _emit(progress, 0.0, "Collecting alleles")
    records = concat_fasta(dbdir, staged_fa, progress=progress, cancel=cancel)
    if records == 0:
        _rmtree(staging)
        raise UpdateError("No alleles were found in %s." % _pubmlst_dir(dbdir))
    _check_cancel(cancel)
    _emit(progress, 0.5, "Indexing %d alleles with makeblastdb" % records)
    # NOT run_tool(makeblastdb ... -in <abs path>) directly: makeblastdb splits
    # both the -in value and the absolute path it re-opens the finished database
    # under on whitespace, so it cannot write an index into a directory such as
    # C:\Program Files\... at all. run_makeblastdb() builds in a whitespace-free
    # scratch directory when needed and moves the (relocatable) index into place.
    completed = blastbin.run_makeblastdb(tools, staged_fa, cancel=cancel)
    if completed.returncode != 0:
        tail = completed.stderr
        if isinstance(tail, bytes):
            tail = tail.decode("utf-8", "replace")
        _rmtree(staging)
        raise UpdateError("makeblastdb failed (exit %d): %s"
                          % (completed.returncode, (tail or "").strip()[-500:]))
    missing = [ext for ext in BLAST_INDEX_EXTENSIONS
               if not os.path.isfile(staged_fa + "." + ext)]
    if missing:
        _rmtree(staging)
        raise UpdateError("makeblastdb produced an incomplete index (missing %s)."
                          % ", ".join(missing))
    _emit(progress, 0.95, "Installing the new search index")
    final_fa = os.path.join(blast_dir, BLAST_FASTA)
    try:
        _swap_dir(blast_dir, staging, None)
    finally:
        _rmtree(staging)
    _emit(progress, 1.0, "Search index ready (%d alleles)" % records)
    return final_fa


def index_is_stale(dbdir: str) -> bool:
    """True when the BLAST index is missing or older than the newest ``.tfa`` (4.9).

    The engine MUST refuse to type with an actionable message when this is True.
    """
    nsq = os.path.join(dbdir, BLAST_DIRNAME, BLAST_FASTA + ".nsq")
    try:
        index_mtime = os.path.getmtime(nsq)
    except OSError:
        return True
    try:
        schemes = _scheme_names_on_disk(dbdir)
    except DatabaseMissingError:
        return True
    for scheme in schemes:
        scheme_dir = os.path.join(_pubmlst_dir(dbdir), scheme)
        try:
            entries = _tfa_names(scheme_dir)
        except OSError:
            continue
        for filename in entries:
            try:
                if os.path.getmtime(os.path.join(scheme_dir, filename)) > index_mtime:
                    return True
            except OSError:
                continue
    return False


# --------------------------------------------------------------------------
# 4.9 export_bundle / import_bundle — air-gapped machines
# --------------------------------------------------------------------------

_BUNDLE_MANIFEST = "wmlst-bundle.json"


def export_bundle(dbdir: str, path: str) -> None:
    """Pack the allele database into one deflated zip for an offline machine (4.9).

    The derived BLAST index is never included: it is rebuilt on import.
    """
    schemes = _scheme_names_on_disk(dbdir)
    meta = {"format": "wmlst-db-bundle/1", "created_at": _now_iso(),
            "wmlst_version": __version__, "db_version": db_version(dbdir),
            "schemes": list(schemes)}
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(_BUNDLE_MANIFEST,
                    json.dumps(meta, indent=2, sort_keys=True) + "\n")
        for extra in (VERSION_FILE, MANIFEST_FILE, SPECIES_MAP_FILE):
            source = os.path.join(dbdir, extra)
            if os.path.isfile(source):
                zf.write(source, extra)
        for scheme in schemes:
            scheme_dir = os.path.join(_pubmlst_dir(dbdir), scheme)
            for filename in sorted(os.listdir(scheme_dir),
                                   key=lambda s: s.encode("utf-8")):
                full = os.path.join(scheme_dir, filename)
                if os.path.isfile(full):
                    zf.write(full, "%s/%s/%s" % (PUBMLST_DIRNAME, scheme, filename))


def import_bundle(dbdir: str, path: str) -> str:
    """Install a bundle written by :func:`export_bundle` (section 4.9).

    Members are extracted into a staging tree and committed per scheme with the
    same atomic swap as an online update. -> the installed ``db/VERSION.txt``.
    """
    if not zipfile.is_zipfile(path):
        raise UpdateError("%s is not a WMLST database bundle." % path)
    root = _pubmlst_dir(dbdir)
    os.makedirs(root, exist_ok=True)
    staging_root = os.path.join(root, STAGING_DIRNAME)
    os.makedirs(staging_root, exist_ok=True)
    work = os.path.join(staging_root, "bundle.%d" % os.getpid())
    _rmtree(work, required=True)
    os.makedirs(work)
    try:
        with zipfile.ZipFile(path) as zf:
            try:
                meta = json.loads(zf.read(_BUNDLE_MANIFEST).decode("utf-8"))
            except KeyError as exc:
                raise UpdateError(
                    "%s is not a WMLST database bundle." % path) from exc
            for member in zf.infolist():
                name = member.filename
                if member.is_dir():
                    continue
                # Normalise FIRST, then judge. Splitting the raw name on "/"
                # leaves a Windows member such as C:\\Windows\\Temp\\evil.txt as a
                # single component, and ntpath.join then returns it unchanged --
                # an arbitrary absolute write outside the staging tree.
                norm = name.replace("\\", "/")
                if (norm.startswith("/") or ".." in norm.split("/")
                        or _DRIVE_RE.match(norm)):
                    raise UpdateError("%s contains an unsafe path (%s)." % (path, name))
                target = os.path.join(work, *norm.split("/"))
                if not _is_within(work, target):
                    raise UpdateError("%s contains an unsafe path (%s)." % (path, name))
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with zf.open(member) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        staged_pubmlst = os.path.join(work, PUBMLST_DIRNAME)
        if not os.path.isdir(staged_pubmlst):
            raise UpdateError("%s carries no scheme directories." % path)
        for scheme in sorted(os.listdir(staged_pubmlst),
                             key=lambda s: s.encode("utf-8")):
            source = os.path.join(staged_pubmlst, scheme)
            if not os.path.isdir(source):
                continue
            _validate_scheme_name(scheme, path)
            _swap_dir(os.path.join(root, scheme), source, None)
        for extra in (VERSION_FILE, MANIFEST_FILE, SPECIES_MAP_FILE):
            staged_extra = os.path.join(work, extra)
            if os.path.isfile(staged_extra):
                with open(staged_extra, "rb") as fh:
                    _write_bytes(os.path.join(dbdir, extra), fh.read())
    finally:
        _rmtree(work)
    version = db_version(dbdir) or str(meta.get("db_version") or _today())
    _Journal(journal_path(dbdir), "import").record(
        "run_end", status="imported", db_version=version, bundle=os.path.basename(path))
    return version


# --------------------------------------------------------------------------
# 7.5 discovery — opt-in, default OFF
# --------------------------------------------------------------------------

def discover_new_schemes(dbdir: str, *, refs=None, progress=None, cancel=None,
                         fetcher=None) -> tuple:
    """Find MLST schemes upstream that are not installed (section 7.5).

    OPT-IN and default OFF everywhere. Keeps only schemes whose description
    matches ``\\bMLST\\b`` and whose ``locus_count`` is in [2, 15]; rejects
    cgMLST / rMLST / wgMLST / ``Core N`` / ``TCS-N``. Adding a scheme that is not
    in tseemann's distribution CHANGES typing results, because a new scheme
    competes in the auto-detect scoring — the caller must say so and must
    require explicit confirmation. -> ``tuple[SchemeRef, ...]`` of proposals.
    """
    refs = tuple(refs) if refs is not None else load_manifest(dbdir)
    known = {r.triple for r in refs if r.resolved}
    taken = {r.name for r in refs} | set(_scheme_names_on_disk(dbdir))
    fetcher = fetcher or _Fetcher(cancel=cancel)
    databases = []
    for source, root in sorted(REST_ROOTS.items()):
        _emit(progress, 0.0, "Listing databases on %s" % root)
        try:
            listing = fetcher.get_json(root + "/db")
        except UpdateError:
            log.warning("could not list databases on %s", root)
            continue
        for group in listing or ():
            for database in group.get("databases", ()) or ():
                name = database.get("name") or ""
                if name.endswith("_seqdef"):
                    databases.append((source, root, name))
    # Never prune for absence from /db: some bundled databases are unlisted yet
    # answer HTTP 200 when addressed directly (6.3).
    for ref in refs:
        if ref.resolved and ref.db in UNLISTED_DATABASES:
            entry = (ref.source, REST_ROOTS[ref.source], ref.db)
            if entry not in databases:
                databases.append(entry)
    proposals = []
    total = max(1, len(databases))
    for position, (source, root, database) in enumerate(sorted(set(databases)), 1):
        _check_cancel(cancel)
        _emit(progress, position / total, "Scanning %s" % database)
        try:
            schemes = fetcher.get_json("%s/db/%s/schemes" % (root, database))
        except UpdateError:
            continue
        for scheme in (schemes or {}).get("schemes", ()) or ():
            description = str(scheme.get("description") or "")
            if not _MLST_DESC_RE.search(description) or _NOT_MLST_RE.search(description):
                continue
            match = _API_RE.match(str(scheme.get("scheme") or ""))
            if not match:
                continue
            scheme_id = match.group(2)
            if (source, database, scheme_id) in known:
                continue
            count = scheme.get("locus_count")
            if not isinstance(count, int) or not (
                    DISCOVERY_LOCUS_RANGE[0] <= count <= DISCOVERY_LOCUS_RANGE[1]):
                continue
            proposed = _derive_dirname(database, scheme_id)
            try:
                _validate_scheme_name(proposed)
            except UpdateError as exc:
                log.warning("skipping upstream scheme %s/%s: %s",
                            database, scheme_id, exc)
                continue
            candidate, suffix = proposed, 1
            while candidate in taken:
                suffix += 1
                candidate = "%s_x%d" % (proposed, suffix)
            taken.add(candidate)
            proposals.append(SchemeRef(candidate, source, database, scheme_id,
                                       count, None))
    return tuple(sorted(proposals, key=lambda r: r.name.encode("utf-8")))


# --------------------------------------------------------------------------
# Maintenance entry point (not the user-facing CLI, which lives in cli.py)
# --------------------------------------------------------------------------

_USAGE = """\
wmlst.updatedb - maintenance helper for the WMLST allele database

  python -m wmlst.updatedb <command> [dbdir]

Commands:
  check [dbdir]            contact PubMLST/Pasteur and report which schemes changed
  manifest [dbdir] [--verify]
                           regenerate db/schemes.manifest.tsv from the live APIs
  blastdb [dbdir]          rebuild the BLAST database from the bundled alleles
  export [dbdir] [out.zip] package the database as a distributable zip

dbdir defaults to the db/ folder shipped alongside this package.

This is a maintenance entry point. The supported user interface for updating
the database is `wmlst --update-db`, or the Database tab in the WMLST GUI.
"""


def _main(argv=None) -> int:
    """``python -m wmlst.updatedb <manifest|check|blastdb|export> <dbdir>``.

    A maintenance helper for regenerating ``db/schemes.manifest.tsv`` and for
    smoke-testing the updater; the supported user interface is ``wmlst --update-db``.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        sys.stderr.write(_USAGE)
        return 2
    command, rest = argv[0], argv[1:]
    if command in ("-h", "--help", "help"):
        sys.stdout.write(_USAGE)
        return 0
    dbdir = rest[0] if rest else os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "db")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    def show(fraction, message):
        sys.stderr.write("\r[%3d%%] %-60s" % (round(fraction * 100), message))
        sys.stderr.flush()

    if command == "manifest":
        refs = build_manifest(dbdir, verify="--verify" in rest, progress=show)
        path = write_manifest(dbdir, refs)
        sys.stderr.write("\n")
        print("wrote %s (%d schemes, %d unresolved)"
              % (path, len(refs), sum(1 for r in refs if not r.resolved)))
        return 0
    if command == "check":
        plan = check(dbdir, progress=show)
        sys.stderr.write("\n")
        for update in plan.updates:
            print("%-24s %-28s %s" % (update.name, update.status, update.detail))
        print("%d changed, %d bytes estimated"
              % (len(plan.changed), plan.total_bytes))
        return 0
    if command == "blastdb":
        path = build_blast_db(dbdir, progress=show)
        sys.stderr.write("\n")
        print(path)
        return 0
    if command == "export":
        target = rest[1] if len(rest) > 1 else "wmlst-db.zip"
        export_bundle(dbdir, target)
        print(target)
        return 0
    sys.stderr.write("unknown command %r\n\n" % command)
    sys.stderr.write(_USAGE)
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
