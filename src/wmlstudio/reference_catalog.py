"""Public, GET-only PubMLST reference discovery and immutable local snapshots.

Official protocol: https://bigsdb.readthedocs.io/en/latest/rest.html (checked
2026-09-12). No sample data, credentials, or submission requests are transmitted.
Unauthenticated PubMLST access can be restricted to an older public dataset; its
access message is retained and must not be presented as a fully current release.
"""

from __future__ import annotations

import csv
import hashlib
import http.client
import io
import json
import re
import shutil
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import UTC, datetime
from html.parser import HTMLParser
from pathlib import Path

from . import __version__, cgmlst_schemes
from .sequence import DNA, SequenceError, check_cancelled
from .typing import SchemeError, load_scheme

API_ROOT = "https://rest.pubmlst.org"
OFFICIAL_HOSTS = {"rest.pubmlst.org", "bigsdb.pasteur.fr", "www.cgmlst.org", "cgmlst.org"}
JSON_LIMIT = 16 * 1024 * 1024
FILE_LIMIT = 1024 * 1024 * 1024
SNAPSHOT_LIMIT = 8 * 1024 * 1024 * 1024
# A public service asked to rate-limit us can name a very long wait; honour it up
# to this bound and then say so, rather than hammering or hanging indefinitely.
MAX_RETRY_AFTER = 120
RESUME_LEDGER = "resume.json"
RESUME_PREFIX = ".resume-"
# Partial downloads live BESIDE the scheme library, never inside it: paths.scheme_locations
# hands every directory under the library to the typing code as a scheme, so a
# half-finished folder kept there would be offered as a reference.
RESUME_DIRNAME = ".wmlstudio-partial-downloads"


class CatalogError(ValueError):
    """A reference could not be fetched or validated without an unsafe assumption."""


def _validated_scheme(path, *, cancelled=None):
    """load_scheme, with its errors reported as catalog errors the caller expects.

    A scheme failure during a download is a download failure: callers documented to
    catch CatalogError must not have to catch SchemeError from a lower layer too.
    """
    try:
        return load_scheme(path, cancelled=cancelled)
    except (SchemeError, SequenceError) as error:
        raise CatalogError(
            f"The downloaded reference could not be read as a scheme: {error} "
            "Nothing was installed.") from error


def screen_allele_file(path, locus, *, cancelled=None) -> list[str]:
    """Remove allele records this application cannot read, naming every one removed.

    Public reference services publish occasional allele sequences containing a
    literal 'X' where a base is unknown. 'X' is not an IUPAC nucleotide code, so a
    single such record used to fail an entire multi-gigabyte scheme download at the
    final validation step. The record is excluded instead. It is never rewritten to
    'N': that would invent a base that was not reported. Each exclusion is returned,
    recorded in the snapshot manifest and shown to the user, because a genome
    carrying an excluded allele is reported as an unmatched sequence and not as that
    allele number.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="ascii")
    except (OSError, UnicodeDecodeError) as error:
        raise CatalogError(f"Allele data for locus {locus} is not readable ASCII FASTA.") from error
    records, header, body = [], None, []
    for line in text.splitlines():
        check_cancelled(cancelled)
        if line.startswith(">"):
            if header is not None:
                records.append((header, body))
            header, body = line, []
        elif header is not None:
            body.append(line)
    if header is not None:
        records.append((header, body))
    if not records:
        raise CatalogError(f"The reference service returned no allele records for locus {locus}.")
    kept, excluded = [], []
    for header, body in records:
        sequence = "".join(body).strip().upper()
        invalid = sorted(set(sequence) - DNA)
        if invalid or not sequence:
            reason = "empty sequence" if not sequence else "non-IUPAC " + ", ".join(invalid[:5])
            excluded.append(f"{locus}:{header[1:].split()[0] if header[1:].split() else '?'} "
                            f"({reason})")
        else:
            kept.append((header, body))
    if not excluded:
        return []
    if not kept:
        raise CatalogError(
            f"Every allele record for locus {locus} carries characters this application cannot "
            f"read ({excluded[0]}); no partial scheme was installed.")
    lines = []
    for header, body in kept:
        lines.append(header)
        lines.extend(body)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return excluded


def _safe_url(url: str) -> str:
    if not isinstance(url, str):
        raise CatalogError("PubMLST returned a non-text resource URL.")
    parts = urllib.parse.urlsplit(url)
    if (parts.scheme != "https" or parts.hostname not in OFFICIAL_HOSTS
            or parts.username or parts.password or parts.fragment or parts.port not in (None, 443)
            or any(part == ".." for part in urllib.parse.unquote(parts.path).split("/"))):
        raise CatalogError("Reference URLs must use an official HTTPS reference service.")
    return url


def _safe_name(name: str) -> str:
    reserved = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                *(f"LPT{i}" for i in range(1, 10))}
    if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,119}", name)
            or name.endswith(".") or name.split(".", 1)[0].upper() in reserved):
        raise CatalogError(f"Reference name {name!r} cannot be stored safely on Windows.")
    return name


def _file_digest(path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def install_destination(library_root, entry, digest="") -> Path:
    """Where a validated snapshot of this scheme is published, by kind.

    A classical seven-locus scheme keeps its content-addressed folder in the
    library the caller named: existing projects store those exact paths. A
    gene-by-gene scheme goes to the ONE cgMLST library instead, into a folder named
    for its organism, provider and target count -- because a 2,358-target scheme
    filed alphabetically among 162 seven-locus schemes under a digest-shaped name
    is a scheme the person who downloaded it cannot find.
    """
    library_root = Path(library_root)
    if cgmlst_schemes.is_gene_by_gene(entry):
        base = cgmlst_schemes.install_root(library_root)
        return cgmlst_schemes.install_folder(base, entry, digest=digest)
    base = _safe_name(f"{entry['database']}_{entry['scheme_id']}") if entry.get("database") \
        else _safe_name(str(entry.get("scheme_id") or "scheme"))
    return library_root / f"{base}_{digest[:16]}"


def _exclusion_notes(excluded) -> list[str]:
    """One plain sentence naming what the snapshot does not contain, and why."""
    if not excluded:
        return []
    shown = ", ".join(excluded[:5]) + (", …" if len(excluded) > 5 else "")
    return [f"{len(excluded)} allele record(s) were excluded because the provider published "
            f"characters that are not IUPAC nucleotide codes: {shown}. They were not rewritten to "
            "'N'. A genome carrying one of these alleles is reported as an unmatched sequence, "
            "never as that allele number."]


def organism_parts(label: str) -> dict:
    """Only parse explicit binomials; organism groups/complexes remain genus-only."""
    words = str(label).strip().split()
    if not words or words[0] in {"Plasmid", "Ribosomal", "Oral", "Sequence", "Candidatus",
                                'Unknown', 'Unidentified', 'Mixed', 'Metagenome'}:
        return {"genus": "", "species": ""}
    genus = words[0] if re.fullmatch(r"[A-Z][a-z]+", words[0]) else ""
    species = words[1] if (genus and len(words) == 2
                           and re.fullmatch(r"[a-z][a-z-]+", words[1])
                           and words[1] not in {"spp", "sp"}) else ""
    return {"genus": genus, "species": species}


class _OfficialRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, newurl):
        _safe_url(newurl)
        return super().redirect_request(request, fp, code, message, headers, newurl)


class PubMLSTCatalog:
    """A cancellable public reference client; instantiate once per background task.

    ``opener`` is injectable for deterministic network-failure tests. Its ``open``
    method follows urllib's context-manager response interface. Timeouts bound a
    blocked socket read; cancellation is checked between every received chunk.
    """

    service_name = "PubMLST"
    RESUMABLE = True

    def __init__(self, cache_dir=None, *, opener=None, timeout=20, retries=2,
                 base_url=API_ROOT):
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.opener = opener or urllib.request.build_opener(_OfficialRedirect())
        self.timeout = timeout
        self.retries = retries
        self.base_url = _safe_url(base_url).rstrip("/")
        self._cache: dict[str, tuple[float, object]] = {}
        self.catalog_from_cache = False

    @staticmethod
    def _retry_after(headers, attempt) -> float:
        """The wait a rate limiter asked for, bounded; otherwise our own backoff.

        Ignoring Retry-After is why a burst of scheme-metadata reads keeps failing:
        the service asks for seconds and the old backoff waited half of one.
        """
        value = (headers or {}).get("Retry-After") if headers is not None else None
        try:
            asked = float(str(value).strip())
        except (TypeError, ValueError):
            asked = 0.0
        return max(0.5 * (attempt + 1), min(asked, MAX_RETRY_AFTER))

    def _wait(self, seconds, cancelled):
        """Interruptible sleep: cancellation is checked ten times a second."""
        for _ in range(max(1, int(round(seconds * 10)))):
            check_cancelled(cancelled)
            threading.Event().wait(0.1)

    def _transfer(self, url, destination=None, *, limit=JSON_LIMIT, cancelled=None,
                  progress=None, label=""):
        url = _safe_url(url)
        for attempt in range(self.retries + 1):
            check_cancelled(cancelled)
            try:
                request = urllib.request.Request(url, headers={
                    "User-Agent": f"WMLSTudio/{__version__} (public reference download)",
                    "Accept": "application/json" if destination is None else "*/*",
                }, method="GET")
                with self.opener.open(request, timeout=self.timeout) as response:
                    _safe_url(response.geturl())
                    length = response.headers.get("Content-Length")
                    if length and int(length) > limit:
                        raise CatalogError("Reference response exceeds the configured size limit.")
                    digest = hashlib.sha256()
                    chunks = []
                    size = 0
                    reported = 0
                    output = Path(destination).open("wb") if destination is not None else None
                    try:
                        while chunk := response.read(65536):
                            check_cancelled(cancelled)
                            size += len(chunk)
                            if size > limit:
                                raise CatalogError("Reference response exceeds the configured size limit.")
                            digest.update(chunk)
                            if output is None:
                                chunks.append(chunk)
                            else:
                                output.write(chunk)
                            # Large bundles arrive chunk-encoded with no declared
                            # length, so without this the progress bar sits still
                            # for minutes and the download looks hung.
                            if progress and size - reported >= 4 * 1024 * 1024:
                                reported = size
                                progress(0, 0, f"{label or 'Downloading'}: {size // (1024 * 1024)} MB received")
                        check_cancelled(cancelled)
                        if length is not None and size != int(length):
                            raise CatalogError("Reference download was truncated; no snapshot was installed.")
                    finally:
                        if output is not None:
                            output.close()
                if destination is None:
                    return b"".join(chunks)
                return {"url": url, "file": Path(destination).name,
                        "sha256": digest.hexdigest(), "bytes": size}
            except urllib.error.HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504} or attempt == self.retries:
                    if error.code == 429:
                        raise CatalogError(
                            "The reference service is rate-limiting this computer (HTTP 429). Wait "
                            "a few minutes and try again, or download one scheme at a time.") from error
                    access = " Public access may require authentication." if error.code in {401, 403} else ""
                    raise CatalogError(f"{self.service_name} returned HTTP {error.code}.{access}") from error
                self._wait(self._retry_after(error.headers, attempt), cancelled)
                continue
            except (urllib.error.URLError, TimeoutError, ConnectionError,
                    http.client.IncompleteRead, http.client.HTTPException) as error:
                if attempt == self.retries:
                    raise CatalogError(f"Could not reach {self.service_name}: {error}") from error
            # Short interruptible backoff; GET-only requests are safe to retry.
            self._wait(0.5 * (attempt + 1), cancelled)
        raise CatalogError(f"{self.service_name} download did not complete.")

    def _json(self, url, *, cancelled=None, refresh=False):
        _safe_url(url)
        cached = self._cache.get(url)
        if not refresh and cached and time.monotonic() - cached[0] < 900:
            check_cancelled(cancelled)
            self.catalog_from_cache = True
            return json.loads(json.dumps(cached[1]))
        self.catalog_from_cache = False
        raw = self._transfer(url, cancelled=cancelled)
        try:
            body = json.loads(raw.decode("utf-8-sig"))
        except (UnicodeError, ValueError) as error:
            raise CatalogError("PubMLST returned invalid JSON reference metadata.") from error
        self._cache[url] = (time.monotonic(), body)
        return json.loads(json.dumps(body))

    def list_organisms(self, *, cancelled=None, refresh=False) -> list[dict]:
        groups = self._json(f"{self.base_url}/db", cancelled=cancelled, refresh=refresh)
        if not isinstance(groups, list):
            raise CatalogError("PubMLST organism catalog has an unsupported response shape.")
        organisms = []
        for group in groups:
            check_cancelled(cancelled)
            if not isinstance(group, dict) or not isinstance(group.get("databases"), list):
                raise CatalogError("PubMLST organism catalog contains an invalid database group.")
            for database in group["databases"]:
                if not isinstance(database, dict):
                    raise CatalogError("PubMLST returned an invalid database record.")
                identifier = database.get("name", "")
                if not identifier.endswith("_seqdef") or identifier == "pubmlst_test_seqdef":
                    continue
                label = str(database.get("description") or group.get("description") or identifier)
                label = label.removesuffix(" sequence/profile definitions")
                if "REST API" in label:
                    label = str(group.get("description") or group.get("name", "")).strip().removesuffix(" REST API group").removesuffix(" REST API group.")
                organisms.append({"id": identifier, "name": label, "organism": label,
                                  "url": _safe_url(database["href"]), **organism_parts(label)})
        return sorted(organisms, key=lambda item: item["name"].casefold())

    @staticmethod
    def _entry(details: dict, url: str, organism: str) -> dict:
        if not isinstance(details, dict) or not isinstance(details.get("loci"), list):
            raise CatalogError("PubMLST scheme metadata does not list its loci.")
        count = details.get("locus_count", len(details["loci"]))
        if isinstance(count, bool) or not isinstance(count, int) or count != len(details["loci"]):
            raise CatalogError("PubMLST scheme locus count does not match the listed loci.")
        name = str(details.get("description") or "Unnamed scheme")
        lower = name.casefold().replace(" ", "")
        kind = ("cgMLST" if "cgmlst" in lower or "coregenome" in lower else
                "wgMLST" if "wgmlst" in lower or "wholegenome" in lower else
                "MLST" if "mlst" in lower else "Other")
        parts = urllib.parse.urlsplit(_safe_url(url)).path.strip("/").split("/")
        if parts and parts[0] == "api":
            parts = parts[1:]
        if len(parts) != 4 or parts[0] != "db" or parts[2] != "schemes":
            raise CatalogError("PubMLST scheme URL has an unexpected path.")
        provider = "BIGSdb-Pasteur" if "pasteur.fr" in url else "PubMLST"
        return {"id": f"{provider}:{parts[1]}:{parts[3]}", "name": name, "organism": organism,
                "locus_count": count, "type": kind, "url": url,
                "provider": provider,
                "database": parts[1], "scheme_id": parts[3],
                "last_updated": details.get("last_updated"),
                "profile_count": details.get("records", 0),
                "has_profiles": bool(details.get("profiles_csv")),
                "access_notice": str(details.get("message") or ""),
                **organism_parts(organism)}

    def list_schemes(self, organism, *, cancelled=None, progress=None, refresh=False,
                     scheme_type=None, min_loci=0, max_loci=None) -> list[dict]:
        if isinstance(organism, str):
            choices = self.list_organisms(cancelled=cancelled)
            matches = [item for item in choices if organism in {item["id"], item["name"], item["url"]}]
            if len(matches) != 1:
                raise CatalogError("Choose one organism from the PubMLST catalog.")
            organism = matches[0]
        if not isinstance(organism, dict) or not organism.get("url"):
            raise CatalogError("Choose one organism from the PubMLST catalog.")
        listing = self._json(organism["url"] + "/schemes", cancelled=cancelled, refresh=refresh)
        if not isinstance(listing, dict) or not isinstance(listing.get("schemes"), list):
            raise CatalogError("PubMLST returned an invalid scheme list.")
        found = []
        for index, row in enumerate(listing["schemes"], 1):
            check_cancelled(cancelled)
            if not isinstance(row, dict) or "scheme" not in row:
                raise CatalogError("PubMLST returned an invalid scheme entry.")
            details = self._json(row["scheme"], cancelled=cancelled, refresh=refresh)
            entry = self._entry(details, row["scheme"], organism["organism"])
            if (entry["locus_count"] >= min_loci
                    and (max_loci is None or entry["locus_count"] <= max_loci)
                    and (not scheme_type or entry["type"].casefold() == scheme_type.casefold())):
                found.append(entry)
            if progress:
                progress(index, len(listing["schemes"]), f"Read scheme: {entry['name']}")
        return sorted(found, key=lambda item: (item["type"], item["locus_count"], item["name"]))

    def search_schemes(self, query, *, cancelled=None, progress=None, **filters) -> dict:
        query = str(query).strip().casefold()
        if not query:
            raise CatalogError("Enter an organism name before searching online schemes.")
        organisms = [item for item in self.list_organisms(cancelled=cancelled)
                     if query in item["name"].casefold() or query in item["id"].casefold()]
        schemes, errors = [], []
        for index, organism in enumerate(organisms, 1):
            check_cancelled(cancelled)
            try:
                schemes.extend(self.list_schemes(organism, cancelled=cancelled, **filters))
            except CatalogError as error:
                errors.append({"organism": organism["name"], "error": str(error)})
            if progress:
                progress(index, len(organisms), f"Searched {organism['name']}")
        return {"schemes": schemes, "errors": errors, "organisms_searched": len(organisms)}

    def download_estimate(self, entry) -> dict:
        """What a download of this scheme actually costs, before a byte is fetched.

        The service supplies no whole-scheme archive, so every target is a separate
        request. For a cgMLST scheme that is thousands of requests and gigabytes:
        the user must be told before starting, not after twenty minutes.
        """
        count = int(entry.get("locus_count") or 0)
        return {"requests": count + (1 if entry.get("has_profiles") else 0),
                "per_target_requests": True, "resumable": True,
                "notice": (f"This scheme has {count} targets and the service offers no whole-scheme "
                           f"archive, so {count} separate downloads are needed. A scheme of this "
                           "size is typically gigabytes on disk and takes many minutes. The "
                           "download can be cancelled and resumed where it stopped."
                           if count > 200 else
                           f"{count} targets are downloaded individually; this is quick.")}

    def _resume_root(self, library_root, current) -> Path:
        base = _safe_name(f"{current['database']}_{current['scheme_id']}")
        library_root = Path(library_root).resolve()
        return library_root.parent / RESUME_DIRNAME / (RESUME_PREFIX + base)

    def resume_state(self, library_root, entry) -> dict | None:
        """What a kept partial download contains, so a person can resume or discard it.

        Returns None when there is nothing partial. The folder is never a usable
        scheme; the count is what would not need downloading again.
        """
        try:
            folder = self._resume_root(library_root, dict(entry))
        except (KeyError, TypeError, CatalogError):
            return None
        ledger = folder / RESUME_LEDGER
        if not ledger.is_file():
            return None
        try:
            stored = json.loads(ledger.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, ValueError):
            return None
        files = stored.get("files") if isinstance(stored, dict) else None
        if not isinstance(files, dict):
            return None
        return {"path": str(folder), "files_held": len(files),
                "expected": int(entry.get("locus_count") or 0),
                "bytes": sum(int(item.get("bytes") or 0) for item in files.values()
                             if isinstance(item, dict)),
                "url": stored.get("url", "")}

    def clear_resume(self, library_root, entry) -> bool:
        """Discard a kept partial download. Installed snapshots are never touched."""
        try:
            folder = self._resume_root(library_root, dict(entry))
        except (KeyError, TypeError, CatalogError):
            return False
        if folder.is_dir() and folder.name.startswith(RESUME_PREFIX):
            shutil.rmtree(folder)
            return True
        return False

    def download_scheme(self, entry, library_root, *, cancelled=None, progress=None,
                        max_bytes=SNAPSHOT_LIMIT, resume=True) -> dict:
        """Install a fully validated snapshot; failures never publish a partial directory.

        A changed remote version is a new directory. Existing snapshots are never
        overwritten, and an existing destination is revalidated before reuse.

        With ``resume`` (the default), target files already fetched are kept in a
        clearly-named partial folder beside the library and reused on the next
        attempt, verified by the SHA-256 recorded when they were received. A
        thousand-target scheme that fails on target 900 then costs one more target,
        not another twenty minutes. The partial folder is never a scheme and is
        never published: publication still requires a complete, validated set. It is
        discarded automatically when the remote scheme has changed underneath it.
        """
        if not isinstance(entry, dict) or not entry.get("url"):
            raise CatalogError("Choose a PubMLST scheme before downloading.")
        check_cancelled(cancelled)
        details = self._json(entry["url"], cancelled=cancelled, refresh=True)
        current = self._entry(details, entry["url"], str(entry.get("organism") or "Unknown organism"))
        if not details["loci"]:
            raise CatalogError("This scheme does not contain any loci.")
        library_root = Path(library_root).resolve()
        library_root.mkdir(parents=True, exist_ok=True)
        resume_root = self._resume_root(library_root, current) if resume else None
        ledger = self._load_ledger(resume_root, details, current)
        if resume_root is not None:
            resume_root.mkdir(parents=True, exist_ok=True)
            staging = resume_root
        else:
            staging = Path(tempfile.mkdtemp(prefix=".pubmlst-download-", dir=library_root))
        keep_staging = resume_root is not None
        sources, names = [], set()
        transferred = 0
        excluded_alleles = []
        try:
            for index, locus_url in enumerate(details["loci"], 1):
                check_cancelled(cancelled)
                _safe_url(locus_url)
                name = _safe_name(urllib.parse.unquote(urllib.parse.urlsplit(locus_url).path.rsplit("/", 1)[-1]))
                if name.casefold() in names:
                    raise CatalogError("The scheme repeats a locus filename, including Windows case folding.")
                names.add(name.casefold())
                target = staging / f"{name}.tfa"
                source = self._reuse(ledger, target, locus_url + "/alleles_fasta")
                if source is None:
                    source = self._transfer(locus_url + "/alleles_fasta", target,
                                            limit=min(FILE_LIMIT, max_bytes - transferred),
                                            cancelled=cancelled, progress=progress,
                                            label=f"Target {index} of {len(details['loci'])}")
                    dropped = screen_allele_file(target, name, cancelled=cancelled)
                    if dropped:
                        excluded_alleles.extend(dropped)
                        source = {**source, "excluded_alleles": dropped,
                                  "stored_sha256": _file_digest(target)}
                    self._record(ledger, resume_root, source)
                else:
                    excluded_alleles.extend(source.get("excluded_alleles", []))
                transferred += source["bytes"]
                sources.append(source)
                if progress:
                    progress(index, len(details["loci"]) + 1, f"Downloaded {name} alleles")
            if details.get("profiles_csv"):
                target = staging / "profiles.tsv"
                source = self._reuse(ledger, target, details["profiles_csv"])
                if source is None:
                    source = self._transfer(details["profiles_csv"], target,
                                            limit=min(FILE_LIMIT, max_bytes - transferred),
                                            cancelled=cancelled, progress=progress,
                                            label="Profile table")
                    self._record(ledger, resume_root, source)
                transferred += source["bytes"]
                sources.append(source)
            metadata = {
                "name": f"{current['organism']} — {current['name']}",
                "organism": current["organism"], "genus": current["genus"],
                "species": current["species"], "description": current["name"],
                "locus_count": current["locus_count"], "type": current["type"],
                "source": current["provider"], "API": current["url"],
                "last_updated": current["last_updated"], "access_notice": current["access_notice"],
                "profile_field": str(details.get("primary_key_field") or "ST").rsplit("/", 1)[-1],
                "authenticated": False,
            }
            (staging / "scheme.json").write_text(json.dumps(metadata, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
            checked = _validated_scheme(staging, cancelled=cancelled)
            if checked.locus_count != current["locus_count"]:
                raise CatalogError("Downloaded reference loci do not match the online scheme.")
            after = self._json(current["url"], cancelled=cancelled, refresh=True)
            for field in ("loci", "last_updated", "records", "message"):
                if details.get(field) != after.get(field):
                    raise CatalogError("PubMLST changed during the download; retry to create a consistent snapshot.")
            source_digest = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
            manifest = {
                "format_version": 1, "provider": current["provider"], "entry": current,
                "retrieved_at": datetime.now(UTC).isoformat(), "sources": sources,
                "source_digest": source_digest, "scheme_digest": checked.digest,
                "access_notice": current["access_notice"], "authenticated": False,
                "downloaded_bytes": transferred,
                "excluded_alleles": excluded_alleles,
                "snapshot_note": "Sources were downloaded separately; the server does not provide a transactional multi-file snapshot.",
            }
            (staging / "reference_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            destination = install_destination(library_root, current, checked.digest)
            check_cancelled(cancelled)
            # A pre-created, labelled cgMLST slot is the destination, not an
            # existing snapshot: only allele files mean a scheme is already there.
            occupied = destination.is_dir() and cgmlst_schemes.has_alleles(destination)
            created = not occupied
            if occupied:
                if _validated_scheme(destination, cancelled=cancelled).digest != checked.digest:
                    raise CatalogError("Existing reference snapshot was modified; it will not be overwritten.")
                shutil.rmtree(staging)
            else:
                (staging / RESUME_LEDGER).unlink(missing_ok=True)
                destination.parent.mkdir(parents=True, exist_ok=True)
                # The partial folder sits beside the library, so a library root
                # that is itself a mount point needs a copying move.
                cgmlst_schemes.install_into(staging, destination)
            staging = None
            notes = list(checked.notes) + _exclusion_notes(excluded_alleles)
            if current["access_notice"] and current["access_notice"] not in notes:
                notes.append(current["access_notice"])
            if progress:
                progress(1, 1, "Public reference snapshot installed")
            return {"path": str(destination), "scheme_digest": checked.digest,
                    "source_digest": source_digest, "created": created,
                    "locus_count": checked.locus_count, "entry": current, "notes": notes,
                    "excluded_alleles": excluded_alleles, "resume_path": None}
        finally:
            if staging is not None:
                if keep_staging and staging.is_dir():
                    self._write_ledger(staging, ledger)
                else:
                    shutil.rmtree(staging, ignore_errors=True)

    # -- resume ledger ---------------------------------------------------------
    # The ledger is bookkeeping, never a scheme: it records which target files were
    # already received and their SHA-256, so an interrupted download continues
    # instead of restarting. A file is reused only when its bytes still hash to the
    # value recorded when they arrived.

    @staticmethod
    def _scheme_fingerprint(details, current) -> str:
        return hashlib.sha256(json.dumps(
            {"loci": details.get("loci"), "last_updated": details.get("last_updated"),
             "records": details.get("records"), "url": current["url"]},
            sort_keys=True).encode()).hexdigest()

    def _load_ledger(self, resume_root, details, current):
        if resume_root is None:
            return None
        fingerprint = self._scheme_fingerprint(details, current)
        path = Path(resume_root) / RESUME_LEDGER
        if path.is_file():
            try:
                stored = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError, ValueError):
                stored = None
            if (isinstance(stored, dict) and stored.get("fingerprint") == fingerprint
                    and isinstance(stored.get("files"), dict)):
                stored["written"] = 0
                return stored
            # The scheme changed upstream: a half-finished copy of the previous
            # definition must never be blended into the new one.
            shutil.rmtree(resume_root, ignore_errors=True)
        return {"format_version": 1, "fingerprint": fingerprint, "url": current["url"],
                "files": {}, "written": 0}

    @staticmethod
    def _write_ledger(resume_root, ledger) -> None:
        if ledger is None or resume_root is None:
            return
        payload = {key: value for key, value in ledger.items() if key != "written"}
        try:
            (Path(resume_root) / RESUME_LEDGER).write_text(
                json.dumps(payload, sort_keys=True), encoding="utf-8")
        except OSError:
            pass

    @staticmethod
    def _reuse(ledger, target, url):
        if ledger is None:
            return None
        recorded = ledger["files"].get(Path(target).name)
        if not isinstance(recorded, dict) or recorded.get("url") != url or not target.is_file():
            return None
        expected = recorded.get("stored_sha256") or recorded.get("sha256")
        if not expected or _file_digest(target) != expected:
            return None
        return {key: value for key, value in recorded.items() if key != "stored_sha256"} | (
            {"stored_sha256": recorded["stored_sha256"]} if recorded.get("stored_sha256") else {})

    def _record(self, ledger, resume_root, source) -> None:
        if ledger is None:
            return
        ledger["files"][source["file"]] = source
        ledger["written"] = ledger.get("written", 0) + 1
        if ledger["written"] % 25 == 0:
            self._write_ledger(resume_root, ledger)

    def update_scheme(self, installed_path, library_root=None, **kwargs) -> dict:
        installed_path = Path(installed_path)
        try:
            manifest = json.loads((installed_path / "reference_manifest.json").read_text(encoding="utf-8"))
            entry = manifest["entry"]
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise CatalogError("This scheme has no PubMLST snapshot manifest; select it in the online catalog.") from error
        return self.download_scheme(entry, library_root or installed_path.parent, **kwargs)


class PasteurCatalog(PubMLSTCatalog):
    """Official BIGSdb-Pasteur API; the same conservative snapshot protocol."""

    service_name = "BIGSdb-Pasteur"

    def __init__(self, cache_dir=None, **kwargs):
        super().__init__(cache_dir, base_url="https://bigsdb.pasteur.fr/api", **kwargs)


class _CGCatalogParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.cells = []
        self.cell = None
        self.href = ""

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.cells, self.href = [], ""
        elif tag == "td":
            self.cell = []
        elif tag == "a" and self.cell is not None:
            link = dict(attrs).get("href", "")
            if "/schema/" in link:
                self.href = link

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag == "td" and self.cell is not None:
            self.cells.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.href and len(self.cells) >= 2:
            self.rows.append((self.href, list(self.cells)))


class CGMLSTOrgCatalog(PubMLSTCatalog):
    """Official cgMLST.org catalog, locus tables and streaming allele ZIP bundles.

    The current public HTML table is parsed; no stale hardcoded organism counts
    or clinical clustering thresholds are inferred from its scheme definitions.
    """

    service_name = "cgMLST.org"
    RESUMABLE = False
    TERMS_URL = 'https://www.cgmlst.org/serverpolicy.html'
    TERMS_NOTICE = ('cgMLST.org nomenclature belongs to Ridom GmbH. Individual downloads are limited to non-commercial use; publication acknowledgement is required. Reuse of database copies in a product or service requires permission. Confirm that your intended use is permitted before downloading. This software does not grant data redistribution rights.')

    def download_estimate(self, entry) -> dict:
        count = int(entry.get("locus_count") or 0)
        return {"requests": 2, "per_target_requests": False, "resumable": False,
                "notice": (f"One archive of all {count} target allele sets, plus a locus table to "
                           "verify it. The archive is tens of megabytes, but the installed scheme "
                           "is much larger: a scheme of this size can occupy several gigabytes on "
                           "disk. Check free space before starting.")}

    def list_organisms(self, *, cancelled=None, refresh=False):
        parser = _CGCatalogParser()
        raw = self._transfer("https://www.cgmlst.org/ncs", cancelled=cancelled)
        if b"Illegal Typing ID" in raw or b"ERROR Occured" in raw:
            raise CatalogError("cgMLST.org reported an error page instead of its scheme catalogue.")
        parser.feed(raw.decode("utf-8"))
        entries = []
        for href, cells in parser.rows:
            _safe_url(urllib.parse.urljoin("https://www.cgmlst.org/ncs/", href))
            slug = _safe_name(urllib.parse.urlsplit(href).path.rstrip("/").rsplit("/", 1)[-1])
            label = cells[0].removesuffix(" cgMLST").strip()
            try:
                count = int(cells[1].replace(",", "").strip())
            except ValueError as error:
                raise CatalogError("cgMLST.org returned an invalid target count.") from error
            entries.append({"id": f"cgMLST.org:{slug}", "name": label, "organism": label,
                            "url": f"https://www.cgmlst.org/ncs/schema/{slug}/",
                            "locus_count": count, "slug": slug, **organism_parts(label)})
        if not entries:
            raise CatalogError("The cgMLST.org catalog layout could not be read; no schemes were assumed.")
        self._cg_entries = entries
        return entries

    def list_schemes(self, organism, *, cancelled=None, progress=None, refresh=False,
                     scheme_type=None, min_loci=0, max_loci=None):
        check_cancelled(cancelled)
        if isinstance(organism, str):
            matches = [entry for entry in self.list_organisms(cancelled=cancelled)
                       if organism in {entry["id"], entry["name"], entry["url"]}]
            if len(matches) != 1:
                raise CatalogError("Choose one organism from the cgMLST.org catalog.")
            organism = matches[0]
        count = organism["locus_count"]
        if count < min_loci or (max_loci is not None and count > max_loci) or (
            scheme_type and scheme_type.casefold() != "cgmlst"
        ):
            return []
        return [{**organism, "name": f"{organism['organism']} cgMLST", "type": "cgMLST",
                 "provider": "cgMLST.org", "last_updated": None, "profile_count": 0,
                 "has_profiles": False,
                 "terms_url": self.TERMS_URL, "terms_notice": self.TERMS_NOTICE,
                 "access_notice": "Allele nomenclature is public; no central ST/profile table is supplied by this download."}]

    def download_scheme(self, entry, library_root, *, cancelled=None, progress=None,
                        max_bytes=SNAPSHOT_LIMIT, terms_acknowledged=False, resume=True):
        check_cancelled(cancelled)
        if terms_acknowledged is not True:
            raise CatalogError('Review the cgMLST.org server policy and explicitly acknowledge that your intended use is permitted before downloading: ' + self.TERMS_URL)
        slug = _safe_name(str(entry.get("slug") or ""))
        base = f"https://www.cgmlst.org/ncs/schema/{slug}/"
        # Every cgMLST.org scheme is gene-by-gene, so it installs into the cgMLST
        # library whichever library the caller named. Staging inside it keeps the
        # final publication a rename on the same filesystem.
        library_root = cgmlst_schemes.install_root(Path(library_root)).resolve()
        library_root.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".cgmlst-download-", dir=library_root))
        excluded_alleles = []
        try:
            table_path = staging / "source_loci.tsv"
            table_source = self._transfer(base + "locus/?content-type=csv", table_path,
                                          limit=min(JSON_LIMIT, max_bytes), cancelled=cancelled,
                                          progress=progress, label="Locus table")
            table_text = table_path.read_text(encoding="utf-8-sig")
            # The service answers an unknown scheme with HTTP 200 and an HTML error
            # body, so a wrong identifier must be named rather than reported as a
            # malformed table.
            if "Illegal Typing ID" in table_text or "ERROR Occured" in table_text:
                raise CatalogError(
                    f"cgMLST.org does not recognise the scheme identifier {slug!r}; refresh the "
                    "catalogue and choose the scheme again.")
            rows = list(csv.reader(io.StringIO(table_text), delimiter="\t"))
            if not rows or not rows[0] or rows[0][0].strip() != "Locus":
                raise CatalogError("cgMLST.org did not return its expected locus table.")
            expected = set()
            for row in rows[1:]:
                while row and not row[0].strip():
                    row = row[1:]
                if row:
                    locus = _safe_name(row[0].strip())
                    if locus.casefold() in {name.casefold() for name in expected}:
                        raise CatalogError("cgMLST.org locus table contains duplicate names.")
                    expected.add(locus)
            if not expected or len(expected) != entry["locus_count"]:
                raise CatalogError(
                    f"The scheme now lists {len(expected)} targets but the catalogue entry says "
                    f"{entry['locus_count']}; refresh the catalogue and select the scheme again. A "
                    "cutoff published for one target set does not carry over to another.")
            archive = staging / "source_alleles.zip"
            archive_source = self._transfer(base + "alleles/", archive,
                                            limit=min(FILE_LIMIT, max_bytes), cancelled=cancelled,
                                            progress=progress,
                                            label=f"{entry['organism']} allele archive")
            sources, seen, unpacked = [table_source, archive_source], set(), table_source['bytes']
            with zipfile.ZipFile(archive) as bundle:
                for member in bundle.infolist():
                    check_cancelled(cancelled)
                    if member.is_dir():
                        continue
                    path_parts = member.filename.replace("\\", "/").split("/")
                    if (member.filename.startswith(("/", "\\")) or ".." in path_parts
                            or re.match(r"^[A-Za-z]:", member.filename)
                            or (member.external_attr >> 16) & 0o170000 == 0o120000):
                        raise CatalogError("Unsafe path or symbolic link in cgMLST.org archive.")
                    filename = _safe_name(path_parts[-1])
                    if Path(filename).suffix.casefold() not in {".fasta", ".fa", ".fna", ".tfa"}:
                        continue
                    locus = Path(filename).stem
                    if locus not in expected:
                        raise CatalogError(f"Archive locus {locus} is not in the downloaded locus table.")
                    if locus.casefold() in seen:
                        raise CatalogError("Archive contains colliding locus filenames.")
                    seen.add(locus.casefold())
                    unpacked += member.file_size
                    if unpacked > max_bytes:
                        raise CatalogError("The allele archive exceeds the configured extraction limit.")
                    digest, size = hashlib.sha256(), 0
                    target = staging / f"{locus}.fasta"
                    with bundle.open(member) as source, target.open("wb") as output:
                        while chunk := source.read(65536):
                            check_cancelled(cancelled)
                            size += len(chunk)
                            if size > member.file_size or size > FILE_LIMIT:
                                raise CatalogError("An allele archive entry exceeds its declared size.")
                            output.write(chunk)
                            digest.update(chunk)
                    source = {"url": base + "alleles/", "file": target.name,
                              "sha256": digest.hexdigest(), "bytes": size}
                    dropped = screen_allele_file(target, locus, cancelled=cancelled)
                    if dropped:
                        excluded_alleles.extend(dropped)
                        source["excluded_alleles"] = dropped
                        source["stored_sha256"] = _file_digest(target)
                    sources.append(source)
                    if progress:
                        progress(len(seen), len(expected), f"Validated archive locus {locus}")
            if seen != {name.casefold() for name in expected}:
                raise CatalogError("The allele archive is missing loci; no partial scheme was installed.")
            metadata = {"name": entry["name"], "organism": entry["organism"], "type": "cgMLST",
                        "source": "cgMLST.org", "API": base, "locus_count": len(expected),
                        "access_notice": entry["access_notice"], 'terms_url': self.TERMS_URL,
                        'terms_notice': self.TERMS_NOTICE, **organism_parts(entry["organism"])}
            (staging / "scheme.json").write_text(json.dumps(metadata, sort_keys=True, indent=2), encoding="utf-8")
            # The original archive is redundant after extraction; all its bytes
            # remain auditable through the source hash in the snapshot manifest.
            archive.unlink()
            checked = _validated_scheme(staging, cancelled=cancelled)
            source_digest = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
            manifest = {"format_version": 1, "provider": "cgMLST.org", "entry": entry,
                        'terms_url': self.TERMS_URL, 'terms_notice': self.TERMS_NOTICE,
                        'user_acknowledged_permitted_use': True,
                        "retrieved_at": datetime.now(UTC).isoformat(), "sources": sources,
                        "source_digest": source_digest, "scheme_digest": checked.digest,
                        "excluded_alleles": excluded_alleles,
                        "downloaded_bytes": archive_source["bytes"] + table_source["bytes"]}
            (staging / "reference_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            destination = cgmlst_schemes.install_folder(
                library_root, {**entry, "type": "cgMLST", "slug": slug}, digest=checked.digest)
            check_cancelled(cancelled)
            # An empty labelled slot is where this scheme belongs; only allele files
            # mean a snapshot is already installed there.
            occupied = destination.is_dir() and cgmlst_schemes.has_alleles(destination)
            created = not occupied
            if occupied:
                if _validated_scheme(destination, cancelled=cancelled).digest != checked.digest:
                    raise CatalogError("The existing reference snapshot was modified and will not be overwritten.")
            else:
                cgmlst_schemes.install_into(staging, destination)
                staging = None
            return {"path": str(destination), "scheme_digest": checked.digest,
                    "source_digest": source_digest, "created": created,
                    "locus_count": checked.locus_count, "entry": entry,
                    "notes": list(checked.notes) + _exclusion_notes(excluded_alleles),
                    "excluded_alleles": excluded_alleles, "resume_path": None}
        finally:
            if staging is not None:
                shutil.rmtree(staging)
