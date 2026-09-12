"""Public, GET-only PubMLST reference discovery and immutable local snapshots.

Official protocol: https://bigsdb.readthedocs.io/en/latest/rest.html (checked
2026-09-12). No sample data, credentials, or submission requests are transmitted.
Unauthenticated PubMLST access can be restricted to an older public dataset; its
access message is retained and must not be presented as a fully current release.
"""

from __future__ import annotations

import csv
import hashlib
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

from . import __version__
from .sequence import check_cancelled
from .typing import load_scheme

API_ROOT = "https://rest.pubmlst.org"
OFFICIAL_HOSTS = {"rest.pubmlst.org", "bigsdb.pasteur.fr", "www.cgmlst.org", "cgmlst.org"}
JSON_LIMIT = 16 * 1024 * 1024
FILE_LIMIT = 1024 * 1024 * 1024
SNAPSHOT_LIMIT = 8 * 1024 * 1024 * 1024


class CatalogError(ValueError):
    """A reference could not be fetched or validated without an unsafe assumption."""


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

    def __init__(self, cache_dir=None, *, opener=None, timeout=20, retries=2,
                 base_url=API_ROOT):
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.opener = opener or urllib.request.build_opener(_OfficialRedirect())
        self.timeout = timeout
        self.retries = retries
        self.base_url = _safe_url(base_url).rstrip("/")
        self._cache: dict[str, tuple[float, object]] = {}
        self.catalog_from_cache = False

    def _transfer(self, url, destination=None, *, limit=JSON_LIMIT, cancelled=None):
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
                    access = " Public access may require authentication." if error.code in {401, 403} else ""
                    raise CatalogError(f"PubMLST returned HTTP {error.code}.{access}") from error
            except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
                if attempt == self.retries:
                    raise CatalogError(f"Could not reach PubMLST: {error}") from error
            # Short interruptible backoff; GET-only requests are safe to retry.
            for _ in range(5 * (attempt + 1)):
                check_cancelled(cancelled)
                threading.Event().wait(0.1)
        raise CatalogError("PubMLST download did not complete.")

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

    def download_scheme(self, entry, library_root, *, cancelled=None, progress=None,
                        max_bytes=SNAPSHOT_LIMIT) -> dict:
        """Install a fully validated snapshot; failures never publish a partial directory.

        A changed remote version is a new directory. Existing snapshots are never
        overwritten, and an existing destination is revalidated before reuse.
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
        staging = Path(tempfile.mkdtemp(prefix=".pubmlst-download-", dir=library_root))
        sources, names = [], set()
        transferred = 0
        try:
            for index, locus_url in enumerate(details["loci"], 1):
                check_cancelled(cancelled)
                _safe_url(locus_url)
                name = _safe_name(urllib.parse.unquote(urllib.parse.urlsplit(locus_url).path.rsplit("/", 1)[-1]))
                if name.casefold() in names:
                    raise CatalogError("The scheme repeats a locus filename, including Windows case folding.")
                names.add(name.casefold())
                source = self._transfer(locus_url + "/alleles_fasta", staging / f"{name}.tfa",
                                        limit=min(FILE_LIMIT, max_bytes - transferred), cancelled=cancelled)
                transferred += source["bytes"]
                sources.append(source)
                if progress:
                    progress(index, len(details["loci"]) + 1, f"Downloaded {name} alleles")
            if details.get("profiles_csv"):
                source = self._transfer(details["profiles_csv"], staging / "profiles.tsv",
                                        limit=min(FILE_LIMIT, max_bytes - transferred), cancelled=cancelled)
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
            checked = load_scheme(staging, cancelled=cancelled)
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
                "snapshot_note": "Sources were downloaded separately; the server does not provide a transactional multi-file snapshot.",
            }
            (staging / "reference_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            base = _safe_name(current["database"] + "_" + current["scheme_id"])
            destination = library_root / f"{base}_{checked.digest[:16]}"
            check_cancelled(cancelled)
            created = not destination.exists()
            if destination.exists():
                if load_scheme(destination, cancelled=cancelled).digest != checked.digest:
                    raise CatalogError("Existing reference snapshot was modified; it will not be overwritten.")
            else:
                staging.rename(destination)
                staging = None
            notes = list(checked.notes)
            if current["access_notice"]:
                notes.append(current["access_notice"])
            if progress:
                progress(1, 1, "Public reference snapshot installed")
            return {"path": str(destination), "scheme_digest": checked.digest,
                    "source_digest": source_digest, "created": created,
                    "locus_count": checked.locus_count, "entry": current, "notes": notes}
        finally:
            if staging is not None:
                shutil.rmtree(staging)

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

    TERMS_URL = 'https://www.cgmlst.org/serverpolicy.html'
    TERMS_NOTICE = ('cgMLST.org nomenclature belongs to Ridom GmbH. Individual downloads are limited to non-commercial use; publication acknowledgement is required. Reuse of database copies in a product or service requires permission. Confirm that your intended use is permitted before downloading. This software does not grant data redistribution rights.')

    def list_organisms(self, *, cancelled=None, refresh=False):
        parser = _CGCatalogParser()
        raw = self._transfer("https://www.cgmlst.org/ncs", cancelled=cancelled)
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
                        max_bytes=SNAPSHOT_LIMIT, terms_acknowledged=False):
        check_cancelled(cancelled)
        if terms_acknowledged is not True:
            raise CatalogError('Review the cgMLST.org server policy and explicitly acknowledge that your intended use is permitted before downloading: ' + self.TERMS_URL)
        slug = _safe_name(str(entry.get("slug") or ""))
        base = f"https://www.cgmlst.org/ncs/schema/{slug}/"
        library_root = Path(library_root).resolve()
        library_root.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".cgmlst-download-", dir=library_root))
        try:
            table_path = staging / "source_loci.tsv"
            table_source = self._transfer(base + "locus/?content-type=csv", table_path,
                                          limit=min(JSON_LIMIT, max_bytes), cancelled=cancelled)
            rows = list(csv.reader(io.StringIO(table_path.read_text(encoding="utf-8-sig")), delimiter="\t"))
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
                raise CatalogError("The current locus table and catalog target count disagree; refresh the catalog.")
            archive = staging / "source_alleles.zip"
            archive_source = self._transfer(base + "alleles/", archive,
                                            limit=min(FILE_LIMIT, max_bytes), cancelled=cancelled)
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
                    sources.append({"url": base + "alleles/", "file": target.name,
                                    "sha256": digest.hexdigest(), "bytes": size})
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
            checked = load_scheme(staging, cancelled=cancelled)
            source_digest = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
            manifest = {"format_version": 1, "provider": "cgMLST.org", "entry": entry,
                        'terms_url': self.TERMS_URL, 'terms_notice': self.TERMS_NOTICE,
                        'user_acknowledged_permitted_use': True,
                        "retrieved_at": datetime.now(UTC).isoformat(), "sources": sources,
                        "source_digest": source_digest, "scheme_digest": checked.digest,
                        "downloaded_bytes": archive_source["bytes"] + table_source["bytes"]}
            (staging / "reference_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            destination = library_root / f"cgmlst_org_{slug}_{checked.digest[:16]}"
            check_cancelled(cancelled)
            created = not destination.exists()
            if destination.exists():
                if load_scheme(destination, cancelled=cancelled).digest != checked.digest:
                    raise CatalogError("The existing reference snapshot was modified and will not be overwritten.")
            else:
                staging.rename(destination)
                staging = None
            return {"path": str(destination), "scheme_digest": checked.digest,
                    "source_digest": source_digest, "created": created,
                    "locus_count": checked.locus_count, "entry": entry, "notes": checked.notes}
        finally:
            if staging is not None:
                shutil.rmtree(staging)
