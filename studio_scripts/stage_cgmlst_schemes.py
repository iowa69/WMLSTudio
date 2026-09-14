"""Stage the cgMLST scheme library: pinned definitions always, alleles only if licensed.

Two separate jobs, deliberately kept apart:

1. THE DEFINITION PACK (default, kilobytes, safe to ship). One folder per catalogued
   scheme, pre-created so a first run never asks a microbiologist to make a
   directory, each carrying a README naming the scheme, its provider, its terms and
   how to install it, plus the pinned target list. No allele sequence is involved,
   so nothing here is restricted data and nothing here is committed to the
   repository (src/wmlstudio/resources/ is gitignored).

2. THE ALLELE DATA (--with-alleles, gigabytes, licence-gated). Only schemes whose
   provider grants redistribution may be staged into a release, and the script
   refuses every other key by name. Measured on 2026-09-14: one cgMLST allele
   database occupies 0.8 GB (E. faecium, cgMLST.org) to 4.6 GB (S. aureus,
   cgMLST.org); a PubMLST cgMLST scheme is roughly 0.9 GB fetched as one request
   per target. A portable Windows release cannot carry several of these, so the
   default build stages definitions only and the user downloads what they need.

--verify re-reads every pinned target list from the live services and reports any
scheme whose target set has changed since the pins were taken. A changed target set
is a different quantity: the bound threshold no longer applies to it.

Never invoked by the application or at startup.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from wmlstudio import cgmlst_schemes  # noqa: E402
from wmlstudio.reference_catalog import (  # noqa: E402
    CatalogError,
    CGMLSTOrgCatalog,
    PasteurCatalog,
    PubMLSTCatalog,
)

DEFAULT_DESTINATION = ROOT / "src/wmlstudio/resources" / cgmlst_schemes.LIBRARY_DIRNAME
MANIFEST_NAME = "manifest.json"
TABLE_LIMIT = 8 * 1024 * 1024
USER_AGENT = "WMLSTudio-reference-staging"
# A courteous fixed pause between catalogue reads. PubMLST answers a burst with
# HTTP 429, which is how a naive "read every scheme" loop fails.
REQUEST_PAUSE = 2.0


def _read(url: str, *, accept: str = "*/*", attempts: int = 5) -> bytes:
    for attempt in range(attempts):
        request = urllib.request.Request(
            url, headers={"User-Agent": USER_AGENT, "Accept": accept}, method="GET")
        try:
            with urllib.request.urlopen(request, timeout=40) as response:
                return response.read(TABLE_LIMIT + 1)
        except urllib.error.HTTPError as error:
            if error.code != 429 or attempt == attempts - 1:
                raise
            try:
                wait = min(120.0, float(str(error.headers.get("Retry-After")).strip()))
            except (TypeError, ValueError):
                wait = 5.0 * (attempt + 1)
            time.sleep(wait)
    raise RuntimeError(f"Could not read {url}")


def fetch_targets(entry: dict) -> list[str]:
    """The current target names of one catalogued scheme, from its own provider."""
    if entry["provider"] == "cgmlst.org":
        raw = _read(entry["source_url"] + "locus/?content-type=csv")
        text = raw.decode("utf-8-sig")
        if "Illegal Typing ID" in text or "ERROR Occured" in text:
            raise RuntimeError(f"{entry['key']}: the provider does not recognise this scheme id.")
        rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
        if not rows or not rows[0] or rows[0][0].strip() != "Locus":
            raise RuntimeError(f"{entry['key']}: unexpected locus table layout.")
        names = []
        for row in rows[1:]:
            while row and not row[0].strip():
                row = row[1:]
            if row:
                names.append(row[0].strip())
        return names
    body = json.loads(_read(entry["source_url"], accept="application/json").decode("utf-8-sig"))
    if not isinstance(body.get("loci"), list):
        raise RuntimeError(f"{entry['key']}: the scheme record does not list its targets.")
    return [urllib.parse.unquote(str(url).rsplit("/", 1)[-1]) for url in body["loci"]]


def target_digest(names) -> str:
    """SHA-256 of the sorted target names: the pin, insensitive to listing order."""
    return hashlib.sha256("\n".join(sorted(names)).encode()).hexdigest()


def verify(keys=None, *, pause: float = REQUEST_PAUSE) -> dict:
    """Compare each pinned target set against the live service, and say what moved."""
    report = {"checked": 0, "matched": [], "changed": [], "unavailable": [],
              "catalog_version": cgmlst_schemes.CATALOG_VERSION,
              "catalog_digest": cgmlst_schemes.catalog_digest(),
              "checked_utc": datetime.now(UTC).isoformat()}
    for entry in cgmlst_schemes.catalog_entries():
        if keys and entry["key"] not in keys:
            continue
        time.sleep(pause)
        report["checked"] += 1
        try:
            names = fetch_targets(entry)
        except (OSError, RuntimeError, ValueError) as error:
            report["unavailable"].append({"key": entry["key"], "reason": str(error)})
            continue
        digest = target_digest(names)
        if len(names) == entry["locus_count"] and digest == entry["target_list_sha256"]:
            report["matched"].append(entry["key"])
        else:
            report["changed"].append({
                "key": entry["key"], "pinned_count": entry["locus_count"],
                "current_count": len(names), "pinned_sha256": entry["target_list_sha256"],
                "current_sha256": digest,
                "consequence": ("The target set is no longer the one the pin describes. Any cutoff "
                                "bound to this scheme was published for the pinned set and does not "
                                "carry over. Re-pin deliberately, do not auto-adopt.")})
    return report


def stage_definitions(destination: Path, *, keys=None, fetch: bool = False,
                      pause: float = REQUEST_PAUSE) -> dict:
    """Create the library layout and its pinned descriptions under `destination`.

    With `fetch`, each scheme's current target list is downloaded and written as
    targets.txt -- but only when it still matches the pin, so a silently changed
    target set is reported instead of shipped.
    """
    targets, drift = {}, []
    if fetch:
        for entry in cgmlst_schemes.catalog_entries():
            if keys and entry["key"] not in keys:
                continue
            time.sleep(pause)
            try:
                names = fetch_targets(entry)
            except (OSError, RuntimeError, ValueError) as error:
                drift.append({"key": entry["key"], "reason": str(error)})
                continue
            if len(names) != entry["locus_count"] or target_digest(names) != entry["target_list_sha256"]:
                drift.append({"key": entry["key"],
                              "reason": (f"current target set ({len(names)} targets, "
                                         f"{target_digest(names)[:16]}) does not match the pin "
                                         f"({entry['locus_count']} targets, "
                                         f"{entry['target_list_sha256'][:16]}); not staged")})
                continue
            targets[entry["key"]] = names
    report = cgmlst_schemes.prepare_library(destination.parent if destination.name ==
                                            cgmlst_schemes.LIBRARY_DIRNAME else destination,
                                            keys=keys, targets=targets)
    report["target_lists_staged"] = sorted(targets)
    report["drift"] = drift
    manifest = {
        "format_version": 1,
        "created_utc": datetime.now(UTC).isoformat(),
        "catalog_version": cgmlst_schemes.CATALOG_VERSION,
        "catalog_digest": cgmlst_schemes.catalog_digest(),
        "licence_reviewed_on": cgmlst_schemes.LICENCE_REVIEWED_ON,
        "contains_sequence_data": False,
        "slots": [{"key": entry["key"], "slot": entry["slot"], "organism": entry["organism"],
                   "provider": entry["provider"], "locus_count": entry["locus_count"],
                   "bundled": entry["bundled"], "redistribution": entry["redistribution"],
                   "target_list_sha256": entry["target_list_sha256"],
                   "threshold_scheme_key": entry["threshold_scheme_key"]}
                  for entry in cgmlst_schemes.catalog_entries()
                  if not keys or entry["key"] in keys],
        "drift": drift,
    }
    Path(report["root"], MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def stage_alleles(key: str, destination: Path, *, max_bytes=None, terms_acknowledged=False) -> dict:
    """Download one scheme's allele data, refusing any the provider does not license.

    This is the licence gate for the build. A scheme whose provider does not grant
    redistribution is never staged into a release, no matter what flag is passed:
    the only way to obtain it is for the person who will use it to download it
    themselves, having seen and accepted that provider's terms.
    """
    entry = cgmlst_schemes.entry_for(key)
    if not entry["bundled"]:
        raise SystemExit(
            f"{key} may not be staged into a release. {entry['provider_name']} does not grant "
            f"redistribution: \"{entry['licence_restriction']}\" Pre-create the folder with the "
            "definition pack instead; the user downloads this scheme themselves.")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if entry["provider"] == "cgmlst.org":
        client, extra = CGMLSTOrgCatalog(), {"terms_acknowledged": bool(terms_acknowledged)}
        online = {"slug": entry["scheme_id"], "organism": entry["organism"],
                  "name": entry["scheme_name"], "locus_count": entry["locus_count"],
                  "url": entry["source_url"], "access_notice": ""}
    else:
        client = PasteurCatalog() if entry["provider"] == "pasteur" else PubMLSTCatalog()
        extra = {}
        online = {"url": entry["source_url"], "organism": entry["organism"]}
    started = time.monotonic()

    def progress(current, total, message):
        if total and current and current % 100 == 0:
            print(f"    {current}/{total} {message}", flush=True)

    try:
        result = client.download_scheme(online, destination, progress=progress,
                                        **({"max_bytes": max_bytes} if max_bytes else {}), **extra)
    except CatalogError as error:
        raise SystemExit(f"{key}: {error}") from error
    if result["locus_count"] != entry["locus_count"]:
        raise SystemExit(
            f"{key}: the provider now serves {result['locus_count']} targets but the catalogue pins "
            f"{entry['locus_count']}. Staged nothing: a different target set is a different scheme.")
    return {"key": key, "path": result["path"], "locus_count": result["locus_count"],
            "scheme_digest": result["scheme_digest"],
            "excluded_alleles": result.get("excluded_alleles", []),
            "seconds": round(time.monotonic() - started, 1),
            "stored_bytes": sum(item.stat().st_size for item in Path(result["path"]).iterdir()
                                if item.is_file()),
            "notes": result["notes"]}


def notices() -> str:
    """Attribution text for studio_packaging/THIRD_PARTY_NOTICES.md."""
    lines = ["## cgMLST scheme catalogue", "",
             f"Catalogue version {cgmlst_schemes.CATALOG_VERSION}; provider terms read "
             f"{cgmlst_schemes.LICENCE_REVIEWED_ON}.", "",
             "WMLSTudio ships pinned scheme *definitions* (target lists and identity), never "
             "allele sequence data it is not licensed to redistribute.", ""]
    for key, provider in cgmlst_schemes.PROVIDERS.items():
        schemes = [entry["key"] for entry in cgmlst_schemes.catalog_entries()
                   if entry["provider"] == key]
        lines.append(f"### {provider['name']}")
        lines.append("")
        lines.append(f"- Terms: {provider['terms_url']} (read {provider['reviewed_on']})")
        lines.append(f"- Redistribution: {provider['redistribution']}")
        if provider["quote"]:
            lines.append(f"- Stated: \"{provider['quote']}\"")
        lines.append(f"- {provider['restriction']}")
        lines.append(f"- May be packed into a WMLSTudio release: "
                     f"{'yes' if provider['may_bundle'] else 'NO'}")
        lines.append(f"- Catalogued schemes: {', '.join(schemes) if schemes else 'none'}")
        lines.append("")
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION,
                        help="Where the cgMLST library layout is created")
    parser.add_argument("--key", action="append", dest="keys", default=None,
                        help="Limit to one catalogue key; repeatable")
    parser.add_argument("--fetch-targets", action="store_true",
                        help="Download each scheme's target list and stage it as targets.txt")
    parser.add_argument("--verify", action="store_true",
                        help="Only re-check the pins against the live services")
    parser.add_argument("--with-alleles", metavar="KEY",
                        help="Also download allele data for one redistributable scheme")
    parser.add_argument("--allele-destination", type=Path, default=None,
                        help="Library root for --with-alleles (default: the staged layout)")
    parser.add_argument("--max-bytes", type=int, default=None,
                        help="Refuse an allele download larger than this")
    parser.add_argument("--accept-provider-terms", action="store_true",
                        help="Record that the operator accepted the provider's terms")
    parser.add_argument("--notices", action="store_true",
                        help="Print the third-party attribution block and exit")
    parser.add_argument("--pause", type=float, default=REQUEST_PAUSE,
                        help="Seconds between catalogue reads (rate-limit courtesy)")
    args = parser.parse_args(argv)

    if args.notices:
        print(notices())
        return 0
    keys = set(args.keys) if args.keys else None
    if keys:
        for key in sorted(keys):
            cgmlst_schemes.entry_for(key)
    if args.verify:
        report = verify(keys, pause=args.pause)
        print(json.dumps(report, indent=2))
        return 1 if report["changed"] else 0
    report = stage_definitions(args.destination, keys=keys, fetch=args.fetch_targets,
                               pause=args.pause)
    summary = {key: report[key] for key in
               ("root", "catalog_version", "catalog_digest", "created", "target_lists_staged",
                "drift")}
    if args.with_alleles:
        summary["alleles"] = stage_alleles(
            args.with_alleles, args.allele_destination or Path(report["root"]),
            max_bytes=args.max_bytes, terms_acknowledged=args.accept_provider_terms)
    print(json.dumps(summary, indent=2))
    return 1 if report["drift"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
