"""A derived, regenerable organism view over the installed reference folders.

Every reference store in this application is content-addressed: a scheme folder's
name carries its digest, and stored scheme_path values in existing projects point
at those exact folders. Moving reference bytes to organise them by organism would
break scheme digests, snapshot re-validation and saved projects alike.

So organisation is delivered as a derived index instead: small JSON pointers in
genus/species folders, built from the labels each reference already carries,
rebuilt on demand and safe to delete. An organism folder here reflects what a
reference says about itself; it is not an independent check of its contents.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import threading
import urllib.parse
from collections import OrderedDict
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

from .identification import scheme_organism
from .project import CGMLST_LOCUS_FLOOR
from .sequence import check_cancelled, file_signature
from .storage import safe_component
from .typing import ALLELE_SUFFIXES

INDEX_DIRNAME = "_by_organism"
UNRESOLVED_NODE = "_Organism_not_recorded"
UNRESOLVED_LABEL = "Organism not recorded in this reference"
# One floor, defined once, in project.CGMLST_LOCUS_FLOOR: a stored profile and the
# scheme it was called against must never end up on opposite sides of it.
CGMLST_LOCUS_THRESHOLD = CGMLST_LOCUS_FLOOR
KINDS = ("mlst", "cgmlst", "unknown")
# A core target set is expected in every isolate of the organism; an accessory or
# whole-genome set is not. They are different quantities and never share a
# threshold, so which one a scheme is gets recorded rather than assumed.
TARGET_SETS = ("core", "accessory")
_ACCESSORY_TYPES = {"wgmlst", "accessory", "wgmlst-accessory", "accessory genome"}
_CORE_TYPES = {"cgmlst", "core", "cgmlst-core", "core genome"}
_PROVIDER_LABELS = {
    "pubmlst": "PubMLST", "bigsdb-pasteur": "BIGSdb-Pasteur", "pasteur": "BIGSdb-Pasteur",
    "cgmlst.org": "cgMLST.org", "cgmlst_org": "cgMLST.org", "ridom": "cgMLST.org",
    "enterobase": "EnteroBase", "chewie-ns": "chewie-NS", "local-ad-hoc": "Local ad-hoc build",
}
_PROVIDER_HOSTS = {
    "rest.pubmlst.org": "PubMLST", "pubmlst.org": "PubMLST",
    "bigsdb.pasteur.fr": "BIGSdb-Pasteur",
    "www.cgmlst.org": "cgMLST.org", "cgmlst.org": "cgMLST.org",
}
_METADATA_BYTES = 4 * 1024 * 1024
_CACHE: OrderedDict = OrderedDict()
_CACHE_LIMIT = 512
_LOCK = threading.Lock()
INDEX_README = """Derived index — safe to delete.

These folders are generated from the organism labels recorded in each installed
reference. They contain pointers, not sequence data, and they are rebuilt on
demand. An organism folder here reflects the label the reference itself carries;
it is not an independent verification of what that reference contains.
"""


def filter_scheme_locations(paths):
    """Drop derived and hidden directories from a list of candidate scheme folders.

    The index lives beside the schemes it describes, and it is not a scheme. Names
    beginning with '_' or '.' are reserved: downloaded scheme folders can never
    start with either, so nothing real is hidden by this filter.
    """
    return [path for path in map(Path, paths) if not path.name.startswith(("_", "."))]


def _signature(path: Path) -> tuple:
    if not path.is_dir():
        return (str(path), "missing")
    return (str(path), tuple((p.name, file_signature(p)) for p in sorted(path.iterdir())
                             if p.is_file()))


def _directory_digest(signature: tuple) -> str:
    """A cheap identity for a folder's current contents.

    This is deliberately NOT the scheme digest: Scheme.digest hashes every allele
    file and is what typing results are bound to. This only says whether the
    folder changed since the index was built.
    """
    return hashlib.sha256(repr(signature).encode()).hexdigest()


def _metadata(path: Path) -> dict:
    files = sorted(p for p in path.iterdir()
                   if p.is_file() and (p.name == "scheme.json" or p.name.endswith("_info.json")))
    for candidate in files:
        if candidate.stat().st_size > _METADATA_BYTES:
            continue
        try:
            value = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, ValueError):
            continue
        if isinstance(value, dict):
            return value
    return {}


def _locus_count(path: Path) -> int:
    return sum((p.with_suffix("") if p.suffix.casefold() in {".gz", ".bz2"} else p)
               .suffix.casefold() in ALLELE_SUFFIXES for p in path.iterdir() if p.is_file())


def _organism_basis(metadata: dict, organism: dict) -> str:
    if not organism.get("genus"):
        return "unresolved"
    label = metadata.get("organism")
    if isinstance(label, dict) or (isinstance(label, str) and label) or metadata.get("genus"):
        return "metadata"
    return "api_slug"


def _declared_type(metadata: dict) -> str:
    """The typing kind a scheme records for itself, normalised and never invented."""
    return str(metadata.get("type") or "").strip().casefold().replace(" ", "")


def classify_scheme(metadata: dict, locus_count: int) -> dict:
    """Say whether an installed folder is a classical MLST or a gene-by-gene scheme.

    The basis travels with the answer because the evidence differs in strength. The
    installed target count is a measured fact about this folder and decides first
    whenever it clears the floor: offering a two-thousand-target scheme where a
    seven-locus one is expected produces nonsense, whatever the folder calls
    itself. Below the floor the declared type decides, because a partial snapshot
    of a gene-by-gene scheme is still a gene-by-gene scheme. A folder holding no
    allele files is 'unknown'; it is never guessed into a kind.

    A declared type that contradicts the count is not discarded — it is reported in
    ``kind_conflict`` so a person can see both facts and fix the scheme metadata.
    """
    declared = _declared_type(metadata)
    if not locus_count:
        return {"kind": "unknown", "kind_basis": "the folder holds no allele FASTA files",
                "kind_conflict": ""}
    if locus_count > CGMLST_LOCUS_THRESHOLD:
        conflict = ("" if declared != "mlst" else
                    f"the scheme metadata declares type 'MLST' but the folder holds "
                    f"{locus_count} targets, far above the {CGMLST_LOCUS_THRESHOLD}-locus floor; "
                    "the installed target count decides")
        return {"kind": "cgmlst",
                "kind_basis": f"{locus_count} targets, above the {CGMLST_LOCUS_THRESHOLD}-locus "
                              "core-genome floor",
                "kind_conflict": conflict}
    if declared in _CORE_TYPES | _ACCESSORY_TYPES:
        return {"kind": "cgmlst",
                "kind_basis": f"scheme metadata declares type {declared!r}",
                "kind_conflict": (f"only {locus_count} of its targets are installed; this is a "
                                  "partial snapshot, not the published target set")}
    if declared == "mlst":
        return {"kind": "mlst", "kind_basis": "scheme metadata declares type 'mlst'",
                "kind_conflict": ""}
    return {"kind": "mlst", "kind_basis": f"{locus_count} loci, at classical scheme size",
            "kind_conflict": ""}


def _target_set(metadata: dict, kind: str) -> str:
    """'core', 'accessory' or '' — read from the scheme, never inferred from size."""
    if kind != "cgmlst":
        return ""
    declared = _declared_type(metadata)
    if declared in _ACCESSORY_TYPES:
        return "accessory"
    if declared in _CORE_TYPES:
        return "core"
    return ""


def _provider(metadata: dict) -> str:
    """Who published this scheme, as the scheme itself records it."""
    for key in ("source", "provider", "provider_name"):
        value = str(metadata.get(key) or "").strip()
        if value:
            return _PROVIDER_LABELS.get(value.casefold(), value)
    host = urllib.parse.urlsplit(str(metadata.get("API") or "")).hostname or ""
    return _PROVIDER_HOSTS.get(host.casefold(), "")


def _version(metadata: dict) -> str:
    """The revision or date the scheme records, so two snapshots can be told apart."""
    for key in ("last_updated", "revision", "download_date", "retrieved_at"):
        value = str(metadata.get(key) or "").strip()
        if value:
            return value.split("T")[0]
    return ""


def _scheme_label(metadata: dict, path: Path, organism_label: str) -> str:
    """The scheme's own name, with the organism removed so a title never repeats it."""
    for key in ("description", "name"):
        value = str(metadata.get(key) or "").strip()
        if not value or value == path.name:
            continue
        if organism_label and value.casefold().startswith(organism_label.casefold()):
            trimmed = value[len(organism_label):].strip(" -–—:·")
            return trimmed or value
        return value
    declared = str(metadata.get("type") or "").strip()
    return declared or path.name


def _title(organism_label: str, scheme_label: str, count: int, kind: str,
           provider: str, version: str) -> str:
    """One readable row title: organism, scheme, size, provider, version.

    The unit is named because it is the point: 'targets' for a gene-by-gene scheme
    and 'loci' for a classical one are different quantities that never share a
    scale, and a row that shows only a number invites them to be compared.
    """
    unit = "targets" if kind == "cgmlst" else "loci"
    size = f"{count} {unit}" if count else "no allele files installed"
    # A scheme that records neither an organism nor a name falls back to its folder
    # name for both; saying it twice reads as a mistake rather than as a title.
    named = "" if scheme_label == organism_label else scheme_label
    parts = [organism_label or UNRESOLVED_LABEL, named, size, provider,
             f"updated {version}" if version else ""]
    return " · ".join(part for part in parts if part)


def _entry(path: Path) -> dict:
    metadata = _metadata(path)
    count = _locus_count(path)
    label, organism = scheme_organism(SimpleNamespace(metadata=metadata, name=path.name))
    classified = classify_scheme(metadata, count)
    scheme_label = _scheme_label(metadata, path, label)
    provider = _provider(metadata)
    version = _version(metadata)
    return {"path": str(path), "id": path.name, "name": str(metadata.get("name") or path.name),
            "genus": organism.get("genus", ""), "species": organism.get("species", ""),
            "label": label, "locus_count": count, "kind": classified["kind"],
            "kind_basis": classified["kind_basis"], "kind_conflict": classified["kind_conflict"],
            "target_set": _target_set(metadata, classified["kind"]),
            "organism_label": label, "scheme_label": scheme_label, "provider": provider,
            "version": version,
            "title": _title(label, scheme_label, count, classified["kind"], provider, version),
            "organism_basis": _organism_basis(metadata, organism)}


def scheme_entries(scheme_paths, *, cancelled=None, kind=None) -> list[dict]:
    """One row per installed scheme folder, with the organism it records for itself.

    Allele files are counted but never parsed, so indexing 160 references stays
    cheap. Rows are cached on the folder's file signature and recomputed the
    moment anything inside it changes.

    ``kind`` restricts the rows to 'mlst', 'cgmlst' or 'unknown'. Ask for one kind
    wherever a widget offers a scheme to be typed against: a classical seven-locus
    picker that also lists a 2,000-target scheme is offering a different quantity
    under the same label.
    """
    if kind is not None and kind not in KINDS:
        raise ValueError(f"A scheme kind is one of {KINDS}, not {kind!r}.")
    entries = []
    for value in scheme_paths:
        check_cancelled(cancelled)
        path = Path(value).expanduser().resolve()
        if not path.is_dir():
            continue
        signature = _signature(path)
        with _LOCK:
            cached = _CACHE.get(signature)
            if cached is not None:
                _CACHE.move_to_end(signature)
        if cached is None:
            cached = {**_entry(path), "directory_digest": _directory_digest(signature)}
            with _LOCK:
                _CACHE[signature] = cached
                while len(_CACHE) > _CACHE_LIMIT:
                    _CACHE.popitem(last=False)
        if kind is None or cached["kind"] == kind:
            entries.append(dict(cached))
    entries.sort(key=lambda entry: (entry["genus"].casefold(), entry["species"].casefold(),
                                    entry["id"].casefold()))
    return entries


def scheme_library(scheme_paths, *, cancelled=None) -> dict[str, list[dict]]:
    """Installed schemes grouped into the two libraries a person browses, plus unknown.

    Every key is always present, so a view that renders a cgMLST tab shows an empty
    library rather than disappearing when nothing is installed yet. A folder whose
    kind could not be read stays in 'unknown' and is offered by neither tab.
    """
    grouped: dict[str, list[dict]] = {name: [] for name in KINDS}
    for entry in scheme_entries(scheme_paths, cancelled=cancelled):
        grouped.setdefault(entry["kind"], []).append(entry)
    return grouped


def clear_index_cache() -> None:
    with _LOCK:
        _CACHE.clear()


def panel_entries(root) -> list[dict]:
    """One row per taxon covered by an installed ANI panel, read from its manifest."""
    base = Path(root).expanduser()
    if not base.is_dir():
        return []
    manifests = [base / "manifest.json"] if (base / "manifest.json").is_file() else sorted(
        child / "manifest.json" for child in base.iterdir()
        if child.is_dir() and (child / "manifest.json").is_file())
    entries = []
    for manifest_path in manifests:
        if manifest_path.stat().st_size > 2 * 1024 * 1024:
            continue
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, ValueError):
            continue
        panel = manifest_path.parent
        for item in manifest.get("species", []):
            label = " ".join(filter(None, (item.get("genus"), item.get("species"))))
            entries.append({
                "path": str(panel), "id": str(item.get("id", "")),
                "name": str(item.get("id", "")), "genus": str(item.get("genus", "")),
                "species": str(item.get("species", "")),
                "label": label, "organism_label": label,
                "scheme_label": "ANI reference panel", "provider": "", "version": "",
                "title": f"{label or UNRESOLVED_LABEL} · ANI reference panel",
                "kind_basis": "listed in the panel manifest", "kind_conflict": "",
                "target_set": "",
                "locus_count": 0, "kind": "ani_panel", "organism_basis": "panel_manifest",
                "directory_digest": str(manifest.get("reference_digest", "")),
                "taxonomy_basis": str(item.get("taxonomy_basis", "")),
                "panel_kind": str(manifest.get("panel_kind", "characterization_starter"))})
    entries.sort(key=lambda entry: (entry["genus"].casefold(), entry["species"].casefold(),
                                    entry["id"].casefold()))
    return entries


def organism_tree(entries) -> dict:
    """Group rows as {genus: {species: [row, ...]}}; an unlabelled reference keeps
    its own explicit node rather than being dropped or guessed at."""
    tree: dict[str, dict[str, list]] = {}
    for entry in entries:
        genus = entry.get("genus") or ""
        species = entry.get("species") or ""
        tree.setdefault(genus, {}).setdefault(species, []).append(entry)
    return {genus: {species: rows for species, rows in sorted(children.items())}
            for genus, children in sorted(tree.items())}


def _pointer(entry: dict) -> dict:
    return {"scheme_path": entry["path"], "id": entry["id"], "name": entry["name"],
            "title": entry.get("title", ""),
            "directory_digest": entry.get("directory_digest", ""),
            "locus_count": entry.get("locus_count", 0), "kind": entry.get("kind", "unknown"),
            "kind_basis": entry.get("kind_basis", ""),
            "target_set": entry.get("target_set", ""),
            "organism_basis": entry.get("organism_basis", "unresolved"),
            "organism": {"genus": entry.get("genus", ""), "species": entry.get("species", "")},
            "generated_utc": datetime.now(UTC).isoformat()}


def write_organism_index(schemes_root, entries) -> Path:
    """Build <schemes_root>/_by_organism/ atomically from already-computed rows.

    Only JSON pointers and a README are written: no sequence bytes are copied, no
    symlink or junction is created, and no scheme folder is read from or touched.
    """
    root = Path(schemes_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    target = root / INDEX_DIRNAME
    with tempfile.TemporaryDirectory(prefix=".by-organism-", dir=root) as temporary:
        stage = Path(temporary) / INDEX_DIRNAME
        stage.mkdir()
        (stage / "README.txt").write_text(INDEX_README, encoding="utf-8", newline="\r\n")
        for genus, children in organism_tree(entries).items():
            for species, rows in children.items():
                folder = (stage / UNRESOLVED_NODE if not genus else
                          stage / safe_component(genus) / safe_component(species, "Unknown_species"))
                folder.mkdir(parents=True, exist_ok=True)
                for row in rows:
                    name = safe_component(row["id"], "reference")
                    pointer = folder / f"{name}.json"
                    if pointer.exists():
                        pointer = folder / f"{name}_{row.get('directory_digest', '')[:8]}.json"
                    pointer.write_text(json.dumps(_pointer(row), indent=2) + "\n", encoding="utf-8")
        if target.exists():
            os.replace(target, Path(temporary) / "previous")
        os.replace(stage, target)
    return target


def clear_organism_index(schemes_root) -> None:
    """Remove the derived index and nothing else."""
    target = Path(schemes_root).expanduser().resolve() / INDEX_DIRNAME
    if not target.is_dir() or target.is_symlink():
        return
    for path in sorted(target.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            path.rmdir()
    target.rmdir()
