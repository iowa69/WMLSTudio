"""Build gates for the pinned reference panels a portable build does and does not carry.

Two separate obligations meet here. The characterization starter travels inside
the ZIP, so it must be provably intact and self-describing before it is packaged
and again after it is frozen. The practice cohorts and the broad species panel
are deliberately *not* bundled: they are sequence data the user downloads to the
adjacent Data root on an explicit request, so a copy that leaked into the source
tree has to fail the build instead of shipping.

Nothing here downloads anything. It reads what staging already produced.
"""

from __future__ import annotations

import os
from pathlib import Path

from wmlstudio import organism_panel, practice_cohorts
from wmlstudio.characterization_refs import validate_characterization_references

# The three manifest sections the organism-specific typing modules read. A
# format 1 snapshot predates all of them; that is a diagnosable state the assays
# report by name, never a negative result.
MODULE_SECTIONS = ("locus_profiles", "sccmec", "capsule")
# Derived, not spelled out, so a change to the download layout cannot leave this
# gate looking for a directory name that no longer exists.
COHORT_DIRECTORY = practice_cohorts.default_destination(Path("."), practice_cohorts.cohort_names()[0]).parent.name
UNBUNDLED_NOTE = ("Practice cohorts and the broad species panel are sequence data downloaded to the "
                  "user's own Data root on an explicit request. They never travel in the portable build.")


def assay_module_imports():
    """Return the assay modules PyInstaller's bytecode scan cannot discover.

    ``organism_modules._load`` reaches its assay modules through
    ``__import__(f"{__package__}.{name}")``, which compiles to a call of the
    builtin rather than an IMPORT_NAME opcode. Nothing else in the package
    imports them, so without an explicit hidden import the frozen build omits
    them and the first characterization run raises ModuleNotFoundError.
    """
    from wmlstudio.organism_modules import registered_modules
    return sorted({module.runner.__module__ for module in registered_modules().values()})


def panel_summary(root):
    """Re-verify a staged characterization snapshot and describe what it carries."""
    root = Path(root)
    manifest = validate_characterization_references(root)
    sccmec = manifest.get("sccmec") or {}
    capsule = manifest.get("capsule") or {}
    missing = [section for section in MODULE_SECTIONS if not manifest.get(section)]
    summary = {
        "path": str(root), "format_version": manifest["format_version"],
        "reference_digest": manifest["reference_digest"],
        "sources": {key: f"{value['revision']} ({value['license']})"
                    for key, value in (manifest.get("sources") or {}).items()},
        "species_references": len(manifest["species"]),
        "virulence_loci": sorted(manifest["virulence"]),
        "locus_profiles": sorted(manifest.get("locus_profiles") or {}),
        "sccmec_targets": len(sccmec.get("targets") or []),
        "sccmec_regions": len(sccmec.get("regions") or []),
        "sccmec_types": [definition["name"] for definition in (sccmec.get("rules") or {}).get("types") or []],
        "sccmec_not_assayed": list(sccmec.get("not_assayed") or []),
        "capsule_loci": {entry["gene"]: entry.get("allele_count") for entry in capsule.get("loci") or []},
        "stored_bytes": sum(entry["bytes"] for entry in manifest["files"]),
        "organism_modules": "absent" if missing else "staged",
        "organism_modules_reason": "" if not missing else (
            f"This snapshot (format {manifest['format_version']}) carries no "
            f"{', '.join(missing)} section. Re-stage it with studio_scripts/stage_characterization.py. "
            "Until then the organism-specific assays report not_run naming the missing section, which "
            "is a diagnosable reference state and never a negative finding."),
    }
    return summary


def sccmec_region_reference(root, subtype="IVa"):
    """Locate one staged SCCmec cassette reference, or None when the panel is absent."""
    root = Path(root)
    manifest = validate_characterization_references(root)
    entries = (manifest.get("sccmec") or {}).get("regions") or []
    entry = next((entry for entry in entries if entry["gene"] == subtype), None)
    return None if entry is None else root / entry["path"]


def unbundled_payload(paths):
    """Report practice-cohort genomes or a staged species panel found under these paths.

    Returns ``[(path, what)]``. The caller decides whether that is a build
    failure (it is) or a report line (it is that too).
    """
    genomes = {practice_cohorts.genome_filename(entry) for name in practice_cohorts.cohort_names()
               for entry in practice_cohorts.cohort_genomes(name)}
    found = []
    for base in paths:
        base = Path(base)
        if base.is_file():
            if base.name in genomes:
                found.append((str(base), "practice-cohort genome"))
            continue
        for directory, names, files in os.walk(base):
            for name in names:
                if name == COHORT_DIRECTORY:
                    found.append((str(Path(directory) / name), "practice-cohort directory"))
                elif name == organism_panel.PANEL_DIRECTORY or name.startswith(organism_panel.PANEL_PREFIX):
                    found.append((str(Path(directory) / name), "broad species panel"))
            found += [(str(Path(directory) / name), "practice-cohort genome")
                      for name in files if name in genomes]
    return sorted(found)


def verify_bundle_payload(sources):
    """Refuse to package anything the portable build is defined not to contain."""
    sources = list(sources)
    found = unbundled_payload(sources)
    if found:
        listing = "\n".join(f"  {what}: {path}" for path, what in found)
        raise SystemExit(f"{UNBUNDLED_NOTE}\nRemove these before packaging:\n{listing}")
    return {"checked": len(sources), "unbundled_payload": [], "note": UNBUNDLED_NOTE}
