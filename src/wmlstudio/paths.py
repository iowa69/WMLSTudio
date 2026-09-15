"""Resolve resources independently of working directory and installed Python."""

import sys
from pathlib import Path

from PySide6.QtCore import QStandardPaths

SCHEME_DIRNAME = "schemes"
CGMLST_DIRNAME = "cgmlst"
# Counted, never parsed: a cgMLST library folder exists before anything is
# installed into it, so "has allele files" is what separates a labelled empty
# slot from an installed scheme. Repeated here rather than imported to keep this
# module free of the typing import chain.
_ALLELE_SUFFIXES = {".tfa", ".fa", ".fasta", ".fna"}


def resource_root() -> Path:
    return Path(__file__).resolve().parent / "resources"


def data_root() -> Path:
    if getattr(sys, "frozen", False):
        root = Path(sys.executable).resolve().parent / "Data"
    else:
        root = Path(QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation))
    root.mkdir(parents=True, exist_ok=True)
    return root


def mlst_library_roots(root: Path) -> list[Path]:
    """Where classical seven-locus schemes live: the bundled snapshot and the user's."""
    return [resource_root() / SCHEME_DIRNAME, Path(root) / SCHEME_DIRNAME]


def cgmlst_library_roots(root: Path) -> list[Path]:
    """Where gene-by-gene schemes live. Separate folder, separate library, separate tab."""
    return [resource_root() / CGMLST_DIRNAME, Path(root) / CGMLST_DIRNAME]


def _has_alleles(path: Path) -> bool:
    try:
        for item in path.iterdir():
            if not item.is_file():
                continue
            name = item.with_suffix("") if item.suffix.casefold() in {".gz", ".bz2"} else item
            if name.suffix.casefold() in _ALLELE_SUFFIXES:
                return True
    except OSError:
        return False
    return False


def scheme_locations(root: Path) -> list[Path]:
    """Every installed scheme directory, and nothing that merely sits beside one.

    Both libraries are searched: <data root>/schemes for the classical seven-locus
    schemes and <data root>/cgmlst for the gene-by-gene ones. Which kind a folder
    holds is decided from the folder itself by reference_index, never from which
    library it happens to sit in, so a cgMLST scheme a user imported by hand into
    the classical folder is still reported as cgMLST.

    Derived folders live here too — reference_index writes `_by_organism` into the
    same place — and handing one to the typing code would offer the user a scheme
    that does not exist. The cgMLST library additionally pre-creates a labelled,
    empty folder for every catalogued scheme; an empty slot is a description, not
    an installed scheme, so only folders that actually hold allele files are
    returned from it.
    """
    from wmlstudio.reference_index import filter_scheme_locations
    classical = {p for base in mlst_library_roots(root) if base.is_dir()
                 for p in base.iterdir() if p.is_dir()}
    gene_by_gene = {p for base in cgmlst_library_roots(root) if base.is_dir()
                    for p in base.iterdir() if p.is_dir() and _has_alleles(p)}
    return filter_scheme_locations(sorted(classical | gene_by_gene,
                                          key=lambda p: (p.name.casefold(), str(p))))


def locations_by_kind(root: Path, *, cancelled=None) -> dict[str, list[Path]]:
    """Installed scheme folders split into 'mlst', 'cgmlst' and 'unknown'.

    A seven-locus MLST distance and a two-thousand-target cgMLST distance are
    different quantities, so the two libraries are offered separately and a folder
    whose kind cannot be read from what it records is listed as unknown rather
    than sorted into whichever tab is closest.
    """
    from wmlstudio.reference_index import scheme_entries
    grouped: dict[str, list[Path]] = {"mlst": [], "cgmlst": [], "unknown": []}
    for entry in scheme_entries(scheme_locations(root), cancelled=cancelled):
        grouped.setdefault(entry["kind"], []).append(Path(entry["path"]))
    return grouped


def mlst_locations(root: Path, *, cancelled=None) -> list[Path]:
    """Installed classical seven-locus scheme folders only."""
    return locations_by_kind(root, cancelled=cancelled)["mlst"]


def cgmlst_locations(root: Path, *, cancelled=None) -> list[Path]:
    """Installed core-genome / whole-genome scheme folders only."""
    return locations_by_kind(root, cancelled=cancelled)["cgmlst"]
