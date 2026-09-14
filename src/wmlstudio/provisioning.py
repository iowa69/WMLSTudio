"""What this installation actually has, what it is missing, and how to get the rest.

Almost every capability here depends on something provisioned separately from the
code: the broad species panel is downloaded on request, BLAST+ is staged beside
the executable, the AMR database store is chosen by the user. When one of those
is absent the feature that needs it does not fail loudly. It returns an empty,
technically defensible and completely unexplained result: every sample
"unresolved", every organism folder never created, every HYDRA run refused. This
module exists so that the reason is a first-class readable object instead of a
silence, and so the fix is an action the caller can carry out.

Nothing is downloaded, installed, created or repaired by reading this report. A
probe looks at what is on disk and says so; the install helpers at the end run
only when a person asks for them. A requirement is never called ready because it
is expected to be there: the panel is validated, the tools are resolved, the
database manifest is read. "I did not check" is reported as such.

The module is free of Qt at import time, so a support call can run it in a plain
Python process, or through studio_scripts/check_setup.py, with no desktop session.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from . import __version__

REPORT_FORMAT_VERSION = 1
ORGANIZATION_NAME = "IOWA-BioTech"
APPLICATION_NAME = "WMLSTudio"
# ready: probed and usable. partial: usable, but a named piece of it is absent.
# missing: not installed at all. unusable: installed and broken or unsupported.
STATES = ("ready", "partial", "missing", "unusable")
# The same four states said in the words a non-expert reads, for one-line summaries.
STATE_WORDS = {"ready": "ready", "partial": "only partly installed",
               "missing": "not installed", "unusable": "installed but not usable"}
# The order a person meets these: import and file samples, type them, screen them.
# The first unmet required item in this order is the single next thing to do.
ORDER = ("species_panel", "assembly_runtime", "mlst_schemes", "cgmlst_schemes", "blast_tools",
         "hydra_engine", "hydra_database", "characterization")
# Which menu stops working when which requirement is unmet. A capability is ready
# only when every requirement it names is ready; nothing here is a partial credit.
CAPABILITIES = {
    "organism_filing": ("Automatic organism folders", ("species_panel",)),
    "assembly": ("Read assembly (FASTQ to contigs)", ("assembly_runtime",)),
    "mlst": ("Classical 7-locus MLST typing", ("mlst_schemes",)),
    "cgmlst": ("cgMLST / wgMLST typing", ("cgmlst_schemes", "blast_tools")),
    "hydra": ("HYDRA resistance screening", ("hydra_engine", "blast_tools", "hydra_database")),
    "organism_assays": ("Organism-specific assays (SCCmec, Klebsiella locus STs, capsule markers)",
                        ("characterization", "blast_tools")),
}
# The two stores HYDRA runs with by default; the protein store also carries the
# point-mutation catalogs and the organism names those catalogs are curated for.
HYDRA_DEFAULT_DATABASES = ("ncbi", "protein")
MANIFEST_BYTES = 2 * 1024 * 1024
TOOL_TIMEOUT = 30


@dataclass(frozen=True)
class Action:
    """The one thing that would fix a requirement, in the words a user will read."""

    label: str
    ui_route: str = ""
    detail: str = ""
    # True when the application can carry this out itself, from inside the app.
    automatic: bool = False
    # "module:function" a caller may invoke, and the single argument it needs.
    entry_point: str = ""
    argument: str = ""
    command: str = ""


@dataclass(frozen=True)
class Requirement:
    """One prerequisite: what it powers, what breaks without it, and its real state."""

    key: str
    title: str
    capability: str
    # A clause, not a sentence: "every imported sample is filed as unresolved".
    blocks: str
    required: bool
    state: str
    reason: str
    probe: str
    action: Action | None = None
    detail: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.state not in STATES:
            raise ValueError(f"Unknown provisioning state: {self.state!r}")

    @property
    def ready(self) -> bool:
        return self.state == "ready"

    @property
    def consequence(self) -> str:
        return f"Without this, {self.blocks}."

    def as_dict(self) -> dict:
        return {"key": self.key, "title": self.title, "capability": self.capability,
                "consequence": self.consequence, "blocks": self.blocks,
                "required": self.required, "state": self.state, "ready": self.ready,
                "reason": self.reason, "probe": self.probe,
                "action": asdict(self.action) if self.action else None,
                "detail": dict(self.detail)}


def _staging_platform() -> str:
    """The platform token the staging scripts take, for a copyable command line."""
    return "windows-x64" if os.name == "nt" else "linux-x64"


def resource_root() -> Path:
    """Mirror of paths.resource_root, repeated so this module never imports Qt."""
    return Path(__file__).resolve().parent / "resources"


def _fallback_data_root() -> Path:
    """The per-user folder Qt would choose, computed without Qt present."""
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData/Local")
    elif sys.platform == "darwin":
        base = Path.home() / "Library/Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local/share")
    return base / ORGANIZATION_NAME / APPLICATION_NAME


def default_data_root() -> Path:
    """The data folder the desktop application itself uses, resolved read-only.

    The application resolves this through Qt (paths.data_root) so the frozen and
    installed builds agree. Qt is imported here, not at module import, and the
    organisation/application names are set only when nothing has set them, so a
    command-line check inspects exactly the folder the desktop would. Unlike
    paths.data_root this never creates the folder: a diagnostic must not write.
    """
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent / "Data"
    try:
        from PySide6.QtCore import QCoreApplication, QStandardPaths
    except ImportError:  # pragma: no cover - the desktop dependency is always present
        return _fallback_data_root()
    if not QCoreApplication.organizationName():
        QCoreApplication.setOrganizationName(ORGANIZATION_NAME)
    if not QCoreApplication.applicationName():
        QCoreApplication.setApplicationName(APPLICATION_NAME)
    location = QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation)
    return Path(location) if location else _fallback_data_root()


def default_scheme_paths(data_root=None) -> list[Path]:
    """The folders the desktop application itself scans for typing schemes.

    Mirrors paths.scheme_locations, then drops the derived _by_organism index the
    way the workbench does, so a count here counts real schemes and nothing else.
    """
    from .reference_index import filter_scheme_locations
    root = Path(data_root) if data_root is not None else default_data_root()
    locations = [resource_root() / "schemes", root / "schemes"]
    found = {path for base in locations if base.is_dir()
             for path in base.iterdir() if path.is_dir()}
    return filter_scheme_locations(sorted(found, key=lambda path: path.name.casefold()))


def hydra_database_candidates(data_root=None, selected=None) -> list[tuple[str, Path]]:
    """Every AMR store the application would consider, most authoritative first.

    Mirrors ui_workbench.active_amr_database: an explicit project selection wins,
    then the snapshot staged inside the application, then the managed download
    folder. They are returned rather than collapsed, so a report can say which
    store would work when the selected one does not.
    """
    from .hydra_runtime import bundled_database_root
    root = Path(data_root) if data_root is not None else default_data_root()
    candidates = []
    if selected:
        candidates.append(("selected in this project", Path(selected)))
    bundled = bundled_database_root()
    if bundled is not None:
        candidates.append(("staged inside the application", Path(bundled)))
    candidates.append(("managed download folder", root / "references" / "hydra"))
    seen, unique = set(), []
    for origin, path in candidates:
        resolved = Path(path).expanduser().resolve()
        if str(resolved) not in seen:
            seen.add(str(resolved))
            unique.append((origin, resolved))
    return unique


def _snapshot_manifest(path, *, verify=True, cancelled=None) -> tuple[dict, str]:
    """Read a reference snapshot's manifest as deeply as the caller asked for.

    verify=True runs the application's own validator, which re-reads and re-hashes
    every file the manifest claims and refuses an unsupported format version.
    verify=False is the refresh-cheap probe: the manifest is parsed and every file
    it names is confirmed present at its recorded size. Neither one guesses, and
    the returned sentence says which of the two actually ran.
    """
    from .characterization_refs import validate_characterization_references
    path = Path(path)
    if verify:
        manifest = validate_characterization_references(path, cancelled=cancelled)
        return manifest, "every reference file was re-read and re-hashed against the manifest"
    manifest_path = path / "manifest.json"
    if manifest_path.stat().st_size > MANIFEST_BYTES:
        raise ValueError("The reference manifest exceeds the 2 MiB limit.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or not isinstance(manifest.get("files"), list):
        raise ValueError("The reference manifest does not list its files.")
    for entry in manifest["files"]:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise ValueError("The reference manifest contains a malformed file entry.")
        target = path / entry["path"]
        if not target.is_file() or target.stat().st_size != entry.get("bytes"):
            raise ValueError(f"A reference file is missing or a different size: {entry['path']}")
    return manifest, ("the manifest was parsed and every file it names was found at its "
                      "recorded size")


def species_panel_requirement(data_root=None, *, verify=True, cancelled=None) -> Requirement:
    """Probe the broad ANI panel that automatic organism identification depends on."""
    from .organism_panel import (
        PANEL_DIRECTORY,
        PANEL_REVISION,
        SPECIES_PANEL,
        installed_species_panel,
    )
    root = Path(data_root) if data_root is not None else default_data_root()
    store = root / PANEL_DIRECTORY
    megabytes = sum(row[5] for row in SPECIES_PANEL) / (1024 * 1024)
    action = Action(
        label="Install the species reference panel",
        ui_route="Samples > Reference data > Species panel",
        detail=(f"Downloads {len(SPECIES_PANEL)} pinned NCBI RefSeq genomes "
                f"(about {megabytes:.0f} MB) into {store}, checking every file against the "
                "checksum NCBI publishes for it. Your own sequences are never uploaded."),
        automatic=True,
        entry_point="wmlstudio.provisioning:install_species_panel",
        argument=str(root),
        command=f'python studio_scripts/stage_species_panel.py --root "{store}"')
    common = {"key": "species_panel", "title": "Species reference panel",
              "capability": "Automatic organism folders",
              "blocks": ("every imported sample is filed as unresolved instead of into an "
                         "organism folder, because nothing is installed to compare it against"),
              "required": True, "action": action}
    try:
        installed = installed_species_panel(root)
    except OSError as error:
        return Requirement(**common, state="unusable", probe=f"Listed {store}.",
                           reason=f"The species panel folder could not be read: {error}",
                           detail={"root": str(store)})
    if installed is None:
        return Requirement(
            **common, state="missing", probe=f"Looked for an installed panel under {store}.",
            reason=("No species reference panel is installed, so organism identification attempts "
                    "no genomic comparison at all: it returns 'unresolved' for every sample and no "
                    "organism folder is ever created. The panel is a download, not part of the "
                    "application, and nothing installs it automatically."),
            detail={"root": str(store), "expected_revision": PANEL_REVISION,
                    "reference_count": len(SPECIES_PANEL)})
    try:
        manifest, probe = _snapshot_manifest(installed, verify=verify, cancelled=cancelled)
    except (OSError, ValueError) as error:
        return Requirement(
            **common, state="unusable", probe=f"Validated {installed}.",
            reason=(f"A species panel is installed at {installed} but it did not pass validation, "
                    f"so it is not used: {error}"),
            detail={"root": str(store), "path": str(installed)})
    taxa = sorted({" ".join(filter(None, (entry.get("genus"), entry.get("species")))).strip()
                   for entry in manifest.get("species", [])})
    return Requirement(
        **common, state="ready", probe=f"Read {installed}: {probe}.",
        reason=(f"{len(manifest.get('species', []))} reference genomes covering {len(taxa)} taxa "
                "are installed, so an imported assembly is compared against them and can be "
                "proposed for an organism folder."),
        detail={"root": str(store), "path": str(installed), "taxa": taxa,
                "reference_count": len(manifest.get("species", [])),
                "reference_digest": manifest.get("reference_digest", ""),
                "source_revision": manifest.get("source_revision", ""),
                "expected_revision": PANEL_REVISION,
                "limitations": list(manifest.get("limitations", []))})


def assembly_runtime_requirement() -> Requirement:
    """Probe the native assembler, which the Assembly submenu needs for raw reads."""
    from .assembly import AssemblyError, resolve_skesa
    action = Action(
        label="Restore the SKESA assembler that ships beside the application",
        ui_route="Settings > Runtime tools",
        detail=("The portable build carries SKESA in Tools/skesa beside the executable; extract "
                "the complete folder rather than the executable alone. Nothing downloads an "
                "assembler on your behalf, and no read is assembled by an emulated engine."),
        automatic=False,
        command=("python studio_packaging/stage_bio_tools.py "
                 f"--platform {_staging_platform()} --require-skesa"))
    common = {"key": "assembly_runtime", "title": "SKESA read assembler",
              "capability": "Read assembly (FASTQ to contigs)",
              "blocks": ("raw FASTQ reads cannot be assembled here, so they can be quality-checked "
                         "but never typed; assemblies you import elsewhere are unaffected"),
              "required": True, "action": action}
    try:
        executable = resolve_skesa()
    except AssemblyError as error:
        return Requirement(**common, state="missing", detail={},
                           probe="Looked for the assembler beside the application, then on PATH.",
                           reason=str(error))
    return Requirement(**common, state="ready", detail={"path": str(executable)},
                       probe=f"Resolved the assembler at {executable}.",
                       reason=(f"Native SKESA is installed at {executable}. It was not run, so a "
                               "binary that cannot start on this machine is not ruled out."))


def installed_scheme_entries(scheme_paths=None, *, data_root=None, cancelled=None) -> list[dict]:
    """One row per installed scheme folder, classified as MLST or cgMLST by its loci."""
    from .reference_index import filter_scheme_locations, scheme_entries
    paths = (filter_scheme_locations(scheme_paths) if scheme_paths is not None
             else default_scheme_paths(data_root))
    return scheme_entries(paths, cancelled=cancelled)


def _scheme_detail(entries, kind) -> dict:
    bundled = resource_root() / "schemes"
    selected = [entry for entry in entries if entry.get("kind") == kind]
    inside = [entry for entry in selected
              if Path(entry["path"]).is_relative_to(bundled)]
    labelled = sorted({entry["label"] for entry in selected
                       if entry.get("genus") and entry.get("label")})
    return {"count": len(selected), "shipped_with_the_application": len(inside),
            "installed_by_you": len(selected) - len(inside),
            "named_organisms": len(labelled), "examples": labelled[:8],
            "roots": sorted({str(Path(entry["path"]).parent) for entry in entries})}


def mlst_scheme_requirement(entries) -> Requirement:
    """Classical 7-locus MLST. A separate capability from cgMLST, and reported so."""
    from .reference_index import CGMLST_LOCUS_THRESHOLD
    detail = _scheme_detail(entries, "mlst")
    action = Action(
        label="Download a classical MLST scheme",
        ui_route="Schemes > Download reference > PubMLST",
        detail=("Fetches one organism's 7-locus allele panel and profile table from PubMLST "
                "into your scheme library. Only the scheme is downloaded; no sample data is sent."),
        automatic=True,
        command="python studio_scripts/stage_schemes.py --download")
    common = {"key": "mlst_schemes", "title": "Classical MLST scheme library",
              "capability": "Classical 7-locus MLST typing",
              "blocks": "the ST menu has no 7-locus scheme to type an assembly against",
              "required": True, "action": action, "detail": detail}
    probe = (f"Counted allele files in every scheme folder; {CGMLST_LOCUS_THRESHOLD} loci or fewer "
             "and the folder is a classical scheme.")
    if not detail["count"]:
        return Requirement(**common, state="missing", probe=probe,
                           reason=("No classical MLST scheme is installed, so nothing can be typed "
                                   "by 7-locus MLST. This is separate from cgMLST: installing a "
                                   "cgMLST scheme does not give you an ST."))
    return Requirement(
        **common, state="ready", probe=probe,
        reason=(f"{detail['count']} classical schemes are installed "
                f"({detail['shipped_with_the_application']} shipped with the application, "
                f"{detail['installed_by_you']} added here), covering "
                f"{detail['named_organisms']} named organisms."))


def cgmlst_scheme_requirement(entries) -> Requirement:
    """Core-genome MLST. Never merged with the classical library: they are not one thing."""
    from .reference_index import CGMLST_LOCUS_THRESHOLD
    detail = _scheme_detail(entries, "cgmlst")
    action = Action(
        label="Download a cgMLST scheme for your organism",
        ui_route="Schemes > Download reference > cgMLST.org",
        detail=("Fetches one organism's core-genome allele panel. These are large downloads "
                "(hundreds of megabytes) and each provider's terms apply; review them first. "
                "cgMLST calling also needs the BLAST+ tools."),
        automatic=True)
    common = {"key": "cgmlst_schemes", "title": "cgMLST / wgMLST scheme library",
              "capability": "cgMLST / wgMLST typing",
              "blocks": ("the cgMLST menu has nothing to type against; classical 7-locus MLST is "
                         "unaffected and keeps working"),
              # Which organism's core genome you need is your decision, so no
              # cgMLST scheme ships with the application and its absence is not a
              # broken installation.
              "required": False, "action": action, "detail": detail}
    probe = (f"Counted allele files in every scheme folder; more than {CGMLST_LOCUS_THRESHOLD} "
             "loci, or a scheme that declares itself cgMLST, is counted here.")
    if not detail["count"]:
        return Requirement(**common, state="missing", probe=probe,
                           reason=("No cgMLST or wgMLST scheme is installed. Core-genome typing "
                                   "and the cgMLST minimum spanning tree have no scheme to use."))
    return Requirement(**common, state="ready", probe=probe,
                       reason=(f"{detail['count']} core-genome schemes are installed, covering "
                               f"{detail['named_organisms']} named organisms."))


def blast_requirement(tools=None, *, verify=True) -> Requirement:
    """Probe the native BLAST+ executables that cgMLST inference and HYDRA both need."""
    from .hydra_runtime import TOOLS, resolve_tool, tool_directories
    resolved = dict(tools) if tools is not None else {name: resolve_tool(name) for name in TOOLS}
    searched = [str(directory) for directory in tool_directories()]
    missing = sorted(name for name in TOOLS if not resolved.get(name))
    found = {name: resolved[name] for name in TOOLS if resolved.get(name)}
    action = Action(
        label="Restore the BLAST+ tools that ship beside the application",
        ui_route="Settings > Runtime tools",
        detail=("The portable build carries these in Tools/blast/bin beside the executable; "
                "extract the complete folder rather than the executable alone. In a source "
                "checkout, stage them or put an installed NCBI BLAST+ on PATH."),
        # Nothing in the application downloads BLAST+ on a user's behalf.
        automatic=False,
        command=("python studio_packaging/stage_bio_tools.py "
                 f"--platform {_staging_platform()}"))
    common = {"key": "blast_tools", "title": "NCBI BLAST+ tools",
              "capability": "cgMLST typing and HYDRA screening",
              "blocks": ("cgMLST full-CDS calling and every HYDRA run stop before they start, "
                         "because both shell out to these executables"),
              "required": True, "action": action}
    detail = {"expected": list(TOOLS), "found": found, "missing": missing,
              "searched_directories": searched,
              "path_searched": bool(os.environ.get("PATH"))}
    if missing:
        state = "partial" if found else "missing"
        return Requirement(
            **common, state=state, detail=detail,
            probe="Resolved each executable in the staged tool folders, then on PATH.",
            reason=(f"{len(missing)} of {len(TOOLS)} BLAST+ executables were not found "
                    f"({', '.join(missing)}). Looked in: {'; '.join(searched)}; and on PATH."))
    if not verify:
        return Requirement(
            **common, state="ready", detail=detail,
            probe="Resolved each executable in the staged tool folders, then on PATH.",
            reason=(f"All {len(TOOLS)} BLAST+ executables were found. They were not run, so a "
                    "binary that cannot start on this machine would not have been detected."))
    try:
        completed = subprocess.run(
            [resolved["blastn"], "-version"], capture_output=True, text=True,
            timeout=TOOL_TIMEOUT, check=False,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    except (OSError, subprocess.SubprocessError) as error:
        detail["execution_error"] = str(error)
        return Requirement(
            **common, state="unusable", detail=detail,
            probe=f"Ran {resolved['blastn']} -version.",
            reason=(f"The BLAST+ executables are present but blastn could not be run: {error}. "
                    "A partially extracted folder or a missing system library does this."))
    if completed.returncode:
        detail["execution_error"] = (completed.stderr or completed.stdout or "").strip()[:500]
        return Requirement(
            **common, state="unusable", detail=detail,
            probe=f"Ran {resolved['blastn']} -version.",
            reason=(f"blastn exited with code {completed.returncode} instead of reporting its "
                    f"version: {detail['execution_error'] or 'no output'}"))
    printed = (completed.stdout or "").strip().splitlines()
    detail["blastn_version"] = printed[0] if printed else ""
    return Requirement(
        **common, state="ready", detail=detail,
        probe=f"Ran {resolved['blastn']} -version.",
        reason=(f"All {len(TOOLS)} BLAST+ executables were found and blastn runs on this machine"
                + (f" ({detail['blastn_version']})." if detail["blastn_version"] else ".")))


def hydra_worker_path() -> Path | None:
    """The separate worker executable a frozen build launches; None in a checkout.

    Mirrors hydra_runtime._worker_command. Running the application from inside a
    zip viewer, or copying only the .exe out of the folder, leaves this behind and
    every HYDRA run then fails at launch instead of at configuration time.
    """
    if not getattr(sys, "frozen", False):
        return None
    suffix = ".exe" if os.name == "nt" else ""
    return Path(sys.executable).resolve().parent / ("WMLSTudio-HYDRA" + suffix)


def hydra_engine_requirement(capabilities=None) -> Requirement:
    """Probe the pinned HYDRA engine itself, separately from its reference data."""
    from .hydra_runtime import HYDRA_COMMIT, HYDRA_VERSION, runtime_capabilities
    capabilities = capabilities if capabilities is not None else runtime_capabilities()
    version = capabilities.get("hydra_version")
    worker = hydra_worker_path()
    action = Action(
        label="Reinstall the application with its pinned HYDRA engine",
        ui_route="",
        detail=(f"HYDRA {HYDRA_VERSION} (commit {HYDRA_COMMIT[:12]}) is a pinned dependency of "
                "this application, not a separate download. A build that lacks it is incomplete; "
                "re-extract the portable folder, or re-sync the source checkout."),
        automatic=False,
        command="uv sync")
    common = {"key": "hydra_engine", "title": "HYDRA analysis engine",
              "capability": "HYDRA resistance screening",
              "blocks": ("the HYDRA menu cannot screen for resistance determinants at all, "
                         "whatever databases are installed"),
              "required": True, "action": action}
    detail = {"expected_version": HYDRA_VERSION, "installed_version": version or "",
              "source_commit": HYDRA_COMMIT, "worker": str(worker or ""),
              "distribution_origin": capabilities.get("distribution_origin") or {}}
    probe = "Read the installed hydra-amr distribution metadata."
    if worker is not None:
        probe += f" Looked for the worker executable at {worker}."
    if not version:
        return Requirement(**common, state="missing", probe=probe, detail=detail,
                           reason=("The pinned HYDRA engine is not installed in this Python "
                                   "environment, so no resistance screening can run."))
    if version != HYDRA_VERSION:
        return Requirement(
            **common, state="unusable", probe=probe, detail=detail,
            reason=(f"HYDRA {version} is installed but this application is pinned to "
                    f"{HYDRA_VERSION}. It refuses to run an engine it was not validated "
                    "against rather than report results from an unknown version."))
    if worker is not None and not worker.is_file():
        extract = Action(
            label="Extract the complete WMLSTudio folder and run it from there",
            detail=("HYDRA runs in a separate worker program that must sit beside the "
                    f"application: {worker}. Copying the application out of the downloaded "
                    "folder, or running it from inside the ZIP, leaves the worker behind."),
            automatic=False)
        return Requirement(
            **{**common, "action": extract}, state="unusable", probe=probe, detail=detail,
            reason=(f"HYDRA {version} is installed, but its worker program is missing from "
                    f"{worker.parent}, so every run would fail the moment it started."))
    return Requirement(**common, state="ready", probe=probe, detail=detail,
                       reason=f"HYDRA {version} is installed, as pinned.")


def hydra_database_requirement(data_root=None, selected=None, *, cancelled=None) -> Requirement:
    """Probe the AMR database store HYDRA would actually use, and name the alternatives.

    A store that exists but holds no database is the failure the packaged build
    hits: HYDRA is present, BLAST is present, and every run still stops with
    "Select at least one installed HYDRA database". When another candidate store
    would work, this says so and points the fix at it instead of at a download.
    """
    from .hydra_runtime import database_status, installed_databases
    root = Path(data_root) if data_root is not None else default_data_root()
    download_root = root / "references" / "hydra"
    probes = []
    for origin, path in hydra_database_candidates(root, selected):
        try:
            entries = installed_databases(path)
            error = ""
        except (OSError, ValueError) as failure:  # HydraRuntimeError is a ValueError.
            entries, error = {}, str(failure)
        probes.append({"origin": origin, "path": str(path), "exists": path.is_dir(),
                       "databases": sorted(entries), "error": error, "entries": entries})
    download = Action(
        label="Download the HYDRA reference databases",
        ui_route="HYDRA > AMR databases / updates",
        detail=(f"Downloads the {' and '.join(HYDRA_DEFAULT_DATABASES)} reference sets into "
                f"{download_root} and keeps any previous snapshot. This needs the BLAST+ tools, "
                "which build the search indexes."),
        automatic=True,
        entry_point="wmlstudio.provisioning:install_hydra_databases",
        argument=str(download_root))
    update = Action(
        label="Update the HYDRA reference databases to the current NCBI release",
        ui_route="Updates > AMR reference databases",
        detail=("Fetches the release NCBI publishes today, keeps the snapshot you have, and "
                "records the new one for this project in one action. Analyses already recorded "
                "keep the reference snapshot they were run against."),
        automatic=True,
        entry_point="wmlstudio.provisioning:update_hydra_databases",
        argument=str(download_root))
    common = {"key": "hydra_database", "title": "HYDRA AMR database store",
              "capability": "HYDRA resistance screening",
              "blocks": ("every HYDRA run stops with \"Select at least one installed HYDRA "
                         "database\", because the engine has no reference data to search"),
              "required": True}
    active = probes[0]
    others = [probe for probe in probes[1:] if probe["databases"]]
    detail = {"active": {key: active[key] for key in ("origin", "path", "databases", "error")},
              "candidates": [{key: probe[key] for key in ("origin", "path", "exists",
                                                          "databases", "error")}
                             for probe in probes]}
    probe_sentence = (f"Read the database manifest of {len(probes)} candidate store"
                      f"{'s' if len(probes) != 1 else ''}, starting with {active['path']}.")
    if active["error"]:
        return Requirement(
            **common, state="unusable", probe=probe_sentence, detail=detail, action=download,
            reason=(f"The database store this project would use ({active['path']}, "
                    f"{active['origin']}) has a manifest that cannot be used: {active['error']}"))
    if not active["databases"]:
        if others:
            usable = others[0]
            return Requirement(
                **common, state="missing", probe=probe_sentence, detail=detail,
                action=Action(
                    label="Point HYDRA at the database store that is already installed",
                    ui_route="HYDRA > AMR databases / updates > Use selected snapshot",
                    detail=(f"{usable['path']} already holds "
                            f"{', '.join(usable['databases'])}. Selecting it records that folder "
                            "for this project; no files are moved or downloaded."),
                    automatic=True,
                    entry_point="wmlstudio.provisioning:select_hydra_database",
                    argument=usable["path"]),
                reason=(f"The store this project would use ({active['path']}, {active['origin']}) "
                        "holds no database, so every HYDRA run refuses to start. A usable store "
                        f"is already present at {usable['path']} ({usable['origin']}) holding "
                        f"{', '.join(usable['databases'])}; it is simply not the one selected."))
        return Requirement(
            **common, state="missing", probe=probe_sentence, detail=detail, action=download,
            reason=(f"No AMR database is installed in any store this application looks at "
                    f"({'; '.join(item['path'] for item in probes)}), so HYDRA has nothing to "
                    "search and every run refuses to start."))
    entries = active["entries"]
    # The release, its age and the day it was staged are read from the store the
    # run would actually use, so a report can say how stale the evidence is rather
    # than only that some reference data exists.
    status = database_status(active["path"])
    versions = {name: str(entry.get("version") or entry.get("installed") or "unrecorded")
                for name, entry in sorted(entries.items())}
    detail.update({"databases": versions, "organisms": status["organisms"],
                   "point_mutation_organisms": status["point_mutation_organisms"],
                   "root": active["path"], "origin": active["origin"],
                   "release": status["release"], "released_utc": status["released_utc"],
                   "staged": status["staged"], "age_days": status["age_days"],
                   "stale": status["stale"], "bundled": status["bundled"],
                   "label": status["label"],
                   "purposes": {name: item["purpose"]
                                for name, item in status["installed"].items()}})
    absent = [name for name in HYDRA_DEFAULT_DATABASES if name not in entries]
    if absent:
        return Requirement(
            **common, state="partial", probe=probe_sentence, detail=detail, action=download,
            reason=(f"{active['path']} holds {', '.join(sorted(entries))} but not "
                    f"{', '.join(absent)}. HYDRA will run and report what it can search for, and "
                    "will silently report nothing for the element classes it cannot"
                    + (", including every point mutation." if "protein" in absent else ".")))
    age = (" Determinants named after that release are not in it; update before reporting."
           if status["stale"] else "")
    return Requirement(
        **common, state="ready", probe=probe_sentence, detail=detail,
        action=update if status["stale"] else download,
        # The label carries the release, where the copy came from and how old it is,
        # because "a database is installed" does not tell anyone what it can find.
        reason=(f"{status['label']} Store: {active['path']} ({active['origin']})."
                + (f" {len(status['organisms'])} organisms are accepted for point-mutation "
                   f"catalogs, {len(status['point_mutation_organisms'])} of them with a DNA "
                   "catalogue." if status["organisms"] else "") + age))


def hydra_prerequisites(data_root=None, selected=None, *, organism=None, databases=None,
                        verify=False, cancelled=None) -> dict:
    """The gate a Run HYDRA button calls before it starts anything at all.

    The three HYDRA requirements plus the engine's own pre-run check, collapsed
    into one answer: can this run start, and if not, which piece is missing, what
    is it for, and what is the one action that fixes it. `warnings` is what the run
    would still not be able to report — an absent protein set, an organism with no
    catalogue, a reference release that predates a determinant — and it is shown
    before the run, not discovered in the result afterwards.

    verify defaults to False here because this runs on a button press: the tools
    are resolved but not executed. Pass verify=True for a settings page.
    """
    from .hydra_runtime import preflight, runtime_capabilities
    root = Path(data_root) if data_root is not None else default_data_root()
    capabilities = runtime_capabilities()
    items = [hydra_engine_requirement(capabilities),
             blast_requirement(capabilities.get("tools"), verify=verify),
             hydra_database_requirement(root, selected, cancelled=cancelled)]
    # The store the run would really use: the one the requirement probed first,
    # so the refusal names that folder rather than a folder nobody chose.
    detail = items[-1].detail
    store = detail.get("root") or (detail.get("active") or {}).get("path") or ""
    checked = preflight(store or None, databases=databases, organism=organism,
                        capabilities=capabilities)
    blocking = [item for item in items if not item.ready]
    action = next((item.action for item in blocking if item.action is not None), None)
    message = checked["message"]
    if blocking and not message:
        # A requirement can be unmet in a way the engine's own probe cannot see —
        # a frozen build whose separate worker program is absent, for instance.
        message = ("HYDRA cannot start: " + blocking[0].title + " is "
                   + STATE_WORDS[blocking[0].state] + ". " + blocking[0].reason)
    return {"ready": not blocking and checked["ready"], "message": message,
            "missing": checked["missing"], "warnings": checked["warnings"],
            "items": [item.as_dict() for item in items],
            "blocking": [item.key for item in blocking],
            "database": checked["database_status"], "databases": checked["databases"],
            "organism": checked["organism"], "limitations": checked["limitations"],
            "action": ({"key": next(item.key for item in blocking if item.action is action),
                        **asdict(action)} if action is not None else None)}


def characterization_requirement(reference_root=None, *, verify=True,
                                 cancelled=None) -> Requirement:
    """Probe the characterization snapshot and report its manifest format version.

    Format 1 carries the species and virulence panels only. The organism-specific
    assays were added in format 2, so a format 1 snapshot is a working but reduced
    installation, and saying "ready" about it would be false.
    """
    # _section_present is the runtime's own section test, imported rather than
    # copied so this report can never disagree with what organism_modules will
    # actually refuse to run.
    from .characterization_refs import bundled_reference_root
    from .organism_modules import _section_present, registered_modules
    root = reference_root if reference_root is not None else bundled_reference_root()
    action = Action(
        label="Install an updated characterization reference snapshot",
        ui_route="Characterization > Install / update reference panel",
        detail=("Stages the pinned Kleborate, sccmec and Kaptive reference panels, each with its "
                "own upstream licence, into the application. None of these is Kleborate, Kaptive "
                "or SCCmecFinder; they are screens over those projects' published reference data."),
        automatic=True,
        command="python studio_scripts/stage_characterization.py")
    common = {"key": "characterization", "title": "Characterization reference snapshot",
              "capability": "Organism-specific assays",
              "blocks": ("the Klebsiella-complex refinement of organism identification and every "
                         "organism-specific assay (SCCmec, Klebsiella locus STs, capsule markers) "
                         "are unavailable"),
              # Typing, filing and HYDRA do not depend on it, so a missing
              # snapshot narrows the application rather than blocking it.
              "required": False, "action": action}
    if root is None:
        return Requirement(**common, state="missing",
                           probe=("Looked for a staged characterization snapshot inside "
                                  "the application."),
                           reason=("No characterization reference snapshot is installed, so no "
                                   "organism-specific assay can run and Klebsiella isolates are "
                                   "identified by the broad panel alone."),
                           detail={"path": ""})
    try:
        manifest, probe = _snapshot_manifest(root, verify=verify, cancelled=cancelled)
    except (OSError, ValueError) as error:
        return Requirement(**common, state="unusable", probe=f"Validated {root}.",
                           reason=(f"A characterization snapshot is installed at {root} but it did "
                                   f"not pass validation, so it is not used: {error}"),
                           detail={"path": str(root)})
    modules = []
    for key, module in registered_modules().items():
        absent = [section for section in module.manifest_sections
                  if not _section_present(manifest, section)]
        modules.append({"key": key, "title": module.title, "available": not absent,
                        "missing_sections": absent,
                        "taxa": sorted(module.match.genera)})
    detail = {"path": str(root), "manifest_format_version": manifest.get("format_version"),
              "reference_digest": manifest.get("reference_digest", ""),
              "species_count": len(manifest.get("species", [])),
              "virulence_loci": sorted(manifest.get("virulence", {})),
              "modules": modules}
    unavailable = [item for item in modules if not item["available"]]
    if unavailable:
        names = ", ".join(item["title"] for item in unavailable)
        sections = sorted({section for item in unavailable for section in item["missing_sections"]})
        return Requirement(
            **common, state="partial", probe=f"Read {root}: {probe}.", detail=detail,
            reason=(f"The installed snapshot is manifest format "
                    f"{manifest.get('format_version')} and carries no {', '.join(sections)} "
                    f"section, so these assays cannot run: {names}. Species comparison and "
                    "virulence loci from this snapshot still work."))
    return Requirement(
        **common, state="ready", probe=f"Read {root}: {probe}.", detail=detail,
        reason=(f"Manifest format {manifest.get('format_version')} is installed, covering "
                f"{detail['species_count']} species references and {len(modules)} "
                "organism-specific assays."))


def requirements(*, data_root=None, scheme_paths=None, hydra_database_root=None,
                 characterization_root=None, verify=True, cancelled=None) -> list[Requirement]:
    """Probe every prerequisite once, in the order a person meets them."""
    from .hydra_runtime import runtime_capabilities
    root = Path(data_root) if data_root is not None else default_data_root()
    # One call resolves the engine version and all five executables; probing them
    # separately would ask the filesystem the same question twice.
    capabilities = runtime_capabilities()
    entries = installed_scheme_entries(scheme_paths, data_root=root, cancelled=cancelled)
    items = [species_panel_requirement(root, verify=verify, cancelled=cancelled),
             assembly_runtime_requirement(),
             mlst_scheme_requirement(entries),
             cgmlst_scheme_requirement(entries),
             blast_requirement(capabilities.get("tools"), verify=verify),
             hydra_engine_requirement(capabilities),
             hydra_database_requirement(root, hydra_database_root, cancelled=cancelled),
             characterization_requirement(characterization_root, verify=verify,
                                          cancelled=cancelled)]
    return sorted(items, key=lambda item: ORDER.index(item.key))


def _lower_first(text: str) -> str:
    """Lower the first letter, unless that would damage a name like HYDRA or NCBI."""
    if len(text) > 1 and text[1].isupper():
        return text
    return text[:1].lower() + text[1:]


def capability_states(items) -> dict:
    """Per-menu readiness: what a user can actually do right now, and what blocks it."""
    known = {item.key: item for item in items}
    states = {}
    for key, (title, needs) in CAPABILITIES.items():
        blocked = [name for name in needs if name in known and not known[name].ready]
        if blocked:
            first = known[blocked[0]]
            label = first.action.label if first.action else "install the missing reference data"
            sentence = (f"{title} is not available yet ({first.title}: "
                        f"{STATE_WORDS[first.state]}). Next: {label}.")
        else:
            sentence = f"{title} is ready to use."
        states[key] = {"title": title, "ready": not blocked, "blocked_by": blocked,
                       "sentence": sentence}
    return states


def summarize(items) -> str:
    """One sentence: can I run an investigation now, and if not what do I do first."""
    blocking = [item for item in items if item.required and not item.ready]
    optional = [item for item in items if not item.required and not item.ready]
    if blocking:
        first = blocking[0]
        action = first.action
        where = f" ({action.ui_route})" if action and action.ui_route else ""
        step = _lower_first(action.label) if action else "install the missing reference data"
        return (f"Not yet: {first.blocks}; the single next thing to do is {step}{where}.")
    if optional:
        first = optional[0]
        step = _lower_first(first.action.label) if first.action else "install it"
        return ("Yes, you can import samples, file them into organism folders and run an "
                f"investigation now, although {_lower_first(first.capability)} stays unavailable "
                f"until you {step}.")
    return ("Yes, everything this application needs is installed, so you can import samples, "
            "let them file themselves into organism folders, and run an investigation now.")


def report(*, data_root=None, scheme_paths=None, hydra_database_root=None,
           characterization_root=None, verify=True, cancelled=None) -> dict:
    """The whole provisioning picture as plain data: states, reasons, and the next action."""
    root = Path(data_root) if data_root is not None else default_data_root()
    paths = list(scheme_paths) if scheme_paths is not None else default_scheme_paths(root)
    items = requirements(data_root=root, scheme_paths=paths,
                         hydra_database_root=hydra_database_root,
                         characterization_root=characterization_root,
                         verify=verify, cancelled=cancelled)
    blocking = [item for item in items if item.required and not item.ready]
    following = blocking or [item for item in items if not item.ready]
    next_item = following[0] if following else None
    next_action = None
    if next_item is not None and next_item.action is not None:
        next_action = {"key": next_item.key, "title": next_item.title,
                       **asdict(next_item.action)}
    return {"format_version": REPORT_FORMAT_VERSION,
            "application_version": __version__,
            "generated_utc": datetime.now(UTC).isoformat(),
            "data_root": str(root), "scheme_paths": [str(path) for path in paths],
            "verified": bool(verify),
            "items": [item.as_dict() for item in items],
            "ready": not blocking,
            "blocking": [item.key for item in blocking],
            "summary": summarize(items),
            "next_action": next_action,
            "capabilities": capability_states(items)}


def install_species_panel(data_root=None, *, cancelled=None, progress=None) -> dict:
    """Download and verify the broad ANI panel into the application's data folder.

    This is the action the species panel requirement names. It is only ever called
    because a person asked for it; nothing in analysis or startup reaches here.
    """
    from .organism_panel import PANEL_DIRECTORY, provision_species_panel
    root = Path(data_root) if data_root is not None else default_data_root()
    return provision_species_panel(root / PANEL_DIRECTORY, cancelled=cancelled, progress=progress)


def install_hydra_databases(download_root=None, names=HYDRA_DEFAULT_DATABASES, *,
                            cancelled=None, progress=None) -> dict:
    """Download the named AMR reference databases, keeping any existing snapshot.

    The upstream downloader publishes a new snapshot folder beside download_root
    and returns its path; the caller must record that path as the store to use
    (select_hydra_database), because existing analyses keep their own provenance.
    """
    from .hydra_runtime import update_databases
    root = (Path(download_root) if download_root is not None
            else default_data_root() / "references" / "hydra")
    return update_databases(root, list(names), cancelled=cancelled, progress=progress)


def hydra_update_available(data_root=None, selected=None, *, cancelled=None) -> dict:
    """Is there a newer AMR reference release than the one installed? One network call.

    The update menu's "check" action. It asks NCBI for the current release string
    and nothing else — no sample data leaves the machine, and nothing is
    downloaded, installed or changed by asking. A check that cannot reach NCBI
    says so and never reports "up to date" on a failed request.
    """
    from .hydra_runtime import database_status, latest_release
    root = Path(data_root) if data_root is not None else default_data_root()
    candidates = hydra_database_candidates(root, selected)
    status = database_status(candidates[0][1] if candidates else None)
    result = {"installed": status["release"], "installed_label": status["label"],
              "latest": "", "newer_available": False, "error": "",
              "checked_utc": datetime.now(UTC).isoformat(),
              "age_days": status["age_days"], "stale": status["stale"]}
    try:
        result["latest"] = latest_release(cancelled=cancelled)
    except (OSError, ValueError) as error:
        result["error"] = f"Could not reach NCBI to check for a newer release: {error}"
        result["message"] = (result["error"] + " The installed release is unchanged: "
                             + status["label"])
        return result
    result["newer_available"] = bool(status["release"]) and result["latest"] != status["release"]
    if not status["installed"]:
        result["message"] = (f"NCBI publishes release {result['latest']} today, and no reference "
                             "database is installed here yet.")
    elif result["newer_available"]:
        result["message"] = (f"NCBI publishes release {result['latest']} today; this installation "
                             f"has {status['release']}. Determinants named between the two are "
                             "not in what you have.")
    else:
        result["message"] = (f"The installed reference release {status['release']} is the one NCBI "
                             "publishes today.")
    return result


def update_hydra_databases(download_root=None, names=HYDRA_DEFAULT_DATABASES, *, project=None,
                           cancelled=None, progress=None) -> dict:
    """Download the current reference release and start using it, in one action.

    install_hydra_databases publishes a snapshot and stops; a person then has to
    find and select it, which is a second step nobody knew was required. This does
    both, so "update" means what it says. The previous snapshot is kept on disk and
    every analysis already recorded keeps the reference snapshot it was run
    against: updating changes what the next run can find, never what a past run
    reported.
    """
    result = install_hydra_databases(download_root, names, cancelled=cancelled, progress=progress)
    published = result["database_root"]
    selected = False
    if project is not None:
        select_hydra_database(project, published)
        selected = True
    from .hydra_runtime import database_status
    status = database_status(published)
    result.update({"selected": selected, "database_status": status,
                   "summary": (f"{status['label']} "
                               + ("This project now uses it; results already recorded keep the "
                                  "snapshot they were run against."
                                  if selected else
                                  "Select it for a project to run against it."))})
    return result


def select_hydra_database(project, path) -> str:
    """Record an installed store as the one this project uses. No files are moved.

    Refuses a folder with no readable database manifest, so the setting can never
    point at something that will fail later, in the middle of a run.
    """
    from .hydra_runtime import installed_databases
    root = Path(path).expanduser().resolve()
    if not installed_databases(root):
        raise ValueError(f"{root} holds no readable HYDRA database manifest; nothing was changed.")
    project.set_setting("hydra_database_root", str(root))
    return str(root)
