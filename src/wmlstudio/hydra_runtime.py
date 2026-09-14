"""Cancellable, native assembly execution of the pinned upstream HYDRA engine.

The wrapper does not reimplement any scientific calling logic. It runs HYDRA
in a dedicated child process, never auto-downloads databases during analysis,
and imports the engine's JSON with execution and reference provenance.
"""

from __future__ import annotations

import importlib.metadata
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request
import uuid
from datetime import UTC, datetime
from pathlib import Path

from wmlstudio.hydra import load_hydra_report
from wmlstudio.sequence import (
    AnalysisCancelled,
    SequenceReader,
    file_sha256,
    file_signature,
    sample_name,
)

HYDRA_COMMIT = "6d36c109491c16544e8919fe6962b4b62e97d3d7"
HYDRA_VERSION = "1.4.0"
TOOLS = ("blastn", "blastx", "blastp", "tblastn", "makeblastdb")
NCBI_REFERENCE_ROOT = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
# The stores a run uses when the caller names none, in the order a report reads.
DEFAULT_DATABASES = ("ncbi", "protein")
# What each store is for, said in the words a microbiologist reads. A run without
# one of these does not fail: it returns nothing for that whole class of evidence,
# which reads exactly like a negative result. Every refusal below names the store
# and this sentence, so "nothing was found" is never confused with "nothing ran".
DATABASE_PURPOSE = {
    "ncbi": ("acquired resistance, stress and virulence genes, searched with blastn against the "
             "NCBI AMRFinderPlus nucleotide catalogue"),
    "protein": ("the translated protein search and every organism point-mutation catalogue, "
                "searched with blastx against AMRProt"),
}
# What each executable is for, so a missing one is a sentence and not a filename.
TOOL_PURPOSE = {
    "blastn": "the nucleotide gene search",
    "blastx": "the translated protein search, which is where point mutations are read",
    "blastp": "protein-to-protein confirmation of a translated hit",
    "tblastn": "recovering a gene that the translated search split across a contig break",
    "makeblastdb": "building the search index a reference store needs before it can be read",
}
# A snapshot older than this still runs and is still reported as the exact release
# it is; the age is surfaced so a user can decide whether to update before
# reporting. Determinants named after this date are simply not in it.
REFERENCE_AGE_DAYS = 180


class HydraRuntimeError(ValueError):
    """A readable configuration, execution, or provenance failure."""


def tool_directories():
    roots = [Path(sys.executable).resolve().parent]
    if getattr(sys, "frozen", False):
        roots.insert(0, Path(sys._MEIPASS))
    directories = [root / relative for root in roots
                   for relative in ("Tools/blast/bin", "Tools/blast+/bin", "Tools/bin")]
    directories.append(Path(__file__).resolve().parent / "resources/tools/blast/bin")
    return directories


def bundled_database_root():
    """Return the explicitly staged starter snapshot, or None in a source-only checkout."""
    candidates = [Path(__file__).resolve().parent / "resources/hydra/starter"]
    if getattr(sys, "frozen", False):
        candidates.insert(0, Path(sys._MEIPASS) / "wmlstudio/resources/hydra/starter")
    return next((path for path in candidates if (path / "manifest.json").is_file()), None)


def resolve_tool(name, explicit=None):
    if name not in TOOLS:
        raise HydraRuntimeError(f"Unsupported native tool: {name}")
    filename = name + (".exe" if os.name == "nt" else "")
    if explicit is not None:
        path = Path(explicit).resolve()
        if not path.is_file():
            raise HydraRuntimeError(f"Tool was not found: {path}")
        return str(path)
    for directory in tool_directories():
        if (directory / filename).is_file():
            return str(directory / filename)
    return shutil.which(filename)


def installed_databases(db_root):
    """Read the upstream manifest without modifying or creating a database store."""
    root = Path(db_root).resolve()
    path = root / "manifest.json"
    if not path.is_file():
        return {}
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        records = manifest.get("databases", {})
        if not isinstance(records, dict):
            raise ValueError("databases must be an object")
        result = {}
        for name, entry in records.items():
            if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
                raise ValueError(f"invalid entry for {name}")
            # A store staged on Windows records "nucl\ncbi"; on POSIX that is one
            # filename, not a folder, and the whole store would read as empty.
            # Normalising here makes one staged snapshot portable, and tightens
            # the containment check rather than loosening it: "..\\outside" is
            # now recognised as an escape instead of a legal filename.
            target = (root / entry["path"].replace("\\", "/")).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"database path leaves the selected store: {name}")
            if target.exists():
                result[name] = dict(entry)
        return result
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise HydraRuntimeError(f"Cannot read HYDRA database manifest: {exc}") from exc


def _entry_directory(root, entry):
    """The folder one manifest entry names, with the same separator normalisation."""
    return (Path(root) / str(entry.get("path", "")).replace("\\", "/")).resolve()


def _release_date(version):
    """The NCBI release date inside a version like '2026-08-07.1', or None."""
    match = re.match(r"([0-9]{4}-[0-9]{2}-[0-9]{2})", str(version or ""))
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%Y-%m-%d").replace(tzinfo=UTC)
    except ValueError:
        return None


def organism_catalogue(db_root):
    """The organism names the installed store will actually accept, and why.

    Upstream resolves --organism against its own taxgroup table and its DNA
    mutation catalogues, so this reads the same two places rather than inventing a
    list. "accepted" is every name the engine will take; "point_mutations" is the
    smaller set that has a curated catalogue behind it. An organism outside the
    first list stops the run upstream; one inside the first but outside the second
    is accepted and simply has no mutations to report, which is not the same thing
    as having none.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    accepted, catalogued = set(), set()
    if root is None or not root.is_dir():
        return {"accepted": [], "point_mutations": [], "root": str(root or "")}
    try:
        entries = installed_databases(root)
    except HydraRuntimeError:
        entries = {}
    for name in (entries.get("protein", {}).get("organisms") or []):
        if isinstance(name, str) and name.strip():
            accepted.add(name.strip())
            catalogued.add(name.strip())
    mutation = root / "mutation" / "dna"
    if mutation.is_dir():
        for path in sorted(mutation.glob("*.fna")):
            accepted.add(path.stem)
            catalogued.add(path.stem)
    protein = entries.get("protein")
    taxgroup = (_entry_directory(root, protein) / "taxgroup.tsv") if protein else None
    if taxgroup is not None and taxgroup.is_file():
        try:
            with taxgroup.open(encoding="utf-8") as handle:
                handle.readline()
                for line in handle:
                    value = line.split("\t")[0].strip()
                    if value and not value.startswith("#"):
                        accepted.add(value)
        except OSError as exc:
            raise HydraRuntimeError(f"Cannot read the organism table: {exc}") from exc
    return {"accepted": sorted(accepted), "point_mutations": sorted(catalogued),
            "root": str(root)}


def match_organism(name, db_root=None, *, accepted=None):
    """Map an assigned organism onto a name the installed catalogues really use.

    Upstream's taxgroups mix species ("Klebsiella_pneumoniae") with genus-level
    groups ("Escherichia", "Salmonella"), and they are underscore-joined, so an
    assignment of "Escherichia coli" belongs to "Escherichia" and one of
    "Klebsiella pneumoniae" to "Klebsiella_pneumoniae". Nothing is resolved by
    similarity: a name with no catalogue returns None so the caller can say that
    plainly instead of substituting a neighbouring organism's mutation list.
    """
    text = " ".join(str(name or "").replace("_", " ").split())
    if not text:
        return None
    choices = list(accepted) if accepted is not None else organism_catalogue(db_root)["accepted"]
    index = {choice.replace("_", " ").casefold(): choice for choice in choices}
    return index.get(text.casefold()) or index.get(text.split(" ")[0].casefold())


def database_status(db_root=None, databases=None):
    """What a store holds, which release it is, how old that release is, what is absent.

    Plain data only: nothing is downloaded, created or repaired by asking. The
    release and its age are reported separately from the day the copy was staged,
    because it is the reference release — not the download — that decides which
    determinants the run could possibly name.
    """
    root = Path(db_root).resolve() if db_root is not None else None
    wanted = [str(name) for name in (databases if databases is not None else DEFAULT_DATABASES)]
    status = {"root": str(root or ""), "bundled": False, "installed": {}, "requested": wanted,
              "missing": [], "error": "", "release": "", "released_utc": "", "staged": "",
              "age_days": None, "stale": False, "label": "No AMR reference database is installed.",
              "organisms": [], "point_mutation_organisms": []}
    bundled = bundled_database_root()
    status["bundled"] = bool(root is not None and bundled is not None
                             and root == Path(bundled).resolve())
    if root is None:
        return status
    try:
        entries = installed_databases(root)
    except HydraRuntimeError as exc:
        status["error"] = str(exc)
        status["label"] = str(exc)
        return status
    for name, entry in sorted(entries.items()):
        status["installed"][name] = {
            "path": str(entry.get("path", "")), "kind": str(entry.get("kind", "")),
            "version": str(entry.get("version", "") or ""),
            "installed": str(entry.get("installed", "") or ""),
            "sequences": entry.get("sequences"),
            "title": str(entry.get("title", "") or ""),
            "purpose": DATABASE_PURPOSE.get(name, "an additional reference set this engine reads")}
    status["missing"] = [{"name": name, "purpose": DATABASE_PURPOSE.get(
        name, "an additional reference set this engine reads")} for name in wanted
        if name not in entries]
    if not entries:
        return status
    catalogue = organism_catalogue(root)
    status["organisms"] = catalogue["accepted"]
    status["point_mutation_organisms"] = catalogue["point_mutations"]
    versions = sorted({item["version"] for item in status["installed"].values() if item["version"]})
    status["release"] = versions[0] if len(versions) == 1 else ("; ".join(versions) or "unrecorded")
    staged = sorted({item["installed"] for item in status["installed"].values()
                     if item["installed"]})
    status["staged"] = staged[0] if staged else ""
    released = _release_date(versions[0]) if len(versions) == 1 else None
    if released is not None:
        status["released_utc"] = released.isoformat()
        status["age_days"] = max(0, (datetime.now(UTC) - released).days)
        status["stale"] = status["age_days"] > REFERENCE_AGE_DAYS
    where = "bundled with this release" if status["bundled"] else "installed"
    age = (f", {status['age_days']} days old" if status["age_days"] is not None else "")
    status["label"] = (f"NCBI AMRFinderPlus reference release {status['release']} "
                       f"({', '.join(sorted(entries))}), {where}"
                       + (f" on {status['staged'][:10]}" if status["staged"] else "") + age + ".")
    return status


def runtime_capabilities(db_root=None):
    try:
        distribution = importlib.metadata.distribution("hydra-amr")
        version = distribution.version
        origin = json.loads(distribution.read_text("direct_url.json") or "{}")
    except importlib.metadata.PackageNotFoundError:
        version = None
        origin = {}
    tools = {name: resolve_tool(name) for name in TOOLS}
    databases = installed_databases(db_root) if db_root is not None else {}
    available = version == HYDRA_VERSION and all(tools.values())
    missing = [name for name, path in tools.items() if not path]
    if not available:
        message = f"HYDRA {version or 'not installed'}; missing tools: {', '.join(missing) or 'none'}."
    elif db_root is None or not databases:
        # The engine being ready is not the same as the screen being runnable, and
        # reporting only the first was how "HYDRA does not work" stayed unexplained.
        message = ("The HYDRA engine and BLAST+ tools are ready, but no AMR reference database is "
                   "installed, so a run would have nothing to search.")
    else:
        message = ("Native HYDRA assembly runtime is ready with "
                   + ", ".join(sorted(databases)) + ".")
    return {"hydra_version": version, "expected_hydra_version": HYDRA_VERSION,
            "source_commit": HYDRA_COMMIT, "distribution_origin": origin,
            "tools": tools, "databases": databases,
            "assembly_available": available, "available": available, "message": message,
            "blastn": tools["blastn"], "makeblastdb": tools["makeblastdb"],
            "reads_available": False,
            "limitations": ["Assembly-only HYDRA execution; its direct-read and minority-allele pipelines are not enabled. Paired-read assembly is handled separately by native SKESA.",
                            "Gene and mutation evidence is not a validated susceptibility phenotype."]}


def preflight(db_root=None, *, databases=None, organism=None, point_mutations=True,
              capabilities=None):
    """Everything a run needs, checked before a single sequence file is opened.

    The failure this exists to stop is the quiet one: HYDRA is installed, BLAST is
    installed, no reference store is, and the run either refuses with a sentence
    about selecting a database or — worse — completes and reports nothing, which a
    reader cannot tell from a clean isolate. So every piece is checked first, each
    missing one is named together with what it is for and what the run loses
    without it, and `ready` is False before any work starts.

    Nothing here downloads, creates or repairs anything, and `warnings` is kept
    separate from `missing`: a warning narrows what the run can report and is
    recorded with the result, a missing item stops it.
    """
    capabilities = capabilities if capabilities is not None else runtime_capabilities(db_root)
    tools = capabilities.get("tools") or {}
    version = capabilities.get("hydra_version")
    missing, warnings = [], []
    if version != HYDRA_VERSION:
        missing.append({"kind": "engine", "name": f"HYDRA {HYDRA_VERSION}",
                        "purpose": "the analysis engine itself",
                        "reason": (f"{version or 'No HYDRA engine'} is installed and this "
                                   f"application is pinned to {HYDRA_VERSION}; it will not report "
                                   "results from an engine it was not validated against.")})
    for name in TOOLS:
        if not tools.get(name):
            missing.append({"kind": "tool", "name": name, "purpose": TOOL_PURPOSE[name],
                            "reason": (f"{name} was not found beside the application or on PATH, "
                                       "so this search cannot be run at all.")})
    status = database_status(db_root, databases)
    if db_root is None:
        missing.append({"kind": "database", "name": "AMR reference database",
                        "purpose": "every gene and mutation this screen can name",
                        "reason": "No database store has been chosen for this project."})
    elif status["error"]:
        missing.append({"kind": "database", "name": "AMR reference database",
                        "purpose": "every gene and mutation this screen can name",
                        "reason": status["error"]})
    chosen = [name for name in status["requested"] if name in status["installed"]]
    if db_root is not None and not status["error"]:
        if not status["installed"]:
            missing.append({"kind": "database", "name": "AMR reference database",
                            "purpose": "every gene and mutation this screen can name",
                            "reason": (f"{status['root']} holds no reference database, so the "
                                       "screen has nothing to search and would report no "
                                       "determinants for every isolate.")})
        elif databases is not None:
            # A database the caller asked for by name is not optional: running
            # without it would answer a different question than the one asked.
            for entry in status["missing"]:
                missing.append({"kind": "database", **entry,
                                "reason": (f"The database '{entry['name']}' is not installed in "
                                           f"{status['root']}.")})
        else:
            for entry in status["missing"]:
                warnings.append(
                    f"The '{entry['name']}' reference set is not installed, so {entry['purpose']} "
                    "did not run. Nothing was reported for it; that is not a negative result.")
    resolved, organism_reason = None, ""
    if organism:
        accepted = status["organisms"]
        resolved = match_organism(organism, accepted=accepted) if accepted else str(organism)
        if accepted and resolved is None:
            organism_reason = (
                f"The installed reference release has no catalogue for '{organism}', so point "
                "mutations were not assessed for this isolate. Genes were still searched for. "
                "An absent catalogue is not an absence of mutations.")
            warnings.append(organism_reason)
        elif resolved and accepted and resolved not in status["point_mutation_organisms"]:
            organism_reason = (
                f"'{resolved}' is accepted by the installed release but has no DNA point-mutation "
                "catalogue in it, so only protein-level mutations could be reported.")
            warnings.append(organism_reason)
    elif point_mutations:
        organism_reason = ("No organism was given, so no point-mutation catalogue was selected. "
                           "Assign a genus and species to this isolate to have them assessed.")
        warnings.append(organism_reason)
    if point_mutations and "protein" not in chosen and status["installed"]:
        warnings.append("Point mutations were requested but the protein reference set is not "
                        "among the databases being searched, so none can be reported.")
    if status["stale"]:
        warnings.append(f"This reference release is {status['age_days']} days old "
                        f"({status['release']}). Determinants named after it are not in it; use "
                        "the update action to fetch the current NCBI release.")
    message = ""
    if missing:
        lines = [f"  • {item['name']} — needed for {item['purpose']}. {item['reason']}"
                 for item in missing]
        fix = ("Install the BLAST+ tools that ship beside the application"
               if any(item["kind"] == "tool" for item in missing) else
               "Download the reference databases from HYDRA > AMR databases / updates")
        message = ("HYDRA cannot start, so nothing was run and no isolate was marked screened.\n"
                   + "\n".join(lines) + f"\n{fix}. Nothing is ever downloaded during a run.")
    return {"ready": not missing, "message": message, "missing": missing, "warnings": warnings,
            "databases": chosen, "database_status": status,
            "organism": {"requested": str(organism or ""), "resolved": resolved or "",
                         "reason": organism_reason,
                         "point_mutations": bool(point_mutations)},
            "engine": {"installed": version or "", "expected": HYDRA_VERSION},
            "tools": dict(tools), "limitations": list(capabilities.get("limitations") or ())}


def _worker_command(arguments):
    if getattr(sys, "frozen", False):
        suffix = ".exe" if os.name == "nt" else ""
        worker = Path(sys.executable).resolve().parent / ("WMLSTudio-HYDRA" + suffix)
        if not worker.is_file():
            raise HydraRuntimeError("The portable HYDRA worker is missing. Extract the complete WMLSTudio folder.")
        return [str(worker), "--upstream", *arguments]
    return [sys.executable, "-m", "wmlstudio.hydra_runtime", "--upstream", *arguments]


def _child_environment():
    environment = dict(os.environ)
    directories = [str(directory) for directory in tool_directories() if directory.is_dir()]
    environment["PATH"] = os.pathsep.join(directories + [environment.get("PATH", "")])
    environment["PYTHONUTF8"] = "1"
    environment["HYDRA_NO_BANNER"] = "1"
    # BLAST receives its explicit per-sample --threads allocation. Prevent
    # numerical helper libraries from creating an additional machine-wide pool
    # in every concurrently running HYDRA worker.
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[key] = "1"
    # Frozen Linux applications prepend their own shared libraries. BLAST must
    # resolve its native dependencies from its own runtime, not Qt's libraries.
    if "LD_LIBRARY_PATH_ORIG" in environment:
        environment["LD_LIBRARY_PATH"] = environment["LD_LIBRARY_PATH_ORIG"]
    elif getattr(sys, "frozen", False):
        environment.pop("LD_LIBRARY_PATH", None)
    return environment


def _stop_process(process):
    if process.poll() is not None:
        return
    if os.name == "nt":
        # PID belongs to the child we launched, and /T includes its BLAST children.
        subprocess.run(["taskkill.exe", "/PID", str(process.pid), "/T", "/F"],
                       capture_output=True, timeout=15, check=False,
                       creationflags=subprocess.CREATE_NO_WINDOW)
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name != "nt":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.kill()
        process.wait(timeout=5)


def _run_child(arguments, directory, cancelled=None, progress=None):
    if cancelled and cancelled():
        raise AnalysisCancelled("HYDRA cancelled before execution.")
    command = _worker_command(arguments)
    log_path = Path(directory) / "hydra-execution.log"
    settings = ({"creationflags": subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP}
                if os.name == "nt" else {"start_new_session": True})
    previous_size = 0
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   cwd=directory, env=_child_environment(), **settings)
        try:
            while process.poll() is None:
                if cancelled and cancelled():
                    raise AnalysisCancelled("HYDRA cancelled; partial results were not imported.")
                size = log_path.stat().st_size
                if progress and size != previous_size:
                    with log_path.open("rb") as reader:
                        reader.seek(max(0, size - 2048))
                        tail = reader.read().decode("utf-8", errors="replace").splitlines()
                    if tail:
                        progress(0, 0, tail[-1][-300:])
                    previous_size = size
                time.sleep(0.1)
        except BaseException:
            _stop_process(process)
            raise
    with log_path.open("rb") as log:
        log.seek(max(0, log_path.stat().st_size - 128_000))
        log_text = log.read().decode("utf-8", errors="replace")
    if process.returncode:
        raise HydraRuntimeError(f"HYDRA exited with code {process.returncode}.\n{log_text[-5000:]}")
    return command, log_text


def _database_provenance(db_root, names, cancelled=None):
    root = Path(db_root).resolve()
    manifest = installed_databases(root)
    entries = {}
    for name in names:
        if name not in manifest:
            raise HydraRuntimeError(f"HYDRA database '{name}' is not installed. Download or choose a populated database store first.")
        entry = manifest[name]
        directory = _entry_directory(root, entry)
        files = {}
        references = list(directory.rglob("*"))
        if name == "protein" and (root / "mutation").is_dir():
            references.extend((root / "mutation").rglob("*"))
        for path in sorted(set(references)):
            if not path.is_file():
                continue
            if not path.resolve().is_relative_to(root):
                raise HydraRuntimeError(f"Database contains an external file link: {path}")
            # Index files are included: the exact reference/index snapshot is auditable.
            files[path.relative_to(root).as_posix()] = {"size": path.stat().st_size,
                                                       "sha256": file_sha256(path, cancelled)}
        entries[name] = {"manifest_entry": entry, "files": files}
    return {"root": str(root), "manifest_sha256": file_sha256(root / "manifest.json", cancelled),
            "databases": entries}


def run_assemblies(inputs, db_root, databases=None, *, sample_names=None, organism=None,
                   threads=2, protein=True, point_mutations=True, cancelled=None,
                   progress=None, work_root=None, min_identity=80, min_coverage=60,
                   protein_min_identity=90, protein_min_coverage=90):
    """Run the original HYDRA assembly pipeline and return its validated JSON report.

    Every run goes through `preflight` first, so a run that could not finish never
    starts and the caller gets a refusal naming each missing tool or reference set.
    What the run could not look for — an absent reference set, an organism with no
    catalogue, the age of the release — is recorded in the report's provenance and
    repeated in its import warnings, because an empty result and an unperformed
    search read identically otherwise.

    Input paths and reference files are never rewritten by this wrapper. The
    database store should be a managed writable cache because upstream may build
    missing indexes there. The caller must not concurrently update that store.
    """
    capabilities = runtime_capabilities(db_root)
    # One gate for every caller: a run that cannot finish never starts, and the
    # refusal names each missing tool or reference set and what it was for.
    checked = preflight(db_root, databases=databases, organism=organism,
                        point_mutations=point_mutations, capabilities=capabilities)
    if not checked["ready"]:
        raise HydraRuntimeError(checked["message"])
    if not isinstance(threads, int) or isinstance(threads, bool) or not 1 <= threads <= 256:
        raise HydraRuntimeError("Threads must be an integer between 1 and 256.")
    thresholds = {"min-identity": min_identity, "min-coverage": min_coverage,
                  "protein-min-identity": protein_min_identity,
                  "protein-min-coverage": protein_min_coverage}
    for name, value in thresholds.items():
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or not 0 <= value <= 100):
            raise HydraRuntimeError(f"{name} must be a finite percentage between 0 and 100.")
    paths = [Path(path).resolve() for path in inputs]
    if not paths or len(paths) != len(set(paths)):
        raise HydraRuntimeError("Choose one or more distinct assembly files.")
    names = list(sample_names) if sample_names is not None else None
    if names is not None and (len(names) != len(paths)
                              or any(not isinstance(name, str) or not name.strip() or any(c in name for c in "\n\r\t") for name in names)
                              or len(set(names)) != len(names)):
        raise HydraRuntimeError("Sample names must be unique nonempty labels, one per assembly.")
    signatures = {}
    input_provenance = []
    for index, path in enumerate(paths):
        if not path.is_file():
            raise HydraRuntimeError(f"Assembly was not found: {path}")
        with SequenceReader(path, cancelled) as reader:
            if reader.kind != "fasta":
                raise HydraRuntimeError(f"HYDRA execution accepts FASTA assemblies only: {path.name}")
        signatures[path] = file_signature(path)
        input_provenance.append({"path": str(path), "sha256": file_sha256(path, cancelled),
                                 "sample_name": names[index] if names is not None else sample_name(path)})
    names_db = list(databases) if databases is not None else list(checked["databases"])
    if not names_db or any(not isinstance(name, str) or not name or name.startswith("-") for name in names_db):
        raise HydraRuntimeError("Select at least one installed HYDRA database.")
    if progress:
        progress(0, 0, "Recording the selected database snapshot…")
    reference = _database_provenance(db_root, names_db, cancelled)
    arguments = ["run", "--db-dir", str(Path(db_root).resolve()), "--format", "json", "--prefix", "hydra",
                 "--threads", str(threads), "--no-banner", "--no-mlst", "--no-typing",
                 "--no-heteroresistance", "--no-reads-mlst", "--no-reads-variants"]
    for name, value in thresholds.items():
        arguments.extend(["--" + name, str(value)])
    # WMLSTudio's native typing evidence is not silently replaced by another engine.
    for name in names_db:
        arguments.extend(["--db", name])
    for path in paths:
        arguments.extend(["--assembly", str(path)])
    if names is not None:
        for name in names:
            arguments.extend(["--name", name])
    # The engine takes underscore-joined taxgroup names and stops the whole run on
    # one it does not know, so the assigned organism is resolved against the
    # installed catalogue first. A name with no catalogue is not passed and not
    # silently replaced by a neighbouring one: the run proceeds without mutation
    # catalogs and says so in the report.
    selected_organism = checked["organism"]["resolved"]
    if selected_organism:
        arguments.extend(["--organism", selected_organism])
    else:
        # Without native species sketches and an explicit organism, mutation
        # catalogs must not be chosen by an unvalidated inference.
        arguments.append("--no-auto-organism")
    if not protein:
        arguments.append("--no-protein")
    if not point_mutations:
        arguments.append("--no-point-mutations")
    if work_root is not None:
        Path(work_root).mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="wmlstudio-hydra-", dir=work_root) as work:
        arguments.extend(["--outdir", str(Path(work) / "results"), "--tmpdir", str(Path(work) / "temporary")])
        if progress:
            progress(0, 0, "Running pinned HYDRA assembly searches with native BLAST+…")
        command, log = _run_child(arguments, work, cancelled, progress)
        for path, signature in signatures.items():
            if file_signature(path) != signature:
                raise HydraRuntimeError(f"Assembly changed during HYDRA analysis: {path.name}. Run it again.")
        if file_sha256(Path(db_root) / "manifest.json", cancelled) != reference["manifest_sha256"]:
            raise HydraRuntimeError("The HYDRA database manifest changed during analysis. Do not update databases while a run is active.")
        report = load_hydra_report(Path(work) / "results/hydra.json")
    report["execution_provenance"] = {
        "engine": "HYDRA", "version": HYDRA_VERSION, "source_commit": HYDRA_COMMIT,
        "mode": "native-assembly", "command": command, "inputs": input_provenance,
        "distribution_origin": capabilities.get("distribution_origin", {}),
        "reference_snapshot": reference, "tools": capabilities["tools"],
        "completed_utc": datetime.now(UTC).isoformat(), "runtime_seconds": round(time.monotonic() - started, 3),
        "log": log, "scientific_limits": capabilities["limitations"],
        # What the run could and could not look for travels with the result, so a
        # report is never read as a complete screen when part of it never ran.
        "organism": dict(checked["organism"]),
        "reference_release": {key: checked["database_status"][key]
                              for key in ("release", "released_utc", "staged", "age_days",
                                          "stale", "bundled", "label")},
        "coverage_warnings": list(checked["warnings"]),
    }
    if checked["warnings"]:
        report.setdefault("import_warnings", []).extend(checked["warnings"])
    return report


def _ncbi_version():
    request = urllib.request.Request(NCBI_REFERENCE_ROOT + "/version.txt",
                                     headers={"User-Agent": "WMLSTudio-reference-provenance"})
    with urllib.request.urlopen(request, timeout=30) as response:
        value = response.read(256).decode("ascii").strip()
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}\.[0-9]+", value):
        raise HydraRuntimeError("NCBI returned an unrecognized reference version; snapshot was not published.")
    return value


def latest_release(*, cancelled=None):
    """Ask NCBI which reference release is current today, without downloading it.

    An explicit, user-initiated network call: nothing in startup, project open or
    analysis reaches here, and no sample data is sent — only a request for one
    version string. It exists so an update menu can say whether there is anything
    to fetch before asking a person to wait for 25 MB.
    """
    if cancelled and cancelled():
        raise AnalysisCancelled("Cancelled before the version check.")
    return _ncbi_version()


def update_databases(db_root, names, *, cancelled=None, progress=None, work_root=None):
    """Explicit user-initiated provider downloads; never called during analysis.

    Upstream owns normalization, version discovery and index construction. Its
    manifest records the selected providers. A selected-database-only snapshot is
    staged beside the requested store and atomically published under a new name.
    Existing stores and their analysis provenance remain untouched.
    """
    selected = list(names)
    if not selected or any(not isinstance(name, str) or not name or name.startswith("-") for name in selected):
        raise HydraRuntimeError("Select explicit database names for download/update.")
    # Nucleotide normalization consumes protein family annotations when present.
    # Fetch protein first so the original upstream importer can attach those
    # curated classes instead of leaving them absent in a brand-new store.
    selected = sorted(set(selected), key=lambda name: (name != "protein", name))
    capabilities = runtime_capabilities()
    if not capabilities["assembly_available"]:
        raise HydraRuntimeError("Install the complete native HYDRA / BLAST runtime before building reference databases.")
    root = Path(db_root).resolve()
    root.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    destination = root.parent / f"{root.name}-snapshot-{stamp}-{uuid.uuid4().hex[:8]}"
    provider_version = _ncbi_version() if set(selected).intersection({"ncbi", "protein"}) else None
    with tempfile.TemporaryDirectory(prefix=".wmlstudio-hydra-stage-", dir=root.parent) as work:
        staged = Path(work) / "snapshot"
        command, log = _run_child(["db", "download", *selected, "--db-dir", str(staged)],
                                  work, cancelled, progress)
        if provider_version:
            if _ncbi_version() != provider_version:
                raise HydraRuntimeError("NCBI changed its current database version during download; retry to obtain one consistent snapshot.")
            manifest_path = staged / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            for name in set(selected).intersection({"ncbi", "protein"}):
                entry = manifest.get("databases", {}).get(name)
                if entry:
                    entry["upstream_reported_version"] = entry.get("version", "unknown")
                    entry["version"] = provider_version
                    entry["version_source"] = NCBI_REFERENCE_ROOT + "/version.txt"
                    entry["source_url"] = NCBI_REFERENCE_ROOT
            manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        reference = _database_provenance(staged, selected, cancelled)
        if cancelled and cancelled():
            raise AnalysisCancelled("Database download cancelled before publishing; previous snapshots are unchanged.")
        os.replace(staged, destination)
    reference["root"] = str(destination)
    return {"command": command, "log": log, "reference_snapshot": reference,
            "database_root": str(destination), "previous_database_root": str(root),
            "snapshot_scope": "Selected databases only; previous snapshots are retained.",
            "completed_utc": datetime.now(UTC).isoformat(), "source_commit": HYDRA_COMMIT}


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments.pop(0) != "--upstream":
        raise SystemExit("This worker is launched by WMLSTudio; use the desktop to select assemblies and databases.")
    from hydra_amr.cli import main as hydra_main
    sys.argv = ["hydra", *arguments]
    return hydra_main(arguments)


if __name__ == "__main__":
    raise SystemExit(main())
