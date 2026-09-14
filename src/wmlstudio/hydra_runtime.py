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
            target = (root / entry["path"]).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"database path leaves the selected store: {name}")
            if target.exists():
                result[name] = dict(entry)
        return result
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise HydraRuntimeError(f"Cannot read HYDRA database manifest: {exc}") from exc


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
    message = ("Native HYDRA assembly runtime is ready." if available else
               f"HYDRA {version or 'not installed'}; missing tools: {', '.join(missing) or 'none'}.")
    return {"hydra_version": version, "expected_hydra_version": HYDRA_VERSION,
            "source_commit": HYDRA_COMMIT, "distribution_origin": origin,
            "tools": tools, "databases": databases,
            "assembly_available": available, "available": available, "message": message,
            "blastn": tools["blastn"], "makeblastdb": tools["makeblastdb"],
            "reads_available": False,
            "limitations": ["Assembly-only HYDRA execution; its direct-read and minority-allele pipelines are not enabled. Paired-read assembly is handled separately by native SKESA.",
                            "Gene and mutation evidence is not a validated susceptibility phenotype."]}


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
        directory = (root / entry["path"]).resolve()
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

    Input paths and reference files are never rewritten by this wrapper. The
    database store should be a managed writable cache because upstream may build
    missing indexes there. The caller must not concurrently update that store.
    """
    capabilities = runtime_capabilities(db_root)
    if not capabilities["assembly_available"]:
        missing = [name for name, path in capabilities["tools"].items() if not path]
        raise HydraRuntimeError(f"Native HYDRA {HYDRA_VERSION} assembly runtime is incomplete. Missing tools: {', '.join(missing) or 'none'}; installed HYDRA: {capabilities['hydra_version'] or 'not installed'}.")
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
    available = capabilities["databases"]
    names_db = list(databases) if databases is not None else [name for name in ("ncbi", "protein") if name in available]
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
    if organism:
        arguments.extend(["--organism", str(organism)])
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
    }
    if not organism and point_mutations:
        report.setdefault("import_warnings", []).append("No organism was specified: organism-specific point mutation catalogs were not selected.")
    return report


def _ncbi_version():
    request = urllib.request.Request(NCBI_REFERENCE_ROOT + "/version.txt",
                                     headers={"User-Agent": "WMLSTudio-reference-provenance"})
    with urllib.request.urlopen(request, timeout=30) as response:
        value = response.read(256).decode("ascii").strip()
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}\.[0-9]+", value):
        raise HydraRuntimeError("NCBI returned an unrecognized reference version; snapshot was not published.")
    return value


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
