"""Genuine offline FastQC execution with app-local Java and immutable inputs.

FastQC module warnings/failures are reported verbatim as QC evidence, not as
process failures, trimming instructions, or proof of isolate purity.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from wmlstudio.assembly import _stop_process
from wmlstudio.sequence import SequenceReader, check_cancelled, file_sha256, file_signature

VERSION = "0.12.1"


class FastQCError(RuntimeError):
    pass


def runtime_capabilities(root=None):
    suffix = ".exe" if os.name == "nt" else ""
    roots = [Path(root)] if root else [
        Path(sys.executable).resolve().parent / "Tools/fastqc",
        Path(__file__).resolve().parent / "resources/tools/fastqc",
        Path(getattr(sys, "_MEIPASS", Path(__file__).parent)) / "Tools/fastqc",
    ]
    for candidate in roots:
        java = candidate / "jre/bin" / ("java" + suffix)
        required = [java, candidate / "uk/ac/babraham/FastQC/FastQCApplication.class",
                    candidate / "htsjdk.jar", candidate / "jbzip2-0.9.jar", candidate / "manifest.json"]
        if not all(path.is_file() for path in required):
            continue
        try:
            manifest = json.loads((candidate / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        expected = "windows-x64" if os.name == "nt" else "linux-x64"
        if manifest.get("fastqc_version") != VERSION or manifest.get("platform") != expected:
            continue
        return {"available": True, "java": str(java.resolve()), "fastqc_root": str(candidate.resolve()),
                "version": VERSION, "java_version": manifest.get("java_version"),
                "message": f"FastQC {VERSION} with bundled Temurin Java; no installation or network required."}
    return {"available": False, "java": None, "fastqc_root": None, "version": VERSION,
            "message": "Optional FastQC/Java tool package is not installed for this platform. Native WMLSTudio read QC remains available; it is not FastQC."}


def parse_report(path):
    """Read original summary and basic statistics; never execute report content."""
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        summaries = [name for name in names if name.endswith("/summary.txt")]
        datasets = [name for name in names if name.endswith("/fastqc_data.txt")]
        if len(summaries) != 1 or len(datasets) != 1:
            raise FastQCError("FastQC output lacks a unique summary/data report")
        if any(archive.getinfo(name).file_size > 32 * 1024 * 1024 for name in summaries + datasets):
            raise FastQCError("FastQC text report exceeds the safe parsing limit")
        modules = []
        for line in archive.read(summaries[0]).decode("utf-8").splitlines():
            status, name, filename = line.split("\t", 2)
            if status not in {"PASS", "WARN", "FAIL"}:
                raise FastQCError(f"Unknown FastQC module status: {status}")
            modules.append({"module": name, "status": status, "filename": filename})
        metrics, in_basic = {}, False
        data = archive.read(datasets[0]).decode("utf-8")
        if not data.startswith(f"##FastQC\t{VERSION}"):
            raise FastQCError("FastQC report version does not match the bundled engine")
        for line in data.splitlines():
            if line.startswith(">>Basic Statistics\t"):
                in_basic = True
            elif line == ">>END_MODULE":
                in_basic = False
            elif in_basic and line and not line.startswith("#"):
                key, value = line.split("\t", 1)
                metrics[key] = value
        if not modules or "Total Sequences" not in metrics:
            raise FastQCError("Incomplete FastQC summary")
    return {"modules": modules, "metrics": metrics,
            "qc_status": "FAIL" if any(m["status"] == "FAIL" for m in modules)
            else "WARN" if any(m["status"] == "WARN" for m in modules) else "PASS"}


def _run(command, directory, log_path, cancelled, progress):
    environment = os.environ.copy()
    for key in ("JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS", "_JAVA_OPTIONS", "CLASSPATH"):
        environment.pop(key, None)
    if "LD_LIBRARY_PATH_ORIG" in environment:
        environment["LD_LIBRARY_PATH"] = environment["LD_LIBRARY_PATH_ORIG"]
    elif getattr(sys, "frozen", False):
        environment.pop("LD_LIBRARY_PATH", None)
    environment["PATH"] = str(Path(command[0]).parent) + os.pathsep + environment.get("PATH", "")
    options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
    started = last_notice = time.monotonic()
    with Path(log_path).open("wb") as log:
        process = subprocess.Popen(command, cwd=directory, env=environment, shell=False,
                                   stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, **options)
        try:
            while True:
                check_cancelled(cancelled)
                try:
                    return process.wait(timeout=0.2)
                except subprocess.TimeoutExpired:
                    if progress and time.monotonic() - last_notice >= 2:
                        progress(0, 1, f"FastQC is reading the complete file ({int(time.monotonic()-started)} s); cancellation is available.")
                        last_notice = time.monotonic()
        finally:
            _stop_process(process)


def run_fastqc(read_paths, output_dir, *, threads=2, memory_gb=1, root=None,
               cancelled=None, progress=None):
    """Analyze complete FASTQs, one Java task/file; publish only complete output.

    Across-sample concurrency belongs to run_bounded; this adapter does not
    spawn an unbounded second queue. Java GC and visible processors are bounded.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("FastQC thread allocation must be a positive integer")
    if isinstance(memory_gb, bool) or not isinstance(memory_gb, int) or memory_gb < 1:
        raise ValueError("FastQC requires at least 1 GiB RAM reservation")
    check_cancelled(cancelled)
    capabilities = runtime_capabilities(root)
    if not capabilities["available"]:
        raise FastQCError(capabilities["message"])
    paths = [Path(path).expanduser().resolve() for path in read_paths]
    if not paths:
        raise ValueError("Choose at least one FASTQ input")
    if any(not path.is_file() for path in paths):
        raise ValueError("A selected FASTQ input is unavailable")
    if any(path.samefile(previous) for index, path in enumerate(paths) for previous in paths[:index]):
        raise ValueError("A FASTQ input was selected more than once")
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise ValueError("FastQC output already exists; choose a new report directory")
    signatures, inputs = [], []
    for path in paths:
        check_cancelled(cancelled)
        signature = file_signature(path)
        with SequenceReader(path, cancelled) as reader:
            if reader.kind != "fastq":
                raise ValueError("FastQC read analysis accepts FASTQ files only")
            next(iter(reader), None)
        checksum = file_sha256(path, cancelled)
        if signature != file_signature(path):
            raise FastQCError("A FASTQ input changed during validation")
        signatures.append(signature)
        inputs.append({"path": str(path), "sha256": checksum, "bytes": path.stat().st_size})
    tool_root = Path(capabilities["fastqc_root"])
    classpath = os.pathsep.join(str(tool_root / name) for name in ("", "htsjdk.jar", "jbzip2-0.9.jar", "cisd-jhdf5.jar"))
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=".fastqc-run-", dir=output_dir.parent) as temporary:
        staging = Path(temporary) / "reports"
        staging.mkdir()
        reports = []
        for index, path in enumerate(paths):
            check_cancelled(cancelled)
            relative = Path(f"read-{index+1:03d}")
            destination = staging / relative
            destination.mkdir()
            scratch = Path(temporary) / f"java-temp-{index}"
            scratch.mkdir()
            # 256 MiB is left outside the heap for JVM/native overhead. GC is
            # serial; one analysis file at a time avoids hidden nested fan-out.
            heap_mb = min(10000, memory_gb * 1024 - 256)
            command = [capabilities["java"], f"-Xmx{heap_mb}m", "-XX:+UseSerialGC",
                       f"-XX:ActiveProcessorCount={threads}", "-XX:+ExitOnOutOfMemoryError",
                       "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", f"-Djava.io.tmpdir={scratch}",
                       "-Dfastqc.threads=1", "-Dfastqc.sequence_format=fastq", "-Dfastqc.unzip=false",
                       f"-Dfastqc.output_dir={destination}", "-cp", classpath,
                       "uk.ac.babraham.FastQC.FastQCApplication", str(path)]
            log = destination / "execution.log"
            if progress:
                progress(index, len(paths), f"FastQC {index+1}/{len(paths)}: {path.name}")
            code = _run(command, tool_root, log, cancelled, progress)
            archives, html = list(destination.glob("*_fastqc.zip")), list(destination.glob("*_fastqc.html"))
            if code != 0 or len(archives) != 1 or len(html) != 1:
                tail = log.read_text(encoding="utf-8", errors="replace")[-3000:]
                raise FastQCError(f"FastQC did not produce a complete report (exit {code}): {tail}")
            parsed = parse_report(archives[0])
            parsed.update({"input": inputs[index], "report_html": str(output_dir / relative / html[0].name),
                           "report_zip": str(output_dir / relative / archives[0].name),
                           "report_sha256": file_sha256(archives[0], cancelled), "command": command})
            reports.append(parsed)
        for path, signature in zip(paths, signatures, strict=True):
            if file_signature(path) != signature:
                raise FastQCError("A FASTQ input changed while FastQC was reading it; report was not published")
        check_cancelled(cancelled)
        result = {"engine": "FastQC", "version": VERSION, "status": "completed", "inputs": inputs,
                  "reports": reports, "generated_at": datetime.now(timezone.utc).isoformat(),
                  "elapsed_seconds": round(time.monotonic() - started, 3),
                  "provenance": {"java": capabilities["java"], "java_version": capabilities["java_version"],
                                 "tool_manifest_sha256": file_sha256(tool_root / "manifest.json"),
                                 "threads_allocated": threads, "files_parallel": 1, "heap_mb": heap_mb,
                                 "sampling": "None; all reads analyzed", "preprocessing": "None; original files unchanged",
                                 "interpretation": "FastQC flags are QC evidence, not automatic trimming recommendations or isolate purity proof."}}
        (staging / "fastqc-result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        os.replace(staging, output_dir)
    if progress:
        progress(len(paths), len(paths), "FastQC reports complete; original FASTQs unchanged.")
    return result
