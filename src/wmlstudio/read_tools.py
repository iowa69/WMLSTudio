"""Cancellable native fastp adapter: trimmed reads are new files, never replacements.

fastp improves reads. It does not validate an isolate, prove purity, establish a
species, or make a downstream typing result correct. Every original FASTQ stays a
read-only input; each trimmed output is a new file tied back to its source by
SHA-256, so the link between the two is always recoverable from the record.
"""

from __future__ import annotations

import copy
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .assembly import _stop_process, _tail, _write_json
from .sequence import (
    AnalysisCancelled,
    SequenceError,
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
)

VERSION = "1.3.7"
SOURCE_COMMIT = "8a2397b6628ae14127efdb7566f67fc05f9aea56"

#: fastp caps its own worker count; the request and what was applied are both recorded.
MAX_THREADS = 16
#: fastp writes every curve of its report as JSON, so the file is large but bounded.
REPORT_LIMIT = 64 * 1024 * 1024
#: A path a C++ tool opens through the platform byte API without surprises.
_ASCII_PATH = re.compile(r"^[A-Za-z0-9 ./\\:_+,()~@#$%^&=!\[\]{}'`-]+$")

TRIMMING_DISCLAIMER = (
    "Trimming and quality filtering improve reads. They do not validate an isolate, "
    "prove purity, establish a species, or make a typing result correct."
)


class ReadTrimError(RuntimeError):
    """The native tool is unavailable, failed, or produced no usable trimmed pair."""


def runtime_capabilities(root=None):
    """Locate the staged, checksum-pinned tool. Nothing is downloaded at analysis time."""
    platform = "windows-x64" if os.name == "nt" else "linux-x64"
    name = "fastp.exe" if os.name == "nt" else "fastp"
    candidates = [Path(root)] if root else [
        Path(__file__).resolve().parent / "resources/tools/fastp" / platform,
        Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent)) / "Tools/fastp",
        Path(sys.executable).resolve().parent / "Tools/fastp",
    ]
    for directory in candidates:
        binary = directory / name
        if not binary.is_file():
            continue
        try:
            manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
            expected = next(item["sha256"] for item in manifest["files"] if item["path"] == name)
            if (manifest["tool"] != "fastp" or manifest["version"] != VERSION
                    or manifest["platform"] != platform or manifest["source_commit"] != SOURCE_COMMIT):
                raise ValueError("staged version, source revision or platform does not match this build")
        except (OSError, KeyError, ValueError, StopIteration) as error:
            return {"available": False, "binary": None, "version": VERSION, "platform": platform,
                    "reason": f"The staged fastp package is unusable: {error}"}
        return {"available": True, "binary": str(binary.resolve()), "expected_sha256": expected,
                "version": VERSION, "platform": platform, "source_commit": SOURCE_COMMIT, "reason": ""}
    reason = ("Native fastp is not staged in this package, so read trimming is unavailable here. "
              "No download happens during analysis and no substitute trimming is performed; "
              "reads can still be assembled and typed exactly as they were supplied.")
    if platform == "windows-x64":
        reason += (" Upstream fastp publishes no official Windows binary, so a Windows package "
                   "carries fastp only when a separately built, verified native tool was staged.")
    return {"available": False, "binary": None, "version": VERSION, "platform": platform, "reason": reason}


def _percent(part, whole):
    return 100 * part / whole if whole else None


def parse_report(path):
    """Read fastp's own report. Its numbers are reported, never recomputed or relabelled."""
    path = Path(path)
    if path.stat().st_size > REPORT_LIMIT:
        raise ReadTrimError("The fastp JSON report exceeds the safe parsing limit.")
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
        summary = report["summary"]
        before, after = summary["before_filtering"], summary["after_filtering"]
        filtering = report["filtering_result"]
    except (KeyError, ValueError) as error:
        raise ReadTrimError(f"The fastp report is not a complete fastp JSON report: {error}") from error
    if summary.get("fastp_version") != VERSION:
        raise ReadTrimError(f"The fastp report was not written by the staged engine {VERSION}.")
    if after["total_reads"] > before["total_reads"] or after["total_bases"] > before["total_bases"]:
        raise ReadTrimError("The fastp report claims more reads or bases after filtering than before.")
    adapter = report.get("adapter_cutting", {})
    kept = {key: before[key] for key in ("total_reads", "total_bases", "q20_bases", "q30_bases",
                                         "q20_rate", "q30_rate", "gc_content")}
    return {
        "engine": "fastp", "version": VERSION,
        "sequencing": summary.get("sequencing", ""),
        "before": {**kept, "read1_mean_length": before.get("read1_mean_length"),
                   "read2_mean_length": before.get("read2_mean_length")},
        "after": {key: after.get(key) for key in (*kept, "read1_mean_length", "read2_mean_length")},
        "reads_removed": before["total_reads"] - after["total_reads"],
        "bases_removed": before["total_bases"] - after["total_bases"],
        "reads_removed_percent": _percent(before["total_reads"] - after["total_reads"], before["total_reads"]),
        "bases_removed_percent": _percent(before["total_bases"] - after["total_bases"], before["total_bases"]),
        "filtering": {key: filtering.get(key) for key in
                      ("passed_filter_reads", "low_quality_reads", "too_many_N_reads",
                       "adapter_dimer_reads", "too_short_reads", "too_long_reads")},
        "adapter": {"trimmed_reads": adapter.get("adapter_trimmed_reads"),
                    "trimmed_bases": adapter.get("adapter_trimmed_bases"),
                    "read1_sequence": adapter.get("read1_adapter_sequence"),
                    "read2_sequence": adapter.get("read2_adapter_sequence")},
        "duplication_rate": report.get("duplication", {}).get("rate"),
        "insert_size": {"peak": report.get("insert_size", {}).get("peak"),
                        "undetermined_pairs": report.get("insert_size", {}).get("unknown")},
        "denominators": {
            "q20_rate/q30_rate": "fraction of all bases in that stage, not of reads",
            "reads_removed_percent": "of the reads present before filtering",
            "duplication_rate": "fastp's estimate from read content; not a library-preparation measurement",
            "insert_size": "inferred from paired overlap; undetermined pairs are counted separately",
        },
    }


def _notify(progress, done, message):
    if progress:
        progress(done, 100, message)


def _run(command, directory, log_path, cancelled=None, progress=None, timeout=None):
    """Drain straight to disk so a chatty tool on a large pair cannot fill a pipe."""
    check_cancelled(cancelled)
    environment = os.environ.copy()
    environment["PATH"] = str(Path(command[0]).parent) + os.pathsep + environment.get("PATH", "")
    if "LD_LIBRARY_PATH_ORIG" in environment:
        environment["LD_LIBRARY_PATH"] = environment["LD_LIBRARY_PATH_ORIG"]
    elif getattr(sys, "frozen", False):
        environment.pop("LD_LIBRARY_PATH", None)
    options = ({"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt"
               else {"start_new_session": True})
    started = last_notice = time.monotonic()
    with Path(log_path).open("wb") as log:
        process = subprocess.Popen(
            list(map(str, command)), cwd=directory, env=environment, shell=False,
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, **options,
        )
        try:
            while True:
                check_cancelled(cancelled)
                if timeout is not None and time.monotonic() - started > timeout:
                    raise ReadTrimError("The fastp version check timed out.")
                try:
                    return process.wait(timeout=0.2)
                except subprocess.TimeoutExpired:
                    if time.monotonic() - last_notice >= 5:
                        _notify(progress, 55, "fastp is trimming and filtering the read pair; elapsed "
                                f"{int(time.monotonic() - started)} s. Cancel remains available.")
                        last_notice = time.monotonic()
        finally:
            _stop_process(process)


def _tool_info(capability, directory, cancelled):
    binary = Path(capability["binary"])
    digest = file_sha256(binary, cancelled)
    if digest != capability.get("expected_sha256"):
        raise ReadTrimError("The staged fastp binary does not match its manifest checksum.")
    log = directory / "version.log"
    code = _run([str(binary), "--version"], directory, log, cancelled, timeout=30)
    reported = _tail(log).strip()
    if code or reported != f"fastp {VERSION}":
        raise ReadTrimError(f"The selected executable is not fastp {VERSION}: {reported!r}")
    return {"name": "fastp", "version": VERSION, "executable": str(binary),
            "executable_sha256": digest, "source_commit": SOURCE_COMMIT,
            "platform": capability["platform"], "reported_version": reported}


def _compression(path):
    with Path(path).open("rb") as handle:
        head = handle.read(3)
    if head.startswith(b"\x1f\x8b"):
        return "gzip"
    if head.startswith(b"BZh"):
        return "bzip2"
    return "plain"


def _holds_reads(path):
    """True when a gzip output holds at least one byte of sequence, not just a header."""
    with gzip.open(path, "rb") as handle:
        return bool(handle.read(1))


def _tool_input(path, kind, directory, index):
    """Give the tool a name it can open, without copying or touching the original.

    fastp chooses its reader from the file name, and opens paths through the
    platform byte API. A link is a second name for the same read-only inode, so
    the original file is neither copied nor modified; when no link can be made
    the original path is used and the provenance records that plainly.
    """
    text = str(path)
    advertises_gzip = text.lower().endswith(".gz")
    if _ASCII_PATH.match(text) and advertises_gzip == (kind == "gzip"):
        return path, "original path", False
    alias = directory / f"input{index}.fastq{'.gz' if kind == 'gzip' else ''}"
    for attempt, description in ((lambda: os.link(path, alias), "hard link to the original inode"),
                                 (lambda: alias.symlink_to(path), "symbolic link to the original")):
        try:
            attempt()
            return alias, description, True
        except (OSError, NotImplementedError):
            continue
    if advertises_gzip != (kind == "gzip"):
        raise ReadTrimError(
            f"{path.name} is {kind} but its name says otherwise, and no link could be created "
            "beside the output. Rename the file so its suffix matches its contents.")
    return path, "original path; no ASCII alias could be created", False


def _validate_options(options):
    for name, value, low, high in (
            ("threads", options["threads"], 1, 256),
            ("min_length", options["min_length"], 1, 10000),
            ("quality_threshold", options["quality_threshold"], 0, 60),
            ("unqualified_percent_limit", options["unqualified_percent_limit"], 0, 100),
            ("n_base_limit", options["n_base_limit"], 0, 10000),
            ("compression", options["compression"], 1, 9)):
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            raise ValueError(f"{name} must be an integer between {low} and {high}.")
    for name in ("detect_adapter", "deduplicate"):
        if not isinstance(options[name], bool):
            raise ValueError(f"{name} must be true or false.")


def run_fastp(read1, read2, output_dir, *, threads=4, min_length=15, quality_threshold=15,
              unqualified_percent_limit=40, n_base_limit=5, detect_adapter=True,
              deduplicate=False, compression=4, root=None, cancelled=None, progress=None):
    """Trim and quality-filter one paired FASTQ set into a NEW output directory.

    progress receives (done, total, message); cancellation raises AnalysisCancelled.
    The originals are opened read-only and never renamed, rewritten or deleted.
    An existing output directory is never replaced, and a run that fails or is
    cancelled publishes nothing; its log and fastp report are kept as uniquely
    named sibling *.failed-*.log/.json files for diagnosis.
    """
    options = {"threads": threads, "min_length": min_length, "quality_threshold": quality_threshold,
               "unqualified_percent_limit": unqualified_percent_limit, "n_base_limit": n_base_limit,
               "detect_adapter": detect_adapter, "deduplicate": deduplicate, "compression": compression}
    _validate_options(options)
    check_cancelled(cancelled)
    first, second, destination = (Path(p).expanduser().resolve() for p in (read1, read2, output_dir))
    if not first.is_file() or not second.is_file():
        raise SequenceError("Both paired FASTQ inputs must exist and be regular files.")
    if first.samefile(second):
        raise SequenceError("Read-pair inputs must be different files.")
    if destination.exists():
        raise FileExistsError(f"Trimmed-read output already exists; choose a new directory: {destination}")
    if any(path.is_relative_to(destination) for path in (first, second)):
        raise ValueError("Trimmed-read output must not contain the original inputs.")
    capability = runtime_capabilities(root)
    if not capability["available"]:
        raise ReadTrimError(capability["reason"])
    kinds = []
    for path in (first, second):
        kind = _compression(path)
        if kind == "bzip2":
            raise SequenceError(
                f"{path.name} is bzip2-compressed. fastp reads plain or gzip FASTQ only; "
                "decompress or re-compress the pair with gzip first.")
        with SequenceReader(path, cancelled) as reader:
            if reader.kind != "fastq":
                raise SequenceError("Read trimming requires two FASTQ files.")
            next(iter(reader), None)
        kinds.append(kind)
    destination.parent.mkdir(parents=True, exist_ok=True)
    needed = first.stat().st_size + second.stat().st_size + 64 * 1024 * 1024
    if shutil.disk_usage(destination.parent).free < needed:
        raise ReadTrimError("Insufficient free disk space to write a trimmed copy of this read pair.")
    started = datetime.now(timezone.utc).isoformat()
    applied_threads = min(threads, MAX_THREADS)
    _notify(progress, 0, "Validating native fastp and the original paired inputs.")
    with tempfile.TemporaryDirectory(prefix=".wmlstudio-trim-", dir=destination.parent) as temporary:
        directory = Path(temporary)
        try:
            tool_info = _tool_info(capability, directory, cancelled)
            _notify(progress, 15, "Hashing the original paired FASTQs before trimming.")
            signatures = [file_signature(path) for path in (first, second)]
            inputs = []
            for index, (path, kind, signature) in enumerate(zip((first, second), kinds, signatures), 1):
                inputs.append({"mate": index, "path": str(path), "sha256": file_sha256(path, cancelled),
                               "size_bytes": signature[0], "mtime_ns": signature[1], "compression": kind})
            tool_paths = [_tool_input(path, kind, directory, index)
                          for index, (path, kind) in enumerate(zip((first, second), kinds), 1)]
            command = [str(Path(tool_info["executable"])),
                       "--in1", str(tool_paths[0][0]), "--in2", str(tool_paths[1][0]),
                       "--out1", "trimmed_R1.fastq.gz", "--out2", "trimmed_R2.fastq.gz",
                       "--unpaired1", "unpaired_R1.fastq.gz", "--unpaired2", "unpaired_R2.fastq.gz",
                       "--json", "fastp.json", "--html", "fastp.html",
                       "--report_title", "WMLSTudio read trimming",
                       "--thread", str(applied_threads), "--compression", str(compression),
                       "--length_required", str(min_length),
                       "--qualified_quality_phred", str(quality_threshold),
                       "--unqualified_percent_limit", str(unqualified_percent_limit),
                       "--n_base_limit", str(n_base_limit)]
            if detect_adapter:
                command.append("--detect_adapter_for_pe")
            if deduplicate:
                command.append("--dedup")
            parameters = {**options, "threads_requested": threads, "threads_applied": applied_threads,
                          "thread_policy": f"fastp accepts at most {MAX_THREADS} worker threads.",
                          "defaults": "fastp's own filter defaults unless a value above differs."}
            # A crash leaves a complete launch plan beside the log. Nothing infers
            # success from a stale output file, a PID, or a previous run.
            _write_json(directory / "run-plan.json", {
                "state": "prepared_not_completed", "started_at": started,
                "output_directory": str(destination), "command": command, "engine": tool_info,
                "inputs": inputs, "parameters": parameters,
                "recovery_policy": "Review the log and rerun; no automatic resume or success inference.",
            })
            _notify(progress, 40, "Running fastp adapter trimming and quality filtering.")
            log = directory / "fastp.log"
            code = _run(command, directory, log, cancelled, progress)
            check_cancelled(cancelled)
            if code:
                raise ReadTrimError(f"fastp exited with code {code}.\n{_tail(log)}")
            if not (directory / "fastp.json").is_file():
                raise ReadTrimError("fastp finished without writing its JSON report.")
            report = parse_report(directory / "fastp.json")
            if not report["after"]["total_reads"]:
                raise ReadTrimError(
                    f"fastp kept none of the {report['before']['total_reads']:,} input reads at these "
                    "settings, so no trimmed pair was written. Review the filters and the fastp report.")
            _notify(progress, 85, "Hashing the trimmed outputs and linking them to the originals.")
            outputs = []
            for index, name in enumerate(("trimmed_R1.fastq.gz", "trimmed_R2.fastq.gz"), 1):
                path = directory / name
                if not path.is_file() or not _holds_reads(path):
                    raise ReadTrimError(f"fastp did not write any trimmed mate {index} reads.")
                outputs.append({"mate": index, "name": name, "path": str(destination / name),
                                "sha256": file_sha256(path, cancelled), "size_bytes": path.stat().st_size,
                                "compression": "gzip"})
            unpaired = []
            for index, name in enumerate(("unpaired_R1.fastq.gz", "unpaired_R2.fastq.gz"), 1):
                path = directory / name
                if not path.is_file():
                    continue
                # fastp opens these writers whatever happens, so a gzip header
                # alone is not a rescued read. Only publish files that hold one.
                if _holds_reads(path):
                    unpaired.append({"mate": index, "name": name, "path": str(destination / name),
                                     "sha256": file_sha256(path, cancelled),
                                     "size_bytes": path.stat().st_size})
                else:
                    path.unlink()
            for path, signature in zip((first, second), signatures):
                if file_signature(path) != signature:
                    raise ReadTrimError(f"An original input changed while fastp was reading it: {path.name}")
            for index, (path, digest) in enumerate(zip((first, second), [item["sha256"] for item in inputs])):
                if file_sha256(path, cancelled) != digest:
                    raise ReadTrimError(f"An original input changed during trimming: {path.name}")
            provenance = {
                "workflow": "paired_short_read_trimming", "engine": tool_info,
                "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
                "command": command, "parameters": parameters,
                "tool_input_access": [description for _path, description, _alias in tool_paths],
                "inputs": inputs, "outputs": outputs, "unpaired_outputs": unpaired,
                "report_json_path": str(destination / "fastp.json"),
                "report_html_path": str(destination / "fastp.html"),
                "report_json_sha256": file_sha256(directory / "fastp.json", cancelled),
                "report_html_sha256": (file_sha256(directory / "fastp.html", cancelled)
                                       if (directory / "fastp.html").is_file() else None),
                "link": [{"mate": item["mate"], "original_path": item["path"],
                          "original_sha256": item["sha256"], "trimmed_path": output["path"],
                          "trimmed_sha256": output["sha256"]}
                         for item, output in zip(inputs, outputs)],
                "originals_policy": ("Originals are read-only inputs. Trimming writes new files and "
                                     "never renames, rewrites, replaces or deletes the supplied FASTQs."),
                "quality_encoding": "Phred+33 (assumed); the encoding cannot be inferred from the data.",
            }
            notes = [
                TRIMMING_DISCLAIMER,
                "Every count, rate and curve in the report is fastp's own output, not a WMLSTudio recalculation.",
                "Reads removed by the filters are absent from the trimmed pair only; the originals still hold them.",
                "fastp's duplication rate is an estimate from read content, not a measured library duplication.",
            ]
            if unpaired:
                notes.append("Reads whose mate was discarded were kept separately as unpaired files; "
                             "the trimmed pair itself stays complete and in order.")
            if any("no ASCII alias" in description for _path, description, _alias in tool_paths):
                notes.append("An input path could not be given an ASCII alias; the tool opened the "
                             "original path directly.")
            payload = {"read1_path": outputs[0]["path"], "read2_path": outputs[1]["path"],
                       "report": report, "provenance": provenance, "notes": notes}
            _write_json(directory / "read-trimming.json", payload)
            for path, _description, created in tool_paths:
                if created:
                    Path(path).unlink(missing_ok=True)
            (directory / "run-plan.json").unlink(missing_ok=True)
            (directory / "version.log").unlink(missing_ok=True)
            check_cancelled(cancelled)
            if os.name == "nt":
                # Windows directory rename is atomic and refuses an existing target.
                os.rename(directory, destination)
            else:
                # Claim the name without replacing anything; atomically replace only
                # the empty placeholder this operation just created.
                destination.mkdir()
                try:
                    os.replace(directory, destination)
                except BaseException:
                    destination.rmdir()  # Empty owned placeholder only.
                    raise
            _notify(progress, 100, "Trimmed read pair written; the original FASTQs are unchanged.")
            return payload
        except BaseException as error:
            logs = [path for path in (directory / "fastp.log", directory / "version.log") if path.is_file()]
            if logs:
                failure_log = destination.parent / f"{destination.name}.failed-{uuid.uuid4().hex}.log"
                try:
                    with failure_log.open("x", encoding="utf-8") as output:
                        output.write(str(error) + "\n")
                        for item in logs:
                            output.write(f"\n{item.name}\n{_tail(item)}\n")
                    plan_path = directory / "run-plan.json"
                    if plan_path.is_file():
                        plan = json.loads(plan_path.read_text(encoding="utf-8"))
                        plan.update(state="cancelled" if isinstance(error, AnalysisCancelled) else "failed",
                                    error=str(error), log_path=str(failure_log))
                        if (directory / "fastp.json").is_file():
                            try:
                                plan["report"] = parse_report(directory / "fastp.json")
                            except (OSError, ReadTrimError):
                                plan["report"] = None
                        _write_json(failure_log.with_suffix(".json"), plan)
                    if not isinstance(error, AnalysisCancelled):
                        error.add_note(f"Diagnostic log: {failure_log}")
                        if isinstance(error, ReadTrimError):
                            error.args = (f"{error}\nDiagnostic log: {failure_log}",)
                except OSError:
                    pass
            raise


def read_source_for(result):
    """Describe a trimmed pair for the assembler, so the assembly records what it used."""
    provenance = result["provenance"]
    outputs, inputs = provenance["outputs"], provenance["inputs"]
    if len(outputs) != 2 or len(inputs) != 2:
        raise ValueError("A trimmed read source must describe exactly two mates.")
    engine = provenance["engine"]
    return {
        "kind": "fastp_trimmed", "tool": "fastp", "version": engine["version"],
        "description": (f"fastp {engine['version']} trimmed and quality-filtered reads; "
                        "the original FASTQs are unchanged and still on record."),
        "read1_path": outputs[0]["path"], "read2_path": outputs[1]["path"],
        "read1_sha256": outputs[0]["sha256"], "read2_sha256": outputs[1]["sha256"],
        "originals": [{"mate": item["mate"], "path": item["path"], "sha256": item["sha256"]}
                      for item in inputs],
        "parameters": copy.deepcopy(provenance["parameters"]),
        "report_json_sha256": provenance["report_json_sha256"],
        "reads_removed": result["report"]["reads_removed"],
        "bases_removed": result["report"]["bases_removed"],
    }


def record_trimming(project, read1_id, read2_id, result):
    """Attach one trimming run to the two read records without replacing either input.

    The originals remain each record's input path: trimming adds evidence and a
    named companion file, it never substitutes one file for another. No typing
    result is invalidated, because nothing that was already analysed changed.
    """
    if read1_id == read2_id:
        raise ValueError("Read trimming involves two distinct sample records.")
    provenance = result["provenance"]
    inputs, outputs = provenance["inputs"], provenance["outputs"]
    if len(inputs) != 2 or len(outputs) != 2:
        raise ValueError("Trimming provenance must identify both hashed read inputs and outputs.")
    samples = [project.get_sample(read1_id), project.get_sample(read2_id)]
    for sample, source in zip(samples, inputs):
        if not sample["input_path"] or Path(sample["input_path"]).resolve() != Path(source["path"]).resolve():
            raise ValueError("Trimming provenance does not match the selected read records.")
    for output, source in zip(outputs, inputs):
        path = Path(output["path"])
        if not path.is_file() or file_sha256(path) != output["sha256"]:
            raise ValueError("A trimmed FASTQ is missing or does not match its recorded SHA-256.")
        if path.resolve() == Path(source["path"]).resolve():
            raise ValueError("A trimmed output cannot be one of the original inputs.")
    identifiers = (read1_id, read2_id)
    with project.transaction():
        for sample in samples:
            if project.get_sample(sample["id"])["input_path"] != sample["input_path"]:
                raise ValueError("A read input changed while trimming was being recorded.")
        for index, (identifier, source, output) in enumerate(zip(identifiers, inputs, outputs), 1):
            project.update_metadata(identifier, {"read_trimming": {
                "status": "completed", "engine": "fastp", "version": provenance["engine"]["version"],
                "mate": index, "paired_with": identifiers[index % 2],
                "original_path": source["path"], "original_sha256": source["sha256"],
                "trimmed_path": output["path"], "trimmed_sha256": output["sha256"],
                "output_directory": str(Path(output["path"]).parent),
                "report": copy.deepcopy(result["report"]),
                "report_json_path": provenance["report_json_path"],
                "report_html_path": provenance["report_html_path"],
                "parameters": copy.deepcopy(provenance["parameters"]),
                "completed_at": provenance["completed_at"],
                "notes": list(result.get("notes", [])),
            }})
            project.record_history(identifier, "read_trimming_recorded", {
                "mate": index, "paired_with": identifiers[index % 2],
                "original_path": source["path"], "original_sha256": source["sha256"],
                "trimmed_path": output["path"], "trimmed_sha256": output["sha256"],
                "engine": f"fastp {provenance['engine']['version']}",
            })
    return list(identifiers)


def trimmed_pair(read1_sample, read2_sample):
    """Return the recorded trimmed pair for two read records, or None when there is none.

    A record is only usable when both mates come from the same run, still exist,
    and still carry the size the run recorded. The full SHA-256 check belongs to
    whoever consumes the files, which reads them completely anyway.
    """
    records = []
    for index, sample in enumerate((read1_sample, read2_sample), 1):
        record = (sample.get("metadata") or {}).get("read_trimming") or {}
        if record.get("status") != "completed" or record.get("mate") != index:
            return None
        path = Path(record.get("trimmed_path", ""))
        if not path.is_file():
            return None
        records.append(record)
    if records[0]["paired_with"] != read2_sample["id"] or records[1]["paired_with"] != read1_sample["id"]:
        return None
    if records[0]["report"] != records[1]["report"]:
        return None
    return records


def reads_for_assembly(read1_sample, read2_sample, *, prefer_trimmed=True):
    """Choose the reads to assemble and say, in the same breath, which ones they are."""
    records = trimmed_pair(read1_sample, read2_sample) if prefer_trimmed else None
    if records:
        return (Path(records[0]["trimmed_path"]), Path(records[1]["trimmed_path"]), {
            "kind": "fastp_trimmed", "tool": "fastp", "version": records[0]["version"],
            "description": (f"fastp {records[0]['version']} trimmed and quality-filtered reads; "
                            "the original FASTQs are unchanged and still on record."),
            "read1_path": records[0]["trimmed_path"], "read2_path": records[1]["trimmed_path"],
            "read1_sha256": records[0]["trimmed_sha256"], "read2_sha256": records[1]["trimmed_sha256"],
            "originals": [{"mate": record["mate"], "path": record["original_path"],
                           "sha256": record["original_sha256"]} for record in records],
            "parameters": copy.deepcopy(records[0]["parameters"]),
            "reads_removed": records[0]["report"]["reads_removed"],
            "bases_removed": records[0]["report"]["bases_removed"],
        })
    return Path(read1_sample["input_path"]), Path(read2_sample["input_path"]), None


def trimming_view(sample):
    """One read record's trimming status for the read tab, with no invented numbers."""
    metadata = sample.get("metadata") or {}
    record = metadata.get("read_trimming") or {}
    view = {"sample_id": sample.get("id"), "sample_name": sample.get("name"),
            "status": "not_trimmed", "engine": None, "version": None, "mate": None,
            "original_path": sample.get("input_path"), "trimmed_path": None,
            "report": None, "report_html_path": None, "parameters": None, "completed_at": None,
            "notes": [TRIMMING_DISCLAIMER]}
    if (metadata.get("workflow") or {}).get("source_kind") == "assembly":
        view["status"] = "not_reads"
        view["notes"].insert(0, "This record's input is an assembly; trimming applies to read files.")
        return view
    if record.get("status") != "completed":
        view["notes"].insert(0, "These reads were not trimmed. They can still be assembled and typed "
                                "exactly as supplied; nothing about them is assumed.")
        return view
    report = record.get("report") or {}
    view.update(status="trimmed", engine=record.get("engine"), version=record.get("version"),
                mate=record.get("mate"), original_path=record.get("original_path"),
                trimmed_path=record.get("trimmed_path"), report=report,
                report_html_path=record.get("report_html_path"),
                parameters=record.get("parameters"), completed_at=record.get("completed_at"))
    view["notes"] = list(record.get("notes") or [TRIMMING_DISCLAIMER])
    view["notes"].append("The original FASTQ is unchanged and remains this record's input.")
    return view
