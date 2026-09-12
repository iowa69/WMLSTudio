"""Native, cancellable SKESA adapter for explicitly paired short-read isolates.

No simulation or automatic download is used. Successful process execution is
not a biological purity/completeness pass. Original FASTQs are never modified.
"""

from __future__ import annotations

import copy
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .sequence import (
    AnalysisCancelled,
    QCAccumulator,
    SequenceError,
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
    inspect_sequence,
    read_pair_identity,
)


class AssemblyError(RuntimeError):
    """The native tool is unavailable, failed, or did not produce valid contigs."""


def resolve_skesa(executable=None) -> Path:
    """Find a separately packaged native executable, never download or emulate."""
    if executable is not None:
        candidate = Path(executable).expanduser().resolve()
        if not candidate.is_file():
            raise AssemblyError(f"SKESA executable does not exist: {candidate}")
        return candidate
    suffix = ".exe" if os.name == "nt" else ""
    roots = [
        Path(sys.executable).resolve().parent / "Tools" / "skesa",
        Path(__file__).resolve().parent / "resources" / "tools" / "skesa",
    ]
    if getattr(sys, "_MEIPASS", None):
        roots.append(Path(sys._MEIPASS) / "wmlstudio" / "resources" / "tools" / "skesa")
    for root in roots:
        candidate = root / ("skesa" + suffix)
        if candidate.is_file():
            return candidate.resolve()
    found = shutil.which("skesa" + suffix)
    if found:
        return Path(found).resolve()
    raise AssemblyError(
        "Native SKESA is not installed. Add the verified SKESA tool package to "
        "Tools/skesa beside WMLSTudio, or choose a SKESA executable. "
        "FASTQ assembly is unavailable until a real native engine is present."
    )


def _notify(progress, done, message):
    if progress:
        progress(done, 100, message)


def _stop_process(process):
    if process.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill.exe", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                creationflags=subprocess.CREATE_NO_WINDOW, timeout=10, check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
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


def _run_process(command, directory, log_path, cancelled=None, progress=None, timeout=None):
    """Drain directly to disk so verbose subprocesses cannot deadlock a pipe."""
    check_cancelled(cancelled)
    environment = os.environ.copy()
    environment["PATH"] = str(Path(command[0]).parent) + os.pathsep + environment.get("PATH", "")
    if "LD_LIBRARY_PATH_ORIG" in environment:
        environment["LD_LIBRARY_PATH"] = environment["LD_LIBRARY_PATH_ORIG"]
    elif getattr(sys, "frozen", False):
        environment.pop("LD_LIBRARY_PATH", None)
    options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
    started, last_notice = time.monotonic(), time.monotonic()
    with Path(log_path).open("wb") as log:
        process = subprocess.Popen(
            list(map(str, command)), cwd=directory, env=environment,
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, shell=False, **options,
        )
        try:
            while True:
                check_cancelled(cancelled)
                if timeout is not None and time.monotonic() - started > timeout:
                    raise AssemblyError("SKESA version check timed out.")
                try:
                    return process.wait(timeout=0.2)
                except subprocess.TimeoutExpired:
                    if time.monotonic() - last_notice >= 5:
                        _notify(progress, 60, "SKESA is assembling reads; elapsed "
                                f"{int(time.monotonic() - started)} s. Cancel remains available.")
                        last_notice = time.monotonic()
        finally:
            _stop_process(process)


def _tail(path, limit=12000):
    with Path(path).open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        handle.seek(max(0, handle.tell() - limit))
        return handle.read().decode("utf-8", errors="replace")


def _tool_info(executable, directory, cancelled):
    log = directory / "version.log"
    code = _run_process([str(executable), "--version"], directory, log, cancelled, timeout=20)
    version = _tail(log).strip()
    if code or not re.search(r"\bSKESA\s+\d+\.\d+\.\d+", version):
        raise AssemblyError(f"The selected executable is not a working SKESA: {version}")
    info = {"name": "SKESA", "version": version, "executable": str(executable),
            "executable_sha256": file_sha256(executable, cancelled), "bundled_manifest": None}
    manifest_path = executable.parent / "manifest.json"
    if manifest_path.is_file():
        if manifest_path.stat().st_size > 8 * 1024 * 1024:
            raise AssemblyError("SKESA tool manifest is unexpectedly large.")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("tool") != "SKESA":
            raise AssemblyError("SKESA tool directory has an unrelated manifest.")
        verified_executable = False
        for item in manifest.get("files", []):
            relative = Path(item["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise AssemblyError("Unsafe file path in SKESA manifest.")
            path = (executable.parent / relative).resolve()
            if not path.is_relative_to(executable.parent):
                raise AssemblyError("SKESA manifest file leaves the tool directory.")
            if path == executable or path.suffix.lower() == ".dll":
                if not path.is_file() or file_sha256(path, cancelled) != item["sha256"]:
                    raise AssemblyError(f"SKESA tool integrity check failed: {relative}")
                verified_executable |= path == executable
        if not verified_executable:
            raise AssemblyError("SKESA manifest does not identify this executable.")
        info["bundled_manifest"] = manifest
        info["manifest_sha256"] = file_sha256(manifest_path, cancelled)
    return info


def _normalise_pairs(first, second, directory, cancelled, progress):
    """Validate every paired record while writing ASCII-named, wrapped-safe copies."""
    signatures = [file_signature(p) for p in (first, second)]
    accumulators = [QCAccumulator(), QCAccumulator()]
    explicit, count, written = True, 0, 0
    with SequenceReader(first, cancelled) as left, SequenceReader(second, cancelled) as right:
        if left.kind != "fastq" or right.kind != "fastq":
            raise SequenceError("Paired-read assembly requires two FASTQ files.")
        with (directory / "mate1.fastq").open("w", encoding="ascii", newline="\n") as first_out, \
                (directory / "mate2.fastq").open("w", encoding="ascii", newline="\n") as second_out:
            left_iter, right_iter = iter(left), iter(right)
            while True:
                check_cancelled(cancelled)
                a, b = next(left_iter, None), next(right_iter, None)
                if a is None and b is None:
                    break
                if a is None or b is None:
                    raise SequenceError("Read-pair files contain different numbers of records.")
                aid, am = read_pair_identity(a.name)
                bid, bm = read_pair_identity(b.name)
                if aid != bid:
                    raise SequenceError(f"Read-pair identifiers differ at record {count + 1}.")
                if am not in {None, 1} or bm not in {None, 2} or (am is None) != (bm is None):
                    raise SequenceError(f"Read-pair mate indicators are reversed or inconsistent at record {count + 1}.")
                explicit = explicit and am == 1 and bm == 2
                count += 1
                for mate, record, output, qc in ((1, a, first_out, accumulators[0]),
                                                (2, b, second_out, accumulators[1])):
                    qc.add(record)
                    text = f"@pair{count}/{mate}\n{record.sequence}\n+\n{record.quality}\n"
                    output.write(text)
                    written += len(text)
                if count % 10000 == 0:
                    if shutil.disk_usage(directory).free < 64 * 1024 * 1024:
                        raise AssemblyError("Insufficient disk space while preparing paired FASTQ input.")
                    _notify(progress, 20, f"Validated and staged {count:,} complete read pairs.")
            for output in (first_out, second_out):
                output.flush()
                os.fsync(output.fileno())
    if not count:
        raise SequenceError("Read-pair input contains no reads.")
    _notify(progress, 35, "Hashing original and normalized paired FASTQs.")
    inputs = []
    for index, path in enumerate((first, second)):
        digest = file_sha256(path, cancelled)
        if file_signature(path) != signatures[index]:
            raise SequenceError(f"Input changed during staging: {path.name}")
        inputs.append({"path": str(path), "sha256": digest, "size_bytes": signatures[index][0],
                       "mtime_ns": signatures[index][1],
                       "normalized_sha256": file_sha256(directory / f"mate{index + 1}.fastq", cancelled)})
    pairing = {"records_checked": count, "complete_file": True, "sampled": False,
               "explicit_mates": explicit, "status": "verified" if explicit else "unmarked",
               "header_normalization": "Validated corresponding records renamed pairN/1 and pairN/2 only in temporary copies.",
               "normalized_bytes": written}
    return inputs, pairing, [qc.result("fastq", True, None) for qc in accumulators]


def run_skesa(read1, read2, output_dir, *, executable=None, threads=4, memory_gb=8,
              min_contig=200, cancelled=None, progress=None):
    """Assemble explicit short-read pairs into a NEW output directory atomically.

    progress receives (done, total, message). Cancellation raises AnalysisCancelled.
    Existing output directories/files are never replaced, even if empty.
    Failure logs are retained as uniquely named sibling *.failed-*.log files.
    """
    for name, value in (("threads", threads), ("memory_gb", memory_gb), ("min_contig", min_contig)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    check_cancelled(cancelled)
    first, second, destination = (Path(p).expanduser().resolve() for p in (read1, read2, output_dir))
    if not first.is_file() or not second.is_file():
        raise SequenceError("Both paired FASTQ inputs must exist and be regular files.")
    if first.samefile(second):
        raise SequenceError("Read-pair inputs must be different files.")
    if destination.exists():
        raise FileExistsError(f"Assembly output already exists; choose a new directory: {destination}")
    if any(path.is_relative_to(destination) for path in (first, second)):
        raise ValueError("Assembly output must not contain the original inputs.")
    tool = resolve_skesa(executable)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(destination.parent).free < first.stat().st_size + second.stat().st_size + 64 * 1024 * 1024:
        raise AssemblyError("Insufficient free disk space to stage paired-read assembly.")
    started = datetime.now(timezone.utc).isoformat()
    _notify(progress, 0, "Validating native assembler and original paired inputs.")
    with tempfile.TemporaryDirectory(prefix=".wmlstudio-assembly-", dir=destination.parent) as temporary:
        directory = Path(temporary)
        try:
            tool_info = _tool_info(tool, directory, cancelled)
            inputs, pairing, read_qc = _normalise_pairs(first, second, directory, cancelled, progress)
            command = [str(tool), "--reads", "mate1.fastq,mate2.fastq", "--cores", str(threads),
                       "--memory", str(memory_gb), "--min_contig", str(min_contig),
                       "--contigs_out", "contigs.fasta"]
            _notify(progress, 50, "Running native SKESA short-read isolate assembly.")
            log = directory / "skesa.log"
            code = _run_process(command, directory, log, cancelled, progress)
            check_cancelled(cancelled)
            if code:
                raise AssemblyError(f"SKESA exited with code {code}.\n{_tail(log)}")
            contigs = directory / "contigs.fasta"
            if not contigs.is_file() or not contigs.stat().st_size:
                raise AssemblyError("SKESA produced no contigs. Review read quality, coverage and the assembly log.")
            _notify(progress, 90, "Validating every assembled contig and computing exact contiguity metrics.")
            result = inspect_sequence(contigs, cancelled=cancelled)
            if result["kind"] != "fasta" or not result["qc"]["records"]:
                raise AssemblyError("SKESA output is not a nonempty FASTA assembly.")
            provenance = {
                "workflow": "paired_short_read_isolate_assembly", "engine": tool_info,
                "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
                "command": command, "working_directory": "temporary ASCII-named inputs; not retained",
                "threads": threads, "memory_gb": memory_gb,
                "memory_policy": "SKESA sorted-counter budget, not an OS-enforced process memory limit.",
                "min_contig": min_contig, "inputs": inputs, "pairing": pairing,
                "assembly_sha256": result["input_sha256"],
                "read_preprocessing": "Format normalization only. No fastp trimming or filtering was performed.",
                "biological_qc": {"completeness": "not assessed", "contamination": "not assessed",
                                  "QUAST": "not run", "QuickClade": "not run; not bundled with native adapter"},
            }
            notes = [
                "Short-read isolate assembly only; not a metagenome or long-read assembly workflow.",
                "Genome completeness and contamination were not assessed; contiguity is not a purity check.",
                "FASTQ quality assumes Phred+33. Read preprocessing was format normalization, not quality filtering.",
                "The memory setting is SKESA's sorted-counter budget, not an operating-system enforced limit.",
            ]
            if not pairing["explicit_mates"]:
                notes.append("Read IDs match completely but lack explicit mate markers; biological pairing is user-assigned.")
            payload = {"assembly_path": str(destination / "contigs.fasta"), "provenance": provenance,
                       "qc": result["qc"], "read_qc": read_qc, "pairing": pairing, "notes": notes}
            with (directory / "assembly.json").open("w", encoding="utf-8") as report:
                json.dump(payload, report, indent=2, allow_nan=False)
                report.write("\n")
                report.flush()
                os.fsync(report.fileno())
            (directory / "mate1.fastq").unlink()
            (directory / "mate2.fastq").unlink()
            check_cancelled(cancelled)
            if os.name == "nt":
                # Windows directory rename is atomic and refuses an existing target.
                os.rename(directory, destination)
            else:
                # Claim the name without replacing anything; atomically replace
                # only the empty placeholder just created by this operation.
                destination.mkdir()
                try:
                    os.replace(directory, destination)
                except BaseException:
                    destination.rmdir()  # Empty owned placeholder only.
                    raise
            return payload
        except BaseException as error:
            logs = [p for p in (directory / "skesa.log", directory / "version.log") if p.is_file()]
            if logs:
                failure_log = destination.parent / f"{destination.name}.failed-{uuid.uuid4().hex}.log"
                try:
                    with failure_log.open("x", encoding="utf-8") as output:
                        output.write(str(error) + "\n")
                        for log in logs:
                            output.write(f"\n{log.name}\n{_tail(log)}\n")
                    if not isinstance(error, AnalysisCancelled):
                        error.add_note(f"Diagnostic log: {failure_log}")
                except OSError:
                    pass
            raise


def associate_assembly(project, primary_id, mate_id, assembly_result):
    """Attach one assembly to two stable read records without deleting either.

    The first record becomes the typing input, while the second remains an
    explicitly linked read-mate record. The caller must exclude source_kind
    read_mate from automatic isolate counting/typing. All changes are one
    transaction; previous typing results and input locations remain in history.
    """
    if primary_id == mate_id:
        raise ValueError("Assembly requires two distinct sample records.")
    primary, mate = project.get_sample(primary_id), project.get_sample(mate_id)
    provenance = assembly_result.get("provenance", {})
    inputs = provenance.get("inputs", [])
    if len(inputs) != 2 or not all(item.get("sha256") for item in inputs):
        raise ValueError("Assembly provenance must identify both hashed read inputs.")
    read_qc = assembly_result.get("read_qc", [None, None])
    if not isinstance(read_qc, list) or len(read_qc) != 2:
        raise ValueError("Assembly evidence requires two read QC entries.")
    for sample, source in zip((primary, mate), inputs):
        if not sample["input_path"] or Path(sample["input_path"]).resolve() != Path(source["path"]).resolve():
            raise ValueError("Assembly read provenance does not match the selected sample records.")
        linked = sample["metadata"].get("workflow", {}).get("paired_with")
        expected = mate_id if sample["id"] == primary_id else primary_id
        if linked and linked != expected:
            raise ValueError("A selected read record already belongs to another read pair.")
    path = Path(assembly_result["assembly_path"]).resolve()
    expected_hash = provenance.get("assembly_sha256")
    if not path.is_file() or not expected_hash or file_sha256(path) != expected_hash:
        raise ValueError("Assembly FASTA does not match its recorded SHA-256.")
    if inspect_sequence(path)["kind"] != "fasta":
        raise ValueError("Assembly evidence must point to a FASTA, not raw reads.")
    if any(path == Path(sample["input_path"]).resolve() for sample in (primary, mate)):
        raise ValueError("Assembly cannot replace either original read input.")
    result = copy.deepcopy(assembly_result)
    result["assembly_path"] = str(path)
    with project.transaction():
        for sample in (primary, mate):
            if project.get_sample(sample["id"])["input_path"] != sample["input_path"]:
                raise ValueError("A paired sample input changed while the assembly was being attached.")
        project.invalidate_result(primary_id, "Paired-read assembly created; type the assembled contigs.")
        project.set_input_path(primary_id, path)
        project.update_metadata(primary_id, {
            "assembly": result,
            "workflow": {"source_kind": "assembly", "paired_with": mate_id,
                         "read1_path": primary["input_path"], "read2_path": mate["input_path"],
                         "managed_sha256": expected_hash},
        })
        project.update_metadata(mate_id, {
            "workflow": {"source_kind": "read_mate", "paired_with": primary_id,
                         "assembly_path": str(path)},
            "assembly": {"primary_sample_id": primary_id, "assembly_sha256": expected_hash,
                         "read_input": copy.deepcopy(inputs[1]),
                         "read_qc": copy.deepcopy(read_qc[1])},
        })
        for identifier, role in ((primary_id, "primary"), (mate_id, "read_mate")):
            project.record_history(identifier, "paired_assembly_linked", {
                "role": role, "primary_sample_id": primary_id, "mate_sample_id": mate_id,
                "assembly_path": str(path), "assembly_sha256": expected_hash,
                "read_inputs": inputs,
            })
    return primary_id
