"""Cancellable sequential analysis queue; no worker thread accesses the GUI/database."""

import shutil
import tempfile
import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal

from wmlstudio.cgtyping import call_cgassembly
from wmlstudio.identification import cached_scheme, identify_assembly
from wmlstudio.sequence import AnalysisCancelled, SequenceReader, inspect_sequence
from wmlstudio.typing import call_assembly, load_scheme


class AnalysisWorker(QThread):
    sample_started = Signal(str)
    sample_finished = Signal(str, dict)
    sample_failed = Signal(str, str)
    sample_cancelled = Signal(str)
    progress = Signal(int, str)

    def __init__(self, samples, scheme_path=None, max_reads=100000, parent=None,
                 installed_scheme_paths=None):
        super().__init__(parent)
        self.samples = samples
        self.scheme_path = scheme_path
        self.max_reads = max_reads
        self.installed_scheme_paths = list(installed_scheme_paths or [])
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        scheme = None
        scheme_error = None
        needs_global = False
        for sample in self.samples:
            metadata = sample.get('metadata') or {}
            workflow = metadata.get('workflow') or {}
            if (workflow.get('typing_mode', metadata.get('typing_mode', 'manual')) == 'manual'
                    and 'scheme_path' not in workflow and 'scheme_path' not in metadata):
                needs_global = True
                break
        if self.scheme_path and needs_global:
            self.progress.emit(0, "Preparing allele scheme…")
            try:
                scheme = load_scheme(self.scheme_path, cancelled=self.cancel_event.is_set)
            except AnalysisCancelled:
                self.progress.emit(100, "Analysis cancelled")
                return
            except Exception as exc:
                scheme_error = str(exc)
        total = len(self.samples)
        for index, sample in enumerate(self.samples):
            if self.cancel_event.is_set():
                break
            sample_id = sample["id"]
            self.sample_started.emit(sample_id)
            self.progress.emit(int(index / total * 100), f"Analysing {sample['name']} · {index + 1} of {total}")
            try:
                path = Path(sample["input_path"])
                metadata = sample.get("metadata") or {}
                workflow = metadata.get("workflow") or {}
                typing_mode = workflow.get("typing_mode", metadata.get("typing_mode"))
                selected_path = workflow.get("scheme_path", metadata.get("scheme_path"))
                selection_overridden = "scheme_path" in workflow or "scheme_path" in metadata
                # Existing callers without per-sample metadata keep the legacy
                # queue-level selection. New workspaces explicitly select auto.
                if typing_mode is None:
                    typing_mode = "manual" if selected_path or self.scheme_path else "unknown"
                if typing_mode not in {"auto", "manual", "unknown"}:
                    raise ValueError(f"Unknown typing assignment mode: {typing_mode}")
                with SequenceReader(path, cancelled=self.cancel_event.is_set) as reader:
                    is_read = reader.kind == "fastq"
                if not is_read and typing_mode == "manual" and not selection_overridden and scheme_error:
                    raise ValueError(f"Scheme could not be loaded: {scheme_error}")
                result = None
                identification = None
                if not is_read and typing_mode == "auto":
                    locations = self.installed_scheme_paths or ([self.scheme_path] if self.scheme_path else [])
                    identification = identify_assembly(
                        path, locations, cancelled=self.cancel_event.is_set,
                        progress=lambda current, count, message: self.progress.emit(
                            int(index / total * 100), message),
                    )
                    result = identification.pop("typing_result", None)
                elif not is_read and typing_mode == "manual":
                    selected_scheme = (cached_scheme(selected_path, self.cancel_event.is_set)
                                       if selected_path else None if selection_overridden else scheme)
                    if selected_scheme is not None:
                        cg_mode = workflow.get("calling_mode", metadata.get("calling_mode", "auto"))
                        if cg_mode not in {'auto', 'exact', 'full_cds'}:
                            raise ValueError(f'Unknown allele calling mode: {cg_mode}')
                        use_cg = cg_mode == "full_cds" or (
                            cg_mode == "auto" and (
                                len(selected_scheme.loci) > 30
                                or str(selected_scheme.metadata.get("type", "")).lower() == "cgmlst"))
                        caller = call_cgassembly if use_cg else call_assembly
                        kwargs = {'genetic_code': workflow.get('genetic_code',
                                   selected_scheme.metadata.get('genetic_code', 11))} if use_cg else {}
                        result = caller(path, selected_scheme, cancelled=self.cancel_event.is_set,
                                        progress=lambda current, count, message: self.progress.emit(
                                            int(index / total * 100), message), **kwargs)
                if result is None:
                    result = inspect_sequence(path, max_reads=self.max_reads, cancelled=self.cancel_event.is_set)
                    result.update(status="qc_only", st=None, alleles={}, calls=[], scheme=None, scheme_digest=None)
                    result.setdefault("notes", []).append(
                        "Read quality only. Raw reads have not been assembled or typed." if is_read
                        else "Assembly quality only. Select an allele scheme to type this assembly.")
                result["typing_mode"] = typing_mode
                if identification is not None:
                    result["identification"] = identification
                    result.setdefault("notes", []).extend(identification["notes"])
                    if identification["identification_status"] == "assigned":
                        result["organism"] = dict(identification["organism"])
                elif metadata.get("organism"):
                    result["organism"] = dict(metadata["organism"])
                    result["organism_assignment"] = "user supplied"
                result["sample_name"] = sample["name"]
                result["sample_id"] = sample_id
                result["software"] = "WMLSTudio"
                from wmlstudio import __version__
                result["software_version"] = __version__
                self.sample_finished.emit(sample_id, result)
            except AnalysisCancelled:
                self.sample_cancelled.emit(sample_id)
                break
            except Exception as exc:
                self.sample_failed.emit(sample_id, str(exc))
        self.progress.emit(100, "Analysis cancelled" if self.cancel_event.is_set() else "Analysis finished")


class SchemeImportWorker(QThread):
    """Validate and copy a scheme without blocking navigation or modifying its source."""

    imported = Signal(str, int)
    failed = Signal(str)
    progress = Signal(int, str)

    def __init__(self, source, root, parent=None):
        super().__init__(parent)
        self.source, self.root = Path(source), Path(root)
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        staging = None
        try:
            self.progress.emit(0, "Checking allele scheme…")
            scheme = load_scheme(self.source, cancelled=self.cancel_event.is_set)
            destination = self.root / "schemes" / (self.source.name + "_" + scheme.digest[:8])
            if not destination.exists():
                staging = Path(tempfile.mkdtemp(prefix="scheme-", dir=self.root))
                files = [p for p in self.source.iterdir() if p.is_file() and p.suffix.lower() in {".tfa", ".fasta", ".fa", ".fna", ".fas", ".txt", ".tsv", ".json", ".gz", ".bz2"}]
                for index, path in enumerate(files):
                    if self.cancel_event.is_set():
                        raise AnalysisCancelled()
                    if path.is_symlink():
                        raise ValueError("Scheme files must be regular files, not symbolic links.")
                    shutil.copy2(path, staging / path.name)
                    self.progress.emit(int((index + 1) / len(files) * 90), "Copying validated scheme files…")
                checked = load_scheme(staging, cancelled=self.cancel_event.is_set)
                if checked.digest != scheme.digest:
                    raise ValueError("The source scheme changed during import. Try again after the files stop changing.")
                destination.parent.mkdir(parents=True, exist_ok=True)
                staging.rename(destination)
                staging = None
            else:
                checked = load_scheme(destination, cancelled=self.cancel_event.is_set)
                if checked.digest != scheme.digest:
                    raise ValueError("The existing imported scheme has changed. Move that modified copy out of the scheme library before importing this snapshot again.")
            self.imported.emit(str(destination), len(scheme.loci))
            self.progress.emit(100, "Scheme ready")
        except AnalysisCancelled:
            self.progress.emit(100, "Scheme import cancelled")
        except Exception as exc:
            self.failed.emit(str(exc))
        finally:
            if staging is not None:
                shutil.rmtree(staging)
