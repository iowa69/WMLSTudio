"""Cancellable sequential analysis queue; no worker thread accesses the GUI/database."""

import shutil
import tempfile
import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal

from wmlstudio.sequence import AnalysisCancelled, SequenceReader, inspect_sequence
from wmlstudio.typing import call_assembly, load_scheme


class AnalysisWorker(QThread):
    sample_started = Signal(str)
    sample_finished = Signal(str, dict)
    sample_failed = Signal(str, str)
    sample_cancelled = Signal(str)
    progress = Signal(int, str)

    def __init__(self, samples, scheme_path=None, max_reads=100000, parent=None):
        super().__init__(parent)
        self.samples = samples
        self.scheme_path = scheme_path
        self.max_reads = max_reads
        self.cancel_event = threading.Event()

    def cancel(self):
        self.cancel_event.set()

    def run(self):
        scheme = None
        scheme_error = None
        if self.scheme_path:
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
                with SequenceReader(path, cancelled=self.cancel_event.is_set) as reader:
                    is_read = reader.kind == "fastq"
                if not is_read and scheme_error:
                    raise ValueError(f"Scheme could not be loaded: {scheme_error}")
                if not is_read and scheme is not None:
                    result = call_assembly(path, scheme, cancelled=self.cancel_event.is_set)
                else:
                    result = inspect_sequence(path, max_reads=self.max_reads, cancelled=self.cancel_event.is_set)
                    result.update(status="qc_only", st=None, alleles={}, calls=[], scheme=None, scheme_digest=None)
                    result.setdefault("notes", []).append(
                        "Read quality only. Raw reads have not been assembled or typed." if is_read
                        else "Assembly quality only. Select an allele scheme to type this assembly.")
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
