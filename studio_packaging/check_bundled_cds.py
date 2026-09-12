"""Test frozen Pyrodigal bytecode/PYDs using a separate same-version interpreter.

This is a bundled-runtime probe, not the actual frozen desktop process. No
installed Pyrodigal package is permitted as an import fallback. The desktop's
own launch is checked separately by check_frozen.py.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.abc
import importlib.machinery
import importlib.util
import json
import os
import random
import sys
from pathlib import Path


class BundledPyrodigal(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def __init__(self, executable, internal):
        from PyInstaller.archive.readers import CArchiveReader
        archive = CArchiveReader(str(executable))
        name = next(name for name, entry in archive.toc.items() if entry[-1] == "z")
        self.archive = archive.open_embedded_archive(name)
        self.internal = Path(internal).resolve()
        self.origins = {}

    def find_spec(self, fullname, path=None, target=None):
        if fullname != "pyrodigal" and not fullname.startswith("pyrodigal."):
            return None
        if fullname in self.archive.toc:
            kind = self.archive.toc[fullname][0]
            self.origins[fullname] = "frozen PYZ: " + fullname
            return importlib.util.spec_from_loader(fullname, self, is_package=kind in (1, 3))
        directory = self.internal.joinpath(*fullname.split(".")[:-1])
        spec = importlib.machinery.PathFinder.find_spec(fullname, [str(directory)])
        if spec is None or not spec.origin or not Path(spec.origin).resolve().is_relative_to(self.internal):
            raise ModuleNotFoundError(f"Required Pyrodigal module is absent from frozen bundle: {fullname}")
        self.origins[fullname] = str(Path(spec.origin).resolve())
        return spec

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        name = module.__name__
        kind = self.archive.toc[name][0]
        location = self.internal.joinpath(*name.split("."))
        module.__file__ = str(location / "__init__.py" if kind in (1, 3) else location.with_suffix(".py"))
        if kind in (1, 3):
            module.__path__ = [str(location)]
        if kind != 3:
            exec(self.archive.extract(name), module.__dict__)


def check(bundle):
    bundle = Path(bundle).resolve()
    suffix = ".exe" if os.name == "nt" else ""
    executable = bundle / f"WMLSTudio{suffix}"
    if any(name == "pyrodigal" or name.startswith("pyrodigal.") for name in sys.modules):
        raise RuntimeError("Run this probe in a fresh interpreter without pre-imported Pyrodigal")
    finder = BundledPyrodigal(executable, bundle / "_internal")
    handles = []
    if os.name == "nt":
        for relative in ("_internal", "_internal/PySide6"):
            handles.append(os.add_dll_directory(str(bundle / relative)))
    sys.meta_path.insert(0, finder)
    try:
        import pyrodigal
        rng = random.Random(12)
        sense = ["GCT", "GCA", "AAA", "GAA", "GAT", "CTG", "ATT", "GGT", "TTC", "TAT"]
        expected = "ATG" + "".join(rng.choice(sense) for _ in range(400)) + "TAA"
        sequence = "TAA" * 100 + expected + "TAA" * 100
        bins = pyrodigal.MetagenomicBins([item for item in pyrodigal.METAGENOMIC_BINS
                                        if item.training_info.translation_table == 11])
        caller = pyrodigal.GeneFinder(meta=True, closed=False, metagenomic_bins=bins)
        genes = caller.find_genes(sequence.encode("ascii"))
        exact = [gene for gene in genes if gene.sequence().upper() == expected
                 and not gene.partial_begin and not gene.partial_end and gene.translation_table == 11]
        if len(exact) != 1:
            raise ValueError("Bundled native CDS caller did not recover the complete synthetic ground-truth CDS")
        if "pyrodigal.lib" not in finder.origins or "pyrodigal.impl.generic" not in finder.origins:
            raise ValueError("Native Pyrodigal implementation provenance was not captured")
        return {"status": "passed", "method": "Separate interpreter loading the desktop PYZ and bundled PYDs; no installed-Pyrodigal fallback",
                "bundle": str(bundle), "pyrodigal_version": pyrodigal.__version__, "genetic_code": 11,
                "input_bases": len(sequence), "expected_cds_bases": len(expected), "matching_complete_cds": len(exact),
                "total_predictions": len(genes), "expected_cds_sha256": hashlib.sha256(expected.encode()).hexdigest(),
                "module_origins": finder.origins}
    finally:
        sys.meta_path.remove(finder)
        for handle in handles:
            handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = check(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
