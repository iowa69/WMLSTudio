"""Exercise a frozen CLI with its bundled reference data and a native GUI demo."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import platform
import random
import subprocess
import tempfile
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.sequence import file_sha256, iter_sequences
from wmlstudio.typing import load_scheme


def check_hydra(bundle, root, suffix):
    """Exercise the actual frozen worker, BLAST binaries and bundled references."""
    database = bundle / "_internal/wmlstudio/resources/hydra/starter"
    with (database / "nucl/ncbi/meta.tsv").open(encoding="utf-8", newline="") as handle:
        control_id = next(row["seqid"] for row in csv.DictReader(handle, delimiter="\t") if row["gene"] == "blaZ")
    reference = next(record for record in iter_sequences(database / "nucl/ncbi/sequences.fna")
                     if record.name == control_id)
    sample = root / "synthetic_amr_positive.fasta"
    sample.write_text(f">synthetic_reference_gene\n{reference.sequence}\n", encoding="utf-8")
    destination = root / "hydra-smoke"
    environment = os.environ.copy()
    environment["PATH"] = str(bundle / "_internal/Tools/blast/bin") + os.pathsep + environment.get("PATH", "")
    command = [str(bundle / f"WMLSTudio-HYDRA{suffix}"), "--upstream", "run",
               "--db-dir", str(database), "--assembly", str(sample), "--db", "ncbi", "--db", "protein",
               "--outdir", str(destination), "--tmpdir", str(root / "hydra-temporary"),
               "--format", "json", "--prefix", "hydra", "--threads", "2", "--no-banner",
               "--no-mlst", "--no-typing", "--no-heteroresistance", "--no-reads-mlst", "--no-reads-variants",
               "--no-auto-organism", "--no-point-mutations"]
    completed = subprocess.run(command, env=environment, text=True, capture_output=True, timeout=180)
    if completed.returncode:
        raise RuntimeError(f"Frozen HYDRA worker failed:\n{completed.stdout[-4000:]}\n{completed.stderr[-8000:]}")
    report = json.loads((destination / "hydra.json").read_text(encoding="utf-8"))
    hits = report["samples"][0]["hits"]
    positives = [hit for hit in hits if float(hit.get("identity_pct", 0)) == 100
                 and float(hit.get("coverage_pct", 0)) == 100]
    if not positives:
        raise ValueError("Frozen HYDRA did not recover its exact synthetic reference gene.")
    if not any(hit.get("method") in {"BLASTX", "EXACTX", "ALLELEX"} for hit in hits):
        raise ValueError("Frozen HYDRA did not return the expected protein-search evidence.")
    return {"status": "passed", "control": "One synthetic sequence copied from the bundled NCBI catalog",
            "reference_id": reference.name, "engine_version": report["hydra_version"],
            "hits": hits, "database_versions": report["parameters"]["databases"]}


def check_fastqc(cli, root):
    source = root / "original reads é.fastq"
    source.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCAACGT\n+\n!!!!!!!!\n", encoding="ascii")
    before = file_sha256(source)
    destination = root / "fastqc reports é"
    process = subprocess.run([str(cli), "fastqc", str(source), "--output", str(destination),
                              "--threads", "2", "--memory-gb", "1"],
                             capture_output=True, text=True, timeout=120)
    if process.returncode:
        raise RuntimeError(f"Frozen FastQC failed: {process.stdout[-2000:]}\n{process.stderr[-4000:]}")
    report = json.loads((destination / "fastqc-result.json").read_text(encoding="utf-8"))
    if report["version"] != "0.12.1" or report["reports"][0]["metrics"]["Total Sequences"] != "2":
        raise ValueError("Frozen FastQC did not complete the original two-read control")
    if report["reports"][0]["metrics"]["Filename"] != source.name or file_sha256(source) != before:
        raise ValueError("Frozen FastQC changed input bytes or mishandled Unicode identity")
    if not Path(report["provenance"]["java"]).resolve().is_relative_to(cli.parent.resolve()):
        raise ValueError("Frozen FastQC used a Java runtime outside the portable bundle")
    return {"status": "passed", "version": report["version"], "input_sha256": before,
            "control": "Two synthetic full reads, Unicode input/output paths; original bytes unchanged",
            "modules": report["reports"][0]["modules"], "provenance": report["provenance"]}


def check_species(cli, bundle, root):
    references = bundle / "_internal/wmlstudio/resources/characterization/starter"
    source = references / "species/GCF_000016305.1.fna.gz"
    before = file_sha256(source)
    destination = root / "species.json"
    process = subprocess.run([str(cli), "characterize", str(source), "--references", str(references),
                              "--no-virulence", "--output", str(destination)],
                             capture_output=True, text=True, timeout=180)
    if process.returncode:
        raise RuntimeError(f"Frozen pyskani characterization failed: {process.stdout[-2000:]}\n{process.stderr[-4000:]}")
    report = json.loads(destination.read_text(encoding="utf-8"))
    result = report["species_evidence"]
    nearest = result.get("nearest") or {}
    if (result["status"] != "completed" or result["species"] != "pneumoniae"
            or nearest.get("ani", 0) < 99.9 or nearest.get("query_fraction", 0) < .99
            or nearest.get("reference_fraction", 0) < .99):
        raise ValueError(f"Frozen native pyskani did not recognize its reference positive control: {result}")
    if file_sha256(source) != before or report["input_sha256"] != before:
        raise ValueError("Frozen species probe altered the reference input or its identity")
    matching = [path for path in bundle.rglob("_skani*") if path.is_file() and path.suffix in {".pyd", ".so"}
                and file_sha256(path) == result["provenance"]["binary_sha256"]]
    if not matching:
        raise ValueError("Executed pyskani native-module hash does not match the portable bundle")
    return {"status": "passed", "control": "Bundled K. pneumoniae reference self-comparison; not an independent biological validation",
            "input_sha256": before, "species_evidence": result, "native_module": str(matching[0])}


def check_ska_cli(cli, root):
    randomizer = random.Random(18091)
    sequence = "".join(randomizer.choices("ACGT", k=24000))
    mutated = list(sequence)
    for position in (500, 3500, 7500, 14500, 21000):
        mutated[position] = next(base for base in "ACGT" if base != sequence[position])
    paths = []
    for name, bases in (("a", sequence), ("b", sequence), ("c", "".join(mutated))):
        path = root / f"SNP control {name} é.fa"
        path.write_text(f">contig\n{bases}\n", encoding="ascii")
        paths.append(path)
    before = [file_sha256(path) for path in paths]
    process = subprocess.run([str(cli), "ska", *map(str, paths), "--output", str(root / "SNP output é"), "--threads", "2"],
                             capture_output=True, text=True, encoding="utf-8", timeout=120)
    if process.returncode:
        raise RuntimeError(f"Frozen SKA2 adapter failed: {process.stdout[-2000:]}\n{process.stderr[-4000:]}")
    result = json.loads(process.stdout)
    observed = sorted(row["distance"] for row in result["rows"])
    if observed != [0, 5, 5] or [file_sha256(path) for path in paths] != before:
        raise ValueError("Frozen SKA2 adapter did not preserve inputs or reproduce known 0/5-SNP control")
    return {"status": "passed", "control": "Three synthetic24kb assemblies: one identical pair and five isolated SNPs",
            "observed_distances": observed, "engine": result["engine"], "version": result["version"],
            "binary_sha256": result["binary_sha256"], "inputs": result["inputs"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--output", type=Path, default=Path("artifacts/frozen-check.json"))
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    suffix = ".exe" if (bundle / "WMLSTudio.exe").exists() else ""
    cli = bundle / f"WMLSTudio-CLI{suffix}"
    gui = bundle / f"WMLSTudio{suffix}"
    manifest = bundle / "_internal/wmlstudio/resources/schemes/manifest.json"
    snapshot = json.loads(manifest.read_text(encoding="utf-8"))
    if snapshot["scheme_count"] < 1:
        raise ValueError("No bundled schemes")
    scheme = load_scheme(manifest.parent / "sepidermidis")
    profile, sequence_types = next(iter(scheme.profiles.items()))
    if len(sequence_types) != 1:
        raise ValueError("The packaged positive-control profile must have a unique ST")
    expected = str(sequence_types[0])
    with tempfile.TemporaryDirectory(prefix="WMLSTudio frozen smoke ") as temporary:
        root = Path(temporary)
        sample = root / "synthetic_reference_profile.fasta"
        sample.write_text("".join(
            f">synthetic_{locus}\n{scheme.alleles[locus][allele]}\n"
            for locus, allele in zip(scheme.loci, profile, strict=True)), encoding="utf-8")
        destination = root / "typed.json"
        subprocess.run([str(cli), "type", str(sample), "--scheme", str(scheme.path),
                        "--output", str(destination)], check=True, timeout=90)
        result = json.loads(destination.read_text(encoding="utf-8"))["samples"][0]
        if result["st"] != expected or result["status"] != "complete":
            raise ValueError(f"Frozen typing expected complete ST {expected}, observed {result['st']}")
        if result.get("engine_version") != __version__:
            raise ValueError(f"Frozen engine version differs from build source: {result.get('engine_version')} != {__version__}")
        hydra = check_hydra(bundle, root, suffix)
        fastqc = check_fastqc(cli, root)
        species = check_species(cli, bundle, root)
        ska = check_ska_cli(cli, root)
        skesa = None
        if suffix:
            module_spec = importlib.util.spec_from_file_location("check_skesa", Path(__file__).with_name("check_skesa.py"))
            module = importlib.util.module_from_spec(module_spec)
            module_spec.loader.exec_module(module)
            skesa = module.check(bundle / "_internal/wmlstudio/resources/tools/skesa/skesa.exe")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    screenshot = args.output.with_suffix(".png").resolve()
    subprocess.run([str(gui), "--smoke-test", "--demo", "--screenshot", str(screenshot)],
                   check=True, timeout=90)
    if not screenshot.is_file() or screenshot.stat().st_size == 0:
        raise ValueError("Frozen desktop did not produce a screenshot")
    report = {
        "platform": platform.platform(), "bundle": str(bundle),
        "application_version": __version__, "engine_version": result["engine_version"],
        "scheme_count": snapshot["scheme_count"], "scheme_snapshot_sha256": snapshot["snapshot_sha256"],
        "control": "Synthetic assembly generated from one complete bundled profile",
        "expected_st": expected, "observed_st": result["st"],
        "typing": "passed", "desktop_demo": "passed", "screenshot": str(screenshot),
        "hydra_frozen_worker": hydra, "native_skesa": skesa,
        "original_fastqc": fastqc, "native_pyskani": species,
        "native_ska2": ska,
    }
    revision = os.environ.get("GITHUB_SHA") or os.environ.get("WMLSTUDIO_SOURCE_REVISION")
    if revision:
        report["source_revision"] = revision
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
