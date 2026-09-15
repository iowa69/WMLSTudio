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
import sys
import tempfile
from pathlib import Path

from wmlstudio import __version__
from wmlstudio.organism_modules import registered_modules
from wmlstudio.sequence import file_sha256, iter_sequences
from wmlstudio.typing import load_scheme

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_size_budget import check as check_budget  # noqa: E402 - this script's own directory
from check_size_budget import format_report as format_budget  # noqa: E402
from check_size_budget import measure_tree  # noqa: E402
from stage_bio_tools import point_mutation_inventory, require_point_mutations  # noqa: E402
from stage_read_tools import verify as verify_fastp  # noqa: E402
from stage_reference_panels import (  # noqa: E402 - resolved from this script's own directory
    UNBUNDLED_NOTE,
    panel_summary,
    sccmec_region_reference,
    unbundled_payload,
)


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


def check_core_reference(bundle):
    """Prove the bundled AMR core reached the package with its point mutations.

    Bundling the core is what makes the first run work with no network, and point
    mutations are part of that core. A package that carried only the acquired-gene
    references would report nothing at all for mutations, which a reader cannot
    distinguish from a negative result — so it fails here rather than shipping.
    """
    from wmlstudio.hydra_runtime import element_counts, organism_catalogue

    store = bundle / "_internal/wmlstudio/resources/hydra/starter"
    inventory = require_point_mutations(point_mutation_inventory(store), store)
    catalogue = organism_catalogue(store)
    counts = element_counts(store)
    if not counts["read"] or counts["total"] <= 0:
        raise ValueError("The frozen bundle's AMR protein reference could not be read at all")
    uncatalogued = sorted(set(catalogue["accepted"]) - set(catalogue["point_mutations"]))
    return {"status": "passed", "store": str(store),
            "dna_point_mutation_catalogues": inventory["dna_catalogues"],
            "protein_point_mutation_organisms": len(catalogue["protein_point_mutations"]),
            "accepted_organisms": len(catalogue["accepted"]),
            "accepted_without_any_catalogue": uncatalogued,
            "elements": counts,
            "boundary": ("An organism with no bundled catalogue is screened for acquired genes "
                         "only. That is an absent catalogue, never a negative mutation result.")}


def check_fastp(bundle, root, suffix):
    """Trim one synthetic pair with the bundle's own fastp, or say why there is none.

    Upstream publishes no Windows binary, so an absent tool is a legitimate
    package state and is reported as `not_bundled` with the reason the application
    itself shows. When the tool is there it is re-verified against its manifest and
    then actually run, because a staged file that cannot execute inside the frozen
    bundle is the failure this check exists to catch. fastp's own JSON report is
    read, never recomputed, and the inputs are hashed before and after.
    """
    directory = bundle / "_internal/Tools/fastp"
    if not (directory / "manifest.json").is_file():
        return {"status": "not_bundled", "reason":
                "This package carries no verified fastp payload, so read trimming is "
                "unavailable in it. Reads remain usable untrimmed, and nothing is "
                "downloaded to substitute for the missing tool."}
    manifest = verify_fastp(directory)
    binary = directory / ("fastp.exe" if suffix else "fastp")
    work = root / "read trimming é"
    work.mkdir()
    reads = []
    for mate, lead in ((1, "ACGT"), (2, "TGCA")):
        path = work / f"original reads é_{mate}.fastq"
        path.write_text("".join(
            f"@pair{index}/{mate}\n{lead}{'ACGTACGTACGTACGTACGTACGTAC'}\n+\n{'I' * 30}\n"
            for index in range(1, 5)), encoding="ascii")
        reads.append(path)
    before = [file_sha256(path) for path in reads]
    report_path = work / "fastp.json"
    command = [str(binary), "--in1", str(reads[0]), "--in2", str(reads[1]),
               "--out1", str(work / "trimmed_1.fastq"), "--out2", str(work / "trimmed_2.fastq"),
               "--json", str(report_path), "--html", str(work / "fastp.html"), "--thread", "1"]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=120)
    if completed.returncode:
        raise RuntimeError(f"Bundled fastp failed: {completed.stdout[-2000:]}\n{completed.stderr[-4000:]}")
    summary = json.loads(report_path.read_text(encoding="utf-8"))["summary"]
    if summary["fastp_version"] != manifest["version"]:
        raise ValueError(f"Bundled fastp reported {summary['fastp_version']}, staged {manifest['version']}")
    kept, supplied = summary["after_filtering"]["total_reads"], summary["before_filtering"]["total_reads"]
    if supplied != 8 or not 0 < kept <= supplied:
        raise ValueError(f"Bundled fastp did not account for its eight synthetic reads: {summary}")
    if [file_sha256(path) for path in reads] != before:
        raise ValueError("Bundled fastp altered the original read files it was given")
    if not Path(binary).resolve().is_relative_to(bundle.resolve()):
        raise ValueError("Bundled fastp resolved to a tool outside the portable package")
    return {"status": "passed", "version": summary["fastp_version"], "platform": manifest["platform"],
            "source_commit": manifest["source_commit"], "input_sha256": before,
            "control": "Four synthetic pairs, Unicode paths; fastp's own report, originals unchanged",
            "reads_supplied": supplied, "reads_kept": kept,
            "boundary": ("Adapter trimming and quality filtering. Not an isolate validation, "
                         "a purity check, a species assignment or a clinical result.")}


def check_size_budget(bundle, output):
    """Measure the package per component and refuse a build that outgrew its budget."""
    report = check_budget(measure_tree(bundle))
    print(format_budget(report))
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if report["status"] == "over_budget":
        raise ValueError(
            f"The portable package predicts {report['total_bytes']:,} bytes compressed, over its "
            f"{report['budget_bytes']:,}-byte download budget. The largest component above is "
            "where the budget went; move it to an explicit on-request download.")
    return report


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
            "input_sha256": before, "species_evidence": result, "native_module": str(matching[0]),
            "characterization_sections": sorted(report)}


def check_reference_panels(bundle):
    """Re-verify the staged reference panels the portable build actually carries.

    The bundled characterization snapshot is re-hashed file by file inside the
    frozen bundle, so a panel that was truncated by packaging fails here rather
    than reading later as a negative assay result. The practice cohorts and the
    broad species panel must be absent: they are the user's own downloads.
    """
    panel = panel_summary(bundle / "_internal/wmlstudio/resources/characterization/starter")
    leaked = unbundled_payload([bundle])
    if leaked:
        raise ValueError(f"{UNBUNDLED_NOTE} This bundle carries {leaked}")
    return {"status": "passed", "characterization": panel,
            "unbundled_payload": [], "unbundled_note": UNBUNDLED_NOTE}


def check_organism_modules(cli, root, panel, sections):
    """Prove the organism-module registry survives freezing, then self-compare SCCmec.

    The assay modules are reached through a computed ``__import__``, which the
    module scan cannot follow, so a bundle that omitted them raises
    ModuleNotFoundError on the first characterization run. ``sections`` is the
    key set of the frozen ``characterize`` report: the registry keys being there
    is the evidence that the import succeeded inside the frozen process. ``panel``
    is the already-verified summary from check_reference_panels.
    """
    expected = sorted(registered_modules())
    absent = [key for key in expected if key not in sections]
    if absent:
        raise ValueError(f"The frozen characterization report carries no evidence block for {absent}; "
                         "the organism-module registry did not load inside the bundle.")
    references = Path(panel["path"])
    registry = {"status": "passed", "registered_modules": expected,
                "control": "Bundled reference self-comparison; not an independent biological validation"}
    if panel["organism_modules"] != "staged":
        return {**registry, "sccmec": {"status": "not_run", "reason": panel["organism_modules_reason"]}}
    usage = subprocess.run([str(cli), "characterize", "--help"], capture_output=True, text=True, timeout=60)
    if "--module" not in usage.stdout:
        return {**registry, "sccmec": {"status": "not_run", "reason":
                "The frozen command line does not offer 'characterize --module', so no frozen SCCmec "
                "control was run here. The assay is reachable from the desktop characterization plan."}}
    source = sccmec_region_reference(references)
    before = file_sha256(source)
    destination = root / "sccmec.json"
    process = subprocess.run([str(cli), "characterize", str(source), "--references", str(references),
                              "--no-species", "--no-virulence", "--module", "sccmec",
                              "--output", str(destination)], capture_output=True, text=True, timeout=300)
    if process.returncode:
        raise RuntimeError(f"Frozen SCCmec typing failed: {process.stdout[-2000:]}\n{process.stderr[-4000:]}")
    evidence = json.loads(destination.read_text(encoding="utf-8"))["sccmec"]
    complexes = {row["name"]: row["state"] for row in (*evidence["ccr_complexes"], *evidence["mec_classes"])}
    if evidence["mecA"] != "detected" or complexes.get("ccr Type 2") != "present":
        raise ValueError(f"Frozen SCCmec typing did not recover its own type IV cassette reference: {evidence}")
    if "IV" not in evidence["candidate_types"]:
        raise ValueError(f"Frozen SCCmec typing proposed {evidence['candidate_types']} for the IVa reference")
    if evidence["official_type"] is not None or evidence["mecC"] != "not_assayed":
        raise ValueError("Frozen SCCmec typing asserted an official type or a mecC result it cannot support")
    if file_sha256(source) != before:
        raise ValueError("Frozen SCCmec typing altered the bundled reference it read")
    return {**registry, "sccmec": {"status": "passed", "reference": str(source), "input_sha256": before,
                                   "complexes": complexes, "candidate_types": evidence["candidate_types"],
                                   "type": evidence["type"], "mecA": evidence["mecA"], "mecC": evidence["mecC"],
                                   "official_type": evidence["official_type"]}}


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
        core = check_core_reference(bundle)
        hydra = check_hydra(bundle, root, suffix)
        fastqc = check_fastqc(cli, root)
        fastp = check_fastp(bundle, root, suffix)
        species = check_species(cli, bundle, root)
        panels = check_reference_panels(bundle)
        modules = check_organism_modules(cli, root, panels["characterization"],
                                         species["characterization_sections"])
        ska = check_ska_cli(cli, root)
        skesa = None
        if suffix:
            module_spec = importlib.util.spec_from_file_location("check_skesa", Path(__file__).with_name("check_skesa.py"))
            module = importlib.util.module_from_spec(module_spec)
            module_spec.loader.exec_module(module)
            skesa = module.check(bundle / "_internal/wmlstudio/resources/tools/skesa/skesa.exe")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    budget = check_size_budget(bundle, args.output.with_name("size-budget.json"))
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
        "bundled_amr_core": core, "hydra_frozen_worker": hydra, "native_skesa": skesa,
        "original_fastqc": fastqc, "original_fastp": fastp, "native_pyskani": species,
        "native_ska2": ska, "reference_panels": panels, "organism_modules": modules,
        "size_budget": budget,
    }
    revision = os.environ.get("GITHUB_SHA") or os.environ.get("WMLSTUDIO_SOURCE_REVISION")
    if revision:
        report["source_revision"] = revision
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
