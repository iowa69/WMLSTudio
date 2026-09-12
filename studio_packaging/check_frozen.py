"""Exercise a frozen CLI with its bundled reference data and a native GUI demo."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import platform
import subprocess
import tempfile
from pathlib import Path

from wmlstudio.sequence import iter_sequences
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
        hydra = check_hydra(bundle, root, suffix)
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
        "scheme_count": snapshot["scheme_count"], "scheme_snapshot_sha256": snapshot["snapshot_sha256"],
        "control": "Synthetic assembly generated from one complete bundled profile",
        "expected_st": expected, "observed_st": result["st"],
        "typing": "passed", "desktop_demo": "passed", "screenshot": str(screenshot),
        "hydra_frozen_worker": hydra, "native_skesa": skesa,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
