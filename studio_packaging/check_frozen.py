"""Exercise a frozen CLI with its bundled reference data and a native GUI demo."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import tempfile
from pathlib import Path

from wmlstudio.typing import load_scheme


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
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
