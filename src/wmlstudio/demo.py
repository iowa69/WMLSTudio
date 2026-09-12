"""Deterministic synthetic exercise: generated DNA, not clinical specimens."""

import random
from pathlib import Path


def create_demo(root: Path) -> tuple[Path, list[Path]]:
    rng = random.Random(42)
    scheme = root / "schemes" / "practice_7"
    samples = root / "practice_inputs"
    scheme.mkdir(parents=True, exist_ok=True)
    samples.mkdir(parents=True, exist_ok=True)
    loci = ["abcZ", "adk", "aroE", "fumC", "gdh", "pdhC", "pgm"]
    sequences = {}
    for locus in loci:
        seq = "".join(rng.choice("ACGT") for _ in range(420))
        variant = seq[:210] + ({"A": "C", "C": "G", "G": "T", "T": "A"}[seq[210]]) + seq[211:]
        sequences[locus] = (seq, variant)
        (scheme / f"{locus}.tfa").write_text(f">{locus}_1\n{seq}\n>{locus}_2\n{variant}\n", encoding="utf-8")
    profiles = [[1] * 7, [2] + [1] * 6, [2, 2] + [1] * 5, [2] * 7]
    (scheme / "profiles.tsv").write_text(
        "ST\t" + "\t".join(loci) + "\n" + "".join(
            f"{i + 1}\t" + "\t".join(map(str, profile)) + "\n"
            for i, profile in enumerate(profiles)), encoding="utf-8")
    names = ["Practice_A01", "Practice_A02", "Practice_A03", "Practice_B01", "Practice_B02", "Practice_C01", "Practice_partial"]
    profile_indices = [0, 0, 1, 2, 2, 3, 0]
    paths = []
    for name, profile_index in zip(names, profile_indices, strict=True):
        profile = profiles[profile_index]
        used = loci[:-2] if name.endswith("partial") else loci
        text = "".join(f">{name}_{locus}\n{sequences[locus][profile[loci.index(locus)] - 1]}\n" for locus in used)
        path = samples / f"{name}.fasta"
        path.write_text(text, encoding="utf-8")
        paths.append(path)
    return scheme, paths
