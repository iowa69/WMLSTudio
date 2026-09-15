"""Measure the portable package per component and refuse a build that outgrows it.

The promise on the download page is one portable ZIP a microbiologist extracts on
a laboratory desktop, and that promise has a number: the archive stays under one
gigabyte. Everything bundled rather than downloaded spends part of that budget, so
it is measured here, attributed to the component that spent it, and gated — rather
than noticed after a release when a 1.4 GB download is already published.

Two ways to measure, both real:

* An archive is read from its own directory table. A ZIP records the compressed
  size of every member, so the number reported is the number the user downloads.
* A frozen bundle directory is deflated file by file. That is a prediction of the
  archive, not the archive itself, and it says so; the authoritative gate is the
  archive step that follows.

Nothing here writes to or modifies the bundle it measures.
"""

from __future__ import annotations

import argparse
import json
import tarfile
import zipfile
import zlib
from pathlib import Path, PurePosixPath

#: The portable download must stay below this. Decimal bytes, so "1 GB" on the
#: download page means what a user's browser will show, not a larger power-of-two.
BUDGET_BYTES = 1_000_000_000
#: Warn while there is still room to act. A build over this is not failed, but the
#: component table is printed with the warning so the growth has a name.
HEADROOM_WARNING = 0.80
#: ZIP per-entry overhead, measured rather than derived from the format: a real
#: 147,572,214-byte Windows package whose members summed to 146,932,272 carried
#: 3,051 entries, which is 210 bytes each for local headers, central directory
#: records, data descriptors and Zip64 extras. Counted so a directory prediction
#: is not quietly optimistic about a bundle with thousands of small files.
ENTRY_OVERHEAD = 210
#: Per-file zlib level 6 came out 1.5% below what PowerShell's Compress-Archive
#: actually produced for that same package, so a directory prediction that ignored
#: the difference would pass a bundle the real archive fails. The margin is larger
#: than the observed gap, in the direction that protects the download promise, and
#: it applies only to the prediction: an archive is measured, never adjusted.
PREDICTION_MARGIN = 0.03

#: Longest prefix wins, so "Tools/fastqc/jre" would be attributed to FastQC even if
#: a JRE row were added later. Every label says what the bytes bought, because a
#: size report that names only directories cannot be used to decide what to cut.
COMPONENTS = (
    ("Tools/blast", "BLAST+ — the nucleotide and translated search"),
    ("Tools/fastqc", "FastQC and its private Java runtime"),
    ("Tools/ska2", "SKA2 — split k-mer SNP distances"),
    ("Tools/fastp", "fastp — read trimming"),
    ("wmlstudio/resources/tools/skesa", "SKESA — paired short-read assembly"),
    ("wmlstudio/resources/tools", "Other staged native tools"),
    ("wmlstudio/resources/schemes", "MLST and cgMLST scheme snapshot"),
    ("wmlstudio/resources/hydra", "AMRFinderPlus core and point-mutation catalogues"),
    ("wmlstudio/resources/characterization", "Species, virulence and organism-module panels"),
    ("wmlstudio/resources", "Other bundled application resources"),
    ("notices", "Licences and corresponding source archives"),
    ("docs", "Offline documentation"),
    ("PySide6", "Qt for Python"),
    ("shiboken6", "Qt for Python"),
)
OTHER = "Python runtime, application code and dependencies"


def component_for(relative):
    """Which reported component a bundle-relative path belongs to."""
    text = PurePosixPath(relative).as_posix()
    best = ""
    label = OTHER
    for prefix, name in COMPONENTS:
        if (text == prefix or text.startswith(prefix + "/")) and len(prefix) > len(best):
            best, label = prefix, name
    return label


def _strip(names):
    """Drop the archive's single top folder and the PyInstaller `_internal` level.

    A path is reported the way a reader thinks about it — `Tools/blast/bin/blastn`
    — instead of `WMLSTudio/_internal/Tools/blast/bin/blastn`. Only a root shared
    by every member is removed, so an archive with two top-level folders keeps
    both and is still attributed correctly.
    """
    roots = {PurePosixPath(name).parts[0] for name in names if PurePosixPath(name).parts}
    root = next(iter(roots)) if len(roots) == 1 else None

    def strip(name):
        parts = PurePosixPath(name).parts
        if root is not None and parts and parts[0] == root:
            parts = parts[1:]
        if parts and parts[0] == "_internal":
            parts = parts[1:]
        return PurePosixPath(*parts).as_posix() if parts else ""

    return strip


def _accumulate(rows):
    """Fold (relative path, stored bytes, original bytes) rows into a component table."""
    totals = {}
    stored = original = count = 0
    for relative, size, raw in rows:
        label = component_for(relative)
        entry = totals.setdefault(label, {"component": label, "bytes": 0,
                                          "uncompressed_bytes": 0, "files": 0})
        entry["bytes"] += size
        entry["uncompressed_bytes"] += raw
        entry["files"] += 1
        stored += size
        original += raw
        count += 1
    table = sorted(totals.values(), key=lambda item: (-item["bytes"], item["component"]))
    for entry in table:
        entry["share"] = round(entry["bytes"] / stored, 6) if stored else 0.0
    return table, stored, original, count


def measure_archive(path):
    """Read a ZIP or tar archive's own directory table; the bytes are not recompressed."""
    path = Path(path)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members = [item for item in archive.infolist() if not item.is_dir()]
            strip = _strip([item.filename for item in members])
            rows = [(strip(item.filename), item.compress_size, item.file_size) for item in members]
        measured = "compressed bytes recorded in the archive's own directory"
    else:
        with tarfile.open(path) as archive:
            members = [item for item in archive.getmembers() if item.isfile()]
        strip = _strip([item.name for item in members])
        # A .tar.gz compresses the whole stream, so no member carries its own
        # compressed size. Members are attributed by their original size and the
        # archive's real byte count is reported as the total, which is the number
        # that is actually gated.
        rows = [(strip(item.name), item.size, item.size) for item in members]
        measured = "the archive's own byte count; per-component shares are uncompressed"
    table, _stored, original, count = _accumulate(rows)
    # The gate is on the file a user downloads, so the archive's own byte count is
    # the total — never the sum of its members, which omits the directory itself.
    return {"target": str(path), "kind": "archive", "measured": measured,
            "total_bytes": path.stat().st_size, "uncompressed_bytes": original,
            "files": count, "components": table}


def _deflate(path, level=6):
    """Compressed size of one file, streamed, without holding it in memory."""
    compressor = zlib.compressobj(level, zlib.DEFLATED, -zlib.MAX_WBITS)
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            size += len(compressor.compress(chunk))
    return size + len(compressor.flush())


def measure_tree(root, *, compress=True, cancelled=None):
    """Measure a frozen bundle directory, optionally predicting its archive size."""
    root = Path(root)
    rows = []
    for path in sorted(root.rglob("*")):
        if cancelled is not None and cancelled():
            raise KeyboardInterrupt("Size measurement was cancelled")
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        raw = path.stat().st_size
        rows.append((relative, _deflate(path) if compress else raw, raw))
    strip = _strip([relative for relative, _, _ in rows])
    table, stored, original, count = _accumulate(
        [(strip(relative), size, raw) for relative, size, raw in rows])
    report = {"target": str(root), "kind": "directory", "uncompressed_bytes": original,
              "files": count, "components": table}
    if compress:
        report["deflate_bytes"] = stored
        report["entry_overhead_bytes"] = count * ENTRY_OVERHEAD
        report["prediction_margin"] = PREDICTION_MARGIN
        report["total_bytes"] = round((stored + count * ENTRY_OVERHEAD) * (1 + PREDICTION_MARGIN))
        report["measured"] = (
            f"deflate-compressed prediction of the archive: per-file zlib level 6, plus "
            f"{ENTRY_OVERHEAD} bytes of ZIP overhead per entry, plus a "
            f"{PREDICTION_MARGIN * 100:.0f}% margin for the packaging compressor. The archive "
            "itself is the authoritative measurement")
    else:
        report["total_bytes"] = stored
        report["measured"] = "bytes on disk, extracted; this is not the download size"
    return report


def measure(target, *, compress=True, cancelled=None):
    target = Path(target)
    if target.is_dir():
        return measure_tree(target, compress=compress, cancelled=cancelled)
    return measure_archive(target)


def check(report, *, budget=BUDGET_BYTES, gate=True):
    """Add the verdict. `budget=None` measures a quantity the budget does not govern."""
    total = report["total_bytes"]
    report["gated"] = bool(gate and budget)
    if not budget:
        # Bytes on disk are not the download, and comparing them with the download
        # budget would report false headroom. The measurement stands on its own.
        report["status"] = "reported"
        return report
    report["budget_bytes"] = budget
    report["headroom_bytes"] = budget - total
    report["used_fraction"] = round(total / budget, 6) if budget else 0.0
    if not gate:
        report["status"] = "reported"
    elif total > budget:
        report["status"] = "over_budget"
    elif total >= budget * HEADROOM_WARNING:
        report["status"] = "approaching_budget"
    else:
        report["status"] = "within_budget"
    return report


def format_report(report):
    lines = [f"Portable package size · {report['target']}",
             f"  measured: {report['measured']}",
             f"  total:    {report['total_bytes']:,} bytes over {report['files']:,} files"]
    if "budget_bytes" in report:
        lines.append(f"  budget:   {report['budget_bytes']:,} bytes "
                     f"({report['used_fraction'] * 100:.1f}% used, "
                     f"{report['headroom_bytes']:,} bytes of headroom)")
    width = max((len(entry["component"]) for entry in report["components"]), default=0)
    for entry in report["components"]:
        lines.append(f"    {entry['component']:<{width}}  {entry['bytes']:>13,} bytes  "
                     f"{entry['share'] * 100:5.1f}%  {entry['files']:>5,} files")
    return "\n".join(lines)


OVER_BUDGET = (
    "The portable package is over its size budget. It is distributed as one ZIP a "
    "user downloads and extracts, so this is a promise to that user, not a build "
    "preference. Move the largest component above to an explicit on-request "
    "download instead of bundling it, or drop it from this platform's package.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="Portable archive, or a frozen bundle directory")
    parser.add_argument("--budget-bytes", type=int, default=BUDGET_BYTES)
    parser.add_argument("--fast", action="store_true",
                        help="Directory only: report bytes on disk without predicting the archive")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if not args.target.exists():
        raise SystemExit(f"Nothing to measure at {args.target}")
    on_disk = args.fast and args.target.is_dir()
    report = check(measure(args.target, compress=not args.fast),
                   budget=None if on_disk else args.budget_bytes, gate=not on_disk)
    print(format_report(report))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if report["status"] == "over_budget":
        raise SystemExit(OVER_BUDGET)
    if report["status"] == "approaching_budget":
        print(f"WARNING: {report['used_fraction'] * 100:.1f}% of the size budget is used. "
              "Decide what stops being bundled before the next addition, not after it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
