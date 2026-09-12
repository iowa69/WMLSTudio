"""Stage the recursive native DLL closure and corresponding SKESA source."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

PIN = "c1413581e4f37211892d3c4310d01f3d9a9b3490"


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def imported_dlls(path, objdump="objdump"):
    text = subprocess.check_output([objdump, "-p", str(path)], text=True)
    return sorted(set(re.findall(r"DLL Name:\s*(\S+)", text, re.IGNORECASE)), key=str.casefold)


def collect_dlls(executable, runtime_bin, windows_system=None, objdump="objdump"):
    """Copy non-system imports recursively; unresolved imports fail the build."""
    executable, runtime_bin = Path(executable), Path(runtime_bin)
    target = executable.parent
    windows_system = Path(windows_system or Path(os.environ.get("SystemRoot", "C:/Windows")) / "System32")
    available = {p.name.casefold(): p for p in runtime_bin.glob("*.dll")}
    system = {p.name.casefold() for p in windows_system.glob("*.dll")}
    pending, seen, staged = [executable], set(), []
    while pending:
        binary = pending.pop()
        for name in imported_dlls(binary, objdump):
            key = name.casefold()
            if key in seen:
                continue
            seen.add(key)
            if key.startswith(("api-ms-win-", "ext-ms-win-")):
                continue
            source = available.get(key)
            if source:
                destination = target / source.name
                if source.resolve() != destination.resolve():
                    shutil.copy2(source, destination)
                staged.append(source)
                pending.append(destination)
            elif key not in system:
                raise RuntimeError(f"Unresolved native dependency {name} imported by {binary.name}")
    return staged


def stage(executable, runtime_bin, source, recipe_dir, commit):
    if commit != PIN:
        raise ValueError("This reviewed port only supports the pinned SKESA 2.4.0 commit.")
    executable, source, recipe_dir = map(Path, (executable, source, recipe_dir))
    actual = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if actual != PIN:
        raise ValueError("Upstream source commit does not match the reviewed pin.")
    staged = collect_dlls(executable, runtime_bin)
    target = executable.parent
    shutil.copy2(source / "LICENSE", target / "SKESA-LICENSE.txt")
    shutil.copy2(recipe_dir / "SKESA_WINDOWS.md", target / "README.txt")
    with urllib.request.urlopen("https://www.gnu.org/licenses/agpl-3.0.txt", timeout=30) as response:
        license_text = response.read(100000)
    if hashlib.sha256(license_text).hexdigest() != "0d96a4ff68ad6d4b6f1f30f713b18d5184912ba8dd389f86aa7710db079abcb0":
        raise ValueError("AGPL license text checksum changed; review before redistribution.")
    (target / "AGPL-3.0.txt").write_bytes(license_text)
    with tarfile.open(target / "skesa-2.4.0-wmlstudio-source.tar.gz", "w:gz") as archive:
        for path in sorted(source.rglob("*")):
            if path.is_file() and ".git" not in path.relative_to(source).parts:
                archive.add(path, arcname="skesa/" + path.relative_to(source).as_posix())
        for name in ("build_skesa_windows.sh", "skesa_windows.patch", "stage_skesa_runtime.py",
                     "check_skesa.py", "SKESA_WINDOWS.md"):
            archive.add(recipe_dir / name, arcname="wmlstudio-port/" + name)
    packages = subprocess.check_output(["pacman", "-Q"], text=True)
    (target / "build-packages.txt").write_text(packages, encoding="utf-8")
    licenses = Path(runtime_bin).parent / "share" / "licenses"
    if licenses.is_dir():
        shutil.copytree(licenses, target / "dependency-licenses")
    manifest = {
        "tool": "SKESA", "version": "2.4.0", "platform": "windows-x64-ucrt",
        "source_url": "https://github.com/ncbi/SKESA", "source_commit": PIN,
        "source_archive": "skesa-2.4.0-wmlstudio-source.tar.gz",
        "patch_sha256": sha256(recipe_dir / "skesa_windows.patch"),
        "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
        "build_flags": "-std=c++14 -DNO_NGS -DBOOST_ALL_DYN_LINK -O3 -msse4.2 -pthread",
        "cpu_requirement": "x86-64 with SSE4.2",
        "license": "NCBI public-domain portions plus AGPL-3.0-or-later GATB portions; see SKESA-LICENSE.txt",
        "runtime_dependencies": [p.name for p in staged],
        "files": [{"path": p.relative_to(target).as_posix(), "sha256": sha256(p),
                   "size": p.stat().st_size} for p in sorted(target.rglob("*")) if p.is_file()],
    }
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "runtime-bin", "source", "recipe-dir", "commit"):
        parser.add_argument("--" + name, required=True)
    stage(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
