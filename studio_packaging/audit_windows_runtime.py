"""Fail closed on non-system PE imports missing from a portable Windows bundle.

BLAST and SKESA receive isolated directory-level checks: copies elsewhere in
the application or a build runner's installed VC++ redist do not satisfy them.
The Python/Qt portion uses PyInstaller's managed package/DLL search locations.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

# Windows 11 system components, not a scan of this build machine's System32.
# VC++, Qt, GCC and other redistributable runtimes are deliberately excluded.
SYSTEM_DLLS = set("""
advapi32.dll authz.dll avrt.dll bcrypt.dll bcryptprimitives.dll cabinet.dll cfgmgr32.dll
combase.dll comctl32.dll comdlg32.dll crypt32.dll d2d1.dll d3d11.dll d3d12.dll
d3d9.dll dbghelp.dll dnsapi.dll dcomp.dll dinput8.dll dsound.dll dwmapi.dll
dwrite.dll dxgi.dll dxva2.dll glu32.dll gdi32.dll hid.dll icu.dll imagehlp.dll
icuuc.dll imm32.dll iphlpapi.dll kernel32.dll kernelbase.dll mf.dll mfplat.dll mfreadwrite.dll
mfuuid.dll mpr.dll msvcrt.dll mswsock.dll ncrypt.dll netapi32.dll normaliz.dll
ntdll.dll ole32.dll oleacc.dll oleaut32.dll opengl32.dll powrprof.dll propsys.dll
psapi.dll rpcrt4.dll secur32.dll setupapi.dll shell32.dll shlwapi.dll shcore.dll
uiautomationcore.dll user32.dll userenv.dll usp10.dll uxtheme.dll ucrtbase.dll version.dll winhttp.dll
wininet.dll winmm.dll winscard.dll winspool.drv wintrust.dll wlanapi.dll wldap32.dll ws2_32.dll wsock32.dll
wtsapi32.dll xinput1_4.dll
""".split())


def imported_dlls(path):
    # pefile is a standard Windows PyInstaller build dependency. It parses both
    # normal and delay imports; objdump/system-PATH tools are not needed in CI.
    import pefile

    with pefile.PE(str(path), fast_load=True) as binary:
        binary.parse_data_directories(directories=[
            pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"],
            pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_DELAY_IMPORT"],
        ])
        return sorted({entry.dll.decode("ascii").casefold()
                       for kind in ("DIRECTORY_ENTRY_IMPORT", "DIRECTORY_ENTRY_DELAY_IMPORT")
                       for entry in getattr(binary, kind, [])})


def verify_directory(directory, *, recursive=False, excluded_roots=()):
    directory = Path(directory).resolve()
    excluded_roots = [Path(root).resolve() for root in excluded_roots]
    paths = directory.rglob("*") if recursive else directory.iterdir()
    binaries = [path for path in paths if path.is_file() and path.suffix.lower() in {".exe", ".dll", ".pyd"}
                and not any(path.is_relative_to(root) for root in excluded_roots)]
    available = {path.name.casefold() for path in binaries}
    unresolved, records = [], []
    for binary in binaries:
        imports = imported_dlls(binary)
        missing = [name for name in imports if name not in available and name not in SYSTEM_DLLS
                   and not name.startswith(("api-ms-win-", "ext-ms-win-"))]
        records.append({"binary": binary.relative_to(directory).as_posix(), "imports": imports})
        # Name where the binary actually is. Several copies of a runtime library can
        # sit in one bundle, and "awt.dll: jvm.dll" does not say which one is wrong.
        unresolved.extend(f"{binary.relative_to(directory).as_posix()}: {name}" for name in missing)
    if not binaries:
        raise ValueError(f"No native binaries found in {directory}")
    if unresolved:
        raise ValueError("Missing app-local Windows dependencies (runner redistributables do not count): "
                         + "; ".join(unresolved))
    return records


def audit(bundle):
    bundle = Path(bundle).resolve()
    tls = verify_tls_providers(bundle)
    reports = {
        "blast_isolated": verify_directory(bundle / "_internal/Tools/blast/bin"),
        "skesa_isolated": verify_directory(bundle / "_internal/wmlstudio/resources/tools/skesa"),
        "ska2_isolated": verify_directory(bundle / "_internal/Tools/ska2"),
        "pyinstaller_managed": verify_directory(bundle, recursive=True, excluded_roots=[
            bundle / "_internal/Tools", bundle / "_internal/wmlstudio/resources/tools"]),
    }
    fastqc_jre = bundle / "_internal/Tools/fastqc/jre"
    if fastqc_jre.is_dir():
        reports["fastqc_java_isolated"] = verify_java_runtime(fastqc_jre)
    notices = json.loads((bundle / "_internal/notices/licenses/manifest.json").read_text(encoding="utf-8"))
    runtime = notices.get("native_runtime", [])
    if {item["name"] for item in runtime} != {"msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"}:
        raise ValueError("Native VC++ runtime provenance is incomplete")
    import hashlib
    for item in runtime:
        path = bundle / "_internal" / item["bundle_path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Native runtime differs from its official-wheel provenance: {item['name']}")
    return {"bundle": str(bundle), "status": "passed", "policy": "Windows11 system allowlist, never runner DLL inventory",
            "groups": reports, "native_runtime": runtime, "tls_providers": tls}


def verify_java_runtime(root):
    """JVM DLLs resolve only against their own directory and app-local JRE bin.

    Do not let a DLL in Python, Qt, BLAST, or a runner's redist mask a missing
    Java dependency. HotSpot loads server/jvm.dll explicitly from the JRE.
    """
    root = Path(root).resolve()
    java_bin = root / "bin"
    records = []
    for binary in sorted(root.rglob("*")):
        if not binary.is_file() or binary.suffix.lower() not in {".exe", ".dll"}:
            continue
        available = {path.name.casefold() for folder in (binary.parent, java_bin)
                     for path in folder.iterdir() if path.is_file()}
        # java.exe explicitly loads server/jvm.dll before JVM code loads AWT,
        # fontmanager or other modules that import the already-loaded JVM.
        if (java_bin / "server/jvm.dll").is_file():
            available.add("jvm.dll")
        imports = imported_dlls(binary)
        missing = [name for name in imports if name not in available and name not in SYSTEM_DLLS
                   and not name.startswith(("api-ms-win-", "ext-ms-win-"))]
        if missing:
            raise ValueError(f"Missing isolated Java runtime dependency for {binary.name}: {', '.join(missing)}")
        records.append({"binary": binary.relative_to(root).as_posix(), "imports": imports})
    if not records or not (java_bin / "java.exe").is_file() or not (java_bin / "server/jvm.dll").is_file():
        raise ValueError("Incomplete app-local Windows Java runtime")
    return records


def verify_tls_providers(bundle):
    """Reject optional Qt OpenSSL binaries discovered outside pinned providers."""
    bundle = Path(bundle)
    files = {path.name.casefold(): path for path in bundle.rglob("*") if path.is_file()}
    forbidden = {"qopensslbackend.dll", "libcrypto-3-x64.dll", "libssl-3-x64.dll"} & files.keys()
    if forbidden:
        raise ValueError("Unpinned ambient Qt TLS libraries must not be redistributed: " + ", ".join(sorted(forbidden)))
    required = {"qschannelbackend.dll", "libcrypto-3.dll", "libssl-3.dll"}
    if not required.issubset(files):
        raise ValueError("Windows native Schannel or official Python TLS provider is missing: "
                         + ", ".join(sorted(required - files.keys())))
    return {"qt": "Windows Schannel", "python": "Official interpreter OpenSSL",
            "files": [files[name].relative_to(bundle).as_posix() for name in sorted(required)]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Native dependency closure passed ({sum(map(len, report['groups'].values()))} binary checks)")
