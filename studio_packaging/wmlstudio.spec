# Build on the target OS: uv run pyinstaller --noconfirm studio_packaging/wmlstudio.spec
import json
import os
import sys
from importlib.metadata import distribution
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

root = Path(SPECPATH).parent
sys.path.insert(0, str(root / "studio_packaging"))
sys.path.insert(0, str(root / "src"))
from stage_bio_tools import filter_windows_qt_tls, stage_hydra_database, verify_skesa_bundle
from stage_fastqc import verify as verify_fastqc
from stage_reference_panels import assay_module_imports, panel_summary, verify_bundle_payload
from stage_ska import verify as verify_ska
schemes = root / "src/wmlstudio/resources/schemes"
manifest = schemes / "manifest.json"
if not manifest.is_file() or json.loads(manifest.read_text(encoding="utf-8"))["scheme_count"] == 0:
    raise SystemExit("Stage schemes with studio_scripts/stage_schemes.py before packaging")
notices = root / "studio_packaging/generated_notices"
if not (notices / "manifest.json").is_file():
    raise SystemExit("Stage license texts with studio_packaging/stage_notices.py before packaging")

datas = [
    (str(root / "src/wmlstudio/resources/ui"), "wmlstudio/resources/ui"),
    (str(schemes), "wmlstudio/resources/schemes"),
    (str(root / "studio_packaging/THIRD_PARTY_NOTICES.md"), "notices"),
    (str(notices), "notices/licenses"),
    (str(root / "docs/STUDIO_WINDOWS.md"), "docs"),
    (str(root / "docs/STUDIO_CAPABILITIES.md"), "docs"),
    (str(root / "docs/HYDRA_INTEGRATION.md"), "docs"),
    (str(root / "docs/MICROBIOLOGY_WORKFLOWS.md"), "docs"),
    (str(root / "docs/WORKBENCH_DESIGN.md"), "docs"),
    (str(root / "docs/WORKFLOW_GUIDE.md"), "docs"),
    (str(root / "docs/ORGANISM_MODULES.md"), "docs"),
    (str(root / "docs/THRESHOLDS.md"), "docs"),
    (str(root / "docs/TEST_DATASETS.md"), "docs"),
]
if (root / "LICENSE").is_file():
    datas.append((str(root / "LICENSE"), "notices"))
for package in ("wmlstudio", "PySide6-Essentials", "shiboken6", "pyahocorasick"):
    datas += copy_metadata(package)
for package in ("hydra-amr", "numpy", "pandas", "python-dateutil", "six", "pyrodigal", "archspec", "pyskani"):
    datas += copy_metadata(package)
datas += collect_data_files("hydra_amr")
datas += collect_data_files("archspec")
tools = Path(os.environ.get("WMLSTUDIO_BLAST_ROOT", root / "src/wmlstudio/resources/tools/blast"))
if not (tools / "manifest.json").is_file():
    raise SystemExit("Stage native BLAST+ with studio_packaging/stage_bio_tools.py before packaging")
platform = "windows-x64" if sys.platform == "win32" else "linux-x64"
fastqc = Path(os.environ.get("WMLSTUDIO_FASTQC_ROOT", root / "src/wmlstudio/resources/tools/fastqc"))
verify_fastqc(fastqc, platform)
datas.append((str(fastqc), "Tools/fastqc"))
ska = root / "src/wmlstudio/resources/tools/ska2" / platform
verify_ska(ska, platform)
datas.append((str(ska), "Tools/ska2"))
if json.loads((tools / "manifest.json").read_text())["platform"] != platform:
    raise SystemExit("Staged BLAST+ archive does not match the target build platform")
datas.append((str(tools), "Tools/blast"))
hydra_database = root / "src/wmlstudio/resources/hydra/starter"
if not (hydra_database / "manifest.json").is_file():
    raise SystemExit("The all-in-one build requires the verified NCBI HYDRA starter snapshot")
stage_hydra_database(hydra_database, hydra_database)
datas.append((str(hydra_database), "wmlstudio/resources/hydra/starter"))
characterization_database = root / "src/wmlstudio/resources/characterization/starter"
if not (characterization_database / "manifest.json").is_file():
    raise SystemExit("Stage the independent species/virulence starter with studio_scripts/stage_characterization.py before packaging")
# Re-hash the whole snapshot before it is copied: a truncated or edited panel must
# fail the build, not ship and then read as a negative assay result.
characterization_panel = panel_summary(characterization_database)
print("Bundled characterization panel: " + json.dumps(characterization_panel))
datas.append((str(characterization_database), "wmlstudio/resources/characterization/starter"))
if sys.platform == "win32":
    skesa = root / "src/wmlstudio/resources/tools/skesa"
    verify_skesa_bundle(skesa)
    datas.append((str(skesa), "wmlstudio/resources/tools/skesa"))
    # Native BLAST imports VC++ runtime libraries absent from NCBI's archive.
    # Keep them app-local: a runner's installed redist is not a portable dependency.
    qt_runtime = Path(distribution("PySide6-Essentials").locate_file("PySide6"))
    for name in ("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"):
        source = qt_runtime / name
        if not source.is_file():
            raise SystemExit(f"The official Qt wheel is missing required native runtime {name}")
        datas.append((str(source), "Tools/blast/bin"))

# Practice cohorts and the broad species panel are downloaded to the user's Data
# root on request. A copy staged into the source tree must stop the build.
print("Unbundled reference payload: " + json.dumps(verify_bundle_payload(source for source, _ in datas)))

common = dict(
    pathex=[str(root / "src")],
    binaries=[],
    datas=datas,
    # impl is a namespace package: collecting only 'pyrodigal' skips it.
    # The organism-module assays are reached through a computed __import__, which
    # leaves no IMPORT_NAME opcode for the module scan; name them explicitly or
    # every characterization run in the frozen build raises ModuleNotFoundError.
    hiddenimports=["ahocorasick", *collect_submodules("pyrodigal.impl"), *collect_submodules("pyskani"),
                   *assay_module_imports()],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets"],
    noarchive=False,
)
gui = Analysis([str(root / "studio_packaging/gui_entry.py")], **common)
gui_exe = EXE(
    PYZ(gui.pure), gui.scripts, [], exclude_binaries=True,
    name="WMLSTudio", debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=False, disable_windowed_traceback=False,
)
cli = Analysis([str(root / "studio_packaging/cli_entry.py")], **common)
cli_exe = EXE(
    PYZ(cli.pure), cli.scripts, [], exclude_binaries=True,
    name="WMLSTudio-CLI", debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=True,
)
hydra = Analysis([str(root / "studio_packaging/hydra_entry.py")], **common)
hydra_exe = EXE(
    PYZ(hydra.pure), hydra.scripts, [], exclude_binaries=True,
    name="WMLSTudio-HYDRA", debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=True,
)
if sys.platform == "win32":
    for analysis in (gui, cli, hydra):
        analysis.binaries = filter_windows_qt_tls(analysis.binaries)
        analysis.datas = filter_windows_qt_tls(analysis.datas)
COLLECT(
    gui_exe, cli_exe, hydra_exe, gui.binaries, gui.datas, cli.binaries, cli.datas, hydra.binaries, hydra.datas,
    strip=False, upx=False, name="WMLSTudio",
)
