# Build on the target OS: uv run pyinstaller --noconfirm studio_packaging/wmlstudio.spec
import json
from pathlib import Path

from PyInstaller.utils.hooks import copy_metadata

root = Path(SPECPATH).parent
schemes = root / "src/wmlstudio/resources/schemes"
manifest = schemes / "manifest.json"
if not manifest.is_file() or json.loads(manifest.read_text(encoding="utf-8"))["scheme_count"] == 0:
    raise SystemExit("Stage schemes with studio_scripts/stage_schemes.py before packaging")
notices = root / "studio_packaging/generated_notices"
if not (notices / "manifest.json").is_file():
    raise SystemExit("Stage license texts with studio_packaging/stage_notices.py before packaging")

datas = [
    (str(schemes), "wmlstudio/resources/schemes"),
    (str(root / "studio_packaging/THIRD_PARTY_NOTICES.md"), "notices"),
    (str(notices), "notices/licenses"),
    (str(root / "docs/STUDIO_WINDOWS.md"), "docs"),
    (str(root / "docs/STUDIO_CAPABILITIES.md"), "docs"),
    (str(root / "docs/HYDRA_INTEGRATION.md"), "docs"),
]
if (root / "LICENSE").is_file():
    datas.append((str(root / "LICENSE"), "notices"))
for package in ("wmlstudio", "PySide6-Essentials", "shiboken6", "pyahocorasick"):
    datas += copy_metadata(package)

common = dict(
    pathex=[str(root / "src")],
    binaries=[],
    datas=datas,
    hiddenimports=["ahocorasick"],
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
COLLECT(
    gui_exe, cli_exe, gui.binaries, gui.datas, cli.binaries, cli.datas,
    strip=False, upx=False, name="WMLSTudio",
)
