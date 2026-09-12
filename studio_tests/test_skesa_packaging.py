import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def staging_module():
    spec = importlib.util.spec_from_file_location("stage_skesa", REPO / "studio_packaging/stage_skesa_runtime.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_recursive_dll_closure_and_system_contracts(tmp_path, monkeypatch):
    module = staging_module()
    target, runtime, system = [tmp_path / name for name in ("target", "ucrt", "system")]
    for directory in (target, runtime, system):
        directory.mkdir()
    executable = target / "skesa.exe"
    executable.write_bytes(b"PE fixture")
    for name in ("boost.dll", "libgcc.dll"):
        (runtime / name).write_bytes(name.encode())
    (system / "KERNEL32.dll").write_bytes(b"system fixture")
    dependencies = {"skesa.exe": ["BOOST.DLL", "kernel32.dll", "api-ms-win-crt-runtime-l1-1-0.dll"],
                    "boost.dll": ["libgcc.dll"], "libgcc.dll": ["boost.dll"]}
    monkeypatch.setattr(module, "imported_dlls", lambda path, objdump: dependencies[path.name])
    staged = module.collect_dlls(executable, runtime, system)
    assert {path.name for path in staged} == {"boost.dll", "libgcc.dll"}
    assert {path.name for path in target.iterdir()} == {"skesa.exe", "boost.dll", "libgcc.dll"}
    assert (target / "boost.dll").read_bytes() == b"boost.dll"


def test_unresolved_dll_fails_instead_of_releasing_broken_tool(tmp_path, monkeypatch):
    module = staging_module()
    executable = tmp_path / "skesa.exe"
    executable.write_bytes(b"fixture")
    monkeypatch.setattr(module, "imported_dlls", lambda *args: ["missing.dll"])
    with pytest.raises(RuntimeError, match="Unresolved native dependency"):
        module.collect_dlls(executable, tmp_path, tmp_path)


def test_build_is_pinned_native_and_requires_real_smoke_assembly():
    script = (REPO / "studio_packaging/build_skesa_windows.sh").read_text()
    assert staging_module().PIN in script
    assert "-DNO_NGS" in script and "UCRT64" in script
    assert "check_skesa.py" in script
    assert script.index("check_skesa.py") < script.index("mv -T")
    workflow = (REPO / ".github/workflows/skesa-windows.yml").read_text()
    assert "windows-latest" in workflow and "workbench-v0.2" in workflow
    assert "workflow_dispatch" in workflow
