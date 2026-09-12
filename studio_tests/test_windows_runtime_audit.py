"""Portable PE checks cannot be satisfied by a build runner's VC++ redist."""

import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location("audit_windows_runtime", Path(__file__).resolve().parents[1]
                                             / "studio_packaging/audit_windows_runtime.py")
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)


def test_vc_runtime_is_not_mistaken_for_windows_system_component(tmp_path, monkeypatch):
    (tmp_path / "blastn.exe").write_bytes(b"synthetic PE placeholder")
    monkeypatch.setattr(audit, "imported_dlls", lambda path: ["kernel32.dll", "msvcp140.dll"])
    with pytest.raises(ValueError, match="msvcp140.dll"):
        audit.verify_directory(tmp_path)
    (tmp_path / "msvcp140.dll").write_bytes(b"synthetic VC runtime placeholder")
    assert len(audit.verify_directory(tmp_path)) == 2


def test_tool_dlls_elsewhere_do_not_satisfy_isolated_tool_closure(tmp_path, monkeypatch):
    tool = tmp_path / "Tools/blast/bin"
    tool.mkdir(parents=True)
    (tool / "blastn.exe").write_bytes(b"fixture")
    (tmp_path / "msvcp140.dll").write_bytes(b"fixture")
    monkeypatch.setattr(audit, "imported_dlls", lambda path: ["msvcp140.dll"])
    with pytest.raises(ValueError, match="msvcp140.dll"):
        audit.verify_directory(tool)


def test_private_tool_runtime_cannot_satisfy_gui_imports(tmp_path, monkeypatch):
    tool = tmp_path / "Tools"
    tool.mkdir()
    (tmp_path / "desktop.exe").write_bytes(b"fixture")
    (tool / "msvcp140.dll").write_bytes(b"fixture")
    monkeypatch.setattr(audit, "imported_dlls", lambda path: ["msvcp140.dll"])
    with pytest.raises(ValueError, match="msvcp140.dll"):
        audit.verify_directory(tmp_path, recursive=True, excluded_roots=[tool])


def test_api_contracts_and_reviewed_system_libraries_are_allowed(tmp_path, monkeypatch):
    (tmp_path / "tool.exe").write_bytes(b"fixture")
    monkeypatch.setattr(audit, "imported_dlls", lambda path: ["api-ms-win-crt-runtime-l1-1-0.dll",
                                                            "kernel32.dll", "bcrypt.dll"])
    assert len(audit.verify_directory(tmp_path)) == 1


def test_unreviewed_system32_filename_is_not_implicitly_trusted(tmp_path, monkeypatch):
    (tmp_path / "tool.exe").write_bytes(b"fixture")
    monkeypatch.setattr(audit, "imported_dlls", lambda path: ["unreviewed-redist.dll"])
    with pytest.raises(ValueError, match="unreviewed-redist"):
        audit.verify_directory(tmp_path)
