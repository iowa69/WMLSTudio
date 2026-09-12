"""The bundled CDS probe must never silently import an installed science wheel."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

_spec = importlib.util.spec_from_file_location("check_bundled_cds", Path(__file__).resolve().parents[1]
                                             / "studio_packaging/check_bundled_cds.py")
probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe)


def test_installed_pyrodigal_is_not_an_import_fallback(tmp_path):
    finder = probe.BundledPyrodigal.__new__(probe.BundledPyrodigal)
    finder.internal, finder.archive, finder.origins = tmp_path, SimpleNamespace(toc={}), {}
    with pytest.raises(ModuleNotFoundError, match="absent from frozen"):
        finder.find_spec("pyrodigal.lib")
    assert finder.find_spec("json") is None


def test_namespace_is_loaded_only_when_present_in_frozen_archive(tmp_path):
    finder = probe.BundledPyrodigal.__new__(probe.BundledPyrodigal)
    finder.internal, finder.archive, finder.origins = tmp_path, SimpleNamespace(toc={"pyrodigal.impl": (3, 0, 0)}), {}
    spec = finder.find_spec("pyrodigal.impl")
    module = importlib.util.module_from_spec(spec)
    finder.exec_module(module)
    assert module.__path__ == [str(tmp_path / "pyrodigal/impl")]
    assert finder.origins["pyrodigal.impl"] == "frozen PYZ: pyrodigal.impl"
