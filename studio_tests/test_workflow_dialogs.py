from PySide6.QtWidgets import QDialog

from wmlstudio.workflow_dialogs import BatchAssignmentDialog, ImportSamplesDialog


def test_assignment_dialog_applies_to_selected_files_and_returns_copy_options(qtbot, tmp_path):
    paths = [tmp_path / "one.fasta", tmp_path / "two.fasta"]
    dialog = ImportSamplesDialog(paths, [("E. coli MLST", tmp_path / "scheme")])
    qtbot.addWidget(dialog)
    dialog.table.selectRow(0)
    dialog.mode.setCurrentIndex(dialog.mode.findData("manual"))
    dialog.genus.setText("Escherichia")
    dialog.species.setText("coli")
    dialog.scheme.setCurrentIndex(1)
    dialog.apply_selected()
    assert dialog.assignments[0]["typing_mode"] == "manual"
    assert dialog.assignments[0]["scheme_path"] == str(tmp_path / "scheme")
    assert dialog.assignments[1]["typing_mode"] == "auto"
    dialog.storage_root.setText(str(tmp_path / "managed"))
    dialog.append_st.setChecked(True)
    dialog.accept()
    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.options == {"managed": True, "storage_root": str(tmp_path / "managed"), "append_st": True}
    assert not (tmp_path / "managed").exists()
    assert not any(path.exists() for path in paths)


def test_manual_assignment_requires_genus_and_unknown_mode_clears_scheme(qtbot, tmp_path):
    dialog = ImportSamplesDialog([tmp_path / "sample.fasta"], [tmp_path / "scheme"])
    qtbot.addWidget(dialog)
    dialog.mode.setCurrentIndex(dialog.mode.findData("manual"))
    dialog.apply_all()
    assert "Enter a genus" in dialog.feedback.text()
    assert dialog.assignments[0]["typing_mode"] == "auto"
    dialog.genus.setText("Listeria")
    dialog.scheme.setCurrentIndex(1)
    dialog.apply_all()
    dialog.mode.setCurrentIndex(dialog.mode.findData("unknown"))
    dialog.apply_all()
    assert dialog.assignments[0]["typing_mode"] == "unknown"
    assert dialog.assignments[0]["scheme_path"] is None
    assert dialog.assignments[0]["genus"] == ""


def test_batch_assignment_keeps_stable_ids_and_existing_values(qtbot):
    samples = [{"id": "stable", "name": "Isolate", "input_path": "input.fasta", "metadata": {
        "organism": {"genus": "Klebsiella", "species": "pneumoniae"},
        "workflow": {"typing_mode": "manual", "scheme_path": "/scheme"},
    }}]
    dialog = BatchAssignmentDialog(samples)
    qtbot.addWidget(dialog)
    assert dialog.assignments[0]["sample_id"] == "stable"
    assert dialog.assignments[0]["genus"] == "Klebsiella"
    assert not dialog.managed.isChecked()
    assert samples[0]["metadata"]["organism"]["genus"] == "Klebsiella"
