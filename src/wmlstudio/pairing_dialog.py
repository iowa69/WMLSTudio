"""Review paired files explicitly; filename hints are never evidence of valid pairing."""

import re
from pathlib import Path

from PySide6.QtWidgets import QComboBox, QDialog, QDialogButtonBox, QVBoxLayout

from wmlstudio.ui_common import cell, make_table
from wmlstudio.widgets import label


def mate_hint(path):
    name = Path(path).name.lower().removesuffix(".gz").removesuffix(".bz2")
    name = re.sub(r"\.(?:fastq|fq)$", "", name)
    match = re.search(r"([_.-])r?([12])(?:_001)?$", name)
    return (name[:match.start()], int(match[2])) if match else (name, None)


class PairReadsDialog(QDialog):
    def __init__(self, samples, parent=None):
        super().__init__(parent)
        self.samples, self.assignments = samples, []
        self.combos = []
        self.setWindowTitle("Review paired reads for assembly")
        self.resize(880, 590)
        layout = QVBoxLayout(self)
        layout.addWidget(label("Pair files into isolates", "title"))
        layout.addWidget(label("Choose a mate beside each forward-read sample you want to assemble. A mate row stays 'Not a forward file'. Filename matches are suggestions: every read identifier, mate order, record count and quality string is validated before SKESA starts.", "muted", True))
        table = make_table(["Forward sample", "Input file", "Reverse mate"])
        table.setSortingEnabled(False)
        table.setRowCount(len(samples))
        for row, sample in enumerate(samples):
            table.setItem(row, 0, cell(sample["name"]))
            table.setItem(row, 1, cell(Path(sample["input_path"]).name))
            combo = QComboBox()
            combo.addItem("Not a forward file / QC only", None)
            for mate in samples:
                if mate["id"] != sample["id"]:
                    combo.addItem(mate["name"], mate["id"])
            stem, direction = mate_hint(sample["input_path"])
            suggestions = [mate["id"] for mate in samples if mate_hint(mate["input_path"]) == (stem, 2)]
            if direction == 1 and len(suggestions) == 1:
                combo.setCurrentIndex(combo.findData(suggestions[0]))
            table.setCellWidget(row, 2, combo)
            self.combos.append((sample["id"], combo))
        layout.addWidget(table, 1)
        layout.addWidget(label("One assembly is linked to each forward-sample ID. The reverse file keeps its existing record as a linked mate; no input or prior evidence is deleted. Unpaired files receive QC only. This workflow is for paired short reads, not long-read or metagenome assembly.", "small", True))
        self.feedback = label("Confirm the biological pairing, especially when names were assigned manually.", "small", True)
        layout.addWidget(self.feedback)
        controls = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        controls.button(QDialogButtonBox.StandardButton.Ok).setText("Confirm paired-read assembly")
        controls.accepted.connect(self.accept)
        controls.rejected.connect(self.reject)
        layout.addWidget(controls)

    def accept(self):
        pairs = [{"primary_id": identifier, "mate_id": combo.currentData()}
                 for identifier, combo in self.combos if combo.currentData()]
        used = [identifier for pair in pairs for identifier in pair.values()]
        if not pairs:
            self.feedback.setText("Choose at least one forward/reverse pair, or cancel and run quality checks only.")
            return
        if len(set(used)) != len(used):
            self.feedback.setText("A read file can belong to only one pair. Leave reverse-mate rows set to 'Not a forward file'.")
            return
        self.assignments = pairs
        super().accept()
