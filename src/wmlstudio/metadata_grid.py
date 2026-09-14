"""Separate cohort-scoped epidemiology editor with explicit import review."""

import json

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QInputDialog,
    QTableWidgetItem,
    QVBoxLayout,
)

from wmlstudio.metadata_ingest import (
    CLEAR_TOKEN,
    DEFAULT_FIELDS,
    annotation_field,
    annotation_values,
    apply_metadata_preview,
    export_metadata,
    preview_metadata,
    read_metadata_table,
)
from wmlstudio.ui_common import FlowLayout, make_table
from wmlstudio.widgets import button, label


class MetadataPreviewDialog(QDialog):
    def __init__(self, project, rows, parent=None, *, sample_ids=None, enable_clear=False):
        super().__init__(parent)
        self.project, self.source_rows, self.sample_ids = project, rows, sample_ids
        self.changed_ids = []
        self.setWindowTitle('Review epidemiology changes before applying')
        self.resize(1080, 640)
        layout = QVBoxLayout(self)
        layout.addWidget(label('Nothing is changed until the complete preview is valid and you approve it. Stable IDs are authoritative; duplicate or unmapped isolates block the whole import.', 'small', True))
        controls = FlowLayout()
        self.allow_names = QCheckBox('Allow unique exact-name fallback when ID is blank')
        self.allow_clear = QCheckBox(f'Interpret {CLEAR_TOKEN} as an explicit clear')
        self.allow_clear.setChecked(enable_clear)
        controls.addWidget(self.allow_names)
        controls.addWidget(self.allow_clear)
        layout.addLayout(controls)
        self.table = make_table(['Row', 'Isolate ID', 'Isolate', 'Review state', 'Proposed changes', 'Conflicts / errors / warnings'])
        layout.addWidget(self.table, 1)
        self.summary = label('', 'small', True)
        layout.addWidget(self.summary)
        self.approve_conflicts = QCheckBox('I reviewed the conflicting values and approve these overwrites / clears')
        layout.addWidget(self.approve_conflicts)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        self.apply_button = self.buttons.button(QDialogButtonBox.StandardButton.Apply)
        self.apply_button.clicked.connect(self.apply_changes)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.allow_names.toggled.connect(self.refresh_preview)
        self.allow_clear.toggled.connect(self.refresh_preview)
        self.approve_conflicts.toggled.connect(self.refresh_apply_state)
        self.refresh_preview()

    def refresh_preview(self):
        self.preview = preview_metadata(self.project, self.source_rows, allow_name_fallback=self.allow_names.isChecked(),
            clear_token=CLEAR_TOKEN if self.allow_clear.isChecked() else None, sample_ids=self.sample_ids)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.preview['rows']))
        for index, row in enumerate(self.preview['rows']):
            notes = [*row['errors'], *row['warnings']]
            if row['conflicts']:
                notes.append(json.dumps(row['conflicts'], ensure_ascii=False))
            values = [row['row_number'], row['sample_id'] or 'Unmapped', row['sample_name'], row['status'],
                      json.dumps(row['changes'], ensure_ascii=False), '; '.join(notes)]
            for column, value in enumerate(values):
                self.table.setItem(index, column, QTableWidgetItem(str(value)))
        self.table.setSortingEnabled(True)
        self.summary.setText(f"{self.preview['changed_samples']} isolates would change · {len(self.preview['errors'])} errors · {len(self.preview['conflicts'])} conflicting rows. Blank cells preserve saved values.")
        self.refresh_apply_state()

    def refresh_apply_state(self):
        self.apply_button.setEnabled(not self.preview['errors'] and
            (not self.preview['conflicts'] or self.approve_conflicts.isChecked()))

    def apply_changes(self):
        try:
            self.changed_ids = apply_metadata_preview(self.project, self.preview,
                accept_conflicts=self.approve_conflicts.isChecked())
            self.accept()
        except Exception as error:
            self.summary.setText(str(error))
            self.apply_button.setEnabled(False)


class MetadataGridDialog(QDialog):
    def __init__(self, window, sample_ids):
        super().__init__(window)
        self.window, self.project, self.sample_ids = window, window.project, set(sample_ids)
        self.changed_ids = set()
        self.extra_fields = set()
        self.setWindowTitle('Epidemiology metadata · selected isolates')
        self.resize(1120, 720)
        layout = QVBoxLayout(self)
        layout.addWidget(label('Edit annotations directly; isolate IDs and names are read-only. Dates use YYYY-MM-DD. Empty cells do not delete previous values: use Clear selected cells for explicit removal.', 'small', True))
        controls = FlowLayout()
        controls.addWidget(button('Review and save edits…', self.review_edits, True))
        controls.addWidget(button('Import CSV / TSV…', self.import_table))
        controls.addWidget(button('Export saved metadata…', self.export_table))
        controls.addWidget(button('Add column…', self.add_column))
        controls.addWidget(button('Clear selected cells', self.clear_cells))
        layout.addLayout(controls)
        self.table = make_table(['Isolate ID', 'Isolate'])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked | QAbstractItemView.EditTrigger.EditKeyPressed)
        layout.addWidget(self.table, 1)
        self.feedback = label('', 'small', True)
        layout.addWidget(self.feedback)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.refresh_grid()

    def refresh_grid(self):
        self.samples = [sample for sample in self.project.samples() if sample['id'] in self.sample_ids]
        existing = {key for sample in self.samples for key in annotation_values(sample)}
        self.fields = list(dict.fromkeys([*DEFAULT_FIELDS, *sorted(existing | self.extra_fields)]))
        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(self.fields) + 2)
        self.table.setHorizontalHeaderLabels(['Isolate ID', 'Isolate', *self.fields])
        self.table.setRowCount(len(self.samples))
        for index, sample in enumerate(self.samples):
            metadata = annotation_values(sample)
            values = [sample['id'], sample['name'], *[metadata.get(field) for field in self.fields]]
            for column, value in enumerate(values):
                item = QTableWidgetItem('' if value is None else str(value))
                item.setData(Qt.ItemDataRole.UserRole, sample['id'])
                if column < 2:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(index, column, item)
        self.table.setSortingEnabled(True)
        self.feedback.setText(f'{len(self.samples)} explicitly selected isolates. Original sequences, typing and linked AMR evidence are protected.')

    def edited_rows(self):
        rows = []
        for index in range(self.table.rowCount()):
            row = {'sample_id': self.table.item(index, 0).data(Qt.ItemDataRole.UserRole),
                   'sample_name': self.table.item(index, 1).text()}
            row.update({field: self.table.item(index, column + 2).text() for column, field in enumerate(self.fields)})
            rows.append(row)
        return rows

    def review(self, rows, *, enable_clear=False):
        dialog = MetadataPreviewDialog(self.project, rows, self, sample_ids=self.sample_ids, enable_clear=enable_clear)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.changed_ids.update(dialog.changed_ids)
            self.refresh_grid()
            self.window.refresh()
            self.feedback.setText(f'Saved changes for {len(dialog.changed_ids)} isolates in one transaction. Other metadata and evidence were retained.')

    def review_edits(self):
        self.review(self.edited_rows(), enable_clear=True)

    def import_table(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Preview epidemiology table', '', 'Metadata tables (*.csv *.tsv *.tab)')
        if not path:
            return
        try:
            self.review(read_metadata_table(path))
        except Exception as error:
            self.feedback.setText(str(error))

    def export_table(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Export saved annotations (unsaved edits excluded)', 'epidemiology.tsv', 'TSV (*.tsv);;CSV (*.csv)')
        if not path:
            return
        try:
            self.window.check_output(path)
            export_metadata(self.project, path, sample_ids=self.sample_ids)
            self.feedback.setText('Exported saved metadata for this cohort. Review and save grid edits before exporting those changes.')
        except Exception as error:
            self.feedback.setText(str(error))

    def add_column(self):
        field, accepted = QInputDialog.getText(self, 'Add epidemiology column', 'Annotation name (e.g. admission_date or ward):')
        if not accepted:
            return
        try:
            field = annotation_field(field)
            if field in self.fields:
                return
            self.extra_fields.add(field)
            self.fields.append(field)
            column = self.table.columnCount()
            self.table.insertColumn(column)
            self.table.setHorizontalHeaderItem(column, QTableWidgetItem(field))
            for index in range(self.table.rowCount()):
                self.table.setItem(index, column, QTableWidgetItem(''))
        except ValueError as error:
            self.feedback.setText(str(error))

    def clear_cells(self):
        for item in self.table.selectedItems():
            if item.column() >= 2:
                item.setText(CLEAR_TOKEN)
        self.feedback.setText(f'Selected editable cells are marked {CLEAR_TOKEN}. Review and approve clears before anything is removed.')


def launch_metadata_grid(window, sample_ids=None):
    from wmlstudio.cohort_picker import CohortPickerDialog
    picker = CohortPickerDialog(window.project.samples(), window.project, 'Choose isolates for epidemiology editing',
                                selected_ids=sample_ids, parent=window)
    if picker.exec() != QDialog.DialogCode.Accepted or not picker.selected_ids:
        return set()
    dialog = MetadataGridDialog(window, picker.selected_ids)
    dialog.exec()
    return dialog.changed_ids
