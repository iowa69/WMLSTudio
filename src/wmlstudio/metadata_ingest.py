"""Reviewed, transactional epidemiological annotations; never sequence evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from copy import deepcopy
from datetime import date
from pathlib import Path

from wmlstudio.export import (
    _atomic_text,
    _csv_cell,
    ensure_separate_destination,
    protected_input_paths,
)

IDENTITY_FIELDS = {'sample_id', 'sample_name'}
PROTECTED_FIELDS = {'id', 'name', 'input_path', 'input_sha256', 'result', 'results', 'status',
                    'st', 'scheme', 'scheme_digest', 'alleles', 'calls', 'workflow', 'hydra',
                    'characterization', 'assembly', 'organism', 'genus', 'species', 'qc', 'created_at', 'updated_at'}
DEFAULT_FIELDS = ['collection_date', 'ward', 'specimen', 'patient_code', 'location']
CLEAR_TOKEN = '<CLEAR>'


def _digest(metadata):
    return hashlib.sha256(json.dumps(metadata, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def annotation_field(value):
    name = str(value).strip().removeprefix('annotations.').casefold().replace(' ', '_')
    if name in PROTECTED_FIELDS or name in IDENTITY_FIELDS:
        raise ValueError(f'Protected field cannot be edited as epidemiology: {value}')
    if not re.fullmatch(r'[a-z][a-z0-9_-]{0,79}', name):
        raise ValueError(f'Annotation field must be a short name using letters, digits, underscore or hyphen: {value}')
    return name


def annotation_values(sample):
    """Nested annotations override legacy scalar metadata, including explicit nulls."""
    metadata = sample.get('metadata') or {}
    result = {}
    for key, value in metadata.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            try:
                result[annotation_field(key)] = value
            except ValueError:
                pass
    for key, value in (metadata.get('annotations') or {}).items():
        if not isinstance(value, (str, int, float, bool)) and value is not None:
            continue
        try:
            result[annotation_field(key)] = value
        except ValueError:
            pass
    return result


def read_metadata_table(path, *, max_rows=20000, max_columns=100):
    source = Path(path)
    if source.stat().st_size > 20 * 1024 * 1024:
        raise ValueError('Metadata tables must be at most 20 MiB; split larger imports into reviewed batches.')
    with source.open(encoding='utf-8-sig', newline='') as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = '\t' if source.suffix.casefold() in {'.tsv', '.tab'} else ','
        if source.suffix.casefold() not in {'.csv', '.tsv', '.tab'}:
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
            except csv.Error as error:
                raise ValueError('Choose a CSV or TSV metadata table with a header.') from error
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        headers = next(reader, [])
        headers = [value.strip().casefold().replace(' ', '_') for value in headers]
        if not headers or len(headers) > max_columns or len(headers) != len(set(headers)) or any(not header for header in headers):
            raise ValueError('Metadata header must contain unique nonempty columns (at most 100).')
        if not IDENTITY_FIELDS.intersection(headers):
            raise ValueError('Include sample_id; sample_name is an optional explicitly enabled fallback.')
        rows = []
        for number, values in enumerate(reader, 2):
            if not values or not any(values):
                continue
            if len(values) != len(headers):
                raise ValueError(f'Ragged metadata row {number}: expected {len(headers)} cells, found {len(values)}.')
            if len(rows) >= max_rows or any(len(value) > 10000 for value in values):
                raise ValueError('Metadata exceeds the row or cell-size limit.')
            rows.append(dict(zip(headers, values)))
    return rows


def _value(field, value, clear_token):
    text = '' if value is None else str(value).strip()
    if not text:
        return False, None
    if clear_token is not None and text == clear_token:
        return True, None
    if len(text) > 10000:
        raise ValueError(f'{field}: annotation values must be at most 10,000 characters.')
    if field == 'date' or field.endswith('_date'):
        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}', text):
            raise ValueError(f'{field}: use an unambiguous ISO date YYYY-MM-DD, not {text!r}.')
        try:
            date.fromisoformat(text)
        except ValueError as error:
            raise ValueError(f'{field}: invalid calendar date {text!r}.') from error
    return True, text


def preview_metadata(project, rows, *, allow_name_fallback=False, clear_token=None, sample_ids=None):
    samples = project.samples()
    by_id = {sample['id']: sample for sample in samples}
    by_name = {}
    for sample in samples:
        by_name.setdefault(sample['name'], []).append(sample['id'])
    preview = {'format_version': 1, 'rows': [], 'errors': [], 'conflicts': [], 'warnings': [],
               'allow_name_fallback': bool(allow_name_fallback), 'clear_token': clear_token,
               'project_path': str(project.path.resolve())}
    seen = set()
    for number, raw in enumerate(rows, 2):
        row = {str(key).strip().casefold().replace(' ', '_'): value for key, value in raw.items()}
        sid, name = str(row.get('sample_id') or '').strip(), str(row.get('sample_name') or '').strip()
        item = {'row_number': number, 'sample_id': sid or None, 'sample_name': name, 'changes': {},
                'conflicts': {}, 'errors': [], 'warnings': [], 'resolution': 'stable_id'}
        if len(row) != len(raw):
            item['errors'].append('Duplicate normalized column names in this row.')
        if not sid:
            item['resolution'] = 'exact_name'
            matches = by_name.get(name, []) if allow_name_fallback else []
            if not allow_name_fallback:
                item['errors'].append('Missing sample_id; exact-name fallback is not enabled.')
            elif len(matches) != 1:
                item['errors'].append('Exact sample name is ambiguous.' if matches else 'No exact sample-name match.')
            else:
                sid = item['sample_id'] = matches[0]
        elif sid not in by_id:
            item['errors'].append('Unknown sample_id; a supplied ID is never replaced by a guessed name.')
        if sid in by_id:
            if sample_ids is not None and sid not in set(sample_ids):
                item['errors'].append('Isolate is outside this explicitly selected metadata cohort.')
            if sid in seen:
                item['errors'].append('Duplicate target sample_id in this import; review and keep one row per isolate.')
            seen.add(sid)
            sample = by_id[sid]
            item['sample_name'] = sample['name']
            item['expected_metadata_digest'] = _digest(sample.get('metadata') or {})
            if name and name != sample['name']:
                item['warnings'].append('Display name differs; immutable sample_id is authoritative. The name will not be changed.')
            existing = annotation_values(sample)
            fields_seen = set()
            for column, value in row.items():
                if column in IDENTITY_FIELDS:
                    continue
                try:
                    field = annotation_field(column)
                    if field in fields_seen:
                        raise ValueError(f'Duplicate annotation field after normalization: {field}')
                    fields_seen.add(field)
                    supplied, converted = _value(field, value, clear_token)
                    if supplied and existing.get(field) != converted:
                        item['changes'][field] = converted
                        if field in existing and existing[field] not in (None, ''):
                            item['conflicts'][field] = {'old': existing[field], 'new': converted}
                except ValueError as error:
                    item['errors'].append(str(error))
        item['status'] = 'error' if item['errors'] else 'conflict' if item['conflicts'] else 'change' if item['changes'] else 'unchanged'
        preview['rows'].append(item)
        preview['errors'].extend(f'Row {number}: {message}' for message in item['errors'])
        preview['warnings'].extend(f'Row {number}: {message}' for message in item['warnings'])
        if item['conflicts']:
            preview['conflicts'].append({'row_number': number, 'sample_id': sid, 'fields': item['conflicts']})
    preview['changed_samples'] = sum(bool(row['changes']) for row in preview['rows'])
    return preview


def apply_metadata_preview(project, preview, *, accept_conflicts=False):
    if preview.get('format_version') != 1 or preview.get('project_path') != str(project.path.resolve()):
        raise ValueError('This metadata preview does not belong to the current project.')
    if preview.get('errors') or any(row.get('errors') for row in preview.get('rows', [])):
        raise ValueError('Resolve every mapping, date and protected-field error before applying this import.')
    if preview.get('conflicts') and not accept_conflicts:
        raise ValueError('Existing values conflict. Explicitly approve overwrites after reviewing the preview.')
    seen, updates = set(), []
    with project.transaction():
        for row in preview.get('rows', []):
            sid = row['sample_id']
            if sid in seen:
                raise ValueError('Duplicate isolate IDs cannot be applied.')
            seen.add(sid)
            sample = project.get_sample(sid)
            if _digest(sample.get('metadata') or {}) != row['expected_metadata_digest']:
                raise ValueError('Metadata changed after this preview. Preview the import again; nothing was applied.')
            changes = {}
            current_values = annotation_values(sample)
            for field, value in row['changes'].items():
                normalized = annotation_field(field)
                if value is None:
                    if preview.get('clear_token') is None:
                        raise ValueError('Clearing requires an explicitly enabled clear token.')
                    changes[normalized] = None
                else:
                    supplied, converted = _value(normalized, value, None)
                    if not supplied:
                        raise ValueError('Blank values do not authorize clearing an annotation.')
                    changes[normalized] = converted
                if current_values.get(normalized) not in (None, '', changes[normalized]) and not accept_conflicts:
                    raise ValueError('Changing existing annotation values requires explicit conflict approval.')
            if changes:
                metadata = deepcopy(sample.get('metadata') or {})
                metadata['annotations'] = {**metadata.get('annotations', {}), **changes}
                updates.append((sid, metadata, changes))
        for sid, metadata, changes in updates:
            project.set_metadata(sid, metadata)
            project.record_history(sid, 'epidemiology_metadata_reviewed', {'changes': changes, 'source': 'reviewed metadata grid/import'})
    return [sid for sid, _, _ in updates]


def export_metadata(project, path, *, sample_ids=None, delimiter=None):
    records = project.samples()
    ids = {sample['id'] for sample in records} if sample_ids is None else set(sample_ids)
    if not ids <= {sample['id'] for sample in records}:
        raise ValueError('Metadata export contains unknown isolate IDs.')
    ensure_separate_destination(path, [project.path, *protected_input_paths(records)])
    selected = [sample for sample in records if sample['id'] in ids]
    fields = sorted({key for sample in selected for key in annotation_values(sample)})
    delimiter = delimiter or ('\t' if Path(path).suffix.casefold() == '.tsv' else ',')
    with _atomic_text(path, newline='') as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(['sample_id', 'sample_name', *fields])
        for sample in selected:
            values = annotation_values(sample)
            writer.writerow([_csv_cell(sample['id']), _csv_cell(sample['name']), *[_csv_cell(values.get(field)) for field in fields]])
    return Path(path).resolve()
