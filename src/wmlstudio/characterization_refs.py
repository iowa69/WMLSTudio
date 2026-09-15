"""Explicit, immutable local reference provisioning for genomic characterization.

Every panel comes from a named public repository pinned to one commit, and each
source's own licence is fetched and hashed beside its data. Sequence accessions
and taxonomy labels are retained, not presented as exhaustive databases. No
download function is called by analysis or application startup.
"""

from __future__ import annotations

import csv
import hashlib
import http.client
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from .sequence import check_cancelled, file_sha256

#: A reference fetch is retried this many times before the staging run gives up.
FETCH_ATTEMPTS = 4
#: Replies worth retrying: rate limiting and the transient server-side failures.
RETRYABLE_STATUS = frozenset({408, 425, 429, 500, 502, 503, 504})


def _wait_before_retry(attempt, retry_after, cancelled=None):
    """Back off, honouring the server's own Retry-After when it sends one.

    The wait is bounded so a provider asking for an hour cannot stall a build, and
    it is broken into short sleeps so cancelling stays responsive.
    """
    delay = min(2.0 * (2 ** attempt), 30.0)
    if retry_after:
        try:
            delay = max(delay, min(float(retry_after), 60.0))
        except (TypeError, ValueError):
            pass
    deadline = time.monotonic() + delay
    while time.monotonic() < deadline:
        check_cancelled(cancelled)
        time.sleep(0.25)


@dataclass(frozen=True)
class ReferenceSource:
    """One pinned upstream repository, its revision and its licence."""

    key: str
    repository: str
    revision: str
    base: str
    license_id: str
    license_path: str = 'LICENSE'


KLEBORATE_REVISION = '550ce22a2c01c76064f4dabf403704ee2293356e'
SCCMEC_REVISION = 'b901cc618be8eb17284ccb0cf6ef9ee428d909c3'
KAPTIVE_REVISION = 'b3856eac6e76b3017aa993319da2a8ea967a1ba0'
SOURCES = {
    'kleborate': ReferenceSource('kleborate', 'https://github.com/klebgenomics/Kleborate', KLEBORATE_REVISION,
                                 f'https://raw.githubusercontent.com/klebgenomics/Kleborate/{KLEBORATE_REVISION}/',
                                 'GPL-3.0-or-later'),
    # rpetit3/sccmec v1.2.0; the maintained successor to staphopia-sccmec.
    'sccmec': ReferenceSource('sccmec', 'https://github.com/rpetit3/sccmec', SCCMEC_REVISION,
                              f'https://raw.githubusercontent.com/rpetit3/sccmec/{SCCMEC_REVISION}/',
                              'MIT'),
    # Kaptive v2.0.9; only the wzi/wzc marker database is staged, not the K/O loci.
    'kaptive': ReferenceSource('kaptive', 'https://github.com/klebgenomics/Kaptive', KAPTIVE_REVISION,
                               f'https://raw.githubusercontent.com/klebgenomics/Kaptive/{KAPTIVE_REVISION}/',
                               'GPL-3.0-or-later'),
}
SOURCE_BASE = SOURCES['kleborate'].base
TAXONOMY_SOURCE = SOURCE_BASE + 'kleborate/modules/enterobacterales__species/test_enterobacterales__species.py'
GENOMES = (
    ('GCF_000016305.1', 'Klebsiella', 'pneumoniae', ''),
    ('GCF_000492415.1', 'Klebsiella', 'quasipneumoniae', 'quasipneumoniae'),
    ('GCF_000492795.1', 'Klebsiella', 'quasipneumoniae', 'similipneumoniae'),
    ('GCF_000019565.1', 'Klebsiella', 'variicola', 'variicola'),
    ('GCF_002806645.1', 'Klebsiella', 'variicola', 'tropica'),
    ('GCF_000523395.1', 'Klebsiella', 'quasivariicola', ''),
    ('GCF_016804125.1', 'Klebsiella', 'africana', ''),
    ('GCF_000215745.1', 'Klebsiella', 'aerogenes', ''),
    ('GCF_000733495.1', 'Klebsiella', 'grimontii', ''),
    ('GCF_000240325.1', 'Klebsiella', 'michiganensis', ''),
    ('GCF_000247855.1', 'Klebsiella', 'oxytoca', ''),
    ('GCF_000648315.1', 'Klebsiella', 'planticola', ''),
    ('GCF_000005845.2', 'Escherichia', 'coli', ''),
    ('GCF_003937345.1', 'Citrobacter', '', ''),
    ('GCF_004010735.1', 'Salmonella', '', ''),
)
KPSC_SPECIES = {'pneumoniae', 'quasipneumoniae', 'variicola', 'quasivariicola', 'africana'}
VIRULENCE_LOCI = {
    'ybt': ('klebsiella__ybst', ('ybtS', 'ybtX', 'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1', 'ybtU', 'ybtT', 'ybtE', 'fyuA')),
    'clb': ('klebsiella__cbst', ('clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ')),
    'iuc': ('klebsiella__abst', ('iucA', 'iucB', 'iucC', 'iucD', 'iutA')),
    'iro': ('klebsiella__smst', ('iroB', 'iroC', 'iroD', 'iroN')),
    'rmp': ('klebsiella__rmst', ('rmpA', 'rmpD', 'rmpC')),
    'rmpA2': ('klebsiella__rmpa2', ('rmpA2',)),
}
# rmpA2 has no upstream profile table and therefore never acquires a locus ST.
PROFILE_LOCI = ('ybt', 'clb', 'iuc', 'iro', 'rmp')
SCCMEC_ENGINE_DEFAULTS = {'target_min_pident': 90, 'target_min_coverage': 80,
                          'region_min_pident': 85, 'region_min_coverage': 83}
# mecC is verifiably absent from the pinned upstream target panel; it is never
# silently treated as absent, and the panel is not extended with our own record.
SCCMEC_NOT_ASSAYED = ('mecC',)
_SAFE_NAME = re.compile(r'[A-Za-z0-9_.-]+')
_DIGEST_KEYS = {
    1: ('format_version', 'source_revision', 'species', 'virulence', 'files'),
    2: ('format_version', 'sources', 'species', 'virulence', 'locus_profiles', 'sccmec', 'capsule', 'files'),
}
_V2_SECTIONS = ('sources', 'locus_profiles', 'sccmec', 'capsule')


def bundled_reference_root():
    """Locate an installed starter without downloading or validating on UI refresh."""
    roots = []
    if getattr(sys, 'frozen', False):
        roots.append(Path(getattr(sys, '_MEIPASS', Path(sys.executable).parent))
                     / 'wmlstudio' / 'resources' / 'characterization' / 'starter')
    roots.append(Path(__file__).resolve().parent / 'resources' / 'characterization' / 'starter')
    return next((root for root in roots if (root / 'manifest.json').is_file()), None)


def reference_digest(manifest):
    """Fingerprint the sections a given format version is defined to cover.

    Format 1 snapshots keep their original fingerprint byte for byte; format 2
    additionally covers the per-source pins and the organism-module sections, so
    a stripped section can never pass unnoticed.
    """
    keys = _DIGEST_KEYS[manifest['format_version']]
    content = {key: manifest[key] for key in keys}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _fetch(source_key, relative, target, cancelled=None, *, source_roots=None):
    check_cancelled(cancelled)
    source = SOURCES[source_key]
    url = source.base + relative
    local_root = (source_roots or {}).get(source_key)
    if local_root is not None:
        original = Path(local_root) / relative
        if original.is_symlink() or not original.is_file():
            raise ValueError(f'Reference source must be a regular file: {original}')
        handle = original.open('rb')
        return _write(handle, url, source_key, relative, target, cancelled)
    # Staging fetches hundreds of files in a row, so a single dropped connection
    # or rate-limit reply used to lose the whole snapshot. Retry a transient
    # failure a few times; a refusal that will not change is raised immediately.
    request = urllib.request.Request(url, headers={'User-Agent': 'WMLSTudio-reference-provisioning'})
    for attempt in range(FETCH_ATTEMPTS):
        check_cancelled(cancelled)
        try:
            handle = urllib.request.urlopen(request, timeout=30)
        except urllib.error.HTTPError as error:
            if error.code not in RETRYABLE_STATUS or attempt == FETCH_ATTEMPTS - 1:
                raise
            _wait_before_retry(attempt, error.headers.get('Retry-After'), cancelled)
            continue
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError):
            if attempt == FETCH_ATTEMPTS - 1:
                raise
            _wait_before_retry(attempt, None, cancelled)
            continue
        try:
            return _write(handle, url, source_key, relative, target, cancelled)
        except (TimeoutError, ConnectionError, http.client.IncompleteRead) as error:
            # A short read leaves a partial file, and the writer creates the target
            # exclusively, so the attempt must be cleared before the next one.
            target.unlink(missing_ok=True)
            if attempt == FETCH_ATTEMPTS - 1:
                raise ValueError(f'Reference download was cut short: {relative} ({error})') from error
            _wait_before_retry(attempt, None, cancelled)
    raise ValueError(f'Reference file could not be downloaded: {relative}')


def _write(handle, url, source_key, relative, target, cancelled):
    target.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with handle, target.open('xb') as output:
        while chunk := handle.read(1024 * 1024):
            check_cancelled(cancelled)
            total += len(chunk)
            if total > 64 * 1024 * 1024:
                raise ValueError(f'Reference file exceeds the 64 MiB download bound: {relative}')
            output.write(chunk)
    return {'path': target.name, 'bytes': total, 'sha256': file_sha256(target, cancelled),
            'source': source_key, 'source_url': url}


def _derived_entry(stage, relative, cancelled, **extra):
    target = stage / relative
    return dict({'path': relative, 'bytes': target.stat().st_size,
                 'sha256': file_sha256(target, cancelled)}, **extra)


def _fasta_records(path):
    """Read (identifier line, sequence) pairs without collapsing duplicate names."""
    header, lines = None, []
    with Path(path).open(encoding='ascii') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(lines)
                header, lines = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError(f'{Path(path).name}: sequence data appears before any FASTA header.')
                lines.append(line.strip())
    if header is not None:
        yield header, ''.join(lines)


def _split_fasta(stage, upstream_relative, directory, rename, cancelled):
    """Split one multi-record FASTA into a group-per-file panel with unique identifiers.

    SequenceReader refuses duplicate identifiers inside one file, and the upstream
    panels reuse their first token across records, so the identifiers must be
    rewritten. Every derived file is re-read and compared against the records it
    came from, so a rewrite that changed sequence bytes aborts provisioning.
    """
    grouped = {}
    for header, sequence in _fasta_records(stage / upstream_relative):
        group, identifier = rename(header)
        for part in (group, identifier):
            if not _SAFE_NAME.fullmatch(part):
                raise ValueError(f'Derived reference name is not a safe identifier: {part!r}')
        records = grouped.setdefault(group, [])
        if any(name == identifier for name, _ in records):
            raise ValueError(f'Derived reference identifier {identifier!r} is not unique within {group}.')
        records.append((identifier, sequence))
    written = []
    for group, records in grouped.items():
        check_cancelled(cancelled)
        relative = f'{directory}/{group}.fasta'
        target = stage / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('w', encoding='ascii') as out:
            for identifier, sequence in records:
                out.write(f'>{identifier}\n{sequence}\n')
        if list(_fasta_records(target)) != records:
            raise ValueError(f'Splitting {upstream_relative} changed sequence bytes for {group}.')
        written.append((group, relative, records))
    return written


def _sccmec_split_name(header):
    name, _, rest = header.partition(' ')
    accession, _, sccmec_type = rest.partition('|')
    if not name or not accession or not sccmec_type:
        raise ValueError(f'Unexpected SCCmec reference header: {header!r}')
    return name, f'{name}__{accession}__{sccmec_type}'


def _capsule_split_name(header):
    fields = header.split('__')
    if len(fields) != 4 or not all(fields):
        raise ValueError(f'Unexpected capsule marker header: {header!r}')
    return fields[1], f'{fields[1]}_{fields[2]}'


def _load_yaml(text):
    try:
        import yaml
    except ImportError as error:  # pragma: no cover - exercised only without the dev dependency
        raise ValueError('Deriving the SCCmec typing rules requires PyYAML, a development-only '
                         'dependency used at staging time; it never enters the frozen bundle.') from error
    return yaml.safe_load(text)


def sccmec_rules(document, table_text, targets):
    """Derive the IWG ccr/mec typing rules and cross-check them against the upstream table.

    The rules live in the manifest, never in Python, so the definitions actually
    in force are provable from the stored reference digest. Nothing is inferred:
    a name the upstream target panel does not define, or a type whose required
    targets do not match the upstream per-type table exactly, aborts staging.
    """
    targets = list(targets)
    if not isinstance(document, dict):
        raise ValueError('The SCCmec rule document is not a mapping.')
    declared = document.get('targets')
    if not isinstance(declared, list) or sorted(declared) != sorted(targets):
        raise ValueError('The SCCmec rule document does not declare exactly the staged target set.')
    aliases, names = [], set(targets)
    for alias in document.get('aliases') or []:
        name, members = alias.get('name'), alias.get('targets')
        if not name or name in names or not isinstance(members, list) or not members:
            raise ValueError(f'SCCmec alias {name!r} is missing, duplicated or empty.')
        unknown = [member for member in members if member not in targets]
        if unknown:
            raise ValueError(f'SCCmec alias {name!r} names targets absent from the panel: {", ".join(unknown)}.')
        names.add(name)
        aliases.append({'name': name, 'targets': list(members)})
    alias_targets = {alias['name']: alias['targets'] for alias in aliases}
    types, seen = [], set()
    for definition in document.get('types') or []:
        name = definition.get('name')
        required, excluded = definition.get('targets'), definition.get('excludes') or []
        if not name or name in seen or not isinstance(required, list) or not required:
            raise ValueError(f'SCCmec type {name!r} is missing, duplicated or empty.')
        for member in (*required, *excluded):
            if member not in names:
                raise ValueError(f'SCCmec type {name!r} names {member!r}, which is neither a target nor an alias.')
        seen.add(name)
        types.append({'name': name, 'targets': list(required), 'excludes': list(excluded)})
    if not types:
        raise ValueError('The SCCmec rule document declares no type definitions.')
    observed = {}
    for row in csv.DictReader(io.StringIO(table_text), delimiter='\t'):
        observed.setdefault(row['type'], set()).add(row['target'])
    for definition in types:
        expected = set()
        for member in definition['targets']:
            expected |= set(alias_targets.get(member, [member]))
        actual = observed.get(definition['name'])
        if actual is None:
            raise ValueError(f'SCCmec type {definition["name"]} has no rows in the upstream target table.')
        if expected != actual:
            raise ValueError(f'SCCmec type {definition["name"]} does not match the upstream target table: '
                             f'{sorted(expected ^ actual)}.')
        for member in definition['excludes']:
            present = set(alias_targets.get(member, [member])) & actual
            if present:
                raise ValueError(f'SCCmec type {definition["name"]} excludes {member!r}, which the upstream '
                                 f'table records as present: {sorted(present)}.')
    engine = ((document.get('engine') or {}).get('params') or {})
    if (engine.get('min_pident'), engine.get('min_coverage')) != (SCCMEC_ENGINE_DEFAULTS['target_min_pident'],
                                                                 SCCMEC_ENGINE_DEFAULTS['target_min_coverage']):
        raise ValueError('Upstream SCCmec engine parameters differ from the reviewed calling floors.')
    return {'schema_version': str((document.get('metadata') or {}).get('version') or ''),
            'targets': targets, 'aliases': aliases, 'types': types}


def _profile_fields(path, locus, genes):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        header = [value.strip() for value in next(csv.reader(handle, delimiter='\t'), [])]
    if len(set(header)) != len(header) or 'ST' not in header:
        raise ValueError(f'{locus}: the profile table has no single ST column.')
    missing = set(genes) - set(header)
    if missing:
        raise ValueError(f'{locus}: the profile table lacks staged loci: {", ".join(sorted(missing))}.')
    extra = [column for column in header if column != 'ST' and column not in genes]
    if len(extra) != 1:
        raise ValueError(f'{locus}: expected exactly one lineage column beside ST and the loci, found {extra}.')
    return 'ST', extra[0]


def provision_characterization_references(root, *, source_roots=None, sources=None, cancelled=None, progress=None):
    """Explicitly install a versioned panel; source_roots may reuse local checkouts.

    Existing snapshots are never overwritten. Local source bytes are identified
    by their hashes; a local checkout is not falsely certified as unmodified.
    Derived panels (FASTA splits, the SCCmec rule table) are produced here, at
    install time, so the runtime only ever reads verified FASTA and JSON.
    """
    check_cancelled(cancelled)
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    selected = tuple(sources if sources is not None else SOURCES)
    unknown = [key for key in selected if key not in SOURCES]
    if unknown:
        raise ValueError(f'Unknown reference source: {", ".join(unknown)}')
    local = {key: Path(value) for key, value in (source_roots or {}).items()}
    manifest = {'format_version': 2, 'source_revision': KLEBORATE_REVISION,
                'sources': {}, 'species': [], 'virulence': {}, 'locus_profiles': {},
                'sccmec': {}, 'capsule': {}, 'files': [], 'created_utc': datetime.now(UTC).isoformat(),
                'source_repository': SOURCES['kleborate'].repository,
                'taxonomy_source': TAXONOMY_SOURCE,
                'source_mode': 'local_checkout' if local else 'pinned_https',
                'license_notice': 'Each staged panel retains its own upstream licence beside its data; see '
                                  'the source-LICENSE-<source> files and the sources section. Sequence records '
                                  'retain their original accessions. These panels are not complete databases.',
                'limitations': ['Reference-supported classification, not independent clinical validation.',
                                'Escherichia coli/Shigella cannot be distinguished by ANI alone.',
                                'One representative per listed taxon is not exhaustive within-taxon diversity.']}
    tasks = []
    for key in selected:
        source = SOURCES[key]
        manifest['sources'][key] = {'repository': source.repository, 'revision': source.revision,
                                    'license': source.license_id, 'license_file': f'source-LICENSE-{key}',
                                    'mode': 'local_checkout' if key in local else 'pinned_https'}
        tasks.append((key, source.license_path, f'source-LICENSE-{key}'))
    if 'kleborate' in selected:
        for accession, genus, species, subspecies in GENOMES:
            path = f'species/{accession}.fna.gz'
            manifest['species'].append({'id': accession, 'path': path, 'genus': genus, 'species': species,
                                        'subspecies': subspecies, 'outgroup': species not in KPSC_SPECIES or genus != 'Klebsiella',
                                        'taxonomy_basis': 'NCBI assembly GCF_000005845.2 (E. coli K-12 MG1655)' if genus == 'Escherichia' else TAXONOMY_SOURCE})
            tasks.append(('kleborate', f'test/test_genomes/{accession}.fna.gz', path))
        for locus, (module, genes) in VIRULENCE_LOCI.items():
            manifest['virulence'][locus] = []
            for gene in genes:
                path = f'virulence/{locus}/{gene}.fasta'
                manifest['virulence'][locus].append({'gene': gene, 'path': path})
                tasks.append(('kleborate', f'kleborate/modules/{module}/data/{gene}.fasta', path))
        for locus in PROFILE_LOCI:
            module = VIRULENCE_LOCI[locus][0]
            tasks.append(('kleborate', f'kleborate/modules/{module}/data/profiles.tsv', f'virulence/{locus}/profiles.tsv'))
    if 'sccmec' in selected:
        for name in ('sccmec-targets.fasta', 'sccmec-targets.yaml', 'sccmec-targets.tsv',
                     'sccmec-regions.fasta', 'sccmec-regions.tsv'):
            tasks.append(('sccmec', f'data/{name}', f'sources/sccmec/{name}'))
    if 'kaptive' in selected:
        tasks.append(('kaptive', 'reference_database/wzi_wzc_db.fasta', 'sources/kaptive/wzi_wzc_db.fasta'))
    with tempfile.TemporaryDirectory(prefix='.characterization-download-', dir=root) as temporary:
        stage = Path(temporary) / 'snapshot'
        stage.mkdir()
        for index, (source_key, upstream, relative) in enumerate(tasks, 1):
            entry = _fetch(source_key, upstream, stage / relative, cancelled, source_roots=local)
            entry['path'] = relative
            manifest['files'].append(entry)
            if sum(item['bytes'] for item in manifest['files']) > 512 * 1024 * 1024:
                raise ValueError('Characterization reference snapshot exceeds 512 MiB.')
            if progress:
                progress(index, len(tasks), f'Verified characterization reference {relative}')
        if 'kleborate' in selected:
            for locus in PROFILE_LOCI:
                relative = f'virulence/{locus}/profiles.tsv'
                st_field, lineage_field = _profile_fields(stage / relative, locus, VIRULENCE_LOCI[locus][1])
                manifest['locus_profiles'][locus] = {'module': VIRULENCE_LOCI[locus][0], 'path': relative,
                                                     'st_field': st_field, 'lineage_field': lineage_field,
                                                     'source': 'kleborate'}
        if 'sccmec' in selected:
            manifest['sccmec'] = _stage_sccmec(stage, manifest, cancelled)
        if 'kaptive' in selected:
            manifest['capsule'] = _stage_capsule(stage, manifest, cancelled)
        manifest['files'].sort(key=lambda item: item['path'])
        manifest['reference_digest'] = reference_digest(manifest)
        (stage / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        validate_characterization_references(stage, cancelled=cancelled)
        destination = root / ('characterization-' + manifest['reference_digest'][:20])
        check_cancelled(cancelled)
        if destination.exists():
            existing = validate_characterization_references(destination, cancelled=cancelled)
            if existing['reference_digest'] != manifest['reference_digest']:
                raise ValueError('Existing reference destination has a different fingerprint.')
        else:
            os.replace(stage, destination)
    return {'path': str(destination), 'reference_digest': manifest['reference_digest'],
            'species_count': len(manifest['species']), 'virulence_loci': list(manifest['virulence']),
            'locus_profiles': list(manifest['locus_profiles']),
            'sccmec_targets': len(manifest['sccmec'].get('targets', [])),
            'sccmec_regions': len(manifest['sccmec'].get('regions', [])),
            'capsule_loci': [entry['gene'] for entry in manifest['capsule'].get('loci', [])]}


def _stage_sccmec(stage, manifest, cancelled):
    targets, regions = [], []
    for group, relative, records in _split_fasta(stage, 'sources/sccmec/sccmec-targets.fasta', 'sccmec/targets',
                                                 _sccmec_split_name, cancelled):
        targets.append({'gene': group, 'path': relative, 'reference_count': len(records)})
        manifest['files'].append(_derived_entry(stage, relative, cancelled, source='sccmec',
                                                derived_from='sources/sccmec/sccmec-targets.fasta'))
    for group, relative, records in _split_fasta(stage, 'sources/sccmec/sccmec-regions.fasta', 'sccmec/regions',
                                                 _sccmec_split_name, cancelled):
        if len(records) != 1:
            raise ValueError(f'SCCmec region {group} must have exactly one cassette reference.')
        identifier, sequence = records[0]
        regions.append({'gene': group, 'path': relative, 'accession': identifier.split('__')[1],
                        'reference_bp': len(sequence)})
        manifest['files'].append(_derived_entry(stage, relative, cancelled, source='sccmec',
                                                derived_from='sources/sccmec/sccmec-regions.fasta'))
    rules = sccmec_rules(_load_yaml((stage / 'sources/sccmec/sccmec-targets.yaml').read_text(encoding='utf-8')),
                         (stage / 'sources/sccmec/sccmec-targets.tsv').read_text(encoding='utf-8-sig'),
                         [entry['gene'] for entry in targets])
    return {'source': 'sccmec', 'schema_version': rules.pop('schema_version'),
            'targets': sorted(targets, key=lambda entry: entry['gene']),
            'regions': sorted(regions, key=lambda entry: entry['gene']), 'rules': rules,
            'engine_defaults': dict(SCCMEC_ENGINE_DEFAULTS), 'not_assayed': list(SCCMEC_NOT_ASSAYED),
            'upstream_paths': ['sources/sccmec/sccmec-targets.fasta', 'sources/sccmec/sccmec-targets.yaml',
                               'sources/sccmec/sccmec-targets.tsv', 'sources/sccmec/sccmec-regions.fasta',
                               'sources/sccmec/sccmec-regions.tsv'],
            'derivation': 'Multi-FASTA split per target/region with identifiers rewritten to '
                          '<name>__<accession>__<type>; sequence bytes asserted unchanged. Rules parsed from '
                          'sccmec-targets.yaml and cross-checked target-for-target against sccmec-targets.tsv.'}


def _stage_capsule(stage, manifest, cancelled):
    loci = []
    for group, relative, records in _split_fasta(stage, 'sources/kaptive/wzi_wzc_db.fasta', 'capsule',
                                                 _capsule_split_name, cancelled):
        loci.append({'gene': group, 'path': relative, 'allele_count': len(records)})
        manifest['files'].append(_derived_entry(stage, relative, cancelled, source='kaptive',
                                                derived_from='sources/kaptive/wzi_wzc_db.fasta'))
    return {'source': 'kaptive', 'loci': sorted(loci, key=lambda entry: entry['gene']),
            'upstream_paths': ['sources/kaptive/wzi_wzc_db.fasta'],
            'k_locus_typing': 'not_staged',
            'derivation': 'SRST2-style headers <cluster>__<gene>__<allele>__<seq> rewritten to <gene>_<allele>; '
                          'sequence bytes asserted unchanged. Only the wzi/wzc marker database is staged: the '
                          'Kaptive K/O locus references and their match-confidence logic are not.'}


def validate_characterization_references(path, *, cancelled=None):
    path = Path(path).resolve()
    manifest_path = path / 'manifest.json'
    if manifest_path.stat().st_size > 2 * 1024 * 1024:
        raise ValueError('Characterization manifest exceeds the 2 MiB limit.')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('format_version') not in _DIGEST_KEYS:
        raise ValueError('Characterization reference manifest format version is not supported.')
    if manifest['format_version'] >= 2 and any(section not in manifest for section in _V2_SECTIONS):
        raise ValueError('Characterization reference manifest omits a declared organism-module section.')
    try:
        computed = reference_digest(manifest)
    except KeyError as error:
        raise ValueError(f'Characterization reference manifest omits a fingerprinted section: {error}') from error
    if computed != manifest.get('reference_digest'):
        raise ValueError('Characterization reference manifest fingerprint is invalid.')
    seen = set()
    for entry in manifest['files']:
        check_cancelled(cancelled)
        relative = entry['path']
        target = path / relative
        if relative in seen or target.is_symlink() or not target.resolve().is_relative_to(path):
            raise ValueError('Unsafe or duplicated characterization reference path.')
        seen.add(relative)
        if not target.is_file() or target.stat().st_size != entry['bytes'] or file_sha256(target, cancelled) != entry['sha256']:
            raise ValueError(f'Characterization reference was changed or is missing: {relative}')
    referenced = [entry['path'] for entry in manifest['species']]
    referenced += [entry['path'] for locus in manifest['virulence'].values() for entry in locus]
    if manifest['format_version'] >= 2:
        referenced += [entry['path'] for entry in manifest['locus_profiles'].values()]
        referenced += [entry['path'] for key in ('targets', 'regions') for entry in manifest['sccmec'].get(key, [])]
        referenced += [entry['path'] for entry in manifest['capsule'].get('loci', [])]
    if not set(referenced).issubset(seen):
        raise ValueError('Some characterization references have no source hash.')
    if len({entry['id'] for entry in manifest['species']}) != len(manifest['species']):
        raise ValueError('Species reference identifiers are not unique.')
    return manifest
