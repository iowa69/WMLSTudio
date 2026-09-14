"""Explicit, immutable local reference provisioning for genomic characterization.

The source snapshot is the official Kleborate repository; sequence accessions
and taxonomy labels are retained, not presented as an exhaustive taxonomy DB.
No download function is called by analysis or application startup.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from .sequence import check_cancelled, file_sha256

KLEBORATE_REVISION = '550ce22a2c01c76064f4dabf403704ee2293356e'
SOURCE_BASE = f'https://raw.githubusercontent.com/klebgenomics/Kleborate/{KLEBORATE_REVISION}/'
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


def bundled_reference_root():
    """Locate an installed starter without downloading or validating on UI refresh."""
    roots = []
    if getattr(sys, 'frozen', False):
        roots.append(Path(getattr(sys, '_MEIPASS', Path(sys.executable).parent))
                     / 'wmlstudio' / 'resources' / 'characterization' / 'starter')
    roots.append(Path(__file__).resolve().parent / 'resources' / 'characterization' / 'starter')
    return next((root for root in roots if (root / 'manifest.json').is_file()), None)


def reference_digest(manifest):
    content = {key: manifest[key] for key in ('format_version', 'source_revision', 'species', 'virulence', 'files')}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _fetch(source, target, cancelled=None, *, local_root=None):
    check_cancelled(cancelled)
    if local_root is not None:
        original = Path(local_root) / source
        if original.is_symlink() or not original.is_file():
            raise ValueError(f'Reference source must be a regular file: {original}')
        handle = original.open('rb')
    else:
        request = urllib.request.Request(SOURCE_BASE + source, headers={'User-Agent': 'WMLSTudio-reference-provisioning'})
        handle = urllib.request.urlopen(request, timeout=30)
    target.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with handle, target.open('xb') as output:
        while chunk := handle.read(1024 * 1024):
            check_cancelled(cancelled)
            total += len(chunk)
            if total > 64 * 1024 * 1024:
                raise ValueError(f'Reference file exceeds the 64 MiB download bound: {source}')
            output.write(chunk)
    return {'path': target.name, 'bytes': total, 'sha256': file_sha256(target, cancelled),
            'source_url': SOURCE_BASE + source}


def provision_characterization_references(root, *, source_root=None, cancelled=None, progress=None):
    """Explicitly install a versioned panel; source_root may reuse a local checkout.

    Existing snapshots are never overwritten. Local source bytes are identified
    by their hashes; a local checkout is not falsely certified as unmodified.
    """
    check_cancelled(cancelled)
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = {'format_version': 1, 'source_revision': KLEBORATE_REVISION,
                'species': [], 'virulence': {}, 'files': [], 'created_utc': datetime.now(UTC).isoformat(),
                'source_repository': 'https://github.com/klebgenomics/Kleborate',
                'taxonomy_source': TAXONOMY_SOURCE, 'source_mode': 'local_checkout' if source_root else 'pinned_https',
                'license_notice': 'Kleborate source/reference snapshot: GPL-3.0-or-later; sequence records retain their NCBI accessions. See source-LICENSE. This panel is not a complete species database.',
                'limitations': ['Reference-supported classification, not independent clinical validation.',
                                'Escherichia coli/Shigella cannot be distinguished by ANI alone.',
                                'One representative per listed taxon is not exhaustive within-taxon diversity.']}
    tasks = [('LICENSE', 'source-LICENSE')]
    for accession, genus, species, subspecies in GENOMES:
        path = f'species/{accession}.fna.gz'
        manifest['species'].append({'id': accession, 'path': path, 'genus': genus, 'species': species,
                                    'subspecies': subspecies, 'outgroup': species not in KPSC_SPECIES or genus != 'Klebsiella',
                                    'taxonomy_basis': 'NCBI assembly GCF_000005845.2 (E. coli K-12 MG1655)' if genus == 'Escherichia' else TAXONOMY_SOURCE})
        tasks.append((f'test/test_genomes/{accession}.fna.gz', path))
    for locus, (module, genes) in VIRULENCE_LOCI.items():
        manifest['virulence'][locus] = []
        for gene in genes:
            path = f'virulence/{locus}/{gene}.fasta'
            manifest['virulence'][locus].append({'gene': gene, 'path': path})
            tasks.append((f'kleborate/modules/{module}/data/{gene}.fasta', path))
    with tempfile.TemporaryDirectory(prefix='.characterization-download-', dir=root) as temporary:
        stage = Path(temporary) / 'snapshot'
        stage.mkdir()
        for index, (source, relative) in enumerate(tasks, 1):
            entry = _fetch(source, stage / relative, cancelled, local_root=source_root)
            entry['path'] = relative
            manifest['files'].append(entry)
            if sum(item['bytes'] for item in manifest['files']) > 512 * 1024 * 1024:
                raise ValueError('Characterization reference snapshot exceeds 512 MiB.')
            if progress:
                progress(index, len(tasks), f'Verified characterization reference {relative}')
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
            'species_count': len(manifest['species']), 'virulence_loci': list(manifest['virulence'])}


def validate_characterization_references(path, *, cancelled=None):
    path = Path(path).resolve()
    manifest_path = path / 'manifest.json'
    if manifest_path.stat().st_size > 2 * 1024 * 1024:
        raise ValueError('Characterization manifest exceeds the 2 MiB limit.')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('format_version') != 1 or reference_digest(manifest) != manifest.get('reference_digest'):
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
    if not set(referenced).issubset(seen):
        raise ValueError('Some characterization references have no source hash.')
    if len({entry['id'] for entry in manifest['species']}) != len(manifest['species']):
        raise ValueError('Species reference identifiers are not unique.')
    return manifest
