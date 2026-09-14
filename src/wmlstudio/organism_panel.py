"""Explicit, immutable provisioning of a broad genus/species ANI triage panel.

The bundled characterization starter is deliberately Klebsiella-scoped, so an
isolate from any other genus cannot be separated from "no match" by it. This
module installs one pinned NCBI RefSeq reference per common clinical taxon so
that whole-genome comparison has something to compare against.

One genome per taxon is a triage panel. It is not a species database, it does
not represent within-species diversity, and a folder created from its labels is
a filing decision rather than a laboratory identification. Nothing here is
downloaded by analysis or at application startup; a person asks for it.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import tempfile
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from .characterization_refs import reference_digest, validate_characterization_references
from .sequence import check_cancelled, file_sha256

PANEL_FORMAT_VERSION = 1
PANEL_DIRECTORY = 'species-panel'
PANEL_PREFIX = 'species-panel-'
# The snapshot identity: bump this when the pinned accession set changes.
PANEL_REVISION = 'ncbi-refseq-species-panel-2026-09-14.2'
NCBI_BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
MAX_FILE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024
PANEL_LIMITATIONS = (
    'A small reference set is a triage panel, not a representation of within-species diversity.',
    'Most taxa carry one reference. L. monocytogenes carries two, because its lineages straddle the 95% ANI species line and one reference leaves common genomes unresolved. Other species with comparable internal diversity may be under-covered in the same way.',
    'This panel does not distinguish Escherichia coli from Shigella.',
    'An organism folder created from these labels records where a copy is stored; it is not a laboratory identification.',
    'Organism names are the labels NCBI records for these assemblies; they were not independently verified here.',
)
LICENSE_NOTICE = (
    'NCBI RefSeq genome records. NCBI places no restrictions on the use or distribution of '
    'these data, although individual submitters may retain rights in their submissions. '
    'This panel is a small curated set, mostly one assembly per taxon, and is not an exhaustive species database.'
)
# (accession, assembly_dir, genus, species, strain, gz_bytes, gz_sha256, fasta_sha256)
# Every row was fetched from the NCBI FTP path below, its published MD5 recomputed
# and matched, and both hashes computed from the downloaded bytes.
SPECIES_PANEL = (
    ('GCF_000016305.1', 'GCF_000016305.1_ASM1630v1', 'Klebsiella', 'pneumoniae', 'MGH 78578', 1679405,
     'b7592666d5ee63aa83517c94885c25a57079496754a7504cb591086a7ac0e8ae',
     'a62bedcb43eb0cebdb3f1fa6c29e23cc220dd01823d6f4abfc858b79c080889c'),
    ('GCF_000005845.2', 'GCF_000005845.2_ASM584v2', 'Escherichia', 'coli', 'K-12 MG1655', 1379902,
     'a96d3cfa58c88d477013c768f90a11d2386d42a70d566abfb6d9013dcbd24255',
     '53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4'),
    ('GCF_000009645.1', 'GCF_000009645.1_ASM964v1', 'Staphylococcus', 'aureus', 'N315', 825612,
     'aced925b806adf46fb45fea80e29981dcdd001f6030a4da63a67d44d2dfe8509',
     'd7b93eeed77ce08d041f597d28fe3d7310fead42f9f790b6325539fb70c458f9'),
    ('GCF_000011925.1', 'GCF_000011925.1_ASM1192v1', 'Staphylococcus', 'epidermidis', 'RP62A', 767953,
     '0551acafbe377ce80f0b5ff69db54ea96566575acbb12a4abe15cde5184e7a0f',
     'bbd0cd4d1c457d3c70b1c6f58884ef6e686eedbf5314200674e8d6ce29db64bc'),
    ('GCF_000250945.2', 'GCF_000250945.2_ASM25094v2', 'Enterococcus', 'faecium', 'Aus0004', 893620,
     '6176d2bd536dd60279150c1d0238e48125b2f24dc2ca43c047e50c87a26408a5',
     '25a669e12bb5c1a530e62b55f51e5772273a2c6c05cfa78554988a9b5bcf4492'),
    ('GCF_000172575.2', 'GCF_000172575.2_ASM17257v2', 'Enterococcus', 'faecalis', 'OG1RF', 806818,
     '10e17448e466fd9ccaaedeb8b961289d397a00dcfdf5f970231e39bf462412c5',
     'a3f1dd01e99d75596d837266edd95aedaecc0ba0fd967466d2e0952cbb1537ef'),
    ('GCF_000018445.1', 'GCF_000018445.1_ASM1844v1', 'Acinetobacter', 'baumannii', 'ACICU', 1184637,
     '86e0492fb57e88d18b4cfa225324ab173e836f50c6c8e719f769b8af6d66e539',
     'bb363169f345de1cc78f012b058ba8fbb33334e442d791bb7ea1d5b732226ba1'),
    ('GCF_000014625.1', 'GCF_000014625.1_ASM1462v1', 'Pseudomonas', 'aeruginosa', 'UCBPP-PA14', 1868144,
     '1d17779cd3bf5cfb32240d2c9ee5a1842fbda13ea03e385b37adf197e79ab621',
     '035464d72feb51de37217d5989695768f4a7571207bdd023ea7732edffa353e4'),
    ('GCF_000195995.1', 'GCF_000195995.1_ASM19599v1', 'Salmonella', 'enterica', 'Typhi CT18', 1526348,
     'd5157297061c228f9675afcd91f5d5591248f8463dba4a0bc172a9fd754ac9ef',
     'cb6ea336f0b444f2662bd86a43ec15c6e1287c7aba2991dafa8c0589cbdc36f8'),
    ('GCF_000007045.1', 'GCF_000007045.1_ASM704v1', 'Streptococcus', 'pneumoniae', 'R6', 602297,
     'd6b5c69dd3cc1b7fe65bb6293b89ab8da7caa58ddacbe06458a78e87786de089',
     'aeb72b21509185ba45cb1b77a91ef117520965e93806be86665c19f4bccb7be2'),
    ('GCF_000012705.1', 'GCF_000012705.1_ASM1270v1', 'Streptococcus', 'agalactiae', 'A909', 624945,
     '435fc1701254278f2af7a70a7f0cf97049224f19b002386f6b2a7f1727081428',
     '1310eea2c3369aaf8751329f8ee1c57d3f42534a162472ffbdec444f284921b7'),
    ('GCF_000011765.3', 'GCF_000011765.3_ASM1176v2', 'Streptococcus', 'pyogenes', 'MGAS5005', 544025,
     '75feedcf567040f8c770c14aeb27ea3969524a90316d9b1a016e67de6a806198',
     '142f5d4d9a8eb33b3f5f97a5c5f799ed6d66d2aa31d20f50d33eb0141d1cf224'),
    ('GCF_000008285.1', 'GCF_000008285.1_ASM828v1', 'Listeria', 'monocytogenes', '4b F2365', 857567,
     'd2684c7aa24f9373c2ee02510b6c9ef1e0d0b6c1ccaaa6a88063b96d76d244b5',
     '999b109bccfb358ed50a1382ee3839e4449cee04849480456fd7b31d8a98a582'),
    # A second L. monocytogenes reference, deliberately. Lineages I and II sit
    # either side of the 95% ANI species line: a lineage II genome measured only
    # 94.83% against the 4b F2365 reference above and was correctly, but
    # unhelpfully, left unresolved. 10403S is lineage II serotype 1/2a and is not
    # the EGD-e genome in the practice cohort, so the match is not a self-match.
    ('GCF_000168695.2', 'GCF_000168695.2_ASM16869v2', 'Listeria', 'monocytogenes', '1/2a 10403S', 856931,
     '2be4f2809f0f487da4b7673c80c8ada7ac9b99ca4d7e09170d855ed34b3a3251',
     '5d030d46e3104a3e014809ca352e5d9d64724b4566fba87a33239e2e7d23a30f'),
    ('GCF_000009105.1', 'GCF_000009105.1_ASM910v1', 'Neisseria', 'meningitidis', 'Z2491', 629571,
     '7c36f924111f1122d592815fb271fa3e159623bdc0caaa92c4eeed09520b830e',
     '2f138b39ace68b954ec807ae32785224c1d9e6b36013abc9649841f7e8a418dd'),
    ('GCF_000012185.1', 'GCF_000012185.1_ASM1218v1', 'Haemophilus', 'influenzae', '86-028NP', 563945,
     '06cb08fe9e280a55a1550072260695f09e6deeaed7111082721bf4e17deb273f',
     'a1988bf4c0a054e374f17a79575c651e5f8ea6e3bde07b29eaba92c7f9454101'),
    ('GCF_015732555.1', 'GCF_015732555.1_ASM1573255v1', 'Clostridioides', 'difficile', 'R20291', 1185904,
     '0fc5a377bddce66ebd46f88d873837053f8e1310cd03c2df164f2a815fda3cbe',
     'e819a841658059de0f74cf597e40b659038852915e15d392d352fea5c15ccb03'),
    ('GCF_000011865.1', 'GCF_000011865.1_ASM1186v1', 'Campylobacter', 'jejuni', 'RM1221', 504168,
     'af5a19817a3cec9362faabfd94a73fce514ad2fd43b9bed9a9c29765d885f29c',
     '865f502272ed7315a125298d411219c6db40e922b62f64276df6d8d0470328c4'),
)


def panel_digest(manifest):
    """The panel fingerprint, computed exactly as the characterization starter's is."""
    return reference_digest(manifest)


def assembly_source(assembly_dir: str, filename: str = '') -> str:
    """The NCBI FTP path for one pinned assembly directory, relative to NCBI_BASE."""
    accession = assembly_dir.split('_')[1].split('.')[0]
    if len(accession) != 9 or not accession.isdigit():
        raise ValueError(f'Unexpected assembly directory name: {assembly_dir}')
    prefix = assembly_dir.split('_')[0]
    parts = (prefix, accession[0:3], accession[3:6], accession[6:9], assembly_dir)
    return '/'.join(parts) + ('/' + filename if filename else '/')


def _open(source: str, cancelled=None, *, source_root=None):
    check_cancelled(cancelled)
    if source_root is not None:
        original = Path(source_root) / source
        if original.is_symlink() or not original.is_file():
            raise ValueError(f'Panel reference source must be a regular file: {original}')
        return original.open('rb')
    request = urllib.request.Request(NCBI_BASE + source,
                                     headers={'User-Agent': 'WMLSTudio-reference-provisioning'})
    return urllib.request.urlopen(request, timeout=30)


def published_md5(assembly_dir: str, cancelled=None, *, source_root=None) -> dict:
    """Read the checksums NCBI publishes beside an assembly, keyed by file name."""
    with _open(assembly_source(assembly_dir, 'md5checksums.txt'), cancelled, source_root=source_root) as handle:
        text = handle.read(1024 * 1024).decode('utf-8', 'strict')
    checksums = {}
    for line in text.splitlines():
        digest, _, name = line.strip().partition(' ')
        name = name.strip().lstrip('./')
        if digest and name:
            checksums[name] = digest
    if not checksums:
        raise ValueError(f'No published checksums were found for {assembly_dir}.')
    return checksums


def _fetch(row, target: Path, cancelled=None, *, source_root=None) -> dict:
    """Download one pinned genome, refusing any bytes the pins do not describe."""
    accession, assembly_dir, _genus, _species, _strain, size, gz_sha256, fasta_sha256 = row
    filename = f'{assembly_dir}_genomic.fna.gz'
    expected_md5 = published_md5(assembly_dir, cancelled, source_root=source_root).get(filename)
    if not expected_md5:
        raise ValueError(f'{accession}: NCBI publishes no checksum for {filename}.')
    source = assembly_source(assembly_dir, filename)
    target.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    published = hashlib.md5(usedforsecurity=False)
    with _open(source, cancelled, source_root=source_root) as handle, target.open('xb') as output:
        while chunk := handle.read(1024 * 1024):
            check_cancelled(cancelled)
            total += len(chunk)
            if total > MAX_FILE_BYTES:
                raise ValueError(f'{accession}: reference exceeds the 16 MiB per-file bound.')
            published.update(chunk)
            output.write(chunk)
    if total != size:
        raise ValueError(f'{accession}: expected {size} bytes from NCBI, received {total}.')
    digest = file_sha256(target, cancelled)
    if digest != gz_sha256:
        raise ValueError(f'{accession}: downloaded bytes do not match the pinned SHA-256.')
    if published.hexdigest() != expected_md5:
        raise ValueError(f'{accession}: downloaded bytes do not match the checksum NCBI publishes.')
    uncompressed = hashlib.sha256()
    length = 0
    with gzip.open(target, 'rb') as stream:
        while chunk := stream.read(1024 * 1024):
            check_cancelled(cancelled)
            uncompressed.update(chunk)
            length += len(chunk)
    if uncompressed.hexdigest() != fasta_sha256:
        raise ValueError(f'{accession}: the decompressed sequence does not match the pinned SHA-256.')
    return {'path': target.name, 'bytes': total, 'sha256': digest, 'md5': expected_md5,
            'uncompressed_bytes': length, 'uncompressed_sha256': fasta_sha256,
            'source_url': NCBI_BASE + source}


def panel_manifest(source_mode: str) -> dict:
    """The manifest skeleton, shaped so the existing panel validator accepts it."""
    manifest = {'format_version': PANEL_FORMAT_VERSION, 'source_revision': PANEL_REVISION,
                'species': [], 'virulence': {}, 'files': [],
                'created_utc': datetime.now(UTC).isoformat(),
                'panel_kind': 'species_panel', 'source_repository': NCBI_BASE,
                'source_mode': source_mode, 'license_notice': LICENSE_NOTICE,
                'limitations': list(PANEL_LIMITATIONS)}
    for accession, assembly_dir, genus, species, strain, *_ in SPECIES_PANEL:
        manifest['species'].append({
            'id': accession, 'path': f'species/{accession}.fna.gz', 'genus': genus,
            'species': species, 'subspecies': '',
            # Cross-genus triage: every reference separates this query from the
            # others, so each one is an outgroup for every other taxon here.
            'outgroup': True,
            'taxonomy_basis': f'NCBI RefSeq assembly {accession} ({genus} {species} {strain})',
            'assembly_directory': assembly_dir})
    return manifest


def provision_species_panel(root, *, source_root=None, cancelled=None, progress=None) -> dict:
    """Install the pinned panel under root, publishing nothing that failed a check.

    Mirrors the characterization starter's discipline: exclusive creation, bounded
    downloads, staging inside root, full validation before publication, and a
    refusal to replace an existing snapshot that has a different fingerprint.
    """
    check_cancelled(cancelled)
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = panel_manifest('local_mirror' if source_root else 'pinned_https')
    with tempfile.TemporaryDirectory(prefix='.species-panel-download-', dir=root) as temporary:
        stage = Path(temporary) / 'snapshot'
        stage.mkdir()
        for index, row in enumerate(SPECIES_PANEL, 1):
            relative = f'species/{row[0]}.fna.gz'
            entry = _fetch(row, stage / relative, cancelled, source_root=source_root)
            entry['path'] = relative
            manifest['files'].append(entry)
            if sum(item['bytes'] for item in manifest['files']) > MAX_TOTAL_BYTES:
                raise ValueError('Species panel snapshot exceeds the 64 MiB download bound.')
            if progress:
                progress(index, len(SPECIES_PANEL), f'Verified species reference {row[0]}')
        manifest['files'].sort(key=lambda item: item['path'])
        manifest['reference_digest'] = panel_digest(manifest)
        (stage / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        validate_characterization_references(stage, cancelled=cancelled)
        destination = root / (PANEL_PREFIX + manifest['reference_digest'][:20])
        check_cancelled(cancelled)
        if destination.exists():
            existing = validate_characterization_references(destination, cancelled=cancelled)
            if existing['reference_digest'] != manifest['reference_digest']:
                raise ValueError('Existing species panel destination has a different fingerprint.')
        else:
            os.replace(stage, destination)
    return {'path': str(destination), 'reference_digest': manifest['reference_digest'],
            'species_count': len(manifest['species']),
            'stored_bytes': sum(item['bytes'] for item in manifest['files'])}


def bundled_species_panel() -> Path | None:
    """No broad panel ships inside the application; it is always installed on request."""
    return None


def installed_species_panel(data_root) -> Path | None:
    """Locate an installed panel without downloading or re-hashing on a UI refresh."""
    base = Path(data_root).expanduser() / PANEL_DIRECTORY
    if not base.is_dir():
        return None
    snapshots = sorted(path for path in base.iterdir()
                       if path.is_dir() and path.name.startswith(PANEL_PREFIX)
                       and (path / 'manifest.json').is_file())
    return snapshots[-1] if snapshots else None


def verify_pins(cancelled=None, *, source_root=None) -> dict:
    """Compare the pinned table against NCBI's published checksums without downloading.

    A curator can run this cheaply; it detects an upstream re-release before a user
    meets it as a failed install. It confirms the files NCBI still publishes, not
    that our pinned sequence hashes are correct.
    """
    report = {'checked': 0, 'available': [], 'missing': [], 'md5': {}}
    for row in SPECIES_PANEL:
        check_cancelled(cancelled)
        accession, assembly_dir = row[0], row[1]
        filename = f'{assembly_dir}_genomic.fna.gz'
        try:
            checksums = published_md5(assembly_dir, cancelled, source_root=source_root)
        except (OSError, ValueError) as error:
            report['missing'].append({'accession': accession, 'reason': str(error)})
            continue
        report['checked'] += 1
        if filename in checksums:
            report['available'].append(accession)
            report['md5'][accession] = checksums[filename]
        else:
            report['missing'].append({'accession': accession,
                                      'reason': f'{filename} is no longer published.'})
    return report
