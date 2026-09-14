"""Pinned practice cohorts: the tables, the verified download and the honesty text."""

import gzip
import hashlib
import importlib.util
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from wmlstudio import practice_cohorts
from wmlstudio.characterization_refs import GENOMES as BUNDLED_PANEL_GENOMES
from wmlstudio.practice_cohorts import (
    COHORTS,
    cohort_digest,
    cohort_genomes,
    cohort_names,
    describe_cohorts,
    download_cohort,
    genome_filename,
    genome_relative_path,
    verify_cohort,
    verify_pins,
)
from wmlstudio.sequence import AnalysisCancelled

REPOSITORY = Path(__file__).resolve().parents[1]
CLI = REPOSITORY / 'studio_scripts' / 'fetch_practice_cohort.py'
_SPEC = importlib.util.spec_from_file_location('fetch_practice_cohort_test', CLI)
fetch_cli = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fetch_cli)


def mirror_cohort(tmp_path, monkeypatch, rows, name='synthetic-cohort', published_md5=None):
    """Build a local NCBI mirror and register a cohort pinned to exactly those bytes."""
    mirror = tmp_path / 'mirror'
    genomes = []
    for accession, assembly_dir, genus, species, strain, sequence in rows:
        entry = {'accession': accession, 'assembly_dir': assembly_dir}
        relative = genome_relative_path(entry)
        target = mirror / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = gzip.compress(f'>{accession}\n{sequence}\n'.encode(), mtime=0)
        target.write_bytes(payload)
        raw = gzip.decompress(payload)
        digest = hashlib.md5(payload).hexdigest()
        (target.parent / 'md5checksums.txt').write_text(
            f'{published_md5 or digest}  ./{target.name}\n', encoding='utf-8')
        genomes.append((accession, assembly_dir, genus, species, strain, len(payload), digest,
                        hashlib.sha256(payload).hexdigest(), hashlib.sha256(raw).hexdigest()))
    monkeypatch.setitem(COHORTS, name, {
        'title': 'Synthetic fixture cohort', 'purpose': 'Exercise the pinned download path.',
        'expect': 'Nothing; these bytes are not real assemblies.',
        'citation_context': ('Synthetic test fixture.',),
        'caveats': practice_cohorts.SHARED_CAVEATS, 'genomes': tuple(genomes)})
    return mirror, name


def three_rows():
    return [('GCF_900000001.1', 'GCF_900000001.1_FIXTUREv1', 'Testgenus', 'alpha', 'A1',
             'ACGTACGTACGTAAAA'),
            ('GCF_900000002.1', 'GCF_900000002.1_FIXTUREv1', 'Testgenus', 'beta', 'B2',
             'TTTTCCCCGGGGAAAA'),
            ('GCF_900000003.1', 'GCF_900000003.1_FIXTUREv1', 'Othergenus', 'gamma', 'C3',
             'GGGGTTTTAAAACCCC')]


def test_pinned_cohort_tables_are_wellformed_and_uniquely_pinned():
    assert cohort_names() == ('kpneumoniae-10', 'mixed-genus-20')
    assert len(COHORTS['kpneumoniae-10']['genomes']) == 10
    assert len(COHORTS['mixed-genus-20']['genomes']) == 20
    for name in cohort_names():
        entries = cohort_genomes(name)
        assert len({entry['accession'] for entry in entries}) == len(entries)
        assert len({genome_filename(entry) for entry in entries}) == len(entries)
        for entry in entries:
            assert re.fullmatch(r'GC[FA]_\d{9}\.\d+', entry['accession'])
            assert entry['assembly_dir'].startswith(entry['accession'] + '_')
            assert entry['gz_bytes'] > 0
            assert re.fullmatch(r'[0-9a-f]{32}', entry['gz_md5'])
            assert re.fullmatch(r'[0-9a-f]{64}', entry['gz_sha256'])
            assert re.fullmatch(r'[0-9a-f]{64}', entry['fasta_sha256'])
            assert entry['gz_sha256'] != entry['fasta_sha256']
            assert entry['genus'][:1].isupper() and entry['species'].islower()


def test_single_species_cohort_is_one_species_and_the_mixed_cohort_spans_many_genera():
    single = cohort_genomes('kpneumoniae-10')
    assert {(entry['genus'], entry['species']) for entry in single} == {('Klebsiella', 'pneumoniae')}
    mixed = cohort_genomes('mixed-genus-20')
    assert len({entry['genus'] for entry in mixed}) >= 15
    assert len({(entry['genus'], entry['species']) for entry in mixed}) >= 18
    required = {('Klebsiella', 'pneumoniae'), ('Staphylococcus', 'aureus'), ('Escherichia', 'coli'),
                ('Enterococcus', 'faecium'), ('Acinetobacter', 'baumannii'),
                ('Pseudomonas', 'aeruginosa'), ('Salmonella', 'enterica')}
    present = {(entry['genus'], entry['species']) for entry in mixed}
    assert required.issubset(present)
    assert any(genus == 'Streptococcus' for genus, _ in present)


def test_mixed_cohort_includes_a_species_complex_member_the_app_should_find_hard():
    mixed = cohort_genomes('mixed-genus-20')
    complex_members = [entry for entry in mixed if entry['genus'] == 'Klebsiella'
                       and entry['species'] in {'variicola', 'quasipneumoniae', 'quasivariicola',
                                                'africana'}]
    assert complex_members, 'the mixed cohort must contain a K. pneumoniae complex member'
    caveats = ' '.join(COHORTS['mixed-genus-20']['caveats'])
    assert 'species complex' in caveats
    assert 'hard call' in caveats


def test_no_practice_genome_reuses_an_assembly_from_the_bundled_reference_panel():
    bundled = {row[0] for row in BUNDLED_PANEL_GENOMES}
    for name in cohort_names():
        used = {entry['accession'] for entry in cohort_genomes(name)}
        assert not used & bundled, f'{name} would compare an isolate against its own reference'


def test_every_cohort_states_its_source_and_refuses_to_claim_validation():
    for summary in describe_cohorts():
        caveats = list(summary['caveats'])
        assert caveats and 'NCBI RefSeq' in caveats[0]
        assert any('not a validation of this software' in line for line in caveats)
        assert any('teaching cohort, not a validation set' in line for line in caveats)
        assert any('not a laboratory identification' in line for line in caveats)
        assert summary['download_bytes'] > 0 and summary['genomes'] >= 10


def test_download_writes_flat_verified_files_and_a_manifest_that_reverifies(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    destination = tmp_path / 'cohorts' / name
    seen = []
    result = download_cohort(name, destination, source_root=mirror,
                             progress=lambda index, total, message: seen.append(message))
    assert result['installed'] and result['downloaded'] == 3 and result['reused'] == 0
    assert len(seen) == 3
    written = sorted(path.name for path in destination.iterdir())
    assert written == ['Othergenus_gamma_C3_GCF_900000003.1.fna.gz',
                       'Testgenus_alpha_A1_GCF_900000001.1.fna.gz',
                       'Testgenus_beta_B2_GCF_900000002.1.fna.gz', 'manifest.json']
    assert not any(path.is_dir() for path in destination.iterdir())
    manifest = verify_cohort(destination)
    assert manifest['cohort'] == name and manifest['cohort_digest'] == cohort_digest(name)
    assert manifest['source'] == 'NCBI RefSeq' and manifest['source_mode'] == 'local_mirror'
    assert manifest['caveats'] == list(practice_cohorts.SHARED_CAVEATS)
    assert 'NCBI places no restrictions' in manifest['license_notice']
    for item in manifest['genomes']:
        assert item['source_url'].startswith(practice_cohorts.NCBI_BASE)
        assert item['checksum_basis'] == 'pinned_sha256_and_published_md5'


def test_a_genome_whose_bytes_do_not_match_the_pin_publishes_nothing(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    entry = cohort_genomes(name)[1]
    tampered = mirror / genome_relative_path(entry)
    tampered.write_bytes(gzip.compress(b'>tampered\nAAAA\n', mtime=0))
    destination = tmp_path / 'cohorts' / name
    with pytest.raises(ValueError, match='pinned SHA-256'):
        download_cohort(name, destination, source_root=mirror)
    assert not destination.exists()
    cache = destination.parent / '_cache'
    assert sorted(path.name for path in cache.iterdir()) == [
        cohort_genomes(name)[0]['gz_sha256'] + '.fna.gz']


def test_a_stale_published_checksum_is_refused_rather_than_accepted_quietly(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows(), published_md5='0' * 32)
    destination = tmp_path / 'cohorts' / name
    with pytest.raises(ValueError, match='publishes a different MD5'):
        download_cohort(name, destination, source_root=mirror)
    assert not destination.exists()


def test_an_accession_shared_between_cohorts_is_fetched_once(tmp_path, monkeypatch):
    rows = three_rows()
    mirror, first = mirror_cohort(tmp_path, monkeypatch, rows, name='shared-a')
    monkeypatch.setitem(COHORTS, 'shared-b', dict(COHORTS[first],
                                                  genomes=COHORTS[first]['genomes'][:2]))
    fetched = []
    original = practice_cohorts._fetch_genome

    def spy(entry, target, **options):
        fetched.append(entry['accession'])
        return original(entry, target, **options)

    monkeypatch.setattr(practice_cohorts, '_fetch_genome', spy)
    download_cohort(first, tmp_path / 'cohorts' / first, source_root=mirror)
    download_cohort('shared-b', tmp_path / 'cohorts' / 'shared-b', source_root=mirror)
    assert fetched == [row[0] for row in rows]
    assert len(fetched) == len(set(fetched))
    assert verify_cohort(tmp_path / 'cohorts' / 'shared-b')['cohort'] == 'shared-b'


def test_repeated_download_reports_the_installed_cohort_without_refetching(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    destination = tmp_path / 'cohorts' / name
    download_cohort(name, destination, source_root=mirror)
    monkeypatch.setattr(practice_cohorts, '_fetch_genome',
                        lambda *a, **k: pytest.fail('an installed cohort must not be refetched'))
    again = download_cohort(name, destination, source_root=mirror)
    assert again['installed'] is False and again['reused'] == 3


def test_a_different_cohort_never_overwrites_an_installed_one(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    monkeypatch.setitem(COHORTS, 'other-cohort', dict(COHORTS[name],
                                                      genomes=COHORTS[name]['genomes'][:2]))
    destination = tmp_path / 'cohorts' / name
    download_cohort(name, destination, source_root=mirror)
    with pytest.raises(ValueError, match='different practice cohort is already installed'):
        download_cohort('other-cohort', destination, source_root=mirror)
    assert verify_cohort(destination)['cohort'] == name
    monkeypatch.setitem(COHORTS, name, dict(COHORTS[name],
                                            genomes=COHORTS[name]['genomes'][:2]))
    with pytest.raises(ValueError, match='does not match the pinned cohort definition'):
        download_cohort(name, destination, source_root=mirror)


def test_verify_cohort_detects_a_mutated_byte_a_missing_file_and_a_short_manifest(tmp_path,
                                                                                 monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    destination = tmp_path / 'cohorts' / name
    download_cohort(name, destination, source_root=mirror)
    target = destination / genome_filename(cohort_genomes(name)[0])
    payload = bytearray(target.read_bytes())
    payload[-1] ^= 0xFF
    target.write_bytes(bytes(payload))
    with pytest.raises(ValueError, match='changed or is missing'):
        verify_cohort(destination)
    target.unlink()
    with pytest.raises(ValueError, match='changed or is missing'):
        verify_cohort(destination)
    manifest = json.loads((destination / 'manifest.json').read_text(encoding='utf-8'))
    manifest['genomes'] = manifest['genomes'][:1]
    (destination / 'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    with pytest.raises(ValueError, match='incomplete'):
        verify_cohort(destination)


def test_verify_cohort_refuses_a_manifest_that_disagrees_with_the_pinned_table(tmp_path,
                                                                              monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    destination = tmp_path / 'cohorts' / name
    download_cohort(name, destination, source_root=mirror)
    path = destination / 'manifest.json'
    manifest = json.loads(path.read_text(encoding='utf-8'))
    manifest['genomes'][0]['sha256'] = '0' * 64
    path.write_text(json.dumps(manifest), encoding='utf-8')
    with pytest.raises(ValueError, match='disagrees with the pinned accession'):
        verify_cohort(destination)
    manifest['cohort_digest'] = '0' * 64
    path.write_text(json.dumps(manifest), encoding='utf-8')
    with pytest.raises(ValueError, match='does not match the pinned cohort definition'):
        verify_cohort(destination)


def test_cohort_digest_is_stable_and_separates_the_two_shipped_cohorts():
    assert cohort_digest('kpneumoniae-10') == cohort_digest('kpneumoniae-10')
    assert cohort_digest('kpneumoniae-10') != cohort_digest('mixed-genus-20')
    assert re.fullmatch(r'[0-9a-f]{64}', cohort_digest('mixed-genus-20'))


def test_cancelling_a_download_leaves_no_cohort_behind(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    destination = tmp_path / 'cohorts' / name
    seen = []
    original = practice_cohorts._fetch_genome

    def spy(entry, target, **options):
        seen.append(entry['accession'])
        return original(entry, target, **options)

    monkeypatch.setattr(practice_cohorts, '_fetch_genome', spy)
    with pytest.raises(AnalysisCancelled):
        download_cohort(name, destination, source_root=mirror, cancelled=lambda: len(seen) >= 1)
    assert seen and not destination.exists()


def test_verify_pins_reads_only_published_checksums_and_reports_drift(tmp_path, monkeypatch):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    monkeypatch.setattr(practice_cohorts, '_fetch_genome',
                        lambda *a, **k: pytest.fail('--verify-pins must not download genomes'))
    report = verify_pins(name, source_root=mirror)
    assert report['stale'] is False and report['checked'] == 3
    assert report['matched'] == [row[0] for row in three_rows()]
    entry = cohort_genomes(name)[2]
    checksums = mirror / practice_cohorts.checksums_relative_path(entry)
    checksums.write_text(f'{"1" * 32}  ./{entry["assembly_dir"]}_genomic.fna.gz\n',
                         encoding='utf-8')
    drifted = verify_pins(name, source_root=mirror)
    assert drifted['stale'] is True
    assert drifted['mismatched'] == [{'accession': entry['accession'],
                                      'pinned_md5': entry['gz_md5'], 'published_md5': '1' * 32}]
    checksums.write_text('', encoding='utf-8')
    assert verify_pins(name, source_root=mirror)['missing'][0]['accession'] == entry['accession']


def test_unknown_cohorts_and_malformed_accessions_are_refused():
    with pytest.raises(ValueError, match='Unknown practice cohort'):
        cohort_genomes('no-such-cohort')
    with pytest.raises(ValueError, match='Unsupported assembly accession'):
        genome_relative_path({'accession': 'bogus', 'assembly_dir': 'bogus_v1'})


def test_download_directory_names_stay_flat_and_path_safe():
    entry = {'genus': 'Klebsiella', 'species': 'pneumoniae', 'strain': '../evil name',
             'accession': 'GCF_000240185.1'}
    filename = genome_filename(entry)
    assert filename == 'Klebsiella_pneumoniae_evil-name_GCF_000240185.1.fna.gz'
    assert '/' not in filename and '\\' not in filename and '..' not in filename
    assert Path(filename).name == filename and filename.endswith('.fna.gz')


def test_git_never_tracks_a_downloaded_practice_cohort():
    candidates = ['practice-cohorts/mixed-genus-20/manifest.json',
                  'practice-cohorts/mixed-genus-20/Klebsiella_pneumoniae_HS11286_GCF_000240185.1'
                  '.fna.gz',
                  'practice-cohorts/_cache/deadbeef.fna.gz',
                  'src/wmlstudio/resources/testdata/anything.json']
    for relative in candidates:
        checked = subprocess.run(['git', 'check-ignore', '-q', relative], cwd=REPOSITORY,
                                 capture_output=True)
        assert checked.returncode == 0, f'{relative} would be committed'


def test_cli_lists_cohorts_and_fetches_and_verifies_from_a_local_mirror(tmp_path, monkeypatch,
                                                                       capsys):
    mirror, name = mirror_cohort(tmp_path, monkeypatch, three_rows())
    monkeypatch.setattr(fetch_cli, 'cohort_names', cohort_names)
    assert fetch_cli.main(['--list']) == 0
    listed = json.loads(capsys.readouterr().out)
    assert {row['name'] for row in listed} >= {'kpneumoniae-10', 'mixed-genus-20'}
    destination = tmp_path / 'cohorts' / name
    assert fetch_cli.main(['--cohort', name, '--destination', str(destination),
                           '--source-root', str(mirror), '--quiet']) == 0
    assert json.loads(capsys.readouterr().out)['installed'] is True
    assert fetch_cli.main(['--verify', str(destination)]) == 0
    assert json.loads(capsys.readouterr().out)['genomes'] == 3
    assert fetch_cli.main(['--cohort', name, '--verify-pins', '--source-root', str(mirror),
                           '--quiet']) == 0
    entry = cohort_genomes(name)[0]
    (mirror / practice_cohorts.checksums_relative_path(entry)).write_text('', encoding='utf-8')
    capsys.readouterr()
    assert fetch_cli.main(['--cohort', name, '--verify-pins', '--source-root', str(mirror),
                           '--quiet']) == 1


@pytest.mark.realdata
@pytest.mark.parametrize('name', ['kpneumoniae-10', 'mixed-genus-20'])
def test_real_ncbi_download_matches_every_pinned_checksum(tmp_path, name):
    if os.environ.get('WMLSTUDIO_REALDATA') != '1':
        pytest.skip('Set WMLSTUDIO_REALDATA=1 to download the pinned cohorts from NCBI.')
    destination = tmp_path / 'practice-cohorts' / name
    result = download_cohort(name, destination)
    assert result['genomes'] == len(COHORTS[name]['genomes'])
    manifest = verify_cohort(destination)
    for item, entry in zip(sorted(manifest['genomes'], key=lambda row: row['accession']),
                           sorted(cohort_genomes(name), key=lambda row: row['accession']),
                           strict=True):
        payload = (destination / item['file']).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == entry['gz_sha256']
        assert hashlib.sha256(gzip.decompress(payload)).hexdigest() == entry['fasta_sha256']
