"""Independent, native whole-genome reference comparisons using pinned Pyskani.

This is not the MLST reference matcher, and does not claim Kleborate/Mash
equivalence. ANI and aligned fractions are retained with conservative review
gates. Subspecies labels are reference-supported candidates, not validated calls.
"""

from __future__ import annotations

import math
import threading
from pathlib import Path

import pyskani

from .characterization_refs import validate_characterization_references
from .sequence import SequenceReader, check_cancelled, file_sha256, file_signature, sample_name

_CACHE = {}
_LOCK = threading.Lock()


def _sequences(path, cancelled=None):
    sequences, seen = [], set()
    with SequenceReader(path, cancelled) as reader:
        if reader.kind != 'fasta':
            raise ValueError('Independent species characterization requires an assembled FASTA.')
        for record in reader:
            if record.identifier in seen:
                raise ValueError(f'Duplicate contig identifier: {record.identifier}')
            seen.add(record.identifier)
            sequences.append(record.sequence.encode('ascii'))
    return sequences


def runtime_provenance():
    return {'engine': 'pyskani', 'version': pyskani.__version__, 'skani_version': pyskani.SKANI_VERSION,
            'binary_sha256': file_sha256(pyskani._skani.__file__),
            'method': 'whole-genome ANI and bidirectional aligned fractions; independent of MLST',
            'sketch_parameters': {'compression': 125, 'marker_compression': 1000, 'k': 15},
            'query_parameters': {'learned_ani': True, 'median': False, 'robust': False, 'faster_small': False},
            'source': 'https://github.com/althonos/pyskani',
            'citation': 'https://doi.org/10.1038/s41592-023-02018-3'}


def _database(root, cancelled=None, progress=None):
    root = Path(root).resolve()
    # Reference revisions are immutable. Check file signatures before reuse;
    # rebuilding verifies every source hash rather than trusting metadata alone.
    with _LOCK:
        cached = _CACHE.get(str(root))
        if cached and all(file_signature(root / relative) == signature for relative, signature in cached[2]):
            return cached[0], cached[1]
    manifest = validate_characterization_references(root, cancelled=cancelled)
    if not manifest['species'] or not any(item.get('outgroup') for item in manifest['species']):
        raise ValueError('The species panel must include explicit outgroup references.')
    signatures = [(entry['path'], file_signature(root / entry['path'])) for entry in manifest['files']]
    signatures.append(('manifest.json', file_signature(root / 'manifest.json')))
    database = pyskani.Database(compression=125, marker_compression=1000, k=15)
    for index, entry in enumerate(manifest['species'], 1):
        check_cancelled(cancelled)
        database.sketch(entry['id'], *_sequences(root / entry['path'], cancelled))
        if progress:
            progress(index, len(manifest['species']), f'Sketched species reference {entry["id"]}')
    check_cancelled(cancelled)
    if any(file_signature(root / relative) != signature for relative, signature in signatures):
        raise ValueError('Species reference files changed while creating the sketch database.')
    with _LOCK:
        _CACHE.clear()  # One bounded installed panel, shared across a batch.
        _CACHE[str(root)] = (database, manifest, signatures)
    return database, manifest


def interpret_species_hits(hits, references, *, min_ani=95.0, min_fraction=0.65, min_gap=0.5):
    """Pure decision kernel; percentages for ANI, fractions in [0,1] for AF."""
    if not (math.isfinite(min_ani) and 90 <= min_ani <= 100
            and math.isfinite(min_fraction) and 0.5 <= min_fraction <= 1
            and math.isfinite(min_gap) and 0 < min_gap <= 10):
        raise ValueError('Invalid species evidence thresholds.')
    by_id = {entry['id']: entry for entry in references}
    rows = []
    for hit in hits:
        if hit['reference_id'] not in by_id:
            raise ValueError('Species output refers to an unknown reference.')
        ani, qf, rf = (hit[key] for key in ('ani', 'query_fraction', 'reference_fraction'))
        if any(not math.isfinite(value) for value in (ani, qf, rf)) or not (0 <= ani <= 100 and 0 <= qf <= 1 and 0 <= rf <= 1):
            raise ValueError('Species output contains invalid identity or aligned fractions.')
        rows.append(dict(by_id[hit['reference_id']], **hit))
    rows.sort(key=lambda row: (-row['ani'], -min(row['query_fraction'], row['reference_fraction']), row['reference_id']))
    adequate = [row for row in rows if min(row['query_fraction'], row['reference_fraction']) >= min_fraction]
    result = {'status': 'ambiguous', 'genus': '', 'species': '', 'subspecies': '',
              'subspecies_status': 'unresolved', 'confidence': 'unresolved',
              'nearest': rows[0] if rows else None, 'runner_up': None, 'gap_ani': None,
              'hits': rows, 'thresholds': {'min_ani_pct': min_ani, 'min_aligned_fraction': min_fraction, 'min_species_gap_ani_pct': min_gap},
              'reason': 'No reference comparison met both aligned-fraction gates.'}
    if not adequate:
        return result
    best = adequate[0]
    taxon = (best['genus'], best['species'])
    alternatives = [row for row in adequate if (row['genus'], row['species']) != taxon]
    runner = alternatives[0] if alternatives else None
    gap = best['ani'] - runner['ani'] if runner else None
    result.update(nearest=best, runner_up=runner, gap_ani=gap)
    discordant = [row for row in rows if (row['genus'], row['species']) != taxon and row['ani'] >= min_ani
                  and 0.1 <= row['query_fraction'] < min_fraction and row['reference_fraction'] >= min_fraction]
    result['discordant_low_fraction_hits'] = discordant
    if best['ani'] < min_ani:
        result['reason'] = 'Nearest adequate reference is below the configured species ANI threshold.'
    elif runner and gap < min_gap:
        result['reason'] = 'Competing species references are too close; review possible hybrid, mixture or reference-panel limitations.'
    elif discordant:
        result['reason'] = 'A second species reference is covered over a substantial minority of the assembly; review a possible mixture before accepting species identity.'
    else:
        result.update(status='completed', genus=best['genus'], species=best['species'],
                      confidence='genomic_reference_supported', reason='Nearest species meets ANI, aligned-fraction and competing-species separation gates.')
        if not best['species']:
            result.update(status='ambiguous', confidence='genus_only', reason='The closest outgroup reference is labelled only to genus; no species name was inferred.')
        if taxon == ('Escherichia', 'coli'):
            result.update(status='ambiguous', confidence='complex_only',
                          reason='ANI cannot reliably distinguish Escherichia coli from Shigella; report the E. coli/Shigella complex.')
        if best.get('subspecies'):
            competitors = [row for row in adequate if (row['genus'], row['species']) == taxon
                           and row.get('subspecies') and row['subspecies'] != best['subspecies']]
            if competitors and best['ani'] >= 98 and best['ani'] - competitors[0]['ani'] >= 1:
                result['subspecies'] = best['subspecies']
                result['subspecies_status'] = 'provisional_reference_match'
    return result


def identify_species(path, reference_root, cancelled=None, progress=None, *, min_ani=95.0, min_fraction=0.65, min_gap=0.5):
    """Return independent genomic species evidence; never silently fall back to MLST."""
    check_cancelled(cancelled)
    path = Path(path).resolve()
    signature = file_signature(path)
    input_hash = file_sha256(path, cancelled)
    sequences = _sequences(path, cancelled)
    database, manifest = _database(reference_root, cancelled, progress)
    if progress:
        progress(0, 0, 'Comparing whole-genome ANI and aligned fractions against the independent reference panel…')
    check_cancelled(cancelled)
    hits = database.query(sample_name(path), *sequences, learned_ani=True, median=False, robust=False, faster_small=False)
    check_cancelled(cancelled)
    normalized = [{'reference_id': hit.reference_name, 'ani': hit.identity * 100,
                   'query_fraction': hit.query_fraction, 'reference_fraction': hit.reference_fraction} for hit in hits]
    result = interpret_species_hits(normalized, manifest['species'], min_ani=min_ani, min_fraction=min_fraction, min_gap=min_gap)
    if file_signature(path) != signature or file_sha256(path, cancelled) != input_hash:
        raise ValueError('Assembly changed during species characterization; no result was accepted.')
    result.update(input_sha256=input_hash, reference_digest=manifest['reference_digest'],
                  reference_count=len(manifest['species']), provenance=runtime_provenance(),
                  limitations=[*manifest.get('limitations', []),
                               'ANI is independent of MLST but does not exclude minority contamination or establish clinical species confirmation.',
                               'Subspecies ANI separation is a conservative review heuristic, not a universal taxonomic cutoff.',
                               'Pneumoniae subspecies ozaenae/rhinoscleromatis are not inferred from an ST or from this panel.',
                               'The in-process sketch query is cancellable immediately before/after the native call, not during that individual native call.'])
    return result
