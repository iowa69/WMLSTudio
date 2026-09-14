"""Explicit opt-in real local positive control; no data downloads or uploads."""

import json
import os
import time
from pathlib import Path

import pytest

from wmlstudio.read_support import investigate_read_support
from wmlstudio.typing import call_assembly, load_scheme


@pytest.mark.skipif(os.environ.get('WMLSTUDIO_REAL_READ_SUPPORT') != '1', reason='Explicit local real-read positive-control opt-in')
def test_real_srr12343864_arcC4_positive_control(tmp_path):
    root = Path(__file__).resolve().parents[1]
    reads = [Path('/media/iowa/u/tesseract_fastq/saureus') / f'SRR12343864_{mate}.fastq.gz' for mate in (1, 2)]
    assembly = root / 'artifacts/assembly-validation/SRR12343864-windows-skesa240/contigs.fasta'
    tools = Path(os.environ['WMLSTUDIO_TEST_BLAST_DIR'])
    suffix = '.exe' if os.name == 'nt' else ''
    assert all(path.is_file() for path in [*reads, assembly])
    scheme = load_scheme(root / 'src/wmlstudio/resources/schemes/saureus')
    original = call_assembly(assembly, scheme)
    assert original['st'] == '20' and original['alleles']['arcC'] == '4'
    started = time.monotonic()
    result = investigate_read_support(*reads, scheme, ['arcC'], assembly_path=assembly, max_pairs=100_000,
        blastn_path=tools / ('blastn' + suffix), makeblastdb_path=tools / ('makeblastdb' + suffix))
    result['validation_control'] = {'role': 'Real known-positive control, not a genuinely missing locus',
                                    'accession': 'SRR12343864', 'expected_assembly_st': '20', 'expected_arcC': '4',
                                    'elapsed_seconds': time.monotonic() - started}
    artifact = tmp_path / 'real-read-support.json'
    artifact.write_text(json.dumps(result, indent=2))
    print(f'Full separate read-support result: {artifact}')
    locus = result['loci'][0]
    candidate = next(row for row in locus['candidates'] if row['allele'] == '4')
    assert locus['status'] == 'supported' and locus['compatible_candidates'] == ['4']
    assert locus['assigned_allele'] is None
    assert candidate['breadth_min_depth'] == 1 and candidate['observed_identity_pct'] == 100
    assert candidate['min_depth'] >= 3 and candidate['mixed_positions'] == 0
    assert result['sampling']['sampled'] and result['sampling']['pairs_checked'] == 100_000
    assert result['panel']['all_alleles_in_selected_loci_tested']
