"""Cancellable native SKA2 with immutable input identities and separate metrics.

SKA Distance is a split-kmer SNP count; Mismatch count is presence/absence and
must never be relabelled SNPs. Reference-mapped/filter-specific distances are
separate evidence and cannot inherit a cutoff from the reference-free run.
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

from .sequence import SequenceReader, check_cancelled, file_sha256, file_signature

VERSION = '0.5.1'
SOURCE_REVISION = 'fcf9413d2768dc6538d31f664a4bf310651449e7'


def runtime_capabilities(root=None):
    platform = 'windows-x64' if os.name == 'nt' else 'linux-x64'
    binary_name = 'ska.exe' if os.name == 'nt' else 'ska'
    candidates = [Path(root)] if root else [
        Path(__file__).resolve().parent / 'resources/tools/ska2' / platform,
        Path(getattr(sys, '_MEIPASS', Path(sys.executable).parent)) / 'Tools/ska2',
        Path(sys.executable).parent / 'Tools/ska2']
    for directory in candidates:
        binary = directory / binary_name
        if not binary.is_file():
            continue
        try:
            manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))
            expected = next(e['sha256'] for e in manifest['files'] if e['path'] == binary_name)
            if manifest['version'] != VERSION or manifest['source_commit'] != SOURCE_REVISION or manifest['platform'] != platform:
                raise ValueError('Wrong SKA2 version/source/platform')
        except (OSError, KeyError, ValueError, StopIteration) as error:
            return {'available': False, 'reason': f'Invalid SKA2 manifest: {error}'}
        return {'available': True, 'binary': str(binary.resolve()), 'expected_sha256': expected,
                'version': VERSION, 'platform': platform, 'source_revision': SOURCE_REVISION}
    return {'available': False, 'reason': 'Native SKA2 is not staged in this package. No download occurs during analysis.'}


def _run(command, directory, stem, cancelled=None, progress=None):
    check_cancelled(cancelled)
    out, err = directory / f'{stem}.stdout.txt', directory / f'{stem}.stderr.txt'
    with out.open('wb') as stdout, err.open('wb') as stderr:
        process = subprocess.Popen(command, cwd=directory, stdout=stdout, stderr=stderr,
            stdin=subprocess.DEVNULL, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        try:
            while True:
                check_cancelled(cancelled)
                try:
                    code = process.wait(timeout=0.2)
                    break
                except subprocess.TimeoutExpired:
                    pass
        except BaseException:
            process.kill()
            process.wait(timeout=10)
            raise
    if code:
        with err.open('rb') as handle:
            handle.seek(max(0, err.stat().st_size - 8192))
            message = handle.read().decode('utf-8', errors='replace')
        raise ValueError(f'SKA2 {stem} failed ({code}): {message}')
    check_cancelled(cancelled)
    if progress:
        progress(0, 0, f'SKA2 {stem} completed')
    return out


def parse_distances(path, aliases, min_shared_fraction=0.95):
    if not 0 <= min_shared_fraction <= 1 or not math.isfinite(min_shared_fraction):
        raise ValueError('Minimum shared split-kmer fraction must be between 0 and 1.')
    expected = {frozenset(pair) for pair in combinations(aliases, 2)}
    seen, rows = set(), []
    with Path(path).open(encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Sample1', 'Sample2', 'Distance', 'Mismatches (proportion)', 'Match count', 'Mismatch count'}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError('SKA2 distance output has an unsupported header.')
        for row in reader:
            a, b = row['Sample1'], row['Sample2']
            pair = frozenset((a, b))
            if pair not in expected or pair in seen:
                raise ValueError('Unknown, duplicate or self pair in SKA2 output.')
            seen.add(pair)
            distance, matched, missing = float(row['Distance']), int(row['Match count']), int(row['Mismatch count'])
            if not math.isfinite(distance) or distance < 0 or not distance.is_integer() or min(matched, missing) < 0 or distance > matched:
                raise ValueError('Non-integer or invalid unambiguous SKA2 SNP evidence.')
            fraction = matched / (matched + missing) if matched + missing else 0
            comparable = matched > 0 and fraction >= min_shared_fraction
            rows.append({'source': aliases[a]['sample_id'], 'target': aliases[b]['sample_id'],
                'source_name': aliases[a]['sample_name'], 'target_name': aliases[b]['sample_name'],
                'distance': int(distance) if comparable else None, 'observed_snp_count': int(distance),
                'shared_split_kmers': matched, 'unshared_split_kmers': missing,
                'shared_fraction': fraction, 'comparable': comparable,
                'reason': '' if comparable else 'Insufficient shared unambiguous split-kmers; no accepted distance.'})
    if seen != expected:
        raise ValueError('SKA2 output is missing cohort pair comparisons.')
    return rows


def _alignment(path, aliases, reference_length):
    """Mapped alignments may contain gaps, unlike ordinary sequence FASTA."""
    records, identifier, parts = {}, None, []
    with Path(path).open('rb') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(b'>'):
                if identifier is not None:
                    records[identifier] = b''.join(parts)
                identifier, parts = line[1:].decode('ascii'), []
                if identifier not in aliases or identifier in records:
                    raise ValueError('Unknown or duplicate SKA2 mapped sample name.')
            elif line:
                if identifier is None or set(line) - set(b'ACGTN-'):
                    raise ValueError('Unexpected unmapped/ambiguous SKA2 alignment format.')
                parts.append(line)
        if identifier is not None:
            records[identifier] = b''.join(parts)
    if set(records) != set(aliases) or any(len(sequence) != reference_length for sequence in records.values()):
        raise ValueError('SKA2 alignment is not a complete reference-coordinate cohort.')
    return records


def _mapped_evidence(binary, work, aliases, reference, *, threads, repeat_mask, annotation_path,
                     annotation_reference_sha256, mask_bed, coding_only, min_shared_fraction, cancelled, progress):
    from .snp_annotation import annotate_snps
    contigs, bases = [], bytearray()
    with SequenceReader(reference, cancelled) as reader:
        if reader.kind != 'fasta':
            raise ValueError('Reference must be FASTA.')
        for record in reader:
            contigs.append((record.identifier, len(bases), len(record.sequence)))
            bases.extend(record.sequence.encode('ascii'))
            if len(bases) > 20_000_000:
                raise ValueError('Reference exceeds the explicit 20 million-base bound.')
    # Raw and masked outputs both retained; masking never destroys evidence.
    command = [binary, 'map', str(Path(reference).resolve()), 'cohort.skf', '--ambig-mask', '--threads', str(threads)]
    raw_path = _run(command, work, 'mapped-unfiltered', cancelled, progress)
    raw = _alignment(raw_path, aliases, len(bases))
    if repeat_mask:
        filtered_path = _run([*command, '--repeat-mask'], work, 'mapped-repeat-masked', cancelled, progress)
        mapped = _alignment(filtered_path, aliases, len(bases))
    else:
        mapped = raw
    # Keep coordinates for every canonical ALT site, including sites masked later.
    variants, offsets = [], []
    for contig, start, length in contigs:
        for local in range(length):
            index = start + local
            if index % 10000 == 0:
                check_cancelled(cancelled)
            ref = chr(bases[index])
            if ref not in 'ACGT':
                continue
            alternatives = sorted({chr(sequence[index]) for sequence in raw.values()} - {ref, 'N', '-'})
            if alternatives:
                variants.append({'contig': contig, 'pos': local + 1, 'ref': ref, 'alt': alternatives,
                    'sample_bases': {aliases[key]['sample_id']: chr(sequence[index]) for key, sequence in raw.items()}})
                offsets.append(index)
                if len(variants) > 100000:
                    raise ValueError('Over 100,000 variant sites: review outliers or narrow the SNP cohort.')
    annotation = annotate_snps(reference, variants, annotation_path=annotation_path,
        annotation_reference_sha256=annotation_reference_sha256, mask_bed=mask_bed,
        coding_only=coding_only, cancelled=cancelled, progress=progress)
    selected = [i for i, site in enumerate(annotation['sites']) if site['included']]
    eligible_bytes = bytearray(len(bases))
    starts = {name: start for name, start, length in contigs}
    for interval in annotation['protocol']['eligible_reference_intervals']:
        start = starts[interval['contig']] + interval['start']
        end = starts[interval['contig']] + interval['end']
        eligible_bytes[start:end] = b'\x01' * (end - start)
    canonical = bytes(1 if i in b'ACGT' else 0 for i in range(256))
    # One bit per byte in big integers gives C-level pairwise intersection/count,
    # avoiding a Python loop across the complete reference for every pair.
    eligible_bits = int.from_bytes(eligible_bytes, 'little') & int.from_bytes(bytes(bases).translate(canonical), 'little')
    eligible_count = eligible_bits.bit_count()
    callable_bits = {key: int.from_bytes(sequence.translate(canonical), 'little') & eligible_bits
                     for key, sequence in mapped.items()}
    rows = []
    for a, b in combinations(aliases, 2):
        check_cancelled(cancelled)
        left, right = mapped[a], mapped[b]
        observed, shared_variants = 0, 0
        for i in selected:
            position = offsets[i]
            x, y = left[position], right[position]
            if x in b'ACGT' and y in b'ACGT':
                observed += x != y
                shared_variants += 1
        shared_bases = (callable_bits[a] & callable_bits[b]).bit_count()
        shared_fraction = shared_bases / eligible_count if eligible_count else 0
        comparable = shared_bases > 0 and shared_fraction >= min_shared_fraction
        rows.append({'source': aliases[a]['sample_id'], 'target': aliases[b]['sample_id'],
            'source_name': aliases[a]['sample_name'], 'target_name': aliases[b]['sample_name'],
            'observed_snp_count': observed, 'distance': observed if comparable else None,
            'comparable': comparable, 'shared_reference_bases': shared_bases,
            'eligible_reference_bases': eligible_count, 'shared_fraction': shared_fraction,
            'shared_selected_variant_sites': shared_variants,
            'reason': '' if comparable else 'Insufficient callable reference overlap after selected filters.'})
    for site, offset in zip(annotation['sites'], offsets):
        site['repeat_masked_sample_ids'] = [aliases[key]['sample_id'] for key in aliases
            if raw[key][offset] in b'ACGT' and mapped[key][offset] not in b'ACGT']
    return {'method': 'reference-mapped-filtered-SNPs', 'rows': rows,
        'annotation': annotation, 'reference_bases': len(bases), 'reference_contigs': contigs,
        'reference_sha256': file_sha256(reference, cancelled), 'repeat_mask': repeat_mask,
        'distance_status': 'callable_overlap_filtered',
        'minimum_shared_fraction': min_shared_fraction,
        'limitations': ['Filtered reference-mapped distances are distinct from reference-free split-kmer counts and require their own threshold validation.',
            'SKA repeat-mask masks repeated reference split-kmers; it is not complete IS/mobile-element detection.',
            'User-supplied BED masks are separate from repeat masking; recombination was not inferred.']}


def run_ska(samples, output_root, *, threads=4, k=31, min_shared_fraction=0.95,
            reference_path=None, annotation_path=None, annotation_reference_sha256=None,
            mask_bed=None, coding_only=False, repeat_mask=True, root=None, cancelled=None, progress=None):
    samples = list(samples)
    if not 2 <= len(samples) <= 200 or len({s['id'] for s in samples}) != len(samples):
        raise ValueError('Choose 2–200 distinct assembled isolates for a bounded SNP cohort.')
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ValueError('Threads must be an integer between 1 and 256.')
    if isinstance(k, bool) or not isinstance(k, int) or k < 5 or k > 63 or k % 2 != 1:
        raise ValueError('SKA k must be an odd integer from 5 to 63.')
    if not 0 <= min_shared_fraction <= 1 or not math.isfinite(min_shared_fraction):
        raise ValueError('Minimum shared split-kmer fraction must be between 0 and 1.')
    if (annotation_path or mask_bed or coding_only) and not reference_path:
        raise ValueError('Coding annotation and interval masks require an exact reference FASTA.')
    capability = runtime_capabilities(root)
    if not capability['available']:
        raise ValueError(capability['reason'])
    binary = capability['binary']
    if file_sha256(binary, cancelled) != capability['expected_sha256']:
        raise ValueError('Native SKA2 binary checksum mismatch.')
    work = Path(output_root).resolve() / ('ska2-' + uuid.uuid4().hex)
    work.mkdir(parents=True, exist_ok=False)
    aliases, tracked, total = {}, {}, 0
    try:
        for index, sample in enumerate(sorted(samples, key=lambda item: item['id'])):
            check_cancelled(cancelled)
            source = Path(sample['input_path']).resolve()
            before = file_signature(source)
            digest = file_sha256(source, cancelled)
            expected = sample.get('input_sha256') or (sample.get('result') or {}).get('input_sha256')
            if not expected or digest != expected:
                raise ValueError('SNP inputs must match their reviewed analysis SHA-256; analyse changed/unverified inputs first.')
            alias = f'isolate{index:04d}'
            with SequenceReader(source, cancelled) as reader, (work / f'{alias}.fa').open('w', encoding='ascii') as out:
                if reader.kind != 'fasta':
                    raise ValueError('This SKA2 workflow accepts assemblies; assemble and review FASTQ inputs first.')
                for record_index, record in enumerate(reader):
                    total += len(record.sequence)
                    if total > 1_000_000_000:
                        raise ValueError('SNP cohort exceeds the explicit one billion input-base bound.')
                    out.write(f'>contig{record_index}\n{record.sequence}\n')
            if before != file_signature(source) or file_sha256(source, cancelled) != digest:
                raise ValueError('An original input changed while preparing SKA2.')
            tracked[str(source)] = (before, digest)
            aliases[alias] = {'sample_id': sample['id'], 'sample_name': sample['name'],
                'input_path': str(source), 'input_sha256': digest, 'staged_file': f'{alias}.fa',
                'staged_sha256': file_sha256(work / f'{alias}.fa', cancelled)}
            if progress:
                progress(index + 1, len(samples), f'Verified SNP input {sample["name"]}')
        for path in (reference_path, annotation_path, mask_bed):
            if path:
                path = str(Path(path).resolve())
                tracked[path] = (file_signature(path), file_sha256(path, cancelled))
        (work / 'inputs.tsv').write_text(''.join(f'{alias}\t{alias}.fa\n' for alias in aliases), encoding='ascii')
        _run([binary, 'build', '-f', 'inputs.tsv', '-o', 'cohort', '-k', str(k), '--threads', str(threads)], work, 'build', cancelled, progress)
        distance_path = _run([binary, 'distance', 'cohort.skf', '--threads', str(threads)], work, 'distance', cancelled, progress)
        rows = parse_distances(distance_path, aliases, min_shared_fraction)
        mapped = _mapped_evidence(binary, work, aliases, reference_path, threads=threads,
            repeat_mask=repeat_mask, annotation_path=annotation_path, annotation_reference_sha256=annotation_reference_sha256,
            mask_bed=mask_bed, coding_only=coding_only, min_shared_fraction=min_shared_fraction,
            cancelled=cancelled, progress=progress) if reference_path else None
        for path, (signature, digest) in tracked.items():
            if file_signature(path) != signature or file_sha256(path, cancelled) != digest:
                raise ValueError('An input/reference/annotation/mask changed during SNP analysis; result rejected.')
        result = {'format_version': 1, 'status': 'completed', 'engine': 'SKA2', 'version': VERSION,
            'run_id': work.name, 'created_at': datetime.now(timezone.utc).isoformat(),
            'method': 'reference-free-assembly-split-kmer-SNPs', 'parameters': {'k': k, 'threads': threads,
                'ambiguous_bases': 'excluded', 'min_frequency': 0.0, 'minimum_shared_fraction': min_shared_fraction},
            'binary_sha256': capability['expected_sha256'], 'source_revision': SOURCE_REVISION,
            'inputs': list(aliases.values()), 'rows': rows, 'mapped': mapped,
            'output_directory': str(work), 'threshold': None,
            'limitations': ['SNP similarity is not direct transmission or direction of spread.',
                'Reference-free counts are not reference-mapped, coding-only, repeat-masked or recombination-filtered counts.',
                'SKA2 assembly results cannot inherit the Higgs 2022 original-SKA read-based cutoff.',
                'Shared split-kmer fraction is a comparison filter, not genome coverage or an assembly quality pass.']}
        from .export import _atomic_text
        with _atomic_text(work / 'result.json') as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False, allow_nan=False)
        return result
    except BaseException as error:
        # Preserve execution logs for diagnosis, never present partial output as a result.
        (work / 'FAILED.txt').write_text(f'{type(error).__name__}: {error}\n', encoding='utf-8')
        raise
