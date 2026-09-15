"""Bounded read-to-reference investigation, never an assembly allele caller.

This is an explicit second-pass BLASTN assay on selected loci. Both mates are
aligned independently; depths count reads, not independent DNA molecules. Every
allele in the selected, bounded panel remains in the result, including zero-hit
alleles. No assembly call, official ST, or distance profile is changed.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import tempfile
from array import array
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

from . import __version__
from .cgtyping import _run, resolve_blast
from .sample_workflow import _sha256, current_input_sha256
from .sequence import (
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
    read_pair_identity,
)
from .typing import Scheme, reverse_complement

# The explicit panel bound this assay enforces. Named here so a caller that has to
# split a long list of loci into runs uses the assay's own number instead of
# repeating a literal that could drift away from it.
MAX_PANEL_LOCI = 20

LIMITATIONS = [
    'Read support is separate investigative evidence: no assembly allele, ST or cgMLST distance is overwritten.',
    'Depth counts independently aligned reads, including overlapping mates and PCR duplicates; it is not molecule depth.',
    'A read can support several similar alleles/loci. Candidate depth is non-unique and not summed across alternatives.',
    'No support in a sampled prefix or a limited reference panel does not establish locus absence.',
    'Local alignments do not reconstruct or phase a complete allele; mixed-base signals can reflect repeats, contamination or sequencing error.',
    'Reads are screened as supplied, without adapter trimming or pair-aware mapping. Low-quality bases are excluded explicitly.',
]


def _integer(value, name, low, high):
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ValueError(f'{name} must be an integer between {low:,} and {high:,}.')


def preflight_read_support(scheme, loci, *, max_loci=MAX_PANEL_LOCI, max_alleles=5000, max_reference_bases=10_000_000):
    """Inspect an already-loaded scheme without reading FASTQ or spawning tools."""
    if not isinstance(scheme, Scheme):
        raise ValueError('Load and validate a local scheme before read-support preflight.')
    _integer(max_loci, 'Maximum loci', 1, 100)
    _integer(max_alleles, 'Maximum alleles', 1, 10_000)
    _integer(max_reference_bases, 'Maximum reference bases', 1, 10_000_000)
    selected = list(loci)
    if not selected or len(selected) != len(set(selected)) or any(locus not in scheme.alleles for locus in selected):
        raise ValueError('Choose a nonempty, unique list of loci from this scheme.')
    counts = {locus: len(scheme.alleles[locus]) for locus in selected}
    total_bases = sum(len(sequence) for locus in selected for sequence in scheme.alleles[locus].values())
    if len(selected) > max_loci or sum(counts.values()) > max_alleles or total_bases > max_reference_bases:
        raise ValueError(f'Read-support panel exceeds the explicit safety bounds: {len(selected):,} loci, '
                         f'{sum(counts.values()):,} alleles, {total_bases:,} reference bases. '
                         'Select fewer loci; alleles are never silently truncated or selected by proximity.')
    if not total_bases or any(not count for count in counts.values()):
        raise ValueError('Every selected locus must contain at least one nonempty allele.')
    # Digest the actual in-memory sequences used, in addition to the whole scheme.
    digest = hashlib.sha256()
    for locus in sorted(selected):
        for allele, sequence in sorted(scheme.alleles[locus].items()):
            if not sequence or set(sequence) - set('ACGT'):
                raise ValueError('Read-support references must contain only unambiguous A/C/G/T bases.')
            digest.update(json.dumps([locus, allele, sequence], separators=(',', ':')).encode())
            digest.update(b'\n')
    return {'loci': selected, 'locus_count': len(selected), 'allele_count': sum(counts.values()),
            'alleles_per_locus': counts, 'reference_bases': total_bases,
            'panel_sha256': digest.hexdigest(), 'all_alleles_in_selected_loci_tested': True,
            'bounds': {'max_loci': max_loci, 'max_alleles': max_alleles, 'max_reference_bases': max_reference_bases}}


def _read_queries(paths, output, *, max_pairs, max_query_bases, min_read_length, max_read_length,
                  cancelled, progress):
    queries, bases, pairs, explicit = {}, 0, 0, True
    with SequenceReader(paths[0], cancelled) as left, SequenceReader(paths[1], cancelled) as right, output.open('w', encoding='ascii') as out:
        if left.kind != 'fastq' or right.kind != 'fastq':
            raise ValueError('Read support requires two FASTQ mates, not assembly FASTA inputs.')
        first, second = iter(left), iter(right)
        for index in range(max_pairs):
            check_cancelled(cancelled)
            a, b = next(first, None), next(second, None)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise ValueError('FASTQ files have different numbers of records in the investigated prefix.')
            aid, am = read_pair_identity(a.name)
            bid, bm = read_pair_identity(b.name)
            if aid != bid or am not in (None, 1) or bm not in (None, 2) or (am is None) != (bm is None):
                raise ValueError(f'FASTQ mate identity or direction disagrees at read pair {index + 1}.')
            explicit = explicit and am == 1 and bm == 2
            for mate, record in enumerate((a, b), 1):
                if not min_read_length <= len(record.sequence) <= max_read_length:
                    raise ValueError(f'Read pair {index + 1} contains a read outside the explicit '
                                     f'{min_read_length}–{max_read_length} bp short-read assay bounds; no reads were silently dropped.')
                bases += len(record.sequence)
                if bases > max_query_bases:
                    raise ValueError('Selected reads exceed the query-base safety bound. Reduce the explicit pair limit; no partial result was accepted.')
                identifier = f'r{index}_{mate}'
                queries[identifier] = record
                out.write(f'>{identifier}\n{record.sequence}\n')
            pairs += 1
            if progress and pairs % 10_000 == 0:
                progress(pairs, max_pairs, f'Preparing read-support prefix · {pairs:,} pairs')
        complete = left.complete and right.complete
        if left.complete != right.complete:
            raise ValueError('FASTQ mates end at different points in the investigated prefix.')
    if not pairs:
        raise ValueError('The selected FASTQ inputs contain no read pairs.')
    return queries, {'pairs_checked': pairs, 'reads_aligned': pairs * 2, 'query_bases': bases,
                     'max_pairs': max_pairs, 'method': 'first N paired records, not random sampling',
                     'complete_files': complete, 'sampled': not complete, 'explicit_mate_markers': explicit,
                     'unexamined_tail_structurally_validated': complete}


def _aligned_bases(row, record, reference, min_base_quality):
    """Recover high-quality oriented observations, excluding deletion depth."""
    qstart, qend, sstart, send = (int(row[key]) for key in ('qstart', 'qend', 'sstart', 'send'))
    qseq, sseq = row['qseq'].upper(), row['sseq'].upper()
    if not (1 <= qstart <= qend <= len(record.sequence) and 1 <= min(sstart, send) <= max(sstart, send) <= len(reference)
            and len(qseq) == len(sseq) and len(qseq.replace('-', '')) == qend - qstart + 1
            and len(sseq.replace('-', '')) == abs(send - sstart) + 1
            and set(qseq + sseq) <= set('ACGTRYSWKMBDHVN-')):
        raise ValueError('BLAST returned invalid read-support alignment coordinates or strings.')
    if qseq.replace('-', '') != record.sequence[qstart - 1:qend]:
        raise ValueError('BLAST read alignment does not match the actual queried sequence.')
    expected = reference[min(sstart, send) - 1:max(sstart, send)]
    if sstart > send:
        expected = reverse_complement(expected)
    if sseq.replace('-', '') != expected:
        raise ValueError('BLAST reference alignment does not match the tested allele sequence.')
    observations, deletions, insertions = {}, set(), set()
    qi, si, step = qstart - 1, sstart - 1, 1 if send >= sstart else -1
    for qb, sb in zip(qseq, sseq, strict=True):
        if qb == '-' and sb == '-':
            raise ValueError('BLAST alignment contains a double gap.')
        quality = ord(record.quality[qi]) - 33 if qb != '-' else None
        if sb == '-':
            if qb in 'ACGT' and quality >= min_base_quality:
                insertions.add(max(0, min(len(reference) - 1, si)))
        elif qb == '-':
            # An inferred deletion is not an observed aligned base and adds no depth.
            # Require measured high-quality bases on both sides of the gap;
            # low-quality alignments cannot create apparently strong disruption.
            if (0 < qi < len(record.sequence) and record.sequence[qi - 1] in 'ACGT' and record.sequence[qi] in 'ACGT'
                    and ord(record.quality[qi - 1]) - 33 >= min_base_quality
                    and ord(record.quality[qi]) - 33 >= min_base_quality):
                deletions.add(si)
        elif qb in 'ACGT' and quality >= min_base_quality:
            observations[si] = qb if step == 1 else reverse_complement(qb)
        if qb != '-':
            qi += 1
        if sb != '-':
            si += step
    return observations, deletions, insertions


def _summarize_candidate(entry, counts, *, min_depth, min_breadth):
    sequence = entry['sequence']
    summary = {key: value for key, value in entry.items() if key != 'sequence'}
    summary.update(length=len(sequence), sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest())
    if counts is None:
        return dict(summary, status='no_support', supporting_reads=0, uniquely_best_reads=0, tied_best_reads=0,
                    mean_depth=0.0, min_depth=0, breadth_1x=0.0, breadth_min_depth=0.0,
                    observed_identity_pct=None, discordant_positions=0, mixed_positions=0,
                    deletion_positions=0, insertion_anchors=0, review_positions=[], review_positions_total=0)
    length = len(sequence)
    vectors, deletions, insertions = counts['bases'], counts['deletions'], counts['insertions']
    depth = [sum(vector[index] for vector in vectors) for index in range(length)]
    total = sum(depth)
    matching = sum(vectors['ACGT'.index(base)][index] for index, base in enumerate(sequence))
    strong, mixed, review = [], [], []
    for index, (base, d) in enumerate(zip(sequence, depth, strict=True)):
        ref_count = vectors['ACGT'.index(base)][index]
        alternate = max(vectors[k][index] for k in range(4) if 'ACGT'[k] != base)
        is_mixed = ref_count >= min_depth and alternate >= min_depth and ref_count / max(1, d) >= .2 and alternate / max(1, d) >= .2
        is_discordant = d >= min_depth and alternate / max(1, d) >= .8
        if is_mixed:
            mixed.append(index + 1)
        if is_discordant:
            strong.append(index + 1)
        if is_mixed or is_discordant or deletions[index] >= min_depth or insertions.get(index, 0) >= min_depth:
            review.append({'position': index + 1, 'reference_base': base,
                           'base_counts': {letter: vectors[k][index] for k, letter in enumerate('ACGT')},
                           'deletion_reads': deletions[index], 'insertion_reads': insertions.get(index, 0)})
    deletion_positions = sum(value >= min_depth for value in deletions)
    insertion_positions = sum(value >= min_depth for value in insertions.values())
    breadth = sum(d >= min_depth for d in depth) / length
    status = ('mixed_support' if mixed else 'discordant_support' if strong or deletion_positions or insertion_positions
              else 'supported' if breadth >= min_breadth else 'incomplete_support' if total else 'no_support')
    return dict(summary, status=status, supporting_reads=counts['reads'], uniquely_best_reads=counts['unique'],
                tied_best_reads=counts['tied'], mean_depth=total / length, min_depth=min(depth),
                breadth_1x=sum(d > 0 for d in depth) / length, breadth_min_depth=breadth,
                observed_identity_pct=100 * matching / total if total else None,
                discordant_positions=len(strong), mixed_positions=len(mixed),
                deletion_positions=deletion_positions, insertion_anchors=insertion_positions,
                review_positions=review[:50], review_positions_total=len(review))


def _summarize_output(output, queries, references, *, min_identity, min_read_coverage, min_base_quality,
                      min_depth, min_breadth, cancelled):
    fields = 'qseqid sseqid pident qstart qend sstart send bitscore qseq sseq'.split()
    counts, seen_queries, discarded_hsps, low_filter_hsps, aligned_reads = {}, set(), 0, 0, 0

    def consume(qid, group):
        nonlocal discarded_hsps, aligned_reads
        if not group:
            return
        aligned_reads += 1
        maximum = max(row['bitscore'] for row in group.values())
        best = [sid for sid, row in group.items() if math.isclose(row['bitscore'], maximum, abs_tol=1e-6)]
        for sid, row in group.items():
            reference = references[sid]
            obs, dels, ins = _aligned_bases(row, queries[qid], reference['sequence'], min_base_quality)
            if sid not in counts:
                length = len(reference['sequence'])
                counts[sid] = {'bases': [array('I', [0]) * length for _ in range(4)],
                               'deletions': array('I', [0]) * length, 'insertions': defaultdict(int),
                               'reads': 0, 'unique': 0, 'tied': 0}
            state = counts[sid]
            state['reads'] += 1
            state['unique'] += sid in best and len(best) == 1
            state['tied'] += sid in best and len(best) > 1
            for position, base in obs.items():
                state['bases']['ACGT'.index(base)][position] += 1
            for position in dels:
                state['deletions'][position] += 1
            for position in ins:
                state['insertions'][position] += 1

    current, group = None, {}
    with output.open(encoding='ascii') as handle:
        for line in handle:
            check_cancelled(cancelled)
            values = line.rstrip('\n').split('\t')
            if len(values) != len(fields):
                raise ValueError('BLAST returned a malformed read-support row.')
            row = dict(zip(fields, values, strict=True))
            qid, sid = row['qseqid'], row['sseqid']
            if qid not in queries or sid not in references:
                raise ValueError('BLAST output refers to an unknown read or allele.')
            if qid != current:
                consume(current, group)
                if qid in seen_queries:
                    raise ValueError('Unexpected noncontiguous BLAST query output; no duplicated depth was accepted.')
                seen_queries.add(qid)
                current, group = qid, {}
            row['pident'], row['bitscore'] = float(row['pident']), float(row['bitscore'])
            if not math.isfinite(row['pident']) or not math.isfinite(row['bitscore']) or not 0 <= row['pident'] <= 100:
                raise ValueError('Invalid BLAST read-support identity or score.')
            coverage = (int(row['qend']) - int(row['qstart']) + 1) / len(queries[qid].sequence)
            if row['pident'] < min_identity or coverage < min_read_coverage:
                low_filter_hsps += 1
                continue
            if sid in group:
                discarded_hsps += 1
            if sid not in group or row['bitscore'] > group[sid]['bitscore']:
                group[sid] = row
        consume(current, group)
    results = [_summarize_candidate(entry, counts.get(sid), min_depth=min_depth, min_breadth=min_breadth)
               for sid, entry in references.items()]
    return results, {'reads_with_qualifying_alignment': aligned_reads, 'filtered_alignment_rows': low_filter_hsps,
                     'additional_hsps_not_counted_twice': discarded_hsps,
                     'per_read_per_allele_policy': 'One highest-scoring local HSP; additional HSPs are disclosed, not double-counted.'}


def _summarize_loci(candidates, loci):
    results = []
    for locus in loci:
        rows = [row for row in candidates if row['locus'] == locus]
        supported = [row for row in rows if row['status'] == 'supported']
        states = {row['status'] for row in rows}
        cross_locus = (bool(supported) and not any(row['uniquely_best_reads'] for row in supported)
                       and any(row['locus'] != locus and row['status'] == 'supported' for row in candidates))
        # Wrong/divergent reference alternatives can generate local mismatch
        # mixtures. They must not override a fully consistent candidate to label
        # the isolate mixed; their individual review evidence remains available.
        state = ('ambiguous' if len(supported) > 1 or cross_locus else 'supported' if supported
                 else 'mixed_support' if 'mixed_support' in states else
                 'incomplete_or_discordant' if states != {'no_support'} else 'no_support')
        results.append({'locus': locus, 'status': state, 'candidates': rows,
                        'compatible_candidates': [row['allele'] for row in supported], 'assigned_allele': None,
                        'interpretation': 'Support is restricted to tested references and aligned reads; no assembly allele has been assigned.'})
    return results


def investigate_read_support(read1, read2, scheme, loci, cancelled=None, progress=None, *,
                             assembly_path=None, expected_read_sha256=None, max_pairs=100_000,
                             max_loci=MAX_PANEL_LOCI, max_alleles=5000, max_reference_bases=10_000_000,
                             max_query_bases=60_000_000, max_output_bytes=128 * 1024 * 1024,
                             min_read_length=50, max_read_length=1000, min_base_quality=20,
                             min_identity=90.0, min_read_coverage=.8, min_depth=3, min_breadth=.95,
                             threads=2, blastn_path=None, makeblastdb_path=None):
    """Run a local, bounded, both-strands short-read support assay.

    Input hashes cover the complete compressed input files, even when alignment
    uses only the disclosed prefix. The unread tail is not structurally validated.
    Bounds cause explicit failure, never incomplete output interpreted as absence.
    """
    check_cancelled(cancelled)
    panel = preflight_read_support(scheme, loci, max_loci=max_loci, max_alleles=max_alleles,
                                   max_reference_bases=max_reference_bases)
    for value, name, low, high in ((max_pairs, 'Maximum pairs', 1, 1_000_000),
                                  (max_query_bases, 'Maximum query bases', 100, 100_000_000),
                                  (max_output_bytes, 'Maximum output bytes', 1, 512 * 1024 * 1024),
                                  (min_read_length, 'Minimum read length', 30, 1000),
                                  (max_read_length, 'Maximum read length', min_read_length, 10_000),
                                  (min_base_quality, 'Minimum Phred+33 base quality', 0, 60),
                                  (min_depth, 'Minimum read depth', 1, 1000), (threads, 'Threads', 1, 256)):
        _integer(value, name, low, high)
    if not (math.isfinite(min_identity) and 80 <= min_identity <= 100 and
            math.isfinite(min_read_coverage) and .5 <= min_read_coverage <= 1 and
            math.isfinite(min_breadth) and .5 <= min_breadth <= 1):
        raise ValueError('Invalid read-support identity or breadth thresholds.')
    paths = [Path(read1).resolve(), Path(read2).resolve()]
    if paths[0].samefile(paths[1]):
        raise ValueError('Choose two different FASTQ mate files.')
    signatures = [file_signature(path) for path in paths]
    hashes = [file_sha256(path, cancelled) for path in paths]
    if expected_read_sha256 is not None and list(expected_read_sha256) != hashes:
        raise ValueError('Attached read SHA-256 identities no longer match the actual FASTQ files.')
    assembly = Path(assembly_path).resolve() if assembly_path is not None else None
    assembly_signature, assembly_sha = None, None
    if assembly is not None:
        assembly_signature, assembly_sha = file_signature(assembly), file_sha256(assembly, cancelled)
        with SequenceReader(assembly, cancelled) as reader:
            if reader.kind != 'fasta':
                raise ValueError('The linked assembly must be FASTA.')
            for _ in reader:
                pass
    blastn, makeblastdb = resolve_blast('blastn', blastn_path), resolve_blast('makeblastdb', makeblastdb_path)
    references, commands = {}, []
    with tempfile.TemporaryDirectory(prefix='wmlstudio-read-support-') as temporary:
        work = Path(temporary)
        query, reference, output = work / 'reads.fasta', work / 'alleles.fasta', work / 'hits.tsv'
        queries, sampling = _read_queries(paths, query, max_pairs=max_pairs, max_query_bases=max_query_bases,
                                           min_read_length=min_read_length, max_read_length=max_read_length,
                                           cancelled=cancelled, progress=progress)
        with reference.open('w', encoding='ascii') as out:
            for locus in panel['loci']:
                for allele, sequence in sorted(scheme.alleles[locus].items()):
                    sid = f'a{len(references)}'
                    references[sid] = {'locus': locus, 'allele': allele, 'sequence': sequence}
                    out.write(f'>{sid}\n{sequence}\n')
        command = [makeblastdb, '-in', str(reference), '-dbtype', 'nucl', '-out', str(work / 'alleles-db')]
        commands.append(command)
        _run(command, cancelled=cancelled)
        command = [blastn, '-task', 'blastn', '-word_size', '11', '-query', str(query),
                   '-db', str(work / 'alleles-db'), '-out', str(output), '-dust', 'no', '-soft_masking', 'false',
                   '-evalue', '1e-10', '-perc_identity', str(min_identity), '-qcov_hsp_perc', str(100 * min_read_coverage),
                   '-max_target_seqs', str(len(references)), '-num_threads', str(threads),
                   '-outfmt', '6 qseqid sseqid pident qstart qend sstart send bitscore qseq sseq']
        commands.append(command)

        def guard():
            if output.exists() and output.stat().st_size > max_output_bytes:
                raise ValueError('BLAST read-support output exceeds the explicit safety bound; no partial evidence was accepted. Select fewer loci or read pairs.')
            return cancelled is not None and cancelled()

        if progress:
            progress(0, 0, f'Aligning {sampling["reads_aligned"]:,} reads against all {len(references):,} selected alleles')
        _run(command, cancelled=guard)
        guard()
        candidates, alignments = _summarize_output(output, queries, references, min_identity=min_identity,
            min_read_coverage=min_read_coverage, min_base_quality=min_base_quality, min_depth=min_depth,
            min_breadth=min_breadth, cancelled=cancelled)
    if signatures != [file_signature(path) for path in paths]:
        raise ValueError('A read input changed during investigation; no result was accepted.')
    if assembly is not None and file_signature(assembly) != assembly_signature:
        raise ValueError('The assembly changed during read investigation; no result was accepted.')
    loci_results = _summarize_loci(candidates, panel['loci'])
    return {'format_version': 1, 'status': 'completed', 'scheme': scheme.name, 'scheme_digest': scheme.digest,
            'input_path': str(assembly) if assembly else None, 'input_sha256': assembly_sha,
            'assembly_signature': list(assembly_signature) if assembly_signature else None,
            'reads': [{'path': str(path), 'sha256': digest, 'signature': list(signature), 'mate': index}
                      for index, (path, digest, signature) in enumerate(zip(paths, hashes, signatures, strict=True), 1)],
            'panel': panel, 'sampling': sampling, 'alignments': alignments, 'loci': loci_results,
            'parameters': {'min_identity_pct': min_identity, 'min_read_coverage': min_read_coverage,
                           'min_base_quality_phred33': min_base_quality, 'min_depth': min_depth,
                           'min_breadth': min_breadth, 'max_query_bases': max_query_bases,
                           'max_output_bytes': max_output_bytes, 'min_read_length': min_read_length, 'max_read_length': max_read_length},
            'provenance': {'software': 'WMLSTudio', 'version': __version__, 'created_utc': datetime.now(UTC).isoformat(),
                           'engine': 'NCBI BLAST+', 'engine_version': _run([blastn, '-version'], cancelled=cancelled).splitlines()[0],
                           'binary_sha256': file_sha256(blastn, cancelled), 'commands': commands,
                           'manual': 'https://www.ncbi.nlm.nih.gov/books/NBK279684/'}, 'limitations': list(LIMITATIONS)}


def current_read_support(record):
    evidence = (record.get('metadata') or {}).get('read_support') or {}
    attached = (record.get('metadata') or {}).get('reads') or {}
    current = current_input_sha256(record)
    if not evidence:
        return {'status': 'missing', 'reason': 'No read-support investigation is attached.', 'evidence': None}
    actual = [_sha256(row.get('sha256')) for row in attached.get('reads', [])]
    tested = [_sha256(row.get('sha256')) for row in evidence.get('reads', [])]
    if not current or not _sha256(evidence.get('input_sha256')) or len(actual) != 2 or None in actual:
        status, reason = 'unverified', 'Recorded assembly or attached-read identity is unavailable.'
    elif current != evidence['input_sha256'] or actual != tested:
        status, reason = 'stale', 'The assembly or attached reads differ from this investigation.'
    else:
        status, reason = 'current', 'Recorded assembly and both attached-read SHA-256 identities match.'
    return {'status': status, 'reason': reason, 'evidence': evidence if status == 'current' else None}


def persist_read_support(project, sample_id, result):
    """Quick serial commit using worker-computed full hashes and stable signatures.

    No large FASTQ hashing belongs on the GUI thread here. The investigator checks
    all hashes before work, signatures after work, and this function checks those
    same signatures again against the explicit attached-read identity.
    """
    result = copy.deepcopy(result)
    if result.get('format_version') != 1 or not _sha256(result.get('input_sha256')) or not result.get('input_path'):
        raise ValueError('Persistent read support requires a verified linked assembly identity.')
    if len(result.get('reads', [])) != 2 or any(not _sha256(read.get('sha256')) for read in result['reads']):
        raise ValueError('Persistent read support requires both full FASTQ SHA-256 identities.')
    with project.transaction():
        record = project.get_sample(sample_id)
        if Path(record['input_path']).resolve() != Path(result['input_path']).resolve():
            raise ValueError('The sample assembly path changed before read support could be attached.')
        current = current_input_sha256(record)
        if current and current != result['input_sha256']:
            raise ValueError('Read support belongs to a different current assembly SHA-256.')
        metadata = record.get('metadata') or {}
        attached = (metadata.get('reads') or {}).get('reads', [])
        if (len(attached) != 2 or [read.get('sha256') for read in attached] != [read['sha256'] for read in result['reads']]
                or [str(Path(read['path']).resolve()) for read in attached] != [read['path'] for read in result['reads']]):
            raise ValueError('Read support does not match the explicitly attached FASTQ pair.')
        checks = [{'path': result['input_path'], 'signature': result['assembly_signature']}, *result['reads']]
        if any(list(file_signature(row['path'])) != row['signature'] for row in checks):
            raise ValueError('An input changed before read support could be saved. Run the investigation again.')
        if metadata.get('read_support'):
            project.record_history(sample_id, 'read_support_superseded', {'evidence': metadata['read_support']})
        metadata['read_support'] = result
        if not current:
            metadata['input_identity'] = {'sha256': result['input_sha256'], 'path': result['input_path'],
                                          'source': 'full-input hash computed and signature-verified by read-support investigation'}
        project.set_metadata(sample_id, metadata)
        project.record_history(sample_id, 'read_support_attached', {'input_sha256': result['input_sha256'],
                               'scheme_digest': result['scheme_digest'], 'loci': result['panel']['loci']})
    return project.get_sample(sample_id)
