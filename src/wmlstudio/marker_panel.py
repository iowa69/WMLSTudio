"""One shared native BLAST+ execution path for every pinned marker reference panel.

Reference records are the BLAST query and the assembly is the subject database,
so every hit carries both reference and contig coordinates. Nothing here
interprets a panel: each organism module supplies its own decision kernel,
thresholds and limitations, and no phenotype is inferred at this layer.
"""

from __future__ import annotations

import math
import tempfile
from collections import defaultdict
from pathlib import Path

from .cgtyping import _run, full_cds_qc, resolve_blast
from .characterization_refs import validate_characterization_references
from .sequence import SequenceReader, check_cancelled, file_sha256, file_signature, inspect_sequence
from .typing import reverse_complement


def panel_entries(manifest, section):
    """Resolve a dotted manifest section to {group: [entry]} without inventing groups."""
    node, walked = manifest, []
    for part in section.split('.'):
        walked.append(part)
        if not isinstance(node, dict) or part not in node:
            raise ValueError(f'This reference snapshot carries no {section} panel section.')
        node = node[part]
    if isinstance(node, list):
        node = {walked[-1]: node}
    if not isinstance(node, dict) or not node or not all(isinstance(value, list) for value in node.values()):
        raise ValueError(f'Reference panel section {section} is not a group-to-entry mapping.')
    return node


def unique_locations(hits, member_key='gene'):
    """Collapse allele alternatives at one physical locus, not separate copies."""
    selected = []
    for hit in sorted(hits, key=lambda row: (-row['bitscore'], -row['identity_pct'], -row['coverage_pct'], row['reference_allele'])):
        duplicate = False
        for old in selected:
            if (old[member_key], old['contig'], old['strand']) != (hit[member_key], hit['contig'], hit['strand']):
                continue
            overlap = max(0, min(old['end'], hit['end']) - max(old['start'], hit['start']) + 1)
            if overlap / min(old['end'] - old['start'] + 1, hit['end'] - hit['start'] + 1) >= .8:
                old.setdefault('alternative_reference_alleles', []).append(hit['reference_allele'])
                duplicate = True
                break
        if not duplicate:
            selected.append(hit)
    return sorted(selected, key=lambda row: (row[member_key], row['contig'], row['start'], row['end']))


def union_coverage(hits, reference_length):
    """Interval union of alignments in REFERENCE coordinates, summed across contigs.

    A long cassette reference never produces one alignment covering most of itself
    in a draft assembly, so per-alignment coverage is the wrong measure. Coverage
    summed across contigs is consistent with, and never proof of, one intact
    element: ``fragmented`` records that the retained alignments came from more
    than one contig. ``identity_pct`` is the mean percent identity of the retained
    alignments weighted by the novel reference bases each one contributed.
    """
    if isinstance(reference_length, bool) or not isinstance(reference_length, int) or reference_length < 1:
        raise ValueError('Union coverage requires a positive reference length.')
    intervals, contigs, weighted = [], [], 0.0
    for hit in sorted(hits, key=lambda row: (-row['bitscore'], row['reference_start'], row['reference_end'])):
        start, end = hit['reference_start'], hit['reference_end']
        if isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('Union coverage requires integer reference coordinates.')
        if not 1 <= start <= end <= reference_length:
            raise ValueError('Union coverage received reference coordinates outside the reference.')
        # ``intervals`` stays merged and disjoint, so overlapping fragments of the
        # same reference region are subtracted exactly once.
        novel = (end - start + 1) - sum(max(0, min(end, old_end) - max(start, old_start) + 1)
                                        for old_start, old_end in intervals)
        if novel <= 0:
            continue
        intervals = _merged_intervals([*intervals, (start, end)])
        weighted += novel * float(hit['identity_pct'])
        if hit['contig'] not in contigs:
            contigs.append(hit['contig'])
    covered = sum(end - start + 1 for start, end in intervals)
    return {'aligned_bp': covered, 'coverage_pct': 100 * covered / reference_length,
            'identity_pct': (weighted / covered) if covered else 0.0,
            'contigs': sorted(contigs), 'contig_count': len(contigs), 'fragmented': len(contigs) > 1}


def _merged_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def summarize_marker_hits(hits, groups, *, adequate_negative_assay, member_key='gene', group_key='locus',
                          min_identity=90, min_coverage=80):
    """Per-member calls with negatives withheld whenever the assay QC gate failed."""
    by_member = defaultdict(list)
    for hit in hits:
        by_member[hit[member_key]].append(hit)
    summary = []
    for group, entries in groups.items():
        members = []
        for entry in entries:
            member_hits = by_member[entry[member_key]]
            accepted = [hit for hit in member_hits if hit['identity_pct'] >= min_identity and hit['coverage_pct'] >= min_coverage]
            intact = [hit for hit in accepted if hit['cds_qc']['valid']]
            status = 'detected' if accepted else ('ambiguous' if member_hits or not adequate_negative_assay else 'not_detected')
            members.append({member_key: entry[member_key], 'status': status, 'intact_cds_copies': len(intact),
                            'completeness': 'intact_cds_detected' if intact else ('partial_or_disrupted' if accepted else 'unresolved'),
                            'hits': member_hits})
        detected = sum(member['status'] == 'detected' for member in members)
        intact = sum(member['intact_cds_copies'] > 0 for member in members)
        state = 'detected' if detected else ('not_detected' if all(member['status'] == 'not_detected' for member in members) else 'ambiguous')
        summary.append({group_key: group, 'status': state, 'genes_detected': detected, 'genes_total': len(members),
                        'intact_genes_detected': intact, 'completeness': 'all_genes_intact_cds' if intact == len(members) else 'incomplete_or_unresolved',
                        'genes': members, 'phenotype': 'not_inferred'})
    return summary


def screen_marker_panel(path, reference_root, cancelled=None, progress=None, *, section, assay_name,
                        progress_message, limitations, subject='Marker panel', member_key='gene', group_key='locus',
                        min_identity=90, min_coverage=80, review_identity=80, review_coverage=40,
                        blast_task='blastn', evalue='1e-10', qc_gates=(1_000_000, 10_000, 1.0),
                        summarize=None, temp_prefix='wmlstudio-marker-panel-',
                        threads=2, blastn_path=None, makeblastdb_path=None):
    """Screen one pinned reference panel against an assembly with native BLAST+.

    The reference snapshot is re-verified before and after the search and the
    assembly is re-hashed, so a result is never accepted across a changed input.
    ``min_identity``/``min_coverage`` are the calling floors applied by
    ``summarize``; ``review_identity``/``review_coverage`` are the looser BLAST
    reporting gates that keep sub-threshold evidence visible for review.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ValueError(f'{subject} threads must be an integer between 1 and 256.')
    summarize = summarize or summarize_marker_hits
    min_total_bases, min_n50, max_ambiguous = qc_gates
    path, root = Path(path).resolve(), Path(reference_root).resolve()
    check_cancelled(cancelled)
    manifest = validate_characterization_references(root, cancelled=cancelled)
    groups = panel_entries(manifest, section)
    reference_signatures = [(entry['path'], file_signature(root / entry['path'])) for entry in manifest['files']]
    reference_signatures.append(('manifest.json', file_signature(root / 'manifest.json')))
    blastn, makeblastdb = resolve_blast('blastn', blastn_path), resolve_blast('makeblastdb', makeblastdb_path)
    signature = file_signature(path)
    qc = inspect_sequence(path, cancelled=cancelled)
    if qc['kind'] != 'fasta':
        raise ValueError(f'{subject} characterization requires a FASTA assembly.')
    input_hash = qc['input_sha256']
    contigs, identifiers = {}, set()
    with SequenceReader(path, cancelled) as reader:
        for index, record in enumerate(reader):
            if record.identifier in identifiers:
                raise ValueError(f'Duplicate contig identifier: {record.identifier}')
            identifiers.add(record.identifier)
            contigs[f'c{index}'] = record
    queries, commands = {}, []
    with tempfile.TemporaryDirectory(prefix=temp_prefix) as temporary:
        work = Path(temporary)
        assembly = work / 'assembly.fasta'
        with assembly.open('w', encoding='ascii') as out:
            for identifier, record in contigs.items():
                out.write(f'>{identifier}\n{record.sequence}\n')
        query = work / 'panel-references.fasta'
        with query.open('w', encoding='ascii') as out:
            for group, entries in groups.items():
                for entry in entries:
                    check_cancelled(cancelled)
                    with SequenceReader(root / entry['path'], cancelled) as reader:
                        for record in reader:
                            identifier = f'q{len(queries)}'
                            queries[identifier] = {group_key: group, member_key: entry[member_key],
                                                   'reference_allele': record.identifier,
                                                   'reference_length': len(record.sequence)}
                            out.write(f'>{identifier}\n{record.sequence}\n')
        if not queries:
            raise ValueError(f'No {section} reference sequences are installed.')
        command = [makeblastdb, '-in', str(assembly), '-dbtype', 'nucl', '-out', str(work / 'assembly-db')]
        commands.append(command)
        _run(command, cancelled=cancelled)
        output = work / 'hits.tsv'
        command = [blastn, '-query', str(query), '-db', str(work / 'assembly-db'), '-out', str(output),
                   '-task', blast_task, '-dust', 'no', '-evalue', evalue, '-perc_identity', str(review_identity)]
        if review_coverage > 0:
            command += ['-qcov_hsp_perc', str(review_coverage)]
        command += ['-max_target_seqs', str(max(1, len(contigs))), '-num_threads', str(threads),
                    '-outfmt', '6 qseqid sseqid pident qstart qend sstart send bitscore']
        commands.append(command)
        if progress:
            progress(0, 0, progress_message)
        _run(command, cancelled=cancelled)
        if output.stat().st_size > 256 * 1024 * 1024:
            raise ValueError(f'{subject} evidence exceeds the 256 MiB safety bound; no partial report was accepted.')
        hits = []
        with output.open(encoding='ascii') as rows:
            for line in rows:
                check_cancelled(cancelled)
                qid, sid, identity, qstart, qend, sstart, send, score = line.rstrip().split('\t')
                entry, contig = queries[qid], contigs[sid]
                identity, score = float(identity), float(score)
                qstart, qend, sstart, send = map(int, (qstart, qend, sstart, send))
                start, end = sorted((sstart, send))
                if not (math.isfinite(identity) and math.isfinite(score) and 0 <= identity <= 100
                        and 1 <= qstart <= qend <= entry['reference_length'] and 1 <= start <= end <= len(contig.sequence)):
                    raise ValueError(f'BLAST returned invalid {section} evidence coordinates or quantities.')
                strand = '+' if sstart <= send else '-'
                sequence = contig.sequence[start - 1:end]
                if strand == '-':
                    sequence = reverse_complement(sequence)
                hits.append(dict(entry, contig=contig.identifier, start=start, end=end, strand=strand,
                                 reference_start=qstart, reference_end=qend,
                                 identity_pct=identity, coverage_pct=100 * (qend - qstart + 1) / entry['reference_length'],
                                 bitscore=score, cds_qc=full_cds_qc(sequence, partial_begin=qstart != 1,
                                                               partial_end=qend != entry['reference_length'])))
    hits = unique_locations(hits, member_key)
    if any(file_signature(root / relative) != signature for relative, signature in reference_signatures):
        raise ValueError(f'{subject} reference files changed during screening; no result was accepted.')
    if file_signature(path) != signature or file_sha256(path, cancelled) != input_hash:
        raise ValueError(f'Assembly changed during {subject.lower()} screening; no result was accepted.')
    metrics = qc['qc']
    adequate = (metrics['total_bases'] >= min_total_bases and (metrics.get('n50') or 0) >= min_n50
                and metrics['ambiguous_percent'] <= max_ambiguous)
    return {'status': 'completed', 'input_sha256': input_hash, 'reference_digest': manifest['reference_digest'],
            'loci': summarize(hits, groups, adequate_negative_assay=adequate), 'hits': hits,
            'assay': {'name': assay_name, 'gene_count': sum(map(len, groups.values())),
                      'loci': list(groups), 'section': section, 'min_identity_pct': min_identity,
                      'min_coverage_pct': min_coverage, 'review_identity_pct': review_identity,
                      'review_coverage_pct': review_coverage, 'blast_task': blast_task,
                      'adequate_negative_assay': adequate,
                      'negative_qc_gates': {'min_total_bases': min_total_bases, 'min_n50': min_n50,
                                            'max_ambiguous_percent': max_ambiguous}},
            'provenance': {'engine': 'NCBI BLAST+', 'version': _run([blastn, '-version'], cancelled=cancelled).splitlines()[0],
                           'binary_sha256': file_sha256(blastn, cancelled), 'commands': commands,
                           'reference_source': manifest.get('source_repository'),
                           'reference_revision': manifest.get('source_revision')},
            'limitations': list(limitations)}
