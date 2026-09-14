"""Local KpSC-associated virulence gene screening with explicit assay coverage.

Uses curated Kleborate/Pasteur allele sequences, but is a WMLSTudio BLAST+
screen, not a reimplementation or claimed equivalent of Kleborate's typing.
No virulence phenotype, official locus ST or plasmid identity is inferred.
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


def summarize_virulence_hits(hits, loci, *, adequate_negative_assay):
    by_gene = defaultdict(list)
    for hit in hits:
        by_gene[hit['gene']].append(hit)
    groups = []
    for locus, entries in loci.items():
        genes = []
        for entry in entries:
            gene_hits = by_gene[entry['gene']]
            accepted = [hit for hit in gene_hits if hit['identity_pct'] >= 90 and hit['coverage_pct'] >= 80]
            intact = [hit for hit in accepted if hit['cds_qc']['valid']]
            status = 'detected' if accepted else ('ambiguous' if gene_hits or not adequate_negative_assay else 'not_detected')
            genes.append({'gene': entry['gene'], 'status': status, 'intact_cds_copies': len(intact),
                          'completeness': 'intact_cds_detected' if intact else ('partial_or_disrupted' if accepted else 'unresolved'),
                          'hits': gene_hits})
        detected = sum(gene['status'] == 'detected' for gene in genes)
        intact = sum(gene['intact_cds_copies'] > 0 for gene in genes)
        state = 'detected' if detected else ('not_detected' if all(gene['status'] == 'not_detected' for gene in genes) else 'ambiguous')
        groups.append({'locus': locus, 'status': state, 'genes_detected': detected, 'genes_total': len(genes),
                       'intact_genes_detected': intact, 'completeness': 'all_genes_intact_cds' if intact == len(genes) else 'incomplete_or_unresolved',
                       'genes': genes, 'official_locus_st': None, 'phenotype': 'not_inferred'})
    return groups


def _unique_locations(hits):
    """Collapse allele alternatives at one physical locus, not separate copies."""
    selected = []
    for hit in sorted(hits, key=lambda row: (-row['bitscore'], -row['identity_pct'], -row['coverage_pct'], row['reference_allele'])):
        duplicate = False
        for old in selected:
            if (old['gene'], old['contig'], old['strand']) != (hit['gene'], hit['contig'], hit['strand']):
                continue
            overlap = max(0, min(old['end'], hit['end']) - max(old['start'], hit['start']) + 1)
            if overlap / min(old['end'] - old['start'] + 1, hit['end'] - hit['start'] + 1) >= .8:
                old.setdefault('alternative_reference_alleles', []).append(hit['reference_allele'])
                duplicate = True
                break
        if not duplicate:
            selected.append(hit)
    return sorted(selected, key=lambda row: (row['gene'], row['contig'], row['start'], row['end']))


def screen_virulence(path, reference_root, cancelled=None, progress=None, *, threads=2, blastn_path=None, makeblastdb_path=None):
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ValueError('Virulence threads must be an integer between 1 and 256.')
    path, root = Path(path).resolve(), Path(reference_root).resolve()
    check_cancelled(cancelled)
    manifest = validate_characterization_references(root, cancelled=cancelled)
    reference_signatures = [(entry['path'], file_signature(root / entry['path'])) for entry in manifest['files']]
    reference_signatures.append(('manifest.json', file_signature(root / 'manifest.json')))
    blastn, makeblastdb = resolve_blast('blastn', blastn_path), resolve_blast('makeblastdb', makeblastdb_path)
    signature = file_signature(path)
    qc = inspect_sequence(path, cancelled=cancelled)
    if qc['kind'] != 'fasta':
        raise ValueError('Virulence characterization requires a FASTA assembly.')
    input_hash = qc['input_sha256']
    contigs, identifiers = {}, set()
    with SequenceReader(path, cancelled) as reader:
        for index, record in enumerate(reader):
            if record.identifier in identifiers:
                raise ValueError(f'Duplicate contig identifier: {record.identifier}')
            identifiers.add(record.identifier)
            contigs[f'c{index}'] = record
    queries, commands = {}, []
    with tempfile.TemporaryDirectory(prefix='wmlstudio-virulence-') as temporary:
        work = Path(temporary)
        assembly = work / 'assembly.fasta'
        with assembly.open('w', encoding='ascii') as out:
            for identifier, record in contigs.items():
                out.write(f'>{identifier}\n{record.sequence}\n')
        query = work / 'virulence-alleles.fasta'
        with query.open('w', encoding='ascii') as out:
            for locus, entries in manifest['virulence'].items():
                for entry in entries:
                    check_cancelled(cancelled)
                    with SequenceReader(root / entry['path'], cancelled) as reader:
                        for record in reader:
                            identifier = f'q{len(queries)}'
                            queries[identifier] = {'locus': locus, 'gene': entry['gene'], 'reference_allele': record.identifier,
                                                   'reference_length': len(record.sequence)}
                            out.write(f'>{identifier}\n{record.sequence}\n')
        if not queries:
            raise ValueError('No virulence allele references are installed.')
        command = [makeblastdb, '-in', str(assembly), '-dbtype', 'nucl', '-out', str(work / 'assembly-db')]
        commands.append(command)
        _run(command, cancelled=cancelled)
        output = work / 'hits.tsv'
        command = [blastn, '-query', str(query), '-db', str(work / 'assembly-db'), '-out', str(output),
                   '-task', 'blastn', '-dust', 'no', '-evalue', '1e-10', '-perc_identity', '80',
                   '-qcov_hsp_perc', '40', '-max_target_seqs', str(max(1, len(contigs))), '-num_threads', str(threads),
                   '-outfmt', '6 qseqid sseqid pident qstart qend sstart send bitscore']
        commands.append(command)
        if progress:
            progress(0, 0, 'Screening six defined KpSC-associated virulence loci with native BLAST+…')
        _run(command, cancelled=cancelled)
        if output.stat().st_size > 256 * 1024 * 1024:
            raise ValueError('Virulence evidence exceeds the 256 MiB safety bound; no partial report was accepted.')
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
                    raise ValueError('BLAST returned invalid virulence evidence coordinates or quantities.')
                strand = '+' if sstart <= send else '-'
                sequence = contig.sequence[start - 1:end]
                if strand == '-':
                    sequence = reverse_complement(sequence)
                hits.append(dict(entry, contig=contig.identifier, start=start, end=end, strand=strand,
                                 identity_pct=identity, coverage_pct=100 * (qend - qstart + 1) / entry['reference_length'],
                                 bitscore=score, cds_qc=full_cds_qc(sequence, partial_begin=qstart != 1,
                                                               partial_end=qend != entry['reference_length'])))
    hits = _unique_locations(hits)
    if any(file_signature(root / relative) != signature for relative, signature in reference_signatures):
        raise ValueError('Virulence reference files changed during screening; no result was accepted.')
    if file_signature(path) != signature or file_sha256(path, cancelled) != input_hash:
        raise ValueError('Assembly changed during virulence screening; no result was accepted.')
    metrics = qc['qc']
    adequate = metrics['total_bases'] >= 1_000_000 and (metrics.get('n50') or 0) >= 10_000 and metrics['ambiguous_percent'] <= 1
    return {'status': 'completed', 'input_sha256': input_hash, 'reference_digest': manifest['reference_digest'],
            'loci': summarize_virulence_hits(hits, manifest['virulence'], adequate_negative_assay=adequate), 'hits': hits,
            'assay': {'name': 'WMLSTudio KpSC-associated virulence allele screen', 'gene_count': sum(map(len, manifest['virulence'].values())),
                      'loci': list(manifest['virulence']), 'min_identity_pct': 90, 'min_coverage_pct': 80,
                      'review_identity_pct': 80, 'review_coverage_pct': 40, 'adequate_negative_assay': adequate,
                      'negative_qc_gates': {'min_total_bases': 1_000_000, 'min_n50': 10_000, 'max_ambiguous_percent': 1}},
            'provenance': {'engine': 'NCBI BLAST+', 'version': _run([blastn, '-version'], cancelled=cancelled).splitlines()[0],
                           'binary_sha256': file_sha256(blastn, cancelled), 'commands': commands,
                           'reference_source': manifest['source_repository'], 'reference_revision': manifest['source_revision']},
            'limitations': ['Not a Kleborate-equivalent result; no virulence-locus ST or lineage has been assigned.',
                            'Intact coding sequence is not proof of expression or virulence; partial/disrupted matches require review.',
                            'Not detected means no hit meeting this defined assay threshold, not proof of genomic absence.',
                            'Negative findings are withheld when assembly quality fails explicit minimum screening gates.',
                            'These loci may occur outside KpSC. No hypervirulence or plasmid carriage phenotype is inferred.']}
