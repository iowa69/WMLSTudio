"""Build an immutable local, reference-anchored cgMLST scheme from assemblies.

This is cohort-defined ortholog selection, not a validated public nomenclature.
All retained sequences are complete CDSs. Copies and cross-locus homologs at the
declared nucleotide thresholds exclude a seed locus from the entire scheme.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pyrodigal

from .cgtyping import CGTypingError, _run, predict_cds, resolve_blast
from .sequence import check_cancelled, file_sha256
from .typing import load_scheme


def _relations(queries, genes, work, blastn, makeblastdb, *, min_identity,
               min_coverage, threads, cancelled):
    """Return every full-coverage homolog; never choose a single best copy."""
    if not queries or not genes:
        return {}
    query_path, subject_path = work / 'seed.fasta', work / 'genes.fasta'
    with query_path.open('w', encoding='ascii') as handle:
        for gene in queries:
            handle.write(f">{gene['id']}\n{gene['sequence']}\n")
    with subject_path.open('w', encoding='ascii') as handle:
        for gene in genes:
            handle.write(f">{gene['id']}\n{gene['sequence']}\n")
    database, output = work / 'genesdb', work / 'hits.tsv'
    _run([makeblastdb, '-in', str(subject_path), '-dbtype', 'nucl',
          '-out', str(database), '-blastdb_version', '4'], cancelled=cancelled)
    _run([blastn, '-task', 'megablast', '-query', str(query_path), '-db', str(database),
          '-out', str(output), '-outfmt', '6 qseqid sseqid pident qstart qend sstart send qlen slen',
          '-num_threads', str(threads), '-perc_identity', str(min_identity * 100),
          '-qcov_hsp_perc', str(min_coverage * 100), '-max_target_seqs', str(len(genes)),
          '-max_hsps', '1', '-dust', 'no', '-evalue', '1e-10'], cancelled=cancelled)
    query_ids, gene_ids = {g['id'] for g in queries}, {g['id'] for g in genes}
    relations = defaultdict(set)
    with output.open(encoding='ascii') as handle:
        for index, line in enumerate(handle):
            if index % 4096 == 0:
                check_cancelled(cancelled)
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[0] not in query_ids or fields[1] not in gene_ids:
                raise CGTypingError('Unexpected BLAST record during ad-hoc scheme construction.')
            qcoverage = (abs(int(fields[4]) - int(fields[3])) + 1) / int(fields[7])
            scoverage = (abs(int(fields[6]) - int(fields[5])) + 1) / int(fields[8])
            if float(fields[2]) / 100 >= min_identity and min(qcoverage, scoverage) >= min_coverage:
                relations[fields[0]].add(fields[1])
    return dict(relations)


def _select_cohort_loci(seeds, cohort_genes, cohort_relations, min_prevalence):
    """Pure scientific selection kernel, used by real native execution and tests."""
    required = math.ceil(len(cohort_genes) * min_prevalence)
    retained, excluded = {}, {}
    gene_maps = [{gene['id']: gene for gene in genes} for genes in cohort_genes]
    shared_genes = []
    for relations in cohort_relations:
        inverse = defaultdict(set)
        for seed, identifiers in relations.items():
            for identifier in identifiers:
                inverse[identifier].add(seed)
        shared_genes.append(inverse)
    for seed in seeds:
        variants, present, reason = {}, [], None
        for index, (by_id, relations) in enumerate(zip(gene_maps, cohort_relations, strict=True)):
            matches = relations.get(seed['id'], set())
            if len(matches) > 1 or any(len(shared_genes[index][identifier]) > 1 for identifier in matches):
                reason = 'multiple CDS copies or cross-locus homology at the configured thresholds'
                break
            if len(matches) == 1:
                candidate = by_id[next(iter(matches))]
                if candidate['cds_qc']['valid']:
                    variants[candidate['sequence_sha256']] = candidate['sequence']
                    present.append(index)
        if reason is None and len(present) < required:
            reason = f'complete unique CDS present in {len(present)}/{len(cohort_genes)} assemblies; {required} required'
        if reason:
            excluded[seed['id']] = reason
        else:
            retained[seed['id']] = {'seed': seed, 'variants': variants, 'present_in': present}
    return retained, excluded


def create_adhoc_scheme(assemblies, library_root, name, *, organism='',
                        min_prevalence=1.0, min_identity=0.95, min_coverage=0.98,
                        genetic_code=11, threads=2, cancelled=None, progress=None,
                        blastn_path=None, makeblastdb_path=None):
    """Build a versioned local scheme; first input is the reference anchor.

    Returns ``path, scheme_digest, locus_count, excluded_loci, sources, created``.
    Original assemblies are read only. Publication is an atomic directory move.
    """
    if not str(name).strip():
        raise ValueError('An ad-hoc scheme needs a descriptive name.')
    if any(not 0 < value <= 1 for value in (min_prevalence, min_identity, min_coverage)):
        raise ValueError('Prevalence, identity and coverage must be fractions in (0, 1].')
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError('The BLAST thread count must be a positive integer.')
    paths = [Path(path).resolve() for path in assemblies]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError('Supply at least one distinct assembly; duplicate paths are not allowed.')
    check_cancelled(cancelled)
    blastn, makeblastdb = resolve_blast('blastn', blastn_path), resolve_blast('makeblastdb', makeblastdb_path)
    versions = {tool: _run([binary, '-version'], cancelled=cancelled).splitlines()[0]
                for tool, binary in [('blastn', blastn), ('makeblastdb', makeblastdb)]}
    parameters = {'method': 'reference-anchored-full-cds-v1', 'name': str(name).strip(),
                  'organism': str(organism).strip(), 'min_prevalence': min_prevalence,
                  'min_identity': min_identity, 'min_coverage': min_coverage,
                  'genetic_code': genetic_code, 'threads': threads,
                  'pyrodigal': pyrodigal.__version__, **versions,
                  'blastn_sha256': file_sha256(blastn, cancelled),
                  'makeblastdb_sha256': file_sha256(makeblastdb, cancelled)}
    sources, cohort_genes = [], []
    for index, path in enumerate(paths):
        check_cancelled(cancelled)
        if progress:
            progress(index, len(paths), f'Predicting complete CDS: {path.name}')
        digest = file_sha256(path, cancelled)
        genes, prediction = predict_cds(path, genetic_code=genetic_code, cancelled=cancelled)
        if file_sha256(path, cancelled) != digest:
            raise CGTypingError(f'Assembly changed during gene prediction: {path.name}')
        cohort_genes.append(genes)
        sources.append({'path': str(path), 'sha256': digest, 'gene_prediction': prediction})
    if len({source['sha256'] for source in sources}) != len(sources):
        raise ValueError('Duplicate assembly contents would bias locus prevalence; supply unique assemblies.')
    seeds = [gene for gene in cohort_genes[0] if gene['cds_qc']['valid']]
    if not seeds:
        raise CGTypingError('The reference assembly contains no validated complete CDS predictions.')
    root = Path(library_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.adhoc-build-', dir=root) as temporary:
        temporary = Path(temporary)
        relations = []
        for index, genes in enumerate(cohort_genes):
            check_cancelled(cancelled)
            if progress:
                progress(index, len(paths), f'Screening orthologs and paralogs: {paths[index].name}')
            work = temporary / f'assembly-{index}'
            work.mkdir()
            # Partial predictions remain targets so a second homologous copy
            # cannot disappear merely because it failed the full-CDS gate.
            relations.append(_relations(seeds, genes, work, blastn, makeblastdb,
                                         min_identity=min_identity, min_coverage=min_coverage,
                                         threads=threads, cancelled=cancelled))
        retained, excluded = _select_cohort_loci(seeds, cohort_genes, relations, min_prevalence)
        if not retained:
            raise CGTypingError('No unique complete-CDS loci passed the declared cohort thresholds.')
        stage = temporary / 'snapshot'
        stage.mkdir()
        for index, data in enumerate(retained.values()):
            check_cancelled(cancelled)
            locus = 'CDS_' + data['seed']['sequence_sha256']
            with (stage / f'{locus}.fasta').open('w', encoding='ascii') as handle:
                for digest, sequence in sorted(data['variants'].items()):
                    handle.write(f'>{locus}_SHA256_{digest}\n{sequence}\n')
            if progress and index % 100 == 0:
                progress(index, len(retained), 'Writing validated local allele nomenclature')
        metadata = {'name': str(name).strip(), 'organism': str(organism).strip(),
                    'type': 'cgMLST', 'provider': 'local-ad-hoc', 'parameters': parameters,
                    'source_sha256': [source['sha256'] for source in sources],
                    'access_notice': 'Local cohort-defined reference-anchored scheme; not a validated public cgMLST nomenclature. No registered STs or transmission threshold.'}
        (stage / 'scheme.json').write_text(json.dumps(metadata, sort_keys=True, indent=2), encoding='utf-8')
        scheme = load_scheme(stage, cancelled=cancelled, sequences=False)
        manifest = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
                    'scheme_digest': scheme.digest, 'parameters': parameters, 'sources': sources,
                    'reference_anchor': sources[0], 'locus_count': scheme.locus_count,
                    'excluded_loci': excluded,
                    'loci': {data['seed']['sequence_sha256']: {'present_in': data['present_in'],
                             'variants': sorted(data['variants'])} for data in retained.values()}}
        (stage / 'adhoc_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        for path, source in zip(paths, sources, strict=True):
            if file_sha256(path, cancelled) != source['sha256']:
                raise CGTypingError(f'Assembly changed while constructing the scheme: {path.name}')
        destination = root / f'adhoc_{scheme.digest[:24]}'
        created = not destination.exists()
        if created:
            check_cancelled(cancelled)
            os.rename(stage, destination)
        elif destination.is_symlink() or load_scheme(destination, cancelled=cancelled, sequences=False).digest != scheme.digest:
            raise CGTypingError('Existing local snapshot was modified; it will not be overwritten.')
        return {'path': str(destination), 'scheme_digest': scheme.digest,
                'locus_count': scheme.locus_count, 'excluded_loci': excluded,
                'sources': sources, 'created': created}
