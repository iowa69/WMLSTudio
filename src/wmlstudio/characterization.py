"""Evidence synthesis for local outbreak investigations, not clinical prediction.

Chromosomal typing, independent taxonomic evidence and accessory determinants
remain separate. No analysis of antimicrobial susceptibility or transmission is
performed by these functions.
"""

from __future__ import annotations

import copy
from datetime import UTC, datetime
from pathlib import Path

from . import __version__
from .hydra import _validate_sample, load_hydra_report
from .organism_modules import module_tasks, registered_modules, stamp_applicability
from .plasmid_evidence import MOB_SUITE_GAP, contig_plasmid_evidence, marker_location
from .sample_workflow import _execution_input_sha256, _sha256, current_input_sha256
from .sequence import (
    AnalysisCancelled,
    SequenceReader,
    check_cancelled,
    file_sha256,
    file_signature,
    inspect_sequence,
)
from .species_evidence import identify_species
from .virulence_evidence import apply_locus_sts, screen_virulence

LIMITATIONS = [
    'Genomic similarity alone does not establish transmission, direction, source attribution or an outbreak.',
    'AMR determinant annotations are genomic associations, not susceptibility, MIC, treatment advice or validated AST.',
    'A replicon and resistance/virulence marker on one contig support co-location only, not a reconstructed or transmissible plasmid.',
    'Negative results are meaningful only for the explicitly completed assay and its reference coverage/QC gates.',
    'The plasmid view is a contig screen over replicon markers. It is not MOB-suite, MOB-typer or MOB-recon '
    'and is not equivalent to them: no relaxase or MPF type, no oriT, no plasmid reconstruction, no cluster '
    'code and no mobility prediction is produced.',
]

# Explicit NCBI reference labels only, not a class-to-member-drug expansion.
# Table reviewed 2026-09-12; unrecognized/future labels remain in raw subclass.
DRUG_LABEL_SOURCE = 'https://github.com/ncbi/amr/wiki/class-subclass'
EXPLICIT_DRUG_LABELS = frozenset({
    'AMIKACIN', 'GENTAMICIN', 'KANAMYCIN', 'PLAZOMICIN', 'TOBRAMYCIN', 'APRAMYCIN',
    'SPECTINOMYCIN', 'STREPTOMYCIN', 'NEOMYCIN', 'PAROMOMYCIN',
    'AMOXICILLIN-CLAVULANIC ACID', 'PIPERACILLIN-TAZOBACTAM', 'TICARCILLIN-CLAVULANIC ACID',
    'AZTREONAM', 'AZTREONAM-AVIBACTAM', 'CEFEPIME', 'CEFIDEROCOL', 'CEFTAZIDIME-AVIBACTAM',
    'CEFTIBUTEN', 'CEFTIBUTEN-AVIBACTAM', 'CEFTAROLINE', 'CEPHALOTHIN', 'MEROPENEM', 'IMIPENEM',
    'METHICILLIN', 'SULBACTAM-DURLOBACTAM', 'TEMOCILLIN', 'TICARCILLIN',
    'AZITHROMYCIN', 'CLARITHROMYCIN', 'CLINDAMYCIN', 'ERYTHROMYCIN', 'TELITHROMYCIN',
    'DOXYCYCLINE', 'MINOCYCLINE', 'TIGECYCLINE', 'CHLORAMPHENICOL', 'COLISTIN',
    'FOSFOMYCIN', 'FUSIDIC ACID', 'VANCOMYCIN', 'DAPTOMYCIN', 'RIFAMPIN', 'TEICOPLANIN',
})


def _not_run(reason):
    return {'status': 'not_run', 'reason': reason}


def _verified_hydra_sample(report, input_sha256, sample=None):
    if report is None:
        return None, _not_run('No HYDRA assay result was supplied.')
    if isinstance(report, (str, Path)):
        report = load_hydra_report(report)
    rows = report.get('samples') or []
    matches = [row for row in rows if row.get('sample') == sample] if sample is not None else rows
    if len(matches) != 1:
        raise ValueError('Choose exactly one HYDRA sample explicitly; report order is not an identity mapping.')
    chosen = _validate_sample(copy.deepcopy(matches[0]), 1, [])
    source_hash = (_sha256(chosen.get('input_sha256'))
                   or _execution_input_sha256(report.get('execution_provenance') or {}, chosen['sample'], single_sample=len(rows) == 1))
    if source_hash is None:
        return None, {'status': 'ambiguous', 'reason': 'Imported HYDRA evidence has no verifiable matching input SHA-256; explicit user mapping alone is not genomic identity verification.'}
    if source_hash != input_sha256:
        return None, {'status': 'failed', 'reason': 'HYDRA evidence refers to a different assembly SHA-256 and was excluded.'}
    return chosen, {'status': 'completed', 'input_sha256': source_hash,
                    'databases': list(report.get('databases') or []),
                    'version': report.get('hydra_version'),
                    'execution_provenance': copy.deepcopy(report.get('execution_provenance') or {}),
                    'import_provenance': copy.deepcopy(report.get('import_provenance') or {})}


def hydra_report_for_record(record):
    """Reconstitute verified per-isolate evidence, never a global/latest report guess."""
    evidence = (record.get('metadata') or {}).get('hydra') or {}
    current = current_input_sha256(record)
    source_name = evidence.get('source_sample')
    execution = evidence.get('execution_provenance') or {}
    source_hash = (_sha256(evidence.get('evidence_input_sha256'))
                   or _execution_input_sha256(execution, source_name, single_sample=False))
    if not source_name or not current or source_hash != current:
        return None
    upstream = copy.deepcopy(evidence.get('upstream') or {})
    if not upstream:
        upstream = {'sample': source_name, 'hits': copy.deepcopy(evidence.get('hits') or []),
                    'species': copy.deepcopy(evidence.get('species') or {}),
                    'mlst': copy.deepcopy(evidence.get('mlst') or {})}
    # This is the actual report-derived hash, not the link-time baseline.
    upstream['input_sha256'] = source_hash
    return {'hydra_version': execution.get('version'), 'samples': [upstream],
            'databases': copy.deepcopy(evidence.get('databases') or list((execution.get('reference_snapshot') or {}).get('databases') or {})),
            'execution_provenance': copy.deepcopy(execution),
            'import_provenance': copy.deepcopy(evidence.get('provenance') or {})}


def current_characterization(record):
    """Presentation-only recorded input comparison; never hash a file on refresh."""
    evidence = (record.get('metadata') or {}).get('characterization') or {}
    current = current_input_sha256(record)
    source_hash = _sha256(evidence.get('input_sha256'))
    if not evidence:
        status, reason = 'missing', 'No characterization result is attached.'
    elif not current or not source_hash:
        status, reason = 'unverified', 'The saved characterization has no matching recorded input identity.'
    elif current != source_hash:
        status, reason = 'stale', 'Characterization belongs to an earlier or different assembly; archived evidence is not current.'
    else:
        status, reason = 'current', 'Characterization input SHA-256 matches the recorded current assembly.'
    return {'status': status, 'reason': reason, 'input_sha256': current,
            'evidence_input_sha256': source_hash, 'evidence': evidence if status == 'current' else None}


def persist_characterization(project, sample_id, result):
    """Attach metadata atomically without marking MLST typing completed."""
    result = copy.deepcopy(result)
    digest = _sha256(result.get('input_sha256'))
    if not digest or result.get('format_version') != 1:
        raise ValueError('Characterization evidence must have format version 1 and a full input SHA-256.')
    source = Path(result.get('input_path') or '').resolve()
    signature = file_signature(source)
    if file_sha256(source) != digest:
        raise ValueError('Characterization assembly bytes changed before results could be attached.')
    with project.transaction():
        record = project.get_sample(sample_id)
        if not record.get('input_path') or Path(record['input_path']).resolve() != source:
            raise ValueError('Sample input path changed before characterization could be attached.')
        current = current_input_sha256(record)
        if current and current != digest:
            raise ValueError('Characterization input does not match the recorded current assembly SHA-256.')
        if file_signature(source) != signature:
            raise ValueError('Characterization assembly changed while attaching evidence.')
        metadata = record.get('metadata') or {}
        if metadata.get('characterization'):
            project.record_history(sample_id, 'characterization_superseded', {'evidence': metadata['characterization']})
        result['attached_input_path'] = str(source)
        metadata['characterization'] = result
        if not current:
            metadata['input_identity'] = {'sha256': digest, 'path': str(source),
                                          'source': 'full-input hash computed and reverified by characterization'}
        project.set_metadata(sample_id, metadata)
        project.record_history(sample_id, 'characterization_attached', {'input_sha256': digest, 'input_path': str(source),
                               'status': result.get('status'), 'reference_digest': result.get('species_evidence', {}).get('reference_digest')})
    return project.get_sample(sample_id)


def synthesize_accessory_evidence(hydra_sample, source, contig_lengths, *, virulence_hits=(),
                                  contig_headers=None):
    """Join only verified-input evidence. Database classes are not expanded to drugs.

    ``contig_headers`` maps a contig identifier to its full FASTA header, so the
    plasmid view can report what the assembler declared about closure and
    coverage. It is optional: without it those stay unknown rather than absent.
    """
    if hydra_sample is None:
        return {'drug_associations': dict(source, associations=[]),
                'plasmid_hypotheses': dict(source, replicons=[], contig_associations=[])}
    primary = [hit for hit in hydra_sample.get('hits', []) if hit.get('primary') is True]
    amr = [hit for hit in primary if hit.get('element_type') == 'AMR']
    associations, review = [], []
    for hit in amr:
        evidence = {key: copy.deepcopy(hit.get(key)) for key in ('gene', 'accession', 'database', 'class', 'subclass', 'element_subtype', 'method', 'resolution', 'identity_pct', 'coverage_pct', 'sequence', 'start', 'end', 'note')}
        if hit.get('resolution') not in {'COMPLETE', 'POINT'} or 'DISRUPT' in str(hit.get('method', '')).upper() or 'DISRUPT' in str(hit.get('element_subtype', '')).upper():
            review.append(dict(evidence, reason='Partial/disrupted or unclassified resolution: no drug effect was inferred.'))
            continue
        if not hit.get('class') and not hit.get('subclass'):
            review.append(dict(evidence, reason='The reference supplied no drug association; none was invented from the gene name.'))
            continue
        explicit = sorted({label.strip().upper() for label in str(hit.get('subclass') or '').split('/')}
                          & EXPLICIT_DRUG_LABELS)
        associations.append(dict(evidence, affected_drugs=explicit,
                                 annotation_scope='reference-reported class/subclass; specific names only when explicitly present in subclass',
                                 drug_label_source=DRUG_LABEL_SOURCE, drug_label_rule='exact reference subclass token; no class expansion',
                                 interpretation='Genomic resistance association; no susceptible/resistant phenotype inferred.'))
    amr_assayed = bool(amr) or bool(set(source.get('databases', [])) & {'ncbi', 'protein', 'card', 'resfinder', 'argannot', 'megares'})
    drugs = {'status': 'completed' if amr_assayed else 'not_run', 'associations': associations, 'requires_review': review,
             'source': source, 'phenotype': 'not_inferred',
             'reason': 'Specific drug names are retained only when explicitly reported by the reference subclass. Class-level annotations are not expanded to member drugs.',
             'limitations': [LIMITATIONS[1], 'No determinant detected does not mean susceptible. Drug-specific predictions require a validated organism/variant-specific interpretation system and laboratory AST.']}
    replicons = [copy.deepcopy(hit) for hit in primary if hit.get('element_type') == 'PLASMID']
    plasmid_assayed = bool(replicons) or 'plasmidfinder' in source.get('databases', [])

    def location(hit):
        return marker_location(hit, contig_lengths)

    marker_hits = [dict(hit, marker_type='AMR') for hit in amr if hit.get('resolution') in {'COMPLETE', 'POINT'}]
    marker_hits += [dict(hit, marker_type='VIRULENCE') for hit in virulence_hits
                    if hit.get('identity_pct', 0) >= 90 and hit.get('coverage_pct', 0) >= 80]
    links = []
    for replicon in replicons:
        rep_position = location(replicon)
        if rep_position is None or replicon.get('resolution') != 'COMPLETE':
            continue
        for hit in marker_hits:
            position = location(hit)
            if position is None or position[0] != rep_position[0]:
                continue
            gap = max(0, max(position[1], rep_position[1]) - min(position[2], rep_position[2]) - 1)
            links.append({'status': 'hypothesis', 'kind': 'same_assembly_contig', 'contig': position[0],
                          'replicon': replicon['gene'], 'replicon_accession': replicon.get('accession'),
                          'marker': hit['gene'], 'marker_type': hit['marker_type'], 'gap_bp': gap,
                          'replicon_coordinates': list(rep_position[1:]), 'marker_coordinates': list(position[1:]),
                          'interpretation': 'Co-located on one assembled contig; plasmid identity, circularity and transferability remain unproven.'})
    # Contig context for every replicon and every determinant, including the
    # determinants that share no contig with a replicon: their absence from this
    # view would read as absence of the determinant.
    contigs = contig_plasmid_evidence(contig_lengths, replicons, [(hit, hit['marker_type']) for hit in marker_hits],
                                      contig_headers=contig_headers)
    plasmids = {'status': 'completed' if plasmid_assayed else 'not_run', 'replicons': replicons,
                'contig_associations': links, 'source': source, 'reconstruction': 'not_run',
                'mobility': 'not_predicted', 'plasmid_count': 'not_estimated',
                'contig_evidence': contigs['contigs'], 'determinant_placement': contigs['determinant_placement'],
                'depth_basis': contigs['depth_basis'], 'depth_departure_ratio': contigs['depth_departure_ratio'],
                'unplaced_determinants': contigs['unplaced_determinants'], 'mob_suite_gap': list(MOB_SUITE_GAP),
                'reason': 'Replicon assay was supplied.' if plasmid_assayed else 'No plasmid-reference assay was run; the NCBI AMR starter is not a comprehensive plasmid assay.',
                'limitations': [LIMITATIONS[2], 'Identical replicon names in different isolates do not identify the same plasmid.',
                                'Draft assemblies may split, collapse or misassemble repeated plasmid sequence. Long-read/hybrid assembly or validated reconstruction is needed to resolve structure.',
                                *contigs['limitations']]}
    return {'drug_associations': drugs, 'plasmid_hypotheses': plasmids}


def characterize_assembly(path, reference_root=None, hydra_report=None, cancelled=None, progress=None, *, threads=2,
                          hydra_sample_name=None, species=True, virulence=True, modules=None, organism=None,
                          blastn_path=None, makeblastdb_path=None):
    """Run selected independent assays, retaining each result/failure and provenance.

    ``hydra_report`` is a normalized native report (or original report JSON path),
    not metadata merely linked by filename. No report/input identity is invented.
    References must be installed explicitly with provision_characterization_references.
    ``modules`` opts in to registered organism-specific assays by key; unselected
    modules report not_run. ``organism`` is a (genus, species) pair used only to
    stamp how applicable each module's reference panel is to this isolate: it
    never gates execution and no assay writes an organism assignment.
    """
    path = Path(path).resolve()
    qc = inspect_sequence(path, cancelled=cancelled)
    if qc['kind'] != 'fasta':
        raise ValueError('Characterization requires an assembly; assemble FASTQ reads first.')
    result = {'format_version': 1, 'sample_name': qc['sample_name'], 'input_path': str(path),
              'input_sha256': qc['input_sha256'], 'qc': qc['qc'], 'limitations': list(LIMITATIONS),
              'provenance': {'software': 'WMLSTudio', 'version': __version__, 'created_utc': datetime.now(UTC).isoformat()}}
    assays = [('species_evidence', species, identify_species, {}, 'Assay not selected.'),
              ('virulence', virulence, screen_virulence, {'threads': threads, 'blastn_path': blastn_path, 'makeblastdb_path': makeblastdb_path},
               'Assay not selected.')]
    assays += module_tasks(modules, reference_root, threads=threads, blastn_path=blastn_path, makeblastdb_path=makeblastdb_path)
    for key, enabled, function, options, reason in assays:
        check_cancelled(cancelled)
        if not enabled or reference_root is None:
            result[key] = _not_run(reason if not enabled else 'Install or select a verified characterization reference snapshot first.')
            continue
        try:
            result[key] = function(path, reference_root, cancelled=cancelled, progress=progress, **options)
        except AnalysisCancelled:
            raise
        except (ValueError, OSError, RuntimeError) as error:
            result[key] = {'status': 'failed', 'reason': str(error)}
    # A locus ST can only come from the exact-allele engine; the BLAST screen
    # never assigns one, so it is attached here or stays None.
    for key, module in registered_modules().items():
        block = result.get(key)
        if module.locus_st_provider and isinstance(block, dict) and block.get('status') not in {'not_run', 'failed'}:
            apply_locus_sts(result.get('virulence'), module.locus_st_provider(block))
    stamp_applicability(result, organism)
    contig_lengths, contig_headers = {}, {}
    with SequenceReader(path, cancelled) as reader:
        for record in reader:
            if record.identifier in contig_lengths:
                raise ValueError(f'Duplicate contig identifier: {record.identifier}')
            contig_lengths[record.identifier] = len(record.sequence)
            # The full header is kept only so the assembler's own closure and
            # coverage claims can be quoted back; no header text becomes a call.
            contig_headers[record.identifier] = record.name
    chosen, source = _verified_hydra_sample(hydra_report, qc['input_sha256'], hydra_sample_name)
    result.update(synthesize_accessory_evidence(chosen, source, contig_lengths,
                                              virulence_hits=result['virulence'].get('hits', []),
                                              contig_headers=contig_headers))
    if file_sha256(path, cancelled) != qc['input_sha256']:
        raise ValueError('Assembly changed during characterization; results were not accepted.')
    states = {result[key]['status'] for key in ('species_evidence', 'virulence', 'drug_associations',
                                                'plasmid_hypotheses', *registered_modules()) if key in result}
    result['status'] = 'failed' if 'failed' in states else ('ambiguous' if 'ambiguous' in states else ('completed' if 'completed' in states else 'not_run'))
    return result
