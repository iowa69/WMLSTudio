"""SCCmec marker screening for Staphylococcus aureus and coagulase-negative staphylococci.

This is a WMLSTudio BLAST+ screen over the pinned public rpetit3/sccmec target
and cassette reference panel (MIT). It is not staphopia-sccmec, not SCCmecFinder
and not equivalent to either, and it is not a validated MRSA diagnostic: no
methicillin phenotype, susceptibility or infection-control conclusion is
inferred. SCCmec is a large, repetitive, IS-rich element that routinely spans
contig breaks in draft assemblies, so a partial or untypeable result is reported
as exactly that and is never rounded to the nearest type definition.
"""

from __future__ import annotations

import html

from .marker_panel import screen_marker_panel, union_coverage
from .organism_modules import OrganismMatch, OrganismModule, register
from .sequence import AnalysisCancelled

LIMITATIONS = [
    'mecC is not in this reference panel and is not assayed here; a mecC-carrying element (SCCmec XI) '
    'is not detected by mecA. The AMR path reports mecC separately, and the two are deliberately not '
    'joined into a type call.',
    'SCCmec is a 14-82 kb IS-rich element that frequently spans contig breaks in draft assemblies: a '
    'detected target proves genome-wide presence, not cassette membership; a missing target may be '
    'assembly collapse rather than true absence; and IS431 and ccr genes recur elsewhere in '
    'staphylococcal genomes. Required targets on different contigs withhold the type call.',
    'Cassette region support is corroboration, not the call. Coverage is summed across contigs, so a '
    'high value is consistent with - never proof of - one intact cassette, and subtypes are offered '
    'only as a provisional reference match.',
    'Coagulase-negative staphylococci frequently carry novel or non-typeable SCCmec elements. The panel '
    'was curated on S. aureus; a negative result on other staphylococci is weak evidence.',
    'This is a WMLSTudio BLAST+ marker screen against a pinned public reference panel. It is not '
    'staphopia-sccmec or SCCmecFinder output and is not equivalent to them. No methicillin/oxacillin '
    'susceptibility, no MRSA/MSSA designation, no infection-control decision and no phenotype is '
    'inferred. Confirmation requires laboratory AST.',
    'Not detected means no hit met this defined assay threshold, not proof of genomic absence. Negative '
    'findings are withheld when assembly quality fails the explicit minimum gates.',
]
TARGET_GATES = {'min_identity': 90, 'min_coverage': 80, 'review_identity': 80, 'review_coverage': 40}
# A staphylococcal genome is roughly 2.8 Mb, so the inherited 1 Mb gate is too
# weak to make a negative defensible here.
TARGET_QC_GATES = (2_000_000, 10_000, 1.0)
REGION_GATES = {'min_identity': 85, 'min_coverage': 83, 'review_identity': 80}


def _name_state(name, detected, aliases):
    """Resolve one rule name - a bare target or a ccr/mec alias - to a single state."""
    if name not in aliases:
        status = detected.get(name, 'not_detected')
        return {'detected': 'present', 'ambiguous': 'unresolved'}.get(status, 'absent')
    states = [detected.get(member, 'not_detected') for member in aliases[name]]
    if any(state == 'ambiguous' for state in states):
        return 'unresolved'
    if all(state == 'detected' for state in states):
        return 'present'
    if not any(state == 'detected' for state in states):
        return 'absent'
    return 'partial'


def _alias_kind(name):
    """The pinned panel names its complexes 'ccr Type n' and its classes 'mec Class X'."""
    lowered = name.casefold()
    return 'ccr' if lowered.startswith('ccr') else ('mec' if lowered.startswith('mec') else 'other')


def _contigs_for(targets, locations):
    contigs = []
    for target in targets:
        for location in locations.get(target) or []:
            if location[0] not in contigs:
                contigs.append(location[0])
    return sorted(contigs)


def interpret_sccmec(detected, locations, rules, *, region_support=None, qc_adequate=True):
    """Pure ccr-complex/mec-class decision kernel; no BLAST, no reference access.

    ``detected`` maps every panel target to 'detected', 'ambiguous' or
    'not_detected'; ``locations`` maps a target to its accepted (contig, start,
    end, strand) placements. The type definitions come from the manifest, never
    from Python, so the rules in force are provable from the reference digest.
    """
    aliases = {alias['name']: list(alias['targets']) for alias in rules.get('aliases') or []}
    alias_states = {name: _name_state(name, detected, aliases) for name in aliases}
    complexes = [{'name': name, 'state': alias_states[name], 'kind': _alias_kind(name),
                  'targets': {member: detected.get(member, 'not_detected') for member in aliases[name]}}
                 for name in aliases]
    candidates = []
    for definition in rules.get('types') or []:
        required = list(definition.get('targets') or [])
        excluded = list(definition.get('excludes') or [])
        if not required or not all(_name_state(name, detected, aliases) == 'present' for name in required):
            continue
        if any(_name_state(name, detected, aliases) == 'present' for name in excluded):
            continue
        candidates.append(definition)
    mec_a = detected.get('mecA', 'not_detected')
    result = {'status': 'ambiguous', 'type': None, 'official_type': None,
              'candidate_types': [definition['name'] for definition in candidates],
              'ccr_complexes': [row for row in complexes if row['kind'] == 'ccr'],
              'mec_classes': [row for row in complexes if row['kind'] == 'mec'],
              'other_rules': [row for row in complexes if row['kind'] == 'other'],
              'targets': dict(detected), 'mecA': mec_a, 'mecC': 'not_assayed',
              'contig_fragmented': False, 'contigs': [], 'ccr_mec_same_contig': None,
              'region_support': region_support or {'status': 'not_run', 'reason': 'The cassette region pass was not run.'},
              'phenotype': 'not_inferred', 'qc_adequate': bool(qc_adequate), 'reason': ''}
    if len(candidates) == 1:
        definition = candidates[0]
        required = list(definition.get('targets') or [])
        expanded = [member for name in required for member in (aliases.get(name) or [name])]
        ccr_names = [name for name in required if name in aliases and _alias_kind(name) == 'ccr']
        mec_names = [name for name in required if name in aliases and _alias_kind(name) == 'mec']
        ccr_contigs = _contigs_for([member for name in ccr_names for member in aliases[name]], locations)
        mec_contigs = _contigs_for([member for name in mec_names for member in aliases[name]], locations)
        result['contigs'] = _contigs_for(expanded, locations)
        result['contig_fragmented'] = len(result['contigs']) > 1
        if ccr_contigs and mec_contigs:
            result['ccr_mec_same_contig'] = bool(set(ccr_contigs) & set(mec_contigs))
        unresolved = sorted(name for name, state in alias_states.items() if state == 'unresolved')
        complexes_resolved = all(alias_states[name] == 'present' for name in (*ccr_names, *mec_names))
        if unresolved:
            result['reason'] = ('One or more ccr complexes or mec classes could not be resolved '
                                f'({", ".join(unresolved)}); the type call is withheld.')
        elif not complexes_resolved:
            result['reason'] = 'The ccr complex and mec class of the matching definition are not both fully resolved.'
        elif result['contig_fragmented']:
            result['reason'] = ('The required targets lie on more than one contig; SCCmec commonly spans contig '
                                'breaks, so the type call is withheld rather than assumed.')
        else:
            result.update(status='completed', type=definition['name'],
                          reason='Exactly one type definition in this reference panel is satisfied, with both '
                                 'complexes fully resolved on a single contig.')
    elif len(candidates) > 1:
        result['contigs'] = _contigs_for([target for target, state in detected.items() if state == 'detected'], locations)
        result['contig_fragmented'] = len(result['contigs']) > 1
        result['reason'] = ('More than one type definition is satisfied: review a possible composite element, '
                            'more than one SCCmec element, or a mixed assembly.')
    else:
        result['contigs'] = _contigs_for([target for target, state in detected.items() if state == 'detected'], locations)
        result['contig_fragmented'] = len(result['contigs']) > 1
        if mec_a == 'detected':
            result['reason'] = ('mecA is present but the detected ccr/mec combination matches no type definition in '
                                'this reference panel: possible novel, composite or truncated element, or assembly '
                                'fragmentation.')
        elif not qc_adequate:
            result['reason'] = ('Assembly quality failed the explicit minimum gates, so a negative SCCmec result is '
                                'withheld rather than reported as absence.')
        else:
            result['status'] = 'not_detected'
            result['reason'] = ('No mecA and no typeable ccr/mec combination was detected by this assay. mecC is not '
                                'in this panel and was not assayed.')
    return result


def _summarize_regions(hits, groups, *, adequate_negative_assay, regions):
    """Union cassette coverage per reference; a subtype is never promoted to a call."""
    by_gene = {}
    for hit in hits:
        if hit['identity_pct'] >= REGION_GATES['min_identity']:
            by_gene.setdefault(hit['gene'], []).append(hit)
    rows = []
    for entry in regions:
        coverage = union_coverage(by_gene.get(entry['gene']) or [], entry['reference_bp'])
        rows.append(dict(coverage, subtype=entry['gene'], accession=entry.get('accession'),
                         reference_bp=entry['reference_bp'],
                         status='provisional_region_match' if coverage['coverage_pct'] >= REGION_GATES['min_coverage']
                         else 'insufficient_coverage'))
    rows.sort(key=lambda row: (-row['coverage_pct'], -row['identity_pct'], row['subtype']))
    return rows


def _region_support(path, reference_root, regions, cancelled, progress, options):
    if not regions:
        return {'status': 'not_run', 'reason': 'This reference snapshot carries no SCCmec cassette references.'}
    def summarize(hits, groups, *, adequate_negative_assay):
        return _summarize_regions(hits, groups, adequate_negative_assay=adequate_negative_assay, regions=regions)
    try:
        evidence = screen_marker_panel(
            path, reference_root, cancelled, progress, section='sccmec.regions',
            assay_name='WMLSTudio SCCmec cassette region coverage screen',
            progress_message='Comparing assembled contigs against pinned SCCmec cassette references with native BLAST+…',
            limitations=[LIMITATIONS[2]], subject='SCCmec region', blast_task='megablast',
            min_identity=REGION_GATES['min_identity'], min_coverage=REGION_GATES['min_coverage'],
            review_identity=REGION_GATES['review_identity'], review_coverage=0,
            qc_gates=TARGET_QC_GATES, summarize=summarize,
            temp_prefix='wmlstudio-sccmec-regions-', **options)
    except AnalysisCancelled:
        raise
    except (ValueError, OSError, RuntimeError) as error:
        return {'status': 'failed', 'reason': str(error)}
    ranked = evidence['loci']
    best = ranked[0] if ranked else None
    return {'status': best['status'] if best else 'insufficient_coverage',
            'best': best, 'ranked': ranked[:10], 'subtype': None,
            'subtype_interpretation': 'Reference-supported cassette similarity only; never promoted to a type or '
                                      'subtype call.',
            'gates': dict(REGION_GATES), 'assay': evidence['assay'], 'provenance': evidence['provenance']}


def type_sccmec(path, reference_root, cancelled=None, progress=None, *, threads=2,
                blastn_path=None, makeblastdb_path=None):
    """Screen the pinned SCCmec target panel, then corroborate with cassette references."""
    options = {'threads': threads, 'blastn_path': blastn_path, 'makeblastdb_path': makeblastdb_path}
    evidence = screen_marker_panel(
        path, reference_root, cancelled, progress, section='sccmec.targets',
        assay_name='WMLSTudio SCCmec ccr/mec marker screen',
        progress_message='Screening pinned SCCmec ccr, mec and IS targets with native BLAST+…',
        limitations=LIMITATIONS, subject='SCCmec target', qc_gates=TARGET_QC_GATES,
        temp_prefix='wmlstudio-sccmec-targets-', **TARGET_GATES, **options)
    manifest_section = _manifest_section(reference_root)
    group = (evidence['loci'] or [{}])[0]
    detected, locations = {}, {}
    for member in group.get('genes') or []:
        detected[member['gene']] = member['status']
        accepted = [hit for hit in member['hits']
                    if hit['identity_pct'] >= TARGET_GATES['min_identity']
                    and hit['coverage_pct'] >= TARGET_GATES['min_coverage']]
        locations[member['gene']] = [(hit['contig'], hit['start'], hit['end'], hit['strand']) for hit in accepted]
    region_support = _region_support(path, reference_root, manifest_section.get('regions') or [],
                                     cancelled, progress, options)
    result = interpret_sccmec(detected, locations, manifest_section.get('rules') or {},
                              region_support=region_support,
                              qc_adequate=evidence['assay']['adequate_negative_assay'])
    result.update(input_sha256=evidence['input_sha256'], reference_digest=evidence['reference_digest'],
                  target_locations={target: [list(location) for location in places]
                                    for target, places in locations.items() if places},
                  hits=evidence['hits'], assay=evidence['assay'], provenance=evidence['provenance'],
                  rules_schema_version=manifest_section.get('schema_version'),
                  not_assayed=list(manifest_section.get('not_assayed') or []),
                  limitations=list(LIMITATIONS))
    return result


def _manifest_section(reference_root):
    from .characterization_refs import validate_characterization_references
    return validate_characterization_references(reference_root).get('sccmec') or {}


def _short(name):
    lowered = name.casefold()
    if lowered.startswith('ccr type '):
        return 'ccr' + name.split()[-1]
    if lowered.startswith('mec class '):
        return 'mec ' + name.split()[-1]
    return name


def summarize_sccmec(evidence):
    """One table cell; a withheld call always reads as withheld."""
    status = (evidence or {}).get('status')
    if status in (None, 'not_run'):
        return 'not_run'
    if status == 'failed':
        return 'failed'
    if status == 'not_detected':
        return 'no mecA or typeable ccr/mec detected'
    present = [_short(row['name']) for row in (*(evidence.get('ccr_complexes') or []),
                                               *(evidence.get('mec_classes') or [])) if row['state'] == 'present']
    if evidence.get('type'):
        detail = ' + '.join(present)
        return f'{evidence["type"]} (candidate)' + (f' - {detail}' if detail else '')
    candidates = evidence.get('candidate_types') or []
    if len(candidates) > 1:
        return f'ambiguous - {len(candidates)} candidate types'
    if candidates:
        return f'ambiguous - {candidates[0]} withheld'
    return 'ambiguous - untypeable by this panel'


def _escape(value):
    return html.escape(str(value if value is not None else ''))


def sccmec_html(evidence):
    """Drill-down and report fragment; every value is escaped and limitations repeat verbatim."""
    evidence = evidence or {}
    escape = _escape
    status = evidence.get('status', 'not_run')
    parts = ['<h4>SCCmec marker screen (Staphylococcus)</h4>',
             '<p><b>Result:</b> ' + escape(summarize_sccmec(evidence)) + ' &middot; status ' + escape(status) + '</p>']
    if status in {'not_run', 'failed'}:
        parts.append('<p>' + escape(evidence.get('reason') or 'This assay was not run for this isolate.') + '</p>')
    else:
        if evidence.get('reason'):
            parts.append('<p>' + escape(evidence['reason']) + '</p>')
        parts.append('<p><b>Applicability:</b> ' + escape(evidence.get('applicability') or 'unknown_organism') +
                     ' &middot; ' + escape(evidence.get('applicability_reason') or '') + '</p>')
        parts.append('<p><b>mecA:</b> ' + escape(evidence.get('mecA')) + ' &middot; <b>mecC:</b> ' +
                     escape(evidence.get('mecC')) + ' (not in this panel) &middot; <b>official type:</b> not assigned '
                     'by this software</p>')
        rows = ''.join('<tr><td>' + escape(row['name']) + '</td><td>' + escape(row['state']) + '</td><td>' +
                       escape(', '.join(f'{gene}: {state}' for gene, state in sorted(row['targets'].items()))) +
                       '</td></tr>'
                       for row in (*(evidence.get('ccr_complexes') or []), *(evidence.get('mec_classes') or [])))
        if rows:
            parts.append('<table><tr><th>Complex or class</th><th>State</th><th>Targets</th></tr>' + rows + '</table>')
        parts.append('<p><b>Candidate types:</b> ' + escape(', '.join(evidence.get('candidate_types') or []) or 'none') +
                     ' &middot; <b>contigs:</b> ' + escape(', '.join(evidence.get('contigs') or []) or 'none') +
                     ' &middot; <b>fragmented:</b> ' + escape(evidence.get('contig_fragmented')) + '</p>')
        support = evidence.get('region_support') or {}
        best = support.get('best') or {}
        if best:
            parts.append('<p><b>Cassette region support:</b> ' + escape(best.get('subtype')) + ' (' +
                         escape(best.get('accession')) + ') ' + escape(f'{best.get("coverage_pct", 0):.1f}') +
                         '% union coverage over ' + escape(best.get('contig_count')) + ' contig(s), identity ' +
                         escape(f'{best.get("identity_pct", 0):.1f}') + '% &middot; ' + escape(support.get('status')) +
                         '. Never promoted to a type or subtype call.</p>')
    for limitation in evidence.get('limitations') or LIMITATIONS:
        parts.append('<p class="muted">' + escape(limitation) + '</p>')
    return ''.join(parts)


MODULE = register(OrganismModule(
    key='sccmec', title='SCCmec typing (Staphylococcus)', column_title='SCCmec',
    match=OrganismMatch(genera=frozenset({'Staphylococcus'}), species=frozenset({'aureus'})),
    runner=type_sccmec, option_keys=('threads', 'blastn_path', 'makeblastdb_path'),
    manifest_sections=('sccmec.targets', 'sccmec.rules'),
    summary=summarize_sccmec, detail_html=sccmec_html, report_default=True,
    limitations=tuple(LIMITATIONS)))
