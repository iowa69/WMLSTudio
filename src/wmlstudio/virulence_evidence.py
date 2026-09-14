"""Local KpSC-associated virulence gene screening with explicit assay coverage.

Uses curated Kleborate/Pasteur allele sequences, but is a WMLSTudio BLAST+
screen, not a reimplementation or claimed equivalent of Kleborate's typing.
No virulence phenotype or plasmid identity is inferred. A locus ST appears only
when the separate exact-allele locus-ST assay assigned one for that locus.
"""

from __future__ import annotations

from .marker_panel import screen_marker_panel, summarize_marker_hits, unique_locations

LIMITATIONS = ['Not a Kleborate-equivalent result; this BLAST screen assigns no virulence-locus lineage.',
               'Intact coding sequence is not proof of expression or virulence; partial/disrupted matches require review.',
               'Not detected means no hit meeting this defined assay threshold, not proof of genomic absence.',
               'Negative findings are withheld when assembly quality fails explicit minimum screening gates.',
               'These loci may occur outside KpSC. No hypervirulence or plasmid carriage phenotype is inferred.']


def summarize_virulence_hits(hits, loci, *, adequate_negative_assay, locus_sts=None):
    """Per-locus gene calls; a locus ST is carried only when one was assigned elsewhere."""
    assignments = locus_sts or {}
    groups = summarize_marker_hits(hits, loci, adequate_negative_assay=adequate_negative_assay)
    for group in groups:
        assigned = assignments.get(group['locus']) or {}
        group['official_locus_st'] = assigned.get('locus_st')
        group['official_locus_st_source'] = assigned.get('source') if assigned.get('locus_st') else None
    return groups


def _unique_locations(hits):
    """Collapse allele alternatives at one physical locus, not separate copies."""
    return unique_locations(hits, 'gene')


def apply_locus_sts(evidence, assignments):
    """Attach externally assigned locus STs to a completed virulence screen in place.

    The screen itself cannot assign one: only the exact-allele engine can, and
    only for a locus it called completely. Anything else stays ``None``.
    """
    if not isinstance(evidence, dict) or evidence.get('status') != 'completed' or not assignments:
        return evidence
    for group in evidence.get('loci') or []:
        assigned = assignments.get(group.get('locus')) or {}
        if assigned.get('locus_st'):
            group['official_locus_st'] = assigned['locus_st']
            group['official_locus_st_source'] = assigned.get('source')
    return evidence


def screen_virulence(path, reference_root, cancelled=None, progress=None, *, threads=2,
                     blastn_path=None, makeblastdb_path=None, locus_sts=None):
    def summarize(hits, groups, *, adequate_negative_assay):
        return summarize_virulence_hits(hits, groups, adequate_negative_assay=adequate_negative_assay, locus_sts=locus_sts)
    return screen_marker_panel(path, reference_root, cancelled, progress, section='virulence',
                               assay_name='WMLSTudio KpSC-associated virulence allele screen',
                               progress_message='Screening six defined KpSC-associated virulence loci with native BLAST+…',
                               limitations=LIMITATIONS, subject='Virulence', summarize=summarize,
                               temp_prefix='wmlstudio-virulence-', threads=threads,
                               blastn_path=blastn_path, makeblastdb_path=makeblastdb_path)
