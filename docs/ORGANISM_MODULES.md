# Organism-specific typing modules

Some questions only make sense for one genus. SCCmec typing is a staphylococcal
question; *Klebsiella* virulence-locus STs and capsule markers are *Klebsiella*
questions. WMLSTudio keeps those assays in a registry so that the interface can
offer the right ones for the isolate in front of you, run them only when you ask,
and record what they can and cannot support.

Nothing here is a validated diagnostic assay, and no module is equivalent to the
published tool its reference data came from. Read [the claim boundary](#the-claim-boundary)
before reading the results.

## The three modules that ship

| Key | Offered for | Reads | Answers |
| --- | --- | --- | --- |
| `sccmec` | *Staphylococcus* (curated on *S. aureus*) | `sccmec.targets`, `sccmec.rules` | Which ccr complex and mec class were detected, and which SCCmec type definitions they satisfy |
| `klebsiella_locus_st` | *Klebsiella* | `locus_profiles` | ybt / clb / iuc / iro / rmp locus STs from exact allele vectors, with the published lineage label |
| `klebsiella_capsule` | *Klebsiella* | `capsule.loci` | A *wzi* and a *wzc* allele number |

Each is a separate registry entry on purpose: a reference snapshot can carry the
capsule markers but not the locus-ST profile tables, and each module then
degrades on its own instead of taking the other down with it.

## How an organism decides what is offered

A module's match rule has **three** outcomes, not two, because running an assay
outside the taxa its panel was curated on is a legitimate thing to do as long as
the result says so.

| State | Meaning | What happens |
| --- | --- | --- |
| `recommended` | Genus and species match the curated taxa | Offered, and pre-selected in the plan |
| `possible` | Genus matches, species is outside the curated set — a coagulase-negative *Staphylococcus*, for example | Offered, runs normally, and the stored evidence is stamped `possible` |
| `off_panel` | Genus is outside the panel, or the species is explicitly excluded | Still selectable, still runs, stamped `off_panel` |
| `unknown_organism` | Neither genus nor species is assigned or detected | Still selectable, stamped `unknown_organism` |

**The match rule never blocks execution.** It decides what is recommended and it
stamps `applicability` and `applicability_reason` onto the stored evidence, so a
result produced outside its validation taxa carries that fact wherever it travels
— into the drill-down, the report and the export.

The plan dialog adds one layer on top of that engine behaviour: off-panel isolates
are left out of a module run unless you tick the box that includes them, and if
you select a module for which no chosen isolate is even `possible`, it says so
before you start. That is a convenience, not a rule — the engine will run any
module against any isolate and record what it was.

The genus and species come from the existing resolution chain: what you assigned,
otherwise what the identification evidence proposed. No module writes
`metadata['organism']`. A module may make its result depend on the organism; only
you (or an explicit identification decision you accept) can change the organism.

## The claim boundary

Every module carries its own limitations and repeats them in every drill-down and
every report section. They are reproduced here verbatim so that they can be read
without opening the application.

### SCCmec

1. mecC is not in this reference panel and is not assayed here; a mecC-carrying element (SCCmec XI) is not detected by mecA. The AMR path reports mecC separately, and the two are deliberately not joined into a type call.
2. SCCmec is a 14-82 kb IS-rich element that frequently spans contig breaks in draft assemblies: a detected target proves genome-wide presence, not cassette membership; a missing target may be assembly collapse rather than true absence; and IS431 and ccr genes recur elsewhere in staphylococcal genomes. Required targets on different contigs withhold the type call.
3. Cassette region support is corroboration, not the call. Coverage is summed across contigs, so a high value is consistent with - never proof of - one intact cassette, and subtypes are offered only as a provisional reference match.
4. Coagulase-negative staphylococci frequently carry novel or non-typeable SCCmec elements. The panel was curated on S. aureus; a negative result on other staphylococci is weak evidence.
5. This is a WMLSTudio BLAST+ marker screen against a pinned public reference panel. It is not staphopia-sccmec or SCCmecFinder output and is not equivalent to them. No methicillin/oxacillin susceptibility, no MRSA/MSSA designation, no infection-control decision and no phenotype is inferred. Confirmation requires laboratory AST.
6. Not detected means no hit met this defined assay threshold, not proof of genomic absence. Negative findings are withheld when assembly quality fails the explicit minimum gates.

### *Klebsiella* virulence locus STs

1. Exact nucleotide allele matching only: an inexact or novel allele yields an incomplete call and no locus ST, never the nearest ST.
2. The profile tables are a pinned snapshot. A locus ST registered upstream after the pinned revision is absent here, and a profile referencing an allele with no reference sequence cannot be assigned.
3. The lineage value is a published Kleborate profile-table lookup reported verbatim, not an inference by this software.
4. A locus ST is not a hypervirulence phenotype, not plasmid identity and not a Kleborate virulence or resistance score; those aggregate scores are not computed here.
5. This is a WMLSTudio exact-allele call using Kleborate public allele and profile data. It is not Kleborate output and is not equivalent to Kleborate.
6. rmpA2 has no upstream profile table and therefore never acquires a locus ST; it stays allele-only in the virulence screen.

### *Klebsiella* capsule markers

1. The result is a wzi or wzc allele number and nothing else. The published wzi-allele to K-type associations are neither shipped nor applied: no K locus, capsule type or serotype is inferred.
2. This is a WMLSTudio exact-allele screen over the public Kaptive wzi/wzc marker database. It is not Kaptive output and is not equivalent to Kaptive: the whole K/O locus references, Kaptive match-confidence grading and the O-locus special logic are not implemented here.
3. wzc entries are short (115-151 bp) variable-region fragments and are more prone to several exact matches; ambiguous and mixed calls are surfaced unchanged rather than tie-broken.
4. Exact nucleotide matching only: a novel or inexact marker sequence is reported as missing, never as the nearest allele.
5. No capsule marker allele establishes virulence, transmission or a phenotype.

## When SCCmec reports a type, and when it refuses

The IWG ccr/mec type table is **not** written in WMLSTudio source. It is derived
from the pinned upstream rule document at staging time and stored in the
reference manifest, so the rules actually in force are provable from the recorded
`reference_digest`, and updating the panel is a re-stage rather than a code
change.

A ccr complex or mec class resolves to `present` when every member target was
detected, `absent` when none was, `partial` in between, and `unresolved` when any
member was ambiguous. A type is a **candidate** when every name it requires is
present and no name it excludes is.

A single `type` is reported only when **all four** hold:

- exactly one candidate type survives,
- no complex or class is `unresolved`,
- both the ccr complex and the mec class are fully resolved, and
- the required targets are on **one contig** (`contig_fragmented` is false).

Otherwise the candidates are listed and `type` stays `None` with the reason.
`official_type` is unconditionally `None`: this software never issues an official
designation. Two or more candidates, or mecA present with no matching definition,
are reported as `ambiguous` — a composite, novel, truncated or fragmented element
— never as a negative. When assembly QC fails the stated minimum gates, every
`not_detected` becomes `ambiguous` rather than a withheld-quality negative
masquerading as absence. `mecC` is always `not_assayed`.

## Reference sources, pinned

| Source | Revision | Licence | What is staged |
| --- | --- | --- | --- |
| [Kleborate](https://github.com/klebgenomics/Kleborate) | `550ce22a2c01c76064f4dabf403704ee2293356e` | GPL-3.0-or-later | Species references, virulence allele FASTAs, five `profiles.tsv` locus-ST tables |
| [rpetit3/sccmec](https://github.com/rpetit3/sccmec) v1.2.0 | `b901cc618be8eb17284ccb0cf6ef9ee428d909c3` | MIT, © 2024 Robert A. Petit III | 20 target FASTAs, 32 cassette-region FASTAs, the rule document and its cross-check table |
| [Kaptive](https://github.com/klebgenomics/Kaptive) v2.0.9 | `b3856eac6e76b3017aa993319da2a8ea967a1ba0` | GPL-3.0-or-later | `wzi` (484 alleles) and `wzc` (120 alleles) only |

Full attribution, byte counts and SHA-256 values are in
[THIRD_PARTY_NOTICES](../studio_packaging/THIRD_PARTY_NOTICES.md). Derivations —
splitting multi-record FASTAs, rewriting headers, parsing the rule document —
happen at staging time only. The upstream files are retained verbatim beside the
derived artefacts so every transform is auditable against the downloaded bytes.

## Reference snapshots and the migration path

The manifest has two format versions. **Format 1** predates these modules;
**format 2** additionally declares `sources`, `locus_profiles`, `sccmec` and
`capsule`, and a format 2 manifest must declare all four even when empty, so a
section can never be silently stripped.

An installed format 1 snapshot keeps validating byte for byte. On it, the
organism modules report:

> `not_run` — "This reference snapshot (format 1) carries no `sccmec` section for
> SCCmec typing (Staphylococcus). Install an updated snapshot from
> Characterization, Install / update."

That is a **diagnosable reference state, not a negative result**, and the
distinction is deliberate: a missing panel must never read as "no SCCmec
detected". Re-staging produces a new digest and therefore a new snapshot
directory; nothing is ever mutated in place.

```sh
python studio_scripts/stage_characterization.py --destination <a fresh path>
```

The frozen build verifies which state it shipped in: `check_frozen.py` reports
`reference_panels.characterization.organism_modules` as `staged` or `absent` with
the reason, and runs a bundled SCCmec self-comparison against the panel's own
type IVa cassette reference when the panel and the command line both support it.
That control is labelled exactly as the species control is — *bundled reference
self-comparison; not an independent biological validation.*

## Deferred: whole K and O locus typing

This is designed and **not implemented**, with the reason recorded rather than
the capability quietly promised.

Kaptive v2.0.9 (commit `b3856eac6e76b3017aa993319da2a8ea967a1ba0`,
GPL-3.0-or-later) ships `Klebsiella_k_locus_primary_reference.gbk` (8,325,855
bytes), `Klebsiella_o_locus_primary_reference.gbk` (325,387 bytes),
`Klebsiella_k_locus_variant_reference.gbk` (1,303,472 bytes) and
`Klebsiella_o_locus_primary_reference.logic` (591 bytes). Kaptive 3 master no
longer carries the databases in-repo, so v2.0.9 is the pinnable artefact.

Two blockers:

1. **Format.** They are GenBank. The runtime depends on PySide6, pyahocorasick,
   pyrodigal and pyskani — there is no GenBank parser and no justification for
   adding one to the runtime. It would need a build-time converter to whole-locus
   FASTA, per-CDS gene FASTA and a locus-to-gene JSON, with Biopython in the
   development group only and the manifest recording the upstream `.gbk`
   SHA-256, the derived SHA-256 and the converter version.
2. **Algorithm.** Kaptive's output is not a best-BLAST-hit. It is a locus
   assignment with a confidence grade derived from expected-gene coverage,
   missing and extra genes and locus contiguity, plus `.logic` special rules for
   O-locus subtypes. Reimplementing that and calling the output a "K locus type"
   is exactly the overstatement this project refuses.

If it is ever pursued, the honest shape is a `capsule_locus_candidate` chosen by
interval-union coverage with per-expected-gene presence and a contig count,
status capped at `provisional_reference_match`, `official_locus` unconditionally
`None`, and the limitation naming what is not implemented. The bundle cost is
about **+9.9 MB**.

## Adding a module

A module is a frozen record: a key that must not collide with the four core
characterization sections, a title and column title, a match rule, a runner, the
option keys to forward, the manifest sections it requires, a one-cell summary and
a drill-down renderer. `register()` refuses a colliding key. `module_tasks()`
turns a selection into the task shape `characterize_assembly` already consumes,
and a selected module whose manifest sections are absent yields `enabled=False`
with a reason that names the missing section.

Two rules are structural rather than stylistic. The registry must not import
`characterization.py` — summary and drill-down callables operate on plain
evidence dictionaries. And because the registry imports its assay modules by
computed name, every new assay module must be reachable by
`stage_reference_panels.assay_module_imports()`, which is what puts it into the
frozen bundle; the packaging test fails if it is not.
