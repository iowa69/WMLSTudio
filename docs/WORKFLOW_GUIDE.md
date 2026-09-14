# WMLSTudio: from a problem to reviewable evidence

This handbook is available offline inside **Help → Problem → solution guide**.
Start with a research question, keep one stable record per isolate, and reuse
the same evidence in tables, graphs and reports. A project is the evidence store;
an investigation is a named, explicitly selected comparison within that store.

The application supports non-commercial research. It is not a validated
diagnostic device. Genetic proximity does not establish a transmission event;
genomic resistance determinants are not measured antimicrobial susceptibility.

You do not need to be a bioinformatician to use it. You do need to read what each
result says about itself. Where something described here is not yet fully wired in
this build, the text says so in italics at that point. The reach of every
capability in this revision is listed in the
[capability audit](STUDIO_CAPABILITIES.md).

## Finding your way around

Seven tabs run across the top of the window. Each one answers a different
question, and each says so in a plain sentence directly under the tab bar, with a
button for the usual next step and a **?** that opens this guide at the right
place.

| Tab | The question it answers |
| --- | --- |
| **Overview** | What is in this project, and what should I do next? |
| **Isolates** | Which isolates do I have, and is each one trustworthy? |
| **Compare** | How close are these isolates to each other? |
| **Schemes** | Which reference definitions do I have installed? |
| **Evidence** | What genes and markers were found, by which assay? |
| **Reports** | What can I hand to a colleague? |
| **Settings** | Text size, screen size, reference data, and what this version cannot do |

The tabs are connected but they do not silently drag each other around. Selecting
isolates anywhere sets a **focus**: the strip at the top says how many isolates
are focused and where they came from — "Focus · 7 isolates · from Graph selection
at 14:32". Focus alone changes nothing. A tab starts reviewing those isolates
only when you press **Use current focus (N)** on that tab, and the tab then states
what it is reviewing and where that cohort came from: "Reviewing 7 · from Graph
selection". That cohort belongs to that tab only. Two tabs can be looking at two
different sets on purpose, and each will tell you which.

The Compare strip says *similarity is not proof of transmission*; the Evidence
strip says *genotype is not measured susceptibility*. Those limits are on the
screen, not buried in a tooltip.

## I have never done this before: can I practise on real data first?

Yes, and you should. Two practice cohorts of published complete genomes are
pinned by accession and checksum. WMLSTudio ships none of the sequences — it
downloads them from NCBI when you ask, and verifies every byte against the pinned
size, MD5 and SHA-256 as well as against the checksums NCBI publishes today.

- **`kpneumoniae-10`** — ten *K. pneumoniae* genomes, about 16 MB. A single-species
  cohort, so you can see clustering, thresholds and a report end to end.
- **`mixed-genus-20`** — twenty genomes across fifteen genera, about 21 MB. Its
  job is to show automatic identification and filing working, *and failing
  honestly*: a *K. variicola* that is a genuinely hard call inside the *K.
  pneumoniae* complex, *E. coli* which cannot be separated from *Shigella* by
  whole-genome identity alone, and *Enterobacter* and *Serratia* which have no
  species reference in the bundled panel and must come out genus-only or
  unresolved.

The files arrive flat, in one folder, with no organism structure — creating the
folders is the demonstration, so the download must not do it for you. No expected
ST, cluster, threshold or "correct answer" ships with either cohort: agreement
with any published investigation would not validate this software.

Per-accession tables, the true organism quoted from each NCBI assembly report and
the full provenance are in [practice cohorts](TEST_DATASETS.md).

Use **Data → Download practice data…**, read the caveats the dialog shows you,
and choose a cohort. The same cohorts are available without the interface:
`python studio_scripts/fetch_practice_cohort.py --list`, then
`--cohort mixed-genus-20 --destination <folder>`. `--verify-pins` re-reads NCBI's
published checksums without downloading a single genome, so a stale pin is
detected cheaply.

## I have 100 presumptive Klebsiella pneumoniae isolates

Your investigation has six connected stages:

**Define isolates → check identity/QC → establish ST/profiles → inspect
resistance/virulence → investigate relatedness → freeze and communicate.**

1. Create a project named for the investigation, not for a single sequencing run.
2. [Import assemblies or reads](wmlstudio:import). Review file-to-isolate identity
   and organism assignments. Names are suggestions; stable IDs join the evidence.
3. Review quality and independent organism evidence. Quarantine conflicts instead
   of forcing all files into the expected species.
4. Establish classical MLST, then choose the appropriate cgMLST/wgMLST snapshot.
5. Inspect linked gene evidence and its limitations alongside the comparison.
6. Save the investigation, review threshold clusters and export a defined cohort.

Do not start by forcing a visually appealing tree. A trustworthy disconnected
node is preferable to an edge based on too little comparable sequence.

## Are all isolates the same species or subspecies?

Three different things must remain visible: the **organism you assigned**, the
**lineage suggested by MLST**, and **independent genome-reference evidence**.
They may disagree. A species-complex MLST scheme can match multiple related taxa.

Independent characterization uses a versioned reference panel and reports ANI,
aligned fractions and competing matches. Reference scope, assembly quality and
the separation from the nearest alternative matter. “No suitable match” means
unresolved, not proof of a novel species. Species-level evidence must not silently
become a subspecies assignment; a dedicated validated resolution method is needed.

[Review characterization](wmlstudio:characterize), then use the sample inspector
to trace the assembly hash, reference accessions and result state. An assigned
organism is never silently overwritten by a weaker computational result.

Background: [skani methods](https://doi.org/10.1038/s41592-023-02018-3),
[Kleborate's organism-specific analyses](https://github.com/klebgenomics/Kleborate).

## Where did my files go? Automatic folders for inputs and references

When you import a batch, each file is identified *before* anything is copied, and
the proposal is shown to you with its evidence: the nearest reference, the margin
over the runner-up, and a confidence word you can hover for the reason. You accept,
change, or defer each row. Only then are the managed copies written.

Accepted isolates are filed under `Genus/species`. Anything that was not resolved
goes to a **`_Unresolved`** folder split by *why*, not by guesswork:

| Folder | What it means |
| --- | --- |
| `Awaiting_identification` | Proposed, not yet reviewed by you |
| `Low_confidence` | The evidence did not reach the confidence your policy requires |
| `Conflicting_evidence` | Two references were too close together to separate |
| `Reads_not_assembled` | A FASTQ file; identification needs an assembly first |
| `Not_in_reference_panel` | Nothing in the installed panel matched — a statement about the panel, not about the organism |
| `User_deferred` | You chose to decide later |

**A folder is a filing decision, not a laboratory identification.** It records
where a copy is stored. Confirm the organism before any clinical interpretation.
Your original files are never moved, renamed or deleted; only the managed copies
are organised, and the original hash is kept so a copy can always be traced back.

Nothing is auto-confirmed below genus level however you set the policy, and the
default is that nothing is auto-confirmed at all: the proposal waits for you.

Installed reference schemes are indexed by genus and species from their own
metadata, with a plainly named "Organism not recorded in this reference" group for
the few whose metadata carries no taxon. That index is a pointer view: it does not
re-read alleles and it changes no scheme.

If you get it wrong, fix it: assign the correct genus and species and the managed
copy moves, the old folder is pruned if it is now empty, and the change is
recorded in the project history with what it was before. During import review you
can also work in a spreadsheet: export a template, fill in the genus and species
per file, and load it back.

Identification beyond *Klebsiella* and *E. coli* needs the broader panel:
**Data → Install broader species panel…** downloads 17 pinned RefSeq references
(about 15.7 MB) and verifies every byte. It is not in the ZIP, and until it is
installed the honest verdicts for most genera are genus-only or
`Not_in_reference_panel`.

## What can I do with the right mouse button?

Right-click is the same everywhere: on an isolate table, a folder tree, the
graph, a scheme list or an evidence row. The menu is built from the selection you
actually have, and it tells you the truth about it:

- The count is always shown — "Export 7 isolates", not "Export".
- Right-clicking outside your selection replaces the selection first, so the
  action cannot quietly apply to rows you cannot see.
- Right-clicking blank space offers only the actions that need no selection.
- An action that needs exactly one isolate is **disabled with the reason**, never
  silently narrowed to the first row.
- Destructive actions sit last, after a separator.
- While an analysis is running, everything that would write to the project or the
  disk is disabled.

The generic actions are add, open, rename, assign genus/species, re-file, archive,
remove, copy identifiers, copy file path, open the containing folder, and export
the selection. **Delete is not the default.** Archiving hides an isolate from the
working views while keeping the row, its results, its analyses, its history and
your file exactly as they were; a removal writes the whole record — including
every stored analysis — into the project history first, so it can be restored from
**Samples → Recently removed**. Removing a record never deletes your input file.

An archived isolate disappears from the working tables, comparisons, evidence and
report cohorts, while its row, results, analyses and history stay exactly as they
were. A project-level export still contains it, on purpose: archiving is a view
decision, not a redaction.

## The text is too small, or my screen is very high resolution

**Settings → Display & text size** has three separate controls, because they are
three different problems:

- **Interface text size** (80–150%) changes immediately, while you watch.
- **Whole interface size** scales everything, including icons and spacing. Only
  the sizes your screen can actually show the whole window at are offered; larger
  ones are hidden with the reason. This one takes effect the next time you open
  the application, and says so in a banner that stays put: *"Close and reopen
  WMLSTudio to use this size. Your project and results are not affected."*
- **Graph text size** is saved with your interface preferences and is intended for
  the labels drawn in the tree only, never for exported images. *It is stored but
  not yet applied to the drawn graph in this build;* the setting is recorded and
  the redraw hook is still to land.

A preview row shows a realistic line of text at the chosen size so you can judge
it before committing. An advanced section shows the rounding policy, the resolved
state ("Active display scale: 150% (from your saved preference) · rounding:
exact") and the recovery command — start with `--display-scale 100` if a saved
size ever leaves the window unusable.

## What are the ST and cgMLST profiles?

[Review the typing plan](wmlstudio:analyse). Automatic classical scheme discovery
provides provisional lineage evidence and a registered ST only for a complete,
unambiguous known profile. A novel combination of known alleles is not a novel
allele, and neither should receive an invented public ST number.

Additional cgMLST/wgMLST analyses preserve the classical result. Local novel
sequence identifiers retain sequence hashes; they are not public allele numbers.
Inspect missing, mixed, duplicate, partial and ambiguous calls before comparison.

Use one identical scheme snapshot within each distance calculation. An explicit
database update creates new reference evidence; it must not rewrite an old report.
The download catalog is subject to provider access and terms. cgMLST.org database
contents are not distributed in the portable ZIP.

## What resistance and virulence determinants were found?

[Open linked evidence](wmlstudio:features). Read the method, reference accession,
identity, coverage, coordinates, primary/secondary state and input identity—not
just a list of gene symbols. Different panels assay different questions.

Keep these states distinct:

| State | Meaning |
| --- | --- |
| Detected | A specified assay produced supporting sequence evidence |
| Not detected | That assay ran; no call passed its stated rules |
| Not tested | No relevant assay/panel was run |
| Inconclusive | Quality, ambiguity or incomplete evidence prevents a conclusion |
| Stale | Evidence belongs to an earlier/different assembly or analysis context |

An AMR-only database is not a comprehensive virulence assay. A fragmented locus
is not automatically a functional virulence system. In Klebsiella, siderophore
and capsule-related questions need their appropriate reference definitions and
QC; generic gene detection does not reproduce all Kleborate/Kaptive outputs.

Source: [NCBI interpretation guidance](https://github.com/ncbi/amr/wiki/Interpreting-results).

## Is there a tool specific to my organism?

Three, in this revision. Some questions only make sense for one genus, so they
are offered only when the isolate's organism matches — and they still run, with a
recorded reason, if you ask for them anyway.

- **SCCmec typing**, offered for *Staphylococcus*. It reports which ccr complex
  and mec class were detected and which SCCmec type definitions they satisfy. It
  reports a single type only when exactly one definition is satisfied, nothing is
  unresolved, both complexes resolve, and the required targets sit on one contig.
  Otherwise you get candidates and the reason the call was withheld.
- ***Klebsiella* virulence locus STs** — ybt, clb, iuc, iro, rmp — from exact
  allele vectors, with the published lineage label reported verbatim.
- ***Klebsiella* capsule markers** — a *wzi* and a *wzc* allele number, and
  nothing else.

The plan dialog shows, per module, how many of your selected isolates it is
recommended for, how many are *possible* (right genus, species outside the curated
set — coagulase-negative staphylococci, for instance) and how many are off-panel.
Modules recommended for at least one isolate are pre-selected. The stored evidence
is stamped with which of those it was for that isolate, so a result produced
outside its validation taxa carries that fact into the report and the export.

Three limits matter more than the rest. **mecC is not in the SCCmec panel and is
never assayed** — a mecC-carrying element is not detected by mecA, and the AMR
path's separate mecC result is deliberately not joined into a type call. **No K
locus, capsule type or serotype is inferred** from a *wzi* allele; the published
wzi-to-K-type associations are neither shipped nor applied. And **none of these is
equivalent to Kleborate, Kaptive, staphopia-sccmec or SCCmecFinder**; no
methicillin susceptibility, MRSA/MSSA designation or hypervirulence phenotype is
implied by any of them.

Full detail, every limitation verbatim and the pinned reference revisions:
[organism modules](ORGANISM_MODULES.md). A reference snapshot staged before these
panels existed reports `not_run` naming the missing section — a diagnosable
reference state, never a negative result.

Outside the desktop, the same assays run from the command line:
`WMLSTudio-CLI characterize --list-modules`, then
`characterize <assembly> --module sccmec --organism "Staphylococcus aureus"`.

## Which drugs could these determinants affect?

Use **reference-backed drug/class associations** to prioritize review. Keep
phenotypic AST results in separately labelled metadata with their method,
measurement, interpretive standard/version and date.

Never convert “no gene found” into “susceptible.” Expression, porins, target
mutations, incomplete assemblies, intrinsic mechanisms and unrepresented
determinants can affect phenotype. Do not recommend a treatment from this view.
A genotype–phenotype disagreement is a useful investigation flag, not an excuse
to silently change either result.

Source: [NCBI genotype versus phenotype warning](https://github.com/ncbi/amr/wiki/Interpreting-results).

## How are core and accessory relationships different?

Core-locus allele distances ask about differences across a defined shared
backbone. A whole-genome scheme may contain accessory targets with different
presence and callability. Presence/absence of AMR, virulence or replicon markers
is another question again.

Compare these layers side by side. Do not append a handful of AMR genes to
cgMLST and call the resulting number a validated whole-genome distance. Each
view must identify its target set, missing-data policy and denominator.
Accessory similarity can support a hypothesis but can also reflect common
mobile elements shared by otherwise distinct lineages.

[Choose the cohort and scheme](wmlstudio:compare). Changing graph colors, labels,
layout or displayed isolates does not require repeating allele calling.

## Can plasmids explain or refine this cluster?

A replicon marker suggests plasmid-associated material. AMR/virulence markers on
the same assembled contig provide **co-location evidence**. Neither establishes
that two isolates carry an identical complete plasmid, or which isolate infected
which patient. Short-read contigs may join or split repeated mobile sequences.

Use core relatedness, marker content and contig context as separate evidence
layers. For stronger claims, consider reconstruction with a documented method,
long-read/hybrid closure, coverage and epidemiological context. MOB-suite's
reconstructions are predictions, not a substitute for experimental confirmation.

Source: [MOB-suite methods and scope](https://github.com/phac-nml/mob-suite).

## How do I find and review every cluster?

Select the intended cohort and set a threshold appropriate to the organism,
scheme, missing-data policy and study protocol. There is no universal outbreak
threshold. Include background isolates when their sampling is meaningful.

WMLSTudio carries a dated catalogue of published cutoffs with their citations. It
never applies one for you. Eight organisms have a cutoff bound to a named scheme
and target count, nine more have a published number this catalogue will show you
but refuses to adopt because no scheme is bound, and twelve listed organisms have
no curated cutoff at all. Before the software accepts a number it requires the
exact scheme, the full target set, the reference fingerprint, the caller, the
missing-data policy and your own written justification — and it then records the
result as a local adaptation requiring validation, never as validated. Read
[which organisms are actually covered](THRESHOLDS.md) before assuming yours is
one of them; a seven-locus MLST distance is not the quantity any of those papers
measured.

Threshold groups must use **all comparable pairwise edges**, not the positions
of drawn nodes. Single-linkage allows chains: A may be close to B and B to C,
while A and C are not close. Review the cluster's maximum observed distance and
any unavailable pairwise comparisons. “Not comparable” is different from an
adequately profiled singleton.

Select a graph node or cluster to inspect its members. Give a reviewed group a
name, color and note; automatic threshold membership and the group's reviewed
membership are distinct. Halos are presentation, not additional genetic evidence.

## What has changed since this investigation started?

The Compare tab can show two trees side by side: the **baseline** — the first
snapshot of this investigation, or any earlier snapshot you pin — and the
**current** state. The baseline is a pointer into the investigation's existing
append-only snapshot list, so it costs no extra storage, re-analysis cannot
rewrite it, and pinning a different one does not disturb a comparison already
running.

The baseline tree is a replay of exactly what that snapshot displayed. Distances
are never recomputed for it; you are looking at the stored evidence, not at old
isolates re-measured with today's rules.

Beside the two trees is a change summary, and its honesty rules are the point of
the feature:

- If the two snapshots used a different scheme fingerprint or a different distance
  method, they are **not comparable**, and the distance, edge and cluster sections
  read "not assessed" rather than "no change". Not assessed is never zero.
- A pair whose distance number is unchanged but whose **shared-locus denominator
  moved** is a different measurement, and it is listed as such.
- A pair touching an isolate that was added or removed is counted as not assessed,
  never as unchanged.
- Merges and splits are computed over the isolates present in both snapshots, so a
  deleted isolate can never read as a split.
- If you changed the threshold, the cluster changes are attributed to that. If you
  changed the minimum overlap, the distance changes are attributed to that.

Five statements travel with every comparison and appear in every export of it: an
MST is not a phylogeny or a transmission chain; changing the cohort alone
re-routes edges; the choice between equal-distance edges is arbitrary; a changed
denominator is a different measurement; and a cluster change is not transmission.

Turn the second tree on from the Compare tab. For a large cohort the baseline is
drawn only when you ask for it, so switching investigations stays responsive, and
you can pin any earlier snapshot as the baseline from the snapshot list. Pinning a
baseline does not disturb a comparison that is already running.

## Can the graph show ST, genes and sample names clearly?

Use a concise primary label—usually isolate ID—and add only the fields needed
for the question. Color nodes by one metadata field and reserve halos for groups.
For a crowded 100-isolate graph, select a subgroup or temporarily hide secondary
labels. Zoom and Fit affect the view, not the data or analysis.

When identical profiles are collapsed, the node represents multiple isolates;
inspect all members rather than treating the representative name as the cohort.
Keep a legend and the scheme/threshold/overlap rules in exported figures.

## New samples arrived next week: do I restart?

No. Open the same project and [import the new batch](wmlstudio:import). Review
identity and read associations, and analyse pending/new inputs only. Keep old
profiles if the same input and reference snapshot remain applicable.

Open the saved investigation, explicitly add the new isolates and create a new
comparison snapshot. Unchanged pairwise results can be reused; only new/changed
pairs need recalculation. Compare additions, splits, merges and unresolved
isolates with the previous snapshot. A previously exported report remains a
record of that earlier cohort and method, not a live document.

If you update a scheme or replace an assembly, review the new evidence before
merging it into an investigation. Never compare incompatible snapshots silently.

## How do I report one isolate, selected isolates or clusters?

[Choose the report cohort](wmlstudio:reports). A report should identify who was
included/excluded, the investigation snapshot, scheme, target set, missing-data
rules and cluster threshold. Selected isolates and cluster groups must remain
explicit even if the broader project contains more samples.

A proximity report should show nearest neighbours, allele differences, shared
and total loci, and a clear “not comparable” state. Add the selected identity,
AMR, virulence, plasmid-hypothesis and metadata sections relevant to the audience.
Record limitations and reviewer notes. Do not present genomic associations as
phenotypic AST or confirmed transmission links.

## I just need something short a colleague can read

Choose the **simple summary** layout. It is one short document in plain language
with five parts: who is in it; how close these isolates are, with the tree as an
embedded picture and the threshold together with where that threshold came from;
which resistance genes were found; the closest matches; and what the report does
not tell you.

Every heading says in one sentence what it shows *and* what it does not show. The
susceptibility caveat — genes found are not a measured susceptibility result — is
printed twice, at the top and welded onto the end of the resistance table, so no
option can detach it; it is printed even when the resistance section is switched
off. A distance that cannot be computed prints **"Not comparable"**, and an
isolate with no profile in the comparison gets its own "Not in this comparison"
row: a pair that cannot be compared has no distance, and that is not a distance
of zero. The word cgMLST is never printed unless the scheme's own name contains
it.

Where the threshold came from is stated explicitly: which paper, which scheme it
was measured on, and whether it has actually been bound to this comparison. If
nothing is bound, the report lists the candidate published cutoffs for that
organism and labels them "None of them is applied to this report".

It is not one page. Measured through the real export pipeline on A4: two pages for
two to five isolates, three for ten, four for twenty. No caveat was shortened to
make it fit.

Press **Make a simple summary (PDF)…** on the Reports tab. If no comparison has
been built yet, it offers to build one first rather than printing a document with
an empty proximity section.

## How are CPU and memory managed?

The run review should show an automatic resource plan: detected CPU allowance,
available RAM, concurrent samples and threads per sample. Four simultaneous
four-thread jobs are possible on a 16-thread budget only if memory permits;
assembly and allele calling have different memory requirements.

Balanced settings reserve capacity for the desktop. Faster settings may use
more available resources; low-memory settings admit fewer jobs. A declared RAM
reservation is a planning estimate, not an operating-system-enforced memory cap.
Cancellation stops new admissions and requests cancellation of running work;
completed evidence is retained and incomplete outputs are not labelled complete.

## How do FASTQ pairs attach to an existing FASTA?

Use filename matching to propose candidates, then review the association.
Identical stems, common R1/R2 naming and lane suffixes can help, but names cannot
establish biological identity. Ambiguous matches need manual review.

Attaching reads to an assembly must preserve the assembly, its ID and existing
typing evidence. Validate the complete pair, record both hashes, and keep the
read association separate from the primary assembly. A manually confirmed
filename association does not prove the reads generated that assembly.

SKESA can assemble reviewed paired short reads without a mandatory fastp step.
Read validation and quality review still matter. Built-in sampled QC is labelled
as sampled QC; it must not be described as a completed FastQC analysis.

## Missing targets or disagreement: can I investigate further?

Start with the failed call and assembly context. A target may be genuinely
absent, fragmented, divergent, duplicated or mixed. A second search or read
support can distinguish some of these explanations; relaxing thresholds until
an allele appears does not establish a valid call.

A read-based check must show the mapping method, candidate reference, breadth,
depth, ambiguity and quality rules. Preserve assembly and read evidence as
separate layers; do not silently replace a missing locus with the nearest allele.
Insufficient evidence remains unresolved. Read evidence can support a hypothesis
without justifying an official allele or ST assignment.

For independent sequence-level refinement, SKA2 is split k-mer analysis, not
cgMLST or a recombination-filtered transmission model. Report its k-mer length,
filtering, comparable sequence and SNP/distance definition separately.

Source: [SKA2 implementation and documentation](https://github.com/bacpop/ska.rust).

## How do I clear the view without losing analysis?

Use **Clear selection** to deselect, **Clear filters** to reveal hidden samples,
or change the comparison/report cohort. These operations do not delete results.
Create a new investigation for a different question within the same project.

Removing records is a separate, confirmed action; original sequence files must
not be deleted. Save a project copy before major curation. A clear workspace,
an empty selection, an excluded isolate and a deleted record are different states
and should never share an ambiguous “Reset everything” button.

## What is the platform's added value—and where is expertise still needed?

The value is a connected, auditable investigation: import once; maintain stable
isolate identities; reuse versioned results; move directly between a question,
its table, graph and selected report; preserve the earlier state as data grow.

Automation should remove file handling and resource-management burden, not hide
scientific uncertainty. Expert review remains important for species conflicts,
contamination, novel/ambiguous calls, unusual resistance, threshold selection,
mobile-element reconstruction and clinical interpretation. No test count alone
establishes “top accuracy” across organisms, laboratories and sequencing methods.
