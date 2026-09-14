# WMLSTudio: from a problem to reviewable evidence

This handbook is available offline inside **Help → Problem → solution guide**.
Start with a research question, keep one stable record per isolate, and reuse
the same evidence in tables, graphs and reports. A project is the evidence store;
an investigation is a named, explicitly selected comparison within that store.

The application supports non-commercial research. It is not a validated
diagnostic device. Genetic proximity does not establish a transmission event;
genomic resistance determinants are not measured antimicrobial susceptibility.

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

Threshold groups must use **all comparable pairwise edges**, not the positions
of drawn nodes. Single-linkage allows chains: A may be close to B and B to C,
while A and C are not close. Review the cluster's maximum observed distance and
any unavailable pairwise comparisons. “Not comparable” is different from an
adequately profiled singleton.

Select a graph node or cluster to inspect its members. Give a reviewed group a
name, color and note; automatic threshold membership and the group's reviewed
membership are distinct. Halos are presentation, not additional genetic evidence.

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
