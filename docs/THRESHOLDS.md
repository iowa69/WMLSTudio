# Cluster thresholds: what is actually curated, and for which organisms

Catalogue `2026-09-12.1`, reviewed 2026-09-12. Source of truth:
`src/wmlstudio/threshold_guidance.py`. This page is an audit of that file, not a
second copy of it; if the two disagree, the code is right and this page is stale.

**The one rule.** WMLSTudio never applies a published cutoff for you. Every
catalogue entry carries `auto_apply: false`. A number becomes an approved local
threshold only after you bind the exact published scheme, match its full target
count, record the reference fingerprint, caller and missing-data policy, and
write your own justification. Until then the number is a citation you can read,
not a setting the software is using.

A genomic cluster is a signal for epidemiological review. It is not proof of
transmission, of the direction of transmission, or of clinical causation. Single
linkage joins distant endpoints through intermediates. Missing calls are unknown,
not identical. A published cutoff is not validated for WMLSTudio merely because
the organism name or the locus count matches.

## Honest coverage, in one line each

The catalogue lists **31 organisms** and holds **34 entries**
drawn from **13 cited sources**, covering **20 of the 31** with at least one entry.

| Tier | Organisms | What you get |
| --- | --- | --- |
| A — a published cutoff bound to a named scheme | **13** | A number the software will let you adopt, after review |
| B — a published number with no scheme this catalogue can bind | **7** | A citation to read. The software refuses to adopt the number |
| C — listed for surveillance, no curated publication | **11** | An explicit evidence gap, with no number at all |

Tier A is the honest answer to "thresholds pre-set for common human and
paediatric pathogens with references": **thirteen organisms, seventeen entries,
eleven papers.** Six organisms were added by primary-source review after the
first audit reported eight; each was independently re-checked against the paper
before it was admitted, and four candidates were **rejected** in the same pass
(see "What was searched for and refused"). No number was added that a source
does not state.

## Tier A — published cutoff, bound to a named scheme

These eleven entries can supply an approved threshold, because each names the
exact scheme its number was measured on and the full target count it requires.

| Organism | Cutoff | Scheme it is bound to | Loci | Source |
| --- | --- | --- | --- | --- |
| *Acinetobacter baumannii* | ≤ 9 allele differences | `cgmlst.org:abaumannii-2390` | 2,390 | Glasgow 2025 |
| *Citrobacter freundii* | ≤ 10 allele differences | `kieninger2025:cfreundii-3250` | 3,250 | Kieninger 2025 |
| *Enterococcus faecalis* | ≤ 7 allele differences | `cgmlst.org:efaecalis-1972` | 1,972 | Glasgow 2025 |
| *Enterococcus faecium* | ≤ 20 allele differences | `cgmlst.org:efaecium-1423` | 1,423 | Glasgow 2025 |
| *Enterococcus faecium* | ≤ 20 allele differences | `cgmlst.org:efaecium-1423` | 1,423 | de Been 2015 |
| *Enterococcus faecium* | ≤ 25 allele differences | `cgmlst.org:efaecium-1423` | 1,423 | Higgs 2022 (population grouping) |
| *Enterococcus faecium* | ≤ 7 SNPs | `higgs2022:ska-short-reads-k15` | — | Higgs 2022 (SKA k=15 protocol) |
| *Escherichia coli* | ≤ 10 allele differences | `cgmlst.org:ecoli-2513` | 2,513 | Glasgow 2025 |
| *Klebsiella pneumoniae* | ≤ 15 allele differences | `cgmlst.org:kpneumoniae-2358` | 2,358 | Glasgow 2025 |
| *Serratia marcescens* | ≤ 12 allele differences | `cgmlst.org:smarcescens-2692` | 2,692 | Kampmeier 2022 |
| *Staphylococcus aureus* | ≤ 24 allele differences | `cgmlst.org:saureus-1861` | 1,861 | Glasgow 2025 |
| *Listeria monocytogenes* | ≤ 7 allelic mismatches | `pasteur:lmonocytogenes-1748` | 1,748 | Moura 2016 |
| *Listeria monocytogenes* | ≤ 7 allele differences | `cgmlst.org:lmonocytogenes-1701` | 1,701 | Van Walle 2018 (ECDC evaluation) |
| *Clostridioides difficile* | ≤ 6 allele differences | `cgmlst.org:cdifficile-2270` | 2,270 | Bletz 2018 |
| *Mycobacterium tuberculosis* complex | ≤ 5 allele differences | `cgmlst.org:mtbc-2891` | 2,891 | Kohl 2018 |
| *Salmonella enterica* | ≤ 10 allele differences | `enterobase:senterica-cgmlst-3002` | 3,002 | Leeper 2023 (PulseNet USA) |
| *Legionella pneumophila* | ≤ 4 allele differences | `cgmlst.org:lpneumophila-1521` | 1,521 | Moran-Gilad 2015 (preliminary) |

Read the caveats that travel with each entry:

- **Glasgow 2025** compared SeqSphere+ public-scheme pipelines on a hospital
  cohort **selected on prior clustering**. Pipeline agreement is not independent
  proof of transmission, and the exact snapshot and caller settings still need
  review.
- ***E. faecium* has four different numbers on purpose.** 20 (Glasgow, de Been),
  25 (Higgs, a *population grouping* cutoff used before fine-scale analysis) and
  7 SNPs (Higgs, the original SKA short-read protocol at k=15). They answer
  different questions. Picking whichever is most convenient is exactly the
  mistake this catalogue exists to prevent.
- **The 7-SNP entry is not a SKA2 assembly cutoff.** Its own limitation text
  says it is "not validated for SKA2 assembly distance, reference-mapped SNPs,
  coding-only SNPs or changed masks. No automatic transfer."
- **Kampmeier 2022** specifies pairwise-ignore for missing targets and describes
  the cutoff as prompting epidemiological investigation, not diagnosing
  transmission.
- **Kieninger 2025** is the *species-specific* 3,250-target scheme and its
  maximum intracluster distance — **not** the paper's separate combined-species
  scheme or its 8-allele criterion. Single-linkage endpoints can exceed the
  published within-cluster bound; read cluster diameters.

## Tier B — a number you can read, not one the software will adopt

Seventeen entries, covering **seven organisms that have nothing better**. The
other nine organisms in this table also appear in tier A, from a different
laboratory — see the disagreement note below. Sixteen of the seventeen come from
one source: **Siddall, Starkey and Patel (2025)**, a Mayo Clinic *local*
related-isolate criterion run on
SeqSphere+ 10.0.5 with SKESA 2.3.0. The paper's exact scheme identity is not
curated here, so the catalogue records the number, the context and the refusal
together. A local validation criterion is not transferable by species name.

| Organism | Published number | Why it cannot be adopted here |
| --- | --- | --- |
| *Acinetobacter baumannii* | ≤ 9 | Local criterion; no scheme bound |
| *Clostridioides difficile* | ≤ 6 | Local criterion; no scheme bound |
| *Cutibacterium acnes* | ≤ 5 | Local criterion; no scheme bound |
| *Enterobacter cloacae* complex | ≤ 15 | Local criterion; no scheme bound |
| *Enterococcus faecalis* | ≤ 7 | Local criterion; no scheme bound |
| *Enterococcus faecium* | ≤ 7 | Local criterion; no scheme bound |
| *Escherichia coli* | ≤ 10 | Local criterion; no scheme bound |
| *Klebsiella pneumoniae* | ≤ 15 | Local criterion; no scheme bound |
| *Legionella pneumophila* | ≤ 4 | Local criterion; no scheme bound |
| *Pseudomonas aeruginosa* | ≤ 6 | Local criterion; no scheme bound |
| *Serratia marcescens* | ≤ 12 | Local criterion; no scheme bound |
| *Staphylococcus aureus* | ≤ 8 | Local criterion; no scheme bound |
| *Staphylococcus epidermidis* | ≤ 8 | Local criterion; no scheme bound |
| *Staphylococcus lugdunensis* | ≤ 8 | Local criterion; no scheme bound |
| *Streptococcus agalactiae* | ≤ 20 | Local criterion; no scheme bound |
| *Streptococcus pyogenes* | ≤ 20 | Local criterion; no scheme bound |
| *Pseudomonas aeruginosa* | **none** | Tönnies 2021 defines the 3,867-locus scheme but this catalogue curates **no universal numeric rule** from it |

Notice that five organisms appear in both tiers with **different numbers** —
*S. aureus* 24 versus 8, *E. faecium* 20/25 versus 7, *P. aeruginosa* a scheme
with no number versus a local 6. Two laboratories reached different answers for
the same species. That disagreement is the evidence; hiding it behind one
"recommended" value would be the dishonest option.

*P. aeruginosa* deserves its own sentence. Tönnies 2021 establishes the
3,867-locus public scheme, and the catalogue deliberately records a `null`
threshold for it with the note: review lineage, hypermutation, environmental
persistence and the exact protocol. Do not substitute the 2025 comparison's ad
hoc scheme for it.

## Tier C — listed, but no curated cutoff at all

These eleven are on the surveillance work-list and have **no entry**. That is an
evidence gap in this catalogue, not proof that no publication exists.

| Organism | Catalogue status |
| --- | --- |
| Klebsiella variicola | no curated cutoff |
| Klebsiella quasipneumoniae | no curated cutoff |
| Klebsiella oxytoca complex | no curated cutoff |
| Klebsiella aerogenes | no curated cutoff |
| Stenotrophomonas maltophilia | no curated cutoff |
| Burkholderia cepacia complex | no curated cutoff |
| Proteus mirabilis | no curated cutoff |
| Morganella morganii | no curated cutoff |
| Providencia stuartii | no curated cutoff |
| Staphylococcus capitis | no curated cutoff |
| Streptococcus pneumoniae | no curated cutoff — searched, none found |

`guidance_for()` returns status `no_curated_transferable_cutoff` for each, with
the message "No reviewed transferable cutoff in this catalog. This is an evidence
gap, not proof that no publications exist."

## The paediatric question, answered directly

Paediatric and neonatal work is covered **unevenly**, and the uncovered names are
not minor ones.

Covered with a bindable cutoff (tier A): *S. aureus*, *E. coli*,
*K. pneumoniae*, *E. faecium*, *E. faecalis*, *S. marcescens* (a classic NICU
outbreak organism), *A. baumannii*, *C. freundii*, and now *L. monocytogenes*
(maternal-neonatal listeriosis), *S. enterica*, *C. difficile*, *M. tuberculosis*
complex and *L. pneumophila*.

Covered by citation only (tier B): *S. agalactiae* — the leading cause of
early-onset neonatal sepsis — *S. pyogenes*, *S. epidermidis* and
*S. lugdunensis* (CoNS device infections), *E. cloacae* complex.

**Not covered at all:** *Streptococcus pneumoniae*, *Staphylococcus capitis*
(the NICU-adapted NRCS-A clone), *Neisseria meningitidis*, *Haemophilus
influenzae*, *Campylobacter jejuni*. The last three are not even on the
work-list.

*S. agalactiae* is the most painful of these, and it was searched specifically.
The one verifiable published cutoff is SNP-based and computed **separately per
clonal complex** against a different reference genome for each CC, so there is no
single cgMLST scheme and no single locus count to bind it to. It stays in tier B
rather than being dressed up as a cgMLST threshold.

Several of those organisms are in the `mixed-genus-20` practice cohort and in
the broad species panel. WMLSTudio will happily identify, file, type and compare
them — and will then tell you, correctly, that it has **no threshold to offer**.
That is the intended behaviour, and it is why the guidance panel prints an
evidence-gap message instead of a number.

## What was searched for and refused

The tier A additions came from a primary-source review that also produced four
refusals. They are recorded because a documented negative is the useful result
here: it stops the same ground being re-searched, and it stops a weak number
being adopted later by someone who assumes nobody looked.

| Organism | Why no entry |
| --- | --- |
| *Streptococcus pneumoniae* | No cgMLST scheme with a published outbreak cutoff was found. cgmlst.org publishes no pneumococcal scheme; PubMLST's Pneumococcal Genome Library publishes cgLIN **hierarchical barcoding levels**, which are a clustering construct, not a validated outbreak threshold |
| *Streptococcus agalactiae* | The only verifiable cutoff is SNP-based and computed per clonal complex against per-CC references — not bindable to one scheme and one locus count |
| *Enterobacter cloacae* complex | No cgMLST scheme exists for the complex as a taxon. cgmlst.org hosts only a single-species *E. hormaechei* scheme; PubMLST hosts only 7-gene MLST |
| *Salmonella enterica*, European practice | ECDC and EFSA set the cutoff **per outbreak** (5 and 7 allele differences in two different 2023–2024 assessments, on a 3,255-locus chewBBACA scheme, not EnteroBase's 3,002). There is no single European number to curate, so only the PulseNet USA criterion is recorded |

**EnteroBase HC5 is deliberately not catalogued.** It is a single-linkage
hierarchical clustering level, and its own authors write that such clusters "may
not necessarily represent traditional transmission chains", recommending that a
transmission chain "continue to be based on epidemiological criteria in addition
to genetic similarities". Recording HC5 as a threshold would misrepresent it.

## Nothing in the ZIP can bind a tier A threshold today

All 162 bundled reference schemes are classical PubMLST MLST schemes of seven to
ten loci. Every tier A cutoff except the SKA SNP entry is bound to a
cgMLST.org (or study-specific) scheme of 1,423 to 3,250 targets, and **cgMLST.org
database contents are not redistributed in the portable build** — its
[server policy](https://www.cgmlst.org/serverpolicy.html) restricts that. So out
of the box, no published cgMLST cutoff can be adopted: you must install the exact
scheme yourself, under its own terms, before the threshold dialog will accept the
number. A seven-locus MLST distance is not the quantity any of these papers
measured.

## What the software checks before it accepts a number

`threshold_guidance.record_decision()` refuses, with a specific message, unless
every one of these holds. It is worth reading as the definition of "pre-set"
in this application:

1. The entry actually carries a numeric cutoff (Tönnies 2021 does not).
2. Your comparison method matches the entry's method — cgMLST alleles for an
   allele entry, SNPs for a SNP entry.
3. The organism matches exactly. There is no genus fallback and no taxonomic
   inference; `guidance_for("Klebsiella variicola")` does not answer with the
   *K. pneumoniae* entry.
4. The **exact scheme key is bound**. An organism name or a locus count alone is
   refused — this is what stops all seventeen tier B entries.
5. The full published target set is present. There is no automatic scaling for
   missing or accessory loci.
6. The reference fingerprint, the caller and the missing-data policy are all
   recorded.
7. You have attested that you reviewed the scheme, the protocol and the
   epidemiology, and written a justification of at least twenty characters.

A decision that passes is stamped `local_adaptation_requires_validation`, never
"validated", and it records whether your number departs from the published one.
A citation recorded without a threshold is stamped
`citation_only_no_threshold_change`.

## Review scope and staleness

This is a targeted primary-source review including 2025 literature. It is **not**
systematic and **not** continuously updated. `review_age_days()` reports how old
it is; treat a stale catalogue as a prompt to re-read the primary sources, not as
a reason to keep using an ageing number.

Adding an organism to this catalogue means reading the paper, extracting the
exact scheme identity and target count, and recording the study's own scope and
missing-data policy. It does not mean copying a number from a review table. If
you cannot state which scheme a cutoff was measured on, the honest entry is a
`null` threshold with the citation — the shape Tönnies 2021 already has.

## Full citations

| Key | Citation | DOI | Type |
| --- | --- | --- | --- |
| `glasgow2025` | Glasgow et al. (2025). Comparison of core genome multi-locus sequencing typing pipelines for hospital outbreak detection of common bacterial pathogens. | [10.1128/jcm.00646-25](https://doi.org/10.1128/jcm.00646-25) | Pipeline comparison |
| `siddall2025` | Siddall, Starkey and Patel (2025). Automated whole genome sequencing platform for bacterial strain typing in clinical microbiology laboratories. | [10.1128/jcm.00178-25](https://doi.org/10.1128/jcm.00178-25) | Local validation |
| `kampmeier2022` | Kampmeier et al. (2022). Development and Evaluation of a Core Genome Multilocus Sequencing Typing (cgMLST) Scheme for *Serratia marcescens* Molecular Surveillance and Outbreak Investigations. | [10.1128/jcm.01196-22](https://doi.org/10.1128/jcm.01196-22) | Scheme evaluation |
| `debeen2015` | de Been et al. (2015). Core Genome Multilocus Sequence Typing Scheme for High-Resolution Typing of *Enterococcus faecium*. | [10.1128/jcm.01946-15](https://doi.org/10.1128/jcm.01946-15) | Scheme evaluation |
| `higgs2022` | Higgs et al. (2022). Optimising genomic approaches for identifying vancomycin-resistant *Enterococcus faecium* transmission in healthcare settings. | [10.1038/s41467-022-28156-4](https://doi.org/10.1038/s41467-022-28156-4) | Genomic epidemiology |
| `tonnies2021` | Tönnies et al. (2021). Establishment and Evaluation of a Core Genome Multilocus Sequence Typing Scheme for Whole-Genome Sequence-Based Typing of *Pseudomonas aeruginosa*. | [10.1128/jcm.01987-20](https://doi.org/10.1128/jcm.01987-20) | Scheme evaluation |
| `citrobacter2025` | Kieninger et al. (2025). Development and validation of a core genome multilocus sequence typing scheme for *Citrobacter freundii*. | [10.1128/jcm.00860-25](https://doi.org/10.1128/jcm.00860-25) | Scheme evaluation |

Each entry also records the locator within its paper (table, figure or methods
section) so a reviewer can find the number without re-reading the whole study.
