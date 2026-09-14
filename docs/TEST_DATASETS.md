# Practice cohorts: pinned public test datasets

Two public datasets are provided so you can watch WMLSTudio propose a genus and a
species for every file you import and file the copies into the folders those
proposals imply. They are **teaching cohorts, not validation sets**.

No sequence data is committed to this repository or shipped in the portable ZIP.
What is committed is a table of NCBI RefSeq accessions with the exact byte count,
MD5 and SHA-256 of each file (`src/wmlstudio/practice_cohorts.py`). The sequences
are downloaded only when you ask for them, and every byte is checked against those
pinned values plus the MD5 that NCBI publishes today. A mismatch aborts the whole
download and nothing partial is written.

---

## How to fetch a cohort

```
.venv/bin/python studio_scripts/fetch_practice_cohort.py --list
.venv/bin/python studio_scripts/fetch_practice_cohort.py \
    --cohort mixed-genus-20 --destination ~/WMLSTudio-practice/mixed-genus-20
```

| Option | Effect |
| --- | --- |
| `--list` | Describe both cohorts, their genera, their size and their caveats |
| `--cohort NAME --destination DIR` | Create `DIR` and install the verified cohort into it |
| `--cache-root DIR` | Shared content-addressed cache; defaults to `_cache` beside the cohort |
| `--verify DIR` | Re-hash an installed cohort against the pinned table |
| `--verify-pins --cohort NAME` | Compare the pinned MD5s with NCBI's published checksums without downloading any genome; exit status 1 when a pin has gone stale |
| `--source-root DIR` | Read from a local mirror of the NCBI `genomes/all/` tree instead of the network |

Files are written **flat and unsorted**, named
`<Genus>_<species>_<strain>_<accession>.fna.gz`. That is deliberate: creating the
organism folders is the thing being demonstrated, so the download must not do it
for you. A `manifest.json` beside them records every hash, every source URL, the
fetch timestamp and the caveats below.

The download directory is git-ignored (`practice-cohorts/` plus the repository-wide
`*.fna.gz` rule), so a cohort fetched inside a checkout can never be committed.

Re-running a download against an already installed cohort re-verifies it and
downloads nothing. An accession that appears in both cohorts (`GCF_000240185.1`)
is fetched once into the shared cache and copied into each cohort.

---

## Cohort A — `kpneumoniae-10`

Ten *Klebsiella pneumoniae* complete genomes. 16 746 561 bytes compressed.

Purpose: single-species filing, MLST/cgMLST comparison and the minimum spanning
forest on one organism.

| # | Organism recorded by NCBI | Strain | Accession |
| --- | --- | --- | --- |
| 1 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH1 | GCF_000281535.2 |
| 2 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH10 | GCF_000281435.2 |
| 3 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH24 | GCF_000714675.1 |
| 4 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH27 | GCF_000695935.1 |
| 5 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH29 | GCF_000784945.1 |
| 6 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH30 | GCF_000784985.1 |
| 7 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH31 | GCF_000785005.1 |
| 8 | *Klebsiella pneumoniae* subsp. *pneumoniae* | KPNIH33 | GCF_000775375.1 |
| 9 | *Klebsiella pneumoniae* subsp. *pneumoniae* | NTUH-K2044 | GCF_000009885.1 |
| 10 | *Klebsiella pneumoniae* subsp. *pneumoniae* | HS11286 | GCF_000240185.1 |

Rows 1–8 are the NIH Clinical Center *K. pneumoniae* complete-genome series
described by Snitkin *et al.*, *Sci Transl Med* 2012
(doi:10.1126/scitranslmed.3004129) and Conlan *et al.*, *Sci Transl Med* 2014
(doi:10.1126/scitranslmed.3009845). Those papers are **background for the cohort,
not an expected result for this software**. Rows 9–10 are widely used unrelated
reference strains, included as visible non-cluster controls.

### What to expect on import

Every isolate should be proposed as one organism, so the whole cohort should land
in a single genus/species folder. After typing, the managed copies re-file into
`ST_*` subfolders. Isolates that fail to type, or that the app declines to call,
stay visible as such — an unassigned ST is not a typing result and must not be read
as one.

---

## Cohort B — `mixed-genus-20`

Twenty complete genomes across **15 genera and 19 species**. 21 703 421 bytes
compressed.

Purpose: demonstrate automatic genus/species proposals and the folders they create
— *including the calls the app should refuse to make*.

| # | Organism recorded by NCBI | Strain | Accession | Role in the demo |
| --- | --- | --- | --- | --- |
| 1 | *Klebsiella pneumoniae* subsp. *pneumoniae* | HS11286 | GCF_000240185.1 | Also in cohort A; fetched once |
| 2 | *Klebsiella variicola* subsp. *variicola* | F2R9T | GCF_020525545.1 | **Species-complex member — a hard call** |
| 3 | *Escherichia coli* O157:H7 | Sakai | GCF_000008865.2 | **ANI cannot separate it from *Shigella*** |
| 4 | *Staphylococcus aureus* subsp. *aureus* | NCTC 8325 | GCF_000013425.1 | Two isolates share one species folder |
| 5 | *Staphylococcus aureus* subsp. *aureus* | USA300 FPR3757 | GCF_000013465.1 | Two isolates share one species folder |
| 6 | *Enterococcus faecium* | DO | GCF_000174395.2 | Two species in one genus |
| 7 | *Enterococcus faecalis* | V583 | GCF_000007785.1 | Two species in one genus |
| 8 | *Acinetobacter baumannii* | ATCC 17978 | GCF_000015425.1 | |
| 9 | *Pseudomonas aeruginosa* | PAO1 | GCF_000006765.1 | |
| 10 | *Salmonella enterica* subsp. *enterica* ser. Typhimurium | LT2 | GCF_000006945.2 | |
| 11 | *Streptococcus pneumoniae* | TIGR4 | GCF_000006885.1 | Three species in one genus |
| 12 | *Streptococcus agalactiae* | NEM316 | GCF_000196055.1 | Three species in one genus |
| 13 | *Streptococcus pyogenes* M1 GAS | SF370 | GCF_000006785.2 | Three species in one genus |
| 14 | *Listeria monocytogenes* | EGD-e | GCF_000196035.1 | |
| 15 | *Neisseria meningitidis* | MC58 | GCF_000008805.1 | |
| 16 | *Haemophilus influenzae* | Rd KW20 | GCF_000027305.1 | |
| 17 | *Clostridioides difficile* | 630 | GCF_000009205.2 | |
| 18 | *Campylobacter jejuni* subsp. *jejuni* | NCTC 11168 | GCF_000009085.1 | |
| 19 | *Enterobacter cloacae* subsp. *cloacae* | ATCC 13047 | GCF_000025565.1 | **No species reference in the bundled panel** |
| 20 | *Serratia marcescens* subsp. *marcescens* | Db11 | GCF_000513215.1 | **No species reference in the bundled panel** |

None of these twenty assemblies is one of the assemblies in the bundled
*Klebsiella* characterization panel, so no comparison in this cohort is a
self-match of an isolate against its own reference.

### What to expect on import — including the honest failures

This cohort is built so that the uncomfortable outcomes are visible rather than
hidden. What you should look for:

- **A confident species proposal** for most isolates: the whole-genome ANI
  evidence has a nearby reference of the same species, the competing species are
  clearly further away, and enough of both genomes aligned. These file into
  `<Genus>/<species>/`.
- **A complex-level answer for *Klebsiella variicola* (row 2).** *K. variicola*
  belongs to the *K. pneumoniae* species complex. Separating complex members by
  whole-genome ANI against a single representative per taxon is a genuinely hard
  call, and the honest outcome may be the complex rather than the species. If the
  app offers `Klebsiella/variicola/`, that is a proposal to review, not a
  confirmed identification; if it declines to choose, that is the correct
  behaviour, not a defect.
- **A complex-level answer for *Escherichia coli* (row 3).** ANI cannot reliably
  distinguish *Escherichia coli* from *Shigella*. The app reports the
  *E. coli*/*Shigella* complex and keeps that verdict below the confidence needed
  to file without a human decision.
- **A genus-only or unresolved outcome for *Enterobacter cloacae* and *Serratia
  marcescens* (rows 19–20).** Neither species has a genomic reference in the
  bundled panel. The honest result is a genus-level label with no species, or no
  organism at all, and the file waits for you in the needs-review area instead of
  being filed under a guessed species. Installing a broader species panel changes
  this; nothing about the isolates changes.
- **An unresolved outcome is a valid outcome.** A file the app cannot place is
  imported and kept, not discarded and not renamed, and its original is never
  touched.

### What running this cohort actually produced

All twenty genomes were run against the staged species panel during development.
Seventeen were proposed with genomic reference support; three were not, and each
of those three is correct:

| Isolate | Outcome | Why |
| --- | --- | --- |
| *Klebsiella variicola* | unresolved, 94.54% ANI | Species-complex member with no reference of its own. It was **not** called *K. pneumoniae*, which is the point of including it |
| *Enterobacter cloacae* | unresolved | No species reference in the panel |
| *Serratia marcescens* | unresolved | No species reference in the panel |

*Escherichia coli* was proposed at 98.08% ANI but held at complex level, as
described above.

That run also found a real gap. *Listeria monocytogenes* EGD-e first came back
**unresolved at 94.83% ANI**, because the panel held only a lineage I reference
and the two *L. monocytogenes* lineages sit either side of the 95% ANI species
line. The software was right to refuse — but for a listeriosis investigation it
was useless. The fix was to add a second, lineage II reference to the panel
rather than to lower the ANI gate, which would have traded a visible refusal for
an invisible wrong answer. EGD-e now resolves at 99.03%. A species whose lineages
are that divergent may be under-covered by one reference; that limitation is
recorded in the panel's own manifest.

Anything that reaches you as a folder name is a **filing decision**. It is not a
laboratory identification, and it is not a susceptibility or virulence statement.
Confirm the organism before any clinical interpretation.

---

## Caveats shipped with both cohorts

These are recorded verbatim in each cohort's `manifest.json` and are surfaced by
the app wherever the cohort is offered:

- These are public complete genomes from NCBI RefSeq. They are a teaching cohort,
  not a validation set.
- No expected ST, cluster, threshold or "correct answer" is shipped with this
  cohort. WMLSTudio computes relatedness from the scheme you choose; agreement
  with any published investigation is not a validation of this software.
- Strain and organism names are the labels NCBI records for these assemblies. They
  were not independently verified here.
- These are finished reference assemblies. Draft assemblies from your own
  sequencing will have more missing loci and different comparison behaviour.
- A genus or species folder created on import is a filing decision, not a
  laboratory identification. Confirm the organism before any clinical
  interpretation.

`mixed-genus-20` adds:

- *Klebsiella variicola* belongs to the *Klebsiella pneumoniae* species complex.
  Separating it from *K. pneumoniae* is a genuinely hard call; treat a
  complex-level answer as the honest one unless you have independent evidence.

---

## Provenance, licence and attribution

**Source.** Every file is the `*_genomic.fna.gz` published for that assembly under
`https://ftp.ncbi.nlm.nih.gov/genomes/all/` — the same path the app builds from the
accession, for example
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz`.

**Verification performed for these pins.** Each accession's FTP directory was
listed, its `*_assembly_report.txt` read (all 30 distinct assemblies report
`Assembly level: Complete Genome`, and the organism and strain names in the tables
above are quoted from those reports), its `md5checksums.txt` fetched, and the
genomic FASTA downloaded in full. The published MD5 matched the downloaded bytes
for all 30. The byte count, MD5, SHA-256 of the `.gz` and SHA-256 of the
decompressed FASTA were computed from those bytes and are the values pinned in
`src/wmlstudio/practice_cohorts.py`. The decompressed hash is pinned as well
because a gzip container can be regenerated upstream while the sequence stays
identical.

**Licence and attribution.** Assemblies are public NCBI RefSeq records. NCBI places
no restrictions on the use or distribution of the data it hosts, but it does not
hold the copyright of those records and cannot grant rights on behalf of the
depositing submitters; individual submitters may assert terms. Cite the assembly
accession and the originating submitters — not WMLSTudio — when you reuse these
sequences. WMLSTudio redistributes none of them: it only records their accessions
and checksums and fetches them from NCBI on your explicit request.

**If a pin goes stale.** NCBI occasionally re-releases an assembly. When that
happens the download stops with the accession named, rather than quietly accepting
different bytes. `--verify-pins` detects the same condition cheaply, without
downloading any genome, so a curator can refresh the table deliberately.
