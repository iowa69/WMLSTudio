# cgMLST schemes: what ships, what does not, and why

Catalogue version `2026-09-14.1`. Provider terms read 2026-09-14. Every target list
below was fetched from its provider on that date and re-verified with
`studio_scripts/stage_cgmlst_schemes.py --verify` (22 of 22 pins matched).

A cgMLST scheme is a **target set** plus an allele nomenclature. Two profiles are
comparable only when both were called against the same installed scheme. A cgMLST
allele distance and a seven-locus MLST allele distance are different quantities:
they never share a scale, a column, an axis or a threshold.

---

## 1. The download failure the user reported, and what it actually was

`reference_catalog.py` already downloaded cgMLST.org schemes. Exercised against the
live services on 2026-09-14, these are the faults found. All of them are fixed.

### 1.1 One unreadable character destroyed a whole 42 MB download

`https://www.cgmlst.org/ncs/schema/Lmonocytogenes/alleles/` returns 42.7 MB
containing 1,701 locus FASTA files. Ten of those files — `lmo0001`, `lmo0002`,
`lmo0003`, `lmo0004`, `lmo0005`, `lmo0006`, `lmo0007`, `lmo0008`, `lmo0009`,
`lmo0012` — each contain exactly one allele carrying a literal `X` where a base is
unknown. `X` is not an IUPAC nucleotide code, so the application's own sequence
reader rejected it. The download completed, the archive was extracted and checked,
and then validation failed with:

```
SchemeError: lmo0001.fasta, line 1778: invalid DNA characters: 'X'.
```

The entire snapshot was deleted. What the user saw, after a long download, was
`Could not complete: lmo0001.fasta, line 1778: invalid DNA characters: 'X'.` — a
message no microbiologist can act on.

**Fix.** Allele *records* the reader cannot accept are excluded, never the locus and
never the scheme. They are **not** rewritten to `N`: that would invent a base the
provider did not report. Each exclusion is named in the snapshot manifest
(`excluded_alleles`), in the result, and in a note shown to the user, because the
consequence is real — a genome carrying `lmo0001` allele 597 is now reported as an
unmatched sequence, not as allele 597. A locus whose *every* allele is unreadable
still fails the download, loudly.

Verified after the fix: **L. monocytogenes 1,701 targets install in 111 s**, with
10 allele records excluded and named.

### 1.2 The other faults found

| Fault | Evidence | Fix |
|---|---|---|
| Typing errors escaped the documented exception type | `typing.SchemeError` / `sequence.SequenceError` propagated out of `download_scheme`, which is documented to raise `CatalogError` | wrapped in `CatalogError` with a plain sentence |
| A wrong scheme id looked like a corrupt table | cgMLST.org answers an unknown slug with **HTTP 200** and an HTML body `ERROR Occured … Illegal Typing ID!`, so the parser reported "did not return its expected locus table" | the error page is detected and the identifier named |
| Rate limiting was not waited out | PubMLST answers a burst with **HTTP 429**; the old backoff ignored `Retry-After` and waited 0.5 s | `Retry-After` honoured (bounded at 120 s), and a final 429 is reported as "the reference service is rate-limiting this computer" |
| Truncated chunked bodies raised an unhandled type | none of these responses carries `Content-Length`, so the truncation check never fires and a short read surfaces as `http.client.IncompleteRead`, which was not caught | caught and retried like any transport error |
| The progress bar froze for minutes | chunked transfer, no declared length, no per-byte progress: the 42 MB Listeria archive and the 211 s S. aureus archive both showed no movement at all | byte progress is reported every 4 MB |
| A failed download threw away everything | a PubMLST/Pasteur cgMLST scheme is **one request per target**: Listeria cgMLST1748 from Institut Pasteur is 1,748 requests, ~2.5 GB, ~17 minutes (measured). One failure at target 1,700 restarted from zero | resume: fetched targets are kept in a hidden `.resume-*` folder with a SHA-256 ledger and reused. A partial folder is never a scheme, is never published, and is discarded automatically if the remote scheme changed |
| No warning of cost before starting | none | `download_estimate(entry)` states the request count, whether it resumes, and that the installed scheme is gigabytes |

**Measured sizes** (2026-09-14, installed on disk): cgMLST.org *M. gallisepticum*
425 targets = 15 MB; *E. faecium* 1,423 = 818 MB; *S. aureus* 1,861 = **4.6 GB**.
PubMLST *S. aureus* cgMLST 1,716 targets ≈ 0.87 GB fetched as 1,716 requests over
~24 minutes.

---

## 2. Licence review: who permits redistribution

Each provider's terms were read directly on 2026-09-14 and the verdict recorded in
`cgmlst_schemes.PROVIDERS` with the sentence it rests on.

| Provider | Stated terms | May we pack it into the release? |
|---|---|---|
| **PubMLST**, University of Oxford | "All data submitted on or before 31 December 2024 may continue to be downloaded, used, and redistributed without restriction, subject only to proper scientific citation and acknowledgement." Data submitted from 2025-01-01 may not be redistributed or incorporated into a software product without written authorisation. | **YES**, for the pre-2025 subset — which is exactly what unauthenticated API access serves ("you are currently restricted to accessing data that was submitted on or prior to 2024-12-31"). |
| **Institut Pasteur BIGSdb** | "You acknowledge and agree that you may not share or publicly post significant parts of the Data (including alleles and profiles definitions)." Academic use only; commercial tools require a licence contract. | **NO** |
| **cgMLST.org / Ridom GmbH** | "Reuse of database copies in a product or service requires permission." | **NO** |
| **EnteroBase**, University of Warwick | Prohibits "Reverse engineering and Reproducing EnteroBase-like database copies without explicit licencing"; academic use only; redistribution allowed only with attribution and statement of modifications. | **NO** — packing the allele database into a product is reproducing a database copy. |
| **Chewie-NS** (chewbbaca.online) | The *software* is GPLv3. **No licence or terms of use is stated for the hosted schema data.** | **NO** — silence is not permission. No chewie-NS schema is catalogued. |

> **This corrects the assumption the batch started from.** Institut Pasteur,
> EnteroBase and chewie-NS were expected to permit redistribution. Read directly,
> Pasteur forbids it, EnteroBase does not grant it, and chewie-NS grants nothing at
> all for its data. Only PubMLST does. If written permission is obtained from any of
> them for a specific scheme, add the scheme to `_SCHEMES` and flip that provider's
> `may_bundle`; nothing else has to change.

---

## 3. Why no cgMLST allele database ships in the release

Two independent reasons, both measured:

1. **Licence.** Of the twelve threshold-bound organisms, ten are bound to
   cgMLST.org schemes, one to an Institut Pasteur scheme — none of which may be
   redistributed. Only Salmonella and E. coli have a redistributable copy.
2. **Size.** One cgMLST allele database is 0.8–4.6 GB installed. Ten of them would
   be roughly 20 GB. The portable Windows release cannot carry that, and the scheme
   loader holds every allele in memory, so a 4.6 GB scheme is not usable on a
   typical laboratory laptop regardless of licence.

So the release ships the **definition pack**: kilobytes per scheme, and the folders
already in the right place.

---

## 4. What ships

`studio_scripts/stage_cgmlst_schemes.py` creates `<data root>/cgmlst/` with one
clearly-labelled folder per catalogued scheme, each containing:

- `README.txt` — which scheme belongs here, its provider, its terms, its target
  count, how to install it, and, for a scheme that cannot be packed, the sentence
  "This scheme CANNOT be packed into WMLSTudio: its provider does not grant
  redistribution."
- `scheme_slot.json` — the machine-readable pin: provider, database, scheme id,
  revision, target count, SHA-256 of the sorted target list, licence verdict, and
  the bound `threshold_scheme_key`.
- `targets.txt` — the pinned target names, when staged with `--fetch-targets`.

Nothing in a slot is sequence data until a scheme is installed into it. No sequence
data is committed to this repository, and `src/wmlstudio/resources/cgmlst/` is
gitignored like every other staged reference directory.

### Redistributable (PubMLST; the build may stage the alleles too)

| Key | Organism | Scheme | Targets | Threshold bound |
|---|---|---|---|---|
| `pubmlst:senterica-cgmlst-3002` | Salmonella enterica | cgMLST v2 (Enterobase) | 3002 | **yes** — `enterobase:senterica-cgmlst-3002` |
| `pubmlst:ecoli-cgmlst-2513` | Escherichia coli | cgMLST | 2513 | **yes** — `cgmlst.org:ecoli-2513` |
| `pubmlst:saureus-cgmlst-1716` | Staphylococcus aureus | cgMLST | 1716 | no (published cutoff is for 1,861 targets) |
| `pubmlst:abaumannii-cgmlst-2133` | Acinetobacter baumannii | cgMLST v1 | 2133 | no (published cutoff is for 2,390 targets) |
| `pubmlst:smarcescens-cgmlst-2692` | Serratia marcescens | S. marcescens cgMLST | 2692 | no — see §5 |
| `pubmlst:mtbc-cgmlst-1561` | M. tuberculosis complex | MTBC cgMLST | 1561 | no (published cutoff is for 2,891 targets) |
| `pubmlst:campylobacter-cgmlst-1142` | C. jejuni / C. coli | cgMLST v2 | 1142 | no curated cutoff |
| `pubmlst:nmeningitidis-cgmlst-1329` | Neisseria meningitidis | cgMLST v3 | 1329 | no curated cutoff |
| `pubmlst:hinfluenzae-cgmlst-1037` | Haemophilus influenzae | cgMLST v1 | 1037 | no curated cutoff |
| `pubmlst:bcc-cgmlst-1925` | B. cepacia complex | cgMLST | 1925 | no curated cutoff |
| `pubmlst:senterica-cgmlst-2750` | Salmonella enterica | SalmcgMLST v1.0 | 2750 | no (PulseNet cites the 3,002-target scheme) |

A PubMLST snapshot taken without credentials contains only submissions up to
2024-12-31. That is what makes it redistributable and it is also what makes it
**not the current allele set**; the snapshot records this and the app repeats it.

### Download-only (folder pre-created, one click, terms shown first)

| Key | Organism | Provider | Targets | Threshold bound |
|---|---|---|---|---|
| `cgmlst.org:lmonocytogenes-1701` | Listeria monocytogenes | cgMLST.org | 1701 | yes (Van Walle 2018, ≤7; ≤4 stricter option) |
| `pasteur:lmonocytogenes-1748` | Listeria monocytogenes | Institut Pasteur | 1748 | yes (Moura 2016, ≤7) |
| `cgmlst.org:kpneumoniae-2358` | K. pneumoniae complex | cgMLST.org | 2358 | yes (Glasgow 2025, ≤15) |
| `cgmlst.org:efaecium-1423` | Enterococcus faecium | cgMLST.org | 1423 | yes (three entries: ≤7, ≤20, ≤25 — different questions) |
| `cgmlst.org:efaecalis-1972` | Enterococcus faecalis | cgMLST.org | 1972 | yes (Glasgow 2025, ≤7) |
| `cgmlst.org:abaumannii-2390` | Acinetobacter baumannii | cgMLST.org | 2390 | yes (Glasgow 2025, ≤9) |
| `cgmlst.org:saureus-1861` | Staphylococcus aureus | cgMLST.org | 1861 | yes (Glasgow 2025, ≤24) — 4.6 GB installed |
| `cgmlst.org:smarcescens-2692` | Serratia marcescens | cgMLST.org | 2692 | yes (Kampmeier 2022, ≤12) |
| `cgmlst.org:cfreundii-3250` | Citrobacter freundii | cgMLST.org | 3250 | yes (Kieninger 2025, ≤10) |
| `cgmlst.org:paeruginosa-3867` | Pseudomonas aeruginosa | cgMLST.org | 3867 | citation only — no number is curated |
| `cgmlst.org:cdifficile-2147` | Clostridioides difficile | cgMLST.org | 2147 | **no** — see §5 |

---

## 5. Binding a threshold to a scheme

`cgmlst_schemes.threshold_for(key, locus_count=<targets found on disk>)` is the
lookup. A number is offered only when **all** of these hold:

1. the method is cgMLST (an MLST question raises, it is not answered);
2. the catalogue entry carries a `threshold_scheme_key`;
3. a `threshold_guidance` entry names that exact key;
4. the installed target count equals the pinned count **and** the published count.

Otherwise the citation is still returned and the number is withheld, with the
reason in words. When several publications are bound to one scheme (E. faecium has
three), the most recently reviewed source is suggested and the rest are returned in
`alternatives` — the same ranking `threshold_guidance.suggested_threshold` uses, so
the per-organism view and this per-scheme lookup can never show a person two
different numbers for the same comparison. `threshold_citations(key)` returns the same rows shaped for a PDF
report footnote: citation, DOI, URL, publication date, locator, scope, missing-data
policy and limitations.

Programmatic entry points (`from wmlstudio import cgmlst_schemes`):

| Call | Use |
|---|---|
| `prepare_library(data_root)` | create the folder layout on first run |
| `library_status(data_root)` | one row per scheme: folder, installed?, ready?, threshold |
| `threshold_for(key, locus_count=…)` | the cutoff lookup by catalogue key |
| `threshold_for_installed(scheme_path)` | the same, from a folder on disk; counts its targets |
| `identify_installed(scheme_path)` | which catalogued scheme a folder is, or `None` |
| `threshold_citations(key)` | citation rows shaped for a PDF report footnote |
| `download_plan(key)` | what a download button must show before transferring |

Three worked cases, all covered by tests:

- **Same count, different scheme.** PubMLST and cgMLST.org each publish a
  2,692-target *S. marcescens* cgMLST scheme. They share **no** target name at all
  (`SERR*` versus `SMDB11_RS*`, verified 2026-09-14). The Kampmeier ≤12 cutoff is
  offered for the cgMLST.org scheme and refused for the PubMLST one: equal counts
  are a coincidence, not an identity.
- **Different provider, provably the same target set.** The PubMLST "cgMLST v2
  (Enterobase)" Salmonella scheme and the cgMLST.org `Senterica` scheme have
  byte-identical sorted target sets (SHA-256
  `e819029e…`), and PubMLST names it the EnteroBase scheme. Leeper 2023 evaluated
  the EnteroBase-derived 3,002-locus core scheme, so the PulseNet ≤10 criterion is
  bound to the PubMLST copy. The same holds for *E. coli* at 2,513 targets
  (`cadaab19…`). The target set matches; the allele nomenclature and the caller do
  not, and the entry says so.
- **Superseded version.** Bletz 2018 published ≤6 on the 2,270-target
  *C. difficile* scheme, which the provider now lists as deprecated v1. The current
  scheme has 2,147 targets, so no number is offered for it at all.

---

## 6. Running the staging script

```bash
# Definition pack only (default; kilobytes, no sequence data)
python studio_scripts/stage_cgmlst_schemes.py

# …with each scheme's current target list, staged only when it still matches the pin
python studio_scripts/stage_cgmlst_schemes.py --fetch-targets

# Re-check every pin against the live services (exit 1 if any target set moved)
python studio_scripts/stage_cgmlst_schemes.py --verify

# Allele data for one redistributable scheme (gigabytes; refuses every other key)
python studio_scripts/stage_cgmlst_schemes.py \
    --with-alleles pubmlst:senterica-cgmlst-3002 --max-bytes 4000000000

# Attribution block for studio_packaging/THIRD_PARTY_NOTICES.md
python studio_scripts/stage_cgmlst_schemes.py --notices
```

`--verify` is the guard against silent drift. A changed target set is a different
quantity: any cutoff bound to the scheme was published for the pinned set and does
not carry over. The script reports the change and never re-pins by itself.

---

## 7. What a cgMLST result is not

- An allele difference count is not a SNP count and not a percentage identity.
- A distance is only meaningful over targets called in **both** isolates; the
  shared-target denominator travels with the number, and insufficient overlap is
  not zero distance.
- A missing, ambiguous or excluded target is unknown, not identical.
- A published cutoff is evidence to review with its scope, not a setting to apply.
  `threshold_guidance.record_decision` still requires the scheme fingerprint, the
  caller, the missing-data policy and a local justification before any number is
  approved for use.
- A cluster is a signal for epidemiological review. It is not proof of transmission,
  of the direction of transmission, or of clinical causation.
