# WMLSTudio workbench revision contract

## User workflow

Import FASTA/FASTQ → assign automatic, manual or unknown organism per sample or
batch → review the run plan → analyse in the background → explore linked results
in a searchable library → select a cohort and a typing scheme → inspect and style
the minimum spanning forest → export selected samples with explicit highlights.

The desktop is native Qt with a dark palette, ordinary menus and contextual
dialogs. Full-page opacity effects are prohibited: repaint correctness matters
more than page-transition animation. The helix remains animated when enabled.

## Data contract

Stable project sample IDs join all evidence. Metadata namespaces are `organism`,
`workflow`, `hydra`, `cluster`, and user annotations. `workflow.typing_mode` is
`auto`, `manual`, or `unknown`; manual labels do not require an installed scheme.
An unknown organism is a valid state, not an input error.

Primary MLST/QC results and additional cgMLST/wgMLST profile snapshots are distinct.
Comparisons require an identical scheme fingerprint and report shared-locus
denominators; missing loci are never silently treated as matches. A provisional
organism inferred from scheme evidence must not overwrite a user's confirmed label.

Managed storage copies inputs into a project-specific library with source hashes;
optional genus/species/ST organisation and ST filename suffixes affect managed
copies only. No original input is renamed, edited or deleted. Copy cancellation,
name collisions, duplicate basenames and project reopening must not lose data.

HYDRA findings link through explicit sample mapping. Ambiguous name matches require
resolution. Gene evidence is not a susceptibility phenotype. Original HYDRA
parameters, database versions and report hashes stay available after linking.

## Reference updates

The user confirmed non-commercial research as the intended use. This does not
grant redistribution rights: cgMLST.org database contents remain excluded from
portable releases, and individual downloads retain explicit provider terms.

User-initiated online discovery and installation must keep an old snapshot usable
until its replacement is fully downloaded, validated and fingerprinted. Downloads
run in the background, are cancellable and do not send sample sequences. API access
restrictions are displayed; a public snapshot is not represented as a complete
current authenticated database. Existing results stay tied to their original snapshot.

## Acceptance gates

- Rapid page switching, resizing, sorting and worker completion repaint correctly.
- Import dialog supports per-row and batch labels, including unknown and custom taxa.
- Managed-copy integrity and input immutability survive cancellation and collisions.
- Auto-typing distinguishes strong, weak and tied evidence with explicit provenance.
- Cohort and scheme selection controls exactly which profiles enter the MST.
- Manual graph colors, layout and highlights persist without changing distances.
- Linked metadata remains consistent between sample tables, graph and reports.
- Reports export an explicit cohort, including selected highlighted cluster members.
- Reference network failures and restricted responses are reported without corrupting caches.
- Source tests, native UI workflow tests, local real-data regressions and frozen
  Windows tests pass before a retest package is published. Clean Windows 11 user
  acceptance remains a separate requirement, not a claim implied by CI.

This contract is the target, not a declaration that all gates have passed.
