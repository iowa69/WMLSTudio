# Validation hypotheses

These are software validation alternatives, not biological study hypotheses.

| Hypothesis | Test evidence required | Status |
| --- | --- | --- |
| Biological variation changes an allele | Practice ST1→ST2 differs at one locus; exact-caller and GUI tests | Tested |
| Parsing or reverse-strand handling loses loci | Wrapped, compressed, reverse-strand and chunk-boundary tests agree | Tested |
| Null: identical profiles have zero distance | Two identical practice profiles yield 0/7 differences | Tested |
| Sampling overstates read QC | Read-prefix limit tests and real FASTQ report disclose sampled=true | Tested |
| Database mismatch creates false type calls | Changed fingerprint/mixed-scheme comparison tests reject pairing | Tested |
| Missing evidence masquerades as similarity | Practice partial profile is disconnected at 0.95 shared coverage | Tested |

Reflection: real reference data contained incomplete profiles and ambiguous allele
sequences. The loader reports exclusions explicitly; it never matches ambiguity
codes as nucleotide wildcards. All 162 staged schemes load with those disclosures.
Agreement with the older caller on four assemblies is useful regression evidence,
not a sensitivity/specificity study or cgMLST equivalence claim.
