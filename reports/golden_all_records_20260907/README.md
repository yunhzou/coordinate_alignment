# Golden benchmark: all 1,851 records

No records excluded. Accuracy checks the exact reference heavy-atom relation, including unmatched atoms, modulo endpoint chemical symmetry. AAM uses explicit hydrogens; this is not a hydrogen-identity accuracy claim.

| Metric | Result |
|---|---:|
| Fixed-search reference-family recovery | **1,818 / 1,851 (98.22%)** |
| Top-ranked representative correct | **1,435 / 1,851 (77.53%)** |
| Separate union including previous opposite-direction searches | **1,833 / 1,851 (99.03%)** |
| Reference absent from saved fixed-search families | 32 |
| Unresolved recovery due to incomplete original search | 1: case 1793 |

There are 19 incomplete original searches. Their full-sweep top-1 remains unknown and they stay in the denominator; 18 nevertheless have certified reference witnesses in saved cuts. Combined-search coverage is not a single fixed-policy result. This evaluator is not a reproduction of the original paper's CGR scoring implementation, so these numbers should not be described as an identical-protocol comparison.

## The 91 one-sided records

Direct inspection of the original RDF found 190 product heavy atoms with labels absent from the reactant side across these 91 records. These are one-sided atom sources, not automatically annotation errors.

- Exact reference recovered: **87 / 91 (95.60%)**.
- Top-ranked representative correct: **81 / 91 (89.01%)**.
- Saved-family misses: **7, 1358, 1377, 1475** (dataset zero-based indices).

The evaluator previously allowed subset agreement in symbolic queries for this group. It now requires the complete reference relation: additional incorrect heavy-atom pairs are rejected, not ignored. Positive symbolic witnesses must pass the same independent chemical-equivalence certificate as representative witnesses. No reference-guided deletions or atom reassignment were used.

## Search and ranking

The original search remains unchanged: three seeds per cut, branch cap 100, tolerance 1.0, smaller-first orientation, independent cut seed streams and single-edge cut sweep. Ranking first maximizes mapped heavy atoms, then total mapped atoms including H, then minimizes bond events; atom indices break ties. Top-1 is the stored representative, not the best mapping within a compressed family. Top-k was not evaluated here.

The 1,760 previously evaluated cases retain their certified results (1,731 recoveries). All 91 one-sided cases were strictly rescored from saved archives. No AAM searches were rerun. There are **145 passing regression tests**, including rejection of a mapping that contains all reference pairs plus an incorrect extra pair.

## Saved artifacts

Strict one-sided evaluations, frozen source snapshot, process logs, job status and timings:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_all_records_20260907`

Original search archives:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`

Corrected complete-reference verification:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_domain_rescore_v2_20260907`

Slurm jobs: 440479 and 440570. Only unstarted array slot 38 was moved off a node stuck in CONFIGURING; no completed calculation was restarted. One worker CPU per task, 24 GiB memory request, five-minute process watchdog and six-minute Slurm limit.

`cases.csv` records every case and its authoritative evaluation path. `witnesses.json` saves the corrected/new positive evaluations and full available mappings. Original cached positive witnesses remain in the original archives. `manifest.json` records source snapshots; `slurm_accounting.psv` retains allocation accounting.
