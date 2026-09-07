# Golden: top-five and bond-event tolerance

All **1,851** records remain in the denominator. The **1,832 completed search archives** were ranked; the 19 incomplete original searches have unknown full-search ranks. All 1,851 evaluation workers exited successfully. No AAM searches were rerun and no ranking policy or chirality setting was changed.

## Deduplicated representative rankings

| Shortlist | Recovered | Percentage |
|---|---:|---:|
| Top 1 | 1,435 | 77.53% |
| Top 3 | 1,741 | 94.06% |
| Top 5 | **1,768** | **95.52%** |
| Top 10 | 1,781 | 96.22% |
| Top 5 including every tied score at fifth place | 1,777 | 96.00% |

Top-five with ties contains a mean of 6.23 classes per completed case. These classes use the same heavy-atom chemical-equivalence certificates as the alternatives viewer. Each is ranked by its best saved explicit-H representative: most mapped heavy atoms, most mapped total atoms, fewest bond events, then deterministic atom-index tie-break. The reference never selects candidates or breaks ties. Top-1 reproduces the previous all-record result exactly.

## Tolerance around the best bond-event score

Only candidates in the **best heavy-atom and total-atom coverage tier** are eligible. The tolerance adds to the minimum event count within that tier. An event is one broken bond, formed bond or bond-order change under the existing ranker; hydrogens are included. This is a ranking-score tolerance, **not a change to AAM's isomorphism tolerance**.

| Extra allowed bond events | Recovered | Percentage | Mean candidate classes |
|---|---:|---:|---:|
| 0 (all best-score ties) | 1,660 | 89.68% | 1.33 |
| 1 | 1,717 | 92.76% | 1.55 |
| **2** | **1,772** | **95.73%** | **3.19** |
| **3** | **1,783** | **96.33%** | **5.31** |
| 5 | 1,788 | 96.60% | 10.28 |

Candidate-count means use the 1,832 completed archives; accuracy percentages use all 1,851. These descriptive results suggest an event window is useful, but selecting its width using this dataset is tuning on the benchmark, not independent validation. No claim of chemical plausibility follows from event count alone.

## Compressed alternatives behind the selected top-five classes

The displayed top-five representatives recover 1,768 references. Independently querying the full explicit-H compressed paths associated with these selected classes recovers **1,780 / 1,851 (96.16%)**: twelve additional references. There are **52 verified absences from the selected saved families** and **19 unknown ranks due to incomplete original searches**. No completed-archive family query remained unresolved.

The twelve additional cases are 44, 124, 1112, 1118, 1313, 1506, 1522, 1589, 1590, 1595, 1786 and 1817. The verifier checks exact heavy-atom reference correspondence including unmatched atoms, modulo endpoint chemical symmetry, and retains explicit-H feasibility. Full positive witnesses are saved in each case's `family_verification.json` when a symbolic query was necessary.

This is coverage of families selected by representative rank, not optimization of the event score over every family member. The event-window table reports **representative** recovery only. Do not substitute the family-recovery percentage for that table.

## Artifacts and runtime

Full class tables, ranking metrics, compressed-family query results, source snapshot, worker logs and statuses:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_topk_20260907`

Original immutable search archives:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`

Script: `bench/golden_topk.py`. Local `cases.json`, `summary.json`, `manifest.json` and `slurm_accounting.psv` retain all case outcomes, the executed source hashes and allocation accounting.

One CPU per worker, 16 GiB memory requested, at most 128 main-array tasks concurrently. Each process had a five-minute watchdog and each Slurm task a six-minute limit. The maximum measured completed scoring time was **99.01 seconds**; summed scoring elapsed time was **4,500.76 seconds** (not parallel wall time or allocated CPU-hours).

Main arrays: 440574 (indices 0–925), 440575 (indices 926–1850). Unstarted CONFIGURING tasks were rescheduled via 441244, 441245 and 441973. No completed scoring calculation was repeated. All 1,851 workers exited zero. **147 regression tests pass**, including score-window coverage priority and fifth-place tie inclusion.
