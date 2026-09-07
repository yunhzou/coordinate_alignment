# Golden: smaller-first AAM and independent sweep-cut seeds

All 1,851 records were attempted and evaluated. No pending cases or evaluator errors remain. A completed evaluation is not necessarily a completed search: 19 searches are explicitly partial.

## Reference recovery

Same 1,760 fully product-heavy-annotated records and scoring convention as the previous fixed search. All input hydrogens are explicit during AAM; Golden reference correctness is checked on the annotated heavy-atom mapping, allowing chemically equivalent endpoint symmetry. This is not a hydrogen-identity accuracy claim.

| Metric | Previous fixed search | New policy |
|---|---:|---:|
| Reference found in saved AAM result | 1,607 / 1,760 (91.31%) | **1,725 / 1,760 (98.01%)** |
| Top-ranked reference-equivalent mapping | 1,302 / 1,760 (73.98%) | **1,354 / 1,760 (76.93%)** |
| Reference not recovered from completed result | 149 | 24 |
| Unresolved reference recovery | 4 | 11 |

The 91 incompletely annotated records are retained separately in `cases.csv`, not used in these accuracy denominators. The old reference score is the improved verification of the **unchanged old search archives**, not a higher-seed diagnostic result.

Of the 19 incomplete searches, 18 have verified reference-equivalent witnesses in saved cuts; one remains unknown. These witnesses count as recovered, but their full-sweep top-1 is deliberately unknown. Of the 11 unresolved complete-reference cases, ten are reference-verification budget limits on completed searches and one is a partial search. Unknowns remain in the denominator.

127 previous misses became recovered; nine previous recoveries became confirmed misses. Two previous unknowns became recovered, and two previous recoveries are now unresolved. Net verified recovery improved by 118 cases. **This is an improvement, not a no-regression or completeness guarantee.**

The nine confirmed losses are 737, 775, 780, 846, 850, 1230, 1380, 1400 and 1817. Case 737 retains R→P, isolating the changed cut streams as the search-policy change. The others reverse direction and also receive new cut streams. Case 1817 is uncapped: direction/seed-sensitive fragment commitment still matters even without cap stops. The CSV and saved graphs retain the evidence; old successes were not merged into new results.

## Case 1665

**Recovered at top-1**, with three seeds per cut, cap 100 and tolerance 1.0. Full search plus archive persistence took **21.98 seconds on one worker**. The planner chose P→R (62 explicit atoms into 140). The recorded reference witness is from source cut `(9, 10)`; the uncut result had no completed terminal at this cap. This successful run used actual sweep exploration, not a guided reference seed.

## What changed

- Each cut receives its own stable pseudorandom stream derived from its edge set; repeated executions remain reproducible. No-cut retains the old seed 42.
- The endpoint with fewer explicit atoms, including H, is searched into the larger endpoint. Ties retain R→P: 986 records search R→P, 865 search P→R.
- Direction selection is a Python orchestration object. Native matching, fragment growth and symmetry compression remain unchanged. Raw AAM graphs retain actual search orientation; scoring and concrete mapping conversion preserve original R/P identities.
- Three seed orders per cut, branch cap 100, tolerance 1.0 and the complete single-edge sweep are unchanged from the previous full campaign. This is not the separate cap-2000 diagnostic.
- Map labels are removed before canonical input ordering. Reference labels never choose direction, seeds or fragments. No reference-guided reruns are included.

## Runtime and watchdogs

- 128 requested one-worker CPU slots at peak; no within-reaction parallelism. Slurm sometimes allocates two logical CPUs per one-worker request, so allocated CPU-hours are not worker compute-hours.
- 1,832 full search checkpoints; median completed search **5.36 seconds**, maximum **297.62 seconds**. Completed-search process CPU total: **10.15 hours**, excluding interrupted searches. Summed completed-search elapsed time is 36,960 seconds, not parallel wall time.
- Five-minute search watchdog and ten-minute Slurm job limit. Twenty search processes reached the watchdog; one already had its finalized checkpoint, leaving 19 incomplete searches. No OOM search exits were recorded.
- Maximum recorded completed-search process RSS: approximately **5.26 GiB**. This is not a measurement of every interrupted process's peak.
- Submission to final recorded outcome: approximately **26.0 minutes**, including cluster delay and saved-cut evaluation repairs. Seventy-five jobs stalled in CONFIGURING without starting Python. Only unstarted jobs were cancelled/resubmitted away from those nodes; no completed search was restarted. `resubmission_unstarted.json` records the 234 unstarted tasks moved, including ordinary pending tasks.
- `slurm_accounting.psv` retains main and replacement allocation records: **79.19 allocated CPU-hours**, including stalled allocations and process overhead. Separate saved-cut recheck allocations are excluded. Do not present this as 79.19 hours of AAM computation.

The frozen campaign's partial-evaluator constructor had an error, repaired without changing the frozen search engine. Saved cuts were reevaluated with the corrected adapter; later partial checks stop at a verified representative instead of ranking every terminal, since no full-sweep top-1 claim is made. Full-search scoring was not changed or rerun for this repair. All raw search checkpoints remain intact.

## Saved artifacts and reproducibility

Full cluster run:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`

Each numbered case has original `input.json`, `reference.json`, `orientation.json`, per-cut JSON under `cuts/`, a finalized `cuts/aam.pkl.gz` when complete, evaluation and process-status records. Partial evaluations identify the precise witness cut. The immutable `engine/` snapshot and hashes in `manifest.json` identify the executed code; its recorded base revision precedes the policy commit, so the hashes, not that base revision alone, define the run.

Local result folder:

`/h/399/yunhengzou/coordinate_alignment/reports/golden_policy_full_20260907`

`cases.csv` includes every record and archive location. `summary.json` includes all misses, unknowns, gains and losses. `manifest.json` includes input and engine hashes. Implementation regression tests: **134 passed**.
