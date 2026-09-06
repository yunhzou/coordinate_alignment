# Example 5: full-bank scan at retro tolerance 1.0

## Scope and defaults

The retro default is now bond-order tolerance **1.0**. Core AAM and reusable
fragment-matcher defaults remain **0.5**. Assembly ranks complete covers by
fewest retained matched fragment units first, then distinct precursor species
and explicit-atom retention. This benchmark measures **detection and persistence**;
it does not measure full-bank assembly, ranking, or establish the ground-truth
set's blind recommendation rank.

- Bank: `data/mcule/merged_fast_delivery_with_inventory.csv.gz`, **155,305** rows.
- Bank SHA256: `847783ec450d437942c2be5cf4fc10011ce988536185ff79bdb8214ac40f4472`.
- Target: `C/C(NC1=C(C(C)C)C=CC=C1C(C)C)=C/C(C)=N/C2=NC3=C(C=C(C#N)C=C3)S2`.
- Explicit H, all source seeds, full augmented matching, branch cap **100**,
  minimum fragment size 1, no sweep and no candidate limit.
- Native-enabled matcher; no core AAM algorithm changes for this experiment.
- Full typed detection checkpoints and compressed result records are saved.

## Measured runtime

**27 minutes 25 seconds** from the first pilot task starting to the final
recovery task finishing (Slurm timestamps 2026-09-05 23:12:45–23:40:10).
This includes the staged launch, queue/node-start delays, watchdog failures,
and recovery. Subsequent archive integrity validation is separate.

256 logical shards, 40 CPUs and 80 GB per task. The peak simultaneously
allocated useful-job resources, reconstructed from Slurm intervals, were
**2,560 CPUs across 40 nodes**. Useful-job allocations consumed **602.38 CPU-hours**,
including failed attempts and recovery; these are allocated hours, not measured
CPU utilization. Six stuck CONFIGURING tasks were cancelled and replaced on
other nodes; their idle allocations are excluded from those resource totals.

Initial saved-record timings (155,300 records, before five recoveries):

| Per-precursor stage | Median | 95th percentile | 99th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: |
| Detection checkpoint events (155,301) | 3.164 s | 14.764 s | 32.333 s | 564.589 s |
| Full saved pipeline (155,300) | 4.837 s | 22.988 s | 43.346 s | 565.961 s |

These are overlapping per-source wall times, not numbers to add to obtain the
parallel scan time. The detection sample contains one checkpoint whose initial
archive write did not finish. Recovery timings below are reported separately.

## Watchdog recoveries

Five of 256 initial shards reached the operational watchdog, each with only
one source left unsaved. All five sources were recovered without changing
matching settings or dropping candidates. A saved detection was reused where
available; the other four were rerun with 39 parallel seed workers on 40 CPUs.
The operational limit remained ten minutes per job.

| Source | Recovery action | Detection | Fully saved recovery |
| --- | --- | ---: | ---: |
| MCULE-4131422119 | Resume existing detection checkpoint | 442.368 s in original attempt | 46.095 s |
| INVENTORY-000385 | Parallel seed retry | 174.816 s | 201.266 s |
| MCULE-7450586738 | Parallel seed retry | 59.488 s | 95.948 s |
| MCULE-1332591442 | Parallel seed retry | 164.501 s | 183.347 s |
| MCULE-5090262351 | Parallel seed retry | 52.296 s | 76.099 s |

Recovery pipeline times exclude Slurm startup and include archive saving.
The initial failed work is already included in the 27:25 elapsed measurement.

## Saved artifacts and reproduction

Shared cluster run root:

```text
/project/yunhengzou/coordinate_alignment/retro_runs/t05_full_merged_tol1_20260905
```

`parts/` contains the original shards, including saved records from interrupted
shards. `recovery_1/parts` through `recovery_4/parts` contain the five additional
records. Pilot parts 0–7 and their logs are symlinked to the original run under
`/h/399/yunhengzou/coordinate_alignment/data/retro_runs/`.
`checkpoints/`, `pilot_checkpoints/`, and recovery checkpoint folders preserve
typed detection intermediates. The large artifacts remain on the cluster,
not in Git.

- Initial jobs: 432179 (pilot), 432187 (main), 432282 (node-start replacements).
- Recovery jobs: 432427, 432442, 432443, 432444.
- Integrity-validation job: 432446.
- Launch: `hpc/search_retro_bank.sbatch` with 256 shards.
- Checkpoint recovery: `bench/prepare_catalog_recovery.py` and
  `hpc/resume_catalog_tail.sbatch`.
- Record identity/JSON integrity check: `bench/validate_catalog_run.py` against
  the original bank and all five output-part directories.

## Interpretation

The flushed-record progress audit accounts for **155,305 / 155,305** source
identities across original and recovery logs. Summaries for 152,278 sources
report **221,517,917 candidate records**; this is a lower bound, not the exact
full-bank count or a count of unique chemical fragments. Remaining saved rows
are in interrupted shards without final summaries.

The separate full JSON/gzip integrity check (job 432446) hit its 580-second
watchdog and did not produce `validated_records.json`. Thus a completed
archive-content/identity audit is **not** claimed. Its extra 9:46 allocation is
not included in the 27:25 scan timing. The original archives and checkpoints
are preserved; no matching was rerun for this check.

The median source remains a few seconds, but this is **not a sub-ten-minute
full-bank run**: scheduling and rare slow sources determine completion.
No controlled tolerance-0.5 run of this same full bank/target/resource schedule
was performed, so these measurements do not establish a tolerance-only
slowdown factor. Processing every bank row is also distinct from exhaustive
matching: a branch-cap hit must remain visible as incomplete search evidence.
