# Full inventory scan after the core-AAM fix

Follow-up: the [occupation compression fix](occupation_compression_fix_2026-09-05.md)
recovered sources 538 and 594. Combined progress is now 1,918/1,919; only 620 is
unresolved. The measurements below describe the earlier scan and recovery.

**Incomplete: 1,916 of 1,919 source records saved after checkpoint recovery.**
The three unfinished sources are not classified as nonmatches. No assembly or
viewer generation was run in this detection benchmark.

## Scope

The entire `data/inventory/processed/inventory_structure_bank.csv.gz` was
submitted, not a sample. This is the 1,919-compound inventory benchmark bank,
**not** the larger merged fast-delivery catalog. SHA256:
`dfcbee63a0aab6571f3edaf6358ec9e13135682d85ca3d8589142fb6d62e00b9`.

Target: the same BIAN molecule used in the preceding full-bank comparison:

```text
CC(C)C1=CC=CC(C(C)C)=C1/N=C2/C(C3=C4C(C=CC=C42)=CC=C3)=N/C5=C(C(C)C)C=CC=C5C(C)C
```

Core revision `fcf26d3`; checkpoint-enabled scanner `6cb60de`.
Explicit hydrogens, tolerance 0.5, branch limit 100, no sweep, orbit-representative
initial seeds, no candidate/seed-count cap, no target-coverage filter, and all
results saved. These are the preceding exact-workflow settings, not the old
approximate scan. A saved result may still carry a branch-cap/incomplete flag.

## Measurements

| Phase | Resources | Outcome |
| --- | --- | --- |
| Fresh scan, job 432109 | 28 nodes × 48 CPUs | 1,913 saved; 23 of 28 shards completed; five hit the watchdog |
| Checkpoint-only recovery, job 432138 | Two small tasks | Two saved results recovered without any new matching; 92.87 s and 136.46 s inside the worker |
| Fresh missing-source recovery, job 432137 | Four nodes × 48 CPUs | One additional result saved in 505.57 s; three hit the watchdog again |

The initial allocation ran from 18:43:55 through 18:53:41 Slurm local time
(9 min 46 s including launch/termination). Each task used a 580-second external
timeout and a ten-minute Slurm limit. Tail recovery is a separate allocation;
these measurements must not be presented as one complete fresh-scan time.

For recovered source `INVENTORY-000400`, detection took 288.54 s, typed checkpoint
writing took 25.68 s, and the remaining approximately 191.36 s includes archive
construction, encoding, compression, and worker overhead. Its 40,062 candidates
were saved. This is one source, not a bank-wide timing decomposition.

The preceding pre-fix full-bank attempt saved 1,901/1,919 before its watchdog.
Neither fresh scan completed, so there is **no measured complete-bank speedup
ratio**. The independently verified 2.33× core improvement does not imply a
2.33× end-to-end scan improvement.

## Saved-record integrity audit

Job 432151 read all 31 result files and verified gzip/JSON syntax, exact bank
row IDs and SMILES, and absence of duplicate source records. It took 4 min 40 s
on four allocated CPUs; this is validation time, not search time. Its exit code
1 explicitly signals the three missing rows, not corrupt saved files.

- 1,916 valid saved source records; 1,830 contain candidates.
- 3,448,290 fragment candidates, not distinct reactant sets or mechanisms.
- Status counts: 1,200 `capped`, 611 `rough`, 78 `matched`, 27 `no_match`.
- Only 105 records carry `complete=true`; completed persistence is not an
  exhaustive-search guarantee. `rough` is the detector's existing completeness
  classification, not an additional shortcut introduced in this run.
- All saved branch limits are 100. Missing zero-based rows: 537, 593, 619.

Compact counts are versioned in
[`full_bank_core_fixed_2026-09-05.json`](full_bank_core_fixed_2026-09-05.json).
The complete per-record locations and SHA256 hashes remain in the cluster's
`validated_index.json`.

## Remaining bottlenecks

Short worker probes sampled the actual hot child processes using SIGUSR1:

- `INVENTORY-000538`: `materialize_target_coverage_orbit` → native
  `occupation_orbit`, while processing augmented AAM results.
- `INVENTORY-000594`: the same distinct-fragment-occupation expansion.
- `INVENTORY-000620`: native `grow_island` inside augmented AAM, observed in
  repeated samples. Native debugger attachment was unavailable, so the exact
  C++ subroutine is not yet identified.

The branch cap does not cap the number of distinct occupations materialized
after AAM. In native growth, the extension-level branch check follows extension
and deduplication; a cap therefore is not a wall-clock guarantee. Neither
observation justifies truncating occupations or changing search semantics.
No such algorithm changes were made during this benchmark.

Existing fixed per-source seed budgets also leave CPUs idle when other sources
finish. Full-node recovery improved source 400 but did not resolve the other
three; scheduling alone is not sufficient.

The short probe job 432148 was stopped once traces were collected. An earlier
diagnostic-only job 432144 was cancelled because rearming a faulthandler timer
in a fork hook blocked its workers. That instrumentation was removed; its
timings are excluded. Neither diagnostic replaced a saved result.

## Saved evidence and continuation

Cluster root:
`/h/399/yunhengzou/coordinate_alignment/data/retro_runs/full_bank_core_fixed_20260905/`

- `run_manifest.json`: bank identity, target, configuration, allocation.
- `parts/` and `tail/parts/`: compressed records; one complete record per source.
- `checkpoints/` and `tail/checkpoints/`: complete typed detections saved before
  archive construction. An unfinished detection has no completed checkpoint.
- `progress_audit.json`: frozen original six-source recovery list.
- `combined_progress.json`: combined saved-record progress, including the three
  unresolved IDs and their exact source SMILES.
- `worker_probe/parts/part_{2,4,5}.stacks.txt`: sampled worker stacks.
- `logs/`, `tail/logs/`: individual elapsed times and checkpoint milestones.

Use `bench/catalog_scan_progress.py --continuation` to audit recovery without
overwriting the frozen recovery list. `bench/validate_catalog_run.py` validates
gzip/JSON syntax and exact bank row identities, indexes record locations and
hashes, and reports cap statuses. It is not a chemical correctness proof.

The next optimization targets are the two identified boundaries above, with
saved full-AAM inputs and old-core comparisons. Do not rescan the 1,916 saved
sources merely to recover the three unresolved ones.
