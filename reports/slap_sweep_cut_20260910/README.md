# SLAP plus blind single-edge sweep

## Result

**Single-edge sweeping substantially improves SLAP, but this test does not reach 99%.**
All 1,851 Golden reactions were evaluated using all saved alternatives, not top-1.

| Search | Reference mapping recovered | Coverage |
| --- | ---: | ---: |
| Uncut control: one input ordering, both directions and bond modes | 1,661 / 1,851 | 89.74% |
| Previously saved expanded SLAP: ten orderings plus all-atom symmetry run | 1,691 / 1,851 | 91.36% |
| New uncut + single-edge sweep, both directions and bond modes | 1,795 / 1,851 | 96.97% |
| Union of new sweep and previous expanded SLAP results | 1,797 / 1,851 | 97.08% |

The sweep adds 134 cases over its own uncut control and 106 beyond the previous
expanded baseline. Binary sweep alone recovers 1,785; weighted sweep alone 1,762.
The combined 54 misses comprise **50 fully swept, valid cases and four incomplete
cases** (1503, 1665, 1723, 1776). Completing those four alone cannot reach 99%.
This is an observation for this specific one-order, single-edge perturbation
policy, not proof that additional SLAP search policies cannot do better.

Of 355,104 planned calls, 352,210 were recorded (352,208 successful, two upstream
exceptions). There were 65 watchdog-limited case/direction/mode variants out of
7,404; their previous completed cuts remain included. In total, 2,894 planned
calls did not complete. No returned prediction changed endpoint chemistry or
failed the strict mapping-format validation. The two exceptions were upstream
`min() arg is an empty sequence` errors on case 836, reverse binary cuts
(31,67) and (35,71); that case was recovered by other blind searches.

## Timing and validation

- Mapping campaign: **8 min 09 s** from first allocation start to last finish,
  including staggered starts, watchdog waits and persistence; **not pure compute time**.
- Maximum 32 allocations × 32 single-threaded CPU workers = 1,024 workers.
  Individual allocations lasted at most 321 seconds; each variant had a
  300-second watchdog. The initial eight pending allocations received reduced
  memory reservations; no completed work was repeated.
- Separate parallel evaluation: **35 seconds**, 64 shards.
- Completed-call native mapping CPU: **14.479 CPU-hours**, or 28.160 CPU-seconds
  per Golden reaction summed over directions, modes and cuts. This is **not total
  campaign CPU**: interrupted in-flight calls have no final CPU reading.
- Completed-call mapping + graph preparation + export: 14.934 CPU-hours;
  encoding: 0.246 CPU-hours, measured separately. Persistence and process startup
  are excluded from these instrumented CPU measurements.
- Slurm reserved 81.493 CPU-hours for mapping allocations, including idle slots.
  Its `TotalCPU` fields report zero, so they cannot supply the missing actual CPU
  consumption. Do not compare the completed-call CPU figure as an all-cost total.
- Three adapter/archive tests passed: exact uncut parity with native SLAP on
  three reactions in both directions and modes; cut/cache isolation; full native
  checkpoint recovery and unchanged endpoints for all cuts of an unbalanced case.
- The strict evaluator source is unchanged from the saved expanded baseline.

## Saved results

Full native outputs, every cut, reference witnesses, logs, frozen upstream source,
and separate timing records:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_edge_sweep_full_20260910`

Compact Git-tracked files in this folder: `summary.json`, `case_metrics.json`,
`manifest.json`, `slurm_accounting.json`, submission records and artifact paths.
`case_metrics.json` provides native archive pointers for one recovered reference
witness per direction/mode. Every alternative remains in the full archive.

## Protocol

The previous AAM optimization campaign was stopped at the user's request.
Slurm array `466079` and dependent comparison `466114` were cancelled.
Existing AAM checkpoints remain under
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/shared_capacity_full_20260910`.

This separate experiment tests pinned upstream SLAP
(`ea248fd9494f52f4865193e87a98cc92c62b5f9e`) on all 1,851 Golden reactions.
Run uncut plus every single source-edge deletion, independently in both
directions and both native binary/weighted modes. Explicit H and native
element-balancing placeholders are retained. No atom is deleted. Rebuild graph
initialization after each cut; otherwise WL caches would describe the wrong graph.

Cuts perturb SLAP's input graph/objective, not its internal algorithm. This is
not identical to AAM's fragment cut constraint. Restore output labels onto the
original molecular endpoints before evaluation. Keep every returned native
alternative, label pool and LAP record. Reference mappings never select cuts,
order searches or stop cases. Native SLAP pruning is unchanged.

Report uncut versus sweep recovery and union with the saved expanded-SLAP
baseline (1,691/1,851 = 91.36%) separately. Use the existing strict heavy-atom
relation / endpoint-chemical-symmetry evaluator, including unmatched atoms.
Explicit H participates in search; reference-H assignment is not this metric.
Failures and unattempted cuts remain visible with the full 1,851 denominator.

Save mapping CPU/wall time separately from graph preparation, output export,
encoding, persistence and evaluation. Fixed hash/random seeds, single-threaded
workers, 300-second process watchdog per case/direction/mode (across its cuts),
and cluster-parallel independent cases. Partial cut outputs survive interruption.
`native.bin` contains independent gzip/pickle frames; committed JSONL records
identify each frame's offset, length and SHA256. Only load these trusted archives.

The home filesystem is full, so development is isolated on branch
`experiment_slap_sweep_cut` at
`/project/yunhengzou/coordinate_alignment/slap_sweep_worktree_20260910`.
No AAM core algorithm change is involved.
