# 140 elementary steps: matching tolerance 1.0

[Offline AAM/SLAP viewer — all 140 cases](viewer.html). Opens on TS_11 (case135).
The shared original white R/P/TS viewer uses yellow buttons to select saved witnesses.
Case123 explicitly displays its separate cap200 follow-up. No ground-truth
mapping labels exist for this holdout; event counts are not chemical accuracy.

## Result

Using the same cap policy as the preceding comparison (100 generally; separately
tested 200 for case123), **case135 improves from 7 to 5 events, and all other
139 best saved scores are unchanged**. No score regressions under that policy.

| Best observed common-score comparison | Cases |
|---|---:|
| AAM = SLAP | 134 |
| AAM fewer events | 6 |
| SLAP fewer events | 0 |

This does not prove global optimality, equality of solution families, or
chemical correctness. In particular, a wider matching tolerance does not make
the finite greedy/capped search monotonically more complete.

## One-factor search protocol

- Copied the preceding benchmark's **frozen engine**, with byte-identity checks
  on its source files and all 140 input records. No core changes.
- Changed only `iso_tolerance`: **0.5 → 1.0**.
- Kept explicit H, event scoring threshold 0.5, bond floor 0.2, 10 seed orders,
  cap100, both directions, sweep cuts and deterministic cut-seed construction.
- 16 CPUs per search, at most 32 concurrent tasks. Initially the healthy
  `cpunodes` nodes; pending tasks subsequently admitted to `cpunodes_nia`.
- Same saved native XYZ SLAP outputs reused byte-for-byte: changing our
  matching tolerance does not change SLAP. SLAP H-group scoring refinement
  rerun, all 176 assignment-pool minima resolved.
- Full archives, exact mappings, input/engine manifests, task environments,
  timing phases and analysis artifacts retained. Search watchdog 300 seconds,
  Slurm search allocation limit 10 minutes. No search failures or watchdog kills.

All 280 cap100 searches completed. Full mappings were returned in 139/140 R→P
calls and 140/140 P→R calls, covering all 140 cases bidirectionally. Caps were
encountered in 65 R→P and 68 P→R calls; cap hit is not synonymous with failure.

### Case123: keep the cap policies separate

At strict cap100, R→P has no full result and P→R has a best saved score of 19.
The old comparison used a **separate cap200 follow-up** for this case, scoring 3.
Repeating that same follow-up at tolerance1 gives **3 events and the shared
SLAP heavy pattern**, with full outputs in both directions.

The untouched fixed-cap100 comparison is 133 ties, six AAM-lower, one SLAP-lower.
The headline 134/6/0 includes the clearly labeled cap200 follow-up. These are
different protocols, not silently merged fixed-cap results.

### Case135 and the old successful record

The original saved record is:

```text
/h/399/yunhengzou/appendix_final/bgcp_rerank_latest/stages/TS_11/rp_stage.json
```

It used matching tolerance1, scoring threshold0.5, three seeds and cap100, and
already contained the same heavy mapping. Original endpoint coordinates and
element order match the current inputs exactly. Rescoring that saved mapping
on current WBOs gives 5 events. Three within-fragment WBO changes are approximately
0.518, 0.786 and 0.697: allowed at tolerance1, not at tolerance0.5. The new frozen
engine tolerance1 run also recovers this pattern. The earlier strict-run miss
was not evidence that AAM inherently could not find it.

## Family overlap, not just best score

Fixed cap100: 159 of 176 SLAP heavy classes represented; 12 excluded from the
saved families (cases11,59,64,123); five unresolved (cases4,6,7,31). Replacing
case123 with its separately saved cap200 result changes this to **160 represented,
11 excluded, five unresolved**. Thus equal best scores do not mean every SLAP
pattern is represented. These are score-preserving graph equivalence classes,
not a chirality equivalence claim.

At cap100, 97 cases have an extra AAM heavy pattern absent from saved SLAP outputs
at the lowest observed AAM event count. Positive witnesses are independently
checked against exact colored-graph certificates; negative results concern only
the saved search families. Query-budget exhaustion remains unresolved.

## Timing

Search timing below excludes queue wait and has no failed-run substitution:

| 280 fixed-cap100 directional calls | Tolerance0.5 | Tolerance1.0 |
|---|---:|---:|
| Total compute CPU, excluding measured persistence/loading | 17,933.75 s | 9,562.02 s |
| Mean directional elapsed, including saving/loading | 9.18 s | 4.57 s |
| Maximum directional elapsed, including saving/loading | 159.88 s | 141.18 s |

New search compute totals **2.66 CPU-hours**, averaging 68.30 CPU-seconds per
reaction for both directions. This is not a directly measured simultaneous
bidirectional latency. Node types differ: on the 120 calls sharing the same CPU
model across campaigns, summed compute CPU is 10,079.70 → 6,525.43 seconds.
That subset is not a randomized equal-load hardware experiment.

The separate cap200 follow-up costs another 34.06 compute CPU-seconds; directional
elapsed times are 0.26 and 8.90 seconds including saving/loading.
Witness comparison totals 78.24 summed case seconds. Bounded family-overlap
analysis totals 1,083.71 summed case seconds, parallelized by case. These
postprocessing sums are neither campaign wall time nor mapper-search CPU time.

## Reproduction and saved artifacts

Entry points: `bench/elementary_tolerance_ab.py` (prepare, submit, analyze,
summarize, cap_followup), existing comparison/family-query scripts, and
`bench/build_elementary_comparison_data.py` (data export) and
`tools/render_mapping_comparison.py` (shared original R/P/TS presentation). `evidence.tar.gz` contains compact
per-case scores/mappings, summaries and run provenance; full AAM checkpoints
remain on the cluster. Old tolerance0.5 artifacts were not overwritten.

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909/case123_cap200
/project/yunhengzou/coordinate_alignment/slap_sweep_worktree_20260910/reports/elementary140_tol1_20260909/viewer.html
```
