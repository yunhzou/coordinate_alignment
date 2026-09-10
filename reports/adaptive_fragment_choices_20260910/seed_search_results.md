# Adaptive sequential seed choices: recovery and remaining limits

**At 400 growth operations, all four cases match the saved full-sweep minimum
event count and recover all eight best-score heavy-mapping classes under
element-compatible score equivalence.** This is not equality of the complete
labelled result sets or proof of performance on the full benchmark.

The mature production pipeline remains unchanged. This is an isolated scheduler
above the existing fragment matcher; no new chemistry rules, native growth
algorithm changes, representative-based pruning or permutation enumeration were
introduced in this phase.

## Measured result

Same fixed four cases as the preceding closure pilot. One CPU per case on
bosque12, fixed RNG 42, explicit H, matching tolerance 1.0, event threshold 0.5,
per-growth branch cap 100. No sweep cuts or historical seed orders enter the
blind adaptive search. All four use the same 400-operation checkpoint.

| Case | Historical full-sweep best events | Adaptive best events | Historical best heavy classes recovered | Compute CPU s |
|---|---:|---:|---:|---:|
| Pd #25, R→P | 5 | 5 | 4/4 | 6.277 |
| Ni TS11 #76, P→R | 3 | 3 | 1/1 | 4.075 |
| Ni TS14 #77, P→R | 8 | 8 | 1/1 | 1.185 |
| V #114, R→P | 4 | 4 | 2/2 | 0.228 |

Compute includes setup, matching and snapshot symmetry finalization. Saving,
queue time, imports and separate pattern evaluation are excluded. These are
single observations, not a repeated-trial speed distribution or a same-hardware
paired speedup against the old full-sweep run. Cumulative compute CPU totals
11.765 seconds across the four one-CPU searches; do not interpret that sum as
parallel elapsed time.

Pd encountered native growth caps; the other three did not at this checkpoint.
Pending alternatives are recorded. Cap 100 is **per growth**, not a global cap
on this experimental agenda. We do not claim exhaustive search or a proven
globally minimum event count.

The vanadium result also contains both literal historical full-H event
signatures, not just matching scores: their saved 400-operation witness ordinals
are 60 and 4. This is agreement with historical outputs, not independently
annotated chemical ground truth.

## The actual root cause and scheduler change

The historical cut sweep also changed the random seed order per cut. Taking one
saved successful Ni TS14 cut run from each of ten seed-index folders, **8/10
still recovered eight events with their order preserved and the cut removed**.
Those order-guided replays were diagnostic only; their seed orders were not
fed into the blind experiment.

The successful Ni path grows 15-, 30-, then 53-atom fragments. The original
fixed-order adaptive pilot first locked a 56-atom region crossing those
boundaries. Earlier closure choices alone did not efficiently discover the
useful order.

The new scheduler:

1. Grows the normal next fragment under the current locked mapping, islands
   and deferred-edge constraints, using the existing `match_fragment` kernel.
2. Keeps one lazy alternative-seed task per exact search state. Source regions
   already grown guide ordering; no covered-region seed is discarded.
3. Finishes **one** mapping, then yields to another seed choice. Ordinary sibling
   branches remain pending instead of forcing an entire subtree to finish.
4. Shares alternative work across depths and between new/previously explored
   regions. This prevents either category from starving the other.
5. Reuses exact equal states while retaining incoming correlated symmetry paths
   in the AAM search DAG.

The preceding depth-only, novelty-only and whole-subtree-continuation policies
are preserved as negative ablations. They were cheaper than old full sweep but
did not reliably recover the difficult result. More work on those policies was
not sufficient; changing scheduling was necessary in these observations.

## Why the pattern counts need careful interpretation

Two comparison definitions are saved side by side. Neither changes AAM search,
the atom-element matching constraint, the event formula or the 0.5 threshold.

**Original, finer labels:** each bond's symmetry label encoded its response
against every opposite-side bond weight, regardless of endpoint elements.
This can distinguish C-C bonds because of thresholds against C-H weights, even
though element-preserving AAM cannot map C-C to C-H.

**Element-compatible score equivalence:** compare a bond only with opposite-side
bond weights having the same unordered element pair. Keep atom elements,
adjacency and every applicable event response. Endpoint automorphisms of these
colored graphs preserve event classification for every admissible atom map:
an R bond can meet only the corresponding element-pair weights in P, and
adjacency preservation also retains bond absence/formation behavior. The same
definition is applied to both historical and adaptive outputs.

This is an exact symmetry of the **specified scoring problem**, not exact
floating-point WBO identity, geometry, electronic environment or chemical
mechanism identity. A generator-based test checks all three event-count
components under sampled element-compatible mappings. Every recovered class has
a saved full-H witness validated for fragment compatibility and the required
event bound. Correlated family constraints are retained; permutations are not
enumerated.

At the larger saved checkpoints, the original finer comparison recovered:

| Case | Historical finer classes at best score | Certified recovered | Evaluation status |
|---|---:|---:|---|
| Pd | 546 | 64 | 60-second diagnostic budget reached |
| Ni TS11 | 48 | 48 | All found |
| Ni TS14 | 17 | 8 | 60-second diagnostic budget reached |
| V | 8 | 8 | All found |

Those results are **not** being relabelled as full labelled parity. The table at
the top uses the separately justified, coarser scoring symmetry. All finer
counts and witnesses remain saved. Heavy-relation classes do not distinguish
every possible hydrogen-transfer pattern, although explicit H enters both
matching and event scoring. Higher-event alternatives and the full benchmark
remain unvalidated for replacement of the mature pipeline.

## Validation, bounds and persistence

- Final scheduler: **489 tests passed**, 24 multiprocessing/fork deprecation
  warnings, 99.06 seconds; job 457201. A subsequently added score-equivalence
  test also passed, along with all six scheduler tests.
- 24 saved checkpoints from the accepted one-route and balanced runs passed
  DAG, element/injectivity, locked-state extension and preserved-bond checks.
- All workers retain a 300-second external timeout, five-second kill grace and
  ten-minute Slurm limit; automatic requeue is disabled. No watchdog kill
  occurred. Work budgets are finite. The larger Pd run stopped at the 30-second
  compute soft limit (2,732 operations); finalization brought compute to 30.618 s.
- Mapping and full compressed graphs are saved at every budget. Further
  pattern analysis used these checkpoints, not new AAM searches.
- This scheduler is resumable in memory. Archives retain results and pending
  descriptors, not a complete cross-process execution checkpoint.

`seed_results.json` contains the ablations, source/input/native hashes, full
positive pattern witnesses, validations and tests. Rebuild that evidence bundle
with `bench/publish_adaptive_seed_results.py`. The two comparison definitions
are selected explicitly by `--element-pair` in `compare_adaptive_patterns.py`
and `query_adaptive_patterns.py`; `--checkpoint-label` evaluates an earlier
saved snapshot without rerunning mapping.

Full accepted-run paths:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_seed_route_20260910_4bbcer
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_seed_balanced_20260910_Nogs0S
```

Diagnosis paths and earlier negative ablations are listed in
`seed_order_diagnosis.md` and `seed_results.json`. Mapping jobs: 457077, 457086,
457096, 457143, 457192. Final validations: 457200; original pattern comparisons:
457198; family queries: 457199 and 457212; element-pair comparisons/queries:
457219/457220; 400-operation snapshot comparisons/queries: 457231/457232.

The first historical-witness comparison used the wrong orientation for P→R
exports. It was corrected by explicitly inverting historical R→P-normalized
witnesses; the initial diagnostic JSON/logs were retained and are not used in
the reported counts. Element compatibility is now asserted during comparison.

**Conclusion:** minimum-score parity and parity of best-score heavy classes
under the stated scoring symmetry are demonstrated on these four cases. Full
labelled-output parity, all hydrogen alternatives, higher-event spectrum parity
and full-dataset performance are not demonstrated. Production defaults are
unchanged.
