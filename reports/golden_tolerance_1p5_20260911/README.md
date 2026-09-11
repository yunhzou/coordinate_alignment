# Golden: matching tolerance 1.0 versus 1.5

Fresh paired measurements on all 1,851 Golden cases. Frozen engine `98b01b1`, one seed, branch cap 100, explicit H, eight workers per directional call, original uncut/single-edge sweep, both directions. Only `iso_tolerance` changes; this affects fragment compatibility and symmetry grouping. Reference scoring, event thresholds, graph floor, seed orders, and branch cap stay fixed.

| Metric | Tolerance 1.0 | Tolerance 1.5 |
|---|---:|---:|
| Recovered reference | 1833/1851 (99.03%) | 1830/1851 (98.87%) |
| Not recovered from saved output | 17 | 19 |
| Unknown | 1 | 2 |
| Completed bidirectional searches | 1851 | 1851 |
| Mapping CPU on paired completed cases, seconds | 26928.816 | 26691.849 |
| Median paired case mapping CPU, seconds | 4.982 | 4.977 |
| Cases with branch-cap flags | 1370 | 1375 |

Paired completed cases: 1851. CPU ratio (1.5/1.0): 0.9912. CPU is parent plus child compute excluding measured archive persistence/loading; it already includes all eight workers. These are single paired measurements with alternating run order on the same host and CPU set. Offline evaluation and initial scheduler queue time are separate.

Golden input weights are discrete bond orders (1, 1.5, 2, 3). The added edge-weight compatibility is aromatic/triple (1.5 versus 3); both occur across endpoints in 281 cases. This input-level count ignores element and subgraph constraints. It differs from the continuous WBO setting of the TS holdout.

## Case changes

| Golden index (zero-based) | At 1.0 | At 1.5 |
|---|---|---|
| 812 | recovered | not_recovered |
| 1517 | recovered | not_recovered |
| 1599 | recovered | unknown |

### Targeted verification of case 1599

Only the saved tolerance-1.5 P-to-R archive was rechecked with a 240-second verification budget. It returned `not_recovered` after 78.93 seconds. Mapping was not rerun. After this extra check, trial recovery counts are {'recovered': 1830, 'not_recovered': 20, 'unknown': 1}. The primary table retains the original equal-budget results.

A not-recovered result describes the saved search families under this cap; it is not a proof that the algorithm cannot find the reference. Evaluator timeouts remain unknown. Coverage is recovery anywhere in the saved output, not top-ranked mapping accuracy.

The fresh 1.0 control differs from the archived 1.0 recovery status on indices: []. Both the historical control and the fresh paired control are retained in `case_metrics.json`.

All reported recovered witnesses passed element, injectivity, and fixed-reference certificate checks: {'tol_1p0': 3578, 'tol_1p5': 3577} directional witnesses. Frozen search/adapters and input/reference hashes were unchanged. This audit covers recovered witnesses; it does not enumerate every compressed family member.

Directional counts, caps, phase exceptions, and all changed-case witnesses are in the accompanying JSON files. The existing paper baseline and engine defaults are unchanged.

## Reproduction

From the paper checkout:

```bash
/project/yunhengzou/coordinate_alignment/.venv/bin/python bench/golden_tolerance_benchmark.py prepare --run /path/to/new_run
/project/yunhengzou/coordinate_alignment/.venv/bin/python bench/golden_tolerance_benchmark.py submit --run /path/to/new_run
# After both phases finish:
PYTHONPATH=/path/to/new_run/original/src:/path/to/new_run/engine/bench /project/yunhengzou/coordinate_alignment/.venv/bin/python bench/publish_golden_tolerance.py --run /path/to/new_run --destination /path/to/new_report
```

Raw archives and logs: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_tol15_20260911`.
