# Independent cut streams and smaller-first search

Two requested search-policy changes are tested on a fresh Golden campaign. Native fragment matching, fragment commitment, automorphism compression, and mechanism-independent search-graph representation are unchanged.

## Cut randomness

`rxn_core.aam.cut_seed` derives a deterministic 128-bit seed from the canonical cut edge set and the root seed using BLAKE2b. Different cuts no longer restart the same seed-order stream. Edge ordering, worker scheduling and Python's randomized hash do not affect the stream. The no-cut case retains seed 42 for continuity. The existing ordering preference for isolated repeated-element atoms is unchanged.

Cut checkpoint manifests now use `rxn_core.aam_checkpoints/v2` and identify this policy. Old-policy cuts cannot be resumed into a new-policy search. Old archives remain readable.

## Direction is orchestration, not a graph rewrite

```python
from rxn_core import plan_aam_search, search_aam

plan = plan_aam_search(input_problem, config)
result = search_aam(plan.problem, plan.config)
path = next(result.graph.paths())
input_R_to_P = plan.to_input_mapping(path.mapping)
```

`plan_aam_search` selects the endpoint with fewer **explicit atoms, including H**, as source. Ties retain R→P. It inverts supplied anchors when necessary, without mutating the input problem or configuration.

`AAMSearchPlan` retains the original input problem, the oriented search problem/configuration, and the direction. A raw `AAMResult` always describes its actual search source/target. Its DAG, fragment histories and correlated symmetry actions are **not inverted or reinterpreted in place**. Only concrete mappings are converted back after selecting/realizing a path. `search_aam` remains a directional primitive; callers such as augmentation workflows are not silently reversed.

The Golden adapter explicitly opts into this policy and stores `orientation.json` beside every original input and reference. Evaluating compressed groups happens in their native search orientation. Ranking, bond-event counts, reported product coverage and reference-completeness denominators retain the original input R/P interpretation. Reference annotations never influence search direction or seed choice.

## New Golden campaign

Settings retained from the previous full campaign: three seed orders, branch cap 100, tolerance 1.0, explicit H, no-cut plus all single-edge source cuts. The new source may be P, so its sweep edge set changes accordingly.

- Original dataset: 1,851 records; 1,760 have complete product-heavy reference annotations.
- Selection: 986 R→P, 865 P→R; no dataset exclusions or reagent removal.
- Fresh archives only; no old search reuse.
- CPU budget: 128 concurrent one-CPU reaction jobs. Twelve GiB memory per job.
- Search watchdog: 300 seconds. Scoring: 180-second process limit, with 120 seconds for the internal reference-query budget. Partial positive-only audit: at most 90 seconds. Slurm hard limit: ten minutes per reaction job.
- Timeouts, OOM events, partial searches and unresolved reference queries are recorded, not treated as ordinary mapping misses or dropped from denominators.
- Original inputs, references, seed policy, source hashes, per-cut JSON, full typed checkpoints, evaluation witnesses and timings are retained.

Full campaign: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`.

Pilot: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_pilot_20260907`.

The eleven-case pilot finished without interrupted searches. Eight recovered, including case 1665 at top-1 in 21.85 seconds and case 1740 in 0.46 seconds. Cases 602, 986 and 1285 still missed. This targeted pilot is not a benchmark accuracy estimate.
