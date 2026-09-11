# Paper benchmark index

Start here when querying results or preparing the paper. Existing reports and
raw runs are preserved in place; this index does not merge experimental protocols.

## Primary evidence

[Forward-only elementary steps: scores, alternatives and mapping time](holdout_forward_cap1000_seed1_20260910/README.md)
— reactant → product only, as requested for the elementary-step comparison.
One-seed, cap-1000 AAM retains every bidirectional best score. It has four lower
scores and 136 ties against forward SLAP+sweep, with mean mapping CPU of 6.733 s
versus 4.955 s. The report includes direction-specific alternative membership,
completeness flags and an offline bond-edit viewer. Earlier bidirectional
reports below remain separate protocol records.

[140 cases: alternative bond-edit patterns at tied minimum scores](holdout_minimum_event_patterns_20260910/README.md)
— full-H event sets with joint symmetry normalization and an offline bond viewer.
On 138 AAM / SLAP-sweep score ties, all 204 verified AAM patterns are shared;
SLAP sweep has three additional patterns excluded from the saved AAM families
(cases 11, 64 and 101). SLAP catalogues are complete; AAM enumeration is complete
in 96/140 cases and otherwise gives lower bounds. Scores are held at the preceding
benchmark minima. This is candidate bond-edit coverage, not pathway validation.

[140 cases: one-seed cap-1000 AAM versus native SLAP and SLAP with cuts](holdout_slap_xyz_sweep_20260910/README.md)
— fresh mapping and CPU comparison on matched hosts. Against the native XYZ
single-edge sweep, AAM has two lower minimum-change scores and 138 ties;
mean mapping CPU is 15.989 s for AAM and 9.441 s for SLAP with cuts.
All 16,748 sweep calls completed. Scores are not reference-mapping accuracy.

[One-seed cap-1000 AAM versus native XYZ SLAP](holdout_cap1000_seed1_20260910/README.md)
— the paired baseline underlying that sweep comparison: six AAM-lower scores,
134 ties, all 140 mappings available. Separate from the manuscript's fixed
ten-seed/cap-100 panels; the frozen stable engine remains unchanged.

[Exact conditioned symmetry and fragment reuse](conditioned_reuse_20260910/README.md)
— full ten-seed stress ablation: 1.58x compute speedup, 24,840 exact graph checks;
opt-in public backend, unchanged seed/cut policy, rejected slower extension cache.

[Concrete AAM / SLAP comparison on real endpoint-to-TS inputs](real_ts_comparison_20260909/README.md)
— six pairs from three cached cases; actual saved alternatives, full-H event
scores, binary/continuous ablation and offline viewer. Exploratory, not curated TS accuracy.

[Output information and continuous weights for TS structures](ts_weight_information_20260909/README.md)
— interface audit and actual cached TS weights; no TS mapping-accuracy claim.

[Dependency repair and native AAM scheduling experiment](dependency_fragment_repair_20260909/README.md)
— unchanged seed/cut policy, exact compressed-result comparisons and paired
timings; 1.24x overall in the tested stress workload, experimental dispatch only.

[Exact incremental cut-replay experiment](incremental_cut_replay_20260909/README.md)
— same-node paired timings and exact compressed-result comparisons; opt-in only,
not a change to the publication search protocol.

[Exploratory SLAP-guided local-search experiments](slap_guided_pilot_20260909/README.md)
— separate from fixed-policy publication results; best-score agreement is not
full-family or chemical accuracy.

[140-step tolerance1 rerun, matched-cap comparison and full viewer](elementary140_tol1_20260909/README.md).

[AAM versus SLAP: theory and saved-family overlap](aam_slap_theory_overlap_20260909/README.md).
The synthetic exact-oracle section is retained as an internal diagnostic only;
it is excluded from the paper evidence by the user's decision (2026-09-09).

Separate from Golden: [140 elementary-step native-input feasibility holdout](elementary140_feasibility_20260908/README.md).
No reference atom mappings are available; this is not an accuracy benchmark.

| Experiment | Report and machine-readable results | Role |
|---|---|---|
| Fixed-policy AAM search | [Report](golden_publication_20260908/README.md), [per case](golden_publication_20260908/per_case.csv), [per direction](golden_publication_20260908/per_direction.csv), [analytics](golden_publication_20260908/analytics.json) | Original compressed archives; single and bidirectional search; seed 10, cap 100, tolerance 1.0, sweep cut |
| Concrete mapping collection | [Report](golden_pattern_benchmark_20260908/README.md), [metrics](golden_pattern_benchmark_20260908/benchmark_report.json), [per case](golden_pattern_benchmark_20260908/per_case.json) | Reuses the preceding search; concrete reference recovery 1839/1851 (99.35%) within +6 bond events, not top-k |
| Default competitors | [Report](golden_competitors_20260908/README.md), [summary](golden_competitors_20260908/summary.json), [saved outputs](golden_competitors_20260908/saved_outputs.tar.gz) | Released mapper configurations evaluated under the common mapping criterion |
| Expanded SLAP alternatives | [Report](golden_slap_budget_20260908/README.md), [main results](golden_slap_budget_20260908/expanded_summary.json), [all-atom results](golden_slap_budget_20260908/all_atoms_summary.json), [combined analysis](golden_slap_budget_20260908/combined_analysis.json), [case evaluations and provenance](golden_slap_budget_20260908/evaluations_and_provenance.tar.gz) | Both directions, binary/weighted, input-order diversity; combined recovery 1691/1851 (91.36%) |

All percentages above use all 1851 original Golden records. Earlier diagnostic
reports remain useful for root-cause analysis, but are not additional independent
paper test sets and must not be unioned into these fixed-policy results.

## Full raw artifacts on the cluster

Git contains reports, case tables, selected viewers, and compact output/provenance
archives. Large AAM archives and expanded candidate banks remain on the cluster:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_publication_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_pattern_collection_20260908_v2
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_competitors_full_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_slap_expanded_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_slap_all_atoms_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_cap1000_seed1_20260910
/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_slap_xyz_sweep_20260910
/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_minimum_event_patterns_20260910
/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_forward_cap1000_seed1_20260910
```

These cluster copies are not a claim of independent/offsite backup. Do not delete
them after pushing Git. Preserve failed attempts, final retry outcomes, frozen
engines, input plans, manifests, environment/native hashes, and Slurm accounting.

## Find an individual reaction without rerunning search

Case identifiers are zero-based Golden indices, not atom labels or SLAP attempt
indices. Use the saved inputs to resolve chemical structures.

- AAM: filter `per_case.csv` on `index`; use `per_direction.csv` for direction-specific
  timing/status. Original search files are under `directions/<index>/<direction>/`.
- Collected mappings: `results/<index>/<direction>/` in the collection run contains
  `candidates.json`, `summary.json`, `patterns.jsonl.gz`, and `paths.jsonl.gz`.
  Patterns retain physical mappings, events and source paths; paths retain
  traversal/completion information. The original compressed archive preserves
  alternatives not extracted within budget.
- Default competitors: `<method>/<index>.json` in the default run contains raw
  output/status/timing. Use the saved common evaluation, not native cost alone.
- Expanded SLAP: `evaluation_<method>_<shard>.json` contains `cases` with `case`,
  `recovered`, and `attempts`. Resolve each attempt using `plans.json` before
  opening `<method>/<attempt-index>.json`; attempt index is not reaction index.
- Example collected-mapping viewer:
  `/h/399/yunhengzou/coordinate_alignment/reports/golden_pattern_benchmark_20260908/case1033/viewer.html`.
  Git stores the compressed `.html.gz` version; decompress before opening.

## Definitions that must accompany paper claims

- Recovery means the exact reference heavy-atom mapping relation, including
  unmatched atoms, modulo endpoint chemical symmetry. It does not measure hydrogen
  provenance or establish that a reference is the only chemically possible answer.
- AAM includes explicit H; search tolerance is 1.0, event scoring tolerance 0.5.
  Concrete event-window recovery, compressed-family feasibility and first-output
  accuracy are different metrics. Collection is not globally exhaustive.
- SLAP has native alternatives. Binary AAM is its default; weighted mode preserves
  bond orders. Our expanded SLAP union is not its default single-call accuracy.
- Timeouts, missing archives and invalid outputs remain in accuracy denominators.
  Reference checks are evaluation, not part of blind candidate generation.

## Timing: query stages separately

For AAM, use search compute, baseline ranking/merge, then additional collection
compute. Collection reused prior rankings: it is neither a fresh full search nor
a replacement timing for the entire pipeline. Record loading, persistence and
reference verification separately. Parallel retry parent-only CPU fields require
the inclusive CPU totals in `slurm_accounting.tsv` (see collection report).

Default competitor mapping CPU excludes startup and external I/O. Expanded SLAP
per-mode CPU sums omit killed calls without a final CPU reading; whole-worker
Slurm accounting includes their cost, startup and retries. The measured combined
expanded SLAP worker total is 20.735 CPU-hours, including both symmetry experiments.

Queue time, failed attempts, successful retry compute, summed worker wall time,
and actual campaign elapsed time are distinct quantities. Keep both clean
successful-compute and total-resource costs. Do not use summed direction wall
time as bidirectional latency or reconstructed scheduling as measured wall time.
The conversational 6–7x comparison is approximate and mixes accounting scopes;
it is not a publication-ready matched-budget speedup claim.

## Reproduction and historical work

See the [script guide](../bench/README.md). Use frozen scripts and exact commands
from each run's manifests/submission records for historical reproduction; current
repository code can change. All other `golden_*` reports here are exploratory or
stage-specific history unless explicitly identified above. In particular, the
SLAP `golden_slap_budget_20260908` cluster pilot is not the completed expanded run.
