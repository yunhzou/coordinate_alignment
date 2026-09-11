# Golden benchmark script guide

Results and exact cluster locations: [paper index](../reports/README.md).
Scripts remain at their original paths so saved commands/imports keep working.
No benchmark needs rerunning merely to inspect saved mappings or metrics.

Use the [two shared viewer styles](../docs/VIEWERS.md). R/P and TS comparisons
use the restored original white-panel renderer; benchmark-specific skins are
retired. `missing_pattern_display_data.py` prepares the three-case witnesses and
`tools/render_mapping_comparison.py` renders them through the shared module.

## Maintained experiment entry points

| Script | Responsibility |
|---|---|
| `holdout_missing_pattern_seeds.py` / `publish_missing_pattern_seeds.py` | Targeted seed-budget diagnostic on cases 11, 64 and 101; frozen search, exact target recovery and paired CPU |
| `certify_heavy_event_exclusion.py` | Necessary-condition exclusion using invariant broken/formed heavy-bond patterns, checked against positive controls |
| `holdout_forward_analysis.py` | Forward-only elementary-step scores, saved directional timings and alternative-family checks |
| `publish_holdout_forward_analysis.py` | Validate forward provenance; publish score, CPU and alternative comparisons with an offline viewer |
| `holdout_minimum_events.py` | Full-H minimum-event snapshots and joint symmetry canonicalization on saved holdout outputs |
| `holdout_minimum_event_families.py` | Correlated AAM path / SLAP label-family membership, with explicit unresolved status |
| `enumerate_slap_minimum_events.py` / `enumerate_aam_minimum_events.py` | Event-pattern enumeration at recorded benchmark minima; completeness and lower-bound records |
| `complete_minimum_event_checks.py` | Targeted longer SLAP queries and AAM membership checks for newly exposed patterns |
| `validate_minimum_event_analysis.py` | Exhaustive finite checks on small real output families |
| `publish_minimum_event_analysis.py` | Independently rescore witnesses; generate the report, figures and offline bond-edit viewer |
| `conditioned_reuse_pilot.py` | Frozen-engine, same-node conditioning/growth-reuse ablations; complete graph checks, recoverable timers and saved outputs |
| `public_reuse_pilot.py` | Public AAM reference/reused-native integration, parent+worker CPU accounting and full raw/final checkpoints |
| `publish_conditioned_reuse.py` | Verify and package those saved experiments, logs, source hashes and tests without rerunning matching |
| `prepare_golden_benchmark.py` | Dataset preparation and pinned inputs |
| `golden_publication.py` | Fixed-policy search, analysis and mode comparison |
| `collect_golden_patterns.py` | Concrete pattern extraction from saved compressed AAM archives |
| `benchmark_pattern_collection.py` | Collection validation/reporting |
| `golden_evaluation.py` | Shared reference-mapping evaluation |
| `golden_competitors.py` | Released mapper adapters, workers and default report |
| `golden_slap_budget.py` | Fixed input-order plans, expanded SLAP evaluation, summary and retry plans |
| `elementary_feasibility.py` | Native XYZ/WBO and SLAP XYZ holdout; reference-free feasibility, saved mappings, timing and reporting |
| `compare_elementary_outputs.py` | Saved native mapping comparison, common WBO scores, symmetry checks and symbolic SLAP hydrogen refinement |
| `elementary_family_overlap.py` | Saved heavy-family membership and low-edit alternative witnesses |
| `elementary_tolerance_ab.py` | Frozen-engine tolerance1 holdout rerun, separate cap follow-up, comparisons and timing |
| `slap_guided_pilot.py` / `slap_local_workflow.py` | Exploratory cheap proposals and targeted cuts; not the default AAM pipeline |
| `aam_small_graph_oracle.py` | Internal synthetic diagnostic only; retained, excluded from paper evidence |
| `audit_aam_slap_overlap.py` | Independent saved-witness and oracle score verification |
| `cut_replay_pilot.py` | Opt-in exact cut replay, paired same-node ablations, full compressed archives and saved-mechanism checks |
| `cut_replay_kernel_probe.py` | Isolated pre/post native-binary output and timing regression probe |
| `profile_fragment_pipeline.py` | Per-cut saved graphs and function/phase diagnostics for Python/native scheduling |
| `audit_ts_weight_information.py` | Real cached TS fractional-cost and output-information audit; no mapping accuracy benchmark |
| `compare_real_ts_mappings.py` | Blind real endpoint-to-TS mapping, five input/engine configurations, native archives and historical-core membership queries |
| `score_real_ts_families.py` | Common all-H scoring of saved AAM paths / SLAP labels plus LAP fingerprints, without mapper reruns |
| `publish_real_ts_comparison.py` | Verify and package those saved results and generate an offline mapping viewer |

The other Golden scripts record diagnostic experiments (seeds, caps, direction,
anchors, ranking, memory, viewers). They are retained for provenance, not silently
included in the fixed-policy paper protocol. No scripts or evidence were moved
or deleted during organization.

## Reproduce safely

1. Select a run through the paper index and inspect its manifest, frozen engine,
   environment, inputs and submission records. Preserve dataset/source pins and
   random seeds; the exact historical settings take precedence over current defaults.
2. Use a new output directory for any fresh experiment. Never initialize or rerun
   workers into a published run merely to inspect its results.
3. Reuse saved mappings for evaluation/figures. Write new analyses separately and
   retain evaluator/source identity. Include missing/invalid outcomes in denominators.
4. Retain intermediate/final mappings, attempt logs, per-case timing, retries and
   accounting along with the report. Push compact evidence; retain large cluster
   banks at the documented paths.

Inspect current CLI options without starting a run:

```bash
.venv/bin/python bench/golden_publication.py --help
.venv/bin/python bench/collect_golden_patterns.py --help
.venv/bin/python bench/golden_competitors.py --help
.venv/bin/python bench/golden_slap_budget.py --help
```

The competitor execution environment is separately pinned at:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/competitor_env_20260908/bin/python
```

Evaluation/adapter regression checks (not a fresh chemistry benchmark):

```bash
.venv/bin/python -m pytest -q tests/test_golden_evaluation.py tests/test_golden_competitors.py tests/test_golden_slap_budget.py
```
