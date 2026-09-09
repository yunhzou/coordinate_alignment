# Golden benchmark script guide

Results and exact cluster locations: [paper index](../reports/README.md).
Scripts remain at their original paths so saved commands/imports keep working.
No benchmark needs rerunning merely to inspect saved mappings or metrics.

## Maintained experiment entry points

| Script | Responsibility |
|---|---|
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
| `aam_small_graph_oracle.py` | Identical binary-graph exact-oracle control for AAM and SLAP |
| `audit_aam_slap_overlap.py` | Independent saved-witness and oracle score verification |

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
