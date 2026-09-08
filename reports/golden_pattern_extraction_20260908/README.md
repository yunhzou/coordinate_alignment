# Patterns from compressed AAM families

The existing alternatives viewer now retains its saved representatives and adds
distinct explicit-atom mapping patterns extracted from compressed paths. The AAM
search engine and frozen publication archives are unchanged. Extraction is
reference-blind; the reference is used only to label the resulting mappings.

## Case 1033, P → R archive

- Previous viewer: 113 representatives; none reference-equivalent.
- Extended viewer: 119 candidates, including six new strict mapping patterns.
- A new reference-equivalent candidate appears at rank 7: three bonds broken,
  three formed, zero bond-order changes. Scores include explicit hydrogen and
  use event tolerance 0.5; archive matching tolerance was 1.0.
- Initial extraction: 15.1 seconds on one path. Resumed pass: 30.4 seconds,
  27 paths visited cumulatively. This is extraction wall time including small
  checkpoint writes, not AAM search time. The resumed viewer generation took
  52.9 seconds including loading, collection, rendering and file writes.
- Extraction remains **incomplete**. None of these 27 path exhaustiveness checks
  finished within the allotted budgets. Finding witnesses works; proving that
  every remaining exact-symmetry class is excluded is still a bottleneck.
  This experiment does not establish an exhaustive or extremely fast collector.

The original archive retains every unexplored family. The checkpoint stores
concrete witnesses, symmetry action traces, terminal/path origins and per-path
completion status. Equal scores are not merged. Exact supplied endpoint graph
symmetries are quotiented using coupled group actions, not independent atom
choices or enumeration of all bijections. Chemical tags refine graph identity;
this is not a new stereochemical feasibility filter or resonance normalization.

## Open the viewer

Full standalone HTML on the cluster:

`/h/399/yunhengzou/coordinate_alignment/reports/golden_pattern_extraction_20260908/case1033/viewer.html`

For GitHub/download, `case1033/viewer.html.gz` contains the same HTML. Decompress
before opening. It retains the clickable rank/event plot, R/P fragment drawings,
atom correspondence inspector, and explicit-hydrogen toggle. The large plain
HTML is excluded from Git; mappings and resume checkpoints are committed.

## Reproduce / resume

```sh
PYTHONPATH=src:bench OPENBLAS_NUM_THREADS=1 .venv/bin/python bench/view_golden_alternatives.py \
  --source /project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_publication_20260908 \
  --direction P_to_R \
  --audit data/aam_benchmarks/golden_original_20260906/audit.jsonl \
  --indices 1033 --output reports/golden_pattern_extraction_20260908 \
  --pattern-seconds 30 --path-seconds 1
```

These are explicit optional postprocessing time budgets, not AAM branch caps or
candidate-count caps. Unfinished queries remain unresolved, not rejected.
Resume validates archive SHA-256, orientation and endpoint annotations.

Validation: 81 focused tests passed, covering the shared family query/scoring
layer, graph histories, mapping viewers, tiny exhaustive pattern oracles,
whole-class exclusion, same-score distinct patterns, explicit-H invariance,
partial mappings, checkpoint resumption and archive provenance.
