# Alignment Architecture

This document describes the module boundaries for the WBO graph alignment
implementation. The main rule is that each layer owns one concept and exports
only the abstraction needed by the layer above it.

## Public Surface

`rxn_core.__init__` is the convenience API used by scripts and the pipeline.
It re-exports stable entry points such as:

- `align_from_arrays`
- `match_wbo_graphs`
- `find_islands`
- `build_graph`
- `classify_bonds`

Code outside the package should prefer these top-level imports unless it is
debugging internals.

## Molecule Alignment: `rxn_core.alignment`

Purpose: turn two WBO matrices into scored molecule-level alignments.

Files:

- `alignment/api.py`: public dataclasses and entry points:
  `MatchCandidate`, `MatchResult`, `match_wbo_graphs`,
  `align_from_arrays`, `analyze_alignment`, `cut_edges_above_floor`.
- `alignment/branch.py`: molecule-level branch scheduling, branch dedupe,
  seed-order generation, and final symmetry repair.
- `alignment/sweep.py`: R-P sweep-cut mechanism discovery via no-cut plus
  one-edge R cuts, with symmetry-canonical mechanism signatures.  Parallel
  cut-sweep work is grouped by cut, so each worker computes the cut-specific R
  orbit map once and reuses it for all seed orders; product orbits are
  invariant and are computed once in the worker initializer.
- `alignment/ts_core.py`: mechanism-local endpoint->TS/IG core mapping used by
  the ranker after R-P mechanisms are known. The BGCP pipeline runs it from
  both R and P, pulls P-derived mappings back through the R-P witness, and
  scores the union in R-core indexing.

This layer owns whole-molecule decisions: seed orders, branch lifecycle,
mechanism discovery, mechanism-local TS core matching, and final representative
choice. It does not implement single-fragment growth or symmetry-block witness
search.

## Fragment Growth: `rxn_core.growth`

Purpose: grow one connected R fragment against a target graph.

Files:

- `growth/island.py`: the live island-growth loop. It initializes seed
  candidates, pops frontier edges, calls the matcher extension rule, commits
  successful atoms, records deferred boundary edges, and locks the island only
  after heap saturation.
- `growth/frontier.py`: heap/frontier utilities and final uniqueness checks.
- `growth/result.py`: `_IsoResult`, the locked fragment result.
- `growth/trace.py`: diagnostic event formatting and failure explanations.

This layer owns traversal order and locking. It does not decide whether a
candidate target atom is valid; validity is delegated to `rxn_core.matcher`.

## Symmetry Matching: `rxn_core.matcher`

Purpose: represent and extend symmetry-compressed partial mappings.

Files:

- `matcher/state.py`: `_SymCand`, `_SymBlock`, witness/materialization helpers.
- `matcher/extend.py`: one-atom extension for a list of compressed candidates.
- `matcher/support.py`: witness search inside unresolved symmetry blocks.
- `matcher/dedupe.py`: orbit plus one-hop-boundary dedupe signatures.
- `matcher/orbits.py`: exact pynauty orbit grouping on 0.2-tolerance
  WBO-colored graphs, plus the old WL/color-refinement helper for fallback
  and comparison.
- `matcher/primitives.py`: WBO access and tolerance helpers.
- `matcher/chemistry.py`: post-hoc chemistry-relevant symmetry expansion.

This layer owns the central compression model:

```text
_SymCand = witness mapping + interchangeable/correlated _SymBlock pools
```

The witness is one deterministic assignment. The blocks encode the actual
symmetry state. A block with two R atoms and four P atoms represents
`P(4, 2)` injective concrete assignments without enumerating them.

## Chemistry And IO Helpers

- `chemistry_computations/`: file and external-computation utilities.
  - `xyz.py`: XYZ parsing and formatting, including extended XYZ mode output.
  - `xtb.py`: cached xtb single-point/Hessian execution and WBO cache loading.
  - `frames.py`: coordinate reindexing into the R atom frame.
- `frag.py`: WBO graph construction, `classify_bonds`, and local witness
  materialization. It re-exports old XYZ/xtb names only for compatibility.
- `modes.py`: normal-mode parsing and score-feature vectors.
- `align.py`: compatibility facade for cache loading and coordinate reindex
  helpers now owned by `chemistry_computations`.
- `pipeline.py`: BGCP cached orchestration, split into explicit resumable
  stages. `run_rp_stage` is the reusable R-P alignment / mechanism-discovery
  entry point and writes `rp_stage.json`. `run_ts_stage` consumes those
  mechanisms plus GT/IG/TS targets, runs mechanism-local R/P endpoint
  `ts_core_pool`, and writes `ts_stage.json`. `write_view_stage` is
  presentation-only: it writes the multi-mechanism viewer plus slim eval JSON
  from stage records. `write_rp_alignment_files` is a Stage 1 export helper
  for downstream NEB/path setup: it writes per-mechanism `R.xyz`,
  `P_aligned.xyz`, `neb_endpoints.xyz`, mapping CSV, and metadata without
  doing any spatial/Kabsch fitting. `write_ts_alignment_files` exports the
  Stage 2 selected best-S GT/IG/TS core mapping as native target XYZ plus an
  R-frame core-aligned materialization and picked-mode extended XYZ. The
  TS selector prefers exact core mappings seen from both R and P endpoints
  when such endpoint-consensus candidates exist, then falls back to score-only
  selection for one-endpoint pools.
  compatibility `process_step` composes all three.
  The CLI exposes the same split through `--stage rp|ts|view|full`, and
  `--mechanism` restricts Stage 2 verification to selected mechanism IDs.
  The default auto scheduler treats `--workers` as a total CPU budget and
  splits it between concurrent steps and each step's inner worker pool. That
  inner pool is reused after R-P discovery for independent TS/IG endpoint
  core-matching tasks across targets, mechanisms, and R/P endpoints. The
  pipeline prefers existing xtb caches but, in `BGCP_XTB_MODE=auto`, can fill
  missing `wbo` and `g98.out` files from available XYZ inputs before
  alignment. Each individual xtb subprocess is capped by
  `BGCP_XTB_MAX_THREADS=8` through `OMP_NUM_THREADS`; target cache filling has
  its own `BGCP_XTB_WORKERS` pool so multi-threaded xtb jobs can be balanced
  separately from R-P/TS alignment workers. Each cache endpoint is loaded as one complete
  molecule/complex graph; the pipeline does not assemble separate
  reactant/product fragments or merge independent WBO matrices.

These modules should not depend on matcher internals unless they are explicitly
running alignment.

## Dependency Direction

The intended direction is:

```text
rxn_core.__init__
  -> alignment
      -> growth
          -> matcher
      -> frag
  -> chemistry_computations / modes / align / pipeline helpers
```

Lower layers must not import higher layers. In particular:

- `matcher` must not import `growth` or `alignment`.
- `growth` must not import `alignment`.
- `alignment` may import `growth`, `matcher`, `frag`, and
  `chemistry_computations`.
- scripts may import public top-level names, or internal packages only when
  producing debug traces.

## Naming

The old `pq.py` module and related compatibility facades were removed. The
priority queue is now an implementation detail of `growth/island.py`, not a
package-level abstraction.
