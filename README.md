# rxn_core

Symmetry-aware WBO graph alignment and BGCP transition-state ranking.

The package discovers R-P mechanisms with sweep cut, aligns GT/IG transition
states through mechanism-local core atoms, scores normal modes, and writes a
self-contained HTML view per step.

## Install

```bash
pip install -e .
```

The installed CLI is:

```bash
rxn-core-pipeline --steps pr7.V.dodh_ts910 --inner-workers 10
```

## Inputs

The pipeline reads precomputed xtb cache directories. By default it looks in:

```text
data/xtb_frequency_calculations/<step>/
  R/                    reactant xyz + wbo
  P/                    product xyz + wbo
  sp_groundtruth/       GT TS xyz + wbo
  hess_groundtruth/     GT g98.out
  sp_iter<N>/           IG xyz + wbo
  hess_iter<N>/         IG g98.out
```

Use `BGCP_WORK=/path/to/xtb_frequency_calculations` when the cache lives
outside the repo.

## Outputs

```text
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
out/bgcp_alignment_eval.json
```

Generated caches, views, and paper artifacts are intentionally not part of the
main repository.

## Pipeline

1. `cut_sweep(...)` runs R-P mechanism discovery: no-cut plus one-edge R cuts
   above `BGCP_CUT_FLOOR`.
2. Mechanisms are deduped by symmetry-canonical broken/formed bond changes.
3. The reactive core is the atoms touching any broken or formed bond.
4. `ts_core_pool(...)` runs endpoint-to-TS matching from both R and P. P-side
   core mappings are pulled back through the R-P mechanism witness, then merged
   with R-side mappings in R-core indexing.
5. Each GT/IG candidate mapping is scored on the selected imaginary mode:

```text
S = beta * (1 + W_RXN * rho) * (1 + W_CORE * kappa) / n_imag^IMAG_PEN
```

where `beta` is bond-axis overlap, `rho` is reaction-coordinate overlap, and
`kappa` is the core-mode fraction.

## Repository Layout

```text
src/rxn_core/
  alignment/              molecule-level WBO alignment and sweep cut
  growth/                 connected-fragment growth and trace events
  matcher/                symmetry-compressed candidate matching
  chemistry_computations/ xyz, xtb cache, and frame helpers
  frag.py                 WBO graph and bond-change classification
  modes.py                Hessian parsing and mode-score features
  pipeline.py             BGCP full-view pipeline and CLI entry point

tests/                    focused unit tests
docs/ARCHITECTURE.md      module boundary notes
ALGORITHM.md              algorithm details
```

## Key Configuration

| name | default | meaning |
|---|---:|---|
| `BGCP_WORK` | `data/xtb_frequency_calculations` | xtb cache root |
| `BGCP_OUT_ROOT` | `out/bgcp_views` | per-step HTML/eval output |
| `BGCP_EVAL_JSON` | `out/bgcp_alignment_eval.json` | merged eval JSON |
| `BGCP_CUT_FLOOR` | `0.2` | R-P sweep cuts every R edge at or above this WBO |
| `BGCP_ISO_TOL` | `1.0` | WBO tolerance used by matching |
| `BGCP_VIEW_MAX_BRANCHES` | `100` | per-cut branch cap; capped cuts are discarded |
| `BGCP_TS_CORE_MAX_CANDIDATES` | `20000` | cap for mechanism-local TS/IG core mappings |
| `BGCP_TIMING` | `0` | set to `1` for per-target timing prints |

CLI scheduling options:

```text
--workers              total CPU budget in auto mode
--parallel-mode        auto | outer | inner
--inner-workers        explicit per-step inner worker count
--auto-inner-workers   target inner workers per concurrent step in auto mode
--steps                explicit cached step names
--limit                first N cached steps
```

In `auto` or `inner` mode, the per-step inner worker pool is used for both
expensive phases: R-P cut sweep and independent TS/IG endpoint core matching.

## Public API

Use top-level imports for stable pieces:

```python
from rxn_core import align_from_arrays, cut_sweep, ts_core_pool
from rxn_core import build_graph, classify_bonds
from rxn_core import parse_g98_modes, bond_overlap_per_mode
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for module boundaries and
[ALGORITHM.md](ALGORITHM.md) for the matching algorithm.
