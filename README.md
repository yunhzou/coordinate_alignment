# rxn_core

Pipeline and library for symmetry-aware WBO graph alignment.  The current
default CLI builds BGCP full views from precomputed xtb caches: it discovers
R-P mechanisms with sweep cut, scores GT/IG transition-state guesses under
each mechanism, and writes an interactive multi-mechanism HTML view.

The signal throughout is **Wiberg bond order** from xtb GFN2 outputs: WBO
matrices on R, P, GT, and every IG-TS, plus normal modes parsed from xtb
Hessians.

## Pipeline overview

```
Inputs: cached BGCP xtb step                    Outputs:
                                                out/bgcp_views/<step>/
  appendix_perparation/
    xtb_frequency_calculations/<step>/ ─────►     view.html
      R/, P/                                      _eval_v2_slim.json
      sp_groundtruth/, hess_groundtruth/       out/bgcp_alignment_eval_v2.json
      sp_iter<N>/, hess_iter<N>/
```

The pipeline:

1. **R-P mechanism discovery** — `cut_sweep(...)` runs no-cut plus every
   one-edge R cut above `BGCP_CUT_FLOOR`, then dedupes by
   symmetry-canonical broken/formed bond changes.
2. **Reactive core per mechanism** — atoms touching any broken or formed
   bond define the mechanism-local core.
3. **GT/IG core matching** — `ts_core_pool(...)` enumerates distinct core
   atom mappings from both endpoints: R→TS and P→TS. Product-side candidates
   are pulled back through the R-P mechanism witness and unioned in R-core
   indexing. Spectator atoms are not expanded into full bijections.
4. **Per candidate TS/IG mapping**: reindex modes to R-frame and compute:
   - β (bond_overlap): mode displacement projected on the bond-axis
     stretch/contract direction at TS coordinates
   - ρ (rxn_overlap): mode core-atom motion projected on the R→P
     direction
   - κ (core_fraction): fraction of mode energy on core atoms
5. **Pick the imaginary mode** with max β.
6. **Score**:
   ```
   S = β · (1 + w_r · ρ) · (1 + w_c · κ) / n_imag^p
       (w_r = 1.0,  w_c = 0.2,  p = 0.3)
   ```
7. **Choose the core mapping** that gives the highest S for GT and each IG.
8. **Render** one page per step with mechanism buttons, R/P/GT, and all IGs.

## Layout

```
rxn_core/
  pyproject.toml              installable Python package metadata
  pipeline.py                 thin CLI shim (also runnable without install)
  src/
    rxn_core/
      __init__.py             public API re-exports
      alignment/              molecule-level WBO alignment API and branch scheduler
      growth/                 connected-fragment growth and trace events
      matcher/                symmetry-compressed candidate matching
      frag.py                 xtb single-point runner, WBO graph,
                              classify_bonds (broken/formed detector)
      modes.py                per-mode features (beta, rho, kappa) +
                              g98.out parser + mode reindex helpers
      align.py                load_cached_xtb, reindex_to_R_frame
      pipeline.py             BGCP full-view pipeline + main()
      plain_pipeline.py       older source-xyz/xtb execution pipeline
  out/                        pipeline output (gitignored)
```

## Install

```bash
# Python 3.10+, xtb on $PATH (e.g. conda install -c conda-forge xtb)
pip install -e .
```

After install:

- `import rxn_core` works from anywhere — re-exports the alignment,
  Hessian, and feature primitives at the top level (see "Public API"
  below).
- `rxn-core-pipeline ...` is on `$PATH` as the BGCP full-view console script.

The default BGCP pipeline reads existing xtb caches.  The older source-xyz
pipeline that runs xtb is preserved as `rxn_core.plain_pipeline`.
https://github.com/grimme-lab/xtb

## Run the pipeline

```bash
# Installed, one or more cached BGCP steps:
rxn-core-pipeline --steps pr7.V.dodh_ts910 --inner-workers 10

# Or, without install, from a clone:
python pipeline.py --steps pr7.V.dodh_ts910 --inner-workers 10

# Backward-compatible wrapper:
python build_bgcp_views_v2.py --steps pr7.V.dodh_ts910 --inner-workers 10
```

By default the pipeline reads:

```
appendix_perparation/xtb_frequency_calculations/<step>/
  R/                    reactant xyz + wbo
  P/                    product xyz + wbo
  sp_groundtruth/       GT TS xyz + wbo
  hess_groundtruth/     GT g98.out
  sp_iter<N>/           IG xyz + wbo
  hess_iter<N>/         IG g98.out
```

Override paths with `RXN_CORE_PROJECT`, `BGCP_WORK`, `BGCP_OUT_ROOT`, and
`BGCP_EVAL_JSON` if needed.  Output lands in `out/bgcp_views/<step>/`.

## Output viewer

`view.html` is self-contained (loads 3Dmol from CDN; no other deps).
Open in any modern browser. Layout:

- **Top row**: R (left) and P (right), static, with broken-bond
  markers (red dashed) on R and formed-bond markers (green dashed) on
  P. P is reindexed to R atom order so bond pairs draw on the right
  atoms.
- **Grid**: every IG panel, sorted by S descending. Each panel
  auto-animates its picked imaginary mode (cylinder markers wiggle
  with the molecule). IGs with `n_imag = 0` render static with a
  badge.
- Each panel labels `S, β, ρ, κ, n_imag, picked_freq`.

## Configurable parameters

All defaults are sane for the targeted system class (small organics
through ~150-atom organometallics). Override only when you know what
you're doing.

### Pipeline (`src/rxn_core/pipeline.py`)

| name | default | meaning |
|---|---|---|
| `--workers` | CPU count - 1 | total CPU budget in auto mode; outer step workers in outer mode |
| `--parallel-mode` | `auto` | `auto` balances concurrent steps and inner R-P/TS workers; `outer` is legacy serial-inside-step mode; `inner` runs one step at a time with parallel inner work |
| `--inner-workers` | `0` | explicit per-step inner workers; `>1` selects inner mode unless `--parallel-mode` is set |
| `--auto-inner-workers` | `8` | target inner workers per concurrent step in auto mode |
| `--steps` | all cached steps | explicit step names |
| `--limit` | none | first N cached steps |
| `BGCP_CUT_FLOOR` | `0.2` | R-P sweep cuts every R edge at or above this WBO |
| `BGCP_ISO_TOL` | `1.0` | WBO tolerance used by BGCP matching |
| `BGCP_VIEW_MAX_BRANCHES` | `100` | per-cut branch cap; a cut is discarded if any seed order reaches it |
| `BGCP_PARALLEL_MODE` | `auto` | env default for `--parallel-mode` |
| `BGCP_AUTO_INNER_WORKERS` | `8` | env default for `--auto-inner-workers` |
| `BGCP_TIMING` | `0` | set to `1` for per-target timing prints |
| `BGCP_TS_CORE_MAX_CANDIDATES` | `20000` | cap for mechanism-local TS/IG core mappings |

In `auto` or `inner` mode, the per-step inner worker budget is used for both
expensive phases: R-P cut sweep and TS/IG endpoint core matching.  The TS/IG
stage dispatches independent `(target, mechanism, endpoint R/P)` tasks, then
merges the R-derived and P-derived core pools before scoring.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the alignment module
boundaries.

### Alignment (`src/rxn_core/alignment/`)

The low-level pair aligner is `align_from_arrays(...)`; key kwargs:

| name | default | meaning |
|---|---|---|
| `graph_floor` | `0.2` | min WBO to admit an edge into the alignment graph (lower = more weak/partial bonds survive) |
| `iso_tol` | `1.0` | per-edge WBO match tolerance during subgraph iso (looser = more permissive matching at the reactive site, where bonds have changed) |
| `n_seeds` | `3` | number of seed orderings explored for one match |
| `max_branches` | `1_000_000` | low-level branch cap for direct matching; R-P cut sweep defaults to `100` and discards any cut whose seed order reaches the cap |
| `min_lock_size` | `1` | minimum fragment size that can be locked during island growth |
| `chirality` | `True` | score chirality violations as a tiebreaker between equally-mapped branches |
| `dwbo_threshold` | `0.5` | (in `classify_bonds`) min |ΔWBO| to classify an edge as broken or formed; also gates "is wR even a real bond" since wP ≥ 0 |
| `return_all` | `False` | when `True`, returns every distinct mapping in `out['all_scored']` (used internally by the per-IG branch sweep) |

`align_from_arrays` is still available as a pure low-level function.  The
BGCP full-view pipeline uses `cut_sweep(...)` for R-P mechanism discovery and
`ts_core_pool(...)` for mechanism-local GT/IG core matching.

R-P mechanism discovery uses the core `cut_sweep(...)` API in
`src/rxn_core/alignment/sweep.py`: no-cut plus one-edge R cuts above
`cut_floor`, deduped by symmetry-canonical broken/formed bond signatures.
Mechanism-local TS/IG core matching lives in
`src/rxn_core/alignment/ts_core.py`.

### Score formula

```
S = β · (1 + W_RXN · ρ) · (1 + W_CORE · κ) / n_imag^IMAG_PEN
```

The viewer applies **no filter** on `n_imag` or `ρ`. Every IG with at
least one imaginary mode is scored and animated; IGs with `n_imag = 0`
render as static structures and sink to the bottom by S = 0.

## Score formula reference

```
β_k  = |d_k · V̂|                / ‖d_k‖    ∈ [0, 1]
ρ_k  = |d_k · Δ̂_core|           / ‖d_k‖    ∈ [0, 1]
κ_k  = Σ_{i ∈ core} ‖d_k[i]‖² / Σ_i ‖d_k[i]‖²  ∈ [0, 1]

picked = argmax_{k : ω_k < 0} β_k

S = β_picked · (1 + w_r · ρ_picked) · (1 + w_c · κ_picked) / n_imag^p
```

V is built from the IG geometry: for each broken bond (i, j), V[i] -= u_ij
and V[j] += u_ij (atoms moving apart); for each formed bond (i, j), V[i]
+= u_ij and V[j] -= u_ij (atoms moving together). Δ is the reaction-coord
displacement P − R after Kabsch alignment, restricted to core.

## Algorithm details

`ALGORITHM.md` covers the priority-queue alignment in depth (seed
ordering, fragment growth, branching on set-non-uniqueness, chirality
scoring).
