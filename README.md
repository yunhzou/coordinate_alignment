# rxn_core

End-to-end pipeline that takes a reactant + product geometry plus a set
of LLM-generated transition-state guesses (IGs) and produces a single
HTML page where the IGs are sorted by a chemistry-aware score and
animated on their picked imaginary mode.

The signal throughout is **Wiberg bond order** from xtb GFN2 single-point
calculations: WBO matrices on R, P, and every IG-TS, plus normal modes
parsed from xtb Hessians.

## Pipeline overview

```
Inputs:                                        Outputs (per step):
                                               out/ranked_views/<step>/
  source/reactant_*.xyz       ────┐              view.html        (open in any browser)
  source/product_*.xyz             ├──► pipeline ──►   alignment.json
  initial_ts_guesses/*.xyz        ┘              scores.csv
  generation_report.json                         aligned/<...>.xyz
    (charge / multiplicity)                      modes/<label>_picked.xyz
                                                 xtb/{R,P,hess_iter<N>}/
                                                 README.md
```

The pipeline:

1. **R/P single-points** — xtb GFN2 on R and P → WBO matrices.
2. **R↔P alignment** — priority-queue subgraph-iso atom mapper
   (`src/rxn_core_pq.py`). Returns mapping + broken/formed bond
   classification.
3. **Reactive core** — atoms touching any broken or formed bond.
4. **xtb hess on each IG** (parallel) — produces normal modes (g98.out)
   and a WBO matrix.
5. **For each IG**: align IG↔R (every alignment branch), reindex modes
   to R-frame, compute three per-mode features:
   - β (bond_overlap): mode displacement projected on the bond-axis
     stretch/contract direction at TS coordinates
   - ρ (rxn_overlap): mode core-atom motion projected on the R→P
     direction
   - κ (core_fraction): fraction of mode energy on core atoms
6. **Pick the imaginary mode** with max β.
7. **Score**:
   ```
   S = β · (1 + w_r · ρ) · (1 + w_c · κ) / n_imag^p
       (w_r = 1.0,  w_c = 0.2,  p = 0.3)
   ```
8. **Choose the alignment branch** that gives the highest S. Branches
   that differ only on non-core (spectator) atoms are score-equivalent
   and collapsed.
9. **Sort IGs by S descending**, render one HTML page with R, P, and
   all 20 IG panels each animated on its picked mode.

## Layout

```
rxn_core/
  pyproject.toml              installable Python package metadata
  pipeline.py                 thin CLI shim (also runnable without install)
  src/
    rxn_core/
      __init__.py             public API re-exports
      pq.py                   priority-queue subgraph-iso atom mapper
      frag.py                 xtb single-point runner, WBO graph,
                              classify_bonds (broken/formed detector)
      modes.py                per-mode features (beta, rho, kappa) +
                              g98.out parser + mode reindex helpers
      align.py                load_cached_xtb, fill_unmapped_greedy,
                              reindex_to_R_frame
      pipeline.py             end-to-end pipeline + main()
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
- `rxn-core-pipeline <step_dir>...` is on `$PATH` as a console script.

xtb (GFN2) handles all electronic-structure work — single-points for R,
P, and every IG, plus Hessians on the IGs.
https://github.com/grimme-lab/xtb

## Run the pipeline

```bash
# Installed:
rxn-core-pipeline <step_dir> [step_dir ...]

# Or, without install, from a clone:
python pipeline.py <step_dir> [step_dir ...]
```

Each `<step_dir>` should contain:

```
<step_dir>/
  source/reactant_*.xyz       at least one reactant xyz
  source/product_*.xyz        product xyz
  initial_ts_guesses/*.xyz    one or more IG TS xyz files
  generation_report.json      JSON with generation_spec.charge
                              and generation_spec.multiplicity
```

Output lands in `out/ranked_views/<workflow_name>/`. The xtb subtree
under that folder is the cache, so re-running the pipeline on the same
step is a sub-second rebuild of the artifacts.

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
