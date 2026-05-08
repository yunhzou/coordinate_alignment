# rxn_core

Reaction-core finder + the ranking-mechanism backend for **El Agente Disco**,
the agentic transition-state discovery workflow described in the NeurIPS 2026
submission. Given a reactant and a product geometry (with possibly different
atom orderings), the alignment side identifies the atom-to-atom mapping and
which bonds broke / formed; the ranking side scores Hessian-derived
imaginary modes against that mapping to surface a top-2 of 20 LLM-generated
TS guesses for expert review.

The signal throughout is **Wiberg bond order** from xtb GFN2 single-point
calculations: WBO matrices on R, P, and every IG-TS, plus normal modes
from g98.out parsing on the IG-TSs.


## Repository layout

```
rxn_core/
  src/                          alignment + ranking sources (importable)
    rxn_core_pq.py              priority-queue subgraph-iso atom mapper
    rxn_core_frag.py            xtb runner, WBO graph, classify_bonds
    ranker.py                   rk_clean_v2: scores + diversity penalty
    bgcp_io.py                  benchmark-path constants, xyz IO
    trace_html.py               slider-driven 3Dmol HTML template
    align_bgcp_coords.py        align all 20 IGs onto GT (Kabsch)
  viewer/                       drivers that emit HTML viewers / traces
    build_pq_viewer.py          combined static viewer for all steps
    build_pq_regression_traces.py    10 random-seed traces per step (regressions)
    build_large_molecule_traces.py   PQ traces for the largest steps
    build_mode_viewer.py        per-step three-panel mode viewer
    build_pq_regression_viewers.py   per-regression static R/P viewer
  analytic_scripts/             one-shot analyses + paper figures
    figure_pass1_pass2_with_reactot.py
    figure_ig_diversity_rmsd.py
    figure_passk_comparison.py
    ablation_score_simplification.py
    verify_pass2_lift_diversity.py
    compare_alignment_variants.py
    eval_dwbo_overlap.py
    improve_ranker.py
  appendix_perparation/         curated geometries + viewer payloads
    Pure_Geometries_Elementary_Step/
      Benchmark_Guesses_Collective_Package/   155 step folders (xyz)
      Benchmark_Guesses_Coordinate_Aligned_Version/   IGs Kabsch-aligned to GT
    viewer/
      mode_viewer/<step>.html   per-step three-panel mode viewer payloads
      flat_view.html            single-page viewer used in pass@2 evaluation
    analtics/
      final_quality_measurement-humanversion (1).csv   human pass@1/pass@2 labels
  out/                          generated viewers / traces / CSVs (gitignored)
  ts_discovery_paper/           the NeurIPS submission (separate git repo)
  run_pq_bgcp.py                CLI: alignment over the whole benchmark
  build_appendix_final.py       build the anonymized supplementary release
```

The dataset itself (`Benchmark_Guesses_Collective_Package/`) is gitignored; a
symlink is created when first needed (see Setup).


## Setup

Requirements:

- Python 3.10+
- [xtb](https://github.com/grimme-lab/xtb) on `$PATH` (GFN2-xTB; provides WBO
  matrices and Hessians)
- pip packages: `numpy`, `scipy`, `networkx`, `matplotlib`, `pypdf` (pdf
  inspection only)

```bash
# Symlink the benchmark dataset to where bgcp_io expects it
ln -s appendix_perparation/Pure_Geometries_Elementary_Step/Benchmark_Guesses_Collective_Package \
      Benchmark_Guesses_Collective_Package
```

xtb working directories live under `work_bgcp/` and are gitignored. They are
caches: deleting them just forces recomputation.


## Algorithm sketch

Atom mapping (`src/rxn_core_pq.py`):

```
xtb on R, P                     ->  WBO matrices
build_graph (graph_floor 0.2)   ->  g_R, g_P

for each random seed ordering:
    for seed in order:
        grow_island_pq:
            heap of (-wbo, frag_atom, ext_atom)
            each pop tries to extend fragment to ext_atom
            cands = all valid subgraph isos of fragment in g_P
                     (incremental; recomputed at every step => order-independent)
            if n is in mapping (island), force whole-island merge
            consume edge if no cand survives
        lock fragment when heap empties or cands == 1
        branching when set-non-unique
expand_mapping for symmetric H / methyls
classify_bonds: |dWBO| >= 0.5 with WBO >= 0.5 floor
score by (broken+formed, chir_violations, -mapped); pick best branch
```

Ranking (`src/ranker.py`, `rk_clean_v2`):

```
S(g) = beta * (1 + w_r * rho) * (1 + w_c * kappa) / n_imag^p
       beta  = bond-direction overlap of mode m with broken+formed bonds
       rho   = reactive-overlap (mode mass on reactive-core atoms)
       kappa = core-fraction (fraction of mode energy in the core)
filter: 1 <= n_imag <= 2  AND  rho >= 0.10
top-2 selection: greedy, with mass-weighted-cosine diversity penalty
                 alpha = 0.7
defaults: w_r=1.0, w_c=0.2, p=0.3
```

Knobs and defaults live in `src/rxn_core_pq.py` (alignment) and
`src/ranker.py` (ranking).


## How to run things

### 1. Run alignment over the whole benchmark

```bash
python run_pq_bgcp.py
```

Writes a per-step CSV under `out/`. Uses 10 random seed orderings per step
and records the best branch.

### 2. Generate alignment traces for one or many steps

The trace is a slider-driven HTML page (3Dmol on R and P side-by-side, plus
an event log) that lets you step through every PQ decision: seed pop, edge
consumption, fragment lock, branching.

For the largest molecules (six steps, 105--149 atoms; Co/Pd/Ni/carbene
chemistry):

```bash
python viewer/build_large_molecule_traces.py
# outputs: out/large_alignment_traces/<step>/pq_seed_<i>.html
#          out/large_alignment_traces/index.html  (top-level)
```

For the regression-set (the seven steps where the new PQ algorithm changed
behavior vs. the old greedy mapper):

```bash
python viewer/build_pq_regression_traces.py
# outputs: out/regressions/<step>/pq_seed_<i>.html
```

To run on your own choice of steps, copy `viewer/build_large_molecule_traces.py`
and edit the `LARGE_STEPS` list.

### 3. Generate the combined static viewer

```bash
python viewer/build_pq_viewer.py
# outputs: out/bgcp_pq_viewer.html
```

A single HTML page with a dropdown over all 155 steps, R/P 3D side-by-side,
broken (red) / formed (green) bond cylinders, and regression markers vs.
the old baseline.

### 4a. Build the one-step ranked view (R + P + 20 IGs, sorted by score)

End-to-end demo for a single step: alignment + core-atom identification +
xtb Hessian + picked imaginary mode + ranker score + a single HTML page
with R, P, GT, and all 20 IG panels (animated on their picked mode,
sorted by descending ranker score; IGs that fail the
`n_imag <= 2 / rho >= 0.10` filter are shown as static structures).

```bash
python viewer/build_ranked_view_one_step.py [STEP_NAME]
# default STEP_NAME = pr16.carbocation_ts11  (18 atoms, 1/1 broken/formed)
# outputs: out/ranked_views/<step>.html
```

Reuses the per-step mode_viewer payload at
`appendix_perparation/viewer/mode_viewer/<step>.html` (so step 4 below must
have been run for `<step>` first).

### 4. Build the per-step three-panel mode viewer

```bash
python viewer/build_mode_viewer.py
# outputs: appendix_perparation/viewer/mode_viewer/<step>.html
#          appendix_perparation/viewer/flat_view.html  (single-page browsing)
```

This is the viewer used by an expert chemist to score pass@1 / pass@2.

### 5. Reproduce the paper figures

```bash
# Headline pass@1 / pass@2 + ReactOT bar chart by atom-count bin
python analytic_scripts/figure_pass1_pass2_with_reactot.py

# IG-vs-GT geometric diversity (Kabsch RMSD): all 20 IGs vs. ranker's top-2
python analytic_scripts/figure_ig_diversity_rmsd.py

# pass@k comparison
python analytic_scripts/figure_passk_comparison.py
# outputs: appendix_perparation/figures/*.{png,pdf}
```

### 6. Verify the appendix claim that diversity drives the pass@2 lift

```bash
python analytic_scripts/verify_pass2_lift_diversity.py
```

Splits the 155 steps into BOTH_GOOD / P1_ONLY / LIFT / BOTH_BAD by their
human (IG#1, IG#2) labels and reports the median mass-weighted cosine and
Kabsch RMSD between the verifier's top-1 and top-2 picks per group, plus a
Mann--Whitney p-value. (LIFT MWC median 0.160 vs BOTH_GOOD 0.870, p=0.001.)

### 7. Ablation: score = bond_overlap only

```bash
python analytic_scripts/ablation_score_simplification.py
```

Reports how often the simplified score `S(g) = beta` would change the top-2
selection vs. the full multiplicative score.

### 8. Build the supplementary-release package

```bash
python build_appendix_final.py
```

Produces `ts_discovery_paper/appendix_final/`:

- `benchmark.zip` -- 155 step folders, anonymized (Jackie_TS_<n> -> TS_<n>;
  /lustre/.../yunhengz paths stripped from xyz comment lines).
- `flat_view.html` -- the per-step three-panel mode viewer with author
  identifiers stripped.
- `README.md` -- contents and intended use.


## Trace HTML output format

Every trace produced by `viewer/build_*_traces.py` is one self-contained
HTML page (no server, no build step). Open in a browser. The page has:

- Two 3Dmol viewers (reactant left, product right).
- A slider plus prev / play / next buttons.
- An event log (monospace) showing the current step's decision.
- Color legend: yellow = seed atom, green = fragment atoms, red =
  consumed/cut edges, purple = candidate mappings, dark green = locked
  island.

Each event in the underlying JSON is one of:
`seed_start`, `commit`, `seed_end`, `island_locked`, `pass_start`,
`consumed`, `done`. The `consumed` event records `(frag_atom, ext_atom,
wbo, reason)` for an edge the algorithm popped from the heap but could
not extend; this is the most informative event when debugging a failed
mapping.


## NeurIPS submission

The paper itself is in `ts_discovery_paper/` (separate git repo, no remote).
Build with `tectonic`:

```bash
cd ts_discovery_paper
tectonic -X compile neurips_2026.tex
# outputs: neurips_2026.pdf  (29 pages: 9-page body + appendix)
```

The supplementary release lives at `ts_discovery_paper/appendix_final/`.
