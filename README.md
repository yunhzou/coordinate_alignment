# rxn_core

Reaction-core finder. Given a reactant and a product geometry (with possibly
different atom orderings), identifies the atom-to-atom mapping and which
bonds broke / formed.

Signal: Wiberg bond order from xtb GFN2 single-point. The WBO matrix is
turned into a weighted graph; matching is done by **priority-queue grow
with order-independent subgraph isomorphism**.

## Layout

- `rxn_core_pq.py` — the alignment algorithm: priority-queue propagation
  from a seed atom, incremental fragment matching that recomputes all
  valid subgraph isos at every step (order-independent), branching on
  set-non-unique saturation, chirality-aware scoring.
- `rxn_core_frag.py` — low-level utilities (xtb runner, WBO graph build,
  expand_mapping, classify_bonds).
- `bgcp_io.py` — shared BGCP path constants and xyz IO helpers.
- `trace_html.py` — slider-driven 3Dmol HTML template for trace animations.

## Drivers

- `run_pq_bgcp.py` — run alignment on every BGCP step, emit a CSV.
- `build_pq_viewer.py` — combined HTML viewer for all 160 BGCP steps with
  R/P 3D, broken (red) / formed (green) bond cylinders, regression markers
  vs OLD baseline, dropdown + filter.
- `build_pq_regression_viewers.py` — per-regression static R/P viewer.
- `build_pq_regression_traces.py` — per-regression 10-seed trace animations
  (PQ algorithm, slider-driven).

## Algorithm at a glance

```
xtb on R, P                    →  WBO matrices
build_graph (graph_floor 0.2)  →  g_R, g_P
                                  
for each of n_seeds random orderings:
    for seed in order:
        grow_island_pq:
            heap of (−wbo, frag_atom, ext_atom)
            each pop tries to extend fragment to ext_atom
            cands = all valid subgraph isos of fragment in g_P
                     (incremental: extend each cand to all valid v's)
            if n is in mapping (island), force whole-island merge
            consume edge if no cand survives
        lock fragment when heap empties or cands == 1
        branching when set-non-unique
expand_mapping for symmetric H/methyls
classify_bonds: |dWBO| ≥ 0.5 with WBO ≥ 0.5 floor
score by (broken+formed, chir_violations, −mapped); pick best
```

See `ALGORITHM.md` for the full design principles.

## Knobs (defaults in `rxn_core_pq.py`)

| name | default | meaning |
|---|---|---|
| `graph_floor` | 0.2 | edge inclusion in g_R / g_P |
| `iso_tol` | 1.0 | per-edge WBO match tolerance during fragment ISO |
| `min_lock_size` | 1 | smallest fragment that can lock |
| `bond_high` | 0.5 | classify_bonds: WBO is a bond |
| `dwbo_threshold` | 0.5 | classify_bonds: change is real |
| `n_seeds` | 10 | random seed orderings |
| `max_branches` | 8 | concurrent live branches per seed order |

## Usage

```bash
python run_pq_bgcp.py                          # all 160 steps → CSV
python build_pq_viewer.py                      # combined HTML viewer
python build_pq_regression_viewers.py          # per-regression static viewer
python build_pq_regression_traces.py           # per-regression trace HTML
```

Outputs land in `out/`.
