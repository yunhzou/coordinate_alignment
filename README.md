# rxn_core

Reaction-core finder. Given a reactant complex and a product complex
(possibly with different atom orderings), identifies the atom-to-atom
correspondence and the bonds that broke / formed in the reaction.

Signal: Wiberg bond order (WBO) computed from xtb GFN2 single-point.
The WBO matrix is converted into a weighted graph; matching is done by
fragment isomorphism + a soft-merge constraint propagation.

## Layout

- `rxn_core_frag.py` — main library: graph build, fragment matching,
  island finding, classification, top-level `analyze`.
- `trace_run.py` — thin renderer producing per-step HTML animation
  (`animate_<step>.html`).
- `build_combined.py` — runs `analyze` on a chosen subset and emits
  `out/viewer.html` with a dropdown for switching between steps.

The trace and the production functions share the same `find_islands` /
`expand_mapping` / `merge_touching_islands` code paths; `trace_run` just
threads `events` and `atom_island_R/P` through.

## Pipeline (`analyze`)

```
xtb single-point on R, P  ->  WBO matrices
build graphs (edge if WBO >= 0.5; node has element, edge has WBO)
find_islands           # WBO-fragment isomorphism with soft merge
expand_mapping         # element-counted symmetric pairing
classify_bonds         # broken / formed by hysteresis (>= 0.6 vs < 0.3)
remove_phantom_pairs   # post-processing band-aid (see below)
```

## Lock policy in grow_island

A fragment seeded at atom `s` accumulates candidate isos that match WBO
within tolerance (`wbo_tol = 0.5`). Locking the fragment commits one
specific iso. The current rule:

1. **Early lock** if `_set_unique(candidates)` (every cand covers the
   same P-atom set, sequences may differ). Use a **valence-aware
   tiebreaker**: pick the cand whose fragment-atom images have matching
   unmapped-neighbor element multisets.
2. **Forced lock** at growth halt (`no_frontier`, `no_strong_frontier`,
   `all_cut`): same valence-aware tiebreaker.

The valence tiebreaker scores each cand by `Σ |R-unmapped-neighbor-
element-count - P-unmapped-neighbor-element-count|` per fragment atom,
and picks the cand with minimum score.

## BFS extension fallback

When extending the fragment, atoms in the closest shell are tried first
in two tiers:

1. **Top band** (max-WBO ± `top_degen=0.1`).
2. **Lower-WBO same-shell fallback** (still ≥ `growth_min_wbo=0.6`).

If all top-band candidates cut at a chemistry boundary (e.g. a C=C
becoming C-C, |dWBO|=0.87), tier 2 lets the fragment grow through other
strong bonds in the ring before giving up. This fixed pr9 cyclopropanation
(2/4 → 0/2) by accumulating constraints through the rest of the ring
before locking.

## Known bug: pseudo-symmetric atom permutation

Pseudo-symmetric atoms (locally-equivalent C, N, etc. in different parts
of the molecule, e.g. two equivalent phosphine arms on a Co complex) can
swap during grow_island: the cands are all WBO-consistent at the moment
of lock, but they differ in which atom maps to which equivalent position.
Pick-any picks an arbitrary permutation. If that permutation is wrong
relative to the *global* chemistry, downstream classify_bonds emits
phantom (broken, formed) pairs that describe the same chemistry on
different atom labels.

The valence tiebreaker fixes this for atoms whose immediate neighbors
already disambiguate them (e.g. a C with H neighbor vs a C with N
neighbor — different valence multisets). It does **not** fix cases where
the swap propagates to atoms 2+ bonds away from the locking fragment;
those atoms aren't in the fragment yet at lock time, so their valence
doesn't enter the score.

### Post-processing: `remove_phantom_pairs`

This is **not a root-cause fix** — the underlying mapping is still
imperfect. It's a display-time band-aid that cancels (broken, formed)
pairs that share:

- the same unordered element pair (e.g. both N-C, both O-H), and
- WBO within `wbo_tol = 0.5`.

The intuition: if the algorithm flags a bond as broken on the R-side
and another bond with identical element pair and similar WBO as formed
on the P-side, those are most likely the same bond mis-labeled by atom
permutation, not two separate chemistry events. Removing them gives a
cleaner view of the *net* chemistry change.

**Caveat:** this filter cannot distinguish a true H migration (e.g.
TEMPO O-H break, new O-H form on a different O) from a phantom
permutation. It will collapse both to 0/0 if the WBOs are similar. On
the 10-step benchmark this affected one case (pr1.tempo_ts2 went 1/1 →
0/0; the H is mid-migration in a TS, so 0/0 as a "net change"
description is also defensible).

A real fix would be a multi-shell valence lookahead during BFS, or
backtracking once a phantom pair is detected post-lock. Neither is
implemented here.

## 10-step benchmark scores

Random sample with seed=42 over `Benchmark/`. `br/fm` is broken/formed
counts after `remove_phantom_pairs`.

| step | N | mapped | br/fm |
|------|---|--------|-------|
| pr1.tempo_ts2 | 29 | 29 | 0/0 |
| pr7.V.dodh_ts1015 | 38 | 38 | 1/0 |
| pr9.carbene.rearr_ts47a | 40 | 40 | 0/2 |
| pr11.cycloadditions_tsIa | 67 | 67 | 1/2 |
| pr13.Cyclobutane_JOC2023_TS-CD_step1 | 46 | 44 | 3/1 |
| pr15.Fe_crosscoupling_TS7 | 82 | 82 | 0/1 |
| pr15.Fe_crosscoupling_TS-a | 83 | 83 | 1/1 |
| pr15.Fe_crosscoupling_TS8 | 82 | 82 | 2/1 |
| pr14.Pd_hydroamination_TS14 | 133 | 133 | 1/1 |
| pr12.Co_Silylation_TS_D*-E* | 149 | 149 | 0/1 |

## Usage

```bash
python build_combined.py            # 10 random steps -> out/viewer.html
python build_combined.py 20         # 20 random steps
python build_combined.py --steps a b c   # specific steps

python trace_run.py <step_name>     # per-step animation
```

## Tunables (in `grow_island`)

- `wbo_tol = 0.5` — fragment bond WBO match tolerance.
- `growth_min_wbo = 0.6` — frontier atoms must have at least this WBO
  to a fragment atom.
- `top_degen = 0.1` — top-band width within a shell.
- `min_lock_size = 2` — singleton seeds don't lock.
