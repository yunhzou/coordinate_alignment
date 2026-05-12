# R<->P Atom Alignment Algorithm

This is a WBO-graph atom mapper for reactant/product and reactant/TS/IG
alignment.  The current implementation is symmetry-aware during growth: it
does not enumerate a concrete one-to-one bijection for every symmetric atom
permutation.  Instead, it carries hierarchical symmetry blocks internally and
materializes one deterministic witness mapping only at API boundaries.

The implementation lives mainly in `src/rxn_core/pq.py`.

## Inputs

- Element lists, coordinates, and Wiberg bond-order matrices for `R` and `P`
  (or `T` / `IG`).
- Identical composition: `Counter(elR) == Counter(elP)`.
- WBO graphs from `build_graph(..., bond_cut=graph_floor)`.

## Output

- A mapping `R atom index -> P atom index`.
- Broken and formed bonds from `classify_bonds`.
- For BGCP views, mode scores for GT and IGs under each minimal mechanism.

## Core Principles

1. **Graph first.** Candidate growth uses element labels and WBO graph edges.
   Coordinates are used later for scoring, chirality, and visualization.

2. **No concrete symmetry explosion.** Symmetric choices are represented as
   local `_SymBlock(r_atoms, p_atoms)` pools inside `_SymCand`.  A block says:
   these R atoms occupy this P atom pool up to symmetry.  The object keeps a
   deterministic witness mapping for cheap WBO checks, but branch identity is
   the block structure, not every possible permutation.

3. **Hierarchical symmetry centers.** `_color_refine_orbits` runs iterated
   1-WL / Morgan color refinement on both graphs using element labels and
   0.2-wide WBO buckets.  The resulting orbit IDs are the hierarchy used to
   group seed targets, extension targets, and chemistry signatures.

4. **Lock only when resolved.** `_set_unique(cands)` is false if any candidate
   has an open symmetry block.  A single unresolved block is not treated as a
   unique concrete bijection.  It can lock only after enough context closes the
   ambiguity or after saturation returns a representative witness.

5. **Chemistry, not spectator labels, creates branches.** `find_islands_pq`
   deduplicates branches by joint `(R orbit pair, P orbit pair)` broken/formed
   chemistry signatures.  Spectator permutations collapse; chemistry-distinct
   mechanisms remain separate.

6. **Public API stays concrete.** `align_from_arrays` and downstream scoring
   still receive ordinary dict mappings.  Compression is internal to island
   growth and branch deduplication.

## Symmetry-Aware Candidate Growth

`grow_island_pq(g_R, g_P, seed, mapping, inv, ...)` grows one island from a
seed using a priority queue of R edges ordered by descending WBO.

### Initialization

```
seed_targets = all unused P atoms with matching element
seed_groups  = seed_targets grouped by P orbit

for each group:
    if len(group) > 1:
        make _SymCand with a non-extendable seed block:
            {seed} -> group
    else:
        make fixed _SymCand({seed: group[0]})
```

The initial seed block is non-extendable because it represents "which anchor
could this seed be" rather than "these sibling R atoms share one target pool."

### Extension

When the heap pops an edge from fragment atom `u` to outside atom `n`,
`_extend_sym_cands` replaces the old concrete fanout:

```
for each compressed candidate:
    cm = deterministic witness mapping
    bonded = R-neighbors of n already in the fragment
    v_set = intersection of P-neighbors of cm[b] for b in bonded
    filter v_set by:
        - v not globally used
        - element(v) == element(n)
        - every bonded R/P WBO matches within iso_tol
```

The valid `v` targets are grouped before constructing children:

```
group key =
    P element
    P orbit id
    relation of v to the current witness atoms
    relation of v to existing symmetry blocks
```

For each target group:

- If it lies inside an existing extendable block and `n` is R-orbit-compatible
  with that block, extend the block by adding `n`.
- Else if the group has multiple P atoms, create a new symmetry block
  `{n} -> group`.
- Else add a fixed mapping `n -> v`.

This is the main change from concrete matching: a K-way symmetric target group
creates one compressed candidate, not K dicts.

### Forced Island Merge

If `n` is already in the global locked mapping, its image is forced.  If it
belongs to a prior island, the whole island is folded into the candidate with
exact images.  The code verifies all cross-bonds between the existing fragment
and the absorbed island before committing the merge.

### Commit Or Consume

If extension succeeds, the fragment grows and outgoing R edges are pushed into
the heap.  If every candidate fails, the popped edge is consumed and the
fragment is unchanged.

At heap exhaustion:

- no candidates or too-small fragment -> fail this seed
- resolved unique candidate -> return one concrete witness
- unresolved/non-unique saturation -> deduplicate by compressed structural
  signature and return one concrete witness per distinct compressed state

## Multi-Island Branching

`find_islands_pq` drives `grow_island_pq` across seed orderings.

```
precompute p_orbits and r_orbits
branches = [empty branch]

for each seed while progress is possible:
    for each live branch:
        if seed is already mapped:
            carry branch forward
        else:
            isos = grow_island_pq(..., p_orbits, r_orbits)
            if no isos:
                carry branch forward
            else:
                dedup isos by joint orbit chemistry signature
                fork one branch per remaining iso

    dedup live branches by the same chemistry signature
    enforce max_branches
```

The joint chemistry signature records:

```
broken: ((R orbit pair), (P orbit pair))
formed: ((R orbit pair), (P orbit pair))
```

This is stricter than P-only or R-only orbit dedup: it collapses symmetric
spectator swaps without erasing product-side mechanism distinctions.

## Outer Alignment

`align_from_arrays(...)` builds graphs, generates seed orders, runs
`find_islands_pq`, expands any unmapped spectators greedily, classifies bonds,
and scores branches by:

```
(number of broken + formed bonds,
 chirality violations,
 -number of mapped atoms)
```

The best lexicographic score wins.  `return_all=True` returns all scored
branches for view/ranking workflows.

## BGCP View Cut Sweep

`build_bgcp_views_v2.py` uses `cut_sweep` to collect mechanisms:

- baseline graph
- plus one run per strong R bond removed (`WBO >= 0.5`)
- each with `N_SEEDS_PER_RUN` seed orders

This is intentionally broader than a single R<->P alignment because the view
needs multiple possible mechanisms for GT/IG scoring.  `BGCP_VIEW_MAX_BRANCHES`
caps per-alignment branch materialization in the view builder; default is
`5000`.

Parallel cut sweeps dispatch `(cut, seed_order)` work units.  Work is ordered
seed-major and defaults to `BGCP_CUTSWEEP_CHUNKSIZE=1`, so the three seed
orders for one pathological cut are not bundled onto one worker.  This improves
tail utilization on high-symmetry cases where a few cuts dominate runtime.

The BGCP view builder keeps `BGCP_ISO_TOL=1.0`, but the dedupe target depends
on the alignment purpose:

- R<->P mechanism discovery deduplicates by symmetry-canonical broken/formed
  bond changes.  Multiple concrete mappings with the same mechanism under R
  symmetry collapse before GT/IG scoring.
- R<->GT and R<->IG verification is mechanism-local.  For each displayed
  mechanism, the view builder enumerates exact mappings only for that
  mechanism's `core_R` atoms, preserves every distinct core mapping, scores
  them all, and keeps the best `S`.  Spectator atoms are never enumerated;
  they are filled greedily only after a core mapping is chosen.

This matters when symmetry touches a core atom.  If an atom is a spectator,
one arbitrary representative of a symmetric group is fine.  If that same
symmetric group contains a core atom, each possible core representative can
give a different `beta/rho/kappa` because TS coordinates and normal modes live
on concrete target atoms.  The code therefore enumerates symmetry alternatives
only on the mechanism core.  A methyl H core in an 18-H symmetric environment
creates up to 18 core candidates, not 18! full spectator permutations.

Mechanism-local TS/IG core enumeration enforces preserved R-core edges against
the target WBO graph (`BGCP_TS_CORE_EDGE_FLOOR`, default `0.2`) and allows
extra TS partial bonds.  The optional cap `BGCP_TS_CORE_MAX_CANDIDATES`
defaults to `20000`; hitting it is a diagnostic warning, not an expected path
for elementary steps.

This makes ranking symmetry/core based instead of full-bijection based.  Each
R<->P `(cut, seed_order)` work unit also has `BGCP_UNIT_TIMEOUT=10` seconds by
default as a safety guard; set it to `0` to disable.

The view still applies a final mechanism dedupe before rendering:

```
key = (
    broken R bonds canonicalized by R symmetry orbit pairs,
    formed  R bonds canonicalized by R symmetry orbit pairs,
)
```

All concrete alignments with the same key are the same displayed mechanism.
The view records collapsed source mechanism IDs/cuts in the slim JSON and
scores/renders IGs only for the deduped mechanism list.  This removes
degenerate mechanisms caused only by swapping equivalent reactant atoms while
preserving different bond-change patterns.

## Bond-Change Core Logic: 1-1, 1-0, 0-1, 0-0

This section is the chemistry core and should stay stable.

For an R atom pair `(r1, r2)` and mapped P pair `(p1, p2)`:

```
wR = wboR[r1, r2]
wP = wboP[p1, p2]
```

Using `dwbo_threshold` (default `0.5`):

- **1-1: preserved bond.** R has a bond and P has the corresponding bond.
  If `abs(wR - wP) < dwbo_threshold`, this is spectator connectivity, not a
  broken/formed event.

- **1-0: broken bond.** R has a bond and P does not have the corresponding
  bond strongly enough:

  ```
  wR - wP >= dwbo_threshold
  ```

  If one or both R endpoints are unmapped, the missing P counterpart is
  treated as `wP = 0`, so the same rule applies.

- **0-1: formed bond.** P has a bond and R does not have the corresponding
  bond strongly enough:

  ```
  wP - wR >= dwbo_threshold
  ```

  If one or both P endpoints are unmapped by the inverse mapping, the missing R
  counterpart is treated as `wR = 0`.

- **0-0: ignored.** Neither graph has a meaningful bond for the pair; it is not
  part of the mechanism.

`classify_bonds(mapping, wboR, wboP, ...)` implements these rules and returns:

```
broken, formed, core_R, core_P
```

`core_R` and `core_P` are the atoms participating in broken or formed events.
The mode scorer only needs these chemistry-relevant atoms.

## Important Knobs

| name | default | meaning |
|---|---:|---|
| `graph_floor` | `0.2` | WBO edge threshold for graph construction |
| `iso_tol` | `1.0` | WBO tolerance during candidate extension |
| `dwbo_threshold` | `0.5` | WBO delta threshold for 1-0 / 0-1 events |
| `max_branches` | `1_000_000` | live branch cap in core alignment |
| `BGCP_VIEW_MAX_BRANCHES` | `5000` | branch cap used by full BGCP view generation |
| `BGCP_CUTSWEEP_CHUNKSIZE` | `1` | multiprocessing chunk size for cut-sweep work units |
| `BGCP_ISO_TOL` | `1.0` | WBO tolerance used by BGCP view cut-sweeps |
| `BGCP_UNIT_TIMEOUT` | `10` | seconds before one cut-sweep work unit is skipped; set `0` to disable |
| `BGCP_TIMING` | `0` | set to `1` to print per-target cut-sweep timings |
| `BGCP_TS_CORE_EDGE_FLOOR` | `0.2` | minimum target WBO for preserving an R-core edge during TS/IG core matching |
| `BGCP_TS_CORE_MAX_CANDIDATES` | `20000` | cap for mechanism-local TS/IG core mappings before spectator fill |
| `n_seeds` | `10` | seed orders in `align_from_arrays` |
| `N_SEEDS_PER_RUN` | `3` | seed orders per cut-sweep unit in views |

## Verification Notes

Recent single-process checks after the symmetry-block implementation:

- perfect hexamethylethane: orbit hierarchy `[18, 6, 2]`
- hexamethylethane PQ growth: one full 26-atom witness, max traced candidates 4
- `pr12.Co_Silylation_JACS2015_TS_B-CStep1`: one full 123-atom branch, peak
  compressed candidates 10
- `pr7.V.dodh_ts910`: two full branches preserved
- `pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion`: full BGCP view
  generated with 10 inner workers; four concrete `2/3` alignments collapse to
  one symmetry-canonical displayed mechanism
