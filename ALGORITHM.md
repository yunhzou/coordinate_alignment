# R<->P Atom Alignment Algorithm

This document states the WBO-weighted atom alignment algorithm for
reactant/product and reactant/TS/IG alignment.  The implementation is
symmetry-aware during growth: it does not enumerate a concrete one-to-one
bijection for every symmetric atom permutation.  Instead, it carries
hierarchical symmetry blocks internally and materializes one deterministic
witness mapping only at API boundaries.

The implementation is split by abstraction; see
`docs/ARCHITECTURE.md`. The main code paths are `src/rxn_core/alignment/`,
`src/rxn_core/growth/`, and `src/rxn_core/matcher/`.

## Inputs

- Element lists, coordinates, and Wiberg bond-order matrices for `R` and `P`
  (or `T` / `IG`).
- Identical composition: `Counter(elR) == Counter(elP)`.
- The alignment object is a thresholded WBO graph plus the full WBO matrix.
  Graph edges, by default `WBO >= graph_floor = 0.2`, define which R-side
  pairs participate in local island matching.  The full matrix is retained for
  exact WBO values, scoring, traces, and bond-change classification.

## Output

- A symmetry-aware witness mapping `R atom index -> P atom index`. This is not
  a promise to materialize every atom in a total concrete bijection.
- Broken and formed bonds from `classify_bonds`.
- For BGCP views, mode scores for GT and IGs under each minimal mechanism.

## Core Principles

1. **Weighted active-edge graph first.** Candidate growth matches element
   labels and WBO values on active R-side graph pairs.  Coordinates are used
   later for scoring, chirality, and visualization.  Zero/non-frontier R pairs
   are not local iso constraints during island growth; final chemistry is
   still scored from the full WBO matrices.

2. **Validity is active-pair WBO consistency.** A partial mapping `m` is valid
   on its grown fragment when it is injective, element-preserving, and every
   active R-side pair in the grown fragment maps to an active target-side pair
   and satisfies:

   ```
   WBO_P[m[i], m[j]] >= graph_floor
   abs(WBO_R[i, j] - WBO_P[m[i], m[j]]) <= iso_tol
   for every R graph edge (i, j) in the current island
   ```

   Non-edges on the R side are not checked during local extension.  This avoids
   rejecting a valid local match because the target has an additional bond to a
   fragment atom that is not connected to the new R atom; such differences are
   handled later as mechanism-level formed/broken bonds.

3. **No concrete symmetry explosion.** Symmetric choices are represented as
   local `_SymBlock(r_atoms, p_atoms)` pools inside `_SymCand`.  A block says:
   these R atoms occupy this P atom pool up to symmetry.  The object keeps a
   deterministic witness mapping, but branch identity is the block structure,
   not every possible permutation.  Extension is allowed to reshuffle the
   witness inside a block when another correlated assignment satisfies the
   same compressed state.  If a child atom resolves an active parent block to
   one witness while all alternatives remain symmetry-equivalent, `_SymCand`
   carries their represented count as `multiplicity`; when concrete alternate
   witnesses are available, they can be re-used if a later frontier atom cannot
   extend the primary witness.

4. **Hierarchical symmetry centers.** `_nauty_orbits` runs exact pynauty
   automorphism orbit detection on both graphs after encoding atom elements
   and 0.2-tolerance WBO buckets as vertex colors.  The resulting orbit IDs
   are the hierarchy used to group seed targets, extension targets, and
   chemistry signatures.

5. **Lock only after saturation.** `_set_unique(cands)` is false if any
   candidate has an open symmetry block.  Even one fully resolved candidate is
   not allowed to lock while the growth heap still has live edges, because
   pending edges may still extend the fragment or become deferred one-hop
   boundary evidence.  `_set_unique` is only a finalization check after heap
   exhaustion.

6. **Growth outcomes are 0/1/many, not arbitrary choices.** A unique candidate
   can grow to many valid targets when it sits at a symmetry center.  That is a
   first-class state, not an ambiguity to resolve by picking one target.

7. **The popped growth edge is not special.** Extension validity uses the same
   active R-pair weighted-vector rule for every already-fragmented atom `r`
   where `(n, r)` is an R graph edge:
   `WBO_P[v, map[r]] >= graph_floor` and
   `abs(WBO_R[n, r] - WBO_P[v, map[r]]) <= iso_tol`.  The heap only chooses
   traversal order.  It must not add a chemistry-specific anchor rule.

8. **Dedupe must preserve observed future distinguishability.** Two candidates
   are true duplicates only if their internal orbit state and their deferred
   one-hop boundary state are symmetry-equivalent.  A side of a symmetric
   island that already failed to absorb a boundary atom is not duplicate with
   the other side.

9. **Public API returns witnesses, not total bijections.** `align_from_arrays`
   and downstream scoring receive ordinary dict witnesses, but those dicts are
   selected representatives of compressed symmetry states.  Unmapped spectators
   must not be completed by geometry.  Before scoring a finished R<->P witness,
   a bounded local symmetry repair may choose a lower bond-change realization
   inside touched product orbits.

## Symmetry-Aware Candidate Growth

`grow_island(g_R, g_P, seed, mapping, inv, ...)` grows one island from a
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
`_extend_sym_cands` replaces the old concrete fanout.  The important point is
that the deterministic witness is not trusted as the only valid assignment.
For each possible target `v`, `_support_witness_for_value` asks whether there
exists a block-internal assignment that supports the WBO vector from `n` to
the already-grown fragment.

```
for each compressed candidate:
    for each unused matching-element P target v:
        test whether some assignment inside every touched symmetry block can
        satisfy:
            WBO_P[v, m[r]] >= graph_floor
            abs(WBO_R[n, r] - WBO_P[v, m[r]]) <= iso_tol
        for every already-grown fragment atom r where R[n,r] is an active edge
```

For fixed atoms the test is direct:

```
- v must be unused by this candidate
- element_R[n] == element_P[v]
- every active R-pair from n to the grown fragment maps to an active target
  edge and has WBO delta <= iso_tol
- R-side non-edges are ignored by local extension; they are not interpreted as
  required target non-bonds
```

For atoms inside a `_SymBlock`, the test is a small constrained matching over
that block's P pool.  The support question is existential:

```
Does there exist an injective assignment inside the touched symmetry blocks
such that the WBO-vector test passes?
```

This handles correlated symmetry: if two core or shell atoms sit in the same
symmetry block, they are assigned jointly, not shuffled independently.  The
support search is capped by `SYM_SUPPORT_MAX_STATES` (default `4096`) so a
pathological block cannot create unbounded backtracking.

The popped edge anchor is included in this same existential support search. If
`u` is inside a symmetry block, the new atom must attach to one compatible
member of that block's target pool. This keeps a leaf atom, such as an H on a
symmetry-related carbon, tied to the parent symmetry state instead of matching
the H as an isolated same-element atom.

If a new atom is valid only under a particular assignment inside an existing
block, the candidate witness is refined to a correlated assignment.  When
several refined witnesses are the same weighted-symmetry state, dedupe merges
them and sums their `multiplicity`.  When concrete alternate witnesses are
retained, they are not display-only: if the primary witness has no extension,
later growth tries those alternates locally.  This is the hierarchical case:
for `Pd(CH3)4`, growing `Pd-C` represents four carbon choices; growing `C-H`
afterward represents `4 * 3 = 12` correlated assignments even though the
candidate keeps one primary witness mapping.

Targets that pass the support test are grouped before constructing children:

```
group key =
    P element
    P orbit id
    WBO-vector relation of v to the current fragment/witness
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

When a support assignment differs from the current witness but still extends an
open block, `_SymCand` updates only the witness for the touched block and keeps
the block.  When the support assignment collapses an open block into a
correlated child relation, the block may become a concrete witness plus
`multiplicity`; the represented alternatives are still part of the state.

### Growth Transition Policy

For a popped growth proposal, evaluate every live candidate and every unused
same-element target by the WBO-vector rule.  Then classify the transition by
the number of input candidates and valid output states:

- **0 -> 0.** No live candidates exist.  Growth cannot continue.

- **1 -> 0.** The only candidate has no valid extension for this proposed
  atom.  The proposal is deferred as a boundary constraint; the candidate
  itself remains live.

- **1 -> 1.** The only candidate has one valid extension.  Commit it.

- **1 -> many.** The only candidate has multiple valid targets.  This is the
  symmetry-center case.  Do not choose one target.  Compress
  symmetry-equivalent targets into a `_SymBlock`; branch only on genuinely
  distinct weighted-symmetry states.

- **many -> 0.** No live candidate can extend through this proposal.  Defer
  the boundary constraint and keep the candidate set unchanged.

- **many -> 1.** Multiple candidates extend, but the resulting states dedupe
  to one symmetry-equivalent weighted state.  Keep the compressed state.

- **many -> many.** Multiple distinct valid states remain after symmetry-aware
  dedupe.  Keep all distinct states, compressed where possible.

The heap chooses which growth proposal to try next.  It does not define
validity.  Popping edge `(u, n)` means "try adding `n` now"; all active R
pairs from `n` to the current fragment are checked, not only `(u, n)`.

In traces, `cut_all_cands` is the `1 -> 0` or `many -> 0` case for an
`extend_free` proposal: every represented candidate variant failed the complete
active-edge WBO-vector test for the popped atom.  The rejection may be caused
by any already-grown fragment atom with an active R edge to the popped atom,
not necessarily by the popped anchor edge.

### Forced Island Merge

If `n` is already in the global locked mapping, its image is forced.  If it
belongs to a prior island, the whole island is folded into the candidate with
exact images.  The code verifies all pairwise WBO constraints between the
existing fragment and the absorbed island before committing the merge.

### Commit Or Defer

If extension succeeds, the fragment grows and outgoing R edges are pushed into
the heap.  If every candidate fails, the popped proposal is not thrown away.
It is recorded as a deferred boundary constraint and the fragment is unchanged.

A deferred boundary constraint records that this island saw an outside atom
through a specific weighted relation but could not absorb it under the current
candidate family.  It is not part of the locked island's internal fragment,
but it remains part of the candidate's boundary state.

This distinction matters for symmetric islands.  If island A has two
internally symmetric sides, but side 1 already has a deferred relation toward
island B, then side 1 and side 2 are no longer interchangeable in that growth
direction.  Dedupe must see that boundary difference.

At heap exhaustion:

- no candidates or too-small fragment -> fail this seed
- resolved unique candidate -> return one concrete witness
- unresolved/non-unique saturation -> deduplicate by compressed
  deferred-boundary-aware structural signature and return one concrete witness per
  distinct compressed state

## Boundary-Aware Dedupe

Dedupe is allowed only when two candidates represent the same orbit state under
the currently observed constraints.  It is not allowed to collapse candidates
when a deferred boundary relation distinguishes one side from another.

The dedupe key has two parts.

### Internal Signature

The internal signature describes the grown fragment:

- mapped R atoms and their P images, with exact-fixed assignments preserved
- `_SymBlock` pools, including R orbit IDs, P orbit IDs, pool sizes, and
  extendability
- core mappings kept exact when the caller marks a mechanism core

This signature collapses spectator permutations inside true symmetry blocks
but keeps core assignments and deferred-boundary-distinguished assignments
when the exact pairing matters.

### One-Hop Boundary Signature

The boundary signature describes what the candidate can see just outside the
fragment:

- deferred proposals that failed to extend earlier
- WBO values from mapped R atoms to those deferred outside R atoms
- element labels and R orbit IDs of the deferred outside atoms
- the corresponding WBO possibilities from mapped P atoms or P pools to unused
  same-element P atoms
- locked neighboring island IDs, when a boundary points at an already locked
  island

Two candidates are duplicates only if both the internal signature and this
one-hop boundary signature are symmetry-equivalent.  This one-hop check catches
the important case where an internally symmetric island has two sides, but
only one side is already coupled to another island.

Boundary-aware dedupe is still compression, not enumeration.  If all boundary
vectors are symmetry-equivalent, the candidates remain compressed.  If a
boundary vector distinguishes one side, the compressed state is refined or
branched only as far as needed to preserve that distinction.

## Multi-Island Branching

`find_islands` drives `grow_island` across seed orderings.

```
precompute p_orbits and r_orbits with `_nauty_orbits(..., wbo_tol=symmetry_wbo_tol)`
branches = [empty branch]

for each seed while progress is possible:
    for each live branch:
        if seed is already mapped:
            carry branch forward
        else:
            isos = grow_island(..., p_orbits, r_orbits)
            if no isos:
                carry branch forward
            else:
                dedup isos by mechanism-state plus deferred boundary
                fork one branch per remaining iso

    dedup live branches by mechanism-state plus deferred boundary
    enforce max_branches; in R-P cut sweep, hitting the cap discards the
    whole cut rather than keeping a truncated branch set
```

This branch dedupe is a mechanism-state dedupe, not a concrete bijection
dedupe.  It uses the current symmetry-canonical broken/formed WBO-change
signature plus the deferred one-hop boundary.  During growth, the key must
preserve future distinguishability already observed through deferred boundary
constraints, without enumerating spectator permutations that have no mechanism
effect.

After complete mappings are scored, the mechanism-level signature records:

```
broken: ((R orbit pair), (P orbit pair))
formed: ((R orbit pair), (P orbit pair))
```

That final mechanism signature collapses symmetric spectator swaps without
erasing distinct bond-change patterns.  It is intentionally later than
boundary-aware growth dedupe: two partial states that might become different
mechanisms must not be collapsed just because their current internal islands
look symmetric.

## Outer Alignment

`align_from_arrays(...)` builds graphs, generates seed orders, runs
`find_islands`, materializes one justified witness from each compressed branch,
applies final symmetry repair for R<->P mappings, classifies bonds, and scores
branches by:

```
(number of broken + formed bonds,
 chirality violations)
```

The best lexicographic score wins.  `return_all=True` returns all scored
branches for view/ranking workflows.

### Final R<->P Symmetry Repair

Compressed growth may still return one legal witness from a symmetric product
orbit.  If that witness creates many bond changes, the final repair pass
searches only the product symmetry orbits touched by current broken/formed
bond endpoints:

```
affected atoms = endpoints of current broken bonds
               + R-frame endpoints of current formed bonds
touched groups = mapped atoms with same (element, product orbit)
```

Within each touched group, the repair swaps images already assigned to that
same `(element, product orbit)` group.  It never introduces a new spectator
target outside the compressed alignment.  The score is:

```
(number of broken + formed bonds,
 total absolute WBO delta on changed bonds)
```

Groups of size up to 6 may try full within-group permutations; larger groups
use improving pair swaps.  The pass is capped by `SYM_REPAIR_MAX_EVALS`
(default `20000`; `BGCP_SYMMETRY_REPAIR_MAX_EVALS` in the view builder).
This is the pr17 TS6a fix: the O/C shell can reshuffle within product orbits
so equivalent O-C pairs stay paired, while the true mechanism-level bond
breaking/forming remains.

## Sweep-Cut Mechanism Discovery

`rxn_core.alignment.cut_sweep(...)` is the core R-P mechanism discovery API.
The package pipeline in `rxn_core.pipeline` passes runtime parameters and
renders results; it does not implement the sweep algorithm. `cut_sweep`
collects mechanisms from:

- baseline graph
- plus one run per R edge removed (`WBO >= cut_floor`, default `0.2`)
- each with `n_seeds` seed orders

This is intentionally broader than a single R<->P alignment because mechanism
discovery needs multiple possible broken/formed bond patterns for later GT/IG
scoring.

Seed generation is capped: `n_seeds=3` means three seed orderings, not
one ordering per heavy atom.  The chosen anchors are heavy atoms in graph
order, then random full-order shuffles only if more trials are requested than
heavy anchors exist.

Parallel cut sweeps dispatch one work unit per cut.  Each work unit runs all
seed orders for that cut, so the cut-specific R orbit map is computed once and
reused across seeds.  Product orbits are computed once per worker process
because they are invariant across every cut and seed.  This avoids the old
`3 * (E + 1)` exact-orbit recomputation pattern for `E` cuttable R edges.

The dedupe target depends on the alignment purpose:

- R<->P mechanism discovery deduplicates by symmetry-canonical broken/formed
  bond changes.  Multiple concrete mappings with the same mechanism under R
  symmetry collapse before GT/IG scoring.
- R/P<->GT and R/P<->IG verification is mechanism-local.  For each mechanism,
  `rxn_core.alignment.ts_core_pool(...)` enumerates exact mappings for the
  same mechanism core from both endpoints: `R -> TS` using R-core WBO context,
  and `P -> TS` using the R-P witness to pull the same core into product
  indexing.  Product-derived candidates are converted back to R-core indexing,
  unioned with the reactant-derived candidates, deduped by the exact
  `R_core -> TS_core` map, scored, and the best `S` is kept.  Spectator atoms
  are never enumerated and are not filled by geometry after a core mapping is
  chosen.  Mode-score numerators use the mapped/core atoms; denominators use
  the full TS mode norm.

This matters when symmetry touches a core atom.  If an atom is a spectator,
one arbitrary representative of a symmetric group is fine.  If that same
symmetric group contains a core atom, each possible core representative can
give a different `beta/rho/kappa` because TS coordinates and normal modes live
on concrete target atoms.  The code therefore enumerates symmetry alternatives
only on the mechanism core.  A methyl H core in an 18-H symmetric environment
creates up to 18 core candidates, not 18! full spectator permutations.

Mechanism-local TS/IG core enumeration enforces preserved endpoint-core edges
against the target WBO graph (`edge_floor`, default `0.2`) and allows extra TS
partial bonds.  The R endpoint preserves R-core active edges; the P endpoint
preserves P-core active edges after the R-P mechanism witness maps the same
core into product indexing.  The optional `max_candidates` cap defaults to
`20000` in the BGCP script; hitting it is a diagnostic warning, not an
expected path for elementary steps.

This makes ranking symmetry/core based instead of full-bijection based.  R-P
cut-sweep work is never skipped by a wall-clock timeout; slow cuts must finish
or be stopped by the caller.

R<->P work units also apply the bounded final symmetry repair by default
(`BGCP_SYMMETRY_REPAIR=1`).  It can be disabled for debugging with
`BGCP_SYMMETRY_REPAIR=0`, and its local search cap is controlled by
`BGCP_SYMMETRY_REPAIR_MAX_EVALS`.

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
| `graph_floor` | `0.2` | threshold for active R/P graph edges used by frontier growth and local iso validity |
| `iso_tol` | `1.0` | WBO tolerance during candidate extension |
| `dwbo_threshold` | `0.5` | WBO delta threshold for 1-0 / 0-1 events |
| `symmetry_wbo_tol` | `0.2` | WBO tolerance for exact automorphism orbit bucketing |
| `max_branches` | `1_000_000` | live branch cap for direct low-level matching |
| `BGCP_VIEW_MAX_BRANCHES` | `100` | R-P cut-sweep branch cap; a cut is discarded if any seed order reaches it |
| `BGCP_CUT_FLOOR` | `0.2` | R-P mechanism discovery cuts every R edge with WBO at or above this floor |
| `BGCP_CUTSWEEP_CHUNKSIZE` | `1` | multiprocessing chunk size for cut-sweep work units |
| `BGCP_ISO_TOL` | `1.0` | WBO tolerance used by BGCP view cut-sweeps |
| `BGCP_DWBO_THRESHOLD` | `0.5` | WBO delta threshold for BGCP broken/formed bond classification |
| `BGCP_SYMMETRY_WBO_TOL` | `0.2` | WBO tolerance for BGCP symmetry-orbit bucketing |
| `BGCP_W_RXN` | `1.0` | reaction-coordinate overlap score weight |
| `BGCP_W_CORE` | `0.2` | core-mode fraction score weight |
| `BGCP_IMAG_PEN` | `0.3` | imaginary-mode count penalty exponent |
| `BGCP_PARALLEL_MODE` | `auto` | pipeline scheduling mode: `auto`, `outer`, or `inner` |
| `BGCP_AUTO_INNER_WORKERS` | `8` | target inner workers per concurrent step in auto mode |
| `BGCP_TIMING` | `0` | set to `1` to print per-target cut-sweep and TS endpoint timings |
| `BGCP_TS_CORE_EDGE_FLOOR` | `0.2` | minimum target WBO for preserving an R-core edge during TS/IG core matching |
| `BGCP_TS_CORE_MAX_CANDIDATES` | `20000` | cap for mechanism-local TS/IG core mappings |
| `n_seeds` | `3` | seed orders in `align_from_arrays` |
| `N_SEEDS_PER_RUN` | `3` | seed orders per cut-sweep unit in views |

The full BGCP pipeline parallelizes two independent stages inside one step.
R-P mechanism discovery is parallel over cut work units.  GT/IG scoring then
parallelizes over endpoint core-matching tasks: every `(target TS or IG,
mechanism, endpoint R/P)` pool is built independently, and the main process
merges R-derived and P-derived core maps before ranking modes.  There is no
sweep-cut step for R-TS/P-TS; those endpoint matches use the known R-P
mechanism core.

## Verification Notes

Recent single-process checks after the symmetry-block implementation:

- perfect hexamethylethane: orbit hierarchy `[18, 6, 2]`
- hexamethylethane island growth: one full 26-atom witness, max traced candidates 4
- `pr12.Co_Silylation_JACS2015_TS_B-CStep1`: one full 123-atom branch, peak
  compressed candidates 10
- `pr7.V.dodh_ts910`: two full branches preserved
- `pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion`: full BGCP view
  generated with 10 inner workers; four concrete `2/3` alignments collapse to
  one symmetry-canonical displayed mechanism
