# Single-Step Retrosynthesis by Precursor Fragment Assembly

## 1. Objective

Given:

- a target product graph `P_target`; and
- a catalog of purchasable precursor graphs `R_i`;

find small combinations of precursors whose retained fragments can cover
`P_target` in one synthetic step.

The basic search unit is a coherent precursor fragment, not an individual
atom. Atoms of a precursor that do not belong to the target are retained as
explicit leftover fragments instead of being forced to match target atoms.

The first stage establishes structural feasibility. Chemical feasibility and
reaction ranking are later stages.

## 2. Candidate Record

One precursor may produce several candidate records. Each record contains:

```text
FragmentCandidate
|- precursor_id
|- retained_atoms_R
|- retained_fragments_R
|- placement_R_to_P
|- covered_atoms_P
|- leftover_fragments_R
|- broken_boundary_bonds_R
|- augmented_fallback_placements
|- island_size
|- branch_count
|- complete
`- status: matched | no_match | capped | rejected_fragmentation
```

`complete=false` means that a branch cap was reached and the returned
candidates are not an exhaustive result.

## 3. Stage A: Discover a Retained Fragment

For each precursor `R`:

1. Grow a connected fragment `H` from `R` into `P_target` using weighted
   subgraph isomorphism.
2. Stop a growth path when its connected island saturates. Do not force later
   unmatched precursor atoms to become additional singleton islands.
3. Retain the largest clean placements and their exact atom mappings.
4. Apply exact symmetry deduplication before counting distinct branches.
5. Stop when the configured branch cap is reached and return `status=capped`.

A match is defined by both the retained precursor atoms and their placement in
the target. If the same retained fragment can occupy several distinct target
regions, those placements remain competing candidates unless an exact
automorphism proves them equivalent.

The initial implementation may keep only maximum-size placements. A later
implementation may retain bounded near-maximum placements when needed for
complete assembly.

## 4. Stage B: Separate the Leftover Fragments

For one retained fragment `H`:

1. Mark every bond of `R` with one endpoint in `H` and the other outside `H`.
   These are the broken boundary bonds.
2. Remove the broken boundary bonds from a working copy of `R`.
3. The connected components containing no atom of `H` are the leftover
   fragments `F_1, F_2, ...`.
4. Record the number of boundary bonds, leftover components, and leftover
   heavy atoms.

This cleanup makes each residual component a coherent competitive unit. It
may subsequently map into unused target atoms, but it cannot be forced across
the cut boundary or fragmented into unrelated singleton matches.

Candidates with excessive boundary cuts or highly fragmented leftovers may
be rejected or deprioritized. The thresholds are configuration values and are
not hard-coded into the graph model.

## 5. Stage C: Construct the Augmented Product

Construct a temporary augmented product graph:

```text
P_augmented = P_target . F_1 . F_2 . ...
```

The dot denotes disconnected graph union.

Each leftover atom in `R` has an identical copied atom in `P_augmented`.
These copies are fallback positions, not anchors. Residual components remain
free to compete for unused positions inside `P_target`; a component becomes a
spectator only when its final placement uses an appended copy.

During discovery, the retained atoms of `H` are free to compete for positions
inside `P_target`. During candidate-specific validation, their discovered
placement is fixed temporarily so that validation checks that exact candidate
without repeating the placement search.

## 6. Stage D: Validate the Complete Precursor Placement

Validate the candidate using:

```text
query  = R with its recorded boundary bonds removed
target = P_augmented
anchors = initially retained R atoms -> their discovered P_target positions
```

Run weighted subgraph matching with only the candidate's initially retained
placement fixed. Residual components can choose unused `P_target` positions
or their copied fallback positions.

A valid result must satisfy all of the following:

- every atom of the cut precursor query is mapped;
- every residual atom maps either into unused `P_target` or an appended copy;
- atoms mapped into `P_target` become additional retained fragments;
- atoms mapped into appended copies become explicit spectators;
- the retained placement agrees with the discovered fragment placement; and
- no target atom receives more than one precursor atom.

The augmented product may contain unused `P_target` atoms. Those atoms are
available for fragments supplied by other precursors.

Example:

```text
precursor R:       Ar-Br
retained H:        Ar
leftover F:        Br
boundary cut:      Ar-Br
augmented product: P_target . Br

locked placement:  precursor Ar -> its discovered region in P_target
competitive Br:    unused P_target position or appended Br fallback
```

## 7. Stage E: Emit Fragment-Unit Candidates

For every validated placement, emit one `FragmentCandidate`.

The important assembly identity is:

```text
(precursor, retained fragments, target placement, boundary cuts)
```

Coverage alone is insufficient because two candidates covering the same
target atoms may require different precursor cuts or leave different
fragments.

## 8. Stage F: Assemble Complete Product Coverage

Combine `FragmentCandidate` records under the following rules:

1. A selected candidate contributes all retained fragments from one precursor
   as one source unit.
2. Selected candidates must not cover the same target atom.
3. The union of their covered target atoms must cover every atom of
   `P_target`.
4. Product bonds whose endpoints occur in the same retained fragment should
   be supported; bonds crossing recorded fragment cuts are proposed edits.
5. Product bonds whose endpoints come from different candidates are proposed
   newly formed bonds.
6. Recorded precursor boundary bonds are proposed broken bonds.
7. Every selected precursor must contribute at least one retained fragment.
8. The number of material-contributing precursors is bounded.

The assembly search operates on fragment candidates rather than raw catalog
pairs. It may use bounded beam search or weighted exact cover.

## 9. Initial Assembly Ranking

Completed assemblies are ranked using structural terms:

1. complete target coverage;
2. fewer selected precursors;
3. fewer formed product bonds;
4. fewer broken precursor boundary bonds;
5. larger retained precursor fragments;
6. fewer and cleaner leftover fragments; and
7. fewer capped or incomplete precursor searches.

Chemical reaction rules, valence repair, charge balance, stereochemistry,
reagent requirements, and learned reaction likelihood are deliberately
separate later validators.

## 10. Branch-Cap Semantics

A branch cap is a computational limit, not evidence that no match exists.

Every capped search must report:

```text
status = capped
complete = false
cap_stage
branch_limit
maximum_branch_count
best_fragment_size_found
```

The implementation must never silently present truncated candidates as an
exhaustive match set. Assemblies containing capped results may be inspected,
but they must remain distinguishable from assemblies built entirely from
complete searches.

## 11. First Experimental Milestone

The first milestone stops before multi-precursor assembly and reports, for a
given `P_target`:

- the largest validated retained fragments;
- their placements in `P_target`;
- their source precursor identifiers;
- their explicit leftover fragments;
- their broken boundary bonds;
- target coverage by coherent fragment units; and
- branch-cap and fragmentation diagnostics.

This result establishes whether the catalog contains useful structural units
before implementing complete product-cover assembly.
