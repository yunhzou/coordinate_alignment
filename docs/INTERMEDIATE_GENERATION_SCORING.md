# Intermediate Generation Scoring

This note records the current working model for proposing the next intermediate
from a single xTB-analyzed structure. The goal is to generate plausible moves
without using known products, child structures, or reaction labels.

## Inputs

For one intermediate, keep only descriptors computed from that structure:

- xTB WBO matrix
- xTB partial charges
- xTB condensed Fukui indices
- Cartesian geometry and distances

These are not strict ab initio observables. They are semiempirical descriptors,
but the proposal logic is product-blind.

## Move Channels

Do not force every possible next step into one global score too early. Treat
intermediate generation as competition between elementary move types.

### 1. Formation Likelihood

This is the original donor-acceptor matrix.

It asks:

```text
Which absent or very weak contact wants to form?
```

Useful terms:

- donor strength from `f-`, charge, and atom identity only through computed
  descriptors
- acceptor strength from `f+`, charge, and computed descriptors
- distance factor
- penalty if the pair already has substantial WBO

This channel works well when the next step is mainly bond formation, as in the
parent `2` examples.

### 2. Cleavage Likelihood

This asks:

```text
Which existing bond/contact wants to weaken or break?
```

Raw WBO should not be compared globally across all bond classes. A metal-ligand
contact and an organic C-O bond can both have WBO near `0.6-0.8`, but that does
not mean they have the same chemical lability.

Therefore cleavage should be ranked within comparable classes:

- organic framework bonds: C-C, C-O, O-C, etc.
- metal-substrate contacts: C-Au, O-Au, etc.
- metal-ligand contacts: Au-P, Au-L, etc.

Useful terms:

- unusually low WBO for that class
- bond stretch relative to a size scale
- charge polarization across the bond
- endpoint Fukui response and local softness
- optional constrained stretch scan for serious candidates

For intermediate `3`, the exploratory static scorer found that the top global
cleavage contacts were Au coordination contacts, while the first organic
framework cleavage was `C7-O26`.

### 3. Coupled Move Likelihood

Most useful mechanistic proposals may be neither pure formation nor pure
cleavage. They can be coupled moves:

```text
breaking A-B makes forming C-D easier
forming C-D pays for weakening A-B
```

This channel should combine the top formation and cleavage signals when they
are spatially and electronically connected.

For intermediate `3`, the product-blind signal was:

```text
cleavage:  C7-O26
formation: C2 -> C27/O26/O28 cluster
proposal:  C7-O26 weakens while C2 engages the acyl/carboxyl fragment
```

## Proposed Workflow

1. Compute one-body donor and acceptor rankings.
2. Build the formation matrix for nonbonded or weakly bonded pairs.
3. Build cleavage rankings separately by bond class.
4. Build coupled proposals by pairing nearby high-scoring formation and
   cleavage events.
5. Generate candidate geometries from:
   - top pure formation moves
   - top pure cleavage moves
   - top coupled formation-cleavage moves
6. Relax or sanity-check generated candidates with xTB.

The final generator should expose the separate channel scores first, then use a
combined priority score only for selecting which candidate geometries to try.
