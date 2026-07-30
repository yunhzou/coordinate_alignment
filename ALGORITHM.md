# Coordinate Alignment: Current Algorithm and Design Rationale

This document describes the current reactant-to-product atom alignment
algorithm, the data structures produced by AAM, and the post-AAM selection
steps for symmetry, mechanism identity, index chirality, and RMSD. It also
records the main problems found during the TS01/TS04 and full-batch
investigation so the old witness-based design is not accidentally restored.

The central principle is:

> AAM returns analytical families of valid atom bijections. A random witness
> is provenance, not the solution. Mechanism identity, chirality, and RMSD
> selection must operate on the exact family represented by graph symmetry.

The primary implementation is in:

- `src/rxn_core/matcher/`: compressed fragment candidate matching;
- `src/rxn_core/growth/`: weighted island growth;
- `src/rxn_core/alignment/branch.py`: multi-island branch construction;
- `src/rxn_core/alignment/sweep.py`: cut sweep and mechanism grouping;
- `src/rxn_core/alignment/post_aam.py`: typed post-AAM data model;
- `src/rxn_core/alignment/index_chirality.py`: analytical families,
  chirality, and fixed-mapping RMSD selection;
- `src/rxn_core/pipeline.py`: orchestration, parallelism, serialization, and
  viewer generation.

## 1. What the Investigation Changed

The implementation originally mixed several distinct concepts:

1. one arbitrary concrete mapping encountered during branch growth;
2. noisy atom-orbit or symmetry indicators accumulated from multiple paths;
3. fragment automorphisms that actually describe allowed assignments;
4. different growth paths that happen to produce the same mechanism;
5. final geometry-based choices such as chirality and RMSD.

This produced several observable problems:

- the viewer colored atoms using aggregated symmetry indicators rather than
  the symmetry of the selected analytical branch;
- several displayed "mechanisms" were duplicate event realizations;
- individual atom orbits were incorrectly used as if they identified an edge
  orbit or an entire multi-edge mechanism;
- group-level ligand permutations disappeared when branches were deduplicated;
- chirality was checked on arbitrary witnesses, although witnesses are not a
  generating set and have no privileged geometric meaning;
- RMSD was used to rank sampled witnesses instead of the exact
  chirality-valid family;
- exact but redundant paths made post-processing much slower than AAM;
- interpolation artifacts were sometimes mistaken for mapping errors.

The current design separates these responsibilities:

```text
weighted graph matching
        |
        v
compressed fragment candidates + exact candidate automorphisms
        |
        v
completed hierarchical AAM branches
        |
        v
exact mechanism event certificate
        |
        v
analytical mapping cosets, deduplicated by containment
        |
        v
signed chirality relations
        |
        v
exact chirality-valid atom action
        |
        v
minimum fixed-mapping proper-fit RMSD
        |
        v
one selected R -> P_aligned mapping and its selected symmetry metadata
```

No geometry remapping or witness fallback is present in this decision flow.

## 2. Inputs, Graphs, and Tolerances

For each endpoint the algorithm receives:

- an element array;
- Cartesian coordinates;
- a full Wiberg bond-order matrix.

R and P must contain the same multiset of elements. Their atom indices do not
need to agree initially.

Two graph views are used:

1. **Active connectivity graph.** An edge exists when
   `WBO >= graph_floor`, normally `0.2`. This graph drives fragment growth and
   persistent-neighbor tests.
2. **Full WBO matrix.** This remains authoritative for weighted matching,
   event classification, and final scoring.

The current pipeline intentionally uses one WBO tolerance for edge matching
and pynauty WBO colors:

```text
symmetry_wbo_tolerance = iso_tol
```

The default BGCP value is `1.0`. Using a different hidden tolerance for
pynauty was one cause of missing degeneracy groups, so the pipeline now
normalizes both settings to the same value.

Bond events use a separate threshold, normally `dwbo_threshold = 0.5`, with
the configured metal-specific threshold where applicable.

## 3. Complete Pipeline

```text
                                  +-----------------------+
R elements/WBO/XYZ -------------->| build weighted R graph|
                                  +-----------+-----------+
                                              |
P elements/WBO/XYZ -------------->+-----------v-----------+
                                  | build weighted P graph|
                                  +-----------+-----------+
                                              |
                        no cut + selected R edge cuts
                                              |
                         +--------------------v--------------------+
                         | per cut: generate seed orders           |
                         | per seed: grow fragments sequentially   |
                         +--------------------+--------------------+
                                              |
                           compressed `_SymCand` fragment results
                                              |
                         +--------------------v--------------------+
                         | `find_islands`: build `_Branch` states  |
                         | exact live-state dedupe; subtree cap    |
                         +--------------------+--------------------+
                                              |
                              completed paths and hierarchies
                                              |
                         +--------------------v--------------------+
                         | classify broken/formed bonds            |
                         | canonicalize whole event set on R graph |
                         +--------------------+--------------------+
                                              |
                                  mechanism-keyed AAM pool
                                              |
                         +--------------------v--------------------+
                         | compile exact relational mapping cosets |
                         | remove equal/subsumed families          |
                         +--------------------+--------------------+
                                              |
                         +--------------------v--------------------+
                         | add local and group orientation         |
                         | solve exact oriented isomorphism        |
                         +--------------------+--------------------+
                                              |
                         +--------------------v--------------------+
                         | minimize fixed-mapping proper-fit RMSD  |
                         +--------------------+--------------------+
                                              |
                                selected complete R -> P mapping
```

Sweep cut is search orchestration, not part of the definition of a mapping.
For focused debugging the no-cut work unit can be run alone. Full production
mechanism discovery normally runs no-cut plus one work unit per eligible R
edge cut.

## 4. Fragment Candidate State

### 4.1 `_SymBlock`

A symmetry block represents a set-valued assignment domain:

```text
SymBlock
|- r_atoms: R roles participating in this domain
|- p_atoms: interchangeable target pool
`- extendable: whether later growth may add another R role
```

It does not mean that every factorial permutation is valid. Correlations are
resolved by the exact candidate automorphism group.

### 4.2 `_SymCand`

One compressed candidate stores:

```text
_SymCand
|- mapping: one deterministic representative R -> P assignment
|- blocks: open/closed assignment domains
|- exact_fixed: roles individualized by prior constraints
|- multiplicity: number of encountered states represented
|- automorph_blocks: connected display domains
`- automorph_generators: authoritative exact target permutations
```

The representative is needed to continue deterministic code, but it is not
treated as the only mapping represented by the candidate.

### 4.3 Candidate canonicalization

After extension, `_CandidateAutomorphismCanonicalizer` constructs a colored
pynauty graph containing:

- target atom elements and WBO buckets;
- fixed target atoms;
- distinct R-role colors;
- set-valued symmetry pools;
- already locked assignments.

Two candidates are merged only when equal canonical certificates and an exact
transporter prove that one full state maps to the other. The transporter is
retained as an exact automorphism generator. Merely placing two atoms in the
same vertex orbit is not sufficient because independent orbit labels discard
correlations.

## 5. Weighted Island Growth

`grow_island` starts from one R seed and uses a priority queue of R edges,
normally strongest WBO first. Queue order affects traversal, not validity.

For a proposed new R atom `n` and candidate target `v`, all active edges from
`n` to the already grown R fragment are checked:

```text
element_R[n] == element_P[v]
v is unused

for every grown r with WBO_R[n,r] >= graph_floor:
    WBO_P[v,m(r)] >= graph_floor
    abs(WBO_R[n,r] - WBO_P[v,m(r)]) <= iso_tol
```

R non-edges are not local negative constraints. An extra P bond may be a
formed bond and is classified at mechanism level.

If a touched R role belongs to a symmetry block, the support question is
existential and correlated:

```text
Does one injective assignment inside all touched blocks make the complete
active-edge WBO vector valid?
```

This small internal matching is bounded by `SYM_SUPPORT_MAX_STATES`. It is not
a global witness enumeration.

The growth transition is explicitly one of:

```text
0 outputs       -> defer the observed boundary; keep the candidate alive
1 output        -> commit the unique compressed state
many outputs    -> merge automorphic outputs; branch only on distinct states
```

At saturation, the fragment record retains:

- its R atom set;
- its representative assignments;
- exact target automorphism generators;
- symmetry domains;
- deferred boundary edges.

## 6. Deferred Boundaries

A fragment that cannot absorb a frontier atom is not necessarily invalid.
The failed weighted relationship is stored as a deferred boundary. This is
essential when an internally symmetric fragment has two sides but only one
side has already encountered another island.

Candidate equality therefore includes:

```text
exact internal pynauty certificate
              +
deferred one-hop boundary state
```

Without the second term, future-distinguishable branches would be collapsed.

## 7. Multi-Fragment AAM Branches

`find_islands` grows fragments sequentially for each seed order.

### 7.1 Live branch object

```text
_Branch
|- mapping: current concrete representative
|- islands_R: R atom -> island ID
|- islands_P: P atom -> island ID
|- deferred_edges
`- symmetry_paths[]
    `- ordered committed fragment records
```

An AAM branch is therefore a combination of fragment candidates, not a single
fragment and not merely a final mapping.

Hard `anchor_map` pairs are preloaded into the branch as fixed mapping/island
state. An anchored atom may seed growth, but it can also remain outside every
grown fragment. During analytical compilation, such an uncovered anchor is
represented as an individually fixed singleton fragment. Missing non-anchor
atoms remain an error. Anchor colors are also carried into the relational
graph, so family dedupe, chirality, symmetry repair, and RMSD selection cannot
move an anchored pair.

### 7.2 Exact live dedupe

Live branches merge only when these concrete cumulative states are equal:

```text
(mapping, islands_R, islands_P, deferred_edges)
```

Their distinct fragment histories are retained in `symmetry_paths`.

We tested a more aggressive coupled R/P automorphic live-state quotient. It
reduced TS04 from four correct mechanisms to two because the live graph did
not encode the full accumulated hierarchy. That optimization was removed.
Semantic automorphic dedupe is delayed until the completed hierarchy exists.

### 7.3 Branch cap

The configured BGCP cap is normally `max_branches = 100`. It applies to the
post-dedupe live leaves of one parent subtree:

```text
if accepting this parent's descendants would make live leaves > 100:
    discard only this overflowing descendant subtree
    retain siblings and other seed paths
```

Exactly 100 is legal. This is an intentional computational completeness cap,
not a symmetry equivalence rule.

## 8. Completed AAM Data Hierarchy

The raw pool and typed post-AAM model have the following ownership:

```text
CutSweepPool
`- mechanism entry, keyed by exact event certificate
   |- representative_mapping
   |- cuts / has_no_cut / encounter count
   `- completed branches[]
      `- AAMBranch
         |- representative_mapping
         |- encounter_count
         |- cuts
         |- hierarchy: AAMHierarchy
         |  `- fragments[]: FragmentMatch
         |     |- R atom set
         |     |- island ID
         |     |- deferred edges
         |     |- symmetry domains
         |     `- exact target generators
         |- analytical mapping-family record
         `- path provenance[]

PostAAMMechanism
|- exact mechanism key
|- endpoint R automorphism group   [auxiliary graph information]
|- endpoint P automorphism group   [auxiliary graph information]
`- maximal AAMBranch families[]    [authoritative mapping choices]
```

Endpoint automorphism groups are not automatically free mapping candidates.
They describe endpoint graph symmetry. Allowed mapping changes come from the
selected branch's analytical relation/coset.

### 8.1 Meaning of a witness

A witness is one complete bijection inside a family. It is useful for:

- continuing deterministic computation;
- recording provenance;
- initializing an exact isomorphism.

It is not sampled uniformly, is not a generator, and must not be ranked as if
the encountered witness list were the solution space.

## 9. Exact Mechanism Identity

For a complete mapping, bond events are computed in R index order:

```text
broken if WBO_R[i,j] - WBO_P[m(i),m(j)] >= event_threshold(i,j)
formed if WBO_P[m(i),m(j)] - WBO_R[i,j] >= event_threshold(i,j)
```

The entire broken/formed edge set is then attached to the full WBO-colored R
graph using typed event vertices. A pynauty certificate of this decorated
graph is the mechanism key.

```text
R graph atom vertices
      |
      +-- WBO-colored graph edges
      |
      +-- broken-event vertices -- broken type marker
      |
      `-- formed-event vertices -- formed type marker
                         |
                         v
               pynauty canonical certificate
```

This fixes a subtle but important bug: equal endpoint vertex-orbit IDs do not
prove that two edges lie in the same edge orbit, and they certainly do not
prove that two multi-edge event sets are equivalent.

Mechanism selection keeps the classes with the minimum number of bond-breaking
plus bond-forming events. Different exact event certificates at that minimum
remain distinct mechanisms.

## 10. Analytical Mapping Families

Every completed branch is compiled into a colored relational isomorphism
between endpoint A and endpoint B.

The relation includes:

- atom element and fragment-owner colors;
- selected fragment WBO relations at `iso_tol`;
- anchors;
- event-invariant pair colors;
- later, signed orientation relations.

If `g` is one isomorphism and `G` is the target automorphism group of the
relation, the branch represents the coset:

```text
F = gG
```

`AnalyticalMappingFamily` stores:

```text
AnalyticalMappingFamily
|- source_mapping
|- representative_mapping g
|- target_generators of G
|- target_orbits
|- group_order
`- colored relational records for exact membership
```

### 10.1 Membership

For a proposed complete mapping `m`, `contains(m)` directly transports all A
atom colors and relation records through `m` and compares them with B. It does
not invoke a new geometry matching operation.

### 10.2 Family inclusion and dedupe

For finite cosets, `F1` is proven to be a subset of `F2` by checking:

1. `F1`'s representative belongs to `F2`;
2. applying every generator of `F1` to that representative remains in `F2`.

Families are processed into a maximal antichain:

```text
equal family       -> merge provenance
strict subset      -> attach its provenance to the containing family
strict superset    -> replace contained families and inherit provenance
incomparable       -> retain both
```

Thus 82 TS01 growth paths become one maximal analytical family without
pretending the paths themselves are group elements. Large cases may retain
several incomparable families under the same mechanism.

## 11. Index Chirality

Index chirality is a correspondence constraint. It asks whether the ordering
of mapped substituent indices is orientation-consistent between endpoints. It
does not assign chemical R/S labels and does not alter AAM chemistry.

### 11.1 Affine simplex sign

For center `c` and ordered neighbors `(a,b,d)`:

```text
v1 = xyz[a] - xyz[c]
v2 = xyz[b] - xyz[c]
v3 = xyz[d] - xyz[c]

s = sign(det([v1, v2, v3]))
```

The determinant is normalized by vector lengths before applying the
degeneracy threshold. Ordinary local centers use the configured dimensionless
near-planarity tolerance.

### 11.2 Local persistent centers

Persistent neighbor simplices are found from the selected branch relation,
not from display color groups. Pynauty stabilizer orbits determine whether a
center's neighbors are actually movable. Signed ordered-simplex relation
vertices are then added to A and B.

The oriented graph is solved as one simultaneous isomorphism problem. There
is no degree-four-only swap rule and no sequential local shuffle.

### 11.3 Higher-coordinate group orientation

For centers with more than four ligands, orientation is a relation among
ligand triples. Physical endpoint geometries can reconfigure so one dependent
triple crosses coplanarity even while the overall ligand assignment remains
consistent. Requiring every one of `C(k,3)` signs as a hard constraint can
therefore reject a valid family.

The current algorithm builds a maximal feasible signed-frame basis:

1. construct all defined group-level triples;
2. rank them by endpoint-normalized geometric robustness;
3. add each complete ordered-sign relation to a trial relational graph;
4. retain it only if the cumulative graph still admits an exact isomorphism;
5. record incompatible dependent triples as geometric reconfiguration.

This is not a witness fallback. Every retained constraint is solved against
the complete analytical family. The excluded frame is explicitly reported.

PR8 demonstrates the distinction: 19 R48 frames are simultaneously
preserved, while the nearly coplanar `[34,41,43]` frame is recorded as
reconfigured instead of invalidating the entire mapping family.

## 12. RMSD Selection Inside the Exact Family

After orientation relations are added, pynauty returns generators of the
chirality-valid target action. Relation-vertex-only kernel permutations are
discarded by restricting each generator to atom vertices.

For any candidate mapping, RMSD uses immutable correspondence:

```text
P_R_order[r] = xyz_P[m(r)]

center R and P_R_order
compute proper Kabsch rotation, det(rotation) = +1
RMSD = sqrt(mean(||R - rotated(P_R_order)||^2))
```

Kabsch removes only global translation and proper rotation. It never performs
assignment, symmetry matching, or atom remapping.

### 12.1 Small groups

When the chirality-valid atom action contains at most 4096 elements, all
actions are evaluated in one bounded vectorized batch. This is faster than
many small Python/SVD calls and cannot grow past the explicit threshold.

### 12.2 Large factorizable groups

Generator supports are joined when they overlap. Disjoint support components
are exact commuting factors:

```text
G = G1 x G2 x ... x Gk
```

The global Cartesian product is searched without materializing it. A greedy
descent supplies only an initial incumbent. Exact branch-and-bound then uses
the rotation-invariant lower bound for already assigned atoms:

```text
RMSD >= sqrt(sum_(i<j) (distance_R(i,j)
                         - distance_P(m(i),m(j)))^2) / N
```

If this lower bound is worse than the incumbent, the entire remaining coset
subtree is skipped. Tie-breaking remains deterministic by rounded RMSD and
the complete mapping tuple.

The search is exact. Its present worst-case risk is a very large connected
support factor: local factor actions are still closed explicitly. The global
independent-factor product no longer needs to be enumerated, but a future
Schreier-Sims representation would be needed to remove that final worst-case
group-size risk.

## 13. Post-AAM Parallelism and Performance

Post-processing previously recomputed endpoint-only event behavior for every
growth path. On a 133-atom case this performed roughly 830,000 repeated event
comparisons per family.

The current immutable compiler context precomputes once:

- endpoint active graphs;
- element/threshold pair classes;
- R-side event behavior vectors;
- P-side event behavior vectors.

It is shared by family compilation and chirality evaluation. Exact relation
records are cached for containment checks. Families are compiled in process
batches, with up to 32 workers for large branch sets, and the remaining
maximal families are evaluated independently in parallel.

Measured on the 133-atom Pd TS12 case with eight CPUs:

```text
initial exact post-AAM: 29.02 s
current exact post-AAM: 11.5-11.6 s
```

The old witness baseline was faster because it proved less and sampled
witnesses. The current result retains the exact family and gives lower RMSD.

## 14. Bounded Symmetry Repair During AAM Scoring

Before analytical-family post-processing, completed AAM mappings may undergo
a bounded symmetry repair inside exact touched target subgroups. The touched
atoms are derived from current event endpoints. The repair never imports an
unrelated target or performs geometric remapping; it composes the mapping with
exact pynauty subgroup actions and scores bond-event count/WBO change.

`symmetry_repair_max_evals`, normally 20000, is a hard diagnostic cap. This
step normalizes a concrete completed representative for event scoring. It does
not replace analytical family compilation or final chirality/RMSD selection.

## 15. Viewer Semantics

The viewer must display only the selected solution:

- one selected mechanism at a time;
- one selected analytical branch/family;
- its selected complete R mapping and aligned P mapping;
- symmetry/degeneracy derived from that branch's exact fragment
  automorphisms;
- atoms that are actually mutable under the selected relation.

It must not aggregate color groups from rejected paths or treat endpoint
orbits as allowed mapping shuffles.

Viewer interpolation is a validation layer, not an alignment step. It uses the
already selected indices. A collision or path crossing in interpolation does
not by itself prove that AAM is wrong; mapping chirality, endpoint pose, and
interpolation constraints must be diagnosed separately. Self-contained HTML
views and per-mechanism `R.xyz` / `P_aligned.xyz` remain the portable debugging
artifacts.

## 16. Failure Rules and Prohibited Shortcuts

The current code deliberately avoids the following shortcuts:

- selecting a random or first witness;
- treating witnesses as generators;
- combining independent atom orbits as if their swaps were uncorrelated;
- deduplicating mechanisms from endpoint orbit pairs;
- merging live automorphic branches without encoding their hierarchy;
- geometry-based remapping before RMSD;
- accepting chirality through a fallback mapping outside the AAM family;
- silently changing WBO tolerance between edge verification and pynauty.

An analytical branch with no oriented isomorphism is rejected with diagnostics.
If every minimum-event mechanism is rejected, the pipeline raises an
`IndexChiralityConflict`; it does not silently return an unverified mapping.

## 17. Important Parameters

| Parameter | Typical value | Role |
|---|---:|---|
| `graph_floor` | `0.2` | active connectivity threshold |
| `iso_tol` | `1.0` | weighted edge tolerance and pynauty WBO color tolerance |
| `dwbo_threshold` | `0.5` | broken/formed bond threshold |
| `n_seeds` | `3` | seed orders per cut work unit |
| `max_branches` | `100` in BGCP | post-dedupe live leaves per parent subtree |
| `SYM_SUPPORT_MAX_STATES` | `4096` | local correlated block-support cap |
| `symmetry_repair_max_evals` | `20000` | bounded completed-representative repair |
| RMSD vector batch threshold | `4096` | largest explicitly materialized global atom action |
| analytical compile workers | up to `32` | process parallelism for large branch sets |

## 18. Verification Record

The current implementation is covered by 133 automated tests. Important
checks include:

- cached and uncached relational graphs are identical;
- a generated 8192-action group gives the same selected mapping and RMSD under
  branch-and-bound as exhaustive enumeration;
- TS01 retains one mechanism, one maximal family, and all 82 paths;
- TS04 retains all four exact mechanisms and 2187 chirality-valid mappings per
  mechanism;
- PR9 TS41a-endo matches the prior corrected mapping and event;
- PR8 retains 19 compatible higher-coordinate frames and records one
  reconfigured frame;
- 133-atom Pd TS12 retains both mechanisms and all four maximal families;
- 133-atom Pd TS14 retains both concrete mechanisms and seven maximal
  families;
- 95-atom Noyori TS65 completes with two exact mechanisms and zero chirality
  violations.

These checks are a regression sample, not a substitute for rerunning the full
140-case batch after future changes to search equivalence, mechanism
certificates, or chirality constraints.
