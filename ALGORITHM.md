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

## Public computational architecture

The public Python interface is typed and immutable. Dictionary records are
not computational inputs or outputs; they exist only inside legacy artifact
adapters. The dependency direction is deliberately one-way:

```text
MolecularEndpoint(R, P)
        |
        v
 AAMProblem + AAMSearchConfig
        |
        v
 search_aam -------------------------------> AAMResult
                                                |
                         +----------------------+
                         | mechanism classes
                         | completed branches
                         | fragment hierarchy
                         | exact fixed roles
                         | assignment domains
                         | target generators
                         | cuts + provenance
                         | search metrics
                         v
             compile_mapping_families
                         |
                         v
                 AnalyticalAAMResult
                         |
                         v
                 select_rp_mappings
                         |
                         v
                       RPResult
                         |
             +-----------+-----------+
             |                       |
             v                       v
      aligned endpoints      analyze_transition_state
                                     |
                    R->TS partial AAM + P->TS partial AAM
                                     |
                         exact core-tuple consensus
                                     |
                                     v
                                  TSResult
```

`search_aam()` performs AAM graph search, mechanism classification, and no
geometry-based selection.
It also finalizes each retained fragment candidate's exact target generators
after branch-family reduction; this is the terminal AAM operation, and its
request/calculation/cache counts are part of `AAMSearchMetrics`.
`compile_mapping_families()` turns the retained hierarchical relations into
maximal exact cosets. `select_rp_mappings()` applies chirality constraints and
then fixed-mapping RMSD ranking. `analyze_transition_state()` composes two
partial AAM searches with mode scoring. None of these stages reads JSON,
writes files, invokes the viewer, or silently recomputes an alternative atom
mapping.

The main result objects are:

```text
AAMResult
`- AAMMechanism[]
   `- AAMBranch[]
      |- representative: AtomBijection
      |- hierarchy: AAMHierarchy
      |  `- FragmentMatch[]
      |     |- representative_assignments
      |     |- symmetry_domains
      |     |- exact_fixed
      |     |- automorph_domains
      |     `- target_generators
      |- cuts and encounter counts
      `- path provenance

AnalyticalAAMResult
`- AnalyticalMechanism[]
   `- AnalyticalBranch[]
      |- original AAMBranch
      `- exact AnalyticalMappingFamily

RPResult
`- RPMechanism[]
   |- selected AtomBijection
   |- broken/formed bonds and core atoms
   |- chirality audit
   `- exact fixed-mapping RMSD

TSResult
`- TSMechanismResult[]
   |- reactant_core_aam: CoreAAMResult
   |- product_core_aam: CoreAAMResult
   |- exact scored core assignments
   `- selected TSScore
```

The TS core search does not compress a multi-atom assignment to independent
vertex orbits. A `CoreAAMBranch` retains the branch hierarchy and its exact
correlated tuple orbit; `CoreAAMResult.assignments` is the deduplicated union
of complete tuples. This is essential because membership of two target atoms
in vertex orbits does not prove that their joint permutation is an
automorphism.

The primary implementation is in:

- `src/rxn_core/matcher/`: compressed fragment candidate matching;
- `src/rxn_core/growth/`: weighted island growth;
- `src/rxn_core/alignment/branch.py`: multi-island branch construction;
- `src/rxn_core/alignment/sweep.py`: cut sweep and mechanism grouping;
- `src/rxn_core/alignment/post_aam.py`: typed post-AAM data model;
- `src/rxn_core/alignment/index_chirality.py`: analytical families,
  chirality, and fixed-mapping RMSD selection;
- `src/rxn_core/aam.py`: typed full AAM search;
- `src/rxn_core/analytical.py`: exact coset compilation and containment;
- `src/rxn_core/rp.py`: R/P composition;
- `src/rxn_core/core_aam.py`: exact partial AAM for mechanism cores;
- `src/rxn_core/ts.py`: typed TS composition and mode scoring;
- `src/rxn_core/pipeline.py`: legacy artifact/CLI adapter, not a public
  computational contract.

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
                        | reuse finalized AAM generators          |
                        | filter exact local actions by chirality |
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

Two nested relations are compiled from the same completed AAM branch:

- atom element and fragment-owner colors;
- selected fragment WBO relations at `iso_tol`;
- anchors;
- the structural relation stops here and represents the pre-event AAM group;
- the conservative event stabilizer additionally contains event-invariant
  pair colors;
- signed orientation relations are added later during chirality selection.

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
|- structural_target_generators
|- structural_group_order
`- colored relational records for exact membership
```

The event-invariant colors preserve threshold behavior against every endpoint
WBO value. That makes them a valid subgroup stabilizer, but occasionally
stricter than preserving the concrete event realized by the selected mapping.
The algorithm therefore traverses the finite quotient of the structural AAM
group by this conservative stabilizer. It retains one representative for each
quotient coset whose actual broken/formed event equals the mechanism event.
Only these event-family representatives proceed to chirality and RMSD.

This is a Schreier-style traversal of event cosets, not enumeration of group
elements or atom bijections. TS01, for example, has two same-event quotient
cosets even though the conservative stabilizer contains only one of them.

### 10.1 Membership

For a proposed complete mapping `m`, `contains(m)` directly transports all A
atom colors and relation records through `m` and compares them with B. It does
not invoke a new geometry matching operation.

### 10.2 Family inclusion and dedupe

For finite cosets, `F1` is proven to be a subset of `F2` by checking:

1. `F1`'s representative belongs to `F2`;
2. applying every generator of `F1` to that representative remains in `F2`.

Identical mapping/fragment payloads are quotient-collapsed before family
compilation. Compiled families are then reduced in two exact phases:

1. canonical relation-certificate buckets prove and merge equal cosets;
2. families are visited in descending group order, using orbit compatibility
   to reject impossible containment before relation transport.

The result is a maximal antichain:

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

For a four-ligand affine tetrahedron, normalization uses the geometric mean
of all six ligand-ligand edge lengths. This scale is invariant to which ligand
is written first; ordered signs are derived from one canonical determinant and
the ligand permutation parity. Consequently, every ordering of one geometric
tetrahedron has the same degeneracy classification and normalized magnitude.

### 11.2 Local persistent centers

Persistent neighbor simplices are found from the selected branch, not from
display color groups. Mutability comes from two exact AAM sources:

1. point stabilizers of the finalized fragment automorphism generators;
2. correlated assignments represented by the maximal branch families.

A center is constrained from the intersection of its mapped R neighbor shell
and its P neighbor shell. Ligands on broken or formed coordination edges are
excluded, while the surviving ligand assignment remains available as an
index-orientation constraint. The persistent intersection is constrained only
when the exact AAM family permits a nontrivial setwise permutation within it.

For a persistent intersection of three or four ligands, its complete affine
orientation is mandatory even when total coordination changes, such as 5→4
or 4→5. If more than four ligands persist, all defined affine simplices are
constructed for the dependent-basis treatment below. Thus coordination change
does not erase surviving chirality, and it also does not impose orientation on
the departing or arriving ligand.

An ordinary persistent three- or four-coordinate frame is a hard constraint.
If its orientation-colored subgroup is empty, that event coset fails; it can
never be relabeled as geometric reconfiguration. Only dependent simplices at
centers above coordination four participate in the maximal-feasible-basis
rule below.

There is no degree-four-only swap rule, sequential greedy shuffle, or sampled
witness fallback.

### 11.3 Higher-coordinate group orientation

For centers with more than four ligands, orientation is a relation among
ligand triples. Physical endpoint geometries can reconfigure so one dependent
triple crosses coplanarity even while the overall ligand assignment remains
consistent. Requiring every one of `C(k,3)` signs as a hard constraint can
therefore reject a valid family.

The current algorithm builds a maximal feasible signed-frame basis directly
inside the stored AAM action:

1. install every ordinary persistent three/four-ligand frame as a hard
   simultaneous constraint;
2. construct symmetry-closed high-coordinate simplex units;
3. rank dependent units and group-level triples by endpoint-normalized
   geometric robustness;
4. filter the still-valid exact action cumulatively by each unit;
5. retain a unit only when at least one cumulative action remains;
6. record incompatible dependent units as geometric reconfiguration.

This priority is invariant across both execution routes: a dependent
high-coordinate frame can never consume the freedom required by an ordinary
hard frame.

This is not a witness fallback. Every retained constraint is solved against
the finalized AAM group. The excluded frame is explicitly reported.

PR8 demonstrates the distinction: 19 R48 frames are simultaneously
preserved, while the nearly coplanar `[34,41,43]` frame is recorded as
reconfigured instead of invalidating the entire mapping family.

## 12. RMSD Selection Inside the Exact Family

Every maximal `AAMBranch` carries the authoritative compiled mapping family:
its representative, exact target generators, target orbits, group order, and
the in-memory colored relation used to prove family equality/containment.
Post-AAM selection reuses this object; endpoint chemistry and fragment
symmetry are not recomputed.

There are two exact execution strategies. SymPy Schreier-Sims first measures
the group and overlapping-support component orders without enumerating group
elements:

1. when the complete action is at most 1,000,000 and every component is at
   most 4,096, components are closed and filtered directly;
2. for a larger entangled action, chirality colors are added to the already
   compiled AAM relation and pynauty returns the exact chirality subgroup.

The second route is not a fallback and does not rematch atoms. It is the
algebraically appropriate representation for wreath-product-like groups where
materializing one support component would itself be combinatorial. Both
routes operate on the same AAM family and return an exact subgroup/coset.

Factors touching chirality frames, a hard anchor, or a potentially changing
bond event are filtered explicitly. Event sensitivity is tested over every
action of a local factor, not just the supplied generator list, so a generator
product cannot silently change the chosen mechanism. Factors independent of
all constraints remain as generators for exact RMSD minimization.

Completed branch representatives contribute correlated mutability detection,
but they are not treated as random RMSD candidates. Each maximal family is
evaluated independently across its same-event structural quotient cosets;
within each coset, only exact chirality-valid group actions are scored.

For any candidate mapping, RMSD uses immutable correspondence:

```text
P_R_order[r] = xyz_P[m(r)]

center R and P_R_order
compute proper Kabsch rotation, det(rotation) = +1
RMSD = sqrt(mean(||R - rotated(P_R_order)||^2))
```

Kabsch removes only global translation and proper rotation. It never performs
assignment, symmetry matching, or atom remapping.

### 12.1 Covariance representation

For centered coordinates, the coordinate norms are invariant under every
permutation. Proper-fit RMSD therefore depends only on the 3x3 covariance:

```text
C(m) = sum_r outer(P[m(r)], R[r])
RMSD(m)^2 = (||R||^2 + ||P||^2 - 2 proper_score(C(m))) / N
```

No global mapping list is constructed.

### 12.2 Exact factor search

Generator supports are joined when they overlap. Disjoint support components
are exact commuting factors:

```text
G = G1 x G2 x ... x Gk
```

Each local action contributes one additive 3x3 covariance matrix. Local
action matrices are stored in a binary tree; every tree node is enclosed by a
rigorous Frobenius ball. A greedy descent supplies only an incumbent. For a
partial covariance and all remaining action balls:

```text
proper_score(C + remaining)
    <= proper_score(C + sum(ball centers))
       + sqrt(3) * sum(ball radii)
```

This is an upper bound on the best possible Kabsch score and therefore a
lower bound on RMSD. If it cannot improve the incumbent, the complete group
subtree is discarded. Tie-breaking remains deterministic by rounded RMSD and
the complete mapping tuple.

The search is exact. On the direct route, local connected-factor actions are
closed explicitly but the global product is never enumerated as atom
bijections. On the entangled route, orientation restriction first reduces the
compiled family (for example, the Ni TS11 family becomes 324 valid actions),
then the same exact covariance search scores the reduced action.

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
batches, with up to 48 workers for large branch sets.

The structural family `G` can contain several right cosets of the conservative
event stabilizer `K`. Their quotient is constructed directly with a
Schreier--Sims transversal, never by quadratic pairwise relation-membership
tests. Concrete event signatures for the quotient representatives are checked
in bounded vectorized batches. A retained member `K g` reuses the already
compiled target subgroup `K`; changing `g` does not justify recompiling the
pynauty relation. The compiled subgroup also proves that its complete action
preserves the selected event, so only the final selected mapping requires an
independent event assertion.

After quotient construction, `(branch, event-coset)` pairs are independent
exact work units. They are scheduled across up to 48 workers and reduced by
the deterministic global order `(RMSD, mapping tuple, branch, coset)`. This is
only a scheduling transformation: no family, constraint, or RMSD candidate is
discarded.

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
| direct component order | `4096` | direct exact action-closure strategy limit |
| direct total group order | `1,000,000` | direct vs compiled-relation exact strategy |
| `symmetry_repair_max_evals` | `20000` | bounded completed-representative repair |
| analytical compile workers | up to `48` | process parallelism for large branch sets |

## 18. Verification Record

The current implementation is covered by 154 automated tests. Important
checks include:

- cached and uncached relational graphs are identical;
- a generated 8192-action group gives the same selected mapping and RMSD as
  exhaustive enumeration while evaluating one complete mapping and proving
  the other 8191 cannot win;
- TS01 retains one mechanism, one maximal family, and all 82 paths;
- TS04 retains all four exact mechanisms and 2187 chirality-valid actions per
  mechanism, with zero violations;
- PR9 TS41a-endo matches the prior corrected mapping and event;
- PR8 retains 19 compatible higher-coordinate frames and records one
  reconfigured frame;
- 133-atom Pd TS12 retains both mechanisms and all four maximal families;
- 133-atom Pd TS14 retains both concrete mechanisms and seven maximal
  families;
- 95-atom Noyori TS65 completes with two exact mechanisms and zero chirality
  violations.

The direct stored-group regression gates additionally show:

- TS01: one mechanism, one maximal family, all 82 growth paths, and a
  self-contained viewer;
- PR9 TS41a-endo: unchanged event, zero violations, and 3.47 seconds total;
- Fe TS2 (83 atoms, 16 workers): 23.9 seconds total, with chirality/RMSD
  reduced from 792 seconds to 0.98 seconds;
- Fe TS8 (82 atoms, 16 workers): 313 seconds total versus 823 seconds before,
  with post-AAM reduced to 13.7 seconds;
- concrete Fe TS8 events using R15 and R26 have the same exact WBO-colored
  mechanism certificate; regression checks therefore compare canonical event
  certificates rather than arbitrary symmetry-equivalent atom labels.

The final revision `3b1e34c` was rerun over the complete 140-case manifest:

- 140/140 cases succeeded, producing 166 mechanisms;
- every selected mechanism reports zero index-chirality violations;
- all 139 cases shared with the previous analytical batch retain their
  mechanism count;
- two concrete event representatives changed atom labels, and both pairs have
  identical exact mechanism certificates;
- 114 mechanisms used direct stored-group actions and 52 used the exact
  compiled-relation subgroup route;
- summed per-case wall time fell from 10,652 s to 9,165 s on the 139 directly
  comparable cases; median time ratio was 0.988;
- the compact self-contained viewer contains all 140 cases, and the aligned
  package contains 332 XYZ files (`R.xyz` and `P_aligned.xyz` for each of 166
  mechanisms).

These checks are a regression sample, not a substitute for rerunning the full
140-case batch after future changes to search equivalence, mechanism
certificates, or chirality constraints.
