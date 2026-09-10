# Scientific novelty audit: continuous fragment-growth atom mapping

Assessment date: 10 September 2026. The algorithm is called GRAFT in the current
draft; naming is still provisional. This assessment concerns the frozen benchmark
source `98b01b175eeed31f70d13e7cbf178b80bf07c9e0`, rather than subsequent development
changes. No benchmarks were rerun.

**Assessment.** The implementation expresses a specific, custom search design.
This review did not identify an earlier publication or inspected implementation
with the same complete combination of shared-fragment growth, compressed live
placements, saturation-triggered whole-mapping branches, continuation-sensitive
state joins, and correlated completed families. That makes a new algorithmic
variant plausible. It does not establish historical priority or prove that the
combination is absent from the literature. Most ingredients, including several
combinations of them, have substantial precedents.

Calling it simply a reimplementation of SLAP, VF2, or maximum common subgraph
search would miss material differences in its control flow and represented
state. Calling fragment growth, branching, symmetry reduction, edge cutting,
or search-graph sharing individually new would also be unsupported. The strongest
candidate contribution is the precise coupling and representation of these
operations, with empirical evidence for its usefulness.

**What the frozen implementation actually does.** These observations come from
source inspection, rather than the proposed method name or benchmark scores.
Links below pin the source revision; the accompanying
[source manifest](novelty-source-manifest.json) records file hashes.

| Operation | Evidence in the implementation | Implication for comparison |
| --- | --- | --- |
| Shared-fragment growth | [`grow_island`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/growth/island.py#L248) proposes one source extension for the current placement collection, accepts it if any placement supports it, and carries the surviving collection forward. | This is a greedy trajectory over a common fragment with multiple placements. It does not independently backtrack over every possible fragment boundary. |
| Continued extension and saturation | The same loop continues until the frontier is exhausted. Uniqueness is evaluated after the loop; remaining distinct candidates become returned placements. | A unique placement is not an early stopping condition. Local placement branching occurs during extension; subsequent whole-mapping branching occurs on saturated placements. |
| Context-preserving candidate deduplication | [`_CandidateAutomorphismCanonicalizer`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/matcher/canonical.py#L288) colors the full target graph with locked images, candidate roles, and pools; [`dedupe.py`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/matcher/dedupe.py) also records boundary information. | The equivalence relation concerns a contextual candidate, rather than a bare atom set or independent atom orbits. The exactness scope is the configured graph coloring and supported relation. |
| Conservative state sharing | [`_Branch.state_key`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/alignment/branch.py#L102) includes assignments, the normalized island partition, and deferred edges. Admission compares keys at the same scheduler step, and [`merge_exact_paths`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/alignment/branch.py#L98) retains histories. | Equal assignments alone are insufficient; joins are conditional on the continuation context. This is narrower than quotienting every partial map by endpoint symmetry. |
| Completed mapping families | [`AnalyticalMappingFamily`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/alignment/index_chirality.py#L888) represents a completed relation through an isomorphism coset and target generators. [`_maximal_families`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b1/src/rxn_core/analytical.py#L107) merges equal/contained families and retains provenance. | Returned alternatives retain correlations. Completed-family containment is a different operation from live-state deduplication. |

The benchmark path dispatches to the native growth implementation when applicable;
the inspected Python growth path specifies the event-replay behavior. This audit
does not add a new native/Python parity proof. Source comments asserting exactness
are not, by themselves, mathematical proofs of the complete bounded algorithm.

**Closest prior methods and the distinctions that matter.** Descriptions below
are limited to the primary material actually inspected. Differences are this
audit's interpretation, not claims made by those authors about GRAFT.

| Prior work | Documented overlap | Distinction or consequence for GRAFT |
| --- | --- | --- |
| **SMSD, Rahman et al., 2009** | Combines molecular subgraph algorithms and extends approximate common substructures using McGregor search. Multiple matches and chemical filtering are already part of the design. [Primary article](https://link.springer.com/article/10.1186/1758-2946-1-12). | Growing partial molecular matches is established. GRAFT's shared-fragment placement collection and staged contextual representation require a more specific comparison than “we grow fragments.” |
| **McSplit, McCreesh, Prosser and Trimble, 2017** | Incrementally builds a mapping using branch and bound, stores candidate domains compactly as paired vertex classes, and refines those classes after assignments. [Paper and Algorithm 1](https://www.ijcai.org/proceedings/2017/0099.pdf). | Compact alternatives during graph search are established. McSplit targets maximum common induced subgraphs and branches on vertex assignments, including unmatched choices. GRAFT uses weighted source-edge compatibility and greedy shared-fragment extension before whole-mapping branching. |
| **Jaworski et al., 2019** | Extends matches through atom neighborhoods, searches a bounded decision tree, and cuts bond subsets starting with single bonds and increasing up to six when needed. Reaction heuristics guide candidate assignments. [Primary article, “Establishing isomorphic mapping”](https://www.nature.com/articles/s41467-019-09440-2). | This is a particularly close precedent for growth + branching + cuts. The described neighborhood-order procedure differs from GRAFT's continued weighted frontier and compressed placements. The full GRAFT state/family combination was not identified in the inspected main article; its supplementary PDF was unavailable in this audit. |
| **Mann et al., 2014** | Uses constraint programming to find cyclic imaginary-transition-state candidates, removes hydrogen-related symmetries, and extends candidates to full atom mappings by graph matching. [Primary article](https://link.springer.com/article/10.1186/s13015-014-0023-3). | Search domains, chemical constraints, partial-map extension, and symmetry handling are established in atom mapping. Its reaction-center/ITS-first formulation differs from GRAFT's fragment-growth search. |
| **Ali et al., 2025** | Enumerates maximum-common-edge-subgraph-based mappings and clusters complete mappings by reactant/product symmetry, retaining representatives of distinct reaction patterns. [Primary article, Section 2.2](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01871). | Discovering symmetry-distinct alternatives is already an explicit contribution in atom mapping. GRAFT's possible distinction concerns live compression and correlated family output, rather than the general motivation to return alternatives. |
| **SLAPMapper, Koda and Saito, 2025** | The frozen author implementation branches over distinct LAP solutions, tracks visited refined-label states, and filters isomorphic completed outputs. [Author source, `SlapMapper`](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py#L168). | Branching and deduplication are not absent from SLAP. The useful contrast is sequential assignment/label refinement versus continued fragment growth and its contextual placement/family representation. The preprint PDF was inaccessible; the benchmarked source was inspected directly. |
| **Laffitte, Phan and Stadler, 2026; preceding WABI 2025 paper** | Formalizes stable extensions of partial atom maps, completion by anchored VF2 or relabeling/isomorphism, and equivalence/co-extension using graph automorphisms. [Primary article](https://link.springer.com/article/10.1186/s13015-026-00300-5). | A necessary comparison for partial maps, fixed context, and hydrogen extensions. Its anchored completion problem starts with a supplied partial map; GRAFT searches for fragment placements and boundaries. Stable-extension uniqueness is up to the paper's equivalence, not unique atom indices. |
| **SynKit v1.5.0, July 2026** | Combines WL/SLAP mapping with uncertainty-region refinement, orbital branching, symmetry-distinct enumeration, certificates, and canonical keys for full and partial mappings. [Archived release](https://zenodo.org/records/21293638), [pinned branching source](https://github.com/TieuLongPhan/SynKit/blob/b57f2a39f1807a2682913183b22e1d8d61531bba/synkit/Chem/Reaction/Mapper/exact/branching.py), [partial-map canonicalization](https://github.com/TieuLongPhan/SynKit/blob/b57f2a39f1807a2682913183b22e1d8d61531bba/synkit/Chem/Reaction/Mapper/exact/enumerate.py#L849). | Another close software precedent. The inspected solver refines an uncertainty kernel from SLAP rather than using GRAFT's saturation-based growing-fragment scheduler. Its presence rules out treating partial-map canonicalization or certificates alone as a distinctive new feature. This audit does not certify its stated solver guarantees. |

**General techniques also have relevant precedents.** OpenChemLib 2026.2.0's
[`setFragmentSymmetryConstraints`](https://github.com/Actelion/openchemlib/blob/0754c781c23f7570b13866e202836b5449282ad6/src/main/java/com/actelion/research/chem/SSSearcher.java#L244)
explicitly distinguishes otherwise symmetric fragment atoms using reaction
context. That is not the same state representation as GRAFT's locks, pools,
boundaries, and continuations, but context-sensitive symmetry handling is not a
new general idea.

Dechter and Mateescu's
[AND/OR search paper](https://www.ics.uci.edu/~csp/r126.pdf) develops merging and
context-based search graphs, including the preservation of incoming paths through
shared continuations. GRAFT is not thereby an implementation of their complete
algorithm; the connection establishes that replacing redundant search-tree
continuations with shared graph states is an established principle.

[McKay and Piperno's nauty/Traces work](https://arxiv.org/abs/1301.1493) supplies
established graph-canonicalization and automorphism machinery. Likewise, if one
isomorphism is known, representing other isomorphisms through automorphism actions
is a standard group-theoretic construction. GRAFT's use of these tools should be
credited as such. A contribution may lie in the encoded relation and how it is
maintained during fragment decisions, rather than in inventing canonical labeling
or cosets.

Additional screening found [ReactionMap, 2013](https://pubs.acs.org/doi/10.1021/ci400326p)
(common chemical subgraphs plus an assignment cost),
[CLCA, 2014](https://pubmed.ncbi.nlm.nih.gov/25412255/) (canonical-label-based
substructure mapping), and
[Heinonen et al., 2011](https://journals.sagepub.com/doi/10.1089/cmb.2009.0216)
(minimum-edit atom mapping). These were used as background leads at abstract level,
not as full algorithm exclusions. The 2026
[NEBscape article](https://www.nature.com/articles/s41524-026-02298-1) also describes
retaining distinct atom-mapping classes and reducing symmetry-equivalent variants
before transition-state searches. That is relevant to downstream motivation, not
evidence that it implements GRAFT's growth loop.

**The most defensible candidate contribution.** The following description is
supported by the implementation and avoids a claim of universal priority:

> We present an atom-mapping algorithm that couples continuous growth of a shared
> fragment with symmetry-compressed candidate placements, branches on retained
> placements at saturation, and preserves contextual search histories and
> correlated mapping families through staged deduplication.

The distinctive point to develop is the interaction between these choices.
For example, a source extension can eliminate some placements while allowing
others to continue on the same enlarged fragment. Later branches must inherit
the correct locks and deferred boundaries; states may join only when their
continuation is compatible; compiled families must preserve joint atom actions.
An implementation that merely enumerates graph matches and removes duplicate
completed mappings does not automatically provide these semantics.

This difference does not prove novelty. A sufficiently general constraint or
graph-search solver can express many specialized algorithms. To justify a new
algorithmic variant, the paper should explain the concrete representation,
transition rules, tradeoffs, and useful behavior of this specialization.

**Claims to avoid on the current evidence.** The search does not support “first
fragment-growth atom mapper,” “first symmetry-aware mapper,” “first mapper to
return alternatives,” “first use of bond-cut relaxation,” or “first deduplicated
atom-mapping search.” It also does not establish global completeness, global
optimality, a new graph-isomorphism algorithm, or chemical validity of every
returned alternative. The absence of an exact match in this review should not
be rewritten as proof that none exists.

**Evidence that would strengthen the algorithm paper.** These are proposed
follow-ups, not experiments performed by this audit.

1. Specify the common-fragment/placement-set transition precisely, including
   greedy acceptance, incompatible-placement removal, deferral, island absorption,
   and saturation. Compare that transition with vertex-by-vertex branching and
   neighborhood-extension methods.
2. State and justify the invariants needed by each deduplication level. In
   particular, explain why its equivalence preserves permitted continuations,
   and why completed-family containment is not used to prune arbitrary live
   partial states. Separate graph-coloring exactness from continuous tolerance
   checks and finite support budgets.
3. Use small mechanistically interpretable examples to show what fails if
   growth stops at a unique placement, only one saturated placement survives,
   or correlations are replaced by independent atom pools. Do not present the
   previously excluded synthetic six-vertex diagnostic as a paper result.
4. Run matched component ablations in a separate experimental revision, keeping
   cut policy, seed orders, budgets, and evaluator fixed where possible. Report
   recovery, distinct families, live-state counts, memory, and CPU. Disabling
   exact duplicate elimination may preserve unbounded results but change which
   branches survive finite caps; make that distinction explicit.
5. Keep the SLAP sweep as the supporting cut-policy ablation requested by the
   author. Its recovery gain demonstrates usefulness under the tested budget;
   it does not isolate GRAFT's growth, branching, or deduplication contributions.

The existing Golden and 140-case results establish measured behavior in their
documented scopes. High recovery or lower recorded CPU does not establish
historical novelty. The adaptive no-sweep run changes additional policies and
therefore cannot substitute for a matched component ablation.

**Search scope and remaining uncertainty.** This was a targeted literature and
source comparison across classical subgraph matching, atom-mapping searches,
symmetry-aware alternatives, partial-map extension, and search-state sharing,
including 2025–2026 work. Sources were reached through queries for fragment growth,
incremental mapping, branching, symmetry, compressed candidates, deduplication,
partial-map canonicalization, and decision/search graphs, then by following
references and inspecting author code. Primary sources support the comparisons;
review articles were used to find relevant methods.

The audit inspected selected algorithms and functions, not every implementation
path of every competitor. It did not establish the user's date of conception,
inspect proprietary mapper internals, reproduce all prior algorithms, or perform
an exhaustive bibliographic search. An exact equivalent could use different
terminology or appear in inaccessible or unlocated material. Availability limits
and pinned source revisions are recorded in the source manifest.

The current draft's five-reference bibliography is insufficient for a novelty
argument. Its related-work discussion should incorporate the closest precedents
above, especially Jaworski, SMSD/McSplit, Mann, Ali, Laffitte, and the relevant
SLAP/SynKit implementations. This audit supplies the basis for that revision;
it does not turn the current manuscript into a priority-certified claim.
