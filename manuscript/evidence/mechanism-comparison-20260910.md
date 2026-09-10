# Mechanism comparison: what would have to change to reproduce this algorithm?

10 September 2026. This follows the [initial novelty audit](novelty-audit-20260910.md)
and addresses the author's request to compare actual mechanisms, rather than
listing shared ideas. Method naming remains provisional. The examined algorithm
is the frozen benchmark baseline `98b01b175eeed31f70d13e7cbf178b80bf07c9e0`.

The subsequent [review-led assessment](review-led-novelty-20260910.md) broadens
the literature coverage and adds Indigo, RDT, SynKit approximate growth and
AMLGAM comparisons. It records original-method access gaps explicitly.

**Conclusion.** Its growth and branch scheduler is mechanically distinct from
the inspected SLAP, OpenChemLib growth, FMCS, McSplit, and SynKit procedures.
The differences change decisions and retained states, rather than merely names
or programming language. A specific specialized AAM algorithm is therefore a
defensible description. The strongest candidate for originality is the coupled
rule below: commit a common fragment extension on existential support across
compressed placements, defer unsupported frontier edges, continue through
uniqueness, then branch over saturation results with whole-island absorption and
contextual continuation joins. No inspected predecessor implements that complete
rule. This is evidence of distinction from those predecessors, not a proof of
historical priority over every greedy common-subgraph algorithm.

**The precise transition being compared.** A local growth state contains:

- `F`: the common source fragment;
- `C`: its retained compressed placement candidates;
- `Q`: a heap of source frontier edges, ordered by decreasing WBO with atom-index
  tie breaking;
- `L, I`: the already committed mapping and its source-island partition;
- `D`: previously deferred edges; and a set of already processed frontier edges.

For the next eligible edge `(u, x)`, the implemented control flow is:

```text
Cnext = contextual_dedup(union of supported Extend(c, x; F, L, I) for c in C)
if Cnext exceeds the live candidate cap:
    report this growth subtree as capped
else if Cnext is nonempty:
    C = Cnext
    F = F plus x, absorbing x's entire locked island when applicable
    add the newly exposed source frontier edges to Q
else:
    retain F and C; record (u, x) in D; continue with Q

when Q is exhausted:
    deduplicate the saturated placements
    return their compressed records to the whole-mapping branch scheduler
```

There is no successful-extension branch that also retains the old fragment as
a competing local topology. Uniqueness is tested after the heap loop, not as its
termination condition. Local candidate children are created during extension;
only the *whole-mapping* branch split is delayed until saturation. A returned
placement may still represent multiple correlated assignments.

This is directly implemented in [Python growth, lines 248 onward](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/growth/island.py#L248)
and in the [native growth loop, lines 1527 onward](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/native/src/engine.cpp#L1527).
Both loops were inspected. This is a source comparison, not a new full-backend
parity test. A deferred edge is bookkeeping evidence; this operation does not
delete it from the WBO graph. Later extension checks all active source edges from
the new atom to the current fragment, not merely the edge popped from the heap.

**A distinguishing decision.** Suppose two retained placements `A` and `B` both
fit `F`. The next selected source extension `x` has support under `A` but none
under `B`.

| Procedure | What happens at this decision? |
| --- | --- |
| This growth loop | Commit `F + x`; only supported children of `A` remain. `B` is not retained on `F` as an alternative fragment topology. |
| Independent greedy growth of each placement | `A` can grow through `x`; `B` may instead keep its smaller fragment and try another frontier extension. |
| Search over alternative fragment shapes | Explore `F + x`, but potentially also shapes that omit `x`. |
| Single-placement greedy growth | A previous choice of `A` or `B` determines whether `x` succeeds; the unchosen placement is unavailable locally. |

This table is a logical consequence of different transition rules, not a new
molecular experiment or a claim that every member of a broad algorithm family
has identical behavior. It exposes a substantive tradeoff: this algorithm can
lose a placement that would have produced a useful *different* fragment. That
possibility distinguishes it from exhaustive common-subgraph enumeration.

**Comparison with specific implementations.** Each row identifies a concrete
replacement that would change the algorithm. These are code/pseudocode
comparisons; no competitor execution was added in this assessment.

| Comparator and inspected operation | Its actual decision rule | What must change to reproduce this algorithm? |
| --- | --- | --- |
| **OpenChemLib 2026.2.0:** [`mapFromRootAtoms`](https://github.com/Actelion/openchemlib/blob/0754c781c23f7570b13866e202836b5449282ad6/src/main/java/com/actelion/research/chem/reaction/mapping/SimilarityGraphBasedReactionMapper.java#L574) | Builds a neighbor-pair similarity matrix around a concrete mapped pair, chooses its best eligible pair, immediately writes mapping numbers, and queues that pair for further growth. The outer routine tries alternative root-pair sequences and keeps the best scored completed mapping. | Replace immediate atom-pair commitment with a collection of placements sharing one source fragment. Choose the next source frontier edge independently of a best target-pair similarity score, and defer whole-mapping branching until saturation. Its root-sequence search means it must not be described as globally having no alternatives. |
| **FMCS:** [documented subgraph enumeration procedure](https://www.rdkit.org/docs/source/rdkit.Chem.fmcs.fmcs.html) | A seed contains a query subgraph and excluded bonds. Growth enumerates nonempty subsets of available bonds; successful subgraph tests admit competing enlarged seeds. A size objective and bounds govern the search. | Replace branching over query shapes with the single accepted fragment trajectory above; retain and update its placement representations for subsequent AAM decisions. “Grow a query while checking target support” is already present, so that description alone cannot identify this algorithm's originality. This row concerns the documented Python FMCS procedure, not every current `rdFMCS` optimization. |
| **McSplit:** [Algorithm 1, especially lines 8–23](https://www.ijcai.org/proceedings/2017/0099.pdf) | Each recursive state has a concrete partial mapping and paired vertex classes. It tries each image of a chosen source vertex, refines the classes, recurses, and then explores leaving that source vertex unmatched; an incumbent bounds maximum-cardinality search. | Introduce an aggregate placement collection with a common fragment, remove the local omit-vertex branch, and defer the whole-mapping split. Its classes describe admissible future assignments under one partial map; they are not the same object as this algorithm's placements of the already grown fragment. |
| **Jaworski et al. 2019:** [isomorphic-mapping procedure](https://www.nature.com/articles/s41467-019-09440-2) | Uses a bounded isomorphism decision tree; extensions match all immediate neighbors, with environment orders descending from four to one. Unmatched regions trigger enumeration of bond-cut subsets. | Replace neighborhood-level recursive extension with shared source-edge-priority growth, its collective accept/defer decision, and saturation-level placement branching. These differences follow from the main article. Its inaccessible supplementary algorithm details remain an unresolved comparison limit. |
| **SLAP:** [`_solve_all_lap`](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py#L479), `_update_labels`, `_hash_current_labels` | Solves label-conditioned LAPs, creates refined graph-label pairs from proper solutions, pushes unseen label states onto a stack, and prunes against the best accumulated cost. Irreducible label classes terminate refinement. | Change the state from paired label partitions to a growing fragment with placements; replace LAP/refinement transitions by supported frontier extension and saturation decisions. Label-refinement completion is not fragment saturation. This is a different search mechanism even though both programs branch and deduplicate. |
| **SynKit v1.5.0:** [`solve_kernel` / `dfs`](https://github.com/TieuLongPhan/SynKit/blob/b57f2a39f1807a2682913183b22e1d8d61531bba/synkit/Chem/Reaction/Mapper/exact/branching.py#L391) | Selects the next uncertainty-kernel atom, assigns one eligible target position, accumulates cost, recurses, and undoes that assignment. Orbit ordering and a remaining-assignment bound prune alternatives. | Replace assignment-by-assignment kernel search with the common-fragment transition and saturated-placement scheduler. Calling both procedures symmetry-aware branching does not make these control flows equivalent. This assessment does not certify SynKit's exactness claims. |

Graph-pattern growth is another relevant structural precedent. In
[gSpan's Subprocedure 1](https://sites.cs.ucsb.edu/~xyan/papers/gSpan-short.pdf),
one graph pattern is extended through its occurrences; supported child patterns
are explored recursively, with minimum DFS codes preventing redundant pattern
generation. Its deduplicated object is a graph pattern, rather than a placement
state under this AAM scheduler. It reinforces the FMCS comparison: common-pattern
growth and support testing alone are insufficiently specific, while the
accept/defer policy and eventual mapping commitments differ.

**What the compression actually does, and what it does not establish.**
[`_SymCand`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/matcher/state.py#L49)
stores a witness, target pools for source blocks, fixed information, and
automorphism-domain metadata. An extension does not blindly test that witness.
[`_support_witness_for_value`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/matcher/support.py#L97)
constructs allowed target domains for affected source atoms and searches for an
injective supporting assignment within their blocks. When the new assignment
depends on a block rearrangement, [`_children_from_context_group`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/matcher/extend.py#L513)
refines/fixes that support instead of treating the new choice as independent.
The support search is bounded and returns a supporting witness; this inspection
does not prove complete retention of every possible support under every input.

There is a useful conditional equivalence here. Let `E(F,L)` mean *all* concrete
embeddings of fragment `F` respecting locks `L`, with fixed graph constraints.
If an extension operator is complete, then for an ordinary free-atom extension:

```text
union(Extend(m, x) for m in E(F,L)) = E(F union {x}, L)
```

An embedding on the right restricts to an embedding of `F`; every valid extension
on the left is an embedding of the enlarged fragment. Thus maintaining candidates
incrementally can reproduce recomputing matches from scratch. This is an
idealized mathematical observation, not an asserted correctness theorem for the
bounded compressed implementation. It separates two claims: incremental reuse
is an implementation/representation contribution; the aggregate greedy
accept/defer and branching policy defines the specialized search algorithm.

**The deduplication relations are also concretely different.**

| Stage | Object compared and implemented equality |
| --- | --- |
| Local growing placements | [`_dedup_sym_cands`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/matcher/dedupe.py#L208) compares canonical certificates of the full, configured WBO-colored product graph with candidate roles/pools and individually marked locked images. Within matching certificates it also checks boundary evidence. Source identities are carried in the roles. This is not equality of occupied atom sets. |
| Whole-mapping continuations | [`_Branch.state_key`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/alignment/branch.py#L102) uses exact atom pairs, normalized source-island partition and deferred edges. Admission compares them at the same pass and seed position. Equal states share continuation nodes and retain incoming histories. |
| Completed families | The compiled [`AnalyticalMappingFamily`](https://github.com/yunhzou/coordinate_alignment/blob/98b01b175eeed31f70d13e7cbf178b80bf07c9e0/src/rxn_core/alignment/index_chirality.py#L888) uses an isomorphism representative and automorphism generators for joint mapping actions; completed-family equality/containment is handled separately from live-state joins. |

For example, the same atom-pair mapping can have island partition
`{{a,b},{c}}` or `{{a},{b,c}}`. Touching locked atom `b` during a later growth
attempt can absorb different source atoms in these two states. Equality of atom
pairs alone therefore cannot justify merging their continuations. This is a
symbolic explanation of the merge rule, not an added molecular benchmark.

SynKit's [`canonical_partial_mapping_key`](https://github.com/TieuLongPhan/SynKit/blob/b57f2a39f1807a2682913183b22e1d8d61531bba/synkit/Chem/Reaction/Mapper/exact/enumerate.py#L849)
instead canonicalizes supplied atom pairs under reactant and product
automorphisms. It has no island-partition, deferred-edge, or scheduler-position
argument. That function cannot directly replace this continuation key. This
does not imply its key is unsuitable for *its own* search problem.

Canonical labeling and memoization remain established techniques. The potential
contribution is the particular state relation required by fragment absorption
and continuation, and its efficient implementation. Exactness of a graph
certificate is relative to the configured coloring; it does not independently
prove preservation of every continuous-WBO tolerance decision or every bounded
search outcome.

**Existing trace evidence for the actual mechanism.** In the saved
[case 1 trace](growth_trace_1.json), event 17 reduces the placement count from
three to one at fragment size nine. Events 19 and 21 still commit extensions,
increasing the fragment to ten and eleven atoms while retaining one compressed
candidate. One candidate need not mean one concrete bijection. This verifies
that local candidate-count uniqueness does not terminate this recorded growth.
In [case 64](growth_trace.json), consecutive commits 25, 27 and 29 change the
compressed counts from 42 to 18 to 12 to 2 while growing the same source
fragment. Count changes alone do not separate feasibility filtering from
canonical deduplication. These are the existing illustrations, not new runs or
proof of superiority to a comparator.

**What this supports in the paper.** Present the transition and state relations
above as the algorithmic definition. The meaningful contribution claim is a
specialized search over shared fragments and retained placements, with specified
conditions for committing, branching, absorbing and joining. A claim to have
invented an elementary subgraph-matching operation is not supported. Nor does a
control-flow difference alone establish that the difference is useful; matched
ablation or trace comparisons should address that separate question.

The [additional source manifest](mechanism-source-manifest.json) records the
newly cached source revisions; the [earlier manifest](novelty-source-manifest.json)
covers the original Python and competitor snapshots. The paper's compiled PDF
and experimental code were not changed by this assessment.
