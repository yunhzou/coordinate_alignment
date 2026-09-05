# AAM search graph: design and coordinated refactor plan

Status: proposed design; **no implementation authorized or included here**.
Recorded 2026-09-04 after discussion with the user. The next step is to finish
reviewing the contracts and then implement the agreed changes together.

## 1. Agreed direction

Core AAM should expose a mechanism-independent search graph. A branch is a
path of fragment-matching decisions. Shared prefixes and exact reconvergence
are explicit. Symmetry belongs to a fragment decision's compressed placement
representation. Mechanism grouping is optional downstream processing.

```text
                       fragment decision A
                   +------------------------> state 2
                   |
root --> state 1 ---+
                   |
                   +------------------------> state 3
                       fragment decision B
```

A selected path supports fragment-by-fragment animation directly from the
recorded matching. A tree display may unfold the graph visually; storage must
retain shared nodes and references instead of copying every path.

This refactor changes representation and responsibility boundaries. It must
preserve the existing growth rules, WBO tolerances, explicit hydrogens, seed
selection, symmetry compression, and branch-cap behavior. It must not quietly
introduce stronger search pruning, exhaustive bijection expansion, or new
chemistry assumptions.

## 2. What the current code actually does

| Layer | Current representation and behavior | Required change |
|---|---|---|
| `matcher/state.py`, `native/src/engine.cpp` | Compressed candidate with witness, pools, fixed assignments and automorph domains | Preserve the candidate semantics and native execution boundary |
| `growth/island.py`, `growth/native.py` | `grow_island` returns `_IsoResult` objects with fragment and symmetry state | Record each committed result as a graph transition |
| `alignment/branch.py` | `_Branch` owns cumulative mapping, islands, deferred edges and copied `symmetry_paths`; `find_islands` owns a list of live branches | Separate immutable search states, fragment transitions, and the live frontier |
| `alignment/sweep.py` | Runs seeds/cuts, repairs/scores mappings, creates mechanism pools, merges mapping-plus-hierarchy records | Separate search collection from mapping selection and mechanism classification |
| `aam.py` | Returns `AAMResult.mechanisms`; attaches completed fragment generators after mechanism pooling | Return a graph before mechanism classification; finalize symmetry independently |
| `analytical.py` | Compiles and compares families within pre-existing mechanism groups | Accept graph-derived families independently of mechanism groups |
| `alignment/post_aam.py`, `domain.py` | Typed histories, full bijections and mechanism-owned branches | Define graph contracts and explicit downstream result types |

Current live-state merging uses literal cumulative mapping/island/deferred-edge
equality. It preserves alternative histories. Completed seed/cut deduplication
uses mapping plus hierarchy equality. Later analytical compilation performs
broader family equality/containment checks. These are distinct operations and
must remain separately named.

Some existing comments describe stronger merging than the executable code;
implementation and regression evidence take precedence when defining parity.

## 3. Proposed minimal object model

Names below are design names, not existing APIs or finalized signatures.

```text
AAMSearchResult
    problem                 source and target atom spaces / graph data
    graph: AAMSearchGraph
    diagnostics             search policy, counts, timing, limits reached

AAMSearchGraph
    contexts                cut graph, seed order, anchors, stopping rule
    roots                   context -> initial state reference
    states                  state ID -> SearchState
    transitions             transition ID -> FragmentTransition
    stops                   StopRecord entries

SearchState
    id, context reference, continuation position
    representative          injective source -> target assignment
    island partition        cumulative ownership / connected-island state
    deferred edges

FragmentTransition
    id, source state ID, destination state ID
    fragment atoms and preserved-edge evidence
    placement               representative + compressed symmetry state
    island/deferred-edge changes
    provenance              growth seed, local attempt / result identity

StopRecord
    state / attempt reference
    reason                  objective met, stalled, cap, configured stop
    limit details           which limit, observed count, configured bound
```

The state supports continuing search. The incoming transitions retain how it
was reached. The transition's placement representation owns local blocks,
automorph domains and the conditioned exact-group reference when finalized.
There is no separate unscoped "symmetry" sibling beside the graph or fragment.

Source/target graph storage is shared by reference, not copied per state.
Mappings may be cached on active states for speed; replay must not depend on
an arbitrary incoming parent when a state has multiple parents.

A representative is a concrete valid assignment used by the search. It does
not replace the compressed relation. Neither representative equality nor a
union of possible atom destinations establishes family equality.

## 4. State sharing and path semantics

1. Each seed/cut run has an explicit context. An equal mapping under a
   different cut graph, remaining seed schedule, anchor set, tolerance or
   stopping objective is not automatically an interchangeable search state.
2. A transition records one committed compressed fragment result. Equivalent
   permutations inside that result do not become separate graph branches.
3. Multiple distinct surviving results create outgoing transitions. Their
   descendants share the incoming prefix by reference.
4. Exact live-state equality at a compatible continuation position permits
   reconvergence. Keep every incoming transition and its provenance.
5. Current merging excludes histories from its live-state key because search
   continuation reads cumulative state. Audit that assumption before encoding
   it as an invariant; if a continuation reads path-specific information, that
   information must be part of the state or those contexts must remain separate.
6. Do not merge whole live states solely by endpoint automorphism or mechanism
   key. Broader mapping-family deduplication remains optional downstream work.
7. No-progress seed attempts are not self-loop matching transitions. Record
   them as attempt/stop metadata if needed. Keep continuation position explicit
   so progression and termination are unambiguous.
8. Initial implementation should share prefixes/reconvergence within runs and
   collect independent runs in one graph. Cross-run state interning requires a
   separate equality argument and is not needed for the first implementation.

Merging two equivalent continuation states must not invent new histories by
allowing an outgoing decision that is invalid for one incoming context. Tests
must exercise reconvergence with different fragment histories, not just equal
representative mappings.

## 5. Symmetry and preservation contract

There are two guarantees:

- Search coverage: alternatives explored under the chosen seeds, cuts and caps.
- Representation fidelity: all retained alternatives remain represented, with
  their correlations, through graph storage and downstream conversion.

The second is mandatory. The first is bounded and must be reported honestly.

During growth, automorph domains summarize variation after certificate-based
merging. They are not independent assignment permissions or a complete group
description. The existing full-AAM stage derives exact generators using the
target's colored graph, candidate relation, and locked prefix. Move this stage
outside mechanism pooling and attach/cache results by their actual transition
context. An explicitly trivial group and a group not yet finalized are
different states of knowledge.

Later consumers must apply correlated transformations in context. If an earlier
fragment is moved to another equivalent placement, dependent mappings and
group constraints must be transported consistently. Taking Cartesian products
of per-atom destination pools is invalid.

Keep the compressed certificate/domain evidence from native growth. Audit
generator reconstruction against it: reconstructed families must neither lose
retained assignments nor introduce unsupported ones. Do not treat agreement
between two implementations of the same rule as the only correctness evidence.

Encounter counts and merged-candidate multiplicities are provenance counters.
They are not exact mapping cardinalities. Any exact cardinality exposed later
must come from the represented relation/group calculation.

## 6. Balanced, partial and augmented matching

The core graph uses injective atom assignments with explicit source/target
domains. A completed balanced mapping is the special case that covers both
endpoints. Completion of a requested core is a separate objective. Search
termination, objective coverage, and search exhaustiveness are separate fields.

Do not automatically remove current endpoint-composition restrictions while
refactoring data structures. Migrate supported entry points explicitly and
test each existing contract. A generalized container alone does not implement
a new matching algorithm or an incomplete-reaction mechanism definition.

Augmented/residual graphs need explicit atom identity maps:

- source atoms retain original endpoint identities;
- target atoms identify original target versus augmentation-copy ownership;
- subgraph-local/native indices convert through stored index maps;
- repeated precursor copies have distinct instance identities.

Record source/target context changes between augmentation rounds. They must not
appear as ordinary fragment decisions in an unchanged graph. Artificial copies
are competitors; their presence is not evidence of actual chemical side products.

Store preserved fragment bonds separately from source bonds deferred/cut during
matching. A viewer must not reconnect separated pieces just because the original
source has a bond between two retained atoms.

## 7. Optional post-processing boundaries

```text
core search -> AAMSearchResult / AAMSearchGraph
    -> exact mapping-family compilation and optional deduplication
    -> optional mechanism classification/grouping
    -> optional chemical, chirality or geometry selection
```

The scheduling layer can still request seeds and cut sweeps. Graph events and
raw representatives must be saved before repair or scoring changes a selection.
Downstream selections reference their source paths/families and retain any
mapping transformation; they never overwrite the original search evidence.

Mechanism classification requires an explicit policy for incomplete endpoints.
The current `classify_bonds` convention for unmapped bonds and full mechanism
signature assumptions cannot silently define that policy. Balanced mechanism
outputs should reproduce their current results through the downstream stage.

One compiled family may need splitting by event class if its allowed mappings
do not all share a mechanism key. Grouping cannot classify the representative
alone unless event invariance over that family is established.

The current `_score_branch_mapping` also applies a coverage acceptance check
and optional symmetry repair. Audit these separately: preserve all raw search
outcomes in the graph, then reproduce existing selection behavior explicitly
where requested. Do not confuse exposing more raw outcomes with changed growth.

## 8. Consumer migration inventory

| Consumer | Required migration / risk |
|---|---|
| `subgraph.py`, `alignment/api.py` | Replace first-history `symmetry_fragments` reads with explicit graph/path selections; preserve all retained alternatives in full results |
| `core_aam.py`, `_core_mapping_variants` | Derive core assignments from all relevant paths with correlations intact; keep explicit core enumeration a requested downstream operation |
| `fragment_matching/detection.py` | Its initial stage calls growth directly; record those calls as graph transitions and associate each retained family with provenance |
| `fragment_matching/augmentation.py` | Carry nested matching contexts and all required histories through augmented projection |
| `fragment_matching/parallel.py` | Transport graph records or references from workers without enumerating paths |
| `fragment_matching/progressive.py` | Current result retains selected mapping/fragments only; preserve graph references for each residual matching round and record the greedy choice as a selection |
| `fragment_matching/symmetry.py` | Transform mapping and its graph/hierarchy reference together; avoid a transformed witness paired with stale symmetry/index context |
| `retrosynthesis/catalog_index.py`, `compressed_coverage.py`, `ranking.py` | Consume correlated occupations; keep coverage signatures as query/index projections, not replacements for matching evidence |
| `analytical.py`, `rp.py`, `ts.py` | Preserve graph provenance across family compilation and representative selection |
| `artifacts.py`, fragment serialization | Persist graph and selection references with explicit schema/version semantics |
| `tools/build_retro_db_viewer.py` and other viewer consumers | Render recorded paths and original atom identities; add tree/path animation without reconstructing invented search history |
| `__init__.py`, CLI and tests | Update public contracts and every in-repo caller together |

This table identifies migration work, not proof that each listed risk currently
causes a wrong chemistry result. In particular, inspect transformations and
first-history readers with targeted cases before concluding which lose evidence.

## 9. Native boundary, parallel execution and storage

- Keep `grow_island` as the Python entry to the native kernel. Build the
  fragment-level search graph around returned `_IsoResult` values. No Python
  callback for every C++ candidate expansion is required for this design.
- Fragment-level animation shows committed matches/forks, not an invented view
  of all inner candidate operations. Detailed inner-kernel tracing is a separate
  instrumentation task; the current events path disables native growth.
- Use compact typed records and indexed adjacency for storage. This proposal
  does not require NetworkX to own the search history.
- Workers produce local immutable graph chunks with context-qualified IDs.
  The parent joins chunks deterministically and builds references. It should
  not rerun matching or eagerly enumerate root-to-terminal paths during merge.
- Save completed chunks/intermediate search results before expensive optional
  processing. Native handles/caches remain runtime-only; persisted graph data
  must be sufficient for offline path replay and family interpretation.
- Version the format explicitly. Existing flat historical results cannot reveal
  unrecorded forks or timings. A deliberate legacy reader may describe its known
  terminal histories, but must never fabricate a complete graph or silently
  fall back to old first-history semantics.
- A live branch cap limits the active search frontier under existing policy,
  not the number of recorded nodes, paths or symmetry permutations. Record a
  capped attempt with its scope/count; an unaffected sibling can still finish.
- Graph recording has memory cost even with shared prefixes. Benchmark node,
  edge and payload growth; avoid per-path history copies and duplicated endpoint
  arrays. Any future retention limit must be explicit, never silent truncation.

## 10. Verification before one coordinated implementation lands

Capture current intermediate results first. Use them as reusable comparison
inputs rather than repeatedly running inventory scans while adjusting adapters.

Required checks:

1. One chain, a real fork, and reconverging identical states retain precisely
   their valid paths. Different continuation contexts do not cross-connect.
2. Symmetric hydrogens and correlated whole-fragment placements stay compressed;
   sampling/replay gives valid mappings without arbitrary atom mixing.
3. Tiny exhaustive reference cases verify compression/group reconstruction
   against actual allowed assignments; capped seed coverage is assessed separately.
4. Repeated seeds retain provenance and do not multiply shared prefix payloads.
   Same mapping with distinct histories remains distinguishable until justified
   family deduplication.
5. Explicit H, partial/core mappings, unequal composition in supported APIs,
   augmentation ownership, and sparse/native atom indices round-trip correctly.
6. Cap at growth versus cap at combined live frontier produce distinct records;
   unaffected siblings survive. No stalled/capped state is labeled a full match.
7. Native/Python parity checks compare terminal relations, fragment information,
   and diagnostics, allowing execution-specific timing and graph IDs to differ.
8. Existing balanced mechanism results are reproduced by the explicit downstream
   pipeline. Selected/repaired representatives reference their raw provenance.
9. Multiworker and serial runs agree semantically. Serialization round-trip
   supports offline fragment-path animation without any AAM rerun.
10. Reuse the long-tail explicit-H case and T05 saved fixtures; compare terminal
    coverage, symmetry relations, cap behavior, wall time, memory, graph size,
    worker transfer size and parent merge time on the same resources.

No new arbitrary numeric performance target is set by this note. Measure the
baseline first, and investigate material regressions before a full bank scan.

## 11. Decisions to finish before implementation

The direction above is agreed; these details still require analysis, not guesswork:

- Final continuation-state key: schedule/pass position and all future-relevant
  data, plus a monotonic progression rule establishing an acyclic history.
- Where graph collection sits for direct fragment growth versus whole-molecule
  scheduling, and how nested augmentation contexts link without duplicated work.
- The minimal persisted symmetry payload and exact-group finalization boundary,
  including transportation of dependent fragments under prefix alternatives.
- The public break from mechanism-owned `AAMResult` to graph-owned results and
  the corresponding CLI/artifact migration.
- Which existing callers legitimately request one selected path and which must
  receive the full relation. Make these explicit choices, not fallback behavior.
- How mechanism/event policy handles incomplete endpoints; until specified,
  retain partial matching output without asserting a full mechanism.

Implementation should then be one coherent change spanning graph contracts,
scheduler, output conversion, consumers, serialization, and regressions. Do not
ship a wrapper around the current branch list while leaving path loss in adapters.
This document itself changes no search code, defaults, artifacts or viewers.

## 12. Online deduplication cleanup to include

Follow-up requirement: suppress redundant continuation online at each growth /
fragment-decision step. Keep necessary exact symmetry checks; simplify the
branch bookkeeping around them. This is planned work, not an implemented speedup.

Confirmed current behavior:

- `_extend_sym_cands` deduplicates each successful atom-extension result;
  the native implementation does likewise. Fragment saturation performs another
  context-specific deduplication, including single-atom/no-extension cases.
- `_admit_subtree` already admits unique live-state keys. The subsequent pass
  over `new_branches` repeats the same state-key deduplication. Verify that no
  admitted key can change, then consolidate these into one admission operation.
- `_progress_key` uses literal island numbers, so identical partitions with
  different numbering can compare unequal. `_island_partition` exists but is
  unused. Define a canonical partition key and normalize all dependent label
  bookkeeping consistently, including future island allocation, before merging.
- Forks copy history lists; `merge_exact_paths` reconstructs frozen keys for
  existing histories. Graph prefix sharing and stable transition/path references
  should replace that repeated copying and hashing.
- `_run_cut_work` accumulates completed records before `_cs_wrun` adds them to
  its local deduplicated pool. Insert completed outcomes incrementally instead,
  retaining encounter/provenance information and original diagnostics.

Desired admission rule:

```text
produce a compressed matched-fragment decision
    -> identify resulting state within compatible continuation context
    -> existing state: attach history/transition, reuse scheduled continuation
    -> new state: register and schedule it
```

State reuse across different seed schedules/cuts is not justified by a shared
fragment object alone. Distinguish deduplicating stored output from skipping
future computation. Reuse continuation only when its full context agrees.

These changes should reduce redundant bookkeeping and repeated compatible
continuations. Quantify the actual gain with the saved benchmark inputs; do not
claim a speedup before measuring it. Do not replace exact candidate-equivalence
proofs with independent atom-orbit labels to make the deduplicator cheaper.

## 13. Concrete automorphism representation and sampling

A generator is a full target permutation represented as an image array:
`g[p]` is the destination of target atom `p`. A single generator can move many
atoms simultaneously. Its cycles describe simultaneous moves, not independently
selectable swaps. Products of the stored generators describe the generated group.

Example: one allowed generator swaps target chains `(10,11,12)` and `(20,21,22)`:

```text
10 <-> 20, 11 <-> 21, 12 <-> 22     all swaps in one generator
```

Applying it to a representative placement of R atoms `(1,2,3)` at `(10,11,12)`
moves the whole fragment to `(20,21,22)`. It does not authorize `(10,21,12)`.

For a complete compiled family in the current code:

```text
m = family.representative_mapping
g = an element generated by family.target_generators
alternative = {r: g[p] for r, p in m.items()}
family.contains(alternative) verifies membership
```

Use the representative and generators from the same family. The compiled
representative can differ from the original search witness. The current
`AnalyticalMappingFamily` is a downstream compiled relation; its availability
must not be mistaken for every raw growth result already having a standalone
whole-branch group.

Choosing a generator or composing a short random sequence gives a valid member
of that group, but is not guaranteed to sample members uniformly. No group-wide
enumeration is needed for such a witness. Uniform sampling is a separate API
contract and is not required for the intended viewer/recommendation use.

At fragment-history level, exact generators are conditioned on the locked
prefix. Sampling earlier choices requires consistent transportation of later
assignments and constraints. Do not shuffle each fragment independently, union
all atom domains, or assume every raw fragment generator is independently
applicable to a completed branch. The proposed graph-level sampler must own this
dependency handling and return one valid path realization with provenance.

For unequal endpoint compositions the same target action yields an injective
partial mapping, not a full bijection. Structural-family membership does not
by itself imply chirality acceptance; that belongs to the requested selection
policy.

## 14. Reusable conditional fragment-matching function

Formalize the existing fragment-growth operation as an independently usable
package function. Both whole-molecule AAM and retrosynthesis detection should
call this contract, rather than reaching through mechanism/result internals.
This is an API proposal around existing matching semantics, not a new solver.

```text
match_fragment(source, target, *, seed, context, config)
    -> FragmentMatchResult
```

- `source`, `target`: prepared weighted graphs with stable atom identity maps.
  Preparing or caching these graphs must be reusable across calls.
- `seed`: the source atom from which this fragment grows. Selection of multiple
  seeds belongs to the calling search policy, rather than being hidden here.
- `context`: locked source-target assignments, island partition, deferred-edge
  evidence, and the graph/availability context needed to continue matching.
- `config`: explicit matching rules (atom compatibility, graph floor, bond
  tolerance, permitted mapped-seed behavior) and fragment-search limits.

```text
FragmentMatchResult
    matches[]              distinct retained compressed fragment placements
        source fragment and preserved-edge evidence
        representative placement
        local blocks / automorph relation and its conditioning context
        deferred edges / island-merge information
    diagnostics            outcome, cap scope/count, elapsed work
```

Each match is the R-P matched object discussed with the user. Its representation
includes correlated alternatives; the top-level list does not enumerate every
equivalent atom assignment. Results can be attached directly to search-graph
transitions. The match references the context under which it is valid.

Responsibilities:

1. Reuse the existing incremental extension, compatibility checks and online
   symmetry deduplication through the native-enabled growth boundary.
2. Return every surviving compressed alternative within the configured search,
   with original atom indices and explicit cap diagnostics.
3. Retain sufficient conditioned symmetry evidence for exact realization. The
   API must distinguish domains awaiting finalization from available exact
   generators; finalization ownership remains the decision listed in section 11.
4. Report a saturated fragment under this seeded search. Do not claim a
   globally largest common fragment or exhaustive matching over all seeds.

The consumers then own different compositions of this same operation:

```text
AAM scheduler:
    choose seed -> match_fragment -> admit graph transitions -> continue

Retro fragment detection:
    choose initial seeds -> match_fragment -> construct augmentation context
    -> match residual components under that context -> retain target ownership
```

The fixed-query subgraph API remains a distinct higher-level contract: it must
cover the requested query and verify its required edges. It may compose multiple
fragment-growth calls. A single seeded-growth result is not automatically a
complete subgraph embedding. Augmentation likewise remains detection workflow
logic; mechanism scoring, catalog ranking and assembly are outside the matcher.

Existing `grow_island`, `match_weighted_subgraph`, and `detect_fragments` represent
these different granularities. The coordinated refactor should clarify their
public contracts and route common work through the fragment primitive, without
duplicating search implementations or introducing one catch-all pipeline API.
