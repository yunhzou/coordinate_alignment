# AAM: fragment primitive, search graph, optional post-processing

Implemented 2026-09-04. This is a representation/API refactor, not a new growth
algorithm. Explicit hydrogens, seed generation, tolerances, cut policy, native
candidate compression, and subtree-cap behavior are preserved.

## The three boundaries

```text
match_fragment(source, target, seed, context, config)
    -> FragmentMatchResult: compressed FragmentPlacement objects + cap evidence

search_aam(problem, config)
    -> AAMResult: AAMSearchGraph + inputs/configuration/diagnostics

group_mechanisms(aam)                         optional balanced event policy
    -> MechanismResult
    -> compile_mechanism_families(grouped)    optional exact event families
    -> select_rp_mappings(...)               optional chirality/geometry selection
```

`compile_mapping_families(aam)` separately compiles complete raw structural
families without first classifying mechanisms. Partial assignments stay in
the graph; the complete-family compiler explicitly rejects incomplete inputs.

## What owns what

```text
AAMResult
  graph
    contexts / roots        cut graph, seed order, anchors, atom domains
    states                  cumulative assignment, island partition, deferred edges
    transitions             committed fragment decisions, or explicit joins
      placement             typed FragmentMatch view
        witness / blocks / automorph domains / conditioned exact generators
      preserved_bonds       excludes deferred/cut fragment bonds
      seed / step           actual growth attempt
    stops                   objective_met / stalled / capped / incompatible_anchors
```

A live `_Branch` is only frontier state plus a graph-node reference. Forks share
history; they do not copy fragment chains. Exact reconvergence creates a join.
Joins are admitted at the same continuation step in the same search context;
different cuts or seed runs are never connected just because mappings agree.
The transition's `step` records the pass/seed position. No-progress attempts
do not create self-loops.

Raw `AAMSearchMetrics.raw_result_count` and `retained_branch_count` count
context-qualified terminal states, before optional cross-run relation/group
deduplication. They are not mechanism counts or atom-permutation counts.

`graph.paths()` lazily projects recorded root-to-terminal histories. It does
not enumerate atom permutations. `aam.branches` groups identical matched
relations while retaining all discovering paths. This is an explicit, possibly
larger downstream projection; collecting workers and saving the graph do not
need to unfold paths.

Online state dedup uses normalized island partitions, literal assignments,
and deferred edges. Completed relation dedup additionally keeps ordered
matched-fragment/symmetry evidence. Source fragment atom sets alone are **not**
an equality key. Broader exact-family comparisons remain post-processing.

## Reusable conditional fragment matching

```python
from rxn_core import match_fragment, FragmentMatchConfig, FragmentMatchContext

result = match_fragment(
    source_graph, target_graph, seed=source_atom,
    context=FragmentMatchContext(locked_mapping=already_matched),
    config=FragmentMatchConfig(iso_tolerance=0.5, branch_limit=100),
)
for placement in result.matches:
    witness = dict(placement)
    source_fragment = placement.fragment
    preserved_edges = placement.preserved_bonds
    compressed_relation = placement.symmetry
if result.capped:
    print(result.branch_count, result.branch_limit)
```

Inputs may be `WeightedGraph` or prepared weighted NetworkX graphs. Orbit maps
can be supplied in the context for reuse; otherwise the primitive prepares
them. This still calls `growth.island.grow_island`, which dispatches to C++ at
the existing boundary. The C++ implementation was not rewritten.

The primitive returns a seeded saturated fragment, not a proof of a globally
largest common subgraph. It neither chooses multiple seeds nor balances
endpoint compositions. AAM and retro detection compose this same primitive
under their existing, distinct scheduling rules.

`match_weighted_subgraph(...)` remains a higher-level fixed-query operation.
Its `SubgraphSearchResult.matches` contains only validated complete query
placements; `.graph` also preserves unsuccessful/capped search evidence.

## Symmetry and one realization

Domains describe compressed possibilities; they are not independent per-atom
choices. `transition.placement.target_generators is None` means exact groups
have not been finalized. An empty tuple means the finalized group is trivial.

Balanced `search_aam` finalizes these groups before returning. Standalone
fragment graphs can be finalized explicitly with `finalize_graph_symmetry`.
Sparse original atom indices remain original indices; unused indices in a
generator image array are fixed.

```python
import random

path = next(aam.graph.paths())
sample = path.sample(random.Random(7), steps_per_fragment=2)
assignment = dict(sample.mapping)
hierarchy_in_sample_frame = sample.hierarchy
```

This is a nonuniform random walk, not uniform group sampling. It composes
whole generators and transports later decisions under earlier choices. It
does not independently shuffle atoms or fragments. The realization retains
the source path, generator choices, and target action.

Caps limit live canonical search candidates—not stored nodes, permutation
counts, or total histories. A capped sibling can coexist with a successful
path. Success of one path is not a claim of exhaustive search.

## Detection, augmentation, and assembly evidence

- Detection records initial seeded calls as graph chunks, including caps.
  Initial families retain every equivalent discovery and identify the chosen
  witness. They do not turn those discoveries into extra reactants.
- Residual matching retains its own target-availability context and paths.
  Augmentation copies remain competitive ownership bookkeeping, not evidence
  that an actual chemical side product exists.
- Candidates retain `FragmentDerivation` references. Symmetry-related
  occupations transport the hierarchy and conjugate generators together with
  the mapping, recording the target action without rewriting raw evidence.
- Progressive matching retains each actual selection, its detection result,
  and source/target local-to-original atom maps. It remains the existing greedy
  selection policy, not a newly exhaustive assembly method.
- Catalog coverage is an index projection; recommendations carry a detection
  row/candidate reference. Existing assembly/ranking chemistry policies were
  not replaced by this refactor.

Optional mechanism repair likewise records its selection action while leaving
raw search witnesses untouched. A mechanism-selected representative is not
automatically the original search witness; use the graph for raw matching.

## Persistence and offline replay

```python
from rxn_core import search_aam, aam_from_record, write_aam_bundle
import json

aam = search_aam(problem, config, workers=8, intermediate_dir="run/aam_search")
write_aam_bundle(aam, "run/raw_view")

# Later, without invoking the search or post-processing:
saved = aam_from_record(json.load(open("run/raw_view/aam.json")))
write_aam_bundle(saved, "run/replayed_view")
```

The search saves completed cut chunks before exact-group finalization and the
full reusable AAM result afterward. The CLI enables this persistence before
entering optional post-processing. `search.html` is self-contained: it shows
the actual DAG, allows choosing paths and stepping through fragment decisions,
and colors original-index R/P atoms and preserved bonds from recorded matches.
It is a raw-search inspector, not a replacement for the assembled retro viewer.

Schema versions: `rxn_core.aam/v1`, `rxn_core.aam_search_graph/v1`, and
`rxn_core.fragment_detection/v3`. Detection rows store a shared graph archive
and references rather than copying a graph into every candidate. Historical
candidate-only records can still be read as historical records; absent search
history is not reconstructed or presented as newly recorded evidence.

## Verification

- Existing native-versus-Python growth tests and all previous chemistry tests
  remain in the suite. T05 still covers all 58 explicit target atoms with its
  three building blocks, with no duplicate ownership.
- New tests exercise forks, joins, cross-context separation, online admission,
  caps/surviving siblings, sparse IDs, correlated sampling/transport, archive
  round-trips, and serial/worker agreement.
- Exhaustive four-carbon ring reference: all 8 automorphisms, no extra maps.
- Saved tetraphenyl and methyl-rich tetratbu benchmarks reproduce the baseline
  mechanism keys, representative mappings, and fragment histories exactly.
  See `bench/search_graph_regression.py` for the reproducible snapshot harness.

This does not establish exhaustive chemical retrosynthesis or remove the
existing seed/cap approximations. The goal is faithful, reusable storage of
the alternatives the search actually retains.
