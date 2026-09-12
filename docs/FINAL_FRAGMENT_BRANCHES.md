# Final fragment branch deduplication

A final branch is the unordered collection of matched fragment pairs `(R atom set, P atom set)`, under a fixed graph floor and matching tolerance. Atom order, island numbers, seed, and fragment growth order do not create additional final branches. The product atom sets are labelled endpoint atoms: different pairings remain different branches.

`AAMResult.branches` now uses this identity. Each returned branch retains all saved paths. Its `representative` and `hierarchy` describe only the first witness and must not be treated as the whole branch. The legacy analytical adapter explicitly uses `graph.literal_branches()` so that consuming one hierarchy at a time cannot discard alternatives.

For postprocessing, `aam.final_catalogue()` creates a flat union of mapping families. Each family retains its representative, required bond constraints, ordered symmetry actions, tolerance policy, and original path references. Independent actions commute and are normalized; overlapping actions retain their order. Identical family records are shared across branches. This performs no explicit mapping or automorphism-group enumeration. It is conservative family deduplication, not a proof that all remaining family relations are disjoint.

Growth history is absent from branch identity and is unnecessary for decoding the exported catalogue. Constraints and symmetry choices that affect the saved allowed mapping relation remain inside the branch. Keeping only the first representative would lose solutions.

## Intrinsic automorphisms

A given final fragment pair has a shared intrinsic automorphism definition. The optional `rebuild_all_symmetries()` labels product bonds by their responses to all matching tests against reactant bond weights in that pair. Exact automorphisms of this labelled graph preserve those tests. These groups depend on the pair and matching policy, not previous growth assignments.

A tolerance test by itself need not define a transitive equivalence relation. Moreover, two admissible mappings for the same pair can lie in different orbits of its intrinsic group. The tests include an explicit three-carbon example. Sharing the automorphism therefore does not justify deleting all but one representative.

The shared intrinsic groups are metadata in the default catalogue. They do not expand the saved mapping relation. Event thresholds remain decoder parameters; matching-response symmetry is not a promise of event invariance. Experimental broader reconstruction is isolated in `bench/experiments/final_fragment_dedup/rebuild_intrinsic.py`. It can introduce matching-admissible alternatives and was not adopted as the default optimization.

## API and CLI

```python
catalogue = aam.final_catalogue()
print(catalogue.branch_count, len(catalogue.families))
catalogue.rebuild_all_symmetries()  # optional shared pair metadata
record = catalogue.to_record()
restored = FinalBranchCatalogue.from_record(aam.problem, record)
path = restored.families[0].as_path(aam.problem)
# Pass path to extract_path_events or query_path.
```

Import `FinalBranchCatalogue` from `rxn_core.final_branches`. The flat archive fingerprint checks endpoint elements and WBO matrices. Reloading requires the endpoint problem but no original search graph. Endpoints must be balanced. Catalogue construction skips partial mappings and reports their count.

From an environment with this package and its native extensions installed:

```sh
python bench/final_fragment_branches.py --archives run/aam.pkl.gz --output final
python bench/final_fragment_branches.py --catalogue final/final_branches.json.gz \
  --problem input.json --output decoded --decode --max-events 8
```

The CLI has a five-minute watchdog and sampled memory limits. A complete result applies only to the requested event window, and every family must finish. Inspect `complete_window` and `unfinished`; a timeout is never a completeness certificate. Very large JSON catalogues can require substantial memory to reload; the CLI guard stops them if needed. Family count, branch count, explicit mapping count, and event-class count are different quantities.
