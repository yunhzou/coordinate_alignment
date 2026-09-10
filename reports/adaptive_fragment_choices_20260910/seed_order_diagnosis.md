# Why the fixed-order closure pilot missed case 77

Diagnostic run:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_gap_20260910_rkSgcC`

For each of the ten historical seed-index folders, take its first saved
8-event cut result. Replay its **same seed order with no cut**, retaining the
normal sequential matcher, explicit H, tolerance 1 and cap 100.

Eight of ten uncut runs still return an 8-event representative. The other two
return 13 and 18 (the latter capped). Thus the apparent benefit of those eight
cuts was confounded with cut-specific random ordering. This is a diagnostic
selected using historical success, not a blind benchmark or seed policy.
Full replay graphs, orders, successful fragment histories and timings are saved.
Job 457068 completed. The first diagnostic attempt (457066) stopped before
search on an incorrect archive field name; the corrected reader uses state
mappings. Failed attempts are not timing observations.

## Next isolated experiment (recorded before implementation)

Replace global seed-order restarts with lazy **next-fragment seed alternatives
at a shared sequential state**. All grows remain conditioned on existing
mapping, islands and boundary constraints. Merge exact equal states and retain
incoming correlated symmetry histories, using the existing search graph.

At each state, try the normal next seed first. Queue additional untried seeds
lazily. Prefer seeds outside the source regions already grown at that state;
seeds inside those regions remain available, not discarded. Complete ordinary
continuations before alternatives at the same discrepancy depth. This changes
exploration ordering only, not native growth or compressed candidate emission.

This first ablation uses **no early closures and no cuts**, so it isolates
fragment ordering. Fixed RNG 42; no supplied historical orders or mappings;
same four target pairs; operation/time budgets and external watchdogs unchanged.
Require saved complete outputs, quality/CPU measurements and graph validations.
Equal best scores are necessary but not sufficient for full family-coverage
equivalence. Production AAM remains untouched.
