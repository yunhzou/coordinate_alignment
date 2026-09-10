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

## First successful shared-state scheduler

Run `adaptive_seed_route_20260910_4bbcer` (under the same benchmark root),
mapping job 457143, recovered best scores **5, 3, 8, 4** by the 400-operation
checkpoint. CPU seconds were **4.535, 0.756, 0.796, 0.236** on bosque12.
This is a no-cut policy, not a reproduction of the full cut/seed result set.
488 tests passed (job 457156). All 12 saved checkpoints passed graph validation
(457155). Full-H bounded symbolic membership checks (457177) recovered
**47/546, 28/48, 8/17, 8/8** historical best-score heavy-relation classes.
The first and third queries hit their 60-second evaluation budgets. These are
positive recovery counts, not minimum-score accuracy or chemical ground truth.

Tracing a missing Ni TS11 class found another order effect: its saved source
was `seed_02/cut_0075.pkl.gz`, cut (50,69). Keeping that order without the cut
recovered the same 3-event heavy pattern in an uncut representative (job 457187).
Its path starts with seed 106; the adaptive root had tried only 25, 68 and 114
after 1,600 operations. Strong preference for new regions was starving seeds
inside previously explored regions, despite retaining them in the queue.

The next bounded ablation will share alternative work between both novelty
classes, preferring new regions only to break service-count ties. Depth sharing
and one-completed-route yielding remain. No seed/cut/reference is supplied from
this diagnosis to the blind search. This changes priority, not represented
choices or native growth.
