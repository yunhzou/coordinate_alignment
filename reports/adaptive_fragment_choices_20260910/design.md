# Adaptive fragment-choice experiment

Baseline: `98b01b1`. Default AAM and its seed/cut protocol stay unchanged.

Replace external seed-times-bond-cut reruns with a bounded agenda of sequential
fragment decisions. One deterministic seed order guides the initial path.
Native growth retains compressed frontiers before extensions; a lazy close
operation finalizes such a frontier as a shorter fragment, under the same
locked mapping. The source WBO/topology is not globally cut. Subsequent
fragments still see the existing mappings, islands and boundary constraints.

Modules:

1. Native growth observer: retain eligible compressed frontiers without changing
   normal growth results, cap decisions or candidate order.
2. Reusable Python fragment-choice session: normal result plus lazy earlier
   closures, all as existing `FragmentPlacement` objects.
3. Experimental sequential scheduler: prioritize normal continuations, keep
   earlier closure choices in an agenda, and continue from their parent states.
   Merge only exact search states under the same remaining scheduling context;
   preserve path-correlated symmetry and full graph provenance.
4. Saved-output evaluation: compare full-H event scores and represented mapping
   patterns with the original uncut pass and saved sweep results. Never treat a
   representative score as an admissible compressed-family pruning bound.

Initial policy exploration uses explicit work/time budgets, records pending
choices at termination, and retains all completed results even when capped.
It is not exhaustive, and no speed or recovery guarantee is assumed. A useful
result requires a cost-versus-recovery comparison, not simply fewer calls.

Before changing production defaults: exact unobserved/observed growth equality,
valid partial closures and conditioned generators, cap/partial-composition
tests, full regression suite, and bounded real-case experiments. Keep all full
outputs and separate compute, scoring and persistence times. A 300-second
per-case watchdog and ten-minute Slurm limit apply.
