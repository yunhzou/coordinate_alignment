# Three-seed bidirectional AAM ablation

Status: preparing a fresh, separately archived experiment. No accuracy or speedup claim yet.

Both the ten-seed baseline and this variation are **bidirectional**. This is
the frozen native-reuse engine `98b01b1`, not the later adaptive-policy candidate.
The earlier publication engine `bb0a7d7` has certified 99.35% Golden recovery;
the native-reuse baseline's saved evaluation certifies 1,836/1,851 (99.19%),
with five verified misses and ten unresolved cases.

Only search configuration change: `seed_count=10` to `seed_count=3` per cut.
The existing deterministic generator makes three orders exactly the prefix of
the ten-order sequence. Keep full sweep, branch cap 100, matching tolerance 1.0,
explicit H, root seed 42, both directions and eight search CPUs per call.
No original algorithm, native binary, search/profiling/evaluation adapter is changed.
Python branches/fragments, correlated symmetry and complete cut archives remain.

Test all 1,851 Golden reactions and all 140 XYZ/WBO holdout steps. Report them
separately; the holdout has no ground-truth labels. Golden recovery is the strict
heavy-atom mapping relation modulo endpoint chemical symmetry, including unmatched
atoms, among returned compressed families—not top-1 accuracy. The bidirectional
result is the union of the two searches; branches are not spliced.

Independent 300-second search and 240-second analysis watchdogs. Search and
evaluation occupy separate allocations. A shared locked work queue dynamically
assigns fresh tasks within a 1,024-CPU budget instead of leaving static batches
waiting behind a slow molecule. Existing attempts are never silently rerun.

Search CPU uses the unchanged profiler, including child processes and excluding
measured archive writing/loading. Evaluation, raw wall time, queue, failures,
and persistence are separate. Compare paired completed reactions to the saved
ten-seed timings; retain the full denominator and all unknowns for accuracy.

Full artifacts:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_three_seeds_bidirectional_20260910`
