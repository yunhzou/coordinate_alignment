# Three-seed bidirectional AAM ablation

## Result

**Three seeds retain 99.19% verified bidirectional Golden recovery while using
2.97× less search CPU on paired completed cases.** All 3,982 directional searches
and all 3,982 analysis processes completed with exit zero; none reached a hard
process watchdog. Soft verification budgets still leave five Golden cases unknown.

| Native-reuse AAM | Recovered / 1,851 | Coverage | Verified absent | Unknown |
| --- | ---: | ---: | ---: | ---: |
| Saved ten-seed baseline | 1,836 | 99.19% | 5 | 10 |
| Fresh three-seed variant | 1,836 | 99.19% | 10 | 5 |

The equal recovery counts are **not identical case sets**. Case 1314 is no longer
recovered with three seeds; cases 850 and 1568 are now unknown. Previous unknowns
1786, 1799 and 1823 are now recovered. The latter two involved incomplete baseline
searches/analysis; 1786 had unresolved baseline verification. These gains do not
mean fewer seeds explore more possibilities. No union with the ten-seed outputs
was used to obtain the three-seed percentage.

On the 140-step XYZ/WBO holdout, both settings return full explicit-atom mappings
for all 140 cases. Best representative event counts agree in 139/140. Case 123
changes from 19 events (ten seeds) to 25 (three seeds). The holdout has no annotated
ground truth; this is not an accuracy or complete-alternative-equivalence claim.

## Compute and SLAP comparison

On 1,837 Golden reactions where both AAM settings completed both searches,
ten seeds cost 176,323.67 CPU-seconds and three seeds 59,466.76 CPU-seconds:
**2.965× lower compute cost**. On all 140 holdout cases the reduction is **2.544×**
(8,516.42 versus 3,348.29 CPU-seconds).

For the same 1,807 Golden reactions completed by both AAM settings and the valid
SLAP sweep, the saved timings are:

| Search | CPU-seconds per reaction |
| --- | ---: |
| AAM ten seeds | 87.87 |
| AAM three seeds | 29.77 |
| SLAP + single-edge sweeps | 17.62 |

This narrows the recorded compute ratio from **4.99× to 1.69×** relative to
SLAP+sweep, whose standalone full-Golden coverage is 96.97%. Timing excludes
measured saving/loading and uses mutually completed cases; it is not a measured
end-to-end or equal-resource parallel latency ratio. AAM instrumentation includes
IPC/process setup; the SLAP workflow timer excludes worker startup. SLAP's
interrupted in-flight work is not included in this paired successful-work comparison.

Across **all** Golden cases, including the difficult cases missing from the old
ten-seed timing, the new search uses 70,570.73 CPU-seconds (19.603 CPU-hours).
Do not substitute the smaller paired total for that complete workload total.

Observed execution spans, measured separately from first actual worker start to
last worker finish: search **263.69 seconds**, analysis **139.25 seconds**. These
include persistence and dispatch, exclude initial Slurm queue, and are not a
combined end-to-end duration. Search ran on 29 working 32-CPU allocations (928 CPUs
maximum); three unused allocations stuck in startup were cancelled only after
all work completed. Analysis used 32 × 8 CPUs. Each AAM search retained eight CPUs.

## Protocol and artifacts

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
