# Root-cause probes before any AAM modification

Production baseline: `032c7d3`. Diagnostic-only checkpoint: `2f91344`. No files under `src/rxn_core` were changed during this investigation. All original benchmark labels and scores remain unchanged.

## Established: case 986

The original RDF itself swaps adjacent pyridine atom labels 15 and 16. Independent re-parsing of the RDF and the existing input preparation agree. This is not an atom-index conversion error.

The ordinary first-order result preserves the pyridine ring and has three changed heavy-atom bonds. The supplied reference has seven: the same three plus four changes caused by the adjacent-carbon swap. A **diagnostic-only** reversal of that swap gives a mapping chemically equivalent to the ordinary result under the existing certificate check. The benchmark reference was not edited.

The traced exclusion is inside fragment growth. With seed 26, the ordinary trajectory commits source edge 13–12 (WBO 1.5). The reference places its endpoints at target 5–3, which have no bond (WBO 0). That reference cannot survive inside this committed fragment at tolerance 1. A different, ring-preserving candidate can extend, so growth commits the shared extension rather than also retaining a branch that declines it. The reference-constrained diagnostic defers this edge and subsequently maps the adjacent-carbon pair as a separate fragment.

At cap 2000, both the ordinary Python trace and ordinary native run are uncapped and return **identical complete graphs**. Both miss the original reference. The constrained probe represents all 16 reference pairs without hitting the cap. Thus the representation can express this reference; the traced ordinary decomposition does not preserve it. This establishes an exclusion point in one trajectory, not a proof that every possible seed/cut combination must fail.

## Established: case 1285

The independently parsed original RDF also agrees with our prepared reference. This reaction's reference exchanges the oxalic-acid carbonyl oxygen attachments across the retained carbons; the ordinary results retain a larger connected oxygen/carbon fragment instead.

In planned P-to-R search space, the observed ordinary fragment contains source C=O edge 7–8 (WBO 2), while the reference maps it to target 13–10, a nonbond. The reference requires separating those atoms. The ordinary candidates give five or six changed heavy-atom bonds; the reference gives nine. These are changed heavy-atom pair counts, **not** the pipeline's full explicit-H event score.

The Python and native graphs again agree exactly and neither hits the cap in this probe. Reference-constrained matching represents all 14 reference pairs, whereas neither ordinary witness is reference-equivalent. The constrained run uses smaller fragments. This demonstrates the same conflict between a successful larger extension and a reference requiring a different decomposition; it does not establish that the reference is chemically invalid.

## Interpretation and limits

- These probes use the original first random atom ordering, tolerance 1, no sweep, and cap 2000. They explain particular missed trajectories; the preceding blind ten-seed/full-sweep campaign remains the separate coverage experiment.
- Reference constraints are implemented only in the diagnostic through the existing node-policy interface. They change available candidates, so their success is **feasibility evidence, not blind recovery**. No reference mapping is preloaded as a completed result.
- The shared-extension commit rule is present in pre-native commit `6ed3e28` and pre-memory-fix `c2e7efe`; it was not introduced by distance seeding or the recent memory correction. Exact Python/native equality here rules out a backend discrepancy in these probes, not all historical regressions.
- Symmetry compression inside a selected fragment does not make a bond-breaking permutation into an automorphism of that fragment. An alternative decomposition is required in these examples.
- Do not conclude that all remaining failures have this cause. Other cases are still unclassified, and prior verification timeouts remain unknown.
- Do not change core growth to enumerate arbitrary fragment splits based on these two examples. Such a change needs an explicit search-completeness objective and bounded-cost design, with preserved baseline and regression tests first.

## Evidence

`986.json` and `1285.json` retain the probe outcomes, graph equality checks, original-RDF checks, atom indices, WBO conflicts, and diagnostic counterfactual evidence.

Full graphs and event traces:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_cause986_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_cause986_audit_20260908
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_cause1285_audit_20260908
```

Reproduce with `bench/golden_cause_probe.py`; summarize without rerunning search using `bench/report_golden_cause.py`. Each completed individual probe took under 0.04 seconds; execution was guarded by a 290-second timeout. Graph saves and traces are retained. These are diagnostic measurements, not performance benchmarks.
