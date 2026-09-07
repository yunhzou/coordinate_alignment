# Corrected Golden reference-family coverage

**The original saved run reaches 1,731 / 1,760 (98.35%), not 100%.** No AAM search was rerun for this result.

| Measurement | Result |
|---|---:|
| Previously certified reference mappings, reused unchanged | 1,725 |
| Additional references certified by the corrected verifier | 6 |
| Corrected fixed-archive coverage | **1,731 / 1,760 (98.35%)** |
| Reference absent from saved completed families | 28 |
| Unresolved because the original search is incomplete | 1: case 1793 |
| Separate union including prior opposite-direction diagnostics | **1,746 / 1,760 (99.20%)** |

The six recovered cases are **44, 1351, 1522, 1729, 1786 and 1817**. The union row combines different searches and must not be labelled fixed-policy accuracy. The 91 incompletely reference-annotated Golden records are outside this denominator, as before. This measures reference mapping recovery, not atom coverage.

## The fix

`rxn_core.family_query.query_path` consumes the existing Python `SearchPath` and its typed fragment hierarchy. It models native assignment-pool permutations, exact correlated automorphism-group factors, fixed/locked assignments, the chronological transport of later decisions, global injectivity, element identity and saved fragment bond constraints.

`automorph_domains` are summaries, not free atomwise domains. They are not expanded independently. Missing exact groups raise an error rather than silently becoming trivial groups. No complete bijections or group elements are enumerated in production.

Every positive query returns a **full explicit-atom mapping including H**, a replayable list of domain/group actions and the number of checked fragment bonds. The choices are replayed independently of solver expressions and checked for injectivity, element identity and fragment bond support. Golden additionally verifies the resulting heavy mapping with an independent chemical-equivalence certificate. Endpoint symmetry used to normalize reference identity does not modify the physical witness.

The benchmark evaluator no longer skips an assignment domain merely because the stored automorphism generators are chemically exact. Its query deduplication retains full-H feasibility history.

## Exact rejection optimization

Many H variants share an impossible heavy-atom projection. The verifier first checks that projection and caches an **UNSAT** result using its projected domain/group history. This is a necessary-condition rejection, not approximate positive matching. A feasible projection always proceeds to the full explicit-H model. Tests include a projection that is feasible while its full-H realization is impossible.

This eliminated the completed-archive verification timeouts in the final run. The first, unoptimized full-domain pass is retained separately for provenance; it is not mixed into the final result.

## What remains

The remaining 28 negative cases are listed in `summary.json`. A negative result means absence from the **saved families**, not proof that another seed/direction/sweep can never recover it. The original search used three seeds per cut, cap 100, tolerance 1.0 and the smaller-first single-edge sweep policy. Those search settings were not changed here.

Case 1793 lacks a finalized full-search checkpoint because its original search reached the watchdog. Saved completed cuts were finalized and checked within a bounded budget, without finding a certificate. Its outcome remains unknown, not a mapping failure.

Top-1 ranking was not optimized over domain alternatives. The correction is to reference-family coverage, not a claim of improved top-1 ranking.

## Artifacts

Final rescoring and all full witnesses:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_domain_rescore_v2_20260907`

Initial full-domain pass:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_domain_rescore_20260907`

Original immutable AAM archives:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`

Final Slurm array: 440443; case 775 was rescheduled as 440478 because its job never left CONFIGURING. One CPU worker and 24 GiB per job; five-minute process watchdog, six-minute Slurm limit. No original or completed AAM calculation was restarted.

`cases.csv`, `witnesses.json`, `manifest.json` and `slurm_accounting.psv` retain case-level results, full action proofs, the executed source hashes and resource accounting. The implementation has **143 passing regression tests**, including domain alternatives, correlated fragments, locked slots, metadata-only automorphism domains, injectivity and explicit-H feasibility.
