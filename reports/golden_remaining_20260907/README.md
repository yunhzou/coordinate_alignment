# Investigation of the remaining Golden mappings

## Decisive finding: the evaluator ignores assignment domains

The Golden evaluator reads each fragment's `target_generators`, but does not generally consume `symmetry_domains`. **These are not interchangeable.** A domain can contain assignments in different orbits of the conditioned automorphism group.

This produces false negatives even when the saved AAM result contains a valid reference-equivalent assignment. It is a result-consumption problem, not evidence that native AAM failed to generate the choice.

### Case 1817: two identical peroxide reactants

The saved representative takes both product oxygens from one peroxide molecule. The reference takes one oxygen from each of two peroxide molecules.

In the saved P→R path, the final oxygen has:

- Representative: P atom 4 → R atom 27.
- Assignment domain: R atoms **27, 28, 29**.
- Conditioned automorphism generators: fix R atom 27; they cannot move that representative to 28.

Choosing the already-saved alternative **P4 → R28** leaves every other assignment unchanged. The resulting explicit-atom mapping is injective, preserves every non-singleton fragment exactly, and independently matches the reference's chemical mapping certificate. No AAM rerun or atomwise reconstruction was needed.

The saved witness is in `singletons.json` and the cluster file `singletons/1817/result.json`.

### Case 1729

The same positive-only domain check recovered a reference-equivalent mapping from the original saved archive. Its certificate, original-orientation mapping, terminal, transitions and selected domain value are also saved. This was previously an unresolved verification case.

### Minimal reproducible example

Source: two separate O atoms. Target: two O–O pairs. After the first O is fixed, the next saved domain includes both its partner and atoms of the other pair. Its conditioned automorphism group cannot swap the partner with the other pair. A generators-only verifier therefore misses an allowed different-pair assignment.

`tests/test_golden_singleton_domains.py` reproduces this without any Golden-case-specific rules. The diagnostic query can check independently validated singleton domains; it is **not** a complete verifier for correlated multi-atom domains.

## Controlled search tests

All tests retain three seeds per cut, explicit H, tolerance 1.0, independent cut seed streams and ordinary single-edge sweep cuts. No reference-directed seeds or chemistry-specific rules were added.

| Diagnostic | Result |
|---|---|
| Opposite direction, cap 100, on all 24 reported misses | **16 certified recoveries**, 8 without a generator-only recovery |
| Same direction, cap 2,000, on the 22 capped misses | 13 complete searches without a generator-only recovery; 9 searches reached the five-minute watchdog |
| Original saved trailing singleton domains | Certified recoveries for **1817 and 1729** |

The opposite-direction recoveries are 576, 602, 603, 737, 775, 780, 833, 845, 846, 867, 964, 1230, 1380, 1400, 1739 and 1817. The two uncapped original misses, 1739 and 1817, both recover in the opposite direction.

For 1739, the top smaller-first representative preserves the azide's N–N bond-order roles; the reference exchanges the terminal-N roles. The opposite search demonstrates that the reference is reachable at tolerance 1.0. This does not prove the smaller-first archive lacks a domain realization.

## What can and cannot be concluded

- Smaller-first is useful, not universally best. Alternate direction has a much stronger observed effect than simply raising the cap in this diagnostic set.
- The reported “remaining 35” are not all demonstrated AAM search failures. At least two were verification failures on existing archives.
- The initial saved-family tests also used a generator-only target-orbit rejection filter. That filter omitted domains too; its negative classifications are **provisional, not proofs of absence**. The updated diagnostic includes domains in the rejection-only orbit overapproximation, but the frozen results are retained honestly rather than silently rewritten.
- The singleton probe is deliberately positive-only: it preserves all non-singleton assignments and explicit-H assignments, checks injectivity and element identity, then independently certifies the full heavy mapping. It is not a replacement for a general correlated-domain verifier. A negative or timed-out probe means unresolved, not absent.
- No core AAM growth, seed, symmetry-compression or branch-cap implementation was changed in this investigation. The previously published fixed-run score is not silently combined with alternate-direction or cap diagnostics.

## Required next correction

Make verification consume the **complete fragment-family contract**: assignment domains, correlated multi-atom constraints, conditioned automorphism actions, and injectivity across successive fragments. Preserve compression; do not replace this with independent atom masks or enumerate complete bijections. Then rescore the existing checkpoints before spending more CPU on search tuning.

## Artifacts

Cluster investigation:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_remaining_20260907`

The immutable initial diagnostic engine and manifest, all search cut checkpoints, process exit statuses, generator-only checks and independent positive witnesses are retained. Main Slurm array: 440281, at most 32 two-worker jobs, 32 GiB per job, five-minute search watchdog and nine-minute job limit. Saved-only rechecks use separate bounded jobs; stalled unstarted jobs were rescheduled, not completed searches.

Local report:

`/h/399/yunhengzou/coordinate_alignment/reports/golden_remaining_20260907`

`cases.csv` records all 81 diagnostic entries with explicit warnings about negative generator-only results. `singletons.json` contains the independent saved-domain witnesses. `summary.json` records code hashes and recovered case lists. **136 regression tests passed.**
