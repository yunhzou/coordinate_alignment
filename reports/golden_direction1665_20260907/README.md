# Case 1665: search direction comparison

**P→R recovers the reference with three seeds, cap 2,000 and tolerance 1.0. R→P does not.** Core AAM is unchanged. All nine original input components and all explicit H are retained.

| Direction | Reference mapping | Scope |
|---|---|---|
| R→P | Not recovered | Fresh full 144-configuration sweep; 4,801 terminal candidates; zero cap stops |
| P→R | Recovered | Already present in the completed no-cut result; 4,001 terminal candidates in that cut; remaining sweep stopped intentionally |

The no-cut reverse witness is terminal 3, reached by transitions 0, 1, 2. Its P seed sequence is **56 → 34 → 20**: the intact 18-carbon core, then the two intact 10-carbon building blocks. It is also top-ranked **within this no-cut result**, using the same original R→P ranking function as the forward candidates. This is not a full-sweep top-1 claim.

Returned P→R pairs are inverted into original R→P indices before scoring. The resulting mapping has the same joint chemical certificate as the ground truth, not merely the same element composition or atom coverage. The 62 explicit-atom assignments are saved; the benchmark reference annotates 38 heavy-atom correspondences, not individual H identities.

## Viewers

- [Verified P→R witness, displayed in the usual R-left/P-right layout](P_to_R/viewer.html).
- [Fresh R→P result](R_to_P/viewer.html).

Both are standalone downloadable HTML. The same `r` label on each drawing identifies one mapped source atom. Clicking highlights the actual pair. Reverse fragment details explicitly label **P seeds and P cuts**, with mapped R atoms listed separately. Reference-derived regions remain distinct from actual AAM fragments.

## Interpretation

This establishes a directional difference for this case, not a general guarantee that reverse search solves all failures. Reversal changes which endpoint supplies seeds and swept edges: a single product versus the nine-component reactant-side graph. It also changes the full-source stopping objective. Both directions use the same existing AAM code; no guided seeds, supplied mappings, or chemical filtering are used.

## Timing and stopping

- R→P: **21.80 s** search plus saving, **7.55 s** reference evaluation, eight workers. This reproduces the previous run's candidate count, top terminal and failure.
- P→R: stopped after positive reference verification, **5:55 allocation elapsed**; 28 completed cut JSON archives preserved. All 18 completed cuts checked by the positive-only diagnostic contained a reference witness. Peak batch RSS approximately **19.3 GiB**. The batch cleanup finished at 6:02.
- The reverse sweep did **not** finish. Its elapsed time is not a completed-sweep runtime, first-witness latency, or a speed improvement. The selected no-cut archive records **49,222 cap stops**, despite containing the correct answer.
- Each run had an eight-CPU allocation, 64 GiB memory, a 570-second process-group watchdog and a cgroup OOM guard. No core search settings changed between directions except swapping endpoints; their sweep edge sets necessarily differ.

## Evidence and reproducibility

Full data and immutable search source snapshot:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_direction1665_20260907`

`comparison.json` records the scope of each result. `early_recovery.json` preserves the original-orientation certificate-verified witness; `partial_reference_check.json` records the checked reverse cut witnesses. All completed cut files remain on the cluster. No search was rerun to construct the viewers.

Launcher/preparation: `bench/golden_direction_probe.py init`. Scoring/viewers: `bench/golden_direction_probe.py compare --reverse-cut .../cut_00000.json`. The explicit cut argument marks a partial result; no synthetic completed AAM archive is constructed.
