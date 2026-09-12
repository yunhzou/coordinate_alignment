# Final fragment branch deduplication — 140 reactions

The median branch count falls from **596.5 to 397 per reaction**. Across the 140 reactions, the count falls from **237,645 to 124,641** (47.55%). These are unordered matched fragment-pair combinations, not decoded mappings or event classes.

| Quantity | Previous ordered branches | Unordered final branches |
|---|---:|---:|
| Total across 140 reactions | 237,645 | 124,641 |
| Median per reaction | 596.5 | 397 |
| Case 25 | 52,669 | 2,856 |
| Case 64 | 229 | 132 |
| Case 129 | 864 | 514 |

The change groups identical final fragment combinations while retaining their union of saved mapping alternatives. Atom order, island identifiers, seed, and growth order do not affect the branch key. Different product fragment sets and matching policies remain distinct.

There are still **236,653 unique flat family records** inside these branches. Only 992 old literal relations collapse to identical normalized family records. Therefore the 47.55% branch reduction must not be described as a decoding speedup. Further compression would require proving equality or containment of different mapping families; deleting all but the first witness would be incorrect.

## Inputs and cost

This experiment reused all 140 saved one-seed, branch-cap-2,000 baseline searches plus accepted results from the previous bounded fragment-competition experiment (up to 128 repair searches per reaction). It did not rerun the matcher. The old branch total includes baseline and accepted repairs, deduplicated within each reaction. The 99 repair graphs previously rejected as inconsistent remain excluded. Case 25 uses its successful retry archive.

The full catalogue construction, shared intrinsic symmetry computation, representative checks, and export took **91.54 seconds wall time**, **180.48 CPU seconds**, with **3 workers**. Peak combined sampled memory was **1,748.56 MiB** (1.71 GiB). Each case had a 300-second watchdog, a 3-GiB sampled worker limit, an 8-GiB combined limit, and a 6-GiB free/inactive memory reserve. All 140 jobs completed. These timing figures exclude later export validation and fresh pilot decoding.

A subsequent streaming export check validated all 140 catalogues and every family's serialization, required bonds, reference indices, and provenance counts in 165.84 seconds, using 61.22 MiB peak sampled memory. An earlier whole-file JSON validation attempt exceeded the 3-GiB worker guard; streaming resolved the validation memory issue. Full in-memory catalogue loading remains guarded by the CLI.

All **238,335 input paths** are retained as provenance references. No explicit mapping enumeration or automorphism-group enumeration was introduced. The family deduplication is conservative; remaining families can overlap semantically.

## Correctness checks

112 focused tests passed, covering search-graph and analytical compatibility, event decoding, order-independent fragment identity, prefix-lock preservation in flat families, noncommuting alternatives, distinct constraint variants, persistence, and the multiple-orbit counterexample.

Six cases were freshly decoded through the flat representation, without a pattern-count cap, in exactly the previously completed event windows:

| Case | Event classes | Result |
|---|---:|---|
| 6 | 5 | Complete; identical |
| 11 | 2 | Complete; identical |
| 59 | 11 | Complete; identical |
| 64 | 27 | Complete; identical |
| 101 | 2 | Complete; identical |
| 129 | 30 | Complete; identical |

Case 129 initially had four families exceeding the three-second pilot allowance. Retrying only those four with twenty seconds per family completed all four in 49.08 seconds wall time. No class changed. The updated result records retain the initial unfinished list and retry evidence.

An independent CLI check reloaded case 64 from the flat JSON catalogue plus endpoint input only, without its original search graph. All 27 event classes were recovered, with the requested window complete.

All 140 cases were counted and their family representatives checked. **Only the six pilot cases were freshly decoded**; this experiment does not claim a fresh exhaustive decoding of all 140. Completeness applies to each recorded event window, not arbitrarily high event counts.

## Why broader symmetry reconstruction is separate

A pair's intrinsic automorphism is independent of fragment growth order, but multiple admissible mappings for that same pair can lie in different intrinsic symmetry orbits. The test suite demonstrates this using three carbons and broad matching tolerance. The default preserves such alternatives inside one branch.

An additional experiment rebuilt broader intrinsic groups and aligned whole-product symmetries. This can add admissible mappings absent from the saved search, so it is an expansion rather than pure deduplication. It was not adopted: case 11 took about 122 seconds in that version, and case 129 left work unfinished at its soft time budget. In case 64, three additional classes were at 7, 8, and 8 events and did not recover the missing four-event SLAP alternatives. The experimental code is isolated from the default API.

## Evidence and reproduction

The previous working state was committed first as `312fff7` (`Checkpoint exact event decoding and bounded fragment competition experiment`). The new implementation is documented in `docs/FINAL_FRAGMENT_BRANCHES.md`; the portable entry point is `bench/final_fragment_branches.py`.

`summary.json` includes all 140 before/after counts, timing, pilot windows, and current source hashes. Local full artifacts are under `/Users/yunhengz/Documents/Codex/2026-09-11/aa/outputs/final_fragment_dedup/`: `lossless140`, `lossless_pilot`, `cli64`, and the exploratory `pilot` / `pilot_global` runs. Frozen implementation files in each experiment folder preserve the code used at measurement time. Subsequent persistence validation added endpoint fingerprints without altering the saved family content.

The benchmark runners in `bench/experiments/final_fragment_dedup/` preserve local dataset paths and depend on the existing saved archives. They are experiment records, not portable dataset installers. The large catalogues stay outside Git; this report and compact evidence are committed.
