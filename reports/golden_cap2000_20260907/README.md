# Golden: cap 2,000 diagnostic

Controlled search comparison, not a replacement for the fixed-configuration paper benchmark. Explicit H, tolerance 1.0, ordinary no-cut plus single-edge sweep cuts, eight cut workers. No reference-directed seed or cut choices. Core AAM and default settings were not changed.

The ten three-seed cases were selected from the previous investigation: five empty outputs, four capped misses tested at cap 1,000, and unresolved blind case 1740. Cases 1665 and 1740 additionally repeat their previous 100-seed configuration. This is a targeted diagnostic, not a representative accuracy sample.

## Completed three-seed tests

| Case | Reference recovered? | Still capped? | Search + archive time |
|---|---|---|---|
| 11 | No | Yes | 11.27 s |
| 602 | No | No | 2.11 s |
| 745 | Yes | No | 0.77 s |
| 853 | Yes | No | 1.21 s |
| 986 | No | Yes | 108.95 s |
| 1285 | No | No | 1.48 s |
| 1405 | Yes, also top-1 | Yes | 288.49 s |
| 1654 | Yes | Yes | 183.51 s |
| 1665 | No | No | 21.81 s |
| 1740 | No | No | 1.27 s |

These four recoveries were missing from the original three-seed cap-100 results. They had already been recovered by other controlled diagnostics; they are not four newly solved benchmark records. Cap relaxation alone fixes their original three-seed failure.

Case 1740 still misses at 100 seeds with **zero cap stops** (67.12 s search/archive; 33.32 s verification). Thus increasing the cap is not sufficient for all misses. Uncapped does not mean exhaustive: seed-dependent fragment growth/commitment remains. This experiment does not prove any reference intrinsically unrepresentable.

## Inspect actual mappings

- [Case 1665, new cap 2,000 / three seeds](case1665_seeds3/viewer.html)
- [Case 1665, previous cap 100 / 100 seeds](../golden_mapping_diagnosis_20260907/case1665/viewer.html)

Both are standalone offline HTML, with reference, blind top-ranked result, and closest saved representative. Display indices are checked against the archive elements and bond matrix. Colors indicate **source-component origins**, not fragment boundaries; actual fragment decisions are listed separately. Explicit-H toggle is available. The reference does not specify individual H mappings.

The reference uses 18 core carbons and two sets of ten naphthyl carbons. The previous cap-100 top result borrows seven carbons from toluene and three from phosphine. The new cap-2,000 closest saved representative uses 18 + 10 + 9 intended carbons but one carbon from phosphine. “Closest” is an explicitly documented orbit-count diagnostic, not a claim of globally optimal alignment over all symmetry realizations.

The saved closest path (terminal 8977, cut R40–R41) gives a concrete chronology:

1. Seed R138 grows the complete ten-carbon region from R6.
2. Seed R13 grows the complete eighteen-carbon core from R1.
3. Seed R113 commits **R41–R113 (a C–H piece of phosphine)** to P1–P39.
4. Seed R128 grows the remaining nine-carbon region from R5, with the other carbon position already occupied.

These are archived transition decisions, not a constructed mechanism. This path illustrates a competing source claiming an atom before the intended building block; it does not prove every missed path fails identically. No cap stop was recorded anywhere in this three-seed run.

## Resource failure and retry

The first case-1665 / 100-seed run was cancelled after inspecting its cgroup: 37 OOM events and seven worker kills at the 32 GiB job limit. The pool replaced workers without completing the lost work. **This is an invalid/incomplete run, not a reference miss.** Its status and partial files are preserved.

A separate retry used 128 GiB, the same eight workers and identical search settings. Its wrapper checks cgroup OOM events every five seconds and aborts if a worker is killed. The 570-second process-group watchdog leaves cleanup time before Slurm's ten-minute limit.

**Case 1665 recovered at 100 seeds, cap 2,000, tolerance 1.0.** Every one of the first eight completed cut archives contains a representative with the exact reference chemical certificate. The no-cut archive alone contains a witness (terminal 243), using the correct 18 + 10 + 10 source-carbon origins. This is mapping recovery modulo endpoint chemical symmetry, not just atom coverage or independent atom-orbit agreement.

- [View the actual recovered witness](case1665_recovered/viewer.html); select “Verified recovered witness”.
- `partial_reference_check.json` preserves the eight witness mappings and their source cut archives.
- `high_memory_retry/results.json` explicitly marks **partial completed cuts**, not a completed sweep.

The remaining sweep was intentionally stopped after confirming recovery: Slurm allocation elapsed **544 seconds (9:04)**, peak batch RSS **43.05 GiB**, no observed OOM events in the retry. This is elapsed time to the stop after verification, **not full-sweep runtime or a precise first-hit search time**. All completed cut graphs are preserved. We did not rerun matching to build the viewer.

Overall, five of the ten selected cases have a recovery in this experiment: four at the original three seeds, and 1665 at the previously tested 100 seeds. Five remain missing (11, 602, 986, 1285, 1740); this does not establish that their references are impossible for AAM.

## Reproduction and records

- `results.csv` / `results.json`: each original attempt, timings and archive paths; verification time is separate from search/archive time.
- Full immutable-source snapshot and intermediate archives: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_cap2000_20260907`.
- High-memory retry: `high_memory_retry` under that directory.
- `bench/golden_cap_campaign.py` and `hpc/golden_cap_probe.sbatch`: launcher, watchdog and report export.
- `bench/view_golden_mapping.py`: visualization from saved archives only; no matching rerun.

HTML rendering and H/selection controls were checked in headless Chromium. Fifteen benchmark-preparation, evaluation, viewer and watchdog tests passed.
