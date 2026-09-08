# Golden: ten seeds, cap2000

Retest the 23 cases not recovered or unresolved in the ten-seed cap100 campaign. Retain its nine verified recoveries without rerunning them. Use ten seeds per cut, branch cap2000, tolerance1.0, explicit hydrogens, smaller-first default orientation, and the full single-edge cut sweep. No reverse-direction rescue or core algorithm changes.

## Result

**No additional verified recoveries; seven confirmed misses and 16 unresolved.** Demonstrated combined coverage remains 1,827/1,851 (98.70%). All 23 workers finished within their ten-minute job budgets (longest552seconds).

- Confirmed misses: 7, 576, 986, 1285, 1475, 1553, 1739.
- Incomplete searches: 19, 590, 780, 833, 871, 1228, 1358, 1377, 1380, 1568.
- Completed searches but incomplete verification: 603, 845, 865, 867, 1033, 1574.
- Cases7 and1739 completed without hitting the cap (maximum live branches624 and14 respectively), yet still missed. Other completed searches hit cap2000.

This bounded experiment does not establish whether the 16 unresolved searches contain or could discover the reference. Their intermediates are retained for further inspection without discarding the work.

See `summary.json` for outcomes, `cases.json` for per-case timings and phase exit statuses, and `witnesses.json` for any new verified mappings. The combined percentage is a diagnostic union with saved earlier results, not an all-dataset high-cap benchmark. Original unresolved case1793 remains outside the requested 32 misses.

## Resources and reproducibility

- Slurm array442512: 23 reactions, eight CPUs and 128GiB allocated per reaction.
- Search watchdog300seconds; verification watchdog240seconds; overall job limit10minutes.
- Incomplete searches are checked using saved completed cuts. A negative result from incomplete search or verification remains `unknown`, not a confirmed miss.
- Frozen source, exact inputs, cut intermediates, completed checkpoints, logs, and evaluations:
  `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_seeds10_cap2000_20260907`
- `manifest.json` records configuration and frozen-source hashes. `slurm_accounting.psv` records resource accounting.
- 141 regression tests passed, including a check that only the requested seed and branch budgets change.

Higher cap does not guarantee an exhaustive search. Even an uncapped miss concerns the sampled seed/cut search, not every mathematically possible mapping.
