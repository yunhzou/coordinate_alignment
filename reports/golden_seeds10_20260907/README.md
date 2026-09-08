# Golden: more seeds, same default direction

Test the 32 misses from the original fixed-direction Golden run. The only search configuration change is **3 to 10 seeds per cut**. Preserve smaller-first orientation, cap100, tolerance1.0, explicit hydrogens, and the full single-edge cut sweep. No reverse-direction rescue or core algorithm change is included.

The first three random seed orders are preserved; ten seeds extend the original stream. Tests check the configuration-only change and seed-order prefix for uncut and cut graphs.

Result: **9 recovered, 17 not recovered, 6 unresolved**. Recovered cases: 602, 737, 775, 846, 850, 964, 1230, 1314, 1400. Case846 has a verified positive in an incomplete search. Cases590 and833 have incomplete searches; cases865,1033,1358,1568 have completed searches but incomplete verification. Combined with the prior fixed-direction run, demonstrated reference coverage increases from 1,818 to **1,827/1,851 (98.70%)**. Ten seeds do not recover all 32.

See `summary.json` for final counts, `cases.json` for individual outcomes and phase statuses, and `witnesses.json` for recovered mappings. A `not_recovered` result means the verifier excluded the reference from the completed saved search, not that no possible search could find it. `unknown` means search or verification was incomplete.

The union percentage combines new recoveries with the original 1,818 recovered records. It is not an all-1,851 ten-seed benchmark. Original unresolved case1793 is outside these 32 selected misses. Reference recovery is mapping-family coverage, not top-one accuracy.

## Resources and persistence

- Four CPUs per reaction; up to 32 concurrent reactions.
- Search watchdog300seconds, verification watchdog240seconds, Slurm limit10minutes.
- Jobs442479 and442511. The latter replaces one allocation stuck configuring before Python started; completed searches were not rerun.
- Frozen search source, inputs, cut intermediates, complete checkpoints where finished, evaluation results and logs:
  `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_seeds10_20260907`
- `manifest.json` records exact configuration and frozen-source hashes; `slurm_accounting.psv` records allocation and runtime data.

No chemistry heuristics or benchmark-specific search rules were added.
