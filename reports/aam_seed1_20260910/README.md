# One-seed bidirectional AAM ablation

Golden: **1,833/1,851 recovered (99.03%)**, 17 verified absent and one unknown.
Full mappings: **140/140 holdout cases**; there is no annotated holdout ground truth.
Best representative event counts agree with ten seeds in 139/140 holdout cases;
case 123 changes from 19 to 27 events.

All 3,982 searches and all 3,982 analysis processes completed with exit zero.
Frozen engine `98b01b1`, one deterministic seed order per cut, both directions,
full sweep, cap 100, tolerance 1.0, explicit H and eight CPUs per directional call.
No AAM core change. Defaults remain unchanged.

Full Golden search: 27,084.86 CPU-seconds; holdout: 1,475.45 CPU-seconds.
Separate search/analysis execution spans: 273.22 / 143.53 seconds.
CPU excludes measured persistence/loading; spans include persistence and dispatch.

See the [combined report](../aam_seed1_2_20260910/README.md) for matched timing
comparisons, all caveats, scheduler intervention and reproduction instructions.

Full archive root:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_one_seed_bidirectional_20260910`
