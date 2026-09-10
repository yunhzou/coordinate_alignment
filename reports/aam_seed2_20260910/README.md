# Two-seed bidirectional AAM ablation

Golden: **1,834/1,851 recovered (99.08%)**, 14 verified absent and three unknown.
Full mappings: **140/140 holdout cases**; there is no annotated holdout ground truth.
Best representative event counts agree with ten seeds in 139/140 holdout cases;
case 123 changes from 19 to 25 events.

All 3,982 searches and all 3,982 analysis processes completed with exit zero.
Frozen engine `98b01b1`, two deterministic seed orders per cut, both directions,
full sweep, cap 100, tolerance 1.0, explicit H and eight CPUs per directional call.
No AAM core change. Defaults remain unchanged.

Full Golden search: 49,161.91 CPU-seconds; holdout: 2,417.94 CPU-seconds.
Separate search/analysis execution spans: 361.65 / 184.13 seconds.
CPU excludes measured persistence/loading; spans include persistence and dispatch.

See the [combined report](../aam_seed1_2_20260910/README.md) for matched timing
comparisons, all caveats and reproduction instructions.

Full archive root:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_two_seeds_bidirectional_20260910`
