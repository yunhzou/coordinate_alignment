# Conditioned symmetry and fine-grained growth reuse (in progress)

First paired pilot: **1.60x total compute-CPU speedup** for native conditioned
symmetry plus the prior dependency/native-scheduler optimization, compared with
the frozen original Python scheduler/native growth. This is four directional
cases, original seed index zero of ten, every cut, cap 100, explicit H,
matching tolerance 1.0. It is not yet the full ten-seed confirmation.

| Mode | Compute CPU seconds | Search | Symmetry |
|---|---:|---:|---:|
| Frozen original baseline | 118.354 | 71.766 | 46.536 |
| Prior native scheduler + whole-fragment repair | 95.949 | 50.595 | 45.246 |
| Prior repair + native conditioned symmetry | 74.053 | 53.031 | 20.920 |
| Prior repair + individual-extension cache | 101.708 | 56.420 | 45.207 |
| Both new components + prior repair | 80.414 | 58.965 | 21.349 |

All **1,656** full cut-graph comparisons agree, including ordered generators,
every state/transition and cap/stop records. The extension cache currently costs
more than it saves and is not enabled by default. The symmetry prototype also
remains opt-in. No seed, cut, chemistry or branch-admission rule was changed.

The initial profiling observation was 309 group requests on one saved long-tail
cut: only about 13% of its instrumented symmetry time was inside pynauty autgrp;
most was building role labels/partitions. Therefore this prototype reuses those
colors and complete conditioned partitions on a fixed native graph. It does
**not** replace nauty with a general permutation-subgroup stabilizer algorithm.

The extension prototype caches each candidate's raw extension children before
ordinary whole-frontier dedupe/cap admission. Keys include every relevant
source-topology, assignment, occupied-target and support constraint; negative
answers are cached too. Source topology reads are repeated during key construction
so enclosing whole-fragment dependency tracking remains sound.

Raw pilot, frozen original and experimental engines, logs and every cut graph:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/conditioned_reuse_pilot_20260910_KciqQq
```

Slurm array 456603: four single-CPU paired tasks; all modes within each task run
on the same node in counterbalanced order. Each mode has a 300 s hard watchdog;
each allocation has a ten-minute limit. Compute excludes queue, imports, input
loading, audit encoding and persistence. No downstream pattern collection or
reference chemical-accuracy claim is included. All comparisons completed.

The source baseline is commit `f158864`; the frozen native/source copy comes
from the real-TS run, whose production core is the same baseline. Python typed
fragment/search-graph objects and correlated symmetry representation remain.
