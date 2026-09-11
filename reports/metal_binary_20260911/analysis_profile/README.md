# What the event-analysis timer actually measures

The 577.155 / 614.502 CPU seconds in the binary-metal experiment include extracting saved mappings, counting raw-WBO events, assigning symmetry-normalized pattern identities, and independently checking described witnesses. They are not timings of bond-event counting alone.

A read-only instrumented replay of the original arm on two saved cases measures these components. No mapping searches were rerun. Both complete score histograms and minimum pattern IDs match the frozen analysis exactly. CPU measurements exclude archive loading; the profile wrappers introduce small measurement overhead.

| Component, CPU seconds | Case 68: 135 atoms, 250 representatives | Case 101 (PR7): 28 atoms, 34 representatives |
|---|---:|---:|
| All timed analysis | 62.264163 | 0.071523 |
| Extract representative vectors | 0.410358 | 0.006789 |
| Vectorized event counts | 0.053229 | 0.000472 |
| Nauty graph canonicalization | 60.551854 | 0.009462 |
| Independent scalar event verification | 0.191044 | 0.015648 |

For case 68, six graph-canonicalization calls consume 60.552 of 62.264 CPU seconds (97.25%). Actual vectorized 0.5/0.3 bond-event counting takes 0.053 seconds. The profiler verifies the entire histogram, so this short count covers all 250 representatives, not a sample.

The event-equivalence implementation builds a colored graph containing every atom pair, including nonbonded pairs, to preserve exact raw-WBO score responses under symmetry. At 135 atoms that is 9,045 pair vertices plus 135 atom vertices before event markers. Canonicalizing this graph can be expensive. The current loop also describes temporary best-scoring representatives before knowing the final minimum. These are implementation choices of the benchmark pattern comparison; the threshold rule itself does not require this work.

This diagnosis covers two selected cases, including one chosen for high analysis cost. It does not establish the exact counting-versus-canonicalization split of all 577/615 seconds. Another source of avoidable work visible in the implementation is reconstructing a mapping dictionary once per atom while extracting vectors; it costs 0.410 seconds in case 68.

## SLAP comparison and timing boundaries

The archived SLAP single-edge sweep records 47.046 CPU seconds for scoring its deduplicated bidirectional candidate union. That phase scores 1,660 unique label partitions and optimizes H assignments. Its timer does not include the later enumeration and symmetry comparison of all alternative minimum-event patterns. It also uses the older uniform-0.5 / bond-floor rule. Consequently 47.046 versus 577.155 is not a comparison of identical postprocessing tasks.

The forward-only SLAP sweep mapping time is 693.756 CPU seconds for 140 cases; the corresponding historical AAM mapping time is 942.576 CPU seconds (the fresh original arm here is 941.484). The archived forward-only SLAP sweep scoring time is unavailable because scoring occurred after combining both directions. It must not be obtained by halving 47.046 seconds.

See [forward comparison](../../holdout_forward_cap1000_seed1_20260910/README.md) and [bidirectional sweep scoring](../../holdout_slap_xyz_sweep_20260910/README.md). The profiler, raw component times, exact job command and completed scheduler accounting are stored alongside this note.
