# Independent fragment generation + compatible assembly: case 590

Experimental workflow only. Sequential AAM, growth, native matching, and the
existing symmetry representation were not changed.

## Design

- Every requested seed calls the existing fragment matcher with an empty mapping.
- Start with every heavy atom; a second arm also seeds every explicit hydrogen.
- Optional sweep uses the existing no-cut plus each single-edge cut, including X–H.
- Exact duplicate records are removed during detection. Finalized identical
  relations are deduplicated across cuts when the cuts outside the fragment do
  not alter that recorded relation. Different witnesses are not assumed equivalent.
- Assembly chooses compressed families and their correlated target actions
  jointly using SMT. All explicit source atoms need fragment support; mapping is
  injective. Overlapping fragments are allowed when their assignments agree.
- No subfragment trimming, invented singleton completion, bijection expansion,
  chemistry ranking, or reference-guided detection. A returned assembly is also
  checked with the existing AAM family-query implementation.
- The reference-constrained query is a separately labelled diagnostic using
  literal reference identities, not benchmark symmetry-normalized accuracy.

## Measured result

The source has 31 heavy atoms / 55 explicit atoms. There are 58 sweep conditions.

| Detection arm | Seed searches | Cap hits | Detection wall seconds |
|---|---:|---:|---:|
| Heavy seeds, no cut, cap 100 | 31 | 4 | 0.148 |
| Heavy seeds, sweep, cap 100 | 1,798 | 257 | 0.550 |
| All explicit atoms, sweep, cap 100 | 3,190 | 261 | 0.810 |
| All explicit atoms, sweep, cap 2000 | 3,190 | 0 | 0.982 |

Sweep jobs used 16 Slurm CPUs. Times exclude scheduling and interpreter startup;
they include the detection pool and per-cut persistence. Summed cut-task elapsed
time for the cap-2000 arm was 12.31 seconds, not a measured CPU-time counter.

**Cap 2000: 1,455 unique stored families; no compatible complete assembly.**
Blind symbolic assembly returned UNSAT in 6.89 seconds. The separate literal
reference query returned UNSAT in 4.19 seconds. Both use the same saved pool.
These times exclude loading the saved pool. No timeout occurred in these checks.

The no-cut pool did not support source H atoms 34,35,42,43. After full-atom sweep,
all atoms occur in the pool, but their fragment placements cannot be combined
into a full injective mapping. Even an initial stricter partition-only assembly
was tested; allowing compatible overlap did not change the negative outcome.
The stricter intermediate results remain archived as `partition_*.json`.

This is a negative result for **these independently generated families**, not
proof that the reaction is impossible or that every independent strategy fails.
All seeds and single cuts do not guarantee retention of all useful smaller
fragments. Source-only partitions exist, so simply adding up source coverage is
insufficient; target placement compatibility matters.

## Reuse and evidence

Code: `src/rxn_core/independent_aam.py`; runner: `bench/golden_independent.py`.
Detection checkpoints can be reassembled without rerunning matching.

Full archives:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_independent_590_all_cap2000_20260908`

Sibling arms are named `golden_independent_590_nocut_20260908`,
`golden_independent_590_sweep_20260908`, and
`golden_independent_590_all_sweep_20260908` in the same parent directory.

Slurm jobs: 443206, 443207, 443208; all detection phases completed. Each had a
300-second process watchdog and six-minute allocation. Assembly phases had a
120-second internal budget and 150-second external watchdog.

Tests: 60 passed across the independent workflow, existing family queries, and
search-graph tests. Tests include correlated symmetry, rejection of inconsistent
assignments, explicit-H coverage, compatible overlaps, and cross-cut deduplication.
The six experimental tests also pass with native matching disabled. A positive
control mapping all case-590 reactants to themselves produced a certified complete
assembly in 0.72 seconds; see `identity_control.json`.
No full Golden benchmark was run for this experimental workflow.

Slurm allocation elapsed times were 11, 5, and 9 seconds respectively, including
launch/startup overhead. These should not be confused with the internal detection
timings above. Scheduler-reported batch RSS was too low to characterize workers'
actual combined memory, so no peak-memory claim is made.
