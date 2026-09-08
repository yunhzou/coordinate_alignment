# Case 590: forced reverse-direction sequential AAM

The existing random-policy run was already R→P (55 explicit atoms on each side).
This experiment forces P→R, retaining 10 random seeds, branch cap 100, tolerance
1.0, explicit hydrogens, and the ordinary single-edge cut sweep. Each direction
sweeps its own source edges. No core algorithm changes or reference-guided search.

Result: **reference not recovered in the returned reverse-search families**.
Verification completed without a timeout or unknown query. It examined 52,738
terminal representatives, 1,330 distinct heavy representatives, and rejected
1,129 symbolic heavy projections. The compressed-family evaluator reports
`not_recovered`; no full explicit symbolic query was needed after the preceding
representative and projection checks. Reference equivalence includes endpoint
chemical symmetry, not just literal pair comparison.

The search hit cap 100, so this is not an exhaustive impossibility statement.
The previous R→P reference evaluation remained unknown; do not reinterpret that
as a verified miss or a recovered mapping. No new recovery was established here.

Timings: 66.96 seconds for search through final archive persistence on 16 Slurm
CPUs; 148.75 seconds inside reference evaluation. These are separate stages, not
a claim of equal-resource speedup over the older run. All 56 reverse cut graphs
and the complete compressed AAM archive were saved.

Full archive:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden590_reverse_20260908/cuts/aam.pkl.gz`

Slurm job 443209. Search watchdog 300 seconds; verification watchdog 260 seconds;
allocation 10 minutes. Both phases completed normally.

Runner: `bench/golden_forced_direction.py`. Direction/evaluation regression tests:
15 passed. Original R→P evidence is under
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_seed_diversity_20260907/random/590`.
