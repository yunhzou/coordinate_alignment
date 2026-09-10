# One- and two-seed bidirectional AAM ablations

## Results

**One seed reaches 99.03% verified Golden recovery at 7.50× lower search CPU
cost than ten seeds; two seeds reach 99.08% at 4.21× lower cost.**
Both fresh experiments completed all 3,982 directional searches and all 3,982
analysis processes with exit zero. No hard process watchdog was reached.

| Method | Recovered / 1,851 | Recovery | Verified absent | Unknown | Search CPU-seconds / reaction |
| --- | ---: | ---: | ---: | ---: | ---: |
| AAM, one seed | 1,833 | 99.03% | 17 | 1 | 11.72 |
| AAM, two seeds | 1,834 | 99.08% | 14 | 3 | 20.85 |
| AAM, three seeds, saved | 1,836 | 99.19% | 10 | 5 | 29.77 |
| AAM, ten seeds, saved | 1,836 | 99.19% | 5 | 10 | 87.87 |
| SLAP + single-edge sweeps, saved | 1,795 | 96.97% | — | — | 17.62 |

Recovery uses **all 1,851 cases**, with unknowns retained in the denominator.
Timing uses **the same 1,807 mutually completed reactions** for every row,
including both directions and the full sweep. Measured saving/loading is excluded.
The one-seed cost is 0.665× SLAP+sweep; two seeds cost 1.184× SLAP+sweep.
These are instrumented compute comparisons, **not equal-resource wall-time
speedups**: AAM includes IPC/process setup, while the SLAP workflow timer excludes
worker startup. Interrupted work in the saved baselines is not included in these
paired successful-work timings. SLAP accuracy here is its standalone sweep result,
not the union with an earlier additional experiment.

The AAM rows use the native-reuse engine `98b01b1`. They must not be confused with
the earlier publication engine `bb0a7d7`, whose frozen result is 99.35%.

### What is lost with fewer seeds?

Counts alone hide changes in case sets. Compared with saved ten seeds:

- One seed no longer verifies cases 44, 602, 850, 912, 1314, 1384 and 1568.
  Case 912 is unresolved; the other six are verified absent from the saved output.
- Two seeds no longer verifies 44, 602, 850, 1314, 1384 and 1568.
  Cases 850 and 1568 are unresolved; the other four are verified absent.
- Both fresh runs verify 1033, 1786, 1799 and 1823, which were unknown in the saved
  ten-seed results. This reflects completed search/verification, not evidence that
  fewer seed orders explore more possibilities.

Indices above are zero-based benchmark indices. No union with another seed run
was used to obtain either reported accuracy. "Unknown" is not a confirmed miss;
soft symbolic-verification limits can remain even when the process exits normally.

Both seed settings return full explicit-atom mappings for **140/140 XYZ/WBO
holdout cases**. Their best saved representative event counts agree with ten seeds
in **139/140** cases. In case 123, the counts are 27 (one seed), 25 (two or three
seeds) and 19 (ten seeds). The holdout has no annotated ground truth: these are
mapping availability and representative-score comparisons, not accuracy or proof
that all alternatives are preserved.

### Full workload cost and execution spans

| Fresh run | Full Golden search CPU-hours | Holdout search CPU-hours | Search span | Analysis span |
| --- | ---: | ---: | ---: | ---: |
| One seed | 7.524 | 0.410 | 273.22 s | 143.53 s |
| Two seeds | 13.656 | 0.672 | 361.65 s | 184.13 s |

Full totals include all 1,851 Golden and 140 holdout cases, including the difficult
cases excluded from the paired baseline timing. Spans cover Golden and holdout
together, from first worker start to last finish in each phase. They exclude the
initial queue but include saving, loading and dispatch. They are separate phase
measurements, not controlled end-to-end latencies. Each experiment used at most
16 × 32 search CPUs (512) and 16 × 8 analysis CPUs (128). The two experiments ran
concurrently. The longest individual search including startup/I/O was 139.46 s
for one seed and 228.55 s for two seeds.

Slurm left one completed one-seed allocation in `COMPLETING`, blocking the
dependent analysis array despite all 3,982 successful saved search statuses.
After checking these statuses and Slurm's successful exit record, the dependency was
cleared with `scontrol update JobId=467414 Dependency=`. This started verification
of existing archives; no search was repeated. Search jobs: 467413 and 467430;
analysis jobs: 467414 and 467431. Scheduler cleanup delay is excluded from the
separate execution spans. Slurm TotalCPU fields on this cluster are unusable;
the unchanged parent-plus-child SearchProfiler supplies the CPU measurements.

## Protocol

Use the same frozen native-reuse AAM engine `98b01b1`, native binary, input
structures and search/profiling/evaluation adapters as the completed three-seed
and saved ten-seed benchmark. Only `seed_count` changes to one or two per cut.
These orders are prefixes of the same deterministic ten-order sequence.

Both directions, full sweep, cap 100, matching tolerance 1.0, explicit H,
root seed 42 and eight CPUs per directional AAM call remain unchanged.
Test all 1,851 Golden reactions and all 140 XYZ/WBO holdout cases separately.
The holdout has no annotated ground truth. Golden recovery uses the strict
heavy-atom relation modulo endpoint chemical symmetry, including unmatched atoms,
among saved compressed families; it is not top-1 accuracy.

Each experiment has a 512-CPU search budget, for 1,024 CPUs combined. Dynamic
queues distribute independent tasks. Search and evaluation are separate phases,
with 300-second search and 240-second evaluation process watchdogs. Save every
cut and final AAM archive. Keep incomplete outcomes in accuracy denominators.

Compare measured search CPU excluding persistence/loading on common completed
cases against three seeds, ten seeds and SLAP+sweep. Publish actual recovery
counts and case-set differences, not an assumption that fewer seeds preserve
accuracy. Cluster execution spans are distinct from CPU cost and exclude only
the initial queue, not persistence or dispatch.

Artifact roots:

- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_one_seed_bidirectional_20260910`
- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_two_seeds_bidirectional_20260910`

The three-seed run and ten-seed baseline remain untouched. Package defaults and
the AAM core are unchanged.

## Reproduction and saved evidence

The seed-prefix, concurrent task-claiming and comparison-denominator tests pass:
4 tests and 27 subtests. Source/native binary and adapter hashes match the saved
three-seed run; search configurations differ only in `seed_count`.

Use `bench/aam_seed_ablation.py collect` and `publish` on the artifact roots to
rebuild compact reports **without rerunning mapping or evaluation**. Run
`bench/compare_aam_seed_counts.py` to regenerate this directory's comparison JSON
and CSV from the saved reports. The report records source-file hashes and all
1,807 common timing indices. Full per-cut and final compressed AAM archives,
evaluation results, timings, environment information and task logs remain in
each artifact root.

- [Combined machine-readable comparison](comparison.json)
- [Combined timing and accuracy table](comparison.csv)
- [One-seed per-case results](../aam_seed1_20260910/per_case.csv)
- [Two-seed per-case results](../aam_seed2_20260910/per_case.csv)
- [Three-seed report](../aam_seed3_20260910/README.md)
- [SLAP sweep report](../slap_sweep_cut_20260910/README.md)
