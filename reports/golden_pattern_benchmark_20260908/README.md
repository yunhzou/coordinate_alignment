# Full Golden validation of compressed-pattern collection

**Concrete candidates recover the reference in 1,839 / 1,851 reactions (99.35%)
bidirectionally, within +6 bond events of the best collected candidate.**
Smaller-to-larger alone reaches 1,821 / 1,851 (98.38%), within +5 events.
These are displayed/saved physical mappings, not reference feasibility inferred
only from an unexpanded compressed family. No new AAM search was run.

## Scope and results

All 1,851 reactions and both directions were checked: 3,702 directional inputs.
3,689 collection workers completed successfully after repairing two baseline
preparation timeouts. Thirteen directional archives were unavailable from the
previous AAM run; they remain missing, not silently counted as successful searches.
Every accuracy denominator remains 1,851, including missing inputs.

| Extra bond events | Bidirectional reference recovery | Percentage |
|---:|---:|---:|
| 0 | 1,692 | 91.41% |
| 1 | 1,757 | 94.92% |
| 2 | 1,817 | 98.16% |
| 3 | 1,830 | 98.87% |
| 4 | 1,834 | 99.08% |
| 5 | 1,838 | 99.30% |
| 6 | 1,839 | 99.35% |

The event window includes all collected candidates at maximum mapped-heavy and
mapped-total counts, within the stated event difference from their minimum.
It is **not top-k accuracy** and is not an exhaustive bound over unseen patterns.
Scores include explicit H with event tolerance 0.5. The archived matching used
tolerance 1.0, seed count 10 and branch cap 100. Reference identity is the benchmark
heavy-atom relation modulo endpoint chemical symmetry; H identity is not annotated
by the benchmark. No energy or new stereochemical feasibility filter was added.

Using the same available archive rankings, saved representatives recovered 1,835
bidirectional cases before extraction; the collector adds cases 44, 865, 1033 and
1786. Single-direction recovery rises from 1,811 to 1,821. The bidirectional
undetected set is exactly the previous full-family check's nine misses plus
three unknowns: 7, 19, 590, 871, 986, 1228, 1285, 1358, 1377, 1475, 1553, 1799.
Thus no previously recovered bidirectional reference was lost.

## What changed

Exact interchangeable atom groups, such as hydrogens attached to one carbon,
are represented by mapping counts between groups. Correlated motions of whole
groups, including ring symmetries, remain coupled. This removes factorial
within-group permutations from class-exclusion queries without changing pattern
identity. The AAM search engine and saved archives remain unchanged.

Same-input comparisons on three reference-blind top-ranked paths completed in
0.076, 0.108 and 0.180 seconds. The previous collector timed out at 3 seconds on
each; the pattern keys were identical. See `exclusion_comparison.json`.

The full pass visited 4,541,718 families, skipped 306,043 identical query histories,
and saved 1,608,664 additional directional patterns. **Extraction is not globally
exhaustive:** 20,916 visited families remain unresolved, and some archives were
not fully traversed within their budgets. Only 634 single-direction cases and
409 two-direction cases have exhaustive archive extraction certificates. Original
compressed archives retain all unexplored alternatives; do not interpret 99.35%
reference recovery as 100% extraction completeness.

## Parallelism, timings and provenance

- Main pass: 16 nodes, up to 24 archive workers each (384 allocated logical CPUs).
  Two stalled allocations were replaced. Main batch elapsed: 11 minutes 57 seconds.
- Two tail baselines contained 56,667 and 17,499 saved classes. Reuse of their saved
  scores/certificates removed redundant work; exact explicit-mapping certificates
  were then parallelized across 16 processes per archive. Live workers on the
  larger case used approximately 94–95% of a CPU each.
- Final parallel tail jobs took 180 and 100 seconds including input loading and
  persistence. Initial pass through final repairs spanned 20 minutes 37 seconds,
  including diagnosis/restarts. These are not all clean-run timing measurements.
- Soft per-path budget: 3 seconds; per-direction collection budget: 60 seconds.
  Preserving all baseline candidates can exceed that soft budget. Five-minute
  process watchdogs apply; parallel retries kill the whole worker process group.
- Archive loading/hashing, extraction, persistence and reference verification
  are separately recorded. Summed extraction wall time is not cluster elapsed.
  Original parallel retry `worker_cpu_seconds` fields contain parent CPU only;
  use `slurm_accounting.tsv` for their inclusive CPU totals. Future instrumentation
  now separately captures child CPU. Do not publish those parent-only fields as
  total CPU cost.
- Immutable main, repaired and parallel worker snapshots are retained. Their
  source hashes are in the manifest files; all prior partial retry artifacts are
  preserved. Final jobs for the parallel repair were 452148 and 452150; job 452149
  exited without doing work and is also retained in the accounting log.

Validation: **84 focused tests passed**, including exact-pattern oracles, partial
mappings, coupled symmetry, resumption, parallel/serial output identity and order,
shared feasibility/scoring, and the existing graph/viewer tests.

## Artifacts and viewer

Complete machine-readable run and witnesses:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_pattern_collection_20260908_v2`

Per direction: `results/<index>/<direction>/patterns.jsonl.gz` stores physical
mappings, events, action traces and originating paths; `paths.jsonl.gz` records
visited families and completion. `candidates.json` and `summary.json` store scores,
reference labels and timing. All original AAM archives remain at the publication
source path recorded in the manifests.

The updated case-1033 viewer shows all 181 collected candidates, including the
reference-equivalent mapping at rank 7, and is rendered from saved mappings only:

`/h/399/yunhengzou/coordinate_alignment/reports/golden_pattern_benchmark_20260908/case1033/viewer.html`

GitHub stores `case1033/viewer.html.gz`; decompress it before opening. The plain
HTML remains on the cluster. The report and per-case tables are committed here.
