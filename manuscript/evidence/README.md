# Evidence used in the manuscript

Snapshots were copied from existing reports on 10 September 2026. Exact input
paths and SHA-256 hashes are in `sources.json`. Rebuilding figures consumes
these copies; it does not rerun a benchmark or read reference labels during search.

| Snapshot | Use | Original report/workload |
| --- | --- | --- |
| `seed_comparison.json` | Figure 3, Table 1, GRAFT seed recovery and CPU | `reports/aam_seed1_2_20260910/comparison.json` |
| `slap_sweep.json` | Figure S1, supporting SLAP cut ablation/completion | `reports/slap_sweep_cut_20260910/summary.json` |
| `publication.json` | Earlier frozen recovery/ranking | `reports/golden_publication_20260908/analytics.json` |
| `collection.json` | Figure S2, concrete event windows | `reports/golden_pattern_benchmark_20260908/benchmark_report.json` |
| `competitors.json` | Supporting Table 2 | `reports/golden_competitors_20260908/summary.json` |
| `holdout_scores.json` | Figure 4a, strict-cap full-H scores | `elementary140_tol1_20260909/tolerance_per_case.json` |
| `holdout_overlap.json` | Figure 4b/c, membership and extras | `elementary140_tol1_20260909/overlap/summary.json` |
| `holdout_manifest.json` | Native-input protocol | `elementary140_tol1_20260909/manifest.json` |
| `adaptive_no_sweep.json` | Separate adaptive no-sweep Golden/holdout result | `adaptive_full_20260910/comparison/summary.json` |
| `adaptive_manifest.json` | No-sweep version, agenda/cap/budgets | `adaptive_full_20260910/manifest.json` |

GRAFT is the manuscript's working method name. Original report identifiers (AAM)
and all source snapshots remain unchanged. Main benchmark panels use the frozen
`98b01b1` baseline; the earlier-engine and adaptive-policy results remain separate
in supporting information. Figure 2 illustrates recorded growth on case 64.

Report paths are relative to the manuscript's parent worktree. Workload paths
are relative to `/project/yunhengzou/coordinate_alignment/aam_benchmarks`.
The adaptive report is a captured snapshot of a report labeled “live”; its saved
adaptive results contain all 3,702 Golden and 280 holdout directional searches.
Verification unknowns remain explicit. It is not a pure cut-only ablation.

Illustration data:

- `growth_trace.json`: carbocation case 64, 18 atoms, 17 accepted extensions,
  two retained terminal mappings, maximum 42 live candidates.
- `growth_trace_1.json`: TEMPO case 1, 57 atoms, four seeds/fragments,
  53 accepted extensions, eight deferred-boundary observations, one terminal.
- Corresponding `case64_input.json` and `case1_input.json` are original endpoint
  elements, coordinates and WBO matrices. Source hashes and the bookkeeping
  extension hash are embedded in each trace.
- `growth_trace_135.json` / `case135_input.json` preserve an inspected candidate
  example that was not selected for a figure or movie.
- `animation_data.json` is derived from the selected traces. Its view bases
  choose an orthographic camera for readability; atomic coordinates and mapping
  events are unchanged. A seed-end event may display one witness while the
  retained graph has other terminals. The full graph is available in the viewer.

Scientific claims also rely on the human-readable reports for accounting scope,
exceptions, case follow-ups, evaluator definitions and excluded experiments.
Those reports remain in the worktree; consult them before changing a comparison.
The compact bundle does not contain every original candidate or timing log.

The [scientific novelty audit](novelty-audit-20260910.md) compares the frozen
`98b01b1` search design with primary literature and selected author implementations.
It identifies a plausible algorithmic contribution in the combined representation
and transition rules, while documenting substantial precedents for the individual
techniques. It is a targeted assessment, not proof of historical priority.
[The source manifest](novelty-source-manifest.json) pins inspected code and records
availability limits. This assessment adds no new benchmark results; its proposed
component ablations remain future work.

The follow-up [mechanism comparison](mechanism-comparison-20260910.md) makes the
comparison at the level of state and transition rules. It inspects the frozen
native loop as well as Python, distinguishes aggregate fragment commitment from
atom-pair and fragment-shape branching, and compares the actual deduplication
keys. OpenChemLib's concrete growth routine and FMCS are additional close
comparators. Newly cached source hashes are in
[mechanism-source-manifest.json](mechanism-source-manifest.json).

The expanded [literature and mechanism assessment](review-led-novelty-20260910.md)
is also available as a [linked PDF](review-led-novelty-20260910.pdf). It follows
five reviews/overviews and adds direct comparisons with Indigo, contemporary RDT,
SynKit's approximate growth routine, AMLGAM and other original methods.
The [coverage ledger](review-coverage-20260910.csv) contains 44 entries: five
reviews/overviews, 37 methods or implementation leads, and two general foundations.
These are coverage entries, **not 44 fully audited algorithms**. Coverage ranges
from inspected code and pseudocode to unresolved historical leads.
The [query log](review-search-log-20260910.json) records 36 additional discovery
queries, and the [source manifest](review-source-manifest-20260910.json) records
download hashes, versions, and access limits. Retrieved third-party PDFs and
source caches remain in ignored build storage rather than the distributable bundle.

The [SLAP counting audit](slap-sweep-count-audit-20260910/README.md) re-scores
all 572,753 saved candidate outputs and reproduces 1,795/1,851. It distinguishes
52 completed standalone misses from four unresolved cases, records their nine
unfinished variants, separates the older expanded union's 1,797 result, and
checks the case sets behind the quoted AAM CPU figures. No mapper was rerun.
