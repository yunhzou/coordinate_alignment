# Repository branches and preserved baselines

Organized on 10 September 2026. `main` remains the default branch. This cleanup
changes reference names and adds documentation; it does not merge experimental
acceleration into the default implementation or rewrite implementation history.

## Where to work

| Branch | Purpose | Preserved implementation / starting point |
| --- | --- | --- |
| `main` | Default project entry point and branch guide | Implementation at `2d5ac7a`; this organization adds documentation only |
| `paper/continuous-fragment-growth` | Manuscript, figures, movies, SLAP sweep and seed benchmarks | Continues `experiment_slap_sweep_cut` from `2704cd9` |
| `stable/pre-acceleration` | Mature sequential pipeline before the September 9–10 cut-replay, dependency-repair and conditioned-reuse work | **`3e9a70a2ab23b2a9e9a0a8d374d951531e1fb0b0`** |
| `stable/benchmark-baseline` | Frozen native-reuse engine used by the current seed ablation, before adaptive-policy acceleration | **`98b01b175eeed31f70d13e7cbf178b80bf07c9e0`** |
| `stable/legacy-aam` | Earlier exact-symmetry stable backup | `a604d1a2e75112a02bf930825707d108a6987b87` |
| `dev/acceleration` | Adaptive and shared-policy acceleration development | `73d37c7393fc59f29005183b502e519b58a4fe3e` |
| `dev/native-engine` | Separate native-engine line, including its unmerged target-element-table fix | `4088b2739e2ce4dfc3796abfed3cec6cc0198b00` |
| `dev/native-index` | Separate native index/chirality/path work | `fc8eb8b34137213096b9cc87e521a4c1ac7ea497` |

Start reading at the [manuscript folder](https://github.com/yunhzou/coordinate_alignment/tree/paper/continuous-fragment-growth/manuscript).
It contains the compiled PDF, five editable figures, two MP4/GIF movies, a
standalone interactive viewer, evidence snapshots and a downloadable source bundle.
The paper is a first draft; author information and specified provenance details
remain to be finalized.

Keep the stable branches at their recorded commits. Make new changes on a
development branch. No hosting-level branch-protection rule is implied by a
`stable/` name; the annotated baseline tags below pin the original commits too.

## Immutable baseline tags

Several different acceleration efforts exist in the history. All relevant
pre-acceleration boundaries are retained to avoid conflating them.

| Annotated tag | Commit | Meaning |
| --- | --- | --- |
| `baseline/sequential-pre-acceleration-20260902` | `23b1e8bef1d455e01cba69dd3e049d124b0cde7b` | Sequential implementation immediately before the earlier September 2 acceleration series beginning at `4f382aa` |
| `baseline/pre-reuse-20260909` | `3e9a70a2ab23b2a9e9a0a8d374d951531e1fb0b0` | Mature pipeline before the recent replay/reuse experiments |
| `baseline/benchmark-native-reuse-20260910` | `98b01b175eeed31f70d13e7cbf178b80bf07c9e0` | Frozen baseline for 1/2/3/10 seed-order comparisons |
| `baseline/golden-publication-20260908` | `bb0a7d7` | Earlier 1,839/1,851 Golden family-recovery and concrete-collection engine |

The `3e9a70a` baseline already includes older native acceleration. The older
`23b1e8b` tag preserves the implementation before that series. The `98b01b1`
baseline includes exact native reuse but precedes the adaptive-policy work.
These are preserved historical baselines, not a new claim that all versions
have interchangeable outputs or identical validation scopes.

## Old branch names

Every remote branch tip present before organization has an annotated tag named
`archive/branches-20260910/<old-branch-name>`. These tags were pushed before
superseded branches were removed. The table lists their replacements.

| Old name | Replacement |
| --- | --- |
| `experiment_slap_sweep_cut` | `paper/continuous-fragment-growth` |
| `experiment_adaptive_fragment_choices` | `dev/acceleration` |
| `experiment_conditioned_symmetry_reuse` | `stable/benchmark-baseline` |
| `experiment_dependency_fragment_repair` | Archived milestone at `f158864`, retained in development and paper ancestry |
| `experiment_incremental_cut_replay` | Archived milestone at `3a4d459`, retained in development and paper ancestry |
| `fable_native_clean_integration` | `stable/pre-acceleration` |
| `Fable_AAM_Opt` | `dev/native-engine` |
| `agent/native-index-chirality-internal-path` | `dev/native-index` |
| `aam_extreme_acceleration` | Archived milestone at `e77125b`, retained in the mature pipeline's ancestry |
| `organic_single_step_retro_syntehsis` | Baseline tag at `23b1e8b` (original misspelling preserved in its archive tag) |
| `stable_backup`, `repair/exact-aam-symmetry` | `stable/legacy-aam` |
| `prop-queue-align` | Archived reference at `fda540c`, also retained in `main` ancestry |

The local-only `fable_aam_catalog_integration` tip `93c35c2` and the older local
`Fable_AAM_Opt` tip `15e89fb` also have tags under
`archive/local-branches-20260910/`. Existing `pr14-works-no-dedup` is unchanged.
Unmerged native-engine and native-index commits are preserved on their own
development branches; they were not folded into another implementation.

To restore an old branch locally without changing current work:

```bash
git fetch origin --tags
git branch recovered/cut-replay archive/branches-20260910/experiment_incremental_cut_replay
```

Reference inventories and patches of existing uncommitted work were also saved
locally under `.git/branch-organization-20260910/` in both working copies.
Existing uncommitted research/HPC changes and raw benchmark folders were not
included in the manuscript commit or removed by branch cleanup.
