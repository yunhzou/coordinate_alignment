# AAM paper benchmark: preparation and proposed protocol

Status: original data acquired and audited; three small full-search smoke tests completed. **No dataset accuracy or baseline speed comparison has been measured.** Core AAM and retrosynthesis code are unchanged.

## Datasets and baselines

- Primary: Lin et al., *Atom-to-atom Mapping: A Benchmarking Study of Popular Mapping Algorithms and Consensus Strategies*, Molecular Informatics, DOI [10.1002/minf.202100138](https://doi.org/10.1002/minf.202100138). The original Golden archive contains 1,851 reactions. Obtain it from the [authors' repository](https://github.com/Laboratoire-de-Chemoinformatique/Reaction_Data_Cleaning), pin its commit and file hashes, and preserve original records.
- Cross-paper comparison: LocalMapper's [2024 paper](https://www.nature.com/articles/s41467-024-46364-y) evaluates 1,758 Golden records after its own exclusions and discusses reference corrections. Reproduce its actual evaluated record list before comparing percentages. Do not compare our strict full-composition subset to their published headline accuracy.
- Difficulty breakdown: the Jaworski-derived 491-reaction subset discussed in that paper includes patent, typical, and complex reactions. It overlaps Golden: use it as a stratum, not an independent test dataset.
- Candidate extension: [SynRXN](https://www.nature.com/articles/s41597-026-07260-w) supplies additional benchmark tasks/biochemical sets. Audit provenance and overlaps before adopting it as independent validation.
- Primary rerun baselines: [RXNMapper](https://github.com/rxn4chemistry/rxnmapper), [LocalMapper](https://github.com/snu-micc/LocalMapper), and an available traditional graph-based mapper such as Indigo or ReactionDecoder. Pin code, model/checkpoint, dependencies, preprocessing, and licenses. Published baseline numbers are context, not timing measurements on our hardware.

## Input audit

Archive pinned to upstream commit `793475e54d8b2c7f714165a61e4eb439435d7d92`.
ZIP SHA256: `baf10464a6cd92dc4719f8becd5f001ab4b6bdafd363d74c5ae9dee77755f3d6`.
RDKit version: `2026.03.5`.

| Audit of unchanged raw RDF | Records |
|---|---:|
| Parsed and sanitized; both sides present | 1,851 |
| Equal full heavy-element composition | 1,015 |
| Equal full composition after adding explicit H | 733 |
| Duplicate nonzero atom-map labels on either side | 0 |
| Duplicate map-free canonical inputs under this reader | 0 |
| At least one explicitly mapped hydrogen in reference | 11 |

These are our reader's input audit counts, not the published exclusion counts or accuracy denominators. Full endpoint equality is much stricter than mapping product atoms to a subset of supplied reactants. Many legitimate benchmark records omit byproducts or include excess reagents. **Unequal full compositions do not make a reaction invalid.**

The current public `AAMProblem` rejects unequal endpoint counts/compositions. This is a public-interface limitation; it does not establish the capabilities of every lower-level matcher. Before claiming a full Golden benchmark, inspect/reuse the appropriate partial-output path and define how its unmatched atoms are represented. Do not bypass the guard, invent byproducts, supply answer-derived anchors, silently drop reagents, or substitute the beta retrosynthesis recommender. Report unsupported records explicitly if the adapter is not yet available.

Keep explicit H during our search. Primary comparable accuracy must use the annotated heavy-atom correspondence: Golden generally does not specify individual hydrogen identities. Do not invent hydrogen ground truth or report H-identity accuracy from automatically added hydrogens. Preserve original map labels only in the evaluator; strip them from AAM input. The audit preserves every row and records exceptions without repairing it.

## Separate search from selection and evaluation

1. **Input adapter:** unmapped reaction to full explicit-H graphs; preserve provenance, charges, isotopes, stereochemical annotations, reagents, and unmatched-atom policy. The smoke adapter uses element/bond-order graphs and zero placeholder coordinates solely because graph search does not use coordinates. It is not a stereo-aware or quantum-WBO benchmark adapter.
2. **Core AAM:** call the existing algorithm, retain its compressed graph and cap status, and save intermediates immediately. No explicit automorphism/bijection enumeration just for scoring.
3. **Output selection:** freeze the reference-blind ranking and deterministic tie policy before evaluation. Optional mechanism grouping remains outside core search and has its own timing.
4. **Evaluator:** compare mapping-induced reaction graphs with side-specific bond/atom labels and the documented stereo/isotope policy. Numeric map IDs and atom order are not chemical identities. Symmetry equivalence must preserve whole correlated assignments, not independently interchangeable atoms. Explicitly account for atoms appearing on only one recorded side.

Report separately:

- **Top-1 reaction accuracy:** the frozen selected answer agrees with the reference up to the declared chemical symmetry equivalence.
- **Compressed-family reference recovery:** a compatible reference assignment exists in the returned families after heavy-atom projection. This is an oracle diagnostic, not top-1 accuracy; count unresolved evaluator searches separately.
- Secondary atom accuracy and reaction-center precision/recall, with unmatched atoms and bond-order changes defined explicitly.
- Unsupported inputs, errors, timeouts, cap hits, and any incomplete searches in the full denominator. Also report conditional accuracy on supported inputs, labelled as such.
- Raw reference agreement and any expert-adjudicated corrections as separate tables. Fewer bond changes alone is not evidence of chemically correct mapping.

## Performance and ablations

Use identical hardware and CPU budgets for CPU comparisons; separate single-reaction latency from batched throughput and GPU results. Pin software versions and record node/CPU model. Count initialization, parsing/graph preparation, AAM growth, symmetry finalization/merge, optional grouping/selection, scoring, and serialization separately; report end-to-end time too. Report median, p95, maximum, CPU-seconds, peak RSS, cap hits and timeouts, not just the fastest run. Persist each reaction so evaluation can be repeated without remapping.

Initial pilot proposal: a fixed, recorded stratified sample spanning size, balance status, rearrangements, stereo and symmetry, followed by the entire pinned set. Keep tuning cases separate from the final test; audit dataset overlap. Compare branch caps 100/500, three/ten seeds, sweep on/off, Python/native backends where semantics agree, and input atom/component ordering. These are proposed experiments, not measured results. Do not change the search defaults based on test labels.

Use a 300-second per-reaction worker watchdog and a 10-minute supervisory watchdog. Timeouts are results, not deleted records. Start with a modest measured CPU budget; scale independent reactions through Slurm after the pilot establishes RAM/runtime needs.

## Smoke test actually run

Three smallest explicit-H-balanced records by atom count then original zero-based index; chosen for adapter validation, not representative speed or accuracy. Configuration: existing `AAMSearchConfig()` defaults, native backend, one CPU worker, three seeds, cap 100, tolerance 1.0, full one-edge sweep. No mechanism grouping or scoring. Full compressed AAM records and per-cut intermediates saved.

| Original index | Explicit atoms | Sweep searches | Search + intermediate save (s) | Cap stops |
|---|---:|---:|---:|---:|
| 1001 | 9 | 8 | 0.094 | 0 |
| 979 | 10 | 9 | 0.191 | 0 |
| 267 | 12 | 13 | 0.142 | 0 |

These were lightweight smoke checks on `bosque1`, not a Slurm performance campaign. Core implementation commit: `eff2f41` (new work in this turn adds only benchmark preparation/smoke tools and documentation). Do not extrapolate them to the full dataset or claim speedup/accuracy from them.

Local data: `/h/399/yunhengzou/coordinate_alignment/data/aam_benchmarks/golden_original_20260906/`.
Saved search results: `/h/399/yunhengzou/coordinate_alignment/data/aam_benchmarks/golden_smoke_20260906/`.

Reproduction: `bench/prepare_golden_benchmark.py --output NEW_DIRECTORY --commit 793475e54d8b2c7f714165a61e4eb439435d7d92`, then `bench/golden_aam_smoke.py --audit NEW_DIRECTORY/audit.jsonl --output NEW_SMOKE_DIRECTORY` with `PYTHONPATH=src`, `RXN_CORE_NATIVE=1` and BLAS/OMP threads set to 1. Both refuse to overwrite an existing run directory.
