# Big-block / gap-first beta

## Slurm CPU budget

`bench/run_beta_distributed.py --cpu-budget N` now limits simultaneous allocated
CPUs for one run, including its coordinator. Slurm admits independent workers
as capacity becomes available, up to that limit; there is no requirement to
obtain a fixed number of nodes. Use one coordinator per run directory.

The default requests the smallest whole-core worker. Slurm topology is queried
and the requested worker size is rounded up to the largest threads-per-core
value in the eligible partitions. On this cluster's `CR_CORE_MEMORY` policy,
newer nodes reserve two logical CPUs per core; workers request and use both.
`--worker-cpus` explicitly changes the requested size; concurrency is
`floor((budget - actual coordinator allocation) / rounded worker CPUs)`.
For example, a 32-CPU budget with a one-CPU coordinator allows 15 two-CPU
workers on the mixed partitions, or 31 one-CPU workers with
`--partitions cpunodes`. Smaller tasks can fill partially occupied nodes.
Memory is requested per CPU.
The budget applies to this run, not unrelated jobs belonging to the user.

Workers run connected scanning and index construction in the **same allocation**.
They release the allocation when finished. Arrays from different query stages
never overlap. This avoids retaining a large idle CPU allocation during serial
refinement or assembly. Workers are not kept alive between separate gap queries;
that still incurs submission latency.

Each worker has a ten-minute Slurm limit and a 580-second process watchdog.
Scanning stops admitting new sources after 75% of the process budget, leaving
room for in-flight work and indexing. A voluntary checkpoint yield (exit 75)
resubmits only unfinished shards and reuses completed source evidence. Real
failures and hard timeouts are reported rather than silently retried. Pending
queue time is not subject to the computation watchdog. The coordinator's own
Slurm wall limit should allow queue time; request only one coordinator CPU.

For a new run, also supply `--catalog`, `--target-smiles`, and `--shards`.
Shards control workload granularity, independently of the CPU ceiling. Existing
run manifests retain their shard layout and checkpoints. Resume with a different
CPU budget after the previous coordinator and its workers have exited.

```bash
PYTHONPATH=src .venv/bin/python bench/run_beta_distributed.py \
  --run /path/to/new-run --catalog /path/to/bank.csv.gz \
  --target-smiles CO --shards 256 --cpu-budget 32
```

`query_*/resources/*.json` records user/system CPU seconds and process wall time
for each scan/index attempt. The final report distinguishes actual CPU work
from CPU capacity reserved during worker process lifetimes. These measurements
exclude Slurm startup/cleanup; hard-killed workers may lack reports. Use Slurm
accounting to include those allocations. The coordinator is measured separately.

This is a **separate, opt-in recommendation workflow**. The full augmented
bank scan and exact assembly remain available and unchanged. Neither workflow
silently invokes the other on failure.

```text
Bank + explicit-H target
    -> connected-fragment scan only (no bank-wide augmentation)
    -> prioritize large target blocks, retaining alternative choices
    -> full augmented AAM for the selected reactant
    -> use that reactant's additional fragments to cover gaps
    -> if gaps remain: query the bank against the missing target subgraph
    -> add/refine another reactant copy; repeat until full coverage
```

## Boundaries

- `fragment_matching.connected.find_connected_fragments` exposes the existing
  initial matcher independently. It retains compressed symmetry and search
  paths. It does not change the core AAM or run augmentation.
- `retrosynthesis.beta.FragmentQueryBank` caches connected queries by requested
  target region and full augmented results by selected source ID. It considers
  all discovered connected families, ordered by the search policy; “largest”
  is not a proven maximum common subgraph.
- `retrosynthesis.beta.recommend_big_blocks` explores provisional assemblies
  best-first: greater coverage, fewer fragment units, fewer distinct species,
  then direct explicit-atom retention. This is **heuristic search order**, not
  the exact full workflow's global ranking or symmetry-adjusted retention.
- `tools/search_retro_beta.py` is the separate entry point. It saves every typed
  detection checkpoint before occupation expansion, an event log, the input
  manifest, a typed final result, and human-readable original-index mappings.

## What is preserved

Explicit H, tolerance 1.0, cap 100, no sweep, and no hidden candidate/beam/
reactant-count cap. All discovered sizes and correlated alternative placements
remain available for backtracking. Repeated reactant copies and overlapping
target coverage are allowed. A refinement replaces the **whole** selected
copy's assignment; it never joins incompatible atomwise assignments from
different AAM alternatives.

Connected gap searches use an induced missing-target subgraph. The local-to-
original target index map is stored explicitly beside the untouched matching
evidence. An augmented replacement must preserve the provisional block's
target coverage, though its internal atom mapping can change. All selected
copies must receive augmentation before a full cover is returned.

## Beta limitations

This deliberately trades exhaustive discovery/ranking for directed proposals.
A gap-only query lacks the already-covered neighborhood and can miss a useful
overlapping *initial* anchor; full selected-R refinement can add overlap.
Alternative branches remain available, but arbitrary fragment recombination
within one source copy is not permitted. Failure is reported as a partial
cover, not a proof that the bank has no solution. Cap hits remain explicit.
The default collects 20 complete assemblies across four construction patterns.
Each pattern preserves whole source-copy assignments and their matched fragment
partition, modulo target symmetry. One representative per chosen pattern is
reserved; all displayed results are then sorted by the final ranking. These are ranked discovered proposals, not
a certified global top 20. Every discovered complete assembly is checkpointed
by the distributed runner before collection continues.

Final beta ranking is independent of discovery order. It uses Pareto layers
over two objectives: maximize explicit-H retention and minimize source boundary
cuts plus unsupported target connections. A candidate dominates another only
if it is at least as good on both and strictly better on at least one. Layer 1
contains the nondominated alternatives; subsequent layers repeat this operation
after removing earlier layers. At identical objective values, fewer total matched
fragments across all source copies wins, then fewer distinct species. Within a
layer, display order is deterministic but is not a preference between trade-offs.

Retention is unique covered target atoms divided by all input atoms across
physical source copies. Overlaps cannot inflate it. There is no separate penalty
for copy count, fitted weighting, or ground-truth bonus. Cuts/connections are
geometric support metrics, not reaction edit counts. Exact rational retention
is used internally; a prefix-max tree computes 2D layers in O(n log n).
Pattern/output quotas may limit displayed alternatives; the saved rerank report
counts nondominated pool and displayed entries. The full workflow's exact solver
and its ranking bounds are unchanged by this beta-specific change. Connected
discovery and the known-supplier join retain their fragment-oriented traversal
heuristics; these are not final Pareto ranking or global-optimality guarantees.

`bench/rerank_beta_saved.py --run RUN` applies this same ranking to all saved
complete beta candidates without any new matching. Optional
`--validation-stem ground_truth_assemblies` compares validation sets but never
adds them to the blind candidate pool.
Outputs are saved as `assemblies_pareto.json` and `assemblies_pareto.pkl.gz`,
preserving earlier rankings. Rebuild the viewer with `--result-stem assemblies_pareto`.
`--scan-summary RUN/assemblies_viewer.json` explicitly reuses the same run's
previously tallied matched-source count, avoiding a bank-wide metadata reread.
Explicitly requesting one recommendation stops at the first fully refined cover; it does
not prove that cover is globally optimal or chemically feasible.

Per-R connected queries can run in parallel with `--workers`. Only that many
queries are in flight at once. Selected-source augmentation is currently
serial and can still be expensive. No full-bank beta speed claim has been
established. The CLI currently loads the bank graphs in memory; benchmark
memory and runtime before scaling this prototype to the full catalog.

## Run

Input CSV (plain or gzip) has `Bank ID` and `SMILES` columns. Example:

```bash
PYTHONPATH=src timeout --kill-after=10s 580s .venv/bin/python tools/search_retro_beta.py \
  --target-smiles CO --catalog /path/to/bank.csv \
  --output-dir /path/to/new-beta-run --workers 4 --recommendations 20 --patterns 4
```

Use a new output directory. On Slurm, request the same CPU count as `--workers`
and a ten-minute job limit. Partial checkpoints survive a watchdog termination;
pickle artifacts are trusted local data and must not be loaded from untrusted
sources. Automatic checkpoint replay is not implemented in this beta entry
point yet. The full workflow's own checkpoint/recovery tools are unchanged.

## Distributed full-bank execution

`bench/scan_beta_catalog.py` / `hpc/scan_beta_catalog.sbatch` shard each connected
target query across nodes. Each source saves typed connected evidence before
occupation projection, then a placement archive and completion record.
`bench/index_beta_query.py` builds one SQLite proposal index per shard.

`bench/run_beta_distributed.py` runs the same beta selection policy against
these indexes. Sorted proposal streams feed the best-first frontier lazily:
only the next alternative from each expansion needs to be on the heap. No
alternative is pruned and full source AAM objects stay in their archives until
needed. This is a memory/execution change, not a new matching rule or beam.
Only selected source graphs are instantiated in the coordinator.

The distributed runner reuses completed source archives, query indexes, and
selected-source AAM checkpoints after interruption. Each query/index Slurm job
has a ten-minute watchdog. Node-start failures must be distinguished from AAM
timeouts; do not rerun completed matching to recover an index-only failure.

## Known-supplier validation and assembled viewer

`retrosynthesis.beta_assembly.assemble_supplier_copies` joins the allowed whole
occupations of specified reactant copies. Its lazy branch-and-bound traversal
returns complete covers in fragment-count order. Overlap and repeated copies
are allowed; mutually exclusive placements of one copy cannot be combined.
This component is independent of how the occupations were detected.

`bench/assemble_beta_known_sources.py --run RUN --case CASE` reads the known
suppliers' actual full-bank connected checkpoints, applies selected-source
augmentation, and saves four distinct full-cover patterns. This validates
recoverability, not the supplier set's blind recommendation rank.

`bench/view_beta_result.py --run RUN --validation-stem ground_truth_assemblies
--case CASE` renders saved blind assemblies and separately labelled validation
panels in `RUN/assemblies.html`. All panels use the actual selected R-to-P
mappings, include hydrogen, and show a combined target-only assembly.

The viewer shows separate score cards for Pareto rank, bond breaking (source
boundary cuts), bond forming (unsupported target connections), their total,
explicit-H target coverage, and reactant retention. These bond counts are
geometric proxies, not validated reaction events. The sort selector orders the
display by any score without changing Pareto ranks or mixing validation panels
into the blind recommendations. Target coverage uses a union, never a sum of
overlapping claims.

A clickable retention-versus-structural-changes plot sits above the molecular
views. Blue points are blind proposals; green points are validation sets.
Exactly coincident scores share one point labelled with their multiplicity,
without artificial jitter. Clicking opens a saved assembly and exposes buttons
for every coincident mapping, ordered by fragment count then species. Selected
markers and molecular views stay synchronized. No matching is run by the plot.
