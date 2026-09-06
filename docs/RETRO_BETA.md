# Big-block / gap-first beta

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
Requesting one recommendation stops at the first fully refined cover; it does
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
  --output-dir /path/to/new-beta-run --workers 4 --recommendations 1
```

Use a new output directory. On Slurm, request the same CPU count as `--workers`
and a ten-minute job limit. Partial checkpoints survive a watchdog termination;
pickle artifacts are trusted local data and must not be loaded from untrusted
sources. Automatic checkpoint replay is not implemented in this beta entry
point yet. The full workflow's own checkpoint/recovery tools are unchanged.
