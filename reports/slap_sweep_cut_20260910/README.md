# SLAP plus blind single-edge sweep

Status: implementation and validation; no 99% claim.

The previous AAM optimization campaign was stopped at the user's request.
Slurm array `466079` and dependent comparison `466114` were cancelled.
Existing AAM checkpoints remain under
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/shared_capacity_full_20260910`.

This separate experiment tests pinned upstream SLAP
(`ea248fd9494f52f4865193e87a98cc92c62b5f9e`) on all 1,851 Golden reactions.
Run uncut plus every single source-edge deletion, independently in both
directions and both native binary/weighted modes. Explicit H and native
element-balancing placeholders are retained. No atom is deleted. Rebuild graph
initialization after each cut; otherwise WL caches would describe the wrong graph.

Cuts perturb SLAP's input graph/objective, not its internal algorithm. This is
not identical to AAM's fragment cut constraint. Restore output labels onto the
original molecular endpoints before evaluation. Keep every returned native
alternative, label pool and LAP record. Reference mappings never select cuts,
order searches or stop cases. Native SLAP pruning is unchanged.

Report uncut versus sweep recovery and union with the saved expanded-SLAP
baseline (1,691/1,851 = 91.36%) separately. Use the existing strict heavy-atom
relation / endpoint-chemical-symmetry evaluator, including unmatched atoms.
Explicit H participates in search; reference-H assignment is not this metric.
Failures and unattempted cuts remain visible with the full 1,851 denominator.

Save mapping CPU/wall time separately from graph preparation, output export,
encoding, persistence and evaluation. Fixed hash/random seeds, single-threaded
workers, 300-second process watchdog per case/direction/mode (across its cuts),
and cluster-parallel independent cases. Partial cut outputs survive interruption.
`native.bin` contains independent gzip/pickle frames; committed JSONL records
identify each frame's offset, length and SHA256. Only load these trusted archives.

The home filesystem is full, so development is isolated on branch
`experiment_slap_sweep_cut` at
`/project/yunhengzou/coordinate_alignment/slap_sweep_worktree_20260910`.
No AAM core algorithm change is involved.
