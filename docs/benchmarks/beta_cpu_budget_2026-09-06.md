# Budgeted Slurm smoke test

The resource policy was checked on a two-source explicit-H methanol target.
This checks CPU admission, both query stages, complete assembly, and measured
CPU accounting. It is not a full-bank performance benchmark.

- CPU budget: **3**, including a one-CPU coordinator.
- At most one two-CPU worker on the mixed SMT partitions; scan and indexing
  share each allocation. Slurm allocation intervals verify the actual peak
  was three CPUs, including the coordinator.
- Coordinator job: 432868. Worker arrays: 432869 and 432871.
- Target coverage: **6/6 atoms**, methane plus water; no cap hits.
- Four worker attempts reported **5.24691 actual CPU-seconds** total versus
  **13.85181 reserved CPU-seconds during their process lifetimes**.
- Coordinator CPU: **0.58158 seconds**. Its measured elapsed time was 54.771 s,
  including queue waits. Very small tasks expose Slurm submission overhead;
  this test is not evidence of improved full-bank utilization or runtime.

The first submission (432860/432861) exposed conflicting inherited Slurm memory
environment variables. Those workers exited before computation. Nested
submission now clears the coordinator's memory mode. A second test (432863,
432864, 432866) completed chemically, but accounting showed two CPUs allocated
for each one-CPU worker request on SMT nodes, exceeding its requested budget.
The final test above uses topology-aware whole-core workers and counts the
coordinator's actual allocation, fixing that error.
Failed allocations are not included in the process-level CPU
numbers because no worker process started to write a resource report.

Artifacts: `/h/399/yunhengzou/coordinate_alignment/data/retro_runs/beta_budget_smt`.
Per-attempt `resources/*.json` distinguishes CPU work from reserved capacity.
The final `result.json` records the CPU ceiling and measured worker/coordinator
CPU time. Slurm accounting is required to include allocation startup, cleanup,
and any workers killed before writing reports.

Automated tests cover budget arithmetic at small and larger limits, inherited
memory isolation, checkpoint yielding/resume without recomputing saved matching,
and preservation of the existing disk-indexed beta recommendation behavior.
