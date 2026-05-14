# Prepared Step Folder

This folder is the same `pr1.tempo_ts3` example converted into the cached
step layout consumed by `load_step_inputs(...)` and `load_ts_targets(...)`.

```text
prepared_steps/pr1.tempo_ts3/
  R/                    reactant endpoint XYZ + wbo
  P/                    product endpoint XYZ + wbo
  sp_iter1/             IG target XYZ + wbo
  hess_iter1/           IG target XYZ + g98.out
  sp_groundtruth/       GT target XYZ + wbo
  hess_groundtruth/     GT target XYZ + g98.out
```

It is runnable in `cache-only` mode; no xtb execution is needed to load the
stored example.
