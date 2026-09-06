# Eight case studies

Download this directory and open `index.html`, then choose a target. The viewers contain assembled R/P mappings, construction patterns, scores, and separately labelled reference-set checks.

Bank: 155,305 structures in merged_fast_delivery_with_inventory. Current beta workflow, explicit hydrogen, tolerance 1.0, branch cap 100, no sweep. Seven fresh connected bank scans; case 5 reused its saved matching and reran assembly with the corrected reference isomer. Recommendations are not exhaustive or certified globally optimal.

| Case | Blind full-cover assemblies | Reference full-cover assemblies | Capped searches: blind / reference |
|---|---:|---:|---:|
| 1: fluorinated phosphorus ligand | 0: memory failure | 1 | 1 / 0 |
| 2: chloro quinoline alcohol | 20 | 0 | 5 / 0 |
| 3: aryl magnesium bromide | 20 | 1 | 2 / 0 |
| 4: dimethyl azobenzene | 20 | 4 | 0 / 0 |
| 5: cyano thiazole amidine | 20 | 4 | 1 / 0 |
| 6: phenyl cyclooctadiene | 20 | 4 | 18 / 0 |
| 7: acetophenone | 20 | 4 | 22 / 1 |
| 8: cyclohexene | 20 | 1 | 0 / 0 |

Case 1 completed its bank scan, but subsequent selected-result processing exhausted the coordinator's 12 GB allocation. Its viewer contains only an independent 78/78-atom reference assembly. This is **not a blind discovery**. Case 2 did not recover a compatible complete assembly from its specified reference set; that does not prove impossibility. Reference validation is separate from blind ranking throughout.

All displayed selected assemblies passed an independent audit of complete explicit-H coverage, per-copy injectivity, element identity, retained-fragment partitions, and preserved bonds within tolerance. This establishes mapping consistency, not chemical feasibility or a synthetic route.

Resource ceiling: 256 Slurm CPUs per case, 2 coordinator CPUs, workers of 8 CPUs with a 10-minute limit. Multiple node-startup stalls required checkpoint-preserving resubmission. `summary.json` records elapsed time including these delays and restarts; it is **not pure compute time**. Case 5's cached run is not a fresh-scan timing comparison.

Full checkpoints and execution logs remain at:
`/project/yunhengzou/coordinate_alignment/retro_runs/eight_beta_20260906`

Regenerate/audit this bundle with `bench/publish_beta_suite.py`. No core AAM changes were made for this batch. The supplier-copy join gained only a safe remaining-coverage capacity bound.
