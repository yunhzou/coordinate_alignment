# Case 123: branch-cap ablation

`pr9.carbene.rearr_ts41a-endo`, 40 explicit atoms. Only the branch cap was changed
from 100 to 200. Identical input bytes and frozen search engine; tolerance 0.5,
10 seeds, 47 cut contexts per direction (470 seed/cut attempts), explicit H.
No reference mappings or previous predictions were supplied.

| Cap | Direction | Full terminal witnesses | Search elapsed incl. checkpoint I/O | Search CPU excl. checkpoint I/O |
|---|---|---:|---:|---:|
| 100 | R→P | 0 | 0.333 s | 2.355 s |
| 100 | P→R | 0 | 0.338 s | 2.681 s |
| 200 | R→P | 12,258 | 5.578 s | 20.977 s |
| 200 | P→R | 985 | 0.441 s | 3.155 s |

At cap 100, all 470 attempts per direction stopped during fragment growth and
no completed terminals existed. At cap 200, both directions returned complete
element-preserving mappings including H. Some other branches still hit caps.
This demonstrates that a cap-hit flag does not imply absence of a returned mapping.
Terminal witnesses are not globally deduplicated mechanism/pattern counts.

The main fixed-cap-100 benchmark remains **139/140**. With this explicitly targeted
cap-200 follow-up, every one of the 140 cases has now returned a full AAM mapping.
This is not a fresh uniform cap-200 benchmark and is not an accuracy measurement.

Two 16-CPU tasks, 32 GB each, five-minute process watchdog, Slurm array 454401.
Node/model information is in each `environment.json`; this is not a same-hardware
speedup experiment. Actual search code and input files were reused unchanged.
Full archives and raw cut checkpoints:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_case123_cap200_20260908
```

Saved full terminal mappings are included here for inspection, and the original
compressed archives remain on the cluster. `comparison.json` contains stop-stage
counts and timing, with no inferred chemical correctness labels.
