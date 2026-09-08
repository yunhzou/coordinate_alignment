# Golden: full-node continuation

Continue the 16 unresolved cases from the ten-seed cap2000 campaign. Keep the same AAM configuration, default orientation, tolerance1.0, explicit hydrogens and single-edge cut sweep. No new seeds, reverse-direction rescue or core AAM changes.

Final result: **zero new recoveries, eight confirmed misses, eight unresolved**. Confirmed misses: 603,845,865,867,1033,1228,1377,1574. Unresolved: 19,590,780,833,871,1358,1380,1568. Combined demonstrated coverage remains 1,827/1,851 (98.70%). Across the originally requested 32 cases, the accumulated evidence is nine recoveries, 15 high-cap confirmed misses, and eight unresolved. A confirmed miss is limited to this sampled search, not all possible mappings.

## Execution

- Slurm array442547 requested 16 exclusive nodes, with 48 worker CPUs per case and all available node memory. The partition advertises 80 logical CPUs per node; exclusive allocation reserves those nodes, while this driver uses up to 48 processes as configured.
- Six completed AAM searches were reused without searching again. Ten incomplete searches resumed using the existing checkpoint-identity validation.
- Independent cuts are verified in parallel by the existing full-domain verifier. No global top-one ranking is computed. Witness terminal/path IDs are local to the named cut archive.
- A positive cut establishes recovery. A negative conclusion requires every expected cut to finish verification negatively. Missing cuts or unresolved verification remain unknown.
- Search watchdog300seconds and verification watchdog240seconds; overall job limit10minutes.
- Existing immutable checkpoint files are hard-linked into the continuation directory; new writes use atomic replacement and do not alter the previous campaign's files.

The six prior verification timeouts (603,845,865,867,1033,1574) were fully checked in 37–77seconds each, without rerunning their searches. All six were negative. Case833 exhausted its node's 184000MiB allocation before saving a completed cut and remains unknown; it was not blindly rerun.

## Saved evidence

Full inputs, frozen source, saved cuts/checkpoints, per-cut verification outputs and logs:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_fullnode_cap2000_20260907`

Local `summary.json`, `cases.json`, `witnesses.json` and `slurm_accounting.psv` record outcomes and resource use. The combined percentage is a diagnostic union with previous saved recoveries, not a new whole-dataset benchmark. Original unresolved case1793 is outside this continuation.

143 regression tests passed. The driver also handles an empty saved-cut set as an unresolved result without launching an empty process pool; that small cleanup was made after the frozen campaign snapshot.
