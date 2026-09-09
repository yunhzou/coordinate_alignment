# Elementary-step holdout: native-input feasibility

Complete: 140 selected elementary steps, 18–137 explicit atoms per endpoint.
Every R/P pair has identical elemental composition. No ground-truth atom mapping
is available. **These are feasibility results, not accuracy measurements.**
No previous AAM predictions or reference TS structures entered the searches.

| Method | Completed calls | Cases with a full element-preserving assignment | Mean CPU seconds/call | Mean elapsed seconds/call |
|---|---:|---:|---:|---:|
| Our AAM, R→P | 140/140 | 139/140 | 55.555 | 7.455 |
| Our AAM, P→R | 140/140 | 139/140 | 72.543 | 10.910 |
| SLAPMapper, native XYZ | 140/140 | 140/140 | 0.07766 | 0.08214 |

The bidirectional AAM union covers the same 139 cases. Every reported full
assignment includes hydrogen, not just heavy atoms. This proves only that a full
element-preserving mapping was returned or is represented by the native label
groups; it does not establish a plausible mechanism or correct atom provenance.
Expert review is required before interpreting chemical quality.

## Native representations and budgets

- AAM uses the existing continuous WBO matrices and corresponding XYZ coordinates.
  Cached coordinate/element order was checked against the selected input, and
  input files were snapshotted and hashed. No new electronic-structure calculations.
- AAM: isomorphism tolerance 0.5, 10 seeds per cut, branch cap 100, existing sweep
  cut and symmetry handling, root seed 42, Python hash seed 0; 16 CPUs per search.
  The two directions run independently. Both endpoints have equal atom counts,
  so the smaller-first convention resolves to R→P.
- SLAPMapper 1.0.0 uses its native `map_3d` interface on original component XYZ
  files, native `bond_scale=1.2`, binary adjacency, and `break_sym='heavy'`.
  Hydrogens remain in the graphs. ASE 3.26.0 was installed in a run-local
  dependency directory without changing the prior competitor environment.
  One compute thread per call. No core mapper modifications.
- SLAP's native XYZ interface infers adjacency from geometry/covalent radii;
  it does not consume our continuous WBO values. Thus the timing table compares
  native workflows with different graph representations and search budgets,
  not interchangeable implementations of one objective or matched-accuracy speed.

## Caps, failures and alternatives

Follow-up: [case 123 at cap 200](case123_cap200/README.md) returns full mappings
in both directions. With this targeted follow-up, all 140 cases have a returned
full AAM mapping. The fixed-cap-100 table above is intentionally unchanged.

Case **123**, `pr9.carbene.rearr_ts41a-endo` (40 atoms), hits cap 100 in both AAM
directions with zero completed terminals. This is a capped search outcome, not
a claim that mapping the reaction is impossible. Its compressed graph, stops,
metrics and input remain available for diagnosis. No cap relaxation or targeted
retry was folded into this fixed-policy result.

AAM reports cap hits on 77 R→P cases and 82 P→R cases; most still return full
mappings. Caps and returned mappings are separate facts. There were **no compute
watchdog timeouts** or mapping exceptions in the final native calls. SLAP's own
internal search guards are unchanged; absence of an exported cap flag does not
prove exhaustive enumeration.

SLAP saved 176 native candidate families across 140 cases. AAM saved 816,600 R→P
and 1,046,543 P→R terminal witnesses, with the compressed search archives retained.
These terminal totals are not globally deduplicated chemical-pattern counts and
are not directly comparable to SLAP's candidate-family count. Neither output set
is claimed to contain every mathematically possible mapping.

## Timing and scheduler accounting

AAM CPU excludes measured checkpoint persistence and loading. Its elapsed column
includes search/checkpoint I/O, but excludes later terminal-witness export.
SLAP timing includes native geometry reading, graph construction and mapping;
startup and final result serialization are outside its mapping timer. Evaluation
and reporting are excluded. Full per-case fields retain these distinctions.

CPU totals: 7,777.679 seconds R→P; 10,156.072 seconds P→R; 10.872 seconds SLAP.
AAM elapsed medians are 3.214 and 3.629 seconds, with maxima 81.317 and 159.883
seconds. SLAP elapsed median is 0.0711 seconds and maximum 0.3137 seconds.
Do not interpret summed worker elapsed times as campaign wall time or CPU time.

Five-minute process watchdog; ten-minute Slurm job limit. Original arrays:
452980 (AAM), 452981 (SLAP XYZ). Twenty-two unstarted node-setup-stalled tasks
were relocated to available nodes; two needed a second relocation. Completed
calculations were not repeated. All actions and exact commands are in
`relocation.json`. Scheduler delay is not mapping compute. The preflight AAM
case 59 R→P and SLAP case 0 were reused, with original outcomes preserved.

## Query and expert-review artifacts

- `feasibility_per_case.json`: case IDs, per-direction outcomes, cap flags and timing.
- `feasibility_summary.json`: aggregate metrics, explicitly null accuracy.
- `native_slap_outputs.tar.gz`: all native mapping strings, compressed label
  groups, original-index graphs and costs. No arbitrary relabeling for display.
- `manifest.json`: historical input hashes and frozen engine identity;
  `report_provenance.json`: current report-driver hash.
- `slurm_native_accounting.psv`: allocation and step resource records. Avoid
  double counting allocation rows and their `.batch` children.

Full artifacts (including every AAM archive and terminal witness):

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_feasibility_20260908
```

Case lookup: `inputs/<index>/input.json`; AAM:
`directions/<index>/<R_to_P|P_to_R>/cuts/aam.pkl.gz` and
`terminal_mappings.jsonl.gz`; SLAP: `slap_xyz/<index>.json`.
Generate reports from saved data using `bench/elementary_feasibility.py report`;
there is no need to rerun searches for expert inspection.

The initial unrequested MOL/SMILES-conversion runs were cancelled and excluded.
Their partial artifacts survive only as an audit trail; see
`WITHDRAWN_SMILES_RUNS.md`. The original manifest records the historical launch,
not approval of those withdrawn protocols. The active driver now supports only
native inputs for this holdout. **21 regression tests passed.**
