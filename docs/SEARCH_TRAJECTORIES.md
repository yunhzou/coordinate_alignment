# Reusable search-trajectory pipeline

This pipeline animates the **algorithm's search history**: seed selection,
one-atom growth attempts, compatible and rejected targets, candidate branching,
deduplication counts, deferred boundaries, and committed fragments. R and P keep
their original coordinates throughout. It uses the shared white 3D reaction
viewer style.

The stages are independent:

1. Select saved archives and contexts in a JSON manifest.
2. Replay the selected fragment calls with their producing matcher and verify
   the placements, fragments, and deferred boundaries against the archive.
3. Store the recorded events and candidates in `trace.json`.
4. Render `trace.json` into a self-contained `algorithm_trajectory.html`.

No cut sweep is run. One recorded history is followed per selected terminal.
All live compressed candidate representatives and symmetry blocks are included;
the pipeline does not enumerate atom permutations. Candidate previews during
growth can change as refinement proceeds. The selected saved branch determines
which placement is committed before the next fragment starts.

## Build any saved case

From the repository root, using its Python environment:

```bash
python tools/build_search_trajectory.py manifest.json output_directory
```

A minimal manifest is:

```json
{
  "selections": [
    {"archive": "saved/aam.pkl.gz", "context": 0}
  ]
}
```

Archives can be AAM checkpoints, AAM JSON records, or gzipped AAM JSON records.
Paths resolve relative to the manifest. Endpoint coordinates must be present.
Multiple selections must refer to exactly the same R/P endpoints and atom order.
The matching tolerance, graph threshold, branch cap, seed order, cuts, and prior
committed state are read from the selected archive context.

Optional manifest fields:

| Field | Purpose |
| --- | --- |
| `title` | Page title |
| `engine_source` | Producing engine's `src` directory; defaults to this checkout |
| `engine_commit` | Human-supplied reference revision; runtime source hashes are recorded separately |
| `key_atoms` | Reactant atom indices labeled by default |
| `watch_targets` | Product atom indices whose occupancy is shown below the structures |
| `event_tolerance` | Bond-order event display threshold, default 0.5 |
| `selections[].terminals` | Specific terminals; defaults to all terminals in that context |
| `selections[].label` | Run button label |
| `selections[].path_labels` | Object mapping terminal IDs to display labels |
| `selections[].default_terminal` | Branch selected when opening this run |
| `selections[].focus` | Optional jump point: `{"source_edge": [a, b], "target": p}` |

Capture runs in a separate process because its diagnostic observers temporarily
wrap the matcher. The actual matching functions decide compatibility; a replay
that differs from the saved placement fails. The bundle records archive hashes,
runtime source hashes, and validation results. Native archives use the matching
Python trace implementation and must pass the same placement verification.

## Re-render an existing trace

```bash
PYTHONPATH=src python tools/render_search_trajectory.py output_directory/trace.json output_directory/algorithm_trajectory.html
```

This stage does not replay matching or load the original archives. The generated
HTML embeds its molecular renderer, shared stylesheet, and trace, so it works
offline. Playback supports frame stepping, scrubbing, candidate selection,
target-trial inspection, branch selection, focus/reset, and final bond overlays.
Final events use red losses on R and green gains on P. Black/gray dashed bonds
mark the tested edge, sweep cuts, or deferred search boundaries.

```bash
PYTHONPATH=src python tools/check_search_trajectory.py output_directory
```

The browser check needs Playwright and Chromium. `MANUSCRIPT_BROWSER` can select
an installed Chromium executable. It checks playback, candidates, fixed endpoint
coordinates, final event colors, offline operation, and mobile width.

## PR7 example

The committed example on `paper/continuous-fragment-growth` is
`reports/pr7_search_trajectory_20260911/manifest.json`. It selects both branches
of the one-seed R18–R21 cut and all three retained branches for order 15 of the
30-seed run with the same cut. The latter starts at O21 and includes the recovered
C–O-cleavage alternative. All 19 distinct fragment-call replays match their saved
placements exactly.

```bash
python tools/build_search_trajectory.py reports/pr7_search_trajectory_20260911/manifest.json reports/pr7_search_trajectory_20260911
```

Open `reports/pr7_search_trajectory_20260911/algorithm_trajectory.html`. It starts
at the P9 decision. **Replay from start** follows growth atom by atom; the
**O21 first** button opens the corresponding successful saved branch.

`tests/test_search_trajectory.py` also exercises a separate four-atom example,
non-default matching settings, JSON/checkpoint loading, input preservation,
manifest orchestration, render-only behavior, and rejection of mixed endpoints.
