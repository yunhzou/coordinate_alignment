# Ground truth and detected mappings: 21-case comparison

Open `viewer.html`, or download `viewer.zip`, unzip it, and open its `viewer.html`. All SVG drawings and data are embedded: no server, external JavaScript, CDN, or network requests are required. The existing report server also serves `/golden_remaining_21/viewer.html` on port 8765.

Cases: 7, 19, 590, 603, 780, 833, 845, 865, 867, 871, 986, 1228, 1285, 1358, 1377, 1380, 1475, 1553, 1568, 1574, 1793.

## Reading the view

- Ground truth remains on the left; detected assignments remain on the right. Both show all input reactants and product(s), with consistent atom indices.
- Select a case, then one of up to three detected raw-heavy representatives. These are reference-blind ranked representatives pooled from the saved random- and distance-policy runs, not all assignments inside compressed families.
- Click atoms to link the same source atom across both mappings. Switch to target linking to inspect different donors for the same product atom. The exact correspondence table gives every target atom's ground-truth and detected assignment.
- Dashed red circles flag different literal reference pairs, not a proof of chemical inequivalence: symmetry-related assignments may use different raw indices.
- Reference colors are reference-derived conserved regions. Detected colors are actual recorded AAM fragment groups. Colors are reused across columns only for equal source-heavy-atom sets (H-only groups separately).
- Explicit H can be shown. Search and detected ranking include H, but reference H identities are unannotated. Comparable bond counts in the headings are heavy-atom-only counts.
- Unknown outcomes are labelled unknown. Case 1793 has an incomplete search: its representatives are ranked only within the first completed raw cut with terminal states from each run, explicitly disclosed in the interface. They are not presented as a completed-sweep ranking.

No matching searches were run for this viewer, and no reference-constrained diagnostic results were inserted. Core AAM is unchanged. `mappings.json` retains selected assignments, fragment groups, terminal IDs and archive provenance without the large drawings; `summary.json` retains case-level provenance. Full per-case drawing payloads remain at:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_remaining_viewer_20260908
```

## Verification

- Molecule atom elements and complete WBO matrices are asserted equal to the input endpoint order used for archived matching. Reversed searches are converted back to input R/P coordinates.
- Display ranking is tested against the existing benchmark ranker on all partial/full injective assignments of a three-atom test problem.
- Layout packing translates disconnected components only. Tests verify unchanged atom identity, molecular connectivity and within-component distances.
- Seven viewer-related tests passed.
- Headless Chromium checked all 21 cases, every candidate selector, explicit-H toggles and cross-panel atom selection. Zero JavaScript errors and zero external network requests. A 900-pixel viewport overflow check passed.
- `browser_check.json` records the UI checks; `browser_check.png` shows the resulting simultaneous comparison. The served HTML returned HTTP 200.

Generation code: `bench/view_golden_remaining.py`, `bench/golden_remaining_template.html`. Browser QA: `bench/check_golden_remaining_viewer.cjs`. Graph loading/ranking jobs used a 290-second execution watchdog and five-minute Slurm limit; these were visualization tasks, not AAM searches.
