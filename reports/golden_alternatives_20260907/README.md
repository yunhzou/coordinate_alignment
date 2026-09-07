# Inspect alternative atom assignments

Open `index.html` or download an individual `case*/viewer.html`; each case viewer is standalone, with embedded SVG and no external dependencies. The existing server exposes the index at http://localhost:8765/golden_alternatives_20260907/index.html after port forwarding.

| Case (zero-based) | Saved terminals | Raw heavy mappings | Chemical classes | Reference representative class rank |
|---|---:|---:|---:|---:|
| 9 | 591 | 34 | 23 | 2 |
| 15 | 334 | 26 | 20 | 1 |
| 21 | 2,319 | 47 | 47 | 4 |

All chemical classes of the saved representatives are displayed, not a top-k subset. Each class uses its best-ranked saved explicit-H witness and one actual fragment-growth path. The reference is a separate view; reference colors are conserved regions, not invented AAM fragments. Colors for search candidates come from their saved fragment path, converted back to original R/P indices when search direction was reversed. Within each drawing pair, identical r-labels identify the actual assigned atom. Color consistency across candidates refers to target regions, not necessarily the same donor atoms.

Click a plot point, table row or dropdown entry to inspect a class. Click atoms to highlight their mapped partners; use the explicit-H toggle to inspect hydrogen assignments. Green plot points mean the shown representative passes the reference chemical-equivalence check. Blue points are not a proof that the entire compressed family lacks a reference-compatible alternative.

## Why multiple classes remain

The search permits different bond edits and different atom sources; it does not require a unique chemically plausible mechanism. Repetition also occurs across seeds and cuts. The table separates saved terminal counts, distinct raw heavy relations and distinct explicit-atom witnesses. H assignments can differ even when the heavy relation is identical. Endpoint symmetry merges different raw heavy relations into a chemical class.

For case 9, the top class has 8 scored bond events and the reference-equivalent class has 9. Minimizing this graph-only event count does not necessarily select the reference. For case 21, classes 3 and 4 tie at 8 events; the reference is class 4, with the existing atom-index tie-break placing class 3 first. Case 15 provides a contrast: its reference-equivalent class is top-ranked with 4 events. Events include explicit H and count breaking, forming and bond-order changes separately.

These are witness-level comparisons, not proofs that each family has been optimally ranked over its compressed alternatives. No core search, chirality filter or ranking policy was changed. Existing archives were reused; generating the three cases took about four seconds total.

## Verification and provenance

Original archives: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_policy_full_20260907`.

Generation script: `bench/view_golden_alternatives.py`. `mapping.json` stores full mappings, selected terminals, actual fragment steps and colors; `classes.json` stores counts and scores. Drawing inputs are checked against archived elements and bond matrices. Every displayed AAM fragment partition and its mapped target atom set was checked, along with class terminal totals and JavaScript syntax. HTTP serving returned 200. A browser screenshot check could not run because the installed browser lacks `libnspr4.so`; visual rendering has not been independently confirmed in a browser here.
