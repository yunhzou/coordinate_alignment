# Viewer styles

There are **two shared viewer styles**. New reports provide data or a layout
adapter; they do not introduce another stylesheet or benchmark-specific skin.

| Style | Purpose | Authoritative files |
| --- | --- | --- |
| `reaction` | R/P comparison, TS structures/modes, mapping diagnostics and recorded growth | `src/rxn_core/static/reaction_viewer.css`; original R/P/TS layout in `reaction_viewer.html` |
| `catalog` | Catalog search, precursor/results lists and batch navigation | `src/rxn_core/static/catalog_viewer.css` |

The molecular presentation restores the original white panels, yellow mechanism
buttons, thin molecular sticks, R-order/spatial-alignment controls and TS mode
animation from `stable/legacy-aam:src/rxn_core/pipeline.py`. The separate dark
`elementary_comparison_viewer.html` skin and its derived missing-pattern viewer
were rejected and removed. Its reports now use the shared original renderer.

`rxn_core.viewers.reaction_html` renders one reaction; `collection_html` adds a
case selector with one shared copy of the renderer and libraries. Existing
diagnostic layouts use `style_document` / `viewer_style`. Golden 2D diagrams
remain 2D where the evidence has no native geometry; they use the same reaction
style. Catalog layouts retain their existing result-navigation behavior.
Different layouts and chemical highlight colors are not independent skins.

Bond-event overlays use **red for breaking or weakening** and **green for forming
or strengthening**, always in the R → P direction. Do not add separate colors for
bond-order changes. In the R/P comparison, **R shows only red changes** and
**P shows only green changes**; each event is highlighted once across the two
panels. Event names and numerical WBO changes remain in the details.

## Saved comparisons

`tools/render_mapping_comparison.py INPUT.json OUTPUT.html` renders saved full
mapping records. `bench/build_elementary_comparison_data.py` exports the older
elementary-step comparison records. `bench/missing_pattern_display_data.py`
exports and independently rescores the nine witnesses for cases 11, 64 and 101.
None of these commands starts a benchmark search.

The active missing-pattern comparison is
`reports/holdout_missing_pattern_seeds_20260910/viewer.html`.
Product fitting is a rigid display transform of the **recorded mapping**;
switching to native P restores its original coordinates. Explicit WBO connectivity
is used when provided. No search result is replaced by display fitting.

## Maintenance and validation

- Edit a shared stylesheet, then refresh affected generated HTML. Keep benchmark
  inputs, mappings, scores and timing records unchanged.
- `tools/sync_viewer_styles.py` migrates historical HTML using the source hashes
  recorded in `viewer-style-migration.json`; it verifies script/data preservation.
- `viewer-style-output-migration.json` records generated-page changes.
- `viewer-bundle-migration.json` records presentation-only archive changes; all
  other archive members were verified byte-for-byte.
- `tools/check_viewer_consolidation.py` checks the nine saved witnesses, native
  and aligned geometry, representative Golden/catalog/growth layouts and offline
  operation. Results and review images are under `docs/viewer-validation/`.
- Frozen stable branches are retained as historical algorithm baselines. This
  cleanup does not rewrite Git history or alter those baselines.

Static manuscript figures, molecular highlight colors and vendor-library internals
are not separate viewer styles. Standalone report pages embed the shared CSS so
they still work offline.

Validation on 11 September 2026: 20 relevant unit tests passed on the paper and
acceleration checkouts (CLI artifacts, catalog and Golden viewers, shared style
policy and proper rigid fitting). Browser checks passed for all nine saved
witnesses, mobile layout, offline use, and representative diagnostic/catalog/
recorded-growth pages. The tracked HTML audit found no private styles.
