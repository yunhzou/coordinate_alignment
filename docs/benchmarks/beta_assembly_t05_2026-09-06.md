# Example 5: complete beta assemblies

Run: `/project/yunhengzou/coordinate_alignment/retro_runs/beta_full_t05_20260906`.
Bank: `merged_fast_delivery_with_inventory`, 155,305 unique structures.
Matching settings: explicit H, tolerance 1.0, branch cap 100, no sweep.

- Replaying the saved full-bank queries and selected-source augmentation
  produced 20 complete assemblies across four construction patterns in
  **31.999 seconds** on a one-CPU coordinator (29.647 CPU seconds).
  No new bank query jobs were submitted. This is cached assembly time,
  not a fresh full-bank scan or HTML rendering time.
- Every returned assembly covers all 58 target atoms, including H.
- The known three-supplier set is **not in those blind top-20 proposals**.
  Separately joining its actual bank-connected matches after selected-source
  augmentation recovered four complete patterns with **6, 6, 6, and 7 fragments**.
  This validation took **12.894 seconds**, including three new selected-source
  augmentations. One augmentation hit its branch cap; no exhaustive matching
  or global optimality claim is made.

Validation suppliers: `INVENTORY-000283` (acetylacetone), `INVENTORY-000436`
(2,6-diisopropylaniline), and `MCULE-7889053722`
(2-amino-6-cyanobenzothiazole). Their full-bank connected occupation counts
were 1, 1, and 6; selected-source refinement yielded 22, 11, and 2,157
whole-copy occupations. These were joined without mixing incompatible
alternatives of a single copy. This is a known-supplier query, not a blind
ground-truth discovery claim.

Saved artifacts: `assemblies.json`, `assemblies.pkl.gz`, individual
`assembly_candidates/*.pkl.gz`, `ground_truth_assemblies.json`,
`ground_truth_assemblies.pkl.gz`, `assemblies_viewer.json`, and
`assemblies.html`. Existing first-cover artifacts remain unchanged.

## Ranking and overlap audit

The first blind assembly has six fragments, two distinct species, 58/175
explicit-atom retention (33.14%), eight source boundary cuts and four unsupported
target connections. The first known-supplier validation has six fragments,
three species, 58/64 retention (90.625%), eight cuts and four connections.
The ground truth is substantially more atom-efficient, but does not have fewer
displayed cuts/connections in these saved mappings. These counts are not a
balanced chemical reaction's bond-edit count.

At the time of this audit, beta ordering compared species count before retention and did not rank
by cuts/connections. Its output quota also stops discovery before global ranking
can be certified. Changing display order alone cannot fix absent candidates.
No ground-truth-specific ranking bonus or discovery rule has been added.

The former grey atoms in blind pattern 1 were overlapping claims from two
copies of R2 (`MCULE-2301670172`), not different R species. The viewer now gives
them R2's color. For overlaps across distinct species it lists all R alternatives
and offers color-selection buttons; this changes display only, not mappings or
which atoms are covered.

## Final beta ranking correction

This lexicographic correction is historical; the Pareto update below supersedes it.

Final ranking is now separate from provisional discovery priority: fewer
fragment units, greater unique-target/all-input explicit-H retention, fewer
source cuts plus unsupported target connections, then fewer distinct species.
No fitted weights or known-source bonuses are used. All selected outputs are
sorted by this objective even when reserving representatives across patterns.

Re-ranking all 23 saved blind full-cover candidates took 0.47 seconds. The new
first blind result has 59.18% retention, six fragments, five cuts and five
connections. The three six-fragment validation patterns outrank all 23 blind
candidates; the seven-fragment validation pattern does not. Validation sets
were compared separately, never inserted into the blind pool. This fixes
ranking, not the earlier discovery stopping criterion or ground-truth absence.
`assemblies_ranked.json` records the objective, scores, and comparison;
`assemblies_ranked.pkl.gz` preserves selected full mappings. `assemblies.html`
is regenerated from this ranked result with separate validation panels.

## Pareto ranking update

Fragment count is no longer a final ranking objective. Explicit-H retention
and structural cuts + connections define Pareto layers; species count only
breaks identical-objective ties. Layer 1 retains both real blind trade-offs:
59.18% retention with 10 structural changes, and 58.0% with 9. Neither dominates
the other. Both are included in the 20 displayed proposals across four patterns.

All four validation patterns, including the seven-fragment variant, now have
the same objective values: 90.625% retention and 12 structural changes. None is
dominated by the 23 saved blind candidates. Each dominates two of those and
is incomparable with the other 21. This is a separate validation comparison,
not ground-truth injection into blind discovery.

Saved reranking took 0.469 seconds. `assemblies_pareto.json` stores layer labels,
metrics and comparisons; `assemblies_pareto.pkl.gz` retains full mappings.
Earlier lexicographic artifacts remain intact. The viewer displays Pareto layer
labels explicitly and does not claim a preferred order within each layer.
