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
