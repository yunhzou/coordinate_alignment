# Ranker analysis & improvement ideas

Documenting where we stand on the vibrational-mode ranker, why we got
here, and what's worth trying next.

---

## 1. Where we are now

**Live ranker:** `bond_overlap` (project mode onto unit-vector reaction
direction built from broken/formed bond list at TS coords). Top-2 IG
selection is the consumed output.

**Performance across 160 BGCP steps**, evaluated by cosine similarity
between the ranker's picked mode and the GT TS's default reaction
mode (3N-Cartesian, sign-blind, both reindexed to R atom-frame):

| ranker            | k | mean  | median | ≥ 0.7 | ≥ 0.5 | gap to oracle |
|-------------------|---|-------|--------|-------|-------|---------------|
| bond_overlap      | 1 | 0.454 | 0.390  | 29 %  | 41 %  | 0.284         |
| **bond_overlap**  | **2** | **0.538** | **0.537** | **38 %** | **52 %** | **0.200** |
| bond_overlap      | 3 | 0.564 | 0.590  | 42 %  | 56 %  | 0.174         |
| dwbo_overlap      | 2 | 0.517 | 0.498  | 34 %  | 50 %  | 0.221         |
| bond × dwbo       | 2 | 0.527 | 0.501  | 36 %  | 50 %  | 0.211         |
| bond + dwbo       | 2 | 0.522 | 0.492  | 36 %  | 49 %  | 0.216         |
| rxn_overlap       | 2 | 0.460 | 0.409  | 28 %  | 42 %  | 0.278         |
| core_fraction     | 2 | 0.487 | 0.479  | 29 %  | 49 %  | 0.251         |
| **oracle**        | – | 0.738 | 0.742  | 57 %  | 85 %  | –             |
| oracle (imag)     | – | 0.676 | 0.726  | 55 %  | 73 %  | –             |

**The k=1 / k=2 ranking flip**

- At k=1, `dwbo_overlap` / `bond × dwbo` slightly outperform plain
  `bond_overlap` (+0.011 mean). The WBO-weighted metric is more
  chemistry-aware on average.
- At k=2, `bond_overlap` *beats* the WBO-weighted variants by 0.011 –
  0.021. The diversity of bond_overlap's top-2 picks (different IGs
  with different chemistry signals) catches the oracle answer more
  often than the more correlated dwbo top-2.

**Decomposition of the 0.200 mean gap at k=2:**

| population | n | description | gap contribution |
|---|---|---|---|
| Lost-cause steps (oracle imag < 0.5) | 42 / 155 (27 %) | No IG has any mode aligning with GT — all 20 optimizers landed on different saddles. | ~0.16 of total |
| Recoverable but ranker missed | ~30 / 113 of recoverable steps | A good IG exists, ranker picks something else. | ~0.27 per affected step → bulk of remaining gap |

So **about half** the addressable gap is "ranker error on recoverable
cases", which a better metric could potentially close. The other half
is **fundamental** to the IG generation quality.

---

## 2. What we tried

| metric | what it captures | result |
|---|---|---|
| `core_fraction` | fraction of mode's KE on broken/formed-bond endpoint atoms | weak — insensitive to direction; rewards localized wags |
| `rxn_overlap` | cosine of mode with R→P direction at core atoms | OK — but R→P direction includes spectator path |
| `bond_overlap` | cosine with binary broken/formed bond-stretch direction at TS | best at k=2; clean noise gate (threshold) |
| `dwbo_overlap` | cosine with continuous ΔWBO-weighted reaction direction | best at k=1; noise-sensitive for diffuse WBO |
| `bond × dwbo`, `bond + dwbo`, rank-fusion | hybrids | tied or marginally worse than bond_overlap at k=2 |
| `most-negative imag freq` | pure frequency-magnitude prior | worst by far (0.377 mean) |

**Key learnings:**

1. *Direction matters more than magnitude*. Frequency-only ranking is
   bad. Chemistry-aware metrics dominate.
2. *Threshold gates help on noisy WBO inputs*. The binary
   broken/formed list in `bond_overlap` acts as a noise filter that
   `dwbo_overlap` lacks.
3. *Top-k diversity is a real signal*. The flip at k=2 says raw
   chemistry-similarity scoring leaves correlations on the table.

---

## 3. Why the ranker misses (concrete cases)

From the per-step analysis (`ranker pick vs oracle pick`):

### Case A — picks higher-bond_ov mode of wrong chemistry

```
pr12.Co_Silylation_TS_Dstar-Estar
  GT default     idx  0  freq -221.79  bond_ov 0.703  ← real reaction mode
  Ranker picks   IG iter19  idx 14  freq -420.89  bond_ov 0.695  align 0.065
  Oracle would   IG iter1   idx  0  freq -235.74  bond_ov 0.695  align 0.595
  miss = 0.530
```
Two modes have nearly identical bond_ov (0.703 vs 0.695). Ranker
arbitrarily picks the high-freq one in a different IG.

### Case B — passes up an exact-match IG

```
Jackie_TS_06
  GT default     idx 0  freq -269.75  bond_ov 0.529
  Ranker picks   IG iter7  idx 0  freq -276.91  bond_ov 0.584  align 0.714
  Oracle would   IG iter5  idx 0  freq -269.75  bond_ov 0.529  align 1.000
  miss = 0.286
```
iter5 *converged to literally the same TS as GT* (matching frequency).
But iter7 has 0.05 higher bond_ov — ranker can't recognize "exact
match" so it picks iter7.

### Case C — GT's own mode has low bond_ov

```
pr1.tempo_ts2
  GT default     idx 0  freq -108.76  bond_ov 0.049  ← metric blind!
  Ranker picks   IG iter18  idx 0  freq -431.44  bond_ov 0.330  align 0.201
  Oracle would   IG iter7   idx 11 freq +135.11  bond_ov 0.186  align 0.695  ← real-freq mode
  miss = 0.494
```
GT's reaction mode is itself uninformative under bond_overlap (0.049).
The metric provides no signal. The oracle's pick is even a *real*-
frequency mode that happens to align — which the imag-only ranker
would never consider.

### Pattern

Most ranker errors fall into:
- **scoring noise tie-breaking** — multiple IGs have similar high
  bond_ov; ranker arbitrarily picks the one with marginally higher
  score (Jackie_TS_06).
- **chemistry mismatch** — ranker's pick is a different *kind* of
  mode (different freq band, different atoms moving) but happens to
  score well on bond projection (pr12, pr7.V.dodh_ts1314).
- **GT not in the metric's coverage** — proton/hydride transfers,
  mappings the alignment fails on, partial bond changes below the
  threshold (pr1.tempo, pr5.Noyori).

---

## 4. Improvement ideas, ordered by depth

### Tier 1 — quick wins (a day each)

1. **Diversity-penalized top-k.** When picking top-2 by bond_overlap,
   penalize the second pick if its mode is too similar to the first.
   E.g.:
   ```
   final_score(IG_k) = bond_overlap(IG_k) × (1 - max_j<k cos(IG_k mode, IG_j mode))
   ```
   This pushes the second pick toward different chemistry. Likely
   closes 0.01–0.03 of the gap.

2. **Frequency-band tiebreaker.** When two IGs have bond_overlap
   within 0.05, prefer the one whose freq is closest to the *median
   imag frequency across the IG pool*. The median is a noise-robust
   proxy for "where the reaction freq actually is". Likely fixes Case
   A and Case B above. Does not require GT.

3. **Mass weighting.** Replace `‖d‖²` with `Σ mᵢ‖dᵢ‖²` everywhere.
   Heavy atoms (the actual atoms moving along reaction coordinate in
   most reactions) get correctly weighted up; H wags get correctly
   weighted down. This is the *physically correct* normalization;
   the current metric is a Cartesian approximation.

4. **Hybrid filter-then-rank.** Use `bond_overlap` to identify
   "valid reaction-mode candidates" (≥ 0.4), then rank within that
   set by `dwbo_overlap`. Avoids dwbo's noise-aggregation failure
   while keeping its sensitivity inside the safe region.

### Tier 2 — moderate research (week each)

5. **Cross-IG consensus.** Compute the centroid of all imag modes
   across the IG pool (each mode normalized + reindexed). The mode
   closest to this centroid is the "consensus" reaction direction.
   Then pick the IG whose own bond-overlap-picked mode is closest to
   that consensus. Catches "exact-match" cases like Jackie_TS_06
   because the matching IG sits at the centroid.

6. **Internal-coordinate projection.** Decompose each mode into
   bond, angle, dihedral basis. The reaction mode should be dominated
   by reaction-bond stretches. Easier to compare across geometries
   (internal coords are intrinsic, Cartesian needs alignment).
   Likely closes Case C — internal-coord projection survives even
   when the broken/formed list is empty in Cartesian.

7. **Train a small classifier.** We have 160 × 21 = 3360 (step, IG)
   pairs with all metrics + freq + GT alignment. A simple gradient-
   boosted tree on (bond_ov, dwbo_ov, core_frac, freq, |freq|, mode
   index, n_imag) → predict gt_alignment. Out-of-fold MAE could
   directly tell us how much linear feature combinations leave on
   the table. Baseline benchmark for what's achievable from these
   features alone.

8. **Use R and P modes as cross-validators.** At R (a minimum), the
   soft modes pointing toward TS direction are weak proxies for the
   reaction coordinate. Compute the projection of the IG TS imag
   mode onto R's softest modes; high projection means "TS mode
   continues a soft motion that started at R", which a wrong mode
   wouldn't show.

### Tier 3 — deep research (month each)

9. **Core-restricted Hessian eigendecomposition.** Block-diagonalize
   the Hessian into core and spectator subspaces (using
   `core_atoms`). The core-restricted Hessian's lowest eigenvalue
   gives the *true* reaction mode independent of spectator
   coupling. Compare each IG's core-restricted imag mode to GT's;
   this should discount any spectator-amplitude bias the current
   metric has.

10. **Energy-aware ranking.** Use SCF energies or relative
    enthalpies. The "true" TS sits at a barrier consistent with
    expected ΔE. IGs that converged to wrong saddles often have
    wildly different energies. Cross-reference with literature ΔH‡
    if available, otherwise use per-step IG energy distribution to
    flag outliers.

11. **IRC tangent estimation.** From force gradient at TS, project
    the negative gradient onto each mode. The mode pointing along
    -∇E at TS is the locally-best descent direction (= IRC tangent
    at TS). This is the formal reaction coordinate definition; our
    metrics are approximations of this. Implementation cost: moderate
    (xtb gives gradient already).

12. **Mass-weighted reaction-coordinate as a proper basis.**
    Project every TS mode onto the mass-weighted reaction coordinate
    Δ̃ᵢ = √mᵢ · (rᵢᴾ - rᵢᴿ). Currently rxn_overlap uses unweighted
    Cartesian; mass-weighting respects nuclear-motion physics. Likely
    helps for hydride-transfer cases where the H is light.

13. **Multi-IG joint analysis.** Treat the 20 IGs not as independent
    candidates but as samples from a noisy "TS landscape." Cluster
    IGs by Hessian-mode similarity; identify the cluster centroid
    nearest to a consistent reaction direction; report the centroid
    and its closest IG. This explicitly models "different IGs landing
    on different saddles" and lets the ranker say "reject this whole
    step — no cluster looks reactive".

14. **Surrogate model in mode space.** Train an end-to-end mode
    classifier on a held-out training set of (TS, mode) → "is this
    the reaction mode" labels (using GT alignment ≥ 0.7 as positive).
    Features: full mode displacement vector, atomic descriptors,
    WBO context. Model: small graph neural network or transformer
    on atoms × modes. This is essentially "learn the mode-quality
    function from data" rather than designing it by hand. The ceiling
    is unclear but with 3360 examples and a focused loss, ~0.05–0.10
    closure of the gap to oracle is plausible.

### Tier 4 — beyond ranker, attack the upstream

15. **Better IG generation.** ~0.16 of the 0.200 gap is from
    "lost-cause steps" where no IG matches GT. This is a guess-
    quality problem, not a ranker problem. Investigate:
    - Are there step types (proton transfers, metal complexes)
      where the IG generator systematically fails?
    - Does increasing the IG count from 20 to 40 reduce the lost-
      cause fraction? (Linear ROI?)
    - Can we use the ranker's confidence (e.g., bond_overlap of best
      pick) as a triage signal — flag low-confidence steps for human
      attention or for a second IG-generation round?

---

## 5. Practical recommendation summary

**For now (zero work):** stay on `bond_overlap` top-2. It's the best
operating point of the simple metrics tested, gap to oracle 0.200.

**Highest-impact / low-cost next step:** Tier 1 #1 (diversity-
penalized top-k) — directly addresses the k=1→k=2 lift mechanism we
observed.

**Most rigorous improvement:** Tier 2 #7 (small classifier) — gives
us a principled benchmark for "what's recoverable from these features
alone", and surfaces unexpected feature interactions.

**For the deep gap:** Tier 4 #15 (better IG generation) — but that's
a pipeline-level investment, not a ranker tweak.

The dwbo branch was a useful experiment that confirmed the WBO
gradient information *is* informative on average (k=1 win) but
doesn't translate to k=2 due to correlation. Worth keeping the code
in `eval_dwbo_overlap.py` and `dwbo_reaction_vector` for future
hybrid attempts (Tier 1 #4) but not promoting it to the live ranker.
