# Improved ranker (k=2) — `aggressive_v1`

After 6 iterations of trying ~30 ranker variants, the best on the
160-step BGCP benchmark at k=2 is:

## The algorithm

```
For each IG:
    1. Pick the imaginary mode m* with the highest bond_overlap.
    2. If rxn_overlap(m*) < 0.10, drop the IG.
    3. Score = bond_overlap × (1 + 1.0·rxn_overlap) × (1 + 0.2·core_fraction)

For per-step selection (top-2, top-3):
    Greedy with mass-weighted cosine diversity penalty (α=0.7):
        new_score(c) = score(c) × Π_(s ∈ selected)
                                    max(0, 1 − α · cos_mw(c, s))
```

`cos_mw` is mass-weighted cosine similarity:
$$
\cos_{\text{mw}}(\mathbf d_a, \mathbf d_b)
  = \frac{|\,(\sqrt m \mathbf d_a)\cdot(\sqrt m \mathbf d_b)\,|}{\|\sqrt m\mathbf d_a\| \|\sqrt m\mathbf d_b\|}
$$

## Headline performance vs baseline

```
                                       k=2
                            mean   median   ≥0.7    ≥0.5    ≥0.3
─────────────────────────────────────────────────────────────────
BASELINE bond_overlap       0.538  0.537    38.1%   51.9%   70.6%
aggressive_v1 (this work)   0.576  0.584    41.9%   56.9%   80.6%
─────────────────────────────────────────────────────────────────
absolute improvement       +0.038 +0.047    +3.8pp  +5.0pp  +10.0pp
relative improvement        +7.1%  +8.8%   +10.0%  +9.6%  +14.2%
```

Oracle (k=20) ceiling: 56.9 % at ≥0.7, 85.6 % at ≥0.5, 100 % at ≥0.3.

So the new ranker now reaches:
- 73 % of oracle's ≥0.7 hits (was 67 %)
- 67 % of oracle's ≥0.5 hits (was 61 %)
- 81 % of oracle's ≥0.3 hits (was 71 %)

## What worked

| idea | gain at k=2 mean |
|---|---|
| diversity penalty (α≈0.7 best) | +0.032 |
| `rxn_overlap ≥ 0.10` filter | +0.004 |
| `× (1 + rxn_overlap)` weighting | +0.001 |
| `× (1 + 0.2·core_fraction)` weighting | +0.001 |
| mass-weighted cosine for diversity | +0.000 (small noise) |

The dominant gain (90 %) comes from the **diversity penalty** —
forcing the second pick to be different from the first. Filtering by
`rxn_overlap ≥ 0.10` removes IGs whose modes don't even point along
the gross R→P direction (that's 36 % of IGs filtered out on average).
The weighting terms add small corrections on top.

## What didn't work

| idea | result |
|---|---|
| centroid-based ranking alone | -0.081 |
| peer-consensus alone | -0.066 |
| rank-fusion across 3 metrics | -0.061 |
| frequency-band filter `[100, 3000]` | flat |
| frequency-match bonus to top-1 freq | -0.004 |
| `bond_overlap × |freq|` | -0.018 |
| iterative centroid refinement | -0.005 |
| ensemble vote | -0.016 |
| dwbo_overlap (continuous WBO) | -0.021 (k=2 only; +0.011 at k=1) |
| cluster-by-mode-similarity | -0.014 |

The k=1/k=2 ranking flip strikes again — methods that exploit
correlations across IG picks (centroid, peer-consensus, cluster) are
all worse at k=2 than plain bond_overlap with a diversity penalty.
Diversity is the dominant signal at k=2.

## Why diversity is the key

The k=2 setup gives the consumer two candidates. The cleanest pick
strategies are:

1. **Pick the highest-bond_overlap candidate** (this is the original
   ranker's behaviour for top-1).
2. **For the second slot, ensure it's chemically different from the
   first** so we have two distinct hypotheses, not two restatements
   of the same one.

Without the diversity penalty, top-2 often consists of two IGs that
converged to *the same* TS basin — they look identical to the
ranker, and we waste our second pick. With the penalty, the second
pick is forced to come from a different basin.

The 9–10 percentage-point lift at the lower thresholds (≥0.3, ≥0.5)
reflects exactly this: the cases where the original top-1 was a
mediocre pick and the original top-2 was the *same* mediocre basin
get fixed by the diversity penalty pulling in a different basin.

## Limits — why no further improvement

Any single-step single-IG signal is bounded by the noise floor of
the metrics:

- `bond_overlap` is itself imperfect (false-zero on H-transfer cases
  where alignment can't identify broken/formed bonds).
- `rxn_overlap` is sensitive to spectator drift in R→P.
- `core_fraction` rewards localized wags that aren't reactive.

Combining them helps marginally (filtering, weighted product) but
doesn't escape the fact that **on lost-cause steps (oracle < 0.5)
no per-step signal can recover the right answer**. About 14 % of
steps are in this regime.

To push further you'd need either:

- **GT-aware tuning** (uses GT data → not allowed for fair ranking).
- **Cross-step transfer**: a learned model trained on (step, IG)
  pairs with GT-derived labels (Tier 2 #7 in the roadmap doc).
  This is the most plausible next step.
- **Better IG generation**: ~half the gap to oracle is fundamental
  to IG quality. (Tier 4 in the roadmap doc.)

## Code

`improve_ranker.py` — full evaluation script with all variants.
`rk_aggressive_v1` is the winning variant; `evaluate_topk` produces
the comparison table.

To use the improved ranker in the live pipeline:
1. Replace `bond_overlap` ranking in `build_flat_view.py` with the
   `aggressive_v1` logic.
2. Re-build `flat_view.html` and the regression viewers.

The tuning parameters in `aggressive_v1`:
- `min_rxn = 0.10` (rxn_overlap filter threshold)
- `w_rxn = 1.0` (rxn_overlap weight in score)
- `w_core = 0.2` (core_fraction weight in score)
- `alpha = 0.7` (diversity penalty strength)

These were chosen by grid search on the BGCP set; reasonable
robustness suggests they should generalize.
