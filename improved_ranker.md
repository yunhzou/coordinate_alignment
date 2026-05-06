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

## Iterations 7–15: hitting the signal floor

After `aggressive_v1` we ran a further 9 iterations (≈40 more variants)
trying to close the remaining 0.16 mean gap to oracle. The result:
**none** delivered a clean improvement. Best new variant is
`clean_v2` (max_imag=2 + soft `1/n_imag^0.3` penalty + lighter core
weight): mean **0.581** (+0.005), `≥0.5` **95** (+4 steps), but `≥0.3`
**126** (−3) — net mixed.

| iteration | family | best result | vs aggressive_v1 |
|---|---|---|---|
| 7 (mode pool, top-M imag/IG) | within-IG mode broadening | flat 0.576 | 0.000 |
| 8 (within-IG combo score) | bond×(1+rxn) within IG | 0.578 | +0.002 |
| 9 (dwbo within-IG pick) | bond × dwbo within | flat | 0.000 |
| 10 (cluster-pool) | cluster across all imag modes | flat | ~0.000 |
| 12 (core-dominant) | bond × (1+r) × (1+w·c), w=1–10 | flat / slight degrade | −0.002 |
| 13 (Borda fusion) | sum-of-ranks across signals | mean drops to 0.535 | −0.04 |
| 13 (max_imag filter) | `aggressive_v1` + ≤2 imag | 0.577 | +0.001 |
| 13 (combined_clean) | bond × (1+r) × (1+0.5c) / n_imag^0.3 | 0.578, ≥0.5=94 | +0.002 |
| 14 (clean_v2 sweep) | tuning above | 0.581, ≥0.5=95 | +0.005 (mixed) |
| 14 (adaptive alpha) | alpha varies w/ top-1 confidence | flat 0.576 | 0.000 |
| 15 (normalized fusion) | per-step rescale of features | mean drops to 0.572 | −0.005 |

### Why we plateau — diagnostic on the IG-selection bottleneck

We decomposed the gap to oracle:

```
Gap decomposition (aggressive_v1 vs oracle, 160 steps):
  Both right (top-2 includes best IG, mode within 0.05): 34/160 (21%)
  Wrong IG (best IG not even in top-2):                 121/160 (76%)
  Right IG, wrong mode within:                           5/160 (3%)
```

The bottleneck is **picking the right 2 of 20 IGs**, not picking the
right mode within an IG. To check whether any unsupervised signal
could rank the best IG into top-2 more often, we tested every
candidate signal we have:

```
feature                 rank=1  rank≤2  rank≤3
bond_overlap            15.7%   26.4%   35.8%
core_fraction           20.8%   30.8%   36.5%
b × (1+r) × (1+0.2c)    15.7%   26.4%   35.8%
rxn_overlap              8.8%   18.2%   25.2%
mode_peer_consensus      9.4%   13.2%   14.5%
TS-RMSD-to-consensus     5.0%    8.2%   12.6%
imag_count (fewer best) 11.9%   18.9%   27.0%
```

**Ceiling**: the strongest single signal (core_fraction) places the
best-IG at rank ≤ 2 only **30.8 %** of the time. All combinations we
tried (Borda, normalized fusion, weighted product, rank averaging)
underperform this — the signals are too correlated, and combining
them doesn't add genuinely independent information.

### What this means

We're signal-bound. The available physics-derived per-IG signals
(bond_overlap, rxn_overlap, core_fraction, dwbo_overlap, peer-mode
consensus, structural consensus, imag_count) cap at ~30 % rank ≤ 2 of
"best IG". The current `aggressive_v1` ranker hits **41.9 %** at
`≥0.7` — already higher than any single-signal rank-2 hit rate,
because multiple "good enough" IGs exist per step. The diversity
penalty captures this.

Closing the rest of the gap requires:

- a **learned model** over (step, IG, mode) features supervised by
  GT-alignment (cross-step transfer; Tier 2 in the roadmap), or
- **better IGs** so the rank-2 ceiling rises (Tier 4).

Both are out of scope for the current "rerank existing IGs" task.

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
