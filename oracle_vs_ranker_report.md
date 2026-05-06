# Oracle vs ranker — by-step report

Per-step measurement: each metric assigns one number to a (BGCP step,
ranker) pair, computed by 3N-Cartesian cosine similarity between the
selected mode displacement and the GT TS's default reaction-mode
displacement.

Underlying data: 160 BGCP steps, 20 initial guesses per step.

CSV: `out/mode_analysis/oracle_vs_ranker_per_step.csv`
(one row per step, six columns of per-step metrics).

---

## 1. Definitions

| metric | k | meaning |
|---|---|---|
| **Oracle (k=20, any mode)** | 20 | Best `gt_alignment` over **all 20 IGs × all modes** (real + imag). With perfect knowledge of GT, the upper bound on what's achievable. |
| Oracle imag (k=20) | 20 | Same, but restricted to imaginary modes only. |
| Ranker top-1 | 1 | `gt_alignment` of the **single** IG with highest bond_overlap (its bond_overlap-picked mode). |
| **Ranker top-2** | 2 | Best `gt_alignment` over the **top-2** IGs by bond_overlap (their bond_overlap-picked modes). The live consumed output. |
| Ranker top-3 | 3 | Same, top-3 IGs. |

Range `[0, 1]`. `1.0` ↔ modes parallel/anti-parallel (sign-blind).
`0.0` ↔ modes orthogonal.

---

## 2. Aggregate statistics

| metric | mean | median | std | min | max |
|---|---|---|---|---|---|
| Oracle (k=20, any mode) | **0.738** | 0.742 | 0.194 | 0.346 | 1.000 |
| Oracle imag (k=20)      | 0.676 | 0.726 | 0.267 | 0.000 | 1.000 |
| Ranker top-1            | 0.454 | 0.390 | 0.314 | 0.000 | 1.000 |
| **Ranker top-2**        | **0.538** | **0.537** | 0.313 | 0.026 | 1.000 |
| Ranker top-3            | 0.564 | 0.590 | 0.306 | 0.041 | 1.000 |

### Headline

- The "k=20 oracle" — having full knowledge of GT and freedom to
  pick any (IG, mode) pair — averages 0.738. **Even with cheating,
  the IG generation can't do better than ~74% mean alignment** with
  GT reaction modes.
- The live ranker at **k=2** averages 0.538 — recovers about 73 % of
  the oracle's mean.
- The k=1 → k=2 step gives the biggest jump: +0.084 mean. Going
  k=2 → k=3 only adds +0.026.

---

## 3. Cumulative coverage

Fraction of steps where the metric reaches each threshold.

| threshold | Oracle (k=20, any) | Oracle imag (k=20) | Ranker top-1 | **Ranker top-2** | Ranker top-3 |
|---|---|---|---|---|---|
| ≥ 0.95 | 21.2 % | 20.6 % |  6.2 % | **10.6 %** | 13.1 % |
| ≥ 0.90 | 28.7 % | 28.1 % | 13.1 % | **18.8 %** | 21.2 % |
| ≥ 0.80 | 40.6 % | 38.8 % | 19.4 % | **27.5 %** | 29.4 % |
| ≥ 0.70 | 56.9 % | 53.8 % | 28.7 % | **38.1 %** | 41.9 % |
| ≥ 0.50 | **85.6 %** | 71.2 % | 41.2 % | **51.9 %** | 55.6 % |
| ≥ 0.30 | **100.0 %** | 91.2 % | 60.6 % | **70.6 %** | 74.4 % |

**Reading the table:**

- 100 % of steps have *some* IG mode with alignment ≥ 0.3 (no step is
  totally hopeless in principle).
- 86 % of steps have an IG mode with alignment ≥ 0.5 — there's a
  reasonable answer in 86 % of cases.
- The ranker at **top-2 reaches 0.5 alignment 52 % of the time** —
  that is, **the ranker successfully retrieves a "decent" pick in
  about half of all steps**.

### Coverage gap (oracle − ranker)

| threshold | oracle (k=20) | ranker top-2 | absolute gap | recovery rate |
|---|---|---|---|---|
| ≥ 0.95 | 21.2 % | 10.6 % | 10.6 pp | 50 % |
| ≥ 0.9  | 28.7 % | 18.8 % |  9.9 pp | 65 % |
| ≥ 0.8  | 40.6 % | 27.5 % | 13.1 pp | 68 % |
| ≥ 0.7  | 56.9 % | 38.1 % | 18.8 pp | 67 % |
| ≥ 0.5  | 85.6 % | 51.9 % | 33.7 pp | 61 % |

"Recovery rate" = ranker / oracle. The ranker recovers about **two
thirds of the oracle's hits** at every quality threshold from 0.5 to
0.9. There's no severe degradation at any level — the ranker is
proportionally consistent.

---

## 4. Distribution histograms

```
  bin             oracle k=20    ranker top-2   ranker top-1
  [0.0, 0.1)         0 ( 0.0%)    15 ( 9.4%)    27 (16.9%)
  [0.1, 0.2)         0 ( 0.0%)    16 (10.0%)    18 (11.2%)
  [0.2, 0.3)         0 ( 0.0%)    16 (10.0%)    18 (11.2%)
  [0.3, 0.4)         4 ( 2.5%)    16 (10.0%)    19 (11.9%)
  [0.4, 0.5)        19 (11.9%)    14 ( 8.8%)    12 ( 7.5%)
  [0.5, 0.6)        23 (14.4%)    10 ( 6.2%)    10 ( 6.2%)
  [0.6, 0.7)        23 (14.4%)    12 ( 7.5%)    10 ( 6.2%)
  [0.7, 0.8)        26 (16.2%)    17 (10.6%)    15 ( 9.4%)
  [0.8, 0.9)        19 (11.9%)    14 ( 8.8%)    10 ( 6.2%)
  [0.9, 0.95)       12 ( 7.5%)    13 ( 8.1%)    11 ( 6.9%)
  [0.95, 1.0]       34 (21.2%)    17 (10.6%)    10 ( 6.2%)
```

### Shape comparison

- **Oracle (k=20)** is **right-skewed**: 0 % below 0.3, 21 % at or above 0.95. This reflects the IG generator's quality — when it works it works well; the failures pile up in the middle (0.4–0.7 range).

- **Ranker top-2** has a **bimodal-ish shape**: a "tail" of ~30 % of steps below 0.4 (where the ranker picked something far from GT), a fairly even spread in the 0.5–0.95 range, and a high-end peak at 0.95+ (10.6 %). The tail at the bottom is what's costing us — many of those steps have an oracle ≥ 0.5 but the ranker missed.

- **Ranker top-1** has a heavier left tail (39 % below 0.4) and a smaller right peak (6 %). The k=1 → k=2 lift mostly moves probability mass from the [0, 0.4) tail to [0.7, 1.0].

---

## 5. Per-step ranker/oracle ratio

For each step, divide the ranker's score by the oracle's:

| | k=1 | k=2 | k=3 |
|---|---|---|---|
| mean ratio | 0.586 | **0.694** | 0.733 |
| median ratio | — | 0.812 | — |
| min ratio | — | 0.046 | — |

So at **k=2**, the ranker's score is on average **69 % of what's
achievable** (median 81 %). The minimum ratio (0.046) is a
worst-case step where the ranker is essentially blind despite the
oracle finding a 0.7+ alignment elsewhere.

---

## 6. Decomposition by oracle quality

How the ranker performs across the oracle distribution:

| oracle bin | n steps | mean ranker top-2 | mean oracle | ranker / oracle |
|---|---|---|---|---|
| [0.3, 0.5) | 23 | 0.305 | 0.412 | 0.74 |
| [0.5, 0.7) | 46 | 0.435 | 0.598 | 0.73 |
| [0.7, 0.9) | 45 | 0.595 | 0.789 | 0.75 |
| [0.9, 1.0] | 46 | 0.770 | 0.962 | 0.80 |

The ranker is most efficient in the **easy regime** (oracle ≥ 0.9):
it recovers 80 % of the oracle's score because near-perfect IGs are
easy to identify (their bond_overlap is also high). In the
mid-difficulty bins (0.5–0.9), recovery drops to ~73 %. This is
where ranker mistakes hurt most.

---

## 7. Practical interpretation

**For the consumer of the top-2 picks:**

- If you draw ten random steps from the BGCP set, expect:
  - ~5 steps where the ranker top-2 has alignment ≥ 0.5 (chemistry
    captured at least roughly).
  - ~4 of those 5 will reach ≥ 0.7 (clearly the right reaction
    direction).
  - ~3 of those 4 will reach ≥ 0.9 (essentially correct).
  - ~5 of the original 10 are below 0.5 (questionable picks).

- The **gap to the upper bound is roughly half ranker error and
  half upstream IG quality**. The lost-cause steps (oracle < 0.5,
  14 % of dataset) cannot be fixed by any ranker.

- The **top-2 over top-1 always helps** — every threshold level
  benefits by 5–10 percentage points. The marginal value of going
  to top-3 is much smaller (1–4 pp).

---

## 8. Summary table

| | top-1 | **top-2** | top-3 | oracle (k=20) |
|---|---|---|---|---|
| mean alignment | 0.454 | **0.538** | 0.564 | 0.738 |
| ≥ 0.5 hit rate | 41 % | **52 %** | 56 % | 86 % |
| ≥ 0.7 hit rate | 29 % | **38 %** | 42 % | 57 % |
| ranker / oracle ratio | 0.59 | **0.69** | 0.73 | 1.00 |
| relative recovery of oracle's headline ≥0.5 | 48 % | **61 %** | 65 % | 100 % |

The current top-2 setup is operating at ≈ 70 % of the GT-knowing
oracle. About a quarter of the remaining gap is fundamentally
unrecoverable (lost-cause steps); the other three quarters is
addressable through ranker improvements (see
`ranker_analysis_and_improvements.md`).
