# Oracle vs ranker — by-step distribution

Each step is independent (different reaction, different molecule).
Aggregating with mean/std doesn't reflect anything physical, since
the underlying chemistry varies. The right view is **how many steps
reach each alignment threshold** — i.e. the cumulative distribution.

Per-step number = the cosine similarity in 3N-Cartesian displacement
space between a selected mode and the GT TS's default reaction mode.
Range `[0, 1]` (sign-blind, see definitions below).

Underlying data: 160 BGCP steps, 20 initial guesses each, all modes
of every TS available, and the GT's bond_overlap-picked default
mode as the gold reference per step.

CSV: `out/mode_analysis/oracle_vs_ranker_per_step.csv` (one row per
step with all four values for direct inspection).

---

## 1. Per-step measurements

For each step we compute four numbers:

| metric | what's selected | k |
|---|---|---|
| **Oracle (k=20)** | best `gt_alignment` over all 20 IGs × all modes (the upper bound) | 20 |
| Ranker top-1 | `gt_alignment` of the bond_overlap-picked mode of the IG with highest bond_overlap | 1 |
| **Ranker top-2** | best `gt_alignment` over the top-2 IGs by bond_overlap (their picked modes) | 2 |
| Ranker top-3 | same, top-3 IGs | 3 |

Ranker = bond_overlap. Live consumed output is **top-2**.

---

## 2. Cumulative coverage (key result)

Number of steps whose metric reaches each threshold (out of 160):

| threshold | oracle (k=20) | ranker top-1 | **ranker top-2** | ranker top-3 |
|---|---|---|---|---|
| ≥ 0.95 |  34 (21.2 %) |  10 ( 6.2 %) |  **17 (10.6 %)** |  21 (13.1 %) |
| ≥ 0.90 |  46 (28.7 %) |  21 (13.1 %) |  **30 (18.8 %)** |  34 (21.2 %) |
| ≥ 0.85 |  55 (34.4 %) |  27 (16.9 %) |  **38 (23.8 %)** |  39 (24.4 %) |
| ≥ 0.80 |  65 (40.6 %) |  31 (19.4 %) |  **44 (27.5 %)** |  47 (29.4 %) |
| ≥ 0.75 |  75 (46.9 %) |  37 (23.1 %) |  **52 (32.5 %)** |  55 (34.4 %) |
| ≥ 0.70 |  91 (56.9 %) |  46 (28.7 %) |  **61 (38.1 %)** |  67 (41.9 %) |
| ≥ 0.65 | 103 (64.4 %) |  52 (32.5 %) |  **67 (41.9 %)** |  72 (45.0 %) |
| ≥ 0.60 | 114 (71.2 %) |  56 (35.0 %) |  **73 (45.6 %)** |  79 (49.4 %) |
| ≥ 0.55 | 124 (77.5 %) |  63 (39.4 %) |  **79 (49.4 %)** |  83 (51.9 %) |
| ≥ 0.50 | 137 (85.6 %) |  66 (41.2 %) |  **83 (51.9 %)** |  89 (55.6 %) |
| ≥ 0.40 | 156 (97.5 %) |  78 (48.8 %) |  **97 (60.6 %)** | 101 (63.1 %) |
| ≥ 0.30 | 160 (100.0 %) | 97 (60.6 %) | **113 (70.6 %)** | 119 (74.4 %) |
| ≥ 0.20 | 160 (100.0 %) | 115 (71.9 %) | **129 (80.6 %)** | 135 (84.4 %) |
| ≥ 0.10 | 160 (100.0 %) | 133 (83.1 %) | **145 (90.6 %)** | 152 (95.0 %) |

### Reading the table

- **Every step has at least one IG with alignment ≥ 0.3** (oracle 100 % at threshold 0.3). No step is utterly hopeless.
- **86 % of steps have oracle ≥ 0.5**: in 86 % of cases there exists a "decent" IG mode that's at least loosely the right reaction direction.
- **At top-2, the ranker reaches threshold 0.5 on 52 % of steps** — captures slightly over half. Drops to 0.7 → 38 % of steps.
- The **gap at each threshold** between top-2 and oracle:
  - ≥ 0.7 → 91 vs 61 steps (30-step gap)
  - ≥ 0.5 → 137 vs 83 steps (54-step gap)
  - ≥ 0.3 → 160 vs 113 steps (47-step gap)

So at the strict end (≥ 0.7) the gap is 30 steps; at the loose end (≥ 0.3) it's 47 steps. The ranker scales reasonably — it catches the easy cases and misses harder ones at proportional rates.

---

## 3. Histograms

How the per-step values are distributed across alignment bins.

```
          bin           oracle       top-1       top-2       top-3
    [0.00, 0.10)        0  ( 0.0%)   27 (16.9%)   15 ( 9.4%)    8 ( 5.0%)
    [0.10, 0.20)        0  ( 0.0%)   18 (11.2%)   16 (10.0%)   17 (10.6%)
    [0.20, 0.30)        0  ( 0.0%)   18 (11.2%)   16 (10.0%)   16 (10.0%)
    [0.30, 0.40)        4  ( 2.5%)   19 (11.9%)   16 (10.0%)   18 (11.2%)
    [0.40, 0.50)       19  (11.9%)   12 ( 7.5%)   14 ( 8.8%)   12 ( 7.5%)
    [0.50, 0.60)       23  (14.4%)   10 ( 6.2%)   10 ( 6.2%)   10 ( 6.2%)
    [0.60, 0.70)       23  (14.4%)   10 ( 6.2%)   12 ( 7.5%)   12 ( 7.5%)
    [0.70, 0.80)       26  (16.2%)   15 ( 9.4%)   17 (10.6%)   20 (12.5%)
    [0.80, 0.90)       19  (11.9%)   10 ( 6.2%)   14 ( 8.8%)   13 ( 8.1%)
    [0.90, 0.95)       12  ( 7.5%)   11 ( 6.9%)   13 ( 8.1%)   13 ( 8.1%)
    [0.95, 1.00]       34  (21.2%)   10 ( 6.2%)   17 (10.6%)   21 (13.1%)
```

### Shape comparison

- **Oracle (k=20)** is right-skewed: mass in the [0.4, 1.0] range,
  peak at [0.95, 1.0] (21 % of steps). No steps below 0.3.

- **Ranker top-2** is bimodal: a "miss" tail in the [0.0, 0.4) range
  (39 % of steps), a fairly even spread across [0.4, 0.95), and a
  "hit" peak at [0.95, 1.0] (11 % of steps).

- **Going k=1 → k=2** moves about 12 steps out of the bottom three
  bins ([0.0, 0.3)) into higher bins. The improvement is concentrated
  in the lower tail rather than at the top.

---

## 4. Per-step distribution shape — what to look at

The CSV at `out/mode_analysis/oracle_vs_ranker_per_step.csv` has one
row per step:

```
step, oracle_k20, ranker_top1, ranker_top2, ranker_top3, gap_top2
```

`gap_top2 = oracle_k20 − ranker_top2` — flags the steps where the ranker
left the most on the table. The largest gaps would be the most
informative cases for "where could a better ranker do better."

For example: sort the CSV by `gap_top2` descending; look at the top
20 entries — those are the steps with both a high-quality IG present
and a ranker that failed to find it. Those are the addressable cases
where a better ranker would actually help.

---

## 5. Definitions reference

`gt_alignment(d_IG, d_GT) = |d_IG · d_GT| / (||d_IG|| · ||d_GT||)`

where d_IG and d_GT are mode displacement vectors flattened to 3N
dimensions, and modes are reindexed into the R-atom-frame so atom i
in d_IG and d_GT refer to the same chemical atom.

| value | meaning |
|---|---|
| 1.0 | parallel modes (sign-arbitrary). Same atomic motions. |
| 0.95 – 1.0 | nearly identical — likely the IG converged to the same TS as GT. |
| 0.7 – 0.95 | similar reaction direction, geometric jitter. ~25–45° angle. |
| 0.5 – 0.7 | shared major direction, substantial differences. ~45–60° angle. |
| 0.3 – 0.5 | weakly correlated. ~60–73° angle. |
| < 0.3 | essentially orthogonal. Different reactions or wrong-direction modes. |
