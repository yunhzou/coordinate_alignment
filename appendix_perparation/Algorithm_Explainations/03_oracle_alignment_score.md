# Oracle alignment score

## Definition

The **alignment score** between two normal modes is the
sign-blind cosine similarity in the $3 n_R$-dimensional Cartesian
displacement space, after both modes have been re-indexed into the
shared R-frame:

$$
a(\mathbf{d}_a, \mathbf{d}_b) =
   \frac{|\,\mathbf{d}_a \cdot \mathbf{d}_b\,|}
        {\|\mathbf{d}_a\|\,\|\mathbf{d}_b\|}
\;\in\; [0, 1],
$$

where $\mathbf{d}_a, \mathbf{d}_b \in \mathbb{R}^{n_R \times 3}$ are
the per-atom 3-vectors of the two modes flattened to a single vector
of length $3 n_R$. Atom $i$ in $\mathbf{d}_a$ and $\mathbf{d}_b$
refers to the same chemical atom because both displacements are
pre-aligned to the canonical R-frame indexing produced by
`align_bgcp_coords.py`.

**Sign-blind.** The eigenvectors of a Hessian have arbitrary
sign — a mode and its negation describe the same physical motion. The
absolute value collapses both into a single equivalence class.

## Geometric interpretation

Let $\theta$ be the angle between $\mathbf{d}_a$ and $\mathbf{d}_b$
(after sign-fixing to make the dot product non-negative). Then
$a = \cos\theta$ and the value bands have stable physical meanings:

| value range | angle | meaning |
|---|---|---|
| $a = 1.0$    | $0°$    | parallel — same atomic motions, same TS basin |
| $0.95 \le a < 1.0$ | $< 18°$ | nearly identical — the IG converged to the same TS as GT, modulo numerical noise |
| $0.7 \le a < 0.95$ | $18°$–$45°$ | shared major reaction direction with geometric jitter (similar TS, slightly different conformer) |
| $0.5 \le a < 0.7$  | $45°$–$60°$ | shared dominant component but substantial differences |
| $0.3 \le a < 0.5$  | $60°$–$73°$ | weakly correlated; partial directional overlap only |
| $a < 0.3$         | $> 73°$ | essentially orthogonal — wrong-direction or unrelated mode |

Because $\mathbf{d}_a$ and $\mathbf{d}_b$ are full-molecule vectors,
spectator-atom motion shows up in both numerator and denominator. A
mode that moves spectators wastefully gets a smaller alignment than a
mode whose amplitude is concentrated on the same reactive atoms as
the reference.

## The "oracle" qualifier

We call an alignment score the **oracle score** when one of the two
modes is the **ground-truth reaction mode** $\mathbf{d}_{\text{GT}}$
— the mode at the genuine TS structure (computed from a literature- or
DFT-derived reference geometry, then propagated through the same xtb
GFN2 Hessian pipeline that the IGs use). The qualifier "oracle"
indicates that this score *uses* the GT label and so cannot be used
to *rank* IGs in production; it is the evaluation target that the
unsupervised verifier (`02_verifier.md`) is trying to approximate.

For each $(IG, \text{mode})$ pair we compute
$a(\mathbf{d}_{\text{mode}}, \mathbf{d}_{\text{GT}})$. Aggregating
this across an IG's modes and across IGs gives several useful
statistics.

## Per-IG aggregations

For IG $g$ with modes $M^{(g)}$:

$$
\text{oracle\_any}(g) = \max_{m \in M^{(g)}} a(\mathbf{d}_m, \mathbf{d}_{\text{GT}})
$$
$$
\text{oracle\_imag}(g) = \max_{m \in M^{(g)} \,:\, \omega_m < 0} a(\mathbf{d}_m, \mathbf{d}_{\text{GT}}).
$$

`oracle_any` is the best alignment achievable from this IG given an
omniscient choice of mode. `oracle_imag` restricts to imaginary modes
— physically, the GT mode is imaginary, so this restriction matches
the canonical TS-mode definition and is the more meaningful quantity
for downstream chemistry.

## Per-step aggregations

For a step with IG pool $\mathcal{G} = \{g_1, \dots, g_{20}\}$:

$$
\text{oracle\_step}^{(k)} = \max_{g \in \text{top-}k(\mathcal{G})} \text{oracle\_any}(g),
$$

where the top-$k$ is over IGs ranked by `oracle_any` descending. Two
edge cases:

* $\text{oracle\_step}^{(20)}$ is the absolute pass@20 ceiling — the
  best alignment achievable on this step given the IG pool we have.
* $\text{oracle\_step}^{(1)}$ is the best individual IG, equal to the
  pass@$\infty$ since further IGs can only equal it.

## Verifier comparison

The same aggregation can be computed against the verifier's chosen
top-$k$ IGs (`02_verifier.md`):

$$
\text{verifier\_step}^{(k)} = \max_{g \in \text{top-}k_{\text{verifier}}(\mathcal{G})} \text{oracle\_any}(g).
$$

The gap $\text{oracle\_step}^{(k)} - \text{verifier\_step}^{(k)}$ is
the *selection penalty*: the alignment quality lost because the
verifier picked the wrong $k$ IGs out of 20. Note that we credit the
verifier with `oracle_any` of its picked IGs (the best mode in those
IGs), not with the alignment of the specific mode the verifier itself
chose. This separates "did the verifier pick the right basin?" from
"did it pick the right mode within the basin?" — gap diagnostics
(`diagnose_gap.py`) show that 76 % of the gap is the IG-selection
problem, not the mode-selection problem.

## Distribution properties (160-step BGCP set, before duplicate cleanup)

```
                       mean    ≥0.7    ≥0.5    ≥0.3
oracle_step (any)     0.738   56.9%   85.6%  100.0%
oracle_step (imag)    0.676   53.8%   71.2%   91.2%
oracle_top2 (any)     0.697   50.0%   78.1%  100.0%
oracle_top3 (any)     0.667   43.1%   73.1%  100.0%
oracle_top10 (any)    0.572   26.2%   55.6%   95.6%
oracle_top20 (any)    0.391    9.4%   21.2%   61.9%
```

Salient observations:

* **No step is hopeless.** Every step has at least one IG with
  alignment $\ge 0.3$; 86 % have at least one with $\ge 0.5$.
* **Imaginary-only is meaningful.** Restricting to imaginary modes
  drops the ceiling from 100 % to 91 % at the $\ge 0.3$ level — about
  9 % of steps achieve their best `oracle_any` on a *real* mode that
  happens to align with the GT direction, often a low-frequency wag
  whose Cartesian shape coincidentally tracks the reaction mode. This
  is geometrically informative but physically unjustified credit.
* **Multiple good IGs are common.** The 2nd-best IG by oracle still
  hits $\ge 0.7$ on 50 % of steps — only 7 percentage points below
  the absolute oracle. Ten IGs all simultaneously reaching $\ge 0.5$
  on 56 % of steps means most steps have a small *cluster* of IGs
  that converged to the same TS basin. Diversity-aware top-$k$
  selection exploits this redundancy.

## Persisted artefacts

* `out/mode_analysis/oracle_alignment_per_ig.csv` — one row per
  $(step, IG)$ pair: `best_align_any`, `best_align_imag`, plus the
  feature values of the best-aligned mode.
* `out/mode_analysis/oracle_alignment_per_step.csv` — one row per
  step: the absolute ceiling and which IG/mode achieves it.
* `out/mode_analysis/oracle_topk_per_step.csv` — one row per step
  with `top1_any … top20_any` and `top1_imag … top20_imag` columns:
  the $k$-th best IG's alignment, sorted descending.
* `appendix_perparation/analtics/final_quality_measurement.csv` — the
  final consolidated table joining oracle top-$k$ and verifier top-$k$
  per step, used for cross-comparison plots and headline tables.
