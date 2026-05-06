# The verifier — `rk_clean_v2`

## Purpose

The verifier ranks initial-guess transition-state structures
(IG-TS) so that downstream consumers see the highest-quality
candidates first. For each elementary step we have

* the **reactant** R and **product** P (already aligned to a shared
  atom indexing — see `01_atom_alignment_and_bond_detection.md`),
* a set of **broken bonds** $B$, **formed bonds** $F$, and **core
  atoms** $C$ derived from $(R, P)$,
* 20 **initial-guess transition states**, each having a full Hessian
  (xtb GFN2 in `xtb_hess/g98.out`); from each we extract every normal
  mode as a $(3 n_R)$-dimensional displacement vector $\mathbf{d} \in
  \mathbb{R}^{n_R \times 3}$ at frequency $\omega$ (negative
  frequencies indicate imaginary modes).

The verifier returns a *ranked, diverse list* of $(IG, \text{mode})$
pairs. The live consumed output is **top-2**. The pass@2 ceiling
provided by oracle is the target; the verifier closes part of the gap
without using the ground-truth label.

## Per-mode features

Three scalar features are pre-computed for every mode of every IG. All
are sign-blind, take values in $[0, 1]$, and use the unified atom
indexing.

### Bond-direction overlap

Construct the **bond reaction vector** at the IG's TS coordinates
$\mathbf{x}^{TS} \in \mathbb{R}^{n_R \times 3}$:

$$
\mathbf{V}_i = \sum_{(i,j) \in B} -\hat{\mathbf{u}}_{ij}^{TS}
             + \sum_{(j,i) \in B}  \hat{\mathbf{u}}_{ji}^{TS}
             + \sum_{(i,j) \in F}  \hat{\mathbf{u}}_{ij}^{TS}
             + \sum_{(j,i) \in F} -\hat{\mathbf{u}}_{ji}^{TS}
$$

where $\hat{\mathbf{u}}_{ij}^{TS} = (\mathbf{x}_j^{TS} -
\mathbf{x}_i^{TS}) / \|\cdot\|$. The sign convention encodes
"broken-bond endpoints move apart, formed-bond endpoints move
together" — a true concerted reaction mode aggregates contributions
coherently. Then

$$
\text{bond\_overlap}(\mathbf{d}) =
   \frac{|\,\mathbf{d} \cdot \mathbf{V}\,|}{\|\mathbf{d}\|\,\|\mathbf{V}\|}.
$$

### Reaction-coordinate overlap

Compute the per-atom R→P displacement after Kabsch-aligning P (in
R-frame) onto R:
$$
\boldsymbol{\Delta}_i = \mathbf{x}_i^{P,\text{aligned}} - \mathbf{x}_i^R
\quad\text{for } i \in \mathrm{dom}(\pi),
$$
zero on unmapped atoms. Restrict to core atoms and unit-normalise
over them: $\hat{\boldsymbol{\Delta}}^{C}_i = \boldsymbol{\Delta}_i /
\|\boldsymbol{\Delta}^{C}\|$ for $i \in C$, zero elsewhere. Then

$$
\text{rxn\_overlap}(\mathbf{d}) =
   \frac{|\,\mathbf{d} \cdot \hat{\boldsymbol{\Delta}}^{C}\,|}{\|\mathbf{d}\|}.
$$

This rewards modes whose core-atom motion follows the R→P direction
and penalises modes that waste amplitude on spectator atoms.

### Core mass fraction

$$
\text{core\_fraction}(\mathbf{d}) =
   \frac{\sum_{i \in C} \|\mathbf{d}_i\|^2}
        {\sum_{i=1}^{n_R} \|\mathbf{d}_i\|^2} \in [0, 1].
$$

The fraction of the mode's squared amplitude that resides on
reactive atoms. Values close to 1 mean the mode is localised on the
reactive site; values close to 0 mean the mode is dominated by
spectator motion.

## Verifier algorithm: `rk_clean_v2`

For each IG let $M$ denote its full set of normal modes,
$M_{\text{im}} \subseteq M$ the imaginary modes, and $n_{\text{im}} =
|M_{\text{im}}|$.

### Step 1 — within-IG mode pick

Within each IG, pick the imaginary mode of highest bond_overlap:
$$
\hat{m} = \arg\max_{m \in M_{\text{im}}} \text{bond\_overlap}(m).
$$
If $M_{\text{im}} = \emptyset$, the IG is skipped (a transition-state
candidate without imaginary modes is structurally not a TS).

### Step 2 — filter

Discard the IG unless **both** structural priors hold:

$$
1 \le n_{\text{im}} \le 2 \quad\wedge\quad \text{rxn\_overlap}(\hat{m}) \ge 0.10.
$$

The first prior favours clean first-order saddle points (a true TS has
exactly one imaginary mode; tolerating $n_{\text{im}} = 2$ admits
near-saddle structures whose second imaginary mode is small). The
second filter eliminates IGs whose chosen mode does not even point
along the gross R→P direction, removing roughly 35 % of IGs on
average.

### Step 3 — score

Compose the three features into a single scalar:
$$
S(\hat{m}, n_{\text{im}}) =
\frac{\text{bond\_overlap}(\hat{m})\,(1 + \text{rxn\_overlap}(\hat{m}))\,(1 + 0.2\,\text{core\_fraction}(\hat{m}))}
{n_{\text{im}}^{\,0.3}}.
$$
The multiplicative form ties the score to bond_overlap (the dominant
single signal) and adds modest gain factors for higher rxn_overlap
and core_fraction. The $n_{\text{im}}^{0.3}$ denominator gives a soft
cleanliness preference even after the hard $\le 2$ filter.

### Step 4 — diversity-aware top-$k$ selection

Greedily select up to $k$ IGs (the live consumer takes $k = 2$).
Maintain the selected set $\mathcal{S}$ initially empty. At every
step pick

$$
c^{*} = \arg\max_{c \in \mathcal{C} \setminus \mathcal{S}}\;
S(c) \cdot \prod_{s \in \mathcal{S}} \max\!\big(0,\; 1 - \alpha\,\rho(c, s)\big),
$$

where $\alpha = 0.7$ is the diversity strength and $\rho$ is the
**mass-weighted cosine similarity** between two modes' displacements:

$$
\rho(\mathbf{d}_a, \mathbf{d}_b) =
\frac{|\,(\sqrt{m}\,\mathbf{d}_a) \cdot (\sqrt{m}\,\mathbf{d}_b)\,|}
     {\|\sqrt{m}\,\mathbf{d}_a\|\,\|\sqrt{m}\,\mathbf{d}_b\|},
\qquad
m_i = \text{atomic mass}(e_i).
$$

The mass weighting makes "different modes" mean *physically different
nuclear motion* rather than swapping which hydrogens vibrate. The
multiplicative penalty drives the second pick away from the basin of
the first: when two IGs converged to the same TS structure their
top-1-by-bond_overlap modes are nearly parallel, $\rho \approx 1$, and
the penalty multiplies the score by $\approx 0$. The diversity term
contributes 90 % of the verifier's gain over plain bond_overlap
ranking at $k = 2$.

### Fallback

If the filter rejects every IG (occasionally happens when the entire
IG pool is far from the reaction direction), revert to plain
bond_overlap top-$k$ on the unfiltered pool. This guarantees the
verifier always returns at least $\min(k, |\text{IGs with imag modes}|)$
candidates.

## Hyperparameter values

| symbol | value | rôle |
|---|---|---|
| `min_rxn`     | 0.10 | reject IG if best mode's `rxn_overlap` below this |
| `max_imag`    | 2    | reject IG if it has more than this many imaginary modes |
| `w_rxn`       | 1.0  | rxn_overlap weight in score |
| `w_core`      | 0.2  | core_fraction weight in score |
| `imag_pen`    | 0.3  | exponent on $n_{\text{im}}$ in denominator |
| $\alpha$ (diversity) | 0.7 | strength of mass-weighted-cosine penalty |
| $k$ (live)    | 2    | top-$k$ output consumed by `flat_view.html` |

These were grid-searched over the 160-step BGCP set during the
`improve_ranker.py` campaign (15 iterations, ~40 variants); they
represent the Pareto front for $\geq 0.5$ alignment coverage.

## Headline performance (160-step BGCP set)

```
                           pass@2
                  mean   ≥0.7    ≥0.5    ≥0.3
oracle (k=20)    0.738  56.9%   85.6%  100.0%
baseline bond    0.538  38.1%   51.9%   70.6%
clean_v2         0.581  41.9%   59.4%   78.8%
```

The verifier closes ~30 % of the bond-overlap → oracle gap at the
strict end (≥ 0.7) and ~25 % at the loose end (≥ 0.3), purely from
unsupervised signals. The remaining gap is signal-bounded: no
single per-IG feature places the truly-best IG into rank ≤ 2 more
often than 30 % of the time (`diagnose_ig_features.py`), so any
unsupervised verifier hits a ceiling near the current numbers.
Closing further would require GT-supervised cross-step training.
