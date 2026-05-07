# Appendix A — benchmark and pipeline detail

## A.1 Benchmark composition

The benchmark contains $N = 155$ unique elementary reaction steps
after deduplication of five redundant cyclobutane-ring-opening
variants and one duplicate Jackie-series TS. Step sizes range from
18 to 149 atoms (mean 58.9, median 57, IQR [34, 82]). The
distribution is bimodal: a small-molecule cluster at 20–40 atoms
(carbocation rearrangements, organocatalysis, simple pericyclic) and
a transition-metal cluster at 60–90 atoms (Pd / Ni / V / Fe /
Co / Au complexes with phosphine, NHC, or cyclopentadienyl ligands).
Element coverage by atom-count bin grows from C, H only at the
smallest sizes to a 12-element set including 4d/5d transition metals
in the 90+ bin.

For every step the dataset provides:
- one reactant geometry $\mathrm{R}$ (xyz),
- one product geometry $\mathrm{P}$ (xyz),
- one ground-truth TS geometry (xyz, "GT"),
- 20 initial-guess TS geometries (xyz, "IG").

For every geometry we compute a single-point Wiberg bond order matrix
$W \in \mathbb{R}^{n \times n}$ and a Hessian using xtb at the
GFN2-xTB level of theory. Hessians are parsed into per-mode 3$N$-
Cartesian displacement vectors $\mathbf{d}_m \in \mathbb{R}^{n \times 3}$
and per-mode harmonic frequencies $\omega_m$.

## A.2 Per-step bond statistics

Distribution of broken-bond / formed-bond counts:

| count | broken (#steps) | formed (#steps) | total (#steps) |
|---|---|---|---|
| 0 | 26 | 11 | 6 |
| 1 | 70 | 93 | 25 |
| 2 | 50 | 42 | 52 |
| 3 |  5 |  6 | 33 |
| 4 |  4 |  3 | 28 |
| 5+ | 0 |  0 | 11 |

The most common profile is one broken and one formed bond (52 of
155 steps). Steps with `broken = 0 ∧ formed = 0` (six cases) are
discussed in Appendix C as bond-identification edge cases.


# Appendix B — Algorithmic detail

## B.1 Atom mapping (priority-queue alignment)

### Bond graph

For each structure we build an undirected graph $G = (V, E, w)$ from
$(e, W)$ using a single threshold $w_{\text{floor}} = 0.2$:
$$
V = \{1, \dots, n\}, \qquad
E = \{(i, j) : W_{ij} \ge w_{\text{floor}}\}, \qquad
w(i, j) = W_{ij}.
$$
Each node carries its element symbol. The threshold is loose by
design so that weak / partial / dative bonds are admitted; this
makes the graph topology robust to mild WBO fluctuations between R
and P.

### Single-island growth

Given a seed atom $s \in V_R$ that is not yet mapped, we initialise
the candidate mapping set
$$
\mathcal{C} = \{\, \{s \mapsto v\} : v \in V_P \setminus \mathrm{img}(\pi),\;
                                       e_R(s) = e_P(v)\, \}
$$
and a max-heap $H$ of edges $(u, n)$ with $u \in F$ (the fragment so
far) and $n \notin F$, keyed by $w(u, n)$ in descending order. We
iterate until either $|\mathcal{C}| = 1$ and the candidate is
*set-unique* or the heap empties:

1. **Pop** the highest-WBO edge $(u, n)$. Mark it consumed.
2. **Try to extend** every candidate $c \in \mathcal{C}$ to include
   $n$: for each $c$ find P-atoms $v$ that are neighbours in $G_P$
   of every $u' \in F$ that is bonded to $n$ in R, satisfy
   $e_R(n) = e_P(v)$, and have WBO discrepancy
   $|W_R(u', n) - W_P(c(u'), v)| \le \tau$ on every bonded edge
   from $n$ into $F$.
3. **Commit or consume.** If any extension survives, replace
   $\mathcal{C}$ with the extended set and add $n$ to $F$. If none
   survive, the edge is consumed and the fragment does not grow.

A candidate set is **set-unique** when every $c_i$ has the same
image set $\mathrm{img}(c_i)$ — the candidates differ only by
permutations of indistinguishable atoms (e.g. methyl hydrogens). If
the heap empties with multiple non-set-unique candidates the
algorithm **branches**: each remaining candidate becomes an
independent hypothesis to be propagated by later seeds.

### Multi-seed driver

We run the single-island growth procedure over $n_{\text{seeds}} = 10$
random seed orderings and accumulate the resulting branches. After
each seed pass, branches are deduplicated by their full mapping
signature (the sorted tuple of $(R\text{-index} \mapsto P\text{-index})$
pairs) and capped at $\text{max\_branches} = 8$.

### Bond classification and scoring

For each branch we classify bonds via
$$
\text{broken}: W_R(i,j) \ge w_{\text{break}} \;\wedge\; W_R(i,j) - W_P(\pi(i), \pi(j)) \ge \Delta_{\text{thr}},
$$
$$
\text{formed}: W_P(\pi(i), \pi(j)) \ge w_{\text{break}} \;\wedge\; W_P(\pi(i), \pi(j)) - W_R(i,j) \ge \Delta_{\text{thr}},
$$
with $w_{\text{break}} = 0.5$ and $\Delta_{\text{thr}} = 0.5$. We
score each branch by the lexicographic tuple
$$
\big( |B| + |F|,\; v_\chi,\; -|\mathrm{dom}(\pi)| \big)
$$
sorted ascending: prefer the simplest reactive change, then fewest
spectator-stereocentre chirality violations, then most atoms mapped.

The chirality term $v_\chi$ counts non-reactive atoms of degree
$\ge 4$ whose mapped neighbours are all preserved (no incident
broken or formed bond), at which the signed determinant of the
first three neighbour-displacement vectors flips sign between
R-frame and P-frame. Near-coplanar centres ($|\det| < 0.05$) are
skipped. The chirality count is the *tiebreaker* between branches
that have the same number of broken plus formed bonds: a
geometrically valid but chirality-inverted mapping is penalised in
favour of the chirality-preserving alternative.

### Tolerances

| symbol | value | role |
|---|---|---|
| $w_{\text{floor}}$ | 0.2 | minimum WBO to admit an edge |
| $\tau$ | 1.0 | maximum $\|W_R - W_P\|$ on a candidate edge |
| $w_{\text{break}}$ | 0.5 | minimum WBO to consider an edge a "bond" for breaking/forming |
| $\Delta_{\text{thr}}$ | 0.5 | minimum $\|\Delta W\|$ to count as a bond change |
| $n_{\text{seeds}}$ | 10 | random seed orderings per step |
| max_branches | 8 | branch cap per pass |

## B.2 Verifier (clean_v2)

### Per-mode features

For every normal mode of every IG we compute three sign-blind
features in $[0, 1]$, all using the unified atom indexing.

**Bond-direction overlap.** Construct the bond reaction vector at
the IG's TS coordinates $\mathbf{x}^{TS}$:
$$
\mathbf{V}_i = \!\!\sum_{(i,j) \in B} -\hat{\mathbf{u}}_{ij}^{TS}
        \;+\!\sum_{(j,i) \in B}  \hat{\mathbf{u}}_{ji}^{TS}
        \;+\!\sum_{(i,j) \in F}  \hat{\mathbf{u}}_{ij}^{TS}
        \;+\!\sum_{(j,i) \in F} -\hat{\mathbf{u}}_{ji}^{TS},
$$
where $\hat{\mathbf{u}}_{ij}^{TS}$ is the unit vector from atom $i$
to atom $j$ at the TS. The sign convention encodes "broken-bond
endpoints separating, formed-bond endpoints approaching" so that a
true concerted reaction mode aggregates contributions coherently.
$$
\mathrm{bond\_overlap}(\mathbf{d}) = \frac{|\,\mathbf{d} \cdot \mathbf{V}\,|}{\|\mathbf{d}\|\,\|\mathbf{V}\|}.
$$

**Reaction-coordinate overlap.** Compute the per-atom R$\to$P
displacement after Kabsch-aligning P (in R-frame) to R, restrict to
core atoms $C$ (atoms touching any broken or formed bond) and
unit-normalise:
$$
\hat{\boldsymbol{\Delta}}^{C}_i =
\boldsymbol{\Delta}_i / \|\boldsymbol{\Delta}^C\| \quad\text{for } i \in C,\quad 0 \text{ else.}
$$
Then
$$
\mathrm{rxn\_overlap}(\mathbf{d}) =
\frac{|\,\mathbf{d} \cdot \hat{\boldsymbol{\Delta}}^{C}\,|}{\|\mathbf{d}\|}.
$$

**Core mass fraction.**
$$
\mathrm{core\_fraction}(\mathbf{d}) =
\frac{\sum_{i \in C} \|\mathbf{d}_i\|^2}
     {\sum_{i = 1}^{n_R} \|\mathbf{d}_i\|^2} \;\in [0, 1].
$$

### Within-IG mode pick + filter

For each IG let $M_{\text{im}} = \{m : \omega_m < 0\}$ be the
imaginary-mode set. We pick the imaginary mode of highest
bond_overlap:
$$
\hat{m} = \arg\max_{m \in M_{\text{im}}} \mathrm{bond\_overlap}(m).
$$
We then apply two structural priors as a hard filter:
$$
1 \le n_{\text{im}} \le 2 \quad\wedge\quad \mathrm{rxn\_overlap}(\hat{m}) \ge 0.10.
$$
The first reflects the canonical TS condition (a clean first-order
saddle has exactly one imaginary mode; we relax to two for
near-saddle structures whose second imaginary mode is small). The
second eliminates IGs whose chosen mode does not point along the
gross R$\to$P direction. The filter eliminates approximately 35 % of
IGs on average.

### Score

The score of a surviving IG is
$$
S(\hat{m}, n_{\text{im}}) =
\frac{\mathrm{bond\_overlap}(\hat{m}) \cdot
      (1 + w_r \,\mathrm{rxn\_overlap}(\hat{m})) \cdot
      (1 + w_c \,\mathrm{core\_fraction}(\hat{m}))}
     {n_{\text{im}}^{\,p}}
$$
with $w_r = 1.0,\; w_c = 0.2,\; p = 0.3$. The multiplicative form
ties the score to bond_overlap (the dominant single signal), adds
modest gain factors for rxn_overlap and core_fraction, and applies
a soft preference for clean first-order saddles via $n_{\text{im}}^{p}$.

### Diversity-aware top-$k$ selection

Maintain the selected set $\mathcal{S}$ (initially empty). At each
step, pick
$$
c^{\star} = \arg\max_{c \in \mathcal{C} \setminus \mathcal{S}}\;
S(c) \cdot \prod_{s \in \mathcal{S}} \max\!\big(0,\; 1 - \alpha\,\rho(c, s)\big),
$$
where $\alpha = 0.7$ is the diversity strength and $\rho$ is the
mass-weighted cosine similarity between two modes' displacements:
$$
\rho(\mathbf{d}_a, \mathbf{d}_b) =
\frac{|\,(\sqrt{m}\,\mathbf{d}_a) \cdot (\sqrt{m}\,\mathbf{d}_b)\,|}
     {\|\sqrt{m}\,\mathbf{d}_a\|\,\|\sqrt{m}\,\mathbf{d}_b\|}.
$$
Mass weighting makes "different modes" mean physically different
nuclear motion rather than a permutation of which hydrogens vibrate.

If the filter rejects every IG (rare; occurs when no IG has any
mode along the gross R$\to$P direction), the verifier falls back to
plain bond_overlap top-$k$ on the unfiltered IG pool to guarantee a
non-empty output.

### Hyperparameters

| symbol | value | role |
|---|---|---|
| `min_rxn` | 0.10 | reject IG if best-mode rxn_overlap below this |
| `max_imag` | 2 | reject IG with more imaginary modes than this |
| $w_r$ | 1.0 | rxn_overlap weight in score |
| $w_c$ | 0.2 | core_fraction weight in score |
| $p$ | 0.3 | exponent on $n_{\text{im}}$ in denominator |
| $\alpha$ | 0.7 | mass-weighted-cosine diversity strength |
| $k$ | 2 | top-$k$ output consumed by the live pipeline |

Hyperparameter values were grid-searched on the BGCP benchmark over
a 15-iteration ranker-improvement campaign; ablations are in
Appendix C.

## B.3 Evaluation

### Oracle alignment score

The alignment score between two normal modes is the sign-blind
cosine similarity in 3$N$-Cartesian displacement space:
$$
a(\mathbf{d}_a, \mathbf{d}_b) =
\frac{|\,\mathbf{d}_a \cdot \mathbf{d}_b\,|}{\|\mathbf{d}_a\|\,\|\mathbf{d}_b\|} \in [0, 1].
$$
The absolute value collapses the eigenvector sign ambiguity; the
shared R-frame indexing ensures atom $i$ in $\mathbf{d}_a$ and
$\mathbf{d}_b$ refers to the same chemical atom. The alignment is
called *oracle* when $\mathbf{d}_b = \mathbf{d}_{\text{GT}}$, the
bond_overlap-best imaginary mode at the GT TS structure.

### Aggregations

For an IG $g$ with mode set $M^{(g)}$ we report
$$
\mathrm{oracle\_any}(g) = \max_{m \in M^{(g)}} a(\mathbf{d}_m, \mathbf{d}_{\text{GT}}),
$$
the best alignment achievable from that IG. For a step with IG pool
$\mathcal{G} = \{g_1, \dots, g_{20}\}$,
$\mathrm{oracle\_step} = \max_{g \in \mathcal{G}} \mathrm{oracle\_any}(g)$
is the absolute pass@20 ceiling.

### Pass@$k$

The verifier's pass@$k$ is the cumulative max of
$\mathrm{oracle\_any}$ over its first $k$ ranked IGs:
$$
\mathrm{pass@k}(\text{verifier}) = \max_{g \in \mathrm{verifier\_top\text{-}}k} \mathrm{oracle\_any}(g).
$$

### Uniform-random baseline

We compare against the closed-form expected pass@$k$ of a
no-signal verifier sampling without replacement from the 20-IG
pool. Sorting the per-step oracle_any values
$a_{(1)} \ge a_{(2)} \ge \dots \ge a_{(N)}$, the expectation is
$$
E[\mathrm{pass@k}] = \sum_{r = 1}^{N - k + 1} a_{(r)}\;\frac{\binom{N - r}{k - 1}}{\binom{N}{k}},
$$
since the probability that the smallest sampled rank equals $r$ is
$\binom{N - r}{k - 1} / \binom{N}{k}$.

### Human evaluation

For every step the chemist reviews the verifier's top-2 surfaced
IGs and marks each correct (1) or incorrect (0). A step *passes*
under human judgement if at least one of the two surfaced IGs is
correct. This protocol matches the operational meaning of pass@2:
the chemist receives a useful TS without needing to look beyond the
two candidates.


# Appendix C — design decisions and ablations

## C.1 Why imaginary-mode count $\le 2$

A clean first-order saddle has exactly one imaginary mode, but
finite-difference Hessians on geometries that are slightly off the
saddle often produce a second small imaginary mode. The
$n_{\text{im}} \le 2$ filter removes IGs with three or more
imaginary modes (typically high-order saddles or non-stationary
points) without being so strict as to discard near-saddle structures
that the agent layer produced reliably. Setting `max_imag = 1`
loses several percentage points of pass@2 relative to `max_imag = 2`.

## C.2 Why mass weighting in the diversity penalty only

Replacing the cosine in the diversity penalty with a mass-weighted
cosine has two effects: it makes "physically different motion" the
diversity criterion rather than "different vibration of the
hydrogens", and it reduces visual duplication in the surfaced top-2
on metal-bearing systems where heavy-atom motion is small in
Cartesian terms but large in mass-weighted terms. Mass weighting in
the score itself (replacing bond_overlap, etc., with mass-weighted
analogues) was tested but did not move pass@2 measurably.

## C.3 Bond-identification edge cases

Six steps in the benchmark have `broken = 0 ∧ formed = 0` despite
being valid TS geometries: the WBO change between R and P falls
below $\Delta_{\text{thr}} = 0.5$ on every bond. These cluster in
TEMPO radical chemistry (where "bond change" is partial spin
redistribution rather than a clean σ-bond event), small carbocation
rearrangements, transition-metal cleavage of partial-covalent
bonds, and concerted hydride transfers. For these cases bond_overlap
and core_fraction are vacuous and the verifier falls back to plain
bond_overlap top-$k$. A second cluster of 13 steps has IG bond
identification working but the GT mode is dominated by motion
perpendicular to the identified bond axes (e.g. wags,
pyramidalisation at carbene centres); for these the human chemist
agrees with the IG choice while the alignment score under-credits
it. We discuss this metric / chemistry mismatch in [§ XX] of the
main text.

## C.4 Hyperparameter ablations (summary)

| change | $\Delta$ pass@2 mean | $\Delta$ pass@2 ($\ge 0.7$) |
|---|---|---|
| baseline (clean_v2 as-is) | — | — |
| no diversity penalty (greedy bond_ov only) | $-0.038$ | $-3.8$ pp |
| no rxn_overlap filter | $-0.004$ | $-0.5$ pp |
| no $(1 + w_r r)$ score factor | $-0.001$ | $-0.1$ pp |
| no $(1 + w_c c)$ score factor | $-0.001$ | $-0.1$ pp |
| no $1 / n_{\text{im}}^p$ damping | flat | flat |
| `max_imag = 1` | $-0.010$ | $-3.1$ pp |
| frequency-band filter | flat | flat |

The diversity penalty contributes ~90 % of the verifier's mean
improvement over a plain bond_overlap ranking. The rxn_overlap
filter and the multiplicative $(1 + w_r r)$ / $(1 + w_c c)$ score
factors each add modest, additive corrections.
