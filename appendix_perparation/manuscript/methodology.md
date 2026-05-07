# Methodology

We describe the four components of our agentic transition-state
discovery pipeline that are evaluated in this work: (i) the benchmark
of elementary steps (Section 1); (ii) the semi-empirical electronic
structure used to assign atomic correspondence and classify bond
events (Section 2); (iii) the priority-queue atom-mapping algorithm
that brings reactant, product, ground-truth TS, and every initial
guess into a shared atom indexing (Section 3); (iv) the verifier
that ranks initial-guess transition-state structures and surfaces a
top-$k$ set to the human chemist (Section 4). Section 5 defines the
oracle alignment score and the pass@$k$ metrics used for evaluation.

## 1. Benchmark dataset

We evaluate on a curated benchmark of $N = 155$ elementary steps
spanning four chemistry domains: organocatalysis (TEMPO radicals,
Morita–Baylis–Hillman, ketenes, carbenes), transition-metal catalysis
(Pd-Suzuki, Ni-Suzuki, V-DODH dehydration, Fe-crosscoupling, Co
silylation, Au alkyl), pericyclic chemistry (cyclobutane ring opening,
Diels–Alder cycloadditions), and main-group / non-metal cases
(carbocation rearrangements, Noyori H$_2$ activation). Each step is a
single elementary reaction $\mathrm{R} \to \mathrm{P}$ with a
literature- or DFT-derived reference transition state (the
"ground truth", GT). Step sizes range from 18 to 149 atoms (median
57, IQR [34, 82]) and contain elements from the first three rows of
the periodic table plus selected 4d/5d transition metals (Pd, Pt, Ru,
Rh, Re, Ir, Au, W). For every step the agent layer
generates 20 initial-guess transition-state structures (IGs); the
benchmark therefore contains 155 reactant geometries, 155 product
geometries, 155 GT-TS geometries, and 3,100 IG-TS geometries.

After deduplication (5 redundant cyclobutane variants and one
duplicate Jackie-series TS removed), the working set is 155 unique
elementary steps. The bond-count distribution per step is unimodal
with mean 1.30 broken bonds and 1.34 formed bonds; the most common
profile is one bond broken and one formed (52 of 155 steps).

## 2. Electronic structure

For every geometry — R, P, GT, and each IG — we compute a single-
point Wiberg bond order (WBO) matrix and a Hessian using xtb at the
GFN2-xTB level of theory. WBO matrices serve as the primary chemical
fingerprint for atom mapping and bond classification (Sections 3–4);
Hessians provide normal-mode displacements used to identify each
TS's reaction coordinate and to score the ranker. We run xtb in a
caching wrapper so that repeated invocations on identical xyz inputs
read the WBO and Hessian from disk; this makes downstream evaluation
scripts effectively zero-cost after a one-time benchmark sweep.

## 3. Atom mapping via priority-queue alignment

A reaction $\mathrm{R} \to \mathrm{P}$ defines an injective mapping
$\pi : V_R \to V_P$ between atomic indices that identifies which
nucleus in R is which in P. The xyz outputs of xtb do not carry such
a mapping, so we infer it from the WBO matrices.

### 3.1 Bond graph construction

For each structure we build an undirected graph $G = (V, E, w)$ from
the elements $e$ and WBO matrix $W$ using a single threshold
$w_{\text{floor}} = 0.2$:
$$
V = \{1, \dots, n\}, \quad
E = \{(i,j) : W_{ij} \ge w_{\text{floor}} \}, \quad
w(i, j) = W_{ij}.
$$
Each node carries its element symbol. The threshold is intentionally
loose so that weak / partial / dative bonds are admitted; this makes
graph topology robust to small WBO fluctuations between R and P.

### 3.2 Single-island growth

The mapping is grown one **island** at a time. An island is a maximal
connected fragment of $G_R$ that has been chemically matched to a
fragment of $G_P$ with consistent element labels and pairwise WBOs.
Given a seed atom $s \in V_R$ that is not yet mapped, we initialise
the candidate mapping set
$$
\mathcal{C} = \{\, \{s \mapsto v\} : v \in V_P \setminus \mathrm{img}(\pi),\;
                                       e_R(s) = e_P(v)\, \}
$$
and a max-heap $H$ of edges $(u, n)$ with $u \in F$ (the fragment so
far) and $n \notin F$, keyed by $w(u, n)$ in descending order. We
iterate until either $|\mathcal{C}| = 1$ and the candidate is
*set-unique* (defined below) or the heap empties:

1. **Pop** the highest-WBO edge $(u, n)$. Mark it consumed.
2. **Try to extend** every candidate $c \in \mathcal{C}$ to include
   $n$: for each $c$ find P-atoms $v$ that are neighbours in $G_P$ of
   every $u' \in F$ that is bonded to $n$ in R, that satisfy the
   element constraint $e_R(n) = e_P(v)$, and whose WBO discrepancy
   $|W_R(u', n) - W_P(c(u'), v)| \le \tau$ on every bonded edge from
   $n$ into $F$ ($\tau = 1.0$).
3. **Commit or consume.** If any extension survives, replace
   $\mathcal{C}$ with the extended set and add $n$ to $F$. If none
   survive, the edge is consumed and the fragment does not grow.

A candidate set is **set-unique** when every $c_i$ has the same
image $\mathrm{img}(c_i)$ — the candidates differ only by permutations
of indistinguishable atoms (e.g. methyl hydrogens). If the heap
empties with multiple non-set-unique candidates, the algorithm
**branches**: each candidate becomes an independent hypothesis to be
propagated by later seeds; branches are deduplicated by mapping
signature and capped at eight.

### 3.3 Multi-seed driver and scoring

We run the single-island growth procedure over ten random seed
orderings and accumulate the resulting branches. For each candidate
mapping $\pi$ we then identify broken and formed bonds via
$$
\text{broken}: W_R(i,j) \ge w_{\text{break}} \;\wedge\; W_R(i,j) - W_P(\pi(i), \pi(j)) \ge \Delta_{\text{thr}},
$$
$$
\text{formed}: W_P(\pi(i), \pi(j)) \ge w_{\text{break}} \;\wedge\; W_P(\pi(i), \pi(j)) - W_R(i, j) \ge \Delta_{\text{thr}},
$$
with $w_{\text{break}} = 0.5$ and $\Delta_{\text{thr}} = 0.5$. We
score each candidate by the lexicographic tuple
$\big(|B| + |F|,\; v_\chi,\; -|\mathrm{dom}(\pi)|\big)$, sorted ascending:
prefer the simplest reactive change, then fewest spectator-stereocentre
chirality violations $v_\chi$, then most atoms mapped. The best-scored
mapping is taken as $\pi$ for downstream analysis.

The chirality term $v_\chi$ counts non-reactive degree-$\ge 4$ centres
where the signed determinant of the first three neighbour displacement
vectors flips sign between R-frame and P-frame; near-coplanar centres
are skipped to avoid sign noise.

The same algorithm aligns every TS structure (GT and the 20 IGs) to
R-frame, producing a unified atom indexing where index $i$ is the
same chemical atom across all 23 structures of a given step.

## 4. Verifier

For each IG we extract every normal mode as a per-atom Cartesian
displacement vector $\mathbf{d} \in \mathbb{R}^{n_R \times 3}$ at
frequency $\omega$. The verifier ranks the 20 IGs and surfaces the
top $k$ to the chemist; in production $k = 2$.

### 4.1 Per-mode features

Three scalar features are computed for every mode of every IG. All
are sign-blind, take values in $[0, 1]$, and use the unified atom
indexing.

**Bond-direction overlap.** Construct the bond reaction vector at the
IG's TS coordinates $\mathbf{x}^{TS}$:
$$
\mathbf{V}_i = \!\!\sum_{(i,j) \in B} -\hat{\mathbf{u}}_{ij}^{TS}
        \;+\!\sum_{(j,i) \in B}  \hat{\mathbf{u}}_{ji}^{TS}
        \;+\!\sum_{(i,j) \in F}  \hat{\mathbf{u}}_{ij}^{TS}
        \;+\!\sum_{(j,i) \in F} -\hat{\mathbf{u}}_{ji}^{TS},
$$
where $\hat{\mathbf{u}}_{ij}^{TS}$ is the unit vector from atom $i$
to atom $j$ at the TS. The sign convention encodes "broken-bond
endpoints separating, formed-bond endpoints approaching" so that a
true concerted reaction mode aggregates contributions coherently. Then
$$
\mathrm{bond\_overlap}(\mathbf{d}) =
\frac{|\,\mathbf{d} \cdot \mathbf{V}\,|}{\|\mathbf{d}\|\,\|\mathbf{V}\|}.
$$

**Reaction-coordinate overlap.** Compute the per-atom R$\to$P
displacement after Kabsch-aligning P (in R-frame) to R, restrict to
core atoms $C$ (atoms touching any broken or formed bond), and
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
     {\sum_{i = 1}^{n_R} \|\mathbf{d}_i\|^2}\;\in [0, 1].
$$

### 4.2 Ranking algorithm

For each IG we first restrict to the imaginary modes
$M_{\text{im}} = \{\, m : \omega_m < 0\, \}$. We pick the mode of
highest bond_overlap within the IG:
$$
\hat{m} = \arg\max_{m \in M_{\text{im}}} \mathrm{bond\_overlap}(m).
$$
We then apply two structural priors as a hard filter:
$$
1 \le n_{\text{im}} \le 2 \quad\wedge\quad \mathrm{rxn\_overlap}(\hat{m}) \ge 0.10.
$$
The first reflects the canonical TS condition (a true first-order
saddle has exactly one imaginary mode; we relax to two); the second
removes IGs whose chosen mode does not even point along the gross
R$\to$P direction. The filter eliminates approximately 35 % of IGs
on average.

The score of the surviving IG is
$$
S(\hat{m}, n_{\text{im}}) =
\frac{\mathrm{bond\_overlap}(\hat{m})\,
      (1 + w_r \,\mathrm{rxn\_overlap}(\hat{m}))\,
      (1 + w_c \,\mathrm{core\_fraction}(\hat{m}))}
     {n_{\text{im}}^{\,p}}
$$
with $w_r = 1.0$, $w_c = 0.2$, $p = 0.3$. The multiplicative form
ties the score to bond_overlap (the dominant single signal), adds
modest gain factors for rxn_overlap and core_fraction, and applies a
soft preference for clean first-order saddles via $n_{\text{im}}^{p}$.

### 4.3 Diversity-aware top-$k$ selection

We greedily pick up to $k$ IGs. Maintain the selected set
$\mathcal{S}$ (initially empty). At each step, pick
$$
c^{\star} = \arg\max_{c \in \mathcal{C} \setminus \mathcal{S}}
\;S(c) \cdot \prod_{s \in \mathcal{S}} \max\!\big(0,\; 1 - \alpha\,\rho(c, s)\big),
$$
where $\alpha = 0.7$ is the diversity strength and $\rho$ is the
mass-weighted cosine similarity:
$$
\rho(\mathbf{d}_a, \mathbf{d}_b) =
\frac{|\,(\sqrt{m}\,\mathbf{d}_a) \cdot (\sqrt{m}\,\mathbf{d}_b)\,|}
     {\|\sqrt{m}\,\mathbf{d}_a\|\,\|\sqrt{m}\,\mathbf{d}_b\|},
$$
with $m_i$ the atomic mass of atom $i$. The mass weighting makes
"different modes" mean *physically different nuclear motion* rather
than a permutation of which hydrogens vibrate. The penalty drives
later picks away from the basin of earlier picks: when two IGs
converged to the same TS structure their best-bond_overlap modes are
near-parallel, $\rho \approx 1$, and the penalty multiplies the score
toward zero. The diversity term contributes 90 % of the verifier's
mean improvement over a plain bond_overlap ranking at $k = 2$.

If the filter rejects every IG (rare, occurs when no IG has any mode
pointing along the gross reaction direction), the verifier falls back
to plain bond_overlap top-$k$ on the unfiltered IG pool to guarantee
a non-empty output.

### 4.4 Hyperparameters

The verifier hyperparameters are summarised below; values were
selected by grid search on the BGCP benchmark over a 15-iteration
ranker-improvement campaign.

| symbol | value | role |
|---|---|---|
| `min_rxn` | 0.10 | reject IG if best-mode `rxn_overlap` below this |
| `max_imag` | 2 | reject IG if it has more imaginary modes than this |
| `w_r` | 1.0 | rxn_overlap weight in score |
| `w_c` | 0.2 | core_fraction weight in score |
| `p` | 0.3 | exponent on $n_{\text{im}}$ in denominator |
| $\alpha$ | 0.7 | mass-weighted-cosine diversity strength |
| $k$ | 2 | top-$k$ output consumed by the live pipeline |

## 5. Evaluation

### 5.1 Oracle alignment score

The alignment score between two normal modes is the sign-blind
cosine similarity in the $3 n_R$-dimensional Cartesian displacement
space, after both modes have been re-indexed into the shared R-frame:
$$
a(\mathbf{d}_a, \mathbf{d}_b) =
\frac{|\,\mathbf{d}_a \cdot \mathbf{d}_b\,|}{\|\mathbf{d}_a\|\,\|\mathbf{d}_b\|}\;\in [0, 1].
$$
The absolute value collapses the eigenvector sign ambiguity (a
vibration mode and its negation are physically equivalent). Because
both vectors live in the shared R-frame, atom $i$ in $\mathbf{d}_a$
and $\mathbf{d}_b$ refers to the same chemical atom; the cosine
captures whether the two modes describe the same nuclear-displacement
direction at every atom.

We call this score the *oracle* alignment when one of the two modes
is the GT reaction mode $\mathbf{d}_{\text{GT}}$ — the
bond_overlap-best imaginary mode at the GT TS structure, computed
through the same xtb GFN2 Hessian pipeline as the IGs. The qualifier
"oracle" indicates that this score uses the GT label and so cannot
be used to rank IGs in production; it is the evaluation target the
unsupervised verifier is approximating.

### 5.2 Aggregations

For an IG $g$ with modes $M^{(g)}$ we report
$$
\mathrm{oracle\_any}(g) = \max_{m \in M^{(g)}} a(\mathbf{d}_m, \mathbf{d}_{\text{GT}}),
$$
the best alignment achievable from that IG given an omniscient choice
of mode. For a full step with IG pool $\mathcal{G} = \{g_1, \dots, g_{20}\}$,
$$
\mathrm{oracle\_step} = \max_{g \in \mathcal{G}} \mathrm{oracle\_any}(g)
$$
is the absolute pass@20 ceiling on that step.

### 5.3 Pass@$k$ metrics

The verifier's $\mathrm{pass}@k$ is the cumulative max of
$\mathrm{oracle\_any}$ over the verifier's first $k$ ranked IGs. For
a side-by-side comparison, we also report a uniform-random baseline:
the closed-form expected max of $k$ IGs sampled without replacement
from a pool of $N = 20$,
$$
E_{\text{rand}}[\max_k] =
\sum_{r=1}^{N - k + 1} a_{(r)}\;\frac{\binom{N - r}{k - 1}}{\binom{N}{k}},
$$
where $a_{(1)} \ge a_{(2)} \ge \dots \ge a_{(N)}$ are the sorted
$\mathrm{oracle\_any}$ values for that step. This gives an exact
"no-signal" baseline against which any verifier improvement can be
measured.

### 5.4 Human evaluation

The verifier's top-2 surfaced for every step is assessed by a
human chemist. For each surfaced IG, the chemist marks it as
chemically correct (1) or incorrect (0). A step *passes* under
human judgement if at least one of the two surfaced IGs is correct,
matching the operational meaning of pass@2: the chemist receives a
useful TS without needing to look beyond the two candidates.
