# Atom alignment and bond breaking / forming detection

## Problem statement

Given two molecular structures — a reactant **R** with $n_R$ atoms and a
product **P** with $n_P$ atoms (we require $n_R = n_P$ and identical
element multisets) — produce:

1. an injective **chemical mapping** $\pi: \{0,\dots,n_R-1\} \to \{0,\dots,n_P-1\}$
   so that atom $i$ of R is the *same physical atom* as atom $\pi(i)$ of P;
2. the set of **broken bonds** $B \subset \binom{[n_R]}{2}$ (atom pairs
   bonded in R but not in P) and **formed bonds** $F$ (bonded in P but
   not in R).

Atom labels are not provided by the input xyzs (xtb output is element-and-coordinate only),
so the mapping must be inferred from chemistry. The resulting
mapping is the basis of every downstream comparison: TS modes, oracle
alignment, and viewer rendering all assume index $i$ refers to the
same atom across R, P, and every transition-state structure.

## Inputs

For each of R and P we obtain from a single xtb GFN2 single-point:

* element list $e \in \Sigma^{n}$ (atomic symbols),
* Cartesian coordinates $\mathbf{x} \in \mathbb{R}^{n \times 3}$,
* Wiberg bond order matrix $W \in \mathbb{R}^{n \times n}$, symmetric,
  $W_{ij} \in [0, \approx 3]$ for typical organic bonds.

The WBO matrix is the load-bearing object: it captures bond
connectivity *and* bond order in one numerical structure that is
robust to small geometric distortions and consistent across the level
of theory.

## Bond graph

Construct an undirected, attributed graph $G = (V, E, w)$ from a
$(e, W)$ pair using a single threshold $w_{\text{floor}} = 0.2$:

$$
V = \{1,\dots,n\}, \quad E = \{(i,j) : W_{ij} \ge w_{\text{floor}} \},
\quad w(i,j) = W_{ij}.
$$

Each node carries its element. The threshold $w_{\text{floor}} = 0.2$
is loose by design: it admits weak / partial bonds (dative,
hyperconjugative, transient interactions) so that the propagation
algorithm sees the same connectivity in R and P even when WBOs differ
mildly.

## Priority-queue alignment

The mapping is grown one **island** at a time. An island is a maximal
connected fragment of $G_R$ that has been chemically matched to a
fragment of $G_P$ with consistent element labels and pairwise WBOs. We
process islands one seed at a time; when growth saturates without
ambiguity, the island is locked.

### Single-island growth

Given a seed atom $s \in V_R$ that is not yet mapped, initialise the
candidate set
$$
\mathcal{C} = \{ \{s \mapsto v\} : v \in V_P \setminus \mathrm{img}(\pi),\;
  e_R(s) = e_P(v) \},
$$
i.e. all not-yet-used P-atoms of the same element as $s$. Maintain a
**fragment** $F \subseteq V_R$ (initially $\{s\}$) and a **max-heap**
$H$ of edges $(u, n)$ with $u \in F$, $n \notin F$, keyed by $w(u, n)$
descending.

Iterate until either (a) $|\mathcal{C}| = 1$ and the candidate is
*set-unique* (defined below), or (b) the heap empties:

1. **Pop** the highest-WBO edge $(u, n)$. Mark it consumed.
2. Try to **extend** every candidate $c \in \mathcal{C}$ to include
   $n$: for each $c$, find P-atoms $v$ that are
   * neighbours in $G_P$ of $c(u)$ for *every* $u' \in F$ that is
     bonded to $n$ in R (intersection of P-neighbour sets);
   * not already in $\mathrm{img}(c)$ or in the global $\mathrm{img}(\pi)$;
   * element-matched to $e_R(n)$;
   * WBO-matched: $|W_R(u', n) - W_P(c(u'), v)| \le \tau$ on every
     bonded edge from $n$ into $F$.
   The extended candidates form $\mathcal{C}'$.
3. **Commit or consume.** If $\mathcal{C}'$ is non-empty, replace
   $\mathcal{C} \leftarrow \mathcal{C}'$, add $n$ to $F$, and push all
   not-yet-consumed outgoing edges of $n$ into $H$. If $\mathcal{C}'$
   is empty, the edge $(u, n)$ is *consumed* without extending the
   fragment — $n$ is reachable from $F$ via this edge but no chemical
   match exists, so we skip to the next pop.
4. If $n$ is already mapped from a previous island (i.e. $n \in
   \mathrm{dom}(\pi)$), perform a **whole-island merge**: extend each
   candidate to include $n$ *and* every other atom of $n$'s
   already-locked island, checking element + WBO consistency on every
   edge that crosses the merge.

The propagation order — highest WBO first — favours the strongest
bonds, which are also the most chemically informative. Weak bonds are
seen later, when more constraints from neighbouring strong bonds are
available to disambiguate.

### Saturation and locking

A candidate set $\mathcal{C} = \{c_1, \dots, c_k\}$ is **set-unique**
if every $c_i$ has the same image set $\mathrm{img}(c_i)$ — the
candidates differ only by permutations within indistinguishable groups
(e.g. the three Hs of a methyl). Set-unique candidates produce
identical island contents, so we lock the island with one representative
mapping. When $|\mathcal{C}| = 1$ at any point we lock immediately.

If the heap empties with $|\mathcal{C}| > 1$ and not set-unique
(genuine residual ambiguity, e.g. enantiomeric C-H configurations
across symmetric scaffolds), we **branch**: each remaining
$c_i$ becomes a separate hypothesis, propagated independently in
later seeds. Branches are deduplicated by their full mapping signature
and capped at `max_branches = 8`.

### Multi-island driver

Iterate single-island growth over every atom in a randomised seed
order; islands accumulate, and later seeds enjoy more locked
neighbours (more constraints, easier disambiguation). Repeat passes
until no new mapping is produced. We run this with $n_{\text{seeds}} =
10$ random permutations to break order-dependent ties and pick the
best ordering by the post-hoc score below.

### Tolerances

| symbol | value | rôle |
|---|---|---|
| $w_{\text{floor}}$ | 0.2 | minimum WBO to admit an edge into the graph |
| $\tau$ | 1.0 | maximum $|W_R - W_P|$ tolerated on a candidate edge |
| $w_{\text{break}}$ | 0.5 | minimum WBO to consider an edge a "bond" for breaking/forming |
| $\Delta_{\text{thr}}$ | 0.5 | minimum $|\Delta W|$ to count as a bond change |

The wide $\tau = 1.0$ accommodates partial bond order changes (e.g.
$\pi$ contributions shifting by ~0.5 between R and P) without breaking
the propagation. Strict matching ($\tau \ll 1$) loses many real
reactions; lax matching ($\tau \gg 1$) collapses chirality and
connectivity.

### Greedy fallback for unmapped atoms

After every seed has been processed, atoms that were never absorbed
into any locked island are paired by a **greedy nearest-element
match**: each unmapped R-atom is paired with the closest still-free
same-element P-atom by Cartesian distance (after a Kabsch alignment
of the locked atoms). This typically settles the few hydrogens that
remained ambiguous.

## Connectivity expansion

Once islands are committed, `expand_mapping` propagates the mapping
across every mapped atom's unmapped neighbours by an element-multiset
match: for each $(u, v) \in \pi$, partition the unmapped R-neighbours
of $u$ and the unmapped P-neighbours of $v$ by element. If counts
agree element-by-element, pair them in arbitrary intra-group order
(symmetric atoms are interchangeable, so any ordering is correct). If
counts differ, leave the atoms unmapped — *that is the reactive
signature*: an atom whose neighbour-element multiset changed between R
and P is touching a broken or formed bond.

## Bond classification

Given the mapping $\pi$, classify every potential bond by its WBO
change:

$$
\text{broken} \iff W_R(i,j) \ge w_{\text{break}} \;\wedge\; W_R(i,j) - W_P(\pi(i), \pi(j)) \ge \Delta_{\text{thr}}
$$
$$
\text{formed} \iff W_P(\pi(i), \pi(j)) \ge w_{\text{break}} \;\wedge\; W_P(\pi(i), \pi(j)) - W_R(i, j) \ge \Delta_{\text{thr}}
$$

Pairs with at least one unmapped endpoint are flagged as broken (if
they were a bond in R) or formed (if a bond in P), since their image
WBO is undefined. The thresholds capture both full bond breaks
($W \approx 1.0 \to 0$) and partial bond-order changes
($W \approx 1.9 \to 1.0$, $\Delta W \approx 0.9$).

The **core atoms** are the union of broken-bond endpoints (already in
R-frame) and formed-bond endpoints (mapped back via $\pi^{-1}$). They
define the localised reactive region used by every downstream metric.

## Scoring multi-seed runs

The PQ run produces one mapping per seed (typically several distinct
mappings across the $n_{\text{seeds}} = 10$ orderings, plus their
branches). We score each by a lexicographic tuple
$$
(\,|B| + |F|,\; v_{\chi},\; -|\,\mathrm{dom}(\pi)|\,)
$$
sorted ascending: prefer fewer broken+formed bonds, then fewer
chirality violations, then more atoms mapped. Lower is better on every
component.

* $|B| + |F|$ favours the simplest reactive change consistent with
  the data — unnecessary "extra" reactive bonds are usually mapping
  artefacts.
* $v_{\chi}$ — the **spectator chirality violation count**: at every
  R-atom of degree $\ge 4$ that does *not* touch a broken or formed
  bond and whose mapped neighbours are all preserved, compare the
  signed determinant $\det[\,\mathbf{r}_{n_1} - \mathbf{r}_u,\,
  \mathbf{r}_{n_2} - \mathbf{r}_u,\, \mathbf{r}_{n_3} - \mathbf{r}_u\,]$
  in R-frame against the analogous determinant in P-frame. A sign flip
  at a non-reactive stereocentre means the mapping has swapped
  enantiomer-equivalent neighbours — a globally consistent but
  chirality-violating choice. Near-coplanar centres (|det| < 0.05) are
  skipped.
* $-|\mathrm{dom}(\pi)|$ rewards more atoms mapped (acts only when the
  earlier components tie).

## Output

`align_from_arrays` returns the best-scored mapping along with $B$,
$F$, the core atoms, and the chirality-violation count. This output
feeds:

* the coordinate-aligned XYZs in
  `Benchmark_Guesses_Coordinate_Aligned_Version/` (every TS reindexed
  to R-frame),
* the per-step viewer payloads in `out/mode_viewer/`,
* the bond-level metrics (`bond_overlap`, `core_fraction`) used by
  the verifier.
