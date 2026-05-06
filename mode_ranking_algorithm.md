# Mode-ranking algorithm

How we score every vibrational mode of a TS structure and pick the one
that best represents the reaction coordinate.

## Inputs

For one (BGCP step, TS structure) pair:

- `R`, `P`           — reactant / product structures (xtb GFN2 single-point cached in `work_modes/<step>/{R,P}/`)
- `TS`               — transition state structure with cached Hessian at `work_modes/<step>/hess_<label>/g98.out` (Gaussian-format normal modes)
- alignment outputs:
  - `mapping_RP`     — R-atom → P-atom (PQ alignment)
  - `mapping_RT`     — R-atom → TS-atom
  - `core_atoms`     — R-frame atoms that touch a broken or formed bond
  - `broken_R`       — list of `(i, j)` pairs in R-frame for broken bonds
  - `formed_R`       — list of `(i, j)` pairs in R-frame for formed bonds (mapped from P-frame via `inv(mapping_RP)`)

All structures are reindexed into the **R atom-index frame** so atom
`i` refers to the same chemical atom across R, TS, and P (see
`align_bgcp_coords.py`).

The Hessian's normal-mode displacements `modes_TS` (n_modes × n_TS_atoms × 3) are reindexed the same way:

```
modes_R[m, r] = modes_TS[m, mapping_RT[r]]
```

Now `modes_R[m, i]` is the Cartesian displacement of *the chemical atom that's R-atom i* in mode `m`.

## The three candidate metrics

We compute three per-mode scalars; only `bond_overlap` is used to rank
and pick defaults, but all three are kept for diagnostics in CSVs and
the viewer.

### 1. `core_fraction` (legacy, weak)

$$
\text{core\_fraction}_m = \frac{\sum_{i \in \text{core}} \|\mathbf{d}_{m,i}\|^2}{\sum_{i=1}^{N} \|\mathbf{d}_{m,i}\|^2}
$$

Fraction of the mode's total Cartesian "kinetic energy" sitting on core
atoms. Range `[0, 1]`. Insensitive to *direction* — a localized wag of
a single core atom (e.g. an out-of-plane terminal-N flap) scores near
1.0 even though the bond is unchanged. **Discarded as a primary
ranker.**

### 2. `rxn_overlap` (intermediate, finer)

Build the reaction-coordinate displacement per atom by Kabsch-aligning P
to R then taking $\boldsymbol{\Delta}_i = \mathbf{r}^P_{m(i)} - \mathbf{r}^R_i$ for mapped atoms (zero elsewhere). Restrict to core atoms and project:

$$
q_m = \frac{|\,\mathbf{d}_m \cdot \hat{\boldsymbol{\Delta}}^{(\text{core})}\,|}{\|\mathbf{d}_m\|}
$$

where $\hat{\boldsymbol{\Delta}}^{(\text{core})}$ is zero outside the core
atoms and unit-norm over them. Range `[0, 1]`.

Penalizes:

- core motion misaligned with the *integrated* R→P displacement (numerator drops),
- mode amplitude wasted on spectator atoms (denominator inflates).

Rewards modes whose core motion follows the gross R→P trajectory.

Weakness: $\boldsymbol{\Delta}$ is the *finite* R-to-P displacement,
which mixes the reaction coordinate with all the spectator
rearrangements that happen along the path (methyl rotations,
conformational drift). The TS imaginary mode is a *local* eigenvector
at the saddle and need not match the long-range $\boldsymbol{\Delta}$
direction even if it's the correct reaction mode.

### 3. `bond_overlap` (current default — sharper)

Skip the gross R→P direction entirely; project the mode onto the
direction that *changes the breaking/forming bond lengths*.

For TS coordinates `xyz_TS` (in R-atom-index order), build a per-atom
"reaction direction" $\mathbf{V}$:

For each broken bond $(i, j)$ — atoms should move **apart** along the bond:

$$
\mathbf{V}_i \mathrel{-{=}} \frac{\mathbf{r}_j - \mathbf{r}_i}{\|\mathbf{r}_j - \mathbf{r}_i\|}, \quad
\mathbf{V}_j \mathrel{+{=}} \frac{\mathbf{r}_j - \mathbf{r}_i}{\|\mathbf{r}_j - \mathbf{r}_i\|}
$$

For each formed bond $(i, j)$ — atoms should move **together** along the bond:

$$
\mathbf{V}_i \mathrel{+{=}} \frac{\mathbf{r}_j - \mathbf{r}_i}{\|\mathbf{r}_j - \mathbf{r}_i\|}, \quad
\mathbf{V}_j \mathrel{-{=}} \frac{\mathbf{r}_j - \mathbf{r}_i}{\|\mathbf{r}_j - \mathbf{r}_i\|}
$$

Sign convention is chosen so that a true concerted reaction mode
accumulates contributions coherently (every bond simultaneously moving
in the right direction). Then

$$
\boxed{\;\;\text{bond\_overlap}_m \;=\; \frac{|\,\mathbf{d}_m \cdot \hat{\mathbf{V}}\,|}{\|\mathbf{d}_m\|}\;\;}
$$

Range `[0, 1]`. Penalizes:

- modes whose core motion is **perpendicular to bond axes** (a wag that doesn't change bond lengths),
- modes that **anti-correlate** the broken/formed pieces (bond A breaks while bond B fails to form on the same beat),
- modes that waste amplitude on **spectator atoms** (denominator inflates).

The sign-coherent construction means a fully concerted asynchronous
multi-bond reaction can still score `1.0` if its mode points exactly
along $\hat{\mathbf{V}}$; non-concerted modes with cancelling pieces
score lower automatically.

## Worked example: pr13.Cyclobutane_JOC2023_TS-CD_step1, GT

A cyclobutane-ring rearrangement with 3 broken + 1 formed bonds; core
atoms are 2 carbons and 2 nitrogens (R-frame indices 0, 1, 8, 42).

| mode | freq (cm⁻¹) | description (visual) | core_fraction | rxn_overlap | **bond_overlap** |
|---|---|---|---|---|---|
| 0 | −294.71 | concerted ring opening | 0.868 | 0.445 | **0.555** ← picked |
| 1 | −64.82  | one N flapping out of plane | 0.900 ❌ | 0.427 | **0.026** |

`core_fraction` ranks the wag *higher* (0.900 vs 0.868) — wrong, because
88 % of the wag's motion is on a single nitrogen which happens to be a
core atom. `rxn_overlap` correctly orders them but with only a 4 %
margin (1.04× ratio). `bond_overlap` gives a **21× ratio** because the
wag's motion is essentially perpendicular to every reaction bond —
exactly the kind of signal the user is asking the metric to express.

## Default selection

For each TS:

```
imag_indices = {m : freq[m] < 0}
default = argmax over imag_indices of bond_overlap   (if any imag modes)
        else argmax over all modes of bond_overlap   (rare; should not happen at a real TS)
```

Tiebreaks (used in CSV ordering): `bond_overlap`, then `rxn_overlap`, then `core_fraction`.

## Implementation map

| component | file | function |
|---|---|---|
| compute V from broken/formed bonds at TS coords | `analyze_core_modes.py` | `bond_reaction_vector` |
| compute bond_overlap per mode | `analyze_core_modes.py` | `bond_overlap_per_mode` |
| compute rxn_overlap (R→P direction, core-restricted) | `analyze_core_modes.py` | `rxn_overlap_per_mode` |
| Δ vector with Kabsch-aligned P | `analyze_core_modes.py` | `reaction_coord_delta`, `kabsch` |
| reindex modes to R-frame | `analyze_core_modes.py` | `reindex_modes_to_R` |
| reindex TS coords to R-frame | `align_bgcp_coords.py` | `reindex_to_R_frame` |
| consume reaction bonds for V (signed) | `analyze_core_modes.py:process_step` | inline |
| viewer ranking + default selection | `build_mode_viewer.py` | `build_step_data` |
| flat-view ranking (GT + top-2 IG) | `build_flat_view.py` | `build_flat_payload`, `best_imag_mode` |

## Why this rewards "concerted reaction modes" specifically

A vibrational mode is a 3N-Cartesian unit eigenvector. Three desirable
properties of a reaction mode at a TS:

1. **Localized on the core**: the displacement should live mostly on
   the atoms whose bonding is changing.
2. **Aligned with bond axes**: the relative motion of bonded atoms
   should change the bond length, not just rotate the molecule or
   translate spectators.
3. **Concerted phase**: all the broken bonds should be lengthening
   *together* and all the formed bonds shortening *together* in the
   same half-period.

`bond_overlap` is the smallest functional that captures all three:

- (1) is enforced because $\mathbf{V}$ is non-zero only at core atoms,
  so $\hat{\mathbf{V}}$ has no spectator components for $\mathbf{d}_m$
  to project onto. A spectator-heavy mode has large $\|\mathbf{d}_m\|$ in
  the denominator without a matching numerator → low score.
- (2) is enforced because $\hat{\mathbf{V}}$ points *along bond axes*,
  not arbitrary directions in space. A pure wag (perpendicular to
  bonds) has zero numerator → zero score.
- (3) is enforced by the **sign convention** in $\mathbf{V}$:
  contributions add coherently for a concerted mode and cancel for an
  asynchronous one. A mode that breaks one bond but does *not* form
  the partner bond at the same time picks up only the partial signal.

`core_fraction` only addresses (1); `rxn_overlap` addresses (1)+(3) but
not (2); `bond_overlap` addresses all three.

## Limitations / caveats

- **No mass weighting.** `bond_overlap` uses raw Cartesian
  displacements. A "true" reaction-coordinate metric would mass-weight
  ($m_i \cdot \|\mathbf{d}_i\|^2$). This is straightforward to add via
  per-element atomic masses, but in practice GFN2 modes are dominated
  by heavy-atom motion already, so the unweighted score and the
  mass-weighted score correlate strongly.

- **Sign of V depends on broken/formed labelling.** The labelling comes
  from the R↔P alignment + `classify_bonds` policy (broken = WBO_R ≥
  0.5 ∧ ΔWBO ≥ 0.5; formed mirror). If alignment puts a bond on the
  wrong side of the broken/formed boundary, $\mathbf{V}$ encodes the
  wrong direction. In practice this is rare given the aligned shared
  atom-index frame, but worth verifying for cases where the bond list
  looks suspect.

- **TS coordinates used for V.** $\hat{\mathbf{V}}$ is built from TS
  bond directions. For an extreme early/late TS where a "broken" bond
  is barely shorter than its R-equilibrium length (or a "formed" bond
  barely longer than its P-equilibrium length), the TS bond axis can
  differ slightly from the R or P axis. The score is robust to this
  because we only need the *direction*, but at extreme geometries the
  metric may slightly under-reward.

- **No explicit penalty for non-imaginary modes.** A real-frequency
  mode that happens to project well onto $\mathbf{V}$ would score high
  too. We work around this by restricting the *default-pick* logic to
  imaginary modes; the table just shows the ranks.

## Pipeline

```
work_modes/<step>/{R, P, sp_*, hess_*}/   ← cached xtb output (no recompute)
   │
   ├── load (R, P, TS) elements/coords/wbo
   ├── PQ alignment R↔P, R↔TS              ← (alignment recompute, no xtb)
   ├── delta_RP = Kabsch(P-R, mapped)
   ├── modes_R = reindex(modes_TS, mapping_RT)
   ├── core_R   = atoms touching broken/formed
   ├── V        = bond_reaction_vector(TS_coords_in_R, broken_R, formed_R)
   │
   ├── core_fraction  = ||d_core||² / ||d||²
   ├── rxn_overlap    = |d · Δ̂_core| / ||d||
   ├── bond_overlap   = |d · V̂|      / ||d||
   │
   ├── per-step CSV    out/mode_analysis/<step>.csv
   ├── per-step viewer out/mode_viewer/<step>.html
   ├── index           out/mode_viewer/index.html
   └── flat view       out/mode_viewer/flat_view.html
```

The patch scripts (`patch_bond_overlap.py`, `patch_mode_index.py`,
`patch_csv_from_html.py`) make the metric retrofittable without
re-running alignment: every input they need is already embedded in the
per-step HTMLs from a prior build.
