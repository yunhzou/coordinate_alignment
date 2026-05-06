# GT-quality failure modes — the 19 steps with `gt_b < 0.1`

Across 160 BGCP steps, 19 have GT picked-mode `bond_overlap < 0.1`,
i.e. the ground-truth reaction mode itself fails our ranking metric.
These steps drag the apparent ranker performance down (gap to oracle
0.30 on these vs 0.14 on the other 141).

Discriminator: **`gt_core_fraction`** — the fraction of GT mode mass
on identified core atoms. If `gt_c ≥ 0.5`, bonds/core were identified
correctly and the mode genuinely concentrates on those atoms; the
issue is purely that `bond_overlap` (projection on bond axes) cannot
capture the motion's geometry. If `gt_c < 0.3`, the bonds/core
themselves were misidentified.

---

## Failure Mode 1: Bond identification wrong (14 steps)

`gt_c < 0.3` — neither bonds nor the resulting core atoms reflect the
real reaction site.

### 1a — Bonds entirely missed: `broken=0 AND formed=0` (6 steps)

The PQ alignment found no WBO change above threshold between R and P.

| step | freq | n_imag |
|---|---|---|
| `pr1.tempo_ts1` | −75 | 1 |
| `pr1.tempo_ts5` | −62 | 1 |
| `pr16.carbocation_ts2` | +76 | 0 |
| `pr16.carbocation_ts4` | −228 | 1 |
| `pr4.Suzuki-Brcleavage_ts1.Br-cleavage-B` | −279 | 1 |
| `pr5.Noyori_ts5H2` | −992 | 3 |

`pr16.carbocation_ts2` is not a TS at all (positive freq, n_imag=0).
The other five are real TS structures whose R↔P WBO change sits below
the detection threshold. Probable causes:
- Pd–Br / metal–ligand bond cleavage (Suzuki) — partial-covalent bond,
  WBO compressed below the threshold.
- H₂ activation (Noyori) — H–H σ-bond breaking with concurrent metal
  binding, but R and P may both have similar WBO patterns if the
  electronic redistribution is subtle.
- TEMPO low-freq modes (ts1, ts5) — these are spin-redistribution
  modes that don't involve geometric bond change.

### 1b — Bonds detected but wrong (low `gt_c`, 8 steps)

| step | gt_b | gt_c | broken / formed | freq | comment |
|---|---|---|---|---|---|
| `pr7.V.dodh_ts1112-triplet`  | 0.03 | 0.07 | 1 / 1 | −20 | barely a TS, n_imag=4 |
| `pr7.V.dodh_ts32`            | 0.04 | 0.05 | 2 / 1 | −25 | barely a TS |
| `pr7.V.dodh_ts910`           | 0.05 | 0.07 | 2 / 2 | −906 | clean TS, wrong bonds |
| `pr7.V.dodh_ts910-water`     | 0.04 | 0.05 | 2 / 2 | −824 | clean TS, wrong bonds |
| `pr1.tempo_ts2`              | 0.05 | 0.29 | 2 / 2 | −109 | over-detected |
| `pr1.tempo_ts8`              | 0.08 | 0.24 | 2 / 2 | −493 | over-detected |
| `pr16.carbocation_ts9`       | 0.09 | 0.30 | 1 / 1 | −731 | sub-threshold dW |
| `pr3.Suzuki.Ni_ts7.OA`       | 0.08 | 0.26 | 1 / 0 | −37 | OA missing formed bonds |

**Patterns:**
- **V-DODH (4 steps):** dative V–O bonds; xtb WBO is calibrated for
  closed-shell organic bonds and compresses metal–ligand WBO. The
  WBO-difference threshold flags spurious changes and misses real
  ones simultaneously.
- **TEMPO ts2/ts8:** `broken=2 formed=2` is too many for an H-atom
  transfer (should be 1/1). Over-detection from spin-redistribution
  WBO shifts that aren't actual bond events.
- **Suzuki Ni OA:** `formed=0` is wrong — oxidative addition forms
  Ni–C and Ni–Br bonds that the threshold misses (partial covalent).

---

## Failure Mode 2: Bond ID is fine, `bond_overlap` is the wrong metric (5 steps)

`gt_c ≥ 0.5` — the GT mode genuinely concentrates on the
identified core atoms (so bonds and core were found correctly). Yet
`gt_b < 0.1` because the mode displacement is *perpendicular to the
bond axes*. These are reactions where the rate-limiting motion is a
wag, pyramidalization, rotation, or angle change rather than a bond
stretch.

| step | gt_b | gt_c | freq | broken / formed | likely motion |
|---|---|---|---|---|---|
| `pr9.carbene.rearr_ts13a`    | 0.07 | 0.79 | −374 | 1 / 2 | carbene insertion / pyramidalization |
| `Jackie_TS_17`               | 0.04 | 0.72 | −380 | 1 / 2 | rotation / wag |
| `pr1.tempo_ts6`              | 0.04 | 0.66 |  −98 | 2 / 2 | low-freq wag |
| `pr8.ketenes_ts2prime`       | 0.08 | 0.61 | −294 | 0 / 1 | ketene π-cycloaddition |
| `Jackie_TS_18`               | 0.01 | 0.51 | −410 | 1 / 1 | angle bend |

**Why bond_overlap fails here:**
`bond_overlap` projects the mode onto unit vectors along the
broken/formed bond axes. For a stretch-driven reaction (most σ-bond
breaking), the mode moves the two bonded atoms toward/away along that
axis — strong projection. For a pyramidalization at carbene, a
keto-enol-style angle bend, or ring-puckering at a TS, the atoms move
*orthogonal* to the bond they're "on", giving near-zero projection
even though the right atoms are moving.

`pr8.ketenes_ts2prime` is also π-bond-driven: the σ skeleton barely
moves while the π-system reorganizes, so atomic displacements are
small along σ axes.

---

## Implications

For the 14 Mode-1 steps, fixing the ranker is hopeless without first
fixing bond identification. Candidate improvements:
- Lower `dW` threshold for transition-metal / radical systems.
- Augment WBO with bond-distance-change criterion (heuristic R↔P
  Δd cutoff per element pair).
- Spin-aware xtb (open-shell with multiplicity from input).

For the 5 Mode-2 steps, `bond_overlap` is structurally inappropriate.
Options:
- Replace with `core_fraction × |freq|` for low-frequency wag/bend
  modes (already part of our toolkit, just up-weight when bond_overlap
  is degenerate).
- Add an *angle-bend* signal: project mode onto changes in select
  bond angles (not lengths) during R→P alignment.
- Use `dwbo_overlap` (per-atom WBO-environment direction) which
  doesn't require bond-axis projection.

The 5 Mode-2 cases are also evidence that on routine reactions the
existing core_fraction is a perfectly good "right atoms?" signal — it
simply cannot be combined with bond_overlap when the latter is
geometric noise. An adaptive ranker that switches scoring rule based
on `core_fraction / bond_overlap` ratio could recover these.
