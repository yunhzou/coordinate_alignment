# Reactivity Descriptors For Mechanism Analysis

This note lists electronic descriptors worth keeping alongside the existing
WBO caches.  The goal is not to turn the pipeline into a full quantum-chemistry
analysis package; it is to preserve cheap, atom-aligned signals that help
explain why a proposed mechanism is plausible.

The current pipeline already has the most important structural signal:

- Atom labels and coordinates for each R, P, TS, GT, or IG structure.
- xTB Wiberg/Mayer bond orders in `wbo`.
- Per-event WBO changes after R/P/TS mapping.

The proposed extension is to keep atom-level and molecule-level reactivity
properties from the same or closely related xTB calculations.

## Descriptor Set

| Scope | Descriptor | Meaning | xTB source | Keep? |
| --- | --- | --- | --- | --- |
| Atom | Partial charge | Electron-rich vs electron-poor atom proxy. Useful but not enough by itself. | `charges` file / population output | Yes |
| Atom | `f_plus` Fukui index | Susceptibility to nucleophilic attack; local electron-acceptor/electrophilic-site signal. | `xtb --vfukui` stdout | Yes |
| Atom | `f_minus` Fukui index | Susceptibility to electrophilic attack; local electron-donor/nucleophilic-site signal. | `xtb --vfukui` stdout | Yes |
| Atom | `f_zero` Fukui index | Radical attack susceptibility. | `xtb --vfukui` stdout | Optional |
| Bond | WBO | Existing bond-order / partial-bond signal. | `wbo` file | Already kept |
| Bond/event | `delta_wbo` | Bond formation, cleavage, and electronic reorganization between mapped atoms. | derived from R/P/TS WBOs | Already kept |
| Molecule | HOMO-LUMO gap | Global softness / low-energy electronic response proxy. | xTB summary output | Yes |
| Molecule | Vertical IP | Global electron-donor resistance. | `xtb --vip` / `--vipea` | Optional |
| Molecule | Vertical EA | Global electron-acceptor tendency. | `xtb --vea` / `--vipea` | Optional |
| Molecule | Global electrophilicity index | Whole-molecule electrophilicity / Lewis acidity proxy from IP and EA. | `xtb --vomega` | Optional |
| Molecule | Dipole / polarizability | Environment and long-range electrostatic context. | `--dipole`, `--alpha` | Optional |

For the mechanism-ranking use case, the best first-pass addition is:

```text
per atom:
  charge
  f_plus
  f_minus
  f_zero

per molecule:
  homo_lumo_gap_ev

per mapped event:
  WBO_R
  WBO_TS
  WBO_P
  delta_WBO_RP
  TS_progress
```

## Interpretation

Use charges as context, not as the main donor/acceptor score.  A negative atom
is often electron-rich, but charge alone does not tell whether that atom is the
site most willing to donate density in the specific molecular environment.

Use Fukui indices for local reactivity:

- High or otherwise extreme `f_plus` marks atoms that are better electron
  acceptors, i.e. atoms more susceptible to nucleophilic attack.
- High or otherwise extreme `f_minus` marks atoms that are better electron
  donors, i.e. atoms more susceptible to electrophilic attack.
- `f_zero` is useful for radical-like mechanisms but can be ignored for closed
  shell polar mechanisms until needed.

xTB's documentation warns that Fukui indices depend on occupation numbers and
population analysis, so they should be treated as relative descriptors within a
consistent calculation protocol rather than absolute transferable constants.

For bond-making and bond-breaking, combine these signals:

```text
forming bond A-B:
  donor-side evidence   = high f_minus on donor atom
  acceptor-side evidence = high f_plus on acceptor atom
  progress evidence     = WBO_TS between WBO_R and WBO_P

breaking bond A-B:
  electronic strain     = large delta_WBO_RP
  charge separation     = charge change on mapped atoms
  leaving/accepting role = f_minus/f_plus pattern near the event
```

This lets the viewer answer questions such as:

- Does the atom receiving the new bond also look like a local acceptor?
- Does the atom donating into the forming bond look like a local donor?
- Do charge changes agree with the proposed electron flow?
- Are the largest WBO changes located on atoms with strong reactivity signals?

## Cache Layout

Keep these as sidecar files in each existing xTB cache directory:

```text
R/
  molecule.xyz
  wbo
  charges
  xtb.stdout
  reactivity.json
  vfukui.stdout

P/
  molecule.xyz
  wbo
  charges
  xtb.stdout
  reactivity.json
  vfukui.stdout

sp_iter<N>/
  molecule.xyz
  wbo
  charges
  xtb.stdout
  reactivity.json
  vfukui.stdout
```

Suggested `reactivity.json` shape:

```json
{
  "method": "GFN2-xTB",
  "charge": 0,
  "multiplicity": 1,
  "atom_properties": [
    {
      "index": 1,
      "element": "C",
      "charge": -0.123,
      "f_plus": 0.041,
      "f_minus": 0.182,
      "f_zero": 0.112
    }
  ],
  "molecule_properties": {
    "homo_lumo_gap_ev": 4.21
  },
  "warnings": []
}
```

If global IP, EA, or electrophilicity are added later, store them in the same
`molecule_properties` object and record their method separately if they come
from a different xTB route.

## Implementation Notes

xTB run types are mutually exclusive in the command-line interface, so Fukui,
IP/EA, and electrophilicity calculations should be separate cached jobs rather
than extra flags on the current `--sp` command.

Recommended first implementation:

1. Extend the xTB adapter to parse `charges`.
2. Parse HOMO-LUMO gap from `xtb.stdout`.
3. Add an optional `--vfukui` cache fill into a separate subdirectory or with
   namespaced output files.
4. Parse the final Fukui table into atom-aligned `f_plus`, `f_minus`, and
   `f_zero` arrays.
5. Merge the parsed properties into `reactivity.json`.
6. Include mapped atom-property changes in Stage 1/2 JSON and the HTML viewer.

For charged or open-shell systems, all descriptor runs must use the same charge
and multiplicity convention as the WBO-generating run.  The pipeline already
passes total charge and multiplicity into xTB cache filling; the descriptor
adapter should reuse that path.

## What This Is Not

These descriptors are not true NBO donor-acceptor stabilization energies.  xTB
can give WBOs, charges, Fukui indices, orbital gaps, and global IP/EA-style
quantities, but it does not directly produce NBO second-order perturbation
terms such as:

```text
LP(O) -> sigma*(C-X), E(2) = ...
```

If that level of donor-acceptor interpretation is needed later, use the xTB
geometry as the structure source and run a sidecar ORCA, Gaussian, Q-Chem, or
ADF calculation with NBO/GenNBO enabled.  Those results should be stored as a
separate backend because their basis, density, and population model are
different from the xTB descriptors above.

## References

- xTB properties and WBO output:
  <https://xtb-docs.readthedocs.io/en/latest/properties.html>
- xTB single-point descriptors, vertical IP/EA, global electrophilicity, and
  Fukui indices:
  <https://xtb-docs.readthedocs.io/en/latest/sp.html>
- xTB command-line option reference:
  <https://github.com/grimme-lab/xtb/blob/main/man/xtb.1.adoc>
