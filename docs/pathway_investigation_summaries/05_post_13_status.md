# Post-13 And Later Status

Current status: not investigated with the xTB/Fukui/AAM scorer in this session.

The source structures exist, but there are no cached post-`13` metric tables in
the current xTB test folder.

## Available XYZ Inputs

| structure | path exists? | xTB/WBO cache in current run? | scorer summary? | note |
|---|---:|---:|---:|---|
| `13.xyz` | yes | no | no | next likely target |
| `14.xyz` | yes | no | no | later target |
| `15.xyz` | yes | no | no | later branch/target |
| `16.xyz` | yes | no | no | later target |
| `17.xyz` | yes | no | no | later target |

Relevant source directory:

```text
/Users/yunhengz/Downloads/gold-catalyzed-rearrangement-3acyloxypropynyloxiranes/
```

## What To Run Next

For post-`13`, use the same product-blind workflow:

1. Run xTB single point and `--vfukui` on native `13.xyz` with charge `+1`.
2. Produce:
   - `current13_intrinsic_forward.atom_stats.csv`
   - `current13_intrinsic_forward.formation_scores.csv`
   - `current13_intrinsic_forward.cleavage_scores_by_class.csv`
   - `current13_intrinsic_forward.coupled_move_scores.csv`
   - `current13_intrinsic_forward.summary.md`
3. Keep raw formation ranks before any tags/filters.
4. Use graph distance only as annotation except for already-bonded pairs.
5. Verify later with local AAM `13 -> 14` only after making the product-blind
   post-`13` call.

## Placeholder Metric Table

| metric | value |
|---|---|
| raw formation ranking | not run |
| cleavage ranking | not run |
| coupled ranking | not run |
| AAM verification | not run |
| native-index viewer | not generated |

Read: do not infer post-`13` chemistry from `12` or from known downstream files
until the native `13` descriptor tables exist.
