# FusedAlgebra.Topics.ShapeDirections

Theorems: 36. GENERAL 24, SAME 9, EQUIV 2, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS, no low-confidence items.

- `periods_outnumber_curves`: PORT-NOTES claim checked and correct. Original: `∀ g, g < 40 → 2 ≤ g → (6 * g - 6 < g * (g + 1) ↔ 4 ≤ g)` over ℕ by `decide`. Port: `(g : ℕ) (hg : 2 ≤ g) : 6 * g - 6 < g * (g + 1) ↔ 4 ≤ g`. Same ℕ-subtraction form with the bound dropped. True for all g ≥ 2 because (g-2)(g-3) > 0 for g ≥ 4 and both sides are false at g = 2, 3. GENERAL (strengthened).
- `jet_rank`: `det3` became `Matrix.det` of `!![…]`. I checked the rows against det3's argument order: (12,0,0),(0,12,0),(12,3,-1) → -144 and (12,0,0),(12,0,0),(0,12,0) → 0. The same.
- `node_three_points`: the left side is now mathlib's `discrim 1 (-1/3) (-(-1/3)^3)` = b² − 4·1·c. That is the same number as the original expression (-1/27).
- The ℤ statements became CommRing statements; I checked the coefficients in `pseries_*` and `additive_at_*`, and they are unchanged.
- The jets (`Jet`, `comp4`, `comp5`, `E`) stay over ℚ, and their formulas are identical.

Route changes: `glue_commutes` is now the term `multAdd_right_comm`. `periods_outnumber_curves` is proved by omega/nlinarith rather than `decide`. Everything else is `ring`/`linarith` where the original used `grind`, which is not a meaningful change.

New port theorems with no original: 2.

Rung index note: no theorem in this module uses `MultAdd.rung`. Where the docstrings say "rung 1" or "rung s" they mean the fractional height s (s = 0 addition, s = 1 mult-add), the originals' own convention in this topic; that differs from Rungs.lean / `MultAdd.rung`, where mult-add is index 2. Same in original and port, so not a port finding.
