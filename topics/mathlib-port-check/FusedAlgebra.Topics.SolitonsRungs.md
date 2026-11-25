# FusedAlgebra.Topics.SolitonsRungs

Theorems: 26. GENERAL 22, SAME 2, EQUIV 1, MATHLIB 1. No WEAKER / DIFFERENT / MISSING / VACUOUS.

- `four_body_half` (EQUIV, medium): a restated numeric check. Original: `2 * 1 + 4 * 1 + 8 * 1 - 32 * 1 + 16 = (-2 : Int)`, a bare numeral fact standing for 16 × the quartic at s = 1/2 (S = 2s = 1). Port: `2 * (1/2 : ℚ)^4 + 2 * (1/2)^3 + 2 * (1/2)^2 - 4 * (1/2) + 1 = -1/8`. The numbers agree (16 · (−1/8) = −2), and the port states what the docstring meant. The literal ℤ numeral identity is no longer stated.
- `theta2` changed from a match table to a double `Finset` sum with exponents `n * (n - 1) / 2` over ℕ. The exponents 0, 0, 1, 3 and the q = 0 values (0^0 = 1) match the original's r-table (1, 1, q, q³), so it is the same polynomial and `theta2_face`/`theta2_free_face` mean the same.
- The "iff" lemmas need `NoZeroDivisors`, which ℤ has. `kdv_weight_below_one` now holds in any strictly ordered commutative ring.

Route changes: `kdv_weight_nonneg` → mathlib `sq_nonneg` (his case-split proof is gone). `face_iff` now uses mathlib's `add_eq_zero_iff_of_nonneg`. `tau2_scaled` is proved from `mul_wAdd` (the original's note already called it `wAdd_scale`).

New port theorems with no original: 0.

Rung index note: no theorem in this module uses `MultAdd.rung`. Where the docstrings say "rung 1" or "rung s" they mean the fractional height s (s = 0 addition, s = 1 mult-add), the originals' own convention in this topic; that differs from Rungs.lean / `MultAdd.rung`, where mult-add is index 2. Same in original and port, so not a port finding.
