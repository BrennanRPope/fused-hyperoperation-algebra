# FusedAlgebra.Topics.LeeYangCluster

Theorems (10): SAME 7, GENERAL 2, DROPPED-OK 1. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

The polynomials are unchanged (`k*k` is now written `k^2`, `k*k*k` is `k^3`).
- `wAdd_units` and `leverage_is_compounding` are generalised to comm semirings. Leverage's right side `wAdd 1 (2a) (2b)` is now `(2a) ✳ (2b)`, which is equal by `wAdd_one`.
- `small_of_mul` → `Int.natAbs_le_of_dvd_ne_zero`, which has the same content.

Route change: the Beraha proofs now use the mathlib divisor bound + `interval_cases` instead of the original's own `small_of_mul` + omega case split.

New port theorems with no original: 0.
