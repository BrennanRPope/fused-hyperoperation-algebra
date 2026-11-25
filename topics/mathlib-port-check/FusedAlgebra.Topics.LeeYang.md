# FusedAlgebra.Topics.LeeYang

Theorems (39): SAME 32, GENERAL 2, EQUIV 1, MATHLIB 3, DROPPED-OK 1. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

Definition change checked: `Mat2` → `Matrix (Fin 2) (Fin 2) R`. `ringT` is `!![1, z; 1, k * z]`, the same entries in the same layout. mul/tr/det are the same formulas. `Mat2.pow` multiplied on the left while `npow` multiplies on the right; the value is the same. The coefficients in the triangle/square/pentagon counts and in `even_ring_on_wall` are identical.

GENERAL: `triangle_disc_nonpos` (ℤ → linearly ordered comm ring) and `pottsRing_lift` (field → comm ring; the hypothesis has no division). EQUIV: `cayley_hamilton`, where the four explicit entries are now `tr • x - det • 1`.

Route changes:
- `odd_pow_of_tr_zero` / `even_pow_of_tr_zero`: the inductions are replaced by `pow_mul` + `sq_of_trace_eq_zero`.
- Helpers `one_pow'`, `two_pow_ge_four`, `neg_one_pow_cases`: replaced by mathlib.
- `triangle_disc_nonpos`: now proved by nlinarith.

New port theorems with no original: 2.
