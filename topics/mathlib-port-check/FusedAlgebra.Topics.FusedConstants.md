# FusedAlgebra.Topics.FusedConstants

Theorems (20): SAME 12, EQUIV 1, GENERAL 4, MATHLIB 2, DROPPED-OK 1. Definitions checked: H (EQUIV; mathlib's hyperoperation (r+1) is the same function on every case), slide (GENERAL).

There is no theorem verdict of WEAKER, DIFFERENT, MISSING or VACUOUS, and none at low confidence.

Route changes:
- `H` is now mathlib's `hyperoperation (r + 1)`, so `H_succ`, `H_zero`, `H_add`, `H_mul` and `H_pow` come from mathlib lemmas instead of his own definition.
- `all_rungs_two_two` was his own induction through `H_two_one`; it is now mathlib's `hyperoperation_two_two_eq_four`. `H_two_one` itself exists only as mathlib's `hyperoperation_ge_two_eq_self`.
- `lt_two_pow'` is replaced by `Nat.lt_two_pow_self`.

New port theorems with no original: 1 (`rung_succ_eq`).
