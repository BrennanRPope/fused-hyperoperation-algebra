# FusedAlgebra.Topics.LinearAlgebra: review

Theorems (51): GENERAL 23, SAME 17, MATHLIB 9, EQUIV 1, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence items.

Record to Matrix check: M2 ⟨a,b,c,d⟩ = [[a,b],[c,d]] maps to M 0 0, M 0 1, M 1 0, M 1 1; M3 aij maps to M (i-1)(j-1). I checked the record mul (row by column, same order), transpose, trace, and det against det_fin_two (ad - bc) and det_fin_three (all six terms, same signs). e2 is the same with indices shifted. I checked every literal (madd_not_comm, rotation_condition, diagonal_is_multAdd, box_cube, kdv_two_soliton, skew, cayley_eq) and every entry hypothesis (coupling_sets_weight; triangular_eigen with a21, a31, a32 becoming A 1 0, A 2 0, A 2 1). Record inv (entries d/det etc.) agrees with Matrix.inv (Ring.inverse det • adjugate) in a field. Inverse theorems moved from "field, det ≠ 0" to "comm ring, IsUnit det", which is equivalent in a field.

- `rank_one_is_add` (EQUIV): the hypothesis b*c = a*d is restated as det = 0 and the conclusion a + d as trace. Equivalent.
- `msub_undo_W` (DROPPED-OK): it is `msubM_multAdd` at A = w - 1.

Table remarks: `M2.lift`, `M3.lift`, `madd_assoc`, `madd_zero`, `zero_madd` are marked "mathlib", but the replacements are the port's own Basic lemmas (one_add_multAdd etc.). They do cover the matrix case, since those lemmas hold in any (non-unital) semiring (cosmetic).

Route changes worth knowing:
- `M2.det_mul` and `M3.det_mul` (hand expansions, the 3×3 one with raised heartbeats) are now mathlib's `Matrix.det_mul`. `shifted_det_hom` now goes through mathlib det_mul.
- `madd_inverse` and `msub_undo` now go through mathlib `mul_nonsing_inv` / `nonsing_inv_mul`; `M2.mul_inv` is mathlib.
- `M3.Dk_hom`: the original split on k = 0 and divided by k (field); the port proves it by `ring` in any comm ring (GENERAL).
- `pow_mul''` (hand induction) is now mathlib `pow_mul`.

New port theorems with no original: 9.
