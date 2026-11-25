# FusedAlgebra.Topics.LinearAlgebra: port remarks

- The record types `M2`, `M3` and all their arithmetic are gone; everything is on mathlib's
  `Matrix`. Size-independent results (`shifted_det_hom`, `transpose_one_add_smul_mul`,
  `multAdd_inv_sub_one`, `msubM_multAdd`, `fl_*`) are stated for `Matrix n n R`, so each 2×2 / 3×3
  pair of originals maps to one theorem.
- `wmadd` is defined once, for any algebra over a commutative semiring. `MultAdd.wAdd` in
  `Basic.lean` needs a commutative semiring, so it cannot be used for matrices;
  `wmadd_eq_wAdd` links the two. A noncommutative `wAdd` in `Basic.lean` would make `wmadd`
  redundant.
- New, not in the original: `wmadd_one`, `wmadd_eq_wAdd`, `det_one_add_multAdd`,
  `det_one_add_smul_wmadd`, `det_one_add_of_isUpperTriangular` (any-size form of
  `triangular_eigen`), `det_one_add_smul_eq_dialTwo`, `det_one_add_smul_eq_dialThree`,
  `transpose_skew`.
- The dial invariant is only defined for sizes 2 and 3 (`dialTwo`, `dialThree`). For any size it
  would be read off mathlib's `Matrix.charpolyRev` (`det (1 - X • A)`); the general dial
  homomorphism needs cancelling `k` in a polynomial ring and was not attempted.
- `M3.Dk_hom` needed a field in the original (it divided by `k`); here `ring` proves it in any
  commutative ring.
- Inverse statements use `IsUnit (1 + A).det` in a commutative ring where the original had
  `det ≠ 0` in a field (`isUnit_iff_ne_zero` converts).
- Products `x * x` in statements are written `x ^ 2`.
- In `lake env lean` the mathlib style linters (line length) do not run; they do under
  `lake build`, which is clean.
