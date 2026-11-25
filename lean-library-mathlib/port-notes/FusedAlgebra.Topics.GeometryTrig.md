# FusedAlgebra.Topics.GeometryTrig: port remarks

- The two anonymous `example`s of the original (3-4-5 and 5-12-13; `norm 3 1 1 (-1) = -1`) are
  kept as `example`s; they have no names, so they are not in the table.
- New, not in the original: `fall_zero`, `fall_succ`, and three bridges saying that over `ℝ` the
  complex corner `(k, q) = (0, 1)` is mathlib's `ℂ`: `mre_one_eq_re_mul`, `mu_zero_eq_im_mul`,
  `norm_zero_one_eq_normSq`.
- The statements stay polynomial identities (as in the original). The analytic addition laws
  (`Real.cos_add`, `Real.sin_add`, `Real.tan_add`) are cited in the docstrings but not restated.
- Generalised from `ℤ` to commutative (semi)rings; `x * x` is written `x ^ 2`.
