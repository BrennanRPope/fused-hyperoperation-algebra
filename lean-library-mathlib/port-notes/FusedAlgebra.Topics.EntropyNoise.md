# FusedAlgebra.Topics.EntropyNoise: port remarks

- New import: `FusedAlgebra.Topics.Finance`. `cancellation_hides` is the two-point power mean
  inequality `FinanceMultAdd.arith_mean_pow_le` at `(1, x)`. This also showed that the original
  hypothesis `x ≤ 1` is not needed; it was removed (a strengthening).
- The "written without `abs`" statements (`mul_bound`, `octahedron_in_ball`,
  `compose_stays_in_octahedron`, `compose_stays_in_cube`) are restated with `|·|`. They are
  equivalent to the originals through `abs_le'` (`|a| ≤ u ↔ a ≤ u ∧ -a ≤ u`); in
  `octahedron_in_ball` the auxiliary bounds `u, v, w` are replaced by `|x|, |y|, |z|`.
  `mul_bound` is renamed `abs_mul_le_of_abs_le_one`; it is a general lemma (two lines from
  `abs_mul` and `mul_le_of_le_one_right`) that could live in a shared file.
- `x * x` is written `x ^ 2` throughout (`len2`, `coll` and every statement using them).
- The hand-rolled order and power facts (`sq_nn`, `pow_nn`, `pow_le_one'`, `mul_ne_zero'`,
  `div_pow'`, `pow_ne_zero'`, `lsum_div`) are mathlib lemmas.
- `Pauli` points are over any commutative ring; statements dividing by 2, 3 or 4 need a field of
  characteristic zero.
