# FusedAlgebra.Topics.SolitonsRungs — port remarks

- Everything was over `ℤ`; identities now hold in any commutative ring, the "exactly when"
  statements in any commutative ring without zero divisors, `kdv_weight_below_one` in any
  linearly ordered commutative ring. `tri_nonneg` and `face_iff` stay on `ℤ`.
- Repeated products are written as powers (`k ^ 3` for `k * k * k`).
- `theta2` is now a double `Finset` sum over `range 4` (same polynomial).
- `four_body_half` was a bare numeral identity (`2 + 4 + 8 - 32 + 16 = -2`); it is restated as
  the value `-1 / 8` of the quartic of `four_body_roads` at `1 / 2` in `ℚ`, which is what the
  docstring meant.
