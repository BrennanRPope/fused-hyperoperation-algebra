# FusedAlgebra.Topics.ClosingRate

- New in the port (not in the original): `log_two_lt`, `lt_log_two`, `two_mul_log_two_sq_lt_one`,
  `inv_log_four_bounds`. The original could only state integer stand-ins for these (`two_lt`,
  `e_pow_lt`, `sq_bound`, `inv_lo`, `inv_hi`) and left the step to `log 2` to the reader; with
  mathlib's `Real.log_two_lt_d9` / `Real.log_two_gt_d9` the real statements are two lines. The
  integer facts are kept so the names still exist, but they are no longer needed for anything.
- `grow` is generalised to an ordered commutative semiring, which costs the explicit hypotheses
  `0 ≤ p` and `0 ≤ q`. On `ℕ` pass `(Nat.zero_le _)`.
- `two_lt` and `e_pow_lt` need `set_option exponentiation.threshold 1000 in`, as in the original.
