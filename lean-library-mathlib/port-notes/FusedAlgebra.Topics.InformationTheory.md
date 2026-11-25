# FusedAlgebra.Topics.InformationTheory: port remarks

- `lsum` is `List.sum`. `outer` is kept (mathlib has no "list of pairwise products"), defined as
  `(xs ×ˢ ys).map fun p ↦ p.1 * p.2`; `outer_cons` recovers the original recursion. `shan` is kept
  as a `List.sum` of a `map`.
- Several originals state the same identity under different topical names (`kraft_slack` =
  `gini_simpson`, `campbell` = `variable_node` = `support_count`, `bec_capacity` = `survival`,
  `parallel_capacity` = `hartley`, `power_sums_multiply` = `outer_sum`, `lengths_add` =
  `shannon_weight`). All names are kept, each proved by the first of its group.
- `union_rule`, `union_is_flipped_madd`, `series_system`, `shannon_weight`, `check_node` are
  one-line corollaries of `MultAdd` lemmas in `Basic`; `climb` was dropped because it is exactly
  `MultAdd.one_add_mul_wpow`.
- Generality: identities in a commutative (semi)ring; `piling_up` and `ladder_four` in a field of
  characteristic zero; inequalities in a linearly ordered commutative ring or field.
- `tsallis`: the variable `R` of the original is renamed `Q` (clash with the ring `R`).
- `parity3` and `sdpi_multiplies` write `x * x` as `x ^ 2`.
- Nothing is connected to mathlib's measure-theoretic entropy; the statements are algebraic.
