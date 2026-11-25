# FusedAlgebra.Topics.QAcross.Partitions — port remarks

- All the series theorems remain finite checks on truncated coefficient lists, as in the
  original (`decide +kernel`). They are not connected to mathlib's `PowerSeries`.
- `partList` (p(0..50)) is not connected to mathlib's `Fintype.card (Nat.Partition n)`: that
  cardinality does not reduce by `decide`, and a proof through `Nat.Partition.hasProd_genFun`
  would be real work. This is the main gap between this file and "mathlib form".
- Replaced by mathlib: `isSq` → `IsSquare`; `sigmaStar` is now a `Finset` sum over
  `Nat.divisors`.
- New general statements next to the finite checks: `one_add_multAddList_neg_pow` (Euler's
  product is `1 +` a mult-add, any ring), `r4_prime` (`r4 p = 8 * (p + 1)` for every prime),
  `r4_pos` (every `n`). Note `r4` is Jacobi's formula, not a count of representations, so
  `r4_pos` is not Lagrange's theorem by itself; mathlib has that as `Nat.sum_four_squares`.
