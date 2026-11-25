# Remarks on the port of `FusedAlgebra.Topics.OsgoodInverse`

- The original imported `Topics.Pi` (for `div_le_div_of`) and `Topics.NSContinuum` (for
  `mpow_one_nat`). Both are mathlib / `FusedAlgebra.Cube` facts now, so the port imports only
  `FusedAlgebra.Cube`.
- Reciprocals `1 / (1 + y)` in the Riccati clock statements are written `(1 + y)⁻¹`, and
  `steps h n` is `n • h`.
- Cauchy condensation is in mathlib (`Finset.le_sum_condensed'`, `Finset.sum_condensed_le`,
  `summable_condensed_iff_of_nonneg`). The port keeps the block-by-block derivation
  (`blk_le`, `le_blk`, `dsum_succ`) because the blocks are the point of the section; the lower
  squeeze `csum_le_dsum` is normalised differently from mathlib's (`Ico 1 (2 ^ n)` instead of
  `Ico 2 (2 ^ n + 1)`), so it is not literally a mathlib corollary.
- `const_block` (any start `s`) is now `blk_const` (start `2 ^ j` only), its only use.
- The instances `harm`, `isq` stay over `ℚ`; the general theory (`blk`, `dsum`, `csum`,
  the squeezes, `condensation`) is stated in any linearly ordered field.
- `lipschitz_unique` holds in any field, including characteristic 2 (handled by a case split).
