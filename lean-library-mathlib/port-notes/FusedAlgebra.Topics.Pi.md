# Remarks on the port of `FusedAlgebra.Topics.Pi`

- The hand-made Gaussian integers `G` are replaced by mathlib's `GaussianInt`; all the ring-law
  lemmas of the original (`gma_*`, `gpow_lift`, `norm_mul`, …) become instances of the general
  `MultAdd` lemmas or of `Zsqrtd` lemmas.
- `cross_10`, `cross_100`, `cross_1000` and `machin_50_digits` are kernel computations and are
  still proved by `decide +kernel` (as in the original; this is not `native_decide`, and adds no
  axiom).
- New, not in the original: `greg_succ` (the old defining recursion) and `hasSum_gregTerm`, a
  one-line link to mathlib's `Real.hasSum_arctan`.
- Not done: that the Machin intervals contain `Real.pi`
  (`machinLo N ≤ π ≤ machinHi N`). It would follow from `hasSum_gregTerm`, the monotonicity
  lemmas here and `Real.four_mul_arctan_inv_5_sub_arctan_inv_239`, but it is more than a
  one-liner.
