# FusedAlgebra.Topics.TwistedBoosts

Theorems (36): GENERAL 31, MATHLIB 4, DROPPED-OK 1. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

The change of carrier was checked:
- C.mul (ac − bd, ad + bc) is complex multiplication.
- C.conj is conj, C.smul is real •, C.nsq is normSq.
- `boost` = `!![1, u; -q * conj u, 1]` and `shape` = `!![α, β; -q * conj β, conj α]`.
- By hand: the entries of `boost_mul` (1 − q·u·v̄ = X + iY, 1 − q·ū·v = X − iY) and of `shape_mul` agree with the stated results.
- `oplus` and `twist` are now defined with complex division. `oplus_eq_smul` and `twist_eq_smul` prove, with no conditions, that they equal the original division-free formulas, because `Complex.inv_def` also gives 0 at 0.

Why the C-statements are rated medium rather than high: the original statements are about ℚ-pairs with q : ℚ, and the port's are about ℂ with q : ℝ. The original is the rational-coordinate special case. The ring identities, the iff and the two inequalities all transfer through the injective, order-preserving Rat cast (which also respects division, including by 0). But the port has no bridge lemma doing this. So "specialising" here means a cast argument, not plain instantiation. If the user wants the literal ℚ statements, one transfer lemma per statement (or a `C ≃ ℚ × ℚ ↪ ℂ` simp set) would settle it.

Other:
- `twist_half_angle` dropped the hypothesis X²+Y² ≠ 0. It follows from X ≠ 0 over ℚ as well, so the port is stronger.
- There are no numeric checks in this module.

Route changes:
- `sdet_mul`: the original grind expansion → `Matrix.det_mul` via `det_shape`.
- `sphere_no_wall`: the hand case split on signs → `Complex.normSq_nonneg`.
- `order_swap`: now by `field_simp` on complex division.

New port theorems with no original: 9.
