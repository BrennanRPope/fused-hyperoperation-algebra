# FusedAlgebra.Topics.NumberPlane: port remarks

- The number system "a + b·u with u² = k·u − q" is exactly mathlib's `QuadraticAlgebra R (-q) k`
  (ω·ω = −q + k·ω), and the norm is `QuadraticAlgebra.norm`. The coordinate functions `mre`, `mu`
  and `norm` are kept as definitions because GeometryTrig (and the prose pages) quote them with
  four explicit numbers; the new lemmas `mre_eq_re_mul`, `mu_eq_im_mul`, `mk_mul_mk` and
  `norm_eq_norm` are the bridge. `mul_assoc_re/u` and `norm_mul` are now read off `mul_assoc` and
  `map_mul` in `QuadraticAlgebra` instead of being expanded by hand.
- Everything is generalised from `ℤ` to an arbitrary commutative ring (the order statement to a
  linearly ordered commutative ring). Squares `x * x` in statements are written `x ^ 2`.
- New helper lemmas not in the original: `mre_def`, `mu_def`, `norm_def`, `mk_mul_mk`,
  `mre_eq_re_mul`, `mu_eq_im_mul`, `norm_eq_norm`.
