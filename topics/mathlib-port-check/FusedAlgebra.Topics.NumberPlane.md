# FusedAlgebra.Topics.NumberPlane: review

Theorems (26): GENERAL 25, MATHLIB 1. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked:
- `mre`, `mu` and `norm` have the same formulas, over any commutative ring.
- `LandscapeMap.num`/`den` (used by `tilted_*`) are unchanged in the port's Landscape: `num` = `wAdd`, `den = 1 - q*a*b`.
- Hypotheses `u * u = k*u - q` and `k * k ≤ 4*q` are now written with `^ 2`, which is equivalent.
- `norm_nonneg_turning` is generalised to a linearly ordered commutative ring, and ℤ qualifies.
- I also checked the bridge to `QuadraticAlgebra R (-q) k` (ω² = −q + kω): `re_mul` gives ac − qbd, and the norm is a² + kab + qb².

Proof-route changes: `mul_assoc_re/u` and `norm_mul` are now read off `mul_assoc` and `map_mul` in `QuadraticAlgebra`, where the original expanded them by hand. The statements are identical.

New port theorems with no original: 7 (`mre_def`, `mu_def`, `norm_def`, `mk_mul_mk`, `mre_eq_re_mul`, `mu_eq_im_mul`, `norm_eq_norm`).
