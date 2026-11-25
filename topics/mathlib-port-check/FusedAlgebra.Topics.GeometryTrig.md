# FusedAlgebra.Topics.GeometryTrig: review

Theorems (24): GENERAL 21, SAME 3. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked:
- `fall` is identical.
- `fallen` is now `∑ i ∈ range n, fall (i + 1)`, with the same recursion.
- Statements over ℤ are generalised to commutative (semi)rings. `x * x` is written `x ^ 2`, `x*x*x` is `x ^ 3`, and `a ✳ a ✳ a` parses as `(a ✳ a) ✳ a`, as in the original.
- `triangle_lengths_real` (k² ≤ 4) is over a linearly ordered commutative ring, and ℤ qualifies.
- The two examples have the same numbers: 3-4-5 and 5-12-13 triangles, and norm 3 1 1 (−1) = −1.

Proof-route changes: none of substance. These are polynomial identities, closed by `ring` / `linear_combination`. `unit_circle_closed` is no longer derived from `norm_mul`, and `euclid_II_4` is mathlib's `add_sq`.

New port theorems with no original: 5 (`fall_zero`, `fall_succ`, and three bridges to ℂ: `mre_one_eq_re_mul`, `mu_zero_eq_im_mul`, `norm_zero_one_eq_normSq`).
