# FusedAlgebra.Topics.DialSweep

Theorems (28): SAME 14, GENERAL 11, EQUIV 1, MATHLIB 2. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

Definitions checked:
- `wmadd`/`madd` on records → `x + y + k • (x * y)` and `✳` in any algebra.
- `dialU = !![0,-q;1,k]` and `unitU = !![0,κ;1,0]` have the same entries.
- Port `e2` equals the original `M3.e2`.
- `outer3` → `vecMulVec`.
- `charpoly2`/`charpoly3` → `Matrix.charpoly` evaluated at t, with `eval_charpoly_fin_two`/`_three` giving back the explicit polynomials.
- `M2.inv` → mathlib's matrix ⁻¹. Both are adj/det, and both are 0 when singular. `cayley` is defined through the inverse in each library, so `cayley_weighted` says the same thing.

GENERAL means the 2×2/3×3 pairs are now one theorem for any n. Specialising the factor `(-k)^card n` gives `k*k` and `-(k*k*k)`, as in the originals. `rank_one_dial3`'s explicit triple sum is `v ⬝ᵥ (A *ᵥ u)`.

Route changes:
- `sylvester2`/`sylvester3`: the original explicit expansions (sylvester3 needed 2M heartbeats) are replaced by mathlib `Matrix.det_one_add_mul_comm`.
- `rank_one3`/`rank_one_dial3`: now derived from mathlib's matrix determinant lemma.
- `det_charpoly*`: now via `eval_charpoly` + `det_smul`.

New port theorems with no original: 3.
