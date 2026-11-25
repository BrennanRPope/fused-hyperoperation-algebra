# FusedAlgebra.Topics.DialSweep: port remarks

- The hand-written `charpoly2`, `charpoly3` are replaced by `Matrix.charpoly`; the characteristic
  polynomial statements (`det_one_add_neg_one_smul`, `det_one_add_smul_eq_charpoly`,
  `wall_iff_isRoot_charpoly`) are now for square matrices of any size, so each 2×2 / 3×3 pair of
  originals maps to one theorem. `eval_charpoly_fin_two` and `eval_charpoly_fin_three` (new) give
  back the explicit polynomials of the original.
- `sylvester2`, `sylvester3` are exactly `Matrix.det_one_add_mul_comm`; `outer3` is
  `Matrix.vecMulVec`; the rank-one results are the matrix determinant lemma.
- `wall_iff_isRoot_charpoly` says "root of the characteristic polynomial"; the step to
  `Module.End.HasEigenvalue` / `spectrum` was not added.
- `wAdd_neg_div` (was `winv`) is a general fact about `MultAdd.wAdd` and would fit in
  `Basic.lean` next to `wAdd_wall`.
