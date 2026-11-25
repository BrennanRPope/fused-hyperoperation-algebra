# FusedAlgebra.Topics.Toric.ToricCharts

- `P1 K` is now mathlib's `OnePoint K` (definitionally `Option K`, so nothing is lost), which
  gives the exact link to mathlib for free: `P1.equivProjectivization : P1 K ≃ ℙ K (Fin 2 → K)`
  is `OnePoint.equivProjectivization`. New lemmas `equivProjectivization_chart1_false/true` show
  the two charts are `x ↦ [x : 1]` and `x ↦ [1 : x]`.
- No schemes are built. `AlgebraicGeometry.ProjectiveSpace` is only mentioned in the docstring.
- New: `trans1_one_add` (the gluing map in mult-add coordinates is `multAddInv`; the original
  said this only in its header comment), `corner_eq_true_iff` (the corner is `Function.support`),
  `card_P1_prod_P1` (the point count over a finite field, which `fan_count` only hinted at),
  and the small API lemmas `chart1_false`, `chart1_true_zero`, `chart1_true_of_ne_zero`.
- Injectivity statements are now `Function.Injective`.
- Statements keep `[Field K]` as in the original, although `chart1`, `trans1` and the cocycle
  only need a `GroupWithZero`.
