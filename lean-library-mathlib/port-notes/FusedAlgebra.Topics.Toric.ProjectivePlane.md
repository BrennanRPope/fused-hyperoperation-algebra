# FusedAlgebra.Topics.Toric.ProjectivePlane

- The point type `P2` is renamed `Plane` (`MultAdd.Toric.P2.Plane`), because
  `MultAdd.Toric.P2.P2` triggers mathlib's duplicated-namespace linter. The namespace is kept.
- The winding on chart coordinates `pw n v` is the monoid power `v ^ n` on `K × K`; the helper
  lemmas `inv_inv'`, `inv_ne`, `inv_pow'`, `pow_succ_eq_zero` are mathlib lemmas.
- `chart_wind` takes `n ≠ 0` instead of being stated for `n + 1`.
- Link to mathlib, new in the port: `toProjectivization : Plane K → ℙ K (Fin 3 → K)` and
  `toProjectivization_chart_c0/_c1/_c2`, which show the charts are `[a : b : 1]`, `[1 : u : v]`,
  `[s : 1 : t]`. NOT done: the proof that `toProjectivization` is a bijection (an `Equiv` like
  mathlib's `OnePoint.equivProjectivization` for the line). It is true and routine but would
  roughly double the file. No schemes are built.
- New: `card_plane`, the point count `q ^ 2 + q + 1` over a finite field.
