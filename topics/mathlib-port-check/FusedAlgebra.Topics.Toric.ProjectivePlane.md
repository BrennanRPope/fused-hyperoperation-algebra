# FusedAlgebra.Topics.Toric.ProjectivePlane

Theorems (16): SAME 9, EQUIV 2, GENERAL 1, MATHLIB 4.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.
- `chart`, `trans`, `Dom`, `wind` identical case by case; the point type is only renamed (`P2` -> `Plane`).
- `chart_wind`: exponent `n + 1` became a hypothesis `n ≠ 0` (equivalent; needed, since n = 0 is false).
- The four helper lemmas are standard mathlib lemmas (inv_inv, inv_ne_zero, inv_pow, pow_eq_zero_iff).

Route changes:
- `trans01_multAdd`, `trans02_multAdd`: grind on unfolded definitions replaced by the library lemmas `one_add_multAdd` and `inv_one_add` (his own, in Basic).
- Helper lemmas he proved by hand (`inv_pow'`, `pow_succ_eq_zero`, ...) are now mathlib citations.

New port theorems with no original: 4 (card_plane, toProjectivization_chart_c0/c1/c2). As the porter says, `toProjectivization` is not proved to be a bijection, so the link to mathlib's ℙ² is partial.
