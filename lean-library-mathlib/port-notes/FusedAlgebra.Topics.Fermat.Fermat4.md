# FusedAlgebra.Topics.Fermat.Fermat4

- The hand-written infinite descent is not redone. `no_x4_y4_z2`, `fermat4` and
  `fourth_power_ne_mpow_four` are corollaries of mathlib's `not_fermat_42` and
  `fermatLastTheoremFour` (`Mathlib.NumberTheory.FLT.Four`).
- The three headline theorems stay in the root namespace under their original names, because
  `FermatRoute` refers to `fermat4` by that name. They are not mathlib-style names; rename them
  (for example into `namespace Fermat4`) if the project prefers.
- `Fermat4.descent_step` is dropped: with the theorem available its hypotheses are contradictory,
  so a port would be a vacuous statement.
- `Fermat4.pyth_param` is kept over ℕ in its original form, proved from
  `PythagoreanTriple.coprime_classification'`.
