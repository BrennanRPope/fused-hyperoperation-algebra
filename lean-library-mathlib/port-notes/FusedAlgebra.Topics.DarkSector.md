# FusedAlgebra.Topics.DarkSector

- `redshift_compose`, `redshift_compose_int`, `redshift_power` are exactly
  `MultAdd.one_add_multAdd` and `MultAdd.one_add_mpow` (with the sides swapped), so they are
  dropped rather than restated.
- `add_le_multAdd` and `sum_le_multAddList` (were `multAdd_ge_add`, `foldMA_ge_sum` on `ℕ`) are
  stated in an ordered semiring and so carry nonnegativity hypotheses; on `ℕ` these are
  `Nat.zero_le`. The helper `multAddList_nonneg` is new; it and `add_le_multAdd` would fit in
  `FusedAlgebra.Basic`.
- `lcdm_is_weighted_cube` is stated with `mpow z 3` on the right; `mpow_three` converts.
