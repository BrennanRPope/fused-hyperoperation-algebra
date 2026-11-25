# FusedAlgebra.Topics.DarkSector

Theorems: 17. GENERAL 10, SAME 4, DROPPED-OK 3. No findings.

- `add_le_multAdd` (was `multAdd_ge_add`) and `sum_le_multAddList` (was `foldMA_ge_sum`): PORT-NOTES claim checked and correct.
  - Both originals are over ℕ: `(a b : Nat) : a + b ≤ multAdd a b` and `∀ l : List Nat, l.sum ≤ foldMA l`.
  - The port's added hypotheses `0 ≤ a`, `0 ≤ b` and `∀ x ∈ l, 0 ≤ x` are automatic on ℕ (`Nat.zero_le`), so the ℕ statements follow with nothing extra. GENERAL, not WEAKER. The hypotheses are needed in general, e.g. over ℤ with a = −1, b = 2: a + b = 1 but a ✳ b = −1.
- `redshift_compose`, `redshift_power` and `redshift_compose_int` are dropped onto `one_add_multAdd` and `one_add_mpow` (sides swapped; the ℤ one with ✳ unfolded). Correct.
- `lcdm_is_weighted_cube`: the right side is now `1 + Om * mpow z 3` instead of `1 + Om * (3z + 3z² + z³)`. These are equal by the port's `mpow_three`.
- The growth no-go lemmas use the same `fexp` (rung 3, `(a+1)^(b+1) - 1`) on ℕ, as in the original.

Route changes: `corner_le_quarter` is now mathlib `four_mul_le_sq_add` instead of his (a−b)² ≥ 0 argument. `cube_is_mpow3` (now `mpow_three`) is proved by `ring` rather than via `redshift_power`.

New port theorems with no original: 1 (`multAddList_nonneg`).
