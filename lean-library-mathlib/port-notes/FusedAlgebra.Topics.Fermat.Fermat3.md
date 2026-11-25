# FusedAlgebra.Topics.Fermat.Fermat3

- **Headline theorems are corollaries of mathlib.** `flt3`, `flt3_nat`, `flt3_multAdd` and
  `no_sol` follow from `fermatLastTheoremThree` (`Mathlib.NumberTheory.FLT.Three`). Euler's
  descent (`step`, 190 lines, and its bookkeeping `sol_cop`, `sol_swap`, `sol_rot`, `sol_even`,
  `measure_lt`) is not re-proved: its hypotheses are contradictory once the theorem is known.
- **`case1` is proved here**, from `cube_mod9`. Mathlib has `fermatLastTheoremThree_case_1`, but
  that declaration is private to its (module-system) file and cannot be referenced.
- **Eisenstein integers.** The hand-rolled structure `Eis` is now
  `abbrev Eis := QuadraticAlgebra ℤ (-1) (-1)`, so ring structure, norm (`QuadraticAlgebra.norm`)
  and conjugation (`star`) come from mathlib and about 45 coordinate and ring lemmas disappear.
  Mathlib's own FLT3 uses `𝓞 (CyclotomicField 3 ℚ)`; `QuadraticAlgebra` was chosen because it has
  the same `re`/`im` coordinates as the original and the proofs carry over.
- **New relative to the original:** `Fermat3.Eis.instEuclideanDomain` (with `Div`/`Mod` instances
  for the rounding quotient, `norm_mod_lt`, `four_mul_norm_sub_mul`) and
  `Fermat3.Eis.instIsDomain`. With these, the original's hand-written Bezout (`bezout_aux`,
  `bezout`, `Cop`, `cop_bezout`) and `cube_of_cop` are replaced by mathlib
  (`EuclideanDomain.gcd_eq_gcd_ab`, `IsCoprime`, `exists_associated_pow_of_mul_eq_pow'`).
  These instances could move upstream or into a shared file if other modules need `ℤ[ω]`.
- **`key_cube` (Euler's lemma) is kept with a full proof**; it is not in mathlib.
- **Hypothesis forms changed (equivalent, not weaker):** `Int.gcd a b = 1` is `IsCoprime a b`
  (`Int.isCoprime_iff_gcd_eq_one`), `n % 2 = 1` is `Odd n`, `n % 2 = 0` is `Even n`. Divisibility
  by 3 and the `% 9` statements are unchanged.
- **`no_sol` changed shape:** from `∀ n x y z, |xyz| = n → ¬ Sol x y z` to `∀ x y z, ¬ Sol x y z`.
- **Renames:** `norm_one` → `eq_one_of_dvd_sq_add_three_mul_sq`, `Eis.four_N` →
  `Eis.four_mul_norm`, `Eis.N_nonneg` → `Eis.norm_nonneg`, `Eis.N_pos` → `Eis.norm_pos`,
  `Eis.eq_zero_of_N` → `Eis.norm_eq_zero_iff`, `Eis.omega_cube` → `Eis.omega_pow_three`,
  `Eis.div_rem` → `Eis.norm_mod_lt`, `Eis.isUnit_of_N` + `Eis.N_of_isUnit` →
  `Eis.isUnit_iff_norm_eq_one`.
- Inside `namespace Fermat3.Eis` with `open QuadraticAlgebra`, the bare name `norm` is ambiguous
  with the analytic `Norm.norm` outside that scope; write `QuadraticAlgebra.norm` there.
