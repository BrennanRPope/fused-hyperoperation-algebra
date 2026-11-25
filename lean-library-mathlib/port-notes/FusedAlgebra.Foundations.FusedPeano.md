# FusedAlgebra.Foundations.FusedPeano: porting remarks

- **Order of construction changed.** The original defined `add` and `mul` on a model by choice
  (`add x y` = the predecessor of `sadd x y`, `mul x y` = what is left of `madd x y`), proved
  their Peano recursions from the axioms, and only then proved categoricity. The port proves
  categoricity first (`toN_bijective`, `equivNat : N ≃ ℕ`, `toN_sadd`, `toN_madd`, all by
  induction from the axioms as before) and takes the semiring structure of `ℕ` across it
  (`Equiv.commSemiring`). The original characterisations are kept as theorems:
  `succ_add_eq_sadd : succ (x + y) = sadd x y` and `madd_eq : madd x y = x + y + x * y`; since
  `succ` is injective and `N` is cancellative these determine `+` and `*`, so it is the same
  structure. The primed Peano-recursion lemmas (`add_succ'`, `mul_succ'`, …) are now one-line
  consequences of mathlib lemmas plus `succ_eq_add_one` and are not restated.
- **New:** `equivNat : N ≃ ℕ`, `ringEquivNat : N ≃+* ℕ`, `succ_eq_add_one`, `sadd_eq`,
  `natCast_eq_toN`, `zero_eq`, `toN_one`, `toN_zero`, `toN_succ`, `ladN_zero`.
- **The ladder on `ℕ` is `MultAdd.rung`.** The original's `lad`/`hyp` and the lift `lad_lift`
  duplicate `FusedAlgebra.Rungs` (`rung`, mathlib's `hyperoperation`, `rung_add_one`), so they
  are not redefined. Three small general lemmas about `rung` that the schema needs live in this
  file and would belong in `Rungs.lean`: `MAFound.FusedPeano.rung_zero`,
  `rung_one_zero_right`, `rung_add_two_zero_right`.
- The module imports `FusedAlgebra.Rungs` (the original was stand-alone).
- `instCommSemiring` is a global instance `[FusedPeano N] → CommSemiring N`, as in the original.
