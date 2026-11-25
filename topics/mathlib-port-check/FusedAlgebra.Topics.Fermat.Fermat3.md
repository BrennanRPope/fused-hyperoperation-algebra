# FusedAlgebra.Topics.Fermat.Fermat3 — review

Theorems (108): MATHLIB 64, SAME 17, EQUIV 16, DROPPED-OK 6, GENERAL 4, MISSING 1.
Definitions checked: Eis, Eis.N, Eis.conj, Eis.IsUnit, Eis.Cop, Sol.

## Findings
- MISSING (minor), high: measure_lt.
  - Original: `(hx : x ≠ 0) (hy : y ≠ 0) (h : N1^3 < (z^3).natAbs) : N1 < (x * y * z).natAbs`.
  - The table maps it to fermatLastTheoremThree, which does not cover it. This is an unconditional arithmetic fact with satisfiable hypotheses, not a statement about Sol, so it is not vacuous once FLT3 is known.
  - It was only used by `step`, which is not ported, so nothing downstream breaks. The table's claim that it is covered is wrong.
- MATHLIB at medium: primitive_of_sol → fermatLastTheoremWith_of_fermatLastTheoremWith_coprime.
  - The mathlib lemma reduces to gcd({a, b, c}) = 1 and is stated as a contrapositive. The original reduces to gcd(x, y) = 1.
  - So the named lemma is not the same statement. The original's statement still follows (ex falso) from flt3.

## Conventions checked (as asked)
- Eisenstein integers: Eis = QuadraticAlgebra ℤ (-1) (-1), checked against mathlib source.
  - Multiplication: mathlib's mul is ⟨re·re' + a·im·im', re·im' + im·re' + b·im·im'⟩. With a = b = −1 this is the original's ⟨ac − bd, ad + bc − bd⟩, i.e. ω² = −1 − ω.
  - Norm: mathlib norm_def is re² + b·re·im − a·im². At a = b = −1 this is re² − re·im + im², the original N.
  - Conjugate: star has re = re − im and im = −im, the original conj.
  - Units list: the same six units in the same order.
- Int.gcd a b = 1 vs IsCoprime: every swap is through Int.isCoprime_iff_gcd_eq_one (exact over ℤ). `% 2 = 1` / `% 2 = 0` → Odd / Even are exact over ℤ. gcd_factors still concludes with Int.gcd = 1 ∨ 3; Sol uses IsCoprime.
- FLT: flt3 (ℤ, all three nonzero) and flt3_nat (ℕ, x, y ≠ 0, z free) have the same statements. The z = 0 gap between FLT and FermatLastTheoremFor is closed in flt3_nat's proof.
- no_sol changed shape: ∀ n x y z, |xyz| = n → ¬Sol becomes ∀ x y z, ¬Sol. These are equivalent (n was only the descent index).
- Mathlib names I was unsure of were checked against mathlib master: Int.isCoprime_two_left, Prime.coprime_iff_not_dvd, exists_associated_pow_of_mul_eq_pow', QuadraticAlgebra.norm_star / algebraMap_norm_eq_mul_star / re_* lemmas. All exist with the stated meaning. re_ofNat needs n ≥ 2; 0 and 1 are re_zero / re_one.

## Route changes (the big ones for the user)
- flt3, flt3_nat, flt3_multAdd and no_sol are now corollaries of mathlib's fermatLastTheoremThree. Euler's infinite descent is not ported:
  - `step` (190 lines) and its bookkeeping sol_cop, sol_swap, sol_rot, sol_even and measure_lt are gone.
  - The dropped Sol-lemmas are "covered" only ex falso, because Sol is uninhabited.
- key_cube (Euler's lemma) is still proved here, but its core changed:
  - The hand-written Bezout (bezout_aux, bezout, Cop, cop_bezout) and cube_of_cop are replaced by mathlib's EuclideanDomain / PID machinery (EuclideanDomain.gcd_eq_gcd_ab, isCoprime_of_dvd, exists_associated_pow_of_mul_eq_pow').
  - These run through a new instEuclideanDomain. That instance is built from the original's own rounding division (norm_mod_lt keeps the original div_rem argument).
- The integer gcd toolkit (cop_*, cube_of_mul_cop, three_dvd_*) is mathlib's IsCoprime / Prime API.
- The 45 or so Eis coordinate and ring lemmas are replaced by QuadraticAlgebra's ring structure.
- case1 is still proved by hand from cube_mod9; mathlib's version is private.

New port theorems with no original: 4 (Eis.norm_eq, Eis.div_def, Eis.mod_def, Eis.four_mul_norm_sub_mul). Also new instances: Div, Mod, IsDomain, EuclideanDomain on Eis.
