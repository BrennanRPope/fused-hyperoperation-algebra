# FusedAlgebra.Topics.Topology.HopfInvariantOne: review

Theorems (17): SAME 8, EQUIV 1, MATHLIB 7, DROPPED-OK 1. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked:
- `OneMod j x` became `multiplicity 2 (x - 1) = j`. The two agree for every j ≥ 1. They differ only at j = 0, x = 1, because mathlib's ℕ-valued `multiplicity 2 0 = 0`. That case never occurs: `level_three_pow` uses levels 1 and k+2.
- `Trunc` with the new `NonUnitalCommRing` (mul `⟨0, p.α*q.α⟩`): I checked the ring laws, including associativity (x³ = 0).
- Library `✳` on `Trunc` is the old `madd`; `Trunc.multAdd_eq` records this. `mpow Trunc.x` has the same recursion as the old `mpow` (`mpow_succ`: `mpow x k ✳ x`).
- `Elt`, `psi`, `psiT` and `IsOdd` = `Odd`: identical.

Medium-confidence items (no finding; content checked, only the lemma mapping is approximate):
- `oneMod_oddPow`: maps to `emultiplicity_pow_sub_pow_of_prime` at p = 2, y = 1 over ℤ. Its hypotheses (2 ∣ A−1, A odd, n odd) hold, but the statement has to be translated ℕ→ℤ and OneMod→emultiplicity.
- `oneMod_sq`: maps to `Int.two_pow_two_pow_sub_pow_two_pow` with n = 1. Same translation needed.
- `le_of_dvd_two_pow_odd`: the named `FiniteMultiplicity.pow_dvd_iff_le_multiplicity` gives only half. The other half, `multiplicity 2 (2^e * odd) = e`, is proved inline as `hmul` in `level_three_pow`.
- `two_adic_split` (`Nat.exists_eq_two_pow_mul_odd`) and `two_pow_dvd_of_mul_odd` (`Int.isCoprime_two_left` + `IsCoprime.dvd_of_dvd_mul_right`): mathlib names believed correct.

Proof-route changes:
- Part 1 no longer uses your `OneMod` lemmas (`oneMod_oddPow`, `oneMod_sq`, `three_pow_two_pow`). Neither general lemma is restated, even for base 3 only. `level_three_pow` and `only_one_two_four` now come from mathlib's lifting-the-exponent lemma `Nat.two_pow_sub_pow` plus `multiplicity`.
- `two_pow_dvd_of_mul_odd`: the hand induction is replaced by mathlib coprimality.
- Parts 2 and 3 (`commute_equation`, `hopf_invariant_one`, `mpow_formula`, `psi_compose`, `cp2_case`) keep your proofs (`linear_combination` in place of `grind`).

New port theorems with no original: 13 (`multiplicity_two_three_pow_sub_one_of_odd`/`_of_even`, ten `Trunc` simp lemmas, `Trunc.multAdd_eq`).
