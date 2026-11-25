module

public import FusedAlgebra.Basic

/-!
# Fermat's Last Theorem for exponent 4, and its mult-add form

The core-Lean library proved by hand, with Fermat's infinite descent, that `x ^ 4 + y ^ 4` is
never a square for nonzero naturals. That is known mathematics and is already in mathlib:
`not_fermat_42` (no nonzero integer solutions of `a ^ 4 + b ^ 4 = c ^ 2`) and
`fermatLastTheoremFour : FermatLastTheoremFor 4`, both in `Mathlib.NumberTheory.FLT.Four`. So
this file does not redo the descent. It states the library's headline theorems as corollaries of
the mathlib ones, and keeps the two pieces of the descent that are statements in their own right.

## Main results

* `Fermat4.sq_mod_four`: a square is `0` or `1` mod `4`, according to parity.
* `Fermat4.pyth_param`: the parametrisation of primitive Pythagorean triples over `ℕ`, a
  corollary of mathlib's `PythagoreanTriple.coprime_classification'`.
* `no_x4_y4_z2`: `x ^ 4 + y ^ 4 ≠ z ^ 2` for nonzero naturals (mathlib's `not_fermat_42`).
* `fermat4`: `x ^ 4 + y ^ 4 ≠ z ^ 4` for nonzero naturals (mathlib's `fermatLastTheoremFour`).
* `fourth_power_ne_mpow_four`: the mult-add form. Since `1 + mpow s 4 = (1 + s) ^ 4`, the
  statement "`q ^ 4` is the mult-add fourth power of `s`" means `1 + q ^ 4 = (1 + s) ^ 4`; no
  nonzero rational `q` and rational `s` do this.

The three headline theorems sit in the root namespace, as they did in the original.
-/

@[expose] public section

namespace Fermat4

/-- Every square is `0` or `1` mod `4`: `0` for an even number, `1` for an odd one. -/
theorem sq_mod_four (n : ℕ) : (n % 2 = 0 ∧ n ^ 2 % 4 = 0) ∨ (n % 2 = 1 ∧ n ^ 2 % 4 = 1) := by
  obtain ⟨k, rfl | rfl⟩ := n.even_or_odd'
  · have : (2 * k) ^ 2 = 4 * (k * k) := by ring
    lia
  · have : (2 * k + 1) ^ 2 = 4 * (k * k + k) + 1 := by ring
    lia

/-- Primitive Pythagorean triples over `ℕ`: `a ^ 2 + b ^ 2 = c ^ 2` with `a` odd and `a`, `b`
coprime gives coprime `m`, `n` with `a = m ^ 2 - n ^ 2` (written `a + n ^ 2 = m ^ 2`),
`b = 2 * m * n` and `c = m ^ 2 + n ^ 2`. This is mathlib's
`PythagoreanTriple.coprime_classification'`, brought down from `ℤ` to `ℕ`. -/
theorem pyth_param {a b c : ℕ} (h : a ^ 2 + b ^ 2 = c ^ 2) (ha : a % 2 = 1)
    (hab : Nat.Coprime a b) :
    ∃ m n, a + n ^ 2 = m ^ 2 ∧ c = m ^ 2 + n ^ 2 ∧ b = 2 * m * n ∧ Nat.Coprime m n := by
  have hc : 0 < c := by
    rcases Nat.eq_zero_or_pos c with rfl | hc
    · have : a = 0 := by nlinarith
      lia
    · exact hc
  have hpt : PythagoreanTriple (a : ℤ) b c := by
    unfold PythagoreanTriple; exact_mod_cast (by rw [← sq, ← sq, ← sq]; exact h)
  obtain ⟨m, n, h1, h2, h3, h4, -, -⟩ := hpt.coprime_classification'
    (by rw [Int.gcd_natCast_natCast]; exact hab) (by lia) (by exact_mod_cast hc)
  refine ⟨m.natAbs, n.natAbs, ?_, ?_, ?_, h4⟩
  · zify; rw [sq_abs, sq_abs]; lia
  · zify; rw [sq_abs, sq_abs]; exact h3
  · have := congrArg Int.natAbs h2
    simpa [Int.natAbs_mul] using this

end Fermat4

/-- No sum of two nonzero fourth powers is a square. This is mathlib's `not_fermat_42`, stated
over `ℕ`. -/
theorem no_x4_y4_z2 (x y z : ℕ) (hx : x ≠ 0) (hy : y ≠ 0) : x ^ 4 + y ^ 4 ≠ z ^ 2 := fun h ↦
  not_fermat_42 (a := (x : ℤ)) (b := y) (c := z) (mod_cast hx) (mod_cast hy) (mod_cast h)

/-- Fermat's Last Theorem for exponent 4, a corollary of `no_x4_y4_z2` (`z ^ 4` is the square of
`z ^ 2`). In mathlib this is `fermatLastTheoremFour`. -/
theorem fermat4 (x y z : ℕ) (hx : x ≠ 0) (hy : y ≠ 0) : x ^ 4 + y ^ 4 ≠ z ^ 4 := fun h ↦
  no_x4_y4_z2 x y (z ^ 2) hx hy (by rw [h, ← pow_mul])

open MultAdd in
/-- The mult-add form of Fermat's Last Theorem for exponent 4: for rationals `q ≠ 0` and `s`,
`q ^ 4` is never the mult-add fourth power `mpow s 4 = (1 + s) ^ 4 - 1`. A corollary of
mathlib's `fermatLastTheoremFour`, through `fermatLastTheoremFor_iff_rat`. -/
theorem fourth_power_ne_mpow_four (q s : ℚ) (hq : q ≠ 0) : q ^ 4 ≠ mpow s 4 := by
  intro h
  have hlift : 1 ^ 4 + q ^ 4 = (1 + s) ^ 4 := by rw [h, one_pow, one_add_mpow]
  have hs : 1 + s ≠ 0 := by
    intro h0
    rw [h0] at hlift
    nlinarith [sq_nonneg (q ^ 2)]
  exact fermatLastTheoremFor_iff_rat.1 fermatLastTheoremFour 1 q (1 + s) one_ne_zero hq hs hlift
