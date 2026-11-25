module

public import FusedAlgebra.Topics.Fermat.Fermat3
public import FusedAlgebra.Topics.Fermat.Fermat4
public import FusedAlgebra.Topics.Fermat.FreyPieces

/-!
# Fermat's Last Theorem, minus exactly one deep step

This file joins the three Fermat files of the library:

* `Fermat4`: exponent `4` (Fermat's descent; mathlib's `fermatLastTheoremFour`);
* `Fermat3`: exponent `3` (Euler's descent in `ℤ[ω]`; mathlib's `fermatLastTheoremThree`);
* `FreyPieces`: reductions, normalisation and the Frey curve.

## What is proved

* `flt_of_prime_case`: if Fermat's Last Theorem holds for every prime `p ≥ 5`, it holds for every
  `n ≥ 3` (mathlib's `FermatLastTheorem`).
* `flt_of_no_frey_curve`: `FermatLastTheorem` follows from the single statement `NoFreyCurve`:
  no coprime solution `a ^ p + b ^ p = c ^ p` (`p ≥ 5` prime, `abc ≠ 0`) exists in the normalised
  shape `a ≡ 3 mod 4`, `b` even, `32 ∣ b ^ p`. That is the shape in which the Frey curve
  `y² = x(x - a ^ p)(x + b ^ p)` is semistable (`FreyPieces.frey_semistable`,
  `FreyPieces.frey_at_two`) with discriminant a `p`-th power up to `2 ^ 8`
  (`FreyPieces.frey_disc_pth_power`).

## What is not proved

`NoFreyCurve` itself. It is the content of Wiles and Taylor–Wiles (semistable elliptic curves are
modular), Mazur (the mod `p` representation is irreducible) and Ribet (level lowering to `N = 2`),
finished by the genus-`0` count for `X₀(2)`. Here it is an explicit hypothesis: never an axiom,
never `sorry`. So this file does **not** prove Fermat's Last Theorem; it proves the elementary
reduction of Fermat's Last Theorem to `NoFreyCurve`. Mathlib does not prove it either: it has
the exponents `3` and `4` and the reduction to odd primes (`FermatLastTheorem.of_odd_primes`).
-/

@[expose] public section

namespace FermatRoute

open FreyPieces

/-- Fermat's Last Theorem for every `n ≥ 3` follows from the prime case `p ≥ 5` alone, because
the exponents `3` and `4` are proved. -/
theorem flt_of_prime_case (hp : ∀ p : ℕ, p.Prime → 5 ≤ p → FermatLastTheoremFor p) :
    FermatLastTheorem :=
  flt_reduction_frey (fermatLastTheoremFor_iff.2 fermat4)
    (fermatLastTheoremFor_iff.2 Fermat3.flt3_nat) hp

/-- Frey normalisation keeping `gcd a' b' = 1`: the input the Frey curve needs. -/
theorem frey_normalize_cop {p : ℕ} (hp : Odd p) (h5 : 5 ≤ p) {a b c : ℤ}
    (hcop : Int.gcd a b = 1) (habc : a * b * c ≠ 0) (h : a ^ p + b ^ p = c ^ p) :
    ∃ a' b' c' : ℤ, a' ^ p + b' ^ p = c' ^ p ∧ a' % 4 = 3 ∧ b' % 2 = 0 ∧ Int.gcd a' b' = 1 ∧
      a' * b' * c' ≠ 0 ∧ (32 : ℤ) ∣ b' ^ p := by
  obtain ⟨a', b', c', h1, h2, h3, -, h5, h6, h7⟩ := frey_normalize hp h5 hcop habc h
  exact ⟨a', b', c', h1, h2, h3, h5, h6, h7⟩

/-- The one deep statement: there is no normalised coprime Fermat solution for a prime `p ≥ 5`
(equivalently, the semistable Frey curve with these numbers cannot exist). This is a definition,
not a theorem: it is what Wiles, Taylor–Wiles, Mazur and Ribet prove, and it is not proved in
this library. -/
def NoFreyCurve : Prop :=
  ∀ p : ℕ, p.Prime → 5 ≤ p → ∀ a b c : ℤ, Int.gcd a b = 1 → a * b * c ≠ 0 → a % 4 = 3 →
    b % 2 = 0 → (32 : ℤ) ∣ b ^ p → a ^ p + b ^ p ≠ c ^ p

/-- Fermat's Last Theorem for every `n ≥ 3`, given only `NoFreyCurve`. -/
theorem flt_of_no_frey_curve (H : NoFreyCurve) : FermatLastTheorem := by
  refine flt_of_prime_case fun p hp h5 ↦ flt_of_coprime hp.ne_zero fun x y z hx hy hcop hxyz ↦ ?_
  have hz : z ≠ 0 := by
    rintro rfl
    rw [zero_pow hp.ne_zero] at hxyz
    exact pow_ne_zero p hx (Nat.eq_zero_of_add_eq_zero_right hxyz)
  have hI : (x : ℤ) ^ p + (y : ℤ) ^ p = (z : ℤ) ^ p := mod_cast hxyz
  have hg : Int.gcd (x : ℤ) (y : ℤ) = 1 := by rwa [Int.gcd_natCast_natCast]
  have habc : (x : ℤ) * y * z ≠ 0 := mod_cast mul_ne_zero (mul_ne_zero hx hy) hz
  obtain ⟨a', b', c', he, ha4, hb2, hg', hne, h32⟩ :=
    frey_normalize_cop (hp.odd_of_ne_two (by lia)) h5 hg habc hI
  exact H p hp h5 a' b' c' hg' hne ha4 hb2 h32 he

/-- `NoFreyCurve` is not only sufficient but necessary: it is a special case of Fermat's Last
Theorem over `ℤ`. So `NoFreyCurve ↔ FermatLastTheorem`, and nothing is lost in the reduction. -/
theorem noFreyCurve_of_flt (H : FermatLastTheorem) : NoFreyCurve := by
  intro p _ h5 a b c _ habc _ _ _
  have h := fermatLastTheoremFor_iff_int.1 (H p (by lia)) a b c
  simp only [mul_ne_zero_iff] at habc
  exact h habc.1.1 habc.1.2 habc.2

end FermatRoute
