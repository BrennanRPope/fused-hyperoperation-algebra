/-
DRAFT (Nov 27, 2026). The two-base letter product: Hoffman–Ihara's q letter product with q a jet 1 + N,
N² = kN − d. Companion to `two_base_checks.py` (13/13) and `QDeformationDial.lean`.

Main points:
* `pow_eq_ring`: qⁿ = 1 + (q − 1)[n]_q in any commutative ring (so q may be a jet or a matrix).
* `letter_mul`: Bradley's letters q^(a n)/[n]^(a+1) multiply by zᵢ◇zⱼ = z_(i+j) + (1 − q) z_(i+j−1)
  whenever [n]_q is a unit, in any commutative ring.
* `jet_letter_mul`: with q = 1 + N the correction is the jet −N.
* `jet_letter_triple`: the triple product carries d as a plain number: (kN − d) z_(S−2).
* `letter_map`: a ring map (an eigen-shadow N ↦ α) sends jet letters to scalar letters.
* `root_sq`: in AdjoinRoot (X² − kX + d) the root satisfies N² = kN − d.
* `dial_from_bases`, `hold_base_tangent`: k = q₁ + q₂ − 2, d = (q₁ − 1)(q₂ − 1); holding one base at c
  gives a line tangent to the border d = k²/4.
Not yet built; meant for CI's report-only Mathlib step.
-/
import Mathlib

open Finset

namespace TwoBaseLetters

variable {R : Type*} [CommRing R]

/-- qⁿ = 1 + (q − 1)[n]_q, in any commutative ring. -/
theorem pow_eq_ring (q : R) (n : ℕ) : q ^ n = 1 + (q - 1) * ∑ j ∈ range n, q ^ j := by
  linear_combination -(geom_sum_mul q n)

/-- Bradley's q-MZV letter with index a (that is, s = a + 1), at the point n, given the unit u = [n]_q. -/
noncomputable def letter (q : R) (u : Rˣ) (n a : ℕ) : R :=
  q ^ (a * n) * ((u⁻¹ : Rˣ) : R) ^ (a + 1)

/-- Hoffman–Ihara's §6.4 letter product, in any commutative ring. -/
theorem letter_mul (q : R) (n : ℕ) (u : Rˣ) (hu : (u : R) = ∑ j ∈ range n, q ^ j) (a b : ℕ) :
    letter q u n a * letter q u n b = letter q u n (a + b + 1) - (q - 1) * letter q u n (a + b) := by
  have key : q ^ n = 1 + (q - 1) * (u : R) := by rw [hu]; exact pow_eq_ring q n
  have huv : (u : R) * ((u⁻¹ : Rˣ) : R) = 1 := Units.mul_inv u
  unfold letter
  have e1 : q ^ ((a + b + 1) * n) = q ^ ((a + b) * n) * q ^ n := by
    rw [← pow_add]; congr 1; ring
  have e2 : q ^ (a * n) * q ^ (b * n) = q ^ ((a + b) * n) := by
    rw [← pow_add]; congr 1; ring
  rw [e1, show q ^ (a * n) * ((u⁻¹ : Rˣ) : R) ^ (a + 1) * (q ^ (b * n) * ((u⁻¹ : Rˣ) : R) ^ (b + 1))
      = (q ^ (a * n) * q ^ (b * n)) * ((u⁻¹ : Rˣ) : R) ^ (a + b + 2) by ring, e2]
  linear_combination (-(q ^ ((a + b) * n) * ((u⁻¹ : Rˣ) : R) ^ (a + b + 2))) * key
    + (-(q - 1) * q ^ ((a + b) * n) * ((u⁻¹ : Rˣ) : R) ^ (a + b + 1)) * huv

/-- With q = 1 + N, Hoffman–Ihara's (1 − q) becomes the jet −N. -/
theorem jet_letter_mul (N : R) (n : ℕ) (u : Rˣ) (hu : (u : R) = ∑ j ∈ range n, (1 + N) ^ j) (a b : ℕ) :
    letter (1 + N) u n a * letter (1 + N) u n b
      = letter (1 + N) u n (a + b + 1) - N * letter (1 + N) u n (a + b) := by
  have h := letter_mul (1 + N) n u hu a b
  rw [add_sub_cancel_left] at h
  exact h

/-- The triple product: d appears as a plain number at third order. -/
theorem jet_letter_triple (N k d : R) (hN : N ^ 2 = k * N - d) (n : ℕ) (u : Rˣ)
    (hu : (u : R) = ∑ j ∈ range n, (1 + N) ^ j) (a b c : ℕ) :
    letter (1 + N) u n a * letter (1 + N) u n b * letter (1 + N) u n c
      = letter (1 + N) u n (a + b + c + 2) - 2 * N * letter (1 + N) u n (a + b + c + 1)
        + (k * N - d) * letter (1 + N) u n (a + b + c) := by
  have h1 := jet_letter_mul N n u hu a b
  have h2 := jet_letter_mul N n u hu (a + b + 1) c
  have h3 := jet_letter_mul N n u hu (a + b) c
  rw [show a + b + 1 + c + 1 = a + b + c + 2 by ring, show a + b + 1 + c = a + b + c + 1 by ring] at h2
  rw [h1]
  linear_combination h2 - N * h3 + letter (1 + N) u n (a + b + c) * hN

/-- An eigen-shadow (any ring map) sends letters to letters. -/
theorem letter_map {S : Type*} [CommRing S] (φ : R →+* S) (q : R) (u : Rˣ) (n a : ℕ) :
    φ (letter q u n a) = letter (φ q) (Units.map (φ : R →* S) u) n a := by
  simp [letter, map_mul, map_pow, Units.coe_map_inv]

/-- The two-base jet ring exists: in AdjoinRoot (X² − kX + d), the root N satisfies N² = kN − d. -/
theorem root_sq (k d : ℝ) :
    let f : Polynomial ℝ := Polynomial.X ^ 2 - Polynomial.C k * Polynomial.X + Polynomial.C d
    (AdjoinRoot.root f) ^ 2 = AdjoinRoot.of f k * AdjoinRoot.root f - AdjoinRoot.of f d := by
  intro f
  have h := AdjoinRoot.eval₂_root f
  simp only [f, Polynomial.eval₂_add, Polynomial.eval₂_sub, Polynomial.eval₂_mul, Polynomial.eval₂_X_pow,
    Polynomial.eval₂_C, Polynomial.eval₂_X] at h
  linear_combination h

/-- The dial from the two shadow bases q₁ = 1 + α, q₂ = 1 + β. -/
theorem dial_from_bases (α β : R) :
    α + β = (1 + α) + (1 + β) - 2 ∧ α * β = ((1 + α) - 1) * ((1 + β) - 1) := by
  constructor <;> ring

/-- Holding one base at c (β = c − 1 fixed) traces the line d = βk − β², which touches the border
d = k²/4 only at k = 2β. -/
theorem hold_base_tangent (β k : ℝ) : k ^ 2 / 4 - (β * k - β ^ 2) = (k / 2 - β) ^ 2 := by ring

/-- The case Brennan asked about, c = 1/2: the line d = −k/2 − 1/4, tangent at (k, d) = (−1, 1/4). -/
theorem hold_half (α : ℝ) :
    let k := α + (1 / 2 - 1); let d := α * (1 / 2 - 1)
    d = -k / 2 - 1 / 4 ∧ k ^ 2 / 4 - d = (k / 2 + 1 / 2) ^ 2 := by
  intro k d
  constructor <;> (simp only [k, d]; ring)

end TwoBaseLetters
