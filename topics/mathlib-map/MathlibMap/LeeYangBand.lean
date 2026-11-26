/-
  STATEMENT ONLY. This file defines the Lee–Yang p-ball rule and its open band as propositions.
  It proves nothing about them and contains no `sorry`. It compiles against Mathlib's interval
  integral. The proof work belongs to the shape-dial thread (shape-dial-survey, section N3).

  The rule. A spin spread over the n-dimensional p-ball, read along one axis, has weight
  (1 − |x|^p)^a with a = (n − 1)/p. Its field function is
      G(s) = ∫₀¹ (1 − x^p)^a cos(s x) dx.
  Claim: every zero of G is real exactly when a ≤ p (that is n ≤ 1 + p²),
  or p is an even whole number.

  Status (shape-dial thread, unchecked here):
  * breaking half (a > p, p not even): proved
  * holding for a ≤ 1, any p: proved (Pólya 1918)
  * holding at p = 2 (Bessel): known
  * holding at p = 4, whole a: proved
  * open: p not even with 1 < a ≤ p, and even p ≥ 4 outside the cases above
-/
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace MathlibMap.LeeYang

open intervalIntegral

/-- The field function G(s) = ∫₀¹ (1 − x^p)^a cos(s x) dx, at complex field s. -/
noncomputable def fieldFn (p a : ℝ) (s : ℂ) : ℂ :=
  ∫ x in (0 : ℝ)..1, (((1 - x ^ p) ^ a : ℝ) : ℂ) * Complex.cos (s * x)

/-- Lee–Yang holds for the weight: every zero of the field function is real. -/
def Holds (p a : ℝ) : Prop :=
  ∀ s : ℂ, fieldFn p a s = 0 → s.im = 0

/-- p is an even whole number. -/
def EvenWhole (p : ℝ) : Prop := ∃ m : ℕ, p = 2 * m

/-- The full rule, as found numerically. -/
def Rule : Prop :=
  ∀ p a : ℝ, 1 ≤ p → 0 ≤ a → (Holds p a ↔ (a ≤ p ∨ EvenWhole p))

/-- The open band: the holding half where no proof is known yet. -/
def OpenBand : Prop :=
  (∀ p a : ℝ, 1 ≤ p → 1 < a → a ≤ p → ¬ EvenWhole p → Holds p a) ∧
  (∀ p a : ℝ, 4 ≤ p → EvenWhole p → 0 ≤ a → Holds p a)

end MathlibMap.LeeYang
