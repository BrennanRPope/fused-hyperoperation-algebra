/-
  Chemistry 12, organic chemistry
  ===============================

  Companion page and checks: ../index.html, ../checks.py (16 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Rings: characteristic polynomials of the 4- and 6-rings, Hückel and Möbius.
  Hammett: rate ratios multiply, so ratio − 1 combines by the circle operation (called mult-add before 1.0). Stereoselectivity:
  with e(r) = (r − 1)/(r + 1), matched selectivities compose by the landscape
  law at (0, −1), (e₁ + e₂)/(1 + e₁e₂), the velocity-addition law; ee = ±1 are
  its walls, and Horeau's duplication is e ⊕ e.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  LIBRARY NOTE: first compiled in the shared library, after ChemAcids. Two exact twins
  removed (hammett_compose, route_losses); each is marked where it stood. No proof step or
  statement changed.
-/
import FusedAlgebra.Base

namespace ChemOrganic
open CircleOp

open scoped CircleOp

/-- Benzene (Hückel 6-ring): x⁶ − 6x⁴ + 9x² − 4 = (x − 2)(x + 2)(x − 1)²(x + 1)². -/
theorem benzene_charpoly (x : Rat) :
    x^6 - 6 * x^4 + 9 * x^2 - 4 = (x - 2) * (x + 2) * ((x - 1) * (x - 1)) * ((x + 1) * (x + 1)) := by
  grind

/-- Cyclobutadiene (Hückel 4-ring): x⁴ − 4x² = x²(x − 2)(x + 2), a doubly degenerate zero level. -/
theorem cyclobutadiene_charpoly (x : Rat) : x^4 - 4 * x^2 = (x * x) * ((x - 2) * (x + 2)) := by
  grind

/-- Möbius 4-ring: x⁴ − 4x² + 4 = (x² − 2)², levels ±√2, each twice. -/
theorem mobius4_charpoly (x : Rat) : x^4 - 4 * x^2 + 4 = (x * x - 2) * (x * x - 2) := by
  grind

/-- Benzene's π energy (in β): 2·2 + 4·1 = 8, two more than three ethylenes (6). -/
theorem benzene_resonance : 2 * 2 + 4 * 1 - 3 * 2 = (2 : Nat) := by decide

-- LIBRARY: hammett_compose (Hammett: σ adds, so rate ratios multiply and ratio − 1 combines by the circle operation) was an exact twin of InfoFused.campbell and is removed.

/-- Matched selectivities: with e = (r − 1)/(r + 1), multiplying the ratios composes ee by
    (e₁ + e₂)/(1 + e₁e₂). Cleared form, with e₁, e₂ given by their defining equations. -/
theorem ee_compose (r1 r2 e1 e2 : Rat)
    (h1 : e1 * (r1 + 1) = r1 - 1) (h2 : e2 * (r2 + 1) = r2 - 1) :
    (e1 + e2) * (r1 * r2 + 1) = (1 + e1 * e2) * (r1 * r2 - 1) := by
  have a1 : e1 * r1 = r1 - 1 - e1 := by grind
  have a2 : e2 * r2 = r2 - 1 - e2 := by grind
  grind

/-- The ee law is the landscape law at (k, q) = (0, −1): numerator e₁ + e₂, denominator 1 + e₁e₂. -/
theorem ee_is_landscape (e1 e2 : Rat) :
    circK 0 e1 e2 = e1 + e2 ∧ 1 - (-1) * e1 * e2 = 1 + e1 * e2 := by
  unfold circK; constructor <;> grind

/-- Its discriminant k² − 4q is 4: the two-walls country. -/
theorem ee_two_walls : (0 : Rat) * 0 - 4 * (-1) = 4 := by grind

/-- ee = 1 is a wall: a perfectly selective step absorbs any other (cleared: 1 + b = 1·(1 + 1·b)). -/
theorem ee_wall (b : Rat) : 1 + b = 1 * (1 + 1 * b) := by grind

/-- Horeau duplication is e ⊕ e: numerator 2e, denominator 1 + e². -/
theorem horeau (e : Rat) : circK 0 e e = 2 * e ∧ 1 - (-1) * e * e = 1 + e * e := by
  unfold circK; constructor <;> grind

-- LIBRARY: route_losses (Linear routes: losses fold at k = −1) was an exact twin of ChemStoichiometry.loss_two and is removed.

end ChemOrganic
