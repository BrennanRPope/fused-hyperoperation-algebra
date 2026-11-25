/-
  Chemistry 2, conservation of mass and charge
  =============================================

  Companion page and checks: ../index.html, ../checks.py (15 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Mass: the invariant mass of two parts is a law of cosines,
      M² = m₁² + m₂² + 2 m₁ m₂ c.
  c = 1 is Lavoisier (M = m₁ + m₂ exactly: the border, a flat triangle);
  c < 1 is a bound state (turning country, c = cos θ); c > 1 is a free pair
  in motion (two walls, c = cosh of the relative rapidity). Hydrogen's 1s
  state has sin θ = α; the checks file compares this with the measured
  ionisation energy.

  Charge: charges add (k = 0), and the fugacity z^Q multiplies, so its
  excess z^Q − 1 combines by the circle operation (called mult-add before 1.0). Every additive conservation law is the log of a
  circle-operation lift.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged at Topics/Chemistry/ChemConservation; it now uses Base's ◯. One exact twin removed.
-/
import FusedAlgebra.Base

namespace ChemConservation
open CircleOp

open scoped CircleOp

/-- Invariant mass squared of two parts with "cosine" c. -/
def mass2 (m1 m2 c : Rat) : Rat := m1 * m1 + m2 * m2 + 2 * m1 * m2 * c

/-- Lavoisier's border: c = 1 gives exactly (m₁ + m₂)². -/
theorem lavoisier (m1 m2 : Rat) : mass2 m1 m2 1 = (m1 + m2) * (m1 + m2) := by
  unfold mass2; grind

/-- The mass deficit squared is 2 m₁ m₂ (1 − c). -/
theorem deficit (m1 m2 c : Rat) : (m1 + m2) * (m1 + m2) - mass2 m1 m2 c = 2 * m1 * m2 * (1 - c) := by
  unfold mass2; grind

/-- With a binding energy B, M = m₁ + m₂ − B: B(2(m₁+m₂) − B) = 2 m₁ m₂ (1 − c), exactly. -/
theorem binding_exact (m1 m2 c B : Rat)
    (h : (m1 + m2 - B) * (m1 + m2 - B) = mass2 m1 m2 c) :
    B * (2 * (m1 + m2) - B) = 2 * m1 * m2 * (1 - c) := by
  unfold mass2 at h; grind

/-- Bound (c < 1) means the whole weighs less than its parts, for positive masses. -/
theorem bound_lighter (m1 m2 c : Rat) (h1 : 0 < m1) (h2 : 0 < m2) (hc : c < 1) :
    mass2 m1 m2 c < (m1 + m2) * (m1 + m2) := by
  have hpos : 0 < 2 * m1 * m2 * (1 - c) :=
    Rat.mul_pos (Rat.mul_pos (Rat.mul_pos (by decide) h1) h2) (by grind)
  have := deficit m1 m2 c
  grind

/-- Free pairs in motion (c > 1) weigh more than their rest masses. -/
theorem free_heavier (m1 m2 c : Rat) (h1 : 0 < m1) (h2 : 0 < m2) (hc : 1 < c) :
    (m1 + m2) * (m1 + m2) < mass2 m1 m2 c := by
  have hpos : 0 < 2 * m1 * m2 * (c - 1) :=
    Rat.mul_pos (Rat.mul_pos (Rat.mul_pos (by decide) h1) h2) (by grind)
  have := deficit m1 m2 c
  grind

/-- The cosine law reads as the circle operation with dial k in the masses' squares: with
    a = m₁², b = m₂² and the corner taken through √(ab) = m₁m₂, the extra
    piece beyond plain addition is 2c·m₁m₂. At c = 1 it is the perfect square. -/
theorem cosine_extra (m1 m2 c : Rat) : mass2 m1 m2 c - (m1 * m1 + m2 * m2) = 2 * c * (m1 * m2) := by
  unfold mass2; grind

/-- Small-angle form: if c = 1 − t/2 (t = θ²), the squared deficit is m₁m₂·t. -/
theorem small_angle (m1 m2 t : Rat) :
    (m1 + m2) * (m1 + m2) - mass2 m1 m2 (1 - t / 2) = m1 * m2 * t := by
  unfold mass2; grind

/-! ## Charge: additive charges, multiplicative fugacities -/

section Charge
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: fugacity_mul (fugacities multiply as charges add) was an exact twin of
--   CircleOp.pw_add and is removed.

/-- So the excess z^Q − 1 combines by the circle operation. -/
theorem fugacity_circ (z : R) (Q1 Q2 : Nat) :
    pw z (Q1 + Q2) - 1 = (pw z Q1 - 1) ◯ (pw z Q2 - 1) := by
  rw [pw_add]; unfold circ; grind

/-- A neutral system (charge 0) has fugacity 1: excess 0, the circle operation's identity. -/
theorem neutral_identity (z : R) : pw z 0 - 1 = 0 := by simp only [pw]; grind

end Charge

/-- Rotation (the phase e^{iQφ}) in Gaussian-rational coordinates: if both
    factors lie on the unit circle, so does their product. The excess of a
    point of the circle is never the wall −1 unless the point is 0. -/
theorem circle_closed (a b c d : Rat) (h1 : a * a + b * b = 1) (h2 : c * c + d * d = 1) :
    (a * c - b * d) * (a * c - b * d) + (a * d + b * c) * (a * d + b * c) = 1 := by
  have : (a * c - b * d) * (a * c - b * d) + (a * d + b * c) * (a * d + b * c)
      = (a * a + b * b) * (c * c + d * d) := by grind
  rw [this, h1, h2]; grind

end ChemConservation
