/-
  The binary pulsar through "twist = curvature × area"
  ====================================================

  Companion checks: ../pulsar_checks.py (13 checks: the coefficient, the
  test-gyroscope limit, the orbit average 1/(1 − e²), the double pulsar and
  Hulse–Taylor numbers, and the field-strength estimate).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports Base only.
  Compiled clean in the shared library (Topics/BinaryPulsar), unchanged.

  A pulsar of mass m₁ orbiting a companion of mass m₂ (total M) turns its spin
  axis by three twists per orbit, each weighted by the companion's mass m₂:
    * the velocity twist (Thomas-type, q = −1/c²):    ½ · m₂/M
    * the space twist (curvature of the space slice):  1 · m₂/M
    * the companion's drag (gravitomagnetic circulation of the companion
      moving relative to the centre of mass):         2 · m₁/M
  Their sum, times m₂, is the standard spin-orbit coefficient m₂(4m₁ + 3m₂)/(2M).
  In the limit m₁ → 0 this is (3/2)·m₂, the Gravity Probe B geodetic law.

  Numbers are in degrees per year (double pulsar B, Breton et al. 2008).
-/
import FusedAlgebra.Base

namespace BinaryPulsar

/-! ## 1. Three twists make the standard coefficient -/

/-- Twist weights per unit of the companion's mass: velocity ½·m₂/M, space m₂/M, drag 2·m₁/M. -/
def twistSum (m₁ m₂ : Rat) : Rat :=
  1 / 2 * (m₂ / (m₁ + m₂)) + m₂ / (m₁ + m₂) + 2 * (m₁ / (m₁ + m₂))

/-- **The three twists sum to the standard spin-orbit coefficient** m₂(4m₁ + 3m₂)/(2M). -/
theorem three_twists (m₁ m₂ : Rat) (hM : m₁ + m₂ ≠ 0) :
    m₂ * twistSum m₁ m₂ = m₂ * (4 * m₁ + 3 * m₂) / (2 * (m₁ + m₂)) := by
  unfold twistSum
  have e : 1 / 2 * (m₂ / (m₁ + m₂)) + m₂ / (m₁ + m₂) + 2 * (m₁ / (m₁ + m₂))
      = (4 * m₁ + 3 * m₂) / (2 * (m₁ + m₂)) := by grind
  rw [e]; grind

/-- Polynomial form (no division): 2M · (½m₂ + m₂ + 2m₁) = 4m₁ + 3m₂ times M/M. -/
theorem three_twists_poly (m₁ m₂ : Rat) :
    (m₂ + 2 * m₂ + 4 * m₁) = 4 * m₁ + 3 * m₂ := by grind

/-- **Test-gyroscope limit** m₁ = 0: the coefficient is (3/2)·m₂, the Gravity Probe B law. -/
theorem test_limit (m₂ : Rat) (h : m₂ ≠ 0) :
    m₂ * (4 * 0 + 3 * m₂) / (2 * (0 + m₂)) = 3 / 2 * m₂ := by grind

/-- In that limit the drag twist vanishes: a light partner barely moves the heavy companion. -/
theorem drag_vanishes (m₂ : Rat) : 2 * ((0 : Rat) / (0 + m₂)) = 0 := by grind

/-- Velocity and space twists keep their one-to-two ratio, as in Gravity Probe B. -/
theorem velocity_space_ratio (w : Rat) (h : w ≠ 0) : (1 / 2 * w) / w = 1 / 2 := by grind

/-! ## 2. The double pulsar numbers (degrees per year) -/

/-- Twist bookkeeping gives 5.07; measured 4.77 with +0.66 / −0.65. Within one sigma. -/
theorem double_pulsar_within_error :
    (507 / 100 : Rat) - 477 / 100 < 66 / 100 ∧ (477 / 100 : Rat) - 507 / 100 < 65 / 100 := by
  grind

/-- The three parts: velocity 0.75, space 1.51, drag 2.81, which add to 5.07. -/
theorem parts_add : (75 / 100 : Rat) + 151 / 100 + 281 / 100 = 507 / 100 := by grind

/-- **Without the companion's drag** the prediction is 2.26, more than 3σ below 4.77. -/
theorem drag_needed : (477 / 100 : Rat) - 226 / 100 > 3 * (65 / 100) := by grind

/-- Without the velocity twist the prediction is 4.32, still within one sigma: not isolated. -/
theorem velocity_not_isolated : (477 / 100 : Rat) - 432 / 100 < 65 / 100 := by grind

/-! ## 3. Why the exact fusion cannot show here -/

/-- The ◯ corner is the product of two first-order twists, each of size about M/r ≈ 4.35e-6,
    so it is below 1e-10 of a turn, far under the measurement's 14% error. -/
theorem corner_tiny :
    (435 / 100000000 : Rat) * (435 / 100000000) < 1 / 10000000000 := by grind

end BinaryPulsar
