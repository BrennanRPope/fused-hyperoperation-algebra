/-
  Gravity Probe B through "twist = curvature × area"
  ==================================================

  Companion checks: ../gpb_checks.py (21 checks: exact Schwarzschild parallel transport,
  the space slice's curvature and holonomy, the GP-B numbers, frame dragging by Stokes).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports Base only.
  Not compiled in the session that wrote it (no Lean there).
  LIBRARY NOTE: first compiled in the shared library, at FusedAlgebra/Topics/GravityProbeB.lean
  after TwistedBoosts; copied unchanged, every grind step closed. No exact twins: fuse_weight_minus_one
  is a relative of MultAdd.wAdd_either, stated on multAdd of negated fractions.

  A gyroscope in a circular orbit turns, per orbit, by two twists:
    * the space twist: the curvature of the space slice enclosed by the orbit,
      fraction p_s = 1 − √(1 − 2M/r) of a turn (Flamm's surface, exact);
    * the velocity twist: the Thomas-type twist of velocity space (q = −1/c²),
      fraction p_v = 1 − 1/γ with γ² = (r − 2M)/(r − 3M) (exact factor; equals
      Thomas's (γ − 1) to first order).
  They combine by multiplying the factors: 1 − p = (1 − p_s)(1 − p_v), so
  −p = (−p_s) ✳ (−p_v): **the precession fractions mult-add**. The total is the
  known exact geodetic precession 1 − √(1 − 3M/r) per orbit.

  Here x = M/r (in units G = c = 1). Lean works with the squares of the factors,
  since square roots of rationals are not rational; the roots are in the checks.
-/
import FusedAlgebra.Base

namespace GravityProbeB
open MultAdd

local infixl:65 " ✳ " => multAdd

/-! ## 1. The two factors and their product -/

/-- Local orbital speed squared, relative to static observers: v² = x/(1 − 2x). -/
def v2 (x : Rat) : Rat := x / (1 - 2 * x)

/-- γ² = 1/(1 − v²) = (1 − 2x)/(1 − 3x). -/
theorem gamma_sq (x : Rat) (h2 : 1 - 2 * x ≠ 0) (h3 : 1 - 3 * x ≠ 0) :
    1 / (1 - v2 x) = (1 - 2 * x) / (1 - 3 * x) := by
  unfold v2
  have e : 1 - x / (1 - 2 * x) = (1 - 3 * x) / (1 - 2 * x) := by grind
  rw [e]; grind

/-- **The factors multiply to the known total** (squared): (1 − 2x)·(1/γ²) = 1 − 3x. -/
theorem factors_multiply (x : Rat) (h2 : 1 - 2 * x ≠ 0) :
    (1 - 2 * x) * (1 - v2 x) = 1 - 3 * x := by
  unfold v2; grind

/-- **The precession fractions mult-add**: if 1 − p = (1 − a)(1 − b) then −p = (−a) ✳ (−b). -/
theorem fractions_fuse (p a b : Rat) (h : 1 - p = (1 - a) * (1 - b)) :
    -p = (-a) ✳ (-b) := by
  show -p = -a + -b + -a * -b; grind

/-- And back: the mult-add form gives the product of factors. -/
theorem fuse_factors (a b : Rat) : 1 + ((-a) ✳ (-b)) = (1 - a) * (1 - b) := by
  show 1 + (-a + -b + -a * -b) = (1 - a) * (1 - b); grind

/-- The same law is the dial's weight −1: p = a + b − ab ("either twist"). -/
theorem fuse_weight_minus_one (a b : Rat) : -((-a) ✳ (-b)) = a + b + (-1) * a * b := by
  show -(-a + -b + -a * -b) = a + b + (-1) * a * b; grind

/-! ## 2. Weak field: one part velocity, two parts space -/

/-- First-order parts: space 1·x, velocity ½·x, total 3/2·x. -/
theorem weak_split (x : Rat) : x + x / 2 = 3 / 2 * x := by grind

/-- The fused total differs from the plain sum only by the corner a·b (second order). -/
theorem corner_is_second_order (a b : Rat) : (a + b) - (-((-a) ✳ (-b))) = a * b := by
  show (a + b) - (-(-a + -b + -a * -b)) = a * b; grind

/-- In PPN form the geodetic coefficient is γ_PPN + ½; general relativity's γ_PPN = 1 gives 3/2. -/
theorem ppn_geodetic : (1 : Rat) + 1 / 2 = 3 / 2 := by grind

/-- The velocity part alone is one third, the space part alone two thirds. -/
theorem thirds : (1 / 2 : Rat) / (3 / 2) = 1 / 3 ∧ (1 : Rat) / (3 / 2) = 2 / 3 := by grind

/-! ## 3. The Gravity Probe B numbers (milliarcseconds per year) -/

/-- Twist bookkeeping at r = 7027.4 km gives 6603.9; GP-B measured 6601.8 ± 18.3. -/
theorem geodetic_within_error :
    (66039 / 10 : Rat) - 66018 / 10 < 183 / 10 ∧ (66018 / 10 : Rat) - 66039 / 10 < 183 / 10 := by
  grind

/-- Space part 4402.6 + velocity part 2201.3 = 6603.9. -/
theorem geodetic_parts : (44026 / 10 : Rat) + 22013 / 10 = 66039 / 10 := by grind

/-- Either part alone misses the measurement by more than 100 standard errors. -/
theorem parts_alone_ruled_out :
    100 * (183 / 10 : Rat) < 66018 / 10 - 22013 / 10 ∧ 100 * (183 / 10 : Rat) < 66018 / 10 - 44026 / 10 := by
  grind

/-- Frame dragging: about 40.8 against the measured 37.2 ± 7.2. -/
theorem frame_dragging_within_error : (408 / 10 : Rat) - 372 / 10 < 72 / 10 := by grind

/-! ## 4. Frame dragging over a polar orbit -/

/-- 3 sin² − 1 averaged at the quarter points (sin² = 0, 1, 0, 1) is ½. -/
theorem polar_average_four : ((3 * 0 - 1) + (3 * 1 - 1) + (3 * 0 - 1) + (3 * 1 - 1) : Rat) / 4 = 1 / 2 := by
  grind

/-- And at the eighth points (sin² = 0, ½, 1, ½, …) it is ½ too. -/
theorem polar_average_eight :
    ((3 * 0 - 1) + (3 * (1/2) - 1) + (3 * 1 - 1) + (3 * (1/2) - 1)
      + (3 * 0 - 1) + (3 * (1/2) - 1) + (3 * 1 - 1) + (3 * (1/2) - 1) : Rat) / 8 = 1 / 2 := by
  grind

/-- Any orbit point and its quarter-turn partner have sin² + cos² = 1, so a pair averages
    (3·s − 1 + 3·(1 − s) − 1)/2 = ½ whatever s is: the average is ½ for every
    set of points closed under quarter turns. -/
theorem polar_average_pair (s : Rat) : ((3 * s - 1) + (3 * (1 - s) - 1)) / 2 = 1 / 2 := by grind

end GravityProbeB
