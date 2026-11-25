/-
  Chemistry 9, rates and Arrhenius
  ================================

  Companion page and checks: ../index.html, ../checks.py (19 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  First order: survivals multiply over consecutive times, so conversions
  compose by circK (−1). Second order: X/(1 − X) adds, so conversions compose
  by the landscape law (a + b − 2ab)/(1 − ab), the point (−2, 1) on the
  border k² = 4q. Order n osculates (−n, n(n−1)/2), whose discriminant is
  n(2 − n). The Tsallis q-logarithm turns products into the dial with
  k = 1 − q; here it is checked in its cleared polynomial form.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged at Topics/Chemistry/ChemRates; it now uses Base's ◯. Three exact twins removed (see the
  LIBRARY comments).
-/
import FusedAlgebra.Base

namespace ChemRates
open CircleOp

open scoped CircleOp

-- LIBRARY: first_order_compose (First order: conversions over consecutive times compose at k = −1) was an exact twin of ThermoFused.ratio_series and is removed.

/-- Parallel channels at a fixed time: survivals multiply too. -/
theorem parallel_channels (sa sb : Rat) : circK (-1) (1 - sa) (1 - sb) = 1 - sa * sb := by
  rw [circK_either]; grind

/-- The landscape law (a + b + k·ab)/(1 − q·ab), written cleared: numerator and denominator. -/
def landNum (k a b : Rat) : Rat := a + b + k * a * b
def landDen (q a b : Rat) : Rat := 1 - q * a * b

/-- Second order: with u = a/(1−a), v = b/(1−b), the sum u + v equals c/(1−c) where
    c = (a + b − 2ab)/(1 − ab). Cleared: (u + v)(1 − c) = c with all denominators multiplied out. -/
theorem second_order_lens (a b : Rat) :
    (a * (1 - b) + b * (1 - a)) * (landDen 1 a b - landNum (-2) a b)
      = landNum (-2) a b * ((1 - a) * (1 - b)) := by
  unfold landNum landDen; grind

/-- (−2, 1) is on the border: k² = 4q. -/
theorem second_order_border : (-2 : Int) ^ 2 = 4 * 1 := by decide

/-- First order is the point (−1, 0): the landscape law reduces to circK (−1). -/
theorem first_order_point (a b : Rat) : landNum (-1) a b = circK (-1) a b := by
  unfold landNum; rw [circK_def]

/-- The osculating point (−n, n(n−1)/2) has discriminant n(2 − n) (times 1, over Rat). -/
theorem osculating_disc (n : Rat) : (-n) * (-n) - 4 * (n * (n - 1) / 2) = n * (2 - n) := by grind

/-- Order 0, 1, 2 give q = 0, 0, 1 on that curve. -/
theorem osculating_small : (0 : Rat) * (0 - 1) / 2 = 0 ∧ (1 : Rat) * (1 - 1) / 2 = 0 ∧ (2 : Rat) * (2 - 1) / 2 = 1 := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

-- LIBRARY: qlog_two (q-logarithm with q = 2 (order-2 survival): ln₂ x = 1 − 1/x. Cleared form of ln₂(xy) = ln₂ x ◯[−1] ln₂ y: (1 − 1/(xy)) = circK (−1) (1 − 1/x) (1 − 1/y)) was an exact twin of InformationFull.qlog_neg_one and is removed.

-- LIBRARY: qlog_zero (q-logarithm with q = 0: ln₀ x = x − 1, and ln₀(xy) = ln₀ x ◯ ln₀ y (k = 1 − 0 = 1)) was an exact twin of ThermoFused.partition_square (also ThermoFused.hess) and is removed.

/-- Order ½ empties in finite time: s^(1/2) = 1 − kt/2 hits 0 at t = 2/k. Cleared: k·(2/k) = 2. -/
theorem half_order_wall (k : Rat) (hk : k ≠ 0) : 1 - k * (2 / k) / 2 = 0 := by grind

/-- Arrhenius ratios over successive steps multiply (as exponents, they add). -/
theorem arrhenius_steps (r1 r2 r3 : Rat) (h1 : r1 ≠ 0) (h2 : r2 ≠ 0) :
    (r2 / r1) * (r3 / r2) = r3 / r1 := by grind

end ChemRates
