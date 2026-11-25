/-
  Fused Hyperoperation Algebra: fill, linger, spill at the landscape's border.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/LeeYangSpill.lean (after LeeYangCluster).

  Repeat one landscape law with a fixed step x:  a ↦ a ⊕ x = (a + x + k·a·x)/(1 − q·a·x).
  The step is a Möbius map with trace 2 + kx and determinant 1 + kx + qx². Its discriminant is
  x²·(k² − 4q), so for every step size it is
    • hyperbolic (two walls, the orbit settles) when 4q < k²,
    • parabolic (one double wall) on the border 4q = k²,
    • a turn (the orbit lingers, spills through infinity, comes round) when 4q > k².
  How far a moves in one step is x·(1 + k·a + q·a²)/(1 − q·a·x), and
  4q·(1 + k·a + q·a²) = (2q·a + k)² + (4q − k²): on the turning side the slowest step is set by
  the distance to the border. That is the bottleneck: linger ~ 1/√(4q − k²).
-/
import FusedAlgebra.Base

namespace LeeYangSpill

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- The step's discriminant: trace² − 4·det = x²·(k² − 4q), the same sign for every step size. -/
theorem step_disc (k q x : R) :
    (2 + k * x) * (2 + k * x) - 4 * (1 + k * x + q * x * x) = x * x * (k * k - 4 * q) := by grind

/-- One step moves a by x·(1 + k·a + q·a²), over the denominator 1 − q·a·x. -/
theorem step_gap (k q x a : R) :
    (a + x + k * a * x) - a * (1 - q * a * x) = x * (1 + k * a + q * a * a) := by grind

/-- Completing the square: the step's size is a square plus the distance to the border. -/
theorem wall_square (k q a : R) :
    4 * q * (1 + k * a + q * a * a) = (2 * q * a + k) * (2 * q * a + k) + (4 * q - k * k) := by grind

/-- Two walls: when k² − 4q = m², the step size factors through the two walls, where it is 0. -/
theorem walls_factor (k q a m : R) (h : k * k - 4 * q = m * m) :
    4 * q * (1 + k * a + q * a * a) = (2 * q * a + k - m) * (2 * q * a + k + m) := by grind

/-- The border (k = 2u, q = u²): n steps of x from 0 land at n·x / (1 − (n − 1)·u·x).
    Stated for a general t in place of n, cross-multiplied: the step takes the t-th point to the
    (t+1)-th. The denominator is linear in t: it passes infinity once, then fills toward the
    double wall −1/u like 1/t, forever. -/
theorem border_orbit (u x t : R) :
    ((t * x) + x * (1 - (t - 1) * u * x) + (2 * u) * (t * x) * x) * (1 - t * u * x)
      = ((t + 1) * x) * ((1 - (t - 1) * u * x) - (u * u) * (t * x) * x) := by grind

end Ring

/-- A square is never negative. -/
theorem sq_nonneg' (a : Int) : 0 ≤ a * a := by
  rcases Int.le_total 0 a with h | h
  · exact Int.mul_nonneg h h
  · have := Int.mul_nonneg (Int.neg_nonneg_of_nonpos h) (Int.neg_nonneg_of_nonpos h)
    simpa [Int.neg_mul_neg] using this

/-- On the turning side the slowest step is bounded below by the distance to the border:
    4q·(1 + k·a + q·a²) ≥ 4q − k² for every a. -/
theorem slowest_step (k q a : Int) : 4 * q * (1 + k * a + q * a * a) ≥ 4 * q - k * k := by
  have h := wall_square k q a
  have hs := sq_nonneg' (2 * q * a + k)
  omega

/-- On the turning side the step never stalls: if q > 0 and 4q > k², the step size is positive. -/
theorem turning_never_stalls (k q a : Int) (hq : 0 < q) (hd : k * k < 4 * q) :
    0 < 1 + k * a + q * a * a := by
  have h := slowest_step k q a
  apply Int.lt_of_not_ge
  intro hle
  have : 4 * q * (1 + k * a + q * a * a) ≤ 0 :=
    Int.mul_nonpos_of_nonneg_of_nonpos (by omega) hle
  omega

/-- Tangent addition (k = 0, q = 1) with step 1: a ↦ (a + 1)/(1 − a). Its matrix [[1,1],[−1,1]]
    squared is [[0,2],[−2,0]] and to the fourth is −4·I, a scalar: every orbit spills through
    infinity and returns to its start in exactly four steps (tan(θ + π/4)). -/
theorem tan_square : (1 * 1 + 1 * (-1) : Int) = 0 ∧ (1 * 1 + 1 * 1 : Int) = 2 ∧
    ((-1) * 1 + 1 * (-1) : Int) = -2 ∧ ((-1) * 1 + 1 * 1 : Int) = 0 := by decide

theorem tan_fourth : (0 * 0 + 2 * (-2) : Int) = -4 ∧ (0 * 2 + 2 * 0 : Int) = 0 ∧
    ((-2) * 0 + 0 * (-2) : Int) = 0 ∧ ((-2) * 2 + 0 * 0 : Int) = -4 := by decide

end LeeYangSpill
