/-
  Fused Hyperoperation Algebra: the elliptic modulus (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib, no `sorry`), compiled in the shared library on Nov 26 as FusedAlgebra.Topics.EllipticModulus
  (the old name MultAdd.multAdd updated to CircleOp.circ).
  Polynomial identities behind ../checks.py, with denominators cleared (Int).
  Notation: m = Jacobi's modulus (parameter), κ = 1 − 2m the curvature reading of the slope chart,
  ε = the t⁴ coefficient under the square root of the law's logarithm, s = height.
-/
import FusedAlgebra.Base

namespace EllipticModulus

/-- Wick rotation m ↦ 1 − m flips the curvature reading κ = 1 − 2m. -/
theorem wick_flips_kappa (m : Int) : 1 - 2 * (1 - m) = -(1 - 2 * m) := by grind

/-- The slope chart's law: G'² = 1 + 2κ t² + t⁴. With x = sn², t² = x(1 − m x)/(1 − x) and
    G' = (1 − 2m x + m x²)/(1 − x), this is the identity below after clearing (1 − x)². -/
theorem slope_log_quartic (m x : Int) :
    (1 - 2 * m * x + m * x * x) * (1 - 2 * m * x + m * x * x)
      = (1 - x) * (1 - x) + 2 * (1 - 2 * m) * (x * (1 - m * x)) * (1 - x)
        + (x * (1 - m * x)) * (x * (1 - m * x)) := by grind

/-- Our k = 0 conic laws are the perfect squares: (1 + d t²)² = 1 + 2d t² + d² t⁴. -/
theorem conic_is_square (d t : Int) :
    (1 + d * t * t) * (1 + d * t * t) = 1 + 2 * d * t * t + d * d * t * t * t * t := by grind

/-- Fingerprint of an odd chart with log' = (1 + 2κ t² + ε t⁴)^(−1/2): the generator is
    u + (κ/3) u³ + ((3ε + κ²)/30) u⁵ + …, and 90·(3/2·c₃² − c₅) = 12κ² − 9ε while 90·c₃² = 10κ²,
    so R = 3/2 − c₅/c₃² = 6/5 − 9ε/(10κ²). -/
theorem fingerprint_general (κ ε : Int) :
    15 * κ * κ - 3 * (3 * ε + κ * κ) = 12 * κ * κ - 9 * ε := by grind

/-- In the slope chart ε = 1 and κ = 1 − 2m; the numerator 12κ² − 9 is 3(16m² − 16m + 1) ... -/
theorem slope_numerator (m : Int) :
    12 * (1 - 2 * m) * (1 - 2 * m) - 9 = 3 * (16 * m * m - 16 * m + 1) := by grind

/-- ... and it is Wick-symmetric, unchanged by m ↦ 1 − m. -/
theorem slope_numerator_wick (m : Int) :
    16 * (1 - m) * (1 - m) - 16 * (1 - m) + 1 = 16 * m * m - 16 * m + 1 := by grind

/-- Heights scale the fingerprint by 1/s: the height-s map has c₃ = s·a and
    c₅ = s·b + (3/2)s(s − 1)a², and 3(sa)² − 2c₅ = s(3a² − 2b). -/
theorem height_scales_fingerprint (a b s : Int) :
    3 * (s * a) * (s * a) - (2 * s * b + 3 * s * (s - 1) * a * a) = s * (3 * a * a - 2 * b) := by grind

/-- Normal form: if 10κ²·R = 12κ² − 9, then κ²(12 − 10R) = 9, so every fingerprint but 6/5 has an elliptic chart. -/
theorem normal_form_inverse (κ2 R : Int) (h : 10 * κ2 * R = 12 * κ2 - 9) :
    κ2 * (12 - 10 * R) = 9 := by grind

/-- The wall in the slope chart is the only one: G' numerator 1 − 2m x + m x² = (1 − m x)² + m(1 − m)x²,
    a sum of squares (positive) for 0 ≤ m < 1 and 0 ≤ x ≤ 1. -/
theorem slope_monotone_sos (m x : Int) :
    1 - 2 * m * x + m * x * x = (1 - m * x) * (1 - m * x) + m * (1 - m) * x * x := by grind

/-- Pendulum: with c² + s² = 1 (c = cos(θ₀/2), s = sin(θ₀/2)), κ = 1 − 2m = 1 − 2s² = c² − s² = cos θ₀. -/
theorem pendulum_kappa (c s : Int) (h : c * c + s * s = 1) : 1 - 2 * (s * s) = c * c - s * s := by grind

/-- Theta as the lift of a fold: (1 + a)(1 + b)(1 + c) = 1 + ((a ◯ b) ◯ c), three steps of the Tate orbit. -/
theorem theta_fold_three (a b c : Int) :
    (1 + a) * (1 + b) * (1 + c) = 1 + CircleOp.circ (CircleOp.circ a b) c := by
  unfold CircleOp.circ; grind

end EllipticModulus
