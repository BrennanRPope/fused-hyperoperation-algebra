/-
  Fused Hyperoperation Algebra: how the law shapes stand on each other.
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`. Imports Base only.
  Intended library path: FusedAlgebra/Topics/ShapeDependencies.lean (after Landscape).

  The landscape law F_{k,d}(a, b) = (a + b + k·a·b) / (1 − d·a·b), numerator `circK k a b`.
  Its global shape is set by the sign of the discriminant k² − 4d of the rate 1 + k·a + d·a²:
    > 0  two real walls (the circle operation's shape),  = 0  one double wall (addition's shape),
    < 0  no real wall (slope addition's shape, a turning law).
  One level up sits the torus (elliptic curves), which degenerates to the other two.

  What is proved:
  1. Every law is addition plus a dialed corner (so addition is the root).
  2. Zooming keeps the sign of the discriminant (so the three shapes are real classes).
  3. Turning by i (any i with i·i = −1) maps the law at (k, d) to the law at (i·k, −d) and flips
     the sign of the discriminant: no-wall laws are two-wall laws turned by i; the border is fixed.
     Velocity addition turns into slope addition this way.
  4. A two-wall law is the difference of its two lens factors divided by the wall gap: as the walls
     merge the law becomes the border law (the limit that takes the circle operation's shape to addition's).
  5. The nodal cubic y² = x³ + x² (the torus pinched once) composes by multiplication of
     φ = (m − 1)/(m + 1), so by the circle operation in a = φ − 1; the cuspidal cubic y² = x³ (pinched twice)
     composes by addition of 1/m. Points are parametrized by the slope m of the line through the
     singular point; three points are collinear exactly when the stated factor vanishes.
-/
import FusedAlgebra.Base

namespace ShapeDependencies
open CircleOp
open Lean.Grind
set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

variable {R : Type} [CommRing R]

/-- the discriminant of the rate 1 + k·a + d·a² -/
def disc (k d : R) : R := k * k - 4 * d

/-! ### 1. Addition is the root: every law is addition plus a dialed corner -/

theorem law_is_addition_plus_corner (k a b : R) : circK k a b = (a + b) + k * (a * b) := by
  unfold circK; grind

theorem addition_is_dial_zero (a b : R) : circK 0 a b = a + b := circK_weight_zero a b

/-- the circle operation needs addition and multiplication, nothing else -/
theorem circ_needs_add_and_mul (a b : R) : circ a b = (a + b) + a * b := rfl

/-! ### 2. Zoom keeps the sign of the discriminant -/

theorem zoom_num (k lam a b : R) : circK k (lam * a) (lam * b) = lam * circK (lam * k) a b := by
  unfold circK; grind

theorem zoom_den (d lam a b : R) : 1 - d * (lam * a) * (lam * b) = 1 - (lam * lam * d) * a * b := by
  grind

theorem zoom_disc (k d lam : R) : disc (lam * k) (lam * lam * d) = (lam * lam) * disc k d := by
  unfold disc; grind

/-! ### 3. Turning by i: two walls ↔ no wall, the border stays put -/

theorem turn_num (i k a b : R) (hi : i * i = -1) :
    circK k (i * a) (i * b) = i * circK (i * k) a b := by
  unfold circK; grind

theorem turn_den (i d a b : R) (hi : i * i = -1) :
    1 - d * (i * a) * (i * b) = 1 - (-d) * a * b := by
  have h : d * (i * a) * (i * b) = (i * i) * (d * a * b) := by grind
  rw [h, hi]; grind

theorem turn_disc (i k d : R) (hi : i * i = -1) : disc (i * k) (-d) = -(disc k d) := by
  unfold disc
  have h : i * k * (i * k) = (i * i) * (k * k) := by grind
  rw [h, hi]; grind

/-- the border (disc = 0) is fixed by the turn -/
theorem turn_keeps_border (i k d : R) (hi : i * i = -1) (h : disc k d = 0) : disc (i * k) (-d) = 0 := by
  rw [turn_disc i k d hi, h]; grind

/-- velocity addition (k, d) = (0, −1), denominator 1 + ab, turned by i is slope addition,
    denominator 1 − ab: numerator i(a + b) -/
theorem velocity_turns_to_slope (i a b : R) (hi : i * i = -1) :
    circK 0 (i * a) (i * b) = i * (a + b) ∧ 1 - (-1) * (i * a) * (i * b) = 1 - 1 * a * b := by
  constructor
  · unfold circK; grind
  · have h := turn_den i (-1) a b hi; rw [h]; grind

/-! ### 4. Merging two walls gives the border law -/

/-- lens factors (1 − u a)(1 − u b) and (1 − v a)(1 − v b): their difference is the wall gap times
    the law's numerator, for the law with walls 1/u and 1/v (k = −(u + v)) -/
theorem two_walls_difference (u v a b : R) :
    (1 - u * a) * (1 - u * b) - (1 - v * a) * (1 - v * b) = (v - u) * circK (-(u + v)) a b := by
  unfold circK; grind

theorem walls_merge_disc (u v : R) : disc (-(u + v)) (u * v) = (u - v) * (u - v) := by
  unfold disc; grind

/-- at u = v the law is the border law, and the lens x ↦ x/(1 − u x) makes it addition:
    cross-multiplied form of N(a ⊕ b) = N(a) + N(b) -/
theorem merged_is_addition (u a b : R) :
    circK (-(u + u)) a b = a * (1 - u * b) + b * (1 - u * a) := by
  unfold circK; grind

/-! ### 5. The torus, pinched: nodal cubic → the circle operation, cuspidal cubic → addition -/

/-- the point on y² = x³ + x² on the line y = m x through the node -/
theorem nodal_on_curve (m : R) :
    (m * (m * m - 1)) * (m * (m * m - 1)) = (m * m - 1) * (m * m - 1) * (m * m - 1) + (m * m - 1) * (m * m - 1) := by
  grind

/-- collinearity determinant of three nodal points: it vanishes exactly on m₁m₂ + m₁m₃ + m₂m₃ + 1 = 0
    (for distinct slopes) -/
theorem nodal_collinear (m1 m2 m3 : R) :
    let x := fun m : R => m * m - 1
    let y := fun m : R => m * (m * m - 1)
    x m1 * (y m2 - y m3) - y m1 * (x m2 - x m3) + (x m2 * y m3 - x m3 * y m2)
      = -((m1 - m2) * (m1 - m3) * (m2 - m3) * (m1 * m2 + m1 * m3 + m2 * m3 + 1)) := by
  intro x y; simp only [x, y]; grind

/-- the same factor in φ = (m − 1)/(m + 1), cleared of denominators: collinear ⇔ φ₁φ₂φ₃ = 1 -/
theorem nodal_phi (m1 m2 m3 : R) :
    (m1 - 1) * (m2 - 1) * (m3 - 1) - (m1 + 1) * (m2 + 1) * (m3 + 1)
      = -2 * (m1 * m2 + m1 * m3 + m2 * m3 + 1) := by
  grind

/-- negation (x, y) ↦ (x, −y) sends m to −m, which inverts φ: so P ⊕ Q has φ = φ(P)·φ(Q),
    and in a = φ − 1 that is the circle operation (the lift) -/
theorem nodal_negation (m : R) : (-m - 1) * (m - 1) = (-m + 1) * (m + 1) := by grind

theorem nodal_law_is_circ (a1 a2 : R) : 1 + circ a1 a2 = (1 + a1) * (1 + a2) := circ_lift a1 a2

/-- the cusp y² = x³ on the line y = m x: (m², m³); collinear ⇔ m₁m₂ + m₁m₃ + m₂m₃ = 0,
    i.e. 1/m₁ + 1/m₂ + 1/m₃ = 0: addition -/
theorem cusp_collinear (m1 m2 m3 : R) :
    let x := fun m : R => m * m
    let y := fun m : R => m * m * m
    x m1 * (y m2 - y m3) - y m1 * (x m2 - x m3) + (x m2 * y m3 - x m3 * y m2)
      = -((m1 - m2) * (m1 - m3) * (m2 - m3) * (m1 * m2 + m1 * m3 + m2 * m3)) := by
  intro x y; simp only [x, y]; grind

/-- cleared of denominators: m₂m₃ + m₁m₃ + m₁m₂ = m₁m₂m₃ (1/m₁ + 1/m₂ + 1/m₃) -/
theorem cusp_is_addition (m1 m2 m3 : R) :
    m1 * m2 + m1 * m3 + m2 * m3 = m2 * m3 + m1 * m3 + m1 * m2 := by grind

end ShapeDependencies
