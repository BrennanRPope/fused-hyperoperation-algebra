/-
  Fused Hyperoperation Algebra: the elliptic lead.
  Core Lean 4 only (no mathlib), no `sorry`. NOT YET COMPILED: Lean is not installed in the
  workspace where this was written. Style follows the shared library (core `Int`, `grind`).
  To compile: copy to lean-library/FusedAlgebra/Topics/Elliptic.lean, add it to build.sh and
  FusedAlgebra.lean, run ./build.sh. Every identity was also checked with sympy (../checks.py,
  ../rung_test.py).

  Curves are y² + a1·xy + a3·y = x³ + a2·x² + a4·x + a6, parameter z = −x/y.
  The landscape point (k, q) of a curve with a double point at the origin is k = −a1, q = −a2.
-/
import FusedAlgebra.Base

namespace EllipticLead
open MultAdd

/-! ### 1. The landscape plane is the cubics with a double point -/

/-- every line y = t·x through the double point meets the curve again at x = t² + a1 t − a2 -/
theorem node_param (a1 a2 t : Int) :
    let x := t * t + a1 * t - a2
    (t * x) * (t * x) + a1 * x * (t * x) = x * x * x + a2 * x * x := by
  intro x; simp only [x]; grind

/-- the walls are the node's branches: with z = −1/t (so z·t = −1),
    t²·(1 + k z + q z²) − (t² + a1 t − a2) is a multiple of (z t + 1), with k = −a1, q = −a2 -/
theorem walls_are_branches (a1 a2 t z : Int) :
    t * t * (1 + (-a1) * z + (-a2) * z * z) - (t * t + a1 * t - a2)
      = -(a1 * t) * (z * t + 1) - a2 * ((z * t) * (z * t) - 1) := by grind

/-- the two tangent slopes m1, m2 factor the slope polynomial: a1 = −(m1+m2), a2 = −m1·m2,
    so k = m1 + m2 and q = m1·m2 -/
theorem slopes_factor (m1 m2 t : Int) :
    t * t + (-(m1 + m2)) * t - (-(m1 * m2)) = (t - m1) * (t - m2) := by grind

/-- the wall polynomial in z factors through the branches: 1 + k z + q z² = (1 + m1 z)(1 + m2 z) -/
theorem wall_poly_factor (m1 m2 z : Int) :
    1 + (m1 + m2) * z + (m1 * m2) * z * z = (1 + m1 * z) * (1 + m2 * z) := by grind

/-- discriminant of the slopes equals the map's k² − 4q -/
theorem slope_discriminant (a1 a2 : Int) : a1 * a1 + 4 * a2 = (-a1) * (-a1) - 4 * (-a2) := by grind

/-- the nodal curve's formal group: w = −1/y as a function of z is rational,
    w·(1 − a1 z − a2 z²) = z³ is the curve equation w = z³ + a1 z w + a2 z² w -/
theorem nodal_w (a1 a2 z w : Int) (h : w * (1 - a1 * z - a2 * z * z) = z * z * z) :
    w = z * z * z + a1 * z * w + a2 * z * z * w := by grind

/-! ### 2. Smoothing the node: Tate's curve at q = 0 is mult-add -/

/-- with z = 1 − v (v = 1/u, Tate's parameter inverted), the law z1 + z2 − z1 z2 is v1·v2 -/
theorem tate_q0_law (v1 v2 : Int) :
    (1 - v1) + (1 - v2) - (1 - v1) * (1 - v2) = 1 - v1 * v2 := by grind

/-- and a = −z turns it into mult-add exactly -/
theorem tate_q0_is_multAdd (a b : Int) :
    -((-a) + (-b) - (-a) * (-b)) = multAdd a b := by unfold multAdd; grind

/-- weights of the invariant differential (from ../rung_test.py): d_n has weight n,
    so the coefficient a_i first acts on the cube's layer i+1, and nothing new can enter at layer 6 -/
def d1 (a1 : Int) : Int := a1
def d2 (a1 a2 : Int) : Int := a1 * a1 + a2
def d3 (a1 a2 a3 : Int) : Int := a1 * a1 * a1 + 2 * a1 * a2 + 2 * a3
def d4 (a1 a2 a3 a4 : Int) : Int := a1 * a1 * a1 * a1 + 3 * a1 * a1 * a2 + 6 * a1 * a3 + a2 * a2 + 2 * a4
def d5 (a1 a2 a3 a4 : Int) : Int :=
  a1 * a1 * a1 * a1 * a1 + 4 * a1 * a1 * a1 * a2 + 12 * a1 * a1 * a3 + 3 * a1 * a2 * a2 + 6 * a1 * a4 + 6 * a2 * a3

theorem d4_weight (u a1 a2 a3 a4 : Int) :
    d4 (u * a1) (u * u * a2) (u * u * u * a3) (u * u * u * u * a4) = u * u * u * u * d4 a1 a2 a3 a4 := by
  unfold d4; grind

theorem d5_weight (u a1 a2 a3 a4 : Int) :
    d5 (u * a1) (u * u * a2) (u * u * u * a3) (u * u * u * u * a4)
      = u * u * u * u * u * d5 a1 a2 a3 a4 := by
  unfold d5; grind

/-- mult-add's own curve (a1 = −1, a2 = a3 = a4 = 0) has d_n = (−1)^n: the differential dz/(1+z),
    the log of mult-add, through layer 6 -/
theorem multAdd_curve_d :
    d1 (-1) = -1 ∧ d2 (-1) 0 = 1 ∧ d3 (-1) 0 0 = -1 ∧ d4 (-1) 0 0 0 = 1 ∧ d5 (-1) 0 0 0 = -1 := by
  unfold d1 d2 d3 d4 d5; decide

/-! ### 3. Connections -/

/-- Ochanine's elliptic genus at ε = δ² (the nodal case): numerator a√R(b) + b√R(a) with
    √R(x) = 1 − δx² equals (a + b)(1 − δab), so the law is (a + b)/(1 + δab): the signature at δ = 1 -/
theorem ochanine_nodal_num (δ a b : Int) :
    a * (1 - δ * b * b) + b * (1 - δ * a * a) = (a + b) * (1 - δ * a * b) := by grind

theorem ochanine_nodal_den (δ a b : Int) :
    1 - δ * δ * a * a * b * b = (1 - δ * a * b) * (1 + δ * a * b) := by grind

/-- Frey curve y² = x(x − A)(x + B) = x³ + (B − A)x² − AB·x: 16 × discriminant of the cubic -/
theorem frey_discriminant (A B : Int) :
    16 * ((B - A) * (B - A) * ((-(A * B)) * (-(A * B))) - 4 * ((-(A * B)) * (-(A * B)) * (-(A * B))))
      = 16 * (A * A) * (B * B) * ((A + B) * (A + B)) := by grind

/-- c4 = b2² − 24 b4 with b2 = 4(B − A), b4 = −2AB -/
theorem frey_c4 (A B : Int) :
    (4 * (B - A)) * (4 * (B - A)) - 24 * (-2 * A * B) = 16 * (A * A + A * B + B * B) := by grind

/-- semistable: at a prime dividing A (A ≡ 0) c4/16 is B²; dividing A + B (B ≡ −A) it is A².
    Either is nonzero mod the prime when gcd(A, B) = 1, so the reduction is a node, never a cusp. -/
theorem frey_c4_at_A (B : Int) : 0 * 0 + 0 * B + B * B = B * B := by grind
theorem frey_c4_at_C (A : Int) : A * A + A * (-A) + (-A) * (-A) = A * A := by grind

/-- Euler factors are wall polynomials: 1 − a T + p T² = (1 − αT)(1 − βT) when a = α + β, p = αβ -/
theorem euler_factor_walls (α β T : Int) :
    1 - (α + β) * T + (α * β) * T * T = (1 - α * T) * (1 - β * T) := by grind

/-- Hasse (a² < 4p) puts the point (k, q) = (−a, p) in the turning country: k² − 4q < 0 -/
theorem hasse_turning (a p : Int) (h : a * a < 4 * p) : (-a) * (-a) - 4 * p < 0 := by grind

/-- units: scaling T by λ moves (k, q) to (λk, λ²q), the landscape's units action -/
theorem euler_units (k q lam T : Int) :
    1 + k * (lam * T) + q * (lam * T) * (lam * T) = 1 + (lam * k) * T + (lam * lam * q) * T * T := by grind

/-! ### 4. Verdict: the fitted curve is singular exactly at the rungs it can copy
    The cubic fitted to rung s through layer 5 has a1 = −s, a2 = −s(s−1)/4, a3 = s(s−1)/24,
    a4 = s(s−1)(11s−4)/288; its layer-6 gap is −s(s−1)(45s² − 125s + 22)/17280.
    Cleared of denominators (×17280): -/
def gap6 (s : Int) : Int := -(s * (s - 1) * (45 * s * s - 125 * s + 22))

theorem gap6_zero_at_rungs : gap6 0 = 0 ∧ gap6 1 = 0 := by unfold gap6; decide

/-- no other whole rung passes layer 6 -/
theorem gap6_whole (s : Int) (h0 : s ≠ 0) (h1' : s ≠ 1) : gap6 s ≠ 0 := by
  unfold gap6
  have hq : 45 * s * s - 125 * s + 22 ≠ 0 := by
    intro h
    -- 45s² − 125s + 22 has no integer root: its roots are 0.189 and 2.589
    rcases (show s ≤ 0 ∨ s = 1 ∨ s = 2 ∨ 3 ≤ s by omega) with hs | hs | hs | hs
    · have : 45 * s * s ≥ 0 := by
        have := Int.mul_nonneg (show (0:Int) ≤ -s by omega) (show (0:Int) ≤ -s by omega); grind
      omega
    · subst hs; omega
    · subst hs; omega
    · have : 45 * s * s ≥ 135 * s := by
        have := Int.mul_le_mul_of_nonneg_left (show (3:Int) ≤ s by omega) (show (0:Int) ≤ 45 * s by omega); grind
      omega
  intro h
  have : s * (s - 1) * (45 * s * s - 125 * s + 22) = 0 := by omega
  rcases Int.mul_eq_zero.mp this with h' | h'
  · rcases Int.mul_eq_zero.mp h' with h'' | h''
    · exact h0 h''
    · exact h1' (by omega)
  · exact hq h'

end EllipticLead
