/-
  Fused Hyperoperation Algebra: differential geometry through the circle operation (called mult-add before 1.0).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base.

  Intended library path: FusedAlgebra/Topics/DiffGeometry.lean (imports only FusedAlgebra.Base).

  What is here (the algebra under the geometry; the curvature and length statements that
  need real numbers are checked in checks.py):
    1. The three conics t² + Kx² = 1: rotations (c, s) compose and keep the form, and the
       central chart r = x/t composes by the landscape law F(0, K) = (r₁ + r₂)/(1 − K r₁ r₂).
       K > 0 sphere, K = 0 flat (chart values add), K < 0 hyperbolic.
    2. The two kinds of chart differ by doubling: 2a/(1 − qa²) = F(a, a) is an endomorphism of
       F(0, q). At q = −1 it takes the Poincaré radius to the Klein radius, at q = 1 the
       stereographic radius to the gnomonic one.
    3. Curvature from a landscape point: 4q − p² is unchanged by every chart change of the line,
       so K_proj = (4q − p²)/4 and K_conf = 4q − p² are well defined. The circle operation gives −1/4 and −1,
       the dial gives −k²/4 and −k², velocity addition gives −1/c² in the Klein chart.
    4. The walls are the absolute: the cross-ratio of (0, x) with the walls is the lens
       (1 − vx)/(1 − ux), and it multiplies under the law.
    5. The circle operation in the upper half plane: height y = 1 + x, translations are dilations of the
       half plane (isometries of (dx² + dy²)/y²), and the wall x = −1 is the horizon y = 0.
    6. Einstein and Möbius additions agree on a line: both are (a + b)/(1 + ab).
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/DiffGeometry.lean; one exact twin removed (see the LIBRARY comment).
-/
import FusedAlgebra.Base

namespace DiffGeometry
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. The three conics and the central chart -/

/-- A turn of the conic t² + Kx² = 1, as the pair (c, s) = (t, x) of the image of (1, 0). -/
def rmul (K : Int) (u v : Int × Int) : Int × Int := (u.1 * v.1 - K * u.2 * v.2, u.2 * v.1 + u.1 * v.2)

/-- The form t² + Kx². -/
def form (K : Int) (u : Int × Int) : Int := u.1 * u.1 + K * u.2 * u.2

/-- Turns keep the conic: form(uv) = form(u)·form(v). -/
theorem rmul_form (K : Int) (u v : Int × Int) : form K (rmul K u v) = form K u * form K v := by
  simp only [form, rmul]; grind

theorem rmul_comm (K : Int) (u v : Int × Int) : rmul K u v = rmul K v u := by
  simp only [rmul]; ext <;> simp <;> grind

theorem rmul_assoc (K : Int) (u v w : Int × Int) : rmul K (rmul K u v) w = rmul K u (rmul K v w) := by
  simp only [rmul]; ext <;> simp <;> grind

theorem rmul_one (K : Int) (u : Int × Int) : rmul K u (1, 0) = u := by
  simp only [rmul]; ext <;> simp

/-- Central chart, cleared of denominators: if s₁ = r₁c₁ and s₂ = r₂c₂ then the turned point
    has x = c₁c₂(r₁ + r₂) and t = c₁c₂(1 − K r₁ r₂), so x/t = F(0, K)(r₁, r₂). -/
theorem central_chart (K c1 c2 r1 r2 : Int) :
    rmul K (c1, r1 * c1) (c2, r2 * c2) = (c1 * c2 * (1 - K * r1 * r2), c1 * c2 * (r1 + r2)) := by
  simp only [rmul]; ext <;> simp <;> grind

/-- The same with the chart point written as (1, r): (1, r₁)(1, r₂) = (1 − K r₁r₂, r₁ + r₂). -/
theorem central_chart_unit (K r1 r2 : Int) : rmul K (1, r1) (1, r2) = (1 - K * r1 * r2, r1 + r2) := by
  simp only [rmul]; ext <;> simp <;> grind

/-- Flat case K = 0: the line t = 1 is carried to itself and chart values add (Galilean shear). -/
theorem flat_adds (s1 s2 : Int) : rmul 0 (1, s1) (1, s2) = (1, s1 + s2) := by
  simp only [rmul]; ext <;> simp <;> grind

/-- Sphere K = 1: the unit circle is closed under turns (complex multiplication). -/
theorem sphere_closed (c1 s1 c2 s2 : Int) (h1 : c1 * c1 + s1 * s1 = 1) (h2 : c2 * c2 + s2 * s2 = 1) :
    form 1 (rmul 1 (c1, s1) (c2, s2)) = 1 := by
  rw [rmul_form]; simp only [form]; grind

/-- Hyperbola K = −1: t² − x² = 1 is closed under boosts. -/
theorem hyperbola_closed (c1 s1 c2 s2 : Int) (h1 : c1 * c1 - s1 * s1 = 1) (h2 : c2 * c2 - s2 * s2 = 1) :
    form (-1) (rmul (-1) (c1, s1) (c2, s2)) = 1 := by
  rw [rmul_form]; simp only [form]; grind

/-! ## 2. Doubling: from angle-true charts to straight-line charts -/

/-- The doubled point a ⊕ a = 2a/(1 − qa²), numerator and denominator. -/
def dN (a : Int) : Int := 2 * a
def dD (q a : Int) : Int := 1 - q * a * a

/-- Doubling is the law applied to (a, a). -/
theorem dbl_is_self (q a : Int) : dN a = a + a ∧ dD q a = 1 - q * a * a := ⟨by simp only [dN]; grind, rfl⟩

/-- The doubled numerators combine to the doubled sum: 2a(1 − qb²) + 2b(1 − qa²) = 2(a + b)(1 − qab). -/
theorem dbl_num (q a b : Int) :
    dN a * dD q b + dN b * dD q a = dN (a + b) * (1 - q * a * b) := by
  simp only [dN, dD]; grind

/-- ... and the denominators: (1 − qa²)(1 − qb²) − q(2a)(2b) = (1 − qab)² − q(a + b)². -/
theorem dbl_den (q a b : Int) :
    dD q a * dD q b - q * dN a * dN b = (1 - q * a * b) * (1 - q * a * b) - q * (a + b) * (a + b) := by
  simp only [dN, dD]; grind

/-- Together, dbl_num and dbl_den say doubling respects the law: dbl(F(a,b)) has numerator
    2(a + b)(1 − qab) and denominator (1 − qab)² − q(a + b)², the same as F(dbl a, dbl b) after
    clearing (1 − qa²)(1 − qb²). Here is the doubled sum itself, cleared by (1 − qab)². -/
theorem dbl_of_sum (q a b : Int) :
    dN (a + b) * (1 - q * a * b) = 2 * (a + b) * (1 - q * a * b) ∧
    dD q (a + b) = 1 - q * (a + b) * (a + b) := by
  simp only [dN, dD]; constructor <;> grind

/-- Poincaré to Klein (q = −1): r_K = 2r_P/(1 + r_P²) is the doubling of the velocity law. -/
theorem poincare_klein (a : Int) : dD (-1) a = 1 + a * a := by simp only [dD]; grind

/-- Stereographic to gnomonic (q = 1): tan θ = 2 tan(θ/2)/(1 − tan²(θ/2)). -/
theorem stereo_gnomonic (a : Int) : dD 1 a = 1 - a * a := by simp only [dD]; grind

/-! ## 3. Curvature read off a landscape point -/

/-- Four times the projective curvature: 4q − p² (minus the discriminant). -/
def curv4 (p q : Int) : Int := 4 * q - p * p

/-- Chart changes x ↦ x/(1 − cx) move (p, q) to (p + 2c, q + pc + c²) and keep 4q − p². -/
theorem curv_chart_free (p q c : Int) : curv4 (p + 2 * c) (q + p * c + c * c) = curv4 p q := by
  simp only [curv4]; grind

/-- Scaling the line by λ scales the curvature by λ² (so the radius by 1/λ). -/
theorem curv_scale (p q l : Int) : curv4 (l * p) (l * l * q) = l * l * curv4 p q := by
  simp only [curv4]; grind

/-- The circle operation: conformal curvature −1 (the half plane), projective curvature −1/4. -/
theorem circ_curv : curv4 1 0 = -1 := by decide

/-- The dial line: curvature −k², flat only at k = 0. -/
theorem dial_curv (k : Int) : curv4 k 0 = -(k * k) := by simp only [curv4]; grind

theorem dial_flat_iff (k : Int) : curv4 k 0 = 0 ↔ k = 0 := by
  simp only [curv4]; constructor
  · intro h; have : k * k = 0 := by grind
    rcases Int.mul_eq_zero.mp this with h | h <;> exact h
  · intro h; subst h; decide

/-- Tangent addition (0, 1) is a sphere, velocity addition (0, −1) a hyperbolic plane. -/
theorem tan_vel_curv : curv4 0 1 = 4 ∧ curv4 0 (-1) = -4 := by decide

/-- Parity (−2, 0) has the same curvature as velocity addition. -/
theorem parity_curv : curv4 (-2) 0 = curv4 0 (-1) := by decide

/-- With walls 1/u and 1/v (so p = −(u+v), q = uv), 4q − p² = −(u − v)²: real walls mean
    negative curvature. -/
theorem walls_curv (u v : Int) : curv4 (-(u + v)) (u * v) = -((u - v) * (u - v)) := by
  simp only [curv4]; grind

/-! ## 4. The walls are the absolute: cross-ratio -/

/-- Cross-ratio (a, b; P, Q). -/
def CR (a b P Q : Rat) : Rat := ((a - P) * (b - Q)) / ((a - Q) * (b - P))

/-- With 0 as base point: CR(0, x; P, Q) = P(x − Q)/(Q(x − P)). -/
theorem cr_base (x P Q : Rat) (hQ : Q ≠ 0) (hP : x - P ≠ 0) :
    CR 0 x P Q = (P * (x - Q)) / (Q * (x - P)) := by
  unfold CR
  have e1 : (0 - P) * (x - Q) = -(P * (x - Q)) := by grind
  have e2 : (0 - Q) * (x - P) = -(Q * (x - P)) := by grind
  rw [e1, e2]
  have h : Q * (x - P) ≠ 0 := by
    intro h0; rcases (Rat.mul_eq_zero.mp h0) with h1 | h1 <;> simp_all
  grind

/-- With walls P = 1/u and Q = 1/v the two halves are the lens factors:
    P(x − Q) = −PQ(1 − vx) ... -/
theorem cr_num (P Q v x : Int) (hv : v * Q = 1) : P * (x - Q) = -(P * Q) * (1 - v * x) := by
  have : P * Q * v * x = P * x := by
    calc P * Q * v * x = P * x * (v * Q) := by grind
      _ = P * x := by rw [hv]; grind
  grind

/-- ... and Q(x − P) = −PQ(1 − ux), so CR(0, x; 1/u, 1/v) = (1 − vx)/(1 − ux), the lens. -/
theorem cr_den (P Q u x : Int) (hu : u * P = 1) : Q * (x - P) = -(P * Q) * (1 - u * x) := by
  have : P * Q * u * x = Q * x := by
    calc P * Q * u * x = Q * x * (u * P) := by grind
      _ = Q * x := by rw [hu]; grind
  grind

/-- The lens multiplies under the law (numerator form): with p = −(u+v), q = uv,
    1 − qab − v(a + b + pab) = (1 − va)(1 − vb). -/
theorem lens_factor (u v a b : Int) :
    1 - (u * v) * a * b - v * (a + b + (-(u + v)) * a * b) = (1 - v * a) * (1 - v * b) := by grind

/-! ## 5. The circle operation in the upper half plane -/

-- LIBRARY: halfplane_dilation (with height y = 1 + x, translation by a is dilation by 1 + a) was an exact twin of EllipticS3.circ_is_1x1 and is removed.

/-- Dilations keep the half-plane metric (dX² + dY²)/Y² (cleared: λ²(dx² + dy²)·y² = (dx² + dy²)·(λy)²). -/
theorem halfplane_isometry (l dx dy y : Int) :
    (l * dx) * (l * dx) * (y * y) + (l * dy) * (l * dy) * (y * y) =
      (dx * dx + dy * dy) * ((l * y) * (l * y)) := by grind

/-- The wall x = −1 is the horizon y = 0, and dilations never move it. -/
theorem halfplane_horizon (a : Int) : (1 + a) * (1 + (-1)) = 0 := by grind

/-! ## 6. Einstein and Möbius addition on a line -/

/-- Möbius (Poincaré) addition of collinear a, b: ((1 + 2ab + b²)a + (1 − a²)b)/(1 + 2ab + a²b²)
    equals (a + b)/(1 + ab), cleared of denominators. -/
theorem mobius_collinear (a b : Int) :
    ((1 + 2 * a * b + b * b) * a + (1 - a * a) * b) * (1 + a * b) =
      (a + b) * (1 + 2 * a * b + a * a * (b * b)) := by grind

/-- Its denominator is the square (1 + ab)². -/
theorem mobius_den (a b : Int) : 1 + 2 * a * b + a * a * (b * b) = (1 + a * b) * (1 + a * b) := by grind

end DiffGeometry
