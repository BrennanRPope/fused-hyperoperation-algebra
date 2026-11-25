/-
  Twisted boosts: two boosts in different directions, in Fused Hyperoperation Algebra
  ===================================================================================

  Companion page and checks: ../index.html, ../checks.py (35 checks, sympy + numpy,
  against real 4×4 Lorentz matrices, SO(3) rotations and the triangle areas).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports Base, Thermodynamics
  and Relativity (reuses `vel`, `dop`, `dop_vel`). Place after Relativity.
  Compiled in the shared library (clean build, 1,616 theorems); one step in sphere_no_wall fixed there. Written in the
  library's usual `ext <;> simp <;> grind` style for its first compile.

  LIBRARY NOTE: copied from twisted-boosts/lean/TwistedBoosts.lean (Nov 25) and
  compiled here for the first time. One proof step changed, no statement: in
  sphere_no_wall, the two calls `Rat.le_total 0 u.re` and `Rat.le_total 0 u.im`
  became `(Rat.le_total : (0 : Rat) ≤ u.re ∨ u.re ≤ 0)` and the same for u.im,
  because Rat.le_total takes its arguments implicitly. No sorry, standard axioms
  only, no exact twins of library theorems (half_kick_ruler is Relativity.dop_vel
  with u = v).

  Speeds are half-velocities on the disk: a velocity β in the plane becomes the
  complex number t = tanh(φ/2)·(direction), which is "half a kick" (t ⊕ t = β).
  Complex numbers are pairs of rationals. The dial q is the landscape's q, the
  number in the denominator 1 − q·ab of the one-dimensional law:
  q = −1 relativity (c = 1), q = 0 Galilean, q = +1 rotations (tangent addition).

    1. The boost of the dial: B_q(u) = [[1, u], [−q·ū, 1]]. In one dimension at
       q = −1 this is the relativity page's [[1, u], [u, 1]].
    2. The product: B(u)·B(v) = [[1 − q·u·v̄, u + v], [−q·(ū + v̄), 1 − q·ū·v]].
       The two diagonal entries are conjugates X + iY and X − iY with
       X = 1 − q·(u·v) and Y = q·(u × v): the corner ū·v has split into a dot part
       and a cross part.
    3. The twist e^{iω} = (X + iY)/(X − iY): a rational point on the unit circle,
       with tan(ω/2) = Y/X. It is 1 exactly when q·(u × v) = 0: collinear boosts,
       or the Galilean dial q = 0.
    4. Order matters only through the twist: u ⊕ v = twist · (v ⊕ u).
    5. One shape per dial, closed under products, with a determinant that
       multiplies; at q = −1 that is the two-dimensional light wall
       |1 + ū·v|² − |u + v|² = (1 − |u|²)(1 − |v|²), so below light stays below
       light in every direction. At q = +1 there is no wall.
    6. The cross part is twice the signed area of the triangle 0, u, v.

  What Lean does not check (real numbers): that the twist angle equals q times the
  exact area of the velocity triangle (Gauss–Bonnet), the Thomas precession
  2π(γ − 1) per lap, and the match with 4×4 Lorentz matrices. Those are in
  ../checks.py.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Thermodynamics
import FusedAlgebra.Topics.Relativity

namespace TwistedBoosts
open MultAdd ThermoFused Relativity

/-- Brennan's notation for mult-add in this file: a ✳ b = a + b + ab. -/
local infixl:65 " ✳ " => multAdd

theorem star_def (a b : Rat) : a ✳ b = a + b + a * b := rfl

/-! ## 0. Complex numbers as pairs of rationals -/

structure C where
  re : Rat
  im : Rat

@[ext] theorem C.ext' {x y : C} (h1 : x.re = y.re) (h2 : x.im = y.im) : x = y := by
  cases x; cases y; simp_all

def C.one : C := ⟨1, 0⟩
def C.add (x y : C) : C := ⟨x.re + y.re, x.im + y.im⟩
def C.mul (x y : C) : C := ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re⟩
def C.conj (x : C) : C := ⟨x.re, -x.im⟩
def C.smul (r : Rat) (x : C) : C := ⟨r * x.re, r * x.im⟩
def C.nsq (x : C) : Rat := x.re * x.re + x.im * x.im

/-- Dot and cross of two plane vectors u = (a, b), v = (c, d). -/
def dot (u v : C) : Rat := u.re * v.re + u.im * v.im
def cross (u v : C) : Rat := u.re * v.im - u.im * v.re

/-- **The corner splits**: ū·v = (u·v) + i(u × v). -/
theorem corner_split (u v : C) : C.mul (C.conj u) v = ⟨dot u v, cross u v⟩ := by
  ext <;> simp [C.mul, C.conj, dot, cross] <;> grind

/-- The other order is the conjugate: v̄·u = conj(ū·v). -/
theorem corner_swap (u v : C) : C.mul (C.conj v) u = C.conj (C.mul (C.conj u) v) := by
  ext <;> simp [C.mul, C.conj] <;> grind

theorem cross_anti (u v : C) : cross v u = - cross u v := by unfold cross; grind
theorem dot_symm (u v : C) : dot v u = dot u v := by unfold dot; grind

/-- The cross part is twice the signed area of the triangle 0, u, v (shoelace). -/
def area2 (p r s : C) : Rat :=
  (r.re - p.re) * (s.im - p.im) - (r.im - p.im) * (s.re - p.re)

theorem cross_is_area (u v : C) : cross u v = area2 ⟨0, 0⟩ u v := by
  unfold cross area2; grind

/-- Collinear speeds (both on the real axis) have no cross part. -/
theorem cross_collinear (a c : Rat) : cross ⟨a, 0⟩ ⟨c, 0⟩ = 0 := by unfold cross; grind

/-! ## 1. The boost of the dial -/

/-- A 2×2 matrix of complex numbers [[a, b], [c, d]]. -/
structure M where
  a : C
  b : C
  c : C
  d : C

@[ext] theorem M.ext' {x y : M} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c)
    (h4 : x.d = y.d) : x = y := by
  cases x; cases y; simp_all

def M.mul (x y : M) : M :=
  ⟨C.add (C.mul x.a y.a) (C.mul x.b y.c), C.add (C.mul x.a y.b) (C.mul x.b y.d),
   C.add (C.mul x.c y.a) (C.mul x.d y.c), C.add (C.mul x.c y.b) (C.mul x.d y.d)⟩

/-- The dial's boost B_q(u) = [[1, u], [−q·ū, 1]]. -/
def boost (q : Rat) (u : C) : M := ⟨C.one, u, C.smul (-q) (C.conj u), C.one⟩

/-- At q = −1 on the real axis it is the relativity page's [[1, u], [u, 1]]. -/
theorem boost_real (u : Rat) :
    boost (-1) ⟨u, 0⟩ = ⟨⟨1, 0⟩, ⟨u, 0⟩, ⟨u, 0⟩, ⟨1, 0⟩⟩ := by
  ext <;> simp [boost, C.one, C.smul, C.conj] <;> grind

/-- The two entries that carry the twist: X ± iY. -/
def tX (q : Rat) (u v : C) : Rat := 1 - q * dot u v
def tY (q : Rat) (u v : C) : Rat := q * cross u v

/-! ## 2. The product of two boosts -/

/-- **Two boosts**: B(u)·B(v) = [[X + iY, u + v], [−q·(ū + v̄), X − iY]]. -/
theorem boost_mul (q : Rat) (u v : C) :
    M.mul (boost q u) (boost q v)
      = ⟨⟨tX q u v, tY q u v⟩, C.add u v, C.smul (-q) (C.conj (C.add u v)),
         ⟨tX q u v, - tY q u v⟩⟩ := by
  ext <;> simp [M.mul, boost, C.add, C.mul, C.one, C.smul, C.conj, tX, tY, dot, cross] <;> grind

/-- The top-left entry is 1 − q·u·v̄ and the bottom-right is 1 − q·ū·v. -/
theorem entry_tl (q : Rat) (u v : C) :
    C.add C.one (C.smul (-q) (C.mul u (C.conj v))) = ⟨tX q u v, tY q u v⟩ := by
  ext <;> simp [C.add, C.one, C.smul, C.mul, C.conj, tX, tY, dot, cross] <;> grind

theorem entry_br (q : Rat) (u v : C) :
    C.add C.one (C.smul (-q) (C.mul (C.conj u) v)) = ⟨tX q u v, - tY q u v⟩ := by
  ext <;> simp [C.add, C.one, C.smul, C.mul, C.conj, tX, tY, dot, cross] <;> grind

/-- The diagonal entries are conjugates. -/
theorem diag_conj (q : Rat) (u v : C) :
    (M.mul (boost q u) (boost q v)).a = C.conj (M.mul (boost q u) (boost q v)).d := by
  rw [boost_mul]; ext <;> simp [C.conj] <;> grind

/-- In one dimension at q = −1 this is the relativity page's (1 + uv)·B(u ⊕ v). -/
theorem boost_mul_real (u v : Rat) :
    M.mul (boost (-1) ⟨u, 0⟩) (boost (-1) ⟨v, 0⟩)
      = ⟨⟨1 + u * v, 0⟩, ⟨u + v, 0⟩, ⟨u + v, 0⟩, ⟨1 + u * v, 0⟩⟩ := by
  rw [boost_mul]; ext <;> simp [tX, tY, dot, cross, C.add, C.smul, C.conj] <;> grind

/-! ## 3. The twist -/

/-- **No twist exactly when the cross part is switched off**: X + iY = X − iY ⇔ q·(u × v) = 0. -/
theorem twist_trivial_iff (q : Rat) (u v : C) :
    (⟨tX q u v, tY q u v⟩ : C) = ⟨tX q u v, - tY q u v⟩ ↔ q * cross u v = 0 := by
  constructor
  · intro h
    have h2 : tY q u v = - tY q u v := congrArg C.im h
    unfold tY at h2; grind
  · intro h
    ext
    · rfl
    · show tY q u v = - tY q u v
      unfold tY; grind

/-- Collinear boosts never twist. -/
theorem collinear_no_twist (q a c : Rat) : tY q ⟨a, 0⟩ ⟨c, 0⟩ = 0 := by
  unfold tY cross; grind

/-- The Galilean dial q = 0 never twists, in any directions. -/
theorem galilean_no_twist (u v : C) : tY 0 u v = 0 ∧ tX 0 u v = 1 := by
  unfold tY tX; constructor <;> grind

/-- The twist as a rational point on the unit circle: cos ω and sin ω. -/
def cosT (q : Rat) (u v : C) : Rat :=
  (tX q u v * tX q u v - tY q u v * tY q u v) / (tX q u v * tX q u v + tY q u v * tY q u v)
def sinT (q : Rat) (u v : C) : Rat :=
  2 * tX q u v * tY q u v / (tX q u v * tX q u v + tY q u v * tY q u v)

/-- It lies on the unit circle. -/
theorem twist_on_circle (q : Rat) (u v : C)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    cosT q u v * cosT q u v + sinT q u v * sinT q u v = 1 := by
  unfold cosT sinT
  generalize tX q u v = X at *
  generalize tY q u v = Y at *
  have e : (X * X - Y * Y) / (X * X + Y * Y) * ((X * X - Y * Y) / (X * X + Y * Y))
      + 2 * X * Y / (X * X + Y * Y) * (2 * X * Y / (X * X + Y * Y))
      = ((X * X - Y * Y) * (X * X - Y * Y) + (2 * X * Y) * (2 * X * Y))
        / ((X * X + Y * Y) * (X * X + Y * Y)) := by grind
  rw [e]
  have e2 : (X * X - Y * Y) * (X * X - Y * Y) + (2 * X * Y) * (2 * X * Y)
      = (X * X + Y * Y) * (X * X + Y * Y) := by grind
  rw [e2]
  have hh : (X * X + Y * Y) * (X * X + Y * Y) ≠ 0 := by
    intro h0; rcases Rat.mul_eq_zero.mp h0 with h' | h' <;> exact h h'
  grind

/-- The twist turns X − iY into X + iY: e^{iω}·(X − iY) = X + iY. -/
theorem twist_turns (q : Rat) (u v : C)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    C.mul ⟨cosT q u v, sinT q u v⟩ ⟨tX q u v, - tY q u v⟩ = ⟨tX q u v, tY q u v⟩ := by
  unfold cosT sinT
  generalize tX q u v = X at *
  generalize tY q u v = Y at *
  ext <;> simp [C.mul] <;> grind

/-- **The half-angle of the twist**: tan(ω/2) = sin ω/(1 + cos ω) = Y/X = q(u × v)/(1 − q(u·v)). -/
theorem twist_half_angle (q : Rat) (u v : C) (hX : tX q u v ≠ 0)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    sinT q u v / (1 + cosT q u v) = tY q u v / tX q u v := by
  unfold cosT sinT
  generalize tX q u v = X at *
  generalize tY q u v = Y at *
  have e1 : 1 + (X * X - Y * Y) / (X * X + Y * Y) = 2 * X * X / (X * X + Y * Y) := by grind
  rw [e1]
  have h2 : 2 * X * X ≠ 0 := by
    intro h0; apply hX
    have : X * X = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h' | h' <;> exact h'
  grind

/-! ## 4. Order matters only through the twist -/

/-- The fused sum's denominator D = 1 − q·ū·v, and the other order's is its conjugate. -/
def den (q : Rat) (u v : C) : C := ⟨tX q u v, - tY q u v⟩

theorem den_swap (q : Rat) (u v : C) : den q v u = C.conj (den q u v) := by
  ext <;> simp [den, C.conj, tX, tY, dot, cross] <;> grind

/-- u ⊕ v = (u + v)·D̄/|D|² (that is (u + v)/D), written without complex division. -/
def oplus (q : Rat) (u v : C) : C :=
  C.smul (1 / C.nsq (den q u v)) (C.mul (C.add u v) (C.conj (den q u v)))

/-- The twist D̄/D = D̄²/|D|². -/
def twist (q : Rat) (u v : C) : C :=
  C.smul (1 / C.nsq (den q u v)) (C.mul (C.conj (den q u v)) (C.conj (den q u v)))

/-- The twist is exactly the unit-circle point (cos ω, sin ω) of section 3. -/
theorem twist_eq (q : Rat) (u v : C)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    twist q u v = ⟨cosT q u v, sinT q u v⟩ := by
  unfold twist cosT sinT den
  generalize tX q u v = X at *
  generalize tY q u v = Y at *
  ext <;> simp [C.smul, C.mul, C.conj, C.nsq] <;> grind

/-- **Order swap**: u ⊕ v = twist(u, v) · (v ⊕ u). -/
theorem order_swap (q : Rat) (u v : C)
    (h : tX q u v * tX q u v + tY q u v * tY q u v ≠ 0) :
    oplus q u v = C.mul (twist q u v) (oplus q v u) := by
  have hs := den_swap q u v
  unfold oplus twist
  rw [hs]
  have hn : C.nsq (C.conj (den q u v)) = C.nsq (den q u v) := by
    simp only [C.nsq, C.conj] <;> grind
  rw [hn]
  have hd : C.conj (C.conj (den q u v)) = den q u v := by ext <;> simp [C.conj]
  rw [hd]
  unfold den
  have hadd : C.add v u = C.add u v := by ext <;> simp [C.add] <;> grind
  rw [hadd]
  generalize tX q u v = X at *
  generalize tY q u v = Y at *
  generalize C.add u v = s
  ext <;> simp [C.smul, C.mul, C.conj, C.nsq] <;> grind

/-! ## 5. One shape per dial, and the light wall in every direction -/

/-- The dial's shape [[α, β], [−q·β̄, ᾱ]]. -/
def shape (q : Rat) (α β : C) : M := ⟨α, β, C.smul (-q) (C.conj β), C.conj α⟩

/-- Boosts have the shape. -/
theorem boost_shape (q : Rat) (u : C) : boost q u = shape q C.one u := by
  ext <;> simp [boost, shape, C.one, C.conj]

/-- **The shape is closed under products**: one group for each q
    (q = −1 the Lorentz spin group SU(1,1), q = 0 the Galilean one, q = +1 SU(2)). -/
theorem shape_mul (q : Rat) (α β γ δ : C) :
    M.mul (shape q α β) (shape q γ δ)
      = shape q (C.add (C.mul α γ) (C.smul (-q) (C.mul β (C.conj δ))))
                (C.add (C.mul α δ) (C.mul β (C.conj γ))) := by
  ext <;> simp [M.mul, shape, C.add, C.mul, C.smul, C.conj] <;> grind

/-- Determinant of the shape: |α|² + q·|β|². -/
def sdet (q : Rat) (α β : C) : Rat := C.nsq α + q * C.nsq β

/-- **The determinant multiplies** (a two-square identity for each q). -/
theorem sdet_mul (q : Rat) (α β γ δ : C) :
    sdet q (C.add (C.mul α γ) (C.smul (-q) (C.mul β (C.conj δ))))
           (C.add (C.mul α δ) (C.mul β (C.conj γ)))
      = sdet q α β * sdet q γ δ := by
  simp only [sdet, C.nsq, C.add, C.mul, C.smul, C.conj] <;> grind

/-- For two boosts: |1 − q·ū·v|² + q·|u + v|² = (1 + q|u|²)(1 + q|v|²). -/
theorem boost_det (q : Rat) (u v : C) :
    C.nsq (den q u v) + q * C.nsq (C.add u v) = (1 + q * C.nsq u) * (1 + q * C.nsq v) := by
  simp only [den, tX, tY, dot, cross, C.nsq, C.add] <;> grind

/-- **The light wall in two dimensions** (q = −1):
    |1 + ū·v|² − |u + v|² = (1 − |u|²)(1 − |v|²). -/
theorem light_wall_2d (u v : C) :
    C.nsq (den (-1) u v) - C.nsq (C.add u v) = (1 - C.nsq u) * (1 - C.nsq v) := by
  have := boost_det (-1) u v; grind

/-- **Below light stays below light, in any directions**: |u|, |v| < 1 ⇒ |u + v| < |1 + ū·v|,
    so |u ⊕ v| < 1. -/
theorem under_light_2d (u v : C) (hu : C.nsq u < 1) (hv : C.nsq v < 1) :
    C.nsq (C.add u v) < C.nsq (den (-1) u v) := by
  have hp : 0 < (1 - C.nsq u) * (1 - C.nsq v) := Rat.mul_pos (by grind) (by grind)
  have := light_wall_2d u v; grind

/-- At the sphere end q = +1 there is no wall: the determinant of a boost is ≥ 1. -/
theorem sphere_no_wall (u : C) : 1 ≤ sdet 1 C.one u := by
  have h1 : 0 ≤ u.re * u.re := by
    rcases (Rat.le_total : (0 : Rat) ≤ u.re ∨ u.re ≤ 0) with h | h
    · exact Rat.mul_nonneg h h
    · have := Rat.mul_nonneg (by grind : 0 ≤ -u.re) (by grind : 0 ≤ -u.re); grind
  have h2 : 0 ≤ u.im * u.im := by
    rcases (Rat.le_total : (0 : Rat) ≤ u.im ∨ u.im ≤ 0) with h | h
    · exact Rat.mul_nonneg h h
    · have := Rat.mul_nonneg (by grind : 0 ≤ -u.im) (by grind : 0 ≤ -u.im); grind
  simp only [sdet, C.nsq, C.one] <;> grind

/-! ## 6. Half a kick, and the tie to mult-add -/

/-- The half-velocity t composed with itself gives the velocity: t ⊕ t = 2t/(1 + t²). -/
theorem half_kick (t : Rat) : vel t t = 2 * t / (1 + t * t) := by
  unfold vel; grind

/-- On the mult-add ruler, the velocity is the half-velocity mult-added to itself:
    d(t ⊕ t) = d(t) ✳ d(t). The half-velocity is half a kick. -/
theorem half_kick_ruler (t : Rat) (h : 1 - t ≠ 0) (h2 : 1 + t * t ≠ 0) :
    dop (vel t t) = dop t ✳ dop t :=
  dop_vel t t h h h2

/-- The twist's tangent, written with the mult-add corner: Y/X = q(u × v)/(1 − q(u·v)),
    and both parts come from the one corner ū·v. -/
theorem twist_from_corner (q : Rat) (u v : C) :
    tX q u v = 1 - q * (C.mul (C.conj u) v).re ∧ tY q u v = q * (C.mul (C.conj u) v).im := by
  rw [corner_split]; exact ⟨rfl, rfl⟩

end TwistedBoosts
