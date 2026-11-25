/-
  General relativity through the circle operation (called mult-add before 1.0): the elementary parts
  ==========================================================

  Companion checks: ../checks.py (sympy, mpmath, scipy). Core Lean 4.37.0-rc1, no mathlib,
  no `sorry`. Imports Base only. Units G = c = 1; x = M/r.

  What is proved here (the algebra; the calculus is in the checks):
  1. The dial's boosts made local: two boosts close on −q times a rotation, preserve the form
     t² + q(x² + y²), and commute exactly at the Galilean point q = 0. The rotation that
     linear lapses produce in the hypersurface bracket (Hojman–Kuchař–Teitelboim) is the same.
  2. Route (a), the dial alone: the light-bending and perihelion fractions (1+γ)/2 and
     (2+2γ−β)/3 for every dial-only theory (γ = 0) or conformal one (γ = −1), and the
     Cassini bound against them.
  3. Route (c), the river: static clocks and radial rulers carry the same boost factor, so
     the time and space deficits of Schwarzschild are circle inverses; the horizon is the
     river at light speed; isotropic radius; β = γ = 1 as exact polynomial remainders;
     Friedmann as Newton's escape speed; Gravity Probe B's factor as the gyroscope's speed
     through the river; redshifts combine by the circle operation.
  4. Route (b): a R + b R·g is divergence-free only for b = −a/2 (given contracted Bianchi).
-/
import FusedAlgebra.Base

namespace GeneralRelativity
open CircleOp

open scoped CircleOp

/-! ## 1. The dial's boosts on (t, x, y) -/

/-- 3×3 matrices over the rationals, rows t, x, y. -/
structure M3 where
  a : Rat
  b : Rat
  c : Rat
  d : Rat
  e : Rat
  f : Rat
  g : Rat
  h : Rat
  i : Rat

namespace M3
def mul (A B : M3) : M3 :=
  ⟨A.a*B.a + A.b*B.d + A.c*B.g, A.a*B.b + A.b*B.e + A.c*B.h, A.a*B.c + A.b*B.f + A.c*B.i,
   A.d*B.a + A.e*B.d + A.f*B.g, A.d*B.b + A.e*B.e + A.f*B.h, A.d*B.c + A.e*B.f + A.f*B.i,
   A.g*B.a + A.h*B.d + A.i*B.g, A.g*B.b + A.h*B.e + A.i*B.h, A.g*B.c + A.h*B.f + A.i*B.i⟩
def sub (A B : M3) : M3 :=
  ⟨A.a-B.a, A.b-B.b, A.c-B.c, A.d-B.d, A.e-B.e, A.f-B.f, A.g-B.g, A.h-B.h, A.i-B.i⟩
def add (A B : M3) : M3 :=
  ⟨A.a+B.a, A.b+B.b, A.c+B.c, A.d+B.d, A.e+B.e, A.f+B.f, A.g+B.g, A.h+B.h, A.i+B.i⟩
def smul (s : Rat) (A : M3) : M3 := ⟨s*A.a, s*A.b, s*A.c, s*A.d, s*A.e, s*A.f, s*A.g, s*A.h, s*A.i⟩
def tr (A : M3) : M3 := ⟨A.a, A.d, A.g, A.b, A.e, A.h, A.c, A.f, A.i⟩
def comm (A B : M3) : M3 := sub (mul A B) (mul B A)
def zero : M3 := ⟨0, 0, 0, 0, 0, 0, 0, 0, 0⟩
end M3

open M3

/-- Dial boost along x: dt = −q·x, dx = t (q = −1/c² relativity, 0 Galilean, +1 rotation). -/
def Kx (q : Rat) : M3 := ⟨0, -q, 0, 1, 0, 0, 0, 0, 0⟩
/-- Dial boost along y. -/
def Ky (q : Rat) : M3 := ⟨0, 0, -q, 0, 0, 0, 1, 0, 0⟩
/-- Rotation in the x–y plane. -/
def Rxy : M3 := ⟨0, 0, 0, 0, 0, 1, 0, -1, 0⟩
/-- The dial's quadratic form t² + q(x² + y²) (for q = −1/c²: t² − r²/c²). -/
def Gq (q : Rat) : M3 := ⟨1, 0, 0, 0, q, 0, 0, 0, q⟩

/-- **Two boosts close on a rotation, with coefficient −q** (the twist, made infinitesimal). -/
theorem boosts_close (q : Rat) : comm (Kx q) (Ky q) = smul (-q) Rxy := by
  simp only [comm, mul, sub, smul, Kx, Ky, Rxy, M3.mk.injEq]; grind

/-- At the Galilean point the boosts commute: no twist. -/
theorem galilean_boosts_commute : comm (Kx 0) (Ky 0) = zero := by
  simp only [comm, mul, sub, Kx, Ky, zero, M3.mk.injEq]; grind

/-- The boost along x keeps the dial's form: Kᵀ G + G K = 0. -/
theorem boost_x_keeps_form (q : Rat) : add (mul (tr (Kx q)) (Gq q)) (mul (Gq q) (Kx q)) = zero := by
  simp only [add, mul, tr, Gq, Kx, zero, M3.mk.injEq]; grind

theorem boost_y_keeps_form (q : Rat) : add (mul (tr (Ky q)) (Gq q)) (mul (Gq q) (Ky q)) = zero := by
  simp only [add, mul, tr, Gq, Ky, zero, M3.mk.injEq]; grind

/-- The rotation keeps the form too, for every q. -/
theorem rotation_keeps_form (q : Rat) : add (mul (tr Rxy) (Gq q)) (mul (Gq q) Rxy) = zero := by
  simp only [add, mul, tr, Gq, Rxy, zero, M3.mk.injEq]; grind

/-- **The hypersurface bracket with linear lapses** N = b·x, M = c·x: the shift N∇M − M∇N
    is the rigid rotation field x ↦ (b × c) × x … written out by components. -/
theorem lapse_bracket_is_rotation (b1 b2 b3 c1 c2 c3 x y z : Rat) :
    let N := b1*x + b2*y + b3*z
    let M := c1*x + c2*y + c3*z
    let w1 := b2*c3 - b3*c2
    let w2 := b3*c1 - b1*c3
    let w3 := b1*c2 - b2*c1
    (N*c1 - M*b1 = w2*z - w3*y) ∧ (N*c2 - M*b2 = w3*x - w1*z) ∧ (N*c3 - M*b3 = w1*y - w2*x) := by
  intro N M w1 w2 w3
  refine ⟨?_, ?_, ?_⟩ <;> (simp only [N, M, w1, w2, w3]; grind)

/-! ## 2. Route (a): the dial alone as a field -/

/-- Light bending as a fraction of GR's: (1 + γ)/2. -/
def bendFrac (γ : Rat) : Rat := (1 + γ) / 2
/-- Perihelion advance as a fraction of GR's: (2 + 2γ − β)/3. -/
def periFrac (γ β : Rat) : Rat := (2 + 2*γ - β) / 3

theorem river_full_bending : bendFrac 1 = 1 := by
  simp only [bendFrac]; grind
theorem river_full_perihelion : periFrac 1 1 = 1 := by
  simp only [periFrac]; grind

/-- **Any dial-only theory (space left flat, γ = 0) bends light by exactly half, whatever β.** -/
theorem dial_only_half_bending : bendFrac 0 = 1/2 := by
  simp only [bendFrac]; grind
/-- Its perihelion is 2/3 of GR's for A = 1 − 2U (β = 0) and 1/3 for A = e^{−2U} (β = 1). -/
theorem dial_only_perihelion : periFrac 0 0 = 2/3 ∧ periFrac 0 1 = 1/3 := by
  simp only [periFrac]; grind
/-- A conformal dial (Nordström, γ = −1): no bending, perihelion backwards by a third. -/
theorem conformal_dial : bendFrac (-1) = 0 ∧ periFrac (-1) 1 = -1/3 := by
  simp only [bendFrac, periFrac]; grind

/-- Full bending forces γ = 1: space must curve as much as time. -/
theorem full_bending_iff (γ : Rat) : bendFrac γ = 1 ↔ γ = 1 := by
  unfold bendFrac; constructor <;> intro h <;> grind

/-- **Cassini** (γ − 1 = (2.1 ± 2.3)·10⁻⁶·10): γ = 0 sits (1 − 0.000021)/0.000023 ≈ 43,477 σ away. -/
theorem cassini_excludes_dial_only : 40000 * 23 < 1000000 - 21 := by decide

/-! ## 3. Route (c): the river (free fall at Newton's escape speed through flat space) -/

/-- River speed squared at x = M/r: Newton's escape speed, v² = 2x. -/
def riverV2 (x : Rat) : Rat := 2 * x

/-- Static clocks run at 1/γ of the river: 1 − v² = 1 − 2M/r = −g_tt. -/
theorem static_clock (x : Rat) : 1 - riverV2 x = 1 - 2 * x := by unfold riverV2; grind

/-- Static radial rulers stretch by γ of the river: 1/(1 − v²) = g_rr. -/
theorem static_ruler (x : Rat) : 1 / (1 - riverV2 x) = 1 / (1 - 2 * x) := by unfold riverV2; grind

/-- **The two deficits are circle inverses**: (−g_tt − 1) ◯ (g_rr − 1) = 0, i.e.
    (−2x) ◯ (2x/(1 − 2x)) = 0. One boost, the same factor for time and space. -/
theorem deficits_inverse (x : Rat) (h : 1 - 2 * x ≠ 0) : (-2 * x) ◯ (2 * x / (1 - 2 * x)) = 0 := by
  show -2 * x + 2 * x / (1 - 2 * x) + -2 * x * (2 * x / (1 - 2 * x)) = 0
  have e : 2 * x / (1 - 2 * x) * (1 - 2 * x) = 2 * x := by
    rw [Rat.div_mul_cancel h]
  grind

/-- Equivalently g_tt · g_rr = −1: (1 − 2x) · 1/(1 − 2x) = 1. -/
theorem metric_product (x : Rat) (h : 1 - 2 * x ≠ 0) : (1 - 2 * x) * (1 / (1 - 2 * x)) = 1 := by
  rw [Rat.div_def, Rat.one_mul, Rat.mul_inv_cancel _ h]

/-- **The horizon is the river at light speed** (the velocity wall): 2M/r = 1 ↔ r = 2M. -/
theorem horizon_is_wall (M r : Rat) (hr : r ≠ 0) : 2 * M / r = 1 ↔ r = 2 * M := by
  constructor
  · intro h
    have : 2 * M / r * r = 1 * r := by rw [h]
    rw [Rat.div_mul_cancel hr] at this; grind
  · intro h; subst h; rw [Rat.div_def, Rat.mul_inv_cancel _ hr]

/-- Isotropic radius: with r = ρ(1 + u)² and M = 2uρ, 1 − 2M/r = ((1 − u)/(1 + u))², written
    without division: (1 + u)²·ρ − 2M equals (1 − u)²·ρ. -/
theorem isotropic_radius (ρ u : Rat) : ρ * (1 + u)^2 - 2 * (2 * u * ρ) = ρ * (1 - u)^2 := by
  grind

/-- **β = 1**: (1 − U/2)² = (1 + U/2)²(1 − 2U + 2U²) − U³(3 + U)/2, so −g_tt = 1 − 2U + 2U² + O(U³). -/
theorem beta_one (U : Rat) :
    (1 - U/2)^2 = (1 + U/2)^2 * (1 - 2*U + 2*U^2) - U^3 * (3 + U) / 2 := by grind

/-- **γ = 1**: (1 + U/2)⁴ = 1 + 2U + U²(3/2 + U/2 + U²/16), so g_ij = (1 + 2U)δ_ij + O(U²). -/
theorem gamma_one (U : Rat) : (1 + U/2)^4 = 1 + 2*U + U^2 * (3/2 + U/2 + U^2/16) := by grind

/-- Friedmann is Newton's escape speed for a uniform ball (π kept as a symbol p):
    (H R)²/2 = ((4/3) p R³ ρ)/R ↔ H² = (8/3) p ρ, for R ≠ 0. -/
theorem friedmann_newton (H R ρ p : Rat) (hR : R ≠ 0) :
    (H * R)^2 / 2 = (4/3 * p * R^3 * ρ) * R⁻¹ ↔ H^2 = 8/3 * p * ρ := by
  have e : R^3 * R⁻¹ = R^2 := by
    have : R^3 = R^2 * R := by grind
    rw [this, Rat.mul_assoc, Rat.mul_inv_cancel R hR, Rat.mul_one]
  have hR2 : R^2 ≠ 0 := by
    intro h0; apply hR
    have : R * R = 0 := by grind
    rcases Rat.mul_eq_zero.mp this with h | h <;> exact h
  constructor
  · intro h
    have h' : (H^2 - 8/3 * p * ρ) * R^2 = 0 := by
      have : 4/3 * p * R^3 * ρ * R⁻¹ = 4/3 * p * ρ * (R^3 * R⁻¹) := by grind
      rw [this, e] at h; grind
    rcases Rat.mul_eq_zero.mp h' with h1 | h1
    · grind
    · exact absurd h1 hR2
  · intro h
    have : 4/3 * p * R^3 * ρ * R⁻¹ = 4/3 * p * ρ * (R^3 * R⁻¹) := by grind
    rw [this, e]; grind

/-- **Gravity Probe B's factor is the gyroscope's speed through the river**: composing the
    river (radial, v_r² = 2x) with the orbit (tangential, v_o² = x/(1 − 2x)) at right angles,
    1 − w² = (1 − v_r²)(1 − v_o²) = 1 − 3x. -/
theorem gpb_through_river (x : Rat) (h : 1 - 2 * x ≠ 0) :
    (1 - 2 * x) * (1 - x / (1 - 2 * x)) = 1 - 3 * x := by
  have e : x / (1 - 2 * x) * (1 - 2 * x) = x := Rat.div_mul_cancel h
  grind

/-- Perpendicular speeds compose by the circle operation on −v²: −w² = (−u²) ◯ (−v²)
    when w² = u² + v²(1 − u²). -/
theorem perpendicular_fuse (u2 v2 : Rat) : -(u2 + v2 * (1 - u2)) = (-u2) ◯ (-v2) := by
  show -(u2 + v2 * (1 - u2)) = -u2 + -v2 + -u2 * -v2; grind

/-- **Redshifts combine by the circle operation**: if 1 + z = (1 + z₁)(1 + z₂) then z = z₁ ◯ z₂ (any chain of
    gravitational, Doppler and cosmological shifts). -/
theorem redshifts_fuse (z z1 z2 : Rat) (h : 1 + z = (1 + z1) * (1 + z2)) : z = z1 ◯ z2 := by
  show z = z1 + z2 + z1 * z2; grind

/-- Light falling in from far away, seen by the river: blueshift γ² times Doppler² equals
    1/(1 + v)², i.e. γ·√((1 − v)/(1 + v)) = 1/(1 + v). Stated without division. -/
theorem infall_light (v : Rat) : (1 - v) * (1 + v)^2 = (1 - v^2) * (1 + v) := by grind

/-! ## 4. Route (b): why the left side is G = R − ½Rg -/

/-- **Lovelock's step in one line**: if ∇·Ric = ½ dR (contracted Bianchi) and the combination
    a·Ric + b·R·g is to be divergence-free where dR ≠ 0, then b = −a/2. -/
theorem einstein_combination (a b dR : Rat) (hdR : dR ≠ 0) (h : a * (dR / 2) + b * dR = 0) :
    b = -a / 2 := by
  have : (b + a / 2) * dR = 0 := by grind
  rcases Rat.mul_eq_zero.mp this with h1 | h1
  · grind
  · exact absurd h1 hdR

/-- And the Einstein combination is divergence-free: ½dR − ½dR = 0. -/
theorem einstein_divergence_free (dR : Rat) : dR / 2 + (-1/2) * dR = 0 := by grind

end GeneralRelativity
