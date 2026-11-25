/-
  Special relativity's velocity addition in Fused Hyperoperation Algebra
  ======================================================================

  Companion page and checks: ../index.html, ../checks.py (37 checks).
  Core Lean 4.37.0-rc1 only (no mathlib), no `sorry`. Imports Base and
  Topics.Thermodynamics (reuses its `vel` and `lens`). Compiled against the
  shared library in the session that wrote it.

  Speeds are in units of light (c = 1); `velC` keeps c and shows that is only
  a change of units. Velocity addition u ⊕ v = (u + v)/(1 + uv).

    1. Mult-add in disguise: with d(u) = 2u/(1 − u) (Doppler factor squared,
       minus 1), d(u ⊕ v) = d(u) ◯ d(v). The lift 1 + d is (1 + u)/(1 − u).
    2. Light is the wall: 1 ⊕ v = 1, and d sends −1 to mult-add's wall −1.
       Below light stays below light.
    3. n equal kicks: 1 + d(u ⊕ … ⊕ u) = (1 + d(u))ⁿ, the mult-add power.
    4. The geometry: a boost [[1, u], [u, 1]] multiplies as boosts add; on the
       light-cone axes it is the squeeze (1 + u, 1 − u), whose ratio is the lift;
       it keeps the hyperbola T² − X² up to the factor 1 − u².
    5. The landscape: no corner (k = 0), third layer −1/c² (q), walls ±c.
    6. The rungs in natural units: rung s = t² sits at (2t, t² − 1), the
       parabola from velocity addition (t = 0) to mult-add at k = 2 (t = 1).

  What Lean does not check: tanh, artanh and exp (rapidity). Rapidity enters
  as "the ruler that turns the lift's multiplication into addition"; the real
  identities are in ../checks.py.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.Thermodynamics

namespace Relativity
open CircleOp ThermoFused

/-! ## 1. Mult-add in disguise -/

/-- The mult-add variable of a speed: Doppler factor squared minus 1. -/
def dop (u : Rat) : Rat := 2 * u / (1 - u)

/-- Its lift is the Doppler factor squared (Bondi's k²). -/
theorem dop_lift (u : Rat) (h : 1 - u ≠ 0) : 1 + dop u = (1 + u) / (1 - u) := by
  unfold dop; grind

/-- Polynomial core: (1 + u ⊕ v)(1 − u)(1 − v) and the lift agree. -/
theorem lift_num (u v : Rat) : (1 + u * v) + (u + v) = (1 + u) * (1 + v) := by grind
theorem lift_den (u v : Rat) : (1 + u * v) - (u + v) = (1 - u) * (1 - v) := by grind

/-- **Velocity addition is mult-add in disguise**: d(u ⊕ v) = d(u) ◯ d(v). -/
theorem dop_vel (u v : Rat) (hu : 1 - u ≠ 0) (hv : 1 - v ≠ 0) (huv : 1 + u * v ≠ 0) :
    dop (vel u v) = circ (dop u) (dop v) := by
  unfold dop vel circ
  have e : 1 - (u + v) / (1 + u * v) = (1 - u) * (1 - v) / (1 + u * v) := by grind
  rw [e]
  have h1 : (1 - u) * (1 - v) ≠ 0 := by
    intro h0; rcases Rat.mul_eq_zero.mp h0 with h | h <;> contradiction
  grind

/-- The way back: u = d/(2 + d). -/
theorem dop_back (u : Rat) (h : 1 - u ≠ 0) : dop u / (2 + dop u) = u := by
  unfold dop
  have e : 2 + 2 * u / (1 - u) = 2 / (1 - u) := by grind
  rw [e]; grind

/-- c is only units: (u ⊕_c v)/c = (u/c) ⊕ (v/c). -/
def velC (c u v : Rat) : Rat := (u + v) / (1 + u * v / (c * c))

theorem velC_units (c u v : Rat) (hc : c ≠ 0) :
    velC c u v / c = vel (u / c) (v / c) := by
  unfold velC vel
  have e : u / c * (v / c) = u * v / (c * c) := by grind
  rw [e]; grind

/-! ## 2. Light is the wall -/

/-- Light absorbs, as −1 does in mult-add: 1 ⊕ v = 1. -/
theorem light_absorbs (v : Rat) (h : 1 + v ≠ 0) : vel 1 v = 1 := vel_wall_pos v h

theorem circ_wall (b : Rat) : circ (-1) b = -1 := by unfold circ; grind

/-- d carries the light wall −1 onto mult-add's wall −1. -/
theorem dop_wall : dop (-1) = -1 := by unfold dop; grind

/-- How far below light: 1 − u ⊕ v = (1 − u)(1 − v)/(1 + uv). -/
theorem gap_to_light (u v : Rat) (h : 1 + u * v ≠ 0) :
    1 - vel u v = (1 - u) * (1 - v) / (1 + u * v) := by
  unfold vel; grind

/-- Two speeds under light, both forward, add to a speed under light. -/
theorem under_light (u v : Rat) (hu0 : 0 ≤ u) (hu : u < 1) (hv0 : 0 ≤ v) (hv : v < 1) :
    vel u v < 1 := by
  have hd : 0 < 1 + u * v := by
    have := Rat.mul_nonneg hu0 hv0; grind
  have hn : 0 < (1 - u) * (1 - v) := Rat.mul_pos (by grind) (by grind)
  unfold vel
  exact (Rat.div_lt_iff hd).mpr (by grind)

/-! ## 3. n equal kicks: the mult-add power -/

/-- Kick n times with speed u. -/
def kicks : Nat → Rat → Rat
  | 0, _ => 0
  | n + 1, u => vel (kicks n u) u

mutual
/-- Numerator and denominator of n kicks. -/
def kN : Nat → Rat → Rat
  | 0, _ => 0
  | n + 1, u => kN n u + u * kD n u
def kD : Nat → Rat → Rat
  | 0, _ => 1
  | n + 1, u => kD n u + u * kN n u
end

/-- The lift of n kicks splits: D + N = (1 + u)ⁿ and D − N = (1 − u)ⁿ. -/
theorem kicks_plus (n : Nat) (u : Rat) : kD n u + kN n u = (1 + u) ^ n := by
  induction n with
  | zero => simp only [kD, kN, Rat.pow_zero]; grind
  | succ n ih => simp only [kD, kN, Rat.pow_succ]; grind

theorem kicks_minus (n : Nat) (u : Rat) : kD n u - kN n u = (1 - u) ^ n := by
  induction n with
  | zero => simp only [kD, kN, Rat.pow_zero]; grind
  | succ n ih => simp only [kD, kN, Rat.pow_succ]; grind

/-- One more kick in the numerator/denominator form is velocity addition. -/
theorem kick_step (n : Nat) (u : Rat) (hD : kD n u ≠ 0) (h2 : kD n u + u * kN n u ≠ 0) :
    vel (kN n u / kD n u) u = kN (n + 1) u / kD (n + 1) u := by
  unfold vel; simp only [kN, kD]
  have e1 : kN n u / kD n u + u = (kN n u + u * kD n u) / kD n u := by grind
  have e2 : 1 + kN n u / kD n u * u = (kD n u + u * kN n u) / kD n u := by grind
  rw [e1, e2]
  grind

/-- **n kicks lift to the mult-add power**: (D + N)/(D − N) = ((1 + u)/(1 − u))ⁿ. -/
theorem kicks_lift (n : Nat) (u : Rat) :
    (kD n u + kN n u) * (1 - u) ^ n = (kD n u - kN n u) * (1 + u) ^ n := by
  rw [kicks_plus, kicks_minus]; grind

/-! ## 4. The geometry: boosts, squeezes and the hyperbola -/

/-- Unnormalized boost [[1, u], [u, 1]] as four numbers. -/
structure B2 where
  a : Rat
  b : Rat
  c : Rat
  d : Rat

def boost (u : Rat) : B2 := ⟨1, u, u, 1⟩
def B2.mul (x y : B2) : B2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩

/-- Boosts multiply as speeds velocity-add: B(u)B(v) = (1 + uv)·B(u ⊕ v). -/
theorem boost_mul (u v : Rat) :
    B2.mul (boost u) (boost v) = ⟨1 + u * v, u + v, u + v, 1 + u * v⟩ := by
  simp only [B2.mul, boost, B2.mk.injEq]; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

theorem boost_mul_vel (u v : Rat) (h : 1 + u * v ≠ 0) :
    (1 + u * v) * vel u v = u + v := by unfold vel; grind

/-- On the light-cone axes P = T + X, M = T − X a boost is a squeeze. -/
theorem squeeze_plus (u T X : Rat) : (T + u * X) + (X + u * T) = (1 + u) * (T + X) := by grind
theorem squeeze_minus (u T X : Rat) : (T + u * X) - (X + u * T) = (1 - u) * (T - X) := by grind

/-- The squeeze ratio is the lift 1 + d. -/
theorem squeeze_ratio (u : Rat) (h : 1 - u ≠ 0) : (1 + u) / (1 - u) = 1 + dop u :=
  (dop_lift u h).symm

/-- The hyperbola: T'² − X'² = (1 − u²)(T² − X²), and T² − X² = P·M. -/
theorem hyperbola (u T X : Rat) :
    (T + u * X) * (T + u * X) - (X + u * T) * (X + u * T)
      = (1 - u * u) * (T * T - X * X) := by grind

theorem hyperbola_cone (T X : Rat) : T * T - X * X = (T + X) * (T - X) := by grind

/-- The light cone is fixed: a ray with T = X (or T = −X) stays one. -/
theorem cone_fixed (u T X : Rat) (h : T - X = 0) : (T + u * X) - (X + u * T) = 0 := by
  rw [squeeze_minus, h]; grind

theorem cone_fixed' (u T X : Rat) (h : T + X = 0) : (T + u * X) + (X + u * T) = 0 := by
  rw [squeeze_plus, h]; grind

/-! ## 5. The landscape: no corner, a bent third layer, walls at ±c -/

/-- The landscape law (a + b + k·ab)/(1 − q·ab) at k = 0, q = −1/c² is velocity addition. -/
theorem landscape_point (c a b : Rat) :
    (a + b + 0 * a * b) / (1 - (-(1 / (c * c))) * a * b) = (a + b) / (1 + a * b / (c * c)) := by
  grind

/-- Its walls 1 + k·a + q·a² = (1 − a/c)(1 + a/c): the light cone. -/
theorem landscape_walls (c a : Rat) (hc : c ≠ 0) :
    1 + 0 * a + (-(1 / (c * c))) * a * a = (1 - a / c) * (1 + a / c) := by grind

/-- The corner tax: w = u ⊕ v satisfies w = u + v − u·v·w. -/
theorem corner_tax (u v : Rat) (h : 1 + u * v ≠ 0) :
    vel u v = u + v - u * v * vel u v := by
  unfold vel; grind

/-! ## 6. The rungs in natural units -/

/-- Rung s's third-layer dot is (s, s(s − 1)/4). In units λ the dot moves to
    (λk, λ²q). With s = t² and λ = 2/t it lands at (2t, t² − 1). -/
theorem rung_natural_k (t : Rat) (ht : t ≠ 0) : (2 / t) * (t * t) = 2 * t := by grind

theorem rung_natural_q (t : Rat) (ht : t ≠ 0) :
    (2 / t) * (2 / t) * ((t * t) * (t * t - 1) / 4) = t * t - 1 := by grind

/-- That curve is the border parabola moved down by 1: q = k²/4 − 1. -/
theorem rung_parabola (t : Rat) : (2 * t) * (2 * t) / 4 - 1 = t * t - 1 := by grind

/-- Its walls factor: 1 + 2t·a + (t² − 1)a² = (1 + (1 + t)a)(1 − (1 − t)a). -/
theorem rung_walls (t a : Rat) :
    1 + 2 * t * a + (t * t - 1) * a * a = (1 + (1 + t) * a) * (1 - (1 - t) * a) := by grind

/-- t = 0: walls ±1, velocity addition. -/
theorem rung_zero_is_velocity (a : Rat) :
    1 + 2 * 0 * a + (0 * 0 - 1) * a * a = (1 + a) * (1 - a) := by grind

/-- t = 1: one wall at −1/2, the other at infinity: mult-add at weight 2. -/
theorem rung_one_is_multadd (a : Rat) :
    1 + 2 * 1 * a + (1 * 1 - 1) * a * a = 1 + 2 * a := by grind

/-- In raw units rung s = t² has walls −2/(t(1 + t)) and 2/(t(1 − t)): a speed
    limit near 2/√s as s → 0, light speed growing to infinity at plain addition. -/
theorem rung_raw_walls (t a : Rat) :
    4 * (1 + (t * t) * a + (t * t) * (t * t - 1) / 4 * a * a)
      = (2 + t * (1 + t) * a) * (2 - t * (1 - t) * a) := by grind

end Relativity
