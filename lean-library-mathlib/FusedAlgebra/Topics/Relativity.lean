module

public import FusedAlgebra.Topics.Thermodynamics

/-!
# Special relativity's velocity addition as mult-add

Speeds are in units of light (`c = 1`). Velocity addition is
`ThermoFused.vel u v = (u + v) / (1 + u * v)`.

## Main definitions

* `Relativity.dop`: the mult-add variable of a speed, `dop u = 2 * u / (1 - u)` (the Doppler
  factor squared, minus `1`).
* `Relativity.velC`: velocity addition with the speed of light `c` kept.
* `Relativity.boost`: the unnormalised boost matrix `!![1, u; u, 1]`.
* `Relativity.kN`, `Relativity.kD`: numerator and denominator of `n` equal kicks, read off
  `boost u ^ n`.

## Main results

* `dop_vel`: velocity addition is mult-add in disguise, `dop (vel u v) = dop u ✳ dop v`.
* `under_light`: two forward speeds under light add to a speed under light.
* `kD_add_kN`, `kD_sub_kN`, `kicks_lift`: `n` equal kicks lift to the `n`-th power.
* `boost_mul`, `boost_mul_boost`: boosts multiply as speeds velocity-add; `det_boost` and
  `hyperbola` are the invariant `T ^ 2 - X ^ 2`.
* `rung_walls`: in natural units rung `s = t ^ 2` has walls `1 + (1 + t) * a` and
  `1 - (1 - t) * a`.

This is the standard algebra of the Lorentz group in one space dimension (Doppler `k`-calculus,
rapidity); what the file records is how it sits inside mult-add. `tanh`, `artanh` and `exp` are
not used. The core-Lean original worked over `ℚ`; here the identities hold in any commutative
ring or field, and `under_light` in any linearly ordered field.
-/

@[expose] public section

namespace Relativity

open MultAdd ThermoFused Matrix

variable {R : Type*} [CommRing R] {K : Type*} [Field K]

/-! ### Mult-add in disguise -/

/-- The mult-add variable of a speed: the Doppler factor squared, minus `1`. -/
def dop (u : K) : K := 2 * u / (1 - u)

/-- The lift of `dop` is the Doppler factor squared (Bondi's `k ^ 2`). -/
theorem one_add_dop (u : K) (h : 1 - u ≠ 0) : 1 + dop u = (1 + u) / (1 - u) := by
  unfold dop; field_simp; ring

/-- How far below light: `1 - vel u v = (1 - u) * (1 - v) / (1 + u * v)`. -/
theorem gap_to_light (u v : K) (h : 1 + u * v ≠ 0) :
    1 - vel u v = (1 - u) * (1 - v) / (1 + u * v) := by
  unfold vel; field_simp; ring

/-- **Velocity addition is mult-add in disguise**: `dop (vel u v) = dop u ✳ dop v`. -/
theorem dop_vel (u v : K) (hu : 1 - u ≠ 0) (hv : 1 - v ≠ 0) (huv : 1 + u * v ≠ 0) :
    dop (vel u v) = dop u ✳ dop v := by
  rw [dop, gap_to_light u v huv]
  unfold dop vel multAdd
  field_simp; ring

/-- The way back: `u = d / (2 + d)`. -/
theorem dop_back [NeZero (2 : K)] (u : K) (h : 1 - u ≠ 0) : dop u / (2 + dop u) = u := by
  have h2 : (2 : K) ≠ 0 := two_ne_zero
  have e : 2 + 2 * u / (1 - u) = 2 / (1 - u) := by field_simp; ring
  rw [dop, e]
  field_simp

/-- Velocity addition with the speed of light `c` kept. -/
def velC (c u v : K) : K := (u + v) / (1 + u * v / (c * c))

/-- `c` is only units: `velC c u v / c = vel (u / c) (v / c)`. -/
theorem velC_units (c u v : K) : velC c u v / c = vel (u / c) (v / c) := by
  have e : u / c * (v / c) = u * v / (c * c) := by ring
  rw [velC, vel, e]; ring

/-! ### Light is the wall -/

/-- `dop` carries the light wall `-1` onto mult-add's wall `-1`. -/
theorem dop_wall [NeZero (2 : K)] : dop (-1 : K) = -1 := by
  have e : (1 : K) - -1 = 2 := by ring
  rw [dop, e, mul_neg_one, neg_div, div_self two_ne_zero]

/-- Two speeds under light, both forward, add to a speed under light. -/
theorem under_light [LinearOrder K] [IsStrictOrderedRing K] (u v : K) (hu0 : 0 ≤ u) (hu : u < 1)
    (hv0 : 0 ≤ v) (hv : v < 1) : vel u v < 1 := by
  have hd : 0 < 1 + u * v := by positivity
  rw [vel, div_lt_one hd]
  nlinarith [mul_pos (sub_pos.2 hu) (sub_pos.2 hv)]

/-! ### The boost matrix -/

/-- The unnormalised boost `!![1, u; u, 1]`. -/
def boost (u : R) : Matrix (Fin 2) (Fin 2) R := !![1, u; u, 1]

/-- Boosts multiply as speeds velocity-add: `B u * B v = (1 + u * v) • B (vel u v)`, written
without division. -/
theorem boost_mul (u v : R) :
    boost u * boost v = !![1 + u * v, u + v; u + v, 1 + u * v] := by
  rw [boost, boost, mul_fin_two]
  simp only [one_mul, mul_one]
  rw [add_comm v u, add_comm (u * v) 1]

/-- The sum of two speeds is `1 + u * v` times their velocity sum. -/
theorem boost_mul_vel (u v : K) (h : 1 + u * v ≠ 0) : (1 + u * v) * vel u v = u + v := by
  unfold vel; field_simp

/-- Boosts multiply as speeds velocity-add, in matrix form. -/
theorem boost_mul_boost (u v : K) (h : 1 + u * v ≠ 0) :
    boost u * boost v = (1 + u * v) • boost (vel u v) := by
  rw [boost_mul]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [boost, boost_mul_vel u v h]

/-- The determinant of a boost is `1 - u ^ 2`, the factor by which it scales `T ^ 2 - X ^ 2`. -/
theorem det_boost (u : R) : (boost u).det = 1 - u ^ 2 := by
  rw [boost, det_fin_two_of]; ring

/-! ### `n` equal kicks: the mult-add power -/

/-- Kick `n` times with speed `u`. -/
def kicks (n : ℕ) (u : K) : K := (vel · u)^[n] 0

@[simp]
theorem kicks_zero (u : K) : kicks 0 u = 0 := rfl

theorem kicks_succ (n : ℕ) (u : K) : kicks (n + 1) u = vel (kicks n u) u :=
  Function.iterate_succ_apply' _ _ _

/-- Numerator of `n` kicks: the off-diagonal entry of `boost u ^ n`. -/
def kN (n : ℕ) (u : R) : R := (boost u ^ n) 0 1

/-- Denominator of `n` kicks: the diagonal entry of `boost u ^ n`. -/
def kD (n : ℕ) (u : R) : R := (boost u ^ n) 0 0

@[simp]
theorem kN_zero (u : R) : kN 0 u = 0 := by simp [kN]

@[simp]
theorem kD_zero (u : R) : kD 0 u = 1 := by simp [kD]

theorem kN_succ (n : ℕ) (u : R) : kN (n + 1) u = kN n u + u * kD n u := by
  rw [kN, pow_succ, Matrix.mul_apply, Fin.sum_univ_two]
  simp only [boost, kN, kD, of_apply, cons_val', cons_val_zero, cons_val_one, cons_val_fin_one]
  ring

theorem kD_succ (n : ℕ) (u : R) : kD (n + 1) u = kD n u + u * kN n u := by
  rw [kD, pow_succ, Matrix.mul_apply, Fin.sum_univ_two]
  simp only [boost, kN, kD, of_apply, cons_val', cons_val_zero, cons_val_one, cons_val_fin_one]
  ring

/-- The lift of `n` kicks splits: `D + N = (1 + u) ^ n`. -/
theorem kD_add_kN (n : ℕ) (u : R) : kD n u + kN n u = (1 + u) ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [kD_succ, kN_succ, pow_succ, ← ih]; ring

/-- The lift of `n` kicks splits: `D - N = (1 - u) ^ n`. -/
theorem kD_sub_kN (n : ℕ) (u : R) : kD n u - kN n u = (1 - u) ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [kD_succ, kN_succ, pow_succ, ← ih]; ring

/-- One more kick in the numerator/denominator form is velocity addition. -/
theorem kick_step (n : ℕ) (u : K) (hD : kD n u ≠ 0) :
    vel (kN n u / kD n u) u = kN (n + 1) u / kD (n + 1) u := by
  rw [kN_succ, kD_succ, vel]
  field_simp

/-- **`n` kicks lift to the mult-add power**: `(D + N) / (D - N) = ((1 + u) / (1 - u)) ^ n`,
written without division. -/
theorem kicks_lift (n : ℕ) (u : R) :
    (kD n u + kN n u) * (1 - u) ^ n = (kD n u - kN n u) * (1 + u) ^ n := by
  rw [kD_add_kN, kD_sub_kN, mul_comm]

/-! ### The geometry: squeezes and the hyperbola -/

/-- On the light-cone axis `P = T + X` a boost is a squeeze by `1 + u`. -/
theorem squeeze_plus (u T X : R) : (T + u * X) + (X + u * T) = (1 + u) * (T + X) := by ring

/-- On the light-cone axis `M = T - X` a boost is a squeeze by `1 - u`. -/
theorem squeeze_minus (u T X : R) : (T + u * X) - (X + u * T) = (1 - u) * (T - X) := by ring

/-- The squeeze ratio is the lift `1 + dop u`. -/
theorem squeeze_ratio (u : K) (h : 1 - u ≠ 0) : (1 + u) / (1 - u) = 1 + dop u :=
  (one_add_dop u h).symm

/-- The hyperbola: `T' ^ 2 - X' ^ 2 = (1 - u ^ 2) * (T ^ 2 - X ^ 2)`. -/
theorem hyperbola (u T X : R) :
    (T + u * X) * (T + u * X) - (X + u * T) * (X + u * T) = (1 - u * u) * (T * T - X * X) := by
  ring

/-- On the light-cone axes, `T ^ 2 - X ^ 2 = P * M`. -/
theorem hyperbola_cone (T X : R) : T * T - X * X = (T + X) * (T - X) := by ring

/-- The light cone is fixed: a ray with `T = X` stays one. -/
theorem cone_fixed (u T X : R) (h : T - X = 0) : (T + u * X) - (X + u * T) = 0 := by
  rw [squeeze_minus, h, mul_zero]

/-- The light cone is fixed: a ray with `T = -X` stays one. -/
theorem cone_fixed_neg (u T X : R) (h : T + X = 0) : (T + u * X) + (X + u * T) = 0 := by
  rw [squeeze_plus, h, mul_zero]

/-! ### The landscape: no corner, a bent third layer, walls at `±c` -/

/-- The landscape law `(a + b + k * a * b) / (1 - q * a * b)` at `k = 0`, `q = -1 / c ^ 2` is
velocity addition. -/
theorem landscape_point (c a b : K) :
    (a + b + 0 * a * b) / (1 - (-(1 / (c * c))) * a * b) = (a + b) / (1 + a * b / (c * c)) := by
  congr 1 <;> ring

/-- Its walls `1 + k * a + q * a ^ 2 = (1 - a / c) * (1 + a / c)`: the light cone. -/
theorem landscape_walls (c a : K) :
    1 + 0 * a + (-(1 / (c * c))) * a * a = (1 - a / c) * (1 + a / c) := by ring

/-- The corner tax: `w = vel u v` satisfies `w = u + v - u * v * w`. -/
theorem corner_tax (u v : K) (h : 1 + u * v ≠ 0) : vel u v = u + v - u * v * vel u v := by
  unfold vel; field_simp; ring

/-! ### The rungs in natural units -/

/-- Rung `s`'s third-layer dot is `(s, s * (s - 1) / 4)`. In units `λ` the dot moves to
`(λ * k, λ ^ 2 * q)`. With `s = t ^ 2` and `λ = 2 / t` it lands at `(2 * t, t ^ 2 - 1)`: this is
the first coordinate. -/
theorem rung_natural_k (t : K) (ht : t ≠ 0) : (2 / t) * (t * t) = 2 * t := by
  field_simp

/-- The second coordinate of the rung's dot in natural units. -/
theorem rung_natural_q [NeZero (2 : K)] (t : K) (ht : t ≠ 0) :
    (2 / t) * (2 / t) * ((t * t) * (t * t - 1) / 4) = t * t - 1 := by
  have h4 : (4 : K) ≠ 0 := by
    rw [show (4 : K) = 2 * 2 by norm_num]; exact mul_ne_zero two_ne_zero two_ne_zero
  field_simp; ring

/-- That curve is the border parabola moved down by `1`: `q = k ^ 2 / 4 - 1`. -/
theorem rung_parabola [NeZero (2 : K)] (t : K) : (2 * t) * (2 * t) / 4 - 1 = t * t - 1 := by
  have h4 : (4 : K) ≠ 0 := by
    rw [show (4 : K) = 2 * 2 by norm_num]; exact mul_ne_zero two_ne_zero two_ne_zero
  field_simp; ring

/-- Its walls factor into `1 + (1 + t) * a` and `1 - (1 - t) * a`. -/
theorem rung_walls (t a : R) :
    1 + 2 * t * a + (t * t - 1) * a * a = (1 + (1 + t) * a) * (1 - (1 - t) * a) := by ring

/-- `t = 0`: walls `±1`, velocity addition. -/
theorem rung_zero_is_velocity (a : R) :
    1 + 2 * 0 * a + (0 * 0 - 1) * a * a = (1 + a) * (1 - a) := by ring

/-- `t = 1`: one wall at `-1 / 2`, the other at infinity: mult-add at weight `2`. -/
theorem rung_one_is_multAdd (a : R) : 1 + 2 * 1 * a + (1 * 1 - 1) * a * a = 1 + 2 * a := by ring

/-- In raw units rung `s = t ^ 2` has walls `-2 / (t * (1 + t))` and `2 / (t * (1 - t))`: a speed
limit near `2 / √s` as `s → 0`, light speed growing to infinity at plain addition. -/
theorem rung_raw_walls [NeZero (2 : K)] (t a : K) :
    4 * (1 + (t * t) * a + (t * t) * (t * t - 1) / 4 * a * a)
      = (2 + t * (1 + t) * a) * (2 - t * (1 - t) * a) := by
  have h4 : (4 : K) ≠ 0 := by
    rw [show (4 : K) = 2 * 2 by norm_num]; exact mul_ne_zero two_ne_zero two_ne_zero
  field_simp; ring

end Relativity
