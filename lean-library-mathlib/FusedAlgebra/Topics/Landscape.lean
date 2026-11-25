module

public import FusedAlgebra.Basic

/-!
# The landscape map

A point `(p, q)` of the landscape plane is the law

`F_{p,q}(a, b) = (a + b + p * a * b) / (1 - q * a * b)`,

where `p` is the weight of the corner cell `a * b` (the dial `k` of `MultAdd.wAdd`) and `q` is
the weight of the next layer (the bend). The law is kept as a numerator `num p a b` and a
denominator `den q a b`, and every statement is a denominator-free polynomial identity, so it
holds in any commutative ring (the core-Lean original worked over `ℤ`).

## Main results

* `assoc_num`, `assoc_den`: every point of the plane is associative (projective form).
* `riccati_rate`: the growth law behind the point `(p, q)` is `y' = 1 + p * y + q * y ^ 2`.
* `lens_factor_u`, `lens_factor_v`, `lens_lift`: with `p = -(u + v)` and `q = u * v` the lens
  `x ↦ (1 - u * x) / (1 - v * x)` turns the law into multiplication; the walls are `1 / u`, `1 / v`.
* `multAdd_point`, `velocity_point`: mult-add and relativistic velocity addition as points.
* `units_num`, `units_den`, `dial_is_multAdd`: a change of units moves `(p, q)` to
  `(l * p, l ^ 2 * q)`; every `k ≠ 0` on the straight line `q = 0` is mult-add in other units.
* `channel_floor`: the floor of the straight line for a symmetric noisy channel.

All of this is known mathematics: these are the one-dimensional formal group laws with a rational
addition formula, and the lens is their logarithm. Mathlib has no ready-made structure for them
(the numerator is `MultAdd.wAdd`, the weighted circle operation), so the statements stay as
polynomial identities, proved by `ring`.
-/

@[expose] public section

namespace LandscapeMap

open MultAdd

variable {R : Type*} [CommRing R]

/-- The numerator of `F_{p,q}(a, b)`: the weighted mult-add `a + b + p * a * b`. -/
def num (p a b : R) : R := wAdd p a b

/-- The denominator of `F_{p,q}(a, b)`: `1 - q * a * b`. -/
def den (q a b : R) : R := 1 - q * a * b

theorem num_def (p a b : R) : num p a b = a + b + p * a * b := rfl

theorem den_def (q a b : R) : den q a b = 1 - q * a * b := rfl

/-! ### Every point of the plane is a group law -/

theorem num_comm (p a b : R) : num p a b = num p b a := wAdd_comm p a b

theorem den_comm (q a b : R) : den q a b = den q b a := by
  simp only [den_def]; ring

/-- Identity: `F(a, 0) = a / 1`. -/
@[simp]
theorem num_zero (p a : R) : num p a 0 = a := wAdd_zero p a

@[simp]
theorem den_zero (q a : R) : den q a 0 = 1 := by simp [den_def]

/-- Associativity in projective form, numerators. Composing `(a ⊕ b) ⊕ c` with fractions `N / D`
gives numerator `N + c * D + p * N * c` and denominator `D - q * N * c`; the other bracketing is
symmetric. -/
theorem assoc_num (p q a b c : R) :
    num p a b + c * den q a b + p * num p a b * c =
      a * den q b c + num p b c + p * a * num p b c := by
  simp only [num_def, den_def]; ring

/-- Associativity in projective form, denominators. -/
theorem assoc_den (p q a b c : R) :
    den q a b - q * num p a b * c = den q b c - q * a * num p b c := by
  simp only [num_def, den_def]; ring

/-- Inverse: `F(a, -a / (1 + p * a)) = 0`. Cross-multiplied, the numerator vanishes. -/
theorem inverse_num (p a : R) : a * (1 + p * a) + -a + p * a * -a = 0 := by ring

/-! ### The growth law behind each point: `y' = 1 + p * y + q * y ^ 2` -/

/-- By the quotient rule, `∂F/∂b` at `b = 0` is `(N_b * D - N * D_b) / D ^ 2` with
`N_b = 1 + p * a`, `D = 1`, `N = a`, `D_b = -q * a`. -/
theorem riccati_rate (p q a : R) :
    (1 + p * a) * 1 - a * -(q * a) = 1 + p * a + q * a ^ 2 := by ring

/-! ### The lens and the walls

Write `p = -(u + v)` and `q = u * v`, so that `1 + p * a + q * a ^ 2 = (1 - u * a) * (1 - v * a)`:
the walls sit at `a = 1 / u` and `a = 1 / v`. The lens `M x = (1 - u * x) / (1 - v * x)` turns the
law into multiplication. -/

theorem rate_factors (u v a : R) :
    1 + -(u + v) * a + u * v * a ^ 2 = (1 - u * a) * (1 - v * a) := by ring

/-- The lift factors: `D - u * N = (1 - u * a) * (1 - u * b)`. -/
theorem lens_factor_u (u v a b : R) :
    den (u * v) a b - u * num (-(u + v)) a b = (1 - u * a) * (1 - u * b) := by
  simp only [num_def, den_def]; ring

/-- The lift factors: `D - v * N = (1 - v * a) * (1 - v * b)`. -/
theorem lens_factor_v (u v a b : R) :
    den (u * v) a b - v * num (-(u + v)) a b = (1 - v * a) * (1 - v * b) := by
  simp only [num_def, den_def]; ring

/-- The lens lift `M (a ⊕ b) = M a * M b`, cross-multiplied. -/
theorem lens_lift (u v a b : R) :
    (den (u * v) a b - u * num (-(u + v)) a b) * ((1 - v * a) * (1 - v * b)) =
      (den (u * v) a b - v * num (-(u + v)) a b) * ((1 - u * a) * (1 - u * b)) := by
  rw [lens_factor_u, lens_factor_v, mul_comm]

/-- Mult-add is the point `u = -1`, `v = 0`: walls at `-1` and at infinity. -/
theorem multAdd_point (a b : R) : num (-(-1 + 0)) a b = a ✳ b ∧ den (-1 * 0) a b = 1 := by
  simp [num, den]

/-- Velocity addition is `u = 1`, `v = -1`: walls at `1` and `-1` (the speed of light). -/
theorem velocity_point (a b : R) :
    num (-(1 + -1)) a b = a + b ∧ den (1 * -1) a b = 1 + a * b := by
  simp [num, den]

/-- The discriminant of the walls is `(u - v) ^ 2`. -/
theorem discriminant (u v : R) : (-(u + v)) ^ 2 - 4 * (u * v) = (u - v) ^ 2 := by ring

/-- On the parabola (a double wall, `u = v`) the lens `N x = x / (1 - u * x)` turns the law into
addition: the numerator is `a * (1 - u * b) + b * (1 - u * a)`. -/
theorem parabola_additive (u a b : R) :
    num (-(u + u)) a b = a * (1 - u * b) + b * (1 - u * a) := by
  simp only [num_def]; ring

/-- On the parabola the lifted denominator is `(1 - u * a) * (1 - u * b)`. -/
theorem parabola_den (u a b : R) :
    den (u * u) a b - u * num (-(u + u)) a b = (1 - u * a) * (1 - u * b) :=
  lens_factor_u u u a b

/-! ### Units -/

/-- Rescaling by `l` moves a point along a parabola through the origin:
`(p, q) ↦ (l * p, l ^ 2 * q)`. -/
theorem units_num (p l a b : R) : num p (l * a) (l * b) = l * num (l * p) a b := by
  simp only [num_def]; ring

theorem units_den (q l a b : R) : den q (l * a) (l * b) = den (l ^ 2 * q) a b := by
  simp only [den_def]; ring

/-- The dial: `k * wAdd k a b = (k * a) ✳ (k * b)`. Every `k ≠ 0` on the straight line is mult-add
in other units. This is `MultAdd.mul_wAdd`. -/
theorem dial_is_multAdd (k a b : R) : k * wAdd k a b = (k * a) ✳ (k * b) := mul_wAdd k a b

/-- The KdV two-soliton corner: `A * (f - 1)` with `f - 1 = e₁ + e₂ + A * e₁ * e₂` is plain
mult-add of `A * e₁` and `A * e₂`. This is `MultAdd.mul_wAdd` again. -/
theorem kdv_phase_shift (A e₁ e₂ : R) : A * wAdd A e₁ e₂ = (A * e₁) ✳ (A * e₂) :=
  mul_wAdd A e₁ e₂

/-! ### Hirzebruch's `χ_y` line

`(u + v + (1 - y) * u * v) / (1 + y * u * v)` sits at `(1 - y, -y)`, on the line `q = p - 1`. Its
discriminant `(1 - y) ^ 2 + 4 * y = (1 + y) ^ 2` vanishes only at `y = -1`, the Euler
characteristic. -/

theorem chi_y_line (y : R) : -y = 1 - y - 1 := by ring

theorem chi_y_touch (y : R) : (1 - y) ^ 2 - 4 * -y = (1 + y) ^ 2 := by ring

/-! ### The rung curve

Rung `s` sits at `(s, s * (s - 1) / 4)`. Scaled by `4` to avoid the division, its discriminant
(scaled by `16`) is `16 * s`, so its unit-free shape `p ^ 2 / (p ^ 2 - 4 * q)` is `s` itself. -/

theorem rung_discriminant (s : R) : (4 * s) ^ 2 - 4 * (4 * (s * (s - 1))) = 16 * s := by ring

theorem rung_equivalent (s : R) : (4 * s) ^ 2 = s * (16 * s) := by ring

/-! ### The floor of the straight line

An `L`-letter symmetric channel with error chance `n / m` (`0 ≤ n ≤ m`) has weight `-L / (L - 1)`
and eigenvalue `1 - L * (n / m) / (L - 1)`. It stays at or above `-1`, for every `L ≥ 2`:
symmetric noise cannot pass the parity wall. -/

/-- The channel floor, cross-multiplied: `L * n ≤ 2 * (L - 1) * m`. -/
theorem channel_floor [LinearOrder R] [IsStrictOrderedRing R] {L n m : R} (hL : 2 ≤ L)
    (hn : 0 ≤ n) (hnm : n ≤ m) : L * n ≤ 2 * (L - 1) * m := by
  nlinarith [mul_nonneg (sub_nonneg.2 hL) (hn.trans hnm),
    mul_nonneg (zero_le_two.trans hL) (sub_nonneg.2 hnm)]

/-- At `L = 2` (parity) the floor is reached exactly at full error: eigenvalue `-1`. -/
theorem parity_reaches_floor (m : R) : 2 * m = 2 * (2 - 1) * m := by ring

end LandscapeMap
