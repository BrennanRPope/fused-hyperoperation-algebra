/-
  Fused Hyperoperation Algebra: the landscape map.
  Core Lean 4 only (no mathlib), no `sorry`. NOT YET COMPILED: Lean was not installed in the
  workspace where this was written. Style follows the project's compiled files (core `Int`, `grind`).
  Every identity below was also checked symbolically with sympy (../algebra_checks.py).

  The map's plane: a point (p, q) is the law
      F_{p,q}(a, b) = (a + b + p·a·b) / (1 − q·a·b)
  p = weight of the corner cell ab (the dial k), q = weight of the next-layer cells a²b, ab² (the bend).
  To stay inside core Lean the laws are written as numerator/denominator pairs over `Int`,
  and every statement is a denominator-free (cross-multiplied) polynomial identity.
-/
namespace LandscapeMap

/-- numerator and denominator of F_{p,q}(a, b) -/
def num (p a b : Int) : Int := a + b + p * a * b
def den (q a b : Int) : Int := 1 - q * a * b

/-! ### 1. Every point of the plane is a group law -/

/-- commutativity, numerator and denominator separately -/
theorem num_comm (p a b : Int) : num p a b = num p b a := by unfold num; grind
theorem den_comm (q a b : Int) : den q a b = den q b a := by unfold den; grind

/-- identity: F(a, 0) = a / 1 -/
theorem num_zero (p a : Int) : num p a 0 = a := by unfold num; grind
theorem den_zero (q a : Int) : den q a 0 = 1 := by unfold den; grind

/-- Associativity, projective form. Composing (a ⊕ b) ⊕ c with fractions N/D gives
    numerator N + c·D + p·N·c and denominator D − q·N·c; the other bracketing is symmetric.
    The two numerators are equal and the two denominators are equal, as polynomials. -/
theorem assoc_num (p q a b c : Int) :
    num p a b + c * den q a b + p * num p a b * c
      = a * den q b c + num p b c + p * a * num p b c := by
  unfold num den; grind

theorem assoc_den (p q a b c : Int) :
    den q a b - q * num p a b * c = den q b c - q * a * num p b c := by
  unfold num den; grind

/-- inverse: F(a, −a/(1 + p a)) = 0, cross-multiplied: the numerator vanishes -/
theorem inverse_num (p a : Int) : a * (1 + p * a) + (-a) + p * a * (-a) = 0 := by grind

/-! ### 2. The growth law behind each point: y' = 1 + p·y + q·y²
    By the quotient rule, ∂F/∂b at b = 0 is (N_b·D − N·D_b)/D² with N_b = 1 + p a, D = 1,
    N = a, D_b = −q a. -/
theorem riccati_rate (p q a : Int) : (1 + p * a) * 1 - a * (-(q * a)) = 1 + p * a + q * a * a := by grind

/-! ### 3. The lens: walls
    Write p = −(u + v), q = u·v, so 1 + p a + q a² = (1 − u a)(1 − v a): the walls sit at a = 1/u, 1/v.
    The lens M(x) = (1 − u x)/(1 − v x) turns the law into multiplication: M(a ⊕ b) = M(a)·M(b). -/
theorem rate_factors (u v a : Int) : 1 + (-(u + v)) * a + (u * v) * a * a = (1 - u * a) * (1 - v * a) := by grind

/-- the lift, cross-multiplied: (D − u N)(1 − v a)(1 − v b) = (D − v N)(1 − u a)(1 − u b) -/
theorem lens_lift (u v a b : Int) :
    (den (u * v) a b - u * num (-(u + v)) a b) * ((1 - v * a) * (1 - v * b))
      = (den (u * v) a b - v * num (-(u + v)) a b) * ((1 - u * a) * (1 - u * b)) := by
  unfold num den; grind

/-- the lift actually factors: D − u N = (1 − u a)(1 − u b), and the same with v -/
theorem lens_factor_u (u v a b : Int) :
    den (u * v) a b - u * num (-(u + v)) a b = (1 - u * a) * (1 - u * b) := by
  unfold num den; grind

theorem lens_factor_v (u v a b : Int) :
    den (u * v) a b - v * num (-(u + v)) a b = (1 - v * a) * (1 - v * b) := by
  unfold num den; grind

/-- mult-add is the point u = −1, v = 0: walls at −1 and at infinity -/
theorem multadd_point (a b : Int) : num (-((-1) + 0)) a b = a + b + a * b ∧ den ((-1) * 0) a b = 1 := by
  unfold num den; constructor <;> grind

/-- velocity addition is u = 1, v = −1: walls at +1 and −1 (the speed of light) -/
theorem velocity_point (a b : Int) : num (-(1 + (-1))) a b = a + b ∧ den (1 * (-1)) a b = 1 + a * b := by
  unfold num den; constructor <;> grind

/-- the parabola: the discriminant of the walls is (u − v)² -/
theorem discriminant (u v : Int) : (-(u + v)) * (-(u + v)) - 4 * (u * v) = (u - v) * (u - v) := by grind

/-- the parabola itself (a double wall, u = v): the lens N(x) = x/(1 − u x) turns the law into addition.
    Cross-multiplied: N·(1 − u a)(1 − u b) = (a(1 − u b) + b(1 − u a))·(D − u N), with D − u N = (1 − u a)(1 − u b). -/
theorem parabola_additive (u a b : Int) :
    num (-(u + u)) a b = a * (1 - u * b) + b * (1 - u * a) := by
  unfold num; grind

theorem parabola_den (u a b : Int) :
    den (u * u) a b - u * num (-(u + u)) a b = (1 - u * a) * (1 - u * b) := by
  unfold num den; grind

/-! ### 4. Units -/

/-- rescaling by λ moves a point along a parabola through the origin: (p, q) ↦ (λ p, λ² q) -/
theorem units_num (p lam a b : Int) : num p (lam * a) (lam * b) = lam * num (lam * p) a b := by
  unfold num; grind
theorem units_den (q lam a b : Int) : den q (lam * a) (lam * b) = den (lam * lam * q) a b := by
  unfold den; grind

/-- the dial: k·(a ⊕ₖ b) = (k a) ⊕ (k b); every k ≠ 0 on the straight line is mult-add in other units -/
theorem dial_is_multadd (k a b : Int) : k * (a + b + k * a * b) = k * a + k * b + (k * a) * (k * b) := by grind

/-- the KdV soliton corner: A(f − 1) with f − 1 = e₁ + e₂ + A e₁e₂ is plain mult-add of A e₁ and A e₂
    (each wave moved by its phase shift log A) -/
theorem kdv_phase_shift (A e1 e2 : Int) : A * (e1 + e2 + A * e1 * e2) = A * e1 + A * e2 + (A * e1) * (A * e2) := by grind

/-! ### 5. Hirzebruch's χ_y line: (u + v + (1 − y)uv)/(1 + y uv) sits at (1 − y, −y), on c₂ = c₁ − 1;
    its discriminant (1 − y)² + 4y = (1 + y)² vanishes only at y = −1, the Euler characteristic. -/
theorem chi_y_line (y : Int) : -y = (1 - y) - 1 := by grind
theorem chi_y_touch (y : Int) : (1 - y) * (1 - y) - 4 * (-y) = (1 + y) * (1 + y) := by grind

/-! ### 6. The rung curve: rung s sits at (s, s(s − 1)/4). Scaled by 4 to stay in `Int`:
    its discriminant (scaled by 16) is 16 s, so its rung-equivalent height p²/(p² − 4q) is s itself. -/
theorem rung_discriminant (s : Int) : (4 * s) * (4 * s) - 4 * (4 * (s * (s - 1))) = 16 * s := by grind
/-- so the unit-free shape p²/(p² − 4q) of rung s is s: (4s)² = s·(16 s) -/
theorem rung_equivalent (s : Int) : (4 * s) * (4 * s) = s * (16 * s) := by grind

/-! ### 7. The floor of the straight line: symmetric noise cannot pass the parity wall.
    An L-letter symmetric channel with error chance p = n/m (0 ≤ n ≤ m) has weight −L/(L − 1) and
    eigenvalue 1 − L p/(L − 1). It stays ≥ −1, i.e. L·n ≤ 2(L − 1)·m, for every L ≥ 2. -/
theorem channel_floor (L n m : Int) (hL : 2 ≤ L) (hn0 : 0 ≤ n) (hnm : n ≤ m) :
    L * n ≤ 2 * (L - 1) * m := by
  have hm0 : 0 ≤ m := Int.le_trans hn0 hnm
  have h1 : L * n ≤ L * m := Int.mul_le_mul_of_nonneg_left hnm (by omega)
  have h2 : L * m ≤ (2 * (L - 1)) * m := Int.mul_le_mul_of_nonneg_right (by omega) hm0
  exact Int.le_trans h1 h2

/-- at L = 2 (parity) the floor is reached exactly at full error: eigenvalue −1 -/
theorem parity_reaches_floor (m : Int) : 2 * m = 2 * (2 - 1) * m := by grind

end LandscapeMap
