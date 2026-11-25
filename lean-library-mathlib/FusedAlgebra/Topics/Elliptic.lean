module

public import FusedAlgebra.Basic

/-!
# The elliptic lead: cubics with a double point, and the Frey curve

Curves are Weierstrass cubics `y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆` (mathlib's
`WeierstrassCurve`), with parameter `z = -x / y`. A cubic with a double point at the origin has
`a₃ = a₄ = a₆ = 0`; its landscape point is `(k, q) = (-a₁, -a₂)`.

This is the mathlib form of `FusedAlgebra.Topics.Elliptic`. The original worked over `ℤ` with
hand-written polynomial identities; here the identities hold in any commutative ring and the
curves are mathlib `WeierstrassCurve`s, so the discriminant and `c₄` are mathlib's `Δ` and `c₄`.

## Main definitions

* `EllipticLead.nodalCurve a₁ a₂`: the cubic `y² + a₁xy = x³ + a₂x²` with a double point at `0`.
* `EllipticLead.freyCurve A B`: the Frey curve `y² = x(x - A)(x + B)`.
* `EllipticLead.d1`, …, `d5`: the first coefficients of the invariant differential.
* `EllipticLead.gap6`: the layer-6 gap of the cubic fitted to rung `s`.

## Main results

* `nodalCurve_Δ`, `node_param`: the nodal cubic is singular, and every line `y = tx` through the
  double point meets it again at `x = t² + a₁t - a₂`.
* `wall_poly_factor`, `slopes_factor`: the wall polynomial `1 + kz + qz²` factors through the two
  tangent slopes of the node.
* `tate_q0_is_multAdd`: the group law of the node with slopes `0` and `1` is mult-add.
* `freyCurve_Δ`, `freyCurve_c₄`: `Δ = 16A²B²(A + B)²` and `c₄ = 16(A² + AB + B²)`.
* `gap6_ne_zero`: the fitted cubic passes layer 6 only at the rungs `s = 0` and `s = 1`.

All the curve facts are classical (Silverman, *The Arithmetic of Elliptic Curves*, III and
App. A); what the file records is how they read on the mult-add landscape.
-/

@[expose] public section

namespace EllipticLead

open MultAdd

variable {R : Type*} [CommRing R]

/-! ### The landscape plane is the cubics with a double point -/

/-- The Weierstrass cubic `y² + a₁xy = x³ + a₂x²`, with a double point at the origin. -/
def nodalCurve (a₁ a₂ : R) : WeierstrassCurve R := ⟨a₁, a₂, 0, 0, 0⟩

/-- A cubic with a double point is singular: its discriminant vanishes. -/
theorem nodalCurve_Δ (a₁ a₂ : R) : (nodalCurve a₁ a₂).Δ = 0 := by
  simp only [nodalCurve, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- `c₄` of the nodal cubic is the square of the slope discriminant `a₁² + 4a₂`: the double point
is a node when it is nonzero and a cusp when it vanishes. -/
theorem nodalCurve_c₄ (a₁ a₂ : R) : (nodalCurve a₁ a₂).c₄ = (a₁ ^ 2 + 4 * a₂) ^ 2 := by
  simp only [nodalCurve, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  ring

/-- Every line `y = tx` through the double point meets the curve again at
`x = t² + a₁t - a₂`. -/
theorem node_param (a₁ a₂ t : R) :
    (nodalCurve a₁ a₂).toAffine.Equation (t ^ 2 + a₁ * t - a₂) (t * (t ^ 2 + a₁ * t - a₂)) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  simp only [nodalCurve]
  ring

/-- The walls are the node's branches: with `z = -1 / t` (so `zt = -1`),
`t²(1 + kz + qz²) - (t² + a₁t - a₂)` is a multiple of `zt + 1`, where `k = -a₁`, `q = -a₂`. -/
theorem walls_are_branches (a₁ a₂ t z : R) :
    t ^ 2 * (1 + -a₁ * z + -a₂ * z ^ 2) - (t ^ 2 + a₁ * t - a₂) =
      -(a₁ * t) * (z * t + 1) - a₂ * ((z * t) ^ 2 - 1) := by
  ring

/-- The two tangent slopes `m₁`, `m₂` factor the slope polynomial: `a₁ = -(m₁ + m₂)` and
`a₂ = -(m₁m₂)`, so `k = m₁ + m₂` and `q = m₁m₂`. -/
theorem slopes_factor (m₁ m₂ t : R) :
    t ^ 2 + -(m₁ + m₂) * t - -(m₁ * m₂) = (t - m₁) * (t - m₂) := by
  ring

/-- The wall polynomial factors through the branches: `1 + kz + qz² = (1 + m₁z)(1 + m₂z)`. -/
theorem wall_poly_factor (m₁ m₂ z : R) :
    1 + (m₁ + m₂) * z + m₁ * m₂ * z ^ 2 = (1 + m₁ * z) * (1 + m₂ * z) := by
  ring

/-- The discriminant of the slopes, which is `b₂` of the nodal cubic, is the landscape's
`k² - 4q`. -/
theorem slope_discriminant (a₁ a₂ : R) : (nodalCurve a₁ a₂).b₂ = (-a₁) ^ 2 - 4 * -a₂ := by
  simp only [nodalCurve, WeierstrassCurve.b₂]
  ring

/-- The nodal curve's formal group: `w = -1 / y` is a rational function of `z`, because
`w(1 - a₁z - a₂z²) = z³` is the curve equation `w = z³ + a₁zw + a₂z²w`. -/
theorem nodal_w (a₁ a₂ z w : R) (h : w * (1 - a₁ * z - a₂ * z ^ 2) = z ^ 3) :
    w = z ^ 3 + a₁ * z * w + a₂ * z ^ 2 * w := by
  linear_combination h

/-! ### Smoothing the node: Tate's curve at `q = 0` is mult-add -/

/-- With `z = 1 - v` (`v = 1 / u`, Tate's parameter inverted), the law `z₁ + z₂ - z₁z₂` is
multiplication `v₁v₂`. -/
theorem tate_q0_law (v₁ v₂ : R) : (1 - v₁) + (1 - v₂) - (1 - v₁) * (1 - v₂) = 1 - v₁ * v₂ := by
  ring

/-- And `a = -z` turns it into mult-add exactly. -/
theorem tate_q0_is_multAdd (a b : R) : -(-a + -b - -a * -b) = a ✳ b := by
  unfold multAdd; ring

/-- First coefficient of the invariant differential `(1 + d₁z + d₂z² + ⋯) dz`. The coefficient
`dₙ` has weight `n`. -/
def d1 (a₁ : R) : R := a₁

/-- Second coefficient of the invariant differential. -/
def d2 (a₁ a₂ : R) : R := a₁ ^ 2 + a₂

/-- Third coefficient of the invariant differential. -/
def d3 (a₁ a₂ a₃ : R) : R := a₁ ^ 3 + 2 * a₁ * a₂ + 2 * a₃

/-- Fourth coefficient of the invariant differential. -/
def d4 (a₁ a₂ a₃ a₄ : R) : R := a₁ ^ 4 + 3 * a₁ ^ 2 * a₂ + 6 * a₁ * a₃ + a₂ ^ 2 + 2 * a₄

/-- Fifth coefficient of the invariant differential. -/
def d5 (a₁ a₂ a₃ a₄ : R) : R :=
  a₁ ^ 5 + 4 * a₁ ^ 3 * a₂ + 12 * a₁ ^ 2 * a₃ + 3 * a₁ * a₂ ^ 2 + 6 * a₁ * a₄ + 6 * a₂ * a₃

/-- `d4` has weight `4`: scaling `aᵢ` by `uⁱ` scales it by `u⁴`. -/
theorem d4_weight (u a₁ a₂ a₃ a₄ : R) :
    d4 (u * a₁) (u ^ 2 * a₂) (u ^ 3 * a₃) (u ^ 4 * a₄) = u ^ 4 * d4 a₁ a₂ a₃ a₄ := by
  unfold d4; ring

/-- `d5` has weight `5`. -/
theorem d5_weight (u a₁ a₂ a₃ a₄ : R) :
    d5 (u * a₁) (u ^ 2 * a₂) (u ^ 3 * a₃) (u ^ 4 * a₄) = u ^ 5 * d5 a₁ a₂ a₃ a₄ := by
  unfold d5; ring

/-- Mult-add's own curve (`a₁ = -1`, `a₂ = a₃ = a₄ = 0`) has `dₙ = (-1)ⁿ`: the differential
`dz / (1 + z)`, the logarithm of mult-add, through layer 6. -/
theorem multAdd_curve_d :
    d1 (-1 : R) = -1 ∧ d2 (-1 : R) 0 = 1 ∧ d3 (-1 : R) 0 0 = -1 ∧ d4 (-1 : R) 0 0 0 = 1 ∧
      d5 (-1 : R) 0 0 0 = -1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> simp only [d1, d2, d3, d4, d5] <;> ring

/-! ### Connections -/

/-- Ochanine's elliptic genus at `ε = δ²` (the nodal case): the numerator
`a√R(b) + b√R(a)` with `√R(x) = 1 - δx²` is `(a + b)(1 - δab)`, so the law is
`(a + b) / (1 + δab)`, the signature at `δ = 1`. -/
theorem ochanine_nodal_num (δ a b : R) :
    a * (1 - δ * b ^ 2) + b * (1 - δ * a ^ 2) = (a + b) * (1 - δ * a * b) := by
  ring

/-- The denominator of Ochanine's law in the nodal case. -/
theorem ochanine_nodal_den (δ a b : R) :
    1 - δ ^ 2 * a ^ 2 * b ^ 2 = (1 - δ * a * b) * (1 + δ * a * b) := by
  ring

/-- The Frey curve `y² = x(x - A)(x + B) = x³ + (B - A)x² - ABx`. -/
def freyCurve (A B : R) : WeierstrassCurve R := ⟨0, B - A, 0, -(A * B), 0⟩

/-- The Frey curve is `y² = x(x - A)(x + B)`. -/
theorem freyCurve_equation (A B x y : R) :
    (freyCurve A B).toAffine.Equation x y ↔ y ^ 2 = x * (x - A) * (x + B) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  simp only [freyCurve]
  constructor <;> intro h <;> linear_combination h

/-- The discriminant of the Frey curve: `Δ = 16A²B²(A + B)²`. -/
theorem freyCurve_Δ (A B : R) : (freyCurve A B).Δ = 16 * A ^ 2 * B ^ 2 * (A + B) ^ 2 := by
  simp only [freyCurve, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- `c₄` of the Frey curve: `c₄ = 16(A² + AB + B²)`. -/
theorem freyCurve_c₄ (A B : R) : (freyCurve A B).c₄ = 16 * (A ^ 2 + A * B + B ^ 2) := by
  simp only [freyCurve, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  ring

/-- Semistability at a prime dividing `A`: there (`A ≡ 0`) `c₄ / 16` is `B²`, which is nonzero
modulo the prime when `A` and `B` are coprime, so the reduction is a node, never a cusp. -/
theorem freyCurve_c₄_zero_left (B : R) : (freyCurve 0 B).c₄ = 16 * B ^ 2 := by
  rw [freyCurve_c₄]; ring

/-- Semistability at a prime dividing `C = A + B`: there (`B ≡ -A`) `c₄ / 16` is `A²`. -/
theorem freyCurve_c₄_neg (A : R) : (freyCurve A (-A)).c₄ = 16 * A ^ 2 := by
  rw [freyCurve_c₄]; ring

/-- Euler factors are wall polynomials: `1 - aT + pT² = (1 - αT)(1 - βT)` when `a = α + β` and
`p = αβ`. -/
theorem euler_factor_walls (α β T : R) :
    1 - (α + β) * T + α * β * T ^ 2 = (1 - α * T) * (1 - β * T) := by
  ring

/-- Hasse's bound `a² < 4p` puts the point `(k, q) = (-a, p)` in the turning country:
`k² - 4q < 0`. -/
theorem hasse_turning [PartialOrder R] [IsOrderedRing R] (a p : R) (h : a ^ 2 < 4 * p) :
    (-a) ^ 2 - 4 * p < 0 := by
  rw [neg_sq]; exact sub_neg.2 h

/-- Units: scaling `T` by `l` moves `(k, q)` to `(lk, l²q)`, the landscape's units action. -/
theorem euler_units (k q l T : R) :
    1 + k * (l * T) + q * (l * T) ^ 2 = 1 + l * k * T + l ^ 2 * q * T ^ 2 := by
  ring

/-! ### Verdict: the fitted curve is singular exactly at the rungs it can copy

The cubic fitted to rung `s` through layer 5 has `a₁ = -s`, `a₂ = -s(s - 1) / 4`,
`a₃ = s(s - 1) / 24` and `a₄ = s(s - 1)(11s - 4) / 288`; its layer-6 gap is
`-s(s - 1)(45s² - 125s + 22) / 17280`. -/

/-- The layer-6 gap of the cubic fitted to rung `s`, cleared of its denominator `17280`. -/
def gap6 (s : ℤ) : ℤ := -(s * (s - 1) * (45 * s ^ 2 - 125 * s + 22))

theorem gap6_zero : gap6 0 = 0 := by simp [gap6]

theorem gap6_one : gap6 1 = 0 := by simp [gap6]

/-- No other whole rung passes layer 6: `45s² - 125s + 22` has no integer root (its roots are
near `0.19` and `2.59`). -/
theorem gap6_ne_zero {s : ℤ} (h0 : s ≠ 0) (h1 : s ≠ 1) : gap6 s ≠ 0 := by
  have hq : 45 * s ^ 2 - 125 * s + 22 ≠ 0 := by
    intro h
    obtain hs | hs | hs | hs : s ≤ 0 ∨ s = 1 ∨ s = 2 ∨ 3 ≤ s := by lia
    · nlinarith [sq_nonneg s]
    · exact h1 hs
    · subst hs; norm_num at h
    · nlinarith
  simp [gap6, h0, sub_ne_zero.2 h1, hq]

end EllipticLead
