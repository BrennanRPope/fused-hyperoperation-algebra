module

public import FusedAlgebra.Basic

/-!
# Solitons on the rungs, and the gluing number `q`

KdV soliton tau functions (Hirota 1971) are polynomials in the exponentials `a = exp η₁`,
`b = exp η₂`, `c = exp η₃`, with pair weights `A i j = ((k i - k j) / (k i + k j)) ^ 2`. The
soliton formulas are standard; what this file records is how they sit on mult-add.

## Main results

* `tau2_eq_one_add_wAdd`: the two-soliton tau function is the lift of the weighted mult-add;
  `mul_tau2_sub_one` is the same law in a rescaled coordinate.
* `tau3_equal_corner`, `wAdd_three_corner`, `soliton_vs_wAdd_corner`: for three solitons every
  pair carries a weight, so the triple corner is `k ^ 3`, while the weighted mult-add of three
  has `k ^ 2`; the two agree exactly for `k = 0` and `k = 1`.
* `three_body_roads`, `rung_three_body`, `four_body_roads`: the corners of the rung flow against
  the corners of the soliton.
* `face_iff`, `theta2_face`: on the theta lattice only the mult-add cube survives at `q = 0`.
* `glue_is_multAdd`, `orbit_product`: gluing by `q` is a mult-add step by `q - 1`.
* `kdv_weight_below_one`: KdV pair weights lie in `[0, 1)`.

The core-Lean original worked over `ℤ`. Here the identities hold in any commutative ring, the
"exactly when" statements in any commutative ring without zero divisors, the inequalities in any
linearly ordered commutative ring; only the lattice statements `tri_nonneg` and `face_iff` are
about `ℤ`.
-/

@[expose] public section

namespace SolitonsRungs

open MultAdd Finset

variable {R : Type*}

section CommRing

variable [CommRing R]

/-! ### 1. Two solitons: the tau function is weighted mult-add -/

/-- The two-soliton tau function `1 + a + b + A * a * b`. -/
def tau2 (A a b : R) : R := 1 + a + b + A * a * b

theorem tau2_eq_one_add_wAdd (A a b : R) : tau2 A a b = 1 + wAdd A a b := by
  unfold tau2 wAdd; ring

/-- Weight `1`: no interaction, the tau function is the plain mult-add lift
`(1 + a) * (1 + b)`. -/
theorem tau2_free (a b : R) : tau2 1 a b = (1 + a) * (1 + b) := by
  unfold tau2; ring

/-- Weight `0`: equal speeds, the two exponentials just add (rung 0). -/
theorem tau2_merge (a b : R) : tau2 0 a b = 1 + (a + b) := by
  unfold tau2; ring

/-- For two solitons the law is rung 1 in a rescaled coordinate:
`A * (τ - 1) = (A * a) ✳ (A * b)`. This is `MultAdd.mul_wAdd` read through `tau2`. -/
theorem mul_tau2_sub_one (A a b : R) : A * (tau2 A a b - 1) = (A * a) ✳ (A * b) := by
  rw [← mul_wAdd, tau2_eq_one_add_wAdd, add_sub_cancel_left]

/-- Read in the `1 + x` picture, the weight is mult-add plus an excess `A - 1` for the pair. -/
theorem tau2_excess (A a b : R) : tau2 A a b = 1 + (a ✳ b) + (A - 1) * a * b := by
  unfold tau2 multAdd; ring

/-! ### 2. Three solitons: corners multiply over pairs -/

/-- The three-soliton tau function: every corner carries the product of its pair weights. -/
def tau3 (A₁₂ A₁₃ A₂₃ a b c : R) : R :=
  1 + a + b + c + A₁₂ * a * b + A₁₃ * a * c + A₂₃ * b * c + A₁₂ * A₁₃ * A₂₃ * (a * b * c)

theorem tau3_free (a b c : R) : tau3 1 1 1 a b c = (1 + a) * (1 + b) * (1 + c) := by
  unfold tau3; ring

/-- A soliton that interacts with nobody (weights `1`) factors off. -/
theorem tau3_splits (A₁₂ a b c : R) : tau3 A₁₂ 1 1 a b c = tau2 A₁₂ a b * (1 + c) := by
  unfold tau3 tau2; ring

/-- Equal weights `k`: the triple corner is `k ^ 3` (three pairs). -/
theorem tau3_equal_corner (k a b c : R) :
    tau3 k k k a b c = 1 + a + b + c + k * (a * b + a * c + b * c) + k ^ 3 * (a * b * c) := by
  unfold tau3; ring

/-- The weighted mult-add of three has triple corner `k ^ 2`, not `k ^ 3`: it is mult-add
rescaled, while solitons weight every pair. -/
theorem wAdd_three_corner (k a b c : R) :
    wAdd k (wAdd k a b) c = a + b + c + k * (a * b + a * c + b * c) + k ^ 2 * (a * b * c) := by
  unfold wAdd; ring

/-! ### 3. The corners of the rung flow against the corners of the soliton

For rung `s` and `n` inputs, the corner (multilinear) coefficient is `n ! * [y ^ n] E ^ s (y)`:
`s` for `n = 2`, `s * (3 * s - 1) / 2` for `n = 3`, `s * (2 * s - 1) * (3 * s - 1) / 2` for
`n = 4`. Solitons with equal weight `A` have `A`, `A ^ 3`, `A ^ 6` (one `A` for each pair). -/

/-- Three bodies: `2 * (soliton corner - rung corner) = s * (2 * s - 1) * (s - 1)`. The roads
meet at `s = 0`, `1 / 2`, `1`. -/
theorem three_body_roads (s : R) : 2 * s ^ 3 - s * (3 * s - 1) = s * (2 * s - 1) * (s - 1) := by
  ring

/-- The three-body term of the rung flow: the triple corner minus the square of the pair corner,
times `2`, is `s * (s - 1)`. -/
theorem rung_three_body (s : R) : s * (3 * s - 1) - 2 * s ^ 2 = s * (s - 1) := by
  ring

/-- Four bodies: `2 * (soliton corner - rung corner)` factors as `s * (s - 1)` times a
quartic. -/
theorem four_body_roads (s : R) :
    2 * s ^ 6 - s * (2 * s - 1) * (3 * s - 1) =
      s * (s - 1) * (2 * s ^ 4 + 2 * s ^ 3 + 2 * s ^ 2 - 4 * s + 1) := by
  ring

/-- The quartic of `four_body_roads` is not zero at `s = 1 / 2` (its value there is `-1 / 8`),
so the three-body meeting at `s = 1 / 2` does not survive to four bodies. -/
theorem four_body_half :
    2 * (1 / 2 : ℚ) ^ 4 + 2 * (1 / 2) ^ 3 + 2 * (1 / 2) ^ 2 - 4 * (1 / 2) + 1 = -1 / 8 := by
  norm_num

/-! ### 4. The theta lattice: only the mult-add cube survives at `q = 0` -/

/-- The two-wave theta series on the lattice points `0, …, 3` in each direction, with weight
`q₁ ^ (n₁ * (n₁ - 1) / 2) * q₂ ^ (n₂ * (n₂ - 1) / 2) * A ^ (n₁ * n₂) * w₁ ^ n₁ * w₂ ^ n₂`. -/
def theta2 (q₁ q₂ A w₁ w₂ : R) : R :=
  ∑ n₁ ∈ range 4, ∑ n₂ ∈ range 4,
    q₁ ^ (n₁ * (n₁ - 1) / 2) * q₂ ^ (n₂ * (n₂ - 1) / 2) * A ^ (n₁ * n₂) * w₁ ^ n₁ * w₂ ^ n₂

/-- At `q₁ = q₂ = 0` the theta series is the weighted mult-add cube `1 + wAdd A w₁ w₂`. -/
theorem theta2_face (A w₁ w₂ : R) : theta2 0 0 A w₁ w₂ = 1 + wAdd A w₁ w₂ := by
  simp [theta2, sum_range_succ, wAdd_def]; ring

/-- With `A = 1` the two-wave theta splits into two one-wave pieces (no interaction). -/
theorem theta2_free_face (w₁ w₂ : R) : theta2 0 0 1 w₁ w₂ = (1 + w₁) * (1 + w₂) := by
  rw [theta2_face, wAdd_one, one_add_multAdd]

/-! ### 5. Gluing is a rung-1 move -/

/-- Multiplying by `q` in the `w` picture is mult-add by `q - 1` in the `1 + x` picture. -/
theorem glue_is_multAdd (w q : R) : w * q - 1 = (w - 1) ✳ (q - 1) := by
  unfold multAdd; ring

/-- A truncated theta product: `(1 + w) * (1 + w * q) * (1 + v * q)`, with `v` standing for
`1 / w`, is `1 +` the mult-add of three points of the `q`-orbit of `w`. -/
theorem orbit_product (w v q : R) :
    (1 + w) * (1 + w * q) * (1 + v * q) = 1 + (w ✳ (w * q ✳ v * q)) := by
  rw [one_add_multAdd, one_add_multAdd, mul_assoc]

/-- At `q = 0` only the soliton factor `1 + w` is left. -/
theorem orbit_product_q0 (w v : R) : (1 + w) * (1 + w * 0) * (1 + v * 0) = 1 + w := by
  ring

end CommRing

/-! ### Where the roads meet -/

section NoZeroDivisors

variable [CommRing R] [NoZeroDivisors R]

/-- The triple corners of the soliton and of the weighted mult-add agree exactly when
`k ^ 3 = k ^ 2`, that is for `k = 0` or `k = 1` (rung 0 or rung 1). -/
theorem soliton_vs_wAdd_corner (k : R) : k ^ 3 = k ^ 2 ↔ k = 0 ∨ k = 1 := by
  constructor
  · intro e
    have z : k ^ 2 * (k - 1) = 0 := by linear_combination e
    rcases mul_eq_zero.1 z with h | h
    · exact Or.inl (pow_eq_zero_iff two_ne_zero |>.1 h)
    · exact Or.inr (sub_eq_zero.1 h)
  · rintro (rfl | rfl) <;> simp

/-- The rung flow has no three-body term only at rungs `0` and `1`. -/
theorem rung_three_body_zero (s : R) : s * (s - 1) = 0 ↔ s = 0 ∨ s = 1 := by
  rw [mul_eq_zero, sub_eq_zero]

/-- The lattice weight `q ^ (n * (n - 1) / 2)` is free of `q` only for `n = 0` and `n = 1`. -/
theorem tri_zero_iff (n : R) : n * (n - 1) = 0 ↔ n = 0 ∨ n = 1 := rung_three_body_zero n

/-- The KdV pair weight is `0` exactly when the speeds match. -/
theorem kdv_weight_zero (k₁ k₂ : R) : (k₁ - k₂) ^ 2 = 0 ↔ k₁ = k₂ := by
  rw [sq_eq_zero_iff, sub_eq_zero]

end NoZeroDivisors

/-! ### The theta lattice over `ℤ` -/

theorem tri_nonneg (n : ℤ) : 0 ≤ n * (n - 1) := by
  rcases le_or_gt n 0 with h | h
  · exact mul_nonneg_of_nonpos_of_nonpos h (by lia)
  · exact mul_nonneg (by lia) (by lia)

/-- Two waves: the exponent `n₁ * (n₁ - 1) + n₂ * (n₂ - 1)` vanishes exactly on the unit square
`{0, 1} × {0, 1}`, the corners of the mult-add cube. -/
theorem face_iff (n₁ n₂ : ℤ) :
    n₁ * (n₁ - 1) + n₂ * (n₂ - 1) = 0 ↔ (n₁ = 0 ∨ n₁ = 1) ∧ (n₂ = 0 ∨ n₂ = 1) := by
  rw [add_eq_zero_iff_of_nonneg (tri_nonneg n₁) (tri_nonneg n₂), tri_zero_iff, tri_zero_iff]

/-! ### 6. KdV pair weights lie in `[0, 1)` -/

/-- For positive wave numbers, `(k₁ - k₂) ^ 2 < (k₁ + k₂) ^ 2`: the weight `A` is below `1`.
(That it is nonnegative is `sq_nonneg`.) -/
theorem kdv_weight_below_one [CommRing R] [LinearOrder R] [IsStrictOrderedRing R] (k₁ k₂ : R)
    (h₁ : 0 < k₁) (h₂ : 0 < k₂) : (k₁ - k₂) ^ 2 < (k₁ + k₂) ^ 2 := by
  nlinarith [mul_pos h₁ h₂]

end SolitonsRungs
