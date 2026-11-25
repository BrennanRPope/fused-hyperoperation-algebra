import ProofComplexity.FrameworkPlane
import Mathlib.Analysis.Complex.Trigonometric
import Mathlib.Tactic.LinearCombination

/-!
# R6. The tangent addition formula, both ways

Same statement for both: for complex `x, y` with `cos x, cos y, cos (x + y) ≠ 0`,
`tan (x + y) = (tan x + tan y) / (1 − tan x · tan y)`.
-/

open Complex

namespace HCT.Results

/-- Standard route (as in Mathlib's `Complex.tan_add`): sine and cosine addition, then
clear denominators. -/
theorem std_tan_add {x y : ℂ} (hx : cos x ≠ 0) (hy : cos y ≠ 0) (hxy : cos (x + y) ≠ 0) :
    tan (x + y) = (tan x + tan y) / (1 - tan x * tan y) := by
  rw [cos_add] at hxy
  rw [tan_eq_sin_div_cos, tan_eq_sin_div_cos, tan_eq_sin_div_cos, sin_add, cos_add]
  field_simp

/-- Through the operation: the plane at `(k, d) = (0, 1)` with `α = i, β = −i`. The lift
`t ↦ (1 + i t)/(1 − i t)` sends `tan θ` to `exp (iθ)²`, so the tangent of a sum is the
plane law of the tangents. -/
theorem op_tan_add {x y : ℂ} (hx : cos x ≠ 0) (hy : cos y ≠ 0) (hxy : cos (x + y) ≠ 0) :
    tan (x + y) = (tan x + tan y) / (1 - tan x * tan y) := by
  -- the lift of tan θ is exp(iθ)²
  have hL : ∀ θ : ℂ, cos θ ≠ 0 → 1 + I * tan θ = (1 - I * tan θ) * exp (θ * I) ^ 2 := by
    intro θ hθ
    rw [tan_eq_sin_div_cos, exp_mul_I]
    field_simp
    linear_combination (cos θ + I * sin θ) * sin θ ^ 2 * I_sq - (cos θ + I * sin θ) * sin_sq_add_cos_sq θ
  -- the lift of a tangent is never at the pole
  have hpole : ∀ θ : ℂ, cos θ ≠ 0 → 1 - I * tan θ ≠ 0 := by
    intro θ hθ h
    have h1 := hL θ hθ
    rw [h, zero_mul] at h1
    have : (2 : ℂ) = 0 := by linear_combination h1 + h
    norm_num at this
  set a := tan x; set b := tan y; set T := tan (x + y)
  have hD : 1 - a * b ≠ 0 := by
    intro h0; apply hxy
    rw [cos_add]
    have : sin x * sin y = cos x * cos y * (a * b) := by
      simp only [a, b, tan_eq_sin_div_cos]; field_simp
    rw [this]; linear_combination cos x * cos y * h0
  -- the plane law with α = i, β = −i (k = 0, d = 1)
  have hP := HCT.plane_lift I (-I) a b
  dsimp only at hP
  simp only [add_neg_cancel, zero_mul, add_zero, mul_neg, I_mul_I, neg_neg] at hP
  -- exp turns the sum of angles into a product
  have hx' := hL x hx; have hy' := hL y hy; have hT := hL (x + y) hxy
  rw [add_mul, exp_add] at hT
  set A := exp (x * I); set B := exp (y * I)
  have key : (1 - I * a) * (1 - I * b) * (2 * I * (T * (1 - a * b) - (a + b))) = 0 := by
    linear_combination (1 - I * T) * (1 - a * b - I * (a + b)) *
        (-(1 + I * b) * hx' - (1 - I * a) * A ^ 2 * hy') - (1 - I * T) * hP +
      (1 - I * a) * (1 - I * b) * (1 - a * b - I * (a + b)) * hT
  rw [eq_div_iff hD]
  have h2I : (2 : ℂ) * I ≠ 0 := by simp [I_ne_zero]
  rcases mul_eq_zero.1 key with h | h
  · exact absurd h (mul_ne_zero (hpole x hx) (hpole y hy))
  · exact sub_eq_zero.1 ((mul_eq_zero.1 h).resolve_left h2I)

end HCT.Results
