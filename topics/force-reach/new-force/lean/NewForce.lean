/-
  Does the theory predict an unknown force? (DRAFT, Nov 26 2026)
  Stated against Mathlib and not compiled as written (CI builds no Mathlib copy of it). All fourteen statements
  are compiled in core Lean in lean-library/FusedAlgebra/Topics/NewForce.lean.
  Every statement is checked true by check_statements.py.
  Pure algebra over ℝ (ring, field_simp, norm_num).
-/
import Mathlib

namespace FusedAlgebra.Topics.NewForce

/-- The rung tower, every rung coupled equally, sums to coth:
    with q = e^{−r/R} (0 ≤ q < 1), 1 + 2q/(1−q) = (1+q)/(1−q), and (1+q)/(1−q) = coth(r/2R). -/
theorem tower_sum (q : ℝ) (hq : q ≠ 1) :
    1 + 2 * q / (1 - q) = (1 + q) / (1 - q) := by
  have h : (1 - q) ≠ 0 := sub_ne_zero.mpr (Ne.symm hq)
  field_simp
  ring

/-- Far away the tower is one Yukawa term of strength 2: (1+q)/(1−q) − 1 − 2q = 2q²/(1−q). -/
theorem tower_leading_yukawa (q : ℝ) (hq : q ≠ 1) :
    (1 + q) / (1 - q) - 1 - 2 * q = 2 * q ^ 2 / (1 - q) := by
  have h : (1 - q) ≠ 0 := sub_ne_zero.mpr (Ne.symm hq)
  field_simp
  ring

/-- Flux exponent s = 3/2 in D-language is field exponent q = s/(s−1) = 3. -/
theorem dual_exponent_mond : (3 / 2 : ℝ) / (3 / 2 - 1) = 3 := by norm_num

/-- Hölder pairing: 1/s + 1/q = 1 for s = 3/2, q = 3. -/
theorem holder_pair : (1 : ℝ) / (3 / 2) + 1 / 3 = 1 := by norm_num

/-- Asymptotic freedom count for SU(3): 11·3/3 − 2·16/3 > 0 and 11·3/3 − 2·17/3 < 0. -/
theorem su3_flavour_bound :
    (0 : ℝ) < 11 * 3 / 3 - 2 * 16 / 3 ∧ 11 * 3 / 3 - 2 * 17 / 3 < (0 : ℝ) := by
  constructor <;> norm_num

/-- Asymptotic freedom count for SU(2): 11·2/3 − 2·10/3 > 0, and at 11 flavours it is exactly 0. -/
theorem su2_flavour_bound :
    (0 : ℝ) < 11 * 2 / 3 - 2 * 10 / 3 ∧ 11 * 2 / 3 - 2 * 11 / 3 = (0 : ℝ) := by
  constructor <;> norm_num

/-! ### Part G: strength from a closed direction (electromagnetism as gravity's circle) -/

/-- Copies split into 4D pieces: tensor 1/2, vector 2, scalar 3/2 (per E₁E₂).
    The vector is 4 times the tensor and the scalar 3 times. -/
theorem vector_over_tensor : (2 : ℝ) / (1 / 2) = 4 := by norm_num
theorem scalar_over_tensor : (3 / 2 : ℝ) / (1 / 2) = 3 := by norm_num

/-- Like circle charges: tensor + scalar − vector = 0. Opposite: tensor + scalar + vector = 4. -/
theorem like_charges_cancel : (1 / 2 : ℝ) + 3 / 2 - 2 = 0 := by norm_num
theorem opposite_charges_add : (1 / 2 : ℝ) + 3 / 2 + 2 = 4 := by norm_num

/-- Klein's relation in our reading: with G = 1 and circle radius R = 2/√α,
    4/R² = α. Stated with R² = 4/α. -/
theorem klein_strength (α : ℝ) (hα : α ≠ 0) : 4 / (4 / α) = α := by
  field_simp

/-- n circles give a scalar force of strength a = n/(n+2); one circle gives 1/3. -/
theorem one_circle_third : (1 : ℝ) / (1 + 2) = 1 / 3 := by norm_num

/-- A scalar of strength a gives γ = (1 − a)/(1 + a); with a = n/(n+2) that is 1/(n+1). -/
theorem gamma_from_strength (n : ℝ) (hn : 0 ≤ n) :
    (1 - n / (n + 2)) / (1 + n / (n + 2)) = 1 / (n + 1) := by
  have h2 : n + 2 ≠ 0 := by positivity
  have h1 : n + 1 ≠ 0 := by positivity
  have h3 : 2 * n + 2 ≠ 0 := by positivity
  field_simp
  ring

/-- One circle: γ − 1 = −1/2. -/
theorem one_circle_gamma : (1 : ℝ) / (1 + 1) - 1 = -1 / 2 := by norm_num

end FusedAlgebra.Topics.NewForce
