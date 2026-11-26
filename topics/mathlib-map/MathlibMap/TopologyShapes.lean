/-
  Topological shapes (DRAFT): holes and twists as whole-number picks.

  Every convex shape on the shape-dials page is one topological type (a ball). The other closed
  surfaces are polygons with their edges glued in pairs by a word; χ = (corner classes) − (edge pairs) + 1.

  Proved here, with no `sorry`:
    * `chi_concat`: writing two gluing words one after the other (connected sum) gives
      χ = χ₁ + χ₂ − 2, because the corner classes share one point and the edge pairs add.
    * `euler_genus_add`: so the Euler genus 2 − χ (holes counted twice, or twists) adds.
    * `curvature_pos_iff`, `curvature_zero_iff`, `curvature_neg_iff`: Gauss–Bonnet on the glued
      polygon, κ · area = 2π χ with area > 0, forces the sign of the curvature to be the sign of χ.
    * `chi_product`: χ is the face polynomial read at t = −1, and face polynomials multiply on
      products, so χ multiplies (the Künneth formula's shadow on counts).
    * `one_sub_chi_mul`: hence 1 − χ composes as x + y − x y on products.
  Numerical companions: shape-dials/topo_checks.py (9 checks).
-/
import Mathlib

namespace MathlibMap.TopologyShapes

/-- Gluing words concatenated: corner classes `V₁ + V₂ − 1`, edge pairs `E₁ + E₂`, one face. -/
theorem chi_concat (V₁ V₂ E₁ E₂ : ℤ) :
    (V₁ + V₂ - 1) - (E₁ + E₂) + 1 = (V₁ - E₁ + 1) + (V₂ - E₂ + 1) - 2 := by
  ring

/-- Euler genus `2 − χ` adds under connected sum. -/
theorem euler_genus_add (χ₁ χ₂ : ℤ) : 2 - (χ₁ + χ₂ - 2) = (2 - χ₁) + (2 - χ₂) := by
  ring

/-- Gauss–Bonnet: positive curvature exactly when χ > 0. -/
theorem curvature_pos_iff (κ A χ : ℝ) (hA : 0 < A) (h : κ * A = 2 * Real.pi * χ) :
    0 < κ ↔ 0 < χ := by
  have hπ := Real.pi_pos
  constructor
  · intro hk
    have h1 : 0 < κ * A := mul_pos hk hA
    by_contra hc
    push_neg at hc
    nlinarith [mul_nonneg hπ.le (neg_nonneg.mpr hc)]
  · intro hc
    have h1 : 0 < 2 * Real.pi * χ := mul_pos (mul_pos two_pos hπ) hc
    by_contra hk
    push_neg at hk
    nlinarith [mul_nonneg (neg_nonneg.mpr hk) hA.le]

/-- Gauss–Bonnet: flat exactly when χ = 0. -/
theorem curvature_zero_iff (κ A χ : ℝ) (hA : 0 < A) (h : κ * A = 2 * Real.pi * χ) :
    κ = 0 ↔ χ = 0 := by
  constructor
  · intro hk
    rw [hk, zero_mul] at h
    rcases mul_eq_zero.mp h.symm with h2 | h2
    · exact absurd h2 (by positivity)
    · exact h2
  · intro hc
    rw [hc, mul_zero] at h
    rcases mul_eq_zero.mp h with h2 | h2
    · exact h2
    · exact absurd h2 hA.ne'

/-- Gauss–Bonnet: negative curvature exactly when χ < 0. -/
theorem curvature_neg_iff (κ A χ : ℝ) (hA : 0 < A) (h : κ * A = 2 * Real.pi * χ) :
    κ < 0 ↔ χ < 0 := by
  have hπ := Real.pi_pos
  constructor
  · intro hk
    have h1 : κ * A < 0 := mul_neg_of_neg_of_pos hk hA
    by_contra hc
    push_neg at hc
    nlinarith [mul_nonneg hπ.le hc]
  · intro hc
    have h1 : 2 * Real.pi * χ < 0 := mul_neg_of_pos_of_neg (mul_pos two_pos hπ) hc
    by_contra hk
    push_neg at hk
    nlinarith [mul_nonneg hk hA.le]

/-- χ is the face polynomial at t = −1; face polynomials multiply on products, so χ multiplies. -/
theorem chi_product (p q : Polynomial ℤ) :
    (p * q).eval (-1) = p.eval (-1) * q.eval (-1) :=
  Polynomial.eval_mul

/-- So `1 − χ` composes on products as `x + y − x y`. -/
theorem one_sub_chi_mul (a b : ℤ) : 1 - a * b = (1 - a) + (1 - b) - (1 - a) * (1 - b) := by
  ring

end MathlibMap.TopologyShapes
