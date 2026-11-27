/-
DRAFT (Nov 26, 2026). Does the theory predict an unknown force? The algebra behind
topics/force-reach/new-force/README.md, compiled in the shared library on Nov 26 as FusedAlgebra.Topics.NewForce.

The session that wrote this stated it against Mathlib over ℝ (topics/force-reach/new-force/lean/NewForce.lean,
fourteen statements; twenty from Nov 27). Here all of them are proved in core Lean over any field of characteristic 0 or, for the two
flavour bounds and γ from the scalar strength, any ordered field. Statements and hypotheses are as written.
-/
import FusedAlgebra.Base

open Lean.Grind

namespace CircleOp.NewForce

section Field
variable {K : Type} [Field K] [IsCharP K 0]

theorem tower_sum (q : K) (hq : q ≠ 1) :
    1 + 2 * q / (1 - q) = (1 + q) / (1 - q) := by
  have h : (1 - q) ≠ 0 := by intro h; apply hq; grind
  grind

theorem tower_leading_yukawa (q : K) (hq : q ≠ 1) :
    (1 + q) / (1 - q) - 1 - 2 * q = 2 * q ^ 2 / (1 - q) := by
  have h : (1 - q) ≠ 0 := by intro h; apply hq; grind
  grind

theorem dual_exponent_mond : (3 / 2 : K) / (3 / 2 - 1) = 3 := by grind
theorem holder_pair : (1 : K) / (3 / 2) + 1 / 3 = 1 := by grind
omit [IsCharP K 0] in
theorem vector_over_tensor : (2 : K) / (1 / 2) = 4 := by grind
theorem scalar_over_tensor : (3 / 2 : K) / (1 / 2) = 3 := by grind
theorem like_charges_cancel : (1 / 2 : K) + 3 / 2 - 2 = 0 := by grind
theorem opposite_charges_add : (1 / 2 : K) + 3 / 2 + 2 = 4 := by grind
theorem klein_strength (α : K) (hα : α ≠ 0) : 4 / (4 / α) = α := by grind
theorem one_circle_third : (1 : K) / (1 + 2) = 1 / 3 := by grind
theorem one_circle_gamma : (1 : K) / (1 + 1) - 1 = -1 / 2 := by grind
end Field

section Ordered
variable {K : Type} [LE K] [LT K] [Field K] [Std.IsLinearOrder K] [Std.LawfulOrderLT K] [OrderedRing K]

theorem su3_flavour_bound :
    (0 : K) < 11 * 3 / 3 - 2 * 16 / 3 ∧ 11 * 3 / 3 - 2 * 17 / 3 < (0 : K) := by
  constructor <;> grind

theorem su2_flavour_bound :
    (0 : K) < 11 * 2 / 3 - 2 * 10 / 3 ∧ 11 * 2 / 3 - 2 * 11 / 3 = (0 : K) := by
  constructor <;> grind

theorem gamma_from_strength (n : K) (hn : 0 ≤ n) :
    (1 - n / (n + 2)) / (1 + n / (n + 2)) = 1 / (n + 1) := by
  have h2 : n + 2 ≠ 0 := by grind
  have h1 : n + 1 ≠ 0 := by grind
  have h3 : 2 * n + 2 ≠ 0 := by grind
  grind
end Ordered

section Hold
variable {K : Type} [Field K] [IsCharP K 0]

/-! ### Part H: holding the circle (vacuum energy plus a 5D vacuum energy); added Nov 27 with sections 7 and 8 -/

/-- With V = a/R + C R*²/R⁶ and a = −6C/R*³, R* is a turning point: V′(R*) = 0. -/
theorem hold_turning_point (C R : K) (hR : R ≠ 0) :
    -(-6 * C / R ^ 3) / R ^ 2 - 6 * C * R ^ 2 / R ^ 7 = 0 := by
  grind

/-- Its curvature there: V″(R*) = 30C/R*⁶ (a minimum when C > 0). -/
theorem hold_curvature (C R : K) (hR : R ≠ 0) :
    2 * (-6 * C / R ^ 3) / R ^ 3 + 42 * C * R ^ 2 / R ^ 8 = 30 * C / R ^ 6 := by
  grind

/-- The vacuum energy left over: V(R*) = −5C/R*⁴. -/
theorem hold_leftover (C R : K) (hR : R ≠ 0) :
    (-6 * C / R ^ 3) / R + C * R ^ 2 / R ^ 6 = -5 * C / R ^ 4 := by
  grind

/-- Radion mass: (2/3)·30 = 20, so m² = 20C/(M̄²R⁴). -/
theorem hold_mass_factor : (2 / 3 : K) * 30 = 20 := by grind

/-! ### Part S: strength from shape and material -/

/-- Sphere: ⟨|K|²⟩ = 2R²/3, so the shape factor 4R²/⟨|K|²⟩ is 6 (a circle gives 4). -/
theorem sphere_shape_factor : (4 : K) / (2 / 3) = 6 := by grind

omit [IsCharP K 0] in
/-- Strong force with 5 flavours: b₀ = 11·3/3 − 2·5/3 = 23/3. -/
theorem strong_b0_five : (11 * 3 / 3 : K) - 2 * 5 / 3 = 23 / 3 := by grind
end Hold

end CircleOp.NewForce
