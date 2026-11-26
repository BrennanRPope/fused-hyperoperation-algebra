/-
DRAFT (Nov 26, 2026). Does the theory predict an unknown force? The algebra behind
topics/force-reach/new-force/README.md, compiled in the shared library on Nov 26 as FusedAlgebra.Topics.NewForce.

The session that wrote this stated it against Mathlib over ℝ (topics/force-reach/new-force/lean/NewForce.lean,
fourteen statements). Here all fourteen are proved in core Lean over any field of characteristic 0 or, for the two
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

end CircleOp.NewForce
