/-
  Fused Hyperoperation Algebra: cell biology (drug and toxicant combinations, cell kill).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/CellBiology.lean

  1. Bliss independence (survivals multiply) is the "or" law on effects, dial −1, and the plain
     circle operation (called mult-add before 1.0) on effect odds o = E/(1 − E).
  2. Chou–Talalay's combination index v₁ + v₂ + α·v₁v₂ is the weighted circle operation v₁ ◯_α v₂.
     α = 0 is Loewe additivity, α = −1 is the Highest Single Agent model, and the wall
     v = −1/α is where one drug makes the other irrelevant.
  3. Multi-target cell kill: survival is a weighted power at dial −1.
-/
import FusedAlgebra.Base

namespace CellBiology
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

theorem star_def (a b : R) : a ◯ b = a + b + a * b := rfl

/-! ## 1. Bliss independence -/

/-- Survivals multiply, so effects combine by the "or" law, dial −1. -/
theorem bliss_or (a b : R) : 1 - (1 - a) * (1 - b) = circK (-1) a b := by
  unfold circK; grind

/-! ## 2. The Chou–Talalay dial -/

/-- The combination index with interaction α is the weighted circle operation of the dose fractions. -/
def ci (α v₁ v₂ : R) : R := v₁ + v₂ + α * v₁ * v₂

theorem ci_is_circK (α v₁ v₂ : R) : ci α v₁ v₂ = circK α v₁ v₂ := by unfold ci circK; grind

/-- α = 0: Loewe dose additivity. -/
theorem ci_loewe (v₁ v₂ : R) : ci 0 v₁ v₂ = v₁ + v₂ := by unfold ci; grind

/-- α = −1: the index sits on the additivity line exactly when one factor 1 − vᵢ vanishes,
    i.e. one drug alone already gives the effect. This is the Highest Single Agent model. -/
theorem ci_hsa (v₁ v₂ : R) : ci (-1) v₁ v₂ - 1 = -((1 - v₁) * (1 - v₂)) := by
  unfold ci; grind

/-- α = 1: Chou's "mutually nonexclusive" index factors as (1 + v₁)(1 + v₂) − 1. -/
theorem ci_nonexclusive (v₁ v₂ : R) : ci 1 v₁ v₂ = v₁ ◯ v₂ := by unfold ci circ; grind

/-- Moving one drug moves the index by the dial lift of the other. -/
theorem ci_shift (α v₁ v₂ v₂' : R) : ci α v₁ v₂ - ci α v₁ v₂' = (v₂ - v₂') * (1 + α * v₁) := by
  unfold ci; grind

/-! ## 3. Multi-target kill -/

/-- A cell with n targets survives if any target survives: survival 1 − (1 − s)ⁿ is the
    dial −1 power of the per-target survival s. -/
theorem multitarget (s : R) (n : Nat) : 1 - (1 - s) ^ n = wpow (-1) s n := by
  have h := wpow_lift (-1 : R) s n
  have e : (1 : R) + -1 * s = 1 - s := by grind
  rw [e] at h; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- Effect odds. -/
def odds (E : K) : K := E / (1 - E)

theorem one_add_odds (E : K) (h : 1 - E ≠ 0) : 1 + odds E = 1 / (1 - E) := by
  unfold odds; grind

/-- **Bliss is the plain circle operation on effect odds**, for any shape of dose–response curve. -/
theorem bliss_odds (a b : K) (ha : 1 - a ≠ 0) (hb : 1 - b ≠ 0) :
    odds (1 - (1 - a) * (1 - b)) = odds a ◯ odds b := by
  unfold odds circ
  have hab : 1 - (1 - (1 - a) * (1 - b)) = (1 - a) * (1 - b) := by grind
  rw [hab]
  have hp : (1 - a) * (1 - b) ≠ 0 := by
    intro e
    have : (1 - a)⁻¹ * ((1 - a) * (1 - b)) = 1 - b := by
      rw [← Lean.Grind.Semiring.mul_assoc, Lean.Grind.Field.inv_mul_cancel ha]; grind
    rw [e] at this; apply hb; grind
  grind

/-- The Highest Single Agent reading of α = −1, over a field: on the additivity line one of
    the two drugs is at its own iso-effective dose. -/
theorem ci_hsa_field (v₁ v₂ : K) (h : ci (-1) v₁ v₂ = 1) : v₁ = 1 ∨ v₂ = 1 := by
  have e : (1 - v₁) * (1 - v₂) = 0 := by have := ci_hsa v₁ v₂; rw [h] at this; grind
  by_cases h1 : 1 - v₁ = 0
  · left; grind
  · right
    have : (1 - v₁)⁻¹ * ((1 - v₁) * (1 - v₂)) = 1 - v₂ := by
      rw [← Lean.Grind.Semiring.mul_assoc, Lean.Grind.Field.inv_mul_cancel h1]; grind
    rw [e] at this; grind

/-- **The wall makes one drug irrelevant.** At v₂ = −1/α the index is −1/α whatever v₁ is.
    For α = −1 the wall is v₂ = 1: drug 2 alone already gives the effect. -/
theorem ci_wall (α v₁ : K) (hα : α ≠ 0) : ci α v₁ (-1 / α) = -1 / α := by
  rw [ci_is_circK, circK_wall α v₁ hα]

/-- Loewe with Hill slope 1: vᵢ = uᵢ/o, so v₁ + v₂ = 1 exactly when the effect odds o are
    u₁ + u₂, odds adding, dial 0. -/
theorem loewe_slope_one (u₁ u₂ o : K) (ho : o ≠ 0) :
    u₁ / o + u₂ / o = 1 ↔ o = u₁ + u₂ := by
  constructor
  · intro h
    have : (u₁ / o + u₂ / o) * o = u₁ + u₂ := by grind
    rw [h] at this; grind
  · intro h; rw [h]; rw [h] at ho; grind

end Field

end CellBiology
