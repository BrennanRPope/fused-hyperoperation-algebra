/-
  Fused Hyperoperation Algebra: ecology (predator interference, diversity, stressors).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Ecology.lean

  1. A predator that searches, handles, and is drawn into interference from searching (rate y)
     and, with relative propensity k, from handling, spends a share 1/(1 + x ◯ₖ y) of its
     time searching, x = αNh, y = γPτ. k = 0 is Beddington–DeAngelis, k = 1 Crowley–Martin.
  2. Tsallis entropy is the q-logarithm of the Hill number; independent communities combine
     by the dial k = 1 − q (proved here for q = 2, where ln₂(D) = 1 − 1/D).
  3. Multiple-stressor nulls: additive is dial 0, multiplicative is the "or" law.
  LIBRARY NOTE: first compiled in the shared library. 3 exact twins removed (star_def,
  stressors, colonisation); each is marked where it stood. No proof step or statement
  changed.
-/
import FusedAlgebra.Base

namespace Ecology
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

/-- The time-budget denominator: 1 + x + y + k·x·y. -/
def budget (k x y : R) : R := 1 + circK k x y

/-- The four states' weights (search 1, handling x, interfering from search y, interfering
    from handling k·x·y) sum to the budget. -/
theorem budget_states (k x y : R) : 1 + x + y + k * x * y = budget k x y := by
  unfold budget circK; grind

/-- Beddington–DeAngelis is dial 0. -/
theorem bd (x y : R) : budget 0 x y = 1 + x + y := by unfold budget circK; grind

/-- Crowley–Martin is dial 1, the product of the two saturations. -/
theorem cm (x y : R) : budget 1 x y = (1 + x) * (1 + y) := by unfold budget circK; grind

/-- The dial's excess over BD is exactly the handling–interference corner k·x·y. -/
theorem corner (k x y : R) : budget k x y - budget 0 x y = k * x * y := by
  unfold budget circK; grind

-- LIBRARY: stressors (Multiple stressors, independent action: proportional losses combine by the "or" law) was an exact twin of CellBiology.bliss_or and is removed.

-- LIBRARY: colonisation (Island colonisation from n sources) was an exact twin of Metabolism.limitation_pow and is removed.

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- Balance of the predator chain: if the handling share is x times the search share and the
    other two states are y and k·x·y times it, and the four shares sum to 1, the search share
    is 1/(1 + x ◯ₖ y). -/
theorem search_share (k x y s : K) (h : s * (1 + x + y + k * x * y) = 1)
    (hb : budget k x y ≠ 0) : s = 1 / budget k x y := by
  have e := budget_states k x y
  rw [e] at h
  have i := Lean.Grind.Field.mul_inv_cancel hb
  rw [Lean.Grind.Field.div_eq_mul_inv]
  have : s = s * (budget k x y * (budget k x y)⁻¹) := by rw [i]; grind
  rw [this, ← Lean.Grind.Semiring.mul_assoc, h]

/-- q = 2 diversity: ln₂(D) = 1 − 1/D. Hill numbers multiply, so these combine by dial −1. -/
theorem tsallis_two (D₁ D₂ : K) :
    1 - 1 / (D₁ * D₂) = circK (-1) (1 - 1 / D₁) (1 - 1 / D₂) := by
  unfold circK
  simp only [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Field.inv_mul]
  grind

/-- For q = 2 the wall −1/k is 1: an entropy at the ceiling absorbs any other community. -/
theorem tsallis_wall (S : K) : circK (-1) 1 S = 1 := by unfold circK; grind

end Field

end Ecology
