import Mathlib.RingTheory.Jacobson.Ideal
import Mathlib.RingTheory.Nilpotent.Basic
import Mathlib.Tactic.LinearCombination
import FusedAlgebra.Base

/-!
# The Jacobson radical is what never meets a wall

On any commutative ring (no field needed), the core library's way back at dial `k`
(`∃ b, a ∘ₖ b = 0`) is exactly "`1 + k a` is a unit" (`wayBack_iff_isUnit`). With Mathlib's
`Ring.mem_jacobson_iff` this gives:

* `mem_jacobson_iff_wayBack`: `x` is in the Jacobson radical exactly when `x` has a way back at
  every dial `k ∈ R`. The radical is the set of elements that never meet a wall.
* `wayBack_of_isNilpotent`: nilpotent elements have a way back at every dial (Mathlib's
  `IsNilpotent.isUnit_one_add`): a finite Neumann sum.
-/

namespace MathlibMap

variable {R : Type} [CommRing R]

/-- A way back at dial `k`, on any commutative ring, is the unit `1 + k * a`. -/
theorem wayBack_iff_isUnit (k a : R) : (∃ b, CircleOp.circK k a b = 0) ↔ IsUnit (1 + k * a) := by
  constructor
  · rintro ⟨b, hb⟩
    have h := CircleOp.circK_lift k a b
    rw [hb, mul_zero, add_zero] at h
    exact IsUnit.of_mul_eq_one _ h.symm
  · intro hu
    obtain ⟨c, hc⟩ := hu.exists_right_inv
    refine ⟨-a * c, ?_⟩
    rw [CircleOp.circK_def]
    linear_combination (-a) * hc

/-- **The Jacobson radical is the set of elements that never meet a wall**: `x` lies in it
exactly when `x` has a way back at every dial `k`. -/
theorem mem_jacobson_iff_wayBack (x : R) :
    x ∈ Ring.jacobson R ↔ ∀ k : R, ∃ b, CircleOp.circK k x b = 0 := by
  rw [Ring.mem_jacobson_iff]
  refine forall_congr' fun k => ?_
  rw [wayBack_iff_isUnit, mul_comm x k, add_comm]

/-- Nilpotent elements have a way back at every dial. -/
theorem wayBack_of_isNilpotent {x : R} (hx : IsNilpotent x) (k : R) :
    ∃ b, CircleOp.circK k x b = 0 :=
  (wayBack_iff_isUnit k x).mpr ((Commute.all k x).isNilpotent_mul_left hx).isUnit_one_add

/-- In a field the only element that never meets a wall is `0`: the Jacobson radical of a
field is zero. -/
theorem wayBack_every_dial_field {K : Type} [Field K] {x : K}
    (h : ∀ k : K, ∃ b, CircleOp.circK k x b = 0) : x = 0 := by
  by_contra hx
  obtain ⟨b, hb⟩ := h (-x⁻¹)
  have := (wayBack_iff_isUnit (-x⁻¹) x).mp ⟨b, hb⟩
  rw [neg_mul, inv_mul_cancel₀ hx, add_neg_cancel] at this
  exact not_isUnit_zero this

end MathlibMap
