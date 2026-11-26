import Mathlib.Algebra.Algebra.Spectrum.Quasispectrum
import Mathlib.Tactic.LinearCombination
import FusedAlgebra.Base

/-!
# Our law at k = 1 is Mathlib's circle operation

Mathlib's `PreQuasiregular R` wraps a ring `R` with the monoid law `x * y = y + x + x * y` and
identity `0` (the "circle operation" of ring theory). The core library's law is
`CircleOp.circK k a b = a + b + k * a * b`. This file proves:

* `circ_eq_quasiregular_mul`, `circK_one_eq_quasiregular_mul`: at k = 1 the two laws are the same operation.
* `circK_scale_quasiregular_mul`: at any k, multiplying by k carries our law to the circle operation.
* `Dial R k`: our law at dial k as a monoid, and `Dial.equivCircle`: for a unit k it is
  isomorphic to Mathlib's `PreQuasiregular R`. So every dial with k invertible is the same
  monoid as Mathlib's, and everything Mathlib proves about the circle monoid transfers.
* `wayBack_iff_isQuasiregular`: a has a way back at dial k (k a unit) exactly when `k * a` is
  Mathlib-quasiregular.
-/

namespace MathlibMap

open PreQuasiregular

section Comm
variable {R : Type} [CommRing R]

set_option linter.unusedSectionVars false

/-- Our circle operation `circ` (dial k = 1) is Mathlib's circle operation (on any commutative ring). -/
theorem circ_eq_quasiregular_mul (a b : R) :
    CircleOp.circ a b = (equiv a * equiv b : PreQuasiregular R).val := by
  change a + b + a * b = b + a + a * b; ring

/-- Our law at dial k = 1 is Mathlib's circle operation. -/
theorem circK_one_eq_quasiregular_mul (a b : R) :
    CircleOp.circK 1 a b = (equiv a * equiv b : PreQuasiregular R).val := by
  change a + b + 1 * a * b = b + a + a * b; ring

/-- At any dial, multiplying by k turns our law into the circle operation. -/
theorem circK_scale_quasiregular_mul (k a b : R) :
    equiv (k * CircleOp.circK k a b) = (equiv (k * a) * equiv (k * b) : PreQuasiregular R) := by
  apply equiv.symm.injective
  simp only [Equiv.symm_apply_apply]
  show k * CircleOp.circK k a b = (equiv (k * a) * equiv (k * b) : PreQuasiregular R).val
  change k * (a + b + k * a * b) = k * b + k * a + k * a * (k * b); ring

/-- Our law at dial `k`, as a type with its own monoid structure (identity `0`). -/
structure Dial (R : Type) (k : R) where
  /-- The underlying ring element. -/
  val : R

namespace Dial
variable {k : R}

theorem ext_val {x y : Dial R k} (h : x.val = y.val) : x = y := by
  cases x; cases y; cases h; rfl

instance : Monoid (Dial R k) where
  one := ⟨0⟩
  mul x y := ⟨CircleOp.circK k x.val y.val⟩
  mul_assoc x y z := ext_val (CircleOp.circK_assoc k x.val y.val z.val)
  one_mul x := ext_val (by show CircleOp.circK k 0 x.val = x.val; rw [CircleOp.circK_def]; ring)
  mul_one x := ext_val (CircleOp.circK_zero k x.val)

@[simp] theorem val_mul (x y : Dial R k) : (x * y).val = CircleOp.circK k x.val y.val := rfl
@[simp] theorem val_one : (1 : Dial R k).val = 0 := rfl

instance : CommMonoid (Dial R k) where
  mul_comm x y := ext_val (CircleOp.circK_comm k x.val y.val)

/-- For a unit `k`, our dial-`k` monoid is isomorphic to Mathlib's circle monoid
`PreQuasiregular R`, by `a ↦ k * a`. -/
def equivCircle (k : Rˣ) : Dial R (k : R) ≃* PreQuasiregular R where
  toFun x := equiv ((k : R) * x.val)
  invFun y := ⟨((k⁻¹ : Rˣ) : R) * y.val⟩
  left_inv x := ext_val (by simp [← mul_assoc])
  right_inv y := by
    apply equiv.symm.injective
    simp [← mul_assoc]
  map_mul' x y := circK_scale_quasiregular_mul (k : R) x.val y.val

end Dial

/-- A two-sided way back at dial `k`, read in Mathlib's language: for a unit `k`, `a` has
`b` with `a ∘ₖ b = 0` exactly when `k * a` is quasiregular (`IsQuasiregular`). -/
theorem wayBack_iff_isQuasiregular (k : Rˣ) (a : R) :
    (∃ b, CircleOp.circK (k : R) a b = 0) ↔ IsQuasiregular ((k : R) * a) := by
  rw [isQuasiregular_iff]
  constructor
  · rintro ⟨b, hb⟩
    refine ⟨(k : R) * b, ?_, ?_⟩
    · have := congrArg ((k : R) * ·) hb
      simp only [CircleOp.circK_def, mul_zero] at this
      linear_combination this
    · have := congrArg ((k : R) * ·) hb
      simp only [CircleOp.circK_def, mul_zero] at this
      linear_combination this
  · rintro ⟨y, hy, -⟩
    refine ⟨((k⁻¹ : Rˣ) : R) * y, ?_⟩
    have hk : ((k⁻¹ : Rˣ) : R) * (k : R) = 1 := by simp
    rw [CircleOp.circK_def]
    linear_combination ((k⁻¹ : Rˣ) : R) * hy - a * hk

end Comm

end MathlibMap
