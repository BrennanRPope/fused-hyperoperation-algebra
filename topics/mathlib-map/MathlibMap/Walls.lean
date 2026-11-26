import Mathlib.Algebra.Algebra.Spectrum.Quasispectrum
import Mathlib.Tactic.LinearCombination
import FusedAlgebra.Base

/-!
# Where an element hits the wall, across the whole dial, is its quasispectrum

Our law needs no commutativity to make sense: in any non-unital algebra `A` over a field `𝕜`,
`a ∘ₖ b = a + b + k • (a * b)` (`dialAdd`). In a commutative ring it is the core library's
`CircleOp.circK` (`dialAdd_self`).

`HasWayBack k a` says `a` has a two-sided way back at dial `k`. The main theorem,
`not_hasWayBack_iff_mem_quasispectrum`, says that for `k ≠ 0`:

  `a` loses its way back at dial `k`  ⟺  `-1/k` is in Mathlib's `quasispectrum 𝕜 a`.

So the set of walls an element meets as the dial turns is exactly its quasispectrum, flipped by
`λ ↦ -1/λ`. With Mathlib's `quasispectrum_eq_spectrum_union_zero`, in a unital algebra the walls
are the spectrum (`not_hasWayBack_iff_mem_spectrum`). Dial `0` is never a wall
(`hasWayBack_zero`): adding has a way back for everything.
-/

namespace MathlibMap

variable {𝕜 A : Type*} [Field 𝕜] [NonUnitalRing A] [Module 𝕜 A] [IsScalarTower 𝕜 A A]
  [SMulCommClass 𝕜 A A]

set_option linter.unusedSectionVars false

/-- Our law at dial `k` in any non-unital `𝕜`-algebra: `a ∘ₖ b = a + b + k • (a * b)`. -/
def dialAdd (k : 𝕜) (a b : A) : A := a + b + k • (a * b)

/-- `a` has a two-sided way back at dial `k`. -/
def HasWayBack (k : 𝕜) (a : A) : Prop := ∃ b, dialAdd k a b = 0 ∧ dialAdd k b a = 0

/-- In a field (or any commutative ring seen as an algebra over itself) `dialAdd` is the core
library's `CircleOp.circK`. -/
theorem dialAdd_self {K : Type} [Field K] (k a b : K) : dialAdd k a b = CircleOp.circK k a b := by
  rw [dialAdd, CircleOp.circK_def, smul_eq_mul]; ring

/-- Dial zero is never a wall. -/
theorem hasWayBack_zero (a : A) : HasWayBack (0 : 𝕜) a :=
  ⟨-a, by simp [dialAdd], by simp [dialAdd]⟩

/-- Scaling by `k` turns our law into Mathlib's circle operation (the `PreQuasiregular` law). -/
theorem smul_dialAdd (k : 𝕜) (a b : A) :
    k • dialAdd k a b = k • b + k • a + (k • a) * (k • b) := by
  rw [dialAdd, smul_add, smul_add, smul_mul_smul_comm, add_comm (k • a) (k • b), mul_smul]

/-- For `k ≠ 0`, a way back at dial `k` is Mathlib-quasiregularity of `k • a`. -/
theorem hasWayBack_iff_isQuasiregular {k : 𝕜} (hk : k ≠ 0) (a : A) :
    HasWayBack k a ↔ IsQuasiregular (k • a) := by
  rw [isQuasiregular_iff]
  constructor
  · rintro ⟨b, h₁, h₂⟩
    refine ⟨k • b, ?_, ?_⟩
    · rw [← smul_dialAdd, h₁, smul_zero]
    · have := congrArg (k • ·) h₂
      simp only [smul_zero] at this
      rw [smul_dialAdd] at this
      rw [← this]
  · rintro ⟨y, h₁, h₂⟩
    refine ⟨k⁻¹ • y, ?_, ?_⟩
    · apply smul_right_injective A hk
      simp only [smul_zero]
      rw [smul_dialAdd, smul_inv_smul₀ hk]; exact h₁
    · apply smul_right_injective A hk
      simp only [smul_zero]
      rw [smul_dialAdd, smul_inv_smul₀ hk, ← h₂]

/-- **The walls of the dial are the quasispectrum.** For `k ≠ 0`, `a` has no way back at dial
`k` exactly when `-1/k` lies in Mathlib's `quasispectrum 𝕜 a`. -/
theorem not_hasWayBack_iff_mem_quasispectrum {k : 𝕜} (hk : k ≠ 0) (a : A) :
    ¬ HasWayBack k a ↔ -k⁻¹ ∈ quasispectrum 𝕜 a := by
  have hu : -k⁻¹ ≠ 0 := neg_ne_zero.mpr (inv_ne_zero hk)
  have key := quasispectrum.units_mem_iff 𝕜 (Units.mk0 (-k⁻¹) hu) a
  simp only [Units.val_mk0] at key
  rw [key, hasWayBack_iff_isQuasiregular hk]
  have : -(((Units.mk0 (-k⁻¹) hu)⁻¹ : 𝕜ˣ) • a) = k • a := by
    rw [Units.smul_def, Units.val_inv_eq_inv_val, Units.val_mk0, inv_neg, inv_inv, neg_smul,
      neg_neg]
  rw [this]

/-- Read the other way: a nonzero `λ` is in the quasispectrum exactly when `a` meets a wall at
dial `-1/λ`. -/
theorem mem_quasispectrum_iff_wall {r : 𝕜} (hr : r ≠ 0) (a : A) :
    r ∈ quasispectrum 𝕜 a ↔ ¬ HasWayBack (-r⁻¹) a := by
  have hk : -r⁻¹ ≠ 0 := neg_ne_zero.mpr (inv_ne_zero hr)
  rw [not_hasWayBack_iff_mem_quasispectrum hk, inv_neg, inv_inv, neg_neg]

/-- In a unital algebra the walls are the spectrum (Mathlib's `spectrum`). -/
theorem not_hasWayBack_iff_mem_spectrum {B : Type*} [Ring B] [Algebra 𝕜 B] {k : 𝕜} (hk : k ≠ 0)
    (a : B) : ¬ HasWayBack k a ↔ -k⁻¹ ∈ spectrum 𝕜 a := by
  rw [not_hasWayBack_iff_mem_quasispectrum hk, quasispectrum_eq_spectrum_union_zero]
  simp [hk]

end MathlibMap
