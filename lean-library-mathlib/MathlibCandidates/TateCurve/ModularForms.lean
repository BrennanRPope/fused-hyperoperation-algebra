module

public import MathlibCandidates.TateCurve
public import Mathlib.NumberTheory.ModularForms.LevelOne.GradedRing

/-!
# The Tate curve and level one modular forms  (DRAFT, candidate for Mathlib)

The formal Eisenstein series `TateCurve.E₄`, `TateCurve.E₆` are the q-expansions of Mathlib's
normalised Eisenstein series `ModularForm.E₄`, `ModularForm.E₆`, and the discriminant of the
Tate curve over `ℂ⟦q⟧` is the q-expansion of the modular discriminant `Δ`.

## Main statements

* `TateCurve.qExpansion_E₄ : qExpansion 1 ModularForm.E₄ = TateCurve.E₄ ℂ`
* `TateCurve.qExpansion_E₆ : qExpansion 1 ModularForm.E₆ = TateCurve.E₆ ℂ`
* `TateCurve.qExpansion_discriminant : qExpansion 1 CuspForm.discriminant = (curve ℂ).Δ`
* `TateCurve.curve_c₄_eq_qExpansion`, `TateCurve.curve_c₆_eq_qExpansion`
-/

@[expose] public section

open PowerSeries ModularForm EisensteinSeries UpperHalfPlane

namespace TateCurve

theorem qExpansion_E₄ : qExpansion 1 ModularForm.E₄ = TateCurve.E₄ ℂ := by
  ext m
  rw [E_qExpansion_coeff (by norm_num) ⟨2, rfl⟩, show bernoulli 4 = -1 / 30 by decide +kernel,
    TateCurve.E₄, map_add, coeff_numeral_mul, coeff_s]
  rcases eq_or_ne m 0 with rfl | hm
  · simp
  · simp only [hm, ↓reduceIte, coeff_one]
    push_cast
    ring

theorem qExpansion_E₆ : qExpansion 1 ModularForm.E₆ = TateCurve.E₆ ℂ := by
  ext m
  rw [E_qExpansion_coeff (by norm_num) ⟨3, rfl⟩, show bernoulli 6 = 1 / 42 by decide +kernel,
    TateCurve.E₆, map_sub, coeff_numeral_mul, coeff_s]
  rcases eq_or_ne m 0 with rfl | hm
  · simp
  · simp only [hm, ↓reduceIte, coeff_one]
    push_cast
    ring

theorem curve_c₄_eq_qExpansion : (curve ℂ).c₄ = qExpansion 1 ModularForm.E₄ := by
  rw [curve_c₄, qExpansion_E₄]

theorem curve_c₆_eq_qExpansion : (curve ℂ).c₆ = -qExpansion 1 ModularForm.E₆ := by
  rw [curve_c₆, qExpansion_E₆]

theorem qExpansion_discriminant : qExpansion 1 CuspForm.discriminant = (curve ℂ).Δ := by
  have hfun : (CuspForm.discriminant : UpperHalfPlane → ℂ) =
      ((1 / 1728 : ℂ) • (ModularForm.mcast (by decide) (ModularForm.E₄.pow 3) -
        ModularForm.mcast (by decide) (ModularForm.E₆.pow 2) : ModularForm _ 12) :
          UpperHalfPlane → ℂ) := by
    ext z
    simp only [CuspForm.coe_discriminant, FunLike.coe_sub, coe_mcast, ModularForm.coe_pow,
      Pi.smul_apply, Pi.sub_apply, Pi.pow_apply, smul_eq_mul]
    rw [discriminant_eq_E₄_cube_sub_E₆_sq]
    ring
  have h1728 : (1728 : ℂ⟦X⟧) * qExpansion 1 CuspForm.discriminant =
      TateCurve.E₄ ℂ ^ 3 - TateCurve.E₆ ℂ ^ 2 := by
    rw [hfun, ModularForm.qExpansion_smul one_pos one_mem_strictPeriods_SL]
    simp only [FunLike.coe_sub, coe_mcast]
    rw [ModularForm.qExpansion_sub one_pos one_mem_strictPeriods_SL (ModularForm.E₄.pow 3)
        (ModularForm.E₆.pow 2),
      ModularForm.qExpansion_pow one_pos one_mem_strictPeriods_SL,
      ModularForm.qExpansion_pow one_pos one_mem_strictPeriods_SL, qExpansion_E₄, qExpansion_E₆,
      smul_eq_C_mul, ← mul_assoc, ← map_ofNat (C (R := ℂ)) 1728, ← map_mul]
    norm_num
  have h0 : (1728 : ℂ⟦X⟧) ≠ 0 := by
    rw [← map_ofNat (C (R := ℂ)) 1728]
    intro h
    simpa using congrArg (constantCoeff (R := ℂ)) h
  exact mul_left_cancel₀ h0 (h1728.trans (curve_c_relation ℂ).symm)

end TateCurve
