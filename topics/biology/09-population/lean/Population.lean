/-
  Fused Hyperoperation Algebra: population biology and epidemics.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Population.lean

  1. Reed–Frost: the chance a susceptible is infected by I infectives is wpow(−1, p, I);
     independent interventions combine by the "or" law; the herd-immunity threshold
     1 − 1/R₀ is −inv(R₀ − 1).
  2. Beverton–Holt: two seasons compose into one, with growth excesses combining by ◯;
     the room below carrying capacity shrinks by 1/R each generation.
  3. Hassell: (1 + aN)^b = 1 + mpow(aN, b); b = 1 is Beverton–Holt.
  4. Survivorship and Levins colonisation are "or" folds.
  LIBRARY NOTE: first compiled in the shared library. 7 exact twins removed (star_def,
  reed_frost, seasons_growth, hassell, survival, levins, herd_threshold); each is marked
  where it stood. No proof step or statement changed.
-/
import FusedAlgebra.Base

namespace Population
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

-- LIBRARY: reed_frost (Reed–Frost: infection chance from I infectives) was an exact twin of Metabolism.limitation_pow and is removed.

/-- Two interventions: transmission kept is (1 − e₁)(1 − e₂) = 1 − (e₁ "or" e₂). -/
theorem interventions (e₁ e₂ : R) : (1 - e₁) * (1 - e₂) = 1 - circK (-1) e₁ e₂ := by
  unfold circK; grind

-- LIBRARY: seasons_growth (Growth factors multiply across seasons, so growth excesses r = R − 1 combine by ◯) was an exact twin of ChemQuantum.fragments_circ and is removed.

-- LIBRARY: hassell (Hassell's denominator: (1 + aN)^b = 1 + mpow(aN, b)) was an exact twin of MolecularBiology.mwc_R and is removed.

/-- Hassell with b = 1 is Beverton–Holt's 1 + aN. -/
theorem hassell_one (x : R) : 1 + mpow x 1 = 1 + x := by
  rw [mpow_eq]; grind

-- LIBRARY: survival (Survival through two age classes: deaths combine by the "or" law) was an exact twin of CellBiology.bliss_or and is removed.

/-- Three age classes, grouped either way. -/
theorem survival_assoc (q₁ q₂ q₃ : R) :
    circK (-1) (circK (-1) q₁ q₂) q₃ = circK (-1) q₁ (circK (-1) q₂ q₃) := by unfold circK; grind

-- LIBRARY: levins (Levins: colonisation from n occupied neighbours) was an exact twin of Metabolism.limitation_pow and is removed.

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

-- LIBRARY: herd_threshold (The herd-immunity threshold is the circle inverse (called mult-add inverse before 1.0) of the excess reproduction number) was an exact twin of Development.logistic_shortfall and is removed.

/-- Beverton–Holt map x ↦ R·x/(1 + x/M), written with reciprocal capacity m = 1/M. -/
def bh (Rg m x : K) : K := Rg * x / (1 + m * x)

/-- **Two seasons are one season.** BH(R₂, m₂) ∘ BH(R₁, m₁) = BH(R₁R₂, m₁ + R₁m₂). -/
theorem bh_compose (R₁ m₁ R₂ m₂ x : K) (h₁ : 1 + m₁ * x ≠ 0)
    (h : 1 + (m₁ + R₁ * m₂) * x ≠ 0) :
    bh R₂ m₂ (bh R₁ m₁ x) = bh (R₁ * R₂) (m₁ + R₁ * m₂) x := by
  unfold bh
  have i₁ := Lean.Grind.Field.inv_mul_cancel h₁
  have i := Lean.Grind.Field.inv_mul_cancel h
  have hd : 1 + m₂ * (R₁ * x / (1 + m₁ * x)) = (1 + (m₁ + R₁ * m₂) * x) / (1 + m₁ * x) := by
    simp only [Lean.Grind.Field.div_eq_mul_inv]; grind
  rw [hd]
  have hne : (1 + (m₁ + R₁ * m₂) * x) / (1 + m₁ * x) ≠ 0 := by
    intro z
    have : (1 + (m₁ + R₁ * m₂) * x) / (1 + m₁ * x) * (1 + m₁ * x) = 1 + (m₁ + R₁ * m₂) * x := by
      simp only [Lean.Grind.Field.div_eq_mul_inv]; grind
    rw [z] at this; apply h; grind
  have i₂ := Lean.Grind.Field.inv_mul_cancel hne
  simp only [Lean.Grind.Field.div_eq_mul_inv] at *
  grind

/-- **Room shrinks by 1/R.** With K = (R − 1)/m, the room K/N − 1 after one step is (K/N − 1)/R. -/
theorem bh_room (Rg m x : K) (hR : Rg ≠ 0) (hx : x ≠ 0) (hm : m ≠ 0) (h : 1 + m * x ≠ 0) :
    ((Rg - 1) / m) / bh Rg m x - 1 = (((Rg - 1) / m) / x - 1) / Rg := by
  unfold bh
  have iR := Lean.Grind.Field.inv_mul_cancel hR
  have ix := Lean.Grind.Field.inv_mul_cancel hx
  have im := Lean.Grind.Field.inv_mul_cancel hm
  have i1 := Lean.Grind.Field.inv_mul_cancel h
  simp only [Lean.Grind.Field.div_eq_mul_inv, Lean.Grind.Field.inv_mul, Lean.Grind.Field.inv_inv]
  grind

end Field

end Population
