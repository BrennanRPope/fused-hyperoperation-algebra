/-
  Fused Hyperoperation Algebra: physiology (muscle, drug antagonists, organ extraction).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Physiology.lean

  1. Hill's force–velocity law (F + a)(v + b) = (F₀ + a)b is the level set X ✳ Y = X₀ in the
     reduced variables X = F/a, Y = v/b; its asymptotes are mult-add's absorbing −1.
  2. Two antagonists: the agonist dose ratio obeys DR₁₂ − 1 = (DR₁ − 1) ✳_ω (DR₂ − 1);
     same site ω = 0 (dose ratios add), independent sites ω = 1 (they multiply).
  3. Organ extraction with u = CLint/Q: parallel pathways add u (dial 0), organs in series
     combine u by plain ✳, and first-pass losses combine by the "or" law.
  LIBRARY NOTE: first compiled in the shared library. 5 exact twins removed (star_def,
  hill_asymptote, first_pass, series_organs, extraction); each is marked where it stood.
  No proof step or statement changed.
-/
import FusedAlgebra.Base

namespace Physiology
open MultAdd

set_option linter.unusedSectionVars false

/-- Brennan's notation for mult-add in this file: a ✳ b = a + b + ab. -/
local infixl:65 " ✳ " => multAdd

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

/-! ## 1. Hill's muscle law -/

/-- In reduced variables Hill's hyperbola (1 + X)(1 + Y) = 1 + X₀ is X ✳ Y = X₀. -/
theorem hill_level_set (X Y X₀ : R) (h : (1 + X) * (1 + Y) = 1 + X₀) : X ✳ Y = X₀ := by
  unfold multAdd; grind

-- LIBRARY: hill_asymptote (The asymptote X = −1 (force −a): mult-add's absorbing element) was an exact twin of MultAdd.neg_one_absorbs and is removed.

/-! ## 2. Antagonists -/

/-- The antagonist part of the receptor's partition function, with interaction ω. -/
def antag (ω b₁ b₂ : R) : R := 1 + b₁ + b₂ + ω * b₁ * b₂

/-- The dose ratio is that partition function; DR − 1 is the weighted mult-add. -/
theorem dose_ratio (ω b₁ b₂ : R) : antag ω b₁ b₂ - 1 = wAdd ω b₁ b₂ := by
  unfold antag wAdd; grind

/-- Same site: dose ratios add (Paton–Rang additivity). -/
theorem same_site (b₁ b₂ : R) :
    antag 0 b₁ b₂ - 1 = (antag 0 b₁ 0 - 1) + (antag 0 0 b₂ - 1) := by unfold antag; grind

/-- Independent sites: dose ratios multiply. -/
theorem separate_sites (b₁ b₂ : R) : antag 1 b₁ b₂ = antag 1 b₁ 0 * antag 1 0 b₂ := by
  unfold antag; grind

/-! ## 3. Organs -/

-- LIBRARY: first_pass (First-pass loss: escapes multiply, so losses combine by the "or" law) was an exact twin of CellBiology.bliss_or and is removed.

-- LIBRARY: series_organs (Organs in series: 1/F = (1 + u₁)(1 + u₂) = 1 + (u₁ ✳ u₂)) was an exact twin of ChemQuantum.fragments_multAdd and is removed.

/-- Splitting clearance across two organs in series beats pooling it in one:
    (1 + u₁)(1 + u₂) exceeds 1 + u₁ + u₂ by exactly u₁u₂, the corner cell. -/
theorem series_gain (u₁ u₂ : R) : (1 + u₁) * (1 + u₂) - (1 + (u₁ + u₂)) = u₁ * u₂ := by grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- **Velocity is a mult-add difference.** On Hill's curve Y = X₀ ✳ inv(X). -/
theorem hill_velocity (X Y X₀ : K) (hX : 1 + X ≠ 0) (h : X ✳ Y = X₀) :
    Y = X₀ ✳ multAddInv X := by
  unfold multAdd multAddInv at *
  have e : Y * (1 + X) = X₀ - X := by grind
  grind

-- LIBRARY: extraction (Well-stirred extraction is minus the inverse: E = u/(1 + u)) was an exact twin of Metabolism.free_fraction and is removed.

end Field

end Physiology
