/-
  Chemistry 13, biochemistry, with photosynthesis and nitrogenase
  ================================================================

  Companion page and checks: ../index.html, ../checks.py (13 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Two inhibitors (Yonetani–Theorell): the four-state binding polynomial is
  1 + a + b + k·ab with k = 1/α, so v₀/vᵢ − 1 = a ⊞ₖ b. k = 0 is mutual
  exclusion (Loewe additivity), k = 1 is independence (Bliss: remaining
  activities multiply). Steps in series add their times. Proofreading
  multiplies discrimination odds. Nitrogenase and Calvin-cycle bookkeeping.

  LIBRARY NOTE: `local infixl " ✳ "` for multAdd; delete it if Base declares ✳.
  LIBRARY NOTE: first compiled in the shared library, after ChemOrganic. Three exact twins
  removed (exclusive_loewe, bliss_fraction, proofreading); each is marked where it stood.
  No proof step or statement changed.
-/
import FusedAlgebra.Base

namespace ChemBiochem
open MultAdd

local infixl:65 " ✳ " => multAdd

/-- Yonetani–Theorell: with k = 1/α the binding polynomial minus one is a ⊞ₖ b. -/
theorem yonetani_theorell (k a b : Rat) : (1 + a + b + k * a * b) - 1 = wAdd k a b := by
  unfold wAdd; grind

-- LIBRARY: exclusive_loewe (Mutually exclusive inhibitors (k = 0): plain addition, Loewe's straight isobole) was an exact twin of InfoFused.shannon_weight and is removed.

/-- Independent inhibitors (k = 1): mult-add, so 1 + (a ⊞₁ b) = (1 + a)(1 + b). -/
theorem independent_lift (a b : Rat) : 1 + wAdd 1 a b = (1 + a) * (1 + b) := by
  unfold wAdd; grind

/-- …which is Bliss independence for the remaining activity: 1/(1+a) · 1/(1+b) = 1/(1 + a ⊞₁ b). -/
theorem bliss (a b : Rat) (ha : 1 + a ≠ 0) (hb : 1 + b ≠ 0) :
    (1 / (1 + a)) * (1 / (1 + b)) = 1 / (1 + wAdd 1 a b) := by
  rw [independent_lift]
  have hab : (1 + a) * (1 + b) ≠ 0 := by
    intro h; rcases Rat.mul_eq_zero.mp h with h | h <;> contradiction
  grind

-- LIBRARY: bliss_fraction (Bliss in fractional-inhibition form: fractions combine at k = −1) was an exact twin of ChemStoichiometry.loss_two and is removed.

/-- The 50% isobole for interaction k: b(1 + k·a) = 1 − a solves a ⊞ₖ b = 1. -/
theorem isobole (k a b : Rat) (h : b * (1 + k * a) = 1 - a) : wAdd k a b = 1 := by
  unfold wAdd; grind

/-- Irreversible steps in series: times add, so with t = 1/k the turnover time of three steps is t₁ + t₂ + t₃
    and the overall rate is its reciprocal (cleared: kcat·(t₁ + t₂ + t₃) = 1). -/
theorem series_rate (t1 t2 t3 : Rat) (h : t1 + t2 + t3 ≠ 0) : (1 / (t1 + t2 + t3)) * (t1 + t2 + t3) = 1 := by
  grind

-- LIBRARY: proofreading (Proofreading: two discrimination stages multiply the odds, so odds − 1 mult-adds) was an exact twin of InfoFused.campbell and is removed.

/-- Nitrogenase: 6 of 8 electrons reach ammonia, 3/4. -/
theorem nitrogenase_electrons : 6 * 4 = 3 * (8 : Nat) ∧ 2 * 3 + 2 = (8 : Nat) := by decide

/-- 16 ATP per N₂ is 8 per NH₃. -/
theorem nitrogenase_atp : 16 = 2 * (8 : Nat) := by decide

/-- Photosynthesis at the 8-photon limit: 478 kJ stored from 1407 kJ of 680 nm photons lies between 33% and 35%. -/
theorem calvin_efficiency : 33 * 1407 < 100 * 478 ∧ 100 * 478 < 35 * (1407 : Nat) := by decide

end ChemBiochem
