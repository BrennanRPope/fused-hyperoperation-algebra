/-
  Fused Hyperoperation Algebra: genetics and heredity.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Genetics.lean

  1. Linkage without interference (Haldane): recombination fractions of adjacent intervals
     combine by the dial at k = −2, and the wall −1/k = 1/2 is "unlinked".
  2. Interference with coincidence c is the dial at k = −2c (bookkeeping identity).
  3. Lens laws (a + b + k·ab)/(1 − q·ab) keep r = 1/2 absorbing exactly on the line
     2k + q = −4; Haldane is (−2, 0) and Kosambi (0, −4). Kosambi is plain ✳ after w = 2R/(1−R).
  4. Gene–environment risk: excess risks combine on the dial; k = 0 additive, k = 1
     multiplicative, and the interaction contrast RERI is k·ERR₁·ERR₂.
  LIBRARY NOTE: first compiled in the shared library. 3 exact twins removed (star_def,
  risk_multiplicative, affected_children); each is marked where it stood. No proof step or
  statement changed.
-/
import FusedAlgebra.Base

namespace Genetics
open MultAdd

set_option linter.unusedSectionVars false

/-- Brennan's notation for mult-add in this file: a ✳ b = a + b + ab. -/
local infixl:65 " ✳ " => multAdd

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

/-! ## 1. Haldane -/

/-- No interference: r₁₂ = r₁ + r₂ − 2r₁r₂, the dial at k = −2. -/
theorem haldane (r₁ r₂ : R) : r₁ + r₂ - 2 * r₁ * r₂ = wAdd (-2) r₁ r₂ := by unfold wAdd; grind

/-- The lift of dial −2 is 1 − 2r, so 1 − 2r multiplies along the chromosome. -/
theorem haldane_lift (r₁ r₂ : R) : 1 - 2 * wAdd (-2) r₁ r₂ = (1 - 2 * r₁) * (1 - 2 * r₂) := by
  unfold wAdd; grind

/-- In the coordinate u = −2r Haldane is plain mult-add. -/
theorem haldane_star (r₁ r₂ : R) : -2 * wAdd (-2) r₁ r₂ = (-2 * r₁) ✳ (-2 * r₂) := by
  unfold wAdd multAdd; grind

/-- Three intervals: the lift is a product of three factors. -/
theorem haldane_three (r₁ r₂ r₃ : R) :
    1 - 2 * wAdd (-2) (wAdd (-2) r₁ r₂) r₃ = (1 - 2 * r₁) * (1 - 2 * r₂) * (1 - 2 * r₃) := by
  unfold wAdd; grind

/-! ## 2. Interference -/

/-- With double-recombinant frequency d₁₂ = c·r₁·r₂, the outer recombination fraction is the
    dial at k = −2c. -/
theorem interference (c r₁ r₂ : R) : r₁ + r₂ - 2 * (c * r₁ * r₂) = wAdd (-2 * c) r₁ r₂ := by
  unfold wAdd; grind

/-- Complete interference (c = 0) is plain addition of map lengths. -/
theorem complete_interference (r₁ r₂ : R) : wAdd (-2 * 0) r₁ r₂ = r₁ + r₂ := by unfold wAdd; grind

/-! ## 4. Gene–environment risk -/

/-- The interaction contrast (RERI) is the dial times the two excess risks. -/
theorem reri (k e₁ e₂ : R) : wAdd k e₁ e₂ - e₁ - e₂ = k * (e₁ * e₂) := by unfold wAdd; grind

-- LIBRARY: risk_multiplicative (k = 1: relative risks multiply) was an exact twin of LeeYang.edge_gas_free and is removed.

/-- k = 0: excess risks add. -/
theorem risk_additive (e₁ e₂ : R) : wAdd 0 e₁ e₂ = e₁ + e₂ := by unfold wAdd; grind

-- LIBRARY: affected_children (Children of two carriers: at least one affected among n is wpow(−1, 1/4, n), stated for any per-child risk p) was an exact twin of Metabolism.limitation_pow and is removed.

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- **The wall is "unlinked".** An interval at r = 1/2 makes the combined pair unlinked,
    whatever the other interval is. -/
theorem haldane_unlinked (r : K) (h2 : (2 : K) ≠ 0) : wAdd (-2) r (1 / 2) = 1 / 2 := by
  have hk : (-2 : K) ≠ 0 := by intro e; apply h2; grind
  have e : (1 : K) / 2 = -1 / -2 := by grind
  rw [e, wAdd_wall (-2) r hk]

/-- The lens law with dial k and denominator weight q. -/
def lens (k q a b : K) : K := (a + b + k * a * b) / (1 - q * a * b)

/-- **Unlinked absorbs exactly on the line 2k + q = −4.** -/
theorem lens_unlinked (k q a : K) (h2 : (2 : K) ≠ 0) (hline : 2 * k + q = -4)
    (hd : 1 - q * a * (1 / 2) ≠ 0) : lens k q a (1 / 2) = 1 / 2 := by
  unfold lens
  have hq : q = -4 - 2 * k := by grind
  subst hq
  grind

/-- Haldane is the lens law at (−2, 0). -/
theorem lens_haldane (a b : K) : lens (-2) 0 a b = wAdd (-2) a b := by
  unfold lens wAdd; grind

/-- Kosambi's addition rule is the lens law at (0, −4). -/
theorem lens_kosambi (a b : K) : lens 0 (-4) a b = (a + b) / (1 + 4 * a * b) := by
  unfold lens; grind

/-- **Kosambi is plain ✳ after a change of variable.** With R = 2r and
    R₁₂ = (R₁ + R₂)/(1 + R₁R₂), the quantity 1 + w = (1 + R)/(1 − R) multiplies. -/
theorem kosambi_lift (R₁ R₂ : K) (hd : 1 + R₁ * R₂ ≠ 0) :
    (1 + (R₁ + R₂) / (1 + R₁ * R₂)) * ((1 - R₁) * (1 - R₂)) =
      (1 - (R₁ + R₂) / (1 + R₁ * R₂)) * ((1 + R₁) * (1 + R₂)) := by
  have e1 : 1 + (R₁ + R₂) / (1 + R₁ * R₂) = (1 + R₁) * (1 + R₂) / (1 + R₁ * R₂) := by grind
  have e2 : 1 - (R₁ + R₂) / (1 + R₁ * R₂) = (1 - R₁) * (1 - R₂) / (1 + R₁ * R₂) := by grind
  rw [e1, e2]; grind

end Field

end Genetics
