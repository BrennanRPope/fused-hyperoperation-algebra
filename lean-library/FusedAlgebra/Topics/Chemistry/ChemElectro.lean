/-
  Chemistry 10, electrochemistry and Nernst, with photosynthesis
  ===============================================================

  Companion page and checks: ../index.html, ../checks.py (22 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Coupled steps: equilibrium constants (and forward/backward rate ratios)
  multiply, so ratio − 1 mult-adds. Potentials combine by Luther's
  electron-weighted mean, which is idempotent and so not a mult-add.
  Two redox sites: Z − 1 = x ⊞ₖ x, with complex zeros exactly when k > 1,
  and the midpoint of the Nernst curve at x = 1/√k. Photosynthesis: stage
  losses compose at k = −1, and the Kok clock's damping eigenvalue has
  modulus at most 1.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  LIBRARY NOTE: first compiled in the shared library, after ChemRates. Three exact twins
  removed (ratio_series, bv_one_sided, stage_losses); each is marked where it stood. No
  proof step or statement changed.
-/
import FusedAlgebra.Base

namespace ChemElectro
open CircleOp

open scoped CircleOp

-- LIBRARY: ratio_series (Coupled steps: ratios (K, or forward/backward rate ratios) multiply, so ratio − 1 mult-adds) was an exact twin of InfoFused.campbell and is removed.

-- LIBRARY: bv_one_sided (One-sided Butler–Volmer (α = 1): with x = e^(Fη/RT), adding overpotentials multiplies x, so the current i/i₀ = x − 1 mult-adds) was an exact twin of InfoFused.campbell and is removed.

/-- Luther's rule: the combined potential m = (n₁E₁ + n₂E₂)/(n₁ + n₂) balances the two couples. -/
theorem luther_balance (n1 n2 E1 E2 : Rat) (h : n1 + n2 ≠ 0) :
    n1 * ((n1 * E1 + n2 * E2) / (n1 + n2) - E1) + n2 * ((n1 * E1 + n2 * E2) / (n1 + n2) - E2) = 0 := by
  have hm : (n1 * E1 + n2 * E2) / (n1 + n2) * (n1 + n2) = n1 * E1 + n2 * E2 := Rat.div_mul_cancel h
  grind

/-- Luther's mean is idempotent: a couple combined with itself gives the same potential,
    so potentials are averaged, not mult-added. -/
theorem luther_idem (n1 n2 E : Rat) (h : n1 + n2 ≠ 0) : (n1 * E + n2 * E) / (n1 + n2) = E := by
  have : n1 * E + n2 * E = E * (n1 + n2) := by grind
  rw [this, Rat.mul_div_cancel h]

/-- Two one-electron sites with interaction factor k: Z − 1 = x ⊞ₖ x. -/
theorem two_site (k x : Rat) : (1 + 2 * x + k * x * x) - 1 = circK k x x := by
  unfold circK; grind

/-- Its discriminant is 4(1 − k): complex zeros exactly when k > 1 (cooperativity). -/
theorem two_site_disc (k : Rat) : 2 * 2 - 4 * k = 4 * (1 - k) := by grind

/-- Midpoint of the two-electron Nernst curve: with s² = k, at x = 1/s the mean number of
    electrons 2x(1 + kx)/(1 + 2x + kx²) equals 1. Cleared form. -/
theorem nernst_midpoint (s : Rat) (hs : s ≠ 0) :
    2 * (1 / s) + 2 * (s * s) * (1 / s) * (1 / s) = 1 + 2 * (1 / s) + (s * s) * (1 / s) * (1 / s) := by
  grind

-- LIBRARY: stage_losses (Photosynthesis (and any chain of stages): losses compose at k = −1) was an exact twin of ChemStoichiometry.loss_two and is removed.

/-- Three stages fold the same way, in either grouping. -/
theorem stage_losses3 (L1 L2 L3 : Rat) :
    circK (-1) (circK (-1) L1 L2) L3 = 1 - (1 - L1) * (1 - L2) * (1 - L3) := by
  unfold circK; grind

/-- Kok clock: one flash advances 0, 1 or 2 steps with probabilities α, 1 − α − β, β.
    The mean advance is 1 − α + β. -/
theorem kok_mean (a b : Rat) : 0 * a + 1 * (1 - a - b) + 2 * b = 1 - a + b := by grind

/-- The damping eigenvalue α − β + i(1 − α − β) has squared modulus at most 1 when α, β ≥ 0
    and α + β ≤ 1, so the four-flash oscillation never grows. -/
theorem kok_damped (a b : Rat) (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : 0 ≤ 1 - a - b) :
    (a - b) * (a - b) + (1 - a - b) * (1 - a - b) ≤ 1 := by
  have p1 : 0 ≤ a * b := Rat.mul_nonneg ha hb
  have p2 : 0 ≤ (a + b) * (1 - a - b) := Rat.mul_nonneg (by grind) hab
  have e : 1 - ((a - b) * (a - b) + (1 - a - b) * (1 - a - b)) = 4 * (a * b) + 2 * ((a + b) * (1 - a - b)) := by grind
  grind

/-- Z-scheme bookkeeping in meV: 4 electrons × 1135 meV stored from 4 × 1823 + 4 × 1771 meV of photons,
    between 30% and a third. -/
theorem zscheme_bounds : 3 * (4 * 1823 + 4 * 1771) < 10 * (4 * 1135) ∧ 3 * (4 * 1135) < 4 * 1823 + 4 * 1771 := by
  decide

end ChemElectro
