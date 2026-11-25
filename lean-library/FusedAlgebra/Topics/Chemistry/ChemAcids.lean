/-
  Chemistry 11, acids and bases
  =============================

  Companion page and checks: ../index.html, ../checks.py (14 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Independent acid groups combine at k = −1; the isoelectric point is
  where h² = K₁K₂; mixing buffers combines base/acid ratios by the mediant.
  The main test: the dial folded over three sites gives protonation
  coefficients 3, 3k, k², which pins one pKa gap at log 3, while a uniform
  pairwise interaction gives 3, 3k, k³, which makes the two gaps equal.
  Tabulated triprotic acids follow the pairwise pattern.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  LIBRARY NOTE: first compiled in the shared library, after ChemElectro. Two exact twins
  removed (indep_groups, chained_K); each is marked where it stood. No proof step or
  statement changed.
-/
import FusedAlgebra.Base

namespace ChemAcids
open CircleOp

open scoped CircleOp

-- LIBRARY: indep_groups (Independent groups: "at least one deprotonated" is failure-union, k = −1) was an exact twin of ChemStoichiometry.loss_two and is removed.

-- LIBRARY: chained_K (Chained dissociation constants multiply, so K − 1 combines by the circle operation (called mult-add before 1.0)) was an exact twin of InfoFused.campbell and is removed.

/-- Isoelectric point of a two-group amino acid: if h² = K₁K₂, the carboxylate fraction
    K₁/(K₁ + h) equals the ammonium fraction h/(K₂ + h). Cleared form. -/
theorem isoelectric (K1 K2 h : Rat) (hh : h * h = K1 * K2) : K1 * (K2 + h) = h * (K1 + h) := by
  grind

/-- Mixing two buffers: the mediant of b₁/a₁ ≤ b₂/a₂ is at least the smaller ratio (cleared, a's positive). -/
theorem mediant_above (a1 a2 b1 b2 : Rat) (h : b1 * a2 ≤ b2 * a1) : b1 * (a1 + a2) ≤ (b1 + b2) * a1 := by
  grind

/-- …and at most the larger one. -/
theorem mediant_below (a1 a2 b1 b2 : Rat) (h : b1 * a2 ≤ b2 * a1) : (b1 + b2) * a2 ≤ b2 * (a1 + a2) := by
  grind

/-- The dial folded over three equal sites: x ◯[k] x ◯[k] x = 3x + 3k x² + k² x³. -/
theorem fold_three (k x : Rat) : circK k (circK k x x) x = 3 * x + 3 * k * (x * x) + k * k * (x * x * x) := by
  unfold circK; grind

/-- Uniform pairwise interaction: β₁ = 3K, β₂ = 3kK², β₃ = k³K³ satisfy β₁³β₃ = β₂³,
    which says the two stepwise ratios K₁/K₂ and K₂/K₃ are equal (equal pKa gaps). -/
theorem pairwise_equal_gaps (k K : Rat) :
    (3 * K) * (3 * K) * (3 * K) * (k * k * k * (K * K * K)) = (3 * k * (K * K)) * (3 * k * (K * K)) * (3 * k * (K * K)) := by
  grind

/-- The fold's coefficients 3, 3k, k² break that equality by the factor k:
    β₁³β₃ − β₂³ = 27k²(1 − k), zero only at k = 0 or k = 1. -/
theorem fold_unequal_gaps (k : Rat) : 3 * 3 * 3 * (k * k) - (3 * k) * (3 * k) * (3 * k) = 27 * (k * k) * (1 - k) := by
  grind

/-- In the fold, the second and third stepwise constants are 3k/3 = k and k²/(3k) = k/3:
    their ratio is exactly 3 (a pKa gap of log 3) whatever k ≠ 0 is. -/
theorem fold_gap_log3 (k : Rat) (hk : k ≠ 0) : (3 * k / 3) / (k * k / (3 * k)) = 3 := by
  have h3 : (3 : Rat) ≠ 0 := by decide
  grind

end ChemAcids
