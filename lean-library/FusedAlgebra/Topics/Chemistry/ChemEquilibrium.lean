/-
  Chemistry 8, equilibrium and mass action
  ========================================

  Companion page and checks: ../index.html, ../checks.py (12 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Two binding sites with interaction factor ω: P − 1 = a ◯[ω] b.
  A symmetric diacid has ω = 4·10^(−ΔpKa).
  For two equal sites the binding polynomial 1 + 2x + ωx² has real zeros
  exactly when ω ≤ 1 (discriminant 4 − 4ω), so positive cooperativity and
  zeros off the real axis begin together at the circle-operation setting (called mult-add setting before 1.0) ω = 1.
  Folding the dial over n sites gives 1 + ((1 + kx)^n − 1)/k, the MWC model
  with binding only to the relaxed state and L = k − 1.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged unchanged at Topics/Chemistry/ChemEquilibrium; it now uses Base's ◯. No exact twins;
  the Rat helper sq_nonneg_rat stays local like the other topics' Rat helpers.
-/
import FusedAlgebra.Base

namespace ChemEquilibrium
open CircleOp

open scoped CircleOp

/-- Two-site binding polynomial minus one is the dial. -/
theorem binding_dial (a b ω : Rat) : (1 + a + b + ω * a * b) - 1 = circK ω a b := by
  rw [circK_def]; grind

/-- Symmetric diacid, association form: β₁ = 2K, β₂ = ωK², so β₂·4 = ω·β₁². -/
theorem diacid_gap (K ω : Rat) : (ω * K * K) * 4 = ω * ((2 * K) * (2 * K)) := by grind

/-- Two equal sites: discriminant of ωx² + 2x + 1 is 4(1 − ω). -/
theorem two_site_disc (ω : Rat) : 2 * 2 - 4 * ω * 1 = 4 * (1 - ω) := by grind

/-- If ω ≤ 1 the polynomial factors over the reals in principle; for ω > 1 there is no
    real (here: rational) zero, since ω(ωx² + 2x + 1) = (ωx + 1)² + (ω − 1). -/
theorem two_site_square (ω x : Rat) : ω * (ω * x * x + 2 * x + 1) = (ω * x + 1) * (ω * x + 1) + (ω - 1) := by
  grind

theorem sq_nonneg_rat (y : Rat) : 0 ≤ y * y := by
  rcases (Rat.le_total : (0 : Rat) ≤ y ∨ y ≤ 0) with h | h
  · exact Rat.mul_nonneg h h
  · have h2 : 0 ≤ -y := by grind
    have := Rat.mul_nonneg h2 h2
    grind

/-- Cooperative pairs (ω > 1) have no real binding zero. -/
theorem cooperative_no_real_zero (ω x : Rat) (h : 1 < ω) : ω * x * x + 2 * x + 1 ≠ 0 := by
  intro h0
  have e := two_site_square ω x
  have s := sq_nonneg_rat (ω * x + 1)
  rw [h0] at e
  grind

/-- Independent pair (ω = 1): the polynomial is (1 + x)², a double zero at x = −1. -/
theorem independent_double_zero (x : Rat) : 1 * x * x + 2 * x + 1 = (1 + x) * (1 + x) := by grind

/-- The dial folded over n copies of x: 1 + k·(fold) = (1 + kx)^n (library wpow_lift). -/
theorem fold_lift (k x : Rat) (n : Nat) : 1 + k * wpow k x n = (1 + k * x) ^ n := wpow_lift k x n

/-- MWC with exclusive binding: k·(1 + fold) = (k − 1) + (1 + kx)^n, i.e. L = k − 1. -/
theorem fold_is_mwc (k x : Rat) (n : Nat) : k * (1 + wpow k x n) = (k - 1) + (1 + k * x) ^ n := by
  rw [← wpow_lift]; grind

/-- Four sites: the fold's coefficients 4, 6k, 4k², k³ (j ligands carry k^(j−1)). -/
theorem four_site_coeffs (k x : Rat) :
    1 + k * (4 * x + 6 * k * x * x + 4 * k * k * x * x * x + k * k * k * x * x * x * x) = (1 + k * x) ^ 4 := by
  simp only [Rat.pow_succ, Rat.pow_zero]; grind

/-- Mass action: activities multiply into Q, so log Q adds; here Q for N₂ + 3H₂ ⇌ 2NH₃
    scales as p⁻² when every partial pressure is multiplied by p. -/
theorem haber_pressure (a b c p : Rat) (hp : p ≠ 0) (hb : b ≠ 0) (ha : a ≠ 0) :
    (p * c) * (p * c) / ((p * a) * ((p * b) * (p * b) * (p * b))) * (p * p) = c * c / (a * (b * b * b)) := by
  grind

end ChemEquilibrium
