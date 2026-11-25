/-
  Chemistry 7, chemical thermodynamics
  ====================================

  Companion page and checks: ../index.html, ../checks.py (11 checks).
  Core Lean 4.37.0-rc1, no mathlib, no `sorry`.

  Free energies add when reactions add, so equilibrium constants multiply:
  the K side is the lift, and K − 1 combines by the circle operation (called mult-add before 1.0).

  Two sites with occupation weights x, y and an interaction Boltzmann factor
  ω have Z = 1 + x + y + ωxy, so Z − 1 is the circle operation x ◯[ω] y.
  ω = 0 (hard exclusion) is plain addition, ω = 1 (no interaction) is
  the circle operation, and the occupancy covariance is xy(ω − 1)/Z²: the dial's
  distance from 1 is the correlation.

  LIBRARY NOTE: the local infix for circ was replaced by Base's scoped ◯ (`open scoped CircleOp`) in the 1.0 rename.
  Merged unchanged at Topics/Chemistry/ChemGibbs; it now uses Base's ◯. No exact twins.
-/
import FusedAlgebra.Base

namespace ChemGibbs
open CircleOp

open scoped CircleOp

/-- Equilibrium constants multiply; their excesses combine by the circle operation. -/
theorem K_excess (K1 K2 : Rat) : (K1 - 1) ◯ (K2 - 1) = K1 * K2 - 1 := by
  rw [circ_def]; grind

/-- Three coupled reactions: the excesses fold associatively. -/
theorem K_three (K1 K2 K3 : Rat) : (K1 - 1) ◯ (K2 - 1) ◯ (K3 - 1) = K1 * K2 * K3 - 1 := by
  rw [K_excess, K_excess]

/-- Two-site partition function with interaction factor ω. -/
def Z2 (x y ω : Rat) : Rat := 1 + x + y + ω * x * y

/-- Z − 1 is the circle operation with weight ω. -/
theorem Z2_dial (x y ω : Rat) : Z2 x y ω - 1 = circK ω x y := by
  unfold Z2; rw [circK_def]; grind

/-- Hard exclusion (ω = 0): occupations simply add. -/
theorem exclusion_adds (x y : Rat) : Z2 x y 0 = 1 + (x + y) := by
  unfold Z2; grind

/-- No interaction (ω = 1): the sites are independent and Z factors. -/
theorem independent_factors (x y : Rat) : Z2 x y 1 = (1 + x) * (1 + y) := by
  unfold Z2; grind

/-- Independent sites: Z − 1 is the circle operation. -/
theorem independent_circ (x y : Rat) : Z2 x y 1 - 1 = x ◯ y := by
  rw [Z2_dial, circK_one]

/-- Covariance identity, cleared of Z²: Z·(ωxy) − (x + ωxy)(y + ωxy) = xy(ω − 1). -/
theorem covariance (x y ω : Rat) :
    Z2 x y ω * (ω * x * y) - (x + ω * x * y) * (y + ω * x * y) = x * y * (ω - 1) := by
  unfold Z2; grind

/-- With x, y > 0, the covariance is positive exactly when ω > 1. -/
theorem covariance_pos (x y ω : Rat) (hx : 0 < x) (hy : 0 < y) (h : 1 < ω) :
    0 < x * y * (ω - 1) := by
  have hxy : 0 < x * y := Rat.mul_pos hx hy
  have hw : 0 < ω - 1 := by grind
  exact Rat.mul_pos hxy hw

/-- Haber: ΔS° = 2·192.77 − 191.61 − 3·130.68 = −198.11 J/K (as hundredths). -/
theorem haber_entropy : (2 * 19277 - 19161 - 3 * 13068 : Int) = -19811 := by decide

/-- Haber: ΔH° = 2 × (−45.90) = −91.80 kJ. -/
theorem haber_enthalpy : (2 * (-4590) : Int) = -9180 := by decide

/-- Free energies at a temperature add when reactions add (ΔG = ΔH − TΔS is linear). -/
theorem gibbs_adds (H1 S1 H2 S2 T : Rat) :
    (H1 + H2) - T * (S1 + S2) = (H1 - T * S1) + (H2 - T * S2) := by grind

end ChemGibbs
