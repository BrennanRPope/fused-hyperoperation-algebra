/-
DRAFT (Nov 26, 2026). Flux-tube dial: the algebra behind topics/force-reach/flux-tube-dial/README.md,
compiled in the shared library on Nov 26 as FusedAlgebra.Topics.FluxTubeDial.

Field energy density ∝ |D|^s / s in D space dimensions. A pair of charges at distance L has energy
∝ L^(D − (D−1)s); the flux's effective dimension is D_eff = 1 + (D−1)(s−1).

The session that wrote this stated it against Mathlib over ℝ (topics/force-reach/flux-tube-dial/lean/
FluxTubeDial.lean, nine statements). Here they are proved in core Lean over any commutative ring, any field
of characteristic 0, or (for the one inequality) any ordered field. `gilbert_s1` uses the real power 2^(s−1)
and needs Mathlib, so it is left out and checked by that folder's check_statements.py.
-/
import FusedAlgebra.Base

namespace CircleOp.FluxTubeDial

open Lean.Grind

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [CommRing R]

/-- Effective dimension of the flux: D_eff = 2 − (exponent of L) = 1 + (D−1)(s−1). -/
theorem deff_from_exponent (D s : R) :
    2 - (D - (D - 1) * s) = 1 + (D - 1) * (s - 1) := by grind

/-- At s = 2 (Maxwell) the flux spreads into all D directions. -/
theorem deff_maxwell (D : R) : 1 + (D - 1) * ((2 : R) - 1) = D := by grind

/-- At s = 1 the flux has one direction: a tube, in every dimension. -/
theorem deff_tube (D : R) : 1 + (D - 1) * ((1 : R) - 1) = 1 := by grind

end Ring

section Field
variable {K : Type} [Field K] [IsCharP K 0]

/-- The confinement threshold: the pair-energy exponent vanishes (log) exactly at s = D/(D−1). -/
theorem threshold_exponent_zero (D : K) (hD : D ≠ 1) :
    D - (D - 1) * (D / (D - 1)) = 0 := by
  have h : D - 1 ≠ 0 := by intro h0; apply hD; grind
  grind

/-- In three dimensions the threshold is s = 3/2. -/
theorem threshold_3d : (3 : K) / (3 - 1) = 3 / 2 := by grind

/-- Pagels–Tomboulis: field energy ∝ |E|^(2+4δ) gives |D|^s with s = (2+4δ)/(1+4δ); at δ = 1/4, s = 3/2. -/
theorem pt_quarter : (2 + 4 * (1 / 4 : K)) / (1 + 4 * (1 / 4)) = 3 / 2 := by grind

/-- The Pagels–Tomboulis map is s = 1 + 1/(1+4δ), so s → 1 as δ grows. -/
theorem pt_map_form (δ : K) (h : 1 + 4 * δ ≠ 0) :
    (2 + 4 * δ) / (1 + 4 * δ) = 1 + 1 / (1 + 4 * δ) := by grind

end Field

section Ordered
variable {K : Type} [LE K] [LT K] [Field K] [Std.IsLinearOrder K] [Std.LawfulOrderLT K] [OrderedRing K]

/-- Below the threshold the exponent is positive: energy grows with distance (confines). -/
theorem confines_below_threshold (D s : K) (hD : 1 < D) (hs : s < D / (D - 1)) :
    0 < D - (D - 1) * s := by
  have h : 0 < D - 1 := by grind
  have hne : D - 1 ≠ 0 := by intro h0; grind
  have e : (D - 1) * (D / (D - 1)) = D := by grind
  have h' : (D - 1) * s < (D - 1) * (D / (D - 1)) := OrderedRing.mul_lt_mul_of_pos_left hs h
  grind

end Ordered

end CircleOp.FluxTubeDial
