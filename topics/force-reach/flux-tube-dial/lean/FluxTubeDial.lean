/-
  Flux-tube dial (DRAFT, Nov 26 2026). Stated against Mathlib and not compiled (CI builds no Mathlib copy of it).
  Eight of the nine statements are compiled in core Lean in lean-library/FusedAlgebra/Topics/FluxTubeDial.lean;
  gilbert_s1 (a real power) needs Mathlib.
  Every statement is checked true symbolically by check_statements.py.
  Field energy density ∝ |D|^s / s in D space dimensions. A pair of charges at distance L
  has energy ∝ L^(D − (D−1)s); the flux's effective dimension is D_eff = 1 + (D−1)(s−1).
  All statements are pure algebra over ℝ (ring / linarith / field_simp), so they should
  compile with Mathlib's basic tactics only.
-/
import Mathlib

namespace FusedAlgebra.Topics.FluxTubeDial

/-- Effective dimension of the flux: D_eff = 1 + (D−1)(s−1) is the pair-energy exponent's partner,
    D_eff = 2 − (exponent of L). -/
theorem deff_from_exponent (D s : ℝ) :
    2 - (D - (D - 1) * s) = 1 + (D - 1) * (s - 1) := by ring

/-- At s = 2 (Maxwell) the flux spreads into all D directions. -/
theorem deff_maxwell (D : ℝ) : 1 + (D - 1) * ((2 : ℝ) - 1) = D := by ring

/-- At s = 1 the flux has one direction: a tube, in every dimension. -/
theorem deff_tube (D : ℝ) : 1 + (D - 1) * ((1 : ℝ) - 1) = 1 := by ring

/-- The confinement threshold: the pair-energy exponent vanishes (log) exactly at s = D/(D−1). -/
theorem threshold_exponent_zero (D : ℝ) (hD : D ≠ 1) :
    D - (D - 1) * (D / (D - 1)) = 0 := by
  have h : D - 1 ≠ 0 := sub_ne_zero.mpr hD
  field_simp
  ring

/-- In three dimensions the threshold is s = 3/2. -/
theorem threshold_3d : (3 : ℝ) / (3 - 1) = 3 / 2 := by norm_num

/-- Below the threshold the exponent is positive: energy grows with distance (confines). -/
theorem confines_below_threshold (D s : ℝ) (hD : 1 < D) (hs : s < D / (D - 1)) :
    0 < D - (D - 1) * s := by
  have h : 0 < D - 1 := by linarith
  have := (lt_div_iff₀ h).mp hs
  linarith

/-- Pagels–Tomboulis: L ∝ F²(F²)^δ means field energy ∝ |E|^(2+4δ); in terms of D it is
    |D|^s with s = (2+4δ)/(1+4δ). At δ = 1/4 this is s = 3/2, the 3D threshold. -/
theorem pt_quarter : (2 + 4 * (1 / 4 : ℝ)) / (1 + 4 * (1 / 4)) = 3 / 2 := by norm_num

/-- The Pagels–Tomboulis map is s = 1 + 1/(1+4δ), so s → 1 as δ grows. -/
theorem pt_map_form (δ : ℝ) (h : 1 + 4 * δ ≠ 0) :
    (2 + 4 * δ) / (1 + 4 * δ) = 1 + 1 / (1 + 4 * δ) := by
  field_simp
  ring

/-- Gilbert junction for two equal unit branches merging into flux 2: the force balance
    2·cos(θ) = 2^s gives cos θ = 2^(s−1); at s = 1 the half angle is 0 (cos = 1). -/
theorem gilbert_s1 : (2 : ℝ) ^ ((1 : ℝ) - 1) = 1 := by norm_num

end FusedAlgebra.Topics.FluxTubeDial
