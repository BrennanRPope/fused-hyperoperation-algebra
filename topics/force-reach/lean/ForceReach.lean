/-
  Force reach (DRAFT, Nov 26 2026). Stated against Mathlib and not compiled (CI builds no Mathlib copy of it).
  Six of the eight statements are compiled in core Lean in lean-library/FusedAlgebra/Topics/ForceReach.lean;
  rung_mass (√, π) and yukawa_reduced (HasDerivAt) need Mathlib.
  Every statement is checked true numerically/symbolically by check_statements.py.
  The law a ∘ b = (a + b + k·a·b)/(1 − d·a·b) has generator v(x) = 1 + k x + d x², a Riccati flow.
  Its linear form u'' − k u' + d u = 0 has exponents with λ² − kλ + d = 0; Δ = k² − 4d.
-/
import Mathlib

namespace FusedAlgebra.Topics.ForceReach

/-- Δ is the squared spread of the two exponents of the law's linear form. -/
theorem disc_is_spread_sq (r₁ r₂ : ℝ) :
    (r₁ - r₂) ^ 2 = (r₁ + r₂) ^ 2 - 4 * (r₁ * r₂) := by ring

/-- Removing the drift: with λ = k/2 + μ, λ² − kλ + d = μ² − Δ/4, so the carrier mass is m² = Δ/4. -/
theorem drift_removal (k d μ : ℝ) :
    (k / 2 + μ) ^ 2 - k * (k / 2 + μ) + d = μ ^ 2 - (k ^ 2 - 4 * d) / 4 := by ring

/-- Massless exactly on the border: a double exponent (no length scale) happens iff Δ = 0. -/
theorem double_root_iff_border (k d : ℝ) :
    (∃ r : ℝ, r ^ 2 - k * r + d = 0 ∧ 2 * r = k) ↔ k ^ 2 - 4 * d = 0 := by
  constructor
  · rintro ⟨r, h1, h2⟩
    subst h2
    linear_combination (-4 : ℝ) * h1
  · intro h
    exact ⟨k / 2, by linear_combination (-1 / 4 : ℝ) * h, by ring⟩

/-- Raychaudhuri as the law: if x' = 1 + K x², then θ = −K x has θ' = −K − θ². -/
theorem raychaudhuri_is_law (K x : ℝ) :
    -K * (1 + K * x ^ 2) = -K - (-K * x) ^ 2 := by ring

/-- Curvature reading of the law with k = 0, d = K is K itself. -/
theorem curvature_reading_k0 (K : ℝ) : -((0 : ℝ) ^ 2 - 4 * K) / 4 = K := by ring

/-- Yukawa's effective dimension: 1 − r·(d ln F/dr) = 3 + m²r²/(1 + m r),
    for F = e^{−mr}(1 + m r)/r². -/
theorem yukawa_deff (m r : ℝ) (hr : 0 < r) (hm : 0 ≤ m) :
    1 - r * (-m - 2 / r + m / (1 + m * r)) = 3 + m ^ 2 * r ^ 2 / (1 + m * r) := by
  have h1 : (1 + m * r) ≠ 0 := by positivity
  have h2 : r ≠ 0 := hr.ne'
  field_simp
  ring

/-- Rung masses of a closed direction: circumference 2π/√(−Δ) gives spacing √(−Δ). -/
theorem rung_mass (Δ : ℝ) (h : Δ < 0) :
    2 * Real.pi / (2 * Real.pi / Real.sqrt (-Δ)) = Real.sqrt (-Δ) := by
  have hs : Real.sqrt (-Δ) ≠ 0 := (Real.sqrt_pos.mpr (by linarith)).ne'
  have hp : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- The reduced Yukawa profile u = e^{−m r}: its derivative −m e^{−mr} has derivative m² e^{−mr}. -/
theorem yukawa_reduced (m r : ℝ) :
    HasDerivAt (fun s => -m * Real.exp (-m * s)) (m ^ 2 * Real.exp (-m * r)) r := by
  have h : HasDerivAt (fun s => Real.exp (-m * s)) (Real.exp (-m * r) * (-m)) r := by
    have := ((hasDerivAt_id r).const_mul (-m)).exp
    simpa using this
  have h2 := h.const_mul (-m)
  convert h2 using 1
  ring

end FusedAlgebra.Topics.ForceReach
