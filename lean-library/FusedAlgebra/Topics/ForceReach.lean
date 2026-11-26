/-
DRAFT (Nov 26, 2026). Force reach: the algebra behind topics/force-reach/README.md, compiled in the
shared library on Nov 26 as FusedAlgebra.Topics.ForceReach.

The law a ∘ b = (a + b + k·a·b)/(1 − d·a·b) has generator v(x) = 1 + k x + d x², a Riccati flow.
Its linear form u'' − k u' + d u = 0 has exponents with λ² − kλ + d = 0; Δ = k² − 4d.

The session that wrote this stated it against Mathlib over ℝ (topics/force-reach/lean/ForceReach.lean,
eight statements). Six are pure algebra and are proved here in core Lean, over any field of
characteristic 0 (or any commutative ring where no division is needed). The other two need Mathlib:
`rung_mass` uses √ and π, and `yukawa_reduced` uses derivatives (HasDerivAt). They are left out of
the library and checked by topics/force-reach/lean/check_statements.py instead.
-/
import FusedAlgebra.Base

namespace CircleOp.ForceReach

open Lean.Grind

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [CommRing R]

/-- Δ is the squared spread of the two exponents of the law's linear form. -/
theorem disc_is_spread_sq (r₁ r₂ : R) :
    (r₁ - r₂) ^ 2 = (r₁ + r₂) ^ 2 - 4 * (r₁ * r₂) := by grind

/-- Raychaudhuri as the law: if x' = 1 + K x², then θ = −K x has θ' = −K − θ². -/
theorem raychaudhuri_is_law (K x : R) :
    -K * (1 + K * x ^ 2) = -K - (-K * x) ^ 2 := by grind

end Ring

section Field
variable {K : Type} [Field K] [IsCharP K 0]

/-- Removing the drift: with λ = k/2 + μ, λ² − kλ + d = μ² − Δ/4, so the carrier mass is m² = Δ/4. -/
theorem drift_removal (k d μ : K) :
    (k / 2 + μ) ^ 2 - k * (k / 2 + μ) + d = μ ^ 2 - (k ^ 2 - 4 * d) / 4 := by grind

/-- Massless exactly on the border: a double exponent (no length scale) happens iff Δ = 0. -/
theorem double_root_iff_border (k d : K) :
    (∃ r : K, r ^ 2 - k * r + d = 0 ∧ 2 * r = k) ↔ k ^ 2 - 4 * d = 0 := by
  constructor
  · rintro ⟨r, h1, h2⟩
    subst h2
    grind
  · intro h
    exact ⟨k / 2, by grind, by grind⟩

/-- Curvature reading of the law with k = 0, d = K is K itself. -/
theorem curvature_reading_k0 (c : K) : -((0 : K) ^ 2 - 4 * c) / 4 = c := by grind

/-- Yukawa's effective dimension: 1 − r·(d ln F/dr) = 3 + m²r²/(1 + m r),
    for F = e^{−mr}(1 + m r)/r² (the logarithmic derivative written out). -/
theorem yukawa_deff (m r : K) (hr : r ≠ 0) (h1 : 1 + m * r ≠ 0) :
    1 - r * (-m - 2 / r + m / (1 + m * r)) = 3 + m ^ 2 * r ^ 2 / (1 + m * r) := by
  grind

end Field

end CircleOp.ForceReach
