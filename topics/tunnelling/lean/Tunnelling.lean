/-
  Affine Composition Theory · tunnelling (DRAFT, Nov 26 2026).
  First compiled in the shared library on Nov 26 2026, as FusedAlgebra.Topics.Tunnelling (the three rewind
  identities added the same day).
  Core Lean 4 (no mathlib), same style as closed-leads-fresh-look/lean/FreshLook.lean. Polynomial identities over Int only; the averages themselves (Jensen, the Legendre product formula) are
  checked numerically in ../checks.py, not here.

  1. mean_law_curvature: with k = −2K, 1 + k(t₁ ∘ₖ t₂) = (1 + kt₁)(1 + kt₂). At K = −1 this is the barriers'
     mean resistance law ∘₂; at K = +1 the mixers' (and the piling-up lemma's) ∘₋₂.
  2. piling_up: independent flips p₁, p₂ give a flip with probability p₁(1 − p₂) + p₂(1 − p₁) = p₁ ∘₋₂ p₂.
  3. pair_identity: (1 + c₁c₂)² − (c₁² − 1)(c₂² − 1) = (c₁ + c₂)², the one step behind both pair laws at k = 0
     (barriers: 1/⟨T⟩ − 1 = t₁ + t₂; mixers: ⟨1/(1 − s)⟩ = 1/|1 − s₁ − s₂|).
  4. sphere_paths: (1 − s₁)(1 − s₂) − s₁s₂ = 1 − s₁ − s₂: the stay–stay and switch–switch paths tie on the wall
     s₁ + s₂ = 1.
  5. hyperbolic_direct_wins: (1 + t₁)(1 + t₂) − t₁t₂ = 1 + t₁ + t₂ > 0 for t ≥ 0, so the barriers' direct path always
     outweighs the double reflection and the typical law ∘₁ has no wall.
  6–7. continuum_mean_unique, continuum_log_unique: the coefficient conditions from the flow ∂ₜ t(1 + εt) ∂ₜ force
     k(k − 2ε) = 0 for an exact mean law and k(k − ε) = 0 for an exact log law.
  8. rewind_identity (added with the fused-ladder rewind): (a ∘ₖ b) − b = a(1 + kb), so the
     rewind (c − b)/(1 + kb) is exact wherever 1 + kb ≠ 0. For a barrier on rung n, 1 + kb = P_n(1 + 2t₁) ≥ 1.
  9–10. rung_two_barrier, rung_two_mixer: 2P₂(1 + 2t) = 2 + 12(t + t²) and 2P₂(1 − 2s) = 2 − 12(s − s²)
     (with 2P₂(x) = 3x² − 1), the second rung's k = ±6 in polynomial form.
-/
namespace Tunnelling

/-- ∘ₖ on Int, written out so the file stands alone. -/
def circ (k a b : Int) : Int := a + b + k * a * b

theorem mean_law_curvature (K t₁ t₂ : Int) :
    1 + (-2 * K) * circ (-2 * K) t₁ t₂ = (1 + (-2 * K) * t₁) * (1 + (-2 * K) * t₂) := by
  simp only [circ]; grind

theorem piling_up (p₁ p₂ : Int) :
    p₁ * (1 - p₂) + p₂ * (1 - p₁) = circ (-2) p₁ p₂ := by
  simp only [circ]; grind

theorem pair_identity (c₁ c₂ : Int) :
    (1 + c₁ * c₂) * (1 + c₁ * c₂) - (c₁ * c₁ - 1) * (c₂ * c₂ - 1) = (c₁ + c₂) * (c₁ + c₂) := by
  grind

theorem sphere_paths (s₁ s₂ : Int) :
    (1 - s₁) * (1 - s₂) - s₁ * s₂ = 1 - s₁ - s₂ := by
  grind

theorem hyperbolic_direct_wins (t₁ t₂ : Int) (h₁ : 0 ≤ t₁) (h₂ : 0 ≤ t₂) :
    0 < (1 + t₁) * (1 + t₂) - t₁ * t₂ := by
  grind

theorem continuum_mean_unique (k c e : Int) (h₀ : k = c) (h₁ : 2 * e * k = c * k) :
    k * (k - 2 * e) = 0 := by
  subst h₀; grind

theorem continuum_log_unique (k c e : Int) (h₀ : k = c) (h₁ : 2 * e * k = 2 * c * k) :
    k * (k - e) = 0 := by
  subst h₀; grind

theorem rewind_identity (k a b : Int) :
    (a + b + k * a * b) - b = a * (1 + k * b) := by
  grind

theorem rung_two_barrier (t : Int) :
    3 * (1 + 2 * t) * (1 + 2 * t) - 1 = 2 + 12 * (t + t * t) := by
  grind

theorem rung_two_mixer (s : Int) :
    3 * (1 - 2 * s) * (1 - 2 * s) - 1 = 2 - 12 * (s - s * s) := by
  grind

end Tunnelling
