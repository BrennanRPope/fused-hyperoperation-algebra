/-
  Affine Composition Theory · closed-leads-fresh-look (DRAFT, Nov 26 2026).
  First compiled in the shared library on Nov 26 2026, as FusedAlgebra.Topics.FreshLook.
  Core Lean 4 (no mathlib), same style as order-dial/lean/OrderDial.lean.
  Polynomial identities only, over Int:
  1. Tunnelling: the key step of the phase average, A² − B² = (1 + t₁ + t₂)², with
     A = 1 + t₁ ∘₂ t₂ and B² = 4 t₁ t₂ (1 + t₁)(1 + t₂).
  2. Tunnelling: 1 + 2(t₁ ∘₂ t₂) = (1 + 2t₁)(1 + 2t₂): the mean resistance law is ∘₂.
  3. Elliptic, Edwards parameter 0: the conic's norm a x² + y² is multiplicative under the law
     (x₁y₂ + y₁x₂, y₁y₂ − a x₁x₂): our k = 0 law is the norm-one group.
  4. Acids: equal off-diagonals give pair-product corner weights only at c = 0, 1 or −2.
  (The order gap a₁c₂ − a₂c₁ for reforms is already `OrderDial.order_gap`.)
-/
namespace FreshLook

/-- ∘ₖ on Int, written out so the file stands alone. -/
def circ (k a b : Int) : Int := a + b + k * a * b

theorem phase_average_key (t₁ t₂ : Int) :
    (1 + circ 2 t₁ t₂) * (1 + circ 2 t₁ t₂) - 4 * t₁ * t₂ * (1 + t₁) * (1 + t₂)
      = (1 + t₁ + t₂) * (1 + t₁ + t₂) := by
  simp only [circ]; grind

theorem mean_resistance_circ2 (t₁ t₂ : Int) :
    1 + 2 * circ 2 t₁ t₂ = (1 + 2 * t₁) * (1 + 2 * t₂) := by
  simp only [circ]; grind

theorem conic_norm_mul (a x₁ y₁ x₂ y₂ : Int) :
    a * (x₁ * y₂ + y₁ * x₂) * (x₁ * y₂ + y₁ * x₂) + (y₁ * y₂ - a * x₁ * x₂) * (y₁ * y₂ - a * x₁ * x₂)
      = (a * x₁ * x₁ + y₁ * y₁) * (a * x₂ * x₂ + y₂ * y₂) := by
  grind

theorem pair_product_gap (c : Int) :
    (1 - 3 * c * c + 2 * c * c * c) - (1 - c * c) * (1 - c * c) * (1 - c * c)
      = (c * c * c) * ((c - 1) * (c - 1)) * (c + 2) := by
  grind

end FreshLook
