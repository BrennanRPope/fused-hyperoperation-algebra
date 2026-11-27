/-
  Fused Hyperoperation Algebra: the pendulum picture of the elliptic modulus (DRAFT, Nov 27 2026).
  Core Lean 4 only (no mathlib, no `sorry`). Polynomial identities behind pendulum_checks.py, over Int.
  Notation: κ = 1 − 2m the curvature reading, m = E/E_top the pendulum's energy as a fraction of the energy
  that just reaches the top. Square roots are named variables tied by their defining relation.
-/
namespace EllipticPendulum

/-- The curvature reading is the energy: with m = E/E_top, E_top · κ = E_top − 2E. -/
theorem kappa_is_energy (E Etop m : Int) (hm : m * Etop = E) : (1 - 2 * m) * Etop = Etop - 2 * E := by grind

/-- The halfway point. The reading t = sin θ/θ′ equals 1 where sin²θ = 2(cos θ − cos θ₀),
    that is, (with sin² + cos² = 1) where cos²θ + 2 cos θ − 1 − 2 cos θ₀ = 0. -/
theorem halfway_condition (s c c0 : Int) (h : s * s + c * c = 1) :
    (s * s - 2 * (c - c0) = 0) ↔ (c * c + 2 * c - 1 - 2 * c0 = 0) := by
  constructor <;> intro hh <;> grind

/-- ... and its root is cos θ* = 2 cos(θ₀/2) − 1 (write h = cos(θ₀/2), so cos θ₀ = 2h² − 1). -/
theorem halfway_root (h : Int) :
    (2 * h - 1) * (2 * h - 1) + 2 * (2 * h - 1) - 1 - 2 * (2 * h * h - 1) = 0 := by grind

/-- Over the top (m > 1), write r = √m, q = √(m − 1), so r² − q² = 1. The cap t = r − q and its partner
    r + q multiply to 1. -/
theorem cap_partner (r q : Int) (h : r * r - q * q = 1) : (r - q) * (r + q) = 1 := by grind

/-- The cap is a root of the law's quartic 1 + 2κt² + t⁴ with κ = 1 − 2m = 1 − 2r²: where G′ = 0. -/
theorem cap_is_root (r q : Int) (h : r * r - q * q = 1) :
    1 + 2 * (1 - 2 * (r * r)) * ((r - q) * (r - q)) + (r - q) * (r - q) * (r - q) * (r - q) = 0 := by grind

/-- The quartic is palindromic: t⁴ P(1/t) = P(t), written with cleared denominators as P(a, b) = P(b, a)
    for P(a, b) = b⁴ + 2κ a²b² + a⁴ (t = a/b). So walls and caps come in pairs t, 1/t. -/
theorem quartic_palindromic (κ a b : Int) :
    b*b*b*b + 2*κ*a*a*b*b + a*a*a*a = a*a*a*a + 2*κ*b*b*a*a + b*b*b*b := by grind

/-- Complementary time: with sn(K − u) = c/d, cn(K − u) = k′s/d, dn(K − u) = k′/d, the product
    G(u)·G(K − u) = (s d / c)·((c/d)(k′/d)/(k′ s/d)) = 1. Cleared: (s d)(c k′) = (c)(k′ s d). -/
theorem complementary_time (s c d k' : Int) : (s * d) * (c * k') = c * (k' * s * d) := by grind

end EllipticPendulum
