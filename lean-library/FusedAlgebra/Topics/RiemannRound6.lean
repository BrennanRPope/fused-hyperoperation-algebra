/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, sixth look: every pair of dials, with face power (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib). Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/RiemannRound6.lean

  First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound6
  (namespace CircleOp.RiemannRound6, the library's convention; no proof changed).

  Divisions are cleared and √3 is a symbol s with s·s = 3, so every statement is a ring identity
  (one is a finite check on ℕ). Nothing here touches the zeros of zeta.

  1. Odds or ratio with curvature: Δ scales by t², the point stays on its parabola.
  2. Height with curvature: the detector on a degree-2 local factor is (x² − w²)/(2x²).
  3. Face power on the hexagon: round at p = 2 and p = 4 (moments of the three directions), not at p = 6;
     p = 4 is the square of the hexagonal quadratic form.
  4. Light on face power: Wick rotation sends the hexagonal norm to 2a² + 2ac − c² (discriminant 12).
  5. Order with homothety: the leading power is a geometric sum (finite version); necklace step at N = 4.
  6. Roundness with homothety: no Euler product at p = 4, since 2 = 1⁴ + 1⁴ and 17 = 1⁴ + 2⁴
     but 34 is not a sum of two fourth powers.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound6

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. Odds (or ratio) with curvature -/

/-- The local factor 1 + k(tu) + d(tu)² has discriminant t²(k² − 4d): the sign of Δ is kept. -/
theorem disc_scaling (k d t : R) : (t * k) * (t * k) - 4 * (t * t * d) = t * t * (k * k - 4 * d) := by
  grind

/-- The scaled point (tk, t²d) stays on the parabola d = c k² through (k, d), cleared of division. -/
theorem parabola (k d t : R) : (t * t * d) * (k * k) = d * ((t * k) * (t * k)) := by grind

/-! ## 2. Height with curvature -/

/-- L = d(x² − w²), L′ = 2dx, L″ = 2d. Cleared: 2x² · L L″ = L′² · (x² − w²), i.e. h* = (x² − w²)/(2x²)
    = 1/2 − w²/(2x²): at most 1/2 when the roots are real (w² ≥ 0), above 1/2 when they are complex. -/
theorem hstar (d x w2 : R) :
    2 * (x * x) * ((d * (x * x - w2)) * (2 * d)) = ((2 * d * x) * (2 * d * x)) * (x * x - w2) := by grind

/-! ## 3. Face power on the hexagon -/

/-- The hexagon's three face projections of a lattice point a·b₁ + c·b₂ are a + c, a and c. Their fourth
    powers sum to twice the square of the hexagonal form: face power 4 is the round gauge squared. -/
theorem hex_p4 (a c : R) :
    (a + c) * (a + c) * (a + c) * (a + c) + a * a * a * a + c * c * c * c
      = 2 * ((a * a + a * c + c * c) * (a * a + a * c + c * c)) := by grind

/-- Second moments of the three directions (1, 0), (1/2, ±√3/2) are round, ×4 to clear halves. -/
theorem hex_design2 (s x y : R) (h : s * s = 3) :
    4 * (x * x) + (x + s * y) * (x + s * y) + (x - s * y) * (x - s * y) = 6 * (x * x + y * y) := by
  grind

/-- Fourth moments are round too, ×16. -/
theorem hex_design4 (s x y : R) (h : s * s = 3) :
    16 * (x * x * x * x) + (x + s * y) * (x + s * y) * (x + s * y) * (x + s * y)
      + (x - s * y) * (x - s * y) * (x - s * y) * (x - s * y)
      = 18 * ((x * x + y * y) * (x * x + y * y)) := by
  grind

/-- Sixth moments, ×64: along (1, 0) the sum is 66 ... -/
theorem hex_six_along_x (s : R) :
    64 * 1 + (1 + s * 0) * (1 + s * 0) * (1 + s * 0) * (1 + s * 0) * (1 + s * 0) * (1 + s * 0)
      + (1 - s * 0) * (1 - s * 0) * (1 - s * 0) * (1 - s * 0) * (1 - s * 0) * (1 - s * 0) = 66 := by
  grind

/-- ... and along (0, 1) it is 54, so face power 6 is not round. -/
theorem hex_six_along_y (s : R) (h : s * s = 3) :
    64 * 0 + (0 + s) * (0 + s) * (0 + s) * (0 + s) * (0 + s) * (0 + s)
      + (0 - s) * (0 - s) * (0 - s) * (0 - s) * (0 - s) * (0 - s) = 54 := by
  grind

/-! ## 4. Light on face power -/

/-- x = a + c/2, y = (√3/2)c; Wick rotation y ↦ iy turns x² + y² into x² − y². Cleared (×4):
    (2a + c)² − (s c)² = 2(2a² + 2ac − c²), a form of discriminant 2² + 4·2 = 12. -/
theorem wick_hex (s a c : R) (h : s * s = 3) :
    (2 * a + c) * (2 * a + c) - (s * c) * (s * c) = 2 * (2 * (a * a) + 2 * a * c - c * c) := by
  grind

theorem wick_disc (b a c : R) (hb : b = 2) (ha : a = 2) (hc : c = -1) : b * b - 4 * a * c = 12 := by
  subst hb; subst ha; subst hc; grind

/-! ## 5. Order with homothety -/

/-- The m = 0 term of the order dial is ∏_j Odds at k q^j, i.e. ζ to the power k(1 + q + q² + …); after
    the 1/k power the leading exponent is the geometric sum 1/(1 − q). Finite form: -/
theorem geometric4 (q : R) : (1 - q) * (1 + q + q * q + q * q * q) = 1 - q * q * q * q := by grind

/-- Necklace step at N = 4: c⁴ = M(c,1) + 2M(c,2) + 4M(c,4), with M(c,1) = c, 2M(c,2) = c² − c,
    4M(c,4) = c⁴ − c². -/
theorem necklace_four (c : R) : c + (c * c - c) + (c * c * c * c - c * c) = c * c * c * c := by grind

end Ring

/-! ## 6. Roundness with homothety: no Euler product at p = 4 -/

theorem two_is_sum : 1 ^ 4 + 1 ^ 4 = (2 : Nat) := by decide
theorem seventeen_is_sum : 1 ^ 4 + 2 ^ 4 = (17 : Nat) := by decide
/-- Fourth powers above 34 start at 3⁴ = 81, so a, b < 3 covers every case. -/
theorem thirty_four_not_sum : ∀ a < 3, ∀ b < 3, a ^ 4 + b ^ 4 ≠ (34 : Nat) := by decide

end CircleOp.RiemannRound6
