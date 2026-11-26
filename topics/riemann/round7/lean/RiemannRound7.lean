/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, seventh look: the quantum toolkit, and the hypothesis on solved shapes
  (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib). Imports FusedAlgebra.Base.
  Intended library path: FusedAlgebra/Topics/RiemannRound7.lean

  First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound7
  (namespace CircleOp.RiemannRound7, the library's convention; no proof changed).

  Every statement is a ring identity or a finite check. Nothing here touches the zeros of zeta.

  1. Prime gas: Fermi local factor; the order dial at q = −1 (Gaussian binomials at −1);
     the fugacity ladder's exponents (geometric sums).
  2. Trace formula: odds weights t^m give Witt multiplicities (cyclotomic identity at N = 6;
     values for t = 2, 3).
  3. Solved shapes: curvature factor completed to a square (Δ = k² − 4d decides the circle);
     Ihara factor from two roots; Lee–Yang for one ferromagnetic bond.
  4. Hilbert–Pólya: E = i(2s − 1) has imaginary part 2σ − 1.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound7

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. The prime gas -/

/-- Fermi local factor: 1 + u = (1 − u²)/(1 − u), so Z_F = ζ(β)/ζ(2β). Cleared: -/
theorem fermi_local (u : R) : (1 + u) * (1 - u) = 1 - u * u := by grind

/-- Gaussian binomial [2 choose 1]_q = 1 + q vanishes at q = −1: two distinct primes in a word
    cancel at the q = −1 end of the order dial. -/
theorem gauss21_at_minus_one (q : R) (h : q = -1) : 1 + q = 0 := by subst h; grind

/-- [4 choose 1]_q = 1 + q + q² + q³ vanishes at q = −1 (n even, k odd). -/
theorem gauss41_at_minus_one (q : R) (h : q = -1) : 1 + q + q * q + q * q * q = 0 := by subst h; grind

/-- [4 choose 2]_q = 1 + q + 2q² + q³ + q⁴ equals 2 = (2 choose 1) at q = −1: the rule
    [n; m]_{−1} = multinomial(⌊n/2⌋; ⌊m_i/2⌋) when at most one m_i is odd. -/
theorem gauss42_at_minus_one (q : R) (h : q = -1) :
    1 + q + 2 * (q * q) + q * q * q + q * q * q * q = 2 := by subst h; grind

/-- [3 choose 1]_q = 1 + q + q² equals 1 = (1 choose 0) at q = −1 (one odd part). -/
theorem gauss31_at_minus_one (q : R) (h : q = -1) : 1 + q + q * q = 1 := by subst h; grind

/-- Fugacity ladder: the exponent of log ζ in G_k is q^k + q^(k+1) + … = q^k/(1 − q).
    Finite form with p = q^k: -/
theorem ladder_geometric (p q : R) : (1 - q) * (p + p * q + p * q * q) = p - p * q * q * q := by grind

/-- Neighbouring rungs differ by exactly one power of ζ: e₀ − e₁ = (1 − q)/(1 − q). Cleared: -/
theorem rung_gap (q : R) : 1 - q = (1 - q) * 1 := by grind

/-! ## 2. Trace formula: orbit weights t^m -/

/-- Cyclotomic identity at N = 6: Σ_{d | 6} d·M(t, d) = t⁶, with M(t,1) = t, 2M(t,2) = t² − t,
    3M(t,3) = t³ − t, 6M(t,6) = t⁶ − t³ − t² + t. -/
theorem cyclotomic_six (t : R) :
    t + (t * t - t) + (t * t * t - t) + (t * t * t * t * t * t - t * t * t - t * t + t)
      = t * t * t * t * t * t := by grind

/-! ## 3. Solved shapes through the curvature dial -/

/-- 4d(1 + k u + d u²) = (2 d u + k)² + (4d − k²): with Δ = k² − 4d < 0 the factor is a sum of a
    square and a positive constant, and its two roots are conjugate with product 1/d, i.e. on |u| = d^(−1/2). -/
theorem curvature_square (k d u : R) :
    4 * d * (1 + k * u + d * (u * u)) = (2 * d * u + k) * (2 * d * u + k) + (4 * d - k * k) := by grind

/-- Ihara: (1 − a u)(1 − b u) = 1 − λ u + q u² with λ = a + b, q = a b. -/
theorem ihara_factor (a b u : R) : (1 - a * u) * (1 - b * u) = 1 - (a + b) * u + (a * b) * (u * u) := by grind

/-- The roots' spread is the discriminant: (a − b)² = λ² − 4q. Real spread (|λ| > 2√q) is exactly
    a non-Ramanujan eigenvalue. -/
theorem ihara_spread (a b : R) : (a - b) * (a - b) = (a + b) * (a + b) - 4 * (a * b) := by grind

/-- Lee–Yang for one ferromagnetic bond: Z(z) = 1 + 2 a z + z² = (z + a)² + (1 − a²); with 0 < a < 1
    the roots are conjugate with product 1, so they lie on |z| = 1. -/
theorem lee_yang_bond (a z : R) : 1 + 2 * a * z + z * z = (z + a) * (z + a) + (1 - a * a) := by grind

/-! ## 4. Hilbert–Pólya -/

/-- With s = σ + iγ, E = i(2s − 1) = −2γ + i(2σ − 1): multiplying by i = (0, 1) on pairs gives
    real part 0·(2σ − 1) − 1·(2γ) and imaginary part 0·(2γ) + 1·(2σ − 1). -/
theorem bbm_re (σ γ : R) : 0 * (2 * σ - 1) - 1 * (2 * γ) = -(2 * γ) := by grind

theorem bbm_im (σ γ : R) : 0 * (2 * γ) + 1 * (2 * σ - 1) = 2 * σ - 1 := by grind

/-- On the critical line (σ + σ = 1) the imaginary part vanishes: E is real. -/
theorem bbm_real_on_line (σ : R) (h : σ + σ = 1) : 2 * σ - 1 = 0 := by grind

/-- M(1/2, 2) = −1/8 is not whole: with 2t = 1, 8·M(t,2) = 4(t² − t) = −1. -/
theorem witt_half (t : R) (h : 2 * t = 1) : 4 * (t * t - t) = -1 := by grind

end Ring

/-! ## Witt numbers, finite checks on ℤ -/

/-- M(2, N) for N = 2, 3, 4, 6 is 1, 2, 3, 9 (Lyndon words in two letters). -/
theorem witt_two : (2 * 1 = (2:Int)^2 - 2) ∧ (3 * 2 = (2:Int)^3 - 2) ∧ (4 * 3 = (2:Int)^4 - 2^2)
    ∧ (6 * 9 = (2:Int)^6 - 2^3 - 2^2 + 2) := by decide

/-- M(3, N) for N = 2, 4 is 3, 18. -/
theorem witt_three : (2 * 3 = (3:Int)^2 - 3) ∧ (4 * 18 = (3:Int)^4 - 3^2) := by decide

end CircleOp.RiemannRound7
