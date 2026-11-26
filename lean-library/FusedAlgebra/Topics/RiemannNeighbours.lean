/-
  Fused Hyperoperation Algebra: Riemann's neighbours on the dial map (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib), no `sorry`. First compiled in the shared library on Nov 26,
  as FusedAlgebra.Topics.RiemannNeighbours (no proof changed).

  What is here: the finite, algebraic identities behind the neighbours write-up. Each one is
  the "positivity half" of a proved Riemann-type statement, said in our (k, d) plane. None of
  them touches the zeros of zeta.

  1. The norm form of the law.  4(m² + k·m·n + d·n²) = (2m + k·n)² + (4d − k²)·n².
     So the norm a² + k·ab + d·b² is positive definite exactly on the elliptic side k² < 4d.
     For a curve over F_p with dial point (k, d) = (−a_p, p) this norm is Weil's degree form
     deg(m − n·F), and its value at (1, 1) is the point count p + 1 − a_p (Hasse's proof).
  2. Selberg's factor.  s² − s + λ is the dial point (k, d) = (−1, λ), and
     4(s² − s + λ) = (2s − 1)² + (4λ − 1): roots on Re s = 1/2 exactly when λ ≥ 1/4.
  3. Ihara and Bartholdi.  One adjacency eigenvalue λ of a (q+1)-regular graph gives the factor
     1 − λ·u + q·u², the dial point (−λ, q).  Bartholdi's bump weight t only moves d:
     d_eff = (1 − t)((q + 1) − (1 − t)) = (1 − t)(q + t).
  4. Lee–Yang with two spins.  1 + a(z₁ + z₂) + z₁z₂ = 1 + a·(z₁ ◯[k] z₂) when a·k = 1:
     two coupled spins are our fused law with k = 1/a.  The uniform field z² + 2a·z + 1 is the
     dial point (2a, 1), whose discriminant 4a² − 4 is ≤ 0 exactly when |a| ≤ 1.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannNeighbours
open CircleOp

set_option linter.unusedSectionVars false

variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. The norm form: positive definite exactly on the elliptic side -/

/-- The norm of the law at dial point (k, d). -/
def normForm (k d m n : R) : R := m * m + k * m * n + d * n * n

/-- Completing the square: 4·N = (2m + kn)² + (4d − k²)·n². -/
theorem four_normForm (k d m n : R) :
    4 * normForm k d m n = (2 * m + k * n) * (2 * m + k * n) + (4 * d - k * k) * (n * n) := by
  unfold normForm; grind

/-- Weil's degree form for a curve over F_p: deg(1 − F) = p + 1 − a_p, the point count. -/
theorem normForm_frobenius_one (a p : R) : normForm (-a) p 1 1 = p + 1 - a := by
  unfold normForm; grind

/-- deg(m − n·F) with F² = a·F − p is the norm at (k, d) = (−a, p). -/
theorem normForm_frobenius (a p m n : R) :
    normForm (-a) p m n = m * m - a * m * n + p * n * n := by
  unfold normForm; grind

/-! ## 2. Selberg's factor s² − s + λ -/

/-- Selberg's local factor as a dial point: k = −1, d = λ. -/
theorem selberg_factor_dial (s l : R) : s * s - s + l = 1 * (s * s) + (-1) * s + l := by
  grind

/-- 4(s² − s + λ) = (2s − 1)² + (4λ − 1): the line Re s = 1/2 is the elliptic side λ ≥ 1/4. -/
theorem selberg_complete_square (s l : R) :
    4 * (s * s - s + l) = (2 * s - 1) * (2 * s - 1) + (4 * l - 1) := by
  grind

/-! ## 3. Ihara and Bartholdi -/

/-- The discriminant of 1 − λu + qu² is our Δ = k² − 4d at (k, d) = (−λ, q). -/
theorem ihara_discriminant (lam q : R) : (-lam) * (-lam) - 4 * q = lam * lam - 4 * q := by
  grind

/-- Bartholdi's bump weight t is the d dial in disguise: d_eff = (1 − t)(q + t). -/
theorem bartholdi_d_eff (q t : R) : (1 - t) * ((q + 1) - (1 - t)) = (1 - t) * (q + t) := by
  grind

/-! ## 4. Lee–Yang with two spins is the fused law -/

/-- Two coupled spins: 1 + a(z₁ + z₂) + z₁z₂ = 1 + a·(z₁ ◯[k] z₂) whenever a·k = 1. -/
theorem leeYang_two (a k z₁ z₂ : R) (h : a * k = 1) :
    1 + a * (z₁ + z₂) + z₁ * z₂ = 1 + a * circK k z₁ z₂ := by
  unfold circK
  have e : a * (z₁ + z₂ + k * z₁ * z₂) = a * z₁ + a * z₂ + (a * k) * z₁ * z₂ := by grind
  rw [e, h]; grind

/-- The uniform field z² + 2az + 1 is the dial point (2a, 1); its Δ is 4a² − 4. -/
theorem leeYang_uniform_discriminant (a : R) : (2 * a) * (2 * a) - 4 * 1 = 4 * (a * a - 1) := by
  grind

end CircleOp.RiemannNeighbours
