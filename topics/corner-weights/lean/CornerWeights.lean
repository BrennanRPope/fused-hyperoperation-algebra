/-
  CornerWeights: size, pair and Lucas weights from one matrix (Nov 25, 2026)
  =========================================================================
  The corners of the n-cube are the principal minors of one matrix M, and the corner sum
  is det(I + tM). Two independent choices give every weighting we have met:
    * the matrix: diagonal (our operation, no pairs), Cauchy 1/(pᵢ + qⱼ) (Hirota's pair
      weights for KdV, KP, sine-Gordon; Belinski–Zakharov's static solitons), or generic
      (three-body terms);
    * the dial you evaluate at: a number k (our size weight k^(m−1)) or a 2×2 matrix with
      trace k and determinant d (the Lucas weights U_m(k, d)).
  Companion: math-theory/corner-weights/README.md, checks.py (54 checks).
  Every statement here is also checked by identities_sympy_check.py (16 of 16 hold).
  STATUS: compiled Nov 25, 2026 with core Lean 4.37.0-rc1 (no mathlib, no `sorry`); all 12
  theorems check. In the shared library as FusedAlgebra.Topics.CornerWeights.
-/
import FusedAlgebra.Base

namespace CornerWeights

section Comm
variable {R : Type} [Lean.Grind.CommRing R]

/-- 3×3 determinant written out. -/
def det3 (a b c d e f g h i : R) : R :=
  a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

/-- The corner sum is a determinant: det(I + tM) = 1 + t·(rank-1 corners) + t²·(rank-2
    principal minors) + t³·det M. Each corner of the 3-cube is a principal minor. -/
theorem det3_I_add_tM (t a b c d e f g h i : R) :
    det3 (1 + t*a) (t*b) (t*c) (t*d) (1 + t*e) (t*f) (t*g) (t*h) (1 + t*i)
      = 1 + t * (a + e + i)
          + t * t * ((a*e - b*d) + (a*i - c*g) + (e*i - f*h))
          + t * t * t * det3 a b c d e f g h i := by
  unfold det3; grind

/-- 2×2 Cauchy minor, denominators cleared: the pair factor's numerator (pᵢ−pⱼ)(qᵢ−qⱼ). -/
theorem cauchy2_cleared (p1 p2 q1 q2 : R) :
    (p1 + q2) * (p2 + q1) - (p1 + q1) * (p2 + q2) = (p1 - p2) * (q1 - q2) := by
  grind

/-- 3×3 Cauchy determinant, denominators cleared (row i times ∏ₗ (xᵢ + yₗ)):
    the minor is a product over the three pairs, with no three-body factor. -/
theorem cauchy3_cleared (x1 x2 x3 y1 y2 y3 : R) :
    det3 ((x1 + y2) * (x1 + y3)) ((x1 + y1) * (x1 + y3)) ((x1 + y1) * (x1 + y2))
         ((x2 + y2) * (x2 + y3)) ((x2 + y1) * (x2 + y3)) ((x2 + y1) * (x2 + y2))
         ((x3 + y2) * (x3 + y3)) ((x3 + y1) * (x3 + y3)) ((x3 + y1) * (x3 + y2))
      = (x1 - x2) * (x1 - x3) * (x2 - x3) * ((y1 - y2) * (y1 - y3) * (y2 - y3)) := by
  unfold det3; grind

/-- KP: the pair factor's numerator is the gap between the two cross products, so
    A = 1 − (pᵢ+qᵢ)(pⱼ+qⱼ)/((pᵢ+qⱼ)(pⱼ+qᵢ)): one minus a cross-ratio. -/
theorem kp_pair_cross_ratio (pi pj qi qj : R) :
    (pi - pj) * (qi - qj) = (pi + qj) * (pj + qi) - (pi + qi) * (pj + qj) := by
  grind

/-- KdV: with u = (p − 1)/(p + 1) = tanh(θ/2), (pᵢ−pⱼ)/(pᵢ+pⱼ) is the relativistic
    velocity subtraction (uᵢ − uⱼ)/(1 − uᵢuⱼ). Numerator and denominator, cleared. -/
theorem kdv_pair_velocity (pi pj : R) :
    (pi - 1) * (pj + 1) - (pj - 1) * (pi + 1) = 2 * (pi - pj) ∧
    (pi + 1) * (pj + 1) - (pi - 1) * (pj - 1) = 2 * (pi + pj) := by
  constructor <;> grind

/-- Straight blend from a diagonal matrix (s = 0) toward unit-diagonal coupling (s = 1):
    the 3×3 corner is the product of the three pair factors plus a three-body remainder.
    Cauchy couplings are exactly the ones that make the remainder vanish at s = 1. -/
theorem blend_three_body (s c12 c13 c23 : R) :
    det3 1 (s*c12) (s*c13) (s*c12) 1 (s*c23) (s*c13) (s*c23) 1
      = (1 - s*s*c12*c12) * (1 - s*s*c13*c13) * (1 - s*s*c23*c23)
        + (2*s*s*s*c12*c13*c23
           - s*s*s*s*(c12*c12*c13*c13 + c12*c12*c23*c23 + c13*c13*c23*c23)
           + s*s*s*s*s*s*(c12*c12*c13*c13*c23*c23)) := by
  unfold det3; grind

/-- Lucas weights from two determinants (numerator): with k = α + β, d = αβ,
    det(I + αM) − det(I + βM) = (α − β)(e₁ + U₂ e₂ + U₃ e₃), U₂ = k, U₃ = k² − d. -/
theorem lucas_two_dets_num (α β a b c d e f g h i : R) :
    det3 (1 + α*a) (α*b) (α*c) (α*d) (1 + α*e) (α*f) (α*g) (α*h) (1 + α*i)
      - det3 (1 + β*a) (β*b) (β*c) (β*d) (1 + β*e) (β*f) (β*g) (β*h) (1 + β*i)
      = (α - β) * ((a + e + i)
          + (α + β) * ((a*e - b*d) + (a*i - c*g) + (e*i - f*h))
          + (α*α + α*β + β*β) * det3 a b c d e f g h i) := by
  unfold det3; grind

/-- Lucas weights from two determinants (denominator): α·det(I+βM) − β·det(I+αM)
    = (α − β)(1 − d·(U₁ e₂ + U₂ e₃)). -/
theorem lucas_two_dets_den (α β a b c d e f g h i : R) :
    α * det3 (1 + β*a) (β*b) (β*c) (β*d) (1 + β*e) (β*f) (β*g) (β*h) (1 + β*i)
      - β * det3 (1 + α*a) (α*b) (α*c) (α*d) (1 + α*e) (α*f) (α*g) (α*h) (1 + α*i)
      = (α - β) * (1 - α*β * (((a*e - b*d) + (a*i - c*g) + (e*i - f*h))
          + (α + β) * det3 a b c d e f g h i)) := by
  unfold det3; grind

/-- The plane law (x + y + kxy)/(1 − dxy) is our operation with dial α − β in the
    resolvent coordinate x/(1 + βx): denominator and numerator, cleared. -/
theorem plane_resolvent (α β x y : R) :
    (1 - α*β*x*y) + β * (x + y + (α + β)*x*y) = (1 + β*x) * (1 + β*y) ∧
    x + y + (α + β)*x*y = x*(1 + β*y) + y*(1 + β*x) + (α - β)*x*y := by
  constructor <;> grind

/-- Cayley–Hamilton for the companion matrix N = [[k, 1], [−d, 0]] (entries listed
    row by row): N³ = U₃·N − d·U₂·I with U₂ = k, U₃ = k² − d. -/
theorem companion_cube (k d : R) :
    let n11 := k; let n12 := (1 : R); let n21 := -d; let n22 := (0 : R)
    let s11 := n11*n11 + n12*n21; let s12 := n11*n12 + n12*n22
    let s21 := n21*n11 + n22*n21; let s22 := n21*n12 + n22*n22
    (s11*n11 + s12*n21 = (k*k - d)*n11 - d*k) ∧ (s11*n12 + s12*n22 = (k*k - d)*n12) ∧
    (s21*n11 + s22*n21 = (k*k - d)*n21) ∧ (s21*n12 + s22*n22 = (k*k - d)*n22 - d*k) := by
  intro n11 n12 n21 n22 s11 s12 s21 s22
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- The size × pair candidate with one shared pair factor A is not an iterated two-input
    law: three inputs differ from iterating a + b + kA·ab by k²A²(A − 1)·abc. -/
theorem size_pair_not_iterated (k A a b c : R) :
    (a + b + c + k*A*(a*b + a*c + b*c) + k*k*A*A*A*(a*b*c))
      - ((a + b + k*A*a*b) + c + k*A*(a + b + k*A*a*b)*c)
      = k*k*A*A*(A - 1)*(a*b*c) := by
  grind

/-- Hard core (A = 0) is plain addition. -/
theorem hard_core_is_addition (k a b c : R) :
    a + b + c + k*0*(a*b + a*c + b*c) + k*k*0*0*0*(a*b*c) = a + b + c := by
  grind

end Comm
end CornerWeights
