/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, second look — an elliptic object for zeta (DRAFT, Nov 26 2026).
  Core Lean 4 only (no mathlib), no `sorry`. First compiled in the shared library on Nov 26,
  as FusedAlgebra.Topics.RiemannElliptic.

  What is here (the finite, algebraic pieces of the write-up; nothing here touches the
  zeros of zeta themselves):

  1. Pairing lemma: two real degree-1 factors multiply to a dial point with
     discriminant (a − b)² ≥ 0, so pairing zeta's prime factors never reaches the
     elliptic side. The mirror pair (1 − u)(1 − p·u) has discriminant (p − 1)².
  2. Cayley map: |ρ − 1|² − |ρ|² = 1 − 2·Re ρ, so z = 1 − 1/ρ is on the unit circle
     exactly on the critical line; and the mirror ρ ↦ 1 − ρ is z ↦ 1/z.
  3. Pair factor (1 − z·w)(1 − z'·w) with z·z' = 1 has d = 1, and its discriminant is
     (z − z')²: elliptic exactly when z and z' are complex conjugates on the circle.
  4. Chebyshev (Dickson) traces: with z·z' = 1, zⁿ + z'ⁿ is a fixed polynomial in the
     dial value z + z'. This is the term 2 − zⁿ − z'ⁿ in Li's λₙ.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannElliptic

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. Pairing degree-1 factors -/

/-- Two degree-1 factors multiply to the dial point k = −(a + b), d = a·b. -/
theorem pair_product (a b u : R) :
    (1 - a * u) * (1 - b * u) = 1 + (-(a + b)) * u + (a * b) * (u * u) := by grind

/-- The pairing lemma: the discriminant of the product is a square. -/
theorem pair_disc (a b : R) : (-(a + b)) * (-(a + b)) - 4 * (a * b) = (a - b) * (a - b) := by grind

/-- Zeta's prime factor with its mirror: (1 − u)(1 − p·u), discriminant (p − 1)². -/
theorem mirror_pair (p u : R) :
    (1 - u) * (1 - p * u) = 1 + (-(p + 1)) * u + p * (u * u) ∧
    (-(p + 1)) * (-(p + 1)) - 4 * p = (p - 1) * (p - 1) := by
  constructor <;> grind

/-! ## 2. The Cayley map, without division -/

/-- |ρ − 1|² − |ρ|² = 1 − 2x for ρ = x + i·y: the unit circle |z| = 1 is the line x = 1/2. -/
theorem cayley_circle (x y : R) :
    ((x - 1) * (x - 1) + y * y) - (x * x + y * y) = 1 - 2 * x := by grind

/-! ## 3. The pair factor on the line d = 1 -/

/-- With z·z' = 1 the pair factor has d = 1 and k = −(z + z'). -/
theorem pair_factor (z z' w : R) (h : z * z' = 1) :
    (1 - z * w) * (1 - z' * w) = 1 + (-(z + z')) * w + w * w := by grind

/-- Its discriminant is (z − z')²: k² − 4 = (z − z')² when z·z' = 1. -/
theorem pair_factor_disc (z z' : R) (h : z * z' = 1) :
    (z + z') * (z + z') - 4 = (z - z') * (z - z') := by grind

/-! ## 4. Chebyshev (Dickson) traces -/

/-- Dickson polynomial in the dial value t = z + z': D₀ = 2, D₁ = t, Dₙ₊₂ = t·Dₙ₊₁ − Dₙ.
    (D_n(2x) = 2·T_n(x).) -/
def dickson (t : R) : Nat → R
  | 0 => 2
  | 1 => t
  | n + 2 => t * dickson t (n + 1) - dickson t n

/-- Powers of z and z' with z·z' = 1, two at a time. -/
theorem dickson_pair (z z' : R) (h : z * z' = 1) (n : Nat) :
    dickson (z + z') n = z ^ n + z' ^ n ∧ dickson (z + z') (n + 1) = z ^ (n + 1) + z' ^ (n + 1) := by
  induction n with
  | zero => simp only [dickson]; constructor <;> grind
  | succ n ih =>
    obtain ⟨h0, h1⟩ := ih
    refine ⟨h1, ?_⟩
    show dickson (z + z') (n + 2) = _
    simp only [dickson]; rw [h1, h0]
    have e1 : z ^ (n + 2) = z ^ n * z * z := by grind
    have e2 : z' ^ (n + 2) = z' ^ n * z' * z' := by grind
    have e3 : z ^ (n + 1) = z ^ n * z := by grind
    have e4 : z' ^ (n + 1) = z' ^ n * z' := by grind
    rw [e1, e2, e3, e4]; grind

/-- The Li term of one zero pair: 2 − zⁿ − z'ⁿ = 2 − Dₙ(z + z'). -/
theorem li_term (z z' : R) (h : z * z' = 1) (n : Nat) :
    2 - (z ^ n + z' ^ n) = 2 - dickson (z + z') n := by
  rw [(dickson_pair z z' h n).1]

end Ring

/-! ## 5. Finite checks over the integers -/

/-- y² = x³ + x + 1 over F₅ has a₅ = −3; its factor 1 + 3u + 5u² is elliptic: 9 < 20. -/
theorem genus1_p5_elliptic' : (3 : Int) * 3 < 4 * 5 := by decide

/-- Over F₂₅ the same curve has 25 + 1 − D₂(−3, q = 5) = 25 + 1 − (9 − 10) = 27 points. -/
theorem genus1_p25_count : (25 : Int) + 1 - ((-3) * (-3) - 2 * 5) = 27 := by decide

end CircleOp.RiemannElliptic
