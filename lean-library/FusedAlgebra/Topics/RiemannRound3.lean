/-
  Affine Composition Theory (built on Fused Hyperoperation Algebra):
  the Riemann hypothesis, third look: the fused rewind, the Kerr turn and the barrier chain
  (DRAFT, Nov 26 2026). Core Lean 4 only (no mathlib). Imports FusedAlgebra.Base.
  First compiled in the shared library on Nov 26, as FusedAlgebra.Topics.RiemannRound3
  (namespace CircleOp.RiemannRound3, the library's convention; no proof changed).

  What is here (the finite, algebraic pieces of the write-up; nothing here touches the
  zeros of zeta themselves):

  1. Fused pair: two roots a, b are kept as one point (s, p) = (a + b, a·b); a conjugate
     pair x ± i·y has the real fused point (2x, x² + y²) with discriminant −4y².
  2. The heat flow on a pair: the monic quadratic z² − s·z + (p − 2t) solves
     ∂H/∂t = −∂²H/∂z² (its second derivative is 2), and its discriminant is
     (s² − 4p) + 8t: a straight line in t that crosses zero, with nothing singular.
  3. Ring averages of a quadratic: shifting by +v and −v and averaging adds v² to the
     constant; with v turned by i (v·v = −w·w) it subtracts w².
  4. Pointwise steps that fix 0 and are one-to-one there leave every zero where it is.
  5. One barrier with a random phase: |1 − r·e^{iφ}|² = 1 + r² − 2r·cos φ.
-/
import FusedAlgebra.Base

namespace CircleOp.RiemannRound3

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. The fused pair -/

/-- Vieta: a pair of roots is the fused point (s, p) = (a + b, a·b). -/
theorem fused_pair (a b z : R) : (z - a) * (z - b) = z * z - (a + b) * z + a * b := by grind

/-- The discriminant of the fused point is the squared gap. -/
theorem fused_disc (a b : R) : (a + b) * (a + b) - 4 * (a * b) = (a - b) * (a - b) := by grind

/-- A conjugate pair x ± i·y, multiplied out without i: the fused point (2x, x² + y²). -/
theorem conj_pair (x y z : R) : (z - x) * (z - x) + y * y = z * z - (2 * x) * z + (x * x + y * y) := by grind

/-- Its discriminant is −4y², negative off the line. -/
theorem conj_disc (x y : R) : (2 * x) * (2 * x) - 4 * (x * x + y * y) = -(4 * (y * y)) := by grind

/-! ## 2. The heat flow on one pair -/

/-- The flowed pair's discriminant is a straight line in t. -/
theorem heat_pair_disc (s p t : R) : s * s - 4 * (p - 2 * t) = (s * s - 4 * p) + 8 * t := by grind

/-- Starting from a real pair with gap g, the discriminant is g² + 8t. -/
theorem heat_pair_gap (a b t : R) :
    (a + b) * (a + b) - 4 * (a * b - 2 * t) = (a - b) * (a - b) + 8 * t := by grind

/-! ## 3. Ring averages of a quadratic -/

/-- Real shifts ±v: the average adds v² to the constant term (a rewind step). -/
theorem real_ring (s p z v : R) :
    ((z + v) * (z + v) - s * (z + v) + p) + ((z - v) * (z - v) - s * (z - v) + p)
      = 2 * ((z * z - s * z + p) + v * v) := by grind

/-- Shifts turned by i (v·v = −w·w): the average subtracts w² (a forward step). -/
theorem turned_ring (s p z v w : R) (h : v * v = -(w * w)) :
    ((z + v) * (z + v) - s * (z + v) + p) + ((z - v) * (z - v) - s * (z - v) + p)
      = 2 * ((z * z - s * z + p) - w * w) := by
  have e := real_ring s p z v
  rw [e, h]; grind

/-! ## 5. One barrier with a random phase -/

/-- |1 − r·e^{iφ}|² with c = cos φ, q = sin φ, c² + q² = 1. -/
theorem barrier_sq (r c q : R) (h : c * c + q * q = 1) :
    (1 - r * c) * (1 - r * c) + (r * q) * (r * q) = 1 + r * r - 2 * (r * c) := by
  have : (1 - r * c) * (1 - r * c) + (r * q) * (r * q) = 1 - 2 * (r * c) + r * r * (c * c + q * q) := by grind
  rw [this, h]; grind

end Ring

/-! ## 4. Pointwise steps fix zeros -/

/-- If g fixes 0 and only 0 goes to 0, then g ∘ f has exactly the zeros of f. -/
theorem pointwise_fixes_zeros {α β : Type} [Zero β] (g : β → β) (h0 : g 0 = 0)
    (h1 : ∀ y, g y = 0 → y = 0) (f : α → β) (x : α) : g (f x) = 0 ↔ f x = 0 := by
  constructor
  · exact h1 (f x)
  · intro h; rw [h]; exact h0

end CircleOp.RiemannRound3
