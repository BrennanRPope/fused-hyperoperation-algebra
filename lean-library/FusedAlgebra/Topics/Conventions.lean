/-
  Conventions before the 1.0 freeze: the sign of k depends on what is composed
  ============================================================================

  Part of Fused Hyperoperation Algebra (1.0 checklist, item 4, names and
  notation). Companion: math-theory/names-review/README.md.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base only.

  The dial k is defined by a ✳ₖ b = a + b + k·a·b on the quantity being
  composed. Pages compose two kinds of quantity:
  * weights (fugacities, odds, dial values t), where one element's factor is
    1 + k·x, and
  * probabilities (occupations, eigenvalues μ), with μ = x/(1 + k·x).
  Changing from weights to probabilities flips the sign of k:
  1 − kμ = 1/(1 + kx) (`one_sub_prob`), and if the weights compose by ✳ₖ then
  the probabilities compose by ✳₋ₖ (`weight_prob_flip`). This is why the
  α-determinantal row (probabilities) has fermions at k = −1 while the counting
  row (weights) has "at most once" at k = +1: both are right, about different
  quantities. It is Tsallis' additive duality q ↔ 2 − q.
-/
import FusedAlgebra.Base

namespace Conventions

open MultAdd
open Lean.Grind

section Flip
variable {K : Type} [Lean.Grind.Field K]

/-- The probability of a weight x under dial k. -/
def prob (k x : K) : K := x / (1 + k * x)

/-- 1 − kμ = 1/(1 + kx). -/
theorem one_sub_prob (k x : K) (h : 1 + k * x ≠ 0) :
    1 - k * prob k x = (1 + k * x)⁻¹ := by
  unfold prob; grind

/-- **The sign flip.** If weights x, y compose to wAdd k x y, their probabilities
    compose by the opposite dial: prob of the composite = wAdd (−k) of the probs. -/
theorem weight_prob_flip (k x y : K) (hx : 1 + k * x ≠ 0) (hy : 1 + k * y ≠ 0) :
    prob k (wAdd k x y) = wAdd (-k) (prob k x) (prob k y) := by
  have hxy : 1 + k * wAdd k x y = (1 + k * x) * (1 + k * y) := wAdd_lift k x y
  have hne : 1 + k * wAdd k x y ≠ 0 := by
    rw [hxy]; intro h0
    rcases (Field.of_mul_eq_zero h0) with h | h <;> contradiction
  unfold prob wAdd at *
  grind

/-- Flipping twice gets back: prob (−k) (prob k x) ... the inverse map is
    x = μ/(1 − kμ). -/
theorem prob_inverse (k x : K) (h : 1 + k * x ≠ 0) :
    prob (-k) (prob k x) = x := by
  unfold prob; grind

end Flip

/-- Counting check: at k = 1 one weight x = 1 (a fair element) has probability 1/2,
    and two fair elements joined (weights 1 ✳ 1 = 3) have probability 3/4 =
    1/2 + 1/2 − 1/2·1/2, the "either" rule at k = −1. -/
theorem fair_example :
    prob (1 : Rat) 1 = 1/2 ∧ wAdd (1 : Rat) 1 1 = 3 ∧ prob (1 : Rat) 3 = 3/4 ∧
    wAdd (-1 : Rat) (1/2) (1/2) = 3/4 := by
  decide +kernel

end Conventions
