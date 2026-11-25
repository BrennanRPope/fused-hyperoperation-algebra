/-
  The zoom dial: a size for each piece, with dimension as the output
  ==================================================================

  Part of Fused Hyperoperation Algebra (1.0 checklist, item 2). Companion files:
  math-theory/zoom-dial/ (README.md, checks.py, figures). Completes the
  fractal-dials page, which found the circle operation (called mult-add before 1.0) had every dial but zoom.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base only.

  The dial: a piece shrunk by zoom r carries the weight r^s.
  * Pieces side by side add their weights: x(s) = Σ rᵢˢ.
  * Zooming in again and again is a word of pieces; all words of depth n weigh
    x(s)ⁿ in total (`words_weigh`), so all depths together are the
    "any number of times" element 1/(1 − x), the k = −1 end of the repetition
    dial (`word_series`: (1 − x)·Σ_{j<n} xʲ = 1 − xⁿ).
  * Dimension is where that element hits its wall, x(d) = 1 (Moran 1946): at the
    wall every depth weighs exactly 1 (`wall_every_depth`) and the wordPartial sums
    grow without end (`wall_partial`). The five half boxes of the fractal page
    are certified at d = 1/2 (`half_box_walls`).
  * A zoom is a twist: scaling dial values by the natural weight w = r^d turns
    the circle operation into the weighted rule at k = w (`zoom_is_twist`), and two zooms
    in a row multiply (`zoom_compose`), so zoom 1/2 at dimension β is the size
    dial's doubling weight 2^(−β).
  Complex dimensions (the wobble) are the complex s with x(s) = 1; they are
  computed in checks.py, not proved here.
-/
import FusedAlgebra.Base

namespace ZoomDial

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

section Words
variable {R : Type} [Lean.Grind.CommRing R]

/-- All words of depth n over the pieces (each letter is one piece's weight). -/
def words (l : List R) : Nat → List (List R)
  | 0 => [[]]
  | n + 1 => l.flatMap (fun r => (words l n).map (fun w => r :: w))

theorem sumL_flatMap_cons (l : List R) (G : List (List R)) :
    sumL ((l.flatMap (fun r => G.map (fun w => r :: w))).map prodL)
      = sumL l * sumL (G.map prodL) := by
  induction l with
  | nil => simp only [List.flatMap_nil, List.map_nil, sumL]; grind
  | cons a l ih =>
    rw [List.flatMap_cons, List.map_append, sumL_append, sumL_map_cons, ih]
    simp only [sumL]; grind

/-- **Words weigh x(s)ⁿ.** The words of depth n, each weighing the product of
    its pieces, add up to (Σ pieces)ⁿ. -/
theorem words_weigh (l : List R) (n : Nat) :
    sumL ((words l n).map prodL) = sumL l ^ n := by
  induction n with
  | zero => simp only [words, List.map, sumL, prodL]; grind
  | succ n ih => simp only [words]; rw [sumL_flatMap_cons, ih, pow_add' (sumL l) n 1]; grind

/-- Partial sums of the word series: Σ_{j<n} xʲ. -/
def wordPartial (x : R) : Nat → R
  | 0 => 0
  | n + 1 => wordPartial x n + x ^ n

/-- **All depths together are the k = −1 element:** (1 − x)·Σ_{j<n} xʲ = 1 − xⁿ. -/
theorem word_series (x : R) (n : Nat) : (1 - x) * wordPartial x n = 1 - x ^ n := by
  induction n with
  | zero => simp only [wordPartial]; grind
  | succ n ih => simp only [wordPartial]; rw [pow_add' x n 1]; grind

/-- **At the wall every depth weighs 1:** if the pieces' weights add to 1, so do
    the words of every depth. -/
theorem wall_every_depth (l : List R) (h : sumL l = 1) (n : Nat) :
    sumL ((words l n).map prodL) = 1 := by
  rw [words_weigh, h, one_pow']

end Words

/-- At the wall the partial sums are n: the series never settles. -/
theorem wall_partial (n : Nat) : wordPartial (1 : Int) n = (n : Int) := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [wordPartial]; rw [ih, one_pow']; grind

/-! ## Dimension one half: the walls of the five half boxes -/

/-- The fractal page's sets A to E at s = 1/2: the natural weights √rᵢ add to 1,
    so x(1/2) = 1. (A: two at 1/4; B, C: three at 1/9; D: 9/100 and 49/100;
    E: four at 1/16.) -/
theorem half_box_walls :
    ((1/2 : Rat) * (1/2) = 1/4 ∧ sumL [(1/2 : Rat), 1/2] = 1) ∧
    ((1/3 : Rat) * (1/3) = 1/9 ∧ sumL [(1/3 : Rat), 1/3, 1/3] = 1) ∧
    ((3/10 : Rat) * (3/10) = 9/100 ∧ (7/10 : Rat) * (7/10) = 49/100 ∧ sumL [(3/10 : Rat), 7/10] = 1) ∧
    ((1/4 : Rat) * (1/4) = 1/16 ∧ sumL [(1/4 : Rat), 1/4, 1/4, 1/4] = 1) := by
  decide +kernel

/-- Uneven set D at depth 3: the eight words still weigh exactly 1. -/
theorem uneven_depth_three : sumL ((words [(3/10 : Rat), 7/10] 3).map prodL) = 1 := by
  decide +kernel

/-! ## A zoom is a twist -/

section Twist
variable {R : Type} [Lean.Grind.CommRing R]

/-- **Zoom is twist.** Shrinking two dial values by the natural weight w = r^d and
    combining by the circle operation is the twisted circle operation at k = w, scaled by w. -/
theorem zoom_is_twist (w t u : R) : circ (w * t) (w * u) = w * circK w t u :=
  (circK_scale w t u).symm

/-- **Zooms compose by multiplying:** zoom r₁ then r₂ weighs (r₁r₂)ⁿ = r₁ⁿr₂ⁿ
    (n a whole-number exponent; checks.py covers real s). -/
theorem zoom_compose (r₁ r₂ : R) (n : Nat) : (r₁ * r₂) ^ n = r₁ ^ n * r₂ ^ n :=
  mul_pow' r₁ r₂ n

/-- Zoom 1/2 in dimension β (whole β) weighs (1/2)^β, the size dial's 2^(−β). -/
theorem halving_weight (β : Nat) : ((1 : Rat) / 2) ^ β * 2 ^ β = 1 := by
  induction β with
  | zero => decide +kernel
  | succ b ih => rw [Rat.pow_succ, Rat.pow_succ]; grind

end Twist

end ZoomDial
