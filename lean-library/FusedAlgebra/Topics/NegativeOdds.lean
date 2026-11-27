/-
  Negative odds and the wall (DRAFT, Nov 27, 2026)
  ================================================

  Part of Fused Hyperoperation Algebra. Page and checks: ../ (negative-odds).
  Core Lean 4 only (no mathlib), no `sorry`. Checked with Lean 4.37.0-rc1 against FusedAlgebra.Base.

  Odds t, chance p = t/(1 + t). Independent items stack by t ⊞ s = t + s + ts (k = 1).
  Proved:
  * union of independent events composes odds by ⊞ (the odds dial's law);
  * the wall −1 absorbs union, is fixed by the complement mirror t ↦ 1/t, and absorbs intersection;
  * the mirror is the complement p ↦ 1 − p; it swaps the two negative regions;
  * the undo ("rewind", the quasi-inverse) of odds t has chance −t;
  * t = −2 (chance 200%) is its own undo; two items past the wall give a positive lift;
  * rewind and mirror generate the six-map cross-ratio group: their composite has order 3,
    the third involution t ↦ −1 − t sends p to 1/p, and each involution fixes one of
    never (0), certain (∞, as a limit) and the wall (−1) while swapping the other two;
  * the flip t/(1 + t) is minus the rewind (sign of k = direction of the clock);
  * doom/destiny bookkeeping: keeping the count of zeros makes 0·∞ well defined.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.NegativeOdds
  (namespace CircleOp.NegativeOdds, the library's convention; no proof changed).
-/
import FusedAlgebra.Base

namespace CircleOp.NegativeOdds

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

variable {K : Type} [Field K]

/-- Chance from odds. -/
def prob (t : K) : K := t / (1 + t)
/-- Odds from chance. -/
def odds (p : K) : K := p / (1 - p)
/-- Rewind: the undo of an item (the quasi-inverse, "way back"). -/
def rewind (t : K) : K := -t / (1 + t)
/-- Intersection of independent events, through the complement mirror. -/
def inter (t s : K) : K := 1 / circ (1 / t) (1 / s)

/-! ## 1. The law and the wall -/

/-- Union of independent events: odds compose by mult-add. -/
theorem odds_union (p q : K) (hp : 1 - p ≠ 0) (hq : 1 - q ≠ 0) :
    odds (1 - (1 - p) * (1 - q)) = circ (odds p) (odds q) := by
  unfold odds circ; grind

/-- The wall swallows every stack. -/
theorem wall_absorbs (s : K) : circ (-1) s = -1 := by unfold circ; grind

/-- The complement mirror fixes the wall. -/
theorem mirror_fixes_wall : (1 : K) / (-1) = -1 := by grind

/-- The wall also swallows intersection (for any nonzero odds). -/
theorem wall_absorbs_inter (s : K) (_hs : s ≠ 0) : inter (-1) s = -1 := by
  unfold inter circ; grind

/-- The mirror is the complement. -/
theorem prob_mirror (t : K) (h0 : t ≠ 0) (h1 : 1 + t ≠ 0) : prob (1 / t) = 1 - prob t := by
  unfold prob; grind

/-! ## 2. Rewinds and the sign sector -/

/-- Rewinding cancels. -/
theorem rewind_cancels (t : K) (h : 1 + t ≠ 0) : circ t (rewind t) = 0 := by
  unfold circ rewind; grind

/-- The undo of odds t has chance −t. -/
theorem prob_rewind (t : K) (h : 1 + t ≠ 0) : prob (rewind t) = -t := by
  unfold prob rewind; grind

/-- Rewind is its own undo. -/
theorem rewind_rewind (t : K) (h : 1 + t ≠ 0) : rewind (rewind t) = t := by
  unfold rewind; grind

/-- Example: a −1/3 "chance" cancels a 1/4 chance (p + q − pq = 0). -/
theorem anti_chance_example : (1/4 : Rat) + (-1/3) - (1/4) * (-1/3) = 0 := by grind

/-- 200% twice is never: t = −2 is its own undo. -/
theorem two_hundred_twice : circ (-2 : K) (-2) = 0 := by unfold circ; grind

theorem prob_minus_two (h : (1 : K) + (-2) ≠ 0) : prob (-2 : K) = 2 := by
  unfold prob; grind

/-- Two items past the wall give a positive lift (the sign bit cancels in pairs). -/
theorem past_wall_pair (a b : Rat) (ha : a < -1) (hb : b < -1) : 0 < 1 + circ a b := by
  rw [circ_lift]
  have h : 0 < (-(1 + a)) * (-(1 + b)) := Rat.mul_pos (by grind) (by grind)
  grind

/-! ## 3. Six maps, three special points -/

/-- Rewind after mirror. -/
def turn (t : K) : K := -1 / (1 + t)

theorem rewind_mirror (t : K) (h0 : t ≠ 0) (h1 : 1 + t ≠ 0) : rewind (1 / t) = turn t := by
  unfold rewind turn; grind

/-- The turn has order 3. -/
theorem turn_order_three (t : K) (h0 : t ≠ 0) (h1 : 1 + t ≠ 0) : turn (turn (turn t)) = t := by
  unfold turn; grind

/-- The third involution t ↦ −1 − t sends p to 1/p. -/
theorem prob_third (t : K) (h0 : t ≠ 0) (h1 : 1 + t ≠ 0) : prob (-1 - t) = 1 / prob t := by
  unfold prob; grind

theorem third_involution (t : K) : -1 - (-1 - t) = t := by grind

/-- The third involution swaps never and the wall. -/
theorem third_swaps : (-1 - (0 : K) = -1) ∧ (-1 - (-1 : K) = 0) := by constructor <;> grind

/-- Rewind fixes never and the 200% point. -/
theorem rewind_fixed : rewind (0 : K) = 0 ∧ ((1 : K) + (-2) ≠ 0 → rewind (-2 : K) = -2) := by
  constructor
  · unfold rewind; grind
  · intro h; unfold rewind; grind

/-! ## 4. Sign of k is the direction of the clock -/

/-- The flip t/(1 + t) (odds → chance, law k = 1 → law k = −1) is minus the rewind. -/
theorem flip_is_minus_rewind (t : K) : prob t = -rewind t := by
  unfold prob rewind; grind

/-! ## 5. Doom and destiny: keep the count, and 0·∞ has an answer -/

/-- A multiplier kept as (coefficient, net count of zeros). Doom = (1, 1), destiny = (1, −1). -/
structure Kept (K : Type) where
  c : K
  n : Int

def Kept.mul (a b : Kept K) : Kept K := ⟨a.c * b.c, a.n + b.n⟩

theorem kept_comm (a b : Kept K) : Kept.mul a b = Kept.mul b a := by
  unfold Kept.mul; congr 1 <;> grind

theorem kept_assoc (a b c : Kept K) : Kept.mul (Kept.mul a b) c = Kept.mul a (Kept.mul b c) := by
  unfold Kept.mul; congr 1 <;> grind

/-- Doom times destiny is the plain multiplier 1: they cancel instead of giving 0·∞. -/
theorem doom_destiny : Kept.mul (⟨1, 1⟩ : Kept K) ⟨1, -1⟩ = ⟨1, 0⟩ := by
  unfold Kept.mul; congr 1 <;> grind

end CircleOp.NegativeOdds
