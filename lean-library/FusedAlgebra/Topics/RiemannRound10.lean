import FusedAlgebra.Base

/-!
# Round 10 (DRAFT): small exact facts behind the prime lifecycle and the hyper integrals

Core Lean only (`decide`, `omega`, `grind`).
- `wheel6`, `wheel30`, `wheel210`: the open spokes of the primorial wheel (2, 8, 48).
- `jet2`, `jet3`: one step of the rung-2 fold, (1 + k a)(1 + k b)(1 + k c) = 1 + k e1 + k² e2 + k³ e3.
  This is det(I + hS) of the jet-river page and the layer expansion of ∏(1 + k p^-s).
- `fold_at_minus_one`: at k = -1 the layers alternate (the Möbius signs).
- `newton2`: Newton's inequality for two layers, 4 e2 ≤ e1².
- `rung3_log_mult`: rung 3 is multiplication read in log coordinates (exponents add).

First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.RiemannRound10
(namespace CircleOp.RiemannRound10, the library's convention). One proof changed:
`newton2` gets 0 ≤ (a − b)² from Int.mul_nonneg by cases on the sign, since core Lean has no
`Int.mul_self_nonneg`; the statement is unchanged.
-/

namespace CircleOp.RiemannRound10

set_option maxRecDepth 8000

/-- Spokes of the wheel mod 6 left open by 2 and 3. -/
theorem wheel6 : (List.range 6).filter (fun n => n % 2 ≠ 0 ∧ n % 3 ≠ 0) = [1, 5] := by decide

/-- Mod 30 the wheel of 2, 3, 5 leaves 8 spokes: 30·(1/2)(2/3)(4/5). -/
theorem wheel30 :
    ((List.range 30).filter (fun n => n % 2 ≠ 0 ∧ n % 3 ≠ 0 ∧ n % 5 ≠ 0)).length = 8 := by decide

/-- Mod 210 the wheel of 2, 3, 5, 7 leaves 48 spokes. -/
theorem wheel210 :
    ((List.range 210).filter
      (fun n => n % 2 ≠ 0 ∧ n % 3 ≠ 0 ∧ n % 5 ≠ 0 ∧ n % 7 ≠ 0)).length = 48 := by decide

/-- Two-factor rung-2 fold: the layers 1, e1, e2. -/
theorem jet2 (k a b : Int) :
    (1 + k * a) * (1 + k * b) = 1 + k * (a + b) + k * k * (a * b) := by grind

/-- Three-factor rung-2 fold: det(I + hS) = 1 + hP + h²Q + h³R for diagonal S. -/
theorem jet3 (k a b c : Int) :
    (1 + k * a) * (1 + k * b) * (1 + k * c)
      = 1 + k * (a + b + c) + k * k * (a * b + b * c + c * a) + k * k * k * (a * b * c) := by grind

/-- At k = -1 the layers alternate in sign (the Möbius signs of 1/ζ). -/
theorem fold_at_minus_one (a b c : Int) :
    (1 - a) * (1 - b) * (1 - c) = 1 - (a + b + c) + (a * b + b * c + c * a) - a * b * c := by grind

/-- Newton's inequality for two variables: 4 e2 ≤ e1². -/
theorem newton2 (a b : Int) : 4 * (a * b) ≤ (a + b) * (a + b) := by
  have h : 0 ≤ (a - b) * (a - b) := by
    rcases Int.le_total 0 (a - b) with hx | hx
    · exact Int.mul_nonneg hx hx
    · have h' := Int.mul_nonneg (Int.neg_nonneg_of_nonpos hx) (Int.neg_nonneg_of_nonpos hx)
      rw [Int.neg_mul_neg] at h'
      exact h'
  grind

/-- Rung 3 read in log coordinates: exponents add. -/
theorem rung3_log_mult (a m n : Nat) : a ^ (m + n) = a ^ m * a ^ n := Nat.pow_add a m n

end CircleOp.RiemannRound10
