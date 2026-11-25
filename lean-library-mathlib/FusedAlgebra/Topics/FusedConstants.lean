module

public import FusedAlgebra.Rungs

/-!
# Natural constants up the fused ladder: the whole-number skeleton

Rung `r + 1` of the ordinary ladder is `H r` (`H 0` is addition, `H 1` multiplication, `H 2`
power, `H 3` tetration, …). In mathlib this is `hyperoperation (r + 1)`, and `H` is defined as
that. Fused rung `r + 1` is the same operation slid by one, `MultAdd.rung (r + 1) a b =
H r (a + 1) (b + 1) - 1` (`rung_succ_eq`).

## Main results

* `slide`, `slide_fixed`, `slide_gap`, `slide_step`, `slide_kiss`: sliding a map by one
  (`x ↦ 1 + x`) moves its fixed points by exactly one and keeps every gap to the diagonal and
  every step, so a kiss with the diagonal survives the slide: each fused constant is the ordinary
  one minus `1`. Stated in any ring (the original was over `ℤ`).
* `rung1_base0`, `rung2_base1`: adding `0` and multiplying by `1` are the diagonal itself.
* `all_rungs_two_two`, `rung_succ_one_one`: `2 (r) 2 = 4` on every rung, so fused `1 (r) 1 = 3`.
  The first is mathlib's `hyperoperation_two_two_eq_four`.
* `base_two_clears`, `base_two_tower_climbs`: from the power rung up, the base-`2` curve is
  strictly above the diagonal at every whole height, so the base-`2` tower climbs at every step.
  This is the whole-number skeleton of "every constant is below `2`".

The real-height statements (the kiss itself, `e ^ (1 / e)`, the limit of the constants) are not
formalised here, as in the core-Lean original.
-/

@[expose] public section

namespace FusedConstants

open MultAdd

/-! ### The ordinary ladder -/

/-- `H r` is rung `r + 1` of the ordinary ladder: mathlib's `hyperoperation (r + 1)`. -/
def H (r a b : ℕ) : ℕ := hyperoperation (r + 1) a b

theorem H_add (a b : ℕ) : H 0 a b = a + b := by rw [H, hyperoperation_one]

theorem H_mul (a b : ℕ) : H 1 a b = a * b := by rw [H, hyperoperation_two]

theorem H_pow (a b : ℕ) : H 2 a b = a ^ b := by rw [H, hyperoperation_three]

theorem H_zero (k a : ℕ) : H (k + 3) a 0 = 1 := hyperoperation_ge_three_eq_one (k + 1) a

theorem H_succ (k a b : ℕ) : H (k + 3) a (b + 1) = H (k + 2) a (H (k + 3) a b) :=
  hyperoperation_recursion (k + 3) a b

/-- Fused rung `r + 1` is `H r` slid by one. -/
theorem rung_succ_eq (r a b : ℕ) : rung (r + 1) a b = H r (a + 1) (b + 1) - 1 := rfl

/-! ### Sliding by one keeps every kiss

`f` is any map (think `x ↦ H r b x`). Its slide is `slide f y = f (1 + y) - 1`, the same curve
drawn with the unit corner as origin. -/

section Slide

variable {R : Type*} [Ring R]

/-- The slide of a map by one: `slide f y = f (1 + y) - 1`. -/
def slide (f : R → R) (y : R) : R := f (1 + y) - 1

/-- Fixed points move by exactly one. -/
theorem slide_fixed (f : R → R) (x : R) : f x = x ↔ slide f (x - 1) = x - 1 := by
  rw [slide, add_sub_cancel, sub_left_inj]

/-- The gap between the curve and the diagonal is the same before and after. -/
theorem slide_gap (f : R → R) (y : R) : slide f y - y = f (1 + y) - (1 + y) := by
  rw [slide]; abel

/-- Every step (the whole-number stand-in for a slope) is the same before and after. -/
theorem slide_step (f : R → R) (y : R) :
    slide f (y + 1) - slide f y = f (1 + y + 1) - f (1 + y) := by
  rw [slide, slide, add_assoc]; abel

/-- A slid kiss: if the original curve stays on or above the diagonal and touches it at `x`, the
slid curve stays on or above the diagonal and touches it at `x - 1`. -/
theorem slide_kiss [PartialOrder R] [IsOrderedAddMonoid R] (f : R → R) (x : R)
    (above : ∀ t, t ≤ f t) (touch : f x = x) :
    (∀ y, y ≤ slide f y) ∧ slide f (x - 1) = x - 1 :=
  ⟨fun y ↦ le_sub_iff_add_le'.2 (above (1 + y)), (slide_fixed f x).1 touch⟩

end Slide

/-! ### The first two constants are exact -/

/-- Rung 1, base 0: adding `0` is the diagonal itself. -/
theorem rung1_base0 (x : ℕ) : H 0 0 x = x := by rw [H_add, zero_add]

/-- Rung 2, base 1: multiplying by `1` is the diagonal itself. -/
theorem rung2_base1 (x : ℕ) : H 1 1 x = x := by rw [H_mul, one_mul]

/-- Any base above `0` on rung 1 clears the diagonal everywhere. -/
theorem rung1_above (b x : ℕ) (hb : 0 < b) : x < H 0 b x := by rw [H_add]; lia

/-! ### Every rung agrees at 2 -/

/-- `2 (r) 2 = 4` on every rung. This is mathlib's `hyperoperation_two_two_eq_four`. -/
theorem all_rungs_two_two (r : ℕ) : H r 2 2 = 4 := hyperoperation_two_two_eq_four r

/-- Fused form: `1 (r) 1 = 3` on every rung. -/
theorem rung_succ_one_one (r : ℕ) : rung (r + 1) 1 1 = 3 := by
  rw [rung_succ_eq, all_rungs_two_two]

/-! ### Base 2 never touches the diagonal at a whole height -/

/-- On every rung from power up, the base-`2` curve is above the diagonal. -/
theorem base_two_clears (k n : ℕ) : n < H (k + 2) 2 n := by
  induction k generalizing n with
  | zero => rw [H_pow]; exact Nat.lt_two_pow_self
  | succ k ihk =>
    induction n with
    | zero => rw [H_zero]; exact Nat.one_pos
    | succ n ihn =>
      rw [H_succ]
      exact lt_of_le_of_lt (Nat.succ_le_of_lt ihn) (ihk _)

/-- The base-`2` tower on every rung from tetration up climbs at every step, so it never
settles. -/
theorem base_two_tower_climbs (k n : ℕ) : H (k + 3) 2 n < H (k + 3) 2 (n + 1) := by
  rw [H_succ]; exact base_two_clears k _

theorem base_two_tower_unbounded (k n : ℕ) : n < H (k + 3) 2 n := base_two_clears (k + 1) n

/-! ### Small cases -/

example : H 3 2 3 = 16 := by simp [H_succ, H_zero, H_pow]

example : H 4 2 3 = 65536 := by simp [H_succ, H_zero, H_pow]

example : H 3 3 2 = 27 := by simp [H_succ, H_zero, H_pow]

example : (List.range 6).map (fun r ↦ H r 2 2) = [4, 4, 4, 4, 4, 4] := by
  simp [all_rungs_two_two]

end FusedConstants
