/-
  Natural constants up the ladder of Fused Hyperoperation Algebra
  ===============================================================

  Companion page: ../index.html.  Numbers: ../constants.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`.  Run it with
  `lean FusedConstants.lean`.  Checked with Lean 4.37.0-rc1.

  Rung r+1 of the ordinary ladder is `H r`:
    H 0 = add, H 1 = multiply, H 2 = power, H 3 = tetration, ...
  Fused rung r+1 is the same operation slid by one:  a (r) b = H r (1+a) (1+b) - 1.

  Proved:
  * mult_add:        (1+a)(1+b) - 1 = a + b + ab, so fused rung 2 is the circle operation (called mult-add before 1.0).
  * slide_fixed, slide_gap, slide_step:  sliding a map by one (x -> 1+x) moves its
                     fixed points by exactly one and keeps every gap to the diagonal
                     and every step size.  So a kiss with the diagonal survives the
                     slide unchanged, and each fused constant is the ordinary one minus 1.
  * rung1_base0, rung2_base1: the first two constants are exact.  Adding 0 and
                     multiplying by 1 ARE the diagonal (every point is a kiss).
  * all_rungs_two_two: 2 (r) 2 = 4 on every rung; fused: 1 (r) 1 = 3 on every rung.
  * base_two_clears: on every rung from power up, H r 2 n > n for every whole n:
                     the base-2 curve never touches the diagonal at a whole height.
  * base_two_tower_climbs: so the base-2 tower at every rung from tetration up
                     climbs by at least one at each step and never settles.
                     This is the whole-number skeleton of "every constant is below 2".

  What Lean does NOT do here: core Lean has no real numbers, so the real-height
  extension, the kiss itself (value = x and slope = 1), e^(1/e), the 1.6353 base,
  and the limit "the constants climb to 2" are checked numerically in
  constants.py and argued on the page, not in Lean.
-/

namespace FusedConstants

/-- `H r` is rung r+1 of the ordinary ladder. -/
def H : Nat → Nat → Nat → Nat
  | 0, a, b => a + b
  | 1, a, b => a * b
  | 2, a, b => a ^ b
  | _ + 3, _, 0 => 1
  | k + 3, a, b + 1 => H (k + 2) a (H (k + 3) a b)
termination_by r _ b => (r, b)

theorem H_succ (k a b : Nat) : H (k + 3) a (b + 1) = H (k + 2) a (H (k + 3) a b) := by
  rw [H]

theorem H_zero (k a : Nat) : H (k + 3) a 0 = 1 := by
  rw [H]

theorem H_add (a b : Nat) : H 0 a b = a + b := by rw [H]
theorem H_mul (a b : Nat) : H 1 a b = a * b := by rw [H]
theorem H_pow (a b : Nat) : H 2 a b = a ^ b := by rw [H]

/-! ### Fused rung 2 is the circle operation -/

theorem mult_add (a b : Int) : (1 + a) * (1 + b) - 1 = a + b + a * b := by
  rw [Int.add_mul, Int.one_mul, Int.mul_add, Int.mul_one]
  omega

/-! ### Sliding by one keeps every kiss

`f` is any map on whole numbers (think: x ↦ H r b x).  Its slide is
`slide f y = f (1 + y) - 1`, the same curve drawn with the unit corner as origin. -/

def slide (f : Int → Int) (y : Int) : Int := f (1 + y) - 1

/-- Fixed points move by exactly one. -/
theorem slide_fixed (f : Int → Int) (x : Int) : f x = x ↔ slide f (x - 1) = x - 1 := by
  unfold slide
  have : 1 + (x - 1) = x := by omega
  rw [this]
  omega

/-- The gap between the curve and the diagonal is the same before and after. -/
theorem slide_gap (f : Int → Int) (y : Int) : slide f y - y = f (1 + y) - (1 + y) := by
  unfold slide; omega

/-- Every step (the whole-number stand-in for a slope) is the same before and after. -/
theorem slide_step (f : Int → Int) (y : Int) :
    slide f (y + 1) - slide f y = f ((1 + y) + 1) - f (1 + y) := by
  unfold slide
  have : 1 + (y + 1) = (1 + y) + 1 := by omega
  rw [this]; omega

/-- A slid kiss: if the original curve stays on or above the diagonal and touches it
    at x, the slid curve stays on or above the diagonal and touches it at x - 1. -/
theorem slide_kiss (f : Int → Int) (x : Int) (above : ∀ t, t ≤ f t) (touch : f x = x) :
    (∀ y, y ≤ slide f y) ∧ slide f (x - 1) = x - 1 := by
  refine ⟨fun y => ?_, (slide_fixed f x).1 touch⟩
  have := above (1 + y)
  unfold slide; omega

/-! ### The first two constants are exact -/

/-- Rung 1, base 0: adding 0 is the diagonal itself. -/
theorem rung1_base0 (x : Nat) : H 0 0 x = x := by
  rw [H_add]; omega

/-- Rung 2, base 1: multiplying by 1 is the diagonal itself. -/
theorem rung2_base1 (x : Nat) : H 1 1 x = x := by
  rw [H_mul]; omega

/-- Any base above 0 on rung 1 clears the diagonal everywhere. -/
theorem rung1_above (b x : Nat) (hb : 0 < b) : x < H 0 b x := by
  rw [H_add]; omega

/-! ### Every rung agrees at 2 -/

theorem H_two_one : ∀ r, 1 ≤ r → H r 2 1 = 2
  | 0, h => absurd h (by decide)
  | 1, _ => by rw [H_mul]
  | 2, _ => by rw [H_pow]
  | k + 3, _ => by
      rw [H_succ, H_zero]
      exact H_two_one (k + 2) (by omega)

/-- 2 (r) 2 = 4 on every rung. -/
theorem all_rungs_two_two : ∀ r, H r 2 2 = 4
  | 0 => by rw [H_add]
  | 1 => by rw [H_mul]
  | 2 => by rw [H_pow]
  | k + 3 => by
      rw [H_succ, H_two_one (k + 3) (by omega)]
      exact all_rungs_two_two (k + 2)

/-- Fused form: 1 (r) 1 = 3 on every rung. -/
theorem fused_one_one (r : Nat) : H r (1 + 1) (1 + 1) - 1 = 3 := by
  rw [show 1 + 1 = 2 from rfl, all_rungs_two_two]

/-! ### Base 2 never touches the diagonal at a whole height -/

theorem lt_two_pow' : ∀ n : Nat, n < 2 ^ n
  | 0 => by decide
  | n + 1 => by
      have ih := lt_two_pow' n
      rw [Nat.pow_succ]; omega

/-- On every rung from power up (r ≥ 2), the base-2 curve is above the diagonal. -/
theorem base_two_clears : ∀ k n : Nat, n < H (k + 2) 2 n
  | 0, n => by rw [show 0 + 2 = 2 from rfl, H_pow]; exact lt_two_pow' n
  | k + 1, 0 => by
      show 0 < H (k + 3) 2 0
      rw [H_zero]; decide
  | k + 1, n + 1 => by
      show n + 1 < H (k + 3) 2 (n + 1)
      rw [H_succ]
      have ih_n := base_two_clears (k + 1) n
      have ih_k := base_two_clears k (H (k + 3) 2 n)
      change n < H (k + 3) 2 n at ih_n
      omega

/-- The base-2 tower on every rung from tetration up climbs at every step
    (tower step n+1 = previous rung applied to tower step n), so it never settles. -/
theorem base_two_tower_climbs (k n : Nat) : H (k + 3) 2 n < H (k + 3) 2 (n + 1) := by
  rw [H_succ]
  exact base_two_clears k (H (k + 3) 2 n)

theorem base_two_tower_unbounded (k n : Nat) : n < H (k + 3) 2 n := base_two_clears (k + 1) n

/-! ### Kernel checks of small cases -/

example : H 3 2 3 = 16 := by                  -- 2^^3
  simp [H_succ, H_zero, H_pow]
example : H 4 2 3 = 65536 := by               -- 2^^^3 = 2^^4
  simp [H_succ, H_zero, H_pow]
example : H 3 3 2 = 27 := by
  simp [H_succ, H_zero, H_pow]
example : (List.range 6).map (fun r => H r 2 2) = [4, 4, 4, 4, 4, 4] := by
  simp only [List.range, List.range.loop, List.map, all_rungs_two_two]

end FusedConstants
