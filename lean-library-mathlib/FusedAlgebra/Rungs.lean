module

public import FusedAlgebra.Basic

/-!
# The fused hyperoperation ladder

Rung `n` of the fused ladder is the ordinary hyperoperation conjugated by `x ↦ x + 1`, on `ℕ`:

* rung 1: `sadd a b = a + b + 1`, the lift of addition;
* rung 2: `a ✳ b = a + b + a * b`, the lift of multiplication (mult-add);
* rung 3: `fexp a b = (a + 1) ^ (b + 1) - 1`, the fused exponent;
* rung 4: `ftet a b`, fused tetration;
* rung 5: `fpent a b`, fused pentation.

This is the mathlib form of `FusedAlgebra.Rungs`. The whole ladder is one definition,
`MultAdd.rung`, built on mathlib's `hyperoperation`; the named rungs are its first five values
(`rung_one`, …, `rung_five`), and "each rung is the rung below, iterated" is the single theorem
`rung_succ_succ`.
-/

@[expose] public section

namespace MultAdd

/-! ### The ladder -/

/-- Rung `n` of the fused ladder: the hyperoperation `n` seen through `x ↦ x + 1`. -/
def rung (n a b : ℕ) : ℕ := hyperoperation n (a + 1) (b + 1) - 1

/-- From rung 1 up, a hyperoperation of positive numbers is positive. -/
theorem hyperoperation_pos : ∀ n {m k : ℕ}, 0 < m → 0 < k → 0 < hyperoperation n m k
  | 0, m, k, _, _ => by simp [hyperoperation_zero]
  | 1, m, k, hm, _ => by rw [hyperoperation_one]; lia
  | 2, m, k, hm, hk => by rw [hyperoperation_two]; exact Nat.mul_pos hm hk
  | n + 3, m, k, hm, _ => by
    have h : ∀ j, 0 < hyperoperation (n + 3) m j := fun j ↦ by
      induction j with
      | zero => simp [hyperoperation_ge_three_eq_one]
      | succ j ih => rw [hyperoperation_recursion]; exact hyperoperation_pos (n + 2) hm ih
    exact h k

/-- The lift of the whole ladder: `rung n a b + 1` is the plain hyperoperation of the lifts. -/
theorem rung_add_one (n a b : ℕ) : rung n a b + 1 = hyperoperation n (a + 1) (b + 1) :=
  Nat.sub_add_cancel (hyperoperation_pos n a.succ_pos b.succ_pos)

/-- Each rung is the rung below, iterated. -/
theorem rung_succ_succ (n a b : ℕ) : rung (n + 1) a (b + 1) = rung n a (rung (n + 1) a b) := by
  apply Nat.add_right_cancel (m := 1)
  rw [rung_add_one, rung_add_one, rung_add_one, hyperoperation_recursion]

/-! ### The named rungs -/

/-- Rung 1, shifted addition: `a + b + 1`. -/
def sadd (a b : ℕ) : ℕ := a + b + 1

/-- Rung 3, the fused exponent: `(a + 1) ^ (b + 1) - 1`. -/
def fexp (a b : ℕ) : ℕ := (a + 1) ^ (b + 1) - 1

/-- Tetration: `tet x n` is a tower of `n` copies of `x`. This is `hyperoperation 4`. -/
def tet (x n : ℕ) : ℕ := hyperoperation 4 x n

/-- Rung 4, fused tetration. -/
def ftet (a b : ℕ) : ℕ := tet (a + 1) (b + 1) - 1

/-- Pentation: iterated tetration. This is `hyperoperation 5`. -/
def pent (x n : ℕ) : ℕ := hyperoperation 5 x n

/-- Rung 5, fused pentation. -/
def fpent (a b : ℕ) : ℕ := pent (a + 1) (b + 1) - 1

@[simp] theorem tet_zero (x : ℕ) : tet x 0 = 1 := hyperoperation_ge_three_eq_one 1 x

@[simp] theorem tet_succ (x n : ℕ) : tet x (n + 1) = x ^ tet x n := by
  rw [tet, hyperoperation_recursion, hyperoperation_three]; rfl

@[simp] theorem pent_zero (x : ℕ) : pent x 0 = 1 := hyperoperation_ge_three_eq_one 2 x

@[simp] theorem pent_succ (x n : ℕ) : pent x (n + 1) = tet x (pent x n) := by
  rw [pent, hyperoperation_recursion]; rfl

theorem rung_one (a b : ℕ) : rung 1 a b = sadd a b := by
  rw [rung, hyperoperation_one, sadd]; lia

theorem rung_two (a b : ℕ) : rung 2 a b = a ✳ b := by
  rw [rung, hyperoperation_two, ← multAdd_add_one, Nat.add_sub_cancel]

theorem rung_three (a b : ℕ) : rung 3 a b = fexp a b := by
  rw [rung, hyperoperation_three, fexp]

theorem rung_four (a b : ℕ) : rung 4 a b = ftet a b := rfl

theorem rung_five (a b : ℕ) : rung 5 a b = fpent a b := rfl

/-! ### Lifting: every rung is the plain operation seen through `x ↦ x + 1` -/

theorem sadd_add_one (a b : ℕ) : sadd a b + 1 = (a + 1) + (b + 1) := by
  rw [← rung_one, rung_add_one, hyperoperation_one]

theorem fexp_add_one (a b : ℕ) : fexp a b + 1 = (a + 1) ^ (b + 1) := by
  rw [← rung_three, rung_add_one, hyperoperation_three]

theorem tet_pos (a n : ℕ) : 0 < tet (a + 1) n := by
  cases n with
  | zero => simp
  | succ n => rw [tet_succ]; positivity

theorem ftet_add_one (a b : ℕ) : ftet a b + 1 = tet (a + 1) (b + 1) :=
  rung_add_one 4 a b

theorem pent_pos (a n : ℕ) : 0 < pent (a + 1) n := by
  cases n with
  | zero => simp
  | succ n => rw [pent_succ]; exact tet_pos a _

theorem fpent_add_one (a b : ℕ) : fpent a b + 1 = pent (a + 1) (b + 1) :=
  rung_add_one 5 a b

/-! ### Rung 3: the fused exponent -/

@[simp] theorem fexp_zero_right (a : ℕ) : fexp a 0 = a := by simp [fexp]

/-- `0` is a left absorber of the fused exponent. -/
@[simp] theorem fexp_zero_left (b : ℕ) : fexp 0 b = 0 := by simp [fexp]

/-- The fused exponent carries mult-add inside it. -/
theorem fexp_eq_multAdd (a b : ℕ) : fexp a b = a ✳ ((a + 1) ^ b - 1) := by
  apply Nat.add_right_cancel (m := 1)
  have h : (a + 1) ^ b - 1 + 1 = (a + 1) ^ b := Nat.sub_add_cancel (Nat.pow_pos a.succ_pos)
  rw [fexp_add_one, multAdd_add_one, h, pow_succ']

/-- Right distributivity over mult-add: `fexp (a ✳ b) c = fexp a c ✳ fexp b c`. -/
theorem fexp_multAdd_left (a b c : ℕ) : fexp (a ✳ b) c = fexp a c ✳ fexp b c := by
  apply Nat.add_right_cancel (m := 1)
  rw [fexp_add_one, multAdd_add_one, multAdd_add_one, fexp_add_one, fexp_add_one, mul_pow]

/-- Tower law: `fexp (fexp a b) c = fexp a (b ✳ c)`. -/
theorem fexp_fexp (a b c : ℕ) : fexp (fexp a b) c = fexp a (b ✳ c) := by
  apply Nat.add_right_cancel (m := 1)
  rw [fexp_add_one, fexp_add_one, fexp_add_one, multAdd_add_one, pow_mul]

/-- Sum law: `fexp a (sadd b c) = fexp a b ✳ fexp a c`. -/
theorem fexp_sadd (a b c : ℕ) : fexp a (sadd b c) = fexp a b ✳ fexp a c := by
  apply Nat.add_right_cancel (m := 1)
  rw [fexp_add_one, multAdd_add_one, fexp_add_one, fexp_add_one, sadd_add_one, ← pow_add]

theorem fexp_not_assoc : fexp (fexp 1 2) 1 ≠ fexp 1 (fexp 2 1) := by decide

theorem fexp_not_comm : fexp 1 2 ≠ fexp 2 1 := by decide

theorem fexp_example : fexp (fexp 1 2) 1 = 63 ∧ fexp 1 (fexp 2 1) = 511 := by decide

/-! ### Each rung is the rung below, iterated -/

theorem multAdd_succ (a b : ℕ) : a ✳ (b + 1) = sadd a (a ✳ b) := by
  rw [← rung_two, ← rung_two, ← rung_one, rung_succ_succ]

theorem fexp_succ (a b : ℕ) : fexp a (b + 1) = a ✳ fexp a b := by
  rw [← rung_three, ← rung_three, ← rung_two, rung_succ_succ]

/-- Rung 4 is rung 3 iterated on the exponent side. -/
theorem ftet_succ (a b : ℕ) : ftet a (b + 1) = fexp a (ftet a b) := by
  rw [← rung_four, ← rung_four, ← rung_three, rung_succ_succ]

/-- Rung 5 is rung 4 iterated on the exponent side. -/
theorem fpent_succ (a b : ℕ) : fpent a (b + 1) = ftet a (fpent a b) := by
  rw [← rung_five, ← rung_five, ← rung_four, rung_succ_succ]

/-- Rung 3 iterated on the base side with a fixed step `c` lifts to a double exponential, the
discrete flow behind `Y' = λ * Y * log Y`. -/
theorem iterate_fexp_add_one (a c n : ℕ) :
    (fun x ↦ fexp x c)^[n] a + 1 = (a + 1) ^ (c + 1) ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [Function.iterate_succ_apply', fexp_add_one, ih, ← pow_mul, pow_succ]

/-! ### Climbing higher: rungs 4 and 5, and which laws survive -/

theorem tet_one (x : ℕ) : tet x 1 = x := by simp

@[simp] theorem one_tet (n : ℕ) : tet 1 n = 1 := hyperoperation_ge_three_one 1 n

@[simp] theorem one_pent (n : ℕ) : pent 1 n = 1 := hyperoperation_ge_three_one 2 n

@[simp] theorem ftet_zero_right (a : ℕ) : ftet a 0 = a := by simp [ftet]

@[simp] theorem ftet_zero_left (b : ℕ) : ftet 0 b = 0 := by simp [ftet]

@[simp] theorem fpent_zero_right (a : ℕ) : fpent a 0 = a := by simp [fpent]

@[simp] theorem fpent_zero_left (b : ℕ) : fpent 0 b = 0 := by simp [fpent]

theorem ftet_values :
    ftet 1 1 = 3 ∧ ftet 1 2 = 15 ∧ ftet 2 1 = 26 ∧ ftet 3 1 = 255 ∧ ftet 1 3 = 65535 := by
  simp [ftet]

theorem fpent_values : fpent 1 1 = 3 ∧ fpent 1 2 = 65535 := by
  simp [fpent]

theorem ftet_not_assoc : ftet (ftet 1 1) 1 ≠ ftet 1 (ftet 1 1) := by
  simp [ftet]

theorem ftet_not_comm : ftet 1 2 ≠ ftet 2 1 := by
  simp [ftet]

/-- The rung-3 laws stop at rung 4: distributivity over mult-add fails ... -/
theorem ftet_multAdd_left_ne : ftet (1 ✳ 1) 1 ≠ ftet 1 1 ✳ ftet 1 1 := by
  simp [ftet, multAdd]

/-- ... the tower law fails ... -/
theorem ftet_ftet_ne : ftet (ftet 1 1) 1 ≠ ftet 1 (fexp 1 1) := by
  simp [ftet, fexp]

/-- ... and the sum law fails. -/
theorem ftet_sadd_ne : ftet 1 (sadd 1 1) ≠ fexp (ftet 1 1) (ftet 1 1) := by
  simp [ftet, fexp, sadd]

end MultAdd
