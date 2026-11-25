/-
  FusedAlgebra.Rungs: the fused hyperoperation ladder
  ===================================================

  Rung n is the ordinary hyperoperation conjugated by x ↦ 1 + x, on Nat:
    rung 1  sadd a b = a + b + 1          = (a+1) + (b+1) − 1
    rung 2  a ⊞ b    = a + b + a·b        = (a+1)·(b+1) − 1   (mult-add, in Base)
    rung 3  fexp a b = (a+1)^(b+1) − 1                          (fused exponent)
    rung 4  ftet a b = (a+1)↑↑(b+1) − 1                         (fused tetration)
    rung 5  fpent a b                                           (fused pentation)

  From ns-rungs/lean/FusedRungs.lean. Core Lean only, no `sorry`.
-/
import FusedAlgebra.Base

namespace CircleOp

def sadd (a b : Nat) : Nat := a + b + 1
def fexp (a b : Nat) : Nat := (a + 1) ^ (b + 1) - 1
def tet (x : Nat) : Nat → Nat
  | 0 => 1
  | n + 1 => x ^ tet x n
def ftet (a b : Nat) : Nat := tet (a + 1) (b + 1) - 1

/-! ### Lifting: every rung is the plain operation seen through x ↦ 1 + x -/
theorem sadd_lift (a b : Nat) : sadd a b + 1 = (a + 1) + (b + 1) := by unfold sadd; omega
theorem pow_pos' (a n : Nat) : 1 ≤ (a + 1) ^ n := Nat.pow_pos (Nat.succ_pos a)
theorem fexp_lift (a b : Nat) : fexp a b + 1 = (a + 1) ^ (b + 1) := by
  unfold fexp; have := pow_pos' a (b + 1); omega
theorem tet_pos (a n : Nat) : 1 ≤ tet (a + 1) n := by
  cases n with
  | zero => simp [tet]
  | succ n => simp only [tet]; exact pow_pos' a _
theorem ftet_lift (a b : Nat) : ftet a b + 1 = tet (a + 1) (b + 1) := by
  unfold ftet; have := tet_pos a (b + 1); omega
theorem lift_inj {x y : Nat} (h : x + 1 = y + 1) : x = y := by omega

/-! ### Rung 2 (mult-add) lives in FusedAlgebra.Base: commutative, associative, identity 0 -/

/-! ### Rung 3 (fused exponent) -/
theorem fexp_zero_right (a : Nat) : fexp a 0 = a := by unfold fexp; simp
/-- 0 is a left absorber, so there is no left identity of the form 0. -/
theorem fexp_zero_left (b : Nat) : fexp 0 b = 0 := by unfold fexp; simp
/-- The fused exponent carries mult-add inside it. -/
theorem fexp_carries_circ (a b : Nat) : fexp a b = circ a ((a + 1) ^ b - 1) := by
  apply lift_inj
  rw [fexp_lift, circ_lift']
  have := pow_pos' a b
  rw [Nat.sub_add_cancel this, Nat.pow_succ, Nat.mul_comm]
/-- Right distributivity over mult-add:  (a ⊕ b) ⊛ c = (a ⊛ c) ⊕ (b ⊛ c). -/
theorem fexp_circ_right (a b c : Nat) : fexp (circ a b) c = circ (fexp a c) (fexp b c) := by
  apply lift_inj
  rw [fexp_lift, circ_lift', circ_lift', fexp_lift, fexp_lift, Nat.mul_pow]
/-- Tower law:  (a ⊛ b) ⊛ c = a ⊛ (b ⊕ c). -/
theorem fexp_tower (a b c : Nat) : fexp (fexp a b) c = fexp a (circ b c) := by
  apply lift_inj
  rw [fexp_lift, fexp_lift, fexp_lift, circ_lift', Nat.pow_mul]
/-- Sum law:  a ⊛ (b ⊞ c) = (a ⊛ b) ⊕ (a ⊛ c), with ⊞ the rung-1 shifted add. -/
theorem fexp_sadd (a b c : Nat) : fexp a (sadd b c) = circ (fexp a b) (fexp a c) := by
  apply lift_inj
  rw [fexp_lift, circ_lift', fexp_lift, fexp_lift, sadd_lift, ← Nat.pow_add]
theorem fexp_not_assoc : fexp (fexp 1 2) 1 ≠ fexp 1 (fexp 2 1) := by decide
theorem fexp_not_comm : fexp 1 2 ≠ fexp 2 1 := by decide
theorem fexp_example : fexp (fexp 1 2) 1 = 63 ∧ fexp 1 (fexp 2 1) = 511 := by decide

/-! ### Each rung is the rung below, iterated -/
theorem circ_iter (a b : Nat) : circ a (b + 1) = sadd a (circ a b) := by
  apply lift_inj; rw [circ_lift', sadd_lift, circ_lift']; grind
theorem fexp_iter (a b : Nat) : fexp a (b + 1) = circ a (fexp a b) := by
  apply lift_inj; rw [fexp_lift, circ_lift', fexp_lift, Nat.pow_succ, Nat.mul_comm]
/-- Rung 4 is rung 3 iterated on the exponent side. -/
theorem ftet_iter (a b : Nat) : ftet a (b + 1) = fexp a (ftet a b) := by
  apply lift_inj; rw [ftet_lift, fexp_lift, ftet_lift]; rfl
/-- Rung 3 iterated on the base side with a fixed step c: lifts to a double exponential,
    the discrete flow behind Y' = λ·Y·log Y. -/
theorem fexp_base_iter (a c n : Nat) :
    (Nat.repeat (fun x => fexp x c) n a) + 1 = (a + 1) ^ ((c + 1) ^ n) := by
  induction n with
  | zero => simp [Nat.repeat]
  | succ n ih =>
    simp only [Nat.repeat]
    rw [fexp_lift, ih, ← Nat.pow_mul, Nat.pow_succ]

/-! ### Climbing higher: rung 5 (fused pentation), and which laws survive at rung 4 -/
def pent (x : Nat) : Nat → Nat
  | 0 => 1
  | n + 1 => tet x (pent x n)
def fpent (a b : Nat) : Nat := pent (a + 1) (b + 1) - 1
theorem pent_pos (a n : Nat) : 1 ≤ pent (a + 1) n := by
  cases n with
  | zero => simp [pent]
  | succ n => simp only [pent]; exact tet_pos a _
theorem fpent_lift (a b : Nat) : fpent a b + 1 = pent (a + 1) (b + 1) := by
  unfold fpent; have := pent_pos a (b + 1); omega
theorem tet_one (x : Nat) : tet x 1 = x := by simp [tet]
theorem ftet_zero_right (a : Nat) : ftet a 0 = a := by unfold ftet; simp [tet]
theorem ftet_zero_left (b : Nat) : ftet 0 b = 0 := by unfold ftet; simp [tet]
theorem fpent_zero_right (a : Nat) : fpent a 0 = a := by unfold fpent; simp [pent, tet]
theorem tet_one_base (n : Nat) : tet 1 n = 1 := by
  cases n with
  | zero => rfl
  | succ n => simp [tet]
theorem pent_one_base (n : Nat) : pent 1 n = 1 := by
  cases n with
  | zero => rfl
  | succ n => simp [pent, tet_one_base]
theorem fpent_zero_left (b : Nat) : fpent 0 b = 0 := by unfold fpent; simp [pent_one_base]
/-- Rung 5 is rung 4 iterated on the exponent side. -/
theorem fpent_iter (a b : Nat) : fpent a (b + 1) = ftet a (fpent a b) := by
  apply lift_inj; rw [fpent_lift, ftet_lift, fpent_lift]; simp only [pent]
theorem ftet_values : ftet 1 1 = 3 ∧ ftet 1 2 = 15 ∧ ftet 2 1 = 26 ∧ ftet 3 1 = 255 ∧ ftet 1 3 = 65535 := by
  decide
theorem fpent_values : fpent 1 1 = 3 ∧ fpent 1 2 = 65535 := by decide
theorem ftet_not_assoc : ftet (ftet 1 1) 1 ≠ ftet 1 (ftet 1 1) := by decide
theorem ftet_not_comm : ftet 1 2 ≠ ftet 2 1 := by decide
/-- The rung-3 laws stop at rung 4: distributivity over mult-add fails ... -/
theorem ftet_no_circ_dist : ftet (circ 1 1) 1 ≠ circ (ftet 1 1) (ftet 1 1) := by decide
/-- ... the tower law fails ... -/
theorem ftet_no_tower : ftet (ftet 1 1) 1 ≠ ftet 1 (fexp 1 1) := by decide
/-- ... and the sum law fails. -/
theorem ftet_no_sum : ftet 1 (sadd 1 1) ≠ fexp (ftet 1 1) (ftet 1 1) := by decide


end CircleOp
