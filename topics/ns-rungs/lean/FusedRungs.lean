/-
  Fused Hyperoperation Algebra: the rungs, and the algebra used on the Navier–Stokes page.
  Core Lean 4 only (no mathlib), no `sorry`.

  Rung n = ordinary hyperoperation conjugated by x ↦ 1 + x, worked on Nat (so 1 + x ≥ 1):
    rung 1  sadd a b = a + b + 1          = (a+1) + (b+1) - 1
    rung 2  circ a b = a + b + a*b        = (a+1) * (b+1) - 1     (the circle operation (called mult-add before 1.0))
    rung 3  fexp a b = (a+1)^(b+1) - 1                             (fused exponent)
    rung 4  ftet a b = (a+1)↑↑(b+1) - 1                            (fused tetration)
-/
namespace FusedRungs

def sadd (a b : Nat) : Nat := a + b + 1
def circ (a b : Nat) : Nat := a + b + a * b
def fexp (a b : Nat) : Nat := (a + 1) ^ (b + 1) - 1
def tet (x : Nat) : Nat → Nat
  | 0 => 1
  | n + 1 => x ^ tet x n
def ftet (a b : Nat) : Nat := tet (a + 1) (b + 1) - 1

/-! ### Lifting: every rung is the plain operation seen through x ↦ 1 + x -/
theorem sadd_lift (a b : Nat) : sadd a b + 1 = (a + 1) + (b + 1) := by unfold sadd; omega
theorem circ_lift (a b : Nat) : circ a b + 1 = (a + 1) * (b + 1) := by unfold circ; grind
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

/-! ### Rung 2 (the circle operation): a commutative monoid with identity 0 -/
theorem circ_comm (a b : Nat) : circ a b = circ b a := by unfold circ; grind
theorem circ_assoc (a b c : Nat) : circ (circ a b) c = circ a (circ b c) := by unfold circ; grind
theorem circ_zero (a : Nat) : circ a 0 = a := by unfold circ; simp

/-! ### Rung 3 (fused exponent) -/
theorem fexp_zero_right (a : Nat) : fexp a 0 = a := by unfold fexp; simp
/-- 0 is a left absorber, so there is no left identity of the form 0. -/
theorem fexp_zero_left (b : Nat) : fexp 0 b = 0 := by unfold fexp; simp
/-- The fused exponent carries the circle operation inside it. -/
theorem fexp_carries_circ (a b : Nat) : fexp a b = circ a ((a + 1) ^ b - 1) := by
  apply lift_inj
  rw [fexp_lift, circ_lift]
  have := pow_pos' a b
  rw [Nat.sub_add_cancel this, Nat.pow_succ, Nat.mul_comm]
/-- Right distributivity over the circle operation:  fexp (a ◯ b) c = fexp a c ◯ fexp b c. -/
theorem fexp_circ_right (a b c : Nat) : fexp (circ a b) c = circ (fexp a c) (fexp b c) := by
  apply lift_inj
  rw [fexp_lift, circ_lift, circ_lift, fexp_lift, fexp_lift, Nat.mul_pow]
/-- Tower law:  fexp (fexp a b) c = fexp a (b ◯ c). -/
theorem fexp_tower (a b c : Nat) : fexp (fexp a b) c = fexp a (circ b c) := by
  apply lift_inj
  rw [fexp_lift, fexp_lift, fexp_lift, circ_lift, Nat.pow_mul]
/-- Sum law:  fexp a (sadd b c) = fexp a b ◯ fexp a c, with sadd the rung-1 shifted add. -/
theorem fexp_sadd (a b c : Nat) : fexp a (sadd b c) = circ (fexp a b) (fexp a c) := by
  apply lift_inj
  rw [fexp_lift, circ_lift, fexp_lift, fexp_lift, sadd_lift, ← Nat.pow_add]
theorem fexp_not_assoc : fexp (fexp 1 2) 1 ≠ fexp 1 (fexp 2 1) := by decide
theorem fexp_not_comm : fexp 1 2 ≠ fexp 2 1 := by decide
theorem fexp_example : fexp (fexp 1 2) 1 = 63 ∧ fexp 1 (fexp 2 1) = 511 := by decide

/-! ### Each rung is the rung below, iterated -/
theorem circ_iter (a b : Nat) : circ a (b + 1) = sadd a (circ a b) := by
  apply lift_inj; rw [circ_lift, sadd_lift, circ_lift]; grind
theorem fexp_iter (a b : Nat) : fexp a (b + 1) = circ a (fexp a b) := by
  apply lift_inj; rw [fexp_lift, circ_lift, fexp_lift, Nat.pow_succ, Nat.mul_comm]
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

/-! ### Burgers: why the top corner of the Cole–Hopf cube must be empty
  A term e^(a x + b t) solves the heat equation φ_t = ν φ_xx iff b = ν a².
  The corner term of e1·e2 has slope a1+a2 and rate b1+b2; its heat defect is -2 ν a1 a2. -/
theorem burgers_corner_defect (ν a1 a2 : Int) :
    (ν * a1 ^ 2 + ν * a2 ^ 2) - ν * (a1 + a2) ^ 2 = -2 * ν * a1 * a2 := by grind
theorem burgers_corner_ok_iff (ν a1 a2 : Int) (hν : ν ≠ 0) :
    (ν * a1 ^ 2 + ν * a2 ^ 2) - ν * (a1 + a2) ^ 2 = 0 ↔ a1 = 0 ∨ a2 = 0 := by
  rw [burgers_corner_defect]
  constructor
  · intro h
    have h2 : ν * (a1 * a2) = 0 := by grind
    rcases Int.mul_eq_zero.mp h2 with h3 | h3
    · exact absurd h3 hν
    · exact Int.mul_eq_zero.mp h3
  · rintro (h | h) <;> subst h <;> simp

/-! ### Helical shell model: every corner of the triad-sign cube conserves energy and helicity
  c0 = (s1 l - s2 l²) G,  c1 = (s2 l² - s0) G,  c2 = c0 + c1,  G = s0 + s1 l + s2 l².
  Energy change ∝ c0 + c1 - c2; helicity change ∝ s0 c0 + s1 l c1 - s2 l² c2.  -/
def c0 (s0 s1 s2 l : Int) : Int := (s1 * l - s2 * l ^ 2) * (s0 + s1 * l + s2 * l ^ 2)
def c1 (s0 s1 s2 l : Int) : Int := (s2 * l ^ 2 - s0) * (s0 + s1 * l + s2 * l ^ 2)
def c2 (s0 s1 s2 l : Int) : Int := c0 s0 s1 s2 l + c1 s0 s1 s2 l
theorem shell_energy (s0 s1 s2 l : Int) : c0 s0 s1 s2 l + c1 s0 s1 s2 l - c2 s0 s1 s2 l = 0 := by
  unfold c2; omega
/-- Helicity is conserved for every sign corner and every shell ratio l, as a polynomial identity
    (it does not even need s² = 1). -/
theorem shell_helicity (s0 s1 s2 l : Int) :
    s0 * c0 s0 s1 s2 l + s1 * l * c1 s0 s1 s2 l - s2 * l ^ 2 * c2 s0 s1 s2 l = 0 := by
  unfold c2 c0 c1; grind

/-! ### The triad-sign cube has layers 1, 3, 3, 1 (by number of minus signs) -/
def signs : List (Bool × Bool × Bool) :=
  [false, true].flatMap fun a => [false, true].flatMap fun b => [false, true].map fun c => (a, b, c)
def minuses (t : Bool × Bool × Bool) : Nat :=
  (if t.1 then 1 else 0) + (if t.2.1 then 1 else 0) + (if t.2.2 then 1 else 0)
theorem sign_cube_layers :
    (List.range 4).map (fun j => (signs.filter (fun t => minuses t == j)).length) = [1, 3, 3, 1] := by
  decide
/-- The homochiral corners (+++) and (---) are exactly layers 0 and 3. -/
theorem homochiral_layers (t : Bool × Bool × Bool) (ht : t.1 = t.2.1 ∧ t.2.1 = t.2.2) :
    minuses t = 0 ∨ minuses t = 3 := by
  obtain ⟨a, b, c⟩ := t; simp at ht; obtain ⟨rfl, rfl⟩ := ht; cases a <;> simp [minuses]

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
/-- The rung-3 laws stop at rung 4: distributivity over the circle operation fails ... -/
theorem ftet_no_circ_dist : ftet (circ 1 1) 1 ≠ circ (ftet 1 1) (ftet 1 1) := by decide
/-- ... the tower law fails ... -/
theorem ftet_no_tower : ftet (ftet 1 1) 1 ≠ ftet 1 (fexp 1 1) := by decide
/-- ... and the sum law fails. -/
theorem ftet_no_sum : ftet 1 (sadd 1 1) ≠ fexp (ftet 1 1) (ftet 1 1) := by decide

end FusedRungs
