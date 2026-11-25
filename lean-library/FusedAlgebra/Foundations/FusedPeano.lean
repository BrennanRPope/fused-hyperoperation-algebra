/-
  Axioms.FusedPeano: whole numbers with mult-add as a primitive
  ==============================================================

  Peano's axioms with the fused ladder in place of + and ×. On a set N with
  a starting point 0, a next-step map S, rung 1 (a ⊞ b, which is a + b + 1)
  and rung 2 (a ⊕ b, which is a + b + ab):

    P1  S n ≠ 0                    0 is where counting starts
    P2  S m = S n → m = n          different numbers have different next steps
    R1  a ⊞ 0 = S a                joining a rod to the unit rod
        a ⊞ S b = S (a ⊞ b)
    R2  a ⊕ 0 = a                  a box one cell wide is just the row
        a ⊕ S b = (a ⊕ b) ⊞ a      a box grows one row at a time
    I   induction: whatever holds at 0 and passes from n to S n holds everywhere

  Rung 1 is needed in the axioms for the same reason ordinary Peano needs +
  before it can define ×: the recursion for rung 2 is written with rung 1.
  As a matter of what can be *defined*, though, 0, S and ⊕ are enough: + and ×
  are first-order definable from them (`add_definable`, `mul_definable`).

  Proved here from the axioms (core Lean, no `sorry`):
  * rung 1 and rung 2 laws by induction alone: ⊕ commutative, associative,
    identity 0, and ⊕ distributes over ⊞ exactly (`madd_sadd`),
  * ordinary + and × defined from ⊞ and ⊕, with their Peano recursions,
  * categoricity: every model is a copy of Nat (`toN_bijective`, `toN_sadd`,
    `toN_madd`), so these axioms pin the numbers down completely,
  * recursion on any model, and the whole fused ladder from one schema
    (`ladN`), whose rungs 1 and 2 are the primitives,
  * every rung n of the fused ladder is the ordinary hyperoperation seen
    through x ↦ 1 + x, for all n at once (`lad_lift`),
  * a derived commutative semiring on any model (`instCommSemiring`).

  LIBRARY NOTE: copied from axioms/lean/Axioms/FusedPeano.lean (Nov 25) as Foundations/FusedPeano.
-/

namespace MAFound

class FusedPeano (N : Type) where
  zero : N
  succ : N → N
  sadd : N → N → N
  madd : N → N → N
  succ_ne_zero : ∀ n, succ n ≠ zero
  succ_inj : ∀ m n, succ m = succ n → m = n
  sadd_zero : ∀ a, sadd a zero = succ a
  sadd_succ : ∀ a b, sadd a (succ b) = succ (sadd a b)
  madd_zero : ∀ a, madd a zero = a
  madd_succ : ∀ a b, madd a (succ b) = sadd (madd a b) a
  ind : ∀ (P : N → Prop), P zero → (∀ n, P n → P (succ n)) → ∀ n, P n

namespace FusedPeano

variable {N : Type} [FusedPeano N]

/-! ## Rung 1 from the axioms -/

theorem succ_sadd (a b : N) : sadd (succ a) b = succ (sadd a b) := by
  refine ind (fun b => sadd (succ a) b = succ (sadd a b)) ?_ ?_ b
  · rw [sadd_zero, sadd_zero]
  · intro b ih; rw [sadd_succ, sadd_succ, ih]

theorem zero_sadd (b : N) : sadd zero b = succ b := by
  refine ind (fun b => sadd zero b = succ b) ?_ ?_ b
  · rw [sadd_zero]
  · intro b ih; rw [sadd_succ, ih]

/-- Rung 1 is commutative: two rods joined either way round. -/
theorem sadd_comm (a b : N) : sadd a b = sadd b a := by
  refine ind (fun b => sadd a b = sadd b a) ?_ ?_ b
  · rw [sadd_zero, zero_sadd]
  · intro b ih; rw [sadd_succ, succ_sadd, ih]

/-- Rung 1 is associative. -/
theorem sadd_assoc (a b c : N) : sadd (sadd a b) c = sadd a (sadd b c) := by
  refine ind (fun c => sadd (sadd a b) c = sadd a (sadd b c)) ?_ ?_ c
  · rw [sadd_zero, sadd_zero, sadd_succ]
  · intro c ih; rw [sadd_succ, sadd_succ, sadd_succ, ih]

theorem sadd_right_cancel {a b c : N} (h : sadd a c = sadd b c) : a = b := by
  revert h
  refine ind (fun c => sadd a c = sadd b c → a = b) ?_ ?_ c
  · intro h; rw [sadd_zero, sadd_zero] at h; exact succ_inj _ _ h
  · intro c ih h; rw [sadd_succ, sadd_succ] at h; exact ih (succ_inj _ _ h)

theorem sadd_ne_zero (a b : N) : sadd a b ≠ zero := by
  refine ind (fun b => sadd a b ≠ zero) ?_ ?_ b
  · rw [sadd_zero]; exact succ_ne_zero a
  · intro b _; rw [sadd_succ]; exact succ_ne_zero _

/-! ## Rung 2 (mult-add) from the axioms -/

/-- 0 is a left identity too. -/
theorem zero_madd (b : N) : madd zero b = b := by
  refine ind (fun b => madd zero b = b) ?_ ?_ b
  · rw [madd_zero]
  · intro b ih; rw [madd_succ, ih, sadd_zero]

/-- Growing the box along its other side. -/
theorem succ_madd (a b : N) : madd (succ a) b = sadd (madd a b) b := by
  refine ind (fun b => madd (succ a) b = sadd (madd a b) b) ?_ ?_ b
  · rw [madd_zero, madd_zero, sadd_zero]
  · intro b ih
    rw [madd_succ, ih, madd_succ, sadd_succ, sadd_succ, sadd_assoc, sadd_assoc, sadd_comm b a]

/-- Mult-add is commutative: the box turned on its side. -/
theorem madd_comm (a b : N) : madd a b = madd b a := by
  refine ind (fun b => madd a b = madd b a) ?_ ?_ b
  · rw [madd_zero, zero_madd]
  · intro b ih; rw [madd_succ, succ_madd, ih]

/-- Rung 2 distributes over rung 1 exactly, with no correction term:
    a ⊕ (b ⊞ c) = (a ⊕ b) ⊞ (a ⊕ c). Two boxes of the same width, stacked. -/
theorem madd_sadd (a b c : N) : madd a (sadd b c) = sadd (madd a b) (madd a c) := by
  refine ind (fun c => madd a (sadd b c) = sadd (madd a b) (madd a c)) ?_ ?_ c
  · rw [sadd_zero, madd_succ, madd_zero]
  · intro c ih; rw [sadd_succ, madd_succ, ih, madd_succ, sadd_assoc]

/-- Mult-add is associative: a 3-D box built two ways. -/
theorem madd_assoc (a b c : N) : madd (madd a b) c = madd a (madd b c) := by
  refine ind (fun c => madd (madd a b) c = madd a (madd b c)) ?_ ?_ c
  · rw [madd_zero, madd_zero]
  · intro c ih; rw [madd_succ, ih, madd_succ, madd_sadd]

/-! ## Ordinary addition and multiplication, defined -/

theorem sadd_is_succ (x y : N) : ∃ z, succ z = sadd x y := by
  refine ind (fun y => ∃ z, succ z = sadd x y) ?_ ?_ y
  · exact ⟨x, (sadd_zero x).symm⟩
  · intro y ⟨z, hz⟩; exact ⟨sadd x y, (sadd_succ x y).symm⟩

/-- Addition: x + y is the number one step below x ⊞ y. -/
noncomputable def add (x y : N) : N := Classical.choose (sadd_is_succ x y)

theorem succ_add_eq (x y : N) : succ (add x y) = sadd x y :=
  Classical.choose_spec (sadd_is_succ x y)

theorem add_zero' (x : N) : add x zero = x :=
  succ_inj _ _ (by rw [succ_add_eq, sadd_zero])

theorem add_succ' (x y : N) : add x (succ y) = succ (add x y) :=
  succ_inj _ _ (by rw [succ_add_eq, sadd_succ, succ_add_eq])

theorem zero_add' (x : N) : add zero x = x :=
  succ_inj _ _ (by rw [succ_add_eq, zero_sadd])

theorem add_comm' (x y : N) : add x y = add y x :=
  succ_inj _ _ (by rw [succ_add_eq, succ_add_eq, sadd_comm])

theorem succ_add' (x y : N) : add (succ x) y = succ (add x y) := by
  rw [add_comm', add_succ', add_comm']

theorem add_assoc' (x y z : N) : add (add x y) z = add x (add y z) := by
  apply succ_inj; apply succ_inj
  rw [succ_add_eq, ← succ_sadd, succ_add_eq, sadd_assoc, ← succ_add_eq y z, sadd_succ, succ_add_eq]

theorem add_right_cancel' {a b c : N} (h : add a c = add b c) : a = b := by
  apply sadd_right_cancel (c := c)
  rw [← succ_add_eq, ← succ_add_eq, h]

theorem madd_is (x y : N) : ∃ z, madd x y = add (add x y) z := by
  refine ind (fun y => ∃ z, madd x y = add (add x y) z) ?_ ?_ y
  · exact ⟨zero, by rw [madd_zero, add_zero', add_zero']⟩
  · intro y ⟨z, hz⟩
    refine ⟨add z x, ?_⟩
    rw [madd_succ, ← succ_add_eq, hz, add_succ', succ_add', add_assoc' (add x y) z x]

/-- Multiplication: x · y is what is left of x ⊕ y after taking off x and y. -/
noncomputable def mul (x y : N) : N := Classical.choose (madd_is x y)

theorem madd_eq' (x y : N) : madd x y = add (add x y) (mul x y) :=
  Classical.choose_spec (madd_is x y)

theorem mul_zero' (x : N) : mul x zero = zero := by
  have h := madd_eq' x zero
  rw [madd_zero, add_zero'] at h
  apply add_right_cancel' (c := x)
  rw [zero_add', add_comm']; exact h.symm

theorem mul_succ' (x y : N) : mul x (succ y) = add (mul x y) x := by
  have h1 := madd_eq' x (succ y)
  rw [madd_succ, ← succ_add_eq, madd_eq', add_succ', succ_add'] at h1
  have h2 := succ_inj _ _ h1
  -- h2 : x + y + xy + x = x + y + x(y+1)
  apply add_right_cancel' (c := add x y)
  rw [add_comm' (mul x (succ y)), ← h2, add_assoc', add_comm' (add x y), add_assoc']

/-! ## Categoricity: every model is a copy of Nat -/

/-- Counting up in the model: k ↦ S(S(…S 0)). -/
def toN : Nat → N
  | 0 => zero
  | k + 1 => succ (toN k)

theorem toN_injective : ∀ {j k : Nat}, toN (N := N) j = toN k → j = k
  | 0, 0, _ => rfl
  | 0, _ + 1, h => absurd h.symm (succ_ne_zero _)
  | _ + 1, 0, h => absurd h (succ_ne_zero _)
  | _ + 1, _ + 1, h => by rw [toN_injective (succ_inj _ _ h)]

theorem toN_surjective (n : N) : ∃ k, toN k = n := by
  refine ind (fun n => ∃ k, toN k = n) ⟨0, rfl⟩ ?_ n
  intro n ⟨k, hk⟩; exact ⟨k + 1, by rw [toN, hk]⟩

theorem toN_bijective :
    (∀ j k : Nat, toN (N := N) j = toN k → j = k) ∧ (∀ n : N, ∃ k, toN k = n) :=
  ⟨fun _ _ h => toN_injective h, toN_surjective⟩

/-- Counting carries rung 1 to rung 1: (a + b + 1) ↦ a ⊞ b. -/
theorem toN_sadd (a b : Nat) : toN (N := N) (a + b + 1) = sadd (toN a) (toN b) := by
  induction b with
  | zero => show toN (a + 0 + 1) = sadd (toN a) zero; rw [sadd_zero]; rfl
  | succ b ih =>
    show toN (a + (b + 1) + 1) = sadd (toN a) (succ (toN b))
    rw [sadd_succ, ← ih]; rfl

/-- Counting carries mult-add to mult-add: (a + b + ab) ↦ a ⊕ b. -/
theorem toN_madd (a b : Nat) : toN (N := N) (a + b + a * b) = madd (toN a) (toN b) := by
  induction b with
  | zero => show toN (a + 0 + a * 0) = madd (toN a) zero; rw [madd_zero]; simp
  | succ b ih =>
    rw [show a + (b + 1) + a * (b + 1) = (a + b + a * b) + a + 1 by grind, toN_sadd, ih]
    show _ = madd (toN a) (succ (toN b)); rw [madd_succ]

theorem toN_add (a b : Nat) : toN (N := N) (a + b) = add (toN a) (toN b) := by
  apply succ_inj; rw [succ_add_eq, ← toN_sadd]; rfl

theorem toN_mul (a b : Nat) : toN (N := N) (a * b) = mul (toN a) (toN b) := by
  apply add_right_cancel' (c := add (toN a) (toN b))
  rw [add_comm' (mul _ _), ← madd_eq', ← toN_madd, ← toN_add, ← toN_add, Nat.add_comm (a * b)]

/-- The way back from the model to Nat. -/
noncomputable def fromN (n : N) : Nat := Classical.choose (toN_surjective n)

theorem toN_fromN (n : N) : toN (fromN n) = n := Classical.choose_spec (toN_surjective n)
theorem fromN_toN (k : Nat) : fromN (toN (N := N) k) = k := toN_injective (toN_fromN _)

/-! ## Recursion on any model, from categoricity -/

/-- Define a function on the model by "value at 0" and "how to take a step". -/
noncomputable def recN {α : Type} (z : α) (s : N → α → α) (n : N) : α :=
  Nat.rec (motive := fun _ => α) z (fun k acc => s (toN k) acc) (fromN n)

theorem recN_zero {α : Type} (z : α) (s : N → α → α) : recN z s (zero : N) = z := by
  unfold recN; rw [show (zero : N) = toN 0 from rfl, fromN_toN]; rfl

theorem recN_succ {α : Type} (z : α) (s : N → α → α) (n : N) :
    recN z s (succ n) = s n (recN z s n) := by
  unfold recN
  obtain ⟨k, rfl⟩ := toN_surjective n
  rw [show succ (toN k) = toN (N := N) (k + 1) from rfl, fromN_toN, fromN_toN]

/-! ## The fused ladder, from one schema

  Rung 0 is the step, and every rung above it is the rung below, iterated:
    F 0 a b = S b
    F 1 a 0 = S a,   F (n+2) a 0 = a      (from rung 2 up, 0 is the right identity)
    F (n+1) a (S b) = F n a (F (n+1) a b)
-/

/-- Where each fused rung starts: F 1 a 0 = a + 1, and a from rung 2 up. -/
def ladBase : Nat → Nat → Nat
  | 0, a => a + 1
  | _ + 1, a => a

/-- The fused ladder on Nat. -/
def lad : Nat → Nat → Nat → Nat
  | 0 => fun _ b => b + 1
  | n + 1 => fun a b => Nat.rec (motive := fun _ => Nat) (ladBase n a) (fun _ acc => lad n a acc) b

/-- Where each ordinary rung starts: a + 0 = a, a·0 = 0, a^0 = 1, … -/
def hypBase : Nat → Nat → Nat
  | 0, a => a
  | 1, _ => 0
  | _ + 2, _ => 1

/-- The ordinary hyperoperation ladder (successor, +, ×, power, tetration, …). -/
def hyp : Nat → Nat → Nat → Nat
  | 0 => fun _ b => b + 1
  | n + 1 => fun a b => Nat.rec (motive := fun _ => Nat) (hypBase n a) (fun _ acc => hyp n a acc) b

theorem hyp_one : ∀ m x, hyp (m + 2) x 1 = x
  | 0, x => rfl
  | m + 1, x => by
    show hyp (m + 2) x (hyp (m + 3) x 0) = x
    rw [show hyp (m + 3) x 0 = 1 from rfl]
    exact hyp_one m x

/-- Every rung of the fused ladder is the ordinary hyperoperation seen
    through x ↦ 1 + x:  F n a b + 1 = H n (a + 1) (b + 1), for every n. -/
theorem lad_lift : ∀ n a b, lad n a b + 1 = hyp n (a + 1) (b + 1)
  | 0, a, b => by simp [lad, hyp]
  | n + 1, a, 0 => by
    cases n with
    | zero => rfl
    | succ m => rw [hyp_one]; rfl
  | n + 1, a, b + 1 => by
    show lad n a (lad (n + 1) a b) + 1 = hyp n (a + 1) (hyp (n + 1) (a + 1) (b + 1))
    rw [lad_lift n a (lad (n + 1) a b), lad_lift (n + 1) a b]
termination_by n _ b => (n, b)

theorem hyp_add (x y : Nat) : hyp 1 x y = x + y := by
  induction y with
  | zero => rfl
  | succ y ih => show hyp 0 x (hyp 1 x y) = _; rw [ih]; simp [hyp]; omega

theorem hyp_mul (x y : Nat) : hyp 2 x y = x * y := by
  induction y with
  | zero => rfl
  | succ y ih => show hyp 1 x (hyp 2 x y) = _; rw [ih, hyp_add, Nat.mul_succ, Nat.add_comm]

theorem hyp_pow (x y : Nat) : hyp 3 x y = x ^ y := by
  induction y with
  | zero => rfl
  | succ y ih => show hyp 2 x (hyp 3 x y) = _; rw [ih, hyp_mul, Nat.pow_succ, Nat.mul_comm]

theorem lad_one (a b : Nat) : lad 1 a b = a + b + 1 := by
  have := lad_lift 1 a b; rw [hyp_add] at this; omega

theorem lad_two (a b : Nat) : lad 2 a b = a + b + a * b := by
  have := lad_lift 2 a b; rw [hyp_mul] at this; grind

theorem lad_three (a b : Nat) : lad 3 a b + 1 = (a + 1) ^ (b + 1) := by
  rw [lad_lift, hyp_pow]

/-- The fused ladder on any model, defined by the same schema with `recN`. -/
noncomputable def ladN : Nat → N → N → N
  | 0, _, b => succ b
  | 1, a, b => recN (succ a) (fun _ acc => succ acc) b
  | n + 2, a, b => recN a (fun _ acc => ladN (n + 1) a acc) b

theorem ladN_zero_right (n : Nat) (a : N) :
    ladN (n + 1) a zero = if n = 0 then succ a else a := by
  cases n with
  | zero => simp [ladN, recN_zero]
  | succ m => simp [ladN, recN_zero]

theorem ladN_succ (n : Nat) (a b : N) : ladN (n + 1) a (succ b) = ladN n a (ladN (n + 1) a b) := by
  cases n with
  | zero => simp only [ladN, recN_succ]
  | succ m => simp only [ladN]; rw [recN_succ]

/-- The schema's rung 1 is the primitive ⊞. -/
theorem ladN_one (a b : N) : ladN 1 a b = sadd a b := by
  refine ind (fun b => ladN 1 a b = sadd a b) ?_ ?_ b
  · rw [ladN_zero_right, sadd_zero]; rfl
  · intro b ih; rw [ladN_succ, ih, sadd_succ]; rfl

/-- The schema's rung 2 is the primitive ⊕. -/
theorem ladN_two (a b : N) : ladN 2 a b = madd a b := by
  refine ind (fun b => ladN 2 a b = madd a b) ?_ ?_ b
  · rw [ladN_zero_right, madd_zero]; rfl
  · intro b ih; rw [ladN_succ, ih, ladN_one, madd_succ, sadd_comm]

/-- Counting carries the whole ladder over: every model has every rung. -/
theorem toN_lad (n a b : Nat) : toN (N := N) (lad n a b) = ladN n (toN a) (toN b) := by
  induction n generalizing b with
  | zero => simp [lad, ladN, toN]
  | succ n ihn =>
    induction b with
    | zero =>
      rw [show toN (N := N) 0 = zero from rfl, ladN_zero_right]
      cases n with
      | zero => simp [lad, ladBase, toN]
      | succ m => rfl
    | succ b ihb =>
      show toN (lad n a (lad (n + 1) a b)) = _
      rw [ihn, ihb, show toN (N := N) (b + 1) = succ (toN b) from rfl, ladN_succ]

/-! ## A derived commutative semiring on any model -/

theorem add_toN_fromN (x y : N) : add x y = toN (fromN x + fromN y) := by
  rw [toN_add, toN_fromN, toN_fromN]

noncomputable def npow' (a : N) : Nat → N
  | 0 => toN 1
  | n + 1 => mul (npow' a n) a

noncomputable instance : Add N := ⟨add⟩
noncomputable instance : Mul N := ⟨mul⟩
instance : NatCast N := ⟨toN⟩
instance (n : Nat) : OfNat N n := ⟨toN n⟩
noncomputable instance : SMul Nat N := ⟨fun n a => mul (toN n) a⟩
noncomputable instance : HPow N Nat N := ⟨npow'⟩

/-- Pull an equation on the model back to Nat. -/
theorem transfer2 {P : N → N → Prop} (h : ∀ a b : Nat, P (toN a) (toN b)) (x y : N) : P x y := by
  rw [← toN_fromN x, ← toN_fromN y]; exact h _ _

theorem transfer3 {P : N → N → N → Prop} (h : ∀ a b c : Nat, P (toN a) (toN b) (toN c))
    (x y z : N) : P x y z := by
  rw [← toN_fromN x, ← toN_fromN y, ← toN_fromN z]; exact h _ _ _

theorem mul_comm' (x y : N) : mul x y = mul y x :=
  transfer2 (P := fun x y => mul x y = mul y x)
    (fun a b => by rw [← toN_mul, ← toN_mul, Nat.mul_comm]) x y

theorem mul_assoc' (x y z : N) : mul (mul x y) z = mul x (mul y z) :=
  transfer3 (P := fun x y z => mul (mul x y) z = mul x (mul y z))
    (fun a b c => by rw [← toN_mul, ← toN_mul, ← toN_mul, ← toN_mul, Nat.mul_assoc]) x y z

theorem mul_add' (x y z : N) : mul x (add y z) = add (mul x y) (mul x z) :=
  transfer3 (P := fun x y z => mul x (add y z) = add (mul x y) (mul x z))
    (fun a b c => by rw [← toN_add, ← toN_mul, ← toN_mul, ← toN_mul, ← toN_add, Nat.mul_add]) x y z

theorem mul_one' (x : N) : mul x (toN 1) = x := by
  rw [← toN_fromN x, ← toN_mul, Nat.mul_one]

/-- Any model of the fused Peano axioms is a commutative semiring, with
    + and × derived from ⊞ and ⊕. -/
noncomputable instance instCommSemiring : Lean.Grind.CommSemiring N where
  add_zero := add_zero'
  add_comm := add_comm'
  add_assoc := add_assoc'
  mul_assoc := mul_assoc'
  mul_one := mul_one'
  one_mul := fun a => by show mul _ a = a; rw [mul_comm']; exact mul_one' a
  left_distrib := mul_add'
  right_distrib := fun a b c => by
    show mul (add a b) c = add (mul a c) (mul b c)
    rw [mul_comm', mul_add', mul_comm' c, mul_comm' c]
  zero_mul := fun a => by show mul zero a = zero; rw [mul_comm', mul_zero']
  mul_zero := mul_zero'
  pow_zero := fun _ => rfl
  pow_succ := fun _ _ => rfl
  ofNat_succ := fun n => by show toN (n + 1) = add (toN n) (toN 1); rw [← toN_add]
  mul_comm := mul_comm'

/-- In the derived semiring, a + b + a·b is the primitive ⊕. -/
theorem multAdd_eq (a b : N) : a + b + a * b = madd a b := (madd_eq' a b).symm

/-! ## Strength: + and × can be defined from successor and mult-add

  Mult-add alone on the whole numbers is multiplication in disguise
  (a ⊕ b + 1 = (a + 1)(b + 1)), and the first-order theory of multiplication
  alone is decidable (Skolem 1930, Mostowski 1952): too weak to carry
  arithmetic. With the successor added, + and × become definable, so the
  structure (Nat, 0, S, ⊕) carries all of first-order arithmetic. -/

/-- Mult-add on Nat. -/
def ma (a b : Nat) : Nat := a + b + a * b

/-- × from S and ⊕: m·n = k exactly when one side is 0 and k = 0, or
    m = S m', n = S n', k = S k' and m' ⊕ n' = k'. -/
theorem mul_definable (m n k : Nat) :
    m * n = k ↔ (m = 0 ∧ k = 0) ∨ (n = 0 ∧ k = 0) ∨
      ∃ m' n' k', m = m' + 1 ∧ n = n' + 1 ∧ k = k' + 1 ∧ ma m' n' = k' := by
  unfold ma
  constructor
  · intro h
    cases m with
    | zero => left; simp at h; exact ⟨rfl, h.symm⟩
    | succ m' =>
      cases n with
      | zero => right; left; simp at h; exact ⟨rfl, h.symm⟩
      | succ n' =>
        right; right
        refine ⟨m', n', m' + n' + m' * n', rfl, rfl, ?_, rfl⟩
        rw [← h]; grind
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨m', n', k', rfl, rfl, rfl, h⟩)
    · simp
    · simp
    · rw [← h]; grind

/-- + from ×: Julia Robinson's trick (1949). For z ≠ 0,
    x + y = z exactly when (xz + 1)(yz + 1) = z²(xy + 1) + 1. -/
theorem add_definable (x y z : Nat) :
    x + y = z ↔ (z = 0 ∧ x = 0 ∧ y = 0) ∨
      (z ≠ 0 ∧ (x * z + 1) * (y * z + 1) = z * z * (x * y + 1) + 1) := by
  constructor
  · intro h
    by_cases hz : z = 0
    · left; omega
    · right; refine ⟨hz, ?_⟩; subst h; grind
  · rintro (⟨rfl, rfl, rfl⟩ | ⟨hz, h⟩)
    · rfl
    · have e : z * (x + y) = z * z := by
        have : (x * z + 1) * (y * z + 1) = x * y * (z * z) + z * (x + y) + 1 := by grind
        have : z * z * (x * y + 1) + 1 = x * y * (z * z) + z * z + 1 := by grind
        omega
      exact Nat.eq_of_mul_eq_mul_left (Nat.pos_of_ne_zero hz) e

end FusedPeano
end MAFound
