/-
  Axioms.MultAddRing: arithmetic with mult-add as the primitive
  ==============================================================

  Mult-add a ⊕ b = a + b + ab is taken here as a basic operation, not as a
  shorthand. Multiplication is then defined from it:  a · b := (a ⊕ b) − a − b.

  The axioms of a mult-add ring (MA-ring), on a set R with 0, +, −, ⊕ and a
  constant w (the wall):

    A1  (a + b) + c = a + (b + c)          laying rods end to end
    A2  a + b = b + a
    A3  a + 0 = a
    A4  −a + a = 0
    M1  a ⊕ b = b ⊕ a                      a box turned on its side
    M2  (a ⊕ b) ⊕ c = a ⊕ (b ⊕ c)          a 3-D box built two ways
    M3  a ⊕ w = w                          the wall: a side of length 0 flattens the box
    M4  a ⊕ (b + c) = (a ⊕ b) + (a ⊕ c) − a
                                           cutting a box in two along one side:
                                           the shared strip is counted once

  "0 is the identity of ⊕" is NOT an axiom: it follows from M4 (`madd_zero`).
  `Axioms.Independence` shows M1, M2, M3 and M4 are each independent.

  Main results (all proved here, core Lean, no `sorry`):
  * every MA-ring is a commutative ring (`instCommRing`), with 1 = −w,
  * in that ring the library's mult-add a + b + ab IS the primitive ⊕ (`multAdd_eq`),
  * every commutative ring is an MA-ring (`ofCommRing`), and going there and
    back changes nothing (`roundtrip_madd`, `roundtrip_mul`),
  * an MA-field (one more axiom: everything off the wall has a way back) is a
    field (`instField`).
-/

namespace MAFound

/-- A mult-add ring: an abelian group (+, −, 0) with a primitive mult-add ⊕
    and a wall w. -/
class MARing (R : Type) where
  zero : R
  add : R → R → R
  neg : R → R
  madd : R → R → R
  wall : R
  add_assoc : ∀ a b c, add (add a b) c = add a (add b c)
  add_comm : ∀ a b, add a b = add b a
  add_zero : ∀ a, add a zero = a
  neg_add : ∀ a, add (neg a) a = zero
  madd_comm : ∀ a b, madd a b = madd b a
  madd_assoc : ∀ a b c, madd (madd a b) c = madd a (madd b c)
  madd_wall : ∀ a, madd a wall = wall
  madd_add : ∀ a b c, madd a (add b c) = add (add (madd a b) (madd a c)) (neg a)

namespace MARing

variable {R : Type} [MARing R]

/-! ## The abelian group: a few rearrangement lemmas (A1–A4 only) -/

local infixl:65 " +' " => MARing.add
local prefix:max "−" => MARing.neg
local infixl:70 " ⊞ " => MARing.madd

instance : Std.Associative (MARing.add (R := R)) := ⟨add_assoc⟩
instance : Std.Commutative (MARing.add (R := R)) := ⟨add_comm⟩

theorem add_neg (a : R) : a +' −a = zero := by rw [add_comm, neg_add]
theorem zero_add (a : R) : zero +' a = a := by rw [add_comm, add_zero]

theorem add_right_cancel {a b c : R} (h : a +' c = b +' c) : a = b := by
  have := congrArg (· +' −c) h
  simp only [add_assoc, add_neg, add_zero] at this
  exact this

theorem add_left_cancel {a b c : R} (h : c +' a = c +' b) : a = b := by
  rw [add_comm c a, add_comm c b] at h; exact add_right_cancel h

theorem add_left_comm (a b c : R) : a +' (b +' c) = b +' (a +' c) := by
  rw [← add_assoc, add_comm a b, add_assoc]

theorem neg_unique {a b : R} (h : a +' b = zero) : b = −a := by
  apply add_left_cancel (c := a); rw [h, add_neg]

theorem neg_neg (a : R) : −(−a) = a := (neg_unique (neg_add a)).symm

theorem neg_zero : −(zero : R) = zero := (neg_unique (add_zero zero)).symm

theorem neg_add_distrib (a b : R) : −(a +' b) = −a +' −b := by
  symm; apply neg_unique
  rw [add_assoc, add_left_comm b, add_neg, add_zero, add_neg]

/-! ## Mult-add laws that follow -/

/-- 0 is the identity of mult-add. Not an axiom: cut a box along a side of
    length 0 and both pieces are the whole box, so the box adds nothing. -/
theorem madd_zero (a : R) : a ⊞ zero = a := by
  have h := madd_add a zero zero
  rw [add_zero] at h
  -- h : X = X + X − a   with X = a ⊞ 0
  have h2 : (a ⊞ zero) +' zero = (a ⊞ zero) +' ((a ⊞ zero) +' −a) := by
    rw [add_zero, ← add_assoc]; exact h
  have h3 := add_left_cancel h2
  have h4 := congrArg (· +' a) h3
  simp only [zero_add, add_assoc, neg_add, add_zero] at h4
  exact h4.symm

theorem zero_madd (a : R) : zero ⊞ a = a := by rw [madd_comm, madd_zero]

/-- The wall is unique: any element that absorbs everything is w. -/
theorem wall_unique (v : R) (hv : ∀ a, a ⊞ v = v) : v = wall := by
  have h1 := hv wall
  rw [madd_comm, madd_wall] at h1
  exact h1.symm

/-! ## Multiplication, defined from mult-add -/

/-- The product: a · b := (a ⊞ b) − a − b, the box minus its two edges. -/
def mul (a b : R) : R := (a ⊞ b) +' −a +' −b

/-- The unit: 1 := −w. -/
def one : R := −wall

local infixl:70 " ⬝ " => MARing.mul

/-- Mult-add is the two edges plus the corner. -/
theorem madd_eq (a b : R) : a ⊞ b = a +' b +' (a ⬝ b) := by
  unfold mul
  calc a ⊞ b = a ⊞ b +' (a +' −a) +' (b +' −b) := by rw [add_neg, add_neg, add_zero, add_zero]
    _ = a +' b +' (a ⊞ b +' −a +' −b) := by ac_rfl

theorem mul_comm (a b : R) : a ⬝ b = b ⬝ a := by
  unfold mul; rw [madd_comm]; ac_rfl

theorem mul_zero (a : R) : a ⬝ zero = zero := by
  unfold mul; rw [madd_zero, add_neg, zero_add, neg_zero]

theorem zero_mul (a : R) : zero ⬝ a = zero := by rw [mul_comm, mul_zero]

/-- Distributivity, from M4. -/
theorem mul_add (a b c : R) : a ⬝ (b +' c) = a ⬝ b +' a ⬝ c := by
  unfold mul
  rw [madd_add, neg_add_distrib]
  -- both sides are (a⊕b) + (a⊕c) − a − a − b − c, rearranged
  ac_rfl

theorem add_mul (a b c : R) : (a +' b) ⬝ c = a ⬝ c +' b ⬝ c := by
  rw [mul_comm, mul_add, mul_comm c, mul_comm c]

theorem mul_neg (a b : R) : a ⬝ (−b) = −(a ⬝ b) := by
  apply neg_unique; rw [← mul_add, add_neg, mul_zero]

/-- Associativity of the product, from M2 and M4. -/
theorem mul_assoc (a b c : R) : (a ⬝ b) ⬝ c = a ⬝ (b ⬝ c) := by
  have h := madd_assoc a b c
  rw [madd_eq a b, madd_eq b c, madd_eq, madd_eq a, add_mul, add_mul, mul_add, mul_add] at h
  -- h : a+b+ab+c+(ac+bc+(ab)c) = a+(b+c+bc)+(ab+ac+a(bc))
  have h' : (a +' b +' a ⬝ b +' c +' a ⬝ c +' b ⬝ c) +' (a ⬝ b) ⬝ c
          = (a +' b +' a ⬝ b +' c +' a ⬝ c +' b ⬝ c) +' a ⬝ (b ⬝ c) := by
    have e1 : a +' b +' a ⬝ b +' c +' (a ⬝ c +' b ⬝ c +' (a ⬝ b) ⬝ c)
        = (a +' b +' a ⬝ b +' c +' a ⬝ c +' b ⬝ c) +' (a ⬝ b) ⬝ c := by
      ac_rfl
    have e2 : a +' (b +' c +' b ⬝ c) +' (a ⬝ b +' a ⬝ c +' a ⬝ (b ⬝ c))
        = (a +' b +' a ⬝ b +' c +' a ⬝ c +' b ⬝ c) +' a ⬝ (b ⬝ c) := by
      ac_rfl
    rw [← e1, ← e2]; exact h
  exact add_left_cancel h'

/-- The wall times anything: a · w = −a. -/
theorem mul_wall (a : R) : a ⬝ wall = −a := by
  unfold mul; rw [madd_wall]
  calc wall +' −a +' −wall = −a +' (wall +' −wall) := by ac_rfl
    _ = −a := by rw [add_neg, add_zero]

theorem mul_one (a : R) : a ⬝ one = a := by
  unfold one; rw [mul_neg, mul_wall, neg_neg]

/-! ## The commutative ring it carries -/

/-- The counting numbers inside R: 1 := −w, 2 := 1 + 1, and so on. -/
def natCast : Nat → R
  | 0 => zero
  | 1 => one
  | n + 2 => natCast (n + 1) +' one

theorem natCast_succ (n : Nat) : natCast (R := R) (n + 1) = natCast n +' one := by
  cases n with
  | zero => show one = zero +' one; rw [zero_add]
  | succ n => rfl

def intCast : Int → R
  | .ofNat n => natCast n
  | .negSucc n => −(natCast (n + 1))

def npow (a : R) : Nat → R
  | 0 => one
  | n + 1 => npow a n ⬝ a

theorem intCast_neg' (i : Int) : intCast (R := R) (-i) = −(intCast i) := by
  match i with
  | .ofNat 0 => show natCast 0 = −(natCast 0); simp only [natCast]; rw [neg_zero]
  | .ofNat (n + 1) => rfl
  | .negSucc n => show natCast (n + 1) = −(−(natCast (n + 1))); rw [neg_neg]

instance : Add R := ⟨add⟩
instance : Mul R := ⟨mul⟩
instance : Neg R := ⟨neg⟩
instance : Sub R := ⟨fun a b => add a (neg b)⟩
instance : NatCast R := ⟨natCast⟩
instance : IntCast R := ⟨intCast⟩
instance (n : Nat) : OfNat R n := ⟨natCast n⟩
instance : SMul Nat R := ⟨fun n a => mul (natCast n) a⟩
instance : SMul Int R := ⟨fun i a => mul (intCast i) a⟩
instance : HPow R Nat R := ⟨npow⟩

/-- Every MA-ring is a commutative ring (with ⬝ as product and −w as 1). -/
instance instCommRing : Lean.Grind.CommRing R where
  add_zero := add_zero
  add_comm := add_comm
  add_assoc := add_assoc
  mul_assoc := mul_assoc
  mul_one := mul_one
  one_mul := fun a => by show mul one a = a; rw [mul_comm]; exact mul_one a
  left_distrib := mul_add
  right_distrib := add_mul
  zero_mul := zero_mul
  mul_zero := mul_zero
  pow_zero := fun _ => rfl
  pow_succ := fun _ _ => rfl
  ofNat_succ := natCast_succ
  neg_add_cancel := neg_add
  sub_eq_add_neg := fun _ _ => rfl
  neg_zsmul := fun i a => by
    show mul (intCast (-i)) a = neg (mul (intCast i) a)
    rw [intCast_neg', mul_comm, mul_neg, mul_comm]
  intCast_neg := intCast_neg'
  mul_comm := mul_comm

/-- In the derived ring, the library's mult-add a + b + a·b is exactly the
    primitive ⊕. So every theorem about mult-add is a theorem about ⊕. -/
theorem multAdd_eq (a b : R) : a + b + a * b = a ⊞ b := (madd_eq a b).symm

theorem one_eq : (1 : R) = −wall := rfl
theorem wall_eq : (wall : R) = -1 := by show wall = −(−wall); rw [neg_neg]
theorem zero_eq : (zero : R) = 0 := rfl

end MARing

/-! ## Every commutative ring is an MA-ring, and back and forth is the identity -/

open Lean.Grind in
/-- A commutative ring gives an MA-ring: ⊞ := a + b + ab and w := −1. -/
@[instance_reducible] def ofCommRing (R : Type) [CommRing R] : MARing R where
  zero := 0
  add := (· + ·)
  neg := (- ·)
  madd a b := a + b + a * b
  wall := -1
  add_assoc := Semiring.add_assoc
  add_comm := Semiring.add_comm
  add_zero := Semiring.add_zero
  neg_add := Ring.neg_add_cancel
  madd_comm _ _ := by grind
  madd_assoc _ _ _ := by grind
  madd_wall _ := by grind
  madd_add _ _ _ := by grind

open Lean.Grind in
/-- Ring → MA-ring → product: the product comes back unchanged. -/
theorem roundtrip_mul (R : Type) [CommRing R] (a b : R) :
    @MARing.mul R (ofCommRing R) a b = a * b := by
  show a + b + a * b + -a + -b = a * b; grind

/-- MA-ring → ring → mult-add: ⊞ comes back unchanged. -/
theorem roundtrip_madd (R : Type) [MARing R] (a b : R) :
    (ofCommRing R).madd a b = MARing.madd a b := MARing.multAdd_eq a b

/-- MA-ring → ring → MA-ring: the wall comes back unchanged. -/
theorem roundtrip_wall (R : Type) [MARing R] :
    (ofCommRing R).wall = (MARing.wall : R) := MARing.wall_eq.symm

/-! ## Fields: one more axiom, the way back

  MA-field = MA-ring with
    F1  0 ≠ w                          (not everything is flat)
    F2  for every a ≠ w there is b with a ⊞ b = 0
                                       (every box off the wall can be undone)
-/

class MAField (R : Type) extends MARing R where
  zero_ne_wall : zero ≠ wall
  way_back : ∀ a, a ≠ wall → ∃ b, madd a b = zero

namespace MAField

open MARing
variable {R : Type} [MAField R]

theorem add_wall_ne {a : R} (h : a ≠ 0) : (a + wall : R) ≠ wall := by
  intro e
  apply h
  have : (a + wall : R) = 0 + wall := by rw [e]; exact (zero_add wall).symm
  exact add_right_cancel this

theorem zero_ne_one' : (0 : R) ≠ 1 := by
  intro e
  apply zero_ne_wall (R := R)
  have w := wall_eq (R := R)
  rw [zero_eq, w, ← e]; grind

open Classical in
/-- The reciprocal, built from the way back: if (a − 1) ⊕ b = 0 then a·(1 + b) = 1. -/
noncomputable def inv (a : R) : R :=
  if h : a = 0 then 0 else 1 + Classical.choose (way_back (a + wall) (add_wall_ne h))

theorem mul_inv {a : R} (h : a ≠ 0) : a * inv a = 1 := by
  have hb := Classical.choose_spec (way_back (a + wall) (add_wall_ne h))
  have hinv : inv a = 1 + Classical.choose (way_back (a + wall) (add_wall_ne h)) := by
    unfold inv; simp only [h, dite_false]
  rw [hinv]
  generalize Classical.choose (way_back (a + wall) (add_wall_ne h)) = b at hb
  have e := multAdd_eq (a + wall) b
  rw [hb, wall_eq, zero_eq] at e
  grind

theorem inv_zero' : inv (0 : R) = 0 := by unfold inv; simp

theorem inv_ne_zero {a : R} (h : a ≠ 0) : inv a ≠ 0 := by
  intro e
  have h1 := mul_inv h
  rw [e] at h1
  apply zero_ne_one' (R := R)
  rw [← h1]; grind

theorem inv_inv' (a : R) : inv (inv a) = a := by
  by_cases h : a = 0
  · rw [h, inv_zero', inv_zero']
  · have h1 := mul_inv h
    have h2 := mul_inv (inv_ne_zero h)
    calc inv (inv a) = (a * inv a) * inv (inv a) := by rw [h1]; grind
      _ = a * (inv a * inv (inv a)) := by grind
      _ = a := by rw [h2]; grind

theorem inv_one' : inv (1 : R) = 1 := by
  have := mul_inv (R := R) (a := 1) (fun e => zero_ne_one' e.symm)
  grind

noncomputable def zpow (a : R) : Int → R
  | .ofNat n => a ^ n
  | .negSucc n => inv (a ^ (n + 1))

noncomputable instance instField : Lean.Grind.Field R :=
  letI : Inv R := ⟨inv⟩
  letI : Div R := ⟨fun a b => a * inv b⟩
  letI : HPow R Int R := ⟨zpow⟩
  { div_eq_mul_inv := fun _ _ => rfl
    zero_ne_one := zero_ne_one'
    inv_zero := inv_zero'
    mul_inv_cancel := fun h => mul_inv h
    zpow_zero := fun a => show a ^ 0 = 1 from Lean.Grind.Semiring.pow_zero a
    zpow_succ := fun a n => show a ^ (n + 1) = a ^ n * a from Lean.Grind.Semiring.pow_succ a n
    zpow_neg := fun a n => by
      match n with
      | .ofNat 0 =>
        show npow a 0 = inv (npow a 0)
        show (1 : R) = inv 1
        rw [inv_one']
      | .ofNat (n + 1) => rfl
      | .negSucc n => show a ^ (n + 1) = inv (inv (a ^ (n + 1))); rw [inv_inv'] }

/-- The way back, in closed form: the mult-add inverse of u is −u/(1 + u). -/
theorem way_back_formula (u : R) (h : 1 + u ≠ 0) :
    MARing.madd u (-u * inv (1 + u)) = (0 : R) := by
  rw [← multAdd_eq]
  have := mul_inv h
  grind

end MAField

end MAFound
