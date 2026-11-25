module

public import FusedAlgebra.Basic

/-!
# Mult-add rings: arithmetic with mult-add as the primitive

Mult-add is taken here as a basic operation, not as shorthand for `a + b + a * b`. A *mult-add
ring* (`MAFound.MARing`) is an abelian group `(R, +, 0)` with a binary operation `madd` and a
constant `wall` such that

* M1: `madd a b = madd b a` (a box turned on its side);
* M2: `madd (madd a b) c = madd a (madd b c)` (a 3-D box built two ways);
* M3: `madd a wall = wall` (the wall: a side of length `0` flattens the box);
* M4: `madd a (b + c) = madd a b + madd a c - a` (cutting a box in two along one side: the shared
  strip is counted once).

Multiplication is then *defined*: `a * b = madd a b - a - b`, and `1 = -wall`. That `0` is the
identity of `madd` is not an axiom; it follows from M4 (`MARing.madd_zero`).

## Main results

* `MARing.instCommRing`: every mult-add ring is a commutative ring (mathlib's `CommRing`).
* `MARing.multAdd_eq`: in that ring, `a ✳ b = a + b + a * b` is the primitive `madd`.
* `ofCommRing`: every commutative ring is a mult-add ring, with `madd = ✳` and `wall = -1`.
* `instCommRing_ofCommRing`, `ofCommRing_instCommRing`, `equivCommRing`: the two constructions
  are inverse to each other, so mult-add ring structures and commutative ring structures on a
  type are the same thing.
* `MAField`, `MAField.instField`, `ofField`: one more axiom (everything off the wall has a way
  back) gives exactly the fields.

This is the mathlib form of `FusedAlgebra.Foundations.MultAddRing`. The underlying mathematics is
known: `madd` is the circle operation of ring theory (see `FusedAlgebra.Basic` and mathlib's
`PreQuasiregular`), and a ring can be recovered from its circle operation and its additive group.
What the file records is that the axioms M1 to M4 are enough.
-/

@[expose] public section

namespace MAFound

open MultAdd

/-- A mult-add ring: an abelian group with a primitive mult-add `madd` and a wall `wall`. The
abelian group axioms (A1 to A4 in the original) are mathlib's `AddCommGroup`. -/
class MARing (R : Type*) extends AddCommGroup R where
  /-- The primitive mult-add operation. -/
  madd : R → R → R
  /-- The wall: the element that absorbs everything under `madd`. -/
  wall : R
  /-- M1: mult-add is commutative. -/
  madd_comm : ∀ a b, madd a b = madd b a
  /-- M2: mult-add is associative. -/
  madd_assoc : ∀ a b c, madd (madd a b) c = madd a (madd b c)
  /-- M3: the wall absorbs. -/
  madd_wall : ∀ a, madd a wall = wall
  /-- M4: cutting a box in two along one side counts the shared strip once. -/
  madd_add : ∀ a b c, madd a (b + c) = madd a b + madd a c - a

namespace MARing

variable {R : Type*} [MARing R]

/-! ### Mult-add laws that follow -/

/-- `0` is the identity of mult-add. Not an axiom: cut a box along a side of length `0` and both
pieces are the whole box, so the box adds nothing. -/
theorem madd_zero (a : R) : madd a 0 = a := by
  have h := madd_add a 0 0
  rw [add_zero, add_sub_assoc] at h
  exact sub_eq_zero.mp (left_eq_add.mp h)

theorem zero_madd (a : R) : madd 0 a = a := by rw [madd_comm, madd_zero]

/-- The wall is unique: any element that absorbs everything is the wall. -/
theorem wall_unique (v : R) (hv : ∀ a, madd a v = v) : v = wall := by
  rw [← hv wall, madd_comm, madd_wall]

/-! ### Multiplication, defined from mult-add -/

/-- The product: `a * b = madd a b - a - b`, the box minus its two edges. -/
instance instMul : Mul R := ⟨fun a b ↦ madd a b - a - b⟩

/-- The unit: `1 = -wall`. -/
instance instOne : One R := ⟨-wall⟩

theorem mul_def (a b : R) : a * b = madd a b - a - b := rfl

theorem one_eq_neg_wall : (1 : R) = -wall := rfl

/-- Mult-add is the two edges plus the corner. -/
theorem madd_eq (a b : R) : madd a b = a + b + a * b := by
  rw [mul_def]; abel

protected theorem mul_comm (a b : R) : a * b = b * a := by
  rw [mul_def, mul_def, madd_comm]; abel

protected theorem mul_zero (a : R) : a * 0 = 0 := by
  rw [mul_def, madd_zero]; abel

protected theorem zero_mul (a : R) : 0 * a = 0 := by
  rw [MARing.mul_comm, MARing.mul_zero]

/-- Distributivity, from M4. -/
protected theorem mul_add (a b c : R) : a * (b + c) = a * b + a * c := by
  simp only [mul_def, madd_add]; abel

protected theorem add_mul (a b c : R) : (a + b) * c = a * c + b * c := by
  rw [MARing.mul_comm, MARing.mul_add, MARing.mul_comm c, MARing.mul_comm c]

protected theorem mul_neg (a b : R) : a * -b = -(a * b) := by
  apply eq_neg_of_add_eq_zero_left
  rw [← MARing.mul_add, neg_add_cancel, MARing.mul_zero]

/-- Associativity of the product, from M2 and M4. -/
protected theorem mul_assoc (a b c : R) : a * b * c = a * (b * c) := by
  have h := madd_assoc a b c
  simp only [madd_eq, MARing.mul_add, MARing.add_mul] at h
  rw [← sub_eq_zero] at h ⊢
  rw [← h]; abel

/-- The wall times anything: `a * wall = -a`. -/
theorem mul_wall (a : R) : a * wall = -a := by
  rw [mul_def, madd_wall]; abel

protected theorem mul_one (a : R) : a * 1 = a := by
  rw [one_eq_neg_wall, MARing.mul_neg, mul_wall, neg_neg]

/-! ### The commutative ring it carries -/

/-- Every mult-add ring is a commutative ring, with product `madd a b - a - b` and unit
`-wall`. -/
instance instCommRing : CommRing R where
  __ := ‹MARing R›.toAddCommGroup
  mul_assoc := MARing.mul_assoc
  one_mul a := by rw [MARing.mul_comm, MARing.mul_one]
  mul_one := MARing.mul_one
  left_distrib := MARing.mul_add
  right_distrib := MARing.add_mul
  zero_mul := MARing.zero_mul
  mul_zero := MARing.mul_zero
  mul_comm := MARing.mul_comm

/-- In the derived ring, the library's mult-add `a ✳ b = a + b + a * b` is exactly the primitive
`madd`. So every theorem about `✳` is a theorem about `madd`. -/
theorem multAdd_eq (a b : R) : a ✳ b = madd a b := (madd_eq a b).symm

theorem wall_eq_neg_one : (wall : R) = -1 := by rw [one_eq_neg_wall, neg_neg]

end MARing

/-! ### Every commutative ring is a mult-add ring, and back and forth is the identity -/

/-- A commutative ring gives a mult-add ring: `madd a b = a ✳ b` and `wall = -1`. -/
abbrev ofCommRing (R : Type*) [CommRing R] : MARing R where
  madd a b := a ✳ b
  wall := -1
  madd_comm := multAdd_comm
  madd_assoc := multAdd_assoc
  madd_wall := multAdd_neg_one
  madd_add a b c := by simp only [multAdd_def]; ring

/-- Ring → mult-add ring → product: the product comes back unchanged. -/
theorem roundtrip_mul (R : Type*) [CommRing R] (a b : R) :
    (ofCommRing R).instMul.mul a b = a * b := by
  change a ✳ b - a - b = a * b
  rw [multAdd_def]; ring

/-- Mult-add ring → ring → mult-add: `madd` comes back unchanged. -/
theorem roundtrip_madd (R : Type*) [MARing R] (a b : R) :
    (ofCommRing R).madd a b = MARing.madd a b := MARing.multAdd_eq a b

/-- Mult-add ring → ring → mult-add ring: the wall comes back unchanged. -/
theorem roundtrip_wall (R : Type*) [MARing R] : (ofCommRing R).wall = (MARing.wall : R) :=
  MARing.wall_eq_neg_one.symm

/-- Two mult-add ring structures with the same abelian group, mult-add and wall are equal. -/
theorem MARing.ext' {R : Type*} {i j : MARing R} (h : i.toAddCommGroup = j.toAddCommGroup)
    (hm : i.madd = j.madd) (hw : i.wall = j.wall) : i = j := by
  obtain ⟨⟩ := i
  obtain ⟨⟩ := j
  congr

/-- Ring → mult-add ring → ring gives back the ring structure we started with. -/
theorem instCommRing_ofCommRing (R : Type*) [i : CommRing R] :
    (ofCommRing R).instCommRing = i :=
  CommRing.ext rfl (by ext a b; exact roundtrip_mul R a b)

/-- Mult-add ring → ring → mult-add ring gives back the structure we started with. -/
theorem ofCommRing_instCommRing (R : Type*) [i : MARing R] : ofCommRing R = i :=
  MARing.ext' rfl (by ext a b; exact roundtrip_madd R a b) (roundtrip_wall R)

/-- Mult-add ring structures and commutative ring structures on a type are the same thing. -/
def equivCommRing (R : Type*) : MARing R ≃ CommRing R where
  toFun i := i.instCommRing
  invFun i := @ofCommRing R i
  left_inv i := @ofCommRing_instCommRing R i
  right_inv i := @instCommRing_ofCommRing R i

/-! ### Fields: one more axiom, the way back -/

/-- A mult-add field: a mult-add ring that is not flat (`0 ≠ wall`) and in which every box off
the wall can be undone. -/
class MAField (R : Type*) extends MARing R where
  /-- F1: not everything is flat. -/
  zero_ne_wall : (0 : R) ≠ wall
  /-- F2: every element off the wall has a way back. -/
  way_back : ∀ a, a ≠ wall → ∃ b, madd a b = 0

namespace MAField

open MARing

variable {R : Type*} [MAField R]

theorem add_wall_ne {a : R} (h : a ≠ 0) : a + wall ≠ wall := by simpa using h

protected theorem zero_ne_one : (0 : R) ≠ 1 := fun e ↦
  zero_ne_wall (R := R) <| by rw [wall_eq_neg_one, ← e, neg_zero]

open Classical in
/-- The reciprocal, built from the way back: if `madd (a - 1) b = 0` then `a * (1 + b) = 1`. -/
noncomputable def inv (a : R) : R :=
  if h : a = 0 then 0 else 1 + Classical.choose (way_back (a + wall) (add_wall_ne h))

theorem mul_inv {a : R} (h : a ≠ 0) : a * inv a = 1 := by
  have hb := Classical.choose_spec (way_back (a + wall) (add_wall_ne h))
  rw [inv, dite_eq_right h]
  generalize Classical.choose (way_back (a + wall) (add_wall_ne h)) = b at hb
  rw [← multAdd_eq, wall_eq_neg_one, multAdd_def] at hb
  linear_combination hb

theorem inv_zero : inv (0 : R) = 0 := dite_eq_left rfl

/-- A mult-add field is a field. -/
noncomputable instance instField : Field R where
  __ := MARing.instCommRing (R := R)
  inv := inv
  exists_pair_ne := ⟨0, 1, MAField.zero_ne_one⟩
  mul_inv_cancel _ h := mul_inv h
  inv_zero := inv_zero
  nnqsmul := _
  nnqsmul_def _ _ := rfl
  qsmul := _
  qsmul_def _ _ := rfl

/-- The way back, in closed form: the mult-add inverse of `u` is `-u / (1 + u)`. -/
theorem way_back_formula (u : R) (h : 1 + u ≠ 0) : madd u (multAddInv u) = 0 := by
  rw [← multAdd_eq, multAdd_multAddInv u h]

end MAField

/-- A field gives a mult-add field: off the wall `-1` the way back is `multAddInv`. -/
abbrev ofField (K : Type*) [Field K] : MAField K :=
  have h0 : (0 : K) ≠ -1 := (neg_ne_zero.mpr one_ne_zero).symm
  have hb (a : K) (h : a ≠ -1) : ∃ b, a ✳ b = 0 :=
    ⟨multAddInv a, multAdd_multAddInv a fun e ↦ h (eq_neg_of_add_eq_zero_right e)⟩
  { ofCommRing K with
    zero_ne_wall := h0
    way_back := hb }

/-- Field → mult-add field → field gives back the ring structure we started with (and so the
same inverse, since inverses in a field are unique). -/
theorem toCommRing_instField_ofField (K : Type*) [i : Field K] :
    (ofField K).instField.toCommRing = i.toCommRing :=
  instCommRing_ofCommRing K

end MAFound
