/-
  Axioms.Independence: none of M1–M4 follows from the others
  ===========================================================

  For each mult-add axiom we give a structure where every other axiom holds
  and that one fails. So the list in `Axioms.MultAddRing` has no spare part
  among M1–M4. (A1–A4 are the usual abelian-group axioms; their independence
  is classical and not redone here.)

  * no M1 (turning the box): 2×2 integer matrices, A ◯ B = A + B + AB.
  * no M2 (building a 3-D box two ways): a 3-number algebra with x·x = y,
    y·y = x, x·y = 0, A ◯ B = A + B + A·B.
  * no M3 (the wall): whole numbers with a ◯ b = a + b + 2ab, the weighted
    family at k = 2. Its wall would be −1/2, which is not a whole number.
  * no M4 (cutting a box in two): whole numbers with a ◯ b = a·b.
-/

namespace CircFound.Independence

/-- The eight axioms, as statements about given operations. -/
structure Ops (R : Type) where
  zero : R
  add : R → R → R
  neg : R → R
  circ : R → R → R
  wall : R

variable {R : Type}
def A1 (o : Ops R) : Prop := ∀ a b c, o.add (o.add a b) c = o.add a (o.add b c)
def A2 (o : Ops R) : Prop := ∀ a b, o.add a b = o.add b a
def A3 (o : Ops R) : Prop := ∀ a, o.add a o.zero = a
def A4 (o : Ops R) : Prop := ∀ a, o.add (o.neg a) a = o.zero
def M1 (o : Ops R) : Prop := ∀ a b, o.circ a b = o.circ b a
def M2 (o : Ops R) : Prop := ∀ a b c, o.circ (o.circ a b) c = o.circ a (o.circ b c)
def M3 (o : Ops R) : Prop := ∀ a, o.circ a o.wall = o.wall
def M4 (o : Ops R) : Prop :=
  ∀ a b c, o.circ a (o.add b c) = o.add (o.add (o.circ a b) (o.circ a c)) (o.neg a)

/-- The group part holds for whole numbers with ordinary +. -/
def intOps (m : Int → Int → Int) (w : Int) : Ops Int := ⟨0, (· + ·), (- ·), m, w⟩
theorem int_group (m : Int → Int → Int) (w : Int) :
    A1 (intOps m w) ∧ A2 (intOps m w) ∧ A3 (intOps m w) ∧ A4 (intOps m w) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [A1, A2, A3, A4, intOps] <;> intros <;> omega

/-! ## No M4: a ◯ b = a·b -/

theorem no_M4 :
    let o := intOps (· * ·) 0
    A1 o ∧ A2 o ∧ A3 o ∧ A4 o ∧ M1 o ∧ M2 o ∧ M3 o ∧ ¬ M4 o := by
  refine ⟨(int_group _ _).1, (int_group _ _).2.1, (int_group _ _).2.2.1,
    (int_group _ _).2.2.2, ?_, ?_, ?_, ?_⟩
  · intro a b; simp only [intOps]; exact Int.mul_comm a b
  · intro a b c; simp only [intOps]; exact Int.mul_assoc a b c
  · intro a; simp only [intOps]; exact Int.mul_zero a
  · intro h; have := h 1 0 0; simp only [intOps] at this; omega

/-! ## No M3: the weighted family at k = 2 -/

theorem no_M3 : ∀ w : Int,
    let o := intOps (fun a b => a + b + 2 * a * b) w
    A1 o ∧ A2 o ∧ A3 o ∧ A4 o ∧ M1 o ∧ M2 o ∧ ¬ M3 o ∧ M4 o := by
  intro w
  refine ⟨(int_group _ _).1, (int_group _ _).2.1, (int_group _ _).2.2.1,
    (int_group _ _).2.2.2, ?_, ?_, ?_, ?_⟩
  · intro a b; simp only [intOps]; grind
  · intro a b c; simp only [intOps]; grind
  · intro h; have := h 1; simp only [intOps] at this; omega
  · intro a b c; simp only [intOps]; grind

/-! ## No M1: 2×2 integer matrices -/

structure Mat where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq

def Mat.add (x y : Mat) : Mat := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def Mat.mul (x y : Mat) : Mat :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def Mat.neg (x : Mat) : Mat := ⟨-x.a, -x.b, -x.c, -x.d⟩
def matOps : Ops Mat :=
  ⟨⟨0, 0, 0, 0⟩, Mat.add, Mat.neg, fun x y => Mat.add (Mat.add x y) (Mat.mul x y), ⟨-1, 0, 0, -1⟩⟩

theorem Mat.ext' {x y : Mat} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c) (h4 : x.d = y.d) :
    x = y := by
  cases x; cases y; simp_all

theorem no_M1 :
    A1 matOps ∧ A2 matOps ∧ A3 matOps ∧ A4 matOps ∧ ¬ M1 matOps ∧ M2 matOps ∧ M3 matOps ∧ M4 matOps := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x y z; apply Mat.ext' <;> simp only [matOps, Mat.add] <;> omega
  · intro x y; apply Mat.ext' <;> simp only [matOps, Mat.add] <;> omega
  · intro x; apply Mat.ext' <;> simp only [matOps, Mat.add] <;> omega
  · intro x; apply Mat.ext' <;> simp only [matOps, Mat.add, Mat.neg] <;> omega
  · intro h
    have := h ⟨0, 1, 0, 0⟩ ⟨0, 0, 1, 0⟩
    revert this; decide
  · intro x y z; apply Mat.ext' <;> simp only [matOps, Mat.add, Mat.mul] <;> grind
  · intro x; apply Mat.ext' <;> simp only [matOps, Mat.add, Mat.mul] <;> grind
  · intro x y z; apply Mat.ext' <;> simp only [matOps, Mat.add, Mat.mul, Mat.neg] <;> grind

/-! ## No M2: a commutative product that is not associative -/

structure Tri where
  p : Int
  q : Int
  r : Int
deriving DecidableEq

/-- Numbers p + q·x + r·y with x·x = y, y·y = x, x·y = 0, and 1 as the unit. -/
def Tri.mul (u v : Tri) : Tri :=
  ⟨u.p * v.p, u.p * v.q + u.q * v.p + u.r * v.r, u.p * v.r + u.r * v.p + u.q * v.q⟩
def Tri.add (u v : Tri) : Tri := ⟨u.p + v.p, u.q + v.q, u.r + v.r⟩
def Tri.neg (u : Tri) : Tri := ⟨-u.p, -u.q, -u.r⟩
def triOps : Ops Tri :=
  ⟨⟨0, 0, 0⟩, Tri.add, Tri.neg, fun u v => Tri.add (Tri.add u v) (Tri.mul u v), ⟨-1, 0, 0⟩⟩

theorem Tri.ext' {u v : Tri} (h1 : u.p = v.p) (h2 : u.q = v.q) (h3 : u.r = v.r) : u = v := by
  cases u; cases v; simp_all

theorem no_M2 :
    A1 triOps ∧ A2 triOps ∧ A3 triOps ∧ A4 triOps ∧ M1 triOps ∧ ¬ M2 triOps ∧ M3 triOps ∧ M4 triOps := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x y z; apply Tri.ext' <;> simp only [triOps, Tri.add] <;> omega
  · intro x y; apply Tri.ext' <;> simp only [triOps, Tri.add] <;> omega
  · intro x; apply Tri.ext' <;> simp only [triOps, Tri.add] <;> omega
  · intro x; apply Tri.ext' <;> simp only [triOps, Tri.add, Tri.neg] <;> omega
  · intro x y; apply Tri.ext' <;> simp only [triOps, Tri.add, Tri.mul] <;> grind
  · intro h
    have := h ⟨0, 1, 0⟩ ⟨0, 1, 0⟩ ⟨0, 0, 1⟩
    revert this; decide
  · intro x; apply Tri.ext' <;> simp only [triOps, Tri.add, Tri.mul] <;> grind
  · intro x y z; apply Tri.ext' <;> simp only [triOps, Tri.add, Tri.mul, Tri.neg] <;> grind

end CircFound.Independence
