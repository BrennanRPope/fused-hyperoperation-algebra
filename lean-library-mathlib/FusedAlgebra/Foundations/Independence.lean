module

public import FusedAlgebra.Foundations.MultAddRing

/-!
# Independence of the mult-add axioms

None of the four mult-add axioms M1 to M4 of a mult-add ring (`MAFound.MARing`) follows from the
others. For each of them we give an abelian group with an operation `m` and a constant `w` where
the other three hold and that one fails:

* `no_M1` (turning the box): `2 × 2` integer matrices with `m A B = A ✳ B = A + B + A * B`;
* `no_M2` (building a 3-D box two ways): `ℤ³` with a commutative, non-associative product
  (`x * x = y`, `y * y = x`, `x * y = 0`) and `m u v = u + v + u * v`;
* `no_M3` (the wall): `ℤ` with `m a b = a + b + 2 * a * b`, the weighted family at `k = 2`. Its
  wall would be `-1 / 2`, which is not an integer;
* `no_M4` (cutting a box in two): `ℤ` with `m a b = a * b`.

The abelian group axioms A1 to A4 of the original are mathlib's `AddCommGroup`; every model here
is an `AddCommGroup` by a mathlib instance, so they are not restated. (Their independence is
classical.)

The file also records which axioms hold for free: the circle operation `✳` of any ring satisfies
M2, M3 and M4 (`m2_multAdd`, `m3_multAdd`, `m4_multAdd`), and the weighted mult-add `wAdd k` of a
commutative ring satisfies M1, M2 and M4. This is the mathlib form of
`FusedAlgebra.Foundations.Independence`.
-/

@[expose] public section

namespace MAFound.Independence

open MultAdd

variable {R : Type*}

/-! ### The four axioms, as statements about a given operation -/

/-- M1: the operation is commutative. -/
def M1 (m : R → R → R) : Prop := ∀ a b, m a b = m b a

/-- M2: the operation is associative. -/
def M2 (m : R → R → R) : Prop := ∀ a b c, m (m a b) c = m a (m b c)

/-- M3: `w` absorbs on the right. -/
def M3 (m : R → R → R) (w : R) : Prop := ∀ a, m a w = w

/-- M4: `m a (b + c) = m a b + m a c - a`. -/
def M4 [AddCommGroup R] (m : R → R → R) : Prop := ∀ a b c, m a (b + c) = m a b + m a c - a

/-- A mult-add ring satisfies the four axioms. -/
theorem axioms_of_maRing [MARing R] :
    M1 (MARing.madd (R := R)) ∧ M2 (MARing.madd (R := R)) ∧
      M3 (MARing.madd (R := R)) MARing.wall ∧ M4 (MARing.madd (R := R)) :=
  ⟨MARing.madd_comm, MARing.madd_assoc, MARing.madd_wall, MARing.madd_add⟩

/-- An abelian group with an operation satisfying the four axioms is a mult-add ring: M1 to M4
are exactly the fields of `MARing`. -/
abbrev maRingOfAxioms [AddCommGroup R] (m : R → R → R) (w : R) (h1 : M1 m) (h2 : M2 m)
    (h3 : M3 m w) (h4 : M4 m) : MARing R where
  madd := m
  wall := w
  madd_comm := h1
  madd_assoc := h2
  madd_wall := h3
  madd_add := h4

/-! ### What holds for free -/

/-- The circle operation of any (possibly non-commutative) ring is associative. -/
theorem m2_multAdd [NonUnitalSemiring R] : M2 (multAdd (R := R)) := multAdd_assoc

/-- In any ring `-1` is a wall for the circle operation. -/
theorem m3_multAdd [Ring R] : M3 (multAdd (R := R)) (-1) := multAdd_neg_one

/-- M4 for the circle operation only needs distributivity. -/
theorem m4_multAdd [NonUnitalNonAssocRing R] : M4 (multAdd (R := R)) := fun a b c ↦ by
  simp only [multAdd_def, mul_add]; abel

/-- The weighted mult-add of a commutative ring satisfies M1, M2 and M4 for every weight. -/
theorem wAdd_axioms [CommRing R] (k : R) : M1 (wAdd k) ∧ M2 (wAdd k) ∧ M4 (wAdd k) :=
  ⟨wAdd_comm k, wAdd_assoc k, fun a b c ↦ by simp only [wAdd_def]; ring⟩

/-! ### No M4: `m a b = a * b` -/

/-- Plain multiplication on `ℤ`, with wall `0`, satisfies M1, M2, M3 but not M4. -/
theorem no_M4 : M1 (· * · : ℤ → ℤ → ℤ) ∧ M2 (· * · : ℤ → ℤ → ℤ) ∧
    M3 (· * · : ℤ → ℤ → ℤ) 0 ∧ ¬ M4 (· * · : ℤ → ℤ → ℤ) :=
  ⟨mul_comm, mul_assoc, mul_zero, fun h ↦ by simpa using h 1 0 0⟩

/-! ### No M3: the weighted family at `k = 2` -/

/-- `a + b + 2 * a * b` on `ℤ` satisfies M1, M2, M4, and no integer is a wall for it. -/
theorem no_M3 (w : ℤ) : M1 (wAdd (2 : ℤ)) ∧ M2 (wAdd (2 : ℤ)) ∧ ¬ M3 (wAdd (2 : ℤ)) w ∧
    M4 (wAdd (2 : ℤ)) := by
  obtain ⟨h1, h2, h4⟩ := wAdd_axioms (2 : ℤ)
  refine ⟨h1, h2, fun h ↦ ?_, h4⟩
  have := h 1
  rw [wAdd_def] at this
  lia

/-! ### No M1: `2 × 2` integer matrices -/

/-- The circle operation on `2 × 2` integer matrices, with wall `-1`, satisfies M2, M3, M4 but
not M1. -/
theorem no_M1 : ¬ M1 (multAdd (R := Matrix (Fin 2) (Fin 2) ℤ)) ∧
    M2 (multAdd (R := Matrix (Fin 2) (Fin 2) ℤ)) ∧
    M3 (multAdd (R := Matrix (Fin 2) (Fin 2) ℤ)) (-1) ∧
    M4 (multAdd (R := Matrix (Fin 2) (Fin 2) ℤ)) := by
  refine ⟨fun h ↦ ?_, m2_multAdd, m3_multAdd, m4_multAdd⟩
  have := congrFun (congrFun (h !![0, 1; 0, 0] !![0, 0; 1, 0]) 0) 0
  simp [multAdd_def] at this

/-! ### No M2: a commutative product that is not associative -/

/-- The carrier of the M2 countermodel: triples `(p, q, r)` standing for `p + q * x + r * y`. -/
abbrev Tri : Type := ℤ × ℤ × ℤ

/-- The product on `Tri`: `x * x = y`, `y * y = x`, `x * y = 0`, and `(1, 0, 0)` is the unit. It
is commutative but not associative. -/
def Tri.mul (u v : Tri) : Tri :=
  (u.1 * v.1, u.1 * v.2.1 + u.2.1 * v.1 + u.2.2 * v.2.2,
    u.1 * v.2.2 + u.2.2 * v.1 + u.2.1 * v.2.1)

/-- The mult-add of the M2 countermodel: `u + v + u * v` for the product `Tri.mul`. -/
def Tri.madd (u v : Tri) : Tri := u + v + Tri.mul u v

/-- `Tri.madd`, with wall `(-1, 0, 0)`, satisfies M1, M3, M4 but not M2. -/
theorem no_M2 : M1 Tri.madd ∧ ¬ M2 Tri.madd ∧ M3 Tri.madd (-1, 0, 0) ∧ M4 Tri.madd := by
  refine ⟨fun u v ↦ ?_, fun h ↦ ?_, fun u ↦ ?_, fun u v w ↦ ?_⟩
  · ext <;> simp [Tri.madd, Tri.mul] <;> ring
  · have := h (0, 1, 0) (0, 1, 0) (0, 0, 1)
    revert this; decide
  · ext <;> simp [Tri.madd, Tri.mul]
  · ext <;> simp [Tri.madd, Tri.mul] <;> ring

end MAFound.Independence
