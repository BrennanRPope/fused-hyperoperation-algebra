module

public import Mathlib

/-!
# Mult-add: the circle operation `a ✳ b = a + b + a * b`

Mult-add is the operation `a ✳ b = a + b + a * b`. Ring theorists call it the *circle operation*
(or quasi-multiplication); mathlib packages it as the monoid structure on `PreQuasiregular R`.
The geometric picture is a box with sides `1 + a` and `1 + b`, with the unit square removed:
`1 + (a ✳ b) = (1 + a) * (1 + b)`.

This file is the mathlib form of `FusedAlgebra.Base` from the core-Lean library. Everything is
built on mathlib's `PreQuasiregular`:

* `MultAdd.multAdd`, written `a ✳ b` (scoped notation): the operation on the ring itself.
* `PreQuasiregular.instCommMonoid`: `(R, ✳, 0)` is a commutative monoid when `R` is commutative.
* `MultAdd.oneAddHom`: the lift `x ↦ 1 + x` as a monoid homomorphism `(R, ✳) →* (R, *)`.
* `MultAdd.oneAddEquiv`: in a ring the lift is an isomorphism of monoids.
* `MultAdd.unitsEquiv`: the group of `✳`-invertible elements is the unit group `Rˣ`.
* `MultAdd.mpow`: the mult-add power, `1 + mpow a n = (1 + a) ^ n`.
* `MultAdd.wAdd`, `MultAdd.wpow`: the weighted family `a + b + k * a * b`.
* `MultAdd.multAddList`: mult-add of a list, with the lift and the corner expansion.
* `MultAdd.multAddInv`, `MultAdd.msub`: the way back, in a field, away from the wall `a = -1`.

## Notation

* `a ✳ b` for `MultAdd.multAdd a b`, in the scope `MultAdd`. Mathlib uses `⊞` for biproducts, so
  the mult-add symbol is `✳`.
-/

@[expose] public section

namespace MultAdd

variable {R : Type*}

/-! ### The operation -/

section Defs

variable [Add R] [Mul R]

/-- Mult-add: `a ✳ b = a + b + a * b`. This is the circle operation of ring theory. -/
def multAdd (a b : R) : R := a + b + a * b

@[inherit_doc] scoped infixl:65 " ✳ " => multAdd

theorem multAdd_def (a b : R) : a ✳ b = a + b + a * b := rfl

end Defs

section NonUnitalNonAssocSemiring

variable [NonUnitalNonAssocSemiring R]

/-- `0` is the identity of mult-add, because `1` is the identity of multiplication. -/
@[simp]
theorem multAdd_zero (a : R) : a ✳ 0 = a := by simp [multAdd]

@[simp]
theorem zero_multAdd (a : R) : 0 ✳ a = a := by simp [multAdd]

instance : Std.LawfulIdentity (multAdd (R := R)) 0 where
  left_id := zero_multAdd
  right_id := multAdd_zero

end NonUnitalNonAssocSemiring

/-! ### Mult-add is the monoid structure of `PreQuasiregular` -/

section NonUnitalSemiring

variable [NonUnitalSemiring R]

open PreQuasiregular

/-- Multiplication in mathlib's `PreQuasiregular R` is mult-add. -/
theorem val_mul_eq_multAdd (x y : PreQuasiregular R) : (x * y).val = x.val ✳ y.val := by
  rw [val_mul, multAdd, add_comm y.val]

/-- Mult-add on `R` is multiplication in `PreQuasiregular R`. -/
theorem equiv_multAdd (a b : R) : equiv (a ✳ b) = equiv a * equiv b :=
  equiv.symm.injective (val_mul_eq_multAdd (equiv a) (equiv b)).symm

theorem multAdd_assoc (a b c : R) : a ✳ b ✳ c = a ✳ (b ✳ c) :=
  equiv.injective <| by simp only [equiv_multAdd, mul_assoc]

instance : Std.Associative (multAdd (R := R)) := ⟨multAdd_assoc⟩

/-- The mult-add power `mpow a n`: `n` copies of `a`, mult-added. It is the `n`-th power of `a` in
the monoid `PreQuasiregular R`. -/
def mpow (a : R) (n : ℕ) : R := (equiv a ^ n).val

@[simp]
theorem mpow_zero (a : R) : mpow a 0 = 0 := rfl

theorem mpow_succ (a : R) (n : ℕ) : mpow a (n + 1) = mpow a n ✳ a := by
  rw [mpow, pow_succ, val_mul_eq_multAdd]; rfl

@[simp]
theorem mpow_one (a : R) : mpow a 1 = a := by simp [mpow]

/-- Powers add: `mpow a (m + n) = mpow a m ✳ mpow a n`. -/
theorem mpow_add (a : R) (m n : ℕ) : mpow a (m + n) = mpow a m ✳ mpow a n := by
  rw [mpow, pow_add, val_mul_eq_multAdd]; rfl

/-- Powers of powers: `mpow a (m * n) = mpow (mpow a m) n`. -/
theorem mpow_mul (a : R) (m n : ℕ) : mpow a (m * n) = mpow (mpow a m) n := by
  rw [mpow, pow_mul]; rfl

/-- The mult-add of a list, `x₁ ✳ x₂ ✳ ⋯ ✳ xₙ`: the product of the list in `PreQuasiregular R`. -/
def multAddList (l : List R) : R := (l.map equiv).prod.val

@[simp]
theorem multAddList_nil : multAddList ([] : List R) = 0 := rfl

@[simp]
theorem multAddList_cons (a : R) (l : List R) : multAddList (a :: l) = a ✳ multAddList l := by
  rw [multAddList, List.map_cons, List.prod_cons, val_mul_eq_multAdd]; rfl

/-- The whole is the first part mult-added to the second part. -/
theorem multAddList_append (l₁ l₂ : List R) :
    multAddList (l₁ ++ l₂) = multAddList l₁ ✳ multAddList l₂ := by
  rw [multAddList, List.map_append, List.prod_append, val_mul_eq_multAdd]; rfl

theorem multAddList_replicate (a : R) (n : ℕ) : multAddList (List.replicate n a) = mpow a n := by
  rw [multAddList, List.map_replicate, List.prod_replicate]; rfl

end NonUnitalSemiring

section NonUnitalCommSemiring

variable [NonUnitalCommSemiring R]

open PreQuasiregular

theorem multAdd_comm (a b : R) : a ✳ b = b ✳ a := by
  simp only [multAdd, add_comm a b, mul_comm a b]

instance : Std.Commutative (multAdd (R := R)) := ⟨multAdd_comm⟩

/-- In a commutative (non-unital) semiring, `(R, ✳, 0)` is a commutative monoid. -/
instance _root_.PreQuasiregular.instCommMonoid : CommMonoid (PreQuasiregular R) where
  mul_comm x y := equiv.symm.injective <| by
    simpa only [equiv_symm_apply, val_mul_eq_multAdd] using multAdd_comm x.val y.val

theorem multAdd_left_comm (a b c : R) : a ✳ (b ✳ c) = b ✳ (a ✳ c) := by
  rw [← multAdd_assoc, multAdd_comm a b, multAdd_assoc]

theorem multAdd_right_comm (a b c : R) : a ✳ b ✳ c = a ✳ c ✳ b := by
  rw [multAdd_assoc, multAdd_comm b c, ← multAdd_assoc]

/-- Order never matters: any reshuffle of a list gives the same mult-add. -/
theorem multAddList_perm {l₁ l₂ : List R} (h : l₁.Perm l₂) : multAddList l₁ = multAddList l₂ :=
  congrArg PreQuasiregular.val (h.map equiv).prod_eq

end NonUnitalCommSemiring

/-! ### The lift `x ↦ 1 + x` carries mult-add to multiplication -/

section Semiring

variable [Semiring R]

open PreQuasiregular

/-- The geometric heart: mult-add is a box with sides `1 + a` and `1 + b`, minus the unit
square. -/
theorem one_add_multAdd (a b : R) : 1 + (a ✳ b) = (1 + a) * (1 + b) := by
  simp only [multAdd, mul_add, add_mul, one_mul, mul_one]; abel

/-- The lift with the `1` on the right (handy on `ℕ`). -/
theorem multAdd_add_one (a b : R) : a ✳ b + 1 = (a + 1) * (b + 1) := by
  rw [add_comm, one_add_multAdd, add_comm 1 a, add_comm 1 b]

/-- The lift `x ↦ 1 + x` as a monoid homomorphism from `(R, ✳)` to `(R, *)`. -/
@[simps]
def oneAddHom : PreQuasiregular R →* R where
  toFun x := 1 + x.val
  map_one' := by simp
  map_mul' x y := by rw [val_mul_eq_multAdd, one_add_multAdd]

/-- `1 + mpow a n = (1 + a) ^ n`. -/
theorem one_add_mpow (a : R) (n : ℕ) : 1 + mpow a n = (1 + a) ^ n :=
  map_pow oneAddHom (equiv a) n

/-- The list lift: `1 + (x₁ ✳ ⋯ ✳ xₙ) = (1 + x₁) * ⋯ * (1 + xₙ)`. -/
theorem one_add_multAddList (l : List R) : 1 + multAddList l = (l.map (1 + ·)).prod := by
  have h := map_list_prod oneAddHom (l.map equiv)
  rwa [List.map_map] at h

/-- The corner expansion: `1 + (x₁ ✳ ⋯ ✳ xₙ)` is the sum over every corner of the cube (every
sublist) of the product of its entries. -/
theorem one_add_multAddList_eq_sum_sublists (l : List R) :
    1 + multAddList l = (l.sublists'.map List.prod).sum := by
  induction l with
  | nil => simp
  | cons a l ih =>
    rw [multAddList_cons, one_add_multAdd, ih, List.sublists'_cons, List.map_append,
      List.sum_append, List.map_map, add_mul, one_mul, ← List.sum_map_mul_left]
    simp [Function.comp_def]

end Semiring

/-! ### The weighted family `a + b + k * a * b` -/

section CommSemiring

variable [CommSemiring R]

/-- Weighted mult-add: the corner carries a coupling weight `k`. -/
def wAdd (k a b : R) : R := a + b + k * a * b

theorem wAdd_def (k a b : R) : wAdd k a b = a + b + k * a * b := rfl

/-- The weighted lift: `1 + k * wAdd k a b = (1 + k * a) * (1 + k * b)`. -/
theorem one_add_mul_wAdd (k a b : R) : 1 + k * wAdd k a b = (1 + k * a) * (1 + k * b) := by
  unfold wAdd; ring

/-- Weighted mult-add is plain mult-add seen through `x ↦ k * x`. -/
theorem mul_wAdd (k a b : R) : k * wAdd k a b = (k * a) ✳ (k * b) := by
  unfold wAdd multAdd; ring

theorem wAdd_comm (k a b : R) : wAdd k a b = wAdd k b a := by
  unfold wAdd; ring

theorem wAdd_assoc (k a b c : R) : wAdd k (wAdd k a b) c = wAdd k a (wAdd k b c) := by
  unfold wAdd; ring

@[simp]
theorem wAdd_zero (k a : R) : wAdd k a 0 = a := by simp [wAdd]

@[simp]
theorem zero_wAdd (k a : R) : wAdd k 0 a = a := by simp [wAdd]

/-- Weight `0`: no corner, plain addition. -/
@[simp]
theorem wAdd_weight_zero (a b : R) : wAdd 0 a b = a + b := by simp [wAdd]

/-- Weight `1`: plain mult-add. -/
@[simp]
theorem wAdd_one (a b : R) : wAdd 1 a b = a ✳ b := by simp [wAdd, multAdd]

/-- The weighted power: `n` copies of `a` under `wAdd k`. -/
def wpow (k a : R) : ℕ → R
  | 0 => 0
  | n + 1 => wAdd k (wpow k a n) a

@[simp]
theorem wpow_zero (k a : R) : wpow k a 0 = 0 := rfl

theorem wpow_succ (k a : R) (n : ℕ) : wpow k a (n + 1) = wAdd k (wpow k a n) a := rfl

/-- `1 + k * wpow k a n = (1 + k * a) ^ n`. -/
theorem one_add_mul_wpow (k a : R) (n : ℕ) : 1 + k * wpow k a n = (1 + k * a) ^ n := by
  induction n with
  | zero => simp
  | succ n ih => rw [wpow_succ, one_add_mul_wAdd, ih, pow_succ]

theorem wpow_one (a : R) (n : ℕ) : wpow 1 a n = mpow a n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [wpow_succ, mpow_succ, ih, wAdd_one]

end CommSemiring

/-! ### Laws that need subtraction -/

section Ring

variable [Ring R]

open PreQuasiregular

/-- In a ring the lift `x ↦ 1 + x` is an isomorphism of monoids from `(R, ✳)` to `(R, *)`. -/
@[simps! apply]
def oneAddEquiv : PreQuasiregular R ≃* R where
  __ := oneAddHom
  invFun y := equiv (y - 1)
  left_inv x := congrArg PreQuasiregular.mk (add_sub_cancel_left 1 x.val)
  right_inv y := by simp

@[simp]
theorem oneAddEquiv_symm_apply (y : R) : oneAddEquiv.symm y = equiv (y - 1) := rfl

/-- The `✳`-invertible elements form a group, and the lift identifies it with the units of `R`. -/
def unitsEquiv : (PreQuasiregular R)ˣ ≃* Rˣ := Units.mapEquiv oneAddEquiv

@[simp]
theorem val_unitsEquiv (u : (PreQuasiregular R)ˣ) : (unitsEquiv u : R) = 1 + u.val.val := rfl

/-- The wall: `-1` absorbs everything mult-added to it. -/
@[simp]
theorem neg_one_multAdd (b : R) : -1 ✳ b = -1 := by simp [multAdd]

@[simp]
theorem multAdd_neg_one (a : R) : a ✳ -1 = -1 := by simp [multAdd]

/-- The power in closed form: `mpow a n = (1 + a) ^ n - 1`. -/
theorem mpow_eq (a : R) (n : ℕ) : mpow a n = (1 + a) ^ n - 1 :=
  eq_sub_of_add_eq' (one_add_mpow a n)

/-- Mult-add is the box minus the unit square. -/
theorem multAdd_eq_mul_sub_one (a b : R) : a ✳ b = (1 + a) * (1 + b) - 1 :=
  eq_sub_of_add_eq' (one_add_multAdd a b)

end Ring

section CommRing

variable [CommRing R]

/-- Weight `-1` is "at least one of two independent events". -/
theorem wAdd_neg_one (a b : R) : wAdd (-1) a b = 1 - (1 - a) * (1 - b) := by
  unfold wAdd; ring

/-- Weight `-1` is mult-add mirrored: `wAdd (-1) a b = -((-a) ✳ (-b))`. -/
theorem wAdd_neg_one_eq_neg_multAdd (a b : R) : wAdd (-1) a b = -(-a ✳ -b) := by
  unfold wAdd multAdd; ring

/-- Flipping both signs flips the weight. -/
theorem wAdd_neg_neg (k a b : R) : wAdd k (-a) (-b) = -wAdd (-k) a b := by
  unfold wAdd; ring

/-- The weighted power at `k = -1`: `n` independent chances, `1 - (1 - a) ^ n`. -/
theorem wpow_neg_one (a : R) (n : ℕ) : wpow (-1) a n = 1 - (1 - a) ^ n := by
  have h := one_add_mul_wpow (-1) a n
  have e : 1 + -1 * a = 1 - a := by ring
  rw [e] at h
  rw [← h]; ring

end CommRing

/-! ### The way back needs division -/

section Field

variable {K : Type*} [Field K]

/-- The mult-add inverse: `-u / (1 + u)`. -/
def multAddInv (u : K) : K := -u / (1 + u)

/-- Everywhere but the wall there is an exact way back. -/
theorem multAdd_multAddInv (u : K) (h : 1 + u ≠ 0) : u ✳ multAddInv u = 0 := by
  unfold multAdd multAddInv; field_simp; ring

theorem multAddInv_multAdd (u : K) (h : 1 + u ≠ 0) : multAddInv u ✳ u = 0 := by
  rw [multAdd_comm]; exact multAdd_multAddInv u h

/-- Under the lift the mult-add inverse is the ordinary reciprocal. -/
theorem inv_one_add (u : K) (h : 1 + u ≠ 0) : (1 + u)⁻¹ = 1 + multAddInv u := by
  unfold multAddInv; field_simp; ring

/-- In a field, an element has a way back exactly when it is off the wall. This is mathlib's
`IsQuasiregular`. -/
theorem isQuasiregular_iff_one_add_ne_zero {u : K} : IsQuasiregular u ↔ 1 + u ≠ 0 := by
  rw [isQuasiregular_iff_isUnit, isUnit_iff_ne_zero]

/-- At the wall there is no way back: nothing mult-adds with `-1` to give `0`. -/
theorem neg_one_multAdd_ne_zero (v : K) : -1 ✳ v ≠ 0 := by simp

/-- Mult-add subtraction: `msub a b = (a - b) / (1 + b)`. -/
def msub (a b : K) : K := (a - b) / (1 + b)

/-- Subtraction undoes mult-add (Fisher's rule for real rates). -/
theorem msub_multAdd (a b : K) (h : 1 + b ≠ 0) : msub a b ✳ b = a := by
  unfold multAdd msub; field_simp; ring

/-- Mult-add subtraction is mult-adding the inverse. -/
theorem msub_eq_multAdd_multAddInv (a b : K) (h : 1 + b ≠ 0) : msub a b = a ✳ multAddInv b := by
  unfold multAdd msub multAddInv; field_simp; ring

/-- The weighted wall: `-1 / k` absorbs everything under `wAdd k`. -/
theorem wAdd_wall (k a : K) (hk : k ≠ 0) : wAdd k a (-1 / k) = -1 / k := by
  unfold wAdd; field_simp; ring

end Field

end MultAdd
