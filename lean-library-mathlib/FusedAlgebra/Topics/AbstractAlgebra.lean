module

public import FusedAlgebra.Basic

/-!
# Abstract algebra through mult-add

The structures that mult-add builds, stated with mathlib's algebraic hierarchy. The core-Lean
original wrote the monoid, group and isomorphism laws as explicit records (`CommMonoidLaws`,
`CommGroupLaws`, `OpIso`); here they are mathlib's `CommMonoid`, `CommGroup` and `≃*`.

Already in `FusedAlgebra.Basic` (and so not repeated): `(R, ✳, 0)` is the commutative monoid
`PreQuasiregular R` (`PreQuasiregular.instCommMonoid`), the shift `x ↦ 1 + x` is the isomorphism
`MultAdd.oneAddEquiv : PreQuasiregular R ≃* R`, and the `✳`-invertible elements are the units of
`R` (`MultAdd.unitsEquiv`); these are mathlib's quasiregular elements (`IsQuasiregular`).

## Main definitions and results

* `multAdd_unit_iff`: `a` has a `✳`-inverse exactly when `1 + a` is a unit.
* `int_multAdd_units`: on `ℤ` the `✳`-units are `0` and `-2`.
* `Weighted k`: a copy of `R` with the weighted mult-add `wAdd k` as multiplication, a commutative
  monoid. `Weighted.scaleHom : Weighted k →* PreQuasiregular R` is `x ↦ k * x`; in a field with
  `k ≠ 0` it is an isomorphism (`Weighted.scaleEquiv`).
* `nilpotent_weight_inverse`, `Weighted.commGroupOfPowEqZero`: if `k ^ m = 0` then
  `(R, wAdd k)` is a commutative group. This is the fact that a nilpotent element is quasiregular
  (mathlib: `IsNilpotent.isUnit_one_add`), which is why a nil ideal sits in the Jacobson radical.
* `zmod4_not_mulEquiv_add`: `(ZMod 4, wAdd 2)` is the Klein four-group, not the cyclic group.
* `Fused R`: a copy of `R` with addition `fadd a b = a + b + 1` and multiplication `✳`, a
  commutative ring with zero `-1` and one `0`, isomorphic to `R` by the shift
  (`Fused.shiftRingEquiv`).
* `mpow_multAdd`: each mult-add power is an endomorphism of mult-add.
* `tateCon q`: the Tate gluing `1 + a ∼ q * (1 + a)` as a congruence on `PreQuasiregular R`; its
  quotient `TateQ q` is a commutative monoid.
* `exists_slope_wall_iff`: slope addition `(a + b) / (1 - a * b)` has walls mod a prime `p`
  exactly when `p % 4 ≠ 3` (from mathlib's `ZMod.exists_sq_eq_neg_one_iff`).

All of this is known ring theory (the circle operation, quasiregularity, radical rings); the
file records how it reads for mult-add.
-/

@[expose] public section

namespace AbstractAlgebra

open MultAdd PreQuasiregular Finset

/-! ### Units of mult-add -/

section Ring

variable {R : Type*} [Ring R]

/-- The shift `MultAdd.oneAddEquiv` sends the wall `-1` to `0`, the absorbing element of
multiplication. -/
theorem shift_wall : oneAddEquiv (equiv (-1 : R)) = 0 := by simp

end Ring

section CommRing

variable {R : Type*} [CommRing R]

/-- Units of `✳` are the quasiregular elements: `a` has a `✳`-inverse exactly when `1 + a` is a
unit. -/
theorem multAdd_unit_iff (a : R) : (∃ b, a ✳ b = 0) ↔ IsUnit (1 + a) := by
  rw [isUnit_iff_exists_inv]
  constructor
  · rintro ⟨b, hb⟩
    exact ⟨1 + b, by rw [← one_add_multAdd, hb, add_zero]⟩
  · rintro ⟨c, hc⟩
    exact ⟨c - 1, by rw [multAdd_def]; linear_combination hc⟩

/-- Having a `✳`-inverse is mathlib's `IsQuasiregular`. -/
theorem multAdd_unit_iff_isQuasiregular (a : R) : (∃ b, a ✳ b = 0) ↔ IsQuasiregular a := by
  rw [multAdd_unit_iff, isQuasiregular_iff_isUnit]

/-- The same for the weighted law, one direction: if `a` has a `wAdd k`-inverse then `1 + k * a`
is a unit. -/
theorem isUnit_one_add_mul_of_wAdd_eq_zero (k a : R) (h : ∃ b, wAdd k a b = 0) :
    IsUnit (1 + k * a) := by
  obtain ⟨b, hb⟩ := h
  exact isUnit_iff_exists_inv.2 ⟨1 + k * b, by rw [← one_add_mul_wAdd, hb, mul_zero, add_zero]⟩

end CommRing

/-! ### Integers and fields -/

/-- On the integers the `✳`-units are exactly `0` and `-2` (the shifts of the units `1` and
`-1`). -/
theorem int_multAdd_units (a : ℤ) : (∃ b, a ✳ b = 0) ↔ a = 0 ∨ a = -2 := by
  rw [multAdd_unit_iff, Int.isUnit_iff]; lia

/-- `-2` is its own `✳`-inverse: the integer group of `✳`-units is the two-element group. -/
theorem neg_two_self_inverse : (-2 : ℤ) ✳ -2 = 0 := by decide

/-- `(ℤ, ✳)` is not a group: `1` has no inverse. -/
theorem int_not_group : ¬∀ a : ℤ, ∃ b, a ✳ b = 0 := fun h ↦ by
  have := (int_multAdd_units 1).1 (h 1)
  lia

/-- In a field every element off the wall is a `✳`-unit. -/
theorem field_unit {K : Type*} [Field K] (a : K) (h : 1 + a ≠ 0) : ∃ b, a ✳ b = 0 :=
  ⟨multAddInv a, multAdd_multAddInv a h⟩

/-! ### The weighted family as monoids -/

/-- A copy of `R` whose multiplication is the weighted mult-add `wAdd k` and whose identity is
`0`. For `k = 1` this is mathlib's `PreQuasiregular R`. -/
@[ext]
structure Weighted {R : Type*} (k : R) where
  /-- The value wrapped into a term of `Weighted k`. -/
  val : R

namespace Weighted

section CommSemiring

variable {R : Type*} [CommSemiring R] {k : R}

/-- `(R, wAdd k, 0)` is a commutative monoid. -/
instance instCommMonoid : CommMonoid (Weighted k) where
  one := ⟨0⟩
  mul x y := ⟨wAdd k x.val y.val⟩
  mul_assoc x y z := Weighted.ext (wAdd_assoc k x.val y.val z.val)
  one_mul x := Weighted.ext (zero_wAdd k x.val)
  mul_one x := Weighted.ext (wAdd_zero k x.val)
  mul_comm x y := Weighted.ext (wAdd_comm k x.val y.val)

@[simp] theorem val_one : (1 : Weighted k).val = 0 := rfl

@[simp] theorem val_mul (x y : Weighted k) : (x * y).val = wAdd k x.val y.val := rfl

variable (k) in
/-- Scaling `x ↦ k * x` is a monoid homomorphism from `(R, wAdd k)` to `(R, ✳)`: weight `k` is
mult-add in other units. -/
def scaleHom : Weighted k →* PreQuasiregular R where
  toFun x := equiv (k * x.val)
  map_one' := by rw [val_one, mul_zero]; rfl
  map_mul' x y := by rw [val_mul, mul_wAdd, equiv_multAdd]

@[simp] theorem val_scaleHom (x : Weighted k) : (scaleHom k x).val = k * x.val := rfl

/-- Where nonzero factors cancel (for instance on `ℤ`), scaling by `k ≠ 0` is injective: weight
`k` is a submonoid of mult-add. -/
theorem scaleHom_injective [IsLeftCancelMulZero R] (hk : k ≠ 0) :
    Function.Injective (scaleHom k) := fun _ _ h ↦
  Weighted.ext (mul_left_cancel₀ hk (congrArg PreQuasiregular.val h))

end CommSemiring

section Field

variable {K : Type*} [Field K]

/-- In a field, weighted laws with `k ≠ 0` are all mult-add in other units: `x ↦ k * x` is an
isomorphism. -/
def scaleEquiv (k : K) (hk : k ≠ 0) : Weighted k ≃* PreQuasiregular K where
  __ := scaleHom k
  invFun y := ⟨y.val / k⟩
  left_inv x := Weighted.ext (mul_div_cancel_left₀ x.val hk)
  right_inv y := congrArg PreQuasiregular.mk (mul_div_cancel₀ y.val hk)

/-- The scaling carries the weighted wall `-1 / k` to the mult-add wall `-1`. -/
theorem scale_wall (k : K) (hk : k ≠ 0) : scaleEquiv k hk ⟨-1 / k⟩ = equiv (-1) :=
  congrArg PreQuasiregular.mk (mul_div_cancel₀ (-1) hk)

end Field

end Weighted

/-- Weight `2`: `-1` is its own inverse. Weight `-2`: `1` is its own inverse. -/
theorem int_weight_units : wAdd (2 : ℤ) (-1) (-1) = 0 ∧ wAdd (-2 : ℤ) 1 1 = 0 := by decide

/-- Weight `-2` on `{0, 1}` is exclusive or. -/
theorem parity_is_xor (x y : Bool) :
    wAdd (-2 : ℤ) (if x then 1 else 0) (if y then 1 else 0) = if (x ^^ y) then 1 else 0 := by
  cases x <;> cases y <;> decide

/-! ### Radical rings: nilpotent weights wrap every element into a group -/

section Radical

variable {R : Type*} [CommRing R]

/-- If `k ^ m = 0`, every `a` has a `wAdd k`-inverse, namely
`-a * (1 - k * a + (k * a) ^ 2 - ⋯)`. This is why a nil ideal sits inside the Jacobson radical;
compare mathlib's `IsNilpotent.isUnit_one_add`. -/
theorem nilpotent_weight_inverse (k a : R) (m : ℕ) (hk : k ^ m = 0) :
    wAdd k a (-a * ∑ i ∈ range m, (-(k * a)) ^ i) = 0 := by
  have h₁ := mul_neg_geom_sum (-(k * a)) m
  have h₂ : (-(k * a)) ^ m = 0 := by rw [← mul_neg, mul_pow, hk, zero_mul]
  rw [wAdd_def]
  linear_combination -a * h₁ + a * h₂

/-- So with a nilpotent weight, `(R, wAdd k)` is a commutative group. -/
abbrev Weighted.commGroupOfPowEqZero (k : R) (m : ℕ) (hk : k ^ m = 0) :
    CommGroup (Weighted k) where
  inv x := ⟨-x.val * ∑ i ∈ range m, (-(k * x.val)) ^ i⟩
  inv_mul_cancel x := Weighted.ext <| by
    rw [Weighted.val_mul, wAdd_comm]; exact nilpotent_weight_inverse k x.val m hk

end Radical

/-! ### `ZMod n`: the exact criterion and the groups it gives -/

/-- `(ZMod 12, wAdd k)` is a group exactly when `6 = 2 * 3` divides `k`. -/
theorem zmod12_group_criterion :
    ∀ k : ZMod 12, (∀ a : ZMod 12, ∃ b : ZMod 12, wAdd k a b = 0) ↔ k.val % 6 = 0 := by decide

/-- The same criterion mod `8`: `2` must divide `k`. -/
theorem zmod8_group_criterion :
    ∀ k : ZMod 8, (∀ a : ZMod 8, ∃ b : ZMod 8, wAdd k a b = 0) ↔ k.val % 2 = 0 := by decide

/-- The same criterion mod `9`: `3` must divide `k`. -/
theorem zmod9_group_criterion :
    ∀ k : ZMod 9, (∀ a : ZMod 9, ∃ b : ZMod 9, wAdd k a b = 0) ↔ k.val % 3 = 0 := by decide

/-- The units of `(ZMod 12, ✳)` are `0, 4, 6, 10`: the units `1, 5, 7, 11` of `(ZMod 12, *)`
shifted down. -/
theorem zmod12_units :
    ∀ a : ZMod 12, (∃ b : ZMod 12, a ✳ b = 0) ↔ a = 0 ∨ a = 4 ∨ a = 6 ∨ a = 10 := by decide

/-- `(ZMod 12, ✳)` is never a group: the wall `11 = -1` absorbs. -/
theorem zmod12_wall : ∀ b : ZMod 12, 11 ✳ b = 11 := by decide

/-- `(ZMod 4, wAdd 2)`: every element is its own inverse (the Klein four-group). -/
theorem zmod4_klein : ∀ a : ZMod 4, wAdd 2 a a = 0 := by decide

/-- So `(ZMod 4, wAdd 2)` and `(ZMod 4, +)` are not isomorphic, though both are groups of order
`4`. -/
theorem zmod4_not_mulEquiv_add :
    IsEmpty (Weighted (2 : ZMod 4) ≃* Multiplicative (ZMod 4)) := ⟨fun f ↦ by
  have h : f.symm (Multiplicative.ofAdd 1) * f.symm (Multiplicative.ofAdd 1) = 1 :=
    Weighted.ext (zmod4_klein _)
  rw [← map_mul, ← ofAdd_add, EmbeddingLike.map_eq_one_iff, ofAdd_eq_one] at h
  exact absurd h (by decide)⟩

/-- Mod `9` with weight `3`: `1` generates all nine elements (cyclic, like `(ZMod 9, +)`). -/
theorem zmod9_cyclic :
    wpow (3 : ZMod 9) 1 9 = 0 ∧ ∀ j : Fin 9, j.val ≠ 0 → wpow (3 : ZMod 9) 1 j.val ≠ 0 := by
  decide

/-- Mod `8` with weight `4`: cyclic of order `8`. -/
theorem zmod8_weight4_cyclic : wpow (4 : ZMod 8) 1 8 = 0 ∧ wpow (4 : ZMod 8) 1 4 ≠ 0 := by decide

/-- Mod `8` with weight `2`: every element has order at most `4`, so the group is not cyclic. -/
theorem zmod8_weight2_not_cyclic : ∀ a : ZMod 8, wpow (2 : ZMod 8) a 4 = 0 := by decide

/-! ### The fused rung pair is a ring -/

section FusedRing

variable {R : Type*} [CommRing R]

/-- Rung 1 on a ring: shifted addition `a + b + 1`. -/
def fadd (a b : R) : R := a + b + 1

theorem fadd_def (a b : R) : fadd a b = a + b + 1 := rfl

/-- The lift of shifted addition: `1 + fadd a b = (1 + a) + (1 + b)`. -/
theorem one_add_fadd (a b : R) : 1 + fadd a b = 1 + a + (1 + b) := by
  rw [fadd_def]; ring

/-- Mult-add distributes over shifted addition. -/
theorem multAdd_fadd (a b c : R) : a ✳ fadd b c = fadd (a ✳ b) (a ✳ c) := by
  simp only [multAdd_def, fadd_def]; ring

/-- A copy of `R` carrying the fused ring structure: addition is `fadd`, multiplication is
mult-add, zero is `-1` and one is `0`. -/
@[ext]
structure Fused (R : Type*) where
  /-- The value wrapped into a term of `Fused R`. -/
  val : R

namespace Fused

/-- The shift `x ↦ 1 + x` from `Fused R` to `R`, as an equivalence. -/
def shift : Fused R ≃ R where
  toFun x := 1 + x.val
  invFun y := ⟨y - 1⟩
  left_inv x := Fused.ext (add_sub_cancel_left 1 x.val)
  right_inv y := add_sub_cancel 1 y

/-- The fused rung pair `(fadd, ✳)` is a commutative ring: the ring structure of `R` carried back
along the shift. -/
instance instCommRing : CommRing (Fused R) := shift.commRing

/-- The shift is a ring isomorphism `(R, fadd, ✳) ≃+* (R, +, *)`. -/
def shiftRingEquiv : Fused R ≃+* R := shift.ringEquiv

@[simp] theorem shiftRingEquiv_apply (x : Fused R) : shiftRingEquiv x = 1 + x.val := rfl

theorem val_add (x y : Fused R) : (x + y).val = fadd x.val y.val := by
  change 1 + x.val + (1 + y.val) - 1 = _
  rw [fadd_def]; ring

theorem val_mul (x y : Fused R) : (x * y).val = x.val ✳ y.val := by
  change (1 + x.val) * (1 + y.val) - 1 = _
  rw [multAdd_eq_mul_sub_one]

/-- The fused zero is the wall `-1`. -/
theorem val_zero : (0 : Fused R).val = -1 := zero_sub 1

/-- The fused one is `0`, the identity of mult-add. -/
theorem val_one : (1 : Fused R).val = 0 := sub_self 1

/-- The fused negative: `fadd a (-a - 2) = -1`. -/
theorem val_neg (x : Fused R) : (-x).val = -x.val - 2 := by
  change -(1 + x.val) - 1 = _
  ring

end Fused

/-- Plain addition and mult-add do not form a ring: `✳` does not distribute over `+`. -/
theorem multAdd_not_distrib_add : (1 : ℤ) ✳ (1 + 1) ≠ 1 ✳ 1 + 1 ✳ 1 := by decide

end FusedRing

/-! ### Mult-add powers: the natural numbers act on `(R, ✳)` -/

section Powers

variable {R : Type*} [NonUnitalCommSemiring R]

/-- `mpow (a ✳ b) n = mpow a n ✳ mpow b n`: each power map is an endomorphism of mult-add. This
is `mul_pow` in the commutative monoid `PreQuasiregular R` (the endomorphism is mathlib's
`powMonoidHom n`). -/
theorem mpow_multAdd (a b : R) (n : ℕ) : mpow (a ✳ b) n = mpow a n ✳ mpow b n := by
  rw [mpow, equiv_multAdd, mul_pow, val_mul_eq_multAdd]; rfl

/-- The identity `0` has all its powers equal to `0`. -/
@[simp]
theorem zero_mpow (n : ℕ) : mpow (0 : R) n = 0 :=
  congrArg PreQuasiregular.val (one_pow (M := PreQuasiregular R) n)

end Powers

/-! ### The Tate gluing is a congruence -/

section Tate

variable {R : Type*} [CommSemiring R]

/-- `a ∼ b` when their lifts differ by powers of `q`: `(1 + a) * q ^ m = (1 + b) * q ^ n`. -/
def tateRel (q a b : R) : Prop := ∃ m n : ℕ, (1 + a) * q ^ m = (1 + b) * q ^ n

theorem tateRel_refl (q a : R) : tateRel q a a := ⟨0, 0, rfl⟩

theorem tateRel_symm (q : R) {a b : R} (h : tateRel q a b) : tateRel q b a := by
  obtain ⟨m, n, h⟩ := h; exact ⟨n, m, h.symm⟩

theorem tateRel_trans (q : R) {a b c : R} (h₁ : tateRel q a b) (h₂ : tateRel q b c) :
    tateRel q a c := by
  obtain ⟨m, n, h₁⟩ := h₁
  obtain ⟨r, s, h₂⟩ := h₂
  refine ⟨m + r, s + n, ?_⟩
  calc (1 + a) * q ^ (m + r) = (1 + a) * q ^ m * q ^ r := by ring
    _ = (1 + b) * q ^ r * q ^ n := by rw [h₁]; ring
    _ = (1 + c) * q ^ (s + n) := by rw [h₂]; ring

/-- The glued point: `q - 1` sits with `0`, that is, multiplying the lift by `q` changes
nothing. -/
theorem tate_glue (q a : R) (h : 1 + a = q) : tateRel q a 0 := ⟨0, 1, by simp [h]⟩

/-- Mult-add respects the gluing on each side. -/
theorem tateRel_multAdd (q : R) {a a' b b' : R} (ha : tateRel q a a') (hb : tateRel q b b') :
    tateRel q (a ✳ b) (a' ✳ b') := by
  obtain ⟨m, n, ha⟩ := ha
  obtain ⟨r, s, hb⟩ := hb
  refine ⟨m + r, n + s, ?_⟩
  rw [one_add_multAdd, one_add_multAdd]
  calc (1 + a) * (1 + b) * q ^ (m + r) = (1 + a) * q ^ m * ((1 + b) * q ^ r) := by ring
    _ = (1 + a') * (1 + b') * q ^ (n + s) := by rw [ha, hb]; ring

/-- The Tate gluing as a congruence on the monoid `(R, ✳)`. -/
def tateCon (q : R) : Con (PreQuasiregular R) where
  r x y := tateRel q x.val y.val
  iseqv := ⟨fun x ↦ tateRel_refl q x.val, tateRel_symm q, tateRel_trans q⟩
  mul' {w x y z} h₁ h₂ := by
    rw [val_mul_eq_multAdd, val_mul_eq_multAdd]; exact tateRel_multAdd q h₁ h₂

theorem tateCon_apply (q : R) (x y : PreQuasiregular R) :
    tateCon q x y ↔ tateRel q x.val y.val := Iff.rfl

/-- The Tate quotient: the points of mult-add with `1 + a` and `q * (1 + a)` glued. It is a
commutative monoid under mult-add (mathlib's `Con.commMonoid`). -/
abbrev TateQ (q : R) : Type _ := (tateCon q).Quotient

/-- In the quotient the class of `q - 1` is the identity. -/
theorem tateQ_glued (q a : R) (h : 1 + a = q) : ((equiv a : PreQuasiregular R) : TateQ q) = 1 :=
  (Con.eq _).2 (tate_glue q a h)

/-- In a field, every class off the wall has an inverse in the Tate quotient. -/
theorem tateQ_inverse {K : Type*} [Field K] (q a : K) (h : 1 + a ≠ 0) :
    ((equiv a : PreQuasiregular K) : TateQ q) * (equiv (multAddInv a) : PreQuasiregular K) = 1 := by
  rw [← Con.coe_mul, ← equiv_multAdd, multAdd_multAddInv a h]; rfl

end Tate

/-! ### Landscape countries change with the prime

The walls of slope addition `(a + b) / (1 - a * b)` are the roots of `1 + a ^ 2`. Mod `5` and `13`
there are two (two-walled country, a split torus with `p - 1` points); mod `3`, `7`, `11` none
(turning country, a non-split torus with `p + 1` points). -/

theorem slope_walls_mod5 : ∀ a : ZMod 5, 1 + a ^ 2 = 0 ↔ a = 2 ∨ a = 3 := by decide

theorem slope_walls_mod13 : ∀ a : ZMod 13, 1 + a ^ 2 = 0 ↔ a = 5 ∨ a = 8 := by decide

theorem slope_no_walls_mod3 : ∀ a : ZMod 3, 1 + a ^ 2 ≠ 0 := by decide

theorem slope_no_walls_mod7 : ∀ a : ZMod 7, 1 + a ^ 2 ≠ 0 := by decide

theorem slope_no_walls_mod11 : ∀ a : ZMod 11, 1 + a ^ 2 ≠ 0 := by decide

/-- The general rule: slope addition has walls mod a prime `p` exactly when `p % 4 ≠ 3` (that is,
`p = 2` or `p ≡ 1 mod 4`). This is mathlib's `ZMod.exists_sq_eq_neg_one_iff`, Fermat's criterion
for `-1` to be a square. -/
theorem exists_slope_wall_iff (p : ℕ) [Fact p.Prime] :
    (∃ a : ZMod p, 1 + a ^ 2 = 0) ↔ p % 4 ≠ 3 := by
  rw [← ZMod.exists_sq_eq_neg_one_iff]
  constructor <;> rintro ⟨a, ha⟩ <;> exact ⟨a, by linear_combination -ha⟩

/-- Mult-add mod `p` always has exactly one finite wall, `-1` (the other wall is at infinity). -/
theorem multAdd_one_wall_mod7 : ∀ a : ZMod 7, 1 + a = 0 ↔ a = 6 := by decide

end AbstractAlgebra
