/-
  Fused Hyperoperation Algebra: abstract algebra through the circle operation (called mult-add before 1.0).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base.

  To fold into the library: copy to lean-library/FusedAlgebra/Topics/AbstractAlgebra.lean,
  add `FusedAlgebra.Topics.AbstractAlgebra` to build.sh and FusedAlgebra.lean, run ./build.sh.

  What is here (the structures the circle operation builds):
    1. Laws as explicit records: commutative monoid, group, isomorphism.
    2. (R, ◯) is a commutative monoid on any commutative semiring; on a ring the shift
       x ↦ 1 + x is an isomorphism (R, ◯, 0) ≅ (R, ·, 1) that sends the wall −1 to 0.
       Its units are exactly the a with 1 + a invertible (quasi-regular elements).
    3. On Int the only units are 0 and −2. On a field everything off the wall is a unit.
    4. The weighted family: for k ≠ 0 in a field, x ↦ k·x is an isomorphism ◯[k] ≅ ◯.
       On Int, x ↦ k·x embeds ◯[k] into ◯. Weight −2 on {0, 1} is XOR.
    5. Radical rings: if kᵐ = 0 then every element has a ◯[k] inverse (the Jacobson-radical
       fact), and Fin n examples: the exact group criterion for n = 12, the units of
       (Fin 12, ◯), the Klein four-group (Fin 4, ◯[2]) that is not (Fin 4, +), cyclic and
       non-cyclic cases mod 8 and 9.
    6. The fused rung pair (shifted add, the circle operation) is a commutative ring, zero −1, one 0,
       isomorphic to (+, ·) by the shift.
    7. Circle powers are a module action: (a ◯ b)^◯n = a^◯n ◯ b^◯n, (a^◯m)^◯n = a^◯(mn).
    8. The Tate gluing (1 + a) ~ q(1 + a) is a congruence for ◯; the quotient is a
       commutative monoid, a group off the wall on a field.
    9. The landscape countries mod p: slope addition has walls mod 5 and 13, none mod 3, 7, 11.

  LIBRARY NOTE: copied from abstract-algebra/lean/AbstractAlgebra.lean (Nov 25) as Topics/AbstractAlgebra.
-/
import FusedAlgebra.Base

namespace AbstractAlgebra
open CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The laws, written out as records -/

/-- A commutative monoid law `op` with identity `e`. -/
structure CommMonoidLaws {α : Type} (op : α → α → α) (e : α) : Prop where
  comm : ∀ a b, op a b = op b a
  assoc : ∀ a b c, op (op a b) c = op a (op b c)
  ident : ∀ a, op a e = a

/-- A commutative group law: a commutative monoid where every element has an inverse. -/
structure CommGroupLaws {α : Type} (op : α → α → α) (e : α) : Prop extends CommMonoidLaws op e where
  inv : ∀ a, ∃ b, op a b = e

/-- An isomorphism of operations: a two-sided inverse pair that carries one law to the other
    and identity to identity. -/
structure OpIso (α β : Type) (op₁ : α → α → α) (e₁ : α) (op₂ : β → β → β) (e₂ : β) where
  toFun : α → β
  invFun : β → α
  left : ∀ x, invFun (toFun x) = x
  right : ∀ y, toFun (invFun y) = y
  map_op : ∀ a b, toFun (op₁ a b) = op₂ (toFun a) (toFun b)
  map_id : toFun e₁ = e₂

/-- An isomorphism carries "x op x = e" across, so it preserves elements of order ≤ 2. -/
theorem OpIso.self_op {α β : Type} {op₁ : α → α → α} {e₁ : α} {op₂ : β → β → β} {e₂ : β}
    (f : OpIso α β op₁ e₁ op₂ e₂) (a : α) (h : op₁ a a = e₁) : op₂ (f.toFun a) (f.toFun a) = e₂ := by
  rw [← f.map_op, h, f.map_id]

/-! ## 2. (R, ◯) is a commutative monoid; the shift makes it (R, ·) -/

section Semiring
variable {R : Type} [Lean.Grind.CommSemiring R]

theorem circ_commMonoid : CommMonoidLaws (circ (R := R)) 0 :=
  ⟨circ_comm, circ_assoc, circ_zero⟩

theorem circK_commMonoid (k : R) : CommMonoidLaws (circK k) 0 :=
  ⟨circK_comm k, circK_assoc k, circK_zero k⟩

theorem mul_commMonoid : CommMonoidLaws (fun a b : R => a * b) 1 :=
  ⟨fun a b => by grind, fun a b c => by grind, fun a => by grind⟩

end Semiring

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- The shift x ↦ 1 + x is an isomorphism (R, ◯, 0) ≅ (R, ·, 1). -/
def shiftIso : OpIso R R circ 0 (fun a b => a * b) 1 where
  toFun x := 1 + x
  invFun y := y - 1
  left x := by grind
  right y := by grind
  map_op a b := circ_lift a b
  map_id := by grind

/-- The shift sends the wall −1 to 0, the absorbing element of multiplication. -/
theorem shift_wall : (shiftIso (R := R)).toFun (-1) = 0 := by
  show (1 : R) + -1 = 0; grind

/-- Units of ◯ are the quasi-regular elements: a has a ◯-inverse exactly when 1 + a is
    invertible for multiplication. -/
theorem circ_unit_iff (a : R) : (∃ b, circ a b = 0) ↔ (∃ c, (1 + a) * c = 1) := by
  constructor
  · rintro ⟨b, hb⟩; exact ⟨1 + b, by rw [← circ_lift, hb]; grind⟩
  · rintro ⟨c, hc⟩; exact ⟨c - 1, by unfold circ; grind⟩

/-- The same for the weighted law: a is a ◯[k]-unit exactly when 1 + k·a is invertible. -/
theorem circK_unit_iff (k a : R) : (∃ b, circK k a b = 0) → (∃ c, (1 + k * a) * c = 1) := by
  rintro ⟨b, hb⟩; exact ⟨1 + k * b, by rw [← circK_lift, hb]; grind⟩

end Ring

/-! ## 3. Integers and fields -/

/-- On the integers the ◯-units are exactly 0 and −2 (the shifts of the units 1 and −1). -/
theorem int_circ_units (a : Int) : (∃ b, circ a b = 0) ↔ a = 0 ∨ a = -2 := by
  rw [circ_unit_iff]
  constructor
  · rintro ⟨c, hc⟩
    have h2 := congrArg Int.natAbs hc
    rw [Int.natAbs_mul] at h2
    have h3 : (1 + a).natAbs = 1 := Nat.eq_one_of_mul_eq_one_right h2
    have := Int.natAbs_eq (1 + a)
    omega
  · rintro (h | h)
    · exact ⟨1, by rw [h]; decide⟩
    · exact ⟨-1, by rw [h]; decide⟩

/-- −2 is its own ◯-inverse: the integer group of ◯-units is the two-element group. -/
theorem neg_two_self_inverse : circ (-2 : Int) (-2) = 0 := by decide

/-- (Int, ◯) has no group law: 1 has no inverse. -/
theorem int_not_group : ¬ ∀ a : Int, ∃ b, circ a b = 0 := by
  intro h
  rcases (int_circ_units 1).mp (h 1) with h1 | h1 <;> omega

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- On a field every element off the wall is a ◯-unit. -/
theorem field_unit (a : K) (h : 1 + a ≠ 0) : ∃ b, circ a b = 0 :=
  ⟨circInv a, circ_circInv a h⟩

/-- Weighted laws with k ≠ 0 are all the circle operation in other units: x ↦ k·x is an isomorphism. -/
def scaleIso (k : K) (hk : k ≠ 0) : OpIso K K (circK k) 0 circ 0 where
  toFun x := k * x
  invFun y := y / k
  left x := by grind
  right y := by grind
  map_op a b := circK_scale k a b
  map_id := by grind

/-- The scaling carries the weighted wall −1/k to the circle-operation wall −1. -/
theorem scale_wall (k : K) (hk : k ≠ 0) : (scaleIso k hk).toFun (-1 / k) = -1 := by
  show k * (-1 / k) = -1; grind

end Field

/-- On the integers x ↦ k·x is an injective homomorphism (Int, ◯[k]) → (Int, ◯): weight k is
    the sub-monoid kℤ of the circle operation. -/
theorem int_scale_hom (k a b : Int) : k * circK k a b = circ (k * a) (k * b) := circK_scale k a b

theorem int_scale_inj (k : Int) (hk : k ≠ 0) (a b : Int) (h : k * a = k * b) : a = b :=
  Int.eq_of_mul_eq_mul_left hk h

/-- Weight 2: −1 is a unit (its own inverse). Weight −2: 1 is its own inverse. -/
theorem int_weight_units : circK (2 : Int) (-1) (-1) = 0 ∧ circK (-2 : Int) 1 1 = 0 := by decide

/-- Weight −2 on {0, 1} is exclusive or. -/
theorem parity_is_xor (x y : Bool) :
    circK (-2 : Int) (if x then 1 else 0) (if y then 1 else 0) = if (x ^^ y) then 1 else 0 := by
  cases x <;> cases y <;> decide

/-! ## 4. Radical rings: nilpotent weights wrap every element into a group -/

section Radical
variable {R : Type} [Lean.Grind.CommRing R]

/-- The geometric sum 1 + x + … + x^(m−1). -/
def geom (x : R) : Nat → R
  | 0 => 0
  | m + 1 => geom x m + x ^ m

theorem geom_mul (x : R) (m : Nat) : (1 - x) * geom x m = 1 - x ^ m := by
  induction m with
  | zero => simp only [geom]; rw [Lean.Grind.Semiring.pow_zero]; grind
  | succ m ih => simp only [geom]; rw [Lean.Grind.Semiring.pow_succ]; grind

/-- If kᵐ = 0, every a has a ◯[k]-inverse, namely −a·(1 − ka + (ka)² − …).
    This is why a nilpotent ideal sits inside the Jacobson radical. -/
theorem nilpotent_weight_inverse (k a : R) (m : Nat) (hk : k ^ m = 0) :
    circK k a (-a * geom (-(k * a)) m) = 0 := by
  have h1 := geom_mul (-(k * a)) m
  have h2 : (-(k * a)) ^ m = (-a) ^ m * k ^ m := by
    rw [← mul_pow']; congr 1; grind
  rw [h2, hk] at h1
  unfold circK
  grind

/-- So with a nilpotent weight, (R, ◯[k]) is a commutative group. -/
theorem nilpotent_weight_group (k : R) (m : Nat) (hk : k ^ m = 0) : CommGroupLaws (circK k) 0 :=
  { toCommMonoidLaws := circK_commMonoid k
    inv := fun a => ⟨_, nilpotent_weight_inverse k a m hk⟩ }

end Radical

/-! ### Fin n: the exact criterion and the groups it gives -/

/-- (Fin 12, ◯[k]) is a group exactly when 6 = 2·3 divides k. -/
theorem fin12_group_criterion :
    ∀ k : Fin 12, (∀ a : Fin 12, ∃ b : Fin 12, circK k a b = 0) ↔ k.val % 6 = 0 := by decide

/-- The same criterion mod 8 (2 must divide k) and mod 9 (3 must divide k). -/
theorem fin8_group_criterion :
    ∀ k : Fin 8, (∀ a : Fin 8, ∃ b : Fin 8, circK k a b = 0) ↔ k.val % 2 = 0 := by decide

theorem fin9_group_criterion :
    ∀ k : Fin 9, (∀ a : Fin 9, ∃ b : Fin 9, circK k a b = 0) ↔ k.val % 3 = 0 := by decide

/-- Units of (Fin 12, ◯) are 0, 4, 6, 10: the units 1, 5, 7, 11 of (Fin 12, ·) shifted down. -/
theorem fin12_units : ∀ a : Fin 12,
    (∃ b : Fin 12, circ a b = 0) ↔ (a = 0 ∨ a = 4 ∨ a = 6 ∨ a = 10) := by decide

/-- (Fin 12, ◯) is never a group: the wall 11 = −1 absorbs. -/
theorem fin12_wall : ∀ b : Fin 12, circ 11 b = 11 := by decide

/-- (Fin 4, ◯[2]): every element is its own inverse (the Klein four-group). -/
theorem fin4_klein : ∀ a : Fin 4, circK 2 a a = 0 := by decide

/-- So (Fin 4, ◯[2]) and (Fin 4, +) are not isomorphic, though both are groups of order 4. -/
theorem fin4_not_iso_add : OpIso (Fin 4) (Fin 4) (circK 2) 0 (· + ·) 0 → False := by
  intro f
  have h := f.self_op (f.invFun 1) (fin4_klein _)
  rw [f.right] at h
  exact absurd h (by decide)

/-- Mod 9 with weight 3: 1 generates all nine elements (cyclic, like (Fin 9, +)). -/
theorem fin9_cyclic : wpow (3 : Fin 9) 1 9 = 0 ∧ ∀ j : Fin 9, j.val ≠ 0 → wpow (3 : Fin 9) 1 j.val ≠ 0 := by
  decide

/-- Mod 8 with weight 4: cyclic of order 8. With weight 2: every element has order ≤ 4. -/
theorem fin8_weight4_cyclic : wpow (4 : Fin 8) 1 8 = 0 ∧ wpow (4 : Fin 8) 1 4 ≠ 0 := by decide

theorem fin8_weight2_not_cyclic : ∀ a : Fin 8, wpow (2 : Fin 8) a 4 = 0 := by decide

/-! ## 5. The fused rung pair is a ring -/

section FusedRing
variable {R : Type} [Lean.Grind.CommRing R]

/-- Rung 1 on a ring: shifted addition a + b + 1.
  Lineage: one of five earlier names for the rung-1 "or"; the unified, standalone
    version is `RungOne.sadd` (Foundations.RungOneSemiring), and
    Foundations.RungOneLineage proves the two agree. Standard name: shifted add, the
    "or" of the ring (R, a+b+1, a+b+ab) obtained by transport along x ↦ x − 1. -/
def fadd (a b : R) : R := a + b + 1

theorem fadd_lift (a b : R) : 1 + fadd a b = (1 + a) + (1 + b) := by unfold fadd; grind

theorem fadd_group : CommGroupLaws (fadd (R := R)) (-1) :=
  { comm := fun a b => by unfold fadd; grind
    assoc := fun a b c => by unfold fadd; grind
    ident := fun a => by unfold fadd; grind
    inv := fun a => ⟨-a - 2, by unfold fadd; grind⟩ }

/-- The circle operation distributes over shifted addition.
  Lineage: one of four earlier proofs of the rung-1 distributive law; the unified,
    standalone proof is `RungOne.circ_sadd_distrib` (Foundations.RungOneSemiring), and
    Foundations.RungOneLineage rederives this one from it. Standard name: the distributive
    law of the transported ring (generated semiring with g(x) = 1 + x). -/
theorem circ_fadd_distrib (a b c : R) : circ a (fadd b c) = fadd (circ a b) (circ a c) := by
  unfold circ fadd; grind

/-- The fused zero −1 absorbs under the fused product (it is the wall). -/
theorem fused_zero_absorbs (a : R) : circ a (-1) = -1 := by unfold circ; grind

/-- The shift is a ring isomorphism (R, fadd, ◯) ≅ (R, +, ·): it carries both operations. -/
theorem shift_ring_hom (a b : R) :
    1 + fadd a b = (1 + a) + (1 + b) ∧ 1 + circ a b = (1 + a) * (1 + b) ∧
    (1 : R) + -1 = 0 ∧ (1 : R) + 0 = 1 :=
  ⟨fadd_lift a b, circ_lift a b, by grind, by grind⟩

/-- Plain addition and the circle operation do not form a ring: ◯ does not distribute over +. -/
theorem circ_not_distrib_add : circ (1 : Int) (1 + 1) ≠ circ 1 1 + circ 1 1 := by decide

end FusedRing

/-! ## 6. Circle powers: the integers act on (R, ◯) -/

section Powers
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- (a ◯ b)^◯n = a^◯n ◯ b^◯n: each power map is an endomorphism of the circle operation. -/
theorem mpow_circ (a b : R) (n : Nat) : mpow (circ a b) n = circ (mpow a n) (mpow b n) := by
  induction n with
  | zero => simp only [mpow]; rw [circ_zero]
  | succ n ih => simp only [mpow]; rw [ih]; unfold circ; grind

/-- (a^◯m)^◯n = a^◯(m·n): composing power maps multiplies the exponents. -/
theorem mpow_mul (a : R) (m n : Nat) : mpow (mpow a m) n = mpow a (m * n) := by
  induction n with
  | zero => simp only [mpow, Nat.mul_zero]
  | succ n ih => rw [Nat.mul_succ, mpow_add a (m * n) m, ← ih]; rfl

theorem mpow_one (a : R) : mpow a 1 = a := by simp only [mpow]; rw [zero_circ]

theorem mpow_zero_base (n : Nat) : mpow (0 : R) n = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [mpow]; rw [ih, circ_zero]

end Powers

/-! ## 7. The Tate gluing is a congruence -/

section Tate
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- a ~ b when their lifts differ by powers of q: (1 + a)·qᵐ = (1 + b)·qⁿ. -/
def tateRel (q : R) (a b : R) : Prop := ∃ m n : Nat, (1 + a) * q ^ m = (1 + b) * q ^ n

theorem tateRel_refl (q a : R) : tateRel q a a := ⟨0, 0, rfl⟩

theorem tateRel_symm (q : R) {a b : R} (h : tateRel q a b) : tateRel q b a := by
  obtain ⟨m, n, h⟩ := h; exact ⟨n, m, h.symm⟩

theorem tateRel_trans (q : R) {a b c : R} (h₁ : tateRel q a b) (h₂ : tateRel q b c) : tateRel q a c := by
  obtain ⟨m, n, h₁⟩ := h₁
  obtain ⟨r, s, h₂⟩ := h₂
  refine ⟨m + r, s + n, ?_⟩
  rw [pow_add', pow_add']
  calc (1 + a) * (q ^ m * q ^ r) = ((1 + a) * q ^ m) * q ^ r := by grind
    _ = ((1 + b) * q ^ r) * q ^ n := by rw [h₁]; grind
    _ = (1 + c) * (q ^ s * q ^ n) := by rw [h₂]; grind

/-- The glued point: q − 1 sits with 0, i.e. multiplying the lift by q changes nothing.
    (Stated as 1 + 0 multiplied by q equals 1 + a whenever 1 + a = q.) -/
theorem tate_glue (q a : R) (h : 1 + a = q) : tateRel q a 0 :=
  ⟨0, 1, by rw [h, Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.pow_succ, Lean.Grind.Semiring.pow_zero]; rw [Lean.Grind.Semiring.add_zero, Lean.Grind.Semiring.mul_one, Lean.Grind.Semiring.one_mul, Lean.Grind.Semiring.one_mul]⟩

/-- The circle operation respects the gluing on each side. -/
theorem tateRel_circ (q : R) {a a' b b' : R} (ha : tateRel q a a') (hb : tateRel q b b') :
    tateRel q (circ a b) (circ a' b') := by
  obtain ⟨m, n, ha⟩ := ha
  obtain ⟨r, s, hb⟩ := hb
  refine ⟨m + r, n + s, ?_⟩
  rw [circ_lift, circ_lift, pow_add', pow_add']
  calc (1 + a) * (1 + b) * (q ^ m * q ^ r) = ((1 + a) * q ^ m) * ((1 + b) * q ^ r) := by grind
    _ = ((1 + a') * q ^ n) * ((1 + b') * q ^ s) := by rw [ha, hb]
    _ = (1 + a') * (1 + b') * (q ^ n * q ^ s) := by grind

def tateSetoid (q : R) : Setoid R where
  r := tateRel q
  iseqv := ⟨tateRel_refl q, tateRel_symm q, tateRel_trans q⟩

/-- The Tate quotient: the points of the circle operation with 1 + a and q·(1 + a) glued. -/
def TateQ (q : R) : Type := Quotient (tateSetoid q)

/-- The circle operation on the Tate quotient. -/
def TateQ.op (q : R) : TateQ q → TateQ q → TateQ q :=
  Quotient.lift₂ (s₁ := tateSetoid q) (s₂ := tateSetoid q)
    (fun a b => Quotient.mk (tateSetoid q) (circ a b))
    (fun _ _ _ _ ha hb => Quotient.sound (tateRel_circ q ha hb))

def TateQ.zero (q : R) : TateQ q := Quotient.mk (tateSetoid q) 0

/-- The Tate quotient is a commutative monoid under the circle operation. -/
theorem tateQ_commMonoid (q : R) : CommMonoidLaws (TateQ.op q) (TateQ.zero q) where
  comm x y := by
    induction x using Quotient.inductionOn; induction y using Quotient.inductionOn
    exact congrArg (Quotient.mk _) (circ_comm _ _)
  assoc x y z := by
    induction x using Quotient.inductionOn; induction y using Quotient.inductionOn
    induction z using Quotient.inductionOn
    exact congrArg (Quotient.mk _) (circ_assoc _ _ _)
  ident x := by
    induction x using Quotient.inductionOn
    exact congrArg (Quotient.mk _) (circ_zero _)

/-- In the quotient the class of q − 1 is the identity. -/
theorem tateQ_glued (q a : R) (h : 1 + a = q) :
    (Quotient.mk (tateSetoid q) a : TateQ q) = TateQ.zero q :=
  Quotient.sound (tate_glue q a h)

end Tate

section TateField
variable {K : Type} [Lean.Grind.Field K]

/-- On a field, every class off the wall has an inverse in the Tate quotient. -/
theorem tateQ_inverse (q a : K) (h : 1 + a ≠ 0) :
    TateQ.op q (Quotient.mk (tateSetoid q) a) (Quotient.mk (tateSetoid q) (circInv a)) = TateQ.zero q :=
  congrArg (Quotient.mk _) (circ_circInv a h)

end TateField

/-! ## 8. Landscape countries change with the prime -/

/-- Walls of slope addition (a + b)/(1 − ab) are the roots of 1 + a². Mod 5 and 13 there are
    two (two-walled country, a split torus with p − 1 points); mod 3, 7, 11 none (turning
    country, a non-split torus with p + 1 points). Split exactly at p ≡ 1 mod 4. -/
theorem slope_walls_mod5 : ∀ a : Fin 5, (1 + a * a = 0 ↔ (a = 2 ∨ a = 3)) := by decide
theorem slope_walls_mod13 : ∀ a : Fin 13, (1 + a * a = 0 ↔ (a = 5 ∨ a = 8)) := by decide
theorem slope_no_walls_mod3 : ∀ a : Fin 3, 1 + a * a ≠ 0 := by decide
theorem slope_no_walls_mod7 : ∀ a : Fin 7, 1 + a * a ≠ 0 := by decide
theorem slope_no_walls_mod11 : ∀ a : Fin 11, 1 + a * a ≠ 0 := by decide

/-- The circle operation mod p always has exactly one finite wall, −1 (the other wall is at infinity). -/
theorem circ_one_wall_mod7 : ∀ a : Fin 7, (1 + a = 0 ↔ a = 6) := by decide

end AbstractAlgebra
