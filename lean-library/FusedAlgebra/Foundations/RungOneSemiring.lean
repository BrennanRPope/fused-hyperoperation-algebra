/-
  Foundations.RungOneSemiring: the rung-1 semiring, stated once and proved on its own
  ===================================================================================

  Core Lean 4 only (no mathlib, no `sorry`). This file imports nothing from the
  library: every law below is proved here from the definitions, so the unified
  result stands on its own. (Nov 25, 2026; asked for by Brennan after the semiring
  audit, math-theory/semiring-audit/.)

  The pair. On any commutative semiring R:
    "or"   sadd a b = a + b + 1          (rung 1, the shifted add)
    "and"  circ a b = a + b + a·b        (the circle operation ∘₁, mult-add)
  and for a dial k ≠ 0 in a field:
    "or"   pseudoAdd k a b = a + b + 1/k
    "and"  circK k a b     = a + b + k·a·b   (∘ₖ)

  What is proved:
    1. `circ` distributes over `sadd` on both sides, on any commutative semiring
       (so on Nat, Int, Rat at once): `circ_sadd_distrib`, `sadd_circ_distrib`.
    2. On a commutative ring the pair satisfies every commutative-ring law, with
       zero −1 (the wall, which absorbs under `circ`) and one 0: `isCommRing`.
       On a commutative semiring without minus (Nat) every law but the zero holds:
       `isCommSemiringNoZero`.
    3. The shift x ↦ 1 + x carries the pair to (+, ×), −1 to 0 and 0 to 1, and on a
       ring it has the inverse x ↦ x − 1: `shift_hom`, `shift_bij`. This is why all
       of 1 and 2 hold: the pair is ordinary arithmetic moved by a translation.
    4. For every dial k ≠ 0 in a field, (pseudoAdd k, circK k) is again a commutative
       ring, zero −1/k, one 0, carried to (+, ×) by g(x) = 1 + kx: `isCommRingK`,
       `gen_hom`. At k = 1 it is the pair above: `pseudoAdd_one`, `circK_one`.

  Standard names (see semiring-audit/README.md for sources):
    * the ring (R, a+b+1, a+b+ab) is a transport of structure along x ↦ x − 1, a
      textbook exercise; Mathlib builds it with `Equiv.commRing`;
    * `circ` is Jacobson's circle composition (Mathlib: `PreQuasiregular`);
    * (pseudoAdd k, circK k) is a generated semiring with generator g(x) = 1 + kx,
      in E. Pap's pseudo-analysis (g-calculus);
    * read on the height ruler a = eˢ − 1, `sadd` is log-sum-exp and `circ` is +:
      the log semiring, whose cold limit is the tropical semiring (Maslov). That
      reading needs exp, so it is checked in semiring-audit/checks.py, not here.

  Lineage: before this file the same pair was defined five times and its
  distributive law proved four times. FusedAlgebra.Foundations.RungOneLineage
  ties each of those to this file.
-/

namespace RungOne

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## The law records -/

/-- Commutative-semiring laws for an "or" `add` and an "and" `mul`, one `one`,
    and no zero required (Nat has no zero for `sadd`). -/
structure CommSemiringLawsNoZero {α : Type} (add mul : α → α → α) (one : α) : Prop where
  add_comm : ∀ a b, add a b = add b a
  add_assoc : ∀ a b c, add (add a b) c = add a (add b c)
  mul_comm : ∀ a b, mul a b = mul b a
  mul_assoc : ∀ a b c, mul (mul a b) c = mul a (mul b c)
  mul_one : ∀ a, mul a one = a
  left_distrib : ∀ a b c, mul a (add b c) = add (mul a b) (mul a c)
  right_distrib : ∀ a b c, mul (add a b) c = add (mul a c) (mul b c)

/-- Commutative-ring laws: the above, plus a zero that is the identity of `add`,
    absorbs under `mul`, and additive inverses. -/
structure CommRingLaws {α : Type} (add mul : α → α → α) (zero one : α) : Prop
    extends CommSemiringLawsNoZero add mul one where
  add_zero : ∀ a, add a zero = a
  mul_zero : ∀ a, mul a zero = zero
  add_inv : ∀ a, ∃ b, add a b = zero

/-! ## 1. The pair on any commutative semiring -/

section Semiring
variable {R : Type} [CommSemiring R]

/-- Rung 1, the shifted add: a + b + 1. The "or" of the rung-1 semiring. -/
def sadd (a b : R) : R := a + b + 1

/-- The circle operation ∘₁ (mult-add): a + b + a·b. The "and" of the rung-1 semiring. -/
def circ (a b : R) : R := a + b + a * b

theorem sadd_comm (a b : R) : sadd a b = sadd b a := by unfold sadd; grind
theorem sadd_assoc (a b c : R) : sadd (sadd a b) c = sadd a (sadd b c) := by unfold sadd; grind
theorem circ_comm (a b : R) : circ a b = circ b a := by unfold circ; grind
theorem circ_assoc (a b c : R) : circ (circ a b) c = circ a (circ b c) := by unfold circ; grind
theorem circ_zero (a : R) : circ a 0 = a := by unfold circ; grind

/-- **The unified distributive law.** a ∘₁ (b ⊕ c) = (a ∘₁ b) ⊕ (a ∘₁ c),
    with ⊕ the shifted add. Proved directly from the definitions. -/
theorem circ_sadd_distrib (a b c : R) : circ a (sadd b c) = sadd (circ a b) (circ a c) := by
  unfold circ sadd; grind

/-- The same law from the right. -/
theorem sadd_circ_distrib (a b c : R) : circ (sadd a b) c = sadd (circ a c) (circ b c) := by
  unfold circ sadd; grind

/-- The shift x ↦ 1 + x carries `sadd` to +. -/
theorem shift_sadd (a b : R) : 1 + sadd a b = (1 + a) + (1 + b) := by unfold sadd; grind

/-- The shift x ↦ 1 + x carries `circ` to ×. -/
theorem shift_circ (a b : R) : 1 + circ a b = (1 + a) * (1 + b) := by unfold circ; grind

/-- Every commutative semiring gives the rung-1 laws (no zero needed). On Nat this is
    the positive integers' (+, ×) moved down by one. -/
theorem isCommSemiringNoZero : CommSemiringLawsNoZero (sadd (R := R)) circ 0 :=
  { add_comm := sadd_comm, add_assoc := sadd_assoc, mul_comm := circ_comm,
    mul_assoc := circ_assoc, mul_one := circ_zero,
    left_distrib := circ_sadd_distrib, right_distrib := sadd_circ_distrib }

end Semiring

/-! ## 2. On a commutative ring: zero −1, the wall -/

section Ring
variable {R : Type} [CommRing R]

/-- −1 is the identity of the shifted add. -/
theorem sadd_neg_one (a : R) : sadd a (-1) = a := by unfold sadd; grind

/-- −1 (the wall) absorbs under ∘₁. -/
theorem circ_neg_one (a : R) : circ a (-1) = -1 := by unfold circ; grind

/-- The shifted-add inverse of a is −a − 2. -/
theorem sadd_inv (a : R) : sadd a (-a - 2) = -1 := by unfold sadd; grind

/-- **The rung-1 pair is a commutative ring**, zero −1, one 0. -/
theorem isCommRing : CommRingLaws (sadd (R := R)) circ (-1) 0 :=
  { toCommSemiringLawsNoZero := isCommSemiringNoZero,
    add_zero := sadd_neg_one, mul_zero := circ_neg_one,
    add_inv := fun a => ⟨-a - 2, sadd_inv a⟩ }

/-- The shift is a homomorphism of the whole structure: both operations, zero and one. -/
theorem shift_hom (a b : R) :
    1 + sadd a b = (1 + a) + (1 + b) ∧ 1 + circ a b = (1 + a) * (1 + b) ∧
    (1 : R) + -1 = 0 ∧ (1 : R) + 0 = 1 :=
  ⟨shift_sadd a b, shift_circ a b, by grind, by grind⟩

/-- … and it is a bijection, with inverse x ↦ x − 1. -/
theorem shift_bij (x y : R) : (1 + x) - 1 = x ∧ 1 + (y - 1) = y := ⟨by grind, by grind⟩

end Ring

/-! ## 3. Every dial k ≠ 0: a generated semiring with g(x) = 1 + kx -/

section Field
variable {K : Type} [Field K]

/-- The "or" for dial k: a + b + 1/k (Pap's pseudo-addition for g(x) = 1 + kx). -/
def pseudoAdd (k a b : K) : K := a + b + 1 / k

/-- The "and" for dial k: a + b + k·a·b (∘ₖ). -/
def circK (k a b : K) : K := a + b + k * a * b

theorem pseudoAdd_one (a b : K) : pseudoAdd 1 a b = sadd a b := by unfold pseudoAdd sadd; grind
theorem circK_one (a b : K) : circK 1 a b = circ a b := by unfold circK circ; grind

/-- The generator g(x) = 1 + kx carries the pair to (+, ×). -/
theorem gen_hom (k a b : K) (hk : k ≠ 0) :
    1 + k * pseudoAdd k a b = (1 + k * a) + (1 + k * b) ∧
    1 + k * circK k a b = (1 + k * a) * (1 + k * b) := by
  unfold pseudoAdd circK
  have h : k * (1 / k) = 1 := by
    rw [Field.div_eq_mul_inv, ← Semiring.mul_assoc, Semiring.mul_one, Field.mul_inv_cancel hk]
  constructor <;> grind

/-- The dial-k distributive law, for every k ≠ 0. -/
theorem circK_pseudoAdd_distrib (k a b c : K) (hk : k ≠ 0) :
    circK k a (pseudoAdd k b c) = pseudoAdd k (circK k a b) (circK k a c) := by
  unfold circK pseudoAdd
  have h : k * (1 / k) = 1 := by
    rw [Field.div_eq_mul_inv, ← Semiring.mul_assoc, Semiring.mul_one, Field.mul_inv_cancel hk]
  grind

/-- **For every k ≠ 0, (pseudoAdd k, circK k) is a commutative ring**, zero −1/k, one 0. -/
theorem isCommRingK (k : K) (hk : k ≠ 0) : CommRingLaws (pseudoAdd k) (circK k) (-1 / k) 0 := by
  have h : k * (1 / k) = 1 := by
    rw [Field.div_eq_mul_inv, ← Semiring.mul_assoc, Semiring.mul_one, Field.mul_inv_cancel hk]
  have hn : -1 / k = -(1 / k) := by
    rw [Field.div_eq_mul_inv, Field.div_eq_mul_inv]; grind
  exact
  { add_comm := fun a b => by unfold pseudoAdd; grind
    add_assoc := fun a b c => by unfold pseudoAdd; grind
    mul_comm := fun a b => by unfold circK; grind
    mul_assoc := fun a b c => by unfold circK; grind
    mul_one := fun a => by unfold circK; grind
    left_distrib := fun a b c => circK_pseudoAdd_distrib k a b c hk
    right_distrib := fun a b c => by unfold circK pseudoAdd; grind
    add_zero := fun a => by unfold pseudoAdd; rw [hn]; grind
    mul_zero := fun a => by unfold circK; rw [hn]; grind
    add_inv := fun a => ⟨-a - 1 / k - 1 / k, by unfold pseudoAdd; rw [hn]; grind⟩ }

end Field

end RungOne
