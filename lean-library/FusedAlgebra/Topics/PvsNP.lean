/-
  Fused Hyperoperation Algebra: where mult-add touches P vs NP
  ===========================================================
  Core Lean 4 only (no mathlib), no `sorry`. Written for Lean 4.37.

  LIBRARY NOTE: written here Nov 25 and compiled on the first try (20 theorems). Near
  neighbours, not twins: StatFused.de_morgan is the same law over Rat with x ↦ 1 − x;
  gf2_lift_de_morgan is the GF(2) lift itself. QDial.qnum_at_zero/at_one/at_neg_one give
  the dial values of q-integers; np_point/sharpP_point/parityP_point give them for a list
  of witnesses.

  Nothing here proves or refutes P ≠ NP. These are the elementary facts behind the
  honest map in math-theory/p-vs-np/README.md:

  1. The shifted ruler changes nothing about computation. With shifted add
     a ⊕ b = a + b + 1 and mult-add a ⊞ b = a + b + ab, the map x ↦ 1 + x carries
     (⊕, ⊞) gate for gate onto (+, ×). Any arithmetic circuit in one world becomes a
     circuit of exactly the same size in the other (`evalS_lift`, `eval_from_shift`,
     `size_mapC`). So no circuit lower bound, and no fast algorithm, can appear or
     disappear by switching to mult-add.

  2. The dial, folded over a list of yes/no witnesses, gives the q-integer of the number
     of yes answers (`witness_fold`). Its settings are the classic counting classes:
       k = −1 (q = 0):  "is there a witness?"        the NP question   (`np_point`)
       k =  0 (q = 1):  "how many witnesses?"         the #P question   (`sharpP_point`)
       k = −2 (q = −1): "odd or even many?"           the ⊕P question   (`parityP_point`)
       q with qᵖ = 1:   the count mod p               the ModₚP family  (`modp_period`)
     q = 2 keeps the whole count (`qnum_two_strictMono`); q = 0 throws all but one bit
     away (`np_point_forgets`). These classes are standard; the dial is one picture of them.

  3. Over the two-element field, mult-add IS the OR gate and the shift x ↦ 1 + x IS
     NOT, so the lift is De Morgan's law (`gf2_multAdd_is_or`, `gf2_lift_de_morgan`).
     That is exactly the arithmetized OR used in IP = PSPACE and in the
     Razborov–Smolensky lower bounds, which is why the algebrization and natural-proofs
     barriers cover any argument built from it.

  4. The n-input mult-add fold has 2ⁿ − 1 monomials but is built from n − 1 gates
     (`foldMA_ones`): it is easy because it factors. The hard counting problems
     (permanent, #SAT) are the ones whose monomials do not factor.

  5. Min-plus (tropical) laws on Nat, the algebra of dynamic programming
     (`minplus_distrib`).
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.QAcross.QDial

namespace PvsNP
open MultAdd
open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Circuits on the shifted ruler -/

/-- An arithmetic circuit (formula): inputs, constants, and two kinds of gate. -/
inductive Circ (R : Type) where
  | var : Nat → Circ R
  | const : R → Circ R
  | add : Circ R → Circ R → Circ R
  | mul : Circ R → Circ R → Circ R

namespace Circ
variable {R : Type}

/-- Size: number of leaves and gates. -/
def size : Circ R → Nat
  | var _ => 1
  | const _ => 1
  | add a b => size a + size b + 1
  | mul a b => size a + size b + 1

/-- Change every constant by f; the wiring stays the same. -/
def mapC (f : R → R) : Circ R → Circ R
  | var i => var i
  | const c => const (f c)
  | add a b => add (mapC f a) (mapC f b)
  | mul a b => mul (mapC f a) (mapC f b)

theorem size_mapC (f : R → R) (C : Circ R) : size (mapC f C) = size C := by
  induction C with
  | var _ => rfl
  | const _ => rfl
  | add a b iha ihb => simp only [mapC, size, iha, ihb]
  | mul a b iha ihb => simp only [mapC, size, iha, ihb]

end Circ

section Semiring
variable {R : Type} [CommSemiring R]
open Circ

/-- Shifted add: a ⊕ b = a + b + 1 (rung 1 seen on the shifted ruler). -/
def sAdd (a b : R) : R := a + b + 1

/-- The shift carries ⊕ to +. -/
theorem sAdd_lift (a b : R) : 1 + sAdd a b = (1 + a) + (1 + b) := by
  unfold sAdd; grind

/-- On the shifted ruler, ⊞ distributes over ⊕ (the semiring law, moved by x ↦ 1 + x). -/
theorem multAdd_sAdd_distrib (a b c : R) :
    multAdd a (sAdd b c) = sAdd (multAdd a b) (multAdd a c) := by
  unfold multAdd sAdd; grind

/-- Ordinary evaluation with + and ×. -/
def eval (env : Nat → R) : Circ R → R
  | var i => env i
  | const c => c
  | add a b => eval env a + eval env b
  | mul a b => eval env a * eval env b

/-- Shifted evaluation: the same wiring, with ⊕ and ⊞ as the gates. -/
def evalS (env : Nat → R) : Circ R → R
  | var i => env i
  | const c => c
  | add a b => sAdd (evalS env a) (evalS env b)
  | mul a b => multAdd (evalS env a) (evalS env b)

/-- Transport: a ⊕/⊞ circuit is a +/× circuit of the same wiring, seen through x ↦ 1 + x. -/
theorem evalS_lift (env : Nat → R) (C : Circ R) :
    1 + evalS env C = eval (fun i => 1 + env i) (mapC (fun c => 1 + c) C) := by
  induction C with
  | var _ => rfl
  | const _ => rfl
  | add a b iha ihb =>
    simp only [evalS, eval, mapC]; rw [sAdd_lift, iha, ihb]
  | mul a b iha ihb =>
    simp only [evalS, eval, mapC]; rw [multAdd_lift, iha, ihb]

end Semiring

section Ring
variable {R : Type} [CommRing R]
open Circ

/-- And back: every +/× circuit is a ⊕/⊞ circuit of the same wiring and size. -/
theorem eval_from_shift (env : Nat → R) (C : Circ R) :
    eval env C = 1 + evalS (fun i => env i - 1) (mapC (fun c => c - 1) C) := by
  induction C with
  | var _ => simp only [eval, evalS, mapC]; grind
  | const _ => simp only [eval, evalS, mapC]; grind
  | add a b iha ihb =>
    simp only [eval, evalS, mapC]; rw [iha, ihb]; unfold sAdd; grind
  | mul a b iha ihb =>
    simp only [eval, evalS, mapC]; rw [iha, ihb]; unfold multAdd; grind

/-- The identities of the shifted semiring: −1 for ⊕ and 0 for ⊞; −1 absorbs under ⊞. -/
theorem shifted_identities (a : R) :
    sAdd a (-1) = a ∧ multAdd a 0 = a ∧ multAdd a (-1) = -1 := by
  unfold sAdd multAdd; refine ⟨?_, ?_, ?_⟩ <;> grind

end Ring

/-! ## 2. The dial over witnesses: NP, #P, ⊕P, ModₚP -/

section Witness
variable {R : Type} [CommRing R]

/-- A yes/no answer as a number. -/
def ind (b : Bool) : R := if b then 1 else 0

/-- Fold a list with the weighted mult-add a + b + k·ab. -/
def foldW (k : R) : List R → R
  | [] => 0
  | a :: l => wAdd k a (foldW k l)

/-- Number of yes answers. -/
def countT : List Bool → Nat
  | [] => 0
  | b :: l => (if b then 1 else 0) + countT l

/-- The dial folded over witnesses is the q-integer of the witness count, q = k + 1. -/
theorem witness_fold (k : R) (bs : List Bool) :
    foldW k (bs.map ind) = QDial.qnum (k + 1) (countT bs) := by
  induction bs with
  | nil => rfl
  | cons b l ih =>
    cases b
    · simp only [List.map, foldW, countT, ind, ih]
      unfold wAdd; simp only [ite_false, Bool.false_eq_true]; grind
    · simp only [List.map, foldW, countT, ind, ih, ite_true]
      rw [Nat.add_comm, QDial.qnum_succ]; unfold wAdd; grind

/-- The ModₚP family: if qᵖ = 1, the count only matters mod p (through the lift). -/
theorem modp_period (q : R) (p m : Nat) (hq : q ^ p = 1) :
    (q - 1) * QDial.qnum q (m + p) = (q - 1) * QDial.qnum q m := by
  have h1 := QDial.qnum_lift q (m + p)
  have h2 := QDial.qnum_lift q m
  rw [pow_add', hq] at h1
  grind

end Witness

/-- "Is there a yes?" = the NP question, from the dial at k = −1 (q = 0). -/
def anyT : List Bool → Bool
  | [] => false
  | b :: l => b || anyT l

/-- Parity of the yes answers = the ⊕P question. -/
def parT : List Bool → Bool
  | [] => false
  | b :: l => xor b (parT l)

theorem np_point (bs : List Bool) :
    foldW (-1 : Int) (bs.map ind) = ind (anyT bs) := by
  induction bs with
  | nil => rfl
  | cons b l ih =>
    cases b <;> simp only [List.map, foldW, anyT, ih, Bool.false_or, Bool.true_or] <;>
      (unfold wAdd ind; split <;> simp_all <;> omega)

theorem sharpP_point (bs : List Bool) :
    foldW (0 : Int) (bs.map ind) = (countT bs : Int) := by
  induction bs with
  | nil => rfl
  | cons b l ih =>
    cases b <;> simp only [List.map, foldW, countT, ih] <;> (unfold wAdd ind; simp; try omega)

theorem parityP_point (bs : List Bool) :
    foldW (-2 : Int) (bs.map ind) = ind (parT bs) := by
  induction bs with
  | nil => rfl
  | cons b l ih =>
    cases b <;> simp only [List.map, foldW, parT, ih, Bool.false_xor, Bool.true_xor] <;>
      (unfold wAdd ind; cases parT l <;> simp)

/-- q = 0 forgets: one witness and five witnesses give the same answer. -/
theorem np_point_forgets : QDial.qnum (0 : Int) 1 = QDial.qnum (0 : Int) 5 := by decide

theorem qnum_two_nonneg (n : Nat) : 0 ≤ QDial.qnum (2 : Int) n := by
  induction n with
  | zero => decide
  | succ n ih => rw [QDial.qnum_succ]; omega

/-- q = 2 keeps the whole count: different counts give different answers. -/
theorem qnum_two_strictMono {m n : Nat} (h : m < n) :
    QDial.qnum (2 : Int) m < QDial.qnum (2 : Int) n := by
  induction n with
  | zero => omega
  | succ n ih =>
    rw [QDial.qnum_succ]
    have := qnum_two_nonneg n
    rcases Nat.lt_succ_iff_lt_or_eq.mp h with h' | h'
    · have := ih h'; omega
    · subst h'; omega

/-! ## 3. The two-element field: mult-add is OR, the shift is NOT -/

/-- Mult-add over GF(2): + is xor, × is and. -/
def gf2MultAdd (a b : Bool) : Bool := xor (xor a b) (a && b)

theorem gf2_multAdd_is_or (a b : Bool) : gf2MultAdd a b = (a || b) := by
  cases a <;> cases b <;> rfl

theorem gf2_shift_is_not (a : Bool) : xor true a = !a := by cases a <;> rfl

/-- The lift 1 + (a ⊞ b) = (1 + a)(1 + b) over GF(2) is De Morgan: ¬(a ∨ b) = ¬a ∧ ¬b. -/
theorem gf2_lift_de_morgan (a b : Bool) :
    xor true (gf2MultAdd a b) = (xor true a && xor true b) := by
  cases a <;> cases b <;> rfl

/-! ## 4. Easy because it factors -/

/-- n gates, 2ⁿ − 1 monomials: 1 + (1 ⊞ 1 ⊞ … ⊞ 1) = 2ⁿ. -/
theorem foldMA_ones (n : Nat) : 1 + foldMA (List.replicate n (1 : Nat)) = 2 ^ n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [List.replicate_succ, foldMA, Nat.pow_succ, ← ih]; unfold multAdd; omega

/-! ## 5. Min-plus: the algebra of dynamic programming -/

theorem minplus_distrib (a b c : Nat) : a + min b c = min (a + b) (a + c) := by omega

theorem minplus_assoc (a b c : Nat) : min (min a b) c = min a (min b c) := by omega

end PvsNP
