/-
  Fused Hyperoperation Algebra: the q dial (one q for the gluing, the counting and the dial)
  ==========================================================================================
  Core Lean 4 only (no mathlib), no `sorry`.

  LIBRARY NOTE: copied from q-across/lean/QDial.lean (Nov 25) and compiled here for the
  first time: no sorry, standard axioms only; no changes were needed. No statement changed.

  The standard facts used: q-integers [n]_q = 1 + q + … + q^(n-1) and the rule
  [m+n]_q = [m]_q + q^m·[n]_q (standard since Gauss and Heine). What this file shows is that
  they are the weighted mult-add with weight k = q − 1:

      [n]_q          = n copies of 1 under ⊞ₖ           (qnum_eq_wpow)
      [m+n]_q        = [m]_q ⊞ₖ [n]_q                    (qnum_add)
      (q − 1)·[n]_q  = (q − 1)^(⊞n), the mult-add power  (qnum_mpow)
      (1+a)·qⁿ       = 1 + (a ⊞ (q−1)·[n]_q)             (tate_glue_n)

  The last line ties this to the Tate gluing (1+a) ~ q·(1+a) of the elliptic work: going n
  times round the torus's q-direction is a mult-add step of exactly (q − 1)·[n]_q.

  The dial values: q = 1 (k = 0) counts, q = 2 (k = 1) is mult-add, q = 0 (k = −1) is the
  failure union, q = −1 (k = −2) is parity.
-/
import FusedAlgebra.Base

namespace QDial
open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [CommRing R]

/-- The q-integer [n]_q = 1 + q + q² + … + q^(n-1), built as [n+1] = 1 + q·[n]. -/
def qnum (q : R) : Nat → R
  | 0 => 0
  | n + 1 => 1 + q * qnum q n

theorem qnum_zero (q : R) : qnum q 0 = 0 := rfl
theorem qnum_succ (q : R) (n : Nat) : qnum q (n + 1) = 1 + q * qnum q n := rfl

/-- The q-integer is the weighted mult-add power of 1 with weight k = q − 1. -/
theorem qnum_eq_wpow (q : R) (n : Nat) : qnum q n = wpow (q - 1) 1 n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [qnum_succ, ih]
    show 1 + q * wpow (q - 1) 1 n = circK (q - 1) (wpow (q - 1) 1 n) 1
    unfold circK; grind

/-- The lift: 1 + (q − 1)·[n]_q = qⁿ. (So q-counting is mult-add counting through x ↦ 1 + x.) -/
theorem qnum_lift (q : R) (n : Nat) : 1 + (q - 1) * qnum q n = q ^ n := by
  rw [qnum_eq_wpow, wpow_lift]
  congr 1; grind

/-- q-integers add by the weighted mult-add, weight q − 1: [m+n] = [m] ⊞ₖ [n]. -/
theorem qnum_add (q : R) (m n : Nat) :
    qnum q (m + n) = circK (q - 1) (qnum q m) (qnum q n) := by
  induction n with
  | zero =>
    show qnum q m = circK (q - 1) (qnum q m) 0
    unfold circK; grind
  | succ n ih =>
    rw [← Nat.add_assoc, qnum_succ, ih, qnum_succ]
    unfold circK; grind

/-- The standard form of the same rule: [m+n] = [m] + qᵐ·[n]. -/
theorem qnum_add_std (q : R) (m n : Nat) :
    qnum q (m + n) = qnum q m + q ^ m * qnum q n := by
  rw [qnum_add, ← qnum_lift q m]
  unfold circK; grind

/-- (q − 1)·[n]_q is the mult-add power (q − 1)^(⊞n). -/
theorem qnum_mpow (q : R) (n : Nat) : (q - 1) * qnum q n = mpow (q - 1) n := by
  have h1 := qnum_lift q n
  have h2 := mpow_lift (q - 1) n
  have h3 : (1 : R) + (q - 1) = q := by grind
  rw [h3] at h2
  grind

/-- One Tate step: multiplying the lift 1 + a by q is a mult-add step by q − 1. -/
theorem tate_glue_step (q a : R) : (1 + a) * q = 1 + circ a (q - 1) := by
  unfold circ; grind

/-- n Tate steps: (1 + a)·qⁿ = 1 + (a ⊞ (q − 1)·[n]_q). Going n times round the torus's
    q-direction is one mult-add step whose length is (q − 1) times the q-integer [n]_q. -/
theorem tate_glue_n (q a : R) (n : Nat) :
    (1 + a) * q ^ n = 1 + circ a ((q - 1) * qnum q n) := by
  rw [← qnum_lift q n]
  unfold circ; grind

end Ring

/-! ## The dial values, on Int (proved for every n) -/

/-- q = 1, k = 0: the q-integer is plain counting. -/
theorem qnum_at_one (n : Nat) : qnum (1 : Int) n = n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [qnum_succ, ih]; omega

/-- q = 2, k = 1: [n]₂ = 2ⁿ − 1, the mult-add power of 1 (and the number of nonzero vectors
    in an n-dimensional space over the two-element field). -/
theorem qnum_at_two (n : Nat) : qnum (2 : Int) n = mpow 1 n := by
  rw [qnum_eq_wpow]
  have h : (2 : Int) - 1 = 1 := by decide
  rw [h]
  exact wpow_one 1 n

/-- q = 0, k = −1: every nonempty count collapses to 1 (the failure union 1 ⊞₋₁ 1 = 1). -/
theorem qnum_at_zero (n : Nat) : qnum (0 : Int) (n + 1) = 1 := by
  rw [qnum_succ]; grind

/-- q = −1, k = −2: parity. -/
theorem qnum_at_neg_one (n : Nat) :
    qnum (-1 : Int) (2 * n) = 0 ∧ qnum (-1 : Int) (2 * n + 1) = 1 := by
  induction n with
  | zero => exact ⟨rfl, by decide⟩
  | succ n ih =>
    obtain ⟨h0, h1⟩ := ih
    have e : 2 * (n + 1) = (2 * n + 1) + 1 := by omega
    constructor
    · rw [e, qnum_succ, h1]; decide
    · rw [e, qnum_succ, qnum_succ, h1]; decide

/-- The parity law is ⊞ at k = −2: on {0, 1} it is "exclusive or".
    Same truth table as StatFused.xor_table, which is over Rat; this one is over Int. -/
theorem parity_xor :
    circK (-2 : Int) 0 0 = 0 ∧ circK (-2 : Int) 0 1 = 1 ∧
    circK (-2 : Int) 1 0 = 1 ∧ circK (-2 : Int) 1 1 = 0 := by
  unfold circK; decide

/-- The failure union at k = −1 on {0, 1} is "or".
    Same truth table as StatFused.or_table, which is over Rat; this one is over Int. -/
theorem failure_or :
    circK (-1 : Int) 0 0 = 0 ∧ circK (-1 : Int) 0 1 = 1 ∧
    circK (-1 : Int) 1 0 = 1 ∧ circK (-1 : Int) 1 1 = 1 := by
  unfold circK; decide

/-- Small values: [3]_q at the four dial points (counting 3, volume 7, union 1, parity 1). -/
theorem qnum_three_at_dial :
    qnum (1 : Int) 3 = 3 ∧ qnum (2 : Int) 3 = 7 ∧ qnum (0 : Int) 3 = 1 ∧ qnum (-1 : Int) 3 = 1 := by
  decide

end QDial
