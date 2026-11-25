/-
  FusedAlgebra.Base: the shared lemmas of Fused Hyperoperation Algebra
  ====================================================================

  Every topic module imports this file, so each law is proved once.
  Core Lean 4 only (no mathlib, no `sorry`). Written for Lean 4.37.

  Mult-add:            a ◯ b = a + b + a·b = (1 + a)(1 + b) − 1
  Weighted family:     a ◯[k] b = a + b + k·a·b,  1 + k(a ◯[k] b) = (1 + ka)(1 + kb)
  Mult-add power:      a^(◯n) = (1 + a)^n − 1
  Way back (inverse):  a ◯ (−a/(1 + a)) = 0, everywhere but the wall a = −1

  The laws are stated once for any commutative semiring (so they hold on
  Nat, Int and Rat at the same time); the ones that need minus are stated
  for rings, and the inverse for fields. `grind` does the ring arithmetic.
-/

namespace CircleOp

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Mult-add on any commutative semiring -/

section Semiring
variable {R : Type} [CommSemiring R]

/-- The circle operation: a ◯ b = a + b + a·b (ring theory's circle operation, plus convention; called mult-add before 1.0). -/
def circ (a b : R) : R := a + b + a * b

@[inherit_doc] scoped infixl:65 " ◯ " => circ

theorem circ_def (a b : R) : circ a b = a + b + a * b := rfl

/-- The lift: mult-add is multiplication seen through x ↦ 1 + x. -/
theorem circ_lift (a b : R) : 1 + circ a b = (1 + a) * (1 + b) := by
  unfold circ; grind

/-- The lift with the 1 on the right (handy on Nat, where we add 1 to both sides). -/
theorem circ_lift' (a b : R) : circ a b + 1 = (a + 1) * (b + 1) := by
  unfold circ; grind

theorem circ_comm (a b : R) : circ a b = circ b a := by
  unfold circ; grind

theorem circ_assoc (a b c : R) :
    circ (circ a b) c = circ a (circ b c) := by
  unfold circ; grind

/-- 0 is the identity of mult-add (because 1 is the identity of ×). -/
theorem circ_zero (a : R) : circ a 0 = a := by
  unfold circ; grind

theorem zero_circ (a : R) : circ 0 a = a := by
  unfold circ; grind

/-- Mult-add is never less than its two pieces added: the corner a·b is extra. -/
theorem circ_excess (a b : R) : circ a b = (a + b) + a * b := rfl

/-! ### The weighted family a ◯[k] b = a + b + k·a·b -/

/-- The circle operation with dial k: a ◯[k] b = a + b + k·a·b. -/
def circK (k a b : R) : R := a + b + k * a * b

@[inherit_doc] scoped notation:65 a:65 " ◯[" k "] " b:66 => circK k a b

theorem circK_def (k a b : R) : circK k a b = a + b + k * a * b := rfl

/-- The weighted lift: 1 + k·(a ◯[k] b) = (1 + k·a)(1 + k·b). -/
theorem circK_lift (k a b : R) : 1 + k * circK k a b = (1 + k * a) * (1 + k * b) := by
  unfold circK; grind

/-- Weighted mult-add is plain mult-add seen through x ↦ k·x. -/
theorem circK_scale (k a b : R) : k * circK k a b = circ (k * a) (k * b) := by
  unfold circK circ; grind

theorem circK_comm (k a b : R) : circK k a b = circK k b a := by
  unfold circK; grind

theorem circK_assoc (k a b c : R) : circK k (circK k a b) c = circK k a (circK k b c) := by
  unfold circK; grind

theorem circK_zero (k a : R) : circK k a 0 = a := by
  unfold circK; grind

/-- k = 0: no corner, plain addition. -/
theorem circK_weight_zero (a b : R) : circK 0 a b = a + b := by
  unfold circK; grind

/-- k = 1: plain mult-add. -/
theorem circK_one (a b : R) : circK 1 a b = circ a b := by
  unfold circK circ; grind

/-! ### Plain powers by repeated multiplication

  `pw x n` is x·x·…·x, built step by step so every power law a topic uses is
  proved here (it equals the built-in `x ^ n`, see `pw_eq_pow`). -/

def pw (x : R) : Nat → R
  | 0 => 1
  | n + 1 => pw x n * x

theorem pw_eq_pow (x : R) (n : Nat) : pw x n = x ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; rfl
  | succ n ih => rw [Semiring.pow_succ, pw, ih]

theorem pw_add (x : R) (m n : Nat) : pw x (m + n) = pw x m * pw x n := by
  induction n with
  | zero => simp only [Nat.add_zero, pw]; grind
  | succ n ih => rw [Nat.add_succ, pw, pw, ih]; grind

theorem pw_mul (x : R) (m n : Nat) : pw x (m * n) = pw (pw x m) n := by
  induction n with
  | zero => simp only [Nat.mul_zero, pw]
  | succ n ih => rw [Nat.mul_succ, pw_add, ih, pw]

theorem pw_one (n : Nat) : pw (1 : R) n = 1 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [pw, ih]; grind

theorem pw_mul_base (x y : R) (n : Nat) : pw (x * y) n = pw x n * pw y n := by
  induction n with
  | zero => simp only [pw]; grind
  | succ n ih => rw [pw, pw, pw, ih]; grind

/-! ### Three power laws core Lean does not state for every number type -/

theorem pow_add' (a : R) (m n : Nat) : a ^ (m + n) = a ^ m * a ^ n := by
  induction n with
  | zero => rw [Nat.add_zero, Semiring.pow_zero]; grind
  | succ n ih => rw [← Nat.add_assoc, Semiring.pow_succ, Semiring.pow_succ, ih]; grind

theorem mul_pow' (a b : R) (n : Nat) : (a * b) ^ n = a ^ n * b ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero, Semiring.pow_zero, Semiring.pow_zero]; grind
  | succ n ih => rw [Semiring.pow_succ, Semiring.pow_succ, Semiring.pow_succ, ih]; grind

theorem one_pow' (n : Nat) : (1 : R) ^ n = 1 := by
  induction n with
  | zero => rw [Semiring.pow_zero]
  | succ n ih => rw [Semiring.pow_succ, ih]; grind

/-! ### Powers: a ◯ a ◯ … ◯ a climbs the cube -/

/-- The mult-add power a^(◯n): n copies of a, mult-added. -/
def mpow (a : R) : Nat → R
  | 0 => 0
  | n + 1 => circ (mpow a n) a

/-- 1 + a^(◯n) = (1 + a)^n. -/
theorem mpow_lift (a : R) (n : Nat) : 1 + mpow a n = (1 + a) ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; simp only [mpow]; grind
  | succ n ih => rw [Semiring.pow_succ, ← ih, mpow, circ_lift]

/-- Powers add: a^(◯(m+n)) = a^(◯m) ◯ a^(◯n). -/
theorem mpow_add (a : R) (m n : Nat) : mpow a (m + n) = circ (mpow a m) (mpow a n) := by
  induction n with
  | zero => simp only [Nat.add_zero, mpow]; rw [circ_zero]
  | succ n ih => rw [← Nat.add_assoc, mpow, mpow, ih, circ_assoc]

/-- The weighted power: n copies of a under ◯[k]. -/
def wpow (k a : R) : Nat → R
  | 0 => 0
  | n + 1 => circK k (wpow k a n) a

/-- 1 + k·a^(◯[k]n) = (1 + k·a)^n. -/
theorem wpow_lift (k a : R) (n : Nat) : 1 + k * wpow k a n = (1 + k * a) ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; simp only [wpow]; grind
  | succ n ih => rw [Semiring.pow_succ, ← ih, wpow, circK_lift]

theorem wpow_one (a : R) (n : Nat) : wpow 1 a n = mpow a n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [wpow, mpow, ih, circK_one]

/-! ### Lists: x₁ ◯ x₂ ◯ … ◯ xₙ -/

/-- The mult-add of a list. -/
def foldCirc : List R → R
  | [] => 0
  | a :: l => circ a (foldCirc l)

/-- Product of a list. -/
def prodL : List R → R
  | [] => 1
  | a :: l => a * prodL l

/-- Sum of a list. -/
def sumL : List R → R
  | [] => 0
  | a :: l => a + sumL l

/-- The whole is the first part ◯ the second part. -/
theorem foldCirc_append (S T : List R) :
    foldCirc (S ++ T) = circ (foldCirc S) (foldCirc T) := by
  induction S with
  | nil => simp only [List.nil_append, foldCirc]; rw [zero_circ]
  | cons a S ih =>
    simp only [List.cons_append, foldCirc]
    rw [ih, circ_assoc]

/-- The list lift: 1 + (x₁ ◯ … ◯ xₙ) = (1 + x₁)⋯(1 + xₙ). -/
theorem foldCirc_lift (l : List R) : 1 + foldCirc l = prodL (l.map (fun x => 1 + x)) := by
  induction l with
  | nil => simp only [foldCirc, List.map, prodL]; grind
  | cons a l ih => simp only [foldCirc, List.map, prodL]; rw [circ_lift, ih]

/-- Order never matters: any reshuffle gives the same total. -/
theorem foldCirc_perm {l l' : List R} (h : List.Perm l l') : foldCirc l = foldCirc l' := by
  induction h with
  | nil => rfl
  | cons x _ ih => simp only [foldCirc]; rw [ih]
  | swap x y l =>
    simp only [foldCirc]
    unfold circ; grind
  | trans _ _ ih₁ ih₂ => exact ih₁.trans ih₂

/-- All sub-lists: the corners of the cube. -/
def groups {α : Type} : List α → List (List α)
  | [] => [[]]
  | a :: l => groups l ++ (groups l).map (fun g => a :: g)

theorem sumL_append (l m : List R) : sumL (l ++ m) = sumL l + sumL m := by
  induction l with
  | nil => simp only [List.nil_append, sumL]; grind
  | cons a l ih => simp only [List.cons_append, sumL]; rw [ih]; grind

theorem sumL_map_cons (a : R) (gs : List (List R)) :
    sumL ((gs.map (fun g => a :: g)).map prodL) = a * sumL (gs.map prodL) := by
  induction gs with
  | nil => simp only [List.map, sumL]; grind
  | cons g gs ih => simp only [List.map, sumL, prodL] at ih ⊢; rw [ih]; grind

/-- The corner expansion: 1 + (x₁ ◯ … ◯ xₙ) is the sum over every corner of
    the cube (every sub-list) of the product of its entries. -/
theorem one_add_foldCirc_eq_corners (l : List R) :
    1 + foldCirc l = sumL ((groups l).map prodL) := by
  induction l with
  | nil => simp only [foldCirc, groups, List.map, prodL, sumL] <;> grind
  | cons a l ih =>
    simp only [foldCirc, groups, List.map_append]
    rw [sumL_append, sumL_map_cons, ← ih]
    unfold circ; grind

end Semiring

/-! ## 2. Laws that need minus (commutative rings) -/

section Ring
variable {R : Type} [CommRing R]

/-- The wall: −1 absorbs everything mult-added to it. -/
theorem neg_one_absorbs (b : R) : circ (-1) b = -1 := by
  unfold circ; grind

/-- The power in closed form: a^(◯n) = (1 + a)^n − 1. -/
theorem mpow_eq (a : R) (n : Nat) : mpow a n = (1 + a) ^ n - 1 := by
  have h := mpow_lift a n
  grind

/-- k = −1 is "at least one of two independent events". -/
theorem circK_either (a b : R) : circK (-1) a b = 1 - (1 - a) * (1 - b) := by
  unfold circK; grind

/-- k = −1 is mult-add mirrored: a ◯[−1] b = −((−a) ◯ (−b)). -/
theorem circK_neg_one_mirror (a b : R) : circK (-1) a b = -(circ (-a) (-b)) := by
  unfold circK circ; grind

/-- Flipping both signs flips the weight. -/
theorem circK_flip (k a b : R) : circK k (-a) (-b) = -(circK (-k) a b) := by
  unfold circK; grind

/-- The weighted power at k = −1: n independent chances, 1 − (1 − a)^n. -/
theorem wpow_neg_one (a : R) (n : Nat) : wpow (-1) a n = 1 - (1 - a) ^ n := by
  have h := wpow_lift (-1) a n
  have e : (1 + -1 * a) = 1 - a := by grind
  rw [e] at h; grind

end Ring

/-! ## 3. The way back needs division (fields) -/

section Field
variable {K : Type} [Field K]

/-- The mult-add inverse: −u/(1 + u). -/
def circInv (u : K) : K := -u / (1 + u)

/-- Everywhere but the wall there is an exact way back. -/
theorem circ_circInv (u : K) (h : 1 + u ≠ 0) : circ u (circInv u) = 0 := by
  unfold circ circInv; grind

theorem circInv_circ (u : K) (h : 1 + u ≠ 0) : circ (circInv u) u = 0 := by
  rw [circ_comm]; exact circ_circInv u h

/-- Under the lift the inverse is the ordinary reciprocal. -/
theorem inv_one_add (u : K) (h : 1 + u ≠ 0) : (1 + u)⁻¹ = 1 + circInv u := by
  unfold circInv; grind

/-- At the wall there is none: nothing mult-adds with −1 to give 0. -/
theorem no_circInv_at_wall (v : K) : circ (-1) v ≠ 0 := by
  unfold circ; grind

/-- Mult-add subtraction: a ⊟ b = (a − b)/(1 + b). -/
def msub (a b : K) : K := (a - b) / (1 + b)

/-- Subtraction undoes mult-add (Fisher's rule for real rates). -/
theorem msub_circ (a b : K) (h : 1 + b ≠ 0) : circ (msub a b) b = a := by
  unfold circ msub; grind

/-- The weighted wall: −1/k absorbs everything under ◯[k]. -/
theorem circK_wall (k a : K) (hk : k ≠ 0) : circK k a (-1 / k) = -1 / k := by
  unfold circK; grind

end Field

end CircleOp
