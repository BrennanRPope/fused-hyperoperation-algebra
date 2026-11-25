/-
  FusedAlgebra.Base: the shared lemmas of Fused Hyperoperation Algebra
  ====================================================================

  Every topic module imports this file, so each law is proved once.
  Core Lean 4 only (no mathlib, no `sorry`). Written for Lean 4.37.

  Mult-add:            a ⊞ b = a + b + a·b = (1 + a)(1 + b) − 1
  Weighted family:     a ⊞ₖ b = a + b + k·a·b,  1 + k(a ⊞ₖ b) = (1 + ka)(1 + kb)
  Mult-add power:      a^(⊞n) = (1 + a)^n − 1
  Way back (inverse):  a ⊞ (−a/(1 + a)) = 0, everywhere but the wall a = −1

  The laws are stated once for any commutative semiring (so they hold on
  Nat, Int and Rat at the same time); the ones that need minus are stated
  for rings, and the inverse for fields. `grind` does the ring arithmetic.
-/

namespace MultAdd

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Mult-add on any commutative semiring -/

section Semiring
variable {R : Type} [CommSemiring R]

/-- Mult-add: a ⊞ b = a + b + a·b. -/
def multAdd (a b : R) : R := a + b + a * b

@[inherit_doc] scoped infixl:65 " ⊞ " => multAdd

theorem multAdd_def (a b : R) : multAdd a b = a + b + a * b := rfl

/-- The lift: mult-add is multiplication seen through x ↦ 1 + x. -/
theorem multAdd_lift (a b : R) : 1 + multAdd a b = (1 + a) * (1 + b) := by
  unfold multAdd; grind

/-- The lift with the 1 on the right (handy on Nat, where we add 1 to both sides). -/
theorem multAdd_lift' (a b : R) : multAdd a b + 1 = (a + 1) * (b + 1) := by
  unfold multAdd; grind

theorem multAdd_comm (a b : R) : multAdd a b = multAdd b a := by
  unfold multAdd; grind

theorem multAdd_assoc (a b c : R) :
    multAdd (multAdd a b) c = multAdd a (multAdd b c) := by
  unfold multAdd; grind

/-- 0 is the identity of mult-add (because 1 is the identity of ×). -/
theorem multAdd_zero (a : R) : multAdd a 0 = a := by
  unfold multAdd; grind

theorem zero_multAdd (a : R) : multAdd 0 a = a := by
  unfold multAdd; grind

/-- Mult-add is never less than its two pieces added: the corner a·b is extra. -/
theorem multAdd_excess (a b : R) : multAdd a b = (a + b) + a * b := rfl

/-! ### The weighted family a ⊞ₖ b = a + b + k·a·b -/

/-- Weighted mult-add: the corner carries a coupling weight k. -/
def wAdd (k a b : R) : R := a + b + k * a * b

theorem wAdd_def (k a b : R) : wAdd k a b = a + b + k * a * b := rfl

/-- The weighted lift: 1 + k·(a ⊞ₖ b) = (1 + k·a)(1 + k·b). -/
theorem wAdd_lift (k a b : R) : 1 + k * wAdd k a b = (1 + k * a) * (1 + k * b) := by
  unfold wAdd; grind

/-- Weighted mult-add is plain mult-add seen through x ↦ k·x. -/
theorem wAdd_scale (k a b : R) : k * wAdd k a b = multAdd (k * a) (k * b) := by
  unfold wAdd multAdd; grind

theorem wAdd_comm (k a b : R) : wAdd k a b = wAdd k b a := by
  unfold wAdd; grind

theorem wAdd_assoc (k a b c : R) : wAdd k (wAdd k a b) c = wAdd k a (wAdd k b c) := by
  unfold wAdd; grind

theorem wAdd_zero (k a : R) : wAdd k a 0 = a := by
  unfold wAdd; grind

/-- k = 0: no corner, plain addition. -/
theorem wAdd_weight_zero (a b : R) : wAdd 0 a b = a + b := by
  unfold wAdd; grind

/-- k = 1: plain mult-add. -/
theorem wAdd_one (a b : R) : wAdd 1 a b = multAdd a b := by
  unfold wAdd multAdd; grind

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

/-! ### Powers: a ⊞ a ⊞ … ⊞ a climbs the cube -/

/-- The mult-add power a^(⊞n): n copies of a, mult-added. -/
def mpow (a : R) : Nat → R
  | 0 => 0
  | n + 1 => multAdd (mpow a n) a

/-- 1 + a^(⊞n) = (1 + a)^n. -/
theorem mpow_lift (a : R) (n : Nat) : 1 + mpow a n = (1 + a) ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; simp only [mpow]; grind
  | succ n ih => rw [Semiring.pow_succ, ← ih, mpow, multAdd_lift]

/-- Powers add: a^(⊞(m+n)) = a^(⊞m) ⊞ a^(⊞n). -/
theorem mpow_add (a : R) (m n : Nat) : mpow a (m + n) = multAdd (mpow a m) (mpow a n) := by
  induction n with
  | zero => simp only [Nat.add_zero, mpow]; rw [multAdd_zero]
  | succ n ih => rw [← Nat.add_assoc, mpow, mpow, ih, multAdd_assoc]

/-- The weighted power: n copies of a under ⊞ₖ. -/
def wpow (k a : R) : Nat → R
  | 0 => 0
  | n + 1 => wAdd k (wpow k a n) a

/-- 1 + k·a^(⊞ₖn) = (1 + k·a)^n. -/
theorem wpow_lift (k a : R) (n : Nat) : 1 + k * wpow k a n = (1 + k * a) ^ n := by
  induction n with
  | zero => rw [Semiring.pow_zero]; simp only [wpow]; grind
  | succ n ih => rw [Semiring.pow_succ, ← ih, wpow, wAdd_lift]

theorem wpow_one (a : R) (n : Nat) : wpow 1 a n = mpow a n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [wpow, mpow, ih, wAdd_one]

/-! ### Lists: x₁ ⊞ x₂ ⊞ … ⊞ xₙ -/

/-- The mult-add of a list. -/
def foldMA : List R → R
  | [] => 0
  | a :: l => multAdd a (foldMA l)

/-- Product of a list. -/
def prodL : List R → R
  | [] => 1
  | a :: l => a * prodL l

/-- Sum of a list. -/
def sumL : List R → R
  | [] => 0
  | a :: l => a + sumL l

/-- The whole is the first part ⊞ the second part. -/
theorem foldMA_append (S T : List R) :
    foldMA (S ++ T) = multAdd (foldMA S) (foldMA T) := by
  induction S with
  | nil => simp only [List.nil_append, foldMA]; rw [zero_multAdd]
  | cons a S ih =>
    simp only [List.cons_append, foldMA]
    rw [ih, multAdd_assoc]

/-- The list lift: 1 + (x₁ ⊞ … ⊞ xₙ) = (1 + x₁)⋯(1 + xₙ). -/
theorem foldMA_lift (l : List R) : 1 + foldMA l = prodL (l.map (fun x => 1 + x)) := by
  induction l with
  | nil => simp only [foldMA, List.map, prodL]; grind
  | cons a l ih => simp only [foldMA, List.map, prodL]; rw [multAdd_lift, ih]

/-- Order never matters: any reshuffle gives the same total. -/
theorem foldMA_perm {l l' : List R} (h : List.Perm l l') : foldMA l = foldMA l' := by
  induction h with
  | nil => rfl
  | cons x _ ih => simp only [foldMA]; rw [ih]
  | swap x y l =>
    simp only [foldMA]
    unfold multAdd; grind
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

/-- The corner expansion: 1 + (x₁ ⊞ … ⊞ xₙ) is the sum over every corner of
    the cube (every sub-list) of the product of its entries. -/
theorem one_add_foldMA_eq_corners (l : List R) :
    1 + foldMA l = sumL ((groups l).map prodL) := by
  induction l with
  | nil => simp only [foldMA, groups, List.map, prodL, sumL] <;> grind
  | cons a l ih =>
    simp only [foldMA, groups, List.map_append]
    rw [sumL_append, sumL_map_cons, ← ih]
    unfold multAdd; grind

end Semiring

/-! ## 2. Laws that need minus (commutative rings) -/

section Ring
variable {R : Type} [CommRing R]

/-- The wall: −1 absorbs everything mult-added to it. -/
theorem neg_one_absorbs (b : R) : multAdd (-1) b = -1 := by
  unfold multAdd; grind

/-- The power in closed form: a^(⊞n) = (1 + a)^n − 1. -/
theorem mpow_eq (a : R) (n : Nat) : mpow a n = (1 + a) ^ n - 1 := by
  have h := mpow_lift a n
  grind

/-- k = −1 is "at least one of two independent events". -/
theorem wAdd_either (a b : R) : wAdd (-1) a b = 1 - (1 - a) * (1 - b) := by
  unfold wAdd; grind

/-- k = −1 is mult-add mirrored: a ⊞₋₁ b = −((−a) ⊞ (−b)). -/
theorem wAdd_neg_one_mirror (a b : R) : wAdd (-1) a b = -(multAdd (-a) (-b)) := by
  unfold wAdd multAdd; grind

/-- Flipping both signs flips the weight. -/
theorem wAdd_flip (k a b : R) : wAdd k (-a) (-b) = -(wAdd (-k) a b) := by
  unfold wAdd; grind

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
def multAddInv (u : K) : K := -u / (1 + u)

/-- Everywhere but the wall there is an exact way back. -/
theorem multAdd_multAddInv (u : K) (h : 1 + u ≠ 0) : multAdd u (multAddInv u) = 0 := by
  unfold multAdd multAddInv; grind

theorem multAddInv_multAdd (u : K) (h : 1 + u ≠ 0) : multAdd (multAddInv u) u = 0 := by
  rw [multAdd_comm]; exact multAdd_multAddInv u h

/-- Under the lift the inverse is the ordinary reciprocal. -/
theorem inv_one_add (u : K) (h : 1 + u ≠ 0) : (1 + u)⁻¹ = 1 + multAddInv u := by
  unfold multAddInv; grind

/-- At the wall there is none: nothing mult-adds with −1 to give 0. -/
theorem no_multAddInv_at_wall (v : K) : multAdd (-1) v ≠ 0 := by
  unfold multAdd; grind

/-- Mult-add subtraction: a ⊟ b = (a − b)/(1 + b). -/
def msub (a b : K) : K := (a - b) / (1 + b)

/-- Subtraction undoes mult-add (Fisher's rule for real rates). -/
theorem msub_multAdd (a b : K) (h : 1 + b ≠ 0) : multAdd (msub a b) b = a := by
  unfold multAdd msub; grind

/-- The weighted wall: −1/k absorbs everything under ⊞ₖ. -/
theorem wAdd_wall (k a : K) (hk : k ≠ 0) : wAdd k a (-1 / k) = -1 / k := by
  unfold wAdd; grind

end Field

end MultAdd
