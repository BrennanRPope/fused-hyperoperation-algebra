/-
  Fused Hyperoperation Algebra: functional analysis through mult-add.
  Core Lean 4 only (no mathlib), no `sorry`. Imports the shared library's Base.

  Operators do not commute, so most of this file works over a general (noncommutative)
  ring `[Lean.Grind.Ring R]`: an algebra of operators is such a ring. Operator mult-add is
  a ◯ₒ b = a + b + ab; on a commutative ring it is the library's ◯ (`opCirc_comm_eq`).

    1. The lift 1 + a ◯ₒ b = (1 + a)(1 + b), associativity, identity 0, the wall −1 on both sides.
    2. Quasi-inverses (Perlis 1942: r + s + rs = 0): a ◯ₒ b = 0 iff (1 + a)(1 + b) = 1, left and
       right quasi-inverses agree, the reverse-order rule (a ◯ₒ b)' = b' ◯ₒ a', Jacobson's circle
       a + b − ab as the k = −1 dial, and the resolvent identity a' − b' = (1 + a')(b − a)(1 + b').
    3. Neumann series: the partial sum s_n = Σ_{j=1}^{n} (−a)^j misses by exactly one term,
       a ◯ₒ s_n = −(−a)^{n+1}; so a nilpotent a has a finite quasi-inverse.
    4. Commuting operators mult-add in either order; conjugation and characters respect ◯ₒ;
       the operator dial with a central weight.
    5. Two noncommuting nilpotent 2×2 matrices: their mult-add has det(1 + ·) = 1 and trace −t
       (eigenvalues on the circle |1 + z| = 1), and at t = 1 it is a sixth mult-add root of 0.
    6. Gap probabilities: the k = −1 dial fold of a list is 1 − ∏(1 − λᵢ).
-/
import FusedAlgebra.Base

namespace FunctionalAnalysis
open CircleOp

open scoped CircleOp

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.Ring R]

/-- Operator mult-add: a ◯ₒ b = a + b + ab (order matters). -/
def opCirc (a b : R) : R := a + b + a * b

local infixl:65 " ◯ₒ " => opCirc

/-! ## 1. The lift and the wall -/

theorem op_lift (a b : R) : 1 + (a ◯ₒ b) = (1 + a) * (1 + b) := by unfold opCirc; grind

theorem op_assoc (a b c : R) : (a ◯ₒ b) ◯ₒ c = a ◯ₒ (b ◯ₒ c) := by unfold opCirc; grind

theorem op_zero (a : R) : a ◯ₒ 0 = a := by unfold opCirc; grind

theorem zero_op (a : R) : 0 ◯ₒ a = a := by unfold opCirc; grind

/-- −1 absorbs from the left… -/
theorem wall_left (b : R) : (-1) ◯ₒ b = -1 := by unfold opCirc; grind

/-- …and from the right. -/
theorem wall_right (a : R) : a ◯ₒ (-1) = -1 := by unfold opCirc; grind

/-- Nothing undoes the wall, unless the ring is the zero ring. -/
theorem wall_no_quasi_inverse (b : R) (h : (-1) ◯ₒ b = 0) : (1 : R) = 0 := by
  rw [wall_left] at h; grind

/-! ## 2. Quasi-inverses -/

/-- Perlis's quasi-regularity is invertibility of the lift. -/
theorem perlis_iff (a b : R) : a ◯ₒ b = 0 ↔ (1 + a) * (1 + b) = 1 := by
  constructor
  · intro h; rw [← op_lift, h]; grind
  · intro h
    have e : a ◯ₒ b = (1 + a) * (1 + b) - 1 := by unfold opCirc; grind
    rw [e, h]; grind

/-- A left quasi-inverse and a right quasi-inverse are the same element. -/
theorem quasi_inverse_unique (a b c : R) (hb : a ◯ₒ b = 0) (hc : c ◯ₒ a = 0) : b = c := by
  have e : (c ◯ₒ a) ◯ₒ b = c ◯ₒ (a ◯ₒ b) := op_assoc c a b
  rw [hb, hc, zero_op, op_zero] at e; exact e

/-- An inverse of 1 + a gives a quasi-inverse c − 1. -/
theorem qi_of_inverse (a c : R) (h : (1 + a) * c = 1) : a ◯ₒ (c - 1) = 0 := by
  have e : a ◯ₒ (c - 1) = (1 + a) * c - 1 := by unfold opCirc; grind
  rw [e, h]; grind

/-- Jacobson's circle a + b − ab is the k = −1 dial of mult-add. -/
theorem jacobson_circle (a b : R) : a + b - a * b = -((-a) ◯ₒ (-b)) := by unfold opCirc; grind

/-- Socks and shoes: the quasi-inverse of a ◯ₒ b is b' ◯ₒ a'. -/
theorem qi_reverse (a a' b b' : R) (ha : a ◯ₒ a' = 0) (hb : b ◯ₒ b' = 0) :
    (a ◯ₒ b) ◯ₒ (b' ◯ₒ a') = 0 := by
  rw [op_assoc, ← op_assoc b b' a', hb, zero_op, ha]

theorem qi_reverse_left (a a' b b' : R) (ha : a' ◯ₒ a = 0) (hb : b' ◯ₒ b = 0) :
    (b' ◯ₒ a') ◯ₒ (a ◯ₒ b) = 0 := by
  rw [op_assoc, ← op_assoc a' a b, ha, zero_op, hb]

/-- The resolvent identity in mult-add form: a' − b' = (1 + a')(b − a)(1 + b'). -/
theorem resolvent_identity (a a' b b' : R) (ha : a' ◯ₒ a = 0) (hb : b ◯ₒ b' = 0) :
    a' - b' = (1 + a') * (b - a) * (1 + b') := by
  have h1 : (1 + a') * (1 + a) = 1 := (perlis_iff a' a).mp ha
  have h2 : (1 + b) * (1 + b') = 1 := (perlis_iff b b').mp hb
  have e : (1 + a') * (b - a) * (1 + b')
      = (1 + a') * ((1 + b) * (1 + b')) - ((1 + a') * (1 + a)) * (1 + b') := by grind
  rw [e, h1, h2]; grind

/-! ## 3. Neumann series -/

/-- Powers built by multiplying on the right. -/
def rpow (x : R) : Nat → R
  | 0 => 1
  | n + 1 => rpow x n * x

theorem rpow_comm (x : R) (n : Nat) : x * rpow x n = rpow x n * x := by
  induction n with
  | zero => simp only [rpow]; grind
  | succ n ih => simp only [rpow]; rw [← Lean.Grind.Semiring.mul_assoc, ih]

/-- Partial Neumann sum s_n = (−a) + (−a)² + … + (−a)ⁿ. -/
def nsum (a : R) : Nat → R
  | 0 => 0
  | n + 1 => nsum a n + rpow (-a) (n + 1)

/-- Each partial sum misses by exactly the next term: a ◯ₒ s_n = −(−a)^{n+1}. -/
theorem neumann_miss (a : R) (n : Nat) : a ◯ₒ nsum a n = -(rpow (-a) (n + 1)) := by
  induction n with
  | zero => simp only [nsum, rpow]; unfold opCirc; grind
  | succ n ih =>
    have c := rpow_comm (-a) (n + 1)
    show a ◯ₒ (nsum a n + rpow (-a) (n + 1)) = -(rpow (-a) (n + 1) * (-a))
    have e1 : a ◯ₒ (nsum a n + rpow (-a) (n + 1))
        = (a ◯ₒ nsum a n) + (rpow (-a) (n + 1) + a * rpow (-a) (n + 1)) := by unfold opCirc; grind
    have e2 : -(rpow (-a) (n + 1)) + (rpow (-a) (n + 1) + a * rpow (-a) (n + 1))
        = -((-a) * rpow (-a) (n + 1)) := by grind
    rw [e1, ih, e2, c]

/-- The same from the other side: s_n ◯ₒ a = −(−a)^{n+1}. -/
theorem neumann_miss_right (a : R) (n : Nat) : nsum a n ◯ₒ a = -(rpow (-a) (n + 1)) := by
  induction n with
  | zero => simp only [nsum, rpow]; unfold opCirc; grind
  | succ n ih =>
    show (nsum a n + rpow (-a) (n + 1)) ◯ₒ a = -(rpow (-a) (n + 1) * (-a))
    have e1 : (nsum a n + rpow (-a) (n + 1)) ◯ₒ a
        = (nsum a n ◯ₒ a) + (rpow (-a) (n + 1) + rpow (-a) (n + 1) * a) := by unfold opCirc; grind
    rw [e1, ih]; grind

theorem rpow_neg (a : R) (n : Nat) : rpow (-a) n = rpow (-1 : R) n * rpow a n := by
  induction n with
  | zero => simp only [rpow]; grind
  | succ n ih =>
    simp only [rpow]; rw [ih]
    have c := rpow_comm (-1 : R) n
    grind

/-- A nilpotent operator (aⁿ⁺¹ = 0) has the finite Neumann sum as its quasi-inverse. -/
theorem nilpotent_quasi_inverse (a : R) (n : Nat) (h : rpow a (n + 1) = 0) :
    a ◯ₒ nsum a n = 0 ∧ nsum a n ◯ₒ a = 0 := by
  have z : rpow (-a) (n + 1) = 0 := by rw [rpow_neg, h]; grind
  rw [neumann_miss, neumann_miss_right, z]; grind

/-- Square zero: the quasi-inverse is −a. -/
theorem square_zero_qi (a : R) (h : a * a = 0) : a ◯ₒ (-a) = 0 := by
  have e : a ◯ₒ (-a) = -(a * a) := by unfold opCirc; grind
  rw [e, h]; grind

/-- Cube zero: the quasi-inverse is −a + a². -/
theorem cube_zero_qi (a : R) (h : a * a * a = 0) : a ◯ₒ (-a + a * a) = 0 ∧ (-a + a * a) ◯ₒ a = 0 := by
  have e1 : a ◯ₒ (-a + a * a) = a * a * a := by clear h; unfold opCirc; grind
  have e2 : (-a + a * a) ◯ₒ a = a * a * a := by clear h; unfold opCirc; grind
  rw [e1, e2, h]; exact ⟨rfl, rfl⟩

/-! ## 4. Commuting, conjugating, characters, the dial -/

theorem op_comm_of_commute (a b : R) (h : a * b = b * a) : a ◯ₒ b = b ◯ₒ a := by
  unfold opCirc; grind

/-- Changing basis (conjugating by an invertible u) respects mult-add. -/
theorem op_conj (u v a b : R) (h : v * u = 1) :
    (u * a * v) ◯ₒ (u * b * v) = u * (a ◯ₒ b) * v := by
  unfold opCirc
  have e : u * a * v * (u * b * v) = u * a * (v * u) * b * v := by grind
  rw [e, h]; grind

/-- A character (additive, multiplicative map) carries ◯ₒ to ◯ₒ. -/
theorem character_op {S : Type} [Lean.Grind.Ring S] (f : R → S)
    (hadd : ∀ x y, f (x + y) = f x + f y) (hmul : ∀ x y, f (x * y) = f x * f y) (a b : R) :
    f (a ◯ₒ b) = f a ◯ₒ f b := by
  unfold opCirc; rw [hadd, hadd, hmul]

/-- The operator dial with a weight k that commutes with a: 1 + k(a ◯ₒ[k] b) = (1 + ka)(1 + kb). -/
theorem op_dial_lift (k a b : R) (hk : a * k = k * a) :
    1 + k * (a + b + k * a * b) = (1 + k * a) * (1 + k * b) := by
  have e : (1 + k * a) * (1 + k * b) = 1 + k * a + k * b + k * (a * k) * b := by grind
  rw [e, hk]; grind

end Ring

/-- On a commutative ring operator mult-add is the library's mult-add. -/
theorem opCirc_comm_eq {R : Type} [Lean.Grind.CommRing R] (a b : R) : opCirc a b = a ◯ b := rfl

/-! ## 5. Two noncommuting nilpotents -/

structure M2 (R : Type) where
  a : R
  b : R
  c : R
  d : R
deriving DecidableEq

theorem M2.ext' {R : Type} {x y : M2 R} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c)
    (h4 : x.d = y.d) : x = y := by
  cases x; cases y; dsimp only at h1 h2 h3 h4; subst h1; subst h2; subst h3; subst h4; rfl

section M2Ring
variable {R : Type} [Lean.Grind.CommRing R]
namespace M2
def one : M2 R := ⟨1, 0, 0, 1⟩
def add (x y : M2 R) : M2 R := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def mul (x y : M2 R) : M2 R :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def circ (x y : M2 R) : M2 R := add (add x y) (mul x y)
def det (x : M2 R) : R := x.a * x.d - x.b * x.c
def tr (x : M2 R) : R := x.a + x.d
end M2

/-- Trace does not see the order of a mult-add. -/
theorem tr_circ_swap (x y : M2 R) : M2.tr (M2.circ x y) = M2.tr (M2.circ y x) := by
  simp only [M2.tr, M2.circ, M2.add, M2.mul]; grind

/-- Neither does the shifted determinant. -/
theorem det_circ_swap (x y : M2 R) :
    M2.det (M2.add M2.one (M2.circ x y)) = M2.det (M2.add M2.one (M2.circ y x)) := by
  simp only [M2.det, M2.add, M2.circ, M2.mul, M2.one]; grind

/-- N = [[0,1],[0,0]] and t·N' = [[0,0],[−t,0]], both nilpotent: their mult-add has
    det(1 + ·) = 1 and trace −t, so its eigenvalues λ satisfy (1 + λ)(1 + λ̄) = 1:
    they ride the circle |1 + z| = 1 and never reach the wall. -/
theorem nilpotent_pair (t : R) :
    M2.det (M2.add M2.one (M2.circ ⟨0, 1, 0, 0⟩ ⟨0, 0, -t, 0⟩)) = 1 ∧
    M2.tr (M2.circ ⟨0, 1, 0, 0⟩ ⟨0, 0, -t, 0⟩) = -t := by
  simp only [M2.det, M2.tr, M2.add, M2.circ, M2.mul, M2.one]; constructor <;> grind

/-- At t = 1 the order matters… -/
theorem nilpotent_pair_noncomm :
    M2.circ (⟨0, 1, 0, 0⟩ : M2 Int) ⟨0, 0, -1, 0⟩ ≠ M2.circ ⟨0, 0, -1, 0⟩ ⟨0, 1, 0, 0⟩ := by decide

/-- …and X = N ◯ₒ N' satisfies X² + X + 1 = 0, so its lift cubes to −1: X is a sixth
    mult-add root of 0, though N and N' each have spectrum {0}. -/
theorem nilpotent_pair_sixth_root :
    let X : M2 Int := M2.circ ⟨0, 1, 0, 0⟩ ⟨0, 0, -1, 0⟩
    let Y := M2.add M2.one X
    M2.mul (M2.mul Y Y) Y = ⟨-1, 0, 0, -1⟩ ∧
    M2.mul (M2.mul (M2.mul Y Y) Y) (M2.mul (M2.mul Y Y) Y) = M2.one := by
  decide

end M2Ring

/-! ## 6. Gap probabilities: the k = −1 dial fold -/

section Gap
variable {R : Type} [Lean.Grind.CommRing R]

/-- The k-dial fold of a list. -/
def foldW (k : R) : List R → R
  | [] => 0
  | x :: l => circK k x (foldW k l)

theorem foldW_lift (k : R) (l : List R) :
    1 + k * foldW k l = prodL (l.map (fun x => 1 + k * x)) := by
  induction l with
  | nil => simp only [foldW, List.map, prodL]; grind
  | cons x l ih => simp only [foldW, List.map, prodL]; rw [circK_lift, ih]

/-- Chance of at least one point = λ₁ ◯₋₁ λ₂ ◯₋₁ … = 1 − ∏(1 − λᵢ) = 1 − det(I − K). -/
theorem gap_fold (l : List R) : foldW (-1) l = 1 - prodL (l.map (fun x => 1 - x)) := by
  have h := foldW_lift (-1) l
  have e : (fun x : R => 1 + -1 * x) = (fun x => 1 - x) := by funext x; grind
  rw [e] at h; grind

/-- At k = 1 the same fold is the shifted Fredholm determinant det(I + K) − 1. -/
theorem fredholm_fold (l : List R) : foldW 1 l = prodL (l.map (fun x => 1 + x)) - 1 := by
  have h := foldW_lift 1 l
  have e : (fun x : R => 1 + 1 * x) = (fun x => 1 + x) := by funext x; grind
  rw [e] at h; grind

end Gap

end FunctionalAnalysis
