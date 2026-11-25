/-
  Infinities between countable and the continuum, through the circle operation (called mult-add before 1.0)
  ==================================================================

  Part of Fused Hyperoperation Algebra. Companion page: ../index.html,
  numeric checks: ../checks.py.

  STATUS: core Lean 4 only (no mathlib), no `sorry`. UNCOMPILED: this
  workspace has no Lean toolchain, so nothing here has been run through
  Lean yet. Every identity below is also checked in ../checks.py.

  The circle operation is a ◯ b = a + b + ab = (1+a)(1+b) − 1.

  What is stated and proved here:
  * sumOnes_eq      n ones added make n.                       (the countable side)
  * circOnes_succ   n ones under the circle operation make 2ⁿ − 1.  (the corner side)
  * layers_eval     the circle-operation cube's layers, weighted by t, fill (1+t)ⁿ
                    (layers 3 = [1,3,3,1]).
  * wOnes_spec      the zoomed circle operation a + b + t·ab, folded n times on 1,
                    gives t·Wₙ + 1 = (1+t)ⁿ; t = 0 is plain addition.
  * cantor          no α-indexed list of 0/1 sequences contains them all:
                    the infinite cube has more corners than axes.
  * single_inj      but every axis gives its own corner, so axes ≤ corners.
  * half_iter_upper / half_iter_lower / half_not_poly
                    a "half rung" f with f ∘ f = g is squeezed: if f were
                    polynomial, g would be polynomial. Used for the claim that
                    fractional rungs grow faster than every polynomial.

  LIBRARY NOTE: copied from the finished infinities file (Nov 25); the circle operation and the weighted family come from FusedAlgebra.Base.
-/
import FusedAlgebra.Base

namespace Infinities

open CircleOp

/-! ## Sum of ones versus the circle operation of ones -/

def sumOnes : Nat → Nat
  | 0 => 0
  | n + 1 => sumOnes n + 1

def circOnes : Nat → Nat
  | 0 => 0
  | n + 1 => circ (circOnes n) 1

theorem sumOnes_eq (n : Nat) : sumOnes n = n := by
  induction n with
  | zero => rfl
  | succ n ih => simp [sumOnes, ih]

theorem circOnes_succ (n : Nat) : circOnes n + 1 = 2 ^ n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    show circOnes n + 1 + circOnes n * 1 + 1 = 2 ^ (n + 1)
    rw [Nat.pow_succ, ← ih, Nat.mul_one]
    omega

/-! ## The circle-operation cube's layers (Pascal), weighted by t -/

def addL : List Nat → List Nat → List Nat
  | [], l => l
  | a :: as, [] => a :: as
  | a :: as, b :: bs => (a + b) :: addL as bs

def evalL (t : Nat) : List Nat → Nat
  | [] => 0
  | c :: cs => c + t * evalL t cs

theorem evalL_addL (t : Nat) : ∀ p q : List Nat, evalL t (addL p q) = evalL t p + evalL t q
  | [], q => by simp [addL, evalL]
  | a :: as, [] => by simp [addL, evalL]
  | a :: as, b :: bs => by
      simp only [addL, evalL, evalL_addL t as bs, Nat.mul_add]
      omega

/-- layer k of the n-cube: corners with k switches on. One more axis: each corner
    either leaves the new switch off (same layer) or turns it on (next layer). -/
def layers : Nat → List Nat
  | 0 => [1]
  | n + 1 => addL (layers n) (0 :: layers n)

example : layers 3 = [1, 3, 3, 1] := rfl

theorem layers_eval (t n : Nat) : evalL t (layers n) = (1 + t) ^ n := by
  induction n with
  | zero => simp [layers, evalL]
  | succ n ih =>
    show evalL t (addL (layers n) (0 :: layers n)) = (1 + t) ^ (n + 1)
    rw [evalL_addL]
    simp only [evalL, ih, Nat.zero_add]
    rw [Nat.pow_succ, Nat.mul_add, Nat.mul_one, Nat.mul_comm t]

/-! ## The zoomed circle operation: a + b + t·ab (t = 0 addition, t = 1 the circle operation) -/

def wOnes (t : Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => circK t (wOnes t n) 1

theorem wkey (t w : Nat) : t * (w + 1 + t * w * 1) + 1 = (t * w + 1) * (1 + t) := by
  grind

theorem wOnes_spec (t n : Nat) : t * wOnes t n + 1 = (1 + t) ^ n := by
  induction n with
  | zero => simp [wOnes]
  | succ n ih =>
    show t * (wOnes t n + 1 + t * wOnes t n * 1) + 1 = (1 + t) ^ (n + 1)
    rw [Nat.pow_succ, ← ih]
    exact wkey t (wOnes t n)

theorem wOnes_zero (n : Nat) : wOnes 0 n = n := by
  induction n with
  | zero => rfl
  | succ n ih => simp [wOnes, circK, ih]

/-! ## The infinite cube: more corners than axes (Cantor) -/

theorem cantor {α : Type} (f : α → (α → Bool)) : ∃ g : α → Bool, ∀ a, f a ≠ g := by
  refine ⟨fun a => !(f a a), fun a h => ?_⟩
  have h2 := congrFun h a
  cases hfa : f a a <;> simp [hfa] at h2

def single (n : Nat) : Nat → Bool := fun k => k == n

theorem single_inj {m n : Nat} (h : single m = single n) : m = n := by
  have h2 := congrFun h m
  simp only [single, beq_self_eq_true] at h2
  exact beq_iff_eq.mp h2.symm

/-! ## Half rungs are squeezed between their square and the whole rung -/

theorem half_iter_upper (f g h : Nat → Nat) (N : Nat)
    (hfg : ∀ n, f (f n) = g n)
    (hgrow : ∀ n, n ≤ f n)
    (hmono : ∀ a b, a ≤ b → h a ≤ h b)
    (hbound : ∀ n, N ≤ n → f n ≤ h n) :
    ∀ n, N ≤ n → g n ≤ h (h n) := by
  intro n hn
  rw [← hfg n]
  have h1 : f (f n) ≤ h (f n) := hbound (f n) (Nat.le_trans hn (hgrow n))
  have h2 : h (f n) ≤ h (h n) := hmono _ _ (hbound n hn)
  exact Nat.le_trans h1 h2

theorem half_iter_lower (f g h : Nat → Nat) (N : Nat)
    (hfg : ∀ n, f (f n) = g n)
    (hfmono : ∀ a b, a ≤ b → f a ≤ f b)
    (hh : ∀ n, N ≤ n → N ≤ h n)
    (hbound : ∀ n, N ≤ n → h n ≤ f n) :
    ∀ n, N ≤ n → h (h n) ≤ g n := by
  intro n hn
  rw [← hfg n]
  have h1 : h (h n) ≤ f (h n) := hbound (h n) (hh n hn)
  have h2 : f (h n) ≤ f (f n) := hfmono _ _ (hbound n hn)
  exact Nat.le_trans h1 h2

/-- If a half rung stayed under n^k, the whole rung would stay under n^(k·k).
    Since the whole rung (exponential) outgrows every polynomial, no half rung
    is polynomial. -/
theorem half_not_poly (f g : Nat → Nat) (N k : Nat)
    (hfg : ∀ n, f (f n) = g n)
    (hgrow : ∀ n, n ≤ f n)
    (hbound : ∀ n, N ≤ n → f n ≤ n ^ k) :
    ∀ n, N ≤ n → g n ≤ n ^ (k * k) := by
  intro n hn
  have h3 := half_iter_upper f g (fun m => m ^ k) N hfg hgrow
    (fun a b hab => Nat.pow_le_pow_left hab k) hbound n hn
  rw [Nat.pow_mul]
  exact h3

end Infinities
