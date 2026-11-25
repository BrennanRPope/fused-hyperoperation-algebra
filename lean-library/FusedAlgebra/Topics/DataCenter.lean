/-
  Data center energy in Fused Hyperoperation Algebra.
  Core Lean 4 only (no mathlib), no `sorry`. Quantities are in whole units (e.g. milliwatts),
  so Nat/Int stand in for the real numbers; every statement is a polynomial identity or
  inequality, so it holds for the reals by the same algebra.

  Mult-add (rung 2):  a ⊕ b = a + b + a*b,  and  (a ⊕ b) + 1 = (a+1)(b+1).
  A stage's overhead odds r = (input - output)/output, so input = output * (1 + r).

  LIBRARY NOTE: mult-add and the list fold (foldCirc, the power chain) come from FusedAlgebra.Base; the local copies were removed.
-/
import FusedAlgebra.Base

namespace DataCenterFused

open CircleOp


/-- Product of (1 + r) over the stages. -/
def lifts : List Nat → Nat
  | [] => 1
  | r :: rs => (r + 1) * lifts rs

/-! ### Lead 1. A power chain's overhead is the mult-add of its stages, in any order. -/
theorem foldCirc_lift_nat (rs : List Nat) : foldCirc rs + 1 = lifts rs := by
  induction rs with
  | nil => rfl
  | cons r rs ih => simp only [foldCirc, lifts]; rw [circ_lift', ih]

theorem lifts_append (xs ys : List Nat) : lifts (xs ++ ys) = lifts xs * lifts ys := by
  induction xs with
  | nil => simp [lifts]
  | cons x xs ih => simp only [List.cons_append, lifts]; rw [ih, Nat.mul_assoc]

/-- Swapping two neighbouring stages never changes the total: blame by position is not physical. -/
theorem foldCirc_swap (xs ys : List Nat) (a b : Nat) :
    foldCirc (xs ++ a :: b :: ys) = foldCirc (xs ++ b :: a :: ys) := by
  have h : lifts (xs ++ a :: b :: ys) = lifts (xs ++ b :: a :: ys) := by
    rw [lifts_append, lifts_append]; simp only [lifts]; grind
  have h1 := foldCirc_lift_nat (xs ++ a :: b :: ys)
  have h2 := foldCirc_lift_nat (xs ++ b :: a :: ys)
  omega

/-- Splitting the foldCirc at the IT meter: (outside ⊕ inside) is the whole foldCirc (TUE = PUE·ITUE). -/
theorem foldCirc_split (xs ys : List Nat) : foldCirc (xs ++ ys) = circ (foldCirc xs) (foldCirc ys) := by
  have h := foldCirc_lift_nat (xs ++ ys)
  rw [lifts_append, ← foldCirc_lift_nat xs, ← foldCirc_lift_nat ys, ← circ_lift'] at h
  omega

/-- Position independence: scaling one stage's factor scales the facility total by the same
    amount, wherever that stage sits. -/
theorem scale_anywhere (xs ys : List Nat) (f k : Nat) :
    lifts xs * (f * k) * lifts ys = k * (lifts xs * f * lifts ys) := by grind

/-! ### Lead 2. PUE overhead = electrical ⊕ cooling. The corner p*c is cooling the losses. -/
theorem pue_corner (p c : Nat) : (1 + p) * (1 + c) = 1 + p + c + p * c := by grind
/-- The additive breakdown PUE = 1 + PLF + CLF books the corner as cooling: CLF = c + p*c. -/
theorem additive_books_corner (it p c : Nat) :
    it * (1 + p) * (1 + c) = it + it * p + (it * c + it * p * c) := by grind
/-- One watt of electrical loss removed saves (1 + c) watts at the meter. -/
theorem loss_costs_one_plus_c (it L c : Nat) :
    (it + (L + 1)) * (1 + c) = (it + L) * (1 + c) + (1 + c) := by grind

/-! ### Lead 3. Paralleled units with loss a + c x² each (fixed + square; the proportional
    part does not depend on how the load is split, so it drops out).
    k units sharing load X lose  k*a + c*X²/k.  Multiply through by k(k+1) to stay whole. -/
/-- Exact gap: k(k+1)·[cost(k) − cost(k+1)] = c X² − k(k+1) a. -/
theorem unit_gap (k a c X : Int) :
    (k * k * (k + 1) * a + (k + 1) * c * X ^ 2) - (k * (k + 1) * (k + 1) * a + k * c * X ^ 2)
      = c * X ^ 2 - k * (k + 1) * a := by grind
/-- So add a unit exactly when X² passes (a/c)·k(k+1): the switch point is x*·√(k(k+1)),
    the geometric mean of k and k+1 in units of the sweet spot x* = √(a/c). -/
theorem add_unit_iff (k a c X : Int) :
    k * (k + 1) * (k + 1) * a + k * c * X ^ 2 < k * k * (k + 1) * a + (k + 1) * c * X ^ 2
      ↔ k * (k + 1) * a < c * X ^ 2 := by
  have h := unit_gap k a c X
  constructor <;> intro h' <;> omega
/-- Equal sharing is best for a square loss: shifting d from one unit to the other costs 2d². -/
theorem equal_share (y d : Int) : (y + d) ^ 2 + (y - d) ^ 2 = 2 * y ^ 2 + 2 * d ^ 2 := by grind
/-- At the sweet spot the fixed loss equals the square loss (a = c x²); a numeric check of the
    k → k+1 = 1 → 2 switch with a = 2, c = 1: two units win for X² > 4, lose for X² < 4. -/
example : (1 * 2 * 2 * 2 + 1 * 1 * 3 ^ 2 : Int) < 1 * 1 * 2 * 2 + 2 * 1 * 3 ^ 2 := by decide
example : ¬ ((1 * 2 * 2 * 2 + 1 * 1 * 1 ^ 2 : Int) < 1 * 1 * 2 * 2 + 2 * 1 * 1 ^ 2) := by decide

/-! ### Lead 4. Fans: power ∝ flow³ (affinity law). The mult-add cube. -/
/-- Layers 1, 3, 3, 1 of the fused cube (1+d)³ − 1 = 3d + 3d² + d³. -/
theorem cube_layers (d : Int) : (1 + d) ^ 3 = 1 + 3 * d + 3 * d ^ 2 + d ^ 3 := by grind
/-- k fans each moving q: k·q³·k² = (k·q)³, i.e. total power = Q³/k² for total flow Q = k q. -/
theorem spread (k q : Nat) : k * q ^ 3 * k ^ 2 = (k * q) ^ 3 := by grind
/-- N+1 fans at speed N versus N fans at speed N+1 (same total flow N(N+1)):
    power ratio N² : (N+1)², so running the spare always helps. -/
theorem spare_fan (N : Nat) :
    ((N + 1) * N ^ 3) * (N + 1) ^ 2 = (N * (N + 1) ^ 3) * N ^ 2 := by grind
theorem spare_fan_helps (N : Nat) (h : 0 < N) : (N + 1) * N ^ 3 < N * (N + 1) ^ 3 := by
  have h1 : N ^ 2 < (N + 1) ^ 2 := Nat.pow_lt_pow_left (by omega) (by omega)
  have : (N + 1) * N ^ 3 = (N * (N + 1)) * N ^ 2 := by grind
  have : N * (N + 1) ^ 3 = (N * (N + 1)) * (N + 1) ^ 2 := by grind
  have hp : 0 < N * (N + 1) := Nat.mul_pos h (by omega)
  rw [‹(N + 1) * N ^ 3 = _›, ‹N * (N + 1) ^ 3 = _›]
  exact Nat.mul_lt_mul_of_pos_left h1 hp

/-! ### Stacking savings: kept fractions multiply, so savings compose by mult-add of the
    (negative) odds; in whole numbers, keep k1 of K and k2 of K: kept = k1*k2 of K². -/
theorem stack (K s1 s2 : Int) :
    (K - s1) * (K - s2) = K * K - (s1 * K + s2 * K - s1 * s2) := by grind

end DataCenterFused
