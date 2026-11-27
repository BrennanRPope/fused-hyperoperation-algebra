/-
  Hat tile and primes (DRAFT, Nov 27 2026). Core Lean 4 only, no imports, no `sorry`.
  Intended library path: FusedAlgebra/Topics/HatPrimes.lean (needs nothing from the library).
  NOT COMPILED where it was written (no Lean toolchain there). Every proof is `decide +kernel` on small
  natural numbers, except `union_rule`, which uses `grind` as the library's other files do.

  1. The hat's substitution on hat counts (H, T, P, F), from Smith–Myers–Kaplan–Goodman-Strauss (2023)
     as written out by Schirmann–Franca–Flicker–Grushin (2023): hats in a level-n H supertile are
     4, 25, 169, 1156, … = F(2n+3)², checked here for n ≤ 11. (checks.py proves it for every n.)
  2. The counts modulo a prime repeat: the step map is the identity after 10 steps mod 11
     (11 = 1 mod 5, period divides p − 1), after 4 steps mod 7 and 14 steps mod 13 (7, 13 = ±2 mod 5,
     period divides p + 1).
  3. The union rule (k = −1): a + b − ab = 1 − (1 − a)(1 − b), the step that turns local removal
     rates into the sieve's density.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.HatPrimes
  (namespace CircleOp.HatPrimes, the library's convention; no proof changed).
-/
namespace CircleOp.HatPrimes

/-- One supertile step on hat counts (H, T, P, F). -/
def step : Nat × Nat × Nat × Nat → Nat × Nat × Nat × Nat
  | (h, t, p, f) => (3*h + t + 3*p + 3*f, h, 2*h + p + 2*f, 2*h + p + 3*f)

/-- Hats in the level-n supertiles H, T, P, F. -/
def hats : Nat → Nat × Nat × Nat × Nat
  | 0 => (4, 1, 2, 2)
  | n + 1 => step (hats n)

/-- Reflected hats (one per level-0 H) in the level-n supertiles. -/
def refl : Nat → Nat × Nat × Nat × Nat
  | 0 => (1, 0, 0, 0)
  | n + 1 => step (refl n)

/-- Fibonacci numbers, F 0 = 0, F 1 = 1. -/
def fib : Nat → Nat
  | 0 => 0
  | 1 => 1
  | n + 2 => fib n + fib (n + 1)

theorem hats_first : (List.range 5).map (fun n => (hats n).1) = [4, 25, 169, 1156, 7921] := by decide +kernel

theorem hats_fib_sq :
    (List.range 12).map (fun n => (hats n).1) = (List.range 12).map (fun n => fib (2*n + 3) ^ 2) := by decide +kernel

theorem refl_first : (List.range 5).map (fun n => (refl n).1) = [1, 3, 22, 147, 1009] := by decide +kernel

/-- The step map reduced mod m. -/
def stepMod (m : Nat) (v : Nat × Nat × Nat × Nat) : Nat × Nat × Nat × Nat :=
  let w := step v
  (w.1 % m, w.2.1 % m, w.2.2.1 % m, w.2.2.2 % m)

def iter {α : Type} : Nat → (α → α) → α → α
  | 0, _, x => x
  | n + 1, f, x => iter n f (f x)

def basis : List (Nat × Nat × Nat × Nat) := [(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]

/-- 11 is 1 mod 5: the counts mod 11 repeat after 10 = 11 − 1 steps. -/
theorem period_11 : basis.map (iter 10 (stepMod 11)) = basis := by decide +kernel
/-- 7 is 2 mod 5: the counts mod 7 repeat after 4 steps (4 divides 7 + 1). -/
theorem period_7 : basis.map (iter 4 (stepMod 7)) = basis := by decide +kernel
/-- 13 is 3 mod 5: the counts mod 13 repeat after 14 = 13 + 1 steps. -/
theorem period_13 : basis.map (iter 14 (stepMod 13)) = basis := by decide +kernel

/-- The union rule (k = −1 in a + b + k·ab) is the complement of a product. -/
theorem union_rule (a b : Int) : a + b + (-1) * (a * b) = 1 - (1 - a) * (1 - b) := by grind

end CircleOp.HatPrimes
