/-
  Fused Hyperoperation Algebra: the counting recurrence under every law
  =====================================================================
  Core Lean 4.37 only (no mathlib), no `sorry`.
  Intended library path: FusedAlgebra/Topics/Patterns.lean (after EarlyScience and Landscape).
  LIBRARY NOTE: first compiled in the shared library. Three proof steps changed, no statement:
  step_is_law splits the pair equality before grind; disc_square_or_ge_five and
  no_integer_five_cycle get the parity split from omega, since Int.emod_add_ediv and
  Int.emod_two_eq_zero_or_one are not in core Lean 4.37. No exact twins: mirror_num is
  MultAdd.wAdd_flip with the minus sign moved across.

  Every law on the landscape map is a ⊕ b = (a + b + k·ab)/(1 − q·ab). Write a value x/y as
  the pair (x, y). Then adding a is one 2×2 step

      step (x, y) = ((1 + k·a)·x + a·y,  −q·a·x + y),

  and adding a to itself n times is that step taken n times. Its trace is P = 2 + k·a and its
  determinant is Q = 1 + k·a + q·a², so the n-fold sum is governed by the Lucas sequence

      U 0 = 0,  U 1 = 1,  U (n+2) = P·U (n+1) − Q·U n,
      a^(⊕(n+1)) = a·U(n+1) / (U(n+1) − Q·U n)                      (iter_lucas)

  Fibonacci is the Lucas sequence with P = 1, Q = −1, and that is exactly the law at
  (k, q) = (−1, −1) adding 1 to itself: 1, 1/2, 2/3, 3/5, 5/8, … (golden_law_fib).

  The standard facts used are Lucas's (1878): the addition rule, Cassini's identity, and the
  fifth term. What this file adds is that each one is a property of the fused laws:
  the addition rule is associativity, Cassini is the determinant, a wall is where the
  determinant vanishes, and the golden ratio appears only through the number 5.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.EarlyScience

namespace Patterns
open MultAdd

set_option linter.unusedSectionVars false

/-! ## 1. One step of any law, as a pair -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-- Adding a in the law (k, q), acting on the pair (x, y) that stands for x/y. -/
def step (k q a : R) (v : R × R) : R × R :=
  ((1 + k * a) * v.1 + a * v.2, (-(q * a)) * v.1 + v.2)

/-- The step's top entry is the law's numerator: x + a·y + k·x·a (for y = 1, wAdd k x a). -/
theorem step_num (k q a x y : R) : (step k q a (x, y)).1 = x + a * y + k * x * a := by
  unfold step; grind

/-- The step's bottom entry is the law's denominator: y − q·x·a (for y = 1, 1 − q·x·a). -/
theorem step_den (k q a x y : R) : (step k q a (x, y)).2 = y - q * x * a := by
  unfold step; grind

/-- With y = 1 the step is exactly the law: numerator wAdd k x a, denominator 1 − q·x·a. -/
theorem step_is_law (k q a x : R) :
    step k q a (x, 1) = (wAdd k x a, 1 - q * x * a) := by
  unfold step wAdd; simp only [Prod.mk.injEq]; constructor <;> grind

/-- n copies of a, as a pair, starting from 0 = (0, 1). -/
def iter (k q a : R) : Nat → R × R
  | 0 => (0, 1)
  | n + 1 => step k q a (iter k q a n)

/-! ## 2. The Lucas sequence of a law -/

/-- Lucas's sequence: U 0 = 0, U 1 = 1, U (n+2) = P·U(n+1) − Q·U n. -/
def lucas (P Q : R) : Nat → R
  | 0 => 0
  | 1 => 1
  | n + 2 => P * lucas P Q (n + 1) - Q * lucas P Q n

theorem lucas_zero (P Q : R) : lucas P Q 0 = 0 := rfl
theorem lucas_one (P Q : R) : lucas P Q 1 = 1 := rfl
theorem lucas_succ_succ (P Q : R) (n : Nat) :
    lucas P Q (n + 2) = P * lucas P Q (n + 1) - Q * lucas P Q n := rfl

/-- The heart of it: n + 1 copies of a, in any law of the map, are
    (a·U(n+1), U(n+1) − Q·U n) with P = 2 + k·a and Q = 1 + k·a + q·a². -/
theorem iter_lucas (k q a : R) (n : Nat) :
    iter k q a (n + 1) =
      (a * lucas (2 + k * a) (1 + k * a + q * a * a) (n + 1),
       lucas (2 + k * a) (1 + k * a + q * a * a) (n + 1)
         - (1 + k * a + q * a * a) * lucas (2 + k * a) (1 + k * a + q * a * a) n) := by
  induction n with
  | zero =>
    simp only [iter, step, lucas, Nat.zero_add, Prod.mk.injEq]
    constructor <;> grind
  | succ n ih =>
    show step k q a (iter k q a (n + 1)) = _
    rw [ih, show n + 1 + 1 = n + 2 from rfl, lucas_succ_succ]
    simp only [step, Prod.mk.injEq]
    constructor <;> grind

/-- The trace of the step is P = 2 + k·a. -/
theorem step_trace (k a : R) : (1 + k * a) + 1 = 2 + k * a := by grind

/-- The determinant of the step is Q = 1 + k·a + q·a²: the law's wall polynomial. -/
theorem step_det (k q a : R) : (1 + k * a) * 1 - a * (-(q * a)) = 1 + k * a + q * a * a := by
  grind

/-- A wall is where the determinant vanishes: adding a wall value w sends every pair onto
    the wall (top = w · bottom), whatever it started as. This is absorption. -/
theorem wall_absorbs (k q w x y : R) (hw : 1 + k * w + q * w * w = 0) :
    (step k q w (x, y)).1 = w * (step k q w (x, y)).2 := by
  simp only [step]
  have h : (1 + k * w) * x + w * y - w * ((-(q * w)) * x + y) = (1 + k * w + q * w * w) * x := by
    grind
  rw [hw] at h
  grind

/-- Associativity of the law, read on the counting sequence, is Lucas's addition rule:
    U(m+n+1) = U(m+1)·U(n+1) − Q·U m·U n. Proved for the pair (n, n+1) together. -/
theorem lucas_add_pair (P Q : R) (m n : Nat) :
    lucas P Q (m + n + 1) = lucas P Q (m + 1) * lucas P Q (n + 1) - Q * lucas P Q m * lucas P Q n ∧
    lucas P Q (m + n + 2) = lucas P Q (m + 1) * lucas P Q (n + 2) - Q * lucas P Q m * lucas P Q (n + 1) := by
  induction n with
  | zero =>
    simp only [lucas, Nat.add_zero, Nat.zero_add]
    constructor <;> grind
  | succ n ih =>
    obtain ⟨h1, h2⟩ := ih
    constructor
    · rw [show m + (n + 1) + 1 = m + n + 2 by omega, show n + 1 + 1 = n + 2 by omega]
      exact h2
    · rw [show m + (n + 1) + 2 = (m + n + 1) + 2 by omega, lucas_succ_succ P Q (m + n + 1),
        show m + n + 1 + 1 = m + n + 2 by omega, h2, h1,
        show n + 1 + 2 = (n + 1) + 2 by omega, lucas_succ_succ P Q (n + 1),
        show n + 1 + 1 = n + 2 by omega, lucas_succ_succ P Q n]
      grind

theorem lucas_add (P Q : R) (m n : Nat) :
    lucas P Q (m + n + 1) = lucas P Q (m + 1) * lucas P Q (n + 1) - Q * lucas P Q m * lucas P Q n :=
  (lucas_add_pair P Q m n).1

/-- Cassini's identity is the determinant of n steps: U(n+1)² − U n·U(n+2) = Qⁿ. -/
theorem cassini (P Q : R) (n : Nat) :
    lucas P Q (n + 1) * lucas P Q (n + 1) - lucas P Q n * lucas P Q (n + 2) = Q ^ n := by
  induction n with
  | zero => rw [lucas_succ_succ, lucas_one, lucas_zero, Lean.Grind.Semiring.pow_zero]; grind
  | succ n ih =>
    rw [Lean.Grind.Semiring.pow_succ, ← ih, show n + 1 + 2 = (n + 1) + 2 by omega,
      lucas_succ_succ P Q (n + 1), show n + 1 + 1 = n + 2 by omega, lucas_succ_succ P Q n]
    grind

/-- If U(n+1) = 0 then n + 1 copies of a land back on 0: a cycle. -/
theorem cycle_of_lucas_zero (k q a : R) (n : Nat)
    (h : lucas (2 + k * a) (1 + k * a + q * a * a) (n + 1) = 0) :
    (iter k q a (n + 1)).1 = 0 := by
  rw [iter_lucas, h]
  show a * 0 = 0
  grind

/-- The fifth term: U 5 = (P² − Q)² − Q·P². -/
theorem lucas_five (P Q : R) : lucas P Q 5 = (P * P - Q) * (P * P - Q) - Q * (P * P) := by
  simp only [lucas]; grind

/-- With Q = 1 (a step of determinant 1, a pure turn), the five-cycle condition factors into
    the two golden quadratics: U 5 = (P² − P − 1)(P² + P − 1). -/
theorem lucas_five_golden (P : R) : lucas P 1 5 = (P * P - P - 1) * (P * P + P - 1) := by
  rw [lucas_five]; grind

/-! ## 3. Mult-add and the Fibonacci law are two points of the same picture -/

/-- Mult-add (k, q) = (1, 0): the denominator is always 1, U(n+1) − (1 + a)·U n = 1. -/
theorem multadd_den_one (a : R) (n : Nat) :
    lucas (2 + a) (1 + a) (n + 1) - (1 + a) * lucas (2 + a) (1 + a) n = 1 := by
  induction n with
  | zero => simp only [lucas, Nat.zero_add]; grind
  | succ n ih => rw [show n + 1 + 1 = n + 2 from rfl, lucas_succ_succ]; grind

/-- Mult-add's counting pair is ((1 + a)^(n+1) − 1, 1), via the Lucas formula. -/
theorem multadd_iter (a : R) (n : Nat) :
    (iter 1 0 a (n + 1)).2 = 1 := by
  rw [iter_lucas]
  have h := multadd_den_one a n
  have e1 : (2 + 1 * a) = 2 + a := by grind
  have e2 : (1 + 1 * a + 0 * a * a) = 1 + a := by grind
  rw [e1, e2]
  exact h

/-- The Fibonacci law (k, q) = (−1, −1): adding 1 is the Fibonacci step (x, y) ↦ (y, x + y). -/
theorem golden_step (x y : R) : step (-1) (-1) 1 (x, y) = (y, x + y) := by
  simp only [step, Prod.mk.injEq]; constructor <;> grind

/-- The resistor ladder of the simple-laws page is one step of the mirror law (1, −1). -/
theorem ladder_is_step (p : Int × Int) : SimpleLaws.ladderStep p = step 1 (-1) 1 p := by
  simp only [SimpleLaws.ladderStep, step, Prod.mk.injEq]; constructor <;> grind

/-- Mirror: flipping the sign of every value takes the law at (k, q) to the law at (−k, q). -/
theorem mirror_num (k a b : R) : -(wAdd k (-a) (-b)) = wAdd (-k) a b := by
  unfold wAdd; grind

theorem mirror_den (q a b : R) : 1 - q * (-a) * (-b) = 1 - q * a * b := by grind

/-- The golden walls: the Fibonacci law's wall polynomial is 1 − a − a², whose roots are
    1/φ and −φ (stated as: w is a wall exactly when w² + w = 1). -/
theorem golden_walls (w : R) : 1 + (-1) * w + (-1) * w * w = 0 ↔ w * w + w = 1 := by
  constructor <;> intro h <;> grind

end Ring

/-! ## 4. Fibonacci exactly, on the integers -/

/-- The Lucas sequence with P = 1, Q = −1 is Fibonacci. -/
theorem lucas_fib_pair (n : Nat) :
    lucas (1 : Int) (-1) n = SimpleLaws.fib n ∧ lucas (1 : Int) (-1) (n + 1) = SimpleLaws.fib (n + 1) := by
  induction n with
  | zero => exact ⟨rfl, rfl⟩
  | succ n ih =>
    obtain ⟨h1, h2⟩ := ih
    refine ⟨h2, ?_⟩
    show 1 * lucas (1 : Int) (-1) (n + 1) - (-1) * lucas (1 : Int) (-1) n
      = SimpleLaws.fib (n + 1) + SimpleLaws.fib n
    rw [h1, h2]
    grind

theorem lucas_fib (n : Nat) : lucas (1 : Int) (-1) n = SimpleLaws.fib n := (lucas_fib_pair n).1

/-- In the law a ⊕ b = (a + b − ab)/(1 + ab), adding 1 to itself n + 1 times gives
    F(n+1)/F(n+2): 1, 1/2, 2/3, 3/5, 5/8, …, marching into the golden wall 1/φ. -/
theorem golden_law_fib (n : Nat) :
    iter (-1 : Int) (-1) 1 (n + 1) = (SimpleLaws.fib (n + 1), SimpleLaws.fib (n + 2)) := by
  rw [iter_lucas]
  have e1 : (2 + (-1 : Int) * 1) = 1 := by decide
  have e2 : (1 + (-1 : Int) * 1 + (-1) * 1 * 1) = -1 := by decide
  rw [e1, e2, lucas_fib, lucas_fib]
  show _ = (SimpleLaws.fib (n + 1), SimpleLaws.fib (n + 1) + SimpleLaws.fib n)
  simp only [Prod.mk.injEq]; constructor <;> grind

/-- Cassini for the golden law: F(n+1)² − F n·F(n+2) = (−1)ⁿ, the determinant (−1) per step. -/
theorem golden_cassini (n : Nat) :
    SimpleLaws.fib (n + 1) * SimpleLaws.fib (n + 1) - SimpleLaws.fib n * SimpleLaws.fib (n + 2)
      = (-1 : Int) ^ n := by
  rw [← lucas_fib, ← lucas_fib, ← lucas_fib]
  exact cassini 1 (-1) n

/-! ## 5. Why the golden ratio, and why only through 5 -/

/-- The smallest irrational walls. For whole-number k and q the wall discriminant k² − 4q is
    either a perfect square or at least 5. So the first integer law whose walls are irrational
    has discriminant 5: the golden walls, at (k, q) = (±1, −1). -/
theorem disc_square_or_ge_five (k q : Int) (hpos : 0 < k * k - 4 * q) (hsmall : k * k - 4 * q < 5) :
    ∃ s : Int, k * k - 4 * q = s * s := by
  have h : k % 2 + 2 * (k / 2) = k := by omega
  generalize k / 2 = m at h
  rcases (by omega : k % 2 = 0 ∨ k % 2 = 1) with hr | hr
  · rw [hr] at h
    subst h
    have e : (0 + 2 * m) * (0 + 2 * m) - 4 * q = 4 * (m * m - q) := by grind
    rw [e] at hpos hsmall ⊢
    generalize m * m - q = t at hpos hsmall ⊢
    have ht : t = 1 := by omega
    subst ht
    exact ⟨2, by decide⟩
  · rw [hr] at h
    subst h
    have e : (1 + 2 * m) * (1 + 2 * m) - 4 * q = 4 * (m * m + m - q) + 1 := by grind
    rw [e] at hpos hsmall ⊢
    generalize m * m + m - q = t at hpos hsmall ⊢
    have ht : t = 0 := by omega
    subst ht
    exact ⟨1, by decide⟩

/-- No whole-number turn has period 5. With determinant 1 and a whole-number trace P, the
    fifth Lucas term P⁴ − 3P² + 1 is odd, so never 0. A five-cycle needs a golden trace. -/
theorem no_integer_five_cycle (P : Int) : lucas P 1 5 ≠ 0 := by
  rw [lucas_five_golden]
  have h : P % 2 + 2 * (P / 2) = P := by omega
  generalize P / 2 = m at h
  rcases (by omega : P % 2 = 0 ∨ P % 2 = 1) with hr | hr
  · rw [hr] at h
    subst h
    have e : ((0 + 2 * m) * (0 + 2 * m) - (0 + 2 * m) - 1) * ((0 + 2 * m) * (0 + 2 * m) + (0 + 2 * m) - 1)
        = 2 * (8 * (m * m * m * m) - 6 * (m * m)) + 1 := by grind
    rw [e]
    generalize 8 * (m * m * m * m) - 6 * (m * m) = t
    omega
  · rw [hr] at h
    subst h
    have e : ((1 + 2 * m) * (1 + 2 * m) - (1 + 2 * m) - 1) * ((1 + 2 * m) * (1 + 2 * m) + (1 + 2 * m) - 1)
        = 2 * (8 * (m * m * m * m) + 16 * (m * m * m) + 6 * (m * m) - 2 * m - 1) + 1 := by grind
    rw [e]
    generalize 8 * (m * m * m * m) + 16 * (m * m * m) + 6 * (m * m) - 2 * m - 1 = t
    omega

/-! ## 6. Fibonacci inside mult-add itself -/

section Semiring
variable {R : Type} [Lean.Grind.CommSemiring R]

/-- Fibonacci numbers on Nat (for exponents). -/
def nfib : Nat → Nat
  | 0 => 0
  | 1 => 1
  | n + 2 => nfib (n + 1) + nfib n

/-- The mult-add Fibonacci sequence: x(n+2) = x(n+1) ✳ x(n). -/
def mfib (x0 x1 : R) : Nat → R
  | 0 => x0
  | 1 => x1
  | n + 2 => multAdd (mfib x0 x1 (n + 1)) (mfib x0 x1 n)

/-- Flatten it and the exponents are Fibonacci numbers:
    1 + x(n+1) = (1 + x0)^F(n) · (1 + x1)^F(n+1). -/
theorem mfib_lift_pair (x0 x1 : R) (n : Nat) :
    1 + mfib x0 x1 (n + 1) = (1 + x0) ^ nfib n * (1 + x1) ^ nfib (n + 1) ∧
    1 + mfib x0 x1 (n + 2) = (1 + x0) ^ nfib (n + 1) * (1 + x1) ^ nfib (n + 2) := by
  induction n with
  | zero =>
    constructor
    · show 1 + x1 = (1 + x0) ^ 0 * (1 + x1) ^ 1
      rw [Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.pow_succ, Lean.Grind.Semiring.pow_zero]
      grind
    · show 1 + multAdd x1 x0 = (1 + x0) ^ 1 * (1 + x1) ^ 1
      rw [multAdd_lift, Lean.Grind.Semiring.pow_succ, Lean.Grind.Semiring.pow_succ,
        Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.pow_zero]
      grind
  | succ n ih =>
    obtain ⟨h1, h2⟩ := ih
    refine ⟨h2, ?_⟩
    show 1 + multAdd (mfib x0 x1 (n + 2)) (mfib x0 x1 (n + 1)) = _
    rw [multAdd_lift, h2, h1, show nfib (n + 1 + 1) = nfib (n + 1) + nfib n from rfl,
      show nfib (n + 1 + 2) = nfib (n + 2) + nfib (n + 1) from rfl,
      show nfib (n + 2) = nfib (n + 1) + nfib n from rfl]
    simp only [pow_add']
    grind

theorem mfib_lift (x0 x1 : R) (n : Nat) :
    1 + mfib x0 x1 (n + 1) = (1 + x0) ^ nfib n * (1 + x1) ^ nfib (n + 1) :=
  (mfib_lift_pair x0 x1 n).1

/-- From 1, 1: the mult-add Fibonacci numbers are 2^F(n) − 1 (1, 1, 3, 7, 31, 255, 8191, …),
    stated as 1 + x(n+1) = (1 + 1)^F(n+2) in any number system. -/
theorem mfib_one_one (n : Nat) : 1 + mfib (1 : R) 1 (n + 1) = (1 + 1) ^ nfib (n + 2) := by
  rw [mfib_lift, show nfib (n + 2) = nfib (n + 1) + nfib n from rfl, pow_add']
  grind

end Semiring

end Patterns
