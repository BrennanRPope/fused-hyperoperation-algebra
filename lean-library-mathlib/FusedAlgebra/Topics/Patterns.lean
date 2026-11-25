module

public import FusedAlgebra.Topics.EarlyScience

/-!
# The counting recurrence under every law

Every law on the landscape map is `a ⊕ b = (a + b + k * a * b) / (1 - q * a * b)`. Write a value
`x / y` as the pair `(x, y)`. Then adding `a` is one linear step

`step (x, y) = ((1 + k * a) * x + a * y, -(q * a) * x + y)`,

the action of the matrix `stepMatrix k q a = !![1 + k * a, a; -(q * a), 1]`, and adding `a` to
itself `n` times is that step taken `n` times (`iter`). Its trace is `P = 2 + k * a` and its
determinant is `Q = 1 + k * a + q * a ^ 2`, so the `n`-fold sum is governed by the Lucas sequence
`lucas P Q` (`U 0 = 0`, `U 1 = 1`, `U (n + 2) = P * U (n + 1) - Q * U n`).

## Main results

* `iter_lucas`: `n + 1` copies of `a` are `(a * U (n + 1), U (n + 1) - Q * U n)`.
* `lucas_add`: Lucas's addition rule, which is associativity read on the counting sequence.
* `cassini`: Cassini's identity `U (n + 1) ^ 2 - U n * U (n + 2) = Q ^ n`, the determinant of `n`
  steps.
* `wall_absorbs`: a wall is where the determinant vanishes.
* `iter_multAdd`: at mult-add's point `(1, 0)` the counting pair is `(mpow a n, 1)`.
* `lucas_fib`, `golden_law_fib`, `golden_cassini`: Fibonacci is the Lucas sequence with `P = 1`,
  `Q = -1`, the law at `(k, q) = (-1, -1)` adding `1` to itself.
* `disc_square_or_ge_five`, `no_integer_five_cycle`: the golden ratio appears only through `5`.
* `one_add_mfib`: the mult-add Fibonacci sequence has Fibonacci exponents.

The facts about Lucas sequences are Lucas's (1878). Mathlib has the Fibonacci numbers (`Nat.fib`)
but no general Lucas sequence `U n (P, Q)`, so `lucas` is defined here.
-/

@[expose] public section

namespace Patterns

open MultAdd Matrix

/-! ### One step of any law, as a pair -/

section CommRing

variable {R : Type*} [CommRing R]

/-- Adding `a` in the law `(k, q)`, acting on the pair `(x, y)` that stands for `x / y`. -/
def step (k q a : R) (v : R × R) : R × R :=
  ((1 + k * a) * v.1 + a * v.2, -(q * a) * v.1 + v.2)

/-- The matrix of `step`. -/
def stepMatrix (k q a : R) : Matrix (Fin 2) (Fin 2) R := !![1 + k * a, a; -(q * a), 1]

/-- `step` is the action of `stepMatrix` on column vectors. -/
theorem stepMatrix_mulVec (k q a x y : R) :
    stepMatrix k q a *ᵥ ![x, y] = ![(step k q a (x, y)).1, (step k q a (x, y)).2] := by
  ext i
  fin_cases i <;> simp [stepMatrix, step, mulVec, dotProduct, Fin.sum_univ_two]

/-- The step's top entry is the law's numerator: `x + a * y + k * x * a`. -/
theorem step_num (k q a x y : R) : (step k q a (x, y)).1 = x + a * y + k * x * a := by
  simp only [step]; ring

/-- The step's bottom entry is the law's denominator: `y - q * x * a`. -/
theorem step_den (k q a x y : R) : (step k q a (x, y)).2 = y - q * x * a := by
  simp only [step]; ring

/-- With `y = 1` the step is exactly the law: numerator `wAdd k x a`, denominator
`1 - q * x * a`. -/
theorem step_is_law (k q a x : R) : step k q a (x, 1) = (wAdd k x a, 1 - q * x * a) :=
  Prod.ext (by simp only [step, wAdd_def]; ring) (by simp only [step]; ring)

/-- `n` copies of `a`, as a pair, starting from `0 = (0, 1)`. -/
def iter (k q a : R) (n : ℕ) : R × R := (step k q a)^[n] (0, 1)

@[simp] theorem iter_zero (k q a : R) : iter k q a 0 = (0, 1) := rfl

theorem iter_succ (k q a : R) (n : ℕ) : iter k q a (n + 1) = step k q a (iter k q a n) :=
  Function.iterate_succ_apply' _ n (0, 1)

/-! ### The Lucas sequence of a law -/

/-- Lucas's sequence of the first kind: `U 0 = 0`, `U 1 = 1`,
`U (n + 2) = P * U (n + 1) - Q * U n`. -/
def lucas (P Q : R) : ℕ → R
  | 0 => 0
  | 1 => 1
  | n + 2 => P * lucas P Q (n + 1) - Q * lucas P Q n

@[simp] theorem lucas_zero (P Q : R) : lucas P Q 0 = 0 := rfl

@[simp] theorem lucas_one (P Q : R) : lucas P Q 1 = 1 := rfl

theorem lucas_succ_succ (P Q : R) (n : ℕ) :
    lucas P Q (n + 2) = P * lucas P Q (n + 1) - Q * lucas P Q n := rfl

/-- The heart of it: `n + 1` copies of `a`, in any law of the map, are
`(a * U (n + 1), U (n + 1) - Q * U n)` with `P = 2 + k * a` and `Q = 1 + k * a + q * a ^ 2`. -/
theorem iter_lucas (k q a : R) (n : ℕ) :
    iter k q a (n + 1) =
      (a * lucas (2 + k * a) (1 + k * a + q * a ^ 2) (n + 1),
        lucas (2 + k * a) (1 + k * a + q * a ^ 2) (n + 1) -
          (1 + k * a + q * a ^ 2) * lucas (2 + k * a) (1 + k * a + q * a ^ 2) n) := by
  induction n with
  | zero => simp [iter_succ, step]
  | succ n ih =>
    rw [iter_succ, ih, lucas_succ_succ]
    exact Prod.ext (by simp only [step]; ring) (by simp only [step]; ring)

/-- The trace of the step is `P = 2 + k * a`. -/
theorem step_trace (k q a : R) : (stepMatrix k q a).trace = 2 + k * a := by
  simp [stepMatrix, trace_fin_two]; ring

/-- The determinant of the step is `Q = 1 + k * a + q * a ^ 2`: the law's wall polynomial. -/
theorem step_det (k q a : R) : (stepMatrix k q a).det = 1 + k * a + q * a ^ 2 := by
  simp [stepMatrix, det_fin_two]; ring

/-- A wall is where the determinant vanishes: adding a wall value `w` sends every pair onto the
wall (top = `w` times bottom), whatever it started as. This is absorption. -/
theorem wall_absorbs (k q w x y : R) (hw : 1 + k * w + q * w ^ 2 = 0) :
    (step k q w (x, y)).1 = w * (step k q w (x, y)).2 := by
  simp only [step]; linear_combination x * hw

/-- Associativity of the law, read on the counting sequence, is Lucas's addition rule:
`U (m + n + 1) = U (m + 1) * U (n + 1) - Q * U m * U n`. -/
theorem lucas_add (P Q : R) (m n : ℕ) :
    lucas P Q (m + n + 1) =
      lucas P Q (m + 1) * lucas P Q (n + 1) - Q * lucas P Q m * lucas P Q n := by
  induction n using Nat.twoStepInduction with
  | zero => simp
  | one => simp [lucas_succ_succ, mul_comm]
  | more n ih₁ ih₂ =>
    have e₁ : lucas P Q (m + (n + 2) + 1) =
        P * lucas P Q (m + (n + 1) + 1) - Q * lucas P Q (m + n + 1) :=
      lucas_succ_succ P Q (m + n + 1)
    have e₂ : lucas P Q (n + 2 + 1) = P * lucas P Q (n + 2) - Q * lucas P Q (n + 1) :=
      lucas_succ_succ P Q (n + 1)
    rw [e₁, ih₁, ih₂, e₂, lucas_succ_succ P Q n]
    ring

/-- Cassini's identity is the determinant of `n` steps: `U (n + 1) ^ 2 - U n * U (n + 2) = Q ^ n`.
-/
theorem cassini (P Q : R) (n : ℕ) :
    lucas P Q (n + 1) ^ 2 - lucas P Q n * lucas P Q (n + 2) = Q ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
    have e : lucas P Q (n + 1 + 2) = P * lucas P Q (n + 2) - Q * lucas P Q (n + 1) :=
      lucas_succ_succ P Q (n + 1)
    rw [e, pow_succ Q, ← ih, lucas_succ_succ P Q n]
    ring

/-- If `U (n + 1) = 0` then `n + 1` copies of `a` land back on `0`: a cycle. -/
theorem cycle_of_lucas_zero (k q a : R) (n : ℕ)
    (h : lucas (2 + k * a) (1 + k * a + q * a ^ 2) (n + 1) = 0) : (iter k q a (n + 1)).1 = 0 := by
  rw [iter_lucas, h, mul_zero]

/-- The fifth term: `U 5 = (P ^ 2 - Q) ^ 2 - Q * P ^ 2`. -/
theorem lucas_five (P Q : R) : lucas P Q 5 = (P ^ 2 - Q) ^ 2 - Q * P ^ 2 := by
  simp only [lucas]; ring

/-- With `Q = 1` (a step of determinant `1`, a pure turn), the five-cycle condition factors into
the two golden quadratics: `U 5 = (P ^ 2 - P - 1) * (P ^ 2 + P - 1)`. -/
theorem lucas_five_golden (P : R) : lucas P 1 5 = (P ^ 2 - P - 1) * (P ^ 2 + P - 1) := by
  rw [lucas_five]; ring

/-! ### Mult-add and the Fibonacci law are two points of the same picture -/

/-- Mult-add, `(k, q) = (1, 0)`: the denominator is always `1`,
`U (n + 1) - (1 + a) * U n = 1`. -/
theorem multAdd_den_one (a : R) (n : ℕ) :
    lucas (2 + a) (1 + a) (n + 1) - (1 + a) * lucas (2 + a) (1 + a) n = 1 := by
  induction n with
  | zero => simp
  | succ n ih => rw [lucas_succ_succ]; linear_combination ih

/-- Mult-add's counting pair is `(mpow a n, 1)`, that is `((1 + a) ^ n - 1, 1)`. -/
theorem iter_multAdd (a : R) (n : ℕ) : iter 1 0 a n = (mpow a n, 1) := by
  induction n with
  | zero => rfl
  | succ n ih => rw [iter_succ, ih, step_is_law, wAdd_one, mpow_succ]; simp

/-- Mult-add's counting pair has denominator `1`. -/
theorem multAdd_iter (a : R) (n : ℕ) : (iter 1 0 a (n + 1)).2 = 1 := by
  rw [iter_multAdd]

/-- The Fibonacci law `(k, q) = (-1, -1)`: adding `1` is the Fibonacci step
`(x, y) ↦ (y, x + y)`. -/
theorem golden_step (x y : R) : step (-1) (-1) 1 (x, y) = (y, x + y) := by
  simp [step]

/-- The resistor ladder of `SimpleLaws.ladderStep` is one step of the mirror law `(1, -1)`. -/
theorem ladder_is_step (p : ℤ × ℤ) : SimpleLaws.ladderStep p = step 1 (-1) 1 p := by
  simp [SimpleLaws.ladderStep, step, one_add_one_eq_two]

/-- Mirror: flipping the sign of every value takes the law at `(k, q)` to the law at `(-k, q)`.
This is `MultAdd.wAdd_neg_neg`. -/
theorem mirror_num (k a b : R) : -wAdd k (-a) (-b) = wAdd (-k) a b := by
  rw [wAdd_neg_neg, neg_neg]

theorem mirror_den (q a b : R) : 1 - q * -a * -b = 1 - q * a * b := by ring

/-- The golden walls: the Fibonacci law's wall polynomial is `1 - w - w ^ 2`, whose roots are
`1 / φ` and `-φ`: `w` is a wall exactly when `w ^ 2 + w = 1`. -/
theorem golden_walls (w : R) : 1 + -1 * w + -1 * w ^ 2 = 0 ↔ w ^ 2 + w = 1 := by
  constructor <;> intro h <;> linear_combination -h

/-! ### Fibonacci exactly -/

/-- The Lucas sequence with `P = 1`, `Q = -1` is Fibonacci (mathlib's `Nat.fib`). -/
theorem lucas_fib (n : ℕ) : lucas (1 : R) (-1) n = Nat.fib n := by
  induction n using Nat.twoStepInduction with
  | zero => simp
  | one => simp
  | more n ih₁ ih₂ => rw [lucas_succ_succ, ih₁, ih₂, Nat.fib_add_two]; push_cast; ring

/-- In the law `a ⊕ b = (a + b - a * b) / (1 + a * b)`, adding `1` to itself `n + 1` times gives
`F (n + 1) / F (n + 2)`: `1, 1/2, 2/3, 3/5, 5/8, …`, marching into the golden wall `1 / φ`. -/
theorem golden_law_fib (n : ℕ) :
    iter (-1 : R) (-1) 1 (n + 1) = ((Nat.fib (n + 1) : R), (Nat.fib (n + 2) : R)) := by
  have e₁ : 2 + (-1 : R) * 1 = 1 := by ring
  have e₂ : 1 + (-1 : R) * 1 + -1 * 1 ^ 2 = -1 := by ring
  rw [iter_lucas, e₁, e₂, lucas_fib, lucas_fib, Nat.fib_add_two]
  exact Prod.ext (by simp) (by push_cast; ring)

/-- Cassini for the golden law: `F (n + 1) ^ 2 - F n * F (n + 2) = (-1) ^ n`, the determinant
`-1` per step. -/
theorem golden_cassini (n : ℕ) :
    (Nat.fib (n + 1) : R) ^ 2 - Nat.fib n * Nat.fib (n + 2) = (-1) ^ n := by
  rw [← lucas_fib, ← lucas_fib, ← lucas_fib]
  exact cassini 1 (-1) n

end CommRing

/-! ### Why the golden ratio, and why only through `5` -/

/-- The smallest irrational walls. For whole-number `k` and `q` the wall discriminant
`k ^ 2 - 4 * q` is either a perfect square or at least `5`. So the first integer law whose walls
are irrational has discriminant `5`: the golden walls, at `(k, q) = (±1, -1)`. -/
theorem disc_square_or_ge_five (k q : ℤ) (hpos : 0 < k ^ 2 - 4 * q) (hsmall : k ^ 2 - 4 * q < 5) :
    IsSquare (k ^ 2 - 4 * q) := by
  obtain ⟨m, rfl | rfl⟩ := Int.even_or_odd' k
  · have e : (2 * m) ^ 2 - 4 * q = 4 * (m ^ 2 - q) := by ring
    rw [e] at hpos hsmall ⊢
    generalize m ^ 2 - q = t at hpos hsmall ⊢
    obtain rfl : t = 1 := by lia
    exact ⟨2, by norm_num⟩
  · have e : (2 * m + 1) ^ 2 - 4 * q = 4 * (m ^ 2 + m - q) + 1 := by ring
    rw [e] at hpos hsmall ⊢
    generalize m ^ 2 + m - q = t at hpos hsmall ⊢
    obtain rfl : t = 0 := by lia
    exact ⟨1, by norm_num⟩

/-- No whole-number turn has period `5`. With determinant `1` and a whole-number trace `P`, the
fifth Lucas term `P ^ 4 - 3 * P ^ 2 + 1` is odd, so never `0`. A five-cycle needs a golden
trace. -/
theorem no_integer_five_cycle (P : ℤ) : lucas P 1 5 ≠ 0 := by
  obtain ⟨d, hd⟩ : Even (P ^ 2 + P) := by
    rw [sq, ← mul_add_one]; exact Int.even_mul_succ_self P
  have h : lucas P 1 5 = 2 * (2 * d ^ 2 - 2 * d - 2 * P * d + P) + 1 := by
    rw [lucas_five_golden]; linear_combination (P ^ 2 - P + 2 * d - 2) * hd
  rw [h]
  generalize 2 * d ^ 2 - 2 * d - 2 * P * d + P = t
  lia

/-! ### Fibonacci inside mult-add itself -/

section CommSemiring

variable {R : Type*} [CommSemiring R]

/-- The mult-add Fibonacci sequence: `x (n + 2) = x (n + 1) ✳ x n`. -/
def mfib (x₀ x₁ : R) : ℕ → R
  | 0 => x₀
  | 1 => x₁
  | n + 2 => mfib x₀ x₁ (n + 1) ✳ mfib x₀ x₁ n

@[simp] theorem mfib_zero (x₀ x₁ : R) : mfib x₀ x₁ 0 = x₀ := rfl

@[simp] theorem mfib_one (x₀ x₁ : R) : mfib x₀ x₁ 1 = x₁ := rfl

theorem mfib_succ_succ (x₀ x₁ : R) (n : ℕ) :
    mfib x₀ x₁ (n + 2) = mfib x₀ x₁ (n + 1) ✳ mfib x₀ x₁ n := rfl

/-- Flatten it and the exponents are Fibonacci numbers:
`1 + x (n + 1) = (1 + x₀) ^ F n * (1 + x₁) ^ F (n + 1)`. -/
theorem one_add_mfib (x₀ x₁ : R) (n : ℕ) :
    1 + mfib x₀ x₁ (n + 1) = (1 + x₀) ^ Nat.fib n * (1 + x₁) ^ Nat.fib (n + 1) := by
  induction n using Nat.twoStepInduction with
  | zero => simp
  | one => rw [mfib_succ_succ, one_add_multAdd]; simp [mul_comm]
  | more n ih₁ ih₂ =>
    have e : Nat.fib (n + 2 + 1) = Nat.fib (n + 1) + Nat.fib (n + 2) := Nat.fib_add_two
    rw [mfib_succ_succ, one_add_multAdd, ih₂, ih₁, e, Nat.fib_add_two (n := n)]
    ring

/-- From `1, 1`: the mult-add Fibonacci numbers are `2 ^ F n - 1`
(`1, 1, 3, 7, 31, 255, 8191, …`), stated as `1 + x (n + 1) = (1 + 1) ^ F (n + 2)` in any
commutative semiring. -/
theorem mfib_one_one (n : ℕ) : 1 + mfib (1 : R) 1 (n + 1) = (1 + 1) ^ Nat.fib (n + 2) := by
  rw [one_add_mfib, ← pow_add, Nat.fib_add_two]

end CommSemiring

end Patterns
