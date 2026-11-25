module

public import FusedAlgebra.Basic

/-!
# Closing rate of the fused natural constants

The structural lemmas behind the closing-rate estimate, and the numerical facts about `log 2`
that the rate uses.

## Main results

* `ladder`: if each rung is the one below it run after one step of `s`, then rung `j` is the
  bottom rung run after `j` steps of `s`. Iteration is mathlib's `Function.iterate` (`s^[j]`).
* `grow`: if each step multiplies the gap by at least `p / q`, then after `j` steps it is at
  least `(p / q) ^ j` times the start (stated without division, in any ordered commutative
  semiring; the original was over `ℕ`).
* `two_lt`, `e_pow_lt`, `sq_bound`, `inv_lo`, `inv_hi`: the exact integer inequalities the
  core-Lean original used as stand-ins for `log 2 < 0.7`, `0.69 < log 2`, `2 * 0.7 ^ 2 < 1` and
  the bounds on `1 / log 4`.
* `log_two_lt`, `lt_log_two`, `two_mul_log_two_sq_lt_one`, `inv_log_four_bounds`: with mathlib's
  real numbers the facts themselves are available, from `Real.log_two_lt_d9` and
  `Real.log_two_gt_d9`. These are new in the port.

The real-analysis steps of the rate itself (fixed points of `s`, the escape-time estimate, the
two sandwich bounds) are not formalised here, as in the original.
-/

@[expose] public section

namespace FusedRate

/-- Ladder lemma: if each rung is the one below it run after one step of `s`, then rung `j` is
the bottom rung run after `j` steps of `s`. On `[1, 2]` the fused-linear rungs satisfy this with
`f x = b ^ x` and `s x = b ^ (x - 1)`. -/
theorem ladder {α : Type*} (f s : α → α) (G : ℕ → α → α) (h0 : ∀ x, G 0 x = f x)
    (hs : ∀ j x, G (j + 1) x = G j (s x)) (j : ℕ) (x : α) : G j x = f (s^[j] x) := by
  induction j generalizing x with
  | zero => exact h0 x
  | succ j ih => rw [hs, ih, Function.iterate_succ_apply]

/-- Growth template: if each step multiplies the gap by at least `p / q`, then after `j` steps
it is at least `(p / q) ^ j` times the start. -/
theorem grow {R : Type*} [CommSemiring R] [PartialOrder R] [IsOrderedRing R] (d : ℕ → R)
    {p q : R} (hp : 0 ≤ p) (hq : 0 ≤ q) (h : ∀ j, p * d j ≤ q * d (j + 1)) (j : ℕ) :
    p ^ j * d 0 ≤ q ^ j * d j := by
  induction j with
  | zero => simp
  | succ j ih =>
    calc p ^ (j + 1) * d 0 = p * (p ^ j * d 0) := by ring
      _ ≤ p * (q ^ j * d j) := mul_le_mul_of_nonneg_left ih hp
      _ = q ^ j * (p * d j) := by ring
      _ ≤ q ^ j * (q * d (j + 1)) := mul_le_mul_of_nonneg_left (h j) (pow_nonneg hq j)
      _ = q ^ (j + 1) * d (j + 1) := by ring

/-! ### Integer facts behind the numbers in the rate

Each line is an exact integer inequality. With `(1 + 1 / n) ^ n < e < 2.7182819` they give
`0.69 < log 2 < 0.7`; the real statements are proved directly below. -/

set_option exponentiation.threshold 1000 in
/-- `2 < 1.001 ^ 700`, so `log 2 < 0.7`. -/
theorem two_lt : 2 * 1000 ^ 700 < 1001 ^ 700 := by norm_num

set_option exponentiation.threshold 1000 in
/-- `2.7182819 ^ 69 < 2 ^ 100`, so `0.69 < log 2`. -/
theorem e_pow_lt : 27182819 ^ 69 < 2 ^ 100 * 10 ^ 483 := by norm_num

/-- `2 * 0.7 ^ 2 < 1`. -/
theorem sq_bound : 2 * 70 * 70 < 100 * 100 := by norm_num

/-- `0.714 < 1 / 1.4`. -/
theorem inv_lo : 714 * 14 < 10000 := by norm_num

/-- `1 / 1.38 < 0.7247`. -/
theorem inv_hi : 1000000 < 7247 * 138 := by norm_num

/-! ### The same facts about the real logarithm -/

theorem log_two_lt : Real.log 2 < 0.7 := Real.log_two_lt_d9.trans (by norm_num)

theorem lt_log_two : 0.69 < Real.log 2 := lt_trans (by norm_num) Real.log_two_gt_d9

/-- The fact that makes the upper bound of the rate work. -/
theorem two_mul_log_two_sq_lt_one : 2 * Real.log 2 ^ 2 < 1 := by
  nlinarith [log_two_lt, lt_log_two]

/-- `0.714 < 1 / log 4 < 0.7247`. -/
theorem inv_log_four_bounds : 0.714 < (Real.log 4)⁻¹ ∧ (Real.log 4)⁻¹ < 0.7247 := by
  have h4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  have h₁ := log_two_lt
  have h₂ := lt_log_two
  have hpos : 0 < Real.log 4 := by rw [h4]; linarith
  rw [h4] at hpos ⊢
  constructor
  · rw [lt_inv_comm₀ (by norm_num) hpos]; norm_num; linarith
  · rw [inv_lt_comm₀ hpos (by norm_num)]; norm_num; linarith

end FusedRate
