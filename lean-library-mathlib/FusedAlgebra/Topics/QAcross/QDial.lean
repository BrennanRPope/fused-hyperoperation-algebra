module

public import FusedAlgebra.Basic

/-!
# The `q` dial: `q`-integers are weighted mult-add powers

The `q`-integer is `[n]_q = 1 + q + ⋯ + q ^ (n - 1)`, a geometric sum
(`∑ i ∈ Finset.range n, q ^ i`, the sum of mathlib's `geom_sum_eq`, `Finset.geom_sum_mul`).
The rule `[m + n]_q = [m]_q + q ^ m * [n]_q` has been standard since Gauss and Heine. What this
file records is that these facts are the weighted mult-add `wAdd k` with weight `k = q - 1`.

## Main results

* `qnum_eq_wpow`: `[n]_q` is `n` copies of `1` under `wAdd (q - 1)`.
* `qnum_add`: `[m + n]_q = wAdd (q - 1) [m]_q [n]_q`.
* `one_add_mul_qnum`: the lift, `1 + (q - 1) * [n]_q = q ^ n` (mathlib's `Finset.mul_geom_sum`).
* `mul_qnum_eq_mpow`: `(q - 1) * [n]_q` is the mult-add power `mpow (q - 1) n`.
* `tate_glue_n`: `(1 + a) * q ^ n = 1 + (a ✳ (q - 1) * [n]_q)`, the Tate gluing as a mult-add
  step.
* The dial values: `q = 1` counts (`qnum_one`), `q = 2` is mult-add (`qnum_two`), `q = 0` is the
  failure union (`qnum_zero_left`, `failure_or`), `q = -1` is parity (`qnum_neg_one`,
  `parity_xor`).

The core-Lean original proved the dial values over `ℤ`; here they hold in any (commutative)
ring.
-/

@[expose] public section

namespace QDial

open MultAdd Finset

variable {R : Type*}

/-! ### The `q`-integer -/

section Semiring

variable [Semiring R]

/-- The `q`-integer `[n]_q = 1 + q + q ^ 2 + ⋯ + q ^ (n - 1)`. -/
def qnum (q : R) (n : ℕ) : R := ∑ i ∈ range n, q ^ i

@[simp]
theorem qnum_zero (q : R) : qnum q 0 = 0 := rfl

/-- The recursion `[n + 1]_q = 1 + q * [n]_q` (the definition in the core-Lean library). -/
theorem qnum_succ (q : R) (n : ℕ) : qnum q (n + 1) = 1 + q * qnum q n := by
  rw [qnum, sum_range_succ', pow_zero, add_comm, qnum, mul_sum]
  simp only [pow_succ']

/-- `q = 1`, weight `0`: the `q`-integer is plain counting. -/
@[simp]
theorem qnum_one (n : ℕ) : qnum (1 : R) n = n := by simp [qnum]

/-- `q = 0`, weight `-1`: every nonempty count collapses to `1`. -/
@[simp]
theorem qnum_zero_left (n : ℕ) : qnum (0 : R) (n + 1) = 1 := by simp [qnum_succ]

/-- The standard addition rule: `[m + n]_q = [m]_q + q ^ m * [n]_q`. -/
theorem qnum_add_eq_add_pow_mul (q : R) (m n : ℕ) :
    qnum q (m + n) = qnum q m + q ^ m * qnum q n := by
  simp only [qnum, sum_range_add, mul_sum, pow_add]

end Semiring

/-! ### `q`-integers and the weighted mult-add of weight `q - 1` -/

section CommRing

variable [CommRing R]

/-- The lift: `1 + (q - 1) * [n]_q = q ^ n`. So `q`-counting is mult-add counting through
`x ↦ 1 + x`. This is mathlib's `Finset.mul_geom_sum`. -/
theorem one_add_mul_qnum (q : R) (n : ℕ) : 1 + (q - 1) * qnum q n = q ^ n := by
  rw [qnum, mul_geom_sum]; ring

/-- The `q`-integer is the weighted mult-add power of `1` with weight `q - 1`. -/
theorem qnum_eq_wpow (q : R) (n : ℕ) : qnum q n = wpow (q - 1) 1 n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [qnum_succ, wpow_succ, ← ih, wAdd_def]; ring

/-- `q`-integers add by the weighted mult-add of weight `q - 1`. -/
theorem qnum_add (q : R) (m n : ℕ) : qnum q (m + n) = wAdd (q - 1) (qnum q m) (qnum q n) := by
  rw [qnum_add_eq_add_pow_mul, ← one_add_mul_qnum q m, wAdd_def]; ring

/-- `(q - 1) * [n]_q` is the mult-add power of `q - 1`. -/
theorem mul_qnum_eq_mpow (q : R) (n : ℕ) : (q - 1) * qnum q n = mpow (q - 1) n := by
  rw [mpow_eq, show (1 : R) + (q - 1) = q by ring, ← one_add_mul_qnum q n]; ring

/-- One Tate step: multiplying the lift `1 + a` by `q` is a mult-add step by `q - 1`. -/
theorem tate_glue_step (q a : R) : (1 + a) * q = 1 + (a ✳ (q - 1)) := by
  rw [one_add_multAdd]; ring

/-- `n` Tate steps: going `n` times round the `q`-direction of the torus is one mult-add step of
length `(q - 1) * [n]_q`. -/
theorem tate_glue_n (q a : R) (n : ℕ) : (1 + a) * q ^ n = 1 + (a ✳ (q - 1) * qnum q n) := by
  rw [one_add_multAdd, one_add_mul_qnum]

/-! ### The dial values -/

/-- `q = 2`, weight `1`: `[n]_2 = 2 ^ n - 1`, the mult-add power of `1` (and the number of
nonzero vectors in an `n`-dimensional space over the field with two elements). -/
theorem qnum_two (n : ℕ) : qnum (2 : R) n = mpow 1 n := by
  have h := mul_qnum_eq_mpow (2 : R) n
  rwa [show (2 : R) - 1 = 1 by ring, one_mul] at h

/-- `q = -1`, weight `-2`: parity. -/
theorem qnum_neg_one (n : ℕ) : qnum (-1 : R) (2 * n) = 0 ∧ qnum (-1 : R) (2 * n + 1) = 1 := by
  induction n with
  | zero => simp [qnum_succ]
  | succ n ih =>
    have e : 2 * (n + 1) = 2 * n + 1 + 1 := by ring
    rw [e, qnum_succ (-1 : R) (2 * n + 1 + 1), qnum_succ (-1 : R) (2 * n + 1), ih.2]
    simp

/-- The parity law is the weighted mult-add of weight `-2`: on `{0, 1}` it is "exclusive or". -/
theorem parity_xor :
    wAdd (-2 : R) 0 0 = 0 ∧ wAdd (-2 : R) 0 1 = 1 ∧
    wAdd (-2 : R) 1 0 = 1 ∧ wAdd (-2 : R) 1 1 = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [wAdd_def] <;> ring

/-- The failure union of weight `-1` on `{0, 1}` is "or". -/
theorem failure_or :
    wAdd (-1 : R) 0 0 = 0 ∧ wAdd (-1 : R) 0 1 = 1 ∧
    wAdd (-1 : R) 1 0 = 1 ∧ wAdd (-1 : R) 1 1 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [wAdd_def] <;> ring

/-- `[3]_q` at the four dial points: counting `3`, volume `7`, union `1`, parity `1`. -/
theorem qnum_three_at_dial :
    qnum (1 : R) 3 = 3 ∧ qnum (2 : R) 3 = 7 ∧ qnum (0 : R) 3 = 1 ∧ qnum (-1 : R) 3 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp only [qnum, sum_range_succ, sum_range_zero] <;> ring

end CommRing

end QDial
