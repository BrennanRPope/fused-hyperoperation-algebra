module

public import FusedAlgebra.Basic

/-!
# Finance as mult-add

Compounding is mult-add: two returns `a` and `b` in a row give the return `a ✳ b`, because growth
factors `1 + a` and `1 + b` multiply. Its inverse `msub a b = (a - b) / (1 + b)` is the real rate.

## Main results

* `fisher`: Fisher's rule, `msub nominal inflation ✳ inflation = nominal`.
* `multAddList_eq_growth`: mult-adding a list of yearly returns gives the total growth minus `1`.
* `drag`, `drag_le`: volatility drag, `(1 + m + s) * (1 + m - s) = (1 + m) ^ 2 - s ^ 2`.
* `gordon`: the Gordon growth formula as a finite geometric sum; `gsum_le` and `gsum_pole` are
  the ceiling for `j > 0` and the pole at growth rate equal to discount rate.
* `inv_arith_mean_le_arith_mean_inv`, `arith_mean_pow_le`, `average_factors`: averaging two
  discount factors is never less than discounting at the averaged growth factor, at every horizon.
* `lowest_rate_wins`: the averaged factor is at least half the factor of the lowest rate.

All of this is known mathematics (the Fisher equation, the Gordon model, the AM-HM and power
mean inequalities); what the file records is that each one is a statement about mult-add. The
core-Lean original worked over `ℚ`; here everything holds in any linearly ordered field.
-/

@[expose] public section

namespace FinanceMultAdd

open MultAdd Finset

variable {K : Type*} [Field K]

/-- Fisher's rule: the real rate mult-added to inflation is the nominal rate. -/
theorem fisher (a b : K) (h : 1 + b ≠ 0) : msub a b ✳ b = a := msub_multAdd a b h

/-- Up `s` then down `s` around a mean `m`: the `s × s` square is cut out of the corner. -/
theorem drag (m s : K) : (1 + m + s) * (1 + m - s) = (1 + m) ^ 2 - s ^ 2 := by ring

/-- Mult-adding a list of yearly returns gives the total growth minus `1`. -/
theorem multAddList_eq_growth (rs : List K) : multAddList rs = (rs.map (1 + ·)).prod - 1 :=
  eq_sub_of_add_eq' (one_add_multAddList rs)

/-- The `N`-year price of a unit dividend with discount-over-growth factor `q`:
`q + q ^ 2 + ⋯ + q ^ N`. -/
def gsum (q : K) (N : ℕ) : K := ∑ n ∈ range N, q ^ (n + 1)

/-- The Gordon growth formula: with `q * (1 + j) = 1`, `j * (q + ⋯ + q ^ N) = 1 - q ^ N`. A
dividend growing at `g` and discounted at `r` is a flat perpetuity at `j = msub r g`. -/
theorem gordon (q j : K) (h : q * (1 + j) = 1) (N : ℕ) : j * gsum q N = 1 - q ^ N := by
  induction N with
  | zero => simp [gsum]
  | succ N ih =>
    have hj : j * q = 1 - q := by linear_combination h
    rw [gsum, sum_range_succ, ← gsum, mul_add, ih, pow_succ]
    linear_combination q ^ N * hj

variable [LinearOrder K] [IsStrictOrderedRing K]

/-- Volatility drag never helps. -/
theorem drag_le (m s : K) : (1 + m + s) * (1 + m - s) ≤ (1 + m) ^ 2 := by
  rw [drag]; linarith [sq_nonneg s]

/-- For `j > 0` the price never passes `D / j`. -/
theorem gsum_le (q j : K) (h : q * (1 + j) = 1) (hq : 0 ≤ q) (N : ℕ) : j * gsum q N ≤ 1 := by
  rw [gordon q j h N]; linarith [pow_nonneg hq N]

/-- For `q ≥ 1` (growth at or above the discount rate) the `N`-year price is at least `N`: the
price has no ceiling. This is the pole at `g = r`. -/
theorem gsum_pole {q : K} (hq : 1 ≤ q) (N : ℕ) : (N : K) ≤ gsum q N :=
  calc (N : K) = ∑ _n ∈ range N, (1 : K) := by simp
    _ ≤ gsum q N := sum_le_sum fun _ _ ↦ one_le_pow₀ hq

/-- Chebyshev's sum inequality for two points: `x ^ n - y ^ n` and `x - y` have the same sign. -/
theorem pow_sub_pow_mul_sub_nonneg {x y : K} (hx : 0 ≤ x) (hy : 0 ≤ y) (n : ℕ) :
    0 ≤ (x ^ n - y ^ n) * (x - y) := by
  rcases le_total y x with h | h
  · exact mul_nonneg (sub_nonneg.2 (pow_le_pow_left₀ hy h n)) (sub_nonneg.2 h)
  · exact mul_nonneg_of_nonpos_of_nonpos (sub_nonpos.2 (pow_le_pow_left₀ hx h n))
      (sub_nonpos.2 h)

/-- The power mean inequality for two points: the power of the average is at most the average of
the powers. -/
theorem arith_mean_pow_le {x y : K} (hx : 0 ≤ x) (hy : 0 ≤ y) (n : ℕ) :
    ((x + y) / 2) ^ n ≤ (x ^ n + y ^ n) / 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
    calc ((x + y) / 2) ^ (n + 1) = ((x + y) / 2) ^ n * ((x + y) / 2) := pow_succ _ _
      _ ≤ (x ^ n + y ^ n) / 2 * ((x + y) / 2) := by gcongr
      _ ≤ (x ^ (n + 1) + y ^ (n + 1)) / 2 := by
        rw [pow_succ, pow_succ]
        linarith [pow_sub_pow_mul_sub_nonneg hx hy n]

/-- The AM-HM inequality for two points: the reciprocal of the average is at most the average of
the reciprocals. -/
theorem inv_arith_mean_le_arith_mean_inv {u v : K} (hu : 0 < u) (hv : 0 < v) :
    ((u + v) / 2)⁻¹ ≤ (u⁻¹ + v⁻¹) / 2 := by
  rw [← sub_nonneg]
  have h : (u⁻¹ + v⁻¹) / 2 - ((u + v) / 2)⁻¹ = (u - v) ^ 2 / (2 * u * v * (u + v)) := by
    field_simp; ring
  rw [h]; positivity

/-- Averaging two discount factors `u⁻¹ ^ n` and `v⁻¹ ^ n` is never less than discounting at the
averaged growth factor, for every horizon `n`. This is why long-run discount rates fall when the
rate is uncertain. -/
theorem average_factors {u v : K} (hu : 0 < u) (hv : 0 < v) (n : ℕ) :
    (((u + v) / 2) ^ n)⁻¹ ≤ (u⁻¹ ^ n + v⁻¹ ^ n) / 2 :=
  calc (((u + v) / 2) ^ n)⁻¹ = ((u + v) / 2)⁻¹ ^ n := (inv_pow _ _).symm
    _ ≤ ((u⁻¹ + v⁻¹) / 2) ^ n :=
      pow_le_pow_left₀ (by positivity) (inv_arith_mean_le_arith_mean_inv hu hv) n
    _ ≤ (u⁻¹ ^ n + v⁻¹ ^ n) / 2 := arith_mean_pow_le (by positivity) (by positivity) n

/-- The averaged factor is at least half the factor of the lowest rate, so over long horizons the
lowest-rate corner sets the rate. -/
theorem lowest_rate_wins {x y : K} (hy : 0 ≤ y) (hxy : y ≤ x) (n : ℕ) :
    x ^ n / 2 ≤ (x ^ n + y ^ n) / 2 ∧ (x ^ n + y ^ n) / 2 ≤ x ^ n := by
  have h₁ := pow_nonneg hy n
  have h₂ := pow_le_pow_left₀ hy hxy n
  constructor <;> linarith

end FinanceMultAdd
