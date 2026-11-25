module

public import FusedAlgebra.Cube

/-!
# `e` from mult-add

Mult-add `1 / n` to itself `n` times and add back the unit corner: that is the staircase
`stair n = (1 + 1 / n) ^ n`. On the `n`-cube every switch weighs `1 / n`, so the corners with `k`
switches on weigh `wt n k = n.choose k / n ^ k` in total.

This is the mathlib form of `FusedAlgebra.Topics.E`. The core-Lean original worked over `ℚ` with a
hand-made sum `rsum`, falling factorial `fall` and triangular number `tri`; here the sums are
`Finset.sum` over `Finset.range`, `fall n k` is `Nat.descFactorial n k`, `tri k` is `k.choose 2`,
and everything holds in any linearly ordered field `K` (an explicit argument of `stair`, `wt`,
`s`, `u`).

## Main results

* `choose_mul_factorial_le_pow`, `collisions`: `n.choose k * k !` counts ordered picks of `k`
  different axes, so it is at most `n ^ k`, and the picks that repeat an axis number at most
  `k.choose 2 * n ^ (k - 1)`.
* `wt_le`, `wt_ge`: layer `k` weighs between `1 / k ! - k.choose 2 / (n * k !)` and `1 / k !`.
* `stair_lt_succ`: the staircase climbs, `(1 + 1 / n) ^ n < (1 + 1 / (n + 1)) ^ (n + 1)`.
* `stair_le_s`, `s_lt_three`: `(1 + 1 / n) ^ n ≤ 1 / 0! + ⋯ + 1 / n ! < 3`.
* `u_step`, `stair_le_u`: `u N = s N + 1 / (N ! * N)` only comes down, and every step of the
  staircase sits under every `u N`.
* `s_le_stair`: for `n ≥ N` the staircase is within `3 * N.choose 2 / n` of `s N`.
* `tendsto_stair`, `tendsto_s`, `s_le_exp_one`, `exp_one_le_u`, `exp_one_bounds`: over `ℝ` the
  staircase and the partial sums close on `Real.exp 1`, which `s 10` and `u 10` trap between
  `2.7182818` and `2.7182819`.

All of this is classical (Bernoulli, Euler, Cauchy). The binomial theorem is `add_pow`, the bound
`n.choose k ≤ n ^ k / k !` is `Nat.choose_le_pow_div`, the two limits are
`Real.tendsto_one_add_div_pow_exp` and `NormedSpace.expSeries_div_hasSum_exp`, and mathlib already
has sharper digits of `e` (`Real.exp_one_gt_d9`, `Real.exp_one_lt_d9`). What the file records is
the reading on the mult-add cube and the elementary route to the bounds. The core-Lean original
stopped before the limit, because it had no real numbers; the last section adds it from mathlib.
-/

@[expose] public section

namespace MultAdd.E

open MultAdd Finset
open scoped Nat

/-! ### Counting ordered picks -/

/-- The layer bound: `n.choose k * k ! ≤ n ^ k`, because distinct ordered picks
(`Nat.descFactorial`) are at most all ordered picks. -/
theorem choose_mul_factorial_le_pow (n k : ℕ) : n.choose k * k ! ≤ n ^ k := by
  rw [mul_comm, ← Nat.descFactorial_eq_factorial_mul_choose]
  exact Nat.descFactorial_le_pow n k

/-- Collisions are few: `n * n ^ k ≤ n * n.descFactorial k + k.choose 2 * n ^ k`. The ordered
picks of `k` axes that repeat an axis number at most `k.choose 2 * n ^ (k - 1)`. -/
theorem collisions (n : ℕ) : ∀ k, k ≤ n →
    n * n ^ k ≤ n * n.descFactorial k + k.choose 2 * n ^ k
  | 0, _ => by simp
  | k + 1, h => by
    have ih := collisions n k (by lia)
    have hF := Nat.descFactorial_le_pow n k
    obtain ⟨d, rfl⟩ : ∃ d, n = d + k := ⟨n - k, by lia⟩
    have hT : (k + 1).choose 2 = k.choose 2 + k := by
      rw [Nat.choose_succ_succ', Nat.choose_one_right, add_comm]
    rw [Nat.descFactorial_succ, Nat.add_sub_cancel, hT, pow_succ]
    nlinarith [Nat.mul_le_mul_left (d + k) ih, Nat.mul_le_mul_left ((d + k) * k) hF]

/-! ### Bernoulli and the climbing staircase, in `ℕ` -/

/-- Bernoulli, in mean-value form: `(c + 1) ^ (k + 1) ≤ c ^ (k + 1) + (k + 1) * (c + 1) ^ k`. -/
theorem bernoulli (c : ℕ) : ∀ k, (c + 1) ^ (k + 1) ≤ c ^ (k + 1) + (k + 1) * (c + 1) ^ k
  | 0 => by simp
  | k + 1 => by
    have ih := bernoulli c k
    have h := Nat.mul_le_mul_left c ih
    rw [pow_succ (c + 1) (k + 1), pow_succ c (k + 1)]
    rw [pow_succ (c + 1) k] at ih h ⊢
    nlinarith [Nat.zero_le (k * (c + 1) ^ k), Nat.zero_le ((c + 1) ^ k)]

/-- Bernoulli with room to spare:
`(c + 1) ^ (k + 2) + (k + 1) * (c + 1) ^ k ≤ c ^ (k + 2) + (k + 2) * (c + 1) ^ (k + 1)`. -/
theorem bernoulli_strict (c k : ℕ) :
    (c + 1) ^ (k + 2) + (k + 1) * (c + 1) ^ k ≤ c ^ (k + 2) + (k + 2) * (c + 1) ^ (k + 1) := by
  have ih := bernoulli c k
  have h := Nat.mul_le_mul_left c ih
  rw [pow_succ (c + 1) (k + 1), pow_succ c (k + 1)]
  rw [pow_succ (c + 1) k] at ih h ⊢
  nlinarith

/-- The staircase climbs, cleared of fractions:
`(m + 2) ^ (2 * (m + 1) + 1) < (m + 1) ^ (m + 1) * (m + 3) ^ (m + 2)`. -/
theorem stair_nat_lt (m : ℕ) :
    (m + 2) ^ (2 * (m + 1) + 1) < (m + 1) ^ (m + 1) * (m + 3) ^ (m + 2) := by
  have hb := bernoulli_strict ((m + 1) * (m + 3)) m
  have hc : (m + 1) * (m + 3) + 1 = (m + 2) ^ 2 := by ring
  rw [hc, mul_pow] at hb
  have e0 : ((m + 2) ^ 2) ^ m = (m + 2) ^ (2 * m) := (pow_mul _ 2 m).symm
  have e1 : ((m + 2) ^ 2) ^ (m + 2) = (m + 2) ^ (2 * m) * (m + 2) ^ 4 := by
    rw [pow_add, e0]; ring
  have e2 : ((m + 2) ^ 2) ^ (m + 1) = (m + 2) ^ (2 * m) * (m + 2) ^ 2 := by
    rw [pow_add, e0]; ring
  have e3 : (m + 2) ^ (2 * (m + 1) + 1) = (m + 2) ^ (2 * m) * (m + 2) ^ 3 := by ring
  have e4 : (m + 1) ^ (m + 2) = (m + 1) * (m + 1) ^ (m + 1) := by ring
  rw [e0, e1, e2, e4, mul_assoc (m + 1)] at hb
  rw [e3]
  have hQ : 0 < (m + 2) ^ (2 * m) := by positivity
  generalize (m + 2) ^ (2 * m) = Q at *
  generalize (m + 1) ^ (m + 1) * (m + 3) ^ (m + 2) = R at *
  have hb' : (m + 1) * (Q * (m + 2) ^ 3 + Q) ≤ (m + 1) * R := by nlinarith
  have := Nat.le_of_mul_le_mul_left hb' m.succ_pos
  lia

/-! ### The cube's corners, weighted -/

/-- The cube's corner count, the binomial theorem: `(1 + x) ^ n = ∑ k, n.choose k * x ^ k`. This
is mathlib's `add_pow`, through `Cube.one_add_mpow_eq_sum_choose`. -/
theorem binomial {R : Type*} [CommSemiring R] (n : ℕ) (x : R) :
    (1 + x) ^ n = ∑ k ∈ range (n + 1), (n.choose k : R) * x ^ k := by
  rw [← one_add_mpow, Cube.one_add_mpow_eq_sum_choose]

variable (K : Type*) [Field K]

/-- The staircase: `1 / n` mult-added to itself `n` times, plus the unit corner. -/
def stair (n : ℕ) : K := (1 + 1 / (n : K)) ^ n

theorem stair_eq_mpow (n : ℕ) : stair K n = mpow (1 / (n : K)) n + 1 := by
  rw [add_comm, one_add_mpow, stair]

/-- The weight of layer `k` of the `n`-cube when every switch weighs `1 / n`. -/
def wt (n k : ℕ) : K := (n.choose k : K) * (1 / (n : K)) ^ k

/-- The staircase is the cube, layer by layer. -/
theorem stair_eq_cube (n : ℕ) : stair K n = ∑ k ∈ range (n + 1), wt K n k :=
  binomial n _

/-- The partial sums of `e`: `s N = 1 / 0! + 1 / 1! + ⋯ + 1 / N !`. -/
def s (N : ℕ) : K := ∑ k ∈ range (N + 1), 1 / (k ! : K)

/-- The upper partner: `u N = s N + 1 / (N ! * N)`. -/
def u (N : ℕ) : K := s K N + 1 / ((N ! : K) * (N : K))

variable [LinearOrder K] [IsStrictOrderedRing K]

theorem wt_nonneg (n k : ℕ) : 0 ≤ wt K n k := by
  unfold wt; positivity

/-- Each layer is at most `1 / k !`: `n.choose k / n ^ k ≤ 1 / k !`. -/
theorem wt_le (n k : ℕ) (hn : 1 ≤ n) : wt K n k ≤ 1 / (k ! : K) := by
  have hn' : (0 : K) < n := by exact_mod_cast hn
  have h := Nat.choose_le_pow_div (α := K) k n
  rw [wt, one_div_pow, mul_one_div, div_le_div_iff₀ (by positivity) (by positivity), one_mul]
  rwa [le_div_iff₀ (by positivity)] at h

/-- Each layer is at least `1 / k !` minus its collisions:
`1 / k ! ≤ n.choose k / n ^ k + k.choose 2 / (n * k !)`. -/
theorem wt_ge (n k : ℕ) (hn : 1 ≤ n) (hk : k ≤ n) :
    1 / (k ! : K) ≤ wt K n k + (k.choose 2 : K) / ((n : K) * (k ! : K)) := by
  have hn' : (0 : K) < n := by exact_mod_cast hn
  have hc : (n : K) * (n : K) ^ k ≤
      (n : K) * ((k ! : K) * (n.choose k : K)) + (k.choose 2 : K) * (n : K) ^ k := by
    have := collisions n k hk
    rw [Nat.descFactorial_eq_factorial_mul_choose] at this
    exact_mod_cast this
  rw [wt, one_div_pow, mul_one_div, div_add_div _ _ (by positivity) (by positivity),
    div_le_div_iff₀ (by positivity) (by positivity)]
  have hk' : (0 : K) < (k ! : K) := by positivity
  nlinarith [mul_le_mul_of_nonneg_left hc hk'.le]

theorem stair_le_s (n : ℕ) (hn : 1 ≤ n) : stair K n ≤ s K n := by
  rw [stair_eq_cube]
  exact sum_le_sum fun k _ ↦ wt_le K n k hn

theorem s_bound : ∀ N, s K N ≤ 3 - 2 / (2 : K) ^ N
  | 0 => by norm_num [s]
  | N + 1 => by
    have ih := s_bound N
    have hf : 1 / ((N + 1)! : K) ≤ 1 / (2 : K) ^ N := by
      apply one_div_le_one_div_of_le (by positivity)
      have : 2 ^ N ≤ (N + 1)! := by
        simpa [add_comm] using Nat.factorial_mul_pow_le_factorial (m := 1) (n := N)
      exact_mod_cast this
    have e : 2 / (2 : K) ^ (N + 1) = 1 / (2 : K) ^ N := by
      rw [pow_succ]; field_simp
    rw [s, sum_range_succ, ← s, e]
    have h2 : 2 / (2 : K) ^ N = 2 * (1 / (2 : K) ^ N) := by ring
    linarith

/-- The partial sums never reach `3`. -/
theorem s_lt_three (N : ℕ) : s K N < 3 := by
  have h := s_bound K N
  have : (0 : K) < 2 / (2 : K) ^ N := by positivity
  linarith

theorem s_mono : Monotone (s K) := fun _ _ h ↦
  sum_le_sum_of_subset_of_nonneg (range_mono (by lia)) fun _ _ _ ↦ by positivity

theorem s_le_u (N : ℕ) : s K N ≤ u K N :=
  le_add_of_nonneg_right (by positivity)

/-- The upper partner only comes down: `u (N + 1) ≤ u N` for `N ≥ 1`. -/
theorem u_step (N : ℕ) (hN : 1 ≤ N) : u K (N + 1) ≤ u K N := by
  have hN' : (0 : K) < N := by exact_mod_cast hN
  have hF : (0 : K) < (N ! : K) := by positivity
  have key : 1 / ((N + 1)! : K) + 1 / (((N + 1)! : K) * ((N + 1 : ℕ) : K)) ≤
      1 / ((N ! : K) * (N : K)) := by
    rw [Nat.factorial_succ]
    push_cast
    rw [div_add_div _ _ (by positivity) (by positivity),
      div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_pos (mul_pos hF hF) (mul_pos hF hN')]
  rw [u, u, s, sum_range_succ, ← s]
  linarith

theorem u_anti {N M : ℕ} (hN : 1 ≤ N) (h : N ≤ M) : u K M ≤ u K N := by
  induction M, h using Nat.le_induction with
  | base => exact le_rfl
  | succ M hM ih => exact (u_step K M (hN.trans hM)).trans ih

/-- The ceiling: every step of the staircase sits under every `u N`. -/
theorem stair_le_u (n N : ℕ) (hn : 1 ≤ n) (hN : 1 ≤ N) : stair K n ≤ u K N := by
  rcases le_total n N with h | h
  · exact (stair_le_s K n hn).trans ((s_mono K h).trans (s_le_u K N))
  · exact (stair_le_s K n hn).trans ((s_le_u K n).trans (u_anti K hN h))

/-- The floor: for `n ≥ N` the staircase is within `3 * N.choose 2 / n` of `s N`. -/
theorem s_le_stair (n N : ℕ) (hn : 1 ≤ n) (hNn : N ≤ n) :
    s K N - 3 * (N.choose 2 : K) / (n : K) ≤ stair K n := by
  have hn' : (0 : K) < n := by exact_mod_cast hn
  set c : K := (N.choose 2 : K) / (n : K) with hc
  have hc0 : 0 ≤ c := by positivity
  -- drop the layers above `N`
  have h1 : ∑ k ∈ range (N + 1), wt K n k ≤ stair K n := by
    rw [stair_eq_cube]
    exact sum_le_sum_of_subset_of_nonneg (range_mono (by lia)) fun k _ _ ↦ wt_nonneg K n k
  -- each kept layer is at least `1 / k ! - c / k !`
  have h2 : ∑ k ∈ range (N + 1), (1 - c) * (1 / (k ! : K)) ≤
      ∑ k ∈ range (N + 1), wt K n k := by
    refine sum_le_sum fun k hk ↦ ?_
    have hkN : k ≤ N := by simpa [Nat.lt_succ_iff] using hk
    have hg := wt_ge K n k hn (hkN.trans hNn)
    have ht : (k.choose 2 : K) / ((n : K) * (k ! : K)) ≤ c * (1 / (k ! : K)) := by
      rw [hc, div_mul_div_comm, mul_one]
      gcongr
    linarith
  rw [← mul_sum, ← s] at h2
  have h4 : c * s K N ≤ c * 3 := mul_le_mul_of_nonneg_left (s_lt_three K N).le hc0
  have e : 3 * (N.choose 2 : K) / (n : K) = c * 3 := by rw [hc]; ring
  rw [e]
  linarith

/-- The staircase climbs: `(1 + 1 / n) ^ n < (1 + 1 / (n + 1)) ^ (n + 1)` for `n ≥ 1`. -/
theorem stair_lt_succ (n : ℕ) (hn : 1 ≤ n) : stair K n < stair K (n + 1) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by lia⟩
  have hA : stair K (m + 1) = ((m : K) + 2) ^ (m + 1) / ((m : K) + 1) ^ (m + 1) := by
    rw [stair, ← div_pow]; congr 1; push_cast; field_simp; ring
  have hB : stair K (m + 1 + 1) = ((m : K) + 3) ^ (m + 2) / ((m : K) + 2) ^ (m + 2) := by
    rw [stair, ← div_pow]; congr 1; push_cast; field_simp; ring
  have hnat : ((m : K) + 2) ^ (2 * (m + 1) + 1) <
      ((m : K) + 1) ^ (m + 1) * ((m : K) + 3) ^ (m + 2) := by
    exact_mod_cast stair_nat_lt m
  rw [hA, hB, div_lt_div_iff₀ (by positivity) (by positivity)]
  calc ((m : K) + 2) ^ (m + 1) * ((m : K) + 2) ^ (m + 2)
      = ((m : K) + 2) ^ (2 * (m + 1) + 1) := by ring
    _ < ((m : K) + 1) ^ (m + 1) * ((m : K) + 3) ^ (m + 2) := hnat
    _ = ((m : K) + 3) ^ (m + 2) * ((m : K) + 1) ^ (m + 1) := mul_comm _ _

/-! ### Small checks -/

/-- Layer `k` of the `100`-cube against `1 / k !`: `k = 3` gives `0.9702 / 6` versus `1 / 6`. -/
example : wt ℚ 100 3 = 4851 / 30000 := by
  norm_num [wt, Nat.choose]

/-- `s 10 = 2.71828180…` and `u 10 = 2.71828183…`. -/
theorem s_ten_u_ten : (27182818 / 10000000 : K) < s K 10 ∧ u K 10 < 27182819 / 10000000 := by
  norm_num [u, s, sum_range_succ, Nat.factorial]

/-- The staircase at `n = 100` is about `2.7048`. -/
example : (2704 / 1000 : ℚ) < stair ℚ 100 ∧ stair ℚ 100 < 2705 / 1000 := by
  norm_num [stair]

/-! ### The limit: both sequences close on `Real.exp 1`

The core-Lean original stops here, because it has no real numbers. With mathlib the last step is
short: the staircase and the partial sums both tend to `Real.exp 1`, so the bounds above are
bounds on `e`. -/

open Filter Topology

/-- The staircase tends to `e`. This is mathlib's `Real.tendsto_one_add_div_pow_exp`. -/
theorem tendsto_stair : Tendsto (stair ℝ) atTop (𝓝 (Real.exp 1)) :=
  Real.tendsto_one_add_div_pow_exp 1

/-- The partial sums tend to `e`. This is mathlib's `NormedSpace.expSeries_div_hasSum_exp`. -/
theorem tendsto_s : Tendsto (s ℝ) atTop (𝓝 (Real.exp 1)) := by
  have h := (NormedSpace.expSeries_div_hasSum_exp (1 : ℝ)).tendsto_sum_nat
  rw [← Real.exp_eq_exp_ℝ] at h
  refine (h.comp (tendsto_add_atTop_nat 1)).congr fun N ↦ ?_
  simp [s]

/-- Every partial sum is below `e`. -/
theorem s_le_exp_one (N : ℕ) : s ℝ N ≤ Real.exp 1 :=
  (s_mono ℝ).ge_of_tendsto tendsto_s N

/-- Every upper partner is above `e`. -/
theorem exp_one_le_u (N : ℕ) (hN : 1 ≤ N) : Real.exp 1 ≤ u ℝ N :=
  le_of_tendsto tendsto_stair <| eventually_atTop.2 ⟨1, fun n hn ↦ stair_le_u ℝ n N hn hN⟩

/-- `e` is trapped between `2.7182818` and `2.7182819` by `s 10` and `u 10`. Mathlib has sharper
bounds: `Real.exp_one_gt_d9` and `Real.exp_one_lt_d9`. -/
theorem exp_one_bounds : 2.7182818 < Real.exp 1 ∧ Real.exp 1 < 2.7182819 := by
  have h := s_ten_u_ten ℝ
  constructor
  · exact lt_of_lt_of_le (by norm_num at h ⊢; exact h.1) (s_le_exp_one 10)
  · exact lt_of_le_of_lt (exp_one_le_u 10 (by norm_num)) (by norm_num at h ⊢; exact h.2)

end MultAdd.E
