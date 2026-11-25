/-
  Finance as the circle operation (called mult-add before 1.0)
  ===================

  Part of Fused Hyperoperation Algebra. Companion page and checks: ../ (econ-forecasting).

  STATUS: core Lean 4 only (no mathlib), no `sorry`. Run it with
  `lean FinanceMultAdd.lean`. Checked with Lean 4.37.0-rc1.

  The circle operation is a ◯ b = a + b + ab = (1+a)(1+b) − 1, and a ⊖ b = (a − b)/(1 + b) undoes it.

  Proved:
  * fisher:       (nominal ⊖ inflation) ◯ inflation = nominal. The real rate is the
                  circle-operation difference, not the plain difference.
  * mall_growth:  combining a list of yearly returns by the circle operation gives the total growth minus 1
                  (the n-cube of returns, all its corners together).
  * drag, drag_le: up s then down s around a mean m: (1+m+s)(1+m−s) = (1+m)² − s².
                  Volatility drag is the s×s square cut out of the corner.
  * gordon:       with q(1 + j) = 1, j·(q + q² + … + q^N) = 1 − q^N. A dividend growing
                  at g discounted at r is a flat perpetuity at j = r ⊖ g (q = (1+g)/(1+r)).
  * gsum_le:      for j > 0 the price never passes D/j.
  * gsum_pole:    for q ≥ 1 (growth at or above the discount rate) the N-year price is at
                  least N: the price has no ceiling. That is the pole at g = r.
  * am_hm, power_mean, average_factors: averaging two discount factors is never less than
                  discounting at the averaged growth factor, for every horizon n. This is the
                  geometric reason long-run discount rates fall when the rate is uncertain.
  * lowest_rate_wins: the averaged factor is at least half the factor of the lowest rate,
                  so over long horizons the lowest-rate corner sets the rate.

  LIBRARY NOTE: the circle operation and circle-operation subtraction (msub, Fisher's rule) come from FusedAlgebra.Base; wrapped in namespace FinanceCircle.
-/
import FusedAlgebra.Base

namespace FinanceCircle

open CircleOp

/-- Fisher's rule: real rate ◯ inflation = nominal rate (msub and the proof live in Base). -/
theorem fisher (a b : Rat) (h : 1 + b ≠ 0) : circ (msub a b) b = a := msub_circ a b h

theorem drag (m s : Rat) : (1 + m + s) * (1 + m - s) = (1 + m) ^ 2 - s ^ 2 := by grind


def gsum (q : Rat) : Nat → Rat
  | 0 => 0
  | n + 1 => gsum q n + q ^ (n + 1)

theorem gordon (q j : Rat) (h : q * (1 + j) = 1) : ∀ N, j * gsum q N = 1 - q ^ N
  | 0 => by simp [gsum]; exact Rat.sub_self.symm
  | N + 1 => by
    simp only [gsum]
    have ih := gordon q j h N
    rw [Rat.pow_succ]
    grind
theorem sq_nonneg' (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

theorem pow_nonneg' {x : Rat} (hx : 0 ≤ x) : ∀ n : Nat, 0 ≤ x ^ n
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by rw [Rat.pow_succ]; exact Rat.mul_nonneg (pow_nonneg' hx n) hx

theorem one_le_pow {q : Rat} (hq : 1 ≤ q) : ∀ n : Nat, 1 ≤ q ^ n
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by
    rw [Rat.pow_succ]
    have ih := one_le_pow hq n
    have h0 : 0 ≤ q ^ n := by grind
    have := Rat.mul_le_mul_of_nonneg_left hq h0
    grind

theorem pow_le_pow' {x y : Rat} (hy : 0 ≤ y) (hxy : y ≤ x) : ∀ n : Nat, y ^ n ≤ x ^ n
  | 0 => by grind
  | n + 1 => by
    rw [Rat.pow_succ, Rat.pow_succ]
    have ih := pow_le_pow' hy hxy n
    have a := Rat.mul_le_mul_of_nonneg_left hxy (pow_nonneg' hy n)
    have b := Rat.mul_le_mul_of_nonneg_left ih (by grind : (0:Rat) ≤ x)
    rw [Rat.mul_comm x, Rat.mul_comm x] at b
    grind

theorem cheb {x y : Rat} (hx : 0 ≤ x) (hy : 0 ≤ y) (n : Nat) : 0 ≤ (x ^ n - y ^ n) * (x - y) := by
  rcases (@Rat.le_total y x) with h | h
  · have := pow_le_pow' hy h n
    exact Rat.mul_nonneg (by grind) (by grind)
  · have := pow_le_pow' hx h n
    have := Rat.mul_nonneg (by grind : (0:Rat) ≤ y ^ n - x ^ n) (by grind : (0:Rat) ≤ y - x)
    grind

theorem power_mean {x y : Rat} (hx : 0 ≤ x) (hy : 0 ≤ y) : ∀ n : Nat, ((x + y) / 2) ^ n ≤ (x ^ n + y ^ n) / 2
  | 0 => by rw [Rat.pow_zero, Rat.pow_zero, Rat.pow_zero]; grind
  | n + 1 => by
    rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_succ]
    have ih := power_mean hx hy n
    have hA : 0 ≤ (x + y) / 2 := by grind
    have s1 := Rat.mul_le_mul_of_nonneg_left ih hA
    have c := cheb hx hy n
    grind

theorem drag_le (m s : Rat) : (1 + m + s) * (1 + m - s) ≤ (1 + m) ^ 2 := by
  rw [drag]; have h := sq_nonneg' s
  rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_zero, Rat.one_mul, Rat.pow_succ, Rat.pow_succ, Rat.pow_zero, Rat.one_mul]
  grind

theorem am_hm {u v x y : Rat} (hu : 0 < u) (hv : 0 < v) (hx : x * u = 1) (hy : y * v = 1) :
    ∀ c : Rat, c * ((u + v) / 2) = 1 → c ≤ (x + y) / 2 := by
  intro c hc
  have key : x * v + y * u - 2 = x * y * ((u - v) * (u - v)) := by grind
  have hx0 : 0 ≤ x := by
    rcases (@Rat.le_total 0 x) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt hu) (by grind : (0:Rat) ≤ -x); grind
  have hy0 : 0 ≤ y := by
    rcases (@Rat.le_total 0 y) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt hv) (by grind : (0:Rat) ≤ -y); grind
  have p := Rat.mul_nonneg (Rat.mul_nonneg hx0 hy0) (sq_nonneg' (u - v))
  -- (x + y)(u + v) = 2 + xv + yu >= 4, and c (u + v) = 2
  have e : (x + y) * (u + v) = x * u + y * v + (x * v + y * u) := by grind
  have p' : 0 ≤ x * v + y * u - 2 := by rw [key]; exact p
  have h4 : 4 ≤ (x + y) * (u + v) := by rw [e, hx, hy]; grind
  have hc0 : 0 ≤ c := by
    rcases (@Rat.le_total 0 c) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (by grind : (0:Rat) ≤ (u + v) / 2) (by grind : (0:Rat) ≤ -c); grind
  have h5 := Rat.mul_le_mul_of_nonneg_left h4 hc0
  have e2 : c * ((x + y) * (u + v)) = (x + y) * (2 * (c * ((u + v) / 2))) := by grind
  rw [e2, hc] at h5
  grind

theorem average_factors {u v x y d : Rat} (hu : 0 < u) (hv : 0 < v) (hx : x * u = 1) (hy : y * v = 1)
    (n : Nat) (hd : d * ((u + v) / 2) ^ n = 1) : d ≤ (x ^ n + y ^ n) / 2 := by
  have hpos : 0 < (u + v) / 2 := by grind
  let c := 1 / ((u + v) / 2)
  have hc : c * ((u + v) / 2) = 1 := Rat.div_mul_cancel (by grind)
  have hcn : c ^ n * ((u + v) / 2) ^ n = 1 := by rw [← mul_pow', hc]; exact one_pow' n
  have hdc : d = c ^ n := by
    have : d * ((u + v) / 2) ^ n * c ^ n = c ^ n * ((u + v) / 2) ^ n * d := by grind
    grind
  have hx0 : 0 ≤ x := by
    rcases (@Rat.le_total 0 x) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt hu) (by grind : (0:Rat) ≤ -x); grind
  have hy0 : 0 ≤ y := by
    rcases (@Rat.le_total 0 y) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt hv) (by grind : (0:Rat) ≤ -y); grind
  have h1 := am_hm hu hv hx hy c hc
  have hc0 : 0 ≤ c := by
    rcases (@Rat.le_total 0 c) with h | h
    · exact h
    · have := Rat.mul_le_mul_of_nonneg_left (Rat.le_of_lt hpos) (by grind : (0:Rat) ≤ -c); grind
  have h2 := pow_le_pow' hc0 h1 n
  have h3 := power_mean hx0 hy0 n
  grind

def mall : List Rat → Rat
  | [] => 0
  | r :: rs => circ r (mall rs)

def growth : List Rat → Rat
  | [] => 1
  | r :: rs => (1 + r) * growth rs

theorem mall_growth : ∀ rs : List Rat, mall rs = growth rs - 1
  | [] => by simp [mall, growth]; exact Rat.sub_self.symm
  | r :: rs => by simp only [mall, growth, circ]; rw [mall_growth rs]; grind

theorem gsum_le (q j : Rat) (h : q * (1 + j) = 1) (hq : 0 ≤ q) (N : Nat) : j * gsum q N ≤ 1 := by
  rw [gordon q j h N]; have := pow_nonneg' hq N; grind

theorem gsum_pole {q : Rat} (hq : 1 ≤ q) : ∀ N : Nat, (N : Rat) ≤ gsum q N
  | 0 => by simp [gsum]
  | N + 1 => by
    simp only [gsum]
    have := gsum_pole hq N
    have := one_le_pow hq (N + 1)
    have : ((N + 1 : Nat) : Rat) = (N : Rat) + 1 := by simp
    grind

theorem lowest_rate_wins {x y : Rat} (hy : 0 ≤ y) (hxy : y ≤ x) (n : Nat) :
    x ^ n / 2 ≤ (x ^ n + y ^ n) / 2 ∧ (x ^ n + y ^ n) / 2 ≤ x ^ n := by
  have := pow_nonneg' hy n
  have := pow_le_pow' hy hxy n
  constructor <;> grind


end FinanceCircle
