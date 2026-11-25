module

public import FusedAlgebra.Basic

/-!
# π from mult-add

Mult-add over the Gaussian integers. The turning law is slope addition
`x ⊕ y = (x + y) / (1 - x * y)`; keeping a slope as a pair `(num, den)`, the lens
`x ↦ den + i * num` turns slope addition into Gaussian multiplication, so into mult-add after the
lift. Every real slope lands on the circle of radius `1` around the wall `-1`; that circle meets
the real line in the two points `0` (the identity) and `-2` (the parity point, slope `∞`), and π
is how far the turning flow runs from one to the other.

This is the mathlib form of `FusedAlgebra.Topics.Pi`. The hand-made Gaussian integers of the
original are mathlib's `GaussianInt` (`ℤ[i] = Zsqrtd (-1)`), Gaussian mult-add is `MultAdd.multAdd`
on that ring, and its laws (lift, commutativity, associativity, identity, multiplicativity of the
norm) are the general ones of `FusedAlgebra.Basic` and `Zsqrtd`.

## Main results

* `circle_closed`: the circle `(1 + a).norm = 1` is closed under mult-add.
* `one_add_mul_self_eq_one_iff`, `multAdd_self_eq_zero_iff`: on a line without zero divisors the
  circle is exactly `{0, -2}`.
* `lens`, `lens_conj`: slope addition is Gaussian multiplication.
* `quarter_twice`, `quarter_four`: the quarter step `-1 + i`.
* `machin`, `euler`, `hermann`, `hutton`: the two-term arctangent identities as exact Gaussian
  products landing on the 45° ray.
* `cross_10`, `cross_100`, `cross_1000`: the number of mult-add steps of `i / n` needed to pass
  the far side is `32`, `315`, `3142`.
* `greg_even_le_odd`: Gregory's partial sums nest, for `0 ≤ x ≤ 1`.
* `machin_nested`, `machin_50_digits`: the Machin intervals nest, and the 22nd one fixes 50
  decimals.

All of this is classical (Brahmagupta–Fibonacci, Gregory 1671, Machin 1706, Størmer 1899). The
analytic facts that the original could not state are in mathlib: Gregory's series is
`Real.hasSum_arctan` (restated here as `hasSum_gregTerm`), Machin's formula is
`Real.four_mul_arctan_inv_5_sub_arctan_inv_239`, and Leibniz's series is
`Real.tendsto_sum_pi_div_four`. That the nested Machin intervals of this file close on `Real.pi`
is not proved here.
-/

@[expose] public section

namespace PiMultAdd

open MultAdd Finset

/-! ### The circle around the wall -/

/-- The circle `|1 + a| = 1` (centre the wall `-1`, radius `1`) is closed under mult-add,
because norms multiply (`Zsqrtd.norm_mul`, the Brahmagupta–Fibonacci identity). -/
theorem circle_closed (a b : GaussianInt) (ha : (1 + a).norm = 1) (hb : (1 + b).norm = 1) :
    (1 + (a ✳ b)).norm = 1 := by
  rw [one_add_multAdd, Zsqrtd.norm_mul, ha, hb, one_mul]

section Line

variable {R : Type*} [Ring R] [NoZeroDivisors R]

/-- On a line without zero divisors, the circle `(1 + x) ^ 2 = 1` is exactly the two points `0`
and `-2`. -/
theorem one_add_mul_self_eq_one_iff (x : R) : (1 + x) * (1 + x) = 1 ↔ x = 0 ∨ x = -2 := by
  rw [mul_self_eq_one_iff, add_eq_left, ← eq_sub_iff_add_eq']
  norm_num

/-- In mult-add units: `x ✳ x = 0` (a half turn twice is a full turn) exactly at `0` and `-2`. -/
theorem multAdd_self_eq_zero_iff (x : R) : x ✳ x = 0 ↔ x = 0 ∨ x = -2 := by
  rw [← one_add_mul_self_eq_one_iff, ← one_add_multAdd, add_eq_left]

end Line

/-- The parity point is its own way back: `-2 ✳ -2 = 0`. -/
theorem parity_twice {R : Type*} [Ring R] : (-2 : R) ✳ -2 = 0 := by
  rw [multAdd]; noncomm_ring

/-! ### The lens: the turning law is mult-add

A slope `num / den` is kept as a pair, so `∞ = (n, 0)` is an ordinary pair. Slope addition
`(x + y) / (1 - x * y)` on pairs is `(p, q) ⊕ (r, s) = (p * s + r * q, q * s - p * r)`. -/

/-- A slope `num / den`, kept as a pair of integers so that `∞ = (n, 0)` is allowed. -/
@[ext]
structure Slope where
  /-- The numerator (the rise). -/
  num : ℤ
  /-- The denominator (the run). -/
  den : ℤ
deriving DecidableEq, Repr

/-- The turning law: slope addition `(x + y) / (1 - x * y)` on pairs. -/
def tl (x y : Slope) : Slope := ⟨x.num * y.den + y.num * x.den, x.den * y.den - x.num * y.num⟩

/-- The Gaussian integer `den + i * num` that a slope stands for. -/
def toG (x : Slope) : GaussianInt := ⟨x.den, x.num⟩

/-- The lens: slope addition is exactly Gaussian multiplication. -/
theorem lens (x y : Slope) : toG (tl x y) = toG x * toG y := by
  ext <;> simp [toG, tl] <;> ring

/-- The conjugates multiply the same way, so the lift `W / W̄` of a sum is the product of the
lifts. -/
theorem lens_conj (x y : Slope) : star (toG (tl x y)) = star (toG x) * star (toG y) := by
  rw [lens, star_mul']

/-- Every real slope lands on the circle: `W` and `W̄` have the same length. -/
theorem lens_on_circle (x : Slope) : (toG x).norm = (star (toG x)).norm :=
  (Zsqrtd.norm_conj _).symm

/-- Slope `0` lands on the identity: `W = W̄`. -/
theorem lens_zero (d : ℤ) : toG ⟨0, d⟩ = star (toG ⟨0, d⟩) := by
  ext <;> simp [toG]

/-- Slope `∞` lands on the parity point: `W = -W̄`, lift `-1`, mult-add value `-2`. -/
theorem lens_infinity (n : ℤ) : toG ⟨n, 0⟩ = -star (toG ⟨n, 0⟩) := by
  ext <;> simp [toG]

/-! ### The quarter step -/

/-- `-1 + i`, the lens image of slope `1` (a quarter of the way round). -/
def q4 : GaussianInt := ⟨-1, 1⟩

theorem quarter_twice : q4 ✳ q4 = -2 := by decide

theorem quarter_four : mpow q4 4 = 0 := by
  rw [mpow_eq]; decide

theorem quarter_on_circle : (1 + q4).norm = 1 := by decide

/-- The same in slopes: `1 ⊕ 1 = 2 / 0 = ∞`. -/
theorem slope_one_twice : tl ⟨1, 1⟩ ⟨1, 1⟩ = ⟨2, 0⟩ := by decide

/-- Four steps of slope `1` come back to slope `0`. -/
theorem slope_one_four : tl (tl (tl ⟨1, 1⟩ ⟨1, 1⟩) ⟨1, 1⟩) ⟨1, 1⟩ = ⟨0, -4⟩ := by decide

/-! ### Steps that land exactly on the eighth turn

Mult-add steps of `i / m` are, after scaling by `m`, the Gaussian integers `m + i`. A mix of
steps lands on the 45° ray when the product is `k * (1 + i)` with `k > 0`. -/

/-- Machin (1706): four steps of `i / 5`, one step back of `i / 239`. The analytic statement is
`Real.four_mul_arctan_inv_5_sub_arctan_inv_239`. -/
theorem machin : (⟨5, 1⟩ : GaussianInt) ^ 4 * ⟨239, -1⟩ = ⟨114244, 114244⟩ := by decide

/-- Euler: one step of `i / 2`, one of `i / 3`. -/
theorem euler : (⟨2, 1⟩ : GaussianInt) * ⟨3, 1⟩ = ⟨5, 5⟩ := by decide

/-- Hermann: two steps of `i / 2`, one back of `i / 7`. -/
theorem hermann : (⟨2, 1⟩ : GaussianInt) ^ 2 * ⟨7, -1⟩ = ⟨25, 25⟩ := by decide

/-- Hutton: two steps of `i / 3`, one of `i / 7`. -/
theorem hutton : (⟨3, 1⟩ : GaussianInt) ^ 2 * ⟨7, 1⟩ = ⟨50, 50⟩ := by decide

/-- Machin again, in the turning law itself: `1/5 ⊕ 1/5 ⊕ 1/5 ⊕ 1/5 ⊕ (-1/239) = 1`. -/
theorem machin_slopes :
    tl (tl (tl (tl ⟨1, 5⟩ ⟨1, 5⟩) ⟨1, 5⟩) ⟨1, 5⟩) ⟨-1, 239⟩ = ⟨114244, 114244⟩ := by decide

/-! ### Counting tiny turning steps until we pass the far side

`(1 + i / n) ^ m = (n + i) ^ m / n ^ m`, so the sign of the imaginary part of `(n + i) ^ m` tells
which side of the real line `m` mult-add steps of `i / n` have reached. -/

/-- `firstCross n fuel m z` multiplies `z` by `n + i` up to `fuel` times, and returns `m` plus the
number of steps at which the imaginary part first becomes negative (`0` if it never does). -/
def firstCross (n : ℤ) : ℕ → ℕ → GaussianInt → ℕ
  | 0, _, _ => 0
  | fuel + 1, m, z =>
    let z' := z * ⟨n, 1⟩
    if z'.im < 0 then m + 1 else firstCross n fuel (m + 1) z'

/-- 32 steps of `i / 10` pass `-2`; 31 do not. -/
theorem cross_10 : firstCross 10 100 0 1 = 32 := by decide +kernel

/-- 315 steps of `i / 100`. -/
theorem cross_100 : firstCross 100 1000 0 1 = 315 := by decide +kernel

/-- 3142 steps of `i / 1000`. -/
theorem cross_1000 : firstCross 1000 4000 0 1 = 3142 := by decide +kernel

/-! ### Gregory's sums and nested intervals

The way back from mult-add to addition is `L a = a - a ^ 2 / 2 + a ^ 3 / 3 - ⋯`. At `a = i * x`
its imaginary part is Gregory's series `x - x ^ 3 / 3 + x ^ 5 / 5 - ⋯`. -/

section Gregory

variable {K : Type*} [Field K]

/-- Gregory's `k`-th term `x ^ (2 * k + 1) / (2 * k + 1)`. -/
def gregTerm (x : K) (k : ℕ) : K := x ^ (2 * k + 1) / ((2 * k + 1 : ℕ) : K)

/-- Gregory's partial sums `∑ k < n, (-1) ^ k * x ^ (2 * k + 1) / (2 * k + 1)`. -/
def greg (x : K) (n : ℕ) : K := ∑ k ∈ range n, (-1) ^ k * gregTerm x k

theorem greg_succ (x : K) (n : ℕ) : greg x (n + 1) = greg x n + (-1) ^ n * gregTerm x n :=
  sum_range_succ _ _

theorem greg_two (x : K) (m : ℕ) :
    greg x (2 * m + 2) = greg x (2 * m) + gregTerm x (2 * m) - gregTerm x (2 * m + 1) := by
  rw [greg_succ, greg_succ, (even_two_mul m).neg_one_pow, (odd_two_mul_add_one m).neg_one_pow]
  ring

/-- Gregory's series sums to the arctangent: this is mathlib's `Real.hasSum_arctan`. -/
theorem hasSum_gregTerm {x : ℝ} (hx : ‖x‖ < 1) :
    HasSum (fun k ↦ (-1) ^ k * gregTerm x k) (Real.arctan x) := by
  simpa only [gregTerm, mul_div_assoc] using Real.hasSum_arctan hx

variable [LinearOrder K] [IsStrictOrderedRing K]

theorem gregTerm_nonneg {x : K} (h0 : 0 ≤ x) (k : ℕ) : 0 ≤ gregTerm x k := by
  unfold gregTerm; positivity

/-- The terms shrink for `0 ≤ x ≤ 1`. -/
theorem gregTerm_succ_le {x : K} (h0 : 0 ≤ x) (h1 : x ≤ 1) (k : ℕ) :
    gregTerm x (k + 1) ≤ gregTerm x k := by
  unfold gregTerm
  apply div_le_div₀ (by positivity) (pow_le_pow_of_le_one h0 h1 (by lia)) (by positivity)
  exact_mod_cast (by lia : 2 * k + 1 ≤ 2 * (k + 1) + 1)

/-- Even partial sums climb. -/
theorem greg_even_up {x : K} (h0 : 0 ≤ x) (h1 : x ≤ 1) (m : ℕ) :
    greg x (2 * m) ≤ greg x (2 * m + 2) := by
  rw [greg_two]; linarith [gregTerm_succ_le h0 h1 (2 * m)]

/-- Each odd partial sum sits above the even one before it. -/
theorem greg_odd_above {x : K} (h0 : 0 ≤ x) (m : ℕ) : greg x (2 * m) ≤ greg x (2 * m + 1) := by
  rw [greg_succ, (even_two_mul m).neg_one_pow]; linarith [gregTerm_nonneg h0 (2 * m)]

/-- Odd partial sums fall. -/
theorem greg_odd_down {x : K} (h0 : 0 ≤ x) (h1 : x ≤ 1) (m : ℕ) :
    greg x (2 * m + 3) ≤ greg x (2 * m + 1) := by
  have h : (-1 : K) ^ (2 * m + 2) = 1 := (even_two_mul (m + 1)).neg_one_pow
  rw [greg_succ (n := 2 * m + 2), greg_two, greg_succ (n := 2 * m), h,
    (even_two_mul m).neg_one_pow]
  linarith [gregTerm_succ_le h0 h1 (2 * m + 1)]

/-- `[even sum, odd sum]` is a nested interval: every even sum sits below every odd sum. -/
theorem greg_even_le_odd {x : K} (h0 : 0 ≤ x) (h1 : x ≤ 1) (m n : ℕ) :
    greg x (2 * m) ≤ greg x (2 * n + 1) := by
  have up : Monotone fun a ↦ greg x (2 * a) := monotone_nat_of_le_succ (greg_even_up h0 h1)
  have down : Antitone fun a ↦ greg x (2 * a + 1) := antitone_nat_of_succ_le (greg_odd_down h0 h1)
  rcases le_total m n with h | h
  · exact (up h).trans (greg_odd_above h0 n)
  · exact (greg_odd_above h0 m).trans (down h)

end Gregory

/-! ### The Machin interval

Lower end: `16 * (even sum at 1/5) - 4 * (odd sum at 1/239)`; upper end: the other way round. -/

/-- The lower end of the `N`-th Machin interval. -/
def machinLo (N : ℕ) : ℚ := 16 * greg (1 / 5) (2 * N) - 4 * greg (1 / 239) (2 * N + 1)

/-- The upper end of the `N`-th Machin interval. -/
def machinHi (N : ℕ) : ℚ := 16 * greg (1 / 5) (2 * N + 1) - 4 * greg (1 / 239) (2 * N)

theorem machinLo_le_machinHi (N : ℕ) : machinLo N ≤ machinHi N := by
  unfold machinLo machinHi
  linarith [greg_odd_above (x := (1 / 5 : ℚ)) (by norm_num) N,
    greg_odd_above (x := (1 / 239 : ℚ)) (by norm_num) N]

/-- The Machin intervals nest: the lower end climbs and the upper end falls with `N`. -/
theorem machin_nested (N : ℕ) :
    machinLo N ≤ machinLo (N + 1) ∧ machinHi (N + 1) ≤ machinHi N := by
  unfold machinLo machinHi
  have a := greg_even_up (x := (1 / 5 : ℚ)) (by norm_num) (by norm_num) N
  have b := greg_odd_down (x := (1 / 239 : ℚ)) (by norm_num) (by norm_num) N
  have c := greg_odd_down (x := (1 / 5 : ℚ)) (by norm_num) (by norm_num) N
  have d := greg_even_up (x := (1 / 239 : ℚ)) (by norm_num) (by norm_num) N
  rw [show 2 * (N + 1) = 2 * N + 2 by ring]
  constructor <;> linarith

/-- 50 decimals: `3.14159265358979323846264338327950288419716939937510 < lo ≤ hi < …511`. -/
theorem machin_50_digits :
    (314159265358979323846264338327950288419716939937510 / 10 ^ 50 : ℚ) < machinLo 22 ∧
      machinHi 22 < (314159265358979323846264338327950288419716939937511 / 10 ^ 50 : ℚ) := by
  decide +kernel

end PiMultAdd
