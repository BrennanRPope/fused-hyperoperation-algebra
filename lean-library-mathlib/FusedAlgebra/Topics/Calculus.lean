module

public import FusedAlgebra.Cube

/-!
# Calculus from mult-add: the exact, discrete part

A quantity that grows step by step, `f (k + 1) = f k ✳ d k`, has mult-add slope `d k`. Through
the lift `x ↦ 1 + x` every step multiplies a box by `1 + d k`.

This is the mathlib form of `FusedAlgebra.Topics.Calculus`. The mult-add integral `mint d n` is a
`Finset.prod` in the monoid `PreQuasiregular R`, the forward difference is mathlib's
`fwdDiff 1` (notation `Δ_[1]`), its iterates are `(Δ_[1])^[k]`, and the hand-made sums are
`Finset.sum`s.

## Main results

* `one_add_mint`: under the lift the mult-add integral is a product integral.
* `ftc`, `ftc_box`: the fundamental theorem, a path whose mult-add steps are `d` ends at its start
  mult-added to the integral of `d`.
* `mint_const`: straight lines are exponentials, `mint (fun _ ↦ a) n = mpow a n`.
* `product_rule_step`, `mint_multAdd`: in mult-add the product rule is the sum rule.
* `leibniz_corner`: its shadow in ordinary differences is Leibniz's rule with the corner kept.
* `rate_product`, `ts_exp_mul`: with steps of size `h`, relative rates combine by the weighted
  mult-add `wAdd h`; this is the `⊕` of time-scale calculus (Hilger 1988, Bohner–Peterson 2001).
* `mslope`, `ftc_mslope`: in a field the mult-add slope can be read off, away from the wall.
* `newton`: Newton's forward formula, `f n = ∑ k ≤ n, n.choose k * Δᵏ f 0`. This is mathlib's
  `shift_eq_sum_fwdDiff_iter`.
* `jackson_gap`: stepping the input by mult-add moves it by `h * (1 + x)`, the gap of Jackson's
  `q`-derivative with `q = 1 + h`.

Everything here is exact; no limit `h → 0` is taken.
-/

@[expose] public section

namespace MultAdd.Calculus

open MultAdd Finset PreQuasiregular fwdDiff

/-! ### The mult-add integral and the fundamental theorem -/

section CommSemiring

variable {R : Type*} [CommSemiring R]

/-- The mult-add integral: the steps `d 0, d 1, …, d (n - 1)`, mult-added. It is the product of
the steps in the monoid `PreQuasiregular R`. -/
def mint (d : ℕ → R) (n : ℕ) : R := (∏ k ∈ range n, equiv (d k)).val

@[simp]
theorem mint_zero (d : ℕ → R) : mint d 0 = 0 := rfl

/-- Fundamental theorem, part 1: the last step of the integral is the integrand. -/
theorem mint_succ (d : ℕ → R) (n : ℕ) : mint d (n + 1) = mint d n ✳ d n := by
  rw [mint, prod_range_succ, val_mul_eq_multAdd]; rfl

/-- Under the lift the mult-add integral is a product integral. -/
theorem one_add_mint (d : ℕ → R) (n : ℕ) : 1 + mint d n = ∏ k ∈ range n, (1 + d k) :=
  map_prod oneAddHom (fun k ↦ equiv (d k)) (range n)

/-- Fundamental theorem, part 2: a path whose mult-add steps are `d` ends at its start mult-added
to the integral of `d`. -/
theorem ftc (f d : ℕ → R) (h : ∀ k, f (k + 1) = f k ✳ d k) (n : ℕ) : f n = f 0 ✳ mint d n := by
  induction n with
  | zero => simp
  | succ n ih => rw [h n, ih, mint_succ, multAdd_assoc]

/-- The fundamental theorem in box form: the end box is the start box times every step. -/
theorem ftc_box (f d : ℕ → R) (h : ∀ k, f (k + 1) = f k ✳ d k) (n : ℕ) :
    1 + f n = (1 + f 0) * ∏ k ∈ range n, (1 + d k) := by
  rw [ftc f d h n, one_add_multAdd, one_add_mint]

/-- Straight lines are exponentials: a constant mult-add slope `a` gives `mpow a n`. -/
theorem mint_const (a : R) (n : ℕ) : mint (fun _ ↦ a) n = mpow a n := by
  rw [mint, prod_const, card_range]; rfl

/-- The integral over two intervals joins by mult-add. -/
theorem mint_add (d : ℕ → R) (a b : ℕ) :
    mint d (a + b) = mint d a ✳ mint (fun j ↦ d (a + j)) b := by
  rw [mint, prod_range_add, val_mul_eq_multAdd]; rfl

/-- Product rule = sum rule: if `f` steps by `d` and `g` steps by `e`, then `f ✳ g` steps by
`d ✳ e`. -/
theorem product_rule_step (f g d e : ℕ → R) (hf : ∀ k, f (k + 1) = f k ✳ d k)
    (hg : ∀ k, g (k + 1) = g k ✳ e k) (k : ℕ) :
    f (k + 1) ✳ g (k + 1) = f k ✳ g k ✳ (d k ✳ e k) := by
  rw [hf, hg, multAdd_assoc, multAdd_assoc, multAdd_left_comm (d k)]

/-- The integral of `d ✳ e` is the integral of `d` mult-added to the integral of `e`. -/
theorem mint_multAdd (d e : ℕ → R) (n : ℕ) :
    mint (fun k ↦ d k ✳ e k) n = mint d n ✳ mint e n := by
  simp only [mint, equiv_multAdd, prod_mul_distrib, val_mul_eq_multAdd]

/-! ### Step size as corner weight (time-scale calculus) -/

/-- If `F` grows by the factor `1 + h * p` over a step of size `h` and `G` by `1 + h * q`, then
`F * G` grows by `1 + h * wAdd h p q`: relative rates combine by weighted mult-add with weight
`h`. This is the time-scale `⊕`. -/
theorem rate_product (h p q F G : R) :
    F * (1 + h * p) * (G * (1 + h * q)) = F * G * (1 + h * wAdd h p q) := by
  rw [one_add_mul_wAdd]; ring

/-- `e_p * e_q = e_(wAdd h p q)`: the time-scale exponential `1 + h * wpow h p n = (1 + h * p) ^ n`
(`one_add_mul_wpow`) turns weighted mult-add into multiplication. -/
theorem ts_exp_mul (h p q : R) (n : ℕ) :
    (1 + h * wpow h p n) * (1 + h * wpow h q n) = 1 + h * wpow h (wAdd h p q) n := by
  simp only [one_add_mul_wpow, ← mul_pow, one_add_mul_wAdd]

end CommSemiring

/-! ### Ordinary differences seen through mult-add

The ordinary forward difference `f (k + 1) - f k` is mathlib's `Δ_[1] f k`. -/

section CommRing

variable {R : Type*} [CommRing R]

/-- Ordinary slope = height × mult-add slope (the log-derivative rule `f' = (1 + f) * r`). -/
theorem step_iff (f d : ℕ → R) (k : ℕ) :
    f (k + 1) = f k ✳ d k ↔ Δ_[1] f k = (1 + f k) * d k := by
  simp only [fwdDiff, multAdd]
  constructor <;> intro h <;> linear_combination h

/-- Leibniz's rule with the corner kept: the box `f * g` grows by two strips and one corner
square. -/
theorem leibniz_corner (f g : ℕ → R) (k : ℕ) :
    Δ_[1] (fun j ↦ f j * g j) k = Δ_[1] f k * g k + f k * Δ_[1] g k + Δ_[1] f k * Δ_[1] g k := by
  simp only [fwdDiff]; ring

/-- The ordinary difference of `f ✳ g`: two lifted strips and the corner. -/
theorem fwdDiff_multAdd (f g : ℕ → R) (k : ℕ) :
    Δ_[1] (fun j ↦ f j ✳ g j) k =
      Δ_[1] f k * (1 + g k) + (1 + f k) * Δ_[1] g k + Δ_[1] f k * Δ_[1] g k := by
  simp only [fwdDiff, multAdd]; ring

/-- An exponential's slope is `a` times itself, exactly. -/
theorem fwdDiff_one_add_pow (a : R) (n : ℕ) :
    Δ_[1] (fun j ↦ (1 + a) ^ j) n = a * (1 + a) ^ n := by
  simp only [fwdDiff, pow_succ]; ring

/-- Stepping the input by mult-add moves it by `h * (1 + x)`: the gap of Jackson's `q`-derivative
with `q = 1 + h`, read in lifted coordinates. -/
theorem jackson_gap (x h : R) : x ✳ h - x = h * (1 + x) := by
  simp only [multAdd]; ring

end CommRing

/-- The additive fundamental theorem: differences telescope. This is `Finset.sum_range_sub`. -/
theorem ftc_add {G : Type*} [AddCommGroup G] (f : ℕ → G) (n : ℕ) :
    ∑ k ∈ range n, Δ_[1] f k = f n - f 0 :=
  sum_range_sub f n

/-! ### The way back: mult-add slope and its inverse -/

section Field

variable {K : Type*} [Field K]

/-- The mult-add slope of a step: `(f (k + 1) - f k) / (1 + f k)`. -/
def mslope (f : ℕ → K) (k : ℕ) : K := msub (f (k + 1)) (f k)

/-- Away from the wall, the mult-add slope is the step: `f k ✳ mslope f k = f (k + 1)`. -/
theorem mslope_step (f : ℕ → K) (k : ℕ) (h : 1 + f k ≠ 0) : f k ✳ mslope f k = f (k + 1) := by
  rw [mslope, multAdd_comm, msub_multAdd _ _ h]

/-- The fundamental theorem with slopes read off: `f n = f 0 ✳ mint (mslope f) n`. -/
theorem ftc_mslope (f : ℕ → K) (h : ∀ k, 1 + f k ≠ 0) (n : ℕ) :
    f n = f 0 ✳ mint (mslope f) n :=
  ftc f (mslope f) (fun k ↦ (mslope_step f k (h k)).symm) n

/-- The time-scale inverse `-p / (1 + h * p)` exists off the wall. -/
theorem ts_inverse (h p : K) (hp : 1 + h * p ≠ 0) : wAdd h p (-p / (1 + h * p)) = 0 := by
  rw [wAdd, div_eq_mul_inv]
  linear_combination (-p) * mul_inv_cancel₀ hp

end Field

/-! ### Newton's forward formula: the cube's layers weigh the differences

The `k`-th forward difference is `(Δ_[1])^[k] f`. -/

section Newton

variable {R : Type*} [Ring R]

/-- Pascal inside a weighted sum: the layers of the `(n + 1)`-cube split into two `n`-cubes. -/
theorem pascal_sum (n : ℕ) (x : ℕ → R) (m : ℕ) :
    ∑ k ∈ range (m + 1), ((n + 1).choose k : R) * x k =
      ∑ k ∈ range (m + 1), (n.choose k : R) * x k + ∑ k ∈ range m, (n.choose k : R) * x (k + 1) := by
  rw [sum_range_succ', sum_range_succ' _ m]
  simp only [Nat.choose_succ_succ, Nat.cast_add, add_mul, sum_add_distrib, Nat.choose_zero_right]
  abel

/-- Newton's forward formula: `f n = ∑ k ≤ n, n.choose k * Δᵏ f 0`, the `k`-th differences
weighed by the cube's layers. This is mathlib's Gregory–Newton formula
`shift_eq_sum_fwdDiff_iter`. -/
theorem newton (n : ℕ) (f : ℕ → R) :
    f n = ∑ k ∈ range (n + 1), (n.choose k : R) * (Δ_[1])^[k] f 0 := by
  simpa [nsmul_eq_mul] using shift_eq_sum_fwdDiff_iter 1 f n 0

/-- A worked check: the cubes `0, 1, 8, 27, …` have third difference `6` and fourth difference
`0`, so Newton rebuilds `n ^ 3` from the numbers `0, 1, 6, 6`. -/
theorem cubes_example :
    (Δ_[1])^[3] (fun j : ℕ ↦ (j : ℤ) ^ 3) 0 = 6 ∧ (Δ_[1])^[4] (fun j : ℕ ↦ (j : ℤ) ^ 3) 0 = 0 := by
  simp [fwdDiff]

end Newton

end MultAdd.Calculus
