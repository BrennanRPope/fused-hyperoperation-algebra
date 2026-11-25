/-
  FusedAlgebra.Topics.Calculus: calculus from mult-add (the exact, discrete part)
  ==============================================================================

  Core Lean 4.37, no mathlib, no `sorry`. Imports the shared base.

  The picture. A quantity that grows step by step, f(k+1) = f(k) ⊞ d(k), has
  mult-add slope d(k). Through the lift x ↦ 1 + x every step multiplies a box
  by (1 + d(k)). Then:

  * Fundamental theorem: the steps mult-added give back the whole climb
    (`ftc`, `mint_lift`): the box after n steps is the product of the steps.
  * Straight lines are exponentials: a constant mult-add slope a gives the
    mult-add power a^(⊞n) = (1 + a)^n − 1 (`mint_const`).
  * Product rule: in mult-add it is just the sum rule (`product_rule_step`,
    `mint_multAdd`). Its shadow in ordinary differences is Leibniz's rule
    with the corner square kept (`leibniz_corner`). Calculus is what is left
    when the corner is dropped.
  * Step size = corner weight: with steps of size h, relative rates combine
    by the weighted mult-add p ⊞ₕ q = p + q + h·p·q (`rate_product`). This is
    the ⊕ of time-scale calculus (Hilger 1988, Bohner–Peterson 2001), and
    h = 0 gives plain addition, the ordinary derivative.
  * Chain rule for blocks: m steps at a time, n times, is a power of a power
    (`mpow_mul`, `mint_add`): the rung above mult-add.
  * Newton's forward formula: f(n) = Σₖ layer(n,k) · Δᵏf(0), the n-th
    differences weighed by the cube's layers (`newton`). Shrinking the step
    turns the layers into xᵏ/k!, which is Taylor's series (not in Lean:
    it needs limits).
  * Mult-add input steps: x ⊞ h − x = h(1 + x), so stepping the input by
    mult-add gives Jackson's q-derivative with q = 1 + h (`jackson_gap`).

  Limits (h → 0) are not in core Lean; everything here is exact.

  LIBRARY NOTE: copied from calculus/lean/FusedCalculus.lean (Nov 25) as Topics/Calculus.
-/
import FusedAlgebra.Base
import FusedAlgebra.Cube

namespace MultAdd.Calculus

open MultAdd
open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Mult-add integral and the fundamental theorem (any commutative semiring) -/

section Semiring
variable {R : Type} [CommSemiring R]

/-- The mult-add integral: the steps d 0, d 1, …, d (n−1), mult-added. -/
def mint (d : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => multAdd (mint d n) (d n)

/-- The product of the lifted steps (1 + d 0)(1 + d 1)⋯(1 + d (n−1)). -/
def liftProd (d : Nat → R) : Nat → R
  | 0 => 1
  | n + 1 => liftProd d n * (1 + d n)

/-- Under the lift the mult-add integral is a product integral. -/
theorem mint_lift (d : Nat → R) (n : Nat) : 1 + mint d n = liftProd d n := by
  induction n with
  | zero => simp only [mint, liftProd]; grind
  | succ n ih => rw [mint, liftProd, multAdd_lift, ih]

/-- Fundamental theorem, part 1: the last step of the integral is the integrand. -/
theorem ftc_step (d : Nat → R) (n : Nat) : mint d (n + 1) = multAdd (mint d n) (d n) := rfl

/-- Fundamental theorem, part 2: a path whose mult-add steps are d ends at
    its start ⊞ the integral of d. -/
theorem ftc (f d : Nat → R) (h : ∀ k, f (k + 1) = multAdd (f k) (d k)) (n : Nat) :
    f n = multAdd (f 0) (mint d n) := by
  induction n with
  | zero => simp only [mint]; rw [multAdd_zero]
  | succ n ih => rw [h n, ih, mint, multAdd_assoc]

/-- The same in box form: the end box is the start box times every step. -/
theorem ftc_box (f d : Nat → R) (h : ∀ k, f (k + 1) = multAdd (f k) (d k)) (n : Nat) :
    1 + f n = (1 + f 0) * liftProd d n := by
  rw [ftc f d h n, multAdd_lift, mint_lift]

/-- Straight lines are exponentials: constant mult-add slope a gives a^(⊞n). -/
theorem mint_const (a : R) (n : Nat) : mint (fun _ => a) n = mpow a n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [mint, mpow, ih]

/-- The integral over two intervals joins by mult-add. -/
theorem mint_add (d : Nat → R) (a b : Nat) :
    mint d (a + b) = multAdd (mint d a) (mint (fun j => d (a + j)) b) := by
  induction b with
  | zero => simp only [Nat.add_zero, mint]; rw [multAdd_zero]
  | succ b ih => rw [← Nat.add_assoc, mint, mint, ih, multAdd_assoc]

/-- Product rule = sum rule: if f steps by d and g steps by e, then f ⊞ g
    steps by d ⊞ e. -/
theorem product_rule_step (f g d e : Nat → R)
    (hf : ∀ k, f (k + 1) = multAdd (f k) (d k))
    (hg : ∀ k, g (k + 1) = multAdd (g k) (e k)) (k : Nat) :
    multAdd (f (k + 1)) (g (k + 1)) = multAdd (multAdd (f k) (g k)) (multAdd (d k) (e k)) := by
  rw [hf, hg]; unfold multAdd; grind

/-- The integral of d ⊞ e is the integral of d ⊞ the integral of e. -/
theorem mint_multAdd (d e : Nat → R) (n : Nat) :
    mint (fun k => multAdd (d k) (e k)) n = multAdd (mint d n) (mint e n) := by
  induction n with
  | zero => simp only [mint]; rw [multAdd_zero]
  | succ n ih => rw [mint, mint, mint, ih]; unfold multAdd; grind

/-- Chain rule for blocks: n blocks of m equal steps is a power of a power. -/
theorem mpow_mul (a : R) (m n : Nat) : mpow a (m * n) = mpow (mpow a m) n := by
  induction n with
  | zero => rfl
  | succ n ih => rw [Nat.mul_succ, mpow_add, ih, mpow]

/-! ### Step size as corner weight (time-scale calculus) -/

/-- If F grows by the factor (1 + h·p) over a step of size h and G by
    (1 + h·q), then F·G grows by (1 + h·(p ⊞ₕ q)): relative rates combine by
    weighted mult-add with weight h. This is time-scale ⊕. -/
theorem rate_product (h p q F G F' G' : R)
    (hF : F' = F * (1 + h * p)) (hG : G' = G * (1 + h * q)) :
    F' * G' = (F * G) * (1 + h * wAdd h p q) := by
  rw [hF, hG, wAdd_lift]; grind

/-- Weight 0 (step size 0) is ordinary calculus: rates just add. -/
theorem rate_product_limit (p q : R) : wAdd 0 p q = p + q := wAdd_weight_zero p q

/-- The time-scale exponential after n steps: 1 + h·p^(⊞ₕn) = (1 + h·p)ⁿ. -/
theorem ts_exp (h p : R) (n : Nat) : 1 + h * wpow h p n = (1 + h * p) ^ n := wpow_lift h p n

/-- e_p · e_q = e_(p ⊞ₕ q): the exponential turns weighted mult-add into ×. -/
theorem ts_exp_mul (h p q : R) (n : Nat) :
    (1 + h * wpow h p n) * (1 + h * wpow h q n) = 1 + h * wpow h (wAdd h p q) n := by
  rw [ts_exp, ts_exp, ts_exp, ← mul_pow', wAdd_lift]

end Semiring

/-! ## 2. Ordinary differences seen through mult-add (commutative rings) -/

section Ring
variable {R : Type} [CommRing R]

/-- The ordinary forward difference. -/
def diff (f : Nat → R) (k : Nat) : R := f (k + 1) - f k

/-- Ordinary slope = height × mult-add slope (the log-derivative rule f' = (1 + f)·r). -/
theorem step_iff (f d : Nat → R) (k : Nat) :
    f (k + 1) = multAdd (f k) (d k) ↔ diff f k = (1 + f k) * d k := by
  unfold diff multAdd; constructor <;> intro h <;> grind

/-- Leibniz's rule with the corner kept: the box f·g grows by two strips
    and one corner square. -/
theorem leibniz_corner (f g : Nat → R) (k : Nat) :
    diff (fun j => f j * g j) k = diff f k * g k + f k * diff g k + diff f k * diff g k := by
  unfold diff; grind

/-- The ordinary difference of f ⊞ g: two lifted strips and the corner. -/
theorem diff_multAdd (f g : Nat → R) (k : Nat) :
    diff (fun j => multAdd (f j) (g j)) k
      = diff f k * (1 + g k) + (1 + f k) * diff g k + diff f k * diff g k := by
  unfold diff multAdd; grind

/-- An exponential's slope is a times itself, exactly. -/
theorem diff_pow (a : R) (n : Nat) : diff (fun j => (1 + a) ^ j) n = a * (1 + a) ^ n := by
  show (1 + a) ^ (n + 1) - (1 + a) ^ n = a * (1 + a) ^ n
  rw [Semiring.pow_succ]; grind

/-- Sum of the first n values. -/
def sumTo (g : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => sumTo g n + g n

/-- The additive fundamental theorem: differences telescope. -/
theorem ftc_add (f : Nat → R) (n : Nat) : sumTo (diff f) n = f n - f 0 := by
  induction n with
  | zero => simp only [sumTo]; grind
  | succ n ih => rw [sumTo, ih]; unfold diff; grind

/-- Stepping the input by mult-add moves it by h·(1 + x): the gap of
    Jackson's q-derivative with q = 1 + h, read in lifted coordinates. -/
theorem jackson_gap (x h : R) : multAdd x h - x = h * (1 + x) := by
  unfold multAdd; grind

end Ring

/-! ## 3. The way back (fields): mult-add slope and its inverse -/

section Field
variable {K : Type} [Field K]

/-- The mult-add slope of a step: (f(k+1) − f(k)) / (1 + f(k)). -/
def mslope (f : Nat → K) (k : Nat) : K := msub (f (k + 1)) (f k)

/-- Away from the wall, the mult-add slope is the step: f(k) ⊞ slope = f(k+1). -/
theorem mslope_step (f : Nat → K) (k : Nat) (h : 1 + f k ≠ 0) :
    multAdd (f k) (mslope f k) = f (k + 1) := by
  unfold mslope; rw [multAdd_comm]; exact msub_multAdd _ _ h

/-- Fundamental theorem with slopes read off: f(n) = f(0) ⊞ ∫ mslope f. -/
theorem ftc_mslope (f : Nat → K) (h : ∀ k, 1 + f k ≠ 0) (n : Nat) :
    f n = multAdd (f 0) (mint (mslope f) n) :=
  ftc f (mslope f) (fun k => (mslope_step f k (h k)).symm) n

/-- Time-scale inverse ⊖p = −p/(1 + h·p): exists exactly off the wall. -/
theorem ts_inverse (h p : K) (hp : 1 + h * p ≠ 0) : wAdd h p (-p / (1 + h * p)) = 0 := by
  unfold wAdd; grind

end Field

/-! ## 4. Newton's forward formula: the cube's layers weigh the differences -/

section Newton

/-- k-th forward difference on Int. -/
def dif : Nat → (Nat → Int) → Nat → Int
  | 0, f => f
  | k + 1, f => fun j => dif k f (j + 1) - dif k f j

/-- Shifting the input commutes with differencing. -/
theorem dif_shift (k : Nat) (f : Nat → Int) (j : Nat) :
    dif k (fun i => f (i + 1)) j = dif k f (j + 1) := by
  induction k generalizing j with
  | zero => rfl
  | succ k ih => simp only [dif]; rw [ih, ih]

/-- Each difference at 1 is the difference at 0 plus the next difference. -/
theorem dif_one (k : Nat) (f : Nat → Int) : dif k f 1 = dif k f 0 + dif (k + 1) f 0 := by
  simp only [dif]; grind

def isum (g : Nat → Int) : Nat → Int
  | 0 => 0
  | n + 1 => isum g n + g n

theorem isum_add (g h : Nat → Int) (m : Nat) :
    isum (fun k => g k + h k) m = isum g m + isum h m := by
  induction m with
  | zero => simp only [isum]; omega
  | succ m ih => simp only [isum] at *; rw [ih]; omega

theorem isum_congr (g h : Nat → Int) (m : Nat) (e : ∀ k, g k = h k) : isum g m = isum h m := by
  have : g = h := funext e
  rw [this]

/-- Pascal inside a weighted sum: the (n+1)-cube's layers split into two n-cubes. -/
theorem pascal_sum (n : Nat) (x : Nat → Int) (m : Nat) :
    isum (fun k => (Cube.layer (n + 1) k : Int) * x k) (m + 1)
      = isum (fun k => (Cube.layer n k : Int) * x k) (m + 1)
        + isum (fun k => (Cube.layer n k : Int) * x (k + 1)) m := by
  induction m with
  | zero => simp only [isum, Cube.layer_zero_right]; omega
  | succ m ih =>
    rw [isum, ih]
    simp only [isum, Cube.layer_succ_succ, Int.natCast_add]
    grind

/-- Newton's forward formula: f(n) = Σ_{k ≤ n} layer(n,k) · Δᵏf(0). -/
theorem newton (n : Nat) (f : Nat → Int) :
    f n = isum (fun k => (Cube.layer n k : Int) * dif k f 0) (n + 1) := by
  induction n generalizing f with
  | zero => simp only [isum, Cube.layer_zero_right, dif]; omega
  | succ n ih =>
    have h1 := ih (fun i => f (i + 1))
    rw [h1]
    rw [isum_congr _ (fun k => (Cube.layer n k : Int) * dif k f 0
          + (Cube.layer n k : Int) * dif (k + 1) f 0) (n + 1)
          (fun k => by rw [dif_shift, dif_one]; grind)]
    rw [isum_add, pascal_sum n (fun k => dif k f 0) (n + 1),
      isum.eq_2 (fun k => (Cube.layer n k : Int) * dif k f 0) (n + 1)]
    rw [Cube.layer_eq_zero n (n + 1) (by omega)]
    simp only [Int.natCast_zero, Int.zero_mul, Int.add_zero]

/-- A worked check: the cubes 0, 1, 8, 27, … have third difference 6 and
    fourth difference 0, so Newton rebuilds n³ from the numbers 0, 1, 6, 6. -/
theorem cubes_example : dif 3 (fun j => (j : Int) ^ 3) 0 = 6 ∧ dif 4 (fun j => (j : Int) ^ 3) 0 = 0 := by
  decide

end Newton

end MultAdd.Calculus
