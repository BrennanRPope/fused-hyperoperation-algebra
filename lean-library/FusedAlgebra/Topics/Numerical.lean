/-
  Numerical analysis through mult-add
  ===================================

  Part of Fused Hyperoperation Algebra. Page and float64 checks: ../ (numerical-analysis).
  Core Lean 4 only (no mathlib), no `sorry`. Checked with Lean 4.37.0-rc1.

  Rounding is not real-number arithmetic, so floats themselves live in ../checks.py.
  What Lean proves here is the exact algebra underneath:
  * relative errors compose by mult-add (products) and mult-add subtraction (quotients);
  * Higham's θₙ is the mult-add fold of the step errors, bounded by (1+u)ⁿ − 1 ≤ nu/(1 − nu);
  * the lift's way back amplifies a rounding by (1 + a)/a; the inverse is sensitive only at the wall;
  * Sterbenz's window: inside the octave 1/2 ≤ 1 + a ≤ 2 the inverse is within a factor 2 of −a;
  * Newton's division-free inverse is v ↦ v ◯ (−(u ◯ v)), and the residual squares: e ↦ −e²;
  * Goldschmidt's factors fold to −e^(2ⁿ);
  * a chain of mult-adds is a chain of Horner (fma) steps, and Bernoulli: a^(◯n) ≥ n·a.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/Numerical.lean; merged unchanged. circ_bound states the same as ComplexAnalysis.corner_bound but is kept,
  because a later proof here uses it and this file does not import ComplexAnalysis; sq_nonneg_rat
  and pow_nonneg_rat stay as local Rat helpers like the ones in Finance and Thermodynamics.
-/
import FusedAlgebra.Base

namespace Numerical

open CircleOp
open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 1. The three fused shapes (what one fma computes) -/

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

theorem fma_shape_corner (a b : R) : a ◯ b = a * b + (a + b) := by unfold circ; grind
theorem fma_shape_right (a b : R) : a ◯ b = a * (1 + b) + b := by unfold circ; grind
theorem fma_shape_left (a b : R) : a ◯ b = b * (1 + a) + a := by unfold circ; grind

/-! ## 2. Relative errors are mult-add numbers -/

/-- Multiply two perturbed values: the relative errors mult-add. -/
theorem relerr_mul (x y d₁ d₂ : R) :
    (x * (1 + d₁)) * (y * (1 + d₂)) = x * y * (1 + (d₁ ◯ d₂)) := by
  have h := circ_lift d₁ d₂; rw [h]; grind

/-- Three perturbed factors: the errors fold. -/
theorem relerr_mul3 (x y z d₁ d₂ d₃ : R) :
    (x * (1 + d₁)) * (y * (1 + d₂)) * (z * (1 + d₃)) = x * y * z * (1 + (d₁ ◯ d₂ ◯ d₃)) := by
  have h := circ_lift (d₁ ◯ d₂) d₃; have h' := circ_lift d₁ d₂; rw [h, h']; grind

/-- Higham's θₙ: the product of the factors (1 + δᵢ), minus 1, is the mult-add fold. -/
theorem theta_is_fold (ds : List R) : prodL (ds.map (fun d => 1 + d)) - 1 = foldCirc ds := by
  have h := foldCirc_lift ds; grind

/-- A chain of rounded products: each step's error δ multiplies in. -/
def chain : List R → R → R
  | [], p => p
  | d :: ds, p => chain ds (p * (1 + d))

theorem chain_eq (ds : List R) (p : R) : chain ds p = p * (1 + foldCirc ds) := by
  induction ds generalizing p with
  | nil => show p = p * (1 + 0); grind
  | cons d ds ih =>
    show chain ds (p * (1 + d)) = p * (1 + (d ◯ foldCirc ds))
    rw [ih, circ_lift]; grind

/-- Perturbing a value moves a ◯ b by exactly h·(1 + b): the sensitivity is the other lift. -/
theorem circ_sensitivity (a b h : R) : ((a + h) ◯ b) - (a ◯ b) = h * (1 + b) := by
  unfold circ; grind

/-- Newton's step for the inverse, without division. -/
def newton (u v : R) : R := v ◯ (-(u ◯ v))

/-- The residual squares (with a sign): the corner e·e is all that is left. -/
theorem newton_residual (u v : R) : u ◯ newton u v = -((u ◯ v) * (u ◯ v)) := by
  unfold newton circ; grind

/-- n Newton steps. -/
def newtonIter (u : R) : Nat → R → R
  | 0, v => v
  | n + 1, v => newtonIter u n (newton u v)

/-- After n steps the residual is −(−e₀)^(2ⁿ): quadratic convergence, digits double. -/
theorem newton_iter_residual (u : R) : ∀ (n : Nat) (v : R),
    u ◯ newtonIter u n v = -(pw (-(u ◯ v)) (2 ^ n))
  | 0, v => by simp only [newtonIter, Nat.pow_zero, pw]; grind
  | n + 1, v => by
    simp only [newtonIter]
    rw [newton_iter_residual u n (newton u v), newton_residual, Nat.pow_succ, pw_mul]
    generalize u ◯ v = e
    have h1 : -(-(e * e)) = (-e) * (-e) := by grind
    rw [h1, pw_mul_base]
    simp only [pw]; grind

/-- At the exact inverse the residual is 0 and Newton stands still. -/
theorem newton_fixed (u v : R) (h : u ◯ v = 0) : newton u v = v := by
  unfold newton; rw [h]; unfold circ; grind

/-- Goldschmidt's factors: (−e) ◯ e ◯ e² ◯ e⁴ ◯ … -/
def gold (e : R) : Nat → R
  | 0 => -e
  | n + 1 => gold e n ◯ pw e (2 ^ n)

theorem gold_eq (e : R) : ∀ n, gold e n = -(pw e (2 ^ n))
  | 0 => by simp only [gold, Nat.pow_zero, pw]; grind
  | n + 1 => by
    simp only [gold]; rw [gold_eq e n, Nat.pow_succ, pw_mul]
    generalize pw e (2 ^ n) = y
    unfold circ; simp only [pw]; grind

/-- In the lift: (1 − e)(1 + e)(1 + e²)…(1 + e^(2^(n−1))) = 1 − e^(2ⁿ). -/
theorem gold_lift (e : R) (n : Nat) : 1 + gold e n = 1 - pw e (2 ^ n) := by
  rw [gold_eq]; grind

/-! ## 3. Horner: a chain of mult-adds is a chain of fma steps -/

/-- Left-to-right fold with one Horner/fma step per entry: acc ↦ acc·(1 + b) + b. -/
def hornerCirc : R → List R → R
  | acc, [] => acc
  | acc, b :: bs => hornerCirc (acc * (1 + b) + b) bs

theorem hornerCirc_eq (acc : R) (bs : List R) : hornerCirc acc bs = acc ◯ foldCirc bs := by
  induction bs generalizing acc with
  | nil => simp only [hornerCirc, foldCirc]; rw [circ_zero]
  | cons b bs ih =>
    simp only [hornerCirc, foldCirc]; rw [ih]; unfold circ; grind

/-- Starting from 0, the n Horner steps give the whole cube: every corner product once. -/
theorem horner_cube (bs : List R) : 1 + hornerCirc 0 bs = sumL ((groups bs).map prodL) := by
  rw [hornerCirc_eq, zero_circ, one_add_foldCirc_eq_corners]

/-- Ordinary Horner evaluation, coefficients from the top. -/
def horner (x : R) : List R → R
  | [] => 0
  | c :: cs => c + x * horner x cs

/-- Mult-adding the same b n times is Horner on the polynomial b(1 + y + … + yⁿ⁻¹), y = 1 + b. -/
theorem mpow_horner (b : R) : ∀ n, mpow b n = horner (1 + b) (List.replicate n b)
  | 0 => rfl
  | n + 1 => by
    simp only [mpow, List.replicate, horner]
    rw [← mpow_horner b n, fma_shape_left]; grind

/-- Powers by squaring: a^(◯2n) = a^(◯n) ◯ a^(◯n). -/
theorem mpow_double (a : R) (n : Nat) : mpow a (2 * n) = mpow a n ◯ mpow a n := by
  rw [Nat.two_mul, mpow_add]

end Ring

/-! ## 4. Quotients, conditioning and Sterbenz's window (fields, ordered Rat) -/

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- Divide two perturbed values: the relative errors mult-subtract (Fisher's rule). -/
theorem relerr_div (x y d₁ d₂ : K) (hy : y ≠ 0) (h₂ : 1 + d₂ ≠ 0) :
    (x * (1 + d₁)) / (y * (1 + d₂)) = (x / y) * (1 + msub d₁ d₂) := by
  unfold msub; grind

/-- The way back from the lift: a relative rounding η of 1 + a becomes η·(1 + a)/a in a. -/
theorem lift_amplifier (a η : K) (ha : a ≠ 0) :
    ((1 + a) * (1 + η) - 1 - a) / a = η * ((1 + a) / a) := by
  have : a / a = 1 := by grind
  grind

/-- The inverse's exact sensitivity: only 1 + u near 0 (the wall) makes it large. -/
theorem inverse_sensitivity (u h : K) (h1 : 1 + u ≠ 0) (h2 : 1 + u + h ≠ 0) :
    circInv (u + h) - circInv u = -h / ((1 + u) * (1 + u + h)) := by
  unfold circInv; grind

/-- The fused inverse undoes the lift's reciprocal: 1 + inv u = 1/(1 + u). -/
theorem inverse_lift (u : K) (h : 1 + u ≠ 0) : (1 + circInv u) * (1 + u) = 1 := by
  unfold circInv; grind

end Field

theorem sq_nonneg_rat (s : Rat) : 0 ≤ s * s := by
  rcases (Rat.le_total : (0:Rat) ≤ s ∨ s ≤ 0) with h | h
  · exact Rat.mul_nonneg h h
  · have := Rat.mul_nonneg (by grind : (0:Rat) ≤ -s) (by grind : (0:Rat) ≤ -s); grind

theorem nonneg_of_mul_pos {x y : Rat} (hy : 0 < y) (h : 0 ≤ x * y) : 0 ≤ x := by
  rcases (Rat.le_total : (0:Rat) ≤ x ∨ x ≤ 0) with hx | hx
  · exact hx
  · rcases Rat.le_iff_lt_or_eq.mp hx with hlt | heq
    · have : 0 < (-x) * y := Rat.mul_pos (by grind) hy
      grind
    · grind

/-- Sterbenz's window, part 1: for 0 ≤ a ≤ 1 (lift in [1, 2]) the inverse b lies in [−2a, −a/2],
    so a + b is computed exactly in floating point. -/
theorem sterbenz_window_pos {a b : Rat} (h0 : 0 ≤ a) (h1 : a ≤ 1) (hb : b * (1 + a) = -a) :
    -2 * a ≤ b ∧ 2 * b ≤ -a := by
  have hp : (0:Rat) < 1 + a := by grind
  constructor
  · have : 0 ≤ (b + 2 * a) * (1 + a) := by
      have e : (b + 2 * a) * (1 + a) = a + 2 * (a * a) := by grind
      rw [e]; have := Rat.mul_nonneg h0 h0; grind
    have := nonneg_of_mul_pos hp this; grind
  · have : 0 ≤ (-a - 2 * b) * (1 + a) := by
      have e : (-a - 2 * b) * (1 + a) = a * (1 - a) := by grind
      rw [e]; exact Rat.mul_nonneg h0 (by grind)
    have := nonneg_of_mul_pos hp this; grind

/-- Sterbenz's window, part 2: for −1/2 ≤ a ≤ 0 (lift in [1/2, 1]) the inverse b lies in
    [−a/2, −2a]. -/
theorem sterbenz_window_neg {a b : Rat} (h0 : a ≤ 0) (h1 : -1 ≤ 2 * a) (hb : b * (1 + a) = -a) :
    -a ≤ 2 * b ∧ b ≤ -2 * a := by
  have hp : (0:Rat) < 1 + a := by grind
  constructor
  · have : 0 ≤ (2 * b + a) * (1 + a) := by
      have e : (2 * b + a) * (1 + a) = (-a) * (1 - a) := by grind
      rw [e]; exact Rat.mul_nonneg (by grind) (by grind)
    have := nonneg_of_mul_pos hp this; grind
  · have : 0 ≤ (-2 * a - b) * (1 + a) := by
      have e : (-2 * a - b) * (1 + a) = (-a) * (1 + 2 * a) := by grind
      rw [e]; exact Rat.mul_nonneg (by grind) (by grind)
    have := nonneg_of_mul_pos hp this; grind

/-! ## 5. Error bounds: θₙ stays inside ±((1+u)ⁿ − 1), and (1+u)ⁿ − 1 ≤ nu/(1 − nu) -/

/-- |a| ≤ u and |t| ≤ M give |a ◯ t| ≤ u ◯ M (the corner of the bounds bounds the corner). -/
theorem circ_bound {a t u M : Rat} (ha1 : -u ≤ a) (ha2 : a ≤ u) (ht1 : -M ≤ t) (ht2 : t ≤ M) :
    -(u ◯ M) ≤ a ◯ t ∧ a ◯ t ≤ u ◯ M := by
  have p1 := Rat.mul_nonneg (by grind : (0:Rat) ≤ u - a) (by grind : (0:Rat) ≤ M + t)
  have p2 := Rat.mul_nonneg (by grind : (0:Rat) ≤ u + a) (by grind : (0:Rat) ≤ M - t)
  have p3 := Rat.mul_nonneg (by grind : (0:Rat) ≤ u + a) (by grind : (0:Rat) ≤ M + t)
  have p4 := Rat.mul_nonneg (by grind : (0:Rat) ≤ u - a) (by grind : (0:Rat) ≤ M - t)
  unfold circ; constructor <;> grind

theorem mpow_nonneg {u : Rat} (hu : 0 ≤ u) : ∀ n, 0 ≤ mpow u n
  | 0 => by simp only [mpow]; grind
  | n + 1 => by
    simp only [mpow]; unfold circ
    have h := mpow_nonneg hu n
    have := Rat.mul_nonneg h hu; grind

/-- Higham's Lemma 3.1 in mult-add form: n step errors of size at most u fold to at most
    (1 + u)ⁿ − 1 = u^(◯n). -/
theorem theta_bound {u : Rat} (hu : 0 ≤ u) : ∀ (ds : List Rat), (∀ d ∈ ds, -u ≤ d ∧ d ≤ u) →
    -(mpow u ds.length) ≤ foldCirc ds ∧ foldCirc ds ≤ mpow u ds.length
  | [], _ => by simp only [foldCirc, List.length, mpow]; grind
  | d :: ds, h => by
    have hd := h d (by simp)
    have ih := theta_bound hu ds (fun x hx => h x (by simp [hx]))
    simp only [foldCirc, List.length]
    have hm : mpow u (ds.length + 1) = u ◯ mpow u ds.length := by
      rw [mpow, circ_comm]
    rw [hm]
    exact circ_bound hd.1 hd.2 ih.1 ih.2

theorem pow_nonneg_rat {x : Rat} (hx : 0 ≤ x) : ∀ n : Nat, 0 ≤ x ^ n
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by rw [Rat.pow_succ]; exact Rat.mul_nonneg (pow_nonneg_rat hx n) hx

/-- (1 + u)ⁿ(1 − nu) ≤ 1 for u ≥ 0: the γₙ bound without division. -/
theorem gamma_lift {u : Rat} (hu : 0 ≤ u) : ∀ n : Nat, (1 + u) ^ n * (1 - n * u) ≤ 1
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by
    have ih := gamma_lift hu n
    have hp := pow_nonneg_rat (by grind : (0:Rat) ≤ 1 + u) n
    rw [Rat.pow_succ]
    have hsq := Rat.mul_nonneg hu hu
    have step : (1 + u) * (1 - ((n + 1 : Nat) : Rat) * u) ≤ 1 - n * u := by
      have hn : (0:Rat) ≤ (n : Rat) := by exact_mod_cast Nat.zero_le n
      have := Rat.mul_nonneg (Rat.add_nonneg hn (by grind : (0:Rat) ≤ 1)) hsq
      push_cast; grind
    have := Rat.mul_le_mul_of_nonneg_left step hp
    grind

/-- γₙ: u^(◯n)·(1 − nu) ≤ nu, i.e. (1 + u)ⁿ − 1 ≤ nu/(1 − nu) whenever nu < 1. -/
theorem gamma_bound {u : Rat} (hu : 0 ≤ u) (n : Nat) : mpow u n * (1 - n * u) ≤ n * u := by
  have h := gamma_lift hu n
  rw [mpow_eq]; grind

/-- Bernoulli: compounding never loses to simple interest, a^(◯n) ≥ n·a for a ≥ −1. -/
theorem bernoulli {a : Rat} (ha : -1 ≤ a) : ∀ n : Nat, (n : Rat) * a ≤ mpow a n
  | 0 => by simp only [mpow]; grind
  | n + 1 => by
    have ih := bernoulli ha n
    have hstep : mpow a (n + 1) = mpow a n + a + mpow a n * a := rfl
    rw [hstep]
    -- mpow a n * a ≥ n a · a: (mpow a n − n a)·a ≥ … split on the sign of a
    have hn : (0:Rat) ≤ (n : Rat) := by exact_mod_cast Nat.zero_le n
    have hm : 0 ≤ (mpow a n - n * a) * (1 + a) :=
      Rat.mul_nonneg (by grind) (by grind)
    have hsq : 0 ≤ (n : Rat) * (a * a) := Rat.mul_nonneg hn (sq_nonneg_rat a)
    push_cast; grind

end Numerical
