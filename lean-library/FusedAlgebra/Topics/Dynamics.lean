/-
  Fused Hyperoperation Algebra: dynamical systems (the lift is a change of coordinates)
  ===================================================================================
  Core Lean 4 only (no mathlib), no `sorry`.

  Classical facts, stated through the circle operation (called mult-add before 1.0):
    * x ◯ x = (1+x)² − 1 is z ↦ z² through the lift z = 1 + x (Fatou, Julia 1918-19);
      the wall −1 is its superattracting fixed point, the identity 0 repels with multiplier 2.
    * Every dial square x ◯ₖ x is x ◯ x in other units; k = −2 is the logistic map at r = 2.
    * The logistic map at r = 4 is x ◯ x folded by X(w) = w·ι(w)/4, ι the circle inverse
      (Ulam and von Neumann 1947, through Joukowski's map).
    * Chebyshev polynomials are zⁿ seen through Joukowski.
    * Near 0 the rung flow is the parabolic Möbius flow x/(1 − t·x/2), whose Abel ruler −2/x adds t.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/Dynamics.lean; one exact twin of KnotTheory removed (see the LIBRARY comment).
-/
import FusedAlgebra.Base

namespace Dynamics
open CircleOp
open Lean.Grind

open scoped CircleOp

set_option linter.unusedSectionVars false

/-- n-fold iteration. -/
def iter {α : Type} (f : α → α) : Nat → α → α
  | 0, x => x
  | n + 1, x => f (iter f n x)

section Ring
variable {R : Type} [CommRing R]

/-- Circle-operation squaring. -/
def msq (x : R) : R := x ◯ x

/-- The lift conjugates circle-operation squaring to z ↦ z². -/
theorem msq_lift (x : R) : 1 + msq x = (1 + x) * (1 + x) := circ_lift x x

/-- (a^m)^n = a^(mn) for any commutative ring. -/
theorem pow_mul_R (a : R) (m n : Nat) : a ^ (m * n) = (a ^ m) ^ n := by
  rw [← pw_eq_pow, ← pw_eq_pow, ← pw_eq_pow, pw_mul]

/-- n steps: 1 + fⁿ(x) = (1 + x)^(2ⁿ). -/
theorem iter_msq_lift (x : R) (n : Nat) : 1 + iter msq n x = (1 + x) ^ (2 ^ n) := by
  induction n with
  | zero => show 1 + x = (1 + x) ^ (2 ^ 0); rw [Nat.pow_zero, ← pw_eq_pow]; show 1 + x = 1 * (1 + x); grind
  | succ n ih =>
    show 1 + msq (iter msq n x) = _
    rw [msq_lift, ih, Nat.pow_succ, pow_mul_R, Semiring.pow_succ, Semiring.pow_succ, Semiring.pow_zero]
    grind

/-- Circle powers compose like ordinary powers: (x^(◯m))^(◯n) = x^(◯mn). -/
theorem mpow_mpow (x : R) (m n : Nat) : mpow (mpow x m) n = mpow x (m * n) := by
  have h1 := mpow_lift (mpow x m) n
  have h2 := mpow_lift x m
  have h3 := mpow_lift x (m * n)
  rw [h2, ← pw_eq_pow, ← pw_eq_pow, ← pw_mul, pw_eq_pow, ← h3] at h1
  grind

/-- Every dial square is circle-operation squaring in other units: k·(x ◯ₖ x) = (kx) ◯ (kx). -/
theorem dial_square_units (k x : R) : k * circK k x x = msq (k * x) := circK_scale k x x

/-- Dial powers too: k·x^(◯ₖn) = (kx)^(◯n). -/
theorem dial_power_units (k x : R) (n : Nat) : k * wpow k x n = mpow (k * x) n := by
  have h1 := wpow_lift k x n
  have h2 := mpow_lift (k * x) n
  grind

/-- k = −2 (parity): the dial square is the logistic map at r = 2. -/
theorem parity_logistic (x : R) : circK (-2) x x = 2 * x * (1 - x) := by unfold circK; grind

/-- k = −1 (failure union): the dial square is 1 − (1 − x)². -/
theorem failure_square (x : R) : circK (-1) x x = 1 - (1 - x) * (1 - x) := by unfold circK; grind

/-- The wall is superattracting: (−1 + h) ◯ (−1 + h) = −1 + h² (multiplier 0). -/
theorem wall_superattracting (h : R) : msq (-1 + h) = -1 + h * h := by unfold msq circ; grind

/-- The identity repels: h ◯ h = 2h + h² (multiplier 2). -/
theorem identity_repelling (h : R) : msq h = 2 * h + h * h := by unfold msq circ; grind

/-- The wall stays put under every iterate. -/
theorem iter_wall (n : Nat) : iter msq n (-1 : R) = -1 := by
  induction n with
  | zero => rfl
  | succ n ih => show msq (iter msq n (-1 : R)) = -1; rw [ih]; exact neg_one_absorbs (-1)

/-- So does the identity. -/
theorem iter_zero (n : Nat) : iter msq n (0 : R) = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => show msq (iter msq n (0 : R)) = 0; rw [ih]; exact circ_zero 0

/-- Periodic points: fⁿ(x) = x exactly when (1 + x)^(2ⁿ) = 1 + x. -/
theorem periodic_iff (x : R) (n : Nat) : iter msq n x = x ↔ (1 + x) ^ (2 ^ n) = 1 + x := by
  rw [← iter_msq_lift]
  constructor
  · intro h; rw [h]
  · intro h; grind

/-- Higher circle powers: the wall is fixed by x ↦ x^(◯(n+1)). -/
theorem wall_fixed_mpow (n : Nat) : mpow (-1 : R) (n + 1) = -1 := by
  rw [mpow, circ_comm]; exact neg_one_absorbs _

-- LIBRARY: chain_circ (multipliers multiply along a composed orbit, so multiplier − 1 composes by ◯) was an exact twin of KnotTheory.poly_circ and is removed.

/-- Derivatives ride along as jets (value, slope): one step of x ◯ x. -/
def msqJet (p : R × R) : R × R := (msq p.1, 2 * (1 + p.1) * p.2)

theorem iter_jet_fst (p : R × R) (n : Nat) : (iter msqJet n p).1 = iter msq n p.1 := by
  induction n with
  | zero => rfl
  | succ n ih => show msq (iter msqJet n p).1 = msq (iter msq n p.1); rw [ih]

/-- In the log ruler L = log(1+x) every step doubles slopes: the slope of fⁿ divided by
    1 + fⁿ(x) is 2ⁿ times the starting slope divided by 1 + x (Lyapunov exponent log 2). -/
theorem iter_jet_ruler (p : R × R) (n : Nat) :
    (iter msqJet n p).2 * (1 + p.1) = 2 ^ n * p.2 * (1 + (iter msqJet n p).1) := by
  induction n with
  | zero => show p.2 * (1 + p.1) = 2 ^ 0 * p.2 * (1 + p.1); rw [Semiring.pow_zero]; grind
  | succ n ih =>
    show (2 * (1 + (iter msqJet n p).1) * (iter msqJet n p).2) * (1 + p.1) =
      2 ^ (n + 1) * p.2 * (1 + msq (iter msqJet n p).1)
    rw [msq_lift, Semiring.pow_succ]
    have : (iter msqJet n p).2 * (1 + p.1) = 2 ^ n * p.2 * (1 + (iter msqJet n p).1) := ih
    grind

/-- Flows: running y′ = 1 + y for time s then t is combining the two lifts' parts by the circle operation:
    (a ◯ (e₁ − 1)) ◯ (e₂ − 1) = a ◯ (e₁e₂ − 1), with eᵢ = e^(tᵢ). -/
theorem flow_time_law (a e1 e2 : R) : (a ◯ (e1 - 1)) ◯ (e2 - 1) = a ◯ (e1 * e2 - 1) := by
  unfold circ; grind

/-! ### Chebyshev polynomials: zⁿ through Joukowski -/

def cheb (u : R) : Nat → R
  | 0 => 1
  | 1 => u
  | n + 2 => 2 * u * cheb u (n + 1) - cheb u n

/-- If 2u = z + z⁻¹ then 2·Tₙ(u) = zⁿ + z⁻ⁿ. -/
theorem cheb_joukowski (u z zi : R) (hz : z * zi = 1) (hu : 2 * u = z + zi) (n : Nat) :
    2 * cheb u n = z ^ n + zi ^ n ∧ 2 * cheb u (n + 1) = z ^ (n + 1) + zi ^ (n + 1) := by
  induction n with
  | zero => rw [Semiring.pow_zero, Semiring.pow_zero, Nat.zero_add, ← pw_eq_pow, ← pw_eq_pow]
            show 2 * 1 = 1 + 1 ∧ 2 * u = 1 * z + 1 * zi; constructor <;> grind
  | succ n ih =>
    refine ⟨ih.2, ?_⟩
    show 2 * (2 * u * cheb u (n + 1) - cheb u n) = _
    have e1 : z ^ (n + 1 + 1) = z ^ n * z * z := by rw [Semiring.pow_succ, Semiring.pow_succ]
    have e2 : zi ^ (n + 1 + 1) = zi ^ n * zi * zi := by rw [Semiring.pow_succ, Semiring.pow_succ]
    have e3 : z ^ (n + 1) = z ^ n * z := Semiring.pow_succ z n
    have e4 : zi ^ (n + 1) = zi ^ n * zi := Semiring.pow_succ zi n
    rw [e1, e2]; rw [e3, e4] at ih
    have h1 := ih.1; have h2 := ih.2
    grind

/-- The logistic map at r = 4, read in u = 1 − 2x, is T₂. -/
theorem logistic_T2 (x : R) : 1 - 2 * (4 * x * (1 - x)) = cheb (1 - 2 * x) 2 := by
  show _ = 2 * (1 - 2 * x) * (1 - 2 * x) - 1; grind

/-- Chebyshev maps commute like powers: T₂ ∘ T₃ = T₃ ∘ T₂ = T₆. -/
theorem cheb_comp (u : R) :
    cheb (cheb u 3) 2 = cheb u 6 ∧ cheb (cheb u 2) 3 = cheb u 6 := by
  simp only [cheb]; constructor <;> grind

end Ring

/-! ## Fields: the fold, the dial wall, the parabolic model -/

section Field
variable {K : Type} [Field K]

/-- The fold X(w) = w·ι(w)/4, ι the circle inverse. -/
def fold (w : K) : K := w * circInv w / 4

theorem fold_eq (w : K) (h : 1 + w ≠ 0) : fold w = -(w * w) / (4 * (1 + w)) := by
  unfold fold circInv; grind

/-- Since w ◯ ι(w) = 0, the sum and product of w and its inverse are negatives. -/
theorem sum_is_neg_product (w : K) (h : 1 + w ≠ 0) : w + circInv w = -(w * circInv w) := by
  have := circ_circInv w h; unfold circ at this; grind

/-- The circle inverse undoes itself. -/
theorem inv_inv (w : K) (h : 1 + w ≠ 0) : circInv (circInv w) = w := by
  unfold circInv
  have : 1 + -w / (1 + w) = 1 / (1 + w) := by grind
  rw [this]; grind

/-- The fold is 2-to-1: w and ι(w) land on the same point. -/
theorem fold_inv (w : K) (h : 1 + w ≠ 0) : fold (circInv w) = fold w := by
  unfold fold; rw [inv_inv w h]; grind

/-- Squaring commutes with the inverse: ι(w) ◯ ι(w) = ι(w ◯ w). -/
theorem msq_inv (w : K) (h : 1 + w ≠ 0) : msq (circInv w) = circInv (msq w) := by
  unfold msq circ circInv; grind

/-- The logistic map at r = 4 is circle-operation squaring folded: 4X(1 − X) = X(w ◯ w). -/
theorem logistic_fold (w : K) (h : 1 + w ≠ 0) : 4 * fold w * (1 - fold w) = fold (msq w) := by
  have h2 : 1 + msq w ≠ 0 := by
    rw [msq_lift]; intro e
    have : (1 + w) * (1 + w) * ((1 + w)⁻¹ * (1 + w)⁻¹) = 0 := by rw [e]; grind
    grind
  rw [fold_eq w h, fold_eq (msq w) h2]
  unfold msq circ; grind

/-- Joukowski in circle-operation terms: z + 1/z = 2 + (w + ι(w)) for z = 1 + w. -/
theorem joukowski_circ (w : K) (h : 1 + w ≠ 0) :
    (1 + w) + 1 / (1 + w) = 2 + (w + circInv w) := by
  unfold circInv; grind

/-- The dial's wall −1/k is the superattracting fixed point of x ◯ₖ x. -/
theorem dial_wall_superattracting (k h : K) (hk : k ≠ 0) :
    circK k (-1 / k + h) (-1 / k + h) = -1 / k + k * h * h := by
  have e : k * (1 / k) = 1 := by grind
  unfold circK
  have : k * (-1 / k + h) * (-1 / k + h) = (1 / k) - 2 * h + k * h * h := by
    have e2 : -1 / k = -(1 / k) := by grind
    rw [e2]; grind
  grind

/-- The parabolic model flow x ↦ x/(1 − t·x/2) (the flow of x′ = x²/2, the first term of the
    rung flow's vector field) is a one-parameter group. -/
def pflow (t x : K) : K := x / (1 - t * x / 2)

theorem pflow_add (s t x : K) (h1 : 1 - t * x / 2 ≠ 0) (h2 : 1 - (s + t) * x / 2 ≠ 0) :
    pflow s (pflow t x) = pflow (s + t) x := by
  unfold pflow
  have h3 : 1 - s * (x / (1 - t * x / 2)) / 2 = (1 - (s + t) * x / 2) / (1 - t * x / 2) := by grind
  rw [h3]; grind

/-- Its Abel ruler −2/x adds the time: −2/φₜ(x) = −2/x + t (Fatou's coordinate, exactly). -/
theorem pflow_ruler (t x : K) (two : (2 : K) ≠ 0) (hx : x ≠ 0) (h1 : 1 - t * x / 2 ≠ 0) :
    -2 / pflow t x = -2 / x + t := by
  unfold pflow
  have e : -2 / (x / (1 - t * x / 2)) = -2 * (1 - t * x / 2) / x := by grind
  rw [e]
  have e2 : (2 : K) * (1 / 2) = 1 := by grind
  grind

end Field

end Dynamics
