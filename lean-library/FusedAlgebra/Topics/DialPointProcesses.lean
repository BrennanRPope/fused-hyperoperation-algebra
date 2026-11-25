/-
  Point processes on the dial
  ===========================

  Companion page: ../index.html.  Numbers: ../checks.py.
  Core Lean 4.37 only (no mathlib), no `sorry`.  Intended module path:
  FusedAlgebra/Topics/DialPointProcesses.lean (imports Base only).

  The α-determinantal processes of Shirai and Takahashi (2003) with α read as
  the dial weight k: k = −1 fermions (determinantal), k = 0 Poisson, k = +1
  bosons (permanental).  On a region the kernel's eigenvalues μᵢ act as
  independent modes.  What this file proves, by part:

  1. The hole lift.  1 + k·(μ₁ ◯ₖ μ₂ ◯ₖ …) = ∏(1 + k·μᵢ) = det(I + kK) on the
     region, so P(no point) = (1 + k·fold)^(−1/k).  An eigenvalue on the wall
     −1/k makes the product 0: that region is never empty.
  2. Number variance.  Each mode contributes μ(1 + kμ), its value times its
     lift, which vanishes at 0 and at the wall.  Summed: Var = Σμ + k·Σμ², a
     straight line in k (the Fano line).  Second factorial moment
     (Σμ)² + k·Σμ².
  3. Superposition.  m independent copies of (k, K) are (k/m, m·K): the dial
     identity c·(a ◯ₖ b) = (c·a) ◯_{k/c} (c·b), and the variance matches.
  4. Tilting.  Weighting a mode's count by tᴺ moves its mean; in the ruler
     x = (mean − μ)/Var the moves compose by the landscape lens law at
     p = 1 + 2s, q = s(1 + s), s = kμ, which always has p² − 4q = 1.
     s = 0 is the circle operation (called mult-add before 1.0), s = −1 the union of failures, s = −1/2 velocity
     addition with speed limit 2.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/DialPointProcesses.lean; merged unchanged; no exact twins. Its Rat helper sq_nonneg'
  stays local like the other topics' Rat helpers.
-/
import FusedAlgebra.Base

namespace DialPointProcesses

open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! # 1. The dial fold of the eigenvalues and the hole -/

/-- μ₁ ◯ₖ μ₂ ◯ₖ … ◯ₖ μₙ. -/
def fold (k : Rat) : List Rat → Rat
  | [] => 0
  | x :: l => circK k x (fold k l)

/-- Two independent groups of modes: the fold of the whole is the dial sum of the parts. -/
theorem fold_append (k : Rat) (A B : List Rat) : fold k (A ++ B) = circK k (fold k A) (fold k B) := by
  induction A with
  | nil => simp only [List.nil_append, fold]; unfold circK; grind
  | cons x A ih => simp only [List.cons_append, fold]; rw [ih, circK_assoc]

/-- The hole lift: 1 + k·fold = ∏(1 + k·μᵢ), the determinant det(I + kK). -/
theorem fold_lift (k : Rat) (l : List Rat) :
    1 + k * fold k l = prodL (l.map (fun x => 1 + k * x)) := by
  induction l with
  | nil => simp only [fold, List.map, prodL]; grind
  | cons x l ih => simp only [fold, List.map, prodL]; rw [circK_lift, ih]

/-- k·fold is the plain circle operation of the scaled eigenvalues k·μᵢ. -/
theorem fold_scale (k : Rat) (l : List Rat) : k * fold k l = foldCirc (l.map (fun x => k * x)) := by
  induction l with
  | nil => simp only [fold, List.map, foldCirc]; grind
  | cons x l ih => simp only [fold, List.map, foldCirc]; rw [circK_scale, ih]

/-- Poisson (k = 0): the fold is the expected count Σμ. -/
theorem fold_poisson (l : List Rat) : fold 0 l = sumL l := by
  induction l with
  | nil => rfl
  | cons x l ih => simp only [fold, sumL]; rw [circK_weight_zero, ih]

/-- Fermions (k = −1): the fold is P(at least one point) = 1 − ∏(1 − μᵢ). -/
theorem fold_fermion (l : List Rat) : fold (-1) l = 1 - prodL (l.map (fun x => 1 - x)) := by
  induction l with
  | nil => simp only [fold, List.map, prodL]; grind
  | cons x l ih => simp only [fold, List.map, prodL]; rw [ih]; unfold circK; grind

/-- Bosons (k = +1): the fold is the plain circle operation of the eigenvalues. -/
theorem fold_boson (l : List Rat) : fold 1 l = foldCirc l := by
  induction l with
  | nil => rfl
  | cons x l ih => simp only [fold, foldCirc]; rw [circK_one, ih]

/-- Order of the modes never matters. -/
theorem fold_swap (k x y : Rat) (l : List Rat) : fold k (x :: y :: l) = fold k (y :: x :: l) := by
  simp only [fold]; unfold circK; grind

/-- Two modes, written out: the corner carries the dial weight. -/
theorem fold_two (k a b : Rat) : fold k [a, b] = a + b + k * a * b := by
  simp only [fold]; unfold circK; grind

theorem prodL_zero_of_mem (l : List Rat) (y : Rat) (h : y ∈ l) (h0 : y = 0) : prodL l = 0 := by
  induction l with
  | nil => cases h
  | cons a l ih =>
    simp only [prodL]
    cases h with
    | head => rw [h0]; grind
    | tail _ h' => rw [ih h']; grind

/-- An eigenvalue on the wall −1/k (k·μ = −1) kills the lift: 1 + k·fold = 0,
    so the hole probability (1 + k·fold)^(−1/k) is 0 for k < 0. -/
theorem wall_mode_fills (k x : Rat) (l : List Rat) (hx : x ∈ l) (hw : k * x = -1) :
    1 + k * fold k l = 0 := by
  rw [fold_lift]
  apply prodL_zero_of_mem _ (1 + k * x)
  · exact List.mem_map_of_mem hx
  · grind

/-- The wall absorbs the fold: once one mode is at −1/k, so is the whole. -/
theorem wall_absorbs (k x : Rat) (hk : k ≠ 0) : circK k (-1 / k) x = -1 / k := by
  rw [circK_comm]; exact circK_wall k x hk

/-- k = −1/m with eigenvalue m (the wall): each of the m slots is filled for sure. -/
theorem slots_full (m : Rat) (hm : m ≠ 0) : 1 + (-1 / m) * m = 0 := by
  have : (-1 / m) * m = -1 := by
    rw [Rat.div_def, Rat.mul_assoc, Rat.inv_mul_cancel m hm]; grind
  grind

/-! # 2. Number variance: value times lift -/

/-- One mode's variance: μ(1 + kμ). Binomial(m, μ/m) at k = −1/m, Poisson at 0,
    negative binomial of shape 1/k at k > 0 (Morris 1982's quadratic variance). -/
def modeVar (k x : Rat) : Rat := x * (1 + k * x)

theorem modeVar_expand (k x : Rat) : modeVar k x = x + k * (x * x) := by unfold modeVar; grind

theorem modeVar_empty (k : Rat) : modeVar k 0 = 0 := by unfold modeVar; grind

/-- On the wall the count is fixed: no variance. -/
theorem modeVar_wall (k : Rat) (hk : k ≠ 0) : modeVar k (-1 / k) = 0 := by
  unfold modeVar
  have : k * (-1 / k) = -1 := by
    rw [Rat.div_def, ← Rat.mul_assoc, Rat.mul_comm k (-1), Rat.mul_assoc, Rat.mul_inv_cancel k hk]; grind
  rw [this]; grind

/-- Poisson: variance equals mean. -/
theorem modeVar_poisson (x : Rat) : modeVar 0 x = x := by unfold modeVar; grind

/-- Fano factor of a mode, Var/mean = 1 + kμ (the lift). -/
theorem mode_fano (k x : Rat) (hx : x ≠ 0) : modeVar k x / x = 1 + k * x := by
  unfold modeVar
  rw [Rat.div_def, Rat.mul_comm x, Rat.mul_assoc, Rat.mul_inv_cancel x hx]; grind

theorem sq_nonneg' (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

/-- Fermion side (k ≤ 0) is sub-Poisson, boson side (k ≥ 0) super-Poisson. -/
theorem fermion_subpoisson (k x : Rat) (hk : k ≤ 0) : modeVar k x ≤ x := by
  rw [modeVar_expand]
  have h := sq_nonneg' x
  have : k * (x * x) ≤ 0 := by
    have h2 : 0 ≤ (-k) * (x * x) := Rat.mul_nonneg (by grind) h
    grind
  grind

theorem boson_superpoisson (k x : Rat) (hk : 0 ≤ k) : x ≤ modeVar k x := by
  rw [modeVar_expand]
  have : 0 ≤ k * (x * x) := Rat.mul_nonneg hk (sq_nonneg' x)
  grind

/-- Total variance over the modes of a region. -/
def varL (k : Rat) (l : List Rat) : Rat := sumL (l.map (modeVar k))

/-- Sum of squared eigenvalues, tr K². -/
def sqL (l : List Rat) : Rat := sumL (l.map (fun x => x * x))

/-- Var N = Σμ + k·Σμ² = tr K + k·tr K². -/
theorem var_formula (k : Rat) (l : List Rat) : varL k l = sumL l + k * sqL l := by
  induction l with
  | nil => simp only [varL, sqL, List.map, sumL]; grind
  | cons x l ih =>
    simp only [varL, sqL, List.map, sumL] at ih ⊢
    rw [ih, modeVar_expand]; grind

/-- The Fano line: variance is affine in k, through the Poisson value at k = 0. -/
theorem fano_line (k : Rat) (l : List Rat) : varL k l = varL 0 l + k * sqL l := by
  rw [var_formula, var_formula]; grind

/-- Disjoint groups of modes: variances add. -/
theorem var_append (k : Rat) (A B : List Rat) : varL k (A ++ B) = varL k A + varL k B := by
  simp only [varL, List.map_append]; rw [sumL_append]

/-- Second factorial moment of independent modes, built one mode at a time:
    E[N(N−1)] for N = M + N', M one mode with E M(M−1) = (1 + k)μ². -/
def fact2 (k : Rat) : List Rat → Rat
  | [] => 0
  | x :: l => (1 + k) * (x * x) + fact2 k l + 2 * x * sumL l

theorem fact2_formula (k : Rat) (l : List Rat) : fact2 k l = sumL l * sumL l + k * sqL l := by
  induction l with
  | nil => simp only [fact2, sumL, sqL, List.map]; grind
  | cons x l ih =>
    simp only [fact2, sqL, List.map, sumL] at ih ⊢
    rw [ih]; grind

/-- Var = E N(N−1) + E N − (E N)², checked against the two formulas. -/
theorem var_from_fact2 (k : Rat) (l : List Rat) :
    varL k l = fact2 k l + sumL l - sumL l * sumL l := by
  rw [var_formula, fact2_formula]; grind

/-! # 3. Superposition: m copies of (k, K) are (k/m, m·K) -/

/-- The dial identity behind it: scaling values by c and the weight by 1/c. -/
theorem dial_rescale (k a b c : Rat) (hc : c ≠ 0) : c * circK k a b = circK (k / c) (c * a) (c * b) := by
  unfold circK
  have : k / c * (c * a) * (c * b) = c * (k * a * b) := by
    rw [Rat.div_def]
    have e : k * c⁻¹ * (c * a) * (c * b) = (c * c⁻¹) * (c * (k * a * b)) := by grind
    rw [e, Rat.mul_inv_cancel c hc]; grind
  rw [this]; grind

/-- The lifted kernel k·K is what the count sees: unchanged by (k, K) ↦ (k/c, cK). -/
theorem lifted_kernel_same (k x c : Rat) (hc : c ≠ 0) : (k / c) * (c * x) = k * x := by
  rw [Rat.div_def]
  have e : k * c⁻¹ * (c * x) = (c * c⁻¹) * (k * x) := by grind
  rw [e, Rat.mul_inv_cancel c hc]; grind

theorem modeVar_rescale (k x c : Rat) (hc : c ≠ 0) : modeVar (k / c) (c * x) = c * modeVar k x := by
  unfold modeVar; rw [lifted_kernel_same k x c hc]; grind

/-- Variance of m superposed copies: Var(k/m, mK) = m·Var(k, K). -/
theorem var_superpose (k c : Rat) (hc : c ≠ 0) (l : List Rat) :
    varL (k / c) (l.map (fun x => c * x)) = c * varL k l := by
  induction l with
  | nil => simp only [varL, List.map, sumL]; grind
  | cons x l ih =>
    simp only [varL, List.map, sumL] at ih ⊢
    rw [ih, modeVar_rescale k x c hc]; grind

/-- Two fermion copies, each with eigenvalue μ/2, leave a mode empty with chance (1 − μ/2)²:
    the k = −1/2 hole with eigenvalue μ. -/
theorem two_fermion_copies (x : Rat) : (1 - x / 2) * (1 - x / 2) = (1 + (-1 / 2) * x) * (1 + (-1 / 2) * x) := by
  grind

/-! # 4. Tilting a mode: the landscape lens law -/

/-- The landscape law (a + b + p·ab)/(1 − q·ab). -/
def lens (p q a b : Rat) : Rat := (a + b + p * a * b) / (1 - q * a * b)

/-- The tilt ruler: tilt factor t = 1 + u (u composes by the circle operation), x = u/(1 − s·u). -/
def ruler (s u : Rat) : Rat := u / (1 - s * u)

/-- Tilted mean: weighting a mode's count by tᴺ gives mean μt/(1 + kμ(1 − t)).
    Read in the ruler (mean − μ)/Var it is u/(1 − kμ·u), u = t − 1. -/
theorem tilt_in_ruler (k x t : Rat) (hx : x ≠ 0) (hv : 1 + k * x ≠ 0)
    (hd : 1 + k * x * (1 - t) ≠ 0) :
    (x * t / (1 + k * x * (1 - t)) - x) / (x * (1 + k * x)) = ruler (k * x) (t - 1) := by
  unfold ruler
  have e1 : 1 - k * x * (t - 1) = 1 + k * x * (1 - t) := by grind
  rw [e1]
  grind

/-- Tilts compose by the circle operation in u, and by the lens law in the ruler:
    p = 1 + 2s, q = s(1 + s). -/
theorem tilt_lens (s u v : Rat) (hu : 1 - s * u ≠ 0) (hv : 1 - s * v ≠ 0)
    (huv : 1 - s * (u ◯ v) ≠ 0) :
    ruler s (u ◯ v) = lens (1 + 2 * s) (s * (1 + s)) (ruler s u) (ruler s v) := by
  unfold ruler lens circ at *
  grind

/-- Every dial count lands on one parabola: p² − 4q = 1, two real walls. -/
theorem count_parabola (s : Rat) : (1 + 2 * s) * (1 + 2 * s) - 4 * (s * (1 + s)) = 1 := by grind

/-- Its two walls are the empty mode and the full mode: 1 + p·x + q·x² = (1 + s·x)(1 + (1+s)·x). -/
theorem count_walls (s x : Rat) :
    1 + (1 + 2 * s) * x + s * (1 + s) * (x * x) = (1 + s * x) * (1 + (1 + s) * x) := by grind

/-- The mean's rate, re-centred: V(μ + h) = V(μ) + (1 + 2kμ)h + k·h². -/
theorem rate_shift (k x h : Rat) : modeVar k (x + h) = modeVar k x + (1 + 2 * k * x) * h + k * (h * h) := by
  unfold modeVar; grind

/-- Shape p = 2·Fano − 1 (Fano = 1 + s), so the units-free shape p²/(p² − 4q) = (2F − 1)². -/
theorem shape_is_fano (s : Rat) : 1 + 2 * s = 2 * (1 + s) - 1 := by grind

/-- Particle-hole: s and −1 − s give the same shape. -/
theorem particle_hole (s : Rat) : (1 + 2 * (-1 - s)) * (1 + 2 * (-1 - s)) = (1 + 2 * s) * (1 + 2 * s) := by grind

/-- s = 0 (rare points, Poisson): the circle operation. -/
theorem empty_mode_multadd (a b : Rat) : lens (1 + 2 * 0) (0 * (1 + 0)) a b = a ◯ b := by
  unfold lens circ; grind

/-- s = −1 (a full fermion mode): the union of failures, k = −1 on the dial line. -/
theorem full_mode_union (a b : Rat) : lens (1 + 2 * (-1)) ((-1) * (1 + (-1))) a b = circK (-1) a b := by
  unfold lens circK; grind

/-- s = −1/2 (a half-filled fermion mode): velocity addition with speed limit 2. -/
theorem half_mode_velocity (a b : Rat) :
    lens (1 + 2 * (-1 / 2)) ((-1 / 2) * (1 + (-1 / 2))) a b = (a + b) / (1 + a * b / 4) := by
  unfold lens
  have e1 : a + b + (1 + 2 * (-1 / 2)) * a * b = a + b := by grind
  have e2 : 1 - (-1 / 2) * (1 + (-1 / 2)) * a * b = 1 + a * b / 4 := by grind
  rw [e1, e2]

end DialPointProcesses
