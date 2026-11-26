/-
  Decrement paths (DRAFT, not yet built)
  ======================================

  Part of Fused Hyperoperation Algebra. Page and checks: ../ (decrement-paths).
  Written for the Mathlib library (lean-library-mathlib, Lean 4.37.0-rc1). NOT BUILT: the machine
  that wrote it has no Lean. The statements are the point; a proof may need repair on first build.

  A policy can end by several causes in a year, with absolute rates qᵢ. The point
  p(t) = (chance each cause has struck by time t) runs through the cube from 0 to q.
  The chance of ending at all is the k = −1 fold (the union rule), 1 − ∏ (1 − pᵢ).
  The dependent rate of cause j along a path is ∫₀¹ (∏_{i ≠ j} (1 − pᵢ t)) · pⱼ' t dt.

  Proved here:
  * `path_total`: on every path the dependent rates add up to the same total (path independence).
  * `udd_corner`: uniform timing (UDD, the straight path pᵢ t = t qᵢ) gives cause j the corners of
    the cube that contain j, each shared equally among its members. By Harsanyi (1959) that equal
    split is the Shapley value of the game v(S) = 1 − ∏_S (1 − qᵢ); for 2 and 3 causes it is already
    proved in core Lean as `GameTheory.shapley2_union` and `GameTheory.shapley3_corners`.
  * `const_force_split`: constant force (pᵢ t = 1 − e^{−μᵢ t}, the straight path in log coordinates)
    gives cause j the share μⱼ / Σ μ of the total.

  Known: the UDD and constant-force formulas (Bowers et al. 1986), Harsanyi 1959, Owen 1972.
  Ours: reading a timing assumption as a path through the cube.
-/
import Mathlib

open Finset intervalIntegral

namespace DecrementPaths

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The union rule (the k = −1 fold): the chance that some cause strikes. -/
noncomputable def total (p : ι → ℝ) : ℝ := 1 - ∏ i, (1 - p i)

/-- The dependent rate of cause `j` along the path `p` with velocity `p'`. -/
noncomputable def rate (p p' : ι → ℝ → ℝ) (j : ι) : ℝ :=
  ∫ t in (0 : ℝ)..1, (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t

/-- The derivative of the total along a path is the sum of the causes' instantaneous rates. -/
theorem hasDerivAt_total (p p' : ι → ℝ → ℝ) (t : ℝ)
    (hp : ∀ i, HasDerivAt (p i) (p' i t) t) :
    HasDerivAt (fun s => total (fun i => p i s))
      (∑ j, (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) t := by
  have h1 : ∀ i ∈ (univ : Finset ι), HasDerivAt (fun s => 1 - p i s) (-p' i t) t :=
    fun i _ => (hp i).const_sub 1
  have h2 := (HasDerivAt.finset_prod h1).const_sub 1
  unfold total
  convert h2 using 1
  simp only [smul_eq_mul, mul_neg, sum_neg_distrib, neg_neg]

/-- Path independence: whatever the timing, the dependent rates add up to the total. -/
theorem path_total (p p' : ι → ℝ → ℝ)
    (hp : ∀ i t, HasDerivAt (p i) (p' i t) t)
    (hc : ∀ i, Continuous (p' i)) :
    ∑ j, rate p p' j = total (fun i => p i 1) - total (fun i => p i 0) := by
  have hcp : ∀ i, Continuous (p i) := fun i =>
    continuous_iff_continuousAt.2 fun t => (hp i t).continuousAt
  have hint : ∀ j ∈ (univ : Finset ι), IntervalIntegrable
      (fun t => (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) MeasureTheory.volume 0 1 := by
    intro j _
    apply Continuous.intervalIntegrable
    exact (continuous_finset_prod _ fun i _ => continuous_const.sub (hcp i)).mul (hc j)
  unfold rate
  rw [← intervalIntegral.integral_finset_sum hint]
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt
  · intro t _; exact hasDerivAt_total p p' t (fun i => hp i t)
  · apply Continuous.intervalIntegrable
    exact continuous_finset_sum _ fun j _ =>
      (continuous_finset_prod _ fun i _ => continuous_const.sub (hcp i)).mul (hc j)

/-- UDD: on the straight path pᵢ t = t qᵢ, cause j gets every corner S ∪ {j} of the cube,
    worth (−1)^|S| qⱼ ∏_S qᵢ, divided equally among its |S| + 1 members. -/
theorem udd_corner (q : ι → ℝ) (j : ι) :
    rate (fun i t => t * q i) (fun i _ => q i) j
      = ∑ S ∈ (univ.erase j).powerset,
          (-1) ^ S.card * q j * (∏ i ∈ S, q i) / (S.card + 1) := by
  unfold rate
  -- expand the product over the other causes into corners
  have hexp : ∀ t : ℝ, (∏ i ∈ univ.erase j, (1 - t * q i)) * q j
      = ∑ S ∈ (univ.erase j).powerset, ((-1) ^ S.card * q j * ∏ i ∈ S, q i) * t ^ S.card := by
    intro t
    have : ∀ i, (1 : ℝ) - t * q i = (-(t * q i)) + 1 := fun i => by ring
    simp_rw [this, prod_add, prod_const_one, mul_one, sum_mul]
    refine sum_congr rfl fun S _ => ?_
    rw [prod_neg, prod_mul_distrib, prod_const]
    ring
  simp_rw [hexp]
  rw [intervalIntegral.integral_finset_sum (fun S _ => by
    apply Continuous.intervalIntegrable; fun_prop)]
  refine sum_congr rfl fun S _ => ?_
  rw [intervalIntegral.integral_const_mul, integral_pow]
  push_cast
  field_simp
  ring

/-- Constant force: with pᵢ t = 1 − e^{−μᵢ t}, cause j gets μⱼ (1 − e^{−M}) / M, M = Σ μ,
    which is the share μⱼ / M of the total 1 − e^{−M}. -/
theorem const_force_split (μ : ι → ℝ) (j : ι) (hM : ∑ i, μ i ≠ 0) :
    rate (fun i t => 1 - Real.exp (-(μ i * t))) (fun i t => μ i * Real.exp (-(μ i * t))) j
      = μ j * (1 - Real.exp (-∑ i, μ i)) / ∑ i, μ i := by
  unfold rate
  have hpt : ∀ t : ℝ, (∏ i ∈ univ.erase j, (1 - (1 - Real.exp (-(μ i * t)))))
      * (μ j * Real.exp (-(μ j * t))) = μ j * Real.exp (-((∑ i, μ i) * t)) := by
    intro t
    simp only [sub_sub_cancel]
    rw [← Real.exp_sum, mul_comm (μ j), mul_assoc, ← Real.exp_add]
    congr 2
    rw [← add_sum_erase _ _ (mem_univ j), sum_mul, ← add_sum_erase _ _ (mem_univ j)]
    simp only [sum_neg_distrib]
    ring
  simp_rw [hpt]
  rw [intervalIntegral.integral_const_mul]
  have hd : ∀ t ∈ Set.uIcc (0 : ℝ) 1, HasDerivAt
      (fun s => -Real.exp (-((∑ i, μ i) * s)) / ∑ i, μ i) (Real.exp (-((∑ i, μ i) * t))) t := by
    intro t _
    have := ((hasDerivAt_id t).const_mul (∑ i, μ i)).neg.exp.neg.div_const (∑ i, μ i)
    convert this using 1
    field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp only [mul_zero, neg_zero, Real.exp_zero, mul_one]
  field_simp
  ring

/-- The total on the constant-force path is 1 − e^{−M}: so the split is in proportion to the μ. -/
theorem const_force_total (μ : ι → ℝ) :
    total (fun i => 1 - Real.exp (-(μ i * 1))) = 1 - Real.exp (-∑ i, μ i) := by
  unfold total
  simp only [mul_one, sub_sub_cancel, ← Real.exp_sum, sum_neg_distrib]

end DecrementPaths
