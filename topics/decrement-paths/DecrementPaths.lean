/-
  Decrement paths (DRAFT, version 3)
  ==================================

  Part of Fused Hyperoperation Algebra. Page and checks: ../ (decrement-paths).
  Written for the Mathlib library (Lean 4.37.0-rc1). Built only by CI's report-only Mathlib step; the
  machine that wrote it has no Lean. Version 2: the lemma names that failed in the first CI run
  (finset product rule, finset continuity, integral of a finset sum) are replaced by small induction
  proofs here (hasDerivAt_prod', continuous_prod', continuous_sum', integral_sum').
  Version 3: two proof steps fixed after the second CI run (the fundamental theorem step in
  `path_total`, and a leftover negated-function goal in `const_force_split`). Statements unchanged.

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

/-! ### Helpers, proved by induction so they do not depend on Mathlib's big-operator calculus names -/

omit [Fintype ι] in
/-- Product rule for a finite product. -/
theorem hasDerivAt_prod' (s : Finset ι) (g g' : ι → ℝ → ℝ) (t : ℝ)
    (hg : ∀ i ∈ s, HasDerivAt (g i) (g' i t) t) :
    HasDerivAt (fun x => ∏ i ∈ s, g i x) (∑ j ∈ s, (∏ i ∈ s.erase j, g i t) * g' j t) t := by
  revert hg
  refine Finset.induction_on s ?_ ?_
  · intro _
    simp only [Finset.prod_empty, Finset.sum_empty]
    exact hasDerivAt_const t (1 : ℝ)
  · intro a s ha ih hg
    have hs : ∀ i ∈ s, HasDerivAt (g i) (g' i t) t := fun i hi => hg i (Finset.mem_insert_of_mem hi)
    have h := (hg a (Finset.mem_insert_self a s)).mul (ih hs)
    have hfun : (fun x => ∏ i ∈ insert a s, g i x) = (fun x => g a x * ∏ i ∈ s, g i x) := by
      funext x; exact Finset.prod_insert ha
    have hrest : ∑ j ∈ s, (∏ i ∈ (insert a s).erase j, g i t) * g' j t
        = g a t * ∑ j ∈ s, (∏ i ∈ s.erase j, g i t) * g' j t := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j hj => ?_
      have haj : a ≠ j := fun h => ha (h ▸ hj)
      have hna : a ∉ s.erase j := fun h => ha (Finset.mem_of_mem_erase h)
      rw [Finset.erase_insert_of_ne haj, Finset.prod_insert hna]
      ring
    rw [hfun]
    refine h.congr_deriv ?_
    rw [Finset.sum_insert ha, Finset.erase_insert ha, hrest]
    ring

omit [Fintype ι] in
theorem continuous_prod' (s : Finset ι) (g : ι → ℝ → ℝ) (hg : ∀ i ∈ s, Continuous (g i)) :
    Continuous (fun x => ∏ i ∈ s, g i x) := by
  revert hg
  refine Finset.induction_on s ?_ ?_
  · intro _; simp only [Finset.prod_empty]; exact continuous_const
  · intro a s ha ih hg
    simp only [Finset.prod_insert ha]
    exact (hg a (Finset.mem_insert_self a s)).mul (ih fun i hi => hg i (Finset.mem_insert_of_mem hi))

theorem continuous_sum' {κ : Type*} [DecidableEq κ] (s : Finset κ) (g : κ → ℝ → ℝ)
    (hg : ∀ i ∈ s, Continuous (g i)) : Continuous (fun x => ∑ i ∈ s, g i x) := by
  revert hg
  refine Finset.induction_on s ?_ ?_
  · intro _; simp only [Finset.sum_empty]; exact continuous_const
  · intro a s ha ih hg
    simp only [Finset.sum_insert ha]
    exact (hg a (Finset.mem_insert_self a s)).add (ih fun i hi => hg i (Finset.mem_insert_of_mem hi))

/-- The integral over [0, 1] of a finite sum of continuous functions is the sum of the integrals. -/
theorem integral_sum' {κ : Type*} [DecidableEq κ] (s : Finset κ) (f : κ → ℝ → ℝ)
    (hf : ∀ i ∈ s, Continuous (f i)) :
    ∫ x in (0 : ℝ)..1, ∑ i ∈ s, f i x = ∑ i ∈ s, ∫ x in (0 : ℝ)..1, f i x := by
  revert hf
  refine Finset.induction_on s ?_ ?_
  · intro _; simp
  · intro a s ha ih hf
    have hs : ∀ i ∈ s, Continuous (f i) := fun i hi => hf i (Finset.mem_insert_of_mem hi)
    simp only [Finset.sum_insert ha]
    rw [intervalIntegral.integral_add ((hf a (Finset.mem_insert_self a s)).intervalIntegrable _ _)
      ((continuous_sum' s f hs).intervalIntegrable _ _), ih hs]

/-! ### The theorems -/

/-- The derivative of the total along a path is the sum of the causes' instantaneous rates. -/
theorem hasDerivAt_total (p p' : ι → ℝ → ℝ) (t : ℝ)
    (hp : ∀ i, HasDerivAt (p i) (p' i t) t) :
    HasDerivAt (fun s => total (fun i => p i s))
      (∑ j, (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) t := by
  have h := (hasDerivAt_prod' univ (fun i x => 1 - p i x) (fun i x => -p' i x) t
    (fun i _ => (hp i).const_sub 1)).const_sub 1
  unfold total
  refine h.congr_deriv ?_
  simp only [mul_neg, Finset.sum_neg_distrib, neg_neg]

/-- Path independence: whatever the timing, the dependent rates add up to the total. -/
theorem path_total (p p' : ι → ℝ → ℝ)
    (hp : ∀ i t, HasDerivAt (p i) (p' i t) t)
    (hc : ∀ i, Continuous (p' i)) :
    ∑ j, rate p p' j = total (fun i => p i 1) - total (fun i => p i 0) := by
  have hcp : ∀ i, Continuous (p i) := fun i =>
    continuous_iff_continuousAt.2 fun t => (hp i t).continuousAt
  have hterm : ∀ j ∈ (univ : Finset ι),
      Continuous (fun t => (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) := fun j _ =>
    (continuous_prod' (univ.erase j) (fun i t => 1 - p i t)
      fun i _ => continuous_const.sub (hcp i)).mul (hc j)
  unfold rate
  rw [← integral_sum' univ (fun j t => (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) hterm]
  have hd : ∀ t ∈ Set.uIcc (0 : ℝ) 1, HasDerivAt (fun s => total (fun i => p i s))
      (∑ j, (∏ i ∈ univ.erase j, (1 - p i t)) * p' j t) t :=
    fun t _ => hasDerivAt_total p p' t (fun i => hp i t)
  exact intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    ((continuous_sum' _ _ hterm).intervalIntegrable _ _)

/-- UDD: on the straight path pᵢ t = t qᵢ, cause j gets every corner S ∪ {j} of the cube,
    worth (−1)^|S| qⱼ ∏_S qᵢ, divided equally among its |S| + 1 members. -/
theorem udd_corner (q : ι → ℝ) (j : ι) :
    rate (fun i t => t * q i) (fun i _ => q i) j
      = ∑ S ∈ (univ.erase j).powerset,
          (-1) ^ S.card * q j * (∏ i ∈ S, q i) / (S.card + 1) := by
  unfold rate
  show ∫ t in (0 : ℝ)..1, (∏ i ∈ univ.erase j, (1 - t * q i)) * q j = _
  -- expand the product over the other causes into corners
  have hexp : ∀ t : ℝ, (∏ i ∈ univ.erase j, (1 - t * q i)) * q j
      = ∑ S ∈ (univ.erase j).powerset, ((-1) ^ S.card * q j * ∏ i ∈ S, q i) * t ^ S.card := by
    intro t
    have h1 : ∀ i, (1 : ℝ) - t * q i = (-t) * q i + 1 := fun i => by ring
    simp_rw [h1, Finset.prod_add, Finset.prod_const_one, mul_one, Finset.sum_mul]
    refine Finset.sum_congr rfl fun S _ => ?_
    rw [Finset.prod_mul_distrib, Finset.prod_const, neg_pow]
    ring
  simp_rw [hexp]
  rw [integral_sum' _ (fun S x => ((-1) ^ S.card * q j * ∏ i ∈ S, q i) * x ^ S.card)
    (fun S _ => by fun_prop)]
  refine Finset.sum_congr rfl fun S _ => ?_
  rw [intervalIntegral.integral_const_mul, integral_pow]
  simp only [one_pow, zero_pow (by omega : S.card + 1 ≠ 0), sub_zero]
  ring

/-- Constant force: with pᵢ t = 1 − e^{−μᵢ t}, cause j gets μⱼ (1 − e^{−M}) / M, M = Σ μ,
    which is the share μⱼ / M of the total 1 − e^{−M}. -/
theorem const_force_split (μ : ι → ℝ) (j : ι) (hM : ∑ i, μ i ≠ 0) :
    rate (fun i t => 1 - Real.exp (-(μ i * t))) (fun i t => μ i * Real.exp (-(μ i * t))) j
      = μ j * (1 - Real.exp (-∑ i, μ i)) / ∑ i, μ i := by
  unfold rate
  show ∫ t in (0 : ℝ)..1, (∏ i ∈ univ.erase j, (1 - (1 - Real.exp (-(μ i * t)))))
      * (μ j * Real.exp (-(μ j * t))) = _
  have hpt : ∀ t : ℝ, (∏ i ∈ univ.erase j, (1 - (1 - Real.exp (-(μ i * t)))))
      * (μ j * Real.exp (-(μ j * t))) = μ j * Real.exp (-((∑ i, μ i) * t)) := by
    intro t
    simp only [sub_sub_cancel]
    rw [← Real.exp_sum]
    have hsplit : (∑ i, μ i) * t = μ j * t + ∑ i ∈ univ.erase j, μ i * t := by
      rw [Finset.sum_mul, Finset.add_sum_erase _ (fun i => μ i * t) (Finset.mem_univ j)]
    rw [hsplit, neg_add, Real.exp_add, Finset.sum_neg_distrib]
    ring
  simp_rw [hpt]
  rw [intervalIntegral.integral_const_mul]
  have hd : ∀ t ∈ Set.uIcc (0 : ℝ) 1, HasDerivAt
      (fun s => -Real.exp (-((∑ i, μ i) * s)) / ∑ i, μ i) (Real.exp (-((∑ i, μ i) * t))) t := by
    intro t _
    have h := (((hasDerivAt_id' t).const_mul (∑ i, μ i)).neg.exp.neg).div_const (∑ i, μ i)
    refine h.congr_deriv ?_
    try simp only [Pi.neg_apply, id]
    rw [mul_one, mul_neg, neg_neg, mul_div_assoc, div_self hM, mul_one]
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp only [mul_one, mul_zero, neg_zero, Real.exp_zero]
  ring

/-- The total on the constant-force path is 1 − e^{−M}: so the split is in proportion to the μ. -/
omit [DecidableEq ι] in
theorem const_force_total (μ : ι → ℝ) :
    total (fun i => 1 - Real.exp (-(μ i * 1))) = 1 - Real.exp (-∑ i, μ i) := by
  unfold total
  simp only [mul_one, sub_sub_cancel, ← Real.exp_sum, sum_neg_distrib]

end DecrementPaths
