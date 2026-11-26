/-
  Depth is derived (DRAFT): the steps behind "Wythoff depth is not a new dial".

  A kaleidoscope shape's support function is h(u) = max over the group of ⟨w·seed, u⟩. When the seed and
  the direction share the chamber, the identity wins that maximum for every such seed (Kostant convexity).
  The common-winner lemma below then makes the support function add, so the shape at depth k + f is the
  Minkowski sum (1 - f)·corner_k + f·corner_(k+1): a weighted sum of shapes already in the model.

  Proved here, with no `sorry`:
    * `argmax_add`: one index that wins at x and at y wins at x + y, and the winning values add.
    * `argmax_smul_add`: the same with non-negative weights, the Minkowski-sum step.
    * `chance_complement`: the depth fraction is the chance t/(1 + t) of an odds t, and the odds
      complement 1/t gives 1 - (chance), the mirror of the depth path.
    * `odds_rescale`: resizing the two corner shapes by α and β multiplies the odds by α/β and changes
      nothing else (the weighted sum is the same shape up to size).
    * `face_to_temperature`: the simplex gap is exact in one limit only: (1 + y/a)^(β a) → exp(β y)
      as the offset a grows, which turns face power into the temperature (log-sum-exp) route.
  Numerical companions: shape-dials/depth_checks.py (20 checks).
-/
import Mathlib

open Filter Topology

namespace MathlibMap.DepthDerived

variable {ι E : Type*} [AddCommGroup E] [Module ℝ E]

/-- Common-winner lemma: if one index wins the maximum at `x` and at `y`, it wins at `x + y`,
and the winning value adds. -/
theorem argmax_add (f : ι → E →ₗ[ℝ] ℝ) (x y : E) (i₀ : ι)
    (hx : ∀ i, f i x ≤ f i₀ x) (hy : ∀ i, f i y ≤ f i₀ y) :
    (∀ i, f i (x + y) ≤ f i₀ (x + y)) ∧ f i₀ (x + y) = f i₀ x + f i₀ y := by
  refine ⟨fun i => ?_, map_add _ _ _⟩
  rw [map_add, map_add]
  exact add_le_add (hx i) (hy i)

/-- The Minkowski-sum step: with non-negative weights `l`, `m`, the common winner still wins at
`l • x + m • y`, and the support value is `l * h(x) + m * h(y)`. -/
theorem argmax_smul_add (f : ι → E →ₗ[ℝ] ℝ) (x y : E) (i₀ : ι) {l m : ℝ} (hl : 0 ≤ l) (hm : 0 ≤ m)
    (hx : ∀ i, f i x ≤ f i₀ x) (hy : ∀ i, f i y ≤ f i₀ y) :
    (∀ i, f i (l • x + m • y) ≤ f i₀ (l • x + m • y)) ∧
      f i₀ (l • x + m • y) = l * f i₀ x + m * f i₀ y := by
  refine ⟨fun i => ?_, by simp only [map_add, map_smul, smul_eq_mul]⟩
  simp only [map_add, map_smul, smul_eq_mul]
  exact add_le_add (mul_le_mul_of_nonneg_left (hx i) hl) (mul_le_mul_of_nonneg_left (hy i) hm)

/-- The odds complement `1/t` gives the complementary chance: the mirror of the depth path. -/
theorem chance_complement {t : ℝ} (ht : 0 < t) : (1 / t) / (1 + 1 / t) = 1 - t / (1 + t) := by
  have h0 : t ≠ 0 := ht.ne'
  have h1 : 1 + t ≠ 0 := by positivity
  have h2 : t + 1 ≠ 0 := by positivity
  field_simp
  ring

/-- Resizing the corner shapes by `α`, `β` > 0 is undone by multiplying the odds by `α/β`: the effective
weights on the original corners keep the ratio `t`, so the shape is the same up to size. -/
theorem odds_rescale {t α β : ℝ} (ht : 0 < t) (hα : 0 < α) (hβ : 0 < β) :
    (t * (α / β) / (1 + t * (α / β))) * β / ((1 / (1 + t * (α / β))) * α) = t := by
  have hα' : α ≠ 0 := hα.ne'
  have hβ' : β ≠ 0 := hβ.ne'
  have h1 : 1 + t * (α / β) ≠ 0 := by positivity
  field_simp

/-- The simplex gap is exact in one limit only: as the offset `a` grows with `p = β a`, each face's
power term tends to the exponential, which is the temperature route. -/
theorem face_to_temperature (y β : ℝ) :
    Tendsto (fun a : ℝ => (1 + y / a) ^ (β * a)) atTop (𝓝 (Real.exp (β * y))) := by
  have h0 : Tendsto (fun x : ℝ => (1 + y / x) ^ x) atTop (𝓝 (Real.exp y)) := by
    first
    | exact Real.tendsto_one_add_div_rpow_exp y
    | exact tendsto_one_add_div_rpow_exp y
    | exact Real.tendsto_one_plus_div_rpow_exp y
  have h : Tendsto (fun a : ℝ => ((1 + y / a) ^ a) ^ β) atTop (𝓝 (Real.exp y ^ β)) :=
    h0.rpow_const (Or.inl (Real.exp_pos y).ne')
  rw [← Real.exp_mul, mul_comm y β] at h
  refine h.congr' ?_
  have hlim : Tendsto (fun a : ℝ => 1 + y / a) atTop (𝓝 1) := by
    simpa using (tendsto_const_nhds (x := (1 : ℝ))).add
      ((tendsto_const_nhds (x := y)).div_atTop tendsto_id)
  have hev : ∀ᶠ a : ℝ in atTop, 0 ≤ 1 + y / a := hlim.eventually (le_mem_nhds zero_lt_one)
  filter_upwards [hev] with a ha
  rw [← Real.rpow_mul ha, mul_comm a β]

end MathlibMap.DepthDerived
