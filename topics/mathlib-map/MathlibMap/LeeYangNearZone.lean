/-
  Lee–Yang near zone on the shared core (DRAFT, stage 1 of moving the near-zone check into Lean).

  What this file proves, in Lean, with no `sorry`:
    1. `taylor_cell`: the Taylor-type bound the box check uses on every cell of width h in s. If the (n + 1)-st
       derivative of g stays under M on [x₀, x₀ + h], then |g t − g x₀| is at most the derivatives at x₀ times
       powers of (t − x₀) over factorials, plus M (t − x₀)^(n+1)/(n+1)!. Proved from the comparison principle
       (`image_norm_le_of_norm_deriv_right_le_deriv_boundary`) by induction on n.
    2. `RunCert.nearAlt`: a "run certificate" (what the box check finds: runs where f changes sign with a
       one-signed slope, gaps where f has no zero, and a sign of g at the run's zero that flips from run to run)
       gives alternation of g at the zeros of f, up to the end of the last run.
    3. `alternates_of_near_far`: near-zone alternation up to Y plus far-zone alternation past T, with a zero of f
       between T and Y, give `Alternates p b` (the core's per-step input) at every positive zero.
    4. `rule_range`, `holds_range`: with the far-zone and near-zone certificates for 1.1 ≤ p ≤ 5, the rule holds at
       every such p (and the holding half at p = 4, which is even).

  What stays cited (as named hypotheses, not proved here):
    * `Toolkit` (Pólya, Krein, the break), as in LeeYangCore.lean;
    * `FarCert`: far_bounds_p.py, the far zone past T*(p) for 1.1 ≤ p ≤ 8;
    * `NearCert`: near_zone_p.py, a run certificate for every (p, b) with 1.1 ≤ p ≤ 5, 0 < b ≤ p − 1.
  Stage 2 replaces `NearCert` with a Lean computation that builds the run certificates from exact numbers.
-/
import MathlibMap.LeeYangCore
import Mathlib

namespace MathlibMap.LeeYang

open Set

/-! ### 1. The Taylor bound on one cell -/

/-- `g` has derivatives of every order. -/
def AllDiff (g : ℝ → ℝ) : Prop := ∀ k : ℕ, Differentiable ℝ (iteratedDeriv k g)

lemma AllDiff.deriv' {g : ℝ → ℝ} (hg : AllDiff g) : AllDiff (deriv g) := by
  intro k
  have h := hg (k + 1)
  rwa [iteratedDeriv_succ'] at h

lemma AllDiff.differentiable {g : ℝ → ℝ} (hg : AllDiff g) : Differentiable ℝ g := by
  simpa using hg 0

/-- Derivative of one term c (t − x₀)^(m+1)/(m+1)!. -/
lemma hasDerivAt_term (c x₀ t : ℝ) (m : ℕ) :
    HasDerivAt (fun t => c * (t - x₀) ^ (m + 1) / ((m + 1).factorial : ℝ))
      (c * (t - x₀) ^ m / (m.factorial : ℝ)) t := by
  have h1 : HasDerivAt (fun t => t - x₀) 1 t := (hasDerivAt_id t).sub_const x₀
  have h2 := (h1.pow (m + 1)).const_mul c
  have h3 := h2.div_const ((m + 1).factorial : ℝ)
  convert h3 using 1
  rw [Nat.add_sub_cancel, Nat.factorial_succ]
  push_cast
  rw [div_eq_div_iff (by positivity) (by positivity)]
  ring

/-- The bound polynomial: Σ_{m<n} c_m u^(m+1)/(m+1)! + M u^(n+1)/(n+1)!, with u = t − x₀. -/
noncomputable def tpoly (c : ℕ → ℝ) (M x₀ : ℝ) (n : ℕ) (t : ℝ) : ℝ :=
  (∑ m ∈ Finset.range n, c m * (t - x₀) ^ (m + 1) / ((m + 1).factorial : ℝ))
    + M * (t - x₀) ^ (n + 1) / ((n + 1).factorial : ℝ)

/-- Its derivative: Σ_{m<n} c_m u^m/m! + M u^n/n!. -/
noncomputable def tpoly' (c : ℕ → ℝ) (M x₀ : ℝ) (n : ℕ) (t : ℝ) : ℝ :=
  (∑ m ∈ Finset.range n, c m * (t - x₀) ^ m / (m.factorial : ℝ)) + M * (t - x₀) ^ n / (n.factorial : ℝ)

lemma hasDerivAt_tpoly (c : ℕ → ℝ) (M x₀ : ℝ) (n : ℕ) (t : ℝ) :
    HasDerivAt (tpoly c M x₀ n) (tpoly' c M x₀ n t) t := by
  unfold tpoly tpoly'
  apply HasDerivAt.add
  · have := HasDerivAt.sum (u := Finset.range n)
      (A := fun m t => c m * (t - x₀) ^ (m + 1) / ((m + 1).factorial : ℝ))
      (A' := fun m => c m * (t - x₀) ^ m / (m.factorial : ℝ)) (x := t)
      (fun m _ => hasDerivAt_term (c m) x₀ t m)
    exact this
  · exact hasDerivAt_term M x₀ t n

lemma tpoly_at (c : ℕ → ℝ) (M x₀ : ℝ) (n : ℕ) : tpoly c M x₀ n x₀ = 0 := by
  simp [tpoly]

/-- The step that reindexes the derivative's bound: tpoly' for g at order n + 1 is |g'(x₀)| plus tpoly for g'
at order n. -/
lemma tpoly'_succ (c : ℕ → ℝ) (M x₀ : ℝ) (n : ℕ) (t : ℝ) :
    tpoly' c M x₀ (n + 1) t = c 0 + tpoly (fun m => c (m + 1)) M x₀ n t := by
  unfold tpoly' tpoly
  rw [Finset.sum_range_succ']
  simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, div_one, mul_one]
  ring

/-- **Taylor bound on a cell.** If the (n+1)-st derivative of g stays under M on [x₀, x₀ + h], then
|g t − g x₀| ≤ Σ_{m<n} |g^(m+1)(x₀)| (t − x₀)^(m+1)/(m+1)! + M (t − x₀)^(n+1)/(n+1)! there. -/
theorem taylor_cell (n : ℕ) : ∀ (g : ℝ → ℝ), AllDiff g → ∀ (x₀ h M : ℝ),
    (∀ t ∈ Icc x₀ (x₀ + h), |iteratedDeriv (n + 1) g t| ≤ M) →
    ∀ t ∈ Icc x₀ (x₀ + h),
      |g t - g x₀| ≤ tpoly (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ n t := by
  induction n with
  | zero =>
    intro g hg x₀ h M hM t ht
    have hd : Differentiable ℝ g := hg.differentiable
    have key := image_norm_le_of_norm_deriv_right_le_deriv_boundary
      (f := fun t => g t - g x₀) (f' := fun t => deriv g t) (a := x₀) (b := x₀ + h)
      (B := tpoly (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ 0)
      (B' := tpoly' (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ 0)
      ((hd.continuous.sub continuous_const).continuousOn)
      (fun x _ => ((hd x).hasDerivAt.sub_const (g x₀)).hasDerivWithinAt)
      (by simp [tpoly_at])
      ((fun x _ => (hasDerivAt_tpoly _ M x₀ 0 x).continuousAt.continuousWithinAt))
      (fun x _ => (hasDerivAt_tpoly _ M x₀ 0 x).hasDerivWithinAt)
      (fun x hx => by
        have := hM x (Ico_subset_Icc_self hx)
        simpa [tpoly', Real.norm_eq_abs, iteratedDeriv_one] using this)
    simpa [Real.norm_eq_abs] using key ht
  | succ n ih =>
    intro g hg x₀ h M hM t ht
    have hd : Differentiable ℝ g := hg.differentiable
    -- the bound for g' from the induction hypothesis
    have hg' : ∀ x ∈ Icc x₀ (x₀ + h),
        |deriv g x| ≤ tpoly' (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ (n + 1) x := by
      intro x hx
      have h1 := ih (deriv g) hg.deriv' x₀ h M
        (fun y hy => by simpa [iteratedDeriv_succ'] using hM y hy) x hx
      rw [tpoly'_succ]
      have e : (fun m => |iteratedDeriv (m + 1) (deriv g) x₀|)
          = (fun m => |iteratedDeriv (m + 1 + 1) g x₀|) := by
        ext m; rw [iteratedDeriv_succ' (n := m + 1)]
      rw [e] at h1
      have h2 : |deriv g x₀| = |iteratedDeriv (0 + 1) g x₀| := by simp [iteratedDeriv_one]
      calc |deriv g x| ≤ |deriv g x₀| + |deriv g x - deriv g x₀| := by
            have := abs_sub_abs_le_abs_sub (deriv g x) (deriv g x₀); linarith
        _ ≤ _ := by rw [h2]; linarith
    have key := image_norm_le_of_norm_deriv_right_le_deriv_boundary
      (f := fun t => g t - g x₀) (f' := fun t => deriv g t) (a := x₀) (b := x₀ + h)
      (B := tpoly (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ (n + 1))
      (B' := tpoly' (fun m => |iteratedDeriv (m + 1) g x₀|) M x₀ (n + 1))
      ((hd.continuous.sub continuous_const).continuousOn)
      (fun x _ => ((hd x).hasDerivAt.sub_const (g x₀)).hasDerivWithinAt)
      (by simp [tpoly_at])
      ((fun x _ => (hasDerivAt_tpoly _ M x₀ (n + 1) x).continuousAt.continuousWithinAt))
      (fun x _ => (hasDerivAt_tpoly _ M x₀ (n + 1) x).hasDerivWithinAt)
      (fun x hx => by simpa [Real.norm_eq_abs] using hg' x (Ico_subset_Icc_self hx))
    simpa [Real.norm_eq_abs] using key ht

/-! ### 2. Run certificates: what the box check finds -/

/-- f is nonzero on the closed interval. -/
def NoZero (f : ℝ → ℝ) (a b : ℝ) : Prop := ∀ t ∈ Icc a b, f t ≠ 0

/-- A run certificate for f (G at b) and g (G at b + 1): N runs [x i, y i] in increasing order. Before the first
run and in every gap f has no zero; across each run f changes sign and its slope keeps one sign; at any zero of f in
run i, g has the sign σ i, and the signs flip from run to run. -/
structure RunCert (f g : ℝ → ℝ) where
  N : ℕ
  x : ℕ → ℝ
  y : ℕ → ℝ
  σ : ℕ → ℝ
  hN : 0 < N
  hx0 : 0 < x 0
  hrun : ∀ i < N, x i < y i
  hord : ∀ i, i + 1 < N → y i ≤ x (i + 1)
  hstart : ∀ t, 0 < t → t ≤ x 0 → f t ≠ 0
  hgap : ∀ i, i + 1 < N → NoZero f (y i) (x (i + 1))
  hends : ∀ i < N, f (x i) * f (y i) < 0
  hslope : ∀ i < N, (∀ t ∈ Icc (x i) (y i), 0 < deriv f t) ∨ (∀ t ∈ Icc (x i) (y i), deriv f t < 0)
  hsign : ∀ i < N, ∀ t ∈ Icc (x i) (y i), f t = 0 → 0 < σ i * g t
  halt : ∀ i, i + 1 < N → σ i * σ (i + 1) < 0

/-- Alternation of g at the zeros of f, up to Y. -/
def NearAlt (f g : ℝ → ℝ) (Y : ℝ) : Prop :=
  (∀ t : ℝ, 0 < t → t ≤ Y → f t = 0 → deriv f t ≠ 0) ∧
  ∀ s t : ℝ, 0 < s → s < t → t ≤ Y → f s = 0 → f t = 0 → (∀ u : ℝ, s < u → u < t → f u ≠ 0) →
    g s * g t < 0

/-- Alternation of g at the zeros of f, from T on. -/
def FarAlt (f g : ℝ → ℝ) (T : ℝ) : Prop :=
  (∀ t : ℝ, T ≤ t → f t = 0 → deriv f t ≠ 0) ∧
  ∀ s t : ℝ, T ≤ s → s < t → f s = 0 → f t = 0 → (∀ u : ℝ, s < u → u < t → f u ≠ 0) →
    g s * g t < 0

namespace RunCert

variable {f g : ℝ → ℝ} (R : RunCert f g)

lemma y_le_x {i j : ℕ} (hij : i < j) (hj : j < R.N) : R.y i ≤ R.x j := by
  induction j with
  | zero => omega
  | succ j ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.mp hij with h | h
    · have h1 := ih h (by omega)
      have h2 := R.hrun j (by omega)
      have h3 := R.hord j (by omega)
      linarith
    · subst h; exact R.hord i (by omega)

lemma x_lt_x {i j : ℕ} (hij : i < j) (hj : j < R.N) : R.x i < R.x j :=
  lt_of_lt_of_le (R.hrun i (by omega)) (R.y_le_x hij hj)

/-- f is strictly monotone or strictly antitone on each run. -/
lemma injOn_run (hf : Differentiable ℝ f) {i : ℕ} (hi : i < R.N) : InjOn f (Icc (R.x i) (R.y i)) := by
  rcases R.hslope i hi with h | h
  · exact (strictMonoOn_of_deriv_pos (convex_Icc _ _) hf.continuous.continuousOn
      (fun t ht => h t (interior_subset ht))).injOn
  · exact (strictAntiOn_of_deriv_neg (convex_Icc _ _) hf.continuous.continuousOn
      (fun t ht => h t (interior_subset ht))).injOn

/-- Each run holds a zero of f. -/
lemma exists_zero (hf : Differentiable ℝ f) {i : ℕ} (hi : i < R.N) :
    ∃ z ∈ Ioo (R.x i) (R.y i), f z = 0 := by
  have hlt := R.hrun i hi
  have hc : ContinuousOn f (Icc (R.x i) (R.y i)) := hf.continuous.continuousOn
  rcases mul_neg_iff.mp (R.hends i hi) with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · obtain ⟨z, hz, hz0⟩ := intermediate_value_Ioo' hlt.le hc ⟨h2, h1⟩
    exact ⟨z, hz, hz0⟩
  · obtain ⟨z, hz, hz0⟩ := intermediate_value_Ioo hlt.le hc ⟨h1, h2⟩
    exact ⟨z, hz, hz0⟩

/-- Every zero of f in (0, y (N − 1)] lies in some run. -/
lemma zero_in_run {t : ℝ} (ht0 : 0 < t) (htY : t ≤ R.y (R.N - 1)) (hz : f t = 0) :
    ∃ i < R.N, t ∈ Icc (R.x i) (R.y i) := by
  -- the last run whose left end is at most t
  have hex : ∃ i, i < R.N ∧ R.x i ≤ t := by
    refine ⟨0, R.hN, ?_⟩
    by_contra h; push_neg at h
    exact R.hstart t ht0 h.le hz
  classical
  set S := (Finset.range R.N).filter (fun i => R.x i ≤ t) with hSdef
  have hS : S.Nonempty := by
    obtain ⟨i, hi, hxi⟩ := hex
    exact ⟨i, by simp [hSdef, hi, hxi]⟩
  set i := S.max' hS with hidef
  have hiS : i ∈ S := S.max'_mem hS
  rw [hSdef, Finset.mem_filter, Finset.mem_range] at hiS
  refine ⟨i, hiS.1, hiS.2, ?_⟩
  by_contra hty; push_neg at hty
  by_cases hlast : i + 1 < R.N
  · -- t is past run i; it isn't in the gap, so t ≥ x (i + 1), against the choice of i
    have hnext : R.x (i + 1) ≤ t := by
      by_contra h; push_neg at h
      exact R.hgap i hlast t ⟨hty.le, h.le⟩ hz
    have hmem : i + 1 ∈ S := by simp [hSdef, hlast, hnext]
    have := S.le_max' _ hmem
    rw [← hidef] at this
    omega
  · have : i = R.N - 1 := by have := hiS.1; omega
    rw [this] at hty
    linarith

/-- **A run certificate gives near-zone alternation** up to the end of the last run. -/
theorem nearAlt (hf : Differentiable ℝ f) : NearAlt f g (R.y (R.N - 1)) := by
  constructor
  · intro t ht0 htY hz
    obtain ⟨i, hi, hti⟩ := R.zero_in_run ht0 htY hz
    rcases R.hslope i hi with h | h
    · exact (h t hti).ne'
    · exact (h t hti).ne
  · intro s t hs0 hst htY hfs hft hbetween
    obtain ⟨j, hj, htj⟩ := R.zero_in_run (hs0.trans hst) htY hft
    obtain ⟨i, hi, hsi⟩ := R.zero_in_run hs0 (hst.le.trans htY) hfs
    -- i < j: same run is impossible (f is injective there), and a later run can't hold the smaller zero
    have hij : i < j := by
      rcases lt_trichotomy i j with h | h | h
      · exact h
      · subst h
        exact absurd (R.injOn_run hf hi hsi htj (by rw [hfs, hft])) hst.ne
      · have := R.y_le_x h hi
        have : t ≤ s := le_trans htj.2 (le_trans this hsi.1)
        linarith
    -- j = i + 1: otherwise run i + 1 holds a zero strictly between s and t
    have hj1 : j = i + 1 := by
      by_contra hne
      have hlt : i + 1 < j := by omega
      obtain ⟨z, hz, hz0⟩ := R.exists_zero hf (show i + 1 < R.N by omega)
      have h1 : s < z := lt_of_le_of_lt (le_trans hsi.2 (R.hord i (by omega))) hz.1
      have h2 : z < t := lt_of_lt_of_le hz.2 (le_trans (R.y_le_x hlt hj) htj.1)
      exact hbetween z h1 h2 hz0
    subst hj1
    have hσs := R.hsign i hi s hsi hfs
    have hσt := R.hsign (i + 1) hj t htj hft
    have hσ := R.halt i hj
    -- σ_i g s > 0, σ_{i+1} g t > 0 and σ_i σ_{i+1} < 0 force g s g t < 0
    have hprod : 0 < (R.σ i * g s) * (R.σ (i + 1) * g t) := mul_pos hσs hσt
    have : (R.σ i * g s) * (R.σ (i + 1) * g t) = (R.σ i * R.σ (i + 1)) * (g s * g t) := by ring
    rw [this] at hprod
    by_contra hge; push_neg at hge
    have : (R.σ i * R.σ (i + 1)) * (g s * g t) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg hσ.le hge
    linarith

end RunCert

/-! ### 3. Near zone + far zone = alternation at every zero -/

/-- Near-zone alternation up to Y, far-zone alternation from T, and a zero of f between T and Y give alternation at
every positive zero. -/
theorem alt_of_near_far {f g : ℝ → ℝ} {T Y : ℝ} (hn : NearAlt f g Y) (hfar : FarAlt f g T)
    (hz : ∃ z, T ≤ z ∧ z ≤ Y ∧ f z = 0) :
    (∀ t : ℝ, 0 < t → f t = 0 → deriv f t ≠ 0) ∧
    ∀ s t : ℝ, 0 < s → s < t → f s = 0 → f t = 0 → (∀ u : ℝ, s < u → u < t → f u ≠ 0) →
      g s * g t < 0 := by
  obtain ⟨z, hTz, hzY, hz0⟩ := hz
  constructor
  · intro t ht0 htz
    by_cases hT : T ≤ t
    · exact hfar.1 t hT htz
    · exact hn.1 t ht0 (by linarith) htz
  · intro s t hs0 hst hfs hft hb
    by_cases hT : T ≤ s
    · exact hfar.2 s t hT hst hfs hft hb
    · push_neg at hT
      by_cases hY : t ≤ Y
      · exact hn.2 s t hs0 hst hY hfs hft hb
      · push_neg at hY
        exact absurd hz0 (hb z (by linarith) (by linarith))

/-- On the core's terms: near + far give `Alternates p b`. -/
theorem alternates_of_near_far {p b T Y : ℝ} (hn : NearAlt (G p b) (G p (b + 1)) Y)
    (hfar : FarAlt (G p b) (G p (b + 1)) T) (hz : ∃ z, T ≤ z ∧ z ≤ Y ∧ G p b z = 0) :
    Alternates p b :=
  alt_of_near_far hn hfar hz

/-! ### 4. The rule for every roundness from 1.1 to 5 -/

/-- far_bounds_p.py (cited): past T*(p), the step b → b + 1 alternates, for 1.1 ≤ p ≤ 8 and 0 < b ≤ p − 1. -/
def FarCert (Tstar : ℝ → ℝ) : Prop :=
  ∀ p b : ℝ, 11 / 10 ≤ p → p ≤ 8 → 0 < b → b ≤ p - 1 → FarAlt (G p b) (G p (b + 1)) (Tstar p)

/-- near_zone_p.py (cited): for 1.1 ≤ p ≤ 5 and 0 < b ≤ p − 1, a run certificate whose last run reaches a zero of
G at b past T*(p). G at b is smooth in s (it is entire). -/
def NearCert (Tstar : ℝ → ℝ) : Prop :=
  ∀ p b : ℝ, 11 / 10 ≤ p → p ≤ 5 → 0 < b → b ≤ p - 1 →
    Differentiable ℝ (G p b) ∧
    ∃ R : RunCert (G p b) (G p (b + 1)), ∃ z, Tstar p ≤ z ∧ z ≤ R.y (R.N - 1) ∧ G p b z = 0

theorem alternatesBelow_range {Tstar : ℝ → ℝ} (hF : FarCert Tstar) (hN : NearCert Tstar) {p : ℝ}
    (hp1 : 11 / 10 ≤ p) (hp5 : p ≤ 5) : AlternatesBelow p := by
  intro b hb0 hbp
  obtain ⟨hd, R, z, hz⟩ := hN p b hp1 hp5 hb0 hbp
  exact alternates_of_near_far (R.nearAlt hd) (hF p b hp1 (by linarith) hb0 hbp) ⟨z, hz⟩

/-- The holding half for every roundness 1.1 ≤ p ≤ 5, even p = 2 and 4 included: Lee–Yang for every a ≤ p. -/
theorem holds_range (T : Toolkit) {Tstar : ℝ → ℝ} (hF : FarCert Tstar) (hN : NearCert Tstar) :
    ∀ p a : ℝ, 11 / 10 ≤ p → p ≤ 5 → 0 ≤ a → a ≤ p → Holds p a :=
  fun p a hp1 hp5 ha0 hap =>
    holds_of_alternatesBelow T (by linarith) (alternatesBelow_range hF hN hp1 hp5) a ha0 hap

/-- **The rule for every roundness 1.1 ≤ p ≤ 5** that isn't an even whole number: Lee–Yang holds exactly when
a ≤ p (n ≤ 1 + p² for a uniform spin). -/
theorem rule_range (T : Toolkit) {Tstar : ℝ → ℝ} (hF : FarCert Tstar) (hN : NearCert Tstar) :
    ∀ p : ℝ, 11 / 10 ≤ p → p ≤ 5 → ¬ EvenWhole p → ∀ a : ℝ, 0 ≤ a → (Holds p a ↔ a ≤ p) :=
  fun p hp1 hp5 hne =>
    rule_of_alternatesBelow T (by linarith) hne (alternatesBelow_range hF hN hp1 hp5)

end MathlibMap.LeeYang
