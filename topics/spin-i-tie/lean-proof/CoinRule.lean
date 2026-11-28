/-
DRAFT. Round 15: the coin rule forces the dial to k = 0, written for Lean 4 with Mathlib.
COMPILED: the repo's CI (run 145) compiled this file against the pinned Mathlib with no errors.
The same steps are also checked symbolically in checks.py.
What it proves: the law's lift, that the law is commutative and associative, that Tsallis
entropy composes independent parts by the law with k = 1 - q, and that the coin rule
(an independent part leaves the leftover unchanged) fails for every q ≠ 1.
What it does not prove: Faddeev's theorem (honest-weight grouping gives Shannon), which is cited.
-/
import Mathlib

open Finset

/-- The law: a ∘ₖ b = a + b + k a b. -/
def law (k a b : ℝ) : ℝ := a + b + k * a * b

/-- The lift 1 + k a turns the law into multiplication. -/
theorem law_lift (k a b : ℝ) : 1 + k * law k a b = (1 + k * a) * (1 + k * b) := by
  unfold law; ring

theorem law_comm (k a b : ℝ) : law k a b = law k b a := by
  unfold law; ring

theorem law_assoc (k a b c : ℝ) : law k (law k a b) c = law k a (law k b c) := by
  unfold law; ring

/-- Tsallis entropy S_q(p) = (1 - Σ p_i^q) / (q - 1). -/
noncomputable def tsallis {ι : Type*} (s : Finset ι) (q : ℝ) (p : ι → ℝ) : ℝ :=
  (1 - ∑ i ∈ s, p i ^ q) / (q - 1)

/-- Independent parts compose by the law, with dial k = 1 - q. -/
theorem tsallis_indep {ι κ : Type*} (s : Finset ι) (t : Finset κ) (q : ℝ) (hq : q ≠ 1)
    (p : ι → ℝ) (r : κ → ℝ) (hp : ∀ i, 0 ≤ p i) (hr : ∀ j, 0 ≤ r j) :
    tsallis (s ×ˢ t) q (fun x => p x.1 * r x.2) = law (1 - q) (tsallis s q p) (tsallis t q r) := by
  have hsum : ∑ x ∈ s ×ˢ t, (p x.1 * r x.2) ^ q = (∑ i ∈ s, p i ^ q) * (∑ j ∈ t, r j ^ q) := by
    rw [Finset.sum_product, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => Real.mul_rpow (hp i) (hr j)))
  have h : q - 1 ≠ 0 := sub_ne_zero.mpr hq
  unfold tsallis law
  rw [hsum]
  field_simp
  ring

/-- The coin rule at the level of numbers: if adding a part with nonzero size leaves the
leftover unchanged, the dial is 0. -/
theorem coin_forces_k_zero (k sA sB : ℝ) (hA : sA ≠ 0) (hB : sB ≠ 0)
    (h : law k sA sB - sA = sB) : k = 0 := by
  unfold law at h
  have h0 : k * sA * sB = 0 := by linarith
  rcases mul_eq_zero.mp h0 with h1 | h1
  · rcases mul_eq_zero.mp h1 with h2 | h2
    · exact h2
    · exact absurd h2 hA
  · exact absurd h1 hB

/-- Every Tsallis entropy with q ≠ 1 breaks the coin rule. -/
theorem tsallis_breaks_coin_rule {ι κ : Type*} (s : Finset ι) (t : Finset κ) (q : ℝ) (hq : q ≠ 1)
    (p : ι → ℝ) (r : κ → ℝ) (hp : ∀ i, 0 ≤ p i) (hr : ∀ j, 0 ≤ r j)
    (hA : tsallis s q p ≠ 0) (hB : tsallis t q r ≠ 0) :
    tsallis (s ×ˢ t) q (fun x => p x.1 * r x.2) - tsallis s q p ≠ tsallis t q r := by
  intro h
  rw [tsallis_indep s t q hq p r hp hr] at h
  have hk := coin_forces_k_zero (1 - q) _ _ hA hB h
  exact hq (by linarith)
