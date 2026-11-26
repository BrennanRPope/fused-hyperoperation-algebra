/-
  Klein-bottle flow (DRAFT): what a twist does to the bookkeeping of a 2D fluid.

  On the Klein bottle K = T² / g, with g(x, y) = (x + π, −y) a glide reflection, the streamfunction
  and the vorticity are odd under g (pseudo-scalars). Discretised on a grid, a field on K is a function
  `f` on a finite set with `f (σ x) = − f x` for the grid version `σ` of g. Proved here, with no `sorry`:

    * `sum_eq_zero_of_odd`: the total of a twist-odd quantity over the cover is 0. So the total
      vorticity ("total spin") carries no information on K: it is not a quantity at all.
    * `sum_comp_eq_zero_of_odd`: the same for any odd function of it (ω³, ω⁵, …): every odd vorticity
      moment vanishes, so on K the vorticity distribution is forced to be symmetric.
    * `uniform_flow_klein_iff`: a uniform flow (U, V) survives the twist (dg = diag(1, −1)) only if
      V = 0: one steady sweeping flow on K, against two on a torus.
    * `relation_odd_on_klein`: if a steady flow has vorticity a function of its streamfunction,
      ω = F ∘ ψ, then on K that function is odd on the values ψ takes.
    * `lone_vortex_maximises`: for ψ = A cos x + B sin y (the lowest shell on K), the quartic moment is
      3/8 (A² + B²)² + 3/4 A² B²; at fixed A² + B² the mixed part A² B² is largest when A² = B², the
      round lone vortex, and smallest at a pure shear (A = 0 or B = 0).
  Numerical companions: klein-flow/checks.py.
-/
import Mathlib

namespace MathlibMap.KleinFlow

open Finset

/-- The total of a quantity that flips sign under an involution-like symmetry is zero. -/
theorem sum_eq_zero_of_odd {α : Type*} [Fintype α] (σ : α ≃ α) (f : α → ℝ)
    (h : ∀ x, f (σ x) = -f x) : ∑ x, f x = 0 := by
  have h1 := Equiv.sum_comp σ f
  simp only [h, Finset.sum_neg_distrib] at h1
  linarith

/-- Every odd function of a twist-odd field also sums to zero (odd vorticity moments vanish). -/
theorem sum_comp_eq_zero_of_odd {α : Type*} [Fintype α] (σ : α ≃ α) (f : α → ℝ)
    (h : ∀ x, f (σ x) = -f x) (F : ℝ → ℝ) (hF : ∀ a, F (-a) = -F a) :
    ∑ x, F (f x) = 0 :=
  sum_eq_zero_of_odd σ (fun x => F (f x)) (fun x => by simp only [h, hF])

/-- A uniform flow is unchanged by the twist's derivative diag(1, −1) exactly when V = 0. -/
theorem uniform_flow_klein_iff (U V : ℝ) : ((U, -V) : ℝ × ℝ) = (U, V) ↔ V = 0 := by
  constructor
  · intro h
    have h2 : -V = V := congrArg Prod.snd h
    linarith
  · intro h
    simp [h]

/-- On the Klein bottle a steady vorticity–streamfunction relation is odd on the values ψ takes. -/
theorem relation_odd_on_klein {α : Type*} (σ : α → α) (ψ ω : α → ℝ) (F : ℝ → ℝ)
    (hψ : ∀ x, ψ (σ x) = -ψ x) (hω : ∀ x, ω (σ x) = -ω x) (hF : ∀ x, ω x = F (ψ x)) (x : α) :
    F (-ψ x) = -F (ψ x) := by
  rw [← hψ x, ← hF (σ x), hω x, hF x]

/-- Quartic moment of A cos + B sin from the averages ⟨c⁴⟩ = ⟨s⁴⟩ = 3/8, ⟨c²s²⟩ = 1/4. -/
theorem quartic_moment (A B : ℝ) :
    A ^ 4 * (3 / 8) + 6 * A ^ 2 * B ^ 2 * (1 / 4) + B ^ 4 * (3 / 8)
      = 3 / 8 * (A ^ 2 + B ^ 2) ^ 2 + 3 / 4 * (A ^ 2 * B ^ 2) := by
  ring

/-- At fixed A² + B², the mixed part A² B² is at most a quarter of (A² + B²)², with equality iff
    A² = B² (the round lone vortex). -/
theorem lone_vortex_maximises (A B : ℝ) :
    A ^ 2 * B ^ 2 ≤ (A ^ 2 + B ^ 2) ^ 2 / 4 ∧
      (A ^ 2 * B ^ 2 = (A ^ 2 + B ^ 2) ^ 2 / 4 ↔ A ^ 2 = B ^ 2) := by
  constructor
  · nlinarith [sq_nonneg (A ^ 2 - B ^ 2)]
  · constructor
    · intro h
      have h2 : (A ^ 2 - B ^ 2) * (A ^ 2 - B ^ 2) = 0 := by nlinarith
      have h3 : A ^ 2 - B ^ 2 = 0 := mul_self_eq_zero.mp h2
      linarith
    · intro h
      rw [h]
      ring

end MathlibMap.KleinFlow
