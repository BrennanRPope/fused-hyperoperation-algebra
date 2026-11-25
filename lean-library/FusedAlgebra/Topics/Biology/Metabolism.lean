/-
  Fused Hyperoperation Algebra: metabolism (enzyme mechanisms, nutrient-limited growth).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Metabolism.lean

  1. Two-substrate enzymes. Writing u_A = K_a/A and u_B = K_b/B, the ordered (sequential)
     rate law is V/v = 1 + (u_A ✳_k u_B) with dial k = K_ia/K_a. Ping-pong is dial 0,
     random independent binding is dial 1 (plain ✳).
  2. Every double-reciprocal line passes through u_A = −1/k, the wall of the dial, where the
     rate stops depending on the second substrate.
  3. Two-nutrient growth: multiplicative Monod is dial 1, the additive-resistance law is dial 0.
     Total limitation combines by the "or" law (dial −1).
  4. Saturation: θ = S/(K + S) = 1 + inv(K/S), and the free fraction is −inv(K/S).
  LIBRARY NOTE: first compiled in the shared library. 2 exact twins removed (star_def,
  monod_mult); each is marked where it stood. No proof step or statement changed.
-/
import FusedAlgebra.Base

namespace Metabolism
open MultAdd

set_option linter.unusedSectionVars false

/-- Brennan's notation for mult-add in this file: a ✳ b = a + b + ab. -/
local infixl:65 " ✳ " => multAdd

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

theorem wAdd_comm' (k a b : R) : wAdd k a b = wAdd k b a := by unfold wAdd; grind

/-! ## 1. Two-substrate rate laws, written in reduced concentrations u = K/S -/

/-- The reduced bisubstrate rate law V/v as a polynomial in u_A, u_B with dial k. -/
def bisub (k uA uB : R) : R := 1 + uA + uB + k * uA * uB

theorem bisub_is_wAdd (k uA uB : R) : bisub k uA uB = 1 + wAdd k uA uB := by
  unfold bisub wAdd; grind

/-- Ping-pong (substituted enzyme): no product term, dial 0. Double-reciprocal lines are
    parallel. -/
theorem pingpong (uA uB : R) : bisub 0 uA uB = 1 + uA + uB := by unfold bisub; grind

/-- Random order with independent sites: dial 1, the rate law factors. -/
theorem random_independent (uA uB : R) : bisub 1 uA uB = (1 + uA) * (1 + uB) := by
  unfold bisub; grind

theorem random_independent' (uA uB : R) : bisub 1 uA uB = 1 + (uA ✳ uB) := by
  unfold bisub multAdd; grind

/-- Read as a line in u_A (a double-reciprocal plot at fixed B): intercept 1 + u_B,
    slope 1 + k·u_B. -/
theorem bisub_line (k uA uB : R) : bisub k uA uB = (1 + uB) + uA * (1 + k * uB) := by
  unfold bisub; grind

/-- Changing B moves V/v by (u_B − u_B')·(1 + k·u_A); the factor is the dial lift. -/
theorem bisub_shift (k uA uB uB' : R) :
    bisub k uA uB - bisub k uA uB' = (uB - uB') * (1 + k * uA) := by
  unfold bisub; grind

/-! ## 3. Two-nutrient growth -/

-- LIBRARY: monod_mult (Multiplicative Monod: 1/(θ₁θ₂) = (1 + u₁)(1 + u₂) = 1 + (u₁ ✳ u₂), dial 1) was an exact twin of ChemQuantum.fragments_multAdd and is removed.

/-- Additive resistance: 1/μ = 1/θ₁ + 1/θ₂ − 1 = 1 + u₁ + u₂, dial 0. -/
theorem monod_additive (u₁ u₂ : R) : (1 + u₁) + (1 + u₂) - 1 = 1 + wAdd 0 u₁ u₂ := by
  unfold wAdd; grind

/-- The shortfall from full growth, 1 − θ₁θ₂, is the "or" of the two shortfalls. -/
theorem limitation_or (θ₁ θ₂ : R) : 1 - θ₁ * θ₂ = wAdd (-1) (1 - θ₁) (1 - θ₂) := by
  unfold wAdd; grind

theorem limitation_or_three (θ₁ θ₂ θ₃ : R) :
    1 - θ₁ * θ₂ * θ₃ = wAdd (-1) (wAdd (-1) (1 - θ₁) (1 - θ₂)) (1 - θ₃) := by
  unfold wAdd; grind

/-- n identical limiting steps: the joint shortfall is a weighted power at dial −1. -/
theorem limitation_pow (x : R) (n : Nat) : 1 - (1 - x) ^ n = wpow (-1) x n := by
  have h := wpow_lift (-1 : R) x n
  have e : (1 : R) + -1 * x = 1 - x := by grind
  rw [e] at h; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- **The ordered mechanism is a weighted mult-add.** With u_A = K_a/A, u_B = K_b/B and
    k = K_ia/K_a, the textbook rate law
    V/v = 1 + K_a/A + K_b/B + K_ia·K_b/(A·B) equals 1 + (u_A ✳_k u_B). -/
theorem ordered_is_dial (Ka Kia Kb A B : K) (hKa : Ka ≠ 0) (hA : A ≠ 0) (hB : B ≠ 0) :
    1 + Ka / A + Kb / B + Kia * Kb / (A * B) = 1 + wAdd (Kia / Ka) (Ka / A) (Kb / B) := by
  unfold wAdd
  have e : Kia * Kb / (A * B) = Kia / Ka * (Ka / A) * (Kb / B) := by
    grind
  rw [e]; grind

/-- **The lines meet at the wall.** At u_A = −1/k the rate no longer depends on B:
    every double-reciprocal line passes through the point (−1/k, 1 − 1/k). -/
theorem lines_meet_at_wall (k uB : K) (hk : k ≠ 0) :
    bisub k (-1 / k) uB = 1 - 1 / k := by
  rw [bisub_is_wAdd, wAdd_comm', wAdd_wall k uB hk]; grind

/-- Only the wall does this: if two different B give the same rate at u_A, then
    1 + k·u_A = 0. -/
theorem meeting_point_unique (k uA uB uB' : K) (hne : uB ≠ uB')
    (h : bisub k uA uB = bisub k uA uB') : 1 + k * uA = 0 := by
  have s := bisub_shift k uA uB uB'
  rw [h] at s
  have hd : uB - uB' ≠ 0 := by intro e; apply hne; grind
  have : (uB - uB') * (1 + k * uA) = 0 := by grind
  have hinv : (uB - uB')⁻¹ * (uB - uB') = 1 := Lean.Grind.Field.inv_mul_cancel hd
  calc 1 + k * uA = (uB - uB')⁻¹ * ((uB - uB') * (1 + k * uA)) := by grind
    _ = 0 := by rw [this]; grind

/-! ## 4. Saturation -/

/-- Saturation is the mult-add inverse lifted: 1/(1 + u) = 1 + inv(u). -/
theorem saturation (u : K) (_h : 1 + u ≠ 0) : 1 / (1 + u) = 1 + multAddInv u := by
  unfold multAddInv; grind

/-- Michaelis–Menten: S/(K + S) = 1 + inv(K/S). -/
theorem michaelis_menten (Km S : K) (hS : S ≠ 0) (h : Km + S ≠ 0) :
    S / (Km + S) = 1 + multAddInv (Km / S) := by
  unfold multAddInv
  have h1 : 1 + Km / S ≠ 0 := by
    intro e; apply h
    have : (1 + Km / S) * S = S + Km := by grind
    grind
  grind

/-- The free (unbound) fraction is minus the inverse: u/(1 + u) = −inv(u). -/
theorem free_fraction (u : K) (_h : 1 + u ≠ 0) : u / (1 + u) = -multAddInv u := by
  unfold multAddInv; grind

end Field

end Metabolism
