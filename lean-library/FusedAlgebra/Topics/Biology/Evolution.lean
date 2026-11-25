/-
  Fused Hyperoperation Algebra: evolution and natural selection.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Evolution.lean

  1. Haploid selection multiplies an allele's odds by 1 + s; successive selective stages
     compose by ◯ and t generations by mpow.
  2. The additive and multiplicative "no epistasis" nulls are dial 0 and dial 1.
  3. Path consistency: if each mutation's effect is linear in background fitness,
     Δ_j(x) = a_j + c_j x, then adding mutations i, j in either order agrees exactly when
     c_j a_i = c_i a_j. With a common ratio k = c/a this is the dial x ◯_k a_j, every
     mutation is neutral at the wall x = −1/k, and the dial is associative and commutative,
     so any order of any set of mutations reaches one fitness.
  4. Mutation load at many loci combines by the "or" law, wpow(−1).
  LIBRARY NOTE: first compiled in the shared library. 4 exact twins removed (star_def,
  stages, generations, load); each is marked where it stood. No proof step or statement
  changed.
-/
import FusedAlgebra.Base

namespace Evolution
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

/-! ## 1. Selection on odds -/

-- LIBRARY: stages (Two selective stages multiply odds by (1 + s₁)(1 + s₂) = 1 + (s₁ ◯ s₂)) was an exact twin of ChemQuantum.fragments_circ and is removed.

-- LIBRARY: generations (t generations of constant selection: odds multiply by (1 + s)^t = 1 + mpow s t) was an exact twin of MolecularBiology.mwc_R and is removed.

/-- Stages split anywhere: t₁ + t₂ generations are t₁ then t₂. -/
theorem generations_add (s : R) (t₁ t₂ : Nat) :
    mpow s (t₁ + t₂) = mpow s t₁ ◯ mpow s t₂ := mpow_add s t₁ t₂

/-! ## 2. Epistasis -/

/-- The multiplicative null is dial 1. -/
theorem mult_null (a b : R) : circK 1 a b = a ◯ b := by unfold circK circ; grind

/-- Epistasis against the multiplicative null is (k − 1)·ab. -/
theorem epistasis (k a b : R) : circK k a b - (a ◯ b) = (k - 1) * a * b := by
  unfold circK circ; grind

/-! ## 3. Global epistasis and path consistency -/

/-- Add mutation (a, c) to a background of excess fitness x: x + a + c·x. -/
def addMut (a c x : R) : R := x + a + c * x

/-- Two orders of adding mutations i and j differ by exactly c_j a_i − c_i a_j. -/
theorem order_gap (aᵢ cᵢ aⱼ cⱼ x : R) :
    addMut aⱼ cⱼ (addMut aᵢ cᵢ x) - addMut aᵢ cᵢ (addMut aⱼ cⱼ x) = cⱼ * aᵢ - cᵢ * aⱼ := by
  unfold addMut; grind

/-- **Path consistency.** The two orders agree iff c_j a_i = c_i a_j. -/
theorem path_consistent_iff (aᵢ cᵢ aⱼ cⱼ x : R) :
    addMut aⱼ cⱼ (addMut aᵢ cᵢ x) = addMut aᵢ cᵢ (addMut aⱼ cⱼ x) ↔ cⱼ * aᵢ = cᵢ * aⱼ := by
  constructor
  · intro h; have g := order_gap aᵢ cᵢ aⱼ cⱼ x; grind
  · intro h; have g := order_gap aᵢ cᵢ aⱼ cⱼ x; grind

/-- With slopes proportional to effects, c = k·a, adding a mutation is the dial. -/
theorem add_is_dial (k a x : R) : addMut a (k * a) x = circK k x a := by
  unfold addMut circK; grind

/-- Proportional slopes are always path-consistent. -/
theorem proportional_consistent (k aᵢ aⱼ x : R) :
    addMut aⱼ (k * aⱼ) (addMut aᵢ (k * aᵢ) x) = addMut aᵢ (k * aᵢ) (addMut aⱼ (k * aⱼ) x) := by
  unfold addMut; grind

/-- A common slope c for two mutations is consistent only if c·(aᵢ − aⱼ) = 0. -/
theorem common_slope_gap (aᵢ aⱼ c x : R) :
    addMut aⱼ c (addMut aᵢ c x) - addMut aᵢ c (addMut aⱼ c x) = c * (aᵢ - aⱼ) := by
  unfold addMut; grind

/-- The dial is associative, so three mutations reach one fitness in any grouping. -/
theorem three_mutations (k x a b : R) :
    circK k (circK k x a) b = circK k x (circK k a b) := by unfold circK; grind

/-- In the shifted coordinate 1 + k·x (fitness measured from the wall) mutations multiply. -/
theorem from_the_wall (k a x : R) : 1 + k * circK k x a = (1 + k * x) * (1 + k * a) := by
  unfold circK; grind

/-- At the wall every mutation is neutral: Δ(x) = a(1 + kx) vanishes when 1 + kx = 0. -/
theorem wall_neutral (k a x : R) (h : 1 + k * x = 0) : circK k x a = x := by
  have e : circK k x a = x + a * (1 + k * x) := by unfold circK; grind
  rw [e, h]; grind

/-! ## 4. Load and dominance -/

-- LIBRARY: load (Load across n loci of equal load u: 1 − (1 − u)ⁿ = wpow(−1, u, n)) was an exact twin of Metabolism.limitation_pow and is removed.

/-- Additive dominance (h = ½ written as two equal copies a) is dial 0. -/
theorem additive_dominance (a : R) : circK 0 a a = a + a := by unfold circK; grind

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- **Path consistency forces one dial.** If mutation i has a nonzero effect and every
    mutation j is path-consistent with it, then c_j = k·a_j with the single k = c_i/a_i. -/
theorem consistent_forces_dial (aᵢ cᵢ aⱼ cⱼ : K) (hi : aᵢ ≠ 0) (h : cⱼ * aᵢ = cᵢ * aⱼ) :
    cⱼ = (cᵢ / aᵢ) * aⱼ := by
  have inv := Lean.Grind.Field.inv_mul_cancel hi
  have e : cⱼ = cⱼ * (aᵢ⁻¹ * aᵢ) := by rw [inv]; grind
  rw [e, Lean.Grind.Field.div_eq_mul_inv]
  grind

/-- The wall: for k ≠ 0, the background x = −1/k makes every mutation neutral. -/
theorem wall_at (k a : K) (hk : k ≠ 0) : circK k (-(1 / k)) a = -(1 / k) := by
  apply wall_neutral
  have inv := Lean.Grind.Field.mul_inv_cancel hk
  rw [Lean.Grind.Field.div_eq_mul_inv]; grind

end Field

end Evolution
