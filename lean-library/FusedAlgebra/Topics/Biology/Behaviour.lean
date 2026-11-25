/-
  Fused Hyperoperation Algebra: behaviour (learning, discounting, vigilance).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports FusedAlgebra.Base only.
  Intended library path: FusedAlgebra/Topics/Biology/Behaviour.lean

  1. Rescorla–Wagner: V_n = λ·wpow(−1, αβ, n); training phases compose by the "or" law.
  2. Discounting, in the 'interest' u = 1/d − 1: hyperbolic u adds over split delays (dial 0),
     exponential u composes by ◯ (dial 1), and any family with 1 + κu multiplicative over
     delays composes by the dial κ. Splitting a hyperbolic delay raises 1/d by exactly k²D₁D₂.
  3. Many-eyes vigilance is wpow(−1, p, n).
  LIBRARY NOTE: first compiled in the shared library. 4 exact twins removed (star_def,
  rw_curve, exponential_split, many_eyes); each is marked where it stood. No proof step or
  statement changed.
-/
import FusedAlgebra.Base

namespace Behaviour
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

-- LIBRARY: star_def was an exact twin of MolecularBiology.star_def and is removed.

/-- One Rescorla–Wagner step in units of λ: V ↦ V + a(1 − V) is V ◯₋₁ a. -/
theorem rw_step (a V : R) : V + a * (1 - V) = circK (-1) V a := by unfold circK; grind

-- LIBRARY: rw_curve (The learning curve after n trials) was an exact twin of CircleOp.wpow_neg_one and is removed.

/-- Two training phases compose by the "or" law. -/
theorem rw_phases (a : R) (m n : Nat) :
    wpow (-1) a (m + n) = circK (-1) (wpow (-1) a m) (wpow (-1) a n) := by
  rw [wpow_neg_one, wpow_neg_one, wpow_neg_one]
  unfold circK
  have e : (1 - a) ^ (m + n) = (1 - a) ^ m * (1 - a) ^ n := pow_add' (1 - a) m n
  rw [e]; grind

/-- Hyperbolic interest k·D adds over split delays: dial 0. -/
theorem hyperbolic_split (k D₁ D₂ : R) : circK 0 (k * D₁) (k * D₂) = k * (D₁ + D₂) := by
  unfold circK; grind

/-- **Subadditivity.** Discounting a hyperbolic delay in two pieces raises 1/d by k²D₁D₂. -/
theorem subadditive (k D₁ D₂ : R) :
    (1 + k * D₁) * (1 + k * D₂) - (1 + k * (D₁ + D₂)) = k * k * D₁ * D₂ := by grind

-- LIBRARY: exponential_split (Exponential interest composes by plain ◯: (1 + u₁)(1 + u₂) = 1 + (u₁ ◯ u₂)) was an exact twin of ChemQuantum.fragments_circ and is removed.

/-- **Dial-stationary discounting.** If 1 + κu multiplies over split delays, u composes by ◯_κ. -/
theorem dial_split (κ u₁ u₂ : R) : 1 + κ * circK κ u₁ u₂ = (1 + κ * u₁) * (1 + κ * u₂) := by
  unfold circK; grind

/-- An integer-step version: n equal pieces of interest u compose to wpow(κ, u, n). -/
theorem dial_steps (κ u : R) (n : Nat) : 1 + κ * wpow κ u n = (1 + κ * u) ^ n := wpow_lift κ u n

-- LIBRARY: many_eyes (Many eyes: a group of n vigilant animals) was an exact twin of Metabolism.limitation_pow and is removed.

end Ring

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- The hyperbolic discount factor is minus the circle inverse (called mult-add inverse before 1.0) of the interest. -/
theorem hyperbolic_factor (k D : K) (h : 1 + k * D ≠ 0) : 1 / (1 + k * D) = 1 + circInv (k * D) := by
  unfold circInv
  have i := Lean.Grind.Field.inv_mul_cancel h
  simp only [Lean.Grind.Field.div_eq_mul_inv]; grind

/-- Herrnstein's hyperbola r/(r + rₑ) is −inv(r/rₑ). -/
theorem herrnstein (r re : K) (hre : re ≠ 0) (h : r + re ≠ 0) :
    r / (r + re) = -circInv (r / re) := by
  unfold circInv
  have i := Lean.Grind.Field.inv_mul_cancel hre
  have j := Lean.Grind.Field.inv_mul_cancel h
  have e : 1 + r * re⁻¹ = (r + re) * re⁻¹ := by grind
  simp only [Lean.Grind.Field.div_eq_mul_inv]
  rw [e, Lean.Grind.Field.inv_mul, Lean.Grind.Field.inv_inv]
  grind

end Field

end Behaviour
