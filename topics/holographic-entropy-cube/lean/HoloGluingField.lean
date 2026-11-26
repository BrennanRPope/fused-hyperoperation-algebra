/-
  HoloGluingField: the same laws with division, on any field   (DRAFT, ours)
  First compiled in the shared library on Nov 26. `gluing` needed help: the
  divisions are split into products of inverses before `grind`.
-/
import FusedAlgebra.Base

namespace HoloGluingField

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

variable {K : Type} [Field K]

/-- Cross-ratio of A = (a₁,a₂) and B = (b₁,b₂) on a line. -/
def xr (a₁ a₂ b₁ b₂ : K) : K := (a₂ - a₁) * (b₂ - b₁) / ((b₁ - a₁) * (b₂ - a₂))

/-- The Möbius map that sends a₁ to 0 and a₂ to ∞ (A's geodesic becomes the vertical axis). -/
def mob (a₁ a₂ p : K) : K := (p - a₁) / (p - a₂)

/-- Shadow lemma: 1 − x(A,B) = f(b₂)/f(b₁). -/
theorem shadow (a₁ a₂ b₁ b₂ : K)
    (h1 : b₁ - a₁ ≠ 0) (h2 : b₂ - a₂ ≠ 0) (h3 : b₁ - a₂ ≠ 0) (h4 : b₂ - a₁ ≠ 0) :
    1 - xr a₁ a₂ b₁ b₂ = mob a₁ a₂ b₂ / mob a₁ a₂ b₁ := by
  unfold xr mob; grind

/-- The gluing law: x(A, B₁ ∪ B₂) = x(A, B₁) ∘₋₁ x(A, B₂). -/
theorem gluing (a₁ a₂ b₁ m b₂ : K)
    (h1 : b₁ - a₁ ≠ 0) (h2 : m - a₂ ≠ 0) (h3 : m - a₁ ≠ 0) (h4 : b₂ - a₂ ≠ 0) :
    xr a₁ a₂ b₁ b₂ = circK (-1) (xr a₁ a₂ b₁ m) (xr a₁ a₂ m b₂) := by
  have e1 : xr a₁ a₂ b₁ b₂ = 1 - (b₁ - a₂) * (b₂ - a₁) / ((b₁ - a₁) * (b₂ - a₂)) := by
    unfold xr; grind
  have e2 : xr a₁ a₂ b₁ m = 1 - (b₁ - a₂) * (m - a₁) / ((b₁ - a₁) * (m - a₂)) := by
    unfold xr; grind
  have e3 : xr a₁ a₂ m b₂ = 1 - (m - a₂) * (b₂ - a₁) / ((m - a₁) * (b₂ - a₂)) := by
    unfold xr; grind
  have e4 : (b₁ - a₂) * (b₂ - a₁) / ((b₁ - a₁) * (b₂ - a₂))
      = ((b₁ - a₂) * (m - a₁) / ((b₁ - a₁) * (m - a₂))) * ((m - a₂) * (b₂ - a₁) / ((m - a₁) * (b₂ - a₂))) := by
    have i2 : (m - a₂) * (m - a₂)⁻¹ = 1 := Field.mul_inv_cancel h2
    have i3 : (m - a₁) * (m - a₁)⁻¹ = 1 := Field.mul_inv_cancel h3
    simp only [Field.div_eq_mul_inv, Field.inv_mul]
    grind
  rw [e1, e2, e3, e4]; unfold circK; grind

/-- The flip φ₋₁(x) = x/(1 − x) carries ∘₋₁ to ∘₁. -/
theorem flip_glue (x y : K) (hx : 1 - x ≠ 0) (hy : 1 - y ≠ 0) :
    circK (-1) x y / (1 - circK (-1) x y) = circK 1 (x / (1 - x)) (y / (1 - y)) := by
  unfold circK; grind

end HoloGluingField
