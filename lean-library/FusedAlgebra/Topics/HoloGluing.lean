/-
  HoloGluing: the AdS₃ gluing law, proved by shadows   (DRAFT, ours)
  =====================================================================

  Written in the style of lean-library (core Lean 4.37, no mathlib, no `sorry`,
  `grind` for the ring arithmetic), now in the library as
  FusedAlgebra/Topics/HoloGluing.lean.

  First compiled in the shared library on Nov 26, with names updated to 1.0
  (CircleOp.circK, called MultAdd.wAdd before 1.0).

  Everything here is stated with denominators cleared, so every proof is a
  pure ring identity (the most robust thing `grind` does). The field forms
  with division are in HoloGluingField.lean.

  What is proved (for points on a line, in any commutative ring):
    ptolemy          (a₂−a₁)(b₂−b₁) + (b₁−a₂)(b₂−a₁) = (b₁−a₁)(b₂−a₂)
    oneSub_numer     1 − x(A,B) has numerator (b₁−a₂)(b₂−a₁): the "shadow" ratio
    gluing_cleared   1 − x multiplies when B is split at m  (so x glues by ∘₋₁)
    complement_mul   1 − (x ∘₋₁ y) = (1 − x)(1 − y)
    odds_glue        the odds X = x/(1−x) glue by ∘₁
    height_flow      (U₁U₂ − 1) = (U₁ − 1) ∘₁ (U₂ − 1)   (U = e^t, t = shadow length)
    conformal        conformal factors ρ cancel from the cross-ratio
    monogamy_glued   max(1,X)·max(1,Y) ≤ max(1, X ∘₁ Y) for X, Y ≥ 0 (on Rat)

  What is NOT in Lean (stated in the write-up as assumptions or as geometry):
    the Ryu–Takayanagi formula, and the hyperbolic-geometry fact that the
    shadow of B on A's geodesic has length t with e^{−t} = 1 − x.
-/
import FusedAlgebra.Base

namespace HoloGluing

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [CommRing R]

/-- Numerator and denominator of the cross-ratio of A = (a₁, a₂) and B = (x, y). -/
def numer (a₁ a₂ x y : R) : R := (a₂ - a₁) * (y - x)
def denom (a₁ a₂ x y : R) : R := (x - a₁) * (y - a₂)

/-- Ptolemy's relation for four points on a line. -/
theorem ptolemy (a₁ a₂ b₁ b₂ : R) :
    (a₂ - a₁) * (b₂ - b₁) + (b₁ - a₂) * (b₂ - a₁) = (b₁ - a₁) * (b₂ - a₂) := by
  grind

/-- 1 − x(A,B) = (b₁ − a₂)(b₂ − a₁) / denom: the numerator is the shadow ratio
    f(b₂)/f(b₁) of the Möbius map f(p) = (p − a₁)/(p − a₂) with denominators cleared. -/
theorem oneSub_numer (a₁ a₂ x y : R) :
    denom a₁ a₂ x y - numer a₁ a₂ x y = (x - a₂) * (y - a₁) := by
  unfold denom numer; grind

/-- The gluing law with denominators cleared: (1 − x₁₂) = (1 − x₁)(1 − x₂)
    when B = (b₁, b₂) is cut at m into B₁ = (b₁, m) and B₂ = (m, b₂). -/
theorem gluing_cleared (a₁ a₂ b₁ m b₂ : R) :
    (denom a₁ a₂ b₁ b₂ - numer a₁ a₂ b₁ b₂) * denom a₁ a₂ b₁ m * denom a₁ a₂ m b₂
      = (denom a₁ a₂ b₁ m - numer a₁ a₂ b₁ m) * (denom a₁ a₂ m b₂ - numer a₁ a₂ m b₂)
          * denom a₁ a₂ b₁ b₂ := by
  unfold denom numer; grind

/-- ∘₋₁ is "complements multiply" (like the union of independent events). -/
theorem complement_mul (x y : R) : 1 - circK (-1) x y = (1 - x) * (1 - y) := by
  unfold circK; grind

/-- So a ∘₋₁ law follows from multiplying complements. -/
theorem of_complement_mul (x y z : R) (h : 1 - z = (1 - x) * (1 - y)) : z = circK (-1) x y := by
  unfold circK; grind

/-- The flip: if X·(1−x) = x and Y·(1−y) = y (odds), then the odds of x ∘₋₁ y is X ∘₁ Y.
    Stated with u = 1 − x, v = 1 − y and denominators cleared. -/
theorem odds_glue (X Y u v : R) (hX : X * u = 1 - u) (hY : Y * v = 1 - v) :
    circK 1 X Y * (u * v) = 1 - u * v := by
  have e : circK 1 X Y * (u * v) = (X * u) * v + (Y * v) * u + (X * u) * (Y * v) := by
    unfold circK; grind
  rw [e, hX, hY]; grind

/-- Height flow: with U = e^t (t the shadow length), odds = U − 1 and adding
    shadows multiplies U, so the odds compose by ∘₁. -/
theorem height_flow (U₁ U₂ : R) : U₁ * U₂ - 1 = circK 1 (U₁ - 1) (U₂ - 1) := by
  unfold circK; grind

/-- Conformal factors cancel: with d(u,v) = (w v − w u)·ρ u·ρ v (circle: w = tan(θ/2);
    BTZ: w = e^{2πu/β}), the cross-ratio built from d equals the one built from w. -/
theorem conformal (w₁ w₂ w₃ w₄ r₁ r₂ r₃ r₄ : R) :
    ((w₂ - w₁) * r₁ * r₂) * ((w₄ - w₃) * r₃ * r₄) * ((w₃ - w₁) * (w₄ - w₂))
      = ((w₂ - w₁) * (w₄ - w₃)) * (((w₃ - w₁) * r₁ * r₃) * ((w₄ - w₂) * r₂ * r₄)) := by
  grind

end Ring

/-! ## Monogamy for A next to B₁, B₂ (multiplicative form, on Rat)

  With RT, e^{3I/c} = max(1, X) where X is the odds (the two phases of the
  two-interval minimal surface). For B = B₁ ∪ B₂ glued, X₁₂ = X₁ ∘₁ X₂, and
  I₃ ≤ 0 becomes max(1,X₁)·max(1,X₂) ≤ max(1, X₁ ∘₁ X₂). -/

/-- e^{3I/c} as a function of the odds: 1 (disconnected) or X (connected). -/
def phase (X : Rat) : Rat := if X ≤ 1 then 1 else X

theorem one_le_phase (X : Rat) : 1 ≤ phase X := by
  unfold phase; split <;> grind

theorem le_circK_one (X Y : Rat) (hX : 0 ≤ X) (hY : 0 ≤ Y) : X * Y ≤ circK 1 X Y := by
  unfold circK; grind

theorem monogamy_glued (X Y : Rat) (hX : 0 ≤ X) (hY : 0 ≤ Y) :
    phase X * phase Y ≤ phase (circK 1 X Y) := by
  have hXY := Rat.mul_nonneg hX hY
  unfold phase circK
  split <;> split <;> split <;> grind

end HoloGluing
