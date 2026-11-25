/-
  Fused Hyperoperation Algebra: the shape through the rungs (whole-rung facts).
  Core Lean 4 only (no mathlib). Library style: drop into FusedAlgebra/Topics/.

  STATUS: NOT COMPILED. No Lean toolchain in the workspace where this was written
  (Nov 25, 2026). Every statement is also checked as an exact algebraic identity in
  ../symbolic_checks.py (sympy). Proofs follow the pattern of the compiled library
  (`unfold` then `grind`, whose ring solver handles CommRing and Field goals).

  What the real-height sheet needs (the regular iterate of e^x − 1) lives in the reals
  and is checked numerically in ../rung_shape.py; here are the exact pieces at whole rungs.

  LIBRARY NOTE: copied from rung-shape/lean/RungShape.lean (Nov 25, the 17-theorem version with
  the sheet-flow lens section) and compiled here with core Lean 4.37.0-rc1. It was already
  written on FusedAlgebra.Base.
-/
import FusedAlgebra.Base

namespace RungShape

open MultAdd
open Lean.Grind

section Ring
variable {R : Type} [CommRing R]

/-! ### The saddle: the twist of rung 1 is exactly the corner rectangle -/

/-- Mixed second difference of mult-add = h·k, the area of the corner rectangle.
  (Rung 0 has twist 0: a flat plane. Rung 1 has twist 1: the saddle.) -/
theorem twist_rung1 (a b h k : R) :
    multAdd (a + h) (b + k) - multAdd (a + h) b - multAdd a (b + k) + multAdd a b = h * k := by
  unfold multAdd; grind

/-- The corner dial a + b + t·ab has twist t (the corner weight). -/
theorem twist_dial (t a b h k : R) :
    wAdd t (a + h) (b + k) - wAdd t (a + h) b - wAdd t a (b + k) + wAdd t a b = t * (h * k) := by
  unfold wAdd; grind

/-- Rung 1's graph is ruled: along any fixed a it is a straight line in b. -/
theorem ruled_rung1 (a b : R) : multAdd a b = a + (1 + a) * b := by
  unfold multAdd; grind

/-- ... and along any fixed b it is a straight line in a (doubly ruled saddle). -/
theorem ruled_rung1' (a b : R) : multAdd a b = b + (1 + b) * a := by
  unfold multAdd; grind

/-- Flat along the axes: rung 1 adds nothing on the a-axis (the identity line). -/
theorem axis_rung1 (a : R) : multAdd a 0 - a = 0 := by
  unfold multAdd; grind

/-! ### The near-identity flow: one step of rung height is a Möbius step

  Near 0, rung height s acts as x ↦ 2x/(2 − s·x) (the flow of x²/2, the first term of
  e^x − 1 − x). Heights add, stated without division: -/

theorem flow_heights_add (s t x y z : R)
    (h1 : y * (2 - t * x) = 2 * x) (h2 : z * (2 - s * y) = 2 * y) :
    2 * (z * (2 - (s + t) * x)) = 2 * (2 * x) := by
  grind

/-- Whole rung height 0 is the identity of the near-identity flow. -/
theorem flow_zero (x : R) : x * (2 - 0 * x) = 2 * x := by grind

end Ring

section Field
variable {K : Type} [Field K]

/-! ### The sheet's slices at whole rungs: the level curves solved exactly -/

/-- Rung 0 slice: a + b = c is the straight line b = c − a. -/
theorem sheet_rung0 (a c : K) : a + (c - a) = c := by grind

/-- Rung 1 slice: a ⊞ b = c is the hyperbola b = (c − a)/(1 + a). -/
theorem sheet_rung1 (a c : K) (h : 1 + a ≠ 0) : multAdd a ((c - a) / (1 + a)) = c := by
  unfold multAdd; grind

/-- Corner-dial slice: a ⊞ₜ b = c is b = (c − a)/(1 + t·a). -/
theorem sheet_dial (t a c : K) (h : 1 + t * a ≠ 0) : wAdd t a ((c - a) / (1 + t * a)) = c := by
  unfold wAdd; grind

/-- Every slice passes through (c, 0) and (0, c): the sheet's two fixed edges. -/
theorem sheet_edges (c : K) : multAdd c 0 = c ∧ multAdd 0 c = c := by
  unfold multAdd; constructor <;> grind

/-- The near-identity flow in closed form: φ_t(x) = 2x/(2 − t x) satisfies its defining relation. -/
theorem flow_closed (t x : K) (h : 2 - t * x ≠ 0) : (2 * x / (2 - t * x)) * (2 - t * x) = 2 * x := by
  grind

end Field

/-! ### The two dials at whole numbers (on Nat) -/

/-- d copies of a, added (count dial at rung 0). -/
def addPow (a : Nat) : Nat → Nat
  | 0 => 0
  | n + 1 => addPow a n + a

/-- Rung 0, count d: the line d·a. -/
theorem addPow_eq (a d : Nat) : addPow a d = d * a := by
  induction d with
  | zero => simp [addPow]
  | succ n ih => simp only [addPow, ih, Nat.succ_mul]

/-- Rung 1, count d, edge 1: the d-cube's corners, 2^d (so the mult-add power is 2^d − 1). -/
theorem cube_count (d : Nat) : 1 + mpow (1 : Nat) d = 2 ^ d := by
  simpa using mpow_lift (1 : Nat) d

/-! ### Sheet versus flow (Nov 25 follow-up)

  Every rung is addition seen through a lens: rung s is φ(ψ a + ψ b) with ψ = E⁻ˢ and φ = Eˢ,
  the flow at height s and its inverse. The rung sheet's slice is then exactly the flow's
  pull-back of a straight line, at every height, whole or not. Stated for any lens: -/

section Lens
variable {α β : Type} [Add β]

/-- The rung built from a lens ψ with inverse φ. -/
def lensOp (φ : β → α) (ψ : α → β) (a b : α) : α := φ (ψ a + ψ b)

/-- A point is on the sheet (rung value c) exactly when its lens image is on the straight
  line u + v = ψ c. -/
theorem sheet_iff_line (φ : β → α) (ψ : α → β)
    (hψφ : ∀ u, ψ (φ u) = u) (hφψ : ∀ x, φ (ψ x) = x) (a b c : α) :
    lensOp φ ψ a b = c ↔ ψ a + ψ b = ψ c := by
  unfold lensOp
  constructor
  · intro h; rw [← h, hψφ]
  · intro h; rw [h, hφψ]

/-- The diagonal point of the sheet: a ⊕ a = c exactly when ψ a + ψ a = ψ c
  (so a is the flow carrying half of ψ c back out). -/
theorem sheet_diagonal (φ : β → α) (ψ : α → β)
    (hψφ : ∀ u, ψ (φ u) = u) (hφψ : ∀ x, φ (ψ x) = x) (a c : α) :
    lensOp φ ψ a a = c ↔ ψ a + ψ a = ψ c :=
  sheet_iff_line φ ψ hψφ hφψ a a c

end Lens

section Ring
variable {R : Type} [CommRing R]

/-- Rung 1 with no logarithms: its lens is x ↦ 1 + x into multiplication, and the slice
  a ⊞ b = c is the "straight line" (1 + a)(1 + b) = 1 + c of that lens. -/
theorem sheet_rung1_lens (a b c : R) : multAdd a b = c ↔ (1 + a) * (1 + b) = 1 + c := by
  constructor
  · intro h; rw [← multAdd_lift, h]
  · intro h; have e := multAdd_lift a b; grind

end Ring

end RungShape
