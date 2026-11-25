module

public import FusedAlgebra.Basic

/-!
# The shape through the rungs (whole-rung facts)

The graph of mult-add `(a, b) ↦ a ✳ b` is a doubly ruled saddle whose twist is the corner
rectangle. This file collects the exact algebraic pieces of that picture at whole rungs (the
real-height sheet, built on the regular iterate of `exp x - 1`, is not treated here).

## Main results

* `twist_rung1`, `twist_dial`: the mixed second difference of mult-add is `h * k`, and of the
  weighted mult-add `wAdd t` it is `t * (h * k)`.
* `ruled_rung1`, `ruled_rung1'`: the graph of mult-add is doubly ruled.
* `flow_heights_add`: the near-identity flow `x ↦ 2 * x / (2 - s * x)` (a Möbius map): heights
  add.
* `sheet_rung1`, `sheet_dial`: the level curves `a ✳ b = c` and `wAdd t a b = c`, solved for
  `b`.
* `lensOp`, `sheet_iff_line`: every rung is addition seen through a lens, and its level curve is
  the pull-back of a straight line; `sheet_rung1_lens` is rung 1, with lens `x ↦ 1 + x` into
  multiplication.

These are elementary identities; none is a named mathlib theorem, and several are one-line
consequences of the core API in `FusedAlgebra.Basic`.
-/

@[expose] public section

namespace RungShape

open MultAdd

section CommRing

variable {R : Type*} [CommRing R]

/-! ### The saddle: the twist of rung 1 is exactly the corner rectangle -/

/-- The mixed second difference of mult-add is `h * k`, the area of the corner rectangle.
(Rung 0 has twist `0`: a flat plane. Rung 1 has twist `1`: the saddle.) -/
theorem twist_rung1 (a b h k : R) :
    ((a + h) ✳ (b + k)) - ((a + h) ✳ b) - (a ✳ (b + k)) + (a ✳ b) = h * k := by
  unfold multAdd; ring

/-- The corner dial `a + b + t * a * b` has twist `t` (the corner weight). -/
theorem twist_dial (t a b h k : R) :
    wAdd t (a + h) (b + k) - wAdd t (a + h) b - wAdd t a (b + k) + wAdd t a b = t * (h * k) := by
  unfold wAdd; ring

/-- The graph of rung 1 is ruled: along any fixed `a` it is a straight line in `b`. -/
theorem ruled_rung1 (a b : R) : a ✳ b = a + (1 + a) * b := by
  unfold multAdd; ring

/-- ... and along any fixed `b` it is a straight line in `a` (a doubly ruled saddle). -/
theorem ruled_rung1' (a b : R) : a ✳ b = b + (1 + b) * a := by
  unfold multAdd; ring

/-- Flat along the axes: rung 1 adds nothing on the `a`-axis (the identity line). -/
theorem axis_rung1 (a : R) : (a ✳ 0) - a = 0 := by simp

/-! ### The near-identity flow: one step of rung height is a Möbius step

Near `0`, rung height `s` acts as `x ↦ 2 * x / (2 - s * x)` (the flow of `x ^ 2 / 2`, the first
term of `exp x - 1 - x`). Heights add, stated without division. -/

theorem flow_heights_add (s t x y z : R)
    (h₁ : y * (2 - t * x) = 2 * x) (h₂ : z * (2 - s * y) = 2 * y) :
    2 * (z * (2 - (s + t) * x)) = 2 * (2 * x) := by
  linear_combination (2 + s * z) * h₁ + (2 - t * x) * h₂

/-- Whole rung height `0` is the identity of the near-identity flow. -/
theorem flow_zero (x : R) : x * (2 - 0 * x) = 2 * x := by ring

/-- Rung 1 with no logarithms: its lens is `x ↦ 1 + x` into multiplication, and the slice
`a ✳ b = c` is the "straight line" `(1 + a) * (1 + b) = 1 + c` of that lens. -/
theorem sheet_rung1_lens (a b c : R) : a ✳ b = c ↔ (1 + a) * (1 + b) = 1 + c := by
  rw [← one_add_multAdd, add_right_inj]

end CommRing

section Field

variable {K : Type*} [Field K]

/-! ### The slices of the sheet at whole rungs: the level curves solved exactly

The rung 0 slice `a + b = c` is the straight line `b = c - a`: `add_sub_cancel`. -/

/-- Rung 1 slice: `a ✳ b = c` is the hyperbola `b = (c - a) / (1 + a)`. -/
theorem sheet_rung1 (a c : K) (h : 1 + a ≠ 0) : a ✳ ((c - a) / (1 + a)) = c := by
  rw [multAdd_comm]; exact msub_multAdd c a h

/-- Corner-dial slice: `wAdd t a b = c` is `b = (c - a) / (1 + t * a)`. -/
theorem sheet_dial (t a c : K) (h : 1 + t * a ≠ 0) : wAdd t a ((c - a) / (1 + t * a)) = c := by
  have e : wAdd t a ((c - a) / (1 + t * a)) = a + (c - a) / (1 + t * a) * (1 + t * a) := by
    unfold wAdd; ring
  rw [e, div_mul_cancel₀ _ h, add_sub_cancel]

/-- Every slice passes through `(c, 0)` and `(0, c)`: the two fixed edges of the sheet. -/
theorem sheet_edges (c : K) : c ✳ 0 = c ∧ 0 ✳ c = c := ⟨multAdd_zero c, zero_multAdd c⟩

/-- The near-identity flow in closed form: `φ t x = 2 * x / (2 - t * x)` satisfies its defining
relation. -/
theorem flow_closed (t x : K) (h : 2 - t * x ≠ 0) :
    2 * x / (2 - t * x) * (2 - t * x) = 2 * x :=
  div_mul_cancel₀ _ h

end Field

/-! ### The two dials at whole numbers

At rung 0, `d` copies of `a` added give the line `d * a` (`smul_eq_mul`). -/

/-- Rung 1, count `d`, edge `1`: the corners of the `d`-cube, `2 ^ d` (so the mult-add power is
`2 ^ d - 1`). -/
theorem cube_count (d : ℕ) : 1 + mpow (1 : ℕ) d = 2 ^ d := by
  rw [one_add_mpow]

/-! ### Sheet versus flow

Every rung is addition seen through a lens: rung `s` is `φ (ψ a + ψ b)` with `ψ` the flow at
height `-s` and `φ` its inverse. The slice of the rung sheet is then exactly the pull-back by the
flow of a straight line, at every height, whole or not. This is stated for any lens. -/

section Lens

variable {α β : Type*} [Add β]

/-- The rung built from a lens `ψ` with inverse `φ`. -/
def lensOp (φ : β → α) (ψ : α → β) (a b : α) : α := φ (ψ a + ψ b)

/-- A point is on the sheet (rung value `c`) exactly when its lens image is on the straight line
`u + v = ψ c`. -/
theorem sheet_iff_line {φ : β → α} {ψ : α → β} (hψφ : Function.LeftInverse ψ φ)
    (hφψ : Function.LeftInverse φ ψ) (a b c : α) : lensOp φ ψ a b = c ↔ ψ a + ψ b = ψ c := by
  unfold lensOp
  constructor
  · rintro rfl; exact (hψφ _).symm
  · intro h; rw [h, hφψ]

/-- The diagonal point of the sheet: `lensOp φ ψ a a = c` exactly when `ψ a + ψ a = ψ c` (so `a`
is the flow carrying half of `ψ c` back out). -/
theorem sheet_diagonal {φ : β → α} {ψ : α → β} (hψφ : Function.LeftInverse ψ φ)
    (hφψ : Function.LeftInverse φ ψ) (a c : α) : lensOp φ ψ a a = c ↔ ψ a + ψ a = ψ c :=
  sheet_iff_line hψφ hφψ a a c

end Lens

end RungShape
