module

public import FusedAlgebra.Rungs

/-!
# The dark sector: the exact statements behind the page

## Main results

* Redshifts compose by mult-add: stretching by `1 + z₁` and then by `1 + z₂` is stretching by
  `1 + (z₁ ✳ z₂)`. This is `MultAdd.one_add_multAdd`, and `n` equal stretches are
  `MultAdd.one_add_mpow`; the core-Lean statements `redshift_compose`, `redshift_compose_int`
  and `redshift_power` are those two lemmas and are not restated.
* `cube_layers`, `mpow_three`, `lcdm_is_weighted_cube`: matter dilutes as the mult-add cube of
  redshift, and flat ΛCDM is `Ωm * (1 + z) ^ 3 + (1 - Ωm) = 1 + Ωm * mpow z 3`.
* `multAdd_growth_is_exponential`: growth fractions that compose by mult-add grow exponentially
  (de Sitter expansion).
* `fexp_growth_trivial`, `fexp_growth_trivial'`, `pow_growth_trivial`: a no-go one rung up;
  growth cannot compose at the fused-exponent rung.
* `add_le_multAdd`, `sum_le_multAddList`, `corner_le_quarter`: combining non-negative pulls by
  mult-add never gives less than their sum, and the corner is at most a quarter of the square on
  the sum (mathlib's `four_mul_le_sq_add`).
* `one_add_mul_one_add_of_multAdd_eq_zero`, `mu_is_neg_minv`, `simple_mond_square`,
  `deep_mond_square`: MOND's simple interpolation `μ x = x / (1 + x)` is minus the mult-add
  inverse, stated without division.

The original stated these over `ℕ`, `ℤ` and `ℚ`; they are generalised to (ordered) commutative
semirings and rings where the same proof works. The growth statements stay on `ℕ`.
-/

@[expose] public section

namespace DarkSector

open MultAdd

/-! ### Matter dilutes as the mult-add cube of redshift -/

section Cube

variable {R : Type*} [CommSemiring R]

/-- `(1 + z) ^ 3`: unit cube, three slabs, three bars, one small cube (layers `1, 3, 3, 1`). -/
theorem cube_layers (z : R) : (1 + z) ^ 3 = 1 + 3 * z + 3 * z ^ 2 + z ^ 3 := by ring

/-- The mult-add cube of `z` is the shell of the `(1 + z)` cube around the unit cube. -/
theorem mpow_three (z : R) : mpow z 3 = 3 * z + 3 * z ^ 2 + z ^ 3 := by
  simp only [mpow_succ, mpow_zero, multAdd]; ring

/-- Flat ΛCDM, exactly: `Ωm * (1 + z) ^ 3 + (1 - Ωm) = 1 + Ωm * (z ✳ z ✳ z)`. The cosmological
constant is the unit cube; matter is the weighted shell of the `(1 + z)` cube around it. -/
theorem lcdm_is_weighted_cube {R : Type*} [CommRing R] (Om z : R) :
    Om * (1 + z) ^ 3 + (1 - Om) = 1 + Om * mpow z 3 := by
  rw [mpow_three]; ring

end Cube

/-! ### Constant expansion rate: growth fractions that compose by mult-add

If the fractional growth `A t` over a time `t` satisfies `A (t + s) = A t ✳ A s`, the scale
factor `1 + A` grows exponentially. That is de Sitter expansion, the cosmological constant. -/

theorem multAdd_growth_is_exponential (A : ℕ → ℕ) (hA : ∀ t s, A (t + s) = A t ✳ A s) (t : ℕ) :
    A t + 1 = (A 1 + 1) ^ t := by
  have h0 : A 0 = 0 := by
    have h := hA 0 0
    rw [multAdd_def] at h
    lia
  induction t with
  | zero => simp [h0]
  | succ n ih => rw [hA n 1, multAdd_add_one, ih, pow_succ]

/-! ### No-go one rung up

Expansion cannot compose at the fused-exponent rung. If growth over `t + s` were `fexp` of the
growths over `t` and `s`, in either order, and nothing has grown at `t = 0`, then nothing ever
grows. So the ladder tops out at rung 2 for expansion. -/

theorem fexp_growth_trivial (A : ℕ → ℕ) (h0 : A 0 = 0)
    (hA : ∀ t s, A (t + s) = fexp (A t) (A s)) (t : ℕ) : A t = 0 := by
  simpa [h0] using hA 0 t

theorem fexp_growth_trivial' (A : ℕ → ℕ) (h0 : A 0 = 0)
    (hA : ∀ t s, A (t + s) = fexp (A s) (A t)) (t : ℕ) : A t = 0 := by
  simpa [h0] using hA t 0

/-- The same in the unshifted form: `a (t + s) = a t ^ a s` with `a 0 = 1` forces `a = 1`. -/
theorem pow_growth_trivial (a : ℕ → ℕ) (h0 : a 0 = 1) (ha : ∀ t s, a (t + s) = a t ^ a s)
    (t : ℕ) : a t = 1 := by
  simpa [h0] using ha 0 t

/-! ### Why mult-add of gravity points the wrong way

Combining non-negative pulls by mult-add never gives less than their sum, and the extra is
exactly the corner rectangle `a * b`. -/

section Order

variable {R : Type*} [Semiring R] [PartialOrder R] [IsOrderedRing R]

theorem add_le_multAdd {a b : R} (ha : 0 ≤ a) (hb : 0 ≤ b) : a + b ≤ a ✳ b :=
  le_add_of_nonneg_right (mul_nonneg ha hb)

theorem multAddList_nonneg {l : List R} (hl : ∀ x ∈ l, 0 ≤ x) : 0 ≤ multAddList l := by
  induction l with
  | nil => simp
  | cons x xs ih =>
    have hx := hl x List.mem_cons_self
    have hxs := ih fun y hy ↦ hl y (List.mem_cons_of_mem x hy)
    rw [multAddList_cons, multAdd_def]
    positivity

theorem sum_le_multAddList {l : List R} (hl : ∀ x ∈ l, 0 ≤ x) : l.sum ≤ multAddList l := by
  induction l with
  | nil => simp
  | cons x xs ih =>
    have hx := hl x List.mem_cons_self
    have hxs : ∀ y ∈ xs, 0 ≤ y := fun y hy ↦ hl y (List.mem_cons_of_mem x hy)
    rw [List.sum_cons, multAddList_cons]
    exact (add_le_add_right (ih hxs) x).trans (add_le_multAdd hx (multAddList_nonneg hxs))

end Order

/-- The corner is at most a quarter of the square on the sum: `4 * (a * b) ≤ (a + b) ^ 2`. This
is mathlib's `four_mul_le_sq_add`. Relative to the plain sum the extra pull fades as the fields
weaken, the opposite of where dark matter shows up. -/
theorem corner_le_quarter {R : Type*} [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]
    (a b : R) : 4 * (a * b) ≤ (a + b) ^ 2 := by
  rw [← mul_assoc]; exact four_mul_le_sq_add a b

/-! ### MOND's simple interpolation is the mult-add inverse (up to sign)

The mult-add inverse of `x` is `y = -x / (1 + x)`, with `x ✳ y = 0`. MOND's "simple"
`μ x = x / (1 + x)` is exactly `-y`. Stated without division. -/

theorem one_add_mul_one_add_of_multAdd_eq_zero {R : Type*} [Semiring R] {x y : R}
    (h : x ✳ y = 0) : (1 + x) * (1 + y) = 1 := by
  rw [← one_add_multAdd, h, add_zero]

theorem mu_is_neg_minv {R : Type*} [CommRing R] {x y : R} (h : x ✳ y = 0) :
    (1 + x) * (-y) = x := by
  rw [multAdd_def] at h
  linear_combination -h

/-- The simple-MOND law `g_N = g * μ (g / a₀)` written without division: the square on the
observed pull equals two rectangles, `g_N` by `g` and `g_N` by `a₀`. -/
theorem simple_mond_square {R : Type*} [Semiring R] {g gN a₀ : R} (h : gN * (g + a₀) = g * g) :
    g * g = gN * g + gN * a₀ := by
  rw [← h, mul_add]

/-- Deep-MOND limit: the observed pull is the side of the square whose area is the rectangle
`g_N × a₀`. -/
theorem deep_mond_square {M : Type*} [Monoid M] {g gN a₀ : M} (h : g * g = gN * a₀) :
    g ^ 2 = gN * a₀ := by
  rw [sq, h]

/-! ### Small checks -/

example : (1 + 2) * (1 + 3) = 1 + ((2 : ℕ) ✳ 3) := by decide

example : mpow (1 : ℕ) 3 = 7 := by simp [mpow_three]

example : mpow (2 : ℕ) 3 = 26 := by simp [mpow_three]

end DarkSector
