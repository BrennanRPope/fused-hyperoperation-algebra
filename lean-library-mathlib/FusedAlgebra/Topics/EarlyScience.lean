module

public import FusedAlgebra.Basic

/-!
# Early science on the landscape

Two kinds of law appear in early science.

* Combination laws (two things of one kind make one thing of that kind) are points `(k, q)` of
  the landscape: lenses a distance `d` apart are weighted mult-add at `k = -d` (`lenses_apart`),
  Newton's cooling is the weighted power at `k = -1` (`cooling`), parallel resistors are addition
  through the lens `x ↦ x⁻¹` (`parallel_is_inverted_addition`).
* Product laws `y = c * x ^ m * z ^ n` (Ohm, Boyle, Kepler III, gravity) are not points. They are
  mult-add itself, read on percent changes: relative errors and percent changes of a product
  combine by mult-add (`product_error`, `constant_product`), and the exponents are mult-add
  powers (`kepler_third`, `inverse_square`).

The resistor ladder (`ladder`) produces the Fibonacci numbers (`ladder_fib`) and walks into the
golden wall (`ladder_walls`, `ladder_hyperbola`).

All of this is known science and known algebra; the file records that each law is a statement
about mult-add. Ray-transfer matrices are mathlib's `Matrix (Fin 2) (Fin 2) R`, and the Fibonacci
numbers are `Nat.fib`. Statements are generalised from `ℤ` to commutative rings where that is
free.
-/

@[expose] public section

namespace SimpleLaws

open MultAdd Matrix

variable {R : Type*} [CommRing R] {K : Type*} [Field K]

/-! ### Counting -/

/-- Non-empty choices: a set of `m` things has `2 ^ m - 1` non-empty subsets; joining two separate
sets combines those counts by mult-add (pick from either, or both, but not nothing). -/
theorem nonempty_subsets (m n : ℕ) : ((2 : R) ^ m - 1) ✳ (2 ^ n - 1) = 2 ^ (m + n) - 1 := by
  simp only [multAdd_def]; ring

/-! ### Measurement -/

/-- Relative errors in a product combine by mult-add, exactly. The textbook rule "they add" is the
`k = 0` shadow that drops the corner `e₁ * e₂`. -/
theorem product_error (x y e₁ e₂ : R) :
    x * (1 + e₁) * (y * (1 + e₂)) = x * y * (1 + (e₁ ✳ e₂)) := by
  rw [one_add_multAdd]; ring

/-- The corner that the first-order rule drops. -/
theorem product_error_corner (e₁ e₂ : R) : e₁ ✳ e₂ - (e₁ + e₂) = e₁ * e₂ := by
  simp only [multAdd_def]; ring

/-! ### Archimedes' lever and Boyle's law: products held constant -/

/-- If one factor of a constant product grows by `a`, the other must change by `b` with
`a ✳ b = 0`: the lever (weight × arm) and Boyle (pressure × volume) are the same picture. -/
theorem constant_product [NoZeroDivisors R] (P V a b : R) (hPV : P * V ≠ 0) :
    P * (1 + a) * (V * (1 + b)) = P * V ↔ a ✳ b = 0 := by
  rw [product_error, mul_add, mul_one, add_eq_left, mul_eq_zero, or_iff_right hPV]

/-- In a field the partner change is the mult-add inverse `-a / (1 + a)`. -/
theorem boyle_partner (a b : K) (ha : 1 + a ≠ 0) (h : a ✳ b = 0) : b = multAddInv a := by
  rw [multAddInv, eq_div_iff ha]
  rw [multAdd_def] at h
  linear_combination h

/-- Pressure gains combine at `k = 1` while the matching volume losses combine at `k = -1`: if
`(1 + a) * (1 - d₁) = 1` and `(1 + b) * (1 - d₂) = 1` then the totals match the same way. -/
theorem boyle_gain_loss (a b d₁ d₂ : R) (h₁ : (1 + a) * (1 - d₁) = 1)
    (h₂ : (1 + b) * (1 - d₂) = 1) : (1 + (a ✳ b)) * (1 - wAdd (-1) d₁ d₂) = 1 := by
  rw [one_add_multAdd, wAdd_neg_one]
  linear_combination (1 + b) * (1 - d₂) * h₁ + h₂

/-! ### Kepler's third law and Newton's inverse square: mult-add powers -/

/-- `T ^ 2 ∝ a ^ 3` on percent changes: `(1 + t) ^ 2 = (1 + α) ^ 3` is `t ✳ t = α ✳ α ✳ α` (the
square's gnomon equals the cube's layers). -/
theorem kepler_third (t α : R) : (1 + t) ^ 2 = (1 + α) ^ 3 ↔ t ✳ t = α ✳ α ✳ α := by
  have e : (1 + t) ^ 2 - (1 + α) ^ 3 = t ✳ t - (α ✳ α ✳ α) := by
    simp only [multAdd_def]; ring
  rw [← sub_eq_zero, e, sub_eq_zero]

/-- Four times as far out, eight times as long a year. -/
example : (1 + 7) ^ 2 = (1 + 3) ^ 3 := by decide

/-- Gravity: when the distance grows by `ρ` the force changes by the mult-add inverse of
`ρ ✳ ρ`. -/
theorem inverse_square (ρ : K) (h : 1 + (ρ ✳ ρ) ≠ 0) :
    ((1 + ρ) ^ 2)⁻¹ = 1 + multAddInv (ρ ✳ ρ) := by
  rw [sq, ← one_add_multAdd]; exact inv_one_add _ h

/-! ### Newton's law of cooling: the gap left behind compounds at `k = -1` -/

/-- Each step closes a fraction `a` of the temperature gap; after `n` steps the part closed is
`1 - (1 - a) ^ n`. This is `MultAdd.wpow_neg_one`. -/
theorem cooling (a : R) (n : ℕ) : wpow (-1) a n = 1 - (1 - a) ^ n := wpow_neg_one a n

/-! ### Optics: lenses in contact add; lenses apart are weighted mult-add with `k = -d`

A ray-transfer matrix is a `2 × 2` matrix `!![A, B; C, D]`; the power of the system is `-C`, the
entry in row `1`, column `0`. -/

/-- The ray-transfer matrix of a thin lens of power `P`. -/
def lens (P : R) : Matrix (Fin 2) (Fin 2) R := !![1, 0; -P, 1]

/-- The ray-transfer matrix of a gap of length `d`. -/
def gap (d : R) : Matrix (Fin 2) (Fin 2) R := !![1, d; 0, 1]

/-- Two lenses in contact: powers add (the point `(0, 0)`). -/
theorem lenses_in_contact (P₁ P₂ : R) : (lens P₂ * lens P₁) 1 0 = -(P₁ + P₂) := by
  simp [lens, add_comm]

/-- Two lenses a distance `d` apart: total power `P₁ + P₂ - d * P₁ * P₂`, weighted mult-add at
`k = -d` (Gullstrand's equation). -/
theorem lenses_apart (P₁ d P₂ : R) :
    (lens P₂ * (gap d * lens P₁)) 1 0 = -wAdd (-d) P₁ P₂ := by
  simp [lens, gap, wAdd_def]; ring

/-- Each lens's squeeze factor `1 - d * P` (how much a parallel beam narrows over the gap)
multiplies. This is `MultAdd.one_add_mul_wAdd` at `k = -d`. -/
theorem squeeze_multiplies (P₁ d P₂ : R) :
    1 - d * wAdd (-d) P₁ P₂ = (1 - d * P₁) * (1 - d * P₂) := by
  simp only [wAdd_def]; ring

/-- The wall: a first lens that focuses exactly on the second (`d * P₁ = 1`) makes the second lens
change nothing: the field lens of telescope eyepieces. -/
theorem field_lens (P₁ d P₂ : R) (h : d * P₁ = 1) : wAdd (-d) P₁ P₂ = P₁ := by
  simp only [wAdd_def]; linear_combination -P₂ * h

/-! ### Circuits and springs: series and parallel -/

/-- Parallel resistors (or springs in series, capacitors in series) are addition read through the
inversion lens `x ↦ x⁻¹`: `(R * S / (R + S)) * (R⁻¹ + S⁻¹) = 1`. -/
theorem parallel_is_inverted_addition (R S : K) (hR : R ≠ 0) (hS : S ≠ 0) (hRS : R + S ≠ 0) :
    R * S / (R + S) * (R⁻¹ + S⁻¹) = 1 := by
  field_simp; ring

/-- One section of the resistor ladder (`1 Ω` in series, then `1 Ω` across): a load `Z = n / d`
goes to `1 + Z / (1 + Z) = (2 * n + d) / (n + d)`. -/
def ladderStep (p : ℤ × ℤ) : ℤ × ℤ := (2 * p.1 + p.2, p.1 + p.2)

/-- The ladder with `n` sections, ending in a short circuit, as a numerator-denominator pair. -/
def ladder (n : ℕ) : ℤ × ℤ := ladderStep^[n] (0, 1)

@[simp] theorem ladder_zero : ladder 0 = (0, 1) := rfl

theorem ladder_succ (n : ℕ) : ladder (n + 1) = ladderStep (ladder n) :=
  Function.iterate_succ_apply' ladderStep n (0, 1)

/-- `n + 1` sections ending in a short circuit give `F (2 * n + 2) / F (2 * n + 1)`:
`1, 3/2, 8/5, 21/13, …`. -/
theorem ladder_fib (n : ℕ) :
    ladder (n + 1) = ((Nat.fib (2 * n + 2) : ℤ), (Nat.fib (2 * n + 1) : ℤ)) := by
  induction n with
  | zero => rfl
  | succ n ih =>
    have h₃ : Nat.fib (2 * n + 3) = Nat.fib (2 * n + 1) + Nat.fib (2 * n + 2) := Nat.fib_add_two
    have h₄ : Nat.fib (2 * n + 4) = Nat.fib (2 * n + 2) + Nat.fib (2 * n + 3) := Nat.fib_add_two
    rw [ladder_succ, ih, ladderStep, show 2 * (n + 1) + 2 = 2 * n + 4 by ring,
      show 2 * (n + 1) + 1 = 2 * n + 3 by ring, h₄, h₃]
    push_cast
    exact Prod.ext (by ring) (by ring)

/-- The walls of one section are its fixed points: `Z = (1 + 2 * Z) / (1 + Z)` exactly when
`Z ^ 2 = Z + 1`, the golden ratio and minus its inverse (two real walls). -/
theorem ladder_walls (n d : R) : (2 * n + d) * d = n * (n + d) ↔ n ^ 2 = n * d + d ^ 2 := by
  constructor <;> intro h <;> linear_combination -h

/-- Every ladder value sits on the hyperbola `n ^ 2 - n * d - d ^ 2 = -1`, which hugs the golden
wall: `Z ^ 2 - Z - 1 = -1 / d ^ 2`, and the ladder walks into the wall. -/
theorem ladder_hyperbola (k : ℕ) :
    (ladder k).1 ^ 2 - (ladder k).1 * (ladder k).2 - (ladder k).2 ^ 2 = -1 := by
  induction k with
  | zero => rfl
  | succ k ih => rw [ladder_succ, ladderStep]; linear_combination ih

end SimpleLaws
