module

public import FusedAlgebra.Topics.NumberPlane

/-!
# Geometry and trigonometry out of mult-add

* Geometry: the corner square `a * b` of `(1 + a) * (1 + b)` is Euclid's gnomon (*Elements* II.4).
  Mult-add of `a` with itself is the L-shaped strip that grows a square (`gnomon`), mult-add three
  times gives the cube's layers `3, 3, 1` (`cube_layers`), and Galileo's falling-body distances
  `1, 3, 5, 7, …` are exactly those gnomons (`fall_is_gnomon`, `fallen_square`).
* Pythagoras: `1 + x * i` is the right triangle with legs `1` and `x`. Multiplying two of them
  adds their angles; the slopes combine by tangent addition (the landscape point `(0, 1)`) and the
  squared hypotenuses multiply (`triangles_compose`). Squaring one gives Euclid's triples
  (`euclid_triple`).
* Trigonometry: the sine and cosine addition laws are the multiplication at `(0, 1)`
  (`angle_addition`); cubing `1 + x * i` sends the cube's layers `1, 3, 3, 1` alternately to the
  two axes (`triple_angle`); the law of cosines is the norm on the line `q = 1`
  (`law_of_cosines_norm`), and real triangles exist exactly in the turning country of that line
  (`triangle_lengths_real`).

All of this is classical. The statements are polynomial identities, valid in any commutative
(semi)ring; the original worked over `ℤ`. Over `ℝ` the complex corner is mathlib's `ℂ`:
`mre_one_eq_re_mul`, `mu_zero_eq_im_mul` and `norm_zero_one_eq_normSq` say so, and the analytic
versions of the addition laws are `Real.cos_add`, `Real.sin_add` and `Real.tan_add`.
-/

@[expose] public section

namespace SimpleLaws

open MultAdd Finset

/-! ### Geometry -/

section CommSemiring

variable {R : Type*} [CommSemiring R]

/-- The rectangle: `(1 + a) * (1 + b)` is the unit square, two strips and the corner square. -/
theorem rectangle_corner (a b : R) : (1 + a) * (1 + b) = 1 + a + b + a * b := by ring

/-- Euclid II.4, the gnomon: mult-add of `a` with itself is the L-shaped strip `2 * a + a ^ 2`. -/
theorem gnomon (a : R) : a ✳ a = 2 * a + a ^ 2 := by
  simp only [multAdd_def]; ring

/-- Euclid II.4 in full: `(a + b) ^ 2 = a ^ 2 + 2 * a * b + b ^ 2`. This is mathlib's `add_sq`. -/
theorem euclid_II_4 (a b : R) : (a + b) ^ 2 = a ^ 2 + 2 * a * b + b ^ 2 := add_sq a b

/-- The cube: mult-add three times gives the layers `3, 3, 1` (with the unit cube,
`1, 3, 3, 1`). -/
theorem cube_layers (a : R) : a ✳ a ✳ a = 3 * a + 3 * a ^ 2 + a ^ 3 := by
  simp only [multAdd_def]; ring

/-- The box: three different sides give the corners `a + b + c`, `a * b + b * c + c * a`,
`a * b * c` (the elementary symmetric functions). -/
theorem box_corners (a b c : R) :
    a ✳ b ✳ c = a + b + c + (a * b + b * c + c * a) + a * b * c := by
  simp only [multAdd_def]; ring

/-- Galileo's square-cube law (1638): grow every length by `a`; area grows by `A = a ✳ a` and
volume by `A ✳ a`, so weight per unit of strength grows by the factor `1 + a`. -/
theorem square_cube (a : R) : 1 + (a ✳ a ✳ a) = (1 + (a ✳ a)) * (1 + a) :=
  one_add_multAdd _ _

end CommSemiring

/-- Galileo's odd-number rule: the distance fallen in the `n`-th of a run of equal times. -/
def fall : ℕ → ℕ
  | 0 => 0
  | n + 1 => 2 * n + 1

@[simp] theorem fall_zero : fall 0 = 0 := rfl

@[simp] theorem fall_succ (n : ℕ) : fall (n + 1) = 2 * n + 1 := rfl

/-- The total distance fallen after `n` equal times. -/
def fallen (n : ℕ) : ℕ := ∑ i ∈ range n, fall (i + 1)

/-- Each step's distance is the gnomon that grows an `n`-square into an `(n + 1)`-square. -/
theorem fall_is_gnomon (n : ℕ) : (n + 1) ^ 2 = n ^ 2 + fall (n + 1) := by
  rw [fall_succ]; ring

/-- So the distance fallen is the square of the time. -/
theorem fallen_square (n : ℕ) : fallen n = n ^ 2 := by
  induction n with
  | zero => rfl
  | succ n ih => rw [fallen, sum_range_succ, ← fallen, ih, ← fall_is_gnomon]

/-- The gnomon step is mult-add's lift on `ℕ`: `n ✳ n + 1 = (n + 1) ^ 2`. -/
theorem gnomon_lift (n : ℕ) : n ✳ n + 1 = (n + 1) ^ 2 := by
  rw [multAdd_add_one, sq]

section CommRing

variable {R : Type*} [CommRing R]

/-- Pythagoras by the "bride's chair": the square on `a + b` holds the tilted square on `c` and
four right triangles of area `a * b / 2`. Given that area count, `c ^ 2 = a ^ 2 + b ^ 2`. -/
theorem pythagoras_from_areas (a b c : R) (h : (a + b) ^ 2 = c ^ 2 + 2 * a * b) :
    c ^ 2 = a ^ 2 + b ^ 2 := by
  linear_combination -h

/-- The right triangle `1 + x * i`: multiplying two of them multiplies the squared hypotenuses,
`(1 + x ^ 2) * (1 + y ^ 2) = (1 - x * y) ^ 2 + (x + y) ^ 2`; the new slope is
`(x + y) / (1 - x * y)`. -/
theorem triangles_compose (x y : R) :
    (1 + x ^ 2) * (1 + y ^ 2) = (1 - x * y) ^ 2 + (x + y) ^ 2 := by ring

/-- The same, as an instance of `norm_mul` at the complex corner. -/
theorem triangles_compose_norm (x y : R) :
    norm 0 1 (mre 1 1 x 1 y) (mu 0 1 x 1 y) = norm 0 1 1 x * norm 0 1 1 y :=
  norm_mul 0 1 1 x 1 y

/-- Euclid's triples: square the triangle `q + p * i`. -/
theorem euclid_triple_parts (p q : R) :
    mre 1 q p q p = q ^ 2 - p ^ 2 ∧ mu 0 q p q p = 2 * p * q := by
  constructor <;> simp only [mre_def, mu_def] <;> ring

/-- Euclid's triples are Pythagorean: `(q ^ 2 - p ^ 2, 2 * p * q, q ^ 2 + p ^ 2)`. -/
theorem euclid_triple (p q : R) :
    (q ^ 2 - p ^ 2) ^ 2 + (2 * p * q) ^ 2 = (q ^ 2 + p ^ 2) ^ 2 := by ring

example : 3 ^ 2 + 4 ^ 2 = 5 ^ 2 ∧ 5 ^ 2 + 12 ^ 2 = 13 ^ 2 := by decide

/-! ### Trigonometry -/

/-- Angle addition: `(cos α + i * sin α) * (cos β + i * sin β)` at the complex point `(0, 1)`. -/
theorem angle_addition (c₁ s₁ c₂ s₂ : R) :
    mre 1 c₁ s₁ c₂ s₂ = c₁ * c₂ - s₁ * s₂ ∧ mu 0 c₁ s₁ c₂ s₂ = c₁ * s₂ + s₁ * c₂ := by
  constructor <;> simp only [mre_def, mu_def] <;> ring

/-- `cos ^ 2 + sin ^ 2 = 1` survives angle addition: the unit circle is closed under the
product. -/
theorem unit_circle_closed (c₁ s₁ c₂ s₂ : R) (h₁ : c₁ ^ 2 + s₁ ^ 2 = 1)
    (h₂ : c₂ ^ 2 + s₂ ^ 2 = 1) :
    (c₁ * c₂ - s₁ * s₂) ^ 2 + (c₁ * s₂ + s₁ * c₂) ^ 2 = 1 := by
  linear_combination (c₂ ^ 2 + s₂ ^ 2) * h₁ + h₂

/-- Tangent addition is the landscape law at `(0, 1)`: numerator `x + y`, denominator
`1 - x * y`. -/
theorem tangent_addition (x y : R) :
    LandscapeMap.num 0 x y = x + y ∧ LandscapeMap.den 1 x y = 1 - x * y := by
  simp [LandscapeMap.num, LandscapeMap.den]

/-- Double angle: `(1 + x * i) ^ 2 = (1 - x ^ 2) + 2 * x * i`, so `tan 2θ = 2 * t / (1 - t ^ 2)`. -/
theorem double_angle (x : R) : mre 1 1 x 1 x = 1 - x ^ 2 ∧ mu 0 1 x 1 x = 2 * x := by
  constructor <;> simp only [mre_def, mu_def] <;> ring

/-- Triple angle: `(1 + x * i) ^ 3 = (1 - 3 * x ^ 2) + (3 * x - x ^ 3) * i`. The cube's layers
`1, 3, 3, 1` land alternately on the real and the `i` axis, with the signs of `i ^ 2`, `i ^ 3`. -/
theorem triple_angle (x : R) :
    mre 1 (mre 1 1 x 1 x) (mu 0 1 x 1 x) 1 x = 1 - 3 * x ^ 2 ∧
      mu 0 (mre 1 1 x 1 x) (mu 0 1 x 1 x) 1 x = 3 * x - x ^ 3 := by
  constructor <;> simp only [mre_def, mu_def] <;> ring

/-- The law of cosines is the norm on the line `q = 1`, with `k = -2 * cos C`:
`c ^ 2 = a ^ 2 + k * a * b + b ^ 2`. At `k = 0` it is Pythagoras. -/
theorem law_of_cosines_norm (k a b : R) : norm k 1 a b = a ^ 2 + k * a * b + b ^ 2 := by
  simp only [norm_def]; ring

theorem law_of_cosines_right_angle (a b : R) : norm 0 1 a b = a ^ 2 + b ^ 2 := norm_complex a b

/-- On the border `k = 2` (angle `180°`) the triangle is flat: `c = a + b`. -/
theorem flat_triangle_open (a b : R) : norm 2 1 a b = (a + b) ^ 2 := by
  simp only [norm_def]; ring

/-- On the border `k = -2` (angle `0°`) the triangle is flat: `c = a - b`. -/
theorem flat_triangle_shut (a b : R) : norm (-2) 1 a b = (a - b) ^ 2 := by
  simp only [norm_def]; ring

/-- Inside the turning country (`k ^ 2 ≤ 4`) every side-squared is a true length, never
negative. -/
theorem triangle_lengths_real [LinearOrder R] [IsStrictOrderedRing R] (k a b : R)
    (h : k ^ 2 ≤ 4) : 0 ≤ 4 * norm k 1 a b :=
  norm_nonneg_turning k 1 a b (by rwa [mul_one])

/-- Past the border (`k ^ 2 > 4`, "`cos C`" beyond `±1`) a side-squared can go negative: the
triangles of relativity, where cosine becomes `cosh`. -/
example : norm (3 : ℤ) 1 1 (-1) = -1 := by decide

end CommRing

/-! ### The complex corner over `ℝ` is mathlib's `ℂ` -/

/-- At `(0, 1)` over `ℝ`, `mre` is the real part of the product in `ℂ`. -/
theorem mre_one_eq_re_mul (a b c d : ℝ) : mre 1 a b c d = ((⟨a, b⟩ * ⟨c, d⟩ : ℂ)).re := by
  simp [mre_def]

/-- At `(0, 1)` over `ℝ`, `mu` is the imaginary part of the product in `ℂ`. -/
theorem mu_zero_eq_im_mul (a b c d : ℝ) : mu 0 a b c d = ((⟨a, b⟩ * ⟨c, d⟩ : ℂ)).im := by
  simp [mu_def]

/-- At `(0, 1)` over `ℝ`, the norm is `Complex.normSq`. -/
theorem norm_zero_one_eq_normSq (a b : ℝ) : norm 0 1 a b = Complex.normSq ⟨a, b⟩ := by
  simp [norm_def, Complex.normSq_apply]

end SimpleLaws
