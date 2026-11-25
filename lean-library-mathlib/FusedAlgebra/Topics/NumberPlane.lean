module

public import FusedAlgebra.Topics.Landscape

/-!
# Every point of the landscape is mult-add in a number system

Pick a new unit `u` with `u ^ 2 = k * u - q`. Numbers are `a + b * u`. Multiplying the two
"tilted ones" `1 + x * u` and `1 + y * u` gives

`(1 + x * u) * (1 + y * u) = (1 - q * x * y) + (x + y + k * x * y) * u`,

and read as a slope (`u`-part over real part) this is exactly the landscape law at `(k, q)`.

The number system is mathlib's `QuadraticAlgebra R (-q) k` (elements `re + im • ω` with
`ω * ω = -q + k • ω`). The core-Lean original wrote the two coordinates of a product and the norm
as functions of four numbers; they are kept here as `mre`, `mu` and `norm` because the later files
quote them, and `mre_eq_re_mul`, `mu_eq_im_mul` and `norm_eq_norm` identify them with the
coordinates of the product and the norm `QuadraticAlgebra.norm` of mathlib.

## Main results

* `mul_rule`: the multiplication really is "`u ^ 2 = k * u - q`".
* `tilted_re`, `tilted_u`: the product of tilted ones is the landscape law.
* `norm_mul`: the norm multiplies (the generalised Brahmagupta–Fibonacci identity). This is
  `map_mul` for mathlib's `QuadraticAlgebra.norm`.
* `four_norm`, `norm_nonneg_turning`, `norm_border`, `norm_walls`: the sign of `4 * q - k ^ 2`
  picks the geometry: Euclid (`k ^ 2 < 4 * q`), Galileo (`k ^ 2 = 4 * q`), Einstein
  (`k ^ 2 > 4 * q`). `four_norm` is mathlib's `QuadraticAlgebra.four_mul_norm_eq` in coordinates.
* `norm_multAdd_lift`, `unit_circle_is_hyperbola`: at mult-add's own point `(1, 0)` the norm of a
  tilted one is its lift, and the unit circle is the hyperbola `x ✳ y = 0`.

This is classical (Clifford 1873, Study, Yaglom's *A Simple Non-Euclidean Geometry*): the three
corners `(0, 1)`, `(0, 0)`, `(0, -1)` are the complex, dual and split-complex numbers. Everything
holds in any commutative ring; the original worked over `ℤ`.
-/

@[expose] public section

namespace SimpleLaws

open MultAdd

variable {R : Type*} [CommRing R]

/-! ### The number system at `(k, q)`: `a + b * u` with `u ^ 2 = k * u - q` -/

/-- The real part of `(a + b * u) * (c + d * u)`. -/
def mre (q a b c d : R) : R := a * c - q * b * d

/-- The `u`-part of `(a + b * u) * (c + d * u)`. -/
def mu (k a b c d : R) : R := a * d + b * c + k * b * d

theorem mre_def (q a b c d : R) : mre q a b c d = a * c - q * b * d := rfl

theorem mu_def (k a b c d : R) : mu k a b c d = a * d + b * c + k * b * d := rfl

/-- `mre` is the real part of the product in mathlib's `QuadraticAlgebra R (-q) k`. -/
theorem mre_eq_re_mul (k q a b c d : R) :
    mre q a b c d = ((⟨a, b⟩ * ⟨c, d⟩ : QuadraticAlgebra R (-q) k)).re := by
  simp only [mre_def, QuadraticAlgebra.re_mul]; ring

/-- `mu` is the `ω`-part of the product in mathlib's `QuadraticAlgebra R (-q) k`. -/
theorem mu_eq_im_mul (k q a b c d : R) :
    mu k a b c d = ((⟨a, b⟩ * ⟨c, d⟩ : QuadraticAlgebra R (-q) k)).im := by
  simp only [mu_def, QuadraticAlgebra.im_mul]

/-- The product in `QuadraticAlgebra R (-q) k`, in the coordinates `mre` and `mu`. -/
theorem mk_mul_mk (k q a b c d : R) :
    (⟨a, b⟩ * ⟨c, d⟩ : QuadraticAlgebra R (-q) k) = ⟨mre q a b c d, mu k a b c d⟩ := by
  ext
  · exact (mre_eq_re_mul k q a b c d).symm
  · exact (mu_eq_im_mul k q a b c d).symm

/-- The multiplication really is "`u ^ 2 = k * u - q`": expand `(a + b * u) * (c + d * u)` and
replace `u ^ 2`. -/
theorem mul_rule (k q a b c d u : R) (hu : u ^ 2 = k * u - q) :
    (a + b * u) * (c + d * u) = mre q a b c d + mu k a b c d * u := by
  simp only [mre_def, mu_def]; linear_combination b * d * hu

theorem mul_comm_re (q a b c d : R) : mre q a b c d = mre q c d a b := by
  simp only [mre_def]; ring

theorem mul_comm_u (k a b c d : R) : mu k a b c d = mu k c d a b := by
  simp only [mu_def]; ring

/-- Associativity, real part. -/
theorem mul_assoc_re (k q a b c d e f : R) :
    mre q (mre q a b c d) (mu k a b c d) e f = mre q a b (mre q c d e f) (mu k c d e f) := by
  have h := congrArg QuadraticAlgebra.re
    (mul_assoc (⟨a, b⟩ : QuadraticAlgebra R (-q) k) ⟨c, d⟩ ⟨e, f⟩)
  rwa [mk_mul_mk, mk_mul_mk, mk_mul_mk, mk_mul_mk] at h

/-- Associativity, `u`-part. -/
theorem mul_assoc_u (k q a b c d e f : R) :
    mu k (mre q a b c d) (mu k a b c d) e f = mu k a b (mre q c d e f) (mu k c d e f) := by
  have h := congrArg QuadraticAlgebra.im
    (mul_assoc (⟨a, b⟩ : QuadraticAlgebra R (-q) k) ⟨c, d⟩ ⟨e, f⟩)
  rwa [mk_mul_mk, mk_mul_mk, mk_mul_mk, mk_mul_mk] at h

/-! ### The main fact: the product of tilted ones is the landscape law -/

/-- The real part of `(1 + x * u) * (1 + y * u)` is the landscape denominator. -/
theorem tilted_re (q x y : R) : mre q 1 x 1 y = LandscapeMap.den q x y := by
  simp only [mre_def, LandscapeMap.den_def]; ring

/-- The `u`-part of `(1 + x * u) * (1 + y * u)` is the landscape numerator, the weighted mult-add
of `x` and `y`. -/
theorem tilted_u (k x y : R) : mu k 1 x 1 y = LandscapeMap.num k x y := by
  simp only [mu_def, LandscapeMap.num_def]; ring

/-- At mult-add's own point `(k, q) = (1, 0)` the unit is idempotent (`u ^ 2 = u`) and the
`u`-part is plain mult-add. -/
theorem multAdd_point_u (x y : R) : mu 1 1 x 1 y = x ✳ y := by
  simp only [mu_def, multAdd_def]; ring

/-- At mult-add's own point the real part is `1`: there is nothing to divide by. -/
theorem multAdd_point_re (x y : R) : mre 0 1 x 1 y = 1 := by simp [mre_def]

/-- The complex corner `(0, 1)`, `u = i`: slopes add like tangents, `(x + y) / (1 - x * y)`. -/
theorem complex_point (x y : R) : mu 0 1 x 1 y = x + y ∧ mre 1 1 x 1 y = 1 - x * y := by
  simp [mre_def, mu_def, add_comm]

/-- The split corner `(0, -1)`, `u ^ 2 = 1`: velocity addition, `(x + y) / (1 + x * y)`. -/
theorem split_point (x y : R) : mu 0 1 x 1 y = x + y ∧ mre (-1) 1 x 1 y = 1 + x * y := by
  simp [mre_def, mu_def, add_comm]

/-- The dual corner `(0, 0)`, `u ^ 2 = 0`: Galileo's velocity addition, plain `x + y` over `1`. -/
theorem dual_point (x y : R) : mu 0 1 x 1 y = x + y ∧ mre 0 1 x 1 y = 1 := by
  simp [mre_def, mu_def, add_comm]

/-! ### The norm: each point's own Pythagoras -/

/-- The norm of `a + b * u` at the point `(k, q)`: `a ^ 2 + k * a * b + q * b ^ 2`. -/
def norm (k q a b : R) : R := a * a + k * a * b + q * b * b

theorem norm_def (k q a b : R) : norm k q a b = a * a + k * a * b + q * b * b := rfl

/-- `norm` is the norm of mathlib's `QuadraticAlgebra R (-q) k`. -/
theorem norm_eq_norm (k q a b : R) :
    norm k q a b = QuadraticAlgebra.norm (⟨a, b⟩ : QuadraticAlgebra R (-q) k) := by
  simp only [norm_def, QuadraticAlgebra.norm_def]; ring

/-- The norm multiplies (the generalised Brahmagupta–Fibonacci identity). This is `map_mul` for
`QuadraticAlgebra.norm`. -/
theorem norm_mul (k q a b c d : R) :
    norm k q (mre q a b c d) (mu k a b c d) = norm k q a b * norm k q c d := by
  rw [norm_eq_norm, norm_eq_norm, norm_eq_norm, ← mk_mul_mk, map_mul]

/-- The norm is the number times its conjugate: with `ū = k - u`,
`(a + b * u) * (a + b * ū) = N`. -/
theorem norm_conj (k q a b u : R) (hu : u ^ 2 = k * u - q) :
    (a + b * u) * (a + b * (k - u)) = norm k q a b := by
  simp only [norm_def]; linear_combination -(b ^ 2) * hu

/-- `4 * N` splits as a square plus `(4 * q - k ^ 2) * b ^ 2`: the sign of `4 * q - k ^ 2` picks
the country. In mathlib this is `QuadraticAlgebra.four_mul_norm_eq`. -/
theorem four_norm (k q a b : R) :
    4 * norm k q a b = (2 * a + k * b) ^ 2 + (4 * q - k ^ 2) * b ^ 2 := by
  simp only [norm_def]; ring

/-- Border: when `k ^ 2 = 4 * q` the norm is one square, `(2 * a + k * b) ^ 2 / 4` (Galileo:
lengths only one way). -/
theorem norm_border (k q a b : R) (h : k ^ 2 = 4 * q) :
    4 * norm k q a b = (2 * a + k * b) ^ 2 := by
  rw [four_norm, h, sub_self, zero_mul, add_zero]

/-- Two walls: write `k = -(r + s)`, `q = r * s`; the norm splits into two light-cone factors. -/
theorem norm_walls (r s a b : R) : norm (-(r + s)) (r * s) a b = (a - r * b) * (a - s * b) := by
  simp only [norm_def]; ring

/-- Pythagoras at the complex corner: `a ^ 2 + b ^ 2`. -/
theorem norm_complex (a b : R) : norm 0 1 a b = a ^ 2 + b ^ 2 := by
  simp only [norm_def]; ring

/-- The spacetime interval at the split corner: `a ^ 2 - b ^ 2`. -/
theorem norm_split (a b : R) : norm 0 (-1) a b = a ^ 2 - b ^ 2 := by
  simp only [norm_def]; ring

/-- Galileo at the dual corner: only `a ^ 2` survives. -/
theorem norm_dual (a b : R) : norm 0 0 a b = a ^ 2 := by
  simp only [norm_def]; ring

/-- Turning country: when `k ^ 2 ≤ 4 * q` the norm is never negative (a true length). -/
theorem norm_nonneg_turning [LinearOrder R] [IsStrictOrderedRing R] (k q a b : R)
    (h : k ^ 2 ≤ 4 * q) : 0 ≤ 4 * norm k q a b := by
  rw [four_norm]
  exact add_nonneg (sq_nonneg _) (mul_nonneg (sub_nonneg.2 h) (sq_nonneg b))

/-! ### Mult-add's own geometry: its unit circle is Boyle's hyperbola -/

/-- At `(1, 0)` the norm is `a * (a + b)`: in the coordinates `X = a`, `Y = a + b` it is
`X * Y`. -/
theorem norm_multAdd (a b : R) : norm 1 0 a b = a * (a + b) := by
  simp only [norm_def]; ring

/-- The tilted one `1 + x * u` at mult-add's point has norm `1 + x`: its lift. -/
theorem norm_multAdd_tilted (x : R) : norm 1 0 1 x = 1 + x := by
  simp only [norm_def]; ring

/-- The product of two tilted ones has norm `(1 + x) * (1 + y) = 1 + (x ✳ y)`. -/
theorem norm_multAdd_lift (x y : R) : norm 1 0 (mre 0 1 x 1 y) (mu 1 1 x 1 y) = 1 + (x ✳ y) := by
  rw [norm_mul, norm_multAdd_tilted, norm_multAdd_tilted, one_add_multAdd]

/-- The unit circle `N = 1` at mult-add's point is the curve `x ✳ y = 0`, that is
`(1 + x) * (1 + y) = 1`: a hyperbola whose asymptotes are the walls `x = -1` and `y = -1`. -/
theorem unit_circle_is_hyperbola (x y : R) : x ✳ y = 0 ↔ (1 + x) * (1 + y) = 1 := by
  rw [← one_add_multAdd, add_eq_left]

end SimpleLaws
