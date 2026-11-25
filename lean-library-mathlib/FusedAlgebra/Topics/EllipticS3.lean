module

public import FusedAlgebra.Basic

/-!
# Why the primes' turning angles live on the 3-sphere

Every law on the landscape `(k, q)` is mult-add's lift in `2 × 2` matrices: the number `1 + a`
becomes the matrix `1 + a • G` with `G = !![k, 1; -q, 0]` (trace `k`, determinant `q`). In the
turning country `k² < 4q` the matrix has no eigenvalue in the base field, and the rotations with
a given trace sit on the 3-sphere of unit quaternions.

This is the mathlib form of `FusedAlgebra.Topics.EllipticS3`. The original wrote every matrix
entry as a separate identity over `ℤ`; here the matrices are mathlib's `Matrix (Fin 2) (Fin 2) R`
over any commutative ring and the quaternions are mathlib's `ℍ[ℝ]`.

## Main definitions

* `EllipticS3.gen k q`: the companion matrix `!![k, 1; -q, 0]`.
* `EllipticS3.quatMatrix q`: a real quaternion as a `2 × 2` complex matrix.
* `EllipticS3.walks`, `EllipticS3.dyck`: counts of `±1` walks (corners of the cube) ending at
  `0`, and of those that never dip below `0`.

## Main results

* `one_add_smul_gen_mul`: `(1 + a • G) * (1 + b • G) = (1 - q * a * b) • 1 + wAdd k a b • G`.
* `gen_mul_self`: Cayley–Hamilton for `G` (a special case of mathlib's
  `Matrix.aeval_self_charpoly`).
* `turning_no_wall`: when `k² < 4q` the polynomial `l² - kl + q` has no root.
* `det_quatMatrix`, `trace_quatMatrix`, `quat_discriminant`: the determinant is the quaternion
  norm, the trace is `2w`, and the wall discriminant is `-4` times the squared radius of the
  latitude 2-sphere.
* `circle_moments`, `sphere_moments`: the first moments are the central binomial coefficients
  (circle) and the Catalan numbers (3-sphere, the Sato–Tate moments).

All of this is known mathematics; the file records how it reads through mult-add.
-/

@[expose] public section

namespace EllipticS3

open MultAdd Matrix

variable {R : Type*} [CommRing R]

/-! ### Every law on the map is mult-add's lift in `2 × 2` matrices -/

/-- The generator `G = !![k, 1; -q, 0]`, the companion matrix of `X² - kX + q`. -/
def gen (k q : R) : Matrix (Fin 2) (Fin 2) R := !![k, 1; -q, 0]

@[simp]
theorem trace_gen (k q : R) : (gen k q).trace = k := by
  simp [gen, Matrix.trace_fin_two]

@[simp]
theorem det_gen (k q : R) : (gen k q).det = q := by
  simp [gen, Matrix.det_fin_two]

/-- Cayley–Hamilton for `G`: `G² = k • G - q • 1`. -/
theorem gen_mul_self (k q : R) : gen k q * gen k q = k • gen k q - q • 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [gen, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- The matrix lift: `(1 + a • G) * (1 + b • G) = (1 - q * a * b) • 1 + N • G` with
`N = wAdd k a b = a + b + k * a * b`. Dividing by `D = 1 - q * a * b` gives the law
`F = N / D` of the landscape point `(k, q)`. -/
theorem one_add_smul_gen_mul (k q a b : R) :
    (1 + a • gen k q) * (1 + b • gen k q) = (1 - q * a * b) • 1 + wAdd k a b • gen k q := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [gen, wAdd, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-! ### The turning country has no wall -/

/-- In the turning country `k² < 4q` the matrix `G` has no eigenvalue: `l² - kl + q ≠ 0`. -/
theorem turning_no_wall [LinearOrder R] [IsStrictOrderedRing R] (k q l : R)
    (h : k ^ 2 < 4 * q) : l ^ 2 - k * l + q ≠ 0 := by
  intro h0
  nlinarith [sq_nonneg (2 * l - k)]

/-! ### Bad primes: Frobenius on the lift is the `p`-th power -/

/-- Turning (non-split) torus: `u ^ (p + 1) = 1`, so `v = u ^ p` satisfies `v * u = 1` and the
mult-add of `v - 1` and `u - 1` is `0`. Frobenius sends the point to its mult-add inverse
(`a_p = -1`). -/
theorem nonsplit_frob_is_inverse (u v : R) (h : v * u = 1) : (v - 1) ✳ (u - 1) = 0 := by
  unfold multAdd; linear_combination h

/-- Split torus: `u ^ p = u`, so Frobenius fixes the lift (`a_p = 1`). -/
theorem split_frob_fixes (u v : R) (h : v = u) : (v - 1) ✳ 0 = u - 1 := by
  rw [multAdd_zero, h]

/-! ### Points are fixed points: the wall polynomial at `1` -/

/-- Lefschetz on the torus's cube of holes (layers `1, 2, 1`):
`1 - (α + β) + αβ = (1 - α)(1 - β)`. -/
theorem lefschetz_count (α β : R) : 1 - (α + β) + α * β = (1 - α) * (1 - β) := by
  ring

/-! ### The 3-sphere: unit quaternions as `2 × 2` matrices -/

open Quaternion in
/-- The quaternion `w + xi + yj + zk` as the complex matrix
`!![w + xi, y + zi; -y + zi, w - xi]`. -/
def quatMatrix (q : ℍ[ℝ]) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![⟨q.re, q.imI⟩, ⟨q.imJ, q.imK⟩; ⟨-q.imJ, q.imK⟩, ⟨q.re, -q.imI⟩]

open Quaternion in
/-- The determinant of the matrix of a quaternion is its squared norm
`w² + x² + y² + z²` (in particular it is real). -/
theorem det_quatMatrix (q : ℍ[ℝ]) : (quatMatrix q).det = ((Quaternion.normSq q : ℝ) : ℂ) := by
  rw [Quaternion.normSq_def']
  apply Complex.ext <;>
    simp [quatMatrix, Matrix.det_fin_two, -Complex.ofReal_pow] <;> ring

open Quaternion in
/-- The trace of the matrix of a quaternion is `2w`. -/
theorem trace_quatMatrix (q : ℍ[ℝ]) : (quatMatrix q).trace = ((2 * q.re : ℝ) : ℂ) := by
  apply Complex.ext
  · simp [quatMatrix, Matrix.trace_fin_two]; ring
  · simp [quatMatrix, Matrix.trace_fin_two]

open Quaternion in
/-- The wall discriminant `k² - 4q` of a quaternion (trace `k = 2w`, determinant `q` the squared
norm) is `-4` times the squared radius of its latitude 2-sphere: how deep a prime sits in the
turning country is how big the sphere of rotations with its angle is. -/
theorem quat_discriminant (q : ℍ[ℝ]) :
    (2 * q.re) ^ 2 - 4 * Quaternion.normSq q = -4 * (q.imI ^ 2 + q.imJ ^ 2 + q.imK ^ 2) := by
  rw [Quaternion.normSq_def']; ring

/-! ### Moments through the cube

Corners of the `m`-cube read as `±1` walks. The circle keeps every balanced corner (the middle
Pascal layer), the 3-sphere keeps only the Dyck ones (Catalan numbers). -/

/-- `walks m h` counts the `±1` walks of `m` steps from height `h` that end at `0`. -/
def walks : ℕ → ℤ → ℕ
  | 0, h => if h = 0 then 1 else 0
  | m + 1, h => walks m (h + 1) + walks m (h - 1)

/-- `dyck m h` counts the `±1` walks of `m` steps from height `h` that end at `0` and never dip
below `0`. -/
def dyck : ℕ → ℕ → ℕ
  | 0, h => if h = 0 then 1 else 0
  | m + 1, h => dyck m (h + 1) + (if h = 0 then 0 else dyck m (h - 1))

/-- The circle's moments are the central binomial coefficients `2, 6, 20, 70`. -/
theorem circle_moments : ∀ n < 5, walks (2 * n) 0 = n.centralBinom := by
  decide

/-- The 3-sphere's moments are the Catalan numbers `1, 2, 5, 14`. -/
theorem sphere_moments : ∀ n < 5, dyck (2 * n) 0 = catalan n := by
  simp_rw [catalan_eq_centralBinom_div]
  decide

/-- The Dyck count is the middle layer minus its neighbour: `20 - 15 = 5` and `70 - 56 = 14`. -/
theorem middle_minus_neighbour :
    Nat.choose 6 3 - Nat.choose 6 2 = dyck 6 0 ∧ Nat.choose 8 4 - Nat.choose 8 3 = dyck 8 0 := by
  decide

end EllipticS3
