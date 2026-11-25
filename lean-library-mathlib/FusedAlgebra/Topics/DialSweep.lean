module

public import FusedAlgebra.Topics.LinearAlgebra

/-!
# Turning the dial in linear algebra

The weight `k` in `wmadd k A B = A + B + k • (A * B)` is a dial. This file is the mathlib form of
`FusedAlgebra.Topics.DialSweep`; it follows the dial through determinants, characteristic
polynomials and a family of 2×2 matrices.

## Main results

* `smul_wmadd`: every weight is one law in other units, `k • wmadd k x y = (k • x) ✳ (k • y)`.
* `det_one_add_smul_eq_charpoly`: `det (1 + k • A) = (-k) ^ n * A.charpoly.eval (-1 / k)`, the
  reversed characteristic polynomial, so the matrix meets the wall of weight `k` exactly when
  `-1 / k` is a root of the characteristic polynomial (`wall_iff_isRoot_charpoly`);
  `det_one_add_neg_one_smul` is the case `k = -1`.
* `dialU k q = !![0, -q; 1, k]`, a matrix on the dial: `dial_product`, `dial_det`, `dial_lens`
  lift the lens law `(a + b + k * a * b) / (1 - q * a * b)`; `unitU κ` gives rotations
  (`κ = -1`), shears (`κ = 0`) and boosts (`κ = 1`).
* `kmadd K x y = x + y + x * K * y`, a matrix weight inside the law, with its lift
  `one_add_kmadd_mul`.
* `det_one_add_vecMulVec`, `det_one_add_mul_vecMulVec`: a rank-one weight leaves no corner
  terms (the matrix determinant lemma).

All of this is known mathematics. In mathlib: `Matrix.eval_charpoly`, `Matrix.charpoly_fin_two`,
`Matrix.det_one_add_mul_comm` (the Weinstein–Aronszajn / Sylvester identity, the original's
`sylvester2` and `sylvester3`) and `Matrix.det_one_add_replicateCol_mul_replicateRow` (the matrix
determinant lemma).
-/

@[expose] public section

namespace DialSweep

open MultAdd LinearAlgebraMA Matrix Polynomial

/-! ### One law in other units -/

section Algebra

variable {R A : Type*} [CommSemiring R] [Semiring A] [Algebra R A]

/-- Every weight is one law in other units: scaling by `k` turns weighted mult-add into plain
mult-add. -/
theorem smul_wmadd (k : R) (x y : A) : k • wmadd k x y = (k • x) ✳ (k • y) := by
  simp only [wmadd, multAdd, smul_add, smul_mul_smul_comm, smul_smul]

end Algebra

/-! ### A matrix weight inside the law -/

section Semiring

variable {A : Type*} [Semiring A]

/-- Mult-add with a weight `K` inside: `x + y + x * K * y`. -/
def kmadd (K x y : A) : A := x + y + x * K * y

/-- The lift of `kmadd`: `1 + kmadd K x y * K = (1 + x * K) * (1 + y * K)`. -/
theorem one_add_kmadd_mul (K x y : A) :
    1 + kmadd K x y * K = (1 + x * K) * (1 + y * K) := by
  simp only [kmadd, add_mul, mul_add, one_mul, mul_one, mul_assoc]
  abel

end Semiring

/-! ### The dial is the characteristic polynomial -/

section Matrix

variable {n R : Type*} [Fintype n] [DecidableEq n] [CommRing R]

/-- `k = -1` is the characteristic polynomial at `1` (Fredholm's `λ` is `-k`). -/
theorem det_one_add_neg_one_smul (A : Matrix n n R) :
    det (1 + (-1 : R) • A) = A.charpoly.eval 1 := by
  rw [eval_charpoly, map_one, neg_one_smul, sub_eq_add_neg]

/-- At `k = -1` the rotation condition is `Xᵀ + X - Xᵀ * X = 0`: `1 - X` is orthogonal. -/
theorem transpose_one_sub_mul (X : Matrix n n R) :
    (1 - X)ᵀ * (1 - X) = 1 - wmadd (-1 : R) Xᵀ X := by
  have h := transpose_one_add_smul_mul (-1 : R) X
  rwa [neg_one_smul, neg_one_smul, ← sub_eq_add_neg, ← sub_eq_add_neg] at h

/-- A rank-one weight leaves no corner terms: `det (1 + u * vᵀ) = 1 + v ⬝ᵥ u`. -/
theorem det_one_add_vecMulVec (u v : n → R) : det (1 + vecMulVec u v) = 1 + v ⬝ᵥ u := by
  rw [vecMulVec_eq Unit, det_one_add_replicateCol_mul_replicateRow]

/-- The matrix determinant lemma: `det (1 + A * (u * vᵀ)) = 1 + v ⬝ᵥ A *ᵥ u`. -/
theorem det_one_add_mul_vecMulVec (A : Matrix n n R) (u v : n → R) :
    det (1 + A * vecMulVec u v) = 1 + v ⬝ᵥ A *ᵥ u := by
  have e : A * vecMulVec u v = vecMulVec (A *ᵥ u) v := by
    ext i j
    simp [vecMulVec_apply, Matrix.mul_apply, mulVec, dotProduct, Finset.sum_mul, mul_assoc]
  rw [e, det_one_add_vecMulVec]

end Matrix

section Small

variable {R : Type*} [CommRing R]

/-- The characteristic polynomial of a 2×2 matrix: `t ^ 2 - trace A * t + det A`. -/
theorem eval_charpoly_fin_two (A : Matrix (Fin 2) (Fin 2) R) (t : R) :
    A.charpoly.eval t = t ^ 2 - trace A * t + det A := by
  rw [eval_charpoly]
  simp [det_fin_two, trace_fin_two]
  ring

/-- The characteristic polynomial of a 3×3 matrix: `t ^ 3 - e₁ * t ^ 2 + e₂ * t - e₃`. -/
theorem eval_charpoly_fin_three (A : Matrix (Fin 3) (Fin 3) R) (t : R) :
    A.charpoly.eval t = t ^ 3 - trace A * t ^ 2 + e2 A * t - det A := by
  rw [eval_charpoly]
  simp [det_fin_three, trace_fin_three, e2]
  ring

/-- The top of the dial: `det (1 + k • A) - k ^ 3 * det A` has only lower powers of `k`. -/
theorem det_top (k : R) (A : Matrix (Fin 3) (Fin 3) R) :
    det (1 + k • A) - k ^ 3 * det A = 1 + k * trace A + k ^ 2 * e2 A := by
  rw [det_one_add_smul_fin_three]; ring

/-! ### A matrix on the dial -/

/-- The dial matrix `U = !![0, -q; 1, k]`, with `U * U = k • U - q • 1`. -/
def dialU (k q : R) : Matrix (Fin 2) (Fin 2) R := !![0, -q; 1, k]

theorem dialU_mul_self (k q : R) : dialU k q * dialU k q = k • dialU k q - q • 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [dialU, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- `(1 + a • U) * (1 + b • U) = (1 - q * a * b) • 1 + (a + b + k * a * b) • U`: the landscape
law, lifted. -/
theorem dial_product (k q a b : R) :
    (1 + a • dialU k q) * (1 + b • dialU k q) =
      (1 - q * a * b) • 1 + (a + b + k * a * b) • dialU k q := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [dialU, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- `det (1 + a • U) = 1 + k * a + q * a ^ 2`, the landscape's wall polynomial. -/
theorem dial_det (k q a : R) : det (1 + a • dialU k q) = 1 + k * a + q * a ^ 2 := by
  simp [dialU, det_fin_two]; ring

/-- With eigenvalues `l₁`, `l₂` (`k = l₁ + l₂`, `q = l₁ * l₂`) the walls are `a = -1 / lᵢ`. -/
theorem dial_det_eigen (l₁ l₂ a : R) :
    det (1 + a • dialU (l₁ + l₂) (l₁ * l₂)) = (1 + a * l₁) * (1 + a * l₂) := by
  rw [dial_det]; ring

/-- The three countries are the sign of the discriminant, the squared eigenvalue gap. -/
theorem dial_discriminant (l₁ l₂ : R) : (l₁ + l₂) ^ 2 - 4 * (l₁ * l₂) = (l₁ - l₂) ^ 2 := by
  ring

/-- The "or" line: a wall at certainty `a = 1` exactly when `1 + k + q = 0`. -/
theorem or_line (k q : R) : det (1 + (1 : R) • dialU k q) = 1 + k + q := by
  rw [dial_det]; ring

/-- The pure unit `U = !![0, κ; 1, 0]`, with `U * U = κ • 1`. -/
def unitU (κ : R) : Matrix (Fin 2) (Fin 2) R := !![0, κ; 1, 0]

theorem unitU_mul_self (κ : R) : unitU κ * unitU κ = κ • 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [unitU, Matrix.mul_apply, Fin.sum_univ_two]

/-- `(1 + a • U) * (1 + b • U) = (1 + κ * a * b) • 1 + (a + b) • U`: projected,
`(a + b) / (1 + κ * a * b)`. -/
theorem unit_product (κ a b : R) :
    (1 + a • unitU κ) * (1 + b • unitU κ) = (1 + κ * a * b) • 1 + (a + b) • unitU κ := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [unitU, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- `κ = 1` (boosts): `det (1 + a • U) = (1 - a) * (1 + a)`, walls at the speed of light `±1`. -/
theorem boost_walls (a : R) : det (1 + a • unitU 1) = (1 - a) * (1 + a) := by
  simp [unitU, det_fin_two]; ring

/-- `κ = -1` (rotations): `det (1 + a • U) = 1 + a ^ 2`, no real wall. -/
theorem rotation_no_wall (a : R) : det (1 + a • unitU (-1)) = 1 + a ^ 2 := by
  simp [unitU, det_fin_two]; ring

/-- `κ = 0` (shears): the law is plain addition. -/
theorem shear_is_add (a b : R) :
    (1 + a • unitU 0) * (1 + b • unitU 0) = 1 + (a + b) • unitU (0 : R) := by
  rw [unit_product]; simp

end Small

/-! ### Fields: the walls, the reversed polynomial, the projected laws -/

section Field

variable {n K : Type*} [Fintype n] [DecidableEq n] [Field K]

/-- `det (1 + k • A)` is the reversed characteristic polynomial:
`det (1 + k • A) = (-k) ^ n * A.charpoly.eval (-1 / k)` for `k ≠ 0`. -/
theorem det_one_add_smul_eq_charpoly (k : K) (hk : k ≠ 0) (A : Matrix n n K) :
    det (1 + k • A) = (-k) ^ Fintype.card n * A.charpoly.eval (-1 / k) := by
  have e : 1 + k • A = (-k) • (scalar n (-1 / k) - A) := by
    rw [scalar_apply, ← smul_one_eq_diagonal, smul_sub, smul_smul]
    have h : -k * (-1 / k) = 1 := by field_simp
    rw [h, one_smul, neg_smul, sub_neg_eq_add]
  rw [eval_charpoly, e, det_smul]

/-- The matrix sits on the wall of weight `k` exactly when `-1 / k` is a root of its
characteristic polynomial, that is, an eigenvalue. -/
theorem wall_iff_isRoot_charpoly (k : K) (hk : k ≠ 0) (A : Matrix n n K) :
    det (1 + k • A) = 0 ↔ A.charpoly.IsRoot (-1 / k) := by
  rw [det_one_add_smul_eq_charpoly k hk, mul_eq_zero, IsRoot.def,
    or_iff_right (pow_ne_zero _ (neg_ne_zero.2 hk))]

/-- The landscape lens law comes out of the dial matrix: off the wall `1 - q * a * b ≠ 0`,
`(1 + a • U) * (1 + b • U) = (1 - q * a * b) • (1 + c • U)` with
`c = (a + b + k * a * b) / (1 - q * a * b)`. -/
theorem dial_lens (k q a b : K) (h : 1 - q * a * b ≠ 0) :
    (1 + a • dialU k q) * (1 + b • dialU k q) =
      (1 - q * a * b) • (1 + ((a + b + k * a * b) / (1 - q * a * b)) • dialU k q) := by
  rw [dial_product, smul_add, smul_smul, mul_div_cancel₀ _ h]

theorem smul_skew (k s : K) : k • skew s = skew (k * s) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [skew]

/-- The weighted Cayley transform `(1 + k • S)⁻¹ * (1 - k • S)` is the Cayley transform of
`k * s`. -/
theorem cayley_weighted (k s : K) :
    (1 + k • skew s)⁻¹ * (1 - k • skew s) = cayley (k * s) := by
  rw [smul_skew, cayley]

/-- The weighted inverse: off the weighted wall, `wAdd k a (-a / (1 + k * a)) = 0`. -/
theorem wAdd_neg_div (k a : K) (h : 1 + k * a ≠ 0) : wAdd k a (-a / (1 + k * a)) = 0 := by
  have e : wAdd k a (-a / (1 + k * a)) = a * (1 - (1 + k * a) / (1 + k * a)) := by
    unfold wAdd; ring
  rw [e, div_self h, sub_self, mul_zero]

end Field

end DialSweep
