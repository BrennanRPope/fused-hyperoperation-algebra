module

public import FusedAlgebra.Basic

/-!
# Linear algebra through mult-add

Matrices form a (noncommutative) ring, so mult-add `A ✳ B = A + B + A * B` makes sense for them
and `1 + (A ✳ B) = (1 + A) * (1 + B)` is `MultAdd.one_add_multAdd`. This file is the mathlib form
of `FusedAlgebra.Topics.LinearAlgebra`. The core-Lean original wrote 2×2 and 3×3 matrices as
records; here they are mathlib's `Matrix (Fin 2) (Fin 2) R` and `Matrix (Fin 3) (Fin 3) R`, and
whatever does not depend on the size is stated for `Matrix n n R` or for any algebra.

## Main definitions

* `wmadd k x y = x + y + k • (x * y)`: weighted mult-add in an `R`-algebra (for matrices,
  `A + B + k • (A * B)`). In a commutative ring this is `MultAdd.wAdd`.
* `dialTwo k A = trace A + k * det A` and `dialThree k A = trace A + k * e2 A + k ^ 2 * det A`:
  the dial invariant, with `det (1 + k • A) = 1 + k * dial k A`. Here `e2` is the sum of the three
  principal 2×2 minors.
* `msubM X A = (X - A) * (1 + A)⁻¹`: matrix mult-add subtraction.
* `skew`, `cayley`: the 2×2 skew matrix and its Cayley transform.
* `fl M v = M *ᵥ (1 + v) - 1`: a fused-linear map, a linear map about the wall corner `-1`.
* `mono p q x₁ x₂ = (1 + x₁) ^ p * (1 + x₂) ^ q - 1`: a monomial map in mult-add coordinates.

## Main results

* `one_add_smul_wmadd`: the weighted lift `1 + k • wmadd k x y = (1 + k • x) * (1 + k • y)`.
* `shifted_det_hom`: `det (1 + (A ✳ B)) - 1 = (det (1 + A) - 1) ✳ (det (1 + B) - 1)`.
* `det_one_add_smul_fin_two`, `det_one_add_smul_fin_three`: `det (1 + k • A)` by layers.
* `dialTwo_wmadd`, `dialThree_wmadd`: the dial homomorphism; `k = 0` is "trace is additive",
  `k = 1` is "`det (1 + A)` is multiplicative".
* `det_one_add_of_isUpperTriangular`: for a triangular matrix `det (1 + A) = ∏ i, (1 + A i i)`,
  so the shifted determinant is the mult-add of the eigenvalues (`triangular_eigen`).
* `coupling_sets_weight`: the off-diagonal entries of a 2×2 matrix set the weight.
* `transpose_one_add_smul_mul`: `(1 + k • X)ᵀ * (1 + k • X) = 1 + k • wmadd k Xᵀ X`.
* `multAdd_inv_sub_one`, `msubM_multAdd`: the way back, off the wall `IsUnit (1 + A).det`.
* `cayley_eq`, `cayley_eq_one_add_msubM`: the Cayley transform is `1 + msubM Sᵀ S`.
* `fl_sadd`, `fl_multAdd`, `mono_multAdd`, `mono_compose`: fused-linear and monomial maps.

All of this is known mathematics. In mathlib: `Matrix.det_mul`, `Matrix.trace_add`,
`Matrix.det_of_upperTriangular`, `Matrix.mul_nonsing_inv`; the layer expansion of
`det (1 + k • A)` is the reversed characteristic polynomial `Matrix.charpolyRev`. What the file
records is that each statement is a statement about mult-add.
-/

@[expose] public section

namespace LinearAlgebraMA

open MultAdd Matrix

/-! ### Weighted mult-add in an algebra -/

section Algebra

variable {R A : Type*} [CommSemiring R] [Semiring A] [Algebra R A]

/-- Weighted mult-add in an `R`-algebra: `x + y + k • (x * y)`. For matrices this is
`A + B + k • (A * B)`. -/
def wmadd (k : R) (x y : A) : A := x + y + k • (x * y)

theorem wmadd_def (k : R) (x y : A) : wmadd k x y = x + y + k • (x * y) := rfl

/-- The weighted lift: `1 + k • wmadd k x y = (1 + k • x) * (1 + k • y)`. -/
theorem one_add_smul_wmadd (k : R) (x y : A) :
    1 + k • wmadd k x y = (1 + k • x) * (1 + k • y) := by
  simp only [wmadd, smul_add, mul_add, add_mul, one_mul, mul_one, smul_mul_smul_comm, smul_smul]
  abel

/-- At weight `0` the corner is gone: plain addition. -/
@[simp]
theorem wmadd_zero_weight (x y : A) : wmadd (0 : R) x y = x + y := by simp [wmadd]

/-- At weight `1`: plain mult-add. -/
@[simp]
theorem wmadd_one (x y : A) : wmadd (1 : R) x y = x ✳ y := by simp [wmadd, multAdd]

/-- In a commutative semiring, `wmadd` is the weighted mult-add `MultAdd.wAdd`. -/
theorem wmadd_eq_wAdd (k a b : R) : wmadd k a b = wAdd k a b := by
  simp [wmadd, wAdd, mul_assoc]

end Algebra

/-! ### Square matrices of any size -/

section Matrix

variable {n R : Type*} [Fintype n] [DecidableEq n] [CommRing R]

/-- `det (1 + ·)` turns matrix mult-add into multiplication. -/
theorem det_one_add_multAdd (A B : Matrix n n R) :
    det (1 + (A ✳ B)) = det (1 + A) * det (1 + B) := by
  rw [one_add_multAdd, det_mul]

/-- The shifted determinant `det (1 + A) - 1` turns matrix mult-add into scalar mult-add. -/
theorem shifted_det_hom (A B : Matrix n n R) :
    det (1 + (A ✳ B)) - 1 = (det (1 + A) - 1) ✳ (det (1 + B) - 1) := by
  rw [det_one_add_multAdd, multAdd_eq_mul_sub_one]; ring

/-- `det (1 + k • ·)` turns weighted matrix mult-add into multiplication. -/
theorem det_one_add_smul_wmadd (k : R) (A B : Matrix n n R) :
    det (1 + k • wmadd k A B) = det (1 + k • A) * det (1 + k • B) := by
  rw [one_add_smul_wmadd, det_mul]

/-- Rotations: `1 + k • X` is orthogonal exactly when `wmadd k Xᵀ X = 0`, because
`(1 + k • X)ᵀ * (1 + k • X) = 1 + k • wmadd k Xᵀ X`. -/
theorem transpose_one_add_smul_mul (k : R) (X : Matrix n n R) :
    (1 + k • X)ᵀ * (1 + k • X) = 1 + k • wmadd k Xᵀ X := by
  rw [one_add_smul_wmadd, transpose_add, transpose_one, transpose_smul]

/-- At `k = 0` the rotation condition drops its corner and becomes the skew condition
`Xᵀ + X = 0` (the Lie algebra). -/
theorem wmadd_zero_transpose (X : Matrix n n R) : wmadd (0 : R) Xᵀ X = Xᵀ + X :=
  wmadd_zero_weight _ _

/-- For a triangular matrix, `det (1 + A)` is the product of the lifted diagonal entries. -/
theorem det_one_add_of_isUpperTriangular [LinearOrder n] {A : Matrix n n R}
    (h : A.IsUpperTriangular) : det (1 + A) = ∏ i, (1 + A i i) := by
  rw [det_of_isUpperTriangular (blockTriangular_one.add h)]
  simp

/-- Off the wall (`det (1 + A)` a unit) every matrix has a mult-add inverse: `(1 + A)⁻¹ - 1`. -/
theorem multAdd_inv_sub_one (A : Matrix n n R) (h : IsUnit (1 + A).det) :
    A ✳ ((1 + A)⁻¹ - 1) = 0 := by
  have e : A ✳ ((1 + A)⁻¹ - 1) = (1 + A) * (1 + A)⁻¹ - 1 := by
    rw [multAdd]; noncomm_ring
  rw [e, mul_nonsing_inv _ h, sub_self]

/-- Matrix mult-add subtraction: `msubM X A = (X - A) * (1 + A)⁻¹`. -/
noncomputable def msubM (X A : Matrix n n R) : Matrix n n R := (X - A) * (1 + A)⁻¹

/-- Subtraction undoes mult-add, off the wall. -/
theorem msubM_multAdd (X A : Matrix n n R) (h : IsUnit (1 + A).det) : msubM X A ✳ A = X := by
  have e : msubM X A ✳ A = (X - A) * ((1 + A)⁻¹ * (1 + A)) + A := by
    rw [multAdd, msubM]; noncomm_ring
  rw [e, nonsing_inv_mul _ h, mul_one, sub_add_cancel]

end Matrix

/-! ### Two by two -/

section FinTwo

variable {R : Type*} [CommRing R]

/-- The dial invariant of a 2×2 matrix: `dialTwo k A = trace A + k * det A`, so that
`det (1 + k • A) = 1 + k * dialTwo k A`. -/
def dialTwo (k : R) (A : Matrix (Fin 2) (Fin 2) R) : R := trace A + k * det A

/-- `det (1 + k • A) = 1 + k * trace A + k ^ 2 * det A`. -/
theorem det_one_add_smul_fin_two (k : R) (A : Matrix (Fin 2) (Fin 2) R) :
    det (1 + k • A) = 1 + k * trace A + k ^ 2 * det A := by
  simp [det_fin_two, trace_fin_two]; ring

theorem det_one_add_smul_eq_dialTwo (k : R) (A : Matrix (Fin 2) (Fin 2) R) :
    det (1 + k • A) = 1 + k * dialTwo k A := by
  rw [det_one_add_smul_fin_two, dialTwo]; ring

/-- The dial homomorphism, for every weight `k` at once. -/
theorem dialTwo_wmadd (k : R) (A B : Matrix (Fin 2) (Fin 2) R) :
    dialTwo k (wmadd k A B) = wAdd k (dialTwo k A) (dialTwo k B) := by
  simp [dialTwo, wmadd, wAdd, det_fin_two, trace_fin_two, Matrix.mul_apply, Fin.sum_univ_two]; ring

/-- `det (1 + M) - 1 = wAdd w m₁₁ m₂₂` whenever `m₁₂ * m₂₁ = (1 - w) * m₁₁ * m₂₂`: the
off-diagonal entries set the dial. -/
theorem coupling_sets_weight (w : R) (M : Matrix (Fin 2) (Fin 2) R)
    (h : M 0 1 * M 1 0 = (1 - w) * M 0 0 * M 1 1) :
    det (1 + M) - 1 = wAdd w (M 0 0) (M 1 1) := by
  simp only [det_fin_two, wAdd, Matrix.add_apply, one_apply_eq, one_apply_ne, ne_eq, zero_ne_one,
    one_ne_zero, not_false_eq_true, zero_add]
  linear_combination -h

/-- No coupling (diagonal): mult-add, weight `1`. -/
theorem diagonal_is_multAdd (p q : R) : det (1 + !![p, 0; 0, q]) - 1 = p ✳ q := by
  simp [det_fin_two, multAdd]; ring

/-- Rank one (determinant `0`): plain addition, weight `0`. -/
theorem rank_one_is_add (M : Matrix (Fin 2) (Fin 2) R) (h : M.det = 0) :
    det (1 + M) - 1 = M.trace := by
  have e := det_one_add_smul_fin_two 1 M
  rw [one_smul, h] at e
  rw [e]; ring

/-- Order matters for matrix mult-add. -/
theorem multAdd_not_comm :
    (!![0, 1; 0, 0] : Matrix (Fin 2) (Fin 2) ℤ) ✳ !![0, 0; 1, 0] ≠
      !![0, 0; 1, 0] ✳ !![0, 1; 0, 0] := by
  intro h
  have h1 := congrFun (congrFun h 0) 0
  simp [multAdd] at h1

/-- For `X = !![x, -y; y, x]`, `wmadd k Xᵀ X` is `(2 * x + k * (x ^ 2 + y ^ 2)) • 1`. -/
theorem rotation_condition (k x y : R) :
    wmadd k !![x, -y; y, x]ᵀ !![x, -y; y, x] =
      !![2 * x + k * (x ^ 2 + y ^ 2), 0; 0, 2 * x + k * (x ^ 2 + y ^ 2)] := by
  have e : !![x, -y; y, x]ᵀ = !![x, y; -y, x] := by
    ext i j
    rw [transpose_apply]
    fin_cases i <;> fin_cases j <;> rfl
  rw [e]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [wmadd] <;> ring

/-- The 2×2 skew matrix `!![0, -s; s, 0]`. -/
def skew (s : R) : Matrix (Fin 2) (Fin 2) R := !![0, -s; s, 0]

/-- The Cayley transform `(1 + S)⁻¹ * (1 - S)` of the skew matrix `S = skew s`. -/
noncomputable def cayley (s : R) : Matrix (Fin 2) (Fin 2) R := (1 + skew s)⁻¹ * (1 - skew s)

theorem transpose_skew (s : R) : (skew s)ᵀ = skew (-s) := by
  ext i j
  rw [transpose_apply]
  fin_cases i <;> fin_cases j <;> simp [skew]

theorem det_one_add_skew (s : R) : det (1 + skew s) = 1 + s ^ 2 := by
  simp [det_fin_two, skew]; ring

end FinTwo

/-! ### Three by three -/

section FinThree

variable {R : Type*} [CommRing R]

/-- The second layer of a 3×3 matrix: the sum of its three principal 2×2 minors. -/
def e2 (A : Matrix (Fin 3) (Fin 3) R) : R :=
  (A 0 0 * A 1 1 - A 0 1 * A 1 0) + (A 0 0 * A 2 2 - A 0 2 * A 2 0)
    + (A 1 1 * A 2 2 - A 1 2 * A 2 1)

/-- The dial invariant of a 3×3 matrix: `trace A + k * e2 A + k ^ 2 * det A`, so that
`det (1 + k • A) = 1 + k * dialThree k A`. -/
def dialThree (k : R) (A : Matrix (Fin 3) (Fin 3) R) : R := trace A + k * e2 A + k ^ 2 * det A

/-- The mult-add cube of a matrix: `det (1 + k • A) = 1 + k * e₁ + k ^ 2 * e₂ + k ^ 3 * e₃`, with
`e₁` the trace, `e₂` the principal 2×2 minors and `e₃` the determinant (counts 1, 3, 3, 1). -/
theorem det_one_add_smul_fin_three (k : R) (A : Matrix (Fin 3) (Fin 3) R) :
    det (1 + k • A) = 1 + k * trace A + k ^ 2 * e2 A + k ^ 3 * det A := by
  simp [det_fin_three, trace_fin_three, e2]; ring

/-- `k = 1`: `det (1 + A) = 1 + e₁ + e₂ + e₃`. -/
theorem det_one_add_fin_three (A : Matrix (Fin 3) (Fin 3) R) :
    det (1 + A) = 1 + trace A + e2 A + det A := by
  have e := det_one_add_smul_fin_three 1 A
  rw [one_smul] at e
  rw [e]; ring

theorem det_one_add_smul_eq_dialThree (k : R) (A : Matrix (Fin 3) (Fin 3) R) :
    det (1 + k • A) = 1 + k * dialThree k A := by
  rw [det_one_add_smul_fin_three, dialThree]; ring

/-- The box picture: stretching the unit cube to `(1 + a) * (1 + b) * (1 + c)` adds 3 slabs,
3 beams and 1 corner block. -/
theorem box_cube (a b c : R) :
    det (1 + !![a, 0, 0; 0, b, 0; 0, 0, c]) =
      1 + (a + b + c) + (a * b + a * c + b * c) + a * b * c := by
  simp [det_fin_three]; ring

/-- For a triangular 3×3 matrix (eigenvalues on the diagonal) the shifted determinant is the
mult-add of the eigenvalues, as the trace is their sum. See `det_one_add_of_isUpperTriangular` for
any size. -/
theorem triangular_eigen (A : Matrix (Fin 3) (Fin 3) R) (h10 : A 1 0 = 0) (h20 : A 2 0 = 0)
    (h21 : A 2 1 = 0) : det (1 + A) - 1 = multAddList [A 0 0, A 1 1, A 2 2] := by
  simp [det_fin_three, multAdd, h10, h20, h21]; ring

/-- The 3×3 dial homomorphism for every weight `k`: `k = 0` is trace additivity, `k = 1` is
`det (1 + A)` multiplicative. -/
theorem dialThree_wmadd (k : R) (A B : Matrix (Fin 3) (Fin 3) R) :
    dialThree k (wmadd k A B) = wAdd k (dialThree k A) (dialThree k B) := by
  simp [dialThree, wmadd, wAdd, e2, det_fin_three, trace_fin_three, Matrix.mul_apply,
    Fin.sum_univ_three]
  ring

end FinThree

/-! ### Fields: Cayley, the rotation circle, the KdV pair -/

section Field

variable {K : Type*} [Field K]

/-- The Cayley transform of a skew matrix is the rational rotation `!![p, q; -q, p]` with
`p = (1 - s ^ 2) / (1 + s ^ 2)` and `q = 2 * s / (1 + s ^ 2)`. -/
theorem cayley_eq (s : K) (h : 1 + s ^ 2 ≠ 0) :
    cayley s = !![(1 - s ^ 2) / (1 + s ^ 2), 2 * s / (1 + s ^ 2);
                  -(2 * s) / (1 + s ^ 2), (1 - s ^ 2) / (1 + s ^ 2)] := by
  rw [cayley, inv_def, det_one_add_skew, Ring.inverse_eq_inv', adjugate_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [skew, Matrix.mul_apply, Fin.sum_univ_two] <;>
    field_simp <;> ring

/-- The Cayley transform is mult-add subtraction of `S` from its transpose: `1 + msubM Sᵀ S`. -/
theorem cayley_eq_one_add_msubM (s : K) (h : 1 + s ^ 2 ≠ 0) :
    cayley s = 1 + msubM (skew s)ᵀ (skew s) := by
  rw [cayley_eq s h, msubM, transpose_skew, inv_def, det_one_add_skew, Ring.inverse_eq_inv',
    adjugate_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [skew] <;> field_simp <;> ring

/-- The Cayley point lies on the unit circle (the matrix is a rotation). -/
theorem cayley_on_circle (s : K) (h : 1 + s ^ 2 ≠ 0) :
    ((1 - s ^ 2) / (1 + s ^ 2)) ^ 2 + (2 * s / (1 + s ^ 2)) ^ 2 = 1 := by
  field_simp; ring

/-- The same point, shifted by `-1`, satisfies the `k = 1` rotation condition
`2 * x + x ^ 2 + y ^ 2 = 0`: it lies on the circle centred at the wall `-1`. -/
theorem cayley_on_wall_circle (s : K) (h : 1 + s ^ 2 ≠ 0) :
    2 * ((1 - s ^ 2) / (1 + s ^ 2) - 1) + ((1 - s ^ 2) / (1 + s ^ 2) - 1) ^ 2
      + (2 * s / (1 + s ^ 2)) ^ 2 = 0 := by
  field_simp; ring

/-- The weighted rotation condition is a circle through `0` centred at the weighted wall
`-1 / k`. -/
theorem rotation_circle (k x y : K) (hk : k ≠ 0) :
    2 * x + k * (x ^ 2 + y ^ 2) = 0 ↔ (x + 1 / k) ^ 2 + y ^ 2 = (1 / k) ^ 2 := by
  have e : (x + 1 / k) ^ 2 + y ^ 2 - (1 / k) ^ 2 = (2 * x + k * (x ^ 2 + y ^ 2)) / k := by
    field_simp; ring
  rw [← sub_eq_zero (a := (x + 1 / k) ^ 2 + y ^ 2), e, div_eq_zero_iff, or_iff_left hk]

/-- The KdV two-soliton: with the Cauchy-type coupling
`m ^ 2 = 4 * k₁ * k₂ / (k₁ + k₂) ^ 2 * E₁ * E₂`, `det (1 + M)` is Hirota's
`τ = 1 + E₁ + E₂ + ((k₁ - k₂) / (k₁ + k₂)) ^ 2 * E₁ * E₂`. -/
theorem kdv_two_soliton (k₁ k₂ E₁ E₂ m : K) (hs : k₁ + k₂ ≠ 0)
    (hm : m ^ 2 = 4 * k₁ * k₂ / (k₁ + k₂) ^ 2 * E₁ * E₂) :
    det (1 + !![E₁, m; m, E₂]) = 1 + E₁ + E₂ + ((k₁ - k₂) / (k₁ + k₂)) ^ 2 * E₁ * E₂ := by
  have e : ((k₁ - k₂) / (k₁ + k₂)) ^ 2 = 1 - 4 * k₁ * k₂ / (k₁ + k₂) ^ 2 := by
    field_simp; ring
  simp only [det_fin_two, Matrix.add_apply, one_apply_eq, one_apply_ne, ne_eq, zero_ne_one,
    one_ne_zero, not_false_eq_true, zero_add, of_apply, cons_val', cons_val_zero, cons_val_one,
    cons_val_fin_one]
  rw [e]
  linear_combination -hm

end Field

/-! ### Fused-linear maps: linear maps about the wall corner -/

section Fused

variable {n R : Type*} [Fintype n] [Ring R]

/-- The fused addition (the rung below mult-add): `u + v + 1`, with zero `-1`. -/
def sadd {α : Type*} [Add α] [One α] (u v : α) : α := u + v + 1

/-- A fused-linear map: shift to the box corner, apply `M`, shift back. -/
def fl (M : Matrix n n R) (v : n → R) : n → R := M *ᵥ (1 + v) - 1

/-- A fused-linear map fixes the wall corner. -/
theorem fl_neg_one (M : Matrix n n R) : fl M (-1) = -1 := by simp [fl]

/-- A fused-linear map respects the fused addition. -/
theorem fl_sadd (M : Matrix n n R) (u v : n → R) : fl M (sadd u v) = sadd (fl M u) (fl M v) := by
  have e : 1 + (u + v + 1) = (1 + u) + (1 + v) := by abel
  rw [sadd, sadd, fl, fl, fl, e, mulVec_add]
  abel

/-- A fused-linear map respects mult-add scaling: `fl M (c ✳ v) = c ✳ fl M v`, entry by entry. -/
theorem fl_multAdd {R : Type*} [CommRing R] (M : Matrix n n R) (c : R) (v : n → R) :
    fl M (fun i ↦ c ✳ v i) = fun i ↦ c ✳ fl M v i := by
  have e : (1 + fun i ↦ c ✳ v i) = (1 + c) • (1 + v) := by
    ext i; simp [one_add_multAdd]
  rw [fl, e, mulVec_smul]
  ext i
  simp only [fl, Pi.sub_apply, Pi.smul_apply, Pi.one_apply, smul_eq_mul, multAdd]
  ring

/-- Mult-adding `a` is zooming by `1 + a` about the wall `-1`. -/
theorem multAdd_is_zoom {R : Type*} [Ring R] (a x : R) :
    a ✳ x - (-1) = (1 + a) * (x - (-1)) := by
  rw [sub_neg_eq_add, sub_neg_eq_add, multAdd_add_one, add_comm a, add_comm x]

end Fused

/-! ### Rate space: monomial maps and exponent matrices -/

section Monomial

/-- A monomial map, first coordinate, in mult-add coordinates:
`(1 + x₁) ^ p * (1 + x₂) ^ q - 1`. -/
def mono {R : Type*} [CommRing R] (p q : ℕ) (x₁ x₂ : R) : R := (1 + x₁) ^ p * (1 + x₂) ^ q - 1

/-- Monomial maps are the linear maps of rate space: they turn mult-add into mult-add. -/
theorem mono_multAdd {R : Type*} [CommRing R] (p q : ℕ) (x₁ x₂ y₁ y₂ : R) :
    mono p q (x₁ ✳ y₁) (x₂ ✳ y₂) = mono p q x₁ x₂ ✳ mono p q y₁ y₂ := by
  rw [mono, one_add_multAdd, one_add_multAdd, mul_pow, mul_pow, multAdd_eq_mul_sub_one, mono, mono]
  ring

/-- Composing monomial maps multiplies exponent matrices (first row shown; the second row is the
same statement with the other row of exponents). -/
theorem mono_compose {M : Type*} [CommMonoid M] (u₁ u₂ : M) (a b c d e f : ℕ) :
    (u₁ ^ a * u₂ ^ b) ^ c * (u₁ ^ d * u₂ ^ e) ^ f =
      u₁ ^ (a * c + d * f) * u₂ ^ (b * c + e * f) := by
  rw [mul_pow, mul_pow, ← pow_mul, ← pow_mul, ← pow_mul, ← pow_mul, pow_add, pow_add,
    mul_mul_mul_comm]

end Monomial

end LinearAlgebraMA
