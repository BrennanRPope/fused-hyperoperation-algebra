/-
  Fused Hyperoperation Algebra: linear algebra through the circle operation (called mult-add before 1.0).
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base.

  LIBRARY NOTE: copied from linear-algebra/lean/LinearAlgebra.lean (Nov 25). It
  compiled against the library on its first build with no changes beyond this
  header note: no sorry, standard axioms only. (Its identities were also checked
  in linear-algebra/checks.py.)

  What is here (2×2 and 3×3 matrices written out as records, so `grind` does the algebra):
    1. The matrix circle operation A ◯ B = A + B + AB: the lift I + (A ◯ B) = (I + A)(I + B), the weighted
       lift with k, associativity, identity 0, and a 2×2 example where order matters.
    2. Determinants: det(I + kA) = 1 + k·e₁ + k²·e₂ + k³·e₃ (principal minors, counts 1,3,3,1),
       det(AB) = det A · det B, the box picture det(I + diag(a,b,c)), and for a triangular
       matrix det(I + A) − 1 = a₁₁ ◯ a₂₂ ◯ a₃₃ (the shifted determinant is the circle operation of the
       eigenvalues).
    3. The dial homomorphism D_k(A ◯[k] B) = D_k(A) ◯[k] D_k(B) with D_k = e₁ + k·e₂ + k²·e₃:
       k = 0 is "trace is additive", k = 1 is "det(I + A) is multiplicative".
    4. 2×2 coupling sets the dial: det(I + M) − 1 = m₁₁ ◯[w] m₂₂ when m₁₂m₂₁ = (1 − w)m₁₁m₂₂;
       diagonal gives w = 1 (the circle operation), rank one gives w = 0 (plain addition). The KdV
       two-soliton τ is the case w = ((k₁ − k₂)/(k₁ + k₂))².
    5. Rotations: (I + kX)ᵀ(I + kX) = I + k(Xᵀ ◯[k] X); for X = [[x, −y],[y, x]] the condition is
       2x + k(x² + y²) = 0, the circle through 0 centred at the weighted wall −1/k; at k = 0 it
       is the skew condition Xᵀ + X = 0.
    6. Inverses on a field: the circle inverse (I + A)⁻¹ − I, matrix circle-operation subtraction
       X ⊟ A = (X − A)(I + A)⁻¹ undoes ◯, and the Cayley transform of a 2×2 skew matrix is
       I + (Sᵀ ⊟ S), the rational point ((1 − s²)/(1 + s²), 2s/(1 + s²)) of the unit circle.
    7. Fused-linear maps (linear maps about the wall corner (−1, −1)): they fix the corner,
       respect the fused addition u + v + 1 and the circle-operation scaling c ◯ v.
    8. Rate space: monomial maps u ↦ (u₁^m₁₁ u₂^m₁₂, u₁^m₂₁ u₂^m₂₂) turn ◯ into ◯, and
       composing two of them multiplies their exponent matrices.
-/
import FusedAlgebra.Base

namespace LinearAlgebraCirc
open CircleOp

set_option linter.unusedSectionVars false

/-! ## 0. Matrices as records -/

/-- A 2×2 matrix [[a, b], [c, d]]. -/
structure M2 (R : Type) where
  a : R
  b : R
  c : R
  d : R

/-- A 3×3 matrix with entries aᵢⱼ. -/
structure M3 (R : Type) where
  a11 : R
  a12 : R
  a13 : R
  a21 : R
  a22 : R
  a23 : R
  a31 : R
  a32 : R
  a33 : R

theorem M2.ext' {R : Type} {x y : M2 R} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c)
    (h4 : x.d = y.d) : x = y := by
  cases x; cases y
  dsimp only at h1 h2 h3 h4
  subst h1; subst h2; subst h3; subst h4; rfl

theorem M3.ext' {R : Type} {x y : M3 R} (h11 : x.a11 = y.a11) (h12 : x.a12 = y.a12)
    (h13 : x.a13 = y.a13) (h21 : x.a21 = y.a21) (h22 : x.a22 = y.a22) (h23 : x.a23 = y.a23)
    (h31 : x.a31 = y.a31) (h32 : x.a32 = y.a32) (h33 : x.a33 = y.a33) : x = y := by
  cases x; cases y
  dsimp only at h11 h12 h13 h21 h22 h23 h31 h32 h33
  subst h11; subst h12; subst h13; subst h21; subst h22; subst h23
  subst h31; subst h32; subst h33; rfl

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

namespace M2

def one : M2 R := ⟨1, 0, 0, 1⟩
def zero : M2 R := ⟨0, 0, 0, 0⟩
def add (x y : M2 R) : M2 R := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def sub (x y : M2 R) : M2 R := ⟨x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d⟩
def mul (x y : M2 R) : M2 R :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def smul (k : R) (x : M2 R) : M2 R := ⟨k * x.a, k * x.b, k * x.c, k * x.d⟩
def transpose (x : M2 R) : M2 R := ⟨x.a, x.c, x.b, x.d⟩
def tr (x : M2 R) : R := x.a + x.d
def det (x : M2 R) : R := x.a * x.d - x.b * x.c

/-- The matrix circle operation: A ◯ B = A + B + AB. -/
def circ (x y : M2 R) : M2 R := add (add x y) (mul x y)
/-- The matrix circle operation with dial k: A ◯[k] B = A + B + k·AB. -/
def circK (k : R) (x y : M2 R) : M2 R := add (add x y) (smul k (mul x y))
/-- The dial invariant: D_k(A) = tr A + k·det A, so det(I + kA) = 1 + k·D_k(A). -/
def Dk (k : R) (x : M2 R) : R := tr x + k * det x

/-! ## 1. The lift: the matrix circle operation is matrix multiplication through A ↦ I + A -/

theorem lift (x y : M2 R) : add one (circ x y) = mul (add one x) (add one y) := by
  apply M2.ext' <;> simp only [circ, add, mul, one] <;> grind

theorem wlift (k : R) (x y : M2 R) :
    add one (smul k (circK k x y)) = mul (add one (smul k x)) (add one (smul k y)) := by
  apply M2.ext' <;> simp only [circK, add, mul, smul, one] <;> grind

theorem circ_assoc (x y z : M2 R) : circ (circ x y) z = circ x (circ y z) := by
  apply M2.ext' <;> simp only [circ, add, mul] <;> grind

theorem circ_zero (x : M2 R) : circ x zero = x := by
  apply M2.ext' <;> simp only [circ, add, mul, zero] <;> grind

theorem zero_circ (x : M2 R) : circ zero x = x := by
  apply M2.ext' <;> simp only [circ, add, mul, zero] <;> grind

/-- At weight 0 the corner is gone: plain matrix addition. -/
theorem circK_zero_weight (x y : M2 R) : circK 0 x y = add x y := by
  apply M2.ext' <;> simp only [circK, add, mul, smul] <;> grind

/-! ## 2. Determinant and trace -/

/-- det(I + kA) = 1 + k·tr A + k²·det A. -/
theorem det_one_add_smul (k : R) (x : M2 R) :
    det (add one (smul k x)) = 1 + k * tr x + k * k * det x := by
  simp only [det, tr, add, smul, one]; grind

theorem det_mul (x y : M2 R) : det (mul x y) = det x * det y := by
  simp only [det, mul]; grind

theorem tr_add (x y : M2 R) : tr (add x y) = tr x + tr y := by
  simp only [tr, add]; grind

/-- The dial homomorphism, for every weight k at once (any commutative ring). -/
theorem Dk_hom (k : R) (x y : M2 R) : Dk k (circK k x y) = CircleOp.circK k (Dk k x) (Dk k y) := by
  simp only [Dk, tr, det, circK, add, smul, mul, CircleOp.circK]; grind

/-- k = 1: the shifted determinant det(I + A) − 1 turns matrix ◯ into scalar ◯. -/
theorem shifted_det_hom (x y : M2 R) :
    det (add one (circ x y)) - 1 = CircleOp.circ (det (add one x) - 1) (det (add one y) - 1) := by
  rw [lift, det_mul]; unfold CircleOp.circ; grind

/-! ## 3. The off-diagonal entries set the dial -/

/-- det(I + M) − 1 = m₁₁ ◯[w] m₂₂ whenever m₁₂·m₂₁ = (1 − w)·m₁₁·m₂₂. -/
theorem coupling_sets_weight (w : R) (m : M2 R) (h : m.b * m.c = (1 - w) * m.a * m.d) :
    det (add one m) - 1 = CircleOp.circK w m.a m.d := by
  simp only [det, add, one, CircleOp.circK]; grind

/-- No coupling (diagonal): the circle operation, weight 1. -/
theorem diagonal_is_circ (p q : R) : det (add one ⟨p, 0, 0, q⟩) - 1 = CircleOp.circ p q := by
  simp only [det, add, one, CircleOp.circ]; grind

/-- Rank one (m₁₂m₂₁ = m₁₁m₂₂, determinant 0): plain addition, weight 0. -/
theorem rank_one_is_add (m : M2 R) (h : m.b * m.c = m.a * m.d) :
    det (add one m) - 1 = m.a + m.d := by
  simp only [det, add, one]; grind

/-- Two by two, a matrix never commutes for free: order matters for matrix ◯. -/
theorem circ_not_comm :
    circ (⟨0, 1, 0, 0⟩ : M2 Int) ⟨0, 0, 1, 0⟩ ≠ circ ⟨0, 0, 1, 0⟩ ⟨0, 1, 0, 0⟩ := by
  intro h
  have h1 := congrArg M2.a h
  revert h1
  first
    | decide
    | (simp only [circ, add, mul]; grind)

/-! ## 4. Rotations: the orthogonal condition is Xᵀ ◯[k] X = 0 -/

theorem orth_lift (k : R) (x : M2 R) :
    mul (transpose (add one (smul k x))) (add one (smul k x))
      = add one (smul k (circK k (transpose x) x)) := by
  apply M2.ext' <;> simp only [transpose, circK, add, mul, smul, one] <;> grind

/-- For X = [[x, −y], [y, x]], Xᵀ ◯[k] X is (2x + k(x² + y²))·I. -/
theorem rotation_condition (k x y : R) :
    circK k (transpose ⟨x, -y, y, x⟩) ⟨x, -y, y, x⟩
      = ⟨2 * x + k * (x * x + y * y), 0, 0, 2 * x + k * (x * x + y * y)⟩ := by
  apply M2.ext' <;> simp only [transpose, circK, add, mul, smul] <;> grind

/-- k = 0 drops the corner: the condition becomes Xᵀ + X = 0 (skew, the Lie algebra). -/
theorem rotation_condition_flat (x : M2 R) : circK 0 (transpose x) x = add (transpose x) x :=
  circK_zero_weight _ _

end M2

/-! ## 5. Three by three -/

namespace M3

def one : M3 R := ⟨1, 0, 0, 0, 1, 0, 0, 0, 1⟩
def zero : M3 R := ⟨0, 0, 0, 0, 0, 0, 0, 0, 0⟩
def add (x y : M3 R) : M3 R :=
  ⟨x.a11 + y.a11, x.a12 + y.a12, x.a13 + y.a13,
   x.a21 + y.a21, x.a22 + y.a22, x.a23 + y.a23,
   x.a31 + y.a31, x.a32 + y.a32, x.a33 + y.a33⟩
def smul (k : R) (x : M3 R) : M3 R :=
  ⟨k * x.a11, k * x.a12, k * x.a13, k * x.a21, k * x.a22, k * x.a23,
   k * x.a31, k * x.a32, k * x.a33⟩
def mul (x y : M3 R) : M3 R :=
  ⟨x.a11 * y.a11 + x.a12 * y.a21 + x.a13 * y.a31,
   x.a11 * y.a12 + x.a12 * y.a22 + x.a13 * y.a32,
   x.a11 * y.a13 + x.a12 * y.a23 + x.a13 * y.a33,
   x.a21 * y.a11 + x.a22 * y.a21 + x.a23 * y.a31,
   x.a21 * y.a12 + x.a22 * y.a22 + x.a23 * y.a32,
   x.a21 * y.a13 + x.a22 * y.a23 + x.a23 * y.a33,
   x.a31 * y.a11 + x.a32 * y.a21 + x.a33 * y.a31,
   x.a31 * y.a12 + x.a32 * y.a22 + x.a33 * y.a32,
   x.a31 * y.a13 + x.a32 * y.a23 + x.a33 * y.a33⟩
def transpose (x : M3 R) : M3 R :=
  ⟨x.a11, x.a21, x.a31, x.a12, x.a22, x.a32, x.a13, x.a23, x.a33⟩
/-- e₁: the trace (the 3 one-by-one principal minors). -/
def tr (x : M3 R) : R := x.a11 + x.a22 + x.a33
/-- e₂: the 3 two-by-two principal minors. -/
def e2 (x : M3 R) : R :=
  (x.a11 * x.a22 - x.a12 * x.a21) + (x.a11 * x.a33 - x.a13 * x.a31)
    + (x.a22 * x.a33 - x.a23 * x.a32)
/-- e₃: the determinant (the 1 three-by-three minor). -/
def det (x : M3 R) : R :=
  x.a11 * (x.a22 * x.a33 - x.a23 * x.a32) - x.a12 * (x.a21 * x.a33 - x.a23 * x.a31)
    + x.a13 * (x.a21 * x.a32 - x.a22 * x.a31)

def circ (x y : M3 R) : M3 R := add (add x y) (mul x y)
def circK (k : R) (x y : M3 R) : M3 R := add (add x y) (smul k (mul x y))
/-- D_k(A) = e₁ + k·e₂ + k²·e₃, so det(I + kA) = 1 + k·D_k(A). -/
def Dk (k : R) (x : M3 R) : R := tr x + k * e2 x + k * k * det x

theorem lift (x y : M3 R) : add one (circ x y) = mul (add one x) (add one y) := by
  apply M3.ext' <;> simp only [circ, add, mul, one] <;> grind

theorem wlift (k : R) (x y : M3 R) :
    add one (smul k (circK k x y)) = mul (add one (smul k x)) (add one (smul k y)) := by
  apply M3.ext' <;> simp only [circK, add, mul, smul, one] <;> grind

theorem circ_assoc (x y z : M3 R) : circ (circ x y) z = circ x (circ y z) := by
  apply M3.ext' <;> simp only [circ, add, mul] <;> grind

theorem circ_zero (x : M3 R) : circ x zero = x := by
  apply M3.ext' <;> simp only [circ, add, mul, zero] <;> grind

/-- The circle-operation cube of a matrix: det(I + kA) = 1 + k·e₁ + k²·e₂ + k³·e₃. -/
theorem det_one_add_smul (k : R) (x : M3 R) :
    det (add one (smul k x)) = 1 + k * tr x + k * k * e2 x + k * k * k * det x := by
  simp only [det, tr, e2, add, smul, one]; grind

/-- k = 1: det(I + A) = 1 + e₁ + e₂ + e₃. -/
theorem det_one_add (x : M3 R) : det (add one x) = 1 + tr x + e2 x + det x := by
  simp only [det, tr, e2, add, one]; grind

set_option maxHeartbeats 2000000 in
theorem det_mul (x y : M3 R) : det (mul x y) = det x * det y := by
  simp only [det, mul]; grind

/-- The box picture: stretching the unit cube to (1 + a)(1 + b)(1 + c) adds 3 slabs,
    3 beams and 1 corner block: 1, 3, 3, 1. -/
theorem box_cube (a b c : R) :
    det (add one ⟨a, 0, 0, 0, b, 0, 0, 0, c⟩)
      = 1 + (a + b + c) + (a * b + a * c + b * c) + a * b * c := by
  simp only [det, add, one]; grind

/-- For a triangular matrix (eigenvalues on the diagonal) the shifted determinant is the
    circle operation of the eigenvalues, as the trace is their sum. -/
theorem triangular_eigen (x : M3 R) (h21 : x.a21 = 0) (h31 : x.a31 = 0) (h32 : x.a32 = 0) :
    det (add one x) - 1 = foldCirc [x.a11, x.a22, x.a33] := by
  simp only [det, add, one, foldCirc, CircleOp.circ]; rw [h21, h31, h32]; grind

theorem tr_add (x y : M3 R) : tr (add x y) = tr x + tr y := by
  simp only [tr, add]; grind

theorem orth_lift (k : R) (x : M3 R) :
    mul (transpose (add one (smul k x))) (add one (smul k x))
      = add one (smul k (circK k (transpose x) x)) := by
  apply M3.ext' <;> simp only [transpose, circK, add, mul, smul, one] <;> grind

end M3

end Ring

/-! ## 6. Fields: the dial homomorphism in 3×3, inverses, Cayley, the KdV pair -/

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- The 3×3 dial homomorphism for every weight k: k = 0 is trace additivity, k ≠ 0 is
    det(I + kA) multiplicative. -/
theorem M3.Dk_hom (k : K) (x y : M3 K) :
    M3.Dk k (M3.circK k x y) = CircleOp.circK k (M3.Dk k x) (M3.Dk k y) := by
  by_cases hk : k = 0
  · subst hk
    simp only [M3.Dk, M3.tr, M3.circK, M3.add, M3.smul, M3.mul, CircleOp.circK]; grind
  · have h1 := M3.det_one_add_smul k (M3.circK k x y)
    rw [M3.wlift, M3.det_mul, M3.det_one_add_smul, M3.det_one_add_smul] at h1
    have h6 : k * M3.Dk k (M3.circK k x y) = k * CircleOp.circK k (M3.Dk k x) (M3.Dk k y) := by
      unfold M3.Dk CircleOp.circK; grind
    grind

namespace M2

/-- The 2×2 inverse by the adjugate. -/
def inv (x : M2 K) : M2 K := ⟨x.d / det x, -x.b / det x, -x.c / det x, x.a / det x⟩

theorem mul_inv (w : M2 K) (h : det w ≠ 0) : mul w (inv w) = one := by
  have hd : det w = w.a * w.d - w.b * w.c := rfl
  apply M2.ext' <;> simp only [mul, inv, one] <;> grind

theorem add_one_eq_one {w : M2 K} (h : add one w = one) : w = zero := by
  have h1 := congrArg M2.a h
  have h2 := congrArg M2.b h
  have h3 := congrArg M2.c h
  have h4 := congrArg M2.d h
  simp only [add, one] at h1 h2 h3 h4
  apply M2.ext' <;> simp only [zero] <;> grind

/-- Off the wall (det(I + A) ≠ 0) every matrix has a circle inverse: (I + A)⁻¹ − I. -/
theorem circ_inverse (x : M2 K) (h : det (add one x) ≠ 0) :
    circ x (sub (inv (add one x)) one) = zero := by
  have hl := lift x (sub (inv (add one x)) one)
  have hz : add one (sub (inv (add one x)) one) = inv (add one x) := by
    apply M2.ext' <;> simp only [add, sub, one] <;> grind
  rw [hz, mul_inv (add one x) h] at hl
  exact add_one_eq_one hl

/-- Matrix circle-operation subtraction X ⊟ A = (X − A)(I + A)⁻¹. -/
def msubM (x a : M2 K) : M2 K := mul (sub x a) (inv (add one a))

theorem msub_undo_W (x w : M2 K) (h : det w ≠ 0) :
    circ (mul (sub x (sub w one)) (inv w)) (sub w one) = x := by
  have hd : det w = w.a * w.d - w.b * w.c := rfl
  apply M2.ext' <;> simp only [circ, add, sub, mul, inv, one] <;> grind

/-- ⊟ undoes ◯: (X ⊟ A) ◯ A = X, off the wall. -/
theorem msub_undo (x a : M2 K) (h : det (add one a) ≠ 0) : circ (msubM x a) a = x := by
  have hA : sub (add one a) one = a := by
    apply M2.ext' <;> simp only [sub, add, one] <;> grind
  have := msub_undo_W x (add one a) h
  rw [hA] at this
  exact this

/-- The 2×2 skew matrix [[0, −s], [s, 0]]. -/
def skew (s : K) : M2 K := ⟨0, -s, s, 0⟩

/-- The Cayley transform (I + S)⁻¹(I − S). -/
def cayley (s : K) : M2 K := mul (inv (add one (skew s))) (sub one (skew s))

theorem det_one_add_skew (s : K) : det (add one (skew s)) = 1 + s * s := by
  simp only [det, add, one, skew]; grind

/-- The Cayley transform of a skew matrix is the rational rotation
    [[p, q], [−q, p]] with p = (1 − s²)/(1 + s²), q = 2s/(1 + s²). -/
theorem cayley_eq (s : K) (h : 1 + s * s ≠ 0) :
    cayley s = ⟨(1 - s * s) / (1 + s * s), 2 * s / (1 + s * s),
                -(2 * s) / (1 + s * s), (1 - s * s) / (1 + s * s)⟩ := by
  apply M2.ext' <;>
  · simp only [cayley, mul, inv]
    rw [det_one_add_skew]
    simp only [add, sub, one, skew]
    grind

/-- The Cayley transform is circle-operation subtraction of S from its transpose: I + (Sᵀ ⊟ S). -/
theorem cayley_is_msub (s : K) (h : 1 + s * s ≠ 0) :
    cayley s = add one (msubM (transpose (skew s)) (skew s)) := by
  rw [cayley_eq s h]
  apply M2.ext' <;>
  · simp only [msubM, mul, inv]
    rw [det_one_add_skew]
    simp only [add, sub, one, skew, transpose]
    grind

/-- The Cayley point lies on the unit circle (the matrix is a rotation). -/
theorem cayley_on_circle (s : K) (h : 1 + s * s ≠ 0) :
    ((1 - s * s) / (1 + s * s)) * ((1 - s * s) / (1 + s * s))
      + (2 * s / (1 + s * s)) * (2 * s / (1 + s * s)) = 1 := by
  grind

/-- The same point, shifted by −1, satisfies the k = 1 rotation condition 2x + x² + y² = 0:
    it lies on the circle centred at the wall −1. -/
theorem cayley_on_wall_circle (s : K) (h : 1 + s * s ≠ 0) :
    2 * ((1 - s * s) / (1 + s * s) - 1)
      + ((1 - s * s) / (1 + s * s) - 1) * ((1 - s * s) / (1 + s * s) - 1)
      + (2 * s / (1 + s * s)) * (2 * s / (1 + s * s)) = 0 := by
  grind

end M2

/-- The weighted rotation condition is a circle through 0 centred at the weighted wall −1/k. -/
theorem rotation_circle (k x y : K) (hk : k ≠ 0) :
    2 * x + k * (x * x + y * y) = 0 ↔ (x + 1 / k) * (x + 1 / k) + y * y = (1 / k) * (1 / k) := by
  constructor
  · intro h; grind
  · intro h; grind

/-- The KdV two-soliton: with the Cauchy-type coupling m² = 4k₁k₂/(k₁ + k₂)²·E₁E₂,
    det(I + M) is Hirota's τ = 1 + E₁ + E₂ + ((k₁ − k₂)/(k₁ + k₂))²·E₁E₂. -/
theorem kdv_two_soliton (k1 k2 E1 E2 m : K) (hs : k1 + k2 ≠ 0)
    (hm : m * m = 4 * k1 * k2 / ((k1 + k2) * (k1 + k2)) * E1 * E2) :
    M2.det (M2.add M2.one ⟨E1, m, m, E2⟩)
      = 1 + E1 + E2 + ((k1 - k2) / (k1 + k2)) * ((k1 - k2) / (k1 + k2)) * E1 * E2 := by
  simp only [M2.det, M2.add, M2.one]; grind

end Field

/-! ## 7. Fused-linear maps: linear maps about the wall corner (−1, −1) -/

section Fused
variable {R : Type} [Lean.Grind.CommRing R]

/-- The fused addition (the rung below the circle operation): u + v + 1, zero −1.
  Lineage: one of five earlier names for the rung-1 "or"; the unified, standalone
    version is `RungOne.sadd` (Foundations.RungOneSemiring), and
    Foundations.RungOneLineage proves the two agree. Standard name: shifted add, the
    "or" of the ring (R, a+b+1, a+b+ab) obtained by transport along x ↦ x − 1. -/
def sadd (u v : R) : R := u + v + 1

/-- A fused-linear map: shift to the box corner, apply M, shift back. -/
def fl1 (m : M2 R) (v1 v2 : R) : R := -1 + (m.a * (1 + v1) + m.b * (1 + v2))
def fl2 (m : M2 R) (v1 v2 : R) : R := -1 + (m.c * (1 + v1) + m.d * (1 + v2))

theorem fl_fixes_wall (m : M2 R) : fl1 m (-1) (-1) = -1 ∧ fl2 m (-1) (-1) = -1 := by
  unfold fl1 fl2; constructor <;> grind

theorem fl_sadd (m : M2 R) (u1 u2 v1 v2 : R) :
    fl1 m (sadd u1 v1) (sadd u2 v2) = sadd (fl1 m u1 u2) (fl1 m v1 v2) ∧
    fl2 m (sadd u1 v1) (sadd u2 v2) = sadd (fl2 m u1 u2) (fl2 m v1 v2) := by
  unfold fl1 fl2 sadd; constructor <;> grind

theorem fl_scale (m : M2 R) (c v1 v2 : R) :
    fl1 m (circ c v1) (circ c v2) = circ c (fl1 m v1 v2) ∧
    fl2 m (circ c v1) (circ c v2) = circ c (fl2 m v1 v2) := by
  unfold fl1 fl2 circ; constructor <;> grind

/-- Combining with a by the circle operation is zooming by 1 + a about the wall −1. -/
theorem circ_is_zoom (a x : R) : circ a x - (-1) = (1 + a) * (x - (-1)) := by
  unfold circ; grind

end Fused

/-! ## 8. Rate space: monomial maps and exponent matrices -/

section Monomial
variable {R : Type} [Lean.Grind.CommRing R]

theorem pow_mul'' (a : R) (m n : Nat) : (a ^ m) ^ n = a ^ (m * n) := by
  induction n with
  | zero => rw [Nat.mul_zero, Lean.Grind.Semiring.pow_zero, Lean.Grind.Semiring.pow_zero]
  | succ n ih => rw [Lean.Grind.Semiring.pow_succ, ih, Nat.mul_succ, pow_add']

/-- A monomial map, first coordinate, in circle-operation coordinates: (1 + x₁)^p (1 + x₂)^q − 1. -/
def mono (p q : Nat) (x1 x2 : R) : R := (1 + x1) ^ p * (1 + x2) ^ q - 1

/-- Monomial maps are the linear maps of rate space: they turn ◯ into ◯. -/
theorem mono_hom (p q : Nat) (x1 x2 y1 y2 : R) :
    mono p q (circ x1 y1) (circ x2 y2) = circ (mono p q x1 x2) (mono p q y1 y2) := by
  unfold mono
  rw [circ_lift, circ_lift, mul_pow', mul_pow']
  generalize (1 + x1) ^ p = P1
  generalize (1 + y1) ^ p = Q1
  generalize (1 + x2) ^ q = P2
  generalize (1 + y2) ^ q = Q2
  unfold circ; grind

/-- Composing monomial maps multiplies exponent matrices (first row shown; the second row is
    the same statement with the other row of exponents). -/
theorem mono_compose (u1 u2 : R) (a b c d e f : Nat) :
    (u1 ^ a * u2 ^ b) ^ c * (u1 ^ d * u2 ^ e) ^ f = u1 ^ (a * c + d * f) * u2 ^ (b * c + e * f) := by
  rw [mul_pow', mul_pow', pow_mul'', pow_mul'', pow_mul'', pow_mul'', pow_add', pow_add']
  generalize u1 ^ (a * c) = A1
  generalize u1 ^ (d * f) = A2
  generalize u2 ^ (b * c) = B1
  generalize u2 ^ (e * f) = B2
  grind

end Monomial

end LinearAlgebraCirc
