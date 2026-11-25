/-
  Fused Hyperoperation Algebra: turning the dial in linear algebra.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports Base and Topics.LinearAlgebra
  (reuses its 2×2 and 3×3 matrix records).

  LIBRARY NOTE: copied from linear-algebra/lean/DialSweep.lean (Nov 25). It compiled
  against the library on its first build with no changes beyond this header note: no
  sorry, standard axioms only. Checked for twins: none states exactly what a library
  theorem states. orth_neg_one is LinearAlgebraCirc's M3.orth_lift at k = −1, kept as a
  named instance. (Its identities were also checked in linear-algebra/checks_dial.py.)

  What is here:
    1. Every nonzero real weight is one law in other units: k(A ◯[k] B) = (kA) ◯ (kB).
    2. det(I + kA) is the reversed characteristic polynomial: det(I + kA) = (−k)ⁿ p_A(−1/k),
       so the matrix meets the wall exactly when −1/k is an eigenvalue. k = −1 is p_A(1).
    3. Weighted rotations and Cayley: I − X orthogonal at k = −1; the weighted Cayley transform
       is the Cayley transform of k·S.
    4. A matrix on the dial. For U = [[0, −q], [1, k]] (so U² = kU − qI):
       (I + aU)(I + bU) = (1 − q·ab)I + (a + b + k·ab)U, det(I + aU) = 1 + k·a + q·a².
       This is the landscape lens law (a + b + k·ab)/(1 − q·ab), its walls are −1/eigenvalue,
       its countries are the sign of k² − 4q = (λ₁ − λ₂)², and the "or" line k + q = −1 is
       "U has eigenvalue −1". The pure units U² = κI give rotations (κ = −1, tangent
       addition), shears (κ = 0, plain addition) and boosts (κ = 1, velocity addition, walls
       at ±1).
    5. A matrix weight inside the law: A ◯[K] B = A + B + AKB lifts through I + AK;
       Sylvester det(I + AK) = det(I + KA); a rank-one weight leaves no corner terms.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.LinearAlgebra

namespace DialSweep
open CircleOp LinearAlgebraCirc

set_option linter.unusedSectionVars false

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

/-! ## 1. One law in other units -/

theorem scale2 (k : R) (x y : M2 R) :
    M2.smul k (M2.circK k x y) = M2.circ (M2.smul k x) (M2.smul k y) := by
  apply M2.ext' <;> simp only [M2.smul, M2.circK, M2.circ, M2.add, M2.mul] <;> grind

theorem scale3 (k : R) (x y : M3 R) :
    M3.smul k (M3.circK k x y) = M3.circ (M3.smul k x) (M3.smul k y) := by
  apply M3.ext' <;> simp only [M3.smul, M3.circK, M3.circ, M3.add, M3.mul] <;> grind

/-! ## 2. The dial is the characteristic polynomial -/

/-- p_A(t) = t² − tr A·t + det A. -/
def charpoly2 (x : M2 R) (t : R) : R := t * t - M2.tr x * t + M2.det x
/-- p_A(t) = t³ − e₁t² + e₂t − e₃. -/
def charpoly3 (x : M3 R) (t : R) : R :=
  t * t * t - M3.tr x * t * t + M3.e2 x * t - M3.det x

/-- k = −1 is the characteristic polynomial at 1 (Fredholm's λ is −k). -/
theorem det_one_sub_2 (x : M2 R) :
    M2.det (M2.add M2.one (M2.smul (-1) x)) = charpoly2 x 1 := by
  rw [M2.det_one_add_smul]; unfold charpoly2; grind

theorem det_one_sub_3 (x : M3 R) :
    M3.det (M3.add M3.one (M3.smul (-1) x)) = charpoly3 x 1 := by
  rw [M3.det_one_add_smul]; unfold charpoly3; grind

/-- The top of the dial: det(I + kA) − k³·det A has only lower powers of k. -/
theorem det_top (k : R) (x : M3 R) :
    M3.det (M3.add M3.one (M3.smul k x)) - k * k * k * M3.det x
      = 1 + k * M3.tr x + k * k * M3.e2 x := by
  rw [M3.det_one_add_smul]; grind

/-! ## 3. Rotations at k = −1 -/

/-- At k = −1 the rotation condition is X + Xᵀ − XᵀX = 0: I − X is orthogonal. -/
theorem orth_neg_one (x : M3 R) :
    M3.mul (M3.transpose (M3.add M3.one (M3.smul (-1) x))) (M3.add M3.one (M3.smul (-1) x))
      = M3.add M3.one (M3.smul (-1) (M3.circK (-1) (M3.transpose x) x)) :=
  M3.orth_lift (-1) x

/-! ## 4. A matrix on the dial -/

/-- The dial matrix U = [[0, −q], [1, k]]. -/
def dialU (k q : R) : M2 R := ⟨0, -q, 1, k⟩

theorem dialU_sq (k q : R) :
    M2.mul (dialU k q) (dialU k q) = M2.sub (M2.smul k (dialU k q)) (M2.smul q M2.one) := by
  apply M2.ext' <;> simp only [M2.mul, M2.sub, M2.smul, M2.one, dialU] <;> grind

/-- (I + aU)(I + bU) = (1 − q·ab)·I + (a + b + k·ab)·U: the landscape law, lifted. -/
theorem dial_product (k q a b : R) :
    M2.mul (M2.add M2.one (M2.smul a (dialU k q))) (M2.add M2.one (M2.smul b (dialU k q)))
      = M2.add (M2.smul (1 - q * a * b) M2.one) (M2.smul (a + b + k * a * b) (dialU k q)) := by
  apply M2.ext' <;> simp only [M2.mul, M2.add, M2.smul, M2.one, dialU] <;> grind

/-- det(I + aU) = 1 + k·a + q·a², the landscape's wall polynomial. -/
theorem dial_det (k q a : R) :
    M2.det (M2.add M2.one (M2.smul a (dialU k q))) = 1 + k * a + q * a * a := by
  simp only [M2.det, M2.add, M2.smul, M2.one, dialU]; grind

/-- With eigenvalues λ₁, λ₂ (k = λ₁ + λ₂, q = λ₁λ₂) the walls are a = −1/λᵢ. -/
theorem dial_det_eigen (l1 l2 a : R) :
    M2.det (M2.add M2.one (M2.smul a (dialU (l1 + l2) (l1 * l2)))) = (1 + a * l1) * (1 + a * l2) := by
  rw [dial_det]; grind

/-- The three countries are the sign of the discriminant, the squared eigenvalue gap. -/
theorem dial_discriminant (l1 l2 : R) :
    (l1 + l2) * (l1 + l2) - 4 * (l1 * l2) = (l1 - l2) * (l1 - l2) := by
  grind

/-- The "or" line: a wall at certainty a = 1 exactly when 1 + k + q = 0. -/
theorem or_line (k q : R) : M2.det (M2.add M2.one (M2.smul 1 (dialU k q))) = 1 + k + q := by
  rw [dial_det]; grind

/-- The pure unit U = [[0, κ], [1, 0]], U² = κ·I. -/
def unitU (κ : R) : M2 R := ⟨0, κ, 1, 0⟩

theorem unitU_sq (κ : R) : M2.mul (unitU κ) (unitU κ) = M2.smul κ M2.one := by
  apply M2.ext' <;> simp only [M2.mul, M2.smul, M2.one, unitU] <;> grind

/-- (I + aU)(I + bU) = (1 + κ·ab)·I + (a + b)·U: projected, (a + b)/(1 + κ·ab). -/
theorem unit_product (κ a b : R) :
    M2.mul (M2.add M2.one (M2.smul a (unitU κ))) (M2.add M2.one (M2.smul b (unitU κ)))
      = M2.add (M2.smul (1 + κ * a * b) M2.one) (M2.smul (a + b) (unitU κ)) := by
  apply M2.ext' <;> simp only [M2.mul, M2.add, M2.smul, M2.one, unitU] <;> grind

/-- κ = 1 (boosts): det(I + aU) = (1 − a)(1 + a), walls at the speed of light ±1. -/
theorem boost_walls (a : R) :
    M2.det (M2.add M2.one (M2.smul a (unitU 1))) = (1 - a) * (1 + a) := by
  simp only [M2.det, M2.add, M2.smul, M2.one, unitU]; grind

/-- κ = −1 (rotations): det(I + aU) = 1 + a², no real wall. -/
theorem rotation_no_wall (a : R) :
    M2.det (M2.add M2.one (M2.smul a (unitU (-1)))) = 1 + a * a := by
  simp only [M2.det, M2.add, M2.smul, M2.one, unitU]; grind

/-- κ = 0 (shears): det(I + aU) = 1, and the law is plain addition. -/
theorem shear_is_add (a b : R) :
    M2.mul (M2.add M2.one (M2.smul a (unitU 0))) (M2.add M2.one (M2.smul b (unitU 0)))
      = M2.add M2.one (M2.smul (a + b) (unitU 0)) := by
  apply M2.ext' <;> simp only [M2.mul, M2.add, M2.smul, M2.one, unitU] <;> grind

/-! ## 5. A matrix weight inside the law -/

/-- A ◯[K] B = A + B + A·K·B. -/
def kcirc (kk x y : M2 R) : M2 R := M2.add (M2.add x y) (M2.mul (M2.mul x kk) y)

theorem kcirc_lift (kk x y : M2 R) :
    M2.add M2.one (M2.mul (kcirc kk x y) kk)
      = M2.mul (M2.add M2.one (M2.mul x kk)) (M2.add M2.one (M2.mul y kk)) := by
  apply M2.ext' <;> simp only [kcirc, M2.add, M2.mul, M2.one] <;> grind

theorem sylvester2 (x kk : M2 R) :
    M2.det (M2.add M2.one (M2.mul x kk)) = M2.det (M2.add M2.one (M2.mul kk x)) := by
  simp only [M2.det, M2.add, M2.mul, M2.one]; grind

set_option maxHeartbeats 2000000 in
theorem sylvester3 (x kk : M3 R) :
    M3.det (M3.add M3.one (M3.mul x kk)) = M3.det (M3.add M3.one (M3.mul kk x)) := by
  simp only [M3.det, M3.add, M3.mul, M3.one]; grind

/-- The rank-one matrix u·vᵀ. -/
def outer3 (u1 u2 u3 v1 v2 v3 : R) : M3 R :=
  ⟨u1 * v1, u1 * v2, u1 * v3, u2 * v1, u2 * v2, u2 * v3, u3 * v1, u3 * v2, u3 * v3⟩

/-- A rank-one weight leaves no corner terms: det(I + u·vᵀ) = 1 + vᵀu. -/
theorem rank_one3 (u1 u2 u3 v1 v2 v3 : R) :
    M3.det (M3.add M3.one (outer3 u1 u2 u3 v1 v2 v3)) = 1 + (v1 * u1 + v2 * u2 + v3 * u3) := by
  simp only [M3.det, M3.add, M3.one, outer3]; grind

/-- The matrix determinant lemma: det(I + A·u·vᵀ) = 1 + vᵀ·A·u. -/
theorem rank_one_dial3 (x : M3 R) (u1 u2 u3 v1 v2 v3 : R) :
    M3.det (M3.add M3.one (M3.mul x (outer3 u1 u2 u3 v1 v2 v3)))
      = 1 + (v1 * (x.a11 * u1 + x.a12 * u2 + x.a13 * u3)
           + v2 * (x.a21 * u1 + x.a22 * u2 + x.a23 * u3)
           + v3 * (x.a31 * u1 + x.a32 * u2 + x.a33 * u3)) := by
  simp only [M3.det, M3.add, M3.one, M3.mul, outer3]; grind

end Ring

/-! ## 6. Fields: the walls, the reversed polynomial, the projected laws -/

section Field
variable {K : Type} [Lean.Grind.Field K]

/-- det(I + kA) = k²·p_A(−1/k) (2×2), for k ≠ 0. -/
theorem det_charpoly2 (k : K) (hk : k ≠ 0) (x : M2 K) :
    M2.det (M2.add M2.one (M2.smul k x)) = k * k * charpoly2 x (-1 / k) := by
  rw [M2.det_one_add_smul]; unfold charpoly2; grind

/-- det(I + kA) = (−k)³·p_A(−1/k) (3×3), for k ≠ 0. -/
theorem det_charpoly3 (k : K) (hk : k ≠ 0) (x : M3 K) :
    M3.det (M3.add M3.one (M3.smul k x)) = -(k * k * k) * charpoly3 x (-1 / k) := by
  rw [M3.det_one_add_smul]; unfold charpoly3; grind

/-- The matrix sits on the wall of weight k exactly when −1/k is an eigenvalue. -/
theorem wall_iff_eigen2 (k : K) (hk : k ≠ 0) (x : M2 K) :
    M2.det (M2.add M2.one (M2.smul k x)) = 0 ↔ charpoly2 x (-1 / k) = 0 := by
  rw [det_charpoly2 k hk x]
  constructor
  · intro h; grind
  · intro h; rw [h]; grind

/-- The landscape lens law comes out of the dial matrix: off the wall 1 − q·ab ≠ 0,
    (I + aU)(I + bU) = (1 − q·ab)·(I + cU) with c = (a + b + k·ab)/(1 − q·ab). -/
theorem dial_lens (k q a b : K) (h : 1 - q * a * b ≠ 0) :
    M2.mul (M2.add M2.one (M2.smul a (dialU k q))) (M2.add M2.one (M2.smul b (dialU k q)))
      = M2.smul (1 - q * a * b)
          (M2.add M2.one (M2.smul ((a + b + k * a * b) / (1 - q * a * b)) (dialU k q))) := by
  apply M2.ext' <;> simp only [M2.mul, M2.add, M2.smul, M2.one, dialU] <;> grind

/-- The weighted Cayley transform (I + kS)⁻¹(I − kS) is the Cayley transform of k·s. -/
theorem cayley_weighted (k s : K) :
    M2.mul (M2.inv (M2.add M2.one (M2.smul k (M2.skew s)))) (M2.sub M2.one (M2.smul k (M2.skew s)))
      = M2.cayley (k * s) := by
  have e : M2.smul k (M2.skew s) = M2.skew (k * s) := by
    apply M2.ext' <;> simp only [M2.smul, M2.skew] <;> grind
  rw [e, M2.cayley]

/-- The weighted inverse: off the weighted wall, a ◯[k] (−a/(1 + k·a)) = 0. -/
theorem winv (k a : K) (h : 1 + k * a ≠ 0) : circK k a (-a / (1 + k * a)) = 0 := by
  unfold circK; grind

end Field

end DialSweep
