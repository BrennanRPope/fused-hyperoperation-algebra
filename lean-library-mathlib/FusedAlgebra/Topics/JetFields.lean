module

public import FusedAlgebra.Topics.LinearAlgebra

/-!
# Fields of jets, mult-add, and the velocity gradient

This is the mathlib form of `FusedAlgebra.Topics.JetFields`.

* A 1-jet with corner weight `k` is a point of `WDual k = R[ε]/(ε² - k ε)`. With `k = 0` these
  are the dual numbers (`WDual.dualNumberEquiv`); with `k ≠ 0` they are the values at two points
  `k` apart (`WDual.twoPt`, a ring homomorphism to `R × R`, and `WDual.twoPtEquiv`). Units
  `1 + a ε` multiply by weighted mult-add (`WDual.mk_one_mul_mk_one`). The wall `1 + k * a = 0`
  is where the unit vanishes at the second point (`WDual.wall_no_inverse`).
* `Jet2 R = R[ε]/(ε³)` holds 2-jets. Relative jets mult-add with the cube cut after layer 2
  (`Jet2.cube_depth_two`), and the way back is a finite sum (`Jet2.nilpotent_inverse`).
* The step `x ↦ x + h • u x` has volume change `det (1 + h • ∇u)`, the mult-add cube of the
  gradient (`LinearAlgebraMA.det_one_add_smul_fin_three`, `det_diag_cube`). For a traceless
  gradient, layer 2 is rotation against strain (`e2_vorticity_strain`) and layer 3 is strain
  self-amplification plus vortex stretching (`det_vorticity_strain`).
* Restricted Euler (Vieillefosse 1982, Cantwell 1992): with `reStep X = -3 • X² + trace X² • 1`,
  the layers move by `re_dQ`, `re_dR`, and `27 R² + 4 Q³` is conserved (`re_conserved`). In two
  dimensions the flow is frozen (`re2_frozen`).
* On a periodic grid the 2×2 minor of a discrete gradient, read with the corner kept, is an exact
  discrete divergence (`minor_is_divergence`) and sums to zero over the box (`box_minor_zero`).

Known mathematics, and where it is in mathlib: the dual numbers are `DualNumber R`;
Cayley–Hamilton is `Matrix.aeval_self_charpoly` (the 3×3 case is restated here as
`cayley_hamilton` in terms of trace, `e2` and determinant); telescoping sums are
`Finset.sum_range_sub`; forward differences are `fwdDiff`. `WDual k` and `Jet2 R` are the
quotients `R[X]/(X² - k X)` and `R[X]/(X³)`; they are kept as two- and three-field structures
with a `CommRing` instance, since that is what the statements use.

Not here (needs real analysis or integrals): the continuous null-Lagrangian integrals, Betchov's
relation, the blow-up time of the Vieillefosse tail, anything about actual Navier–Stokes
solutions.
-/

@[expose] public section

namespace MultAdd.JetFields

open MultAdd LinearAlgebraMA Matrix Finset

/-! ### Weighted dual numbers: `R[ε]/(ε² - k ε)` -/

/-- An element `re + du * ε` of `R[ε]/(ε² - k ε)`: a value and a slope. The weight `k` is part
of the type, so that the ring structure is an instance. -/
@[ext]
structure WDual {R : Type*} (k : R) where
  /-- The value. -/
  re : R
  /-- The slope: the coefficient of `ε`. -/
  du : R

namespace WDual

section CommRing

variable {R : Type*} [CommRing R] {k : R}

instance : Zero (WDual k) := ⟨⟨0, 0⟩⟩

instance : One (WDual k) := ⟨⟨1, 0⟩⟩

instance : Add (WDual k) := ⟨fun x y ↦ ⟨x.re + y.re, x.du + y.du⟩⟩

instance : Neg (WDual k) := ⟨fun x ↦ ⟨-x.re, -x.du⟩⟩

instance : Sub (WDual k) := ⟨fun x y ↦ ⟨x.re - y.re, x.du - y.du⟩⟩

/-- Multiplication in `R[ε]/(ε² - k ε)`. -/
instance : Mul (WDual k) :=
  ⟨fun x y ↦ ⟨x.re * y.re, x.re * y.du + x.du * y.re + k * x.du * y.du⟩⟩

@[simp] theorem zero_re : (0 : WDual k).re = 0 := rfl
@[simp] theorem zero_du : (0 : WDual k).du = 0 := rfl
@[simp] theorem one_re : (1 : WDual k).re = 1 := rfl
@[simp] theorem one_du : (1 : WDual k).du = 0 := rfl
@[simp] theorem add_re (x y : WDual k) : (x + y).re = x.re + y.re := rfl
@[simp] theorem add_du (x y : WDual k) : (x + y).du = x.du + y.du := rfl
@[simp] theorem neg_re (x : WDual k) : (-x).re = -x.re := rfl
@[simp] theorem neg_du (x : WDual k) : (-x).du = -x.du := rfl
@[simp] theorem sub_re (x y : WDual k) : (x - y).re = x.re - y.re := rfl
@[simp] theorem sub_du (x y : WDual k) : (x - y).du = x.du - y.du := rfl
@[simp] theorem mul_re (x y : WDual k) : (x * y).re = x.re * y.re := rfl
@[simp] theorem mul_du (x y : WDual k) :
    (x * y).du = x.re * y.du + x.du * y.re + k * x.du * y.du := rfl

/-- `R[ε]/(ε² - k ε)` is a commutative ring. -/
instance instCommRing : CommRing (WDual k) where
  add_assoc x y z := by ext <;> simp [add_assoc]
  zero_add x := by ext <;> simp
  add_zero x := by ext <;> simp
  add_comm x y := by ext <;> simp [add_comm]
  nsmul := nsmulRec
  zsmul := zsmulRec
  neg_add_cancel x := by ext <;> simp
  sub_eq_add_neg x y := by ext <;> simp [sub_eq_add_neg]
  left_distrib x y z := by ext <;> simp <;> ring
  right_distrib x y z := by ext <;> simp <;> ring
  zero_mul x := by ext <;> simp
  mul_zero x := by ext <;> simp
  mul_assoc x y z := by ext <;> simp <;> ring
  one_mul x := by ext <;> simp
  mul_one x := by ext <;> simp
  mul_comm x y := by ext <;> simp <;> ring

/-- `ε * ε = k * ε`: the corner weight is the square of the infinitesimal. -/
theorem eps_mul_eps : (⟨0, 1⟩ : WDual k) * ⟨0, 1⟩ = ⟨0, k⟩ := by
  ext <;> simp

/-- Units multiply by weighted mult-add: `(1 + a ε) * (1 + b ε) = 1 + (wAdd k a b) ε`. -/
theorem mk_one_mul_mk_one (a b : R) : (⟨1, a⟩ : WDual k) * ⟨1, b⟩ = ⟨1, wAdd k a b⟩ := by
  ext <;> simp [wAdd]; ring

/-- With weight `1` the units multiply by plain mult-add. -/
theorem units_multAdd (a b : R) : (⟨1, a⟩ : WDual (1 : R)) * ⟨1, b⟩ = ⟨1, a ✳ b⟩ := by
  rw [mk_one_mul_mk_one, wAdd_one]

/-- With weight `0` (dual numbers, 1-jets) the corner is dropped: relative slopes add. This is
the logarithmic derivative of a product. -/
theorem dual_units_add (a b : R) : (⟨1, a⟩ : WDual (0 : R)) * ⟨1, b⟩ = ⟨1, a + b⟩ := by
  rw [mk_one_mul_mk_one, wAdd_weight_zero]

/-- At order 1 only layer 1 survives: relative jets just add. -/
theorem cube_depth_one (a₁ a₂ a₃ : R) :
    (⟨1, a₁⟩ : WDual (0 : R)) * ⟨1, a₂⟩ * ⟨1, a₃⟩ = ⟨1, a₁ + a₂ + a₃⟩ := by
  rw [dual_units_add, dual_units_add]

/-- The two-point reading, a ring homomorphism: the value at `x` and the value at `x + k`
(slope times gap added). It turns the twisted product into the plain product at each point. -/
def twoPt : WDual k →+* R × R where
  toFun x := (x.re, x.re + k * x.du)
  map_one' := by simp
  map_mul' x y := by ext <;> simp; ring
  map_zero' := by simp
  map_add' x y := by ext <;> simp; ring

@[simp] theorem twoPt_apply (x : WDual k) : twoPt x = (x.re, x.re + k * x.du) := rfl

/-- At weight `0` the two points coincide: the secant has become the tangent. -/
theorem twoPt_weight_zero (x : WDual (0 : R)) : twoPt x = (x.re, x.re) := by simp

/-- Discrete Leibniz with the corner. If `f (x + h) = f₀ + h * f'` and `g (x + h) = g₀ + h * g'`
then `(f * g) (x + h) = f₀ * g₀ + h * (the ε-part of the product in R[ε]/(ε² - h ε))`. -/
theorem diff_leibniz (h f₀ f' g₀ g' : R) :
    (f₀ + h * f') * (g₀ + h * g') = f₀ * g₀ + h * ((⟨f₀, f'⟩ : WDual h) * ⟨g₀, g'⟩).du := by
  simp; ring

/-- The dual numbers have no wall: `1 + a ε` always has the inverse `1 - a ε`. -/
theorem dual_inverse (a : R) : (⟨1, a⟩ : WDual (0 : R)) * ⟨1, -a⟩ = 1 := by
  ext <;> simp

/-- At weight `0` the weighted dual numbers are mathlib's dual numbers `R[ε]`. -/
def dualNumberEquiv : WDual (0 : R) ≃+* DualNumber R where
  toFun x := TrivSqZeroExt.inl x.re + TrivSqZeroExt.inr x.du
  invFun y := ⟨y.fst, y.snd⟩
  left_inv x := by ext <;> simp
  right_inv y := by ext <;> simp
  map_mul' x y := by ext <;> simp [mul_comm]
  map_add' x y := by ext <;> simp

end CommRing

section Field

variable {K : Type*} [Field K] {k : K}

/-- From two point values back to value and slope (useful for a nonzero gap `k`). -/
def ofTwoPt (k : K) (p : K × K) : WDual k := ⟨p.1, (p.2 - p.1) / k⟩

theorem twoPt_ofTwoPt (hk : k ≠ 0) (p : K × K) : twoPt (ofTwoPt k p) = p := by
  ext <;> simp [ofTwoPt]
  field_simp; ring

theorem ofTwoPt_twoPt (hk : k ≠ 0) (x : WDual k) : ofTwoPt k (twoPt x) = x := by
  ext <;> simp [ofTwoPt]
  field_simp

/-- For a nonzero gap the weighted dual numbers are two points: `WDual k ≃+* K × K`. -/
def twoPtEquiv (hk : k ≠ 0) : WDual k ≃+* K × K :=
  { twoPt with
    invFun := ofTwoPt k
    left_inv := ofTwoPt_twoPt hk
    right_inv := twoPt_ofTwoPt hk }

/-- The wall is a zero at the second point. If `1 + k * a = 0`, the unit `1 + a ε` has no inverse
in `R[ε]/(ε² - k ε)`: read at `x + k`, it is `0`. -/
theorem wall_no_inverse (a : K) (hw : 1 + k * a = 0) (y : WDual k) :
    (⟨1, a⟩ : WDual k) * y ≠ 1 := by
  intro h
  have h2 := congrArg (fun z ↦ (twoPt z).2) h
  simp only [map_mul, map_one, Prod.snd_mul, Prod.snd_one, twoPt_apply, hw, zero_mul] at h2
  exact zero_ne_one h2

/-- Off the wall, the inverse is the weighted mult-add way back `-a / (1 + k * a)`. -/
theorem off_wall_inverse (a : K) (hw : 1 + k * a ≠ 0) :
    (⟨1, a⟩ : WDual k) * ⟨1, -a / (1 + k * a)⟩ = 1 := by
  ext
  · simp
  · have e : 1 * (-a / (1 + k * a)) + a * 1 + k * a * (-a / (1 + k * a)) =
        a * (1 - (1 + k * a) / (1 + k * a)) := by ring
    rw [mul_du, e, div_self hw, sub_self, mul_zero, one_du]

end Field

end WDual

/-! ### Second-order jets: the cube cut after layer 2 -/

/-- A 2-jet as Taylor coefficients `(f, f', f'' / 2)`: an element of `R[ε]/(ε³)`. -/
@[ext]
structure Jet2 (R : Type*) where
  /-- The value. -/
  c0 : R
  /-- The first Taylor coefficient. -/
  c1 : R
  /-- The second Taylor coefficient. -/
  c2 : R

namespace Jet2

variable {R : Type*} [CommRing R]

instance : Zero (Jet2 R) := ⟨⟨0, 0, 0⟩⟩

instance : One (Jet2 R) := ⟨⟨1, 0, 0⟩⟩

instance : Add (Jet2 R) := ⟨fun x y ↦ ⟨x.c0 + y.c0, x.c1 + y.c1, x.c2 + y.c2⟩⟩

instance : Neg (Jet2 R) := ⟨fun x ↦ ⟨-x.c0, -x.c1, -x.c2⟩⟩

instance : Sub (Jet2 R) := ⟨fun x y ↦ ⟨x.c0 - y.c0, x.c1 - y.c1, x.c2 - y.c2⟩⟩

/-- Multiplication in `R[ε]/(ε³)`. -/
instance : Mul (Jet2 R) :=
  ⟨fun x y ↦ ⟨x.c0 * y.c0, x.c0 * y.c1 + x.c1 * y.c0,
    x.c0 * y.c2 + x.c1 * y.c1 + x.c2 * y.c0⟩⟩

@[simp] theorem zero_c0 : (0 : Jet2 R).c0 = 0 := rfl
@[simp] theorem zero_c1 : (0 : Jet2 R).c1 = 0 := rfl
@[simp] theorem zero_c2 : (0 : Jet2 R).c2 = 0 := rfl
@[simp] theorem one_c0 : (1 : Jet2 R).c0 = 1 := rfl
@[simp] theorem one_c1 : (1 : Jet2 R).c1 = 0 := rfl
@[simp] theorem one_c2 : (1 : Jet2 R).c2 = 0 := rfl
@[simp] theorem add_c0 (x y : Jet2 R) : (x + y).c0 = x.c0 + y.c0 := rfl
@[simp] theorem add_c1 (x y : Jet2 R) : (x + y).c1 = x.c1 + y.c1 := rfl
@[simp] theorem add_c2 (x y : Jet2 R) : (x + y).c2 = x.c2 + y.c2 := rfl
@[simp] theorem neg_c0 (x : Jet2 R) : (-x).c0 = -x.c0 := rfl
@[simp] theorem neg_c1 (x : Jet2 R) : (-x).c1 = -x.c1 := rfl
@[simp] theorem neg_c2 (x : Jet2 R) : (-x).c2 = -x.c2 := rfl
@[simp] theorem sub_c0 (x y : Jet2 R) : (x - y).c0 = x.c0 - y.c0 := rfl
@[simp] theorem sub_c1 (x y : Jet2 R) : (x - y).c1 = x.c1 - y.c1 := rfl
@[simp] theorem sub_c2 (x y : Jet2 R) : (x - y).c2 = x.c2 - y.c2 := rfl
@[simp] theorem mul_c0 (x y : Jet2 R) : (x * y).c0 = x.c0 * y.c0 := rfl
@[simp] theorem mul_c1 (x y : Jet2 R) : (x * y).c1 = x.c0 * y.c1 + x.c1 * y.c0 := rfl
@[simp] theorem mul_c2 (x y : Jet2 R) :
    (x * y).c2 = x.c0 * y.c2 + x.c1 * y.c1 + x.c2 * y.c0 := rfl

/-- `R[ε]/(ε³)` is a commutative ring. -/
instance instCommRing : CommRing (Jet2 R) where
  add_assoc x y z := by ext <;> simp [add_assoc]
  zero_add x := by ext <;> simp
  add_zero x := by ext <;> simp
  add_comm x y := by ext <;> simp [add_comm]
  nsmul := nsmulRec
  zsmul := zsmulRec
  neg_add_cancel x := by ext <;> simp
  sub_eq_add_neg x y := by ext <;> simp [sub_eq_add_neg]
  left_distrib x y z := by ext <;> simp <;> ring
  right_distrib x y z := by ext <;> simp <;> ring
  zero_mul x := by ext <;> simp
  mul_zero x := by ext <;> simp
  mul_assoc x y z := by ext <;> simp <;> ring
  one_mul x := by ext <;> simp
  mul_one x := by ext <;> simp
  mul_comm x y := by ext <;> simp <;> ring

/-- The jet product is the product of functions: second-order Leibniz,
`(f * g)'' = f'' * g + 2 * f' * g' + f * g''` (with `f'' = 2 * c2`). -/
theorem leibniz2 (x y : Jet2 R) :
    2 * (x * y).c2 = x.c0 * (2 * y.c2) + 2 * (x.c1 * y.c1) + (2 * x.c2) * y.c0 := by
  simp; ring

/-- Jet order is cube depth. Three relative 2-jets (no constant term) mult-add to their sum plus
the three pair corners; the triple corner, which would sit at `ε³`, is gone. -/
theorem cube_depth_two (a₁ b₁ a₂ b₂ a₃ b₃ : R) :
    (⟨0, a₁, b₁⟩ : Jet2 R) ✳ (⟨0, a₂, b₂⟩ ✳ ⟨0, a₃, b₃⟩) =
      ⟨0, a₁ + a₂ + a₃, (b₁ + b₂ + b₃) + (a₁ * a₂ + a₁ * a₃ + a₂ * a₃)⟩ := by
  ext <;> simp [multAdd] <;> ring

/-- No wall in jet space. A relative jet `ν` has the finite way back `-ν + ν ^ 2` (the geometric
series stops because `ε³ = 0`). -/
theorem nilpotent_inverse (a b : R) : (⟨0, a, b⟩ : Jet2 R) ✳ ⟨0, -a, -b + a ^ 2⟩ = 0 := by
  ext <;> simp [multAdd]; ring

end Jet2

/-! ### The velocity gradient: the mult-add cube of the step `x ↦ x + h • u x` -/

section Gradient

variable {R : Type*} [CommRing R]

/-- In principal axes the step's volume change is literally mult-add:
`det (1 + h • diagonal ![a, b, c]) = 1 + (h * a ✳ h * b ✳ h * c)`. -/
theorem det_diag_cube (h a b c : R) :
    det (1 + h • diagonal ![a, b, c]) = 1 + ((h * a) ✳ (h * b) ✳ (h * c)) := by
  simp [det_fin_three, multAdd]; ring

/-- The layers of a diagonal gradient are the cube's layers of `(a, b, c)`. -/
theorem diag_layers (a b c : R) :
    trace (diagonal ![a, b, c]) = a + b + c ∧
      e2 (diagonal ![a, b, c]) = a * b + a * c + b * c ∧
      det (diagonal ![a, b, c]) = a * b * c := by
  refine ⟨?_, ?_, ?_⟩
  · simp [Fin.sum_univ_three]
  · simp [e2]
  · simp [Fin.prod_univ_three]

/-- Cayley–Hamilton for 3×3 matrices, with the layers as coefficients:
`X ^ 3 - P • X ^ 2 + Q • X - R • 1 = 0`. In mathlib this is `Matrix.aeval_self_charpoly`. -/
theorem cayley_hamilton (X : Matrix (Fin 3) (Fin 3) R) :
    X ^ 3 - trace X • X ^ 2 + e2 X • X - det X • (1 : Matrix (Fin 3) (Fin 3) R) = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [pow_succ, Matrix.mul_apply, Fin.sum_univ_three, trace_fin_three, e2, det_fin_three] <;>
    ring

/-- Newton's identity at degree 2: `trace (X * X) = P ^ 2 - 2 * Q`. -/
theorem trace_mul_self_fin_three (X : Matrix (Fin 3) (Fin 3) R) :
    trace (X * X) = trace X ^ 2 - 2 * e2 X := by
  simp [Matrix.mul_apply, Fin.sum_univ_three, trace_fin_three, e2]; ring

/-- Newton's identity at degree 3: `trace (X * (X * X)) = P ^ 3 - 3 * P * Q + 3 * R`. -/
theorem trace_mul_mul_self_fin_three (X : Matrix (Fin 3) (Fin 3) R) :
    trace (X * (X * X)) = trace X ^ 3 - 3 * trace X * e2 X + 3 * det X := by
  simp [Matrix.mul_apply, Fin.sum_univ_three, trace_fin_three, e2, det_fin_three]; ring

/-- For an incompressible gradient, `trace (X * X) = -2 * Q`. -/
theorem tr_sq_traceless (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    trace (X * X) = -2 * e2 X := by
  rw [trace_mul_self_fin_three, h0]; ring

/-- For an incompressible gradient, `trace (X * (X * X)) = 3 * R`. -/
theorem tr_cube_traceless (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    trace (X * (X * X)) = 3 * det X := by
  rw [trace_mul_mul_self_fin_three, h0]; ring

/-- The vorticity vector `ω = curl u` read off the gradient. -/
def vort (X : Matrix (Fin 3) (Fin 3) R) : Fin 3 → R :=
  ![X 2 1 - X 1 2, X 0 2 - X 2 0, X 1 0 - X 0 1]

/-- The squared Frobenius norm: the sum of the squares of all entries, `trace (X * Xᵀ)`. -/
def fsq {n : Type*} [Fintype n] (X : Matrix n n R) : R := trace (X * Xᵀ)

/-- Layer 2 against rotation and strain, for any gradient: with `T = X + Xᵀ` (twice the strain)
and `ω` the vorticity, `8 * Q = 2 * |ω|² - |T|² + 4 * P ^ 2`. -/
theorem e2_vorticity_strain_add (X : Matrix (Fin 3) (Fin 3) R) :
    8 * e2 X = 2 * (vort X ⬝ᵥ vort X) - fsq (X + Xᵀ) + 4 * trace X ^ 2 := by
  simp [vort, fsq, e2, Matrix.mul_apply, Fin.sum_univ_three, trace_fin_three, dotProduct]
  ring

/-- Layer 2 is rotation against strain. For `P = 0`: `8 * Q = 2 * |ω|² - |T|²`, that is,
`Q = |ω|² / 4 - |S|² / 2`. -/
theorem e2_vorticity_strain (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    8 * e2 X = 2 * (vort X ⬝ᵥ vort X) - fsq (X + Xᵀ) := by
  rw [e2_vorticity_strain_add, h0]; ring

/-- Layer 3, the far corner, is strain self-amplification plus vortex stretching. For `P = 0`:
`8 * R = det T + ω ⬝ᵥ T *ᵥ ω`, that is, `R = det S + (ω ⬝ᵥ S *ᵥ ω) / 4`. -/
theorem det_vorticity_strain (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    8 * det X = det (X + Xᵀ) + vort X ⬝ᵥ (X + Xᵀ) *ᵥ vort X := by
  have h : X 2 2 = -X 0 0 - X 1 1 := by
    rw [trace_fin_three] at h0; linear_combination h0
  simp [vort, det_fin_three, mulVec, dotProduct, Fin.sum_univ_three, h]
  ring

/-- For a symmetric traceless `S` (a pure strain), `trace (S * (S * S)) = 3 * det S`. -/
theorem sym_tr_cube (a b c e f : R) :
    trace (!![a, b, c; b, e, f; c, f, -a - e] * (!![a, b, c; b, e, f; c, f, -a - e] *
      !![a, b, c; b, e, f; c, f, -a - e])) = 3 * det !![a, b, c; b, e, f; c, f, -a - e] :=
  tr_cube_traceless _ (by simp [trace_fin_three])

/-- For a pure traceless strain, `-2 * Q` is the sum of the squares of its entries. -/
theorem sym_e2_sumsq (a b c e f : R) :
    -2 * e2 !![a, b, c; b, e, f; c, f, -a - e] = fsq !![a, b, c; b, e, f; c, f, -a - e] := by
  have ht : !![a, b, c; b, e, f; c, f, -a - e]ᵀ = !![a, b, c; b, e, f; c, f, -a - e] := by
    ext i j
    rw [transpose_apply]
    fin_cases i <;> fin_cases j <;> rfl
  rw [fsq, ht]
  simp [e2, trace_fin_three]
  ring

/-- In an ordered ring `Q ≤ 0` for a pure strain: the strain's second layer only ever
dissipates. -/
theorem sym_e2_nonpos {R : Type*} [CommRing R] [LinearOrder R] [IsStrictOrderedRing R]
    (a b c e f : R) : e2 !![a, b, c; b, e, f; c, f, -a - e] ≤ 0 := by
  simp [e2]
  nlinarith [sq_nonneg a, sq_nonneg b, sq_nonneg c, sq_nonneg e, sq_nonneg f, sq_nonneg (a + e)]

end Gradient

/-! ### Restricted Euler: the layers drive themselves -/

section RestrictedEuler

variable {R : Type*} [CommRing R]

/-- The linear part of `Q = e2` at `X` in direction `Y`:
`e2 (X + t • Y) = e2 X + t * dE2 X Y + t ^ 2 * e2 Y`. -/
def dE2 (X Y : Matrix (Fin 3) (Fin 3) R) : R := trace X * trace Y - trace (X * Y)

/-- The linear part of `det` at `X` in direction `Y`: one row from `Y` at a time. -/
def dDet {n : Type*} [Fintype n] [DecidableEq n] (X Y : Matrix n n R) : R :=
  ∑ i, det (X.updateRow i (Y i))

/-- `dE2` really is the first-order part of `Q`. -/
theorem e2_expand (X Y : Matrix (Fin 3) (Fin 3) R) (t : R) :
    e2 (X + t • Y) = e2 X + t * dE2 X Y + t ^ 2 * e2 Y := by
  simp [e2, dE2, Matrix.mul_apply, Fin.sum_univ_three, trace_fin_three]; ring

/-- `dDet` really is the first-order part of `det` (the second-order part is `dDet Y X`). -/
theorem det_expand (X Y : Matrix (Fin 3) (Fin 3) R) (t : R) :
    det (X + t • Y) = det X + t * dDet X Y + t ^ 2 * dDet Y X + t ^ 3 * det Y := by
  simp [dDet, det_fin_three, Fin.sum_univ_three, updateRow_apply]; ring

/-- Three times the restricted Euler velocity: `3 * dA/dt = -3 • A ^ 2 + trace (A ^ 2) • 1`
(self-advection plus the isotropic part of the pressure Hessian). -/
def reStep (X : Matrix (Fin 3) (Fin 3) R) : Matrix (Fin 3) (Fin 3) R :=
  (-3 : R) • (X * X) + trace (X * X) • (1 : Matrix (Fin 3) (Fin 3) R)

/-- The flow keeps the gradient incompressible. -/
theorem re_traceless (X : Matrix (Fin 3) (Fin 3) R) : trace (reStep X) = 0 := by
  simp [reStep, trace_fin_three]; ring

/-- Layer 2 is driven by layer 3: `3 * dQ/dt = 9 * R` (so `dQ/dt = 3 * R`). -/
theorem re_dQ (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    dE2 X (reStep X) = 9 * det X := by
  have h : trace (X * reStep X) = -3 * trace (X * (X * X)) + trace (X * X) * trace X := by
    simp [reStep, mul_add, trace_smul]
  rw [dE2, re_traceless, h, tr_cube_traceless X h0, h0]; ring

/-- Layer 3 is driven by layer 2 squared: `3 * dR/dt = -2 * Q ^ 2` (so `dR/dt = -2 / 3 * Q ^ 2`). -/
theorem re_dR (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    dDet X (reStep X) = -2 * e2 X ^ 2 := by
  have h : X 2 2 = -X 0 0 - X 1 1 := by
    rw [trace_fin_three] at h0; linear_combination h0
  simp [dDet, reStep, e2, det_fin_three, Fin.sum_univ_three, updateRow_apply, Matrix.mul_apply,
    trace_fin_three, h]
  ring

/-- The discriminant is conserved: along restricted Euler,
`d/dt (27 * R ^ 2 + 4 * Q ^ 3) = 54 * R * R' + 12 * Q ^ 2 * Q' = 0`. The quantity
`27 * R ^ 2 + 4 * Q ^ 3` is minus the discriminant of `λ ^ 3 + Q * λ - R`; its zero set is the
Vieillefosse tail. -/
theorem re_conserved (X : Matrix (Fin 3) (Fin 3) R) (h0 : trace X = 0) :
    54 * det X * dDet X (reStep X) + 12 * e2 X ^ 2 * dE2 X (reStep X) = 0 := by
  rw [re_dQ X h0, re_dR X h0]; ring

end RestrictedEuler

/-! ### Two dimensions: the square has no far corner -/

section TwoD

variable {R : Type*} [CommRing R]

/-- In two dimensions an incompressible gradient has no third invariant:
`trace (X * (X * X)) = 0`. -/
theorem tr_cube_2d (a b c : R) :
    trace (!![a, b; c, -a] * (!![a, b; c, -a] * !![a, b; c, -a])) = 0 := by
  simp [trace_fin_two]; ring

/-- Two-dimensional restricted Euler is frozen: `-2 • A ^ 2 + trace (A ^ 2) • 1 = 0` for every
traceless 2×2 gradient, so the self-advection is cancelled exactly by the pressure. -/
theorem re2_frozen (a b c : R) :
    (-2 : R) • (!![a, b; c, -a] * !![a, b; c, -a]) +
      trace (!![a, b; c, -a] * !![a, b; c, -a]) • (1 : Matrix (Fin 2) (Fin 2) R) = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [trace_fin_two] <;> ring

end TwoD

/-! ### Null Lagrangians on a periodic grid: the corner keeps the minor a divergence -/

section Grid

variable {R : Type*} [CommRing R]

/-- The forward difference in the first grid direction. This is `fwdDiff 1` in the first
variable. -/
def d1 (u : ℕ → ℕ → R) (i j : ℕ) : R := u (i + 1) j - u i j

/-- The forward difference in the second grid direction. This is `fwdDiff 1` in the second
variable. -/
def d2 (u : ℕ → ℕ → R) (i j : ℕ) : R := u i (j + 1) - u i j

theorem d1_eq_fwdDiff (u : ℕ → ℕ → R) (i j : ℕ) : d1 u i j = fwdDiff 1 (fun i ↦ u i j) i := rfl

theorem d2_eq_fwdDiff (u : ℕ → ℕ → R) (i j : ℕ) : d2 u i j = fwdDiff 1 (u i) j := rfl

/-- The discrete 2×2 minor of the map `(u, v)`, read with the corner kept: each second factor is
taken one step along the first factor's direction. -/
def minor (u v : ℕ → ℕ → R) (i j : ℕ) : R :=
  d1 u i j * d2 v (i + 1) j - d2 u i j * d1 v i (j + 1)

/-- The shift is the corner: `S₁ Δ₂ v = Δ₂ v + Δ₁ Δ₂ v` (shift is `1 + Δ`). -/
theorem shift_is_corner (v : ℕ → ℕ → R) (i j : ℕ) :
    d2 v (i + 1) j = d2 v i j + d1 (d2 v) i j := by
  simp only [d1, d2]; ring

/-- The corner-kept minor is an exact discrete divergence:
`minor = Δ₁ (u * Δ₂ v) - Δ₂ (u * Δ₁ v)`. -/
theorem minor_is_divergence (u v : ℕ → ℕ → R) (i j : ℕ) :
    minor u v i j =
      d1 (fun i j ↦ u i j * d2 v i j) i j - d2 (fun i j ↦ u i j * d1 v i j) i j := by
  simp only [minor, d1, d2]; ring

/-- The minor sums to zero over a periodic box: the grid version of "the integral of layer 2 of
`∇u` vanishes" for a map of the 2-torus. -/
theorem box_minor_zero (u v : ℕ → ℕ → R) (n : ℕ)
    (pu1 : ∀ i j, u (i + n) j = u i j) (pu2 : ∀ i j, u i (j + n) = u i j)
    (pv1 : ∀ i j, v (i + n) j = v i j) (pv2 : ∀ i j, v i (j + n) = v i j) :
    ∑ i ∈ range n, ∑ j ∈ range n, minor u v i j = 0 := by
  have hg : ∀ j, u n j * d2 v n j = u 0 j * d2 v 0 j := fun j ↦ by
    have h1 := pu1 0 j
    have h2 := pv1 0 (j + 1)
    have h3 := pv1 0 j
    simp only [zero_add] at h1 h2 h3
    rw [d2, d2, h1, h2, h3]
  have hg' : ∀ i, u i n * d1 v i n = u i 0 * d1 v i 0 := fun i ↦ by
    have h1 := pu2 i 0
    have h2 := pv2 (i + 1) 0
    have h3 := pv2 i 0
    simp only [zero_add] at h1 h2 h3
    rw [d1, d1, h1, h2, h3]
  have inner : ∀ i, ∑ j ∈ range n, minor u v i j =
      (∑ j ∈ range n, u (i + 1) j * d2 v (i + 1) j - ∑ j ∈ range n, u i j * d2 v i j) -
        (u i n * d1 v i n - u i 0 * d1 v i 0) := fun i ↦ by
    rw [← sum_range_sub (fun j ↦ u i j * d1 v i j) n, ← sum_sub_distrib, ← sum_sub_distrib]
    exact sum_congr rfl fun j _ ↦ minor_is_divergence u v i j
  simp only [inner, hg', sub_self, sub_zero]
  rw [sum_range_sub (fun i ↦ ∑ j ∈ range n, u i j * d2 v i j) n]
  simp only [hg, sub_self]

end Grid

end MultAdd.JetFields
