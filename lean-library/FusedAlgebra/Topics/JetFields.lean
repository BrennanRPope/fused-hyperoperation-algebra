/-
  JetFields: fields of jets, mult-add, and the Navier–Stokes velocity gradient
  ===========================================================================

  Core Lean 4.37, no mathlib, no `sorry`. Imports the shared base.
  Written Nov 25, 2026 for the page "Fields of jets in mult-add".

  LIBRARY NOTE: copied from jet-fields/lean/JetFields.lean (Nov 25), unchanged
  apart from this note.

  The picture.

  * A jet is a point of a truncated ring. The 1-jet of f at x lives in
    R[ε]/(ε² − kε). With k = 0 these are the dual numbers (the derivative);
    with k = h they hold the values at two points h apart (the forward
    difference). Units 1 + aε multiply by the weighted mult-add:
    (1 + aε)(1 + bε) = 1 + (a ⊞ₖ b)ε (`units_circK`). So the corner weight k
    is the gap between two sample points, and k → 0 merges them into a
    tangent vector (`twoPt_mul`, `twoPt_weight_zero`). The wall 1 + ka = 0 is
    where the unit vanishes at the second point (`wall_no_inverse`); the dual
    numbers have no wall (`dual_inverse`).
  * In the 2-jet ring the relative jets ν (no constant term) mult-add with
    the cube cut after layer 2: the triple corner is gone (`cube_depth_two`),
    and the way back is a finite sum (`nilpotent_inverse`).
  * The Lagrangian step x ↦ x + h·u(x) has volume change det(I + h∇u). Its
    expansion in h is the mult-add cube of the gradient's eigen-rates:
    1 + hP + h²Q + h³R (`det_one_add_smul`, `det_diag_cube`), where
    P = trace (divergence), Q = sum of principal 2×2 minors, R = det.
    These are the P, Q, R invariants that turbulence research plots.
  * For an incompressible gradient (P = 0), with T = ∇u + ∇uᵀ (twice the
    strain) and ω the vorticity:  8Q = 2|ω|² − |T|²  (`e2_vorticity_strain`)
    and 8R = det T + ω·Tω  (`det_vorticity_strain`). Layer 2 is rotation
    against strain; layer 3, the far corner, is strain self-amplification
    plus vortex stretching.
  * Restricted Euler (Vieillefosse 1982, Cantwell 1992): with
    3·dA/dt = −3A² + tr(A²)·I, the layers move by 3Q' = 9R, 3R' = −2Q²
    (`re_dQ`, `re_dR`), and 27R² + 4Q³ is conserved (`re_conserved`). In 2D
    the same flow is frozen: −2A² + tr(A²)·I = 0 (`re2_frozen`), since the
    square has no third layer (`tr_cube_2d`).
  * Every minor of a gradient is a divergence (a null Lagrangian). On a
    periodic grid the 2×2 minor, read with the corner kept (one index
    shifted), is an exact discrete divergence (`minor_is_divergence`) and so
    sums to zero over the box (`box_minor_zero`). The shift is the corner:
    S·Δv = Δv + Δ²v (`shift_is_corner`).

  Not here (needs real analysis or integrals): the continuous null-Lagrangian
  integrals ⟨Q⟩ = ⟨R⟩ = 0 and Betchov's relation, the blow-up time of the
  Vieillefosse tail, anything about actual Navier–Stokes solutions.

  STATUS: compiled with core Lean 4.37.0-rc1 against the shared library
  (Base), standard axioms only. Library-ready: copy to
  lean-library/FusedAlgebra/Topics/JetFields.lean, add
  FusedAlgebra.Topics.JetFields to build.sh and FusedAlgebra.lean.
-/
import FusedAlgebra.Base

namespace CircleOp.JetFields

open CircleOp
open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 1. Weighted dual numbers: R[ε]/(ε² − kε) -/

/-- An element re + du·ε of R[ε]/(ε² − kε): a value and a slope. -/
structure WDual (R : Type) where
  re : R
  du : R

section WDualRing
variable {R : Type} [CommRing R]

/-- Multiplication in R[ε]/(ε² − kε). -/
def wmul (k : R) (x y : WDual R) : WDual R :=
  ⟨x.re * y.re, x.re * y.du + x.du * y.re + k * x.du * y.du⟩

theorem wmul_comm (k : R) (x y : WDual R) : wmul k x y = wmul k y x := by
  obtain ⟨a, b⟩ := x; obtain ⟨c, d⟩ := y
  simp only [wmul, WDual.mk.injEq]; constructor <;> grind

theorem wmul_assoc (k : R) (x y z : WDual R) :
    wmul k (wmul k x y) z = wmul k x (wmul k y z) := by
  obtain ⟨a, b⟩ := x; obtain ⟨c, d⟩ := y; obtain ⟨e, f⟩ := z
  simp only [wmul, WDual.mk.injEq]; constructor <;> grind

/-- ε² = kε: the corner weight is the square of the infinitesimal. -/
theorem eps_sq (k : R) : wmul k ⟨0, 1⟩ ⟨0, 1⟩ = ⟨0, k⟩ := by
  simp only [wmul, WDual.mk.injEq]; constructor <;> grind

/-- **Units multiply by weighted mult-add**: (1 + aε)(1 + bε) = 1 + (a ⊞ₖ b)ε. -/
theorem units_circK (k a b : R) : wmul k ⟨1, a⟩ ⟨1, b⟩ = ⟨1, circK k a b⟩ := by
  simp only [wmul, circK, WDual.mk.injEq]; constructor <;> grind

/-- With weight 1 the units multiply by plain mult-add. -/
theorem units_circ (a b : R) : wmul 1 ⟨1, a⟩ ⟨1, b⟩ = ⟨1, circ a b⟩ := by
  rw [units_circK, circK_one]

/-- With weight 0 (dual numbers, 1-jets) the corner is dropped: relative
    slopes add. This is the logarithmic derivative of a product. -/
theorem dual_units_add (a b : R) : wmul 0 ⟨1, a⟩ ⟨1, b⟩ = ⟨1, a + b⟩ := by
  rw [units_circK, circK_weight_zero]

/-- The two-point reading: value at x and value at x + k (slope · gap added). -/
def twoPt (k : R) (x : WDual R) : R × R := (x.re, x.re + k * x.du)

/-- **The weighted dual numbers are two points.** Reading an element at x and at
    x + k turns the twisted product into the plain product at each point. -/
theorem twoPt_mul (k : R) (x y : WDual R) :
    twoPt k (wmul k x y) = ((twoPt k x).1 * (twoPt k y).1, (twoPt k x).2 * (twoPt k y).2) := by
  obtain ⟨a, b⟩ := x; obtain ⟨c, d⟩ := y
  simp only [twoPt, wmul, Prod.mk.injEq]; constructor <;> grind

/-- At weight 0 the two points coincide: the secant has become the tangent. -/
theorem twoPt_weight_zero (x : WDual R) : twoPt 0 x = (x.re, x.re) := by
  simp only [twoPt, Prod.mk.injEq]; constructor <;> grind

/-- **Discrete Leibniz with the corner.** If f(x+h) = f₀ + h·f' and
    g(x+h) = g₀ + h·g', then (fg)(x+h) = f₀g₀ + h·(the ε-part of the product in
    R[ε]/(ε² − hε)). The forward difference of a product is that ε-part. -/
theorem diff_leibniz (h f₀ f' g₀ g' : R) :
    (f₀ + h * f') * (g₀ + h * g') = f₀ * g₀ + h * (wmul h ⟨f₀, f'⟩ ⟨g₀, g'⟩).du := by
  simp only [wmul]; grind

/-- The dual numbers have no wall: 1 + aε always has the inverse 1 − aε. -/
theorem dual_inverse (a : R) : wmul 0 ⟨1, a⟩ ⟨1, -a⟩ = ⟨1, 0⟩ := by
  simp only [wmul, WDual.mk.injEq]; constructor <;> grind

end WDualRing

section WDualField
variable {K : Type} [Field K]

/-- From two point values back to value and slope (needs a nonzero gap). -/
def ofTwoPt (k : K) (p : K × K) : WDual K := ⟨p.1, (p.2 - p.1) / k⟩

theorem twoPt_ofTwoPt (k : K) (hk : k ≠ 0) (p : K × K) : twoPt k (ofTwoPt k p) = p := by
  obtain ⟨a, b⟩ := p
  simp only [twoPt, ofTwoPt, Prod.mk.injEq]
  refine ⟨trivial, ?_⟩
  grind

theorem ofTwoPt_twoPt (k : K) (hk : k ≠ 0) (x : WDual K) : ofTwoPt k (twoPt k x) = x := by
  obtain ⟨a, b⟩ := x
  simp only [twoPt, ofTwoPt, WDual.mk.injEq]
  refine ⟨trivial, ?_⟩
  grind

/-- **The wall is a zero at the second point.** If 1 + ka = 0, the unit 1 + aε
    has no inverse in R[ε]/(ε² − kε): read at x + k, it is 0. -/
theorem wall_no_inverse (k a : K) (hw : 1 + k * a = 0) (y : WDual K) :
    wmul k ⟨1, a⟩ y ≠ ⟨1, 0⟩ := by
  intro h
  have h2 := congrArg (fun z => z.re + k * z.du) h
  simp only [wmul] at h2
  have e : 1 * y.re + k * (1 * y.du + a * y.re + k * a * y.du)
      = (1 + k * a) * (y.re + k * y.du) := by grind
  rw [e, hw] at h2
  grind

/-- Off the wall, the inverse is the weighted mult-add way back −a/(1 + ka). -/
theorem off_wall_inverse (k a : K) (hw : 1 + k * a ≠ 0) :
    wmul k ⟨1, a⟩ ⟨1, -a / (1 + k * a)⟩ = ⟨1, 0⟩ := by
  simp only [wmul, WDual.mk.injEq]
  constructor <;> grind

end WDualField

/-! ## 2. Second-order jets: the cube cut after layer 2 -/

/-- A 2-jet as Taylor coefficients (f, f', f''/2): an element of R[ε]/(ε³). -/
structure Jet2 (R : Type) where
  c0 : R
  c1 : R
  c2 : R

section Jet2Ring
variable {R : Type} [CommRing R]

/-- Multiplication in R[ε]/(ε³). -/
def jmul (x y : Jet2 R) : Jet2 R :=
  ⟨x.c0 * y.c0, x.c0 * y.c1 + x.c1 * y.c0, x.c0 * y.c2 + x.c1 * y.c1 + x.c2 * y.c0⟩

/-- Mult-add in R[ε]/(ε³): x + y + xy. -/
def jma (x y : Jet2 R) : Jet2 R :=
  ⟨x.c0 + y.c0 + x.c0 * y.c0, x.c1 + y.c1 + (x.c0 * y.c1 + x.c1 * y.c0),
   x.c2 + y.c2 + (x.c0 * y.c2 + x.c1 * y.c1 + x.c2 * y.c0)⟩

/-- The jet product is the product of functions: second-order Leibniz,
    (fg)'' = f''g + 2f'g' + fg'' (with f'' = 2·c2). -/
theorem leibniz2 (x y : Jet2 R) :
    2 * (jmul x y).c2 = x.c0 * (2 * y.c2) + 2 * (x.c1 * y.c1) + (2 * x.c2) * y.c0 := by
  simp only [jmul]; grind

/-- Mult-add of jets is multiplication through 1 + x, as everywhere else. -/
theorem jma_lift (x y : Jet2 R) :
    jmul ⟨1 + x.c0, x.c1, x.c2⟩ ⟨1 + y.c0, y.c1, y.c2⟩
      = ⟨1 + (jma x y).c0, (jma x y).c1, (jma x y).c2⟩ := by
  obtain ⟨a, b, c⟩ := x; obtain ⟨d, e, f⟩ := y
  simp only [jmul, jma, Jet2.mk.injEq]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- **Jet order = cube depth.** Three relative 2-jets (no constant term)
    mult-add to their sum plus the three pair corners; the triple corner, which
    would sit at ε³, is gone. -/
theorem cube_depth_two (a₁ b₁ a₂ b₂ a₃ b₃ : R) :
    jma ⟨0, a₁, b₁⟩ (jma ⟨0, a₂, b₂⟩ ⟨0, a₃, b₃⟩)
      = ⟨0, a₁ + a₂ + a₃, (b₁ + b₂ + b₃) + (a₁ * a₂ + a₁ * a₃ + a₂ * a₃)⟩ := by
  simp only [jma, Jet2.mk.injEq]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- At order 1 only layer 1 survives: relative jets just add. -/
theorem cube_depth_one (a₁ a₂ a₃ : R) :
    wmul 0 (wmul 0 ⟨1, a₁⟩ ⟨1, a₂⟩) ⟨1, a₃⟩ = ⟨1, a₁ + a₂ + a₃⟩ := by
  rw [dual_units_add, dual_units_add]

/-- **No wall in jet space.** A relative jet ν has the finite way back
    −ν + ν² (the geometric series stops because ε³ = 0). -/
theorem nilpotent_inverse (a b : R) : jma ⟨0, a, b⟩ ⟨0, -a, -b + a * a⟩ = ⟨0, 0, 0⟩ := by
  simp only [jma, Jet2.mk.injEq]; refine ⟨?_, ?_, ?_⟩ <;> grind

end Jet2Ring

/-! ## 3. The velocity gradient: the mult-add cube of the step x ↦ x + h·u -/

/-- A 3×3 matrix, entry mᵢⱼ. For the velocity gradient, mᵢⱼ = ∂uᵢ/∂xⱼ. -/
structure M3 (R : Type) where
  (m11 m12 m13 m21 m22 m23 m31 m32 m33 : R)

section M3Ring
variable {R : Type} [CommRing R]

def M3.mul (X Y : M3 R) : M3 R :=
  ⟨X.m11*Y.m11 + X.m12*Y.m21 + X.m13*Y.m31, X.m11*Y.m12 + X.m12*Y.m22 + X.m13*Y.m32,
   X.m11*Y.m13 + X.m12*Y.m23 + X.m13*Y.m33,
   X.m21*Y.m11 + X.m22*Y.m21 + X.m23*Y.m31, X.m21*Y.m12 + X.m22*Y.m22 + X.m23*Y.m32,
   X.m21*Y.m13 + X.m22*Y.m23 + X.m23*Y.m33,
   X.m31*Y.m11 + X.m32*Y.m21 + X.m33*Y.m31, X.m31*Y.m12 + X.m32*Y.m22 + X.m33*Y.m32,
   X.m31*Y.m13 + X.m32*Y.m23 + X.m33*Y.m33⟩

/-- I + h·X: the gradient of the step map x ↦ x + h·u(x). -/
def M3.step (h : R) (X : M3 R) : M3 R :=
  ⟨1 + h*X.m11, h*X.m12, h*X.m13, h*X.m21, 1 + h*X.m22, h*X.m23, h*X.m31, h*X.m32, 1 + h*X.m33⟩

/-- Layer 1, P: the trace (for ∇u, the divergence). -/
def M3.tr (X : M3 R) : R := X.m11 + X.m22 + X.m33

/-- Layer 2, Q: the sum of the three principal 2×2 minors. -/
def M3.e2 (X : M3 R) : R :=
  (X.m11*X.m22 - X.m12*X.m21) + (X.m11*X.m33 - X.m13*X.m31) + (X.m22*X.m33 - X.m23*X.m32)

/-- The determinant of three rows. -/
def det9 (a b c d e f g h i : R) : R := a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)

/-- Layer 3, R: the determinant. -/
def M3.det (X : M3 R) : R :=
  det9 X.m11 X.m12 X.m13 X.m21 X.m22 X.m23 X.m31 X.m32 X.m33

/-- **The volume change of a step is the mult-add cube.**
    det(I + h∇u) = 1 + h·P + h²·Q + h³·R. -/
theorem det_one_add_smul (h : R) (X : M3 R) :
    (X.step h).det = 1 + h * X.tr + h^2 * X.e2 + h^3 * X.det := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.step, M3.det, det9, M3.tr, M3.e2]; grind

/-- A diagonal gradient (principal axes): stretch rates a, b, c. -/
def diag3 (a b c : R) : M3 R := ⟨a, 0, 0, 0, b, 0, 0, 0, c⟩

/-- **In principal axes the step's volume change is literally mult-add:**
    det(I + h·diag(a,b,c)) = 1 + (ha ⊞ hb ⊞ hc). -/
theorem det_diag_cube (h a b c : R) :
    ((diag3 a b c).step h).det = 1 + circ (circ (h * a) (h * b)) (h * c) := by
  simp only [diag3, M3.step, M3.det, det9, circ]; grind

/-- The layers of a diagonal gradient are the cube's layers of (a, b, c). -/
theorem diag_layers (a b c : R) :
    (diag3 a b c).tr = a + b + c ∧ (diag3 a b c).e2 = a*b + a*c + b*c ∧
    (diag3 a b c).det = a*b*c := by
  simp only [diag3, M3.tr, M3.e2, M3.det, det9]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- Cayley–Hamilton, entrywise: X³ − P·X² + Q·X − R·I = 0. -/
theorem cayley_hamilton (X : M3 R) :
    let X2 := X.mul X; let X3 := X.mul X2
    (X3.m11 - X.tr*X2.m11 + X.e2*X.m11 - X.det = 0 ∧ X3.m12 - X.tr*X2.m12 + X.e2*X.m12 = 0 ∧
     X3.m13 - X.tr*X2.m13 + X.e2*X.m13 = 0 ∧ X3.m21 - X.tr*X2.m21 + X.e2*X.m21 = 0 ∧
     X3.m22 - X.tr*X2.m22 + X.e2*X.m22 - X.det = 0 ∧ X3.m23 - X.tr*X2.m23 + X.e2*X.m23 = 0 ∧
     X3.m31 - X.tr*X2.m31 + X.e2*X.m31 = 0 ∧ X3.m32 - X.tr*X2.m32 + X.e2*X.m32 = 0 ∧
     X3.m33 - X.tr*X2.m33 + X.e2*X.m33 - X.det = 0) := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.mul, M3.tr, M3.e2, M3.det, det9]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> grind

/-- For an incompressible gradient, tr(X²) = −2Q. -/
theorem tr_sq_traceless (X : M3 R) (h0 : X.tr = 0) : (X.mul X).tr = -2 * X.e2 := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.mul, M3.tr, M3.e2] at h0 ⊢
  have : i = -a - e := by grind
  subst this; grind

/-- For an incompressible gradient, tr(X³) = 3R. -/
theorem tr_cube_traceless (X : M3 R) (h0 : X.tr = 0) : (X.mul (X.mul X)).tr = 3 * X.det := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.mul, M3.tr, M3.det, det9] at h0 ⊢
  have : i = -a - e := by grind
  subst this; grind

/-- The vorticity vector ω = curl u read off the gradient. -/
def M3.vort (X : M3 R) : R × R × R := (X.m32 - X.m23, X.m13 - X.m31, X.m21 - X.m12)

/-- T = X + Xᵀ, twice the strain-rate tensor. -/
def M3.symm2 (X : M3 R) : M3 R :=
  ⟨2*X.m11, X.m12 + X.m21, X.m13 + X.m31, X.m21 + X.m12, 2*X.m22, X.m23 + X.m32,
   X.m31 + X.m13, X.m32 + X.m23, 2*X.m33⟩

/-- |ω|². -/
def vsq (w : R × R × R) : R := w.1*w.1 + w.2.1*w.2.1 + w.2.2*w.2.2

/-- |T|², the sum of the squares of all nine entries. -/
def M3.fsq (X : M3 R) : R :=
  X.m11*X.m11 + X.m12*X.m12 + X.m13*X.m13 + X.m21*X.m21 + X.m22*X.m22 + X.m23*X.m23 +
  X.m31*X.m31 + X.m32*X.m32 + X.m33*X.m33

/-- ω·Tω. -/
def quad (X : M3 R) (w : R × R × R) : R :=
  w.1*(X.m11*w.1 + X.m12*w.2.1 + X.m13*w.2.2) + w.2.1*(X.m21*w.1 + X.m22*w.2.1 + X.m23*w.2.2) +
  w.2.2*(X.m31*w.1 + X.m32*w.2.1 + X.m33*w.2.2)

/-- **Layer 2 is rotation against strain.** For P = 0:
    8Q = 2|ω|² − |T|², that is, Q = ¼|ω|² − ½|S|². -/
theorem e2_vorticity_strain (X : M3 R) (h0 : X.tr = 0) :
    8 * X.e2 = 2 * vsq X.vort - X.symm2.fsq := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.tr] at h0
  have : i = -a - e := by grind
  subst this
  simp only [M3.e2, vsq, M3.vort, M3.symm2, M3.fsq]; grind

/-- **Layer 3, the far corner, is strain self-amplification plus vortex
    stretching.** For P = 0: 8R = det T + ω·Tω, that is,
    R = det S + ¼ ω·Sω. -/
theorem det_vorticity_strain (X : M3 R) (h0 : X.tr = 0) :
    8 * X.det = X.symm2.det + quad X.symm2 X.vort := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.tr] at h0
  have : i = -a - e := by grind
  subst this
  simp only [M3.det, det9, quad, M3.vort, M3.symm2]; grind

/-- For a symmetric traceless S (a pure strain), tr(S³) = 3·det S. -/
theorem sym_tr_cube (a b c e f : R) :
    let S : M3 R := ⟨a, b, c, b, e, f, c, f, -a - e⟩
    (S.mul (S.mul S)).tr = 3 * S.det := by
  simp only [M3.mul, M3.tr, M3.det, det9]; grind

/-- For a pure traceless strain, −2Q is the sum of the squares of its entries. -/
theorem sym_e2_sumsq (a b c e f : R) :
    let S : M3 R := ⟨a, b, c, b, e, f, c, f, -a - e⟩
    (-2) * S.e2 = S.fsq := by
  simp only [M3.e2, M3.fsq]; grind

/-- A square is never negative (core Lean has no ready-made lemma for it). -/
theorem sq_nonneg_int (a : Int) : 0 ≤ a * a := by
  rcases Int.le_total 0 a with h | h
  · exact Int.mul_nonneg h h
  · have := Int.mul_nonneg (Int.neg_nonneg_of_nonpos h) (Int.neg_nonneg_of_nonpos h)
    rwa [Int.neg_mul_neg] at this

/-- So on the integers (and likewise any ordered ring) Q ≤ 0 for a pure strain:
    the strain's second layer only ever dissipates. -/
theorem sym_e2_nonpos (a b c e f : Int) :
    (⟨a, b, c, b, e, f, c, f, -a - e⟩ : M3 Int).e2 ≤ 0 := by
  have h := sym_e2_sumsq a b c e f
  simp only [M3.fsq] at h
  have h1 := sq_nonneg_int a; have h2 := sq_nonneg_int b
  have h3 := sq_nonneg_int c; have h4 := sq_nonneg_int e
  have h5 := sq_nonneg_int f; have h6 := sq_nonneg_int (-a - e)
  generalize (⟨a, b, c, b, e, f, c, f, -a - e⟩ : M3 Int).e2 = q at h ⊢
  generalize a * a = p1 at h h1; generalize b * b = p2 at h h2
  generalize c * c = p3 at h h3; generalize e * e = p4 at h h4
  generalize f * f = p5 at h h5; generalize (-a - e) * (-a - e) = p6 at h h6
  omega

end M3Ring

/-! ## 4. Restricted Euler: the layers drive themselves -/

section RestrictedEuler
variable {R : Type} [CommRing R]

/-- Linear part of Q at X in direction Y: Q(X + tY) = Q(X) + t·dE2 X Y + t²·Q(Y). -/
def dE2 (X Y : M3 R) : R :=
  (X.m11*Y.m22 + Y.m11*X.m22 - X.m12*Y.m21 - Y.m12*X.m21) +
  (X.m11*Y.m33 + Y.m11*X.m33 - X.m13*Y.m31 - Y.m13*X.m31) +
  (X.m22*Y.m33 + Y.m22*X.m33 - X.m23*Y.m32 - Y.m23*X.m32)

/-- Linear part of det at X in direction Y (one row from Y at a time). -/
def dDet (X Y : M3 R) : R :=
  det9 Y.m11 Y.m12 Y.m13 X.m21 X.m22 X.m23 X.m31 X.m32 X.m33 +
  det9 X.m11 X.m12 X.m13 Y.m21 Y.m22 Y.m23 X.m31 X.m32 X.m33 +
  det9 X.m11 X.m12 X.m13 X.m21 X.m22 X.m23 Y.m31 Y.m32 Y.m33

def M3.addS (X Y : M3 R) (t : R) : M3 R :=
  ⟨X.m11 + t*Y.m11, X.m12 + t*Y.m12, X.m13 + t*Y.m13, X.m21 + t*Y.m21, X.m22 + t*Y.m22,
   X.m23 + t*Y.m23, X.m31 + t*Y.m31, X.m32 + t*Y.m32, X.m33 + t*Y.m33⟩

/-- dE2 really is the first-order part of Q. -/
theorem e2_expand (X Y : M3 R) (t : R) :
    (X.addS Y t).e2 = X.e2 + t * dE2 X Y + t^2 * Y.e2 := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X; obtain ⟨a', b', c', d', e', f', g', k', i'⟩ := Y
  simp only [M3.addS, M3.e2, dE2]; grind

/-- dDet really is the first-order part of det (the second-order part is dDet Y X). -/
theorem det_expand (X Y : M3 R) (t : R) :
    (X.addS Y t).det = X.det + t * dDet X Y + t^2 * dDet Y X + t^3 * Y.det := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X; obtain ⟨a', b', c', d', e', f', g', k', i'⟩ := Y
  simp only [M3.addS, M3.det, dDet, det9]; grind

/-- Three times the restricted Euler velocity: 3·dA/dt = −3A² + tr(A²)·I
    (self-advection plus the isotropic part of the pressure Hessian). -/
def reStep (X : M3 R) : M3 R :=
  let X2 := X.mul X
  ⟨-3*X2.m11 + X2.tr, -3*X2.m12, -3*X2.m13, -3*X2.m21, -3*X2.m22 + X2.tr, -3*X2.m23,
   -3*X2.m31, -3*X2.m32, -3*X2.m33 + X2.tr⟩

/-- The flow keeps the gradient incompressible. -/
theorem re_traceless (X : M3 R) : (reStep X).tr = 0 := by
  simp only [reStep, M3.tr]; grind

/-- **Layer 2 is driven by layer 3**: 3·dQ/dt = 9R (so dQ/dt = 3R). -/
theorem re_dQ (X : M3 R) (h0 : X.tr = 0) : dE2 X (reStep X) = 9 * X.det := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.tr] at h0
  have : i = -a - e := by grind
  subst this
  simp only [dE2, reStep, M3.mul, M3.tr, M3.det, det9]; grind

/-- **Layer 3 is driven by layer 2 squared**: 3·dR/dt = −2Q² (so dR/dt = −⅔Q²). -/
theorem re_dR (X : M3 R) (h0 : X.tr = 0) : dDet X (reStep X) = -2 * X.e2^2 := by
  obtain ⟨a, b, c, d, e, f, g, k, i⟩ := X
  simp only [M3.tr] at h0
  have : i = -a - e := by grind
  subst this
  simp only [dDet, reStep, M3.mul, M3.tr, M3.e2, det9]; grind

/-- **The discriminant is conserved**: along restricted Euler,
    d/dt (27R² + 4Q³) = 54R·R' + 12Q²·Q' = 0. 27R² + 4Q³ is minus the
    discriminant of the cube polynomial λ³ + Qλ − R; its zero set is the
    Vieillefosse tail. -/
theorem re_conserved (X : M3 R) (h0 : X.tr = 0) :
    54 * X.det * dDet X (reStep X) + 12 * X.e2^2 * dE2 X (reStep X) = 0 := by
  rw [re_dQ X h0, re_dR X h0]; grind

end RestrictedEuler

/-! ## 5. Two dimensions: the square has no far corner -/

section TwoD
variable {R : Type} [CommRing R]

structure M2 (R : Type) where
  (m11 m12 m21 m22 : R)

def M2.mul (X Y : M2 R) : M2 R :=
  ⟨X.m11*Y.m11 + X.m12*Y.m21, X.m11*Y.m12 + X.m12*Y.m22,
   X.m21*Y.m11 + X.m22*Y.m21, X.m21*Y.m12 + X.m22*Y.m22⟩

def M2.tr (X : M2 R) : R := X.m11 + X.m22
def M2.det (X : M2 R) : R := X.m11*X.m22 - X.m12*X.m21

/-- The 2D step: det(I + hX) = 1 + h·P + h²·Q, the square (1, 2, 1). -/
theorem det2_step (h : R) (X : M2 R) :
    (1 + h*X.m11)*(1 + h*X.m22) - (h*X.m12)*(h*X.m21) = 1 + h * X.tr + h^2 * X.det := by
  simp only [M2.tr, M2.det]; grind

/-- In 2D an incompressible gradient has no third invariant: tr(X³) = 0. -/
theorem tr_cube_2d (a b c : R) :
    ((⟨a, b, c, -a⟩ : M2 R).mul ((⟨a, b, c, -a⟩ : M2 R).mul ⟨a, b, c, -a⟩)).tr = 0 := by
  simp only [M2.mul, M2.tr]; grind

/-- **2D restricted Euler is frozen**: −2A² + tr(A²)·I = 0 for every traceless
    2×2 gradient, so the self-advection is cancelled exactly by the pressure. -/
theorem re2_frozen (a b c : R) :
    let X : M2 R := ⟨a, b, c, -a⟩
    let X2 := X.mul X
    (-2)*X2.m11 + X2.tr = 0 ∧ (-2)*X2.m12 = 0 ∧ (-2)*X2.m21 = 0 ∧ (-2)*X2.m22 + X2.tr = 0 := by
  simp only [M2.mul, M2.tr]; refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

end TwoD

/-! ## 6. Null Lagrangians on a periodic grid: the corner keeps the minor a divergence -/

section Grid
variable {R : Type} [CommRing R]

/-- Forward differences in the two grid directions. -/
def d1 (u : Nat → Nat → R) (i j : Nat) : R := u (i + 1) j - u i j
def d2 (u : Nat → Nat → R) (i j : Nat) : R := u i (j + 1) - u i j

/-- The discrete 2×2 minor of the map (u, v), read with the corner kept:
    each second factor is taken one step along the first factor's direction. -/
def minor (u v : Nat → Nat → R) (i j : Nat) : R :=
  d1 u i j * d2 v (i + 1) j - d2 u i j * d1 v i (j + 1)

/-- The shift is the corner: S₁Δ₂v = Δ₂v + Δ₁Δ₂v (shift = 1 + Δ). -/
theorem shift_is_corner (v : Nat → Nat → R) (i j : Nat) :
    d2 v (i + 1) j = d2 v i j + d1 (d2 v) i j := by
  simp only [d1, d2]; grind

/-- **The corner-kept minor is an exact discrete divergence:**
    minor = Δ₁(u·Δ₂v) − Δ₂(u·Δ₁v). -/
theorem minor_is_divergence (u v : Nat → Nat → R) (i j : Nat) :
    minor u v i j = d1 (fun i j => u i j * d2 v i j) i j - d2 (fun i j => u i j * d1 v i j) i j := by
  simp only [minor, d1, d2]; grind

/-- Finite sum f 0 + ⋯ + f (n−1). -/
def sumTo (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => sumTo f n + f n

theorem sumTo_sub (f g : Nat → R) (n : Nat) :
    sumTo (fun k => f k - g k) n = sumTo f n - sumTo g n := by
  induction n with
  | zero => simp only [sumTo]; grind
  | succ n ih => simp only [sumTo]; rw [ih]; grind

theorem sumTo_telescope (G : Nat → R) (n : Nat) :
    sumTo (fun k => G (k + 1) - G k) n = G n - G 0 := by
  induction n with
  | zero => simp only [sumTo]; grind
  | succ n ih => simp only [sumTo]; rw [ih]; grind

theorem sumTo_congr (f g : Nat → R) (n : Nat) (h : ∀ k, f k = g k) : sumTo f n = sumTo g n := by
  have : f = g := funext h
  rw [this]

/-- **The minor sums to zero over a periodic box**: the grid version of
    "∫ (layer 2 of ∇u) = 0" for a 2D map of the torus. -/
theorem box_minor_zero (u v : Nat → Nat → R) (n : Nat)
    (pu1 : ∀ i j, u (i + n) j = u i j) (pu2 : ∀ i j, u i (j + n) = u i j)
    (pv1 : ∀ i j, v (i + n) j = v i j) (pv2 : ∀ i j, v i (j + n) = v i j) :
    sumTo (fun i => sumTo (fun j => minor u v i j) n) n = 0 := by
  let g : Nat → Nat → R := fun i j => u i j * d2 v i j
  let g' : Nat → Nat → R := fun i j => u i j * d1 v i j
  have inner : ∀ i, sumTo (fun j => minor u v i j) n
      = (sumTo (fun j => g (i + 1) j) n - sumTo (fun j => g i j) n) - (g' i n - g' i 0) := by
    intro i
    have e1 : sumTo (fun j => minor u v i j) n
        = sumTo (fun j => d1 g i j - d2 g' i j) n :=
      sumTo_congr _ _ n (fun j => minor_is_divergence u v i j)
    rw [e1, sumTo_sub]
    have e2 : sumTo (fun j => d1 g i j) n = sumTo (fun j => g (i + 1) j) n - sumTo (fun j => g i j) n := by
      simp only [d1]; rw [sumTo_sub]
    have e3 : sumTo (fun j => d2 g' i j) n = g' i n - g' i 0 := by
      simp only [d2]; rw [sumTo_telescope (fun j => g' i j)]
    rw [e2, e3]
  have wrapj : ∀ i, g' i n - g' i 0 = 0 := by
    intro i
    have h1 := pu2 i 0; have h2 := pv2 (i + 1) 0; have h3 := pv2 i 0
    simp only [Nat.zero_add] at h1 h2 h3
    simp only [g', d1]; rw [h1, h2, h3]; grind
  have e4 : sumTo (fun i => sumTo (fun j => minor u v i j) n) n
      = sumTo (fun i => sumTo (fun j => g (i + 1) j) n - sumTo (fun j => g i j) n) n := by
    apply sumTo_congr; intro i; rw [inner i, wrapj i]; grind
  rw [e4, sumTo_telescope (fun i => sumTo (fun j => g i j) n)]
  have wrapi : sumTo (fun j => g n j) n = sumTo (fun j => g 0 j) n := by
    apply sumTo_congr; intro j
    have h1 := pu1 0 j; have h2 := pv1 0 (j + 1); have h3 := pv1 0 j
    simp only [Nat.zero_add] at h1 h2 h3
    simp only [g, d2]; rw [h1, h2, h3]
  rw [wrapi]; grind

end Grid

end CircleOp.JetFields
