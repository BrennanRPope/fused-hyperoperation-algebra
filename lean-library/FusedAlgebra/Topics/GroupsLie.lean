/-
  Fused Hyperoperation Algebra: group theory and Lie groups through mult-add.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base.

  Intended library path: FusedAlgebra/Topics/GroupsLie.lean (imports only FusedAlgebra.Base).

  What is here:
    1. The landscape law F_{p,q}(a,b) = (a + b + pab)/(1 − qab) is a product of Möbius
       matrices: mob a = I + a·X with X = [[p, 1], [−q, 0]], and
       mob a · mob b = (1 − qab)·I + (a + b + pab)·X. X obeys X² = pX − qI.
    2. The same thing as numbers 1 + xθ with θ² = pθ − q (the lift for every landscape point).
       The norm x² + pxy + qy² is multiplicative and the norm of 1 + aθ is the rate
       V(a) = 1 + pa + qa². Hence V(a)V(b) = (1 − qab)²·V(F) and ∂F/∂b = V(a)/(1 − qab)²:
       the ruler dx/V(x) is carried by every translation (the invariant differential).
    3. Lie's algebra of the line, fields α + βx + γx²: the bracket closes, is antisymmetric,
       obeys Jacobi, and the form β² − 4αγ (the discriminant, the Killing form up to a factor)
       is ad-invariant. The landscape is the plane α = 1.
    4. Re-coordinatizing the line by x ↦ x/(1 − cx) moves (p, q) to (p + 2c, q + pc + c²) and
       keeps p² − 4q; scaling multiplies it by λ². So three classes: line, border, circle.
    5. The dial inside the ax+b group: T_a(x) = (1 + ka)x + a composes by ◯ₖ, is a dilation
       about the wall −1/k, rescales k under conjugation, and at k = 0 is a fixed-point-free
       translation (not conjugate to any dilation).
    6. The twisting corner on 2×2 matrices: A ◯ B − B ◯ A = AB − BA exactly; brackets are
       trace-free; landscape matrices commute.
    7. Finite clocks: (F_p minus the wall, ◯) is cyclic for p = 5, 7, 11, 13 (by decide), the
       mult-add powers of 1 are Mersenne numbers 2ⁿ − 1, and the circle x² + y² = 1 has
       p + 1 or p − 1 points.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/GroupsLie.lean; two exact twins of ShapeDirections removed (see the LIBRARY comments).
-/
import FusedAlgebra.Base

namespace GroupsLie
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. The landscape law as a Möbius product -/

/-- Numerator and denominator of the landscape law F_{p,q}. -/
def N (p a b : Int) : Int := a + b + p * a * b
def Dn (q a b : Int) : Int := 1 - q * a * b
/-- The rate (Riccati field) V(x) = 1 + px + qx². -/
def V (p q x : Int) : Int := 1 + p * x + q * x * x

/-- 2×2 integer matrices written out, so `grind` does the algebra. -/
structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq

theorem M2.ext' {x y : M2} (h1 : x.a = y.a) (h2 : x.b = y.b) (h3 : x.c = y.c) (h4 : x.d = y.d) :
    x = y := by
  cases x; cases y; simp_all

def M2.mul (x y : M2) : M2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d, x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩
def M2.add (x y : M2) : M2 := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def M2.sub (x y : M2) : M2 := ⟨x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d⟩
def M2.smul (k : Int) (x : M2) : M2 := ⟨k * x.a, k * x.b, k * x.c, k * x.d⟩
def M2.one : M2 := ⟨1, 0, 0, 1⟩
def M2.tr (x : M2) : Int := x.a + x.d

/-- The landscape generator X = [[p, 1], [−q, 0]]. -/
def gen (p q : Int) : M2 := ⟨p, 1, -q, 0⟩
/-- The Möbius matrix of the translation by a: [[1 + pa, a], [−qa, 1]]. -/
def mob (p q a : Int) : M2 := ⟨1 + p * a, a, -(q * a), 1⟩

/-- mob a = I + a·X: each translation is one step along the generator. -/
theorem mob_linear (p q a : Int) : mob p q a = M2.add M2.one (M2.smul a (gen p q)) := by
  apply M2.ext' <;> simp only [mob, M2.add, M2.one, M2.smul, gen] <;> grind

theorem mob_zero (p q : Int) : mob p q 0 = M2.one := by
  apply M2.ext' <;> simp only [mob, M2.one] <;> grind

/-- Cayley–Hamilton for the generator: X² = pX − qI. -/
theorem gen_sq (p q : Int) :
    M2.mul (gen p q) (gen p q) = M2.sub (M2.smul p (gen p q)) (M2.smul q M2.one) := by
  apply M2.ext' <;> simp only [M2.mul, M2.sub, M2.smul, gen, M2.one] <;> grind

/-- The law is matrix multiplication up to the scale 1 − qab:
    mob a · mob b = (1 − qab)·I + (a + b + pab)·X. -/
theorem mob_mul (p q a b : Int) :
    M2.mul (mob p q a) (mob p q b) =
      M2.add (M2.smul (Dn q a b) M2.one) (M2.smul (N p a b) (gen p q)) := by
  apply M2.ext' <;> simp only [M2.mul, M2.add, M2.smul, mob, M2.one, gen, Dn, N] <;> grind

/-- Same statement entry by entry: the product is (1 − qab)·mob(F(a,b)) with F = N/Dn. -/
theorem mob_mul_entries (p q a b : Int) :
    M2.mul (mob p q a) (mob p q b) = ⟨Dn q a b + p * N p a b, N p a b, -(q * N p a b), Dn q a b⟩ := by
  apply M2.ext' <;> simp only [M2.mul, mob, Dn, N] <;> grind

/-- On the dial line (q = 0) there is no scale at all: matrices multiply exactly as ◯ₖ. -/
theorem mob_dial (k a b : Int) : M2.mul (mob k 0 a) (mob k 0 b) = mob k 0 (circK k a b) := by
  apply M2.ext' <;> simp only [M2.mul, mob, circK] <;> grind

/-- Mult-add itself as 2×2 matrices [[1 + a, a], [0, 1]]. -/
theorem mob_circ (a b : Int) : M2.mul (mob 1 0 a) (mob 1 0 b) = mob 1 0 (a ◯ b) := by
  apply M2.ext' <;> simp only [M2.mul, mob, circ] <;> grind

/-- Landscape matrices commute: the corner cannot twist inside one landscape law. -/
theorem mob_comm (p q a b : Int) : M2.mul (mob p q a) (mob p q b) = M2.mul (mob p q b) (mob p q a) := by
  apply M2.ext' <;> simp only [M2.mul, mob] <;> grind

/-- Acting on a point: the Möbius map of mob a sends b to F(a,b) (cross-multiplied). -/
theorem mob_acts (p q a b : Int) :
    (mob p q a).a * b + (mob p q a).b = N p a b ∧ (mob p q a).c * b + (mob p q a).d = Dn q a b := by
  simp only [mob, N, Dn]; constructor <;> grind

/-! ## 2. The lift for every point: 1 + xθ with θ² = pθ − q -/

/-- Multiplication of x + yθ in the algebra θ² = pθ − q, as pairs (x, y). -/
def tmul (p q : Int) (u v : Int × Int) : Int × Int :=
  (u.1 * v.1 - q * u.2 * v.2, u.1 * v.2 + u.2 * v.1 + p * u.2 * v.2)

/-- (1 + aθ)(1 + bθ) = (1 − qab) + (a + b + pab)θ: the landscape law is θ-multiplication. -/
theorem theta_lift (p q a b : Int) : tmul p q (1, a) (1, b) = (Dn q a b, N p a b) := by
  simp only [tmul, Dn, N]; ext <;> simp <;> grind

theorem tmul_comm (p q : Int) (u v : Int × Int) : tmul p q u v = tmul p q v u := by
  simp only [tmul]; ext <;> simp <;> grind

theorem tmul_assoc (p q : Int) (u v w : Int × Int) :
    tmul p q (tmul p q u v) w = tmul p q u (tmul p q v w) := by
  simp only [tmul]; ext <;> simp <;> grind

/-- p = 1, q = 0: θ² = θ, and the first slot of (1 + aθ)(1 + bθ) is the lift of a ◯ b. -/
theorem theta_circ (a b : Int) : tmul 1 0 (1, a) (1, b) = (1, a ◯ b) := by
  simp only [tmul, circ]; ext <;> simp <;> grind

/-- p = q = 0: dual numbers, θ² = 0, plain addition. -/
theorem theta_dual (a b : Int) : tmul 0 0 (1, a) (1, b) = (1, a + b) := by
  simp only [tmul]; ext <;> simp <;> grind

/-- p = 0, q = 1: θ² = −1, the complex numbers; the law is tan addition (a + b)/(1 − ab). -/
theorem theta_tan (a b : Int) : tmul 0 1 (1, a) (1, b) = (1 - a * b, a + b) := by
  simp only [tmul]; ext <;> simp <;> grind

/-- p = 0, q = −1: θ² = 1, split-complex numbers; velocity addition (a + b)/(1 + ab). -/
theorem theta_velocity (a b : Int) : tmul 0 (-1) (1, a) (1, b) = (1 + a * b, a + b) := by
  simp only [tmul]; ext <;> simp <;> grind

/-- The norm x² + pxy + qy² of x + yθ. -/
def tnorm (p q : Int) (u : Int × Int) : Int := u.1 * u.1 + p * u.1 * u.2 + q * u.2 * u.2

theorem tnorm_mul (p q : Int) (u v : Int × Int) :
    tnorm p q (tmul p q u v) = tnorm p q u * tnorm p q v := by
  simp only [tnorm, tmul]; grind

/-- The norm of 1 + aθ is the rate V(a). -/
theorem tnorm_one_theta (p q a : Int) : tnorm p q (1, a) = V p q a := by
  simp only [tnorm, V]; grind

/-- Rate law: V(a)·V(b) = (1 − qab)²·V(F(a,b)), cross-multiplied. -/
theorem rate_mul (p q a b : Int) :
    V p q a * V p q b = Dn q a b * Dn q a b + p * N p a b * Dn q a b + q * N p a b * N p a b := by
  simp only [V, Dn, N]; grind

/-- The b-derivative of F has numerator ∂N·Dn − N·∂Dn = V(a). -/
theorem deriv_numerator (p q a b : Int) : (1 + p * a) * Dn q a b - N p a b * (-(q * a)) = V p q a := by
  simp only [Dn, N, V]; grind

/-- Invariant ruler: ∂F/∂b · V(b) = V(F), cross-multiplied by Dn²
    (so the tick spacing dx/V(x) is carried along by every translation). -/
theorem invariant_ruler (p q a b : Int) :
    ((1 + p * a) * Dn q a b - N p a b * (-(q * a))) * V p q b =
      Dn q a b * Dn q a b + p * N p a b * Dn q a b + q * N p a b * N p a b := by
  rw [deriv_numerator, rate_mul]

/-- The rate splits over its walls: 1 − (u+v)x + uvx² = (1 − ux)(1 − vx). -/
theorem rate_walls (u v x : Int) : V (-(u + v)) (u * v) x = (1 - u * x) * (1 - v * x) := by
  simp only [V]; grind

/-! ## 3. Lie's algebra of the line: fields α + βx + γx² -/

/-- A field α + βx + γx² (times ∂/∂x), stored as its three coefficients. -/
structure Fld where
  al : Int
  be : Int
  ga : Int
deriving DecidableEq

theorem Fld.ext' {f g : Fld} (h1 : f.al = g.al) (h2 : f.be = g.be) (h3 : f.ga = g.ga) : f = g := by
  cases f; cases g; simp_all

/-- The bracket [f, g] = f g' − g f', which stays in the span of 1, x, x². -/
def br (f g : Fld) : Fld :=
  ⟨f.al * g.be - f.be * g.al, 2 * (f.al * g.ga - f.ga * g.al), f.be * g.ga - f.ga * g.be⟩

theorem br_antisymm (f g : Fld) : br f g = ⟨-(br g f).al, -(br g f).be, -(br g f).ga⟩ := by
  apply Fld.ext' <;> simp only [br] <;> grind

theorem br_jacobi (f g h : Fld) :
    let s := br f (br g h); let t := br g (br h f); let u := br h (br f g)
    s.al + t.al + u.al = 0 ∧ s.be + t.be + u.be = 0 ∧ s.ga + t.ga + u.ga = 0 := by
  simp only [br]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- The quadratic form β² − 4αγ (proportional to the Killing form) and its bilinear partner. -/
def Q (f : Fld) : Int := f.be * f.be - 4 * f.al * f.ga
def Bl (f g : Fld) : Int := 2 * f.be * g.be - 4 * (f.al * g.ga + f.ga * g.al)

theorem Bl_diag (f : Fld) : Bl f f = 2 * Q f := by simp only [Bl, Q]; grind

/-- Ad-invariance: B([g,f], h) + B(f, [g,h]) = 0, so the light cone Q = 0 is a Lie-theoretic object. -/
theorem Bl_invariant (f g h : Fld) : Bl (br g f) h + Bl f (br g h) = 0 := by
  simp only [Bl, br]; grind

/-- The landscape point (p, q) is the field 1 + px + qx²; its form is the discriminant. -/
theorem landscape_Q (p q : Int) : Q ⟨1, p, q⟩ = p * p - 4 * q := by simp only [Q]; grind

/-- Translation ∂ and dilation x∂ bracket to a translation: [x∂, ∂] = −∂. -/
theorem br_dilation_translation : br ⟨0, 1, 0⟩ ⟨1, 0, 0⟩ = ⟨-1, 0, 0⟩ := by decide

/-- Two dial generators 1 + kx and 1 + k'x bracket to a pure translation (k' − k)·∂. -/
theorem br_dial (k k' : Int) : br ⟨1, k, 0⟩ ⟨1, k', 0⟩ = ⟨k' - k, 0, 0⟩ := by
  apply Fld.ext' <;> simp only [br] <;> grind

/-! ## 4. Changing the chart: three classes -/

/-- In the chart y = x/(1 − cx) the rate becomes 1 + (p + 2c)y + (q + pc + c²)y². -/
theorem shift_rate (p q c y : Int) :
    (1 + c * y) * (1 + c * y) + p * y * (1 + c * y) + q * y * y = V (p + 2 * c) (q + p * c + c * c) y := by
  simp only [V]; grind

-- LIBRARY: shift_disc (the discriminant does not see the chart) was an exact twin of ShapeDirections.ruler_keeps_disc and is removed.

-- LIBRARY: scale_disc (scaling the line multiplies the discriminant by λ²) was an exact twin of ShapeDirections.zoom_disc and is removed.

/-- The chart change carries one law to the other. Write ψ(y) = y/(1 + cy) for the way back
    (ψ undoes φ(x) = x/(1 − cx)). Then F_{p,q}(ψA, ψB) = ψ(F_{p+2c, q+pc+c²}(A, B)):
    after clearing the factors (1 + cA)(1 + cB), the numerators agree ... -/
theorem shift_conj_num (p c A B : Int) :
    A * (1 + c * B) + B * (1 + c * A) + p * A * B = N (p + 2 * c) A B := by
  simp only [N]; grind

/-- ... and so do the denominators: (1 + cA)(1 + cB) − qAB = D' + c·N'. -/
theorem shift_conj_den (p q c A B : Int) :
    (1 + c * A) * (1 + c * B) - q * A * B = Dn (q + p * c + c * c) A B + c * N (p + 2 * c) A B := by
  simp only [Dn, N]; grind

/-! ## 5. The dial inside the ax+b group -/

/-- Translation by a under the dial: the affine map x ↦ (1 + ka)x + a. -/
def T (k a x : Int) : Int := (1 + k * a) * x + a

theorem T_eq_circK (k a x : Int) : T k a x = circK k a x := by simp only [T, circK]; grind

theorem T_comp (k a b x : Int) : T k a (T k b x) = T k (circK k a b) x := by
  simp only [T, circK]; grind

/-- Scaling conjugation moves the dial: T^k_{λa}(λx) = λ·T^{kλ}_a(x). -/
theorem T_scale (k l a x : Int) : T k (l * a) (l * x) = l * T (k * l) a x := by
  simp only [T]; grind

/-- About the wall −1/k every translation is a pure dilation by 1 + ka. -/
theorem T_dilation (k a y : Rat) (hk : k ≠ 0) :
    (1 + k * a) * (y - 1 / k) + a + 1 / k = (1 + k * a) * y := by grind

/-- At k = 0 a nonzero translation has no fixed point, so it is no dilation. -/
theorem T_zero_no_fixed (a x : Int) (ha : a ≠ 0) : T 0 a x ≠ x := by
  simp only [T]; grind

/-! ## 6. The twisting corner on matrices -/

/-- Matrix mult-add A ◯ B = A + B + AB. -/
def maM (x y : M2) : M2 := M2.add (M2.add x y) (M2.mul x y)

/-- The failure of matrix mult-add to commute is exactly the commutator AB − BA. -/
theorem maM_bracket (x y : M2) : M2.sub (maM x y) (maM y x) = M2.sub (M2.mul x y) (M2.mul y x) := by
  apply M2.ext' <;> simp only [maM, M2.sub, M2.add, M2.mul] <;> grind

/-- Brackets are trace-free: they land in sl₂. -/
theorem bracket_trace (x y : M2) : (M2.sub (M2.mul x y) (M2.mul y x)).tr = 0 := by
  simp only [M2.tr, M2.sub, M2.mul]; grind

/-- On the landscape, a·X and b·X mult-add the same in either order. -/
theorem landscape_no_twist (p q a b : Int) :
    maM (M2.smul a (gen p q)) (M2.smul b (gen p q)) = maM (M2.smul b (gen p q)) (M2.smul a (gen p q)) := by
  apply M2.ext' <;> simp only [maM, M2.add, M2.mul, M2.smul, gen] <;> grind

/-- A 2×2 pair where the twist is visible: [[0,1],[0,0]] and [[0,0],[1,0]]. -/
theorem twist_example :
    M2.sub (maM ⟨0, 1, 0, 0⟩ ⟨0, 0, 1, 0⟩) (maM ⟨0, 0, 1, 0⟩ ⟨0, 1, 0, 0⟩) = ⟨1, 0, 0, -1⟩ := by
  decide

/-! ## 7. Finite clocks -/

/-- n copies of 1 under mult-add are the Mersenne number 2ⁿ − 1. -/
theorem mersenne (n : Nat) : mpow (1 : Int) n = 2 ^ n - 1 := by
  rw [mpow_eq]; rfl

/-- (F₅ minus the wall 4, ◯) is a 4-cycle generated by 1. -/
theorem f5_cyclic : mpow (1 : Fin 5) 4 = 0 ∧ ∀ a : Fin 5, a ≠ 4 → ∃ j : Fin 4, mpow (1 : Fin 5) j.val = a := by
  decide

/-- Mod 7, 1 has order 3 (2 is not a primitive root), but 2 generates all six. -/
theorem f7_cyclic : mpow (1 : Fin 7) 3 = 0 ∧ mpow (2 : Fin 7) 6 = 0 ∧
    ∀ a : Fin 7, a ≠ 6 → ∃ j : Fin 6, mpow (2 : Fin 7) j.val = a := by
  decide

theorem f11_cyclic : mpow (1 : Fin 11) 10 = 0 ∧
    ∀ a : Fin 11, a ≠ 10 → ∃ j : Fin 10, mpow (1 : Fin 11) j.val = a := by
  decide

theorem f13_cyclic : mpow (1 : Fin 13) 12 = 0 ∧
    ∀ a : Fin 13, a ≠ 12 → ∃ j : Fin 12, mpow (1 : Fin 13) j.val = a := by
  decide

/-- Points on the circle x² + y² = 1 over F_p. -/
def circleCount (p : Nat) [NeZero p] : Nat :=
  ((List.finRange p).flatMap fun x => (List.finRange p).filter fun y => x * x + y * y = (1 : Fin p)).length

/-- The finite circle: p + 1 points when p ≡ 3 mod 4, p − 1 when p ≡ 1 mod 4. -/
theorem circle_counts : circleCount 3 = 4 ∧ circleCount 5 = 4 ∧ circleCount 7 = 8 ∧
    circleCount 11 = 12 ∧ circleCount 13 = 12 := by decide

/-- The circle over F₇ is closed under complex multiplication (the circle group). -/
theorem circle7_closed : ∀ a b c d : Fin 7, a * a + b * b = 1 → c * c + d * d = 1 →
    (a * c - b * d) * (a * c - b * d) + (a * d + b * c) * (a * d + b * c) = 1 := by decide

end GroupsLie
