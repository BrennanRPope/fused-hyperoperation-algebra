/-
  Fused Hyperoperation Algebra: complex analysis through mult-add.
  Core Lean 4 only (no mathlib), no `sorry`. Imports the shared library's Base.

  Core Lean has no complex numbers, so a complex number is written out as a record
  ⟨re, im⟩ over any commutative ring (Int gives the Gaussian integers), and the cube
  roots of unity as Eisenstein numbers a + bω with ω² = −1 − ω.

    1. Mult-add on the complex plane is a spin and a stretch about the wall −1:
       1 + (z ◯ w) = (1 + w)(1 + z), the step z ↦ z ◯ w is the vector w turned by 1 + z,
       squared distance from the wall multiplies, the circle |1 + z| = 1 is closed.
    2. Affine maps that fix −1 are exactly the maps z ↦ z ◯ c.
    3. Mult-add roots of 0: if ζⁿ = 1 then (ζ − 1)^(◯n) = 0. Worked: the square n = 4 in
       the Gaussian integers and the triangle n = 3 in the Eisenstein integers.
    4. Möbius: the landscape law F(a,b) = (a + b + p·ab)/(1 − q·ab) in projective form.
       Translation matrices compose; det = wall polynomial; tr² − 4·det = a²(p² − 4q);
       the fixed points of every translation are the walls; two walls → multiplication
       (Cayley for the tan law), one double wall → addition.
    5. Products: a mult-add fold hits the wall iff one piece does (Int), and the fold of
       the pieces is bounded by the fold of their sizes (Rat), the inequality behind
       Weierstrass's convergence test.
    6. Tate gluing (1 + a) ~ q(1 + a) respects mult-add.
  LIBRARY NOTE: in the shared library at FusedAlgebra/Topics/ComplexAnalysis.lean; merged unchanged; no exact twins.
-/
import FusedAlgebra.Base

namespace ComplexAnalysis
open CircleOp

open scoped CircleOp

set_option linter.unusedSectionVars false

/-! ## 0. Complex numbers as records -/

structure Cx (R : Type) where
  re : R
  im : R
deriving DecidableEq

theorem Cx.ext' {R : Type} {x y : Cx R} (h1 : x.re = y.re) (h2 : x.im = y.im) : x = y := by
  cases x; cases y; dsimp only at h1 h2; subst h1; subst h2; rfl

section Ring
variable {R : Type} [Lean.Grind.CommRing R]

namespace Cx
def one : Cx R := ⟨1, 0⟩
def zero : Cx R := ⟨0, 0⟩
def wall : Cx R := ⟨-1, 0⟩
def I : Cx R := ⟨0, 1⟩
def add (z w : Cx R) : Cx R := ⟨z.re + w.re, z.im + w.im⟩
def sub (z w : Cx R) : Cx R := ⟨z.re - w.re, z.im - w.im⟩
def mul (z w : Cx R) : Cx R := ⟨z.re * w.re - z.im * w.im, z.re * w.im + z.im * w.re⟩
/-- Squared length |z|². -/
def normSq (z : Cx R) : R := z.re * z.re + z.im * z.im
/-- Complex mult-add z ◯ w = z + w + zw. -/
def ma (z w : Cx R) : Cx R := add (add z w) (mul z w)
/-- n copies of z under complex mult-add. -/
def mpowC (z : Cx R) : Nat → Cx R
  | 0 => zero
  | n + 1 => ma (mpowC z n) z
end Cx

open Cx

/-! ## 1. A spin and a stretch about the wall -/

/-- The lift on the complex plane. -/
theorem ma_lift (z w : Cx R) : add one (ma z w) = mul (add one z) (add one w) := by
  apply Cx.ext' <;> simp only [ma, add, mul, one] <;> grind

theorem ma_comm (z w : Cx R) : ma z w = ma w z := by
  apply Cx.ext' <;> simp only [ma, add, mul] <;> grind

theorem ma_assoc (z w u : Cx R) : ma (ma z w) u = ma z (ma w u) := by
  apply Cx.ext' <;> simp only [ma, add, mul] <;> grind

theorem ma_zero (z : Cx R) : ma z zero = z := by
  apply Cx.ext' <;> simp only [ma, add, mul, zero] <;> grind

/-- The wall −1 absorbs: it is the fixed point of every map z ↦ z ◯ w. -/
theorem wall_fixed (w : Cx R) : ma wall w = wall := by
  apply Cx.ext' <;> simp only [ma, add, mul, wall] <;> grind

/-- Seen from the wall, z ↦ z ◯ w is multiplication by 1 + w: a spin and a stretch. -/
theorem spin_about_wall (z w : Cx R) :
    sub (ma z w) wall = mul (add one w) (sub z wall) := by
  apply Cx.ext' <;> simp only [ma, add, sub, mul, one, wall] <;> grind

/-- The step from z to z ◯ w is the vector w turned and stretched by 1 + z
    (so triangle (−1, 0, w) and triangle (−1, z, z ◯ w) are similar). -/
theorem step_is_turned_w (z w : Cx R) : sub (ma z w) z = mul (add one z) w := by
  apply Cx.ext' <;> simp only [ma, add, sub, mul, one] <;> grind

/-- Lengths multiply (Brahmagupta–Fibonacci two-square identity). -/
theorem normSq_mul (z w : Cx R) : normSq (mul z w) = normSq z * normSq w := by
  simp only [normSq, mul]; grind

/-- Squared distances from the wall multiply under mult-add. -/
theorem wall_distance_mul (z w : Cx R) :
    normSq (add one (ma z w)) = normSq (add one z) * normSq (add one w) := by
  rw [ma_lift, normSq_mul]

/-- The circle |1 + z| = 1 (centre the wall, through 0) is closed under mult-add. -/
theorem circle_closed (z w : Cx R) (hz : normSq (add one z) = 1) (hw : normSq (add one w) = 1) :
    normSq (add one (ma z w)) = 1 := by
  rw [wall_distance_mul, hz, hw]; grind

/-- The same circle as an equation: 2x + x² + y² = 0, the k = 1 rotation condition of
    the linear-algebra page. -/
theorem circle_equation (z : Cx R) :
    normSq (add one z) - 1 = 2 * z.re + (z.re * z.re + z.im * z.im) := by
  simp only [normSq, add, one]; grind

/-! ## 2. Affine maps that fix the wall -/

/-- An affine map z ↦ αz + β that fixes −1 is mult-add by α − 1. -/
theorem affine_fixing_wall (α β z : R) (h : α * (-1) + β = -1) : α * z + β = z ◯ (α - 1) := by
  unfold circ; grind

/-- Composing two of them mult-adds their parameters. -/
theorem affine_compose (z c d : R) : (z ◯ c) ◯ d = z ◯ (c ◯ d) := circ_assoc z c d

/-! ## 3. Mult-add roots of zero -/

/-- If ζⁿ = 1 then ζ − 1 is an n-th mult-add root of 0. -/
theorem mroot_of_unity (ζ : R) (n : Nat) (h : ζ ^ n = 1) : mpow (ζ - 1) n = 0 := by
  rw [mpow_eq]
  have e : (1 + (ζ - 1)) = ζ := by grind
  rw [e, h]; grind

/-- n = 4: i − 1 mult-added four times is 0 (any commutative ring). -/
theorem gauss_root_four : mpowC (⟨-1, 1⟩ : Cx R) 4 = zero := by
  apply Cx.ext' <;> simp only [mpowC, ma, add, mul, zero] <;> grind

/-- …but not before: twice gives −2, three times gives −1 − i (on the integers). -/
theorem gauss_root_steps :
    mpowC (⟨-1, 1⟩ : Cx Int) 2 = ⟨-2, 0⟩ ∧ mpowC (⟨-1, 1⟩ : Cx Int) 3 = ⟨-1, -1⟩ := by
  decide

/-- The four roots of the square mult-add to −2 (even n). -/
theorem gauss_roots_total :
    ma (ma (ma (zero : Cx Int) ⟨-1, 1⟩) ⟨-2, 0⟩) ⟨-1, -1⟩ = ⟨-2, 0⟩ := by decide

/-- All four lie on the circle |1 + z| = 1. -/
theorem gauss_roots_on_circle :
    normSq (add one (⟨-1, 1⟩ : Cx Int)) = 1 ∧ normSq (add one (⟨-2, 0⟩ : Cx Int)) = 1 ∧
    normSq (add one (⟨-1, -1⟩ : Cx Int)) = 1 ∧ normSq (add one (zero : Cx Int)) = 1 := by
  decide

/-- Eisenstein numbers a + bω with ω² = −1 − ω. -/
structure Eis (R : Type) where
  a : R
  b : R

theorem Eis.ext' {x y : Eis R} (h1 : x.a = y.a) (h2 : x.b = y.b) : x = y := by
  cases x; cases y; dsimp only at h1 h2; subst h1; subst h2; rfl

namespace Eis
def zero : Eis R := ⟨0, 0⟩
def add (x y : Eis R) : Eis R := ⟨x.a + y.a, x.b + y.b⟩
def mul (x y : Eis R) : Eis R := ⟨x.a * y.a - x.b * y.b, x.a * y.b + x.b * y.a - x.b * y.b⟩
def ma (x y : Eis R) : Eis R := add (add x y) (mul x y)
end Eis

/-- ω really is a cube root of 1: ω·ω·ω = 1. -/
theorem omega_cubed : Eis.mul (Eis.mul (⟨0, 1⟩ : Eis R) ⟨0, 1⟩) ⟨0, 1⟩ = ⟨1, 0⟩ := by
  apply Eis.ext' <;> simp only [Eis.mul] <;> grind

/-- n = 3: ω − 1 mult-added three times is 0. -/
theorem eis_root_three :
    Eis.ma (Eis.ma (⟨-1, 1⟩ : Eis R) ⟨-1, 1⟩) ⟨-1, 1⟩ = Eis.zero := by
  apply Eis.ext' <;> simp only [Eis.ma, Eis.add, Eis.mul, Eis.zero] <;> grind

/-- The three roots 0, ω − 1, ω² − 1 = −2 − ω mult-add to 0 (odd n). -/
theorem eis_roots_total : Eis.ma (Eis.ma (Eis.zero : Eis R) ⟨-1, 1⟩) ⟨-2, -1⟩ = Eis.zero := by
  apply Eis.ext' <;> simp only [Eis.ma, Eis.add, Eis.mul, Eis.zero] <;> grind

/-! ## 4. Möbius: the landscape laws -/

/-- The fixed points of b ↦ F(a,b) solve a·(1 + p·b + q·b²) = 0: the walls, for every a. -/
theorem fixed_points_are_walls (p q a b : R) :
    (a + b + p * a * b) - b * (1 - q * a * b) = a * (1 + p * b + q * b * b) := by grind

/-- A Möbius map in projective form [[α, β], [γ, δ]]. -/
structure Mob (R : Type) where
  α : R
  β : R
  γ : R
  δ : R

theorem Mob.ext' {x y : Mob R} (h1 : x.α = y.α) (h2 : x.β = y.β) (h3 : x.γ = y.γ)
    (h4 : x.δ = y.δ) : x = y := by
  cases x; cases y; dsimp only at h1 h2 h3 h4; subst h1; subst h2; subst h3; subst h4; rfl

def Mob.comp (x y : Mob R) : Mob R :=
  ⟨x.α * y.α + x.β * y.γ, x.α * y.β + x.β * y.δ, x.γ * y.α + x.δ * y.γ, x.γ * y.β + x.δ * y.δ⟩

/-- Translation by the point N/D of the law (p, q): b ↦ F(N/D, b). -/
def trans (p q N D : R) : Mob R := ⟨D + p * N, N, -(q * N), D⟩

/-- Translations compose by the law itself: T(a)∘T(b) = T(a + b + p·ab over 1 − q·ab). -/
theorem trans_compose (p q a b : R) :
    Mob.comp (trans p q a 1) (trans p q b 1) = trans p q (a + b + p * a * b) (1 - q * a * b) := by
  apply Mob.ext' <;> simp only [Mob.comp, trans] <;> grind

/-- det T(a) is the wall polynomial: the translation degenerates exactly at a wall. -/
theorem trans_det (p q a : R) :
    (trans p q a 1).α * (trans p q a 1).δ - (trans p q a 1).β * (trans p q a 1).γ
      = 1 + p * a + q * a * a := by
  simp only [trans]; grind

/-- tr² − 4·det = N²(p² − 4q): the discriminant alone sets every translation's type. -/
theorem trans_type (p q N D : R) :
    ((trans p q N D).α + (trans p q N D).δ) * ((trans p q N D).α + (trans p q N D).δ)
      - 4 * ((trans p q N D).α * (trans p q N D).δ - (trans p q N D).β * (trans p q N D).γ)
      = N * N * (p * p - 4 * q) := by
  simp only [trans]; grind

/-- Two walls (reciprocals r₁, r₂, so p = −(r₁ + r₂), q = r₁r₂): ζ(x) = (1 − r₁x)/(1 − r₂x)
    turns the law into multiplication. Cross-multiplied, with F = N/D. -/
theorem two_walls_lift (r1 r2 a b : R) :
    let N := a + b + (-(r1 + r2)) * a * b
    let D := 1 - r1 * r2 * a * b
    (D - r1 * N) * ((1 - r2 * a) * (1 - r2 * b)) = ((1 - r1 * a) * (1 - r1 * b)) * (D - r2 * N) := by
  intro N D; simp only [N, D]; grind

/-- Mult-add is the case of one wall at −1 and the other at infinity (r₁ = −1, r₂ = 0). -/
theorem two_walls_circ (a b : R) :
    let N := a + b + (-((-1) + 0)) * a * b
    (1 - (-1) * N) = (1 - (-1) * a) * (1 - (-1) * b) := by
  intro N; simp only [N]; grind

/-- Cayley: the tan law (p, q) = (0, 1) has walls ±i, and (1 + i·x)/(1 − i·x) multiplies. -/
theorem cayley_tan (i a b : R) (hi : i * i = -1) :
    let N := a + b
    let D := 1 - a * b
    (D + i * N) * ((1 - i * a) * (1 - i * b)) = ((1 + i * a) * (1 + i * b)) * (D - i * N) := by
  intro N D; simp only [N, D]; grind

/-- One double wall (p = 2c, q = c²): ψ(x) = x/(1 + cx) turns the law into addition. -/
theorem double_wall_adds (c a b : R) :
    let N := a + b + 2 * c * a * b
    let D := 1 - c * c * a * b
    N * ((1 + c * a) * (1 + c * b)) = (D + c * N) * (a * (1 + c * b) + b * (1 + c * a)) := by
  intro N D; simp only [N, D]; grind

end Ring

/-! ## 5. Products: wall hits and the size bound -/

section Products
variable {R : Type} [Lean.Grind.CommRing R]

/-- If one piece sits on the wall, the whole fold sits on the wall. -/
theorem fold_wall_of_mem (l : List R) (h : (-1 : R) ∈ l) : foldCirc l = -1 := by
  induction l with
  | nil => cases h
  | cons a l ih =>
    simp only [foldCirc]
    cases h with
    | head => exact neg_one_absorbs _
    | tail _ hm => rw [ih hm, circ_comm]; exact neg_one_absorbs _

end Products

theorem prodL_eq_zero_int (l : List Int) (h : prodL l = 0) : (0 : Int) ∈ l := by
  induction l with
  | nil => simp only [prodL] at h; exact absurd h (by decide)
  | cons a l ih =>
    simp only [prodL] at h
    rcases Int.mul_eq_zero.mp h with h0 | h0
    · rw [h0]; exact List.Mem.head _
    · exact List.Mem.tail _ (ih h0)

/-- On the integers the converse holds too: a fold hits the wall only if a piece does
    (the zeros of a product are the zeros of its factors). -/
theorem fold_wall_iff_int (l : List Int) : foldCirc l = -1 ↔ (-1 : Int) ∈ l := by
  constructor
  · intro h
    have hl := foldCirc_lift l
    rw [h] at hl
    have h0 : prodL (l.map (fun x => 1 + x)) = 0 := by rw [← hl]; decide
    have hm := prodL_eq_zero_int _ h0
    rcases List.mem_map.mp hm with ⟨x, hx, hx0⟩
    have : x = -1 := by omega
    rw [← this]; exact hx
  · exact fold_wall_of_mem l

/-- The size bound for one corner: if |a| ≤ A and |b| ≤ B then |a ◯ b| ≤ A ◯ B. -/
theorem corner_bound (a b A B : Rat) (ha1 : -A ≤ a) (ha2 : a ≤ A) (hb1 : -B ≤ b) (hb2 : b ≤ B) :
    -(A ◯ B) ≤ a ◯ b ∧ a ◯ b ≤ A ◯ B := by
  have p1 : 0 ≤ (A - a) * (B + b) := Rat.mul_nonneg (by grind) (by grind)
  have p2 : 0 ≤ (A + a) * (B - b) := Rat.mul_nonneg (by grind) (by grind)
  have p3 : 0 ≤ (A - a) * (B - b) := Rat.mul_nonneg (by grind) (by grind)
  have p4 : 0 ≤ (A + a) * (B + b) := Rat.mul_nonneg (by grind) (by grind)
  unfold circ
  constructor <;> grind

/-- Sizes are never negative under mult-add of sizes. -/
theorem size_nonneg (A B : Rat) (hA : 0 ≤ A) (hB : 0 ≤ B) : 0 ≤ A ◯ B := by
  have := Rat.mul_nonneg hA hB
  unfold circ; grind

/-- Weierstrass's bound: the fold of the pieces is bounded by the fold of their sizes,
    |(a₁ ◯ … ◯ aₙ)| ≤ |a₁| ◯ … ◯ |aₙ| = ∏(1 + |aᵢ|) − 1. -/
theorem fold_bound : ∀ (l : List (Rat × Rat)), (∀ x ∈ l, -x.2 ≤ x.1 ∧ x.1 ≤ x.2) →
    -(foldCirc (l.map Prod.snd)) ≤ foldCirc (l.map Prod.fst) ∧
      foldCirc (l.map Prod.fst) ≤ foldCirc (l.map Prod.snd)
  | [], _ => by simp only [List.map, foldCirc]; grind
  | x :: l, h => by
    simp only [List.map, foldCirc]
    have hx := h x (List.Mem.head _)
    have ih := fold_bound l (fun y hy => h y (List.Mem.tail _ hy))
    exact corner_bound _ _ _ _ hx.1 hx.2 ih.1 ih.2

/-! ## 6. Tate gluing -/

/-- Gluing (1 + a) ~ q(1 + a) respects mult-add, so the punctured plane folds into a torus. -/
theorem tate_glue {R : Type} [Lean.Grind.CommRing R] (q a a' b : R) (h : 1 + a' = q * (1 + a)) :
    1 + (a' ◯ b) = q * (1 + (a ◯ b)) := by
  rw [circ_lift, circ_lift, h]; grind

end ComplexAnalysis
