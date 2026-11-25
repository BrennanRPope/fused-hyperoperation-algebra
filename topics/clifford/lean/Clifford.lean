/-
  Fused Hyperoperation Algebra: W. K. Clifford's work, piece by piece, against mult-add.
  Core Lean 4.37 only (no mathlib), no `sorry`. Imports the shared library's Base and GroupsLie
  (for the θ-lift `tmul`, so the new statements attach to the existing ones).

  Intended library path: FusedAlgebra/Topics/Clifford.lean (after GroupsLie).

  What is here:
    1. The Clifford algebra of the rate form x² + pxy + qy² (generators u, v with u² = 1, v² = q,
       uv + vu = p), as 4-tuples. It is associative, a vector squares to the form, and θ = uv
       obeys θ² = pθ − q: the θ-lift of GroupsLie is its even part. Two lifts multiply to the
       landscape law; at (p, q) = (1, 0) they multiply to mult-add.
    2. Conjugation θ ↦ p − θ on the θ-lift: multiplicative, X·X̄ is the norm (a scalar), the lift
       times its conjugate is the rate V(a), and mult-add's wall −1 is a zero divisor.
    3. A two-part dial K = k₀ + k₁e with e² = κ: K² = tr(K)·K − N(K); the lift product reads as
       the landscape law with p = tr K, q = N K; the class invariant p² − 4q is 4κk₁²; mult-add's
       dial is the null element c(1 + e) (κ = 1, N = 0).
    4. The three 2-D Clifford algebras: no zero divisor when κ = −1 (complex: the norm vanishes
       only at 0), a zero divisor (1 + e)(1 − e) = 0 when κ = 1, a nilpotent e² = 0 when κ = 0.
       Clifford's elliptic ω (κ = 1): (1 + ω)² = 2(1 + ω), and x + yω ↦ (x + y, x − y) is a
       multiplicative split into two copies of the line (left and right halves).
    5. Cl(2,0)-type algebra with squares κ₁, κ₂, as 4-tuples: associative; the multivector
       mult-add X + Y + XKY is associative and has the lift (1 + KX)(1 + KY) = 1 + K(X ◯ Y);
       its commutator is XKY − YKX; for vectors with K = 1 it is twice the wedge (oriented area),
       and the corner of two vectors splits into dot + wedge; reversion turns X ◯_K Y into
       Ỹ ◯_K̃ X̃; two factors in different directions keep their product on its own blade.
    6. Dual numbers: mult-add with ε² = 0 carries the derivative (a₁(1 + b₀) + b₁(1 + a₀)).
    7. The q = +1 twin of the light-wall identity: |1 − ūv|² + |u + v|² = (1 + |u|²)(1 + |v|²)
       (norm of the quaternion lift; no wall on the sphere side).
    8. Bott's count: the number of blades squaring to +1 minus those squaring to −1, sig(p, q),
       obeys sig(p+1, q+1) = 2·sig(p, q) and sig(p+8, q) = 16·sig(p, q) for p, q ≤ 4, and among
       p + q ≤ 8 only (0,0), (0,1), (0,2) (R, C, H) have no blade squaring to +1.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.GroupsLie

namespace CliffordW
open CircleOp

set_option linter.unusedSectionVars false

open scoped CircleOp

/-! ## 1. The Clifford algebra of the rate form x² + pxy + qy² -/

/-- x₀ + x₁u + x₂v + x₃uv. -/
structure B4 where
  c0 : Int
  c1 : Int
  c2 : Int
  c3 : Int
deriving DecidableEq

theorem B4.ext' {x y : B4} (h0 : x.c0 = y.c0) (h1 : x.c1 = y.c1) (h2 : x.c2 = y.c2)
    (h3 : x.c3 = y.c3) : x = y := by
  cases x; cases y; simp_all

/-- Product from u² = 1, v² = q, uv + vu = p. -/
def B4.mul (p q : Int) (x y : B4) : B4 :=
  ⟨x.c0 * y.c0 + x.c1 * y.c1 + p * x.c2 * y.c1 + q * x.c2 * y.c2 - q * x.c3 * y.c3,
   x.c0 * y.c1 + x.c1 * y.c0 - q * x.c2 * y.c3 + p * x.c3 * y.c1 + q * x.c3 * y.c2,
   x.c0 * y.c2 + x.c2 * y.c0 + x.c1 * y.c3 + p * x.c2 * y.c3 - x.c3 * y.c1,
   x.c0 * y.c3 + x.c3 * y.c0 + x.c1 * y.c2 - x.c2 * y.c1 + p * x.c3 * y.c3⟩

theorem B4.mul_assoc (p q : Int) (x y z : B4) :
    B4.mul p q (B4.mul p q x y) z = B4.mul p q x (B4.mul p q y z) := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

/-- The generators obey the defining relations. -/
theorem B4.u_sq (p q : Int) : B4.mul p q ⟨0, 1, 0, 0⟩ ⟨0, 1, 0, 0⟩ = ⟨1, 0, 0, 0⟩ := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

theorem B4.v_sq (p q : Int) : B4.mul p q ⟨0, 0, 1, 0⟩ ⟨0, 0, 1, 0⟩ = ⟨q, 0, 0, 0⟩ := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

theorem B4.anticomm (p q : Int) :
    (B4.mul p q ⟨0, 1, 0, 0⟩ ⟨0, 0, 1, 0⟩).c0 + (B4.mul p q ⟨0, 0, 1, 0⟩ ⟨0, 1, 0, 0⟩).c0 = p ∧
    (B4.mul p q ⟨0, 1, 0, 0⟩ ⟨0, 0, 1, 0⟩).c3 + (B4.mul p q ⟨0, 0, 1, 0⟩ ⟨0, 1, 0, 0⟩).c3 = 0 := by
  simp only [B4.mul]; constructor <;> grind

/-- A vector squares to the rate form: (xu + yv)² = x² + pxy + qy². -/
theorem B4.vector_sq (p q x y : Int) :
    B4.mul p q ⟨0, x, y, 0⟩ ⟨0, x, y, 0⟩ = ⟨x * x + p * x * y + q * y * y, 0, 0, 0⟩ := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

/-- θ = uv squares to pθ − q: the θ-lift is the even part of this Clifford algebra. -/
theorem B4.theta_sq (p q : Int) : B4.mul p q ⟨0, 0, 0, 1⟩ ⟨0, 0, 0, 1⟩ = ⟨-q, 0, 0, p⟩ := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

/-- The even part multiplies exactly like `tmul` (GroupsLie): same algebra. -/
theorem B4.even_is_tmul (p q x y z w : Int) :
    B4.mul p q ⟨x, 0, 0, y⟩ ⟨z, 0, 0, w⟩ =
      ⟨(GroupsLie.tmul p q (x, y) (z, w)).1, 0, 0, (GroupsLie.tmul p q (x, y) (z, w)).2⟩ := by
  apply B4.ext' <;> simp only [B4.mul, GroupsLie.tmul] <;> grind

/-- Two lifts multiply to the landscape law: (1 + aθ)(1 + bθ) = (1 − qab) + (a + b + pab)θ. -/
theorem B4.lift_mul (p q a b : Int) :
    B4.mul p q ⟨1, 0, 0, a⟩ ⟨1, 0, 0, b⟩ = ⟨1 - q * a * b, 0, 0, a + b + p * a * b⟩ := by
  apply B4.ext' <;> simp only [B4.mul] <;> grind

/-- At mult-add's point (1, 0) the geometric product of lifts is the lift of a ◯ b. -/
theorem B4.lift_circ (a b : Int) :
    B4.mul 1 0 ⟨1, 0, 0, a⟩ ⟨1, 0, 0, b⟩ = ⟨1, 0, 0, a ◯ b⟩ := by
  apply B4.ext' <;> simp only [B4.mul, circ] <;> grind

/-- Mult-add's rate form x² + xy factors as x(x + y): two light lines (a split form). -/
theorem circ_form_split (x y : Int) : x * x + 1 * x * y + 0 * y * y = x * (x + y) := by grind

/-! ## 2. Conjugation, norm and the wall -/

def tconj (p : Int) (u : Int × Int) : Int × Int := (u.1 + p * u.2, -u.2)

theorem tconj_mul (p q : Int) (u v : Int × Int) :
    tconj p (GroupsLie.tmul p q u v) = GroupsLie.tmul p q (tconj p u) (tconj p v) := by
  simp only [tconj, GroupsLie.tmul]; ext <;> simp <;> grind

theorem tconj_tconj (p : Int) (u : Int × Int) : tconj p (tconj p u) = u := by
  simp only [tconj]; ext <;> simp <;> grind

/-- X·X̄ is the scalar x² + pxy + qy² (the norm form, = the rate form). -/
theorem mul_tconj (p q : Int) (u : Int × Int) :
    GroupsLie.tmul p q u (tconj p u) = (u.1 * u.1 + p * u.1 * u.2 + q * u.2 * u.2, 0) := by
  simp only [tconj, GroupsLie.tmul]; ext <;> simp <;> grind

/-- The lift times its conjugate is the rate V(a) = 1 + pa + qa². -/
theorem lift_tconj (p q a : Int) :
    GroupsLie.tmul p q (1, a) (tconj p (1, a)) = (1 + p * a + q * a * a, 0) := by
  simp only [tconj, GroupsLie.tmul]; ext <;> simp <;> grind

/-- Mult-add: lift times conjugate = 1 + a, so the conjugate (1 + a) − aθ is (1 + a)·(lift of the inverse). -/
theorem circ_lift_tconj (a : Int) : GroupsLie.tmul 1 0 (1, a) (tconj 1 (1, a)) = (1 + a, 0) := by
  simp only [tconj, GroupsLie.tmul]; ext <;> simp <;> grind

/-- The conjugate of a lift, read in the chart, is the landscape inverse: its θ-slot over its
    scalar slot, −a/(1 + pa), composes with a to 0 (numerator of F vanishes). -/
theorem tconj_chart_inverse (p a : Int) :
    a * (1 + p * a) + (-a) + p * a * (-a) = (1 + p * a) * a - a - p * a * a := by grind

/-- Mult-add's wall is a zero divisor: (1 − θ)·θ = 0 when θ² = θ. -/
theorem wall_zero_divisor : GroupsLie.tmul 1 0 (1, -1) (0, 1) = (0, 0) := by decide

/-- θ is idempotent at mult-add's point. -/
theorem theta_idem : GroupsLie.tmul 1 0 (0, 1) (0, 1) = (0, 1) := by decide

/-! ## 3. A two-part dial K = k₀ + k₁e, e² = κ -/

/-- x + ye with e² = κ. -/
def cmul (κ : Int) (u v : Int × Int) : Int × Int :=
  (u.1 * v.1 + κ * u.2 * v.2, u.1 * v.2 + u.2 * v.1)

def ctr (K : Int × Int) : Int := 2 * K.1
def cnorm (κ : Int) (K : Int × Int) : Int := K.1 * K.1 - κ * K.2 * K.2

theorem cmul_assoc (κ : Int) (u v w : Int × Int) :
    cmul κ (cmul κ u v) w = cmul κ u (cmul κ v w) := by
  simp only [cmul]; ext <;> simp <;> grind

theorem cmul_comm (κ : Int) (u v : Int × Int) : cmul κ u v = cmul κ v u := by
  simp only [cmul]; ext <;> simp <;> grind

/-- Cayley–Hamilton: K² = tr(K)·K − N(K). -/
theorem dial_sq (κ : Int) (K : Int × Int) :
    cmul κ K K = (ctr K * K.1 - cnorm κ K, ctr K * K.2) := by
  simp only [cmul, ctr, cnorm]; ext <;> simp <;> grind

/-- (1 + aK)(1 + bK) = (1 − N·ab) + (a + b + tr·ab)K: a multivector dial with scalar arguments
    gives the landscape law with p = tr K, q = N K. -/
theorem dial_lift (κ a b : Int) (K : Int × Int) :
    cmul κ (1 + a * K.1, a * K.2) (1 + b * K.1, b * K.2) =
      (1 - cnorm κ K * a * b + (a + b + ctr K * a * b) * K.1, (a + b + ctr K * a * b) * K.2) := by
  simp only [cmul, ctr, cnorm]; ext <;> simp <;> grind

/-- The class invariant p² − 4q of the dial is 4κk₁²: the country is the sign of κ. -/
theorem dial_discriminant (κ : Int) (K : Int × Int) :
    ctr K * ctr K - 4 * cnorm κ K = 4 * κ * K.2 * K.2 := by
  simp only [ctr, cnorm]; grind

/-- Mult-add's dial is null: K = c(1 + e) with e² = 1 has N = 0 and tr = 2c, so K² = 2c·K. -/
theorem circ_dial_null (c : Int) :
    cnorm 1 (c, c) = 0 ∧ ctr (c, c) = 2 * c ∧ cmul 1 (c, c) (c, c) = (2 * c * c, 2 * c * c) := by
  simp only [cnorm, ctr, cmul]; refine ⟨?_, ?_, ?_⟩ <;> (try ext) <;> simp <;> grind

/-! ## 4. The three 2-D Clifford algebras, and Clifford's ω -/

/-- Squares are nonnegative (core Int has no named lemma for this). -/
theorem sq_nonneg' (x : Int) : 0 ≤ x * x := by
  rcases Int.le_total 0 x with h | h
  · exact Int.mul_nonneg h h
  · have := Int.mul_nonneg (Int.neg_nonneg_of_nonpos h) (Int.neg_nonneg_of_nonpos h)
    simpa [Int.neg_mul_neg] using this

theorem sq_eq_zero' (x : Int) (h : x * x = 0) : x = 0 := by
  rcases Int.mul_eq_zero.mp h with h | h <;> exact h

/-- κ = −1 (complex numbers): the norm vanishes only at 0, so the wall stays one point. -/
theorem complex_no_zero_norm (x y : Int) (h : x * x + y * y = 0) : x = 0 ∧ y = 0 := by
  have hx : 0 ≤ x * x := sq_nonneg' x
  have hy : 0 ≤ y * y := sq_nonneg' y
  have hx0 : x * x = 0 := by omega
  have hy0 : y * y = 0 := by omega
  exact ⟨sq_eq_zero' x hx0, sq_eq_zero' y hy0⟩

/-- κ = +1 (split-complex): (1 + e)(1 − e) = 0, a zero divisor; the wall grows to two lines. -/
theorem split_zero_divisor : cmul 1 (1, 1) (1, -1) = (0, 0) := by decide

/-- κ = 0 (dual numbers): e² = 0. -/
theorem dual_nilpotent : cmul 0 (0, 1) (0, 1) = (0, 0) := by decide

/-- Clifford's elliptic ω (ω² = 1): (1 + ω)² = 2(1 + ω), so (1 + ω)/2 is idempotent, like θ. -/
theorem omega_projector : cmul 1 (1, 1) (1, 1) = (2, 2) := by decide

/-- x + yω ↦ (x + y, x − y) turns ω-multiplication into slot-by-slot multiplication:
    the elliptic biquaternion splits into a left half and a right half. -/
theorem omega_split (u v : Int × Int) :
    ((cmul 1 u v).1 + (cmul 1 u v).2, (cmul 1 u v).1 - (cmul 1 u v).2) =
      ((u.1 + u.2) * (v.1 + v.2), (u.1 - u.2) * (v.1 - v.2)) := by
  simp only [cmul]; ext <;> simp <;> grind

/-- The quaternion case: the norm of x₀ + x₁i + x₂j + x₃k vanishes only at 0. -/
theorem quaternion_no_zero_norm (a b c d : Int) (h : a * a + b * b + c * c + d * d = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 := by
  have ha : 0 ≤ a * a := sq_nonneg' a
  have hb : 0 ≤ b * b := sq_nonneg' b
  have hc : 0 ≤ c * c := sq_nonneg' c
  have hd : 0 ≤ d * d := sq_nonneg' d
  refine ⟨sq_eq_zero' a (by omega), sq_eq_zero' b (by omega),
    sq_eq_zero' c (by omega), sq_eq_zero' d (by omega)⟩

/-! ## 5. Multivector mult-add in a Clifford algebra with squares κ₁, κ₂ -/

/-- x₀ + x₁e₁ + x₂e₂ + x₃e₁₂. -/
structure G4 where
  s : Int
  e1 : Int
  e2 : Int
  b : Int
deriving DecidableEq

theorem G4.ext' {x y : G4} (h0 : x.s = y.s) (h1 : x.e1 = y.e1) (h2 : x.e2 = y.e2)
    (h3 : x.b = y.b) : x = y := by
  cases x; cases y; simp_all

def G4.mul (k1 k2 : Int) (x y : G4) : G4 :=
  ⟨x.s * y.s + k1 * x.e1 * y.e1 + k2 * x.e2 * y.e2 - k1 * k2 * x.b * y.b,
   x.s * y.e1 + x.e1 * y.s - k2 * x.e2 * y.b + k2 * x.b * y.e2,
   x.s * y.e2 + x.e2 * y.s + k1 * x.e1 * y.b - k1 * x.b * y.e1,
   x.s * y.b + x.b * y.s + x.e1 * y.e2 - x.e2 * y.e1⟩
def G4.add (x y : G4) : G4 := ⟨x.s + y.s, x.e1 + y.e1, x.e2 + y.e2, x.b + y.b⟩
def G4.sub (x y : G4) : G4 := ⟨x.s - y.s, x.e1 - y.e1, x.e2 - y.e2, x.b - y.b⟩
def G4.one : G4 := ⟨1, 0, 0, 0⟩
def G4.rev (x : G4) : G4 := ⟨x.s, x.e1, x.e2, -x.b⟩

/-- Multivector mult-add X ◯_K Y = X + Y + XKY. -/
def G4.star (k1 k2 : Int) (K x y : G4) : G4 :=
  G4.add (G4.add x y) (G4.mul k1 k2 (G4.mul k1 k2 x K) y)

theorem G4.mul_assoc (k1 k2 : Int) (x y z : G4) :
    G4.mul k1 k2 (G4.mul k1 k2 x y) z = G4.mul k1 k2 x (G4.mul k1 k2 y z) := by
  apply G4.ext' <;> simp only [G4.mul] <;> grind

theorem G4.e12_sq (k1 k2 : Int) : G4.mul k1 k2 ⟨0, 0, 0, 1⟩ ⟨0, 0, 0, 1⟩ = ⟨-(k1 * k2), 0, 0, 0⟩ := by
  apply G4.ext' <;> simp only [G4.mul] <;> grind

/-- Kept: the left lift. 1 + K(X ◯_K Y) = (1 + KX)(1 + KY). -/
theorem G4.star_lift (k1 k2 : Int) (K x y : G4) :
    G4.add G4.one (G4.mul k1 k2 K (G4.star k1 k2 K x y)) =
      G4.mul k1 k2 (G4.add G4.one (G4.mul k1 k2 K x)) (G4.add G4.one (G4.mul k1 k2 K y)) := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.one, G4.star] <;> grind

/-- Kept: associativity, for every multivector dial K. -/
theorem G4.star_assoc (k1 k2 : Int) (K x y z : G4) :
    G4.star k1 k2 K (G4.star k1 k2 K x y) z = G4.star k1 k2 K x (G4.star k1 k2 K y z) := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.star] <;> grind

theorem G4.star_zero (k1 k2 : Int) (K x : G4) : G4.star k1 k2 K x ⟨0, 0, 0, 0⟩ = x := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.star] <;> grind

/-- Lost: commutativity. X ◯ Y − Y ◯ X = XKY − YKX. -/
theorem G4.star_comm_defect (k1 k2 : Int) (K x y : G4) :
    G4.sub (G4.star k1 k2 K x y) (G4.star k1 k2 K y x) =
      G4.sub (G4.mul k1 k2 (G4.mul k1 k2 x K) y) (G4.mul k1 k2 (G4.mul k1 k2 y K) x) := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.sub, G4.star] <;> grind

/-- For two vectors and K = 1 the defect is twice the wedge: 2(x₁y₂ − x₂y₁)e₁₂, twice the
    oriented area, whatever the squares κ₁, κ₂. -/
theorem G4.vector_comm_wedge (k1 k2 x1 x2 y1 y2 : Int) :
    G4.sub (G4.star k1 k2 G4.one ⟨0, x1, x2, 0⟩ ⟨0, y1, y2, 0⟩)
        (G4.star k1 k2 G4.one ⟨0, y1, y2, 0⟩ ⟨0, x1, x2, 0⟩) = ⟨0, 0, 0, 2 * (x1 * y2 - x2 * y1)⟩ := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.sub, G4.star, G4.one] <;> grind

/-- The corner of two vectors (K = 1, Euclidean squares) splits into dot + wedge:
    x ◯ y = x + y + (x·y) + (x∧y)e₁₂. -/
theorem G4.vector_corner (x1 x2 y1 y2 : Int) :
    G4.star 1 1 G4.one ⟨0, x1, x2, 0⟩ ⟨0, y1, y2, 0⟩ =
      ⟨x1 * y1 + x2 * y2, x1 + y1, x2 + y2, x1 * y2 - x2 * y1⟩ := by
  apply G4.ext' <;> simp only [G4.mul, G4.add, G4.star, G4.one] <;> grind

/-- Reversion reverses products. -/
theorem G4.rev_mul (k1 k2 : Int) (x y : G4) :
    G4.rev (G4.mul k1 k2 x y) = G4.mul k1 k2 (G4.rev y) (G4.rev x) := by
  apply G4.ext' <;> simp only [G4.mul, G4.rev] <;> grind

/-- Reversion turns X ◯_K Y into Ỹ ◯_K̃ X̃. -/
theorem G4.rev_star (k1 k2 : Int) (K x y : G4) :
    G4.rev (G4.star k1 k2 K x y) = G4.star k1 k2 (G4.rev K) (G4.rev y) (G4.rev x) := by
  apply G4.ext' <;> simp only [G4.mul, G4.rev, G4.add, G4.star] <;> grind

/-- Information kept: two factors in different directions keep their product on its own blade,
    a₁e₁ ◯ a₂e₂ = a₁e₁ + a₂e₂ + a₁a₂e₁₂, whatever the squares; the scalar shadow
    (sum of the parts) is a₁ ◯ a₂. -/
theorem G4.blades_keep (k1 k2 a1 a2 : Int) :
    G4.star k1 k2 G4.one ⟨0, a1, 0, 0⟩ ⟨0, 0, a2, 0⟩ = ⟨0, a1, a2, a1 * a2⟩ ∧
      0 + a1 + a2 + a1 * a2 = a1 ◯ a2 := by
  constructor
  · apply G4.ext' <;> simp only [G4.mul, G4.add, G4.star, G4.one] <;> grind
  · simp only [circ]; grind

/-- Two walls for a vector dial: (1 + a k)(1 − a k) = 1 − a²k² with k = 3e₁ + 4e₂ (k² = 25). -/
theorem G4.vector_dial_walls (a : Int) :
    G4.mul 1 1 ⟨1, 3 * a, 4 * a, 0⟩ ⟨1, -(3 * a), -(4 * a), 0⟩ = ⟨1 - 25 * a * a, 0, 0, 0⟩ := by
  apply G4.ext' <;> simp only [G4.mul] <;> grind

/-! ## 6. Dual numbers carry the derivative -/

theorem dual_circ (a0 a1 b0 b1 : Int) :
    (GroupsLie.tmul 0 0 (a0, a1) (b0, b1)).2 + a1 + b1 = a1 * (1 + b0) + b1 * (1 + a0) ∧
      (GroupsLie.tmul 0 0 (a0, a1) (b0, b1)).1 + a0 + b0 = a0 ◯ b0 := by
  simp only [GroupsLie.tmul, circ]; constructor <;> grind

/-! ## 7. The sphere twin of the light-wall identity (no wall at q = +1) -/

theorem sphere_norm_identity (u1 u2 v1 v2 : Int) :
    (1 - (u1 * v1 + u2 * v2)) * (1 - (u1 * v1 + u2 * v2)) + (u1 * v2 - u2 * v1) * (u1 * v2 - u2 * v1)
      + (u1 + v1) * (u1 + v1) + (u2 + v2) * (u2 + v2) =
      (1 + u1 * u1 + u2 * u2) * (1 + v1 * v1 + v2 * v2) := by grind

/-! ## 8. Bott's eight-fold count -/

def choose' : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose' n k + choose' n (k + 1)

/-- Sign of the square of a blade with r₊ positive and r₋ negative vectors:
    (−1)^(r(r−1)/2) (−1)^(r₋), r = r₊ + r₋. -/
def bsign (rp rm : Nat) : Int :=
  let r := rp + rm
  (if (r * (r - 1) / 2) % 2 = 0 then 1 else -1) * (if rm % 2 = 0 then 1 else -1)

/-- Blades squaring to +1 minus blades squaring to −1 in Cl(p, q). -/
def sig (p q : Nat) : Int :=
  ((List.range (p + 1)).map fun rp =>
    ((List.range (q + 1)).map fun rm =>
      (choose' p rp : Int) * (choose' q rm : Int) * bsign rp rm).sum).sum

/-- Number of non-scalar blades squaring to +1 (each gives a zero divisor 1 ± e_S). -/
def plusBlades (p q : Nat) : Nat :=
  ((List.range (p + 1)).map fun rp =>
    ((List.range (q + 1)).map fun rm =>
      if rp + rm ≥ 1 ∧ bsign rp rm = 1 then choose' p rp * choose' q rm else 0).sum).sum

theorem sig_examples : sig 0 0 = 1 ∧ sig 0 1 = 0 ∧ sig 0 2 = -2 ∧ sig 1 0 = 2 ∧ sig 1 1 = 2 ∧
    sig 3 0 = 0 ∧ sig 0 4 = -4 := by decide

/-- Cl(p+1, q+1) ≅ M₂(Cl(p, q)): the count doubles. -/
theorem sig_shift_11 : ∀ p < 5, ∀ q < 5, sig (p + 1) (q + 1) = 2 * sig p q := by decide

/-- Bott: eight more positive directions multiply the count by 16 = 2⁴. -/
theorem sig_bott_8 : ∀ p < 5, ∀ q < 5, sig (p + 8) q = 16 * sig p q := by decide

/-- The check behind `division_only`, as a Boolean table over p, q < 9. -/
def divisionTable : Bool :=
  (List.range 9).all fun p => (List.range 9).all fun q =>
    !(decide (p + q ≤ 8)) || ((plusBlades p q == 0) == (p == 0 && q ≤ 2))

/-- Only R, C, H have no blade squaring to +1 (p + q ≤ 8): the wall of 1 + X stays a single
    point exactly in the three division algebras Cl(0,0), Cl(0,1), Cl(0,2). -/
theorem division_only : divisionTable = true := by decide

theorem division_cases : plusBlades 0 0 = 0 ∧ plusBlades 0 1 = 0 ∧ plusBlades 0 2 = 0 ∧
    plusBlades 0 3 = 1 ∧ plusBlades 1 0 = 1 ∧ plusBlades 0 4 = 5 := by decide

end CliffordW
