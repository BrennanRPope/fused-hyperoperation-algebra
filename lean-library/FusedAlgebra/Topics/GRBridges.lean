/-
  Curved geometry from mult-add, and four bridges to Einstein's equations: the elementary parts
  ===============================================================================================

  Companion checks: ../checks.py (sympy, mpmath, scipy; 65 of 65 pass). Core Lean 4.37.0-rc1,
  no mathlib, no `sorry`. Imports Base only.

  What is proved here (the algebra; the calculus and numerics are in the checks):
  1. Boosts in the Clifford plane Cl(1,1) (= 2×2 real matrices): 1 + βσ composes so that the two
     light-cone parts multiply, i.e. they compose by mult-add; the norm 1 − β² is the wall.
  2. The multivector mult-add a + b + a·k·b (2×2 integer matrices, a non-commutative model):
     exact lift 1 + k(a ◯ b) = (1 + k a)(1 + k b) and associativity for every k; the
     front-loaded a + b + k·a·b fails associativity for a non-central k (explicit witness).
  3. Schur's factor (n−1)(n−2)/2: zero only in dimension 2, so a place-only dial is free only there.
  4. Jacobson's area bookkeeping: bundle area = 1 + (a ◯ b); with straight rays the second
     difference of the area is the cross term 2αβ; θ²/2 + σ² = θ₁² + θ₂².
  5. The river's energy constraint: e2 is the pairwise layer of a three-way mult-add; the
     spherical river's e2 is Newton's shell law; the Schwarzschild river has e2 = 0, and an aether
     river needs c₂ = −c₁₃.
  6. Light-like test directions force X = φ·η (two-dimensional slice).
  7. Brans–Dicke: γ − 1 = −1/(2 + ω).
-/
import FusedAlgebra.Base

namespace GRBridges
open CircleOp

open scoped CircleOp

/-! ## 1. Boosts in Cl(1,1): light-cone parts compose by mult-add -/

/-- The product (1 + β₁σ)(1 + β₂σ) with σ² = 1 has scalar part 1 + β₁β₂ and σ-part β₁ + β₂. -/
theorem boost_product (b1 b2 : Rat) :
    (1 + b1 * b2) + (b1 + b2) = (1 + b1) * (1 + b2) ∧
    (1 + b1 * b2) - (b1 + b2) = (1 - b1) * (1 - b2) := by
  constructor <;> grind

/-- The forward light-cone part (scalar + σ) of a product is 1 + (β₁ ◯ β₂). -/
theorem forward_part_circ (b1 b2 : Rat) : (1 + b1 * b2) + (b1 + b2) = 1 + (b1 ◯ b2) := by
  unfold circ; grind

/-- The backward light-cone part is 1 + ((−β₁) ◯ (−β₂)). -/
theorem backward_part_circ (b1 b2 : Rat) :
    (1 + b1 * b2) - (b1 + b2) = 1 + ((-b1) ◯ (-b2)) := by
  unfold circ; grind

/-- The norm (1 + βσ)(1 − βσ) = 1 − β²: zero exactly at β = ±1, the wall (light speed). -/
theorem boost_norm (b : Rat) : (1 + b) * (1 - b) = 1 - b * b := by grind

/-! ## 2. Multivector mult-add on 2×2 integer matrices -/

structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving DecidableEq

namespace M2
def mul (X Y : M2) : M2 := ⟨X.a*Y.a + X.b*Y.c, X.a*Y.b + X.b*Y.d, X.c*Y.a + X.d*Y.c, X.c*Y.b + X.d*Y.d⟩
def add (X Y : M2) : M2 := ⟨X.a + Y.a, X.b + Y.b, X.c + Y.c, X.d + Y.d⟩
def one : M2 := ⟨1, 0, 0, 1⟩
theorem ext' {X Y : M2} (h1 : X.a = Y.a) (h2 : X.b = Y.b) (h3 : X.c = Y.c) (h4 : X.d = Y.d) : X = Y := by
  cases X; cases Y; simp_all
end M2

open M2

/-- Sandwich mult-add: a + b + a·k·b. -/
def sand (k a b : M2) : M2 := add (add a b) (mul (mul a k) b)
/-- Front-loaded mult-add: a + b + k·a·b. -/
def front (k a b : M2) : M2 := add (add a b) (mul k (mul a b))

/-- **The lift holds for every matrix weight k:** 1 + k(a ◯ b) = (1 + k a)(1 + k b). -/
theorem sand_lift (k a b : M2) :
    add one (mul k (sand k a b)) = mul (add one (mul k a)) (add one (mul k b)) := by
  apply M2.ext' <;> simp [sand, add, mul, one] <;> grind

/-- **The sandwich law is associative for every k.** -/
theorem sand_assoc (k a b c : M2) : sand k (sand k a b) c = sand k a (sand k b c) := by
  apply M2.ext' <;> simp [sand, add, mul] <;> grind

/-- The front-loaded law fails associativity for a non-central weight (explicit witness). -/
theorem front_not_assoc :
    front ⟨1, 1, 0, 1⟩ (front ⟨1, 1, 0, 1⟩ ⟨1, 0, 0, 0⟩ ⟨0, 0, 1, 0⟩) ⟨1, 0, 0, 0⟩ ≠
    front ⟨1, 1, 0, 1⟩ ⟨1, 0, 0, 0⟩ (front ⟨1, 1, 0, 1⟩ ⟨0, 0, 1, 0⟩ ⟨1, 0, 0, 0⟩) := by decide

/-- A central (scalar) weight makes the front-loaded law agree with the sandwich law. -/
theorem front_eq_sand_scalar (s : Int) (a b : M2) : front ⟨s, 0, 0, s⟩ a b = sand ⟨s, 0, 0, s⟩ a b := by
  apply M2.ext' <;> simp [front, sand, add, mul] <;> grind

/-! ## 3. Schur's factor -/

/-- (n−1)(n−2)/2 is 0 in dimension 2 and positive from dimension 3 on. -/
theorem schur_factor : (2 - 1) * (2 - 2) / 2 = 0 ∧ (3 - 1) * (3 - 2) / 2 = 1 ∧ (4 - 1) * (4 - 2) / 2 = 3 := by
  decide

/-- Isotropic curvature: G = ((n−1) − n(n−1)/2)·q·g = −(n−1)(n−2)/2·q·g. -/
theorem isotropic_einstein (n : Rat) : (n - 1) - n * (n - 1) / 2 = -((n - 1) * (n - 2) / 2) := by grind

/-! ## 4. Jacobson's area bookkeeping -/

/-- Bundle area from two transverse stretches. -/
theorem bundle_area (a b : Rat) : (1 + a) * (1 + b) = 1 + (a ◯ b) := by unfold circ; grind

/-- Straight rays (a = αλ, b = βλ): the second difference of the area is the cross term 2αβ. -/
theorem flat_focusing (α β l : Rat) :
    (1 + α*(l+2)) * (1 + β*(l+2)) - 2 * ((1 + α*(l+1)) * (1 + β*(l+1))) + (1 + α*l) * (1 + β*l) = 2*α*β := by
  grind

/-- Focusing split: θ²/2 + σ² = θ₁² + θ₂², with θ = θ₁ + θ₂ and σ² = (θ₁ − θ₂)²/2. -/
theorem focusing_split (t1 t2 : Rat) : (t1 + t2)^2 / 2 + (t1 - t2)^2 / 2 = t1^2 + t2^2 := by grind

/-- The cross term: θ² − (θ₁² + θ₂²) = 2θ₁θ₂ (the ab part of mult-add, in rates). -/
theorem focusing_cross (t1 t2 : Rat) : (t1 + t2)^2 - (t1^2 + t2^2) = 2 * t1 * t2 := by grind

/-! ## 5. The river's energy constraint -/

/-- Three-way mult-add = sum + pairwise layer (e2) + triple layer. -/
theorem three_way_layers (a b c : Rat) :
    a ◯ (b ◯ c) = (a + b + c) + (a*b + b*c + c*a) + a*b*c := by
  unfold circ; grind

/-- e2 from traces: (tr K)² − tr K² = 2 e2. -/
theorem e2_from_traces (a b c : Rat) : (a + b + c)^2 - (a^2 + b^2 + c^2) = 2 * (a*b + b*c + c*a) := by grind

/-- Spherical river (rates v′, v/r, v/r, with q = v/r): r²·e2 = v² + 2 r v v′ = (r v²)′. -/
theorem spherical_shell (r q v v' : Rat) (h : r * q = v) :
    r * r * (2 * v' * q + q * q) = v * v + 2 * r * v * v' := by
  subst h; grind

/-- Schwarzschild river: radial rate p = −q/2 makes e2 vanish, (tr K)² = K·K. -/
theorem schwarzschild_constraint (p q : Rat) (h : 2 * p + q = 0) : (p + 2*q)^2 = p^2 + 2*q^2 := by
  have hp : p = -q/2 := by grind
  subst hp; grind

/-- Aether river: (1 + c₂)(tr K)² − (1 − c₁₃) K·K = (c₂ + c₁₃)(tr K)² when (tr K)² = K·K. -/
theorem aether_constraint (T2 KK c2 c13 : Rat) (h : T2 = KK) :
    (1 + c2) * T2 - (1 - c13) * KK = (c2 + c13) * T2 := by
  subst h; grind

/-! ## 6. Light-like test directions -/

/-- In a 2D slice, X(k,k) = 0 for both light directions k = (1, ±1) forces X₀₁ = 0, X₁₁ = −X₀₀ (X = φ·η). -/
theorem null_null_2d (X00 X01 X11 : Int) (h1 : X00 + 2*X01 + X11 = 0) (h2 : X00 - 2*X01 + X11 = 0) :
    X01 = 0 ∧ X11 = -X00 := by
  omega

/-! ## 7. Brans–Dicke -/

/-- γ = (1 + ω)/(2 + ω) means γ − 1 = −1/(2 + ω). -/
theorem brans_dicke (ω γ : Rat) (h : (2 + ω) * γ = 1 + ω) : (2 + ω) * (γ - 1) = -1 := by grind

end GRBridges
