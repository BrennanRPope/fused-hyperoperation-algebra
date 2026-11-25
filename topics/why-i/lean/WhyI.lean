/-
  Fused Hyperoperation Algebra: why quantum mechanics needs i (topics/why-i).
  Core Lean 4 only (no mathlib), no `sorry`, no imports.
  Library path: FusedAlgebra/Topics/WhyI.lean (compiled, 21 theorems, standard axioms only).

  Everything here is a small integer identity behind a picture on the page:
  1. the twist (order gap) is antisymmetric, and between two mult-add dials it is (k − k₂)·a·b;
  2. the dial-k elements stay on the line λ = 1 + k·s through the identity (the half-plane picture);
  3. a quarter turn J squares to −1, keeps length, and turns the dot product into the area form;
  4. a mirror flips the area form and turns J into −J (the twist's sign is i's sign);
  5. two real symmetric 2×2 matrices have an antisymmetric commutator (Stueckelberg's point);
  6. path integrals: a chain of weighted-circle steps is the product of (1 + k·aⱼ) (two and three
     steps; over Int the steps commute, the ordered matrix version is check P1), and a path's
     accumulated twist (sum of pairwise wedges) is its shoelace area (three steps; check P8).
-/
namespace WhyI

/-! ## 1. The twist -/

/-- Order gap of (a, c) and (b, d), where (a, c) is x ↦ (1 + a)x + c (topics/order-dial). -/
def wedge (a c b d : Int) : Int := a * d - b * c

theorem wedge_antisymm (a c b d : Int) : wedge a c b d = - wedge b d a c := by
  unfold wedge; grind

theorem wedge_self (a c : Int) : wedge a c a c = 0 := by
  unfold wedge; grind

/-- Mult-add element a at dial k is the pair (k·a, a). Two dials twist by (k − k₂)·a·b. -/
theorem dial_twist (k k₂ a b : Int) : wedge (k * a) a (k₂ * b) b = (k - k₂) * a * b := by
  unfold wedge; grind

theorem same_dial_no_twist (k a b : Int) : wedge (k * a) a (k * b) b = 0 := by
  unfold wedge; grind

/-! ## 2. The dial line in the half-plane: λ = 1 + k·s is closed under composition -/

/-- Composing (λ₁, s₁) after (λ₂, s₂) gives (λ₁λ₂, λ₁s₂ + s₁). On the dial-k line it stays on the line. -/
theorem dial_line_closed (k a b : Int) :
    (1 + k * a) * (1 + k * b) = 1 + k * ((1 + k * a) * b + a) := by
  grind

/-- ... and the shift part of the composite is mult-add. -/
theorem dial_line_multAdd (k a b : Int) : (1 + k * a) * b + a = a + b + k * a * b := by
  grind

/-! ## 3. Twist plus length gives a quarter turn J with J² = −1 -/

def omega (u₁ u₂ v₁ v₂ : Int) : Int := u₁ * v₂ - u₂ * v₁
def dot (u₁ u₂ v₁ v₂ : Int) : Int := u₁ * v₁ + u₂ * v₂
/-- J(x, y) = (−y, x), written by components. -/
def J₁ (_x y : Int) : Int := -y
def J₂ (x _y : Int) : Int := x

theorem J_sq_one (x y : Int) : J₁ (J₁ x y) (J₂ x y) = -x := by
  unfold J₁ J₂; grind
theorem J_sq_two (x y : Int) : J₂ (J₁ x y) (J₂ x y) = -y := by
  unfold J₁ J₂; grind

theorem J_keeps_length (x y : Int) : dot (J₁ x y) (J₂ x y) (J₁ x y) (J₂ x y) = dot x y x y := by
  unfold dot J₁ J₂; grind

/-- The compatible-triple relation ω(u, v) = g(Ju, v). -/
theorem omega_eq_dot_J (u₁ u₂ v₁ v₂ : Int) : omega u₁ u₂ v₁ v₂ = dot (J₁ u₁ u₂) (J₂ u₁ u₂) v₁ v₂ := by
  unfold omega dot J₁ J₂; grind

/-- J is "positive" for the orientation: ω(u, Ju) = |u|² ≥ 0. -/
theorem omega_u_Ju (x y : Int) : omega x y (J₁ x y) (J₂ x y) = x * x + y * y := by
  unfold omega J₁ J₂; grind

/-! ## 4. A mirror flips the twist and turns J into −J -/

/-- Mirror R(x, y) = (x, −y). -/
theorem mirror_flips_omega (u₁ u₂ v₁ v₂ : Int) : omega u₁ (-u₂) v₁ (-v₂) = - omega u₁ u₂ v₁ v₂ := by
  unfold omega; grind

/-- R J R = −J, by components: R(J(R(x, y))) = (−y, −x)·… = −J(x, y). -/
theorem mirror_conj_J_one (x y : Int) : J₁ x (-y) = - J₁ x y := by
  unfold J₁; grind
theorem mirror_conj_J_two (x y : Int) : -(J₂ x (-y)) = - J₂ x y := by
  unfold J₂; grind

/-! ## 5. Stueckelberg: real symmetric observables have an antisymmetric order cost -/

/-- Commutator entries of A = [[p, q], [q, r]] and B = [[s, t], [t, w]]. -/
def comm11 (p q _r s t _w : Int) : Int := (p * s + q * t) - (s * p + t * q)
def comm22 (_p q r _s t w : Int) : Int := (q * t + r * w) - (t * q + w * r)
def comm12 (p q r s t w : Int) : Int := (p * t + q * w) - (s * q + t * r)
def comm21 (p q r s t w : Int) : Int := (q * s + r * t) - (t * p + w * q)

theorem comm_diag_zero (p q r s t w : Int) : comm11 p q r s t w = 0 ∧ comm22 p q r s t w = 0 := by
  unfold comm11 comm22; constructor <;> grind

theorem comm_antisymm (p q r s t w : Int) : comm21 p q r s t w = - comm12 p q r s t w := by
  unfold comm12 comm21; grind

/-- So the commutator is c·J with c = comm12, and (c·J)² = −c²: the product of its two
    off-diagonal entries is −c², the determinant-side fact behind eigenvalues ±ic. -/
theorem comm_square_negative (p q r s t w : Int) :
    comm12 p q r s t w * comm21 p q r s t w = -(comm12 p q r s t w * comm12 p q r s t w) := by
  rw [comm_antisymm]; grind

/-! ## 6. Path integrals: whole-chain composition and the area a path sweeps -/

theorem chain_two (k a b : Int) : 1 + k * (a + b + k * a * b) = (1 + k * a) * (1 + k * b) := by
  grind

theorem chain_three (k a b c : Int) :
    1 + k * (c + (b + a + k * b * a) + k * c * (b + a + k * b * a))
      = (1 + k * c) * (1 + k * b) * (1 + k * a) := by
  grind

/-- k = −1 is survival analysis: the chance of the event over two periods is 1 − (1 − h₁)(1 − h₂). -/
theorem survival_two (h₁ h₂ : Int) : h₁ + h₂ + (-1) * h₁ * h₂ = 1 - (1 - h₁) * (1 - h₂) := by
  grind

/-- Three steps v₁, v₂, v₃ from the origin: the sum of the pairwise twists equals the shoelace
    sum of the polygon 0, v₁, v₁+v₂, v₁+v₂+v₃ (twice its signed area). -/
theorem levy_three (x₁ y₁ x₂ y₂ x₃ y₃ : Int) :
    omega x₁ y₁ x₂ y₂ + omega x₁ y₁ x₃ y₃ + omega x₂ y₂ x₃ y₃
      = omega x₁ y₁ (x₁ + x₂) (y₁ + y₂) + omega (x₁ + x₂) (y₁ + y₂) (x₁ + x₂ + x₃) (y₁ + y₂ + y₃) := by
  unfold omega; grind

end WhyI
