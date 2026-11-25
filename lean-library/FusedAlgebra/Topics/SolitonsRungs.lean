/-
  Fused Hyperoperation Algebra: solitons on the rungs, and the gluing number q.
  Core Lean 4 only (no mathlib), no `sorry`. Style follows the shared library (core `Int`, `grind`).
  LIBRARY NOTE: copied from solitons-rungs/lean/SolitonsRungs.lean (Nov 25) and compiled here for
  the first time, with no changes beyond this note: no sorry, standard axioms only. No theorem is an
  exact twin of a library theorem; tau2_scaled is Landscape's kdv_phase_shift (CircleOp.circK_scale)
  read through tau2. (Its identities were also checked in solitons-rungs/checks.py.)

  Conventions. KdV soliton tau functions (Hirota 1971) as polynomials in the exponentials
  a = e^η₁, b = e^η₂, c = e^η₃, with pair weights A_ij = ((k_i − k_j)/(k_i + k_j))².
  Theta lattice weights use the triangular number n(n−1)/2 (here as n·(n−1), the q = p² convention).
-/
import FusedAlgebra.Base

namespace SolitonsRungs
open CircleOp

/-! ### 1. Two solitons: the tau function is weighted mult-add -/

/-- the two-soliton tau function 1 + a + b + A·ab -/
def tau2 (A a b : Int) : Int := 1 + a + b + A * a * b

theorem tau2_is_circK (A a b : Int) : tau2 A a b = 1 + circK A a b := by
  unfold tau2 circK; grind

/-- weight 1: no interaction, the tau function is the plain mult-add lift (1+a)(1+b) -/
theorem tau2_free (a b : Int) : tau2 1 a b = (1 + a) * (1 + b) := by
  unfold tau2; grind

/-- weight 0: equal speeds, the two exponentials just add (rung 0) -/
theorem tau2_merge (a b : Int) : tau2 0 a b = 1 + (a + b) := by
  unfold tau2; grind

/-- for two solitons the law is rung 1 in a rescaled coordinate: A·(τ − 1) = (A a) ◯ (A b) -/
theorem tau2_scaled (A a b : Int) : A * (tau2 A a b - 1) = circ (A * a) (A * b) := by
  unfold tau2 circ; grind

/-- read in the 1+x picture, the weight is mult-add plus an excess (A − 1) per pair -/
theorem tau2_excess (A a b : Int) : tau2 A a b = 1 + circ a b + (A - 1) * a * b := by
  unfold tau2 circ; grind

/-! ### 2. Three solitons: corners multiply over pairs -/

/-- the three-soliton tau function: every corner carries the product of its pair weights -/
def tau3 (A12 A13 A23 a b c : Int) : Int :=
  1 + a + b + c + A12 * a * b + A13 * a * c + A23 * b * c + (A12 * A13 * A23) * (a * b * c)

theorem tau3_free (a b c : Int) : tau3 1 1 1 a b c = (1 + a) * (1 + b) * (1 + c) := by
  unfold tau3; grind

/-- a soliton that interacts with nobody (weights 1) factors off -/
theorem tau3_splits (A12 a b c : Int) : tau3 A12 1 1 a b c = tau2 A12 a b * (1 + c) := by
  unfold tau3 tau2; grind

/-- equal weights k: the triple corner is k³ (three pairs) -/
theorem tau3_equal_corner (k a b c : Int) :
    tau3 k k k a b c = 1 + a + b + c + k * (a * b + a * c + b * c) + (k * k * k) * (a * b * c) := by
  unfold tau3; grind

/-- weighted mult-add of three has triple corner k², not k³: it is mult-add rescaled,
    while solitons weight every pair -/
theorem circK_three_corner (k a b c : Int) :
    circK k (circK k a b) c = a + b + c + k * (a * b + a * c + b * c) + (k * k) * (a * b * c) := by
  unfold circK; grind

/-- the two agree exactly when k³ = k², i.e. k = 0 or k = 1 (rung 0 or rung 1) -/
theorem soliton_vs_circK_corner (k : Int) : k * k * k = k * k ↔ k = 0 ∨ k = 1 := by
  have h : k * k * k - k * k = (k * k) * (k - 1) := by grind
  constructor
  · intro e
    have z : (k * k) * (k - 1) = 0 := by omega
    cases Int.mul_eq_zero.mp z with
    | inl h0 =>
      left
      cases Int.mul_eq_zero.mp h0 with
      | inl h1 => exact h1
      | inr h1 => exact h1
    | inr h1 => right; omega
  · intro hk
    cases hk with
    | inl h0 => subst h0; decide
    | inr h1 => subst h1; decide

/-! ### 3. The rung flow's corners versus the soliton's corners

  Rung s, n inputs: the corner (multilinear) coefficient is n!·[yⁿ] E^s(y):
  n = 2: s,  n = 3: s(3s−1)/2,  n = 4: s(2s−1)(3s−1)/2  (computed in ../checks.py).
  Solitons with equal weight A: A, A³, A⁶ (one A per pair). -/

/-- three bodies: 2·(soliton corner − rung corner) = s(2s−1)(s−1): the roads meet at s = 0, ½, 1 -/
theorem three_body_roads (s : Int) : 2 * (s * s * s) - s * (3 * s - 1) = s * (2 * s - 1) * (s - 1) := by
  grind

/-- the rung flow's three-body term: triple corner minus (pair)², times 2, is s(s−1) -/
theorem rung_three_body (s : Int) : s * (3 * s - 1) - 2 * (s * s) = s * (s - 1) := by
  grind

/-- on whole numbers the rung flow has no three-body term only at rungs 0 and 1 -/
theorem rung_three_body_zero (s : Int) : s * (s - 1) = 0 ↔ s = 0 ∨ s = 1 := by
  constructor
  · intro e
    cases Int.mul_eq_zero.mp e with
    | inl h => left; exact h
    | inr h => right; omega
  · intro hs
    cases hs with
    | inl h => subst h; decide
    | inr h => subst h; decide

/-- four bodies: 2·(soliton corner A⁶ − rung corner) factors as s(s−1)·(quartic) -/
theorem four_body_roads (s : Int) :
    2 * (s * s * s * s * s * s) - s * (2 * s - 1) * (3 * s - 1)
      = s * (s - 1) * (2 * s * s * s * s + 2 * s * s * s + 2 * s * s - 4 * s + 1) := by
  grind

/-- the quartic is not zero at s = ½ (16 × its value there, written with S = 2s = 1, is −2):
    so the three-body meeting at s = ½ does not survive to four bodies -/
theorem four_body_half : 2 * 1 + 4 * 1 + 8 * 1 - 32 * 1 + 16 = (-2 : Int) := by decide

/-! ### 4. The theta lattice: only the mult-add cube survives at q = 0 -/

theorem tri_nonneg (n : Int) : 0 ≤ n * (n - 1) := by
  by_cases h : 1 ≤ n
  · exact Int.mul_nonneg (by omega) (by omega)
  · have e : n * (n - 1) = (-n) * (1 - n) := by grind
    rw [e]
    exact Int.mul_nonneg (by omega) (by omega)

/-- the lattice weight q^(n(n−1)/2) is free of q only for n = 0 and n = 1 -/
theorem tri_zero_iff (n : Int) : n * (n - 1) = 0 ↔ n = 0 ∨ n = 1 := rung_three_body_zero n

/-- two waves: the exponent n₁(n₁−1) + n₂(n₂−1) vanishes exactly on the unit square {0,1}²,
    the corners of the mult-add cube -/
theorem face_iff (n1 n2 : Int) :
    n1 * (n1 - 1) + n2 * (n2 - 1) = 0 ↔ (n1 = 0 ∨ n1 = 1) ∧ (n2 = 0 ∨ n2 = 1) := by
  have h1 := tri_nonneg n1
  have h2 := tri_nonneg n2
  constructor
  · intro e
    have z1 : n1 * (n1 - 1) = 0 := by omega
    have z2 : n2 * (n2 - 1) = 0 := by omega
    exact ⟨(tri_zero_iff n1).mp z1, (tri_zero_iff n2).mp z2⟩
  · intro hn
    have z1 := (tri_zero_iff n1).mpr hn.1
    have z2 := (tri_zero_iff n2).mpr hn.2
    omega

/-- the two-wave theta series on the lattice points 0..3 in each direction, written out:
    weight q₁^(n₁(n₁−1)/2) q₂^(n₂(n₂−1)/2) A^(n₁n₂) w₁^n₁ w₂^n₂ (exponents 0, 0, 1, 3) -/
def theta2 (q1 q2 A w1 w2 : Int) : Int :=
  let r (q : Int) (n : Nat) : Int := match n with
    | 0 => 1 | 1 => 1 | 2 => q | _ => q * q * q
  let row (n1 : Nat) : Int :=
    r q2 0 * A ^ (n1 * 0) * w2 ^ 0 + r q2 1 * A ^ (n1 * 1) * w2 ^ 1 +
    r q2 2 * A ^ (n1 * 2) * w2 ^ 2 + r q2 3 * A ^ (n1 * 3) * w2 ^ 3
  r q1 0 * w1 ^ 0 * row 0 + r q1 1 * w1 ^ 1 * row 1 + r q1 2 * w1 ^ 2 * row 2 + r q1 3 * w1 ^ 3 * row 3

/-- at q₁ = q₂ = 0 the theta series is the weighted mult-add cube 1 + (w₁ ◯[A] w₂) -/
theorem theta2_face (A w1 w2 : Int) : theta2 0 0 A w1 w2 = 1 + circK A w1 w2 := by
  unfold theta2 circK; simp <;> grind

/-- with A = 1 the two-wave theta splits into two one-wave pieces (no interaction) -/
theorem theta2_free_face (w1 w2 : Int) : theta2 0 0 1 w1 w2 = (1 + w1) * (1 + w2) := by
  rw [theta2_face]; unfold circK; grind

/-! ### 5. Gluing is a rung-1 move -/

/-- multiplying by q in the w picture is mult-add by q − 1 in the 1+x picture -/
theorem glue_is_circ (w q : Int) : w * q - 1 = circ (w - 1) (q - 1) := by
  unfold circ; grind

/-- a truncated theta product: (1+w)(1+wq)(1+vq), v standing for 1/w, is 1 + the mult-add of
    three points of w's q-orbit -/
theorem orbit_product (w v q : Int) :
    (1 + w) * (1 + w * q) * (1 + v * q) = 1 + circ w (circ (w * q) (v * q)) := by
  unfold circ; grind

/-- and at q = 0 only the soliton factor 1 + w is left -/
theorem orbit_product_q0 (w v : Int) : (1 + w) * (1 + w * 0) * (1 + v * 0) = 1 + w := by
  grind

/-! ### 6. KdV pair weights lie in [0, 1) -/

/-- for positive wave numbers, (k₁ − k₂)² < (k₁ + k₂)²: the weight A is below 1 -/
theorem kdv_weight_below_one (k1 k2 : Int) (h1 : 0 < k1) (h2 : 0 < k2) :
    (k1 - k2) * (k1 - k2) < (k1 + k2) * (k1 + k2) := by
  have e : (k1 + k2) * (k1 + k2) - (k1 - k2) * (k1 - k2) = 4 * (k1 * k2) := by grind
  have p : 0 < k1 * k2 := Int.mul_pos h1 h2
  omega

theorem kdv_weight_nonneg (k1 k2 : Int) : 0 ≤ (k1 - k2) * (k1 - k2) := by
  by_cases h : 0 ≤ k1 - k2
  · exact Int.mul_nonneg h h
  · have e : (k1 - k2) * (k1 - k2) = (k2 - k1) * (k2 - k1) := by grind
    rw [e]
    exact Int.mul_nonneg (by omega) (by omega)

/-- the weight is 0 exactly when the speeds match -/
theorem kdv_weight_zero (k1 k2 : Int) : (k1 - k2) * (k1 - k2) = 0 ↔ k1 = k2 := by
  constructor
  · intro e
    cases Int.mul_eq_zero.mp e with
    | inl h => omega
    | inr h => omega
  · intro h; subst h; grind

end SolitonsRungs
