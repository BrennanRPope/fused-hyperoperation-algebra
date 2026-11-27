/-
  Fused Hyperoperation Algebra: the "only energy differences matter" theorem, algebraic core
  (topics/why-i/rotations-force-i, paper Theorem 1).
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/EnergyShift.lean.

  Works in any ring R (not necessarily commutative), e.g. real n×n matrices.
  J := ħA(1), B := ħA(H). Applying the Planck–Einstein rule (A3) to H + c for every real c,
  with the shift rule (A4), gives three equations (powers c², c¹, c⁰):
      J² = −1,   BJ + JB = −2H,   B² = −H².
  From these alone, with no use of energy conservation (A2):
    * J commutes with H, and B commutes with H, so (A2) is a consequence;
    * E := B − JH anticommutes with J, commutes with H, and squares to 0
      (steps (iv) of the proof; step (v), E antisymmetric with E² = 0 ⇒ E = 0, uses the
      transpose and trace and is left to the paper).
  "No 2-torsion" (x + x = 0 → x = 0) holds for real matrices.
-/
namespace EnergyShift
open Lean.Grind

variable {R : Type} [Ring R]

theorem dbl_sub (a b : R) (h : -(a + a) = -(b + b)) : (a - b) + (a - b) = 0 := by grind
theorem eq_of_sub (a b : R) (h : a - b = 0) : a = b := by grind
theorem add_comm' (a b : R) : a + b = b + a := by grind

/-- The c-expansion of (B + cJ)² = −(H + c)² is its three coefficient equations.
    Stated for c = 0, 1, −1, which already determine the three coefficients. -/
theorem shift_square_coeffs (J B H : R) (hJ : J * J = -1) (hB : B * B = -(H * H))
    (hBJ : B * J + J * B = -(H + H)) (c : R) (hc : c * J = J * c) (hcH : c * H = H * c)
    (hcB : c * B = B * c) :
    (B + c * J) * (B + c * J) = -((H + c) * (H + c)) := by
  have e1 : (B + c * J) * (B + c * J) = B * B + c * (B * J + J * B) + c * c * (J * J) := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]
    have : c * J * B = c * (J * B) := Semiring.mul_assoc _ _ _
    have h2 : B * (c * J) = c * (B * J) := by
      rw [← Semiring.mul_assoc, ← hcB, Semiring.mul_assoc]
    have h3 : c * J * (c * J) = c * c * (J * J) := by
      rw [Semiring.mul_assoc, ← Semiring.mul_assoc J c J, ← hc, Semiring.mul_assoc,
        ← Semiring.mul_assoc]
    rw [this, h2, h3, Semiring.left_distrib]; grind
  rw [e1, hB, hBJ, hJ]
  rw [Semiring.left_distrib (H + c), Semiring.right_distrib, Semiring.right_distrib, ← hcH]
  rw [Ring.mul_neg, Ring.mul_neg, Semiring.mul_one, Semiring.left_distrib]
  grind

/-- (A3) at the shifts c = 0, 1, −1 gives its three coefficient equations. -/
theorem sq_plus (B J : R) : (B + J) * (B + J) = B * B + (B * J + J * B) + J * J := by
  rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]; grind
theorem sq_minus (B J : R) : (B - J) * (B - J) = B * B - (B * J + J * B) + J * J := by
  rw [Ring.sub_eq_add_neg, Ring.sub_eq_add_neg, Semiring.left_distrib, Semiring.right_distrib,
    Semiring.right_distrib, Ring.mul_neg, Ring.neg_mul, Ring.mul_neg, Ring.neg_mul]; grind
theorem lin_a (q x a h : R) (hp : -q + x + a = -(q + (h + h) + 1))
    (hm : -q - x + a = -(q - (h + h) + 1)) : (a + 1) + (a + 1) = 0 := by grind
theorem lin_x (q x a h : R) (hp : -q + x + a = -(q + (h + h) + 1))
    (hm : -q - x + a = -(q - (h + h) + 1)) : (x + (h + h)) + (x + (h + h)) = 0 := by grind
theorem eq_neg_of_add (a b : R) (h : a + b = 0) : a = -b := by grind
/-- Reading the c², c¹ coefficients of (B + cJ)² = −(H + c)² off the shifts c = ±1. -/
theorem coeffs_of_shifts (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (h0 : B * B = -(H * H))
    (hp : (B + J) * (B + J) = -((H + 1) * (H + 1)))
    (hm : (B - J) * (B - J) = -((H - 1) * (H - 1))) :
    J * J = -1 ∧ B * J + J * B = -(H + H) := by
  rw [sq_plus, sq_plus, h0] at hp
  rw [sq_minus, sq_minus, h0] at hm
  simp only [Semiring.mul_one, Semiring.one_mul] at hp hm
  have ha := no2 _ (lin_a _ _ _ _ hp hm)
  have hx := no2 _ (lin_x _ _ _ _ hp hm)
  exact ⟨eq_neg_of_add _ _ ha, eq_neg_of_add _ _ hx⟩

/-- Step 1: J commutes with any combination BJ + JB, given J² = −1. -/
theorem J_comm_anticomm (J B : R) (hJ : J * J = -1) :
    J * (B * J + J * B) = (B * J + J * B) * J := by
  have e1 : J * (B * J + J * B) = J * B * J + (J * J) * B := by grind
  have e2 : (B * J + J * B) * J = B * (J * J) + J * B * J := by grind
  rw [e1, e2, hJ]; grind

/-- Energy conservation for the common rotation: [J, H] = 0, from (A3) and (A4) alone. -/
theorem J_comm_H (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (hJ : J * J = -1) (hBJ : B * J + J * B = -(H + H)) : J * H = H * J := by
  have h := J_comm_anticomm J B hJ
  rw [hBJ, Ring.mul_neg, Ring.neg_mul, Semiring.left_distrib, Semiring.right_distrib] at h
  exact eq_of_sub _ _ (no2 _ (dbl_sub _ _ h))

/-- Squares commute with what their roots commute with. -/
theorem comm_sq (J H : R) (h : J * H = H * J) : H * H * J = J * (H * H) := by
  calc H * H * J = H * (H * J) := Semiring.mul_assoc _ _ _
    _ = H * (J * H) := by rw [h]
    _ = (H * J) * H := (Semiring.mul_assoc _ _ _).symm
    _ = (J * H) * H := by rw [h]
    _ = J * (H * H) := Semiring.mul_assoc _ _ _

/-- Energy conservation (A2), [B, H] = 0, from (A3) and (A4) alone. -/
theorem A2_from_A3_A4 (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (hJ : J * J = -1) (hB : B * B = -(H * H)) (hBJ : B * J + J * B = -(H + H)) :
    B * H = H * B := by
  have hJH := J_comm_H J B H no2 hJ hBJ
  have h : B * (B * J + J * B) = (B * J + J * B) * B := by
    have e1 : B * (B * J + J * B) = (B * B) * J + B * J * B := by grind
    have e2 : (B * J + J * B) * B = B * J * B + J * (B * B) := by grind
    rw [e1, e2, hB, Ring.mul_neg, Ring.neg_mul, comm_sq J H hJH, add_comm']
  rw [hBJ, Ring.mul_neg, Ring.neg_mul, Semiring.left_distrib, Semiring.right_distrib] at h
  exact eq_of_sub _ _ (no2 _ (dbl_sub _ _ h))

/-- Step (iv): E := B − JH anticommutes with J. -/
theorem E_anticomm_J (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (hJ : J * J = -1) (hBJ : B * J + J * B = -(H + H)) :
    J * (B - J * H) + (B - J * H) * J = 0 := by
  have hJH := J_comm_H J B H no2 hJ hBJ
  have a : J * (J * H) = -H := by rw [← Semiring.mul_assoc, hJ, Ring.neg_mul, Semiring.one_mul]
  have b : J * H * J = -H := by rw [Semiring.mul_assoc, ← hJH, ← Semiring.mul_assoc, hJ,
    Ring.neg_mul, Semiring.one_mul]
  have e : J * (B - J * H) + (B - J * H) * J = (B * J + J * B) - (J * (J * H) + J * H * J) := by
    rw [Ring.sub_eq_add_neg, Ring.sub_eq_add_neg, Semiring.left_distrib, Semiring.right_distrib,
      Ring.mul_neg, Ring.neg_mul]; grind
  rw [e, hBJ, a, b]; grind

/-- Step (iv): E := B − JH commutes with H. -/
theorem E_comm_H (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (hJ : J * J = -1) (hB : B * B = -(H * H)) (hBJ : B * J + J * B = -(H + H)) :
    (B - J * H) * H = H * (B - J * H) := by
  have hJH := J_comm_H J B H no2 hJ hBJ
  have hBH := A2_from_A3_A4 J B H no2 hJ hB hBJ
  have : J * H * H = H * (J * H) := by
    rw [Semiring.mul_assoc, ← Semiring.mul_assoc H J H, ← hJH, Semiring.mul_assoc]
  rw [Ring.sub_eq_add_neg, Semiring.right_distrib, Semiring.left_distrib, Ring.neg_mul,
    Ring.mul_neg, this, hBH]

/-- Step (iv): E := B − JH squares to 0. With E antisymmetric, tr(EᵀE) = 0 gives E = 0. -/
theorem E_sq_zero (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (hJ : J * J = -1) (hB : B * B = -(H * H)) (hBJ : B * J + J * B = -(H + H)) :
    (B - J * H) * (B - J * H) = 0 := by
  have hJH := J_comm_H J B H no2 hJ hBJ
  have hEJ := E_anticomm_J J B H no2 hJ hBJ
  have hEH := E_comm_H J B H no2 hJ hB hBJ
  -- B = JH + E, so B² = (JH)² + JH·E + E·JH + E² and the cross terms cancel.
  generalize hE : B - J * H = E at hEJ hEH ⊢
  have hBdef : B = J * H + E := by rw [← hE]; grind
  have JHsq : J * H * (J * H) = -(H * H) := by
    rw [Semiring.mul_assoc, ← Semiring.mul_assoc H J H, ← hJH, Semiring.mul_assoc,
      ← Semiring.mul_assoc J J, hJ, Ring.neg_mul, Semiring.one_mul]
  have EJ : E * J = -(J * E) := by grind
  have cross : J * H * E + E * (J * H) = 0 := by
    rw [← Semiring.mul_assoc E J H, EJ, Ring.neg_mul, Semiring.mul_assoc J H E,
      Semiring.mul_assoc J E H, hEH.symm]
    grind
  have expand : B * B = J * H * (J * H) + (J * H * E + E * (J * H)) + E * E := by
    rw [hBdef, Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]; grind
  rw [hB, JHsq, cross] at expand
  grind

/-- Headline: energy conservation (A2) is not an assumption. Planck–Einstein (A3) applied to
    H, H + 1 and H − 1, with the shift rule (A4) giving ħA(H ± 1) = B ± J, forces [B, H] = 0. -/
theorem A2_redundant (J B H : R) (no2 : ∀ x : R, x + x = 0 → x = 0)
    (h0 : B * B = -(H * H))
    (hp : (B + J) * (B + J) = -((H + 1) * (H + 1)))
    (hm : (B - J) * (B - J) = -((H - 1) * (H - 1))) :
    J * H = H * J ∧ B * H = H * B := by
  have ⟨hJ, hBJ⟩ := coeffs_of_shifts J B H no2 h0 hp hm
  exact ⟨J_comm_H J B H no2 hJ hBJ, A2_from_A3_A4 J B H no2 hJ h0 hBJ⟩

end EnergyShift
