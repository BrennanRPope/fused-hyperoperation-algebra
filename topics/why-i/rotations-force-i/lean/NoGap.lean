/-
  Fused Hyperoperation Algebra: composites rule out the gap
  (topics/why-i/rotations-force-i, paper Proposition "Planck–Einstein from the energy zero").
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/NoGap.lean.

  On the joint space, P = J_A ⊗ 1 and Q = 1 ⊗ J_B commute and square to −1; the balanced
  product W is where P w = Q w. The gap terms X = M_A ⊗ 1 and Y = 1 ⊗ M_B anticommute with
  their own side's quarter turn and commute with the other's. Elements of a ring stand for
  operators, and a vector w for the operator it is applied through (products compose).
    * `gap_to_odd`: each gap term sends W into W⁻ (where P = −Q).
    * `gap_sum_zero`: if X + Y maps W into W (it is part of the joint generator, which acts on W),
      then (X + Y) w = 0 on W.
    * `gap_zero`: then X w = 0 on W, so M_A ⊗ 1 vanishes on W, and by faithfulness M_A = 0.
-/
namespace NoGap
open Lean.Grind

variable {R : Type} [Ring R]

theorem two_zero (no2 : ∀ x : R, x + x = 0 → x = 0) (a : R) (h : a = -a) : a = 0 :=
  no2 a (by grind)

/-- A term anticommuting with P and commuting with Q sends W into W⁻. -/
theorem gap_to_odd (P Q X w : R) (hXP : X * P = -(P * X)) (hXQ : X * Q = Q * X)
    (hw : P * w = Q * w) : P * (X * w) = -(Q * (X * w)) := by
  have e1 : P * (X * w) = -(X * (P * w)) := by
    rw [← Semiring.mul_assoc, ← Semiring.mul_assoc]
    have : P * X = -(X * P) := by rw [hXP]; grind
    rw [this, Ring.neg_mul]
  rw [e1, hw, ← Semiring.mul_assoc, hXQ, Semiring.mul_assoc]

/-- A term anticommuting with Q and commuting with P also sends W into W⁻. -/
theorem gap_to_odd' (P Q Y w : R) (hYQ : Y * Q = -(Q * Y)) (hYP : Y * P = P * Y)
    (hw : P * w = Q * w) : P * (Y * w) = -(Q * (Y * w)) := by
  have e1 : Q * (Y * w) = -(Y * (Q * w)) := by
    rw [← Semiring.mul_assoc, ← Semiring.mul_assoc]
    have : Q * Y = -(Y * Q) := by rw [hYQ]; grind
    rw [this, Ring.neg_mul]
  have e2 : Y * (P * w) = P * (Y * w) := by
    rw [← Semiring.mul_assoc, hYP, Semiring.mul_assoc]
  rw [e1, ← hw, e2]; grind

theorem cancel_turn (Q v : R) (hQ : Q * Q = -1) (h : Q * v = 0) : v = 0 := by
  have : Q * (Q * v) = -v := by rw [← Semiring.mul_assoc, hQ, Ring.neg_mul, Semiring.one_mul]
  rw [h, Semiring.mul_zero] at this
  grind

/-- The gap sum vanishes on W when it maps W into W. -/
theorem gap_sum_zero (P Q X Y w : R) (no2 : ∀ x : R, x + x = 0 → x = 0) (hQ : Q * Q = -1)
    (hXP : X * P = -(P * X)) (hXQ : X * Q = Q * X) (hYQ : Y * Q = -(Q * Y)) (hYP : Y * P = P * Y)
    (hw : P * w = Q * w) (inW : P * ((X + Y) * w) = Q * ((X + Y) * w)) : (X + Y) * w = 0 := by
  have a := gap_to_odd P Q X w hXP hXQ hw
  have b := gap_to_odd' P Q Y w hYQ hYP hw
  have odd : P * ((X + Y) * w) = -(Q * ((X + Y) * w)) := by
    rw [Semiring.right_distrib, Semiring.left_distrib, Semiring.left_distrib, a, b]; grind
  rw [inW] at odd
  exact cancel_turn Q _ hQ (two_zero no2 _ odd)

/-- Composites rule out the gap: M_A ⊗ 1 vanishes on every w in W. -/
theorem gap_zero (P Q X Y w : R) (no2 : ∀ x : R, x + x = 0 → x = 0) (hP : P * P = -1)
    (hQ : Q * Q = -1) (hPQ : P * Q = Q * P)
    (hXP : X * P = -(P * X)) (hXQ : X * Q = Q * X) (hYQ : Y * Q = -(Q * Y)) (hYP : Y * P = P * Y)
    (hw : P * w = Q * w)
    (inW : ∀ v : R, P * v = Q * v → P * ((X + Y) * v) = Q * ((X + Y) * v)) : X * w = 0 := by
  have s0 := gap_sum_zero P Q X Y w no2 hQ hXP hXQ hYQ hYP hw (inW w hw)
  have hPw : P * (P * w) = Q * (P * w) := by
    have e : P * (P * w) = P * (Q * w) := by rw [hw]
    rw [e, ← Semiring.mul_assoc, ← Semiring.mul_assoc, hPQ]
  have s1 := gap_sum_zero P Q X Y (P * w) no2 hQ hXP hXQ hYQ hYP hPw (inW _ hPw)
  have yw : Y * w = -(X * w) := by
    rw [Semiring.right_distrib] at s0; grind
  have ypw : Y * (P * w) = -(X * (P * w)) := by
    rw [Semiring.right_distrib] at s1; grind
  -- X P w = −Y P w = −P Y w = P X w, while X P w = −P X w
  have c1 : X * (P * w) = P * (X * w) := by
    have : Y * (P * w) = P * (Y * w) := by
      rw [← Semiring.mul_assoc, hYP, Semiring.mul_assoc]
    rw [this, yw, Ring.mul_neg] at ypw
    grind
  have c2 : X * (P * w) = -(P * (X * w)) := by
    rw [← Semiring.mul_assoc, hXP, Ring.neg_mul, Semiring.mul_assoc]
  rw [c1] at c2
  exact cancel_turn P _ hP (two_zero no2 _ c2)

end NoGap
