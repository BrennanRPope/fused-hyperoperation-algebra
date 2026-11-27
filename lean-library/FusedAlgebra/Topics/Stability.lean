/-
  Fused Hyperoperation Algebra: antisymmetry (A1) from stability
  (topics/why-i/rotations-force-i, paper Proposition "antisymmetry from stability").
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/Stability.lean.

  Without (A1), Planck–Einstein (A3) and the shift rule (A4) still give J² = −1, [J,H] = 0 and
  B = JH + E with E² = 0 (EnergyShift.lean never uses the transpose). Here:
  1. shifting the energy zero to a level λ leaves the generator B − λJ = E + J(H − λ), which is E on
     the λ-eigenspace of H;
  2. E² = 0 makes the flow of E the straight line 1 + tE (a group law in t), which grows without
     bound unless E = 0, so stability forces E = 0 and A(H) = JH/ħ;
  3. if Jᵀ = a + bJ with a, b scalars and (Jᵀ)² = −1, then (a² − b² + 1) + 2ab·J = 0;
     since 1 and J are independent this gives a = 0, b = ±1, and Jᵀ = J is impossible for a real
     symmetric J with J² = −1, so Jᵀ = −J and every A(H) = JH/ħ is antisymmetric.
  Operators live in an abstract ring R; scalars a, b, t, λ are elements of a commutative subring
  that commutes with everything, given here as hypotheses.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.Stability
  (namespace CircleOp.Stability, the library's convention; no proof changed).
-/
namespace CircleOp.Stability
open Lean.Grind

variable {R : Type} [Ring R]

/-- Shifting to level λ: the generator B − λJ is E plus J(H − λ), and J(H − λ) vanishes on the
    λ-eigenspace of H. -/
theorem shift_to_level (J H E lam : R) (hl : lam * J = J * lam) :
    (J * H + E) - lam * J = E + J * (H - lam) := by
  rw [Ring.sub_eq_add_neg H lam, Semiring.left_distrib, Ring.mul_neg, hl]; grind

/-- With E² = 0 the flow of E is t ↦ 1 + tE: it composes as a group in t, so it is the exponential
    of tE, and it grows linearly unless E = 0. -/
theorem nilpotent_flow (E s t : R) (hE : E * E = 0) (ht : t * E = E * t) :
    (1 + s * E) * (1 + t * E) = 1 + (s + t) * E := by
  have : s * E * (t * E) = s * t * (E * E) := by
    rw [Semiring.mul_assoc, ← Semiring.mul_assoc E t E, ← ht, Semiring.mul_assoc,
      ← Semiring.mul_assoc]
  rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib, this, hE]
  rw [Semiring.mul_one, Semiring.one_mul, Semiring.right_distrib]
  have z : s * t * (0 : R) = 0 := by grind
  rw [z]; grind

/-- If Jᵀ = a + bJ with central scalars, squaring (Jᵀ)² = −1 gives (a² − b² + 1) + 2ab·J = 0. -/
theorem transpose_turn (J a b : R) (hJ : J * J = -1) (ha : a * J = J * a) (hb : b * J = J * b)
    (hab : a * b = b * a) :
    (a + b * J) * (a + b * J) + 1 = (a * a - b * b + 1) + (a * b + a * b) * J := by
  have e1 : (a + b * J) * (a + b * J) = a * a + a * (b * J) + b * J * a + b * J * (b * J) := by
    rw [Semiring.left_distrib, Semiring.right_distrib, Semiring.right_distrib]; grind
  have e2 : b * J * a = a * b * J := by
    rw [Semiring.mul_assoc, ← ha, ← Semiring.mul_assoc, hab]
  have e3 : b * J * (b * J) = -(b * b) := by
    rw [Semiring.mul_assoc, ← Semiring.mul_assoc J b J, ← hb, Semiring.mul_assoc,
      ← Semiring.mul_assoc, hJ, Ring.mul_neg, Semiring.mul_one]
  have e4 : a * (b * J) = a * b * J := (Semiring.mul_assoc _ _ _).symm
  rw [e1, e2, e3, e4, Semiring.right_distrib (a * b) (a * b) J]
  grind

/-- With Jᵀ = −J, the generator JH is antisymmetric: (JH)ᵀ = H Jᵀ = −HJ = −JH. -/
theorem JH_antisymm (J H Jt : R) (hJt : Jt = -J) (hJH : J * H = H * J) :
    H * Jt = -(J * H) := by
  rw [hJt, Ring.mul_neg, hJH]

end CircleOp.Stability
