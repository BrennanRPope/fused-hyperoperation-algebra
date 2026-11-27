/-
  Fused Hyperoperation Algebra: the energy-zero shift is unseen, so observables commute with J
  (topics/why-i/rotations-force-i, paper assumption (A0')).
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/Unseen.lean.

  A shift of the energy zero by c turns every state by exp(tcJ/ħ) = cos(tc/ħ) + sin(tc/ħ) J.
  (A0'): no measured average changes, i.e. Rᵀ O R = O for these rotations R. `Jt` stands for Jᵀ.
    * `comm_of_quarter`: invariance under the quarter turn alone (R = J) forces O J = J O.
    * `invariant_of_comm`: conversely, if O J = J O, every rotation c + sJ with c² + s² = 1
      leaves O unchanged, so (A0') says exactly that every observable commutes with J.

  First compiled in the shared library on Nov 27, as FusedAlgebra.Topics.Unseen
  (namespace CircleOp.Unseen, the library's convention; no proof changed).
-/
namespace CircleOp.Unseen
open Lean.Grind

variable {R : Type} [Ring R]

/-- Quarter-turn invariance JᵀOJ = O, with Jᵀ = −J and J² = −1, gives [J, O] = 0. -/
theorem comm_of_quarter (J Jt O : R) (hJ : J * J = -1) (hJt : Jt = -J)
    (hinv : Jt * O * J = O) : O * J = J * O := by
  rw [hJt] at hinv
  have h : J * (-J * O * J) = J * O := by rw [hinv]
  have e : J * (-J * O * J) = O * J := by
    rw [Ring.neg_mul, Ring.neg_mul, Ring.mul_neg, ← Semiring.mul_assoc, ← Semiring.mul_assoc, hJ,
      Ring.neg_mul, Ring.neg_mul, Semiring.one_mul]
    grind
  rw [e] at h; exact h

/-- Conversely, if O commutes with J, the rotation c + sJ (c, s central, c² + s² = 1) leaves O fixed. -/
theorem invariant_of_comm (J O c s : R) (hJ : J * J = -1) (hc : ∀ x : R, c * x = x * c)
    (hs : ∀ x : R, s * x = x * s) (hcs : c * c + s * s = 1) (hOJ : O * J = J * O) :
    (c - s * J) * O * (c + s * J) = O := by
  have e1 : (c - s * J) * O * (c + s * J)
      = c * c * O + c * s * (O * J - J * O) - s * s * (J * O * J) := by
    have a1 : s * J * O = s * (J * O) := Semiring.mul_assoc _ _ _
    have a2 : O * (s * J) = s * (O * J) := by
      rw [← Semiring.mul_assoc, ← hs O, Semiring.mul_assoc]
    have a3 : s * J * O * (s * J) = s * s * (J * O * J) := by
      rw [Semiring.mul_assoc (s * J) O, Semiring.mul_assoc s J, ← Semiring.mul_assoc O s J, ← hs O,
        Semiring.mul_assoc s O J, ← Semiring.mul_assoc J s, ← hs J, Semiring.mul_assoc s J,
        ← Semiring.mul_assoc s s, ← Semiring.mul_assoc J O J]
    have a4 : c * O * c = c * c * O := by rw [Semiring.mul_assoc, ← hc O, Semiring.mul_assoc]
    have a5 : c * O * (s * J) = c * s * (O * J) := by
      rw [Semiring.mul_assoc, a2, ← Semiring.mul_assoc]
    have a6 : s * J * O * c = c * s * (J * O) := by
      rw [a1, Semiring.mul_assoc, ← hc (J * O), ← Semiring.mul_assoc, hs c]
    rw [Ring.sub_eq_add_neg, Semiring.right_distrib, Semiring.right_distrib, Semiring.left_distrib,
      Semiring.left_distrib, Ring.neg_mul, Ring.neg_mul, Ring.neg_mul, a3, a4, a5, a6]
    grind
  have jo : J * O * J = -O := by
    rw [Semiring.mul_assoc, hOJ, ← Semiring.mul_assoc, hJ, Ring.neg_mul, Semiring.one_mul]
  have z : O * J - J * O = 0 := by rw [hOJ]; grind
  rw [e1, z, jo, Semiring.mul_zero, Ring.mul_neg]
  have f : c * c * O + s * s * O = O := by
    rw [← Semiring.right_distrib, hcs, Semiring.one_mul]
  grind

end CircleOp.Unseen
