/-
  Fused Hyperoperation Algebra: the shift rule (A4) from composite systems
  (topics/why-i/rotations-force-i, paper Proposition "shift from composites").
  Core Lean 4 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/Booking.lean.

  Two systems A, B. A constant c added to the total energy of non-interacting systems is the same
  joint Hamiltonian whether it is booked to A or to B:
      (H_A + c) ⊗ 1 + 1 ⊗ H_B  =  H_A ⊗ 1 + 1 ⊗ (H_B + c).
  If the joint generator depends only on the joint Hamiltonian, and non-interacting parts evolve
  independently (generators add: A_AB = A_A ⊗ 1 + 1 ⊗ A_B), and local operators act faithfully on
  the joint space (X ⊗ 1 = Y ⊗ 1 ⇒ X = Y), then the change D(c) := A(H + c) − A(H) does not
  depend on H, and it is additive in c. (Linearity, D(c) = c·A(1), then follows from boundedness:
  (A1) and (A3) give ‖A(c)‖ = |c|/ħ; that analytic step is in the paper.)
  The operator spaces are abstract rings (e.g. real matrices); `L` is X ↦ X ⊗ 1 and `R` is Y ↦ 1 ⊗ Y.
-/
namespace Booking
open Lean.Grind

variable {S HA HB HW OA OB OW : Type} [Ring OA] [Ring OB] [Ring OW]

theorem sub_eq_iff (a b c d : OW) (h : a + b = c + d) : a - c = d - b := by grind

/-- The two bookings give the same joint generator, so the shifts of A and of B agree on W. -/
theorem booking (shA : HA → S → HA) (shB : HB → S → HB) (joint : HA → HB → HW)
    (AA : HA → OA) (AB : HB → OB) (AW : HW → OW) (L : OA → OW) (R : OB → OW)
    (hL : ∀ x y, L (x - y) = L x - L y) (hR : ∀ x y, R (x - y) = R x - R y)
    (book : ∀ hA hB c, joint (shA hA c) hB = joint hA (shB hB c))
    (indep : ∀ hA hB, AW (joint hA hB) = L (AA hA) + R (AB hB)) (hA : HA) (hB : HB) (c : S) :
    L (AA (shA hA c) - AA hA) = R (AB (shB hB c) - AB hB) := by
  have e := congrArg AW (book hA hB c)
  rw [indep, indep] at e
  rw [hL, hR]
  exact sub_eq_iff _ _ _ _ e

/-- The shift of A's generator by c is the same for every Hamiltonian of A. -/
theorem shift_independent (shA : HA → S → HA) (shB : HB → S → HB) (joint : HA → HB → HW)
    (AA : HA → OA) (AB : HB → OB) (AW : HW → OW) (L : OA → OW) (R : OB → OW)
    (hL : ∀ x y, L (x - y) = L x - L y) (hR : ∀ x y, R (x - y) = R x - R y)
    (faithful : ∀ x y, L x = L y → x = y)
    (book : ∀ hA hB c, joint (shA hA c) hB = joint hA (shB hB c))
    (indep : ∀ hA hB, AW (joint hA hB) = L (AA hA) + R (AB hB)) (hB : HB)
    (hA hA' : HA) (c : S) :
    AA (shA hA c) - AA hA = AA (shA hA' c) - AA hA' := by
  apply faithful
  rw [booking shA shB joint AA AB AW L R hL hR book indep hA hB c,
      booking shA shB joint AA AB AW L R hL hR book indep hA' hB c]

/-- Telescoping: an H-independent shift is additive, D(c + c') = D(c) + D(c'). -/
theorem shift_additive [Add S] (shA : HA → S → HA) (AA : HA → OA)
    (shsh : ∀ h c c', shA (shA h c) c' = shA h (c + c'))
    (indepH : ∀ h h' c, AA (shA h c) - AA h = AA (shA h' c) - AA h') (h : HA) (c c' : S) :
    AA (shA h (c + c')) - AA h = (AA (shA h c) - AA h) + (AA (shA h c') - AA h) := by
  have e := indepH (shA h c) h c'
  rw [shsh] at e
  have : AA (shA h (c + c')) - AA h
      = (AA (shA h (c + c')) - AA (shA h c)) + (AA (shA h c) - AA h) := by grind
  rw [this, e]; grind

end Booking
