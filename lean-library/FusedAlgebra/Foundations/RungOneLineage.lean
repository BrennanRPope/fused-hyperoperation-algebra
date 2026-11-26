/-
  Foundations.RungOneLineage: where each older name for the rung-1 semiring came from
  ==================================================================================

  Core Lean 4 only (no mathlib, no `sorry`). Nov 25, 2026.

  The rung-1 semiring (shifted add a + b + 1 as "or", the circle operation
  a + b + ab as "and") is proved once, on its own, in
  FusedAlgebra.Foundations.RungOneSemiring. Before that file existed, the project met
  the same pair five times, in five topics, and gave it a local name each time. All
  five are kept, unchanged, as a record of how the idea grew. This file shows that
  each one IS the unified definition (by `rfl`) and that each older distributive law
  follows from the unified law `RungOne.circ_sadd_distrib`.

  Names are given as they are in the library now; the name before the 1.0 rename
  (Nov 25, 2026) follows in brackets where it changed. NAMING-COLLISIONS.md has the full table.

  | Older name                                | Introduced in (library copy)                                 | Standard name                                  |
  |-------------------------------------------|--------------------------------------------------------------|------------------------------------------------|
  | `CircleOp.sadd` (on Nat) [`MultAdd.sadd`] | ns-rungs/lean/FusedRungs.lean → Rungs (rung 1 of the ladder)  | the hyperoperation + moved by x ↦ x − 1        |
  | `CircFound.FusedPeano.sadd`, `circ_sadd` [`MAFound…`, `madd_sadd`] | axioms/lean/Axioms/FusedPeano.lean → Foundations.FusedPeano | (ℕ⁺, +, ×) by transport; Peano with the shift |
  | `AbstractAlgebra.fadd`, `circ_fadd_distrib` [`multAdd_fadd_distrib`] | abstract-algebra/lean → Topics.AbstractAlgebra ("fused rung pair") | transported ring; Mathlib `Equiv.commRing` |
  | `PvsNP.sAdd`, `circ_sAdd_distrib` [`multAdd_sAdd_distrib`] | p-vs-np (Nov 25, 2026) → Topics.PvsNP ("shifted ruler") | circuits over an isomorphic semiring |
  | `LinearAlgebraCirc.sadd`, `fl_sadd` [`LinearAlgebraMA…`] | linear-algebra/lean → Topics.LinearAlgebra ("fused addition") | linear maps about the corner (−1, −1) |

  Elsewhere the pair is known as: the ring (R, a+b+1, a+b+ab) obtained by transport of
  structure (a textbook exercise); Jacobson's circle composition for the "and"
  (Mathlib `PreQuasiregular`); for a dial k, a generated semiring with generator
  g(x) = 1 + kx in E. Pap's pseudo-analysis; and, read on the height ruler a = eˢ − 1,
  the log semiring (log-sum-exp, +) whose cold limit is the tropical semiring.
  Sources and pictures: math-theory/semiring-audit/README.md.
-/
import FusedAlgebra.Foundations.RungOneSemiring
import FusedAlgebra.Rungs
import FusedAlgebra.Foundations.FusedPeano
import FusedAlgebra.Topics.AbstractAlgebra
import FusedAlgebra.Topics.PvsNP
import FusedAlgebra.Topics.LinearAlgebra

namespace RungOneLineage

open Lean.Grind

set_option linter.unusedSectionVars false

/-! ## 0. The circle operation itself is the unified "and" -/

/-- `CircleOp.circ` (Base, the project's first operation, called `CircleOp.circ` and "mult-add"
    before 1.0) is `RungOne.circ`. -/
theorem multAdd_eq_circ {R : Type} [CommSemiring R] (a b : R) :
    CircleOp.circ a b = RungOne.circ a b := rfl

/-! ## 1. Rungs: rung 1 on Nat -/

/-- `CircleOp.sadd` (Rungs, `MultAdd.sadd` before 1.0, from ns-rungs/lean/FusedRungs.lean): rung 1 of the fused
    hyperoperation ladder on Nat. It is the unified shifted add. -/
theorem rungs_sadd_eq (a b : Nat) : CircleOp.sadd a b = RungOne.sadd a b := rfl

/-! ## 2. FusedPeano: rung 1 and rung 2 as Peano primitives -/

section Peano
open CircFound CircFound.FusedPeano
variable {N : Type} [FusedPeano N]

/-- `FusedPeano.circ_sadd` (`madd_sadd` before 1.0; axioms/lean/Axioms/FusedPeano.lean), "rung 2 distributes
    over rung 1", proved there by induction from the fused Peano axioms. Here it is
    rederived from the unified law: every model is a copy of Nat (`toN_surjective`),
    counting carries the model's sadd and circ to `RungOne.sadd` and `RungOne.circ`
    (`toN_sadd`, `toN_circ`), and on Nat the unified law holds. -/
theorem madd_sadd_from_unified (a b c : N) :
    circ a (sadd b c) = sadd (circ a b) (circ a c) := by
  obtain ⟨x, rfl⟩ := toN_surjective (N := N) a
  obtain ⟨y, rfl⟩ := toN_surjective (N := N) b
  obtain ⟨z, rfl⟩ := toN_surjective (N := N) c
  have h := RungOne.circ_sadd_distrib x y z
  unfold RungOne.circ RungOne.sadd at h
  rw [← toN_sadd, ← toN_circ, ← toN_circ, ← toN_circ, ← toN_sadd, h]

end Peano

/-! ## 3. AbstractAlgebra: "the fused rung pair is a ring" -/

/-- `AbstractAlgebra.fadd` (abstract-algebra/lean, section 5, "rung 1 on a ring") is
    the unified shifted add. -/
theorem fadd_eq {R : Type} [CommRing R] (a b : R) :
    AbstractAlgebra.fadd a b = RungOne.sadd a b := rfl

/-- `AbstractAlgebra.circ_fadd_distrib` (`multAdd_fadd_distrib` before 1.0) is an instance of the unified law. -/
theorem multAdd_fadd_distrib_from_unified {R : Type} [CommRing R] (a b c : R) :
    CircleOp.circ a (AbstractAlgebra.fadd b c)
      = AbstractAlgebra.fadd (CircleOp.circ a b) (CircleOp.circ a c) :=
  RungOne.circ_sadd_distrib a b c

/-- `AbstractAlgebra.fadd_group` and `fused_zero_absorbs` are parts of the unified
    `RungOne.isCommRing` (additive group with zero −1; −1 absorbs). -/
theorem fused_ring_from_unified {R : Type} [CommRing R] :
    RungOne.CommRingLaws (AbstractAlgebra.fadd (R := R)) CircleOp.circ (-1) 0 :=
  RungOne.isCommRing

/-! ## 4. PvsNP: the shifted ruler for circuits -/

/-- `PvsNP.sAdd` (p-vs-np, Nov 25, 2026, "shifted add on the shifted ruler") is the
    unified shifted add. -/
theorem sAdd_eq {R : Type} [CommSemiring R] (a b : R) :
    PvsNP.sAdd a b = RungOne.sadd a b := rfl

/-- `PvsNP.circ_sAdd_distrib` (`multAdd_sAdd_distrib` before 1.0) is an instance of the unified law. -/
theorem multAdd_sAdd_distrib_from_unified {R : Type} [CommSemiring R] (a b c : R) :
    CircleOp.circ a (PvsNP.sAdd b c) = PvsNP.sAdd (CircleOp.circ a b) (CircleOp.circ a c) :=
  RungOne.circ_sadd_distrib a b c

/-! ## 5. LinearAlgebra: fused-linear maps -/

/-- `LinearAlgebraCirc.sadd` (`LinearAlgebraCirc.sadd` before 1.0; linear-algebra/lean, section 7, "the fused addition, zero −1")
    is the unified shifted add. -/
theorem linearAlgebra_sadd_eq {R : Type} [CommRing R] (a b : R) :
    LinearAlgebraCirc.sadd a b = RungOne.sadd a b := rfl

/-- `LinearAlgebraCirc.fl_sadd` says fused-linear maps respect the unified shifted add:
    they are the linear maps of the rung-1 ring (semimodule maps over it). Restated with
    the unified name. -/
theorem fl_sadd_unified {R : Type} [CommRing R] (m : LinearAlgebraCirc.M2 R) (u1 u2 v1 v2 : R) :
    LinearAlgebraCirc.fl1 m (RungOne.sadd u1 v1) (RungOne.sadd u2 v2)
      = RungOne.sadd (LinearAlgebraCirc.fl1 m u1 u2) (LinearAlgebraCirc.fl1 m v1 v2) ∧
    LinearAlgebraCirc.fl2 m (RungOne.sadd u1 v1) (RungOne.sadd u2 v2)
      = RungOne.sadd (LinearAlgebraCirc.fl2 m u1 u2) (LinearAlgebraCirc.fl2 m v1 v2) :=
  LinearAlgebraCirc.fl_sadd m u1 u2 v1 v2

end RungOneLineage
