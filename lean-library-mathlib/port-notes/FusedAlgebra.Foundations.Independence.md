# FusedAlgebra.Foundations.Independence: porting remarks

- **Shape of the four theorems changed.** The original `no_M1` … `no_M4` were eight-part
  conjunctions `A1 ∧ A2 ∧ A3 ∧ A4 ∧ M1 ∧ …` over a hand-rolled `Ops` bundle. In the port the
  group part is mathlib's `AddCommGroup` instance on the model (`ℤ`, `Matrix (Fin 2) (Fin 2) ℤ`,
  `ℤ × ℤ × ℤ`), so each theorem is the four-part conjunction about M1 to M4 only. Nothing
  mathematical is lost: A1 to A4 hold by the instance.
- The module now imports `Foundations.MultAddRing` (the original was stand-alone) so that the
  predicates can be tied to the class: `axioms_of_maRing` and `maRingOfAxioms` are new.
- New general lemmas: `m2_multAdd`, `m3_multAdd`, `m4_multAdd` (the circle operation of any ring,
  commutative or not, satisfies M2, M3, M4) and `wAdd_axioms` (the weighted mult-add of any
  commutative ring satisfies M1, M2, M4). These could live in `FusedAlgebra.Basic`.
- The models are the original ones (over `ℤ`, so concrete but not finite).
