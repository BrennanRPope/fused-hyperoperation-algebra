# FusedAlgebra.Foundations.MultAddRing: porting remarks

- **`MARing` now extends `AddCommGroup`.** The original class carried raw `zero`/`add`/`neg`
  fields with axioms A1 to A4; those are exactly mathlib's `AddCommGroup`, so the port uses it as
  the parent. All the group rearrangement helpers became mathlib lemmas. One side effect: `nsmul`
  and `zsmul` are data of the parent, not `natCast n * a` as in the original (they agree, by
  `nsmul_eq_mul` / `zsmul_eq_mul`).
- **`MARing.instCommRing` and `MAField.instField` are global instances** to mathlib's `CommRing`
  and `Field`, as in the original (which targeted `Lean.Grind.CommRing`/`Field`). `ofCommRing` and
  `ofField` are `abbrev`s, not instances, so there is no instance loop.
- **New, beyond the original** (round trips as structure equalities):
  `MARing.ext'`, `instCommRing_ofCommRing`, `ofCommRing_instCommRing`,
  `equivCommRing : MARing R ≃ CommRing R`, `ofField`, `toCommRing_instField_ofField`.
  For fields only the ring structure is shown to round-trip: mathlib has no extensionality lemma
  for `Field` (it carries extra data such as `ratCast`, `qsmul`, `zpow`), and the inverse is
  determined by the ring anyway. No `MAField R ≃ Field R` is stated for the same reason.
- `MARing.mul`/`MARing.one` are no longer separate definitions; they are the `Mul`/`One`
  instances (`mul_def`, `one_eq_neg_wall`). The ring-axiom lemmas (`MARing.mul_comm`, …) are kept
  as `protected` theorems because they are the content of "M1 to M4 give a ring".
