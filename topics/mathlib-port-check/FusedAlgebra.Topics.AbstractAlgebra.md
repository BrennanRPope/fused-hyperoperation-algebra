# FusedAlgebra.Topics.AbstractAlgebra

Theorems (62): SAME 30, EQUIV 10, MATHLIB 13, GENERAL 5, DROPPED-OK 4. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

Definitions checked (12): the records CommMonoidLaws/CommGroupLaws/OpIso → CommMonoid/CommGroup/≃*; tateSetoid/TateQ → Con on PreQuasiregular R and its quotient (same relation through `equiv`); geom → `∑ i ∈ range m, x ^ i` (same value); Fin n → ZMod n (definitionally Fin n for n ≥ 1, same ring operations, `ZMod.val` = `Fin.val`). All the numbers in the decide checks (6/2/3 criteria; units 0,4,6,10; wall 11; roots 2,3 mod 5 and 5,8 mod 13; wall 6 mod 7) are unchanged.

Medium-confidence item:
- `fadd_group` (`CommGroupLaws fadd (-1)`) → `Fused.instCommRing`. Nothing in the port states the group laws on `fadd` directly. They follow from the additive group of `Fused R` (built by `Equiv.commRing`) through `Fused.val_add` (`(x+y).val = fadd x.val y.val`), `val_zero` (`= -1`) and `val_neg` (`= -x.val - 2`). The content is all there, but you need a small unwrapping step to get the original form back.

Notes on table claims (all correct):
- `wAdd_unit_iff` was an implication, as the table says, and it is ported as the same implication.
- `int_scale_hom` maps to a def (`Weighted.scaleHom`). The identity itself is `MultAdd.mul_wAdd` (any comm semiring).
- `shift_ring_hom`'s last two conjuncts (`1 + -1 = 0`, `1 + 0 = 1`) are trivial ring facts, so the RingEquiv covers the conjunction.

Route changes:
- `geom_mul`: the original's own induction is replaced by mathlib `mul_neg_geom_sum`, and `nilpotent_weight_inverse` now uses it.
- `int_multAdd_units`: the natAbs argument is replaced by `Int.isUnit_iff` + `lia`.
- `mpow_multAdd`, `zero_mpow`: induction is replaced by `mul_pow`/`one_pow` in PreQuasiregular.
- `mpow_mul`, `mpow_one`: the AbstractAlgebra proofs are dropped in favour of the Basic versions.
- `tateQ_commMonoid`: replaced by mathlib `Con.commMonoid`.

New port theorems with no original: 13.
