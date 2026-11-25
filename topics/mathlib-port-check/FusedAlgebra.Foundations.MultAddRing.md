# FusedAlgebra.Foundations.MultAddRing

Theorems (51): SAME 28, MATHLIB 18, DROPPED-OK 3, EQUIV 2. No WEAKER / DIFFERENT / MISSING / VACUOUS, no low-confidence items.

- MARing now extends AddCommGroup. The extra data (sub, nsmul, zsmul) is forced by its laws, so original MARing structures and port MARing structures correspond one to one. A1–A4 and the group helper lemmas (add_neg, neg_unique, neg_add_distrib, …) are mathlib lemmas, each checked to state the original.
- M4 (madd_add): `(madd a b + madd a c) + -a` became `madd a b + madd a c - a`, which are equal.
- The derived product is `madd a b - a - b` and the unit is `-wall`, as before. madd_zero is still derived from M4 alone. mul_comm, mul_assoc, mul_add and the other ring laws are still proved from M1–M4 by the user's route, with abel replacing ac_rfl.
- Default natCast, intCast, npow and zpow have the same values as the original recursions.
- way_back_formula: `multAddInv u = -u/(1+u) = -u * inv (1+u)`, the same statement with the same hypothesis.
- instCommRing and instField are now mathlib CommRing and Field instead of Lean.Grind classes, so they are stronger. Port extras: equivCommRing : MARing R ≃ CommRing R, and ofField.

Route changes: inv_ne_zero, inv_inv' and inv_one' are now mathlib consequences of the Field instance (they were stepping stones to Lean.Grind's zpow laws). way_back_formula now goes through Basic's multAdd_multAddInv.

New port theorems with no original: 5 (mul_def, MARing.ext', instCommRing_ofCommRing, ofCommRing_instCommRing, toCommRing_instField_ofField).
