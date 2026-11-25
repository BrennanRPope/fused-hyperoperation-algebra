# FusedAlgebra.Base (port: FusedAlgebra/Basic.lean)

Theorems (46): SAME 20, GENERAL 15, MATHLIB 9, DROPPED-OK 2. No WEAKER / DIFFERENT / MISSING / VACUOUS, no low-confidence items.

Definitions checked:
- multAdd: same formula a + b + a*b, now for any [Add R] [Mul R]. Notation ✳ is infixl 65, as ⊞ was.
- mpow: now `(equiv a ^ n).val` in mathlib's PreQuasiregular. PreQuasiregular's product is y + x + x*y, so mpow 0 = 0 and mpow (n+1) = mpow n ✳ a, which is the original recursion exactly (mpow_zero, mpow_succ).
- foldMA → multAddList: List.prod in PreQuasiregular, with nil = 0 and cons = a ✳ rest, the same recursion.
- groups → List.sublists': the same recursion (sublists'_cons).
- wAdd, wpow, multAddInv, msub: unchanged. pw, prodL and sumL are mathlib's ^, List.prod and List.sum, with the same recursions.

GENERAL items only drop typeclass assumptions (CommSemiring becomes Semiring, NonUnitalSemiring, NonUnitalNonAssocSemiring or Add/Mul; CommRing becomes Ring). Nothing was added.

Medium confidence: `sumL_map_cons` → `List.sum_map_mul_left`. The original form also needs List.map_map and List.prod_cons; the content is covered.

Route changes worth knowing:
- The pw power laws and pow_add', mul_pow', one_pow' are now mathlib lemmas. The hand inductions are gone (they were routine).
- mpow_lift, foldMA_lift, mpow_add, foldMA_append, foldMA_perm and multAdd_assoc are now proved through PreQuasiregular: map_pow and map_list_prod of the hom x ↦ 1 + x, pow_add, List.prod_append, Perm.prod_eq and mul_assoc. Before, they were proved by induction or by direct algebra. wpow_lift is still proved by induction.

New port theorems with no original: 20.
