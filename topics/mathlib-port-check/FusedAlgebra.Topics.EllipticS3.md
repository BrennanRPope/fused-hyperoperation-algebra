# FusedAlgebra.Topics.EllipticS3 — review

Theorems (17): EQUIV 7 (lift_11/12/21/22, cayley_hamilton, circle_moments, sphere_moments), GENERAL 7, SAME 1 (middle_minus_neighbour), MATHLIB 1 (sq_nonneg_rat → mul_self_nonneg), DROPPED-OK 1 (multAdd_is_1x1 → MultAdd.one_add_multAdd). Definitions walks, dyck SAME.

No WEAKER / DIFFERENT / MISSING / VACUOUS items, no low-confidence items.

Checks made: each of the four entry identities lift_ij is the (i,j) entry of the port's matrix identity (1 + a•G)(1 + b•G) = (1 - qab)•1 + wAdd k a b•G with the same G = !![k,1;-q,0]; quatMatrix is entry-for-entry the original's U = [[w+ix, y+iz],[-y+iz, w-ix]]; the moment values (2,6,20,70 and 1,2,5,14) and binomials (20,15,70,56) are the same numbers; quaternion coordinates moved from ℤ to ℝ (ℤ case follows by casting).

Route changes:
- sq_nonneg_rat: the original's hand proof (cases on the sign of x) is replaced by mathlib's mul_self_nonneg.
- turning_no_wall: the original's explicit identity 4(l² - kl + q) = (2l - k)² + (4q - k²) is replaced by nlinarith with the hint sq_nonneg (2l - k); same idea, not written out.

New port theorems with no original: 3 (trace_gen, det_gen, trace_quatMatrix); new definitions gen, quatMatrix.
