# FusedAlgebra.Rungs

Theorems (36): SAME 31, EQUIV 3, MATHLIB 2. No WEAKER / DIFFERENT / MISSING / VACUOUS, no low-confidence items.

## The rung definition (checked in full)
`rung n a b = hyperoperation n (a+1) (b+1) - 1`. Mathlib's hyperoperation is: n = 0 successor; n = 1 starts at m (addition); n = 2 starts at 0 (multiplication); n ≥ 3 starts at 1 (n = 3 power, n = 4 tetration, n = 5 pentation). So rung 1 = sadd, rung 2 = ✳, rung 3 = fexp, rung 4 = ftet and rung 5 = fpent, which matches the original numbering (rung_one … rung_five). rung 0 a b = b + 1, the step.

Natural subtraction: `hyperoperation_pos` proves hyperoperation n m k > 0 for all n when m, k > 0, and rung only ever evaluates it at a+1 and b+1. So the "- 1" never truncates, and rung n a b + 1 = hyp n (a+1) (b+1) for all n, a, b (rung_add_one). Zeros are safe: hyperoperation (n+4) 0 k can be 0, but rung never uses base 0.

`tet` = hyperoperation 4 and `pent` = hyperoperation 5 have exactly the original recursions (tet x 0 = 1, tet x (n+1) = x ^ tet x n; pent likewise) for every x including 0. So every original law about tet, ftet, pent and fpent means the same thing.

rung_succ_succ (rung (n+1) a (b+1) = rung n a (rung (n+1) a b)) is correct for every n, including n = 0, where it reads a + b + 2 = (a + b + 1) + 1.

All numeric checks have the same numbers, recomputed by hand: 63, 511, 7 ≠ 8, ftet values 3, 15, 26, 255, 65535, fpent values 3, 65535, and the four "laws fail" inequalities.

EQUIV items: tet_pos and pent_pos (1 ≤ x became 0 < x, the same proposition on ℕ), and fexp_base_iter (Nat.repeat became Function.iterate, equal for all n).

Route changes: multAdd_iter, fexp_iter, ftet_iter and fpent_iter are now instances of the single general rung_succ_succ, no longer separate lift computations. tet_one_base and pent_one_base use mathlib's hyperoperation_ge_three_one. ftet_values and its neighbours are checked by simp evaluation, no longer by decide. These are minor.

New port theorems with no original: 12 (hyperoperation_pos, rung_add_one, rung_succ_succ, rung_one…rung_five, tet_zero/succ, pent_zero/succ).
