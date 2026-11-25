# FusedAlgebra.Topics.LeeYangSpill

Theorems (10): SAME 5, GENERAL 2, EQUIV 2, MATHLIB 1. WEAKER/DIFFERENT/MISSING/VACUOUS: none. Low confidence: none.

- `slowest_step` and `turning_never_stalls` are generalised from ℤ to linearly ordered comm rings, with the same hypotheses.
- `tan_square` and `tan_fourth` were four bare ℤ arithmetic facts each (the entry sums). The port states the matrix identities `!![1,1;-1,1]^2 = !![0,2;-2,0]` and `^4 = (-4)•1`. The entries match (0, 2, -2, 0 and -4, 0, 0, -4), so these are EQUIV and, if anything, more meaningful.

Route changes:
- `sq_nonneg'` → `mul_self_nonneg`.
- The inequalities are now proved via `le_add_of_nonneg_left` / `positivity` instead of omega.

New port theorems with no original: 0.
