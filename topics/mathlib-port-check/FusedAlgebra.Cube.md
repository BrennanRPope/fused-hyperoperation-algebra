# FusedAlgebra.Cube

Theorems (7): MATHLIB 7. No WEAKER / DIFFERENT / MISSING / VACUOUS.

layer is Nat.choose (same three-case recursion with Pascal's rule) and fact is Nat.factorial (same recursion). Each lemma matches its mathlib name up to the sides of the equation and mul_comm. layer_step's natural subtraction (n - k) is the same as in Nat.choose_succ_right_eq.

Medium confidence: `absorb` → `Nat.add_one_mul_choose_eq`. The statement fits, but this is the newer name of Nat.succ_mul_choose_eq, and I can't confirm it exists under this name in the pinned mathlib (`#check Nat.add_one_mul_choose_eq` settles it). The mathlib rows of the name tables are not checked by tools/CheckNames.lean, which only checks rows marked "ported".

Route: the whole module's own proofs are gone, notably absorb and the 30-line layer_step, replaced by mathlib lemmas.

New port theorems with no original: 4 (one_add_mpow_eq_sum_choose, mpow_eq_sum_choose, card_corners, mpow_one_add_one).
