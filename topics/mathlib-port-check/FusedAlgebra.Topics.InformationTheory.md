# FusedAlgebra.Topics.InformationTheory

Theorems (52): GENERAL 48, MATHLIB 2, DROPPED-OK 1, SAME 1. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence verdicts.

## Checks
- `outer` is now `(xs ×ˢ ys).map (p.1 * p.2)`. List product runs row by row, so the list and its order are the same, and `outer_cons` recovers the recursion.
- `shan` is now `(xs.map (x * L x)).sum`. `shan_nil`/`shan_cons` give the old equations.
- `lsum` is `List.sum` and `pw` is `^`.
- `climb` is dropped and covered exactly by `MultAdd.one_add_mul_wpow`.
- Numbers unchanged: `ladder_four` uses `-(1+1/3)` and `-4/3`; `piling_up` uses `1/2`. Both are in a characteristic-zero field, which is automatic over ℚ.
- `bec_dpi` is stated in mathlib's ordered-monoid generality, which ℚ satisfies.
- Every table claim checked out. The duplicate-name groups (`kraft_slack` = `gini_simpson`, and so on) are each re-proved from the first member, with the same statements.

## Route changes
- `union_rule`, `union_is_flipped_madd`, `series_system`, `shannon_weight`, `lengths_add`, `check_node`, `uniform_tsallis_is_fused_power`, `cube_corners` are now one-line uses of Basic lemmas, rather than the original's `unfold; grind` or `climb` arguments.
- `hartley` and `parallel_capacity` now use mathlib's `pow_add` instead of his `pw_add`. `bec_dpi` is mathlib's `mul_le_of_le_one_right`.

New port theorems with no original: 4 (`outer_nil`, `outer_cons`, `shan_nil`, `shan_cons`).
