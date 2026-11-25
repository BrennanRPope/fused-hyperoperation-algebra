# FusedAlgebra.Topics.DataCenter

Theorems: 16. GENERAL 15, MATHLIB 1. No findings.

- `lifts` is now `(rs.map (· + 1)).prod`: the same recursion. `lifts_append` → `List.prod_append` (with `map_append`). Correct.
- `scale_anywhere`: `lifts xs` and `lifts ys` are replaced by arbitrary x and y. Strictly more general, and the original is an instance.
- `add_unit_iff` and `spare_fan_helps` now hold in ordered (semi)rings with the same hypotheses.

Route changes:
- `foldMA_swap` (now `multAddList_swap`) is proved from `multAddList_perm` instead of his argument through `lifts` and omega.
- `add_unit_iff` is proved by `sub_pos` and `unit_gap` instead of omega.
- His `lifts_append` induction is replaced by mathlib's `List.prod_append`.

New port theorems with no original: 0.
