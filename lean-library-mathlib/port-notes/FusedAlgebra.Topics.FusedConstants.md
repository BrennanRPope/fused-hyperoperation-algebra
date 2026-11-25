# FusedAlgebra.Topics.FusedConstants

- `H r` is kept as a name but is now `hyperoperation (r + 1)` by definition, so the ladder is
  mathlib's. New lemma `rung_succ_eq` ties it to `MultAdd.rung`.
- `fused_one_one` is restated as `rung_succ_one_one : rung (r + 1) 1 1 = 3`.
- The `slide` lemmas are stated in any ring instead of `ℤ`.
