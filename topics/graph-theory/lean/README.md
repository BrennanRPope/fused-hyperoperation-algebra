# GraphTheory.lean

- 38 theorems (3 of them small division helpers: add_div', ddiv, cancel_left), core Lean 4.37.0-rc1,
  no mathlib, no `sorry`.
- Compiles clean (no errors, no warnings) with
  `LEAN_PATH=/tmp/lib/build lean lean/GraphTheory.lean`.
- `#print axioms` on parallel_tr, potts_is_amh, ising_fixed, c4_indep, p3_strong_p3: only propext,
  Classical.choice, Quot.sound (some none).
- Imports FusedAlgebra.Base and FusedAlgebra.Topics.Statistics (for `StatFused.lens`, `StatFused.amh`
  and `amh_union_is_lens`). Intended module path: FusedAlgebra/Topics/GraphTheory.lean
  (must build after Topics/Statistics).

Main theorems: parallel_join, series_local, series_tr (series multiplies t), parallel_tr (parallel is
the landscape law F_{Q−2,1−Q}), or_line, parallel_lift, q_one (reliability), q_two (velocity addition),
q_zero_cusp, potts_is_amh (θ = 1 − Q), perfect_bond_absorbs, colouring_bond, colouring_wall_absorbs,
series_fails, series_odds, wheatstone, diamond_poly, golden_fixed, golden_slope, diamondQ_one,
ising_fixed, three_state_fixed, strong_entry, strong_count, p3_strong_p3, c4_strong_c4, c5_strong_p2,
union_circ, join_adds, c4_indep, two_k2_indep, c5_indep.
