# Lean: InformationFull.lean

- 31 theorems, core Lean 4.37.0-rc1, no mathlib, no `sorry`.
- Compiles clean: `LEAN_PATH=/tmp/lib/build lean lean/InformationFull.lean` (no errors, no warnings).
- Axioms checked with `#print axioms` on qlogN_mul, qlogNeg_mul, qexpN_mul, tsallis_engine, graft_circ, full_inner, kraft_one, channel_cascade, dary_graft: only propext, Classical.choice, Quot.sound.
- Imports `FusedAlgebra.Base` and `FusedAlgebra.Topics.InformationTheory` (reuses `lsum`, `outer`, `shan`, `outer_sum`).
- Intended module path: `FusedAlgebra/Topics/InformationFull.lean`.

Parts: (1) q-logarithm at whole-number weights turns products into the dial (`qlogN_mul`), its lift, inverse (`qlog_of_root`), k = −1 case, duality q ↔ 2 − q as weight flip (`qlog_dual`, `qlogNeg_mul`), q-exponential at k = 1/n (`qexpN_mul`, `qexp_climb` = compound interest); (2) weighted outer-log lemma `outer_shan_w` and `tsallis_engine` (Tsallis' law as expectation through the dial; Shannon at k = 0); (3) `sm_compose` (Sharma–Mittal), `flip_weight`; (4) binary trees: `leaves_lift`, `graft_leaves`, `graft_inner`, `graft_circ`, `full_inner` (= circle power of 1), `kraft_one`; D-ary counts `dary_graft`, `dary_leaves`, `dary_full`; (5) any binary channel: `channel_cascade` (s = a + b chains by k = −1), `channel_det`, `bsc_from_union`, `z_cascade`.
