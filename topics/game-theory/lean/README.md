# GameTheory.lean

- 36 theorems, core Lean 4.37.0-rc1 (no mathlib), no `sorry`, axioms only propext, Classical.choice, Quot.sound.
- Compiles clean against the shared library build: `LEAN_PATH=/tmp/lib/build lean lean/GameTheory.lean`.
- Imports `FusedAlgebra.Base` only. Namespace `GameTheory`, the library's scoped notation `◯` for `circ` (`open scoped CircleOp`).
  Defines its own dial fold `dfold k l` (with `dfold_lift`, `dfold_append`, `dfold_one`).
  Section 5 works over `Rat` (`abbrev K := Rat`) because the 1/2, 1/3, 1/6 weights need characteristic 0.
- Intended module path: `FusedAlgebra/Topics/GameTheory.lean`.

Main results: `kelly_many`, `kelly_ruin`, `kelly_win_lift`, `kelly_loss_lift`, `kelly_ratio`, `kelly_foc`;
`hedge_step`, `hedge_expert`, `run_total`, `run_expert`, `lifted_regret` (two experts: best expert's
growth ≤ 2 × the learner's); `odds_circ`, `edge_compose` (velocity addition), `odds_is_lift`,
`prob_compose`; `dial_game_union` (Sugeno's rule), `marginal`, `corners3`, `multilinear2/3`,
`shapley2_eq`, `shapley2_union`, `shapley3_eq`, `shapley3_efficient`, `shapley3_corners`,
`banzhaf3_eq`, `bonus_gap`, `bigger_smaller_bonus`.
