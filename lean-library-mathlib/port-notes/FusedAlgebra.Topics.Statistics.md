# FusedAlgebra.Topics.Statistics: port remarks

- Indicators are mathlib's `IsIdempotentElem`; `IsInd`, `ind_cases`, `and_ind`, `not_ind` are
  mathlib declarations. Parts 1 to 4 are stated in any commutative ring, the parts that divide in
  any field (`clt_skew`, `clt_kurt`, `fair_coin_kurt` need characteristic 0; `bern_var_parity`
  and `pool_velocity` need `2 ≠ 0`, and `pool_velocity` is false in characteristic 2).
- The five moment definitions `mass`, `m1`, …, `m4` are one definition `moment r`. The five
  `_append` and five `_shift` lemmas are the single general `moment_append` and `moment_shift`.
  New general theorem `moment_conv`: `E[(X + Y) ^ r] = ∑ j, r.choose j * E[X ^ (r - j)] * E[Y ^ j]`
  for every `r`; the original rows `1 1`, `1 2 1`, `1 3 3 1`, `1 4 6 4 1` are kept as named
  corollaries (`moment_zero_conv`, …, `moment_four_conv`). `scale_m` became `moment_scale` for
  every `r`.
- Conjunction theorems `scale_k` and `bern_stats` are split into one lemma per fact.
- `bp_layers` no longer needs `j ≤ n`; `bp_zero_col` and `bp_above` follow from it.
- `bayes_odds` lost two hypotheses the proof never used (`p ≠ 0`, `lH ≠ 0`).
- Distributions stay finite lists of weighted outcomes, as in the original. They were not moved
  to mathlib's `PMF` / measure theory: the theorems here do not assume nonnegative or normalised
  weights (they are ring identities), and `PMF` would force `ℝ≥0∞` weights. The module docstring
  points to the measure-theoretic counterparts in mathlib.
- New small lemmas: `unionL_nil`, `unionL_cons`, `expect_nil`, `expect_cons`, `gen_nil`,
  `gen_cons`, `nconv_nil`, `nconv_cons`, `conv_nil`, `conv_cons`, `moment_nil`, `moment_cons`,
  `foldW_neg_one`.
