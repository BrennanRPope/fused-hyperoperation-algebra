# FusedAlgebra.Topics.E: port remarks

- `stair`, `wt`, `s`, `u` take the field as an explicit first argument (`stair K n`), so the
  bounds hold in any linearly ordered field and can be used over both `ℚ` and `ℝ`.
- The three kernel-checked `example`s are kept. The `s 10` / `u 10` check became the named theorem
  `s_ten_u_ten` (any ordered field), because the new `exp_one_bounds` uses it over `ℝ`. They are
  proved by `norm_num` instead of `decide +kernel`.
- New, not in the original (which had no real numbers): `tendsto_stair`, `tendsto_s`,
  `s_le_exp_one`, `exp_one_le_u`, `exp_one_bounds`. These are short links to mathlib
  (`Real.tendsto_one_add_div_pow_exp`, `NormedSpace.expSeries_div_hasSum_exp`) and supply the
  step the original header listed as "what Lean does NOT do here". Mathlib already has sharper
  digits (`Real.exp_one_gt_d9`, `Real.exp_one_lt_d9`).
- `bernoulli`, `bernoulli_strict`, `stair_nat_lt` are kept in `ℕ` with the original statements;
  they are the elementary route to `stair_lt_succ`.
