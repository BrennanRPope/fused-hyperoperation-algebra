# FusedAlgebra.Topics.E

Theorems (37): SAME 5, GENERAL 16, MATHLIB 16. Definitions checked:
- `fall` = `Nat.descFactorial n k` (same argument order) and `tri k` = `k.choose 2`: both correct.
- `rsum` = `Finset.sum` over `range m`, with no off-by-one in `s` or `stair_eq_cube`.
- `layer`/`fact` = `choose`/`factorial`.
- `stair`, `wt`, `s`, `u` are the same over any field.

There is no theorem verdict of WEAKER, DIFFERENT, MISSING or VACUOUS, and none at low confidence.

Numeric checks:
- `wt 100 3 = 4851/30000`: the same numbers, over ℚ.
- `s 10` / `u 10` against 27182818/10^7 and 27182819/10^7: the same numbers, now the named theorem `s_ten_u_ten` in any ordered field.
- `stair 100` between 2704/1000 and 2705/1000: the same numbers.
All three went from `decide +kernel` to `norm_num`.

Route changes:
- `binomial`: his Pascal induction over `rsum` is replaced by mathlib's `add_pow`.
- `wt_le`: it no longer goes through his `layer_fact_le`/`layer_mul_fact`; it uses `Nat.choose_le_pow_div`. His `layer_fact_le` (now `choose_mul_factorial_le_pow`) is itself proved from mathlib's descFactorial lemmas.
- `s_bound`: uses mathlib's `factorial_mul_pow_le_factorial` in place of his `pow2_le_fact`.
- Kept as his own proofs: `collisions`, `bernoulli`, `bernoulli_strict`, `stair_nat_lt`, `wt_ge`, `stair_lt_succ`, `s_le_stair`.

New port theorems with no original: 6 (`s_ten_u_ten`, which promotes an example, plus `tendsto_stair`, `tendsto_s`, `s_le_exp_one`, `exp_one_le_u`, `exp_one_bounds`).
