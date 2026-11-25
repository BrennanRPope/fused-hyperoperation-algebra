# FusedAlgebra.Topics.Statistics

Theorems (98): GENERAL 88, MATHLIB 7, EQUIV 2, DROPPED-OK 1. Definitions were checked row by row:
- `IsInd` is `IsIdempotentElem`, the same `x * x = x`.
- `groups` is literally `List.sublists'`, with the same cons recursion.
- `Ex`, `gen` and `litCorners` became list or Finset sums with the same terms.
- `nconv`/`conv` via `flatMap` keep the same cons recursion.
- `mass`/`m1`…`m4` are `moment 0..4`, using `x^0 = 1`.
- `sgn` is `(-1)^length`.
- `pEven`/`pOdd` use `Even` for `% 2 = 0`.
- `k2`/`k3`/`k4`, `skewSq` and `kurt` are the same polynomials.

There is no theorem verdict of WEAKER, DIFFERENT, MISSING or VACUOUS.

The two rows the brief flagged, checked:
- `bayes_odds`. Original hypotheses: `p ≠ 0, p ≠ 1, lH ≠ 0, lN ≠ 0, p*lH + (1-p)*lN ≠ 0`. Port: `p ≠ 1, lN ≠ 0, hden`. The port statement is still true: with D ≠ 0, both sides equal `p*lH / ((1-p)*lN)` whatever p and lH are. So this is GENERAL (strengthened), and the table's claim is correct.
- `bp_layers`. The original needs `j ≤ n`; the port has `∀ n j`. It is still true, because above the diagonal `n.choose j = 0` and `bp = 0`. The port proves it directly by the same Pascal recursion, with `Nat.choose_succ_succ'` and `choose_eq_zero_of_lt`. `bp_zero_col` and `bp_above` keep their exact statements and are now corollaries. GENERAL, and the table is correct.

Other things checked:
- Indicators as `IsIdempotentElem`. Part 1 now holds for all idempotents in a commutative ring. `or_ind`, `xor_ind` and `sq_diff_is_xor` are true there; the port proves them algebraically rather than by 0/1 cases. `ind_cases` holds only with no zero divisors, which is why it is left to mathlib's `IsIdempotentElem.iff_eq_zero_or_one` and not restated.
- Moment specialisations. `moment_shift` at r = 0..4 and `moment_conv` at rows 1/121/1331/14641 give exactly the original coefficients, in the original order (`moment (r-j) d * moment j e`). The `_append`, `_shift` and `scale_m` rows are covered by the general-r lemmas.
- Added typeclass hypotheses (`CharZero K` for `clt_skew`/`clt_kurt`/`fair_coin_kurt`, `NeZero (2:K)` for `bern_var_parity`/`pool_velocity`) are automatic for ℚ, and `0 < n` became `n ≠ 0`, which is the same on ℕ.

Medium-confidence items (mathlib names only; the statements are fine):
- `ind_cases` → `IsIdempotentElem.iff_eq_zero_or_one`. Settle with `#check`.
- `pw_neg_one_parity` → `neg_one_pow_eq_ite`. Settle with `#check`. The port itself uses `Even.neg_one_pow`/`Odd.neg_one_pow`, which certainly exist.

Route changes:
- The moment lemmas, previously five hand inductions each for append, shift and conv, are now derived from one general `moment_shift` (via `add_pow`) and `moment_conv`.
- `bp_zero_col` and `bp_above` are now corollaries of `bp_layers` instead of their own inductions.
- `gen_eq_litCorners` uses `E.binomial`, which is now mathlib's `add_pow`.
- `pw_neg_one_parity`, `groups_map`, `prodL_neg_sgn` and `sgn_map` are replaced by mathlib list and parity lemmas.

New port theorems with no original: 14 (12 nil/cons simp lemmas, `foldW_neg_one`, `moment_conv`).
