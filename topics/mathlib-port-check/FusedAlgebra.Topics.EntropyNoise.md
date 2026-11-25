# FusedAlgebra.Topics.EntropyNoise

Theorems (41): GENERAL 28, MATHLIB 7, EQUIV 4, SAME 2. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence verdicts.

## PORT-NOTES claims, checked
- **Import of Finance: true.** The port imports `FusedAlgebra.Topics.Finance`, used only by `cancellation_hides` through `FinanceMultAdd.arith_mean_pow_le` (`0 ≤ x → 0 ≤ y → ((x+y)/2)^n ≤ (x^n+y^n)/2`). That is the port of his own Finance `power_mean`, still proved by hand induction.
- **`x ≤ 1` removed from `cancellation_hides`: correct strengthening.** The two-point power mean holds for every `x ≥ 0`.
  Original: `0 ≤ x → x ≤ 1 → ((1+x)/2)^n ≤ (1+x^n)/2`.
  Port: the same without `x ≤ 1`, over any linearly ordered field.
- **abs restatements (EQUIV): checked.**
  - `mul_bound` (`a≤u, -a≤u, b≤1, -1≤b ⊢ ab≤u ∧ -(ab)≤u`) ⇔ `|a|≤u → |b|≤1 → |ab|≤u`.
  - `octahedron_in_ball` with auxiliary bounds `u,v,w` (`±x≤u …, u+v+w≤1`) ⇔ `|x|+|y|+|z|≤1 → x²+y²+z²≤1`, in both directions (take `u = |x|`).
  - `compose_stays_in_octahedron` keeps all three conjuncts. `compose_stays_in_cube` is `mul_bound` again.
- `Pauli` is now generic in `R`. The `lx`/`ly`/`lz`/`comp`/`maj3` formulas are unchanged from the original. The numeric point checks (`corner`, `centre`, `edgeMid`, `diagonal_crossings`) use the same numbers, now by `norm_num`/`simp` instead of `decide +kernel`.
- Added `CharZero`/`NeZero 2` assumptions are automatic in ℚ.
- `kw_lens_polar` still has no `1 - p ≠ 0` hypothesis, and it is true at p = 1 (both sides are 0).

## Route changes
- `cancellation_hides`: his inline induction is replaced by a citation of Finance's `arith_mean_pow_le` (his own proof, moved).
- `mul_bound`: his case split is replaced by mathlib's `abs_mul` + `mul_le_of_le_one_right`. `octahedron_in_ball` is now `nlinarith`.
- `escort_escort`: his hand algebra is replaced by mathlib's `div_pow`, `List.sum_map_mul_right` and `div_div_div_cancel_right₀`.
- The seven small helpers (`sq_nn`, `pow_nn`, …) are mathlib lemmas.

New port theorems with no original: 0.
