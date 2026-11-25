# FusedAlgebra.Topics.Toric.ToricCharts

Theorems (16): SAME 11, EQUIV 2, GENERAL 3.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.
- `P1 K` became mathlib's `OnePoint K`, which mathlib defines as `Option K`; `chart1` sends x to `↑x` (= some x) in chart false and 0 to `∞` (= none), x ≠ 0 to `↑x⁻¹` in chart true, exactly as before. `trans1`, `chart2`, `trans2`, `InChart1/2` are identical, so every chart, cover, gluing and cocycle statement says what the original did.
- Injectivity is now `Function.Injective` (EQUIV).
- `corner_*` generalised from a field to a `MulZeroClass` (with `NoZeroDivisors` for `corner_scale`, automatic in a field).

Route changes: injectivity/cover proofs now use the OnePoint API (`coe_injective`, `coe_ne_infty`, `OnePoint.rec`) instead of Option case splits and grind; minor.

New port theorems with no original: 8 (corner_eq_true_iff, chart1_false, chart1_true_zero, chart1_true_of_ne_zero, equivProjectivization_chart1_false/true, trans1_one_add, card_P1_prod_P1). `trans1_one_add` states what the original only said in its header comment.
