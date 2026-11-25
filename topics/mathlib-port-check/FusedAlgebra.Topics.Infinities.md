# FusedAlgebra.Topics.Infinities

Theorems (12): SAME 3, EQUIV 2, GENERAL 5, DROPPED-OK 2 (evalL_addL, wkey). Definitions checked: sumOnes, maddOnes, evalL, layers, wOnes (MATHLIB, they match); addL (DROPPED-OK); single (SAME).

No WEAKER / DIFFERENT / MISSING / VACUOUS items.

Medium-confidence items:
- `layers_eval`. Original: `evalL t (layers n) = (1 + t) ^ n`, where `layers` is his own Pascal recursion. Port: `∑ k ∈ range (n + 1), (n.choose k : R) * t ^ k = (1 + t) ^ n`. The port replaces `layers n` by `Nat.choose` by fiat and never proves that the recursion gives the binomial row. The row really is the binomial row, so the meaning holds, but the original's content (that this particular recursion fills (1+t)^n) is now just the binomial theorem. To settle it in Lean, prove `layers n = (List.range (n+1)).map n.choose`.
- `evalL_addL` / `addL` were dropped and pointed to `Finset.sum_add_distrib`. That lemma does not literally cover sums of padded lists of different lengths (`Polynomial.eval_add` would). This is harmless, since they were only helpers.

Route changes:
- `layers_eval`: his Pascal induction is replaced by mathlib's `add_pow`.
- `sumOnes_eq`: `simp`.
- `maddOnes_succ` and `wOnes_spec` now come from his Base lift lemmas (`one_add_mpow`, `one_add_mul_wpow`), so they are still his proofs.
- New port theorems with no original: 0 (one example on `Nat.choose 3`).
