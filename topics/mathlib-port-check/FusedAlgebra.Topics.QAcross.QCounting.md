# FusedAlgebra.Topics.QAcross.QCounting

Theorems (13): SAME 9, EQUIV 2, GENERAL 1, MATHLIB 1.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.
- `gauss_4_2`, `gauss_at_neg_one` (Fin 9 × Fin 9), `q_binomial_theorem` (Fin 7): bounds, list [1,1,2,1,1] and values 6/35/1/2 identical.
- EQUIV only because `Cube.layer` became `Nat.choose` (identical recursion) in `gauss_at_one` and `gauss_at_neg_one`.
- `uq_sl2_weight` is `mul_assoc q q x`, correctly dropped as MATHLIB.

Route changes:
- `peval_padd`, `peval_pzeros`: original direct list inductions; port derives them from `toPolynomial` and `Polynomial.eval` via the new `peval_eq_eval`.
- `ising_multAdd`: `mpow_lift` -> `mpow_eq` (same fact).

New port theorems with no original: 9 (gauss equation lemmas x3, peval_eq_eval, peval_gauss_succ_succ, ising_zero, ising_succ, transfer_pow, trace_transfer_pow).
