# FusedAlgebra.Topics.QAcross.QCounting — port remarks

- Mathlib has no Gaussian binomials (only a TODO in `RingTheory/Polynomial/Pochhammer`), so
  `gauss` stays a recursive definition on coefficient lists, with equation lemmas. The
  `q`-binomial theorem stays a finite check (`n ≤ 6`) as in the original.
- `peval` is kept (computable) and identified with `Polynomial.eval` by `peval_eq_eval`.
- New: `transfer_pow` and `trace_transfer_pow` state the Ising result for the actual matrix
  `!![1, u; u, 1] ^ n`.
- `uq_sl2_weight` (`q * q * x = q * (q * x)`) is `mul_assoc` and was not restated.
