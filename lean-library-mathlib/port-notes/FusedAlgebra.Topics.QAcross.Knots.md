# FusedAlgebra.Topics.QAcross.Knots — port remarks

- Polynomials stay coefficient lists (`QKnots.Poly := List ℤ`). The knot invariants are kernel
  computations (`decide +kernel`, as in the original; no `native_decide`), and mathlib's
  `Polynomial` is not computable. New: `QKnots.toPolynomial : Poly → ℤ[X]` with
  `toPolynomial_padd/_pscale/_pmul/_ppow/_pzeros/_trim`, so the list operations are certified to
  be the ring operations of `ℤ[X]`. The Kauffman bracket itself is not restated in `ℤ[X]` or
  `LaurentPolynomial`.
- `state_layers` is now proved for every `n` and `k` (the original was a finite check for
  `n, k ≤ 5`), stated with `List.countP` instead of `filter … |>.length`. `length_allStates` is
  new.
- The loop counter (`relabel`, `countLoops`) is kept as written; mathlib has no computable
  connected-components count that would replace it.
- The local rules are written with powers (`A ^ 2` instead of `A * A`).
