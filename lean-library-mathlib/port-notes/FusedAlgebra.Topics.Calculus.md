# Remarks on the port of `FusedAlgebra.Topics.Calculus`

- The forward difference `diff` and its iterates `dif k` are mathlib's `fwdDiff 1` (`Δ_[1]`) and
  `(Δ_[1])^[k]`. With all of `Mathlib` imported, `Δ_[1]^[k]` does not parse (some other notation
  owns the token `]^`), so the iterate is written with parentheses, `(Δ_[1])^[k]`.
- `newton` is mathlib's Gregory–Newton formula `shift_eq_sum_fwdDiff_iter`; `pascal_sum` is kept
  as a statement although the proof of `newton` no longer needs it.
- `mint` is defined as a `Finset.prod` in `PreQuasiregular R`, so its lemmas are one-line
  consequences of `Finset.prod_range_succ`, `prod_range_add`, `prod_mul_distrib`, `map_prod`.
