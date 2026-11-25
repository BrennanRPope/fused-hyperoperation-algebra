# FusedAlgebra.Topics.EarlyScience: port remarks

- `Ray` (a hand-rolled 2 × 2 integer matrix) is replaced by `Matrix (Fin 2) (Fin 2) R`; the
  entry `.C` is `M 1 0`. `lens` and `gap` are now matrices and the two lens theorems are about
  the `1 0` entry of a matrix product.
- `fib : ℕ → ℤ` is replaced by `Nat.fib` with a cast to `ℤ`.
- `ladder` is `ladderStep^[n] (0, 1)`; `ladder_zero` and `ladder_succ` (new) give the old
  recursion equations.
- The anonymous `example` (Kepler: 8² = 4³) is kept as an `example`.
- `constant_product` needs "no zero divisors" (it was stated on `ℤ`); everything else that was on
  `ℤ` holds in any commutative ring, and the `Lean.Grind.Field` statements in any field.
