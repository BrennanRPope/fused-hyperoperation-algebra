import Mathlib.RingTheory.Polynomial.Chebyshev
import FusedAlgebra.Topics.Patterns

/-!
# Our Lucas sequences are Mathlib's Chebyshev (Vieta–Fibonacci) polynomials

Repeating an element under one of our laws counts with a Lucas sequence
`U 0 = 0, U 1 = 1, U (n+2) = P U(n+1) − Q U n` (`Patterns.lucas`, and `Patterns.iter_lucas`).
When `Q = 1` this is Mathlib's rescaled Chebyshev polynomial of the second kind
`Polynomial.Chebyshev.S` (`S_n(2x) = U_n(x)`), evaluated at `P`:

* `lucas_one_eq_S`: `lucas P 1 (n + 1) = (S R n).eval P`.

So every Mathlib fact about `S` and `U` (their roots `2 cos(kπ/(n+1))`, products, compositions)
applies to our repetition counts at `Q = 1`.
-/

namespace MathlibMap

open Polynomial Polynomial.Chebyshev

variable {R : Type} [CommRing R]

theorem lucas_one_eq_S (P : R) (n : ℕ) :
    Patterns.lucas P 1 (n + 1) = (S R (n : ℤ)).eval P := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    match n, ih with
    | 0, _ => simp [Patterns.lucas, S_zero]
    | 1, _ => simp [Patterns.lucas, S_one]
    | m + 2, ih =>
      have h1 := ih (m + 1) (by omega)
      have h0 := ih m (by omega)
      rw [Patterns.lucas_succ_succ, h1, h0]
      have := S_add_two R (m : ℤ)
      push_cast
      rw [show ((m : ℤ) + 2) = (m : ℤ) + 2 from rfl, this]
      simp only [eval_sub, eval_mul, eval_X, one_mul]


end MathlibMap
