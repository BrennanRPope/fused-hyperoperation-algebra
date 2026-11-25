# FusedAlgebra.Topics.EarlyScience: review

Theorems (17): GENERAL 12, SAME 5. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence items.

- Ray record to Matrix (Fin 2) (Fin 2): lens ⟨1,0,-P,1⟩ = !![1,0;-P,1] and gap ⟨1,d,0,1⟩ = !![1,d;0,1]. The field .C is entry (1,0). Ray.mul's C = M.C*N.A + M.D*N.C, which is (M*N) 1 0. The products keep the same order (lens P2 * lens P1, and lens P2 * (gap d * lens P1)).
- fib becomes the cast of Nat.fib (same seeds and recurrence); ladder_fib keeps the indices 2n+2 and 2n+1. ladder becomes ladderStep^[n] (0,1), the same recursion (ladder_succ).
- `constant_product` adds NoZeroDivisors, which holds automatically for ℤ.
- The Lean.Grind.Field statements are now on mathlib Field.
- The Kepler example (8² = 4³) is kept as an example.

Route changes: none of substance.

New port theorems with no original: 2 (ladder_zero, ladder_succ).
