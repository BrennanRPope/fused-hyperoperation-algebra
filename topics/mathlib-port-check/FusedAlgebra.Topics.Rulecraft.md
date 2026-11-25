# FusedAlgebra.Topics.Rulecraft

Theorems: 10. SAME 7, GENERAL 2, MATHLIB 1. No findings.

- `boundary`: `foldMA` is now `multAddList`, the same right fold since `multAddList_cons` gives `a ✳ multAddList l`. It now holds in any non-unital semiring.
- `sq_nonneg'` → `mul_self_nonneg`, which is exactly `0 ≤ d * d`.
- `balance_le`: ℤ → any linearly ordered commutative ring.
- The worked `example`s are kept. The corner-counting one is rewritten as `(2^2 - 1) ✳ (2^3 - 1) = 2^5 - 1` over ℕ, which unfolds to the original sum.

Route changes: `balance` is now the mathlib term `(sq_sub_sq s d).symm` instead of his `grind`. His case-split proof of `sq_nonneg'` is replaced by the mathlib lemma.

New port theorems with no original: 0.
