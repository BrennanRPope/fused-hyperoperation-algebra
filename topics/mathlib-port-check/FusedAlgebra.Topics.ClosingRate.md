# FusedAlgebra.Topics.ClosingRate

Theorems (7): SAME 5, GENERAL 2. (Definition iter: MATHLIB, checked.)

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.

- `grow`: original over ℕ `(h : ∀ j, p * d j ≤ q * d (j + 1)) : ∀ j, p ^ j * d 0 ≤ q ^ j * d j`; port in any ordered commutative semiring with added `0 ≤ p`, `0 ≤ q`, which hold automatically on ℕ, so GENERAL is right. The proof is the same calc chain.
- Route changes: the five numeric facts went from `decide` (checked by the kernel) to `norm_num`. The numbers are identical.
- New port theorems with no original: 4 (`log_two_lt`, `lt_log_two`, `two_mul_log_two_sq_lt_one`, `inv_log_four_bounds`; they come from mathlib's `Real.log_two_lt_d9` / `Real.log_two_gt_d9`).
