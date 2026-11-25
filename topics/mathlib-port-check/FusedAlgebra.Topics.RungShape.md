# FusedAlgebra.Topics.RungShape

Theorems: 17. SAME 13, EQUIV 2, MATHLIB 2. No WEAKER / DIFFERENT / MISSING / VACUOUS.

- `addPow_eq` → `smul_eq_mul` (MATHLIB, medium). The definition `addPow a d` (d copies of a added, on ℕ) was replaced by `d • a`, and the theorem became `d • a = d * a`. That is true, and `nsmul_eq_mul` is the more apt mathlib name. But his induction proof that repeated addition equals multiplication is gone: on ℕ the content is now part of how nsmul is defined. Worth knowing as a lost proof, not as a wrong statement. A one-line Lean check would settle the exact lemma: `example (a d : ℕ) : d • a = d * a := smul_eq_mul ..`.
- `sheet_rung0` → mathlib `add_sub_cancel` (`a + (b - a) = b` in post-2024 naming). Correct.
- `sheet_iff_line` and `sheet_diagonal`: the hypotheses `∀ u, ψ (φ u) = u` and `∀ x, φ (ψ x) = x` are now `Function.LeftInverse ψ φ` and `Function.LeftInverse φ ψ`. These unfold to exactly the same thing.

Route changes: `sheet_rung1` is now proved from `msub_multAdd` (a Basic lemma), and `sheet_edges` and `flow_closed` are one-term proofs from library lemmas. These are minor. The `addPow_eq` induction proof is lost, as above.

New port theorems with no original: 0.

Rung index note: no theorem in this module uses `MultAdd.rung`. Where the docstrings say "rung 1" or "rung s" they mean the fractional height s (s = 0 addition, s = 1 mult-add), the originals' own convention in this topic; that differs from Rungs.lean / `MultAdd.rung`, where mult-add is index 2. Same in original and port, so not a port finding.
