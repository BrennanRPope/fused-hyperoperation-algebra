# FusedAlgebra.Topics.Consensus

Theorems: 27. GENERAL 16, SAME 9, DROPPED-OK 2. No findings.

- `mpow_lift` and `mpow_lift1` are dropped. They are exactly `MultAdd.one_add_mul_wpow`, and `wpow_one` composed with `one_add_mpow`, once `pw` is read as `^`.
- `depth_useless` is strengthened from weight 1 to every weight k.
- `half_is_identity` and `maj3_half` need `NeZero 2`, which is automatic in ℚ.
- `dbl` is now `(fun x ↦ x ✳ x)^[n] a`. `dbl_succ` recovers the original recursion, so `one_add_dbl` (was `dbl_lift`) means the same.
- `additive_split`: f : ℚ → ℚ became f : R → S, with n moved before the colon.
- The numeric checks (`quarter_six`, `coin_age_flip`, `fused_age_flip`, `linear_keeps`, `merge_cut`) keep the same numbers over ℚ. The BFT lemmas keep the same ℕ subtraction.

Route changes:
- `additive_split` now uses mathlib `map_nsmul`; his induction proof is gone.
- `pietrzak_halves` uses `pow_mul`/`pow_add` instead of his `pw_mul`.
- `caught` and `either` are terms from Basic lemmas instead of local proofs.
- The local induction for `mpow_lift` is dropped in favour of the Basic lemma.

New port theorems with no original: 3 (`dbl_zero`, `dbl_succ`, `dbl_eq_mpow`).
