# FusedAlgebra.Topics.Fermat.Fermat4 — review

Theorems (9): SAME 5, MATHLIB 3 (sq_inj, dvd_sq_self, sq_of_coprime_mul), DROPPED-OK 1 (descent_step).

No WEAKER / DIFFERENT / MISSING / VACUOUS items, no low-confidence items. Every headline statement is word-for-word the original (same ℕ / ℚ types, same nonzero hypotheses; fermat4 still allows z = 0).

Route changes (all of them matter to the user, since this file was his own descent):
- no_x4_y4_z2 is a one-line corollary of mathlib's not_fermat_42. The original's strong induction on z with Fermat's infinite descent is not in the port.
- descent_step (original: x odd, y ≠ 0, coprime, x⁴ + y⁴ = z² ⇒ smaller nonzero solution) is dropped. That is correct as a statement (its hypotheses can't hold, by not_fermat_42), but the descent step itself is gone.
- pyth_param keeps the same statement but is now proved from mathlib's PythagoreanTriple.coprime_classification' instead of the original's own derivation through sq_of_coprime_mul.
- sq_of_coprime_mul (proved in the original without prime factorisation) is replaced by mathlib's exists_eq_pow_of_mul_eq_pow.
- fourth_power_ne_mpow_four goes through fermatLastTheoremFor_iff_rat + fermatLastTheoremFour instead of clearing denominators by hand down to fermat4.

New port theorems with no original: 0.
